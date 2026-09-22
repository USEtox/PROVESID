"""
Read the tables out of a ``mysqldump`` file without a MySQL server.

ChEMBL publishes each release three ways: a SQLite database (5.8 GB compressed,
27.7 GiB extracted), a PostgreSQL dump in ``pg_dump``'s binary custom format,
and a plain-text MySQL dump (2.1 GB compressed).  PROVESID reads eight of
ChEMBL's 74 tables, and the MySQL dump is the only one of the three that can
be read *as a stream*: one ``CREATE TABLE`` per table, then its rows as
``INSERT INTO ... VALUES (...),(...);`` lines.  Reading those lines as they are
decompressed, and keeping only the eight tables, builds the same extract as
[`provesid.CheMBL.compact`][provesid.chembl.CheMBL.compact] without the 27.7
GiB release ever existing on disk.  That is what
[`provesid.CheMBL.build_from_mysql_dump`][provesid.chembl.CheMBL.build_from_mysql_dump]
does with this module.

What is understood
------------------

Exactly what ``mysqldump`` writes, and nothing more general:

* ``CREATE TABLE `name` (`` followed by one column definition per line and a
  closing ``) ENGINE=...;`` line.  Key and constraint lines are skipped.
* ``INSERT INTO `name` VALUES (...),(...);`` on a single line, with or without
  a column list after the table name.  ``mysqldump`` never splits a statement
  across lines: a newline inside a value is written as ``\\n``.
* Values that are a quoted string, ``NULL``, or an unquoted number.  Strings
  use MySQL's backslash escapes (``\\'``, ``\\\\``, ``\\n``, ``\\r``, ``\\t``,
  ``\\0``, ``\\Z``, ``\\b``, ``\\"``); any other escaped character stands for
  itself, as MySQL defines.

Anything else inside a wanted table's ``INSERT`` --- a hex literal from
``--hex-blob``, a ``_binary`` introducer, a statement with no closing ``;`` ---
raises [`DumpFormatError`][provesid.mysqldump.DumpFormatError] naming the table
and the position.  A reader that guessed would put wrong values into a database
that is then trusted, so this one refuses instead.  Tables that are not wanted
are skipped by looking at the first few bytes of each line, which is what makes
a 2.1 GB dump cheap to read when most of it is bioactivity data.

Examples:
    >>> dump = [
    ...     "CREATE TABLE `compound` (",
    ...     "  `molregno` bigint NOT NULL,",
    ...     "  `name` varchar(255) DEFAULT NULL,",
    ...     "  PRIMARY KEY (`molregno`)",
    ...     ") ENGINE=InnoDB;",
    ...     "INSERT INTO `compound` VALUES (1,'aspirin'),(2,'it\\\\'s'),(3,NULL);",
    ... ]
    >>> for statement in read_statements(dump, tables={"compound"}):
    ...     print(statement)
    CreateTable(table='compound', columns=(Column(name='molregno', type='bigint'), Column(name='name', type='varchar(255)')))
    Insert(table='compound', columns=None, rows=[(1, 'aspirin'), (2, "it's"), (3, None)])
"""

import re
from typing import Iterable, Iterator, List, NamedTuple, Optional, Set, Tuple, Union

__all__ = [
    "Column",
    "CreateTable",
    "DumpFormatError",
    "Insert",
    "parse_values",
    "read_statements",
    "sqlite_affinity",
    "unescape",
]


class DumpFormatError(ValueError):
    """
    A line of a dump could not be read as ``mysqldump`` writes it.

    Raised rather than skipped, because a value this module does not understand
    would otherwise be dropped or mangled on its way into a database that is
    then used as a reference.  The message names the table and the character
    offset within the line.

    Examples:
        >>> list(read_statements(["INSERT INTO `compound` VALUES (1);"], tables={"compound"}))
        Traceback (most recent call last):
        ...
        provesid.mysqldump.DumpFormatError: compound: rows found before its CREATE TABLE
    """


class Column(NamedTuple):
    """
    One column of a ``CREATE TABLE``: its name, and its MySQL type as written.

    Examples:
        >>> column = Column("molregno", "bigint")
        >>> column.name, sqlite_affinity(column.type)
        ('molregno', 'INT')
    """

    name: str
    type: str


class CreateTable(NamedTuple):
    """
    A ``CREATE TABLE`` statement: the table name and its columns, in order.

    Examples:
        >>> table = CreateTable("compound", (Column("molregno", "bigint"),))
        >>> [column.name for column in table.columns]
        ['molregno']
    """

    table: str
    columns: Tuple[Column, ...]


class Insert(NamedTuple):
    """
    One ``INSERT`` statement's rows.

    ``columns`` is the explicit column list when the statement had one
    (``mysqldump --complete-insert``), and None when the rows follow the
    ``CREATE TABLE`` column order, which is ``mysqldump``'s default.

    Examples:
        >>> insert = Insert("compound", None, parse_values("(1,'aspirin'),(2,NULL);"))
        >>> insert.rows
        [(1, 'aspirin'), (2, None)]
    """

    table: str
    columns: Optional[Tuple[str, ...]]
    rows: List[tuple]


Statement = Union[CreateTable, Insert]


_CREATE_RE = re.compile(r"CREATE TABLE `([^`]+)` \(")
_COLUMN_RE = re.compile(r"\s*`([^`]+)`\s+(\S+)")
_INSERT_RE = re.compile(r"INSERT INTO `([^`]+)`\s*(?:\(([^)]*)\)\s*)?VALUES\s*")

# One value and the delimiter after it.  The string alternative consumes runs
# of ordinary characters with ``[^'\\]+`` rather than one character per
# repetition: a ChEMBL molfile is kilobytes long, and the per-character form is
# several times slower in ``re``.  The quantifiers are possessive (``++``,
# ``*+``) because a run inside a repetition is the classic shape for
# exponential backtracking, and an unterminated string --- a truncated line ---
# is exactly the input that would trigger it.
_VALUE_RE = re.compile(
    r"""\s*(?:
        '((?:[^'\\]++|\\.)*+)'      # 1: quoted string, MySQL escapes inside
      | (NULL)                      # 2: NULL
      | ([-+]?[0-9][0-9.eE+-]*)     # 3: an unquoted number
    )\s*([,)])                      # 4: the delimiter after it
    """,
    re.VERBOSE | re.DOTALL,
)
_INTEGER_RE = re.compile(r"[-+]?[0-9]+")

# MySQL's escape sequences.  An escaped character not listed stands for itself,
# which covers \' \" and \\ as well as anything unusual.
_ESCAPES = {"0": "\0", "b": "\b", "n": "\n", "r": "\r", "t": "\t", "Z": "\x1a"}
_ESCAPE_RE = re.compile(r"\\(.)", re.DOTALL)


def unescape(text: str) -> str:
    """
    Undo MySQL's backslash escaping of a string literal's contents.

    Args:
        text: What stood between the quotes, escapes still in place.

    Returns:
        The string the dump encodes.

    Examples:
        >>> unescape(r"it\\'s a\\ttab and a \\\\ backslash")
        "it's a\\ttab and a \\\\ backslash"
        >>> unescape("nothing to do")
        'nothing to do'
    """
    if "\\" not in text:
        return text
    return _ESCAPE_RE.sub(lambda match: _ESCAPES.get(match.group(1), match.group(1)),
                          text)


def _number(token: str) -> Union[int, float]:
    """An unquoted numeric literal as an ``int`` when it is one, else a ``float``."""
    if _INTEGER_RE.fullmatch(token):
        return int(token)
    return float(token)


def parse_values(text: str, pos: int = 0, *, table: str = "?") -> List[tuple]:
    """
    Parse the row tuples of an ``INSERT ... VALUES`` clause.

    Args:
        text: The statement, or just its ``VALUES`` part.
        pos: Offset of the first ``(`` in ``text``.
        table: Table name, used only in error messages.

    Returns:
        One tuple per row.  Quoted strings become ``str``, ``NULL`` becomes
        ``None``, integers ``int`` and other numbers ``float`` --- the types
        SQLite then applies its column affinity to, exactly as it does to the
        values in ChEMBL's own SQLite release.

    Raises:
        DumpFormatError: If anything other than rows of those three kinds of
            value, separated by commas and ended by ``;``, is found.

    Examples:
        >>> parse_values("(1,'a',NULL),(2,'b\\\\nc',4.50);")
        [(1, 'a', None), (2, 'b\\nc', 4.5)]
    """
    rows: List[tuple] = []
    end = len(text)
    value_match = _VALUE_RE.match

    while True:
        while pos < end and text[pos].isspace():
            pos += 1
        if pos >= end or text[pos] != "(":
            raise DumpFormatError(
                f"{table}: expected '(' at offset {pos}, found {text[pos:pos + 20]!r}"
            )
        pos += 1

        row = []
        while True:
            match = value_match(text, pos)
            if match is None:
                raise DumpFormatError(
                    f"{table}: unreadable value at offset {pos}: "
                    f"{text[pos:pos + 40]!r}"
                )
            string, null, number, delimiter = match.groups()
            if string is not None:
                row.append(unescape(string))
            elif null is not None:
                row.append(None)
            else:
                try:
                    row.append(_number(number))
                except ValueError:
                    raise DumpFormatError(
                        f"{table}: {number!r} at offset {pos} is not a number"
                    ) from None
            pos = match.end()
            if delimiter == ")":
                break
        rows.append(tuple(row))

        while pos < end and text[pos].isspace():
            pos += 1
        if pos < end and text[pos] == ",":
            pos += 1
            continue
        if pos < end and text[pos] == ";":
            return rows
        raise DumpFormatError(
            f"{table}: expected ',' or ';' after a row at offset {pos}, found "
            f"{text[pos:pos + 20]!r} --- a statement split across lines, or a "
            "truncated dump"
        )


def sqlite_affinity(mysql_type: str) -> str:
    """
    The SQLite column type a MySQL column's values should be stored under.

    Applies SQLite's own affinity rules (section 3.1 of its datatype
    documentation) to the MySQL type name, and returns the affinity keyword
    that ``CREATE TABLE ... AS SELECT`` writes.  That is the point: the extract
    [`provesid.CheMBL.compact`][provesid.chembl.CheMBL.compact] builds from
    ChEMBL's SQLite release gets its column types from exactly that statement,
    so an extract built from the MySQL dump with these types stores every value
    the same way --- a ``decimal(9,2)`` of ``180.00`` becomes the integer 180
    in both.

    Args:
        mysql_type: The type as written in the dump, e.g. ``bigint``,
            ``varchar(255)``, ``decimal(9,2)``.

    Returns:
        One of ``INT``, ``TEXT``, ``REAL``, ``NUM``, or ``""`` for a blob.

    Examples:
        >>> [sqlite_affinity(t) for t in
        ...  ("bigint", "varchar(20)", "longtext", "decimal(9,2)", "double", "blob")]
        ['INT', 'TEXT', 'TEXT', 'NUM', 'REAL', '']
    """
    name = mysql_type.upper()
    if "INT" in name:
        return "INT"
    if "CHAR" in name or "CLOB" in name or "TEXT" in name:
        return "TEXT"
    if "BLOB" in name:
        return ""
    if "REAL" in name or "FLOA" in name or "DOUB" in name:
        return "REAL"
    return "NUM"


def read_statements(
    lines: Iterable[Union[str, bytes]],
    tables: Optional[Set[str]] = None,
) -> Iterator[Statement]:
    """
    Yield the ``CREATE TABLE`` and ``INSERT`` statements of a dump, in order.

    Args:
        lines: The dump, line by line.  ``bytes`` are decoded as UTF-8, which
            is what ``mysqldump`` writes under ``SET NAMES utf8mb4``; a binary
            file object from ``tarfile`` or ``gzip`` can be passed directly.
        tables: Only statements for these tables are parsed and yielded.  None
            yields every table.  The others are skipped after a glance at the
            start of the line, which is what keeps a 2.1 GB dump with 66
            unwanted tables affordable.

    Yields:
        [`CreateTable`][provesid.mysqldump.CreateTable] once per wanted table,
        before any of its [`Insert`][provesid.mysqldump.Insert] statements.

    Raises:
        DumpFormatError: If a wanted table's rows arrive before its
            ``CREATE TABLE``, if its column definitions never close, or if an
            ``INSERT`` cannot be parsed
            ([`parse_values`][provesid.mysqldump.parse_values]).

    Examples:
        >>> dump = [
        ...     "CREATE TABLE `other` (", "  `x` int,", ") ENGINE=InnoDB;",
        ...     "INSERT INTO `other` VALUES (1);",
        ...     "CREATE TABLE `compound` (", "  `molregno` bigint NOT NULL,", ") ENGINE=InnoDB;",
        ...     "INSERT INTO `compound` VALUES (1),(2);",
        ... ]
        >>> [type(s).__name__ for s in read_statements(dump, tables={"compound"})]
        ['CreateTable', 'Insert']
    """
    creating: Optional[str] = None
    columns: List[Column] = []
    created: Set[str] = set()

    for line in lines:
        if isinstance(line, bytes):
            # Skipping on bytes first means an unwanted table's rows are never
            # decoded at all.
            if not (line.startswith(b"INSERT") or line.startswith(b"CREATE")
                    or creating is not None):
                continue
            line = line.decode("utf-8")
        line = line.rstrip("\r\n")

        if creating is not None:
            if line.startswith(")"):
                created.add(creating)
                yield CreateTable(creating, tuple(columns))
                creating = None
                continue
            match = _COLUMN_RE.match(line)
            if match:
                columns.append(Column(match.group(1), match.group(2)))
            # Anything else is a key or constraint line: not a column.
            continue

        if line.startswith("INSERT"):
            match = _INSERT_RE.match(line)
            if match is None or (tables is not None and match.group(1) not in tables):
                continue
            table = match.group(1)
            if table not in created:
                raise DumpFormatError(
                    f"{table}: rows found before its CREATE TABLE"
                )
            listed = None
            if match.group(2) is not None:
                listed = tuple(name.strip().strip("`")
                               for name in match.group(2).split(","))
            yield Insert(table, listed, parse_values(line, match.end(), table=table))
            continue

        if line.startswith("CREATE TABLE"):
            match = _CREATE_RE.match(line)
            if match and (tables is None or match.group(1) in tables):
                creating = match.group(1)
                columns = []

    if creating is not None:
        raise DumpFormatError(
            f"{creating}: the dump ends inside its CREATE TABLE --- truncated?"
        )
