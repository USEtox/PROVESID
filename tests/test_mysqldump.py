"""
Tests for ``provesid.mysqldump`` — reading ChEMBL's MySQL dump without MySQL.

The reader feeds a database that is then trusted as a reference, so most of
these tests are about refusing: a value the reader does not understand must
raise, not be dropped or guessed at.
"""

import time

import pytest

from provesid.mysqldump import (
    Column,
    CreateTable,
    DumpFormatError,
    Insert,
    parse_values,
    read_statements,
    sqlite_affinity,
    unescape,
)


def _dump(*inserts, table="t", columns=("id bigint", "name varchar(20)")):
    """A minimal dump in mysqldump's shape: one table, the given INSERT lines."""
    lines = [
        "-- MySQL dump 10.13  Distrib 8.0.36",
        "/*!40101 SET NAMES utf8mb4 */;",
        f"DROP TABLE IF EXISTS `{table}`;",
        f"CREATE TABLE `{table}` (",
    ]
    lines += [f"  `{c.split()[0]}` {c.split()[1]} DEFAULT NULL," for c in columns]
    lines += [
        "  PRIMARY KEY (`id`),",
        "  KEY `ix_name` (`name`),",
        "  CONSTRAINT `fk` FOREIGN KEY (`id`) REFERENCES `other` (`id`)",
        ") ENGINE=InnoDB DEFAULT CHARSET=utf8mb3;",
        f"LOCK TABLES `{table}` WRITE;",
        f"/*!40000 ALTER TABLE `{table}` DISABLE KEYS */;",
    ]
    lines += list(inserts)
    lines += ["UNLOCK TABLES;", "-- Dump completed on 2026-09-21"]
    return lines


def _rows(lines, tables=None):
    return [row for s in read_statements(lines, tables)
            if isinstance(s, Insert) for row in s.rows]


class TestValues:
    """What a single value turns into."""

    @pytest.mark.parametrize("escaped, meaning", [
        (r"it\'s", "it's"),
        (r"a\\b", "a\\b"),
        (r"line\nbreak", "line\nbreak"),
        (r"cr\rlf", "cr\rlf"),
        (r"tab\there", "tab\there"),
        (r"nul\0byte", "nul\0byte"),
        (r"ctrl\Zz", "ctrl\x1az"),
        (r"back\bspace", "back\bspace"),
        (r"dq\"x", 'dq"x'),
        (r"other\qchar", "otherqchar"),
    ])
    def test_every_mysql_escape(self, escaped, meaning):
        assert unescape(escaped) == meaning
        assert parse_values(f"('{escaped}');") == [(meaning,)]

    def test_an_escaped_backslash_before_a_quote_ends_the_string(self):
        assert parse_values(r"('ends with \\','next');") == [("ends with \\", "next")]

    def test_unicode_passes_through(self):
        assert parse_values("('β-carotène','日本');") == [("β-carotène", "日本")]

    def test_null_is_none_and_the_string_null_is_not(self):
        assert parse_values("(NULL,'NULL');") == [(None, "NULL")]

    @pytest.mark.parametrize("literal, value", [
        ("0", 0), ("-1", -1), ("123456789012", 123456789012),
        ("4.0", 4.0), ("180.16", 180.16), ("-0.50", -0.5), ("1e-05", 1e-05),
    ])
    def test_numbers_keep_their_kind(self, literal, value):
        parsed = parse_values(f"({literal});")[0][0]
        assert parsed == value
        assert type(parsed) is type(value)

    def test_many_rows_and_mysqldump_spacing(self):
        assert parse_values("(1,'a') , (2, 'b');") == [(1, "a"), (2, "b")]


class TestRefusals:
    """Anything not in mysqldump's vocabulary raises instead of being guessed."""

    @pytest.mark.parametrize("text", [
        "(0x4142);",                 # --hex-blob
        "(_binary 'ab');",           # a charset introducer
        "(1,'a'),(2,'b')",           # no closing ';' -- truncated
        "(1,'a';",                   # row never closed
        "(1,'unterminated);",
        "1,'a');",                   # no opening '('
        "(1,'a') (2,'b');",          # no ',' between rows
        "(1-2-3);",                  # looks numeric, is not a number
    ])
    def test_malformed_values_raise(self, text):
        with pytest.raises(DumpFormatError):
            parse_values(text)

    def test_an_unterminated_string_fails_fast(self):
        """The possessive regex must not backtrack exponentially."""
        text = "('" + "ab\\n" * 50_000
        started = time.perf_counter()
        with pytest.raises(DumpFormatError):
            parse_values(text)
        assert time.perf_counter() - started < 1.0

    def test_the_error_names_the_table(self):
        with pytest.raises(DumpFormatError, match="molecule_dictionary"):
            parse_values("(0x00);", table="molecule_dictionary")


class TestStatements:
    """Reading whole dumps."""

    def test_create_table_lists_columns_not_keys(self):
        statements = list(read_statements(_dump("INSERT INTO `t` VALUES (1,'a');")))
        assert statements[0] == CreateTable(
            "t", (Column("id", "bigint"), Column("name", "varchar(20)"))
        )

    def test_insert_rows_follow_their_create(self):
        lines = _dump("INSERT INTO `t` VALUES (1,'a'),(2,'b');",
                      "INSERT INTO `t` VALUES (3,NULL);")
        assert _rows(lines) == [(1, "a"), (2, "b"), (3, None)]

    def test_a_column_list_is_reported(self):
        lines = _dump("INSERT INTO `t` (`name`, `id`) VALUES ('a',1);")
        insert = [s for s in read_statements(lines) if isinstance(s, Insert)][0]
        assert insert.columns == ("name", "id")
        assert insert.rows == [("a", 1)]

    def test_unwanted_tables_are_not_parsed_at_all(self):
        """A value the reader would reject is harmless in a skipped table."""
        lines = _dump("INSERT INTO `t` VALUES (0xDEADBEEF);", table="t")
        assert list(read_statements(lines, tables={"wanted"})) == []

    def test_a_wanted_table_with_an_unreadable_value_raises(self):
        lines = _dump("INSERT INTO `t` VALUES (0xDEADBEEF);")
        with pytest.raises(DumpFormatError):
            list(read_statements(lines, tables={"t"}))

    def test_bytes_lines_are_decoded(self):
        lines = [line.encode("utf-8") + b"\n"
                 for line in _dump("INSERT INTO `t` VALUES (1,'β');")]
        assert _rows(lines, {"t"}) == [(1, "β")]

    def test_windows_line_endings_are_tolerated(self):
        lines = [line + "\r\n" for line in _dump("INSERT INTO `t` VALUES (1,'a');")]
        assert _rows(lines) == [(1, "a")]

    def test_rows_before_their_create_table_raise(self):
        with pytest.raises(DumpFormatError, match="before its CREATE TABLE"):
            list(read_statements(["INSERT INTO `t` VALUES (1);"]))

    def test_a_dump_ending_inside_create_table_raises(self):
        with pytest.raises(DumpFormatError, match="truncated"):
            list(read_statements(["CREATE TABLE `t` (", "  `id` bigint,"]))


class TestAffinity:
    """MySQL types map to the affinity CREATE TABLE AS SELECT would write."""

    @pytest.mark.parametrize("mysql_type, affinity", [
        ("bigint", "INT"), ("bigint(20)", "INT"), ("tinyint(1)", "INT"),
        ("smallint", "INT"), ("int", "INT"),
        ("varchar(4000)", "TEXT"), ("char(1)", "TEXT"), ("longtext", "TEXT"),
        ("decimal(9,2)", "NUM"), ("decimal(2,1)", "NUM"), ("date", "NUM"),
        ("double", "REAL"), ("float", "REAL"),
        ("longblob", ""),
    ])
    def test_affinity(self, mysql_type, affinity):
        assert sqlite_affinity(mysql_type) == affinity
