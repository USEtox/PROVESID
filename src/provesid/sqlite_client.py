"""
One connection policy for the four SQLite-backed clients.

:class:`~provesid.PubChemID`, :class:`~provesid.CompToxID`,
:class:`~provesid.ZeroPM` and :class:`~provesid.CheMBL` each open a local
database in ``__init__`` and each used to close it in ``__del__`` and nowhere
else.  That is three defects in one shape:

* **No way to let go.**  A notebook cell that re-runs ``db = PubChemID()``
  leaves the previous connection open until the collector happens to run.  On
  Windows the file stays locked while it is, so the next download cannot
  replace the database the old object is still holding.
* **No ``with``.**  Every other resource in the package is scoped; these four
  were not, so a script that opens one in a ``try`` had nothing to put in the
  ``finally``.
* **One connection, shared by every thread.**  ``sqlite3`` refuses by default
  to use a connection from a thread other than the one that created it, so a
  user who reaches for :class:`~concurrent.futures.ThreadPoolExecutor` over a
  list of CAS numbers --- the obvious thing to do with a 2.2 GB local database
  --- meets ``ProgrammingError: SQLite objects created in a thread can only be
  used in that same thread`` on the first worker.

:class:`SQLiteClient` is the mixin all four now inherit.  It gives them
:meth:`~SQLiteClient.close`, ``with`` support, and **one connection and cursor
per thread**, created on that thread's first query and closed together when
the owner is closed.

Why per-thread connections rather than ``check_same_thread=False``
------------------------------------------------------------------

``check_same_thread=False`` only removes the *check*.  It would leave every
thread sharing one connection and, worse, one ``self.cursor`` --- and these
classes are written as ``self.cursor.execute(...)`` followed by
``self.cursor.fetchone()``, so two threads interleaving those two statements
would read each other's rows.  Silently wrong answers are a far worse failure
than the exception they replace.  A connection per thread makes the existing
code correct as written, and SQLite serialises the writes that
:class:`~provesid.ZeroPM` performs when it builds an index or a view.

What threads buy, and what they do not
--------------------------------------

A pool now *works*; whether it is faster is a separate question, and for a
tight loop of nothing but local lookups the answer is no.  These queries take
tens of microseconds against a warm page cache --- less than the GIL handoff
around each one costs --- so 5 000 ``get_by_cas`` calls measured 0.29 s
serially against 14.65 s on eight threads.  That is ``sqlite3`` under CPython
rather than anything this module adds: a plain ``sqlite3.connect`` per thread
measures the same.  A pool pays when each item also does something slower ---
a request, a file read, an RDKit call: 400 lookups each followed by 20 ms of
waiting took 8.52 s serially and 1.17 s on eight threads.

Examples:
    >>> from provesid import PubChemID
    >>> with PubChemID() as db:
    ...     row = db.get_by_cas("50-00-0")
    >>> db.closed
    True

    >>> from concurrent.futures import ThreadPoolExecutor
    >>> cas_numbers = ["50-00-0", "64-17-5", "50-78-2"]
    >>> with PubChemID() as db:
    ...     with ThreadPoolExecutor(8) as pool:
    ...         rows = list(pool.map(db.get_by_cas, cas_numbers))
    >>> [row["cid"] for row in rows]
    [712, 702, 2244]
"""

import logging
import os
import sqlite3
import threading
from types import TracebackType
from typing import Any, List, Optional, Type

__all__ = ["SQLiteClient", "DatabaseClosedError"]

logger = logging.getLogger(__name__)


class DatabaseClosedError(RuntimeError):
    """Raised when a closed client is asked for its connection.

    Inherits from :class:`RuntimeError` so that ``except RuntimeError`` in
    existing calling code still catches it, and carries the class name and
    the database path so the message says which client was closed and which
    file it was reading.

    Examples:
        >>> client = SQLiteClient()
        >>> _ = client._open_database(":memory:")
        >>> client.close()
        >>> client.conn
        Traceback (most recent call last):
        ...
        provesid.sqlite_client.DatabaseClosedError: SQLiteClient was closed; its connection to :memory: is gone. Construct a new client to query it again.
    """


class SQLiteClient:
    """Connection lifetime and thread affinity for a local SQLite database.

    A mixin, not a base class with behaviour of its own: the four clients keep
    their own constructors, download logic and query methods, and call
    :meth:`_open_database` once the file they want is known to be on disk.
    After that they use :attr:`conn` and :attr:`cursor` exactly as they did
    when both were plain attributes.

    Attributes:
        conn (sqlite3.Connection): This thread's connection.  Opened on first
            access from a thread that does not have one yet.
        cursor (sqlite3.Cursor): This thread's cursor, belonging to
            :attr:`conn`.  Long-lived, as the query methods assume: they
            ``execute`` and then ``fetchone`` as two statements.
        closed (bool): True once :meth:`close` has run.

    Examples:
        >>> class Tiny(SQLiteClient):
        ...     def __init__(self, path):
        ...         self._open_database(path)
        ...     def count(self):
        ...         return self.conn.execute("SELECT COUNT(*) FROM t").fetchone()[0]
        >>> import os, sqlite3, tempfile
        >>> path = os.path.join(tempfile.mkdtemp(), "demo.db")
        >>> sqlite3.connect(path).executescript("CREATE TABLE t (x); INSERT INTO t VALUES (1), (2);")  # doctest: +ELLIPSIS
        <sqlite3.Cursor object at ...>
        >>> with Tiny(path) as tiny:
        ...     tiny.count()
        2
        >>> tiny.closed
        True
    """

    #: Seconds a connection waits for a write lock before raising
    #: ``sqlite3.OperationalError``.  The stdlib default is 5; these databases
    #: are read mostly, but ZeroPM's index build can hold a lock for longer
    #: than that on a slow disk.
    _SQLITE_TIMEOUT = 30.0

    # ── Opening ───────────────────────────────────────────────────────────────

    def _open_database(
        self,
        db_path: str,
        *,
        row_factory: Optional[Any] = sqlite3.Row,
        timeout: Optional[float] = None,
    ) -> sqlite3.Connection:
        """Take ownership of a database file and open it for this thread.

        Call once, from the subclass constructor, after the file exists.  The
        connection is opened eagerly rather than on first query so that an
        unreadable file still fails inside ``__init__``, where it failed
        before this mixin existed.

        Args:
            db_path: Path to the SQLite file.  Stored as given; the callers
                have already made it absolute.
            row_factory: Assigned to every connection this client opens.
                :class:`sqlite3.Row` (the default) gives access by column
                name; pass ``None`` for plain tuples, as :class:`ZeroPM`
                does.
            timeout: Lock timeout in seconds.  Defaults to
                :attr:`_SQLITE_TIMEOUT`.

        Returns:
            sqlite3.Connection: The calling thread's connection, also
            reachable as :attr:`conn`.

        Raises:
            sqlite3.Error: If the file cannot be opened.  Subclasses that
                promise their own exception type catch and translate it ---
                :class:`~provesid.CheMBL` raises
                :class:`~provesid.ChEMBLError`.
        """
        self._db_file = os.fspath(db_path)
        self._row_factory = row_factory
        self._sqlite_timeout = self._SQLITE_TIMEOUT if timeout is None else timeout
        self._thread_state = threading.local()
        self._open_connections: List[sqlite3.Connection] = []
        self._connection_lock = threading.RLock()
        self._closed = False
        return self._connect_this_thread()

    def _adopt_connection(
        self,
        connection: Any,
        *,
        db_path: str = "<adopted connection>",
    ) -> None:
        """Take over a connection that is already open.

        The supported way to give a client a connection it did not open
        itself: an in-memory database, or a stub in a test that builds the
        instance with ``object.__new__`` and never runs ``__init__``.  The
        connection joins the registry, so :meth:`close` closes it like any
        other.

        There is deliberately no setter on :attr:`conn`.  Assigning to it
        would leave the replaced connection open and unreachable --- nothing
        would ever close it --- and this method is the version of that
        assignment which keeps the bookkeeping straight.

        Args:
            connection: An open :class:`sqlite3.Connection`, or anything
                quacking like one.  Its ``row_factory`` is left as it is, and
                it becomes the calling thread's connection; other threads open
                their own from ``db_path``, so this is for single-threaded use
                unless ``db_path`` names a real file.
            db_path: What to report as :attr:`db_file`, and what any other
                thread would open.

        Examples:
            >>> import sqlite3
            >>> client = object.__new__(SQLiteClient)    # skips any __init__
            >>> client._adopt_connection(sqlite3.connect(":memory:"))
            >>> client.conn.execute("SELECT 1 + 1").fetchone()
            (2,)
            >>> client.db_file
            '<adopted connection>'
        """
        self._db_file = db_path
        self._row_factory = None
        self._sqlite_timeout = self._SQLITE_TIMEOUT
        self._thread_state = threading.local()
        self._open_connections = [connection]
        self._connection_lock = threading.RLock()
        self._closed = False
        self._thread_state.conn = connection
        self._thread_state.cursor = connection.cursor()

    def _connect_this_thread(self) -> sqlite3.Connection:
        """Open a connection and cursor for the calling thread.

        Returns:
            sqlite3.Connection: The new connection, already recorded in the
            registry :meth:`close` walks.

        Raises:
            DatabaseClosedError: If the client was closed first.
            sqlite3.Error: If the file cannot be opened.
        """
        self._check_open()
        # check_same_thread=False is not what makes this safe -- one
        # connection per thread is.  It is here so that close(), called from
        # whichever thread owns the object, may close the connections the
        # worker threads opened.
        connection = sqlite3.connect(
            self._db_file,
            timeout=self._sqlite_timeout,
            check_same_thread=False,
        )
        if self._row_factory is not None:
            connection.row_factory = self._row_factory

        with self._connection_lock:
            if self._closed:
                # Closed between the check above and here.  The registry has
                # already been drained, so this connection would never be
                # closed by anyone.
                connection.close()
                self._check_open()
            self._open_connections.append(connection)

        self._thread_state.conn = connection
        self._thread_state.cursor = connection.cursor()
        logger.debug(
            "Opened %s connection to %s for thread %s",
            type(self).__name__, self._db_file, threading.current_thread().name,
        )
        return connection

    def _check_open(self) -> None:
        """Raise if the client has been closed.

        Raises:
            DatabaseClosedError: If :meth:`close` has run.
        """
        if getattr(self, "_closed", False):
            raise DatabaseClosedError(
                f"{type(self).__name__} was closed; its connection to "
                f"{getattr(self, '_db_file', 'the database')} is gone. "
                "Construct a new client to query it again."
            )

    # ── Handles ───────────────────────────────────────────────────────────────

    @property
    def conn(self) -> sqlite3.Connection:
        """This thread's connection, opening one if the thread has none.

        Returns:
            sqlite3.Connection: A connection owned by the calling thread.

        Raises:
            DatabaseClosedError: If the client has been closed.

        Examples:
            >>> import os, sqlite3, tempfile
            >>> path = os.path.join(tempfile.mkdtemp(), "demo.db")
            >>> sqlite3.connect(path).executescript("CREATE TABLE t (x); INSERT INTO t VALUES (1), (2);")  # doctest: +ELLIPSIS
            <sqlite3.Cursor object at ...>
            >>> client = SQLiteClient()
            >>> _ = client._open_database(path)
            >>> client.conn.execute("SELECT COUNT(*) FROM t").fetchone()[0]
            2
        """
        self._check_open()
        connection = getattr(self._thread_state, "conn", None)
        if connection is None:
            connection = self._connect_this_thread()
        return connection

    @property
    def cursor(self) -> sqlite3.Cursor:
        """This thread's cursor, opening a connection if the thread has none.

        A cursor per thread rather than per call, because the query methods
        were written against a long-lived one: they ``execute`` in one
        statement and ``fetchone`` in the next.

        Returns:
            sqlite3.Cursor: A cursor belonging to :attr:`conn`.

        Raises:
            DatabaseClosedError: If the client has been closed.

        Examples:
            >>> import os, sqlite3, tempfile
            >>> path = os.path.join(tempfile.mkdtemp(), "demo.db")
            >>> sqlite3.connect(path).executescript("CREATE TABLE t (x); INSERT INTO t VALUES (1), (2);")  # doctest: +ELLIPSIS
            <sqlite3.Cursor object at ...>
            >>> client = SQLiteClient()
            >>> _ = client._open_database(path)
            >>> cursor = client.cursor
            >>> _ = cursor.execute("SELECT x FROM t ORDER BY x")
            >>> cursor.fetchone()["x"]
            1
        """
        self.conn  # opens this thread's connection and cursor together
        cursor = getattr(self._thread_state, "cursor", None)
        if cursor is None:
            # Only reachable if another thread closed the client between the
            # line above and this one, replacing the thread-local.  Report
            # that, rather than an AttributeError from under it.
            self._check_open()
            self._connect_this_thread()
            cursor = self._thread_state.cursor
        return cursor

    @property
    def closed(self) -> bool:
        """Whether :meth:`close` has run.

        Returns:
            bool: True after :meth:`close`, False while the client is usable.
            True as well for a client whose constructor failed before it
            reached :meth:`_open_database`, since such an object owns nothing.

        Examples:
            >>> client = SQLiteClient()
            >>> client.closed                  # never opened anything
            True
            >>> _ = client._open_database(":memory:")
            >>> client.closed
            False
        """
        return getattr(self, "_closed", True)

    @property
    def db_file(self) -> Optional[str]:
        """The database file this client opened, or None if it never opened one.

        Returns:
            str | None: The path passed to :meth:`_open_database`.

        Examples:
            >>> client = SQLiteClient()
            >>> client.db_file is None
            True
            >>> _ = client._open_database(":memory:")
            >>> client.db_file
            ':memory:'
        """
        return getattr(self, "_db_file", None)

    # ── Closing ───────────────────────────────────────────────────────────────

    def close(self) -> None:
        """Close every connection this client opened, on every thread.

        Idempotent, and safe to call on an object whose constructor raised
        before it opened anything.  After it returns, the database file is
        no longer held open by this client --- which is what lets a
        re-download replace it on Windows --- and any further use raises
        :class:`DatabaseClosedError` rather than the ``Cannot operate on a
        closed database`` that bare ``sqlite3`` would give.

        Threads other than the caller are not consulted.  A query already
        executing on another thread when this runs will fail; closing a
        client while it is being queried is a caller error, and the
        alternative --- refusing to close, or blocking until the workers
        finish --- makes ``with`` unable to guarantee anything.

        Examples:
            >>> client = SQLiteClient()
            >>> _ = client._open_database(":memory:")
            >>> client.close()
            >>> client.closed
            True
            >>> client.close()                 # idempotent
        """
        if getattr(self, "_closed", True):
            return

        with self._connection_lock:
            self._closed = True
            connections, self._open_connections = self._open_connections, []

        # A fresh threading.local drops every thread's handles at once: the
        # old object is unreferenced, so no thread can reach a closed
        # connection through it.
        self._thread_state = threading.local()

        for connection in connections:
            try:
                connection.close()
            except Exception as exc:  # pragma: no cover - sqlite3 rarely fails here
                logger.warning(
                    "Error closing %s connection to %s: %s",
                    type(self).__name__, self._db_file, exc,
                )
        logger.debug(
            "Closed %d %s connection(s) to %s",
            len(connections), type(self).__name__, self._db_file,
        )

    def __enter__(self) -> "SQLiteClient":
        """Return the client, so ``with Client() as db`` binds the client.

        Returns:
            SQLiteClient: ``self``.

        Raises:
            DatabaseClosedError: If the client has already been closed.
        """
        self._check_open()
        return self

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc_value: Optional[BaseException],
        traceback: Optional[TracebackType],
    ) -> bool:
        """Close the client on the way out of a ``with`` block.

        Args:
            exc_type: Exception class, or None.
            exc_value: Exception instance, or None.
            traceback: Traceback, or None.

        Returns:
            bool: False --- an exception raised in the block propagates.
        """
        self.close()
        return False

    def __del__(self):
        """Close as a backstop for callers who never did.

        Interpreter shutdown can have torn down enough of the module for
        :meth:`close` to fail, and an exception in ``__del__`` is printed and
        discarded rather than raised, so it is swallowed here.
        """
        try:
            self.close()
        except Exception:  # pragma: no cover - only reachable during shutdown
            pass
