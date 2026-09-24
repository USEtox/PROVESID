"""
Connection lifetime and thread affinity for the SQLite-backed clients.

Two halves.  The first exercises :class:`provesid.sqlite_client.SQLiteClient`
against a three-row database built in ``tmp_path``, so every property of the
mixin --- closing, re-closing, ``with``, per-thread connections, the backstop
in ``__del__`` --- is tested without needing 34 GB of installed datasets.

The second runs the same contract against the four real clients and
:class:`~provesid.Search`, skipping whichever datasets are not installed.  The
regression they exist for is the one a user meets first: querying a client
from a :class:`~concurrent.futures.ThreadPoolExecutor`, which raised
``sqlite3.ProgrammingError`` on the first worker before this change.
"""

import gc
import os
import sqlite3
import threading
from concurrent.futures import ThreadPoolExecutor

import pytest

from provesid.sqlite_client import SQLiteClient, DatabaseClosedError
from provesid.utils import data_path


# ── A miniature client, so the mixin can be tested on its own ────────────────


class TinyClient(SQLiteClient):
    """Smallest thing that uses the mixin: one table, one query."""

    def __init__(self, db_path, row_factory=sqlite3.Row):
        self.db_path = str(db_path)
        self._open_database(self.db_path, row_factory=row_factory)

    def name_of(self, key):
        """Return the name for ``key`` through this thread's cursor."""
        self.cursor.execute("SELECT name FROM things WHERE id = ?", (key,))
        row = self.cursor.fetchone()
        return None if row is None else row[0]


@pytest.fixture
def tiny_db(tmp_path):
    """A three-row SQLite file on disk."""
    path = tmp_path / "tiny.db"
    conn = sqlite3.connect(path)
    conn.execute("CREATE TABLE things (id INTEGER PRIMARY KEY, name TEXT)")
    conn.executemany(
        "INSERT INTO things VALUES (?, ?)",
        [(1, "formaldehyde"), (2, "ethanol"), (3, "benzene")],
    )
    conn.commit()
    conn.close()
    return path


class TestClosing:
    """close(), the closed flag, and what a closed client does."""

    def test_query_then_close(self, tiny_db):
        client = TinyClient(tiny_db)
        assert client.closed is False
        assert client.name_of(2) == "ethanol"
        client.close()
        assert client.closed is True

    def test_close_releases_every_connection(self, tiny_db):
        client = TinyClient(tiny_db)
        client.name_of(1)
        assert len(client._open_connections) == 1
        client.close()
        assert client._open_connections == []

    def test_close_is_idempotent(self, tiny_db):
        client = TinyClient(tiny_db)
        client.close()
        client.close()
        assert client.closed is True

    def test_use_after_close_names_the_client_and_the_file(self, tiny_db):
        client = TinyClient(tiny_db)
        client.close()
        with pytest.raises(DatabaseClosedError) as excinfo:
            client.name_of(1)
        message = str(excinfo.value)
        assert "TinyClient" in message
        assert str(tiny_db) in message

    def test_use_after_close_is_a_runtime_error(self, tiny_db):
        """``except RuntimeError`` in existing code still catches it."""
        client = TinyClient(tiny_db)
        client.close()
        with pytest.raises(RuntimeError):
            client.conn

    def test_del_closes_as_a_backstop(self, tiny_db):
        """The pre-existing behaviour --- dropping the object --- still works.

        Asserted on the connection itself, since after the collector has run
        there is no client left to ask.
        """
        client = TinyClient(tiny_db)
        client.name_of(1)
        connection = client._open_connections[0]
        del client
        gc.collect()
        with pytest.raises(sqlite3.ProgrammingError):
            connection.execute("SELECT 1")

    def test_closed_is_true_before_any_database_is_opened(self):
        """An object whose constructor failed early owns nothing.

        :meth:`close` is called from ``__del__`` on such an object, so it has
        to survive the attributes never having been set.
        """
        orphan = SQLiteClient()
        assert orphan.closed is True
        assert orphan.db_file is None
        orphan.close()  # must not raise


class TestContextManager:
    """``with`` support."""

    def test_with_closes_on_exit(self, tiny_db):
        with TinyClient(tiny_db) as client:
            assert client.name_of(3) == "benzene"
        assert client.closed is True

    def test_with_returns_the_client_itself(self, tiny_db):
        client = TinyClient(tiny_db)
        with client as bound:
            assert bound is client
        assert client.closed

    def test_exception_propagates_and_still_closes(self, tiny_db):
        client = TinyClient(tiny_db)
        with pytest.raises(ValueError):
            with client:
                raise ValueError("boom")
        assert client.closed is True

    def test_reentering_a_closed_client_raises(self, tiny_db):
        client = TinyClient(tiny_db)
        client.close()
        with pytest.raises(DatabaseClosedError):
            with client:
                pass


class TestThreadAffinity:
    """One connection and one cursor per thread."""

    def test_each_thread_gets_its_own_connection(self, tiny_db):
        """Four threads, four connections, none of them the main thread's.

        The barrier is what makes this a test: without it the pool finishes
        each task before starting the next thread, and all four run on one
        worker.
        """
        workers = 4
        with TinyClient(tiny_db) as client:
            main_connection = client.conn
            barrier = threading.Barrier(workers, timeout=10)

            def record(_):
                barrier.wait()          # hold all four threads open at once
                return client.conn

            with ThreadPoolExecutor(workers) as pool:
                seen = list(pool.map(record, range(workers)))

            assert len({id(c) for c in seen}) == workers
            assert main_connection not in seen
            assert len(client._open_connections) == workers + 1

    def test_the_same_thread_reuses_its_connection(self, tiny_db):
        with TinyClient(tiny_db) as client:
            first = client.conn
            client.name_of(1)
            assert client.conn is first
            assert len(client._open_connections) == 1

    def test_the_cursor_belongs_to_this_threads_connection(self, tiny_db):
        with TinyClient(tiny_db) as client:
            assert client.cursor.connection is client.conn

    def test_querying_from_a_pool_no_longer_raises_programming_error(self, tiny_db):
        """The regression this step exists for.

        Before the mixin, a worker thread touching the single connection
        opened in ``__init__`` raised ``sqlite3.ProgrammingError: SQLite
        objects created in a thread can only be used in that same thread``.
        """
        keys = [1, 2, 3] * 20
        with TinyClient(tiny_db) as client:
            with ThreadPoolExecutor(8) as pool:
                names = list(pool.map(client.name_of, keys))
        assert names.count("formaldehyde") == 20
        assert names.count("ethanol") == 20
        assert names.count("benzene") == 20

    def test_close_reaches_the_connections_worker_threads_opened(self, tiny_db):
        workers = 4
        client = TinyClient(tiny_db)
        barrier = threading.Barrier(workers, timeout=10)

        def query(key):
            barrier.wait()
            return client.name_of(key)

        with ThreadPoolExecutor(workers) as pool:
            list(pool.map(query, [1, 2, 3, 1]))
        opened = list(client._open_connections)
        assert len(opened) == workers + 1  # one per worker, plus __init__'s

        client.close()

        for connection in opened:
            with pytest.raises(sqlite3.ProgrammingError):
                connection.execute("SELECT 1")

    def test_a_worker_thread_sees_the_close(self, tiny_db):
        """A thread that had a connection gets the package error, not sqlite3's."""
        client = TinyClient(tiny_db)
        ready, released = threading.Event(), threading.Event()
        failure = {}

        def worker():
            client.name_of(1)          # opens this thread's connection
            ready.set()
            released.wait(5)
            try:
                client.name_of(1)
            except Exception as exc:   # noqa: BLE001 - the type is the assertion
                failure["exc"] = exc

        thread = threading.Thread(target=worker)
        thread.start()
        ready.wait(5)
        client.close()
        released.set()
        thread.join(5)

        assert isinstance(failure.get("exc"), DatabaseClosedError)


class TestRowFactory:
    """The row shape each client asked for is preserved."""

    def test_row_factory_defaults_to_sqlite3_row(self, tiny_db):
        with TinyClient(tiny_db) as client:
            row = client.conn.execute("SELECT * FROM things WHERE id = 1").fetchone()
            assert row["name"] == "formaldehyde"

    def test_row_factory_none_gives_tuples(self, tiny_db):
        """ZeroPM indexes rows by position, so it opens with no row factory."""
        with TinyClient(tiny_db, row_factory=None) as client:
            row = client.conn.execute("SELECT * FROM things WHERE id = 1").fetchone()
            assert row == (1, "formaldehyde")

    def test_every_thread_gets_the_same_row_factory(self, tiny_db):
        with TinyClient(tiny_db, row_factory=None) as client:
            with ThreadPoolExecutor(2) as pool:
                rows = list(
                    pool.map(
                        lambda _: client.conn.execute(
                            "SELECT * FROM things WHERE id = 1"
                        ).fetchone(),
                        range(2),
                    )
                )
        assert rows == [(1, "formaldehyde"), (1, "formaldehyde")]


class TestAdoptedConnections:
    """A connection the client did not open itself."""

    def test_an_adopted_connection_answers_queries(self, tiny_db):
        client = object.__new__(TinyClient)
        client._adopt_connection(sqlite3.connect(tiny_db))
        try:
            assert client.name_of(2) == "ethanol"
        finally:
            client.close()

    def test_an_adopted_connection_is_closed_by_close(self, tiny_db):
        connection = sqlite3.connect(tiny_db)
        client = object.__new__(TinyClient)
        client._adopt_connection(connection)
        client.close()
        with pytest.raises(sqlite3.ProgrammingError):
            connection.execute("SELECT 1")

    def test_an_in_memory_database_works(self):
        """The case `object.__new__` plus a connection exists for."""
        connection = sqlite3.connect(":memory:")
        connection.execute("CREATE TABLE things (id INTEGER PRIMARY KEY, name TEXT)")
        connection.execute("INSERT INTO things VALUES (1, 'methane')")

        client = object.__new__(TinyClient)
        client._adopt_connection(connection, db_path=":memory:")
        with client:
            assert client.name_of(1) == "methane"
            assert client.db_file == ":memory:"

    def test_conn_has_no_setter(self, tiny_db):
        """Assigning would orphan the replaced connection; _adopt_connection does not."""
        with TinyClient(tiny_db) as client:
            with pytest.raises(AttributeError):
                client.conn = sqlite3.connect(tiny_db)


class TestWrites:
    """ZeroPM writes; per-thread connections must not break that."""

    def test_a_write_on_one_thread_is_visible_on_another(self, tiny_db):
        with TinyClient(tiny_db) as client:
            client.conn.execute("INSERT INTO things VALUES (4, 'toluene')")
            client.conn.commit()

            with ThreadPoolExecutor(1) as pool:
                seen = pool.submit(client.name_of, 4).result()

        assert seen == "toluene"


# ── The four real clients ────────────────────────────────────────────────────


def _installed(filename):
    """Whether a dataset file is present in the repository data directory."""
    return os.path.exists(os.path.join(data_path(), filename))


def _chembl_on_disk():
    """The newest ``chembl_*.db`` in the data directory, or None."""
    import glob

    candidates = sorted(glob.glob(os.path.join(data_path(), "chembl_*.db")))
    return candidates[-1] if candidates else None


def _pubchem():
    from provesid import PubChemID

    return PubChemID(auto_download=False)


def _comptox():
    from provesid import CompToxID

    return CompToxID(auto_download=False)


def _zeropm():
    from provesid import ZeroPM
    import glob

    files = sorted(glob.glob(os.path.join(data_path(), "zeropm-*.sqlite")))
    return ZeroPM(db_path=files[-1], auto_download=False)


def _chembl():
    from provesid import CheMBL

    return CheMBL(db_path=_chembl_on_disk(), auto_download=False)


#: (name, factory, a query that returns something truthy, availability)
REAL_CLIENTS = [
    pytest.param(
        "PubChemID", _pubchem, lambda db: db.get_by_cas("50-00-0")["cid"],
        marks=pytest.mark.skipif(
            not _installed("pubchem_id.db"), reason="pubchem_id.db not installed"
        ),
    ),
    pytest.param(
        "CompToxID", _comptox, lambda db: db.casrn_to_dtxsid("50-78-2"),
        marks=pytest.mark.skipif(
            not _installed("comptox_chemicals.db"),
            reason="comptox_chemicals.db not installed",
        ),
    ),
    pytest.param(
        "ZeroPM", _zeropm, lambda db: db.get_smiles_from_cas("50-00-0"),
        marks=pytest.mark.skipif(
            not _installed("zeropm-v0-0-3.sqlite")
            and not _installed("zeropm-v0-0-4.sqlite"),
            reason="no ZeroPM database installed",
        ),
    ),
    pytest.param(
        "CheMBL", _chembl, lambda db: db.search_by_chembl_id("CHEMBL25")["pref_name"],
        marks=pytest.mark.skipif(
            _chembl_on_disk() is None, reason="no ChEMBL database installed"
        ),
    ),
]


@pytest.mark.parametrize("name, factory, query", REAL_CLIENTS)
class TestTheFourClients:
    """The same contract, against the databases users actually query."""

    def test_is_a_sqlite_client(self, name, factory, query):
        with factory() as db:
            assert isinstance(db, SQLiteClient)

    def test_with_block_closes_it(self, name, factory, query):
        with factory() as db:
            assert query(db)
        assert db.closed is True

    def test_close_then_query_raises(self, name, factory, query):
        db = factory()
        db.close()
        with pytest.raises(DatabaseClosedError):
            query(db)

    def test_answers_the_same_from_eight_threads(self, name, factory, query):
        with factory() as db:
            expected = query(db)
            with ThreadPoolExecutor(8) as pool:
                answers = list(pool.map(lambda _: query(db), range(16)))
        assert answers == [expected] * 16

    def test_db_file_is_the_file_it_opened(self, name, factory, query):
        with factory() as db:
            assert db.db_file == db.db_path


# ── Search ───────────────────────────────────────────────────────────────────


NEEDS_PUBCHEM = pytest.mark.skipif(
    not _installed("pubchem_id.db"), reason="pubchem_id.db not installed"
)


@NEEDS_PUBCHEM
class TestSearchClosesWhatItOpened:
    """``Search`` owns the clients it built, and only those."""

    def test_close_closes_the_clients_it_constructed(self):
        from provesid import Search

        with Search("cas", datasets="present") as search:
            search.search("50-00-0")
            owned = {key: search._clients[key] for key in search._owned_clients}

        assert owned, "Search constructed no clients; the test proves nothing"
        for key, client in owned.items():
            if isinstance(client, SQLiteClient):
                assert client.closed is True, f"{key} was left open"

    def test_a_client_passed_in_is_left_open(self):
        """It belongs to the caller, who may still be using it."""
        from provesid import Search

        db = _pubchem()
        try:
            with Search("cas", datasets="present", pubchem=db) as search:
                search.search("50-00-0")
            assert "pubchem" not in search._owned_clients
            assert db.closed is False
            assert db.get_by_cas("50-00-0")["cid"]
        finally:
            db.close()

    def test_close_is_idempotent(self):
        from provesid import Search

        search = Search("cas", datasets="present")
        search.close()
        search.close()

    def test_searching_after_close_raises(self):
        from provesid import Search

        search = Search("cas", datasets="present")
        search.search("50-00-0")
        search.close()
        with pytest.raises(DatabaseClosedError):
            search.search("64-17-5")

    def test_close_before_any_search_is_harmless(self):
        """Clients are built lazily, so there may be nothing to close."""
        from provesid import Search

        search = Search("cas", datasets="present")
        search.close()
        assert search._owned_clients == []
