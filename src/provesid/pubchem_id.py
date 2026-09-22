"""
The offline PubChem identifier database: :class:`PubChemID`.

``pubchem_id.db`` is a local SQLite file of the ~1.43 M PubChem compounds that
carry a CAS number, with their identifiers, names, synonyms, formula and
masses. :class:`PubChemID` answers lookups against it with no network, and
falls back to PUG-REST (through :class:`~provesid.pubchem.PubChemAPI`) only
for what the file does not hold, the computed descriptors among them.

The file is built from PubChem's FTP site by :mod:`provesid.pubchem_ftp`, or
downloaded prebuilt from Zenodo. The online client lives in
:mod:`provesid.pubchem`.

This module also holds :func:`rdkit_descriptors`, the RDKit computation
behind :meth:`PubChemID.descriptors`, for structures that are not in PubChem.

Examples:
    >>> from provesid import PubChemID
    >>> with PubChemID() as db:                         # doctest: +SKIP
    ...     db.cas_to_cid("50-78-2")
    2244
"""

import logging
import os
import sqlite3
from typing import Any, Dict, List, Optional, Union

import pandas as pd

from .datasets import download_file
from .pubchem import PROPERTY_CHUNK_SIZE, PubChemAPI
from .sqlite_client import SQLiteClient
from .utils import user_dataset_path


#: Descriptors :func:`rdkit_descriptors` computes, in the order it reports
#: them. The names are PubChem's wherever the quantity is the same one ---
#: a polar surface area, a count of donors --- so that a table can switch
#: source without renaming its columns. The logP is the exception: PubChem's
#: ``XLogP`` is the XLogP3 model and RDKit's is Crippen's, a different model
#: with a different number, so it keeps RDKit's name, ``MolLogP``. PubChem's
#: ``Complexity`` has no RDKit counterpart and is not here.
RDKIT_DESCRIPTORS = (
    'MolLogP',
    'TPSA',
    'HBondDonorCount',
    'HBondAcceptorCount',
    'RotatableBondCount',
    'HeavyAtomCount',
    'Charge',
)

#: The computed descriptors PubChem publishes, which ``pubchem_id.db`` no
#: longer stores (see :meth:`PubChemID.descriptors`).
PUBCHEM_DESCRIPTORS = (
    'XLogP',
    'TPSA',
    'Complexity',
    'HBondDonorCount',
    'HBondAcceptorCount',
    'RotatableBondCount',
    'HeavyAtomCount',
    'Charge',
)


def _rdkit_descriptor_functions() -> Dict[str, Any]:
    """
    Map each name in :data:`RDKIT_DESCRIPTORS` to the RDKit function computing it.

    Built on call because RDKit is slow to import, and a session that never
    asks for a descriptor should not pay for it.

    ``TPSA`` counts sulfur and phosphorus. Ertl's original definition, and
    RDKit's default, leave them out; PubChem puts them in, and counting them
    agrees with PubChem's TPSA on 70% of compounds against 60% without.

    Both floats are sums of tabulated per-atom contributions given to at most
    four decimals, so rounding to four drops the floating-point residue
    (aspirin's TPSA is 63.6, not 63.60000000000001) and nothing else.
    """
    from rdkit import Chem
    from rdkit.Chem import Crippen, rdMolDescriptors

    return {
        'MolLogP': lambda mol: round(Crippen.MolLogP(mol), 4),
        'TPSA': lambda mol: round(rdMolDescriptors.CalcTPSA(mol, includeSandP=True), 4),
        'HBondDonorCount': rdMolDescriptors.CalcNumHBD,
        'HBondAcceptorCount': rdMolDescriptors.CalcNumHBA,
        'RotatableBondCount': rdMolDescriptors.CalcNumRotatableBonds,
        'HeavyAtomCount': lambda mol: mol.GetNumHeavyAtoms(),
        'Charge': Chem.GetFormalCharge,
    }


def rdkit_descriptors(smiles: str,
                      descriptors: Optional[List[str]] = None) -> Optional[Dict[str, Any]]:
    """
    Compute molecular descriptors for a structure with RDKit.

    This is what :meth:`PubChemID.descriptors` runs on each stored SMILES,
    exposed for structures that are not in PubChem. No network; about half a
    millisecond per molecule, half of it the logP.

    The numbers are RDKit's, and they are not always PubChem's. PubChem
    computes its descriptors with Cactvs, which counts differently. Against
    PubChem's own values for 20 000 random CAS-bearing compounds, measured on
    2026-09-21:

    ========================  ==================================================
    ``HeavyAtomCount``        identical for all
    ``Charge``                identical for all
    ``HBondDonorCount``       identical for 94%
    ``RotatableBondCount``    identical for 74%: Cactvs counts, for instance,
                              the bond to a CF3 group
    ``TPSA``                  identical for 70%
    ``HBondAcceptorCount``    identical for 63%: Cactvs counts, for instance,
                              fluorine and halide counter-ions
    ``MolLogP``               Crippen's model, not XLogP3: within 0.5 of
                              PubChem's ``XLogP`` for 62%, median gap 0.37
    ========================  ==================================================

    Args:
        smiles: The structure, as SMILES.
        descriptors: Names from :data:`RDKIT_DESCRIPTORS` to compute. Defaults
            to all of them.

    Returns:
        A dict from descriptor name to value, in the order requested, or None
        when RDKit cannot parse ``smiles``. ``MolLogP`` and ``TPSA`` are
        floats, the rest ints.

    Raises:
        ValueError: If a name is not in :data:`RDKIT_DESCRIPTORS`. The message
            says where to get it instead, for PubChem's ``XLogP`` and
            ``Complexity``.

    Examples:
        >>> rdkit_descriptors("CC(=O)OC1=CC=CC=C1C(=O)O", ["TPSA", "HBondDonorCount"])
        {'TPSA': 63.6, 'HBondDonorCount': 1}
        >>> rdkit_descriptors("not a molecule") is None
        True
    """
    from rdkit import Chem, rdBase

    names = _check_descriptor_names(descriptors, 'rdkit')
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return None
    functions = _rdkit_descriptor_functions()
    return {name: functions[name](mol) for name in names}


def _check_descriptor_names(descriptors: Optional[List[str]], source: str) -> List[str]:
    """
    Validate descriptor names against what ``source`` can provide.

    Args:
        descriptors: The names asked for, or None for all of them.
        source: ``'rdkit'`` or ``'pubchem'``.

    Returns:
        The names to compute, defaulting to every descriptor ``source`` has.

    Raises:
        ValueError: If ``source`` is neither, ``descriptors`` is an empty list,
            or it names something ``source`` does not provide. A name the
            *other* source provides gets a message saying so, since that is
            the likely mistake.
    """
    if source not in ('rdkit', 'pubchem'):
        raise ValueError(f"source must be 'rdkit' or 'pubchem', got {source!r}")
    available = RDKIT_DESCRIPTORS if source == 'rdkit' else PUBCHEM_DESCRIPTORS
    if descriptors is None:
        return list(available)
    if not descriptors:
        raise ValueError("descriptors must name at least one descriptor, or be None")

    redirects = {
        ('rdkit', 'XLogP'): "XLogP is PubChem's XLogP3 model; RDKit's logP is "
                            "Crippen's, named 'MolLogP'. For PubChem's value "
                            "use source='pubchem'.",
        ('rdkit', 'Complexity'): "Complexity has no RDKit equivalent (BertzCT is "
                                 "a different number); use source='pubchem'.",
        ('pubchem', 'MolLogP'): "MolLogP is RDKit's Crippen logP; PubChem's logP "
                                "is 'XLogP'. For RDKit's value use source='rdkit'.",
    }
    for name in descriptors:
        if (source, name) in redirects:
            raise ValueError(redirects[(source, name)])
        if name not in available:
            raise ValueError(f"Unknown {source} descriptor {name!r}; "
                             f"choose from {', '.join(available)}")
    return list(descriptors)


class PubChemID(SQLiteClient):
    """
    Interface to PubChem ID SQLite database for fast identifier lookup and conversion.

    This class provides access to a local SQLite database of the ~1.43 M PubChem
    compounds that carry a CAS number, with their identifiers (CID, CAS, InChI,
    InChIKey, SMILES), names, synonyms, formula and masses.

    Where the database comes from is the ``source`` argument, and only matters
    when there is none on disk yet:

    * ``"ftp"`` (default) builds it from a dated monthly snapshot of PubChem's
      FTP site with :func:`provesid.pubchem_ftp.build_pubchem_id_db`. The
      result records its release and the MD5 of every source file (see
      :meth:`provenance`), and carries cross-references to DSSTox, ChEBI,
      ChEMBL, EC and UNII (see :meth:`xrefs`).
    * ``"zenodo"`` downloads a prebuilt copy, refreshed by hand every few
      months. Quicker to fetch, but it is whatever release it was built from.

    Both hold the same tables, so every lookup works on either. They differ in
    the property columns: a database built from FTP has ``MonoisotopicMass``
    and none of the eight computed descriptors (XLogP, TPSA and the like).
    :meth:`descriptors` computes those with RDKit from the stored SMILES, or
    fetches PubChem's own on request; :meth:`properties` fetches PubChem's.
    See :attr:`offline_properties` for what the open database can answer.

    Connection handling comes from
    :class:`~provesid.sqlite_client.SQLiteClient`: use the class as a context
    manager, or call :meth:`~provesid.sqlite_client.SQLiteClient.close` when
    finished, and query it from as many threads as you like --- each gets its
    own connection.

    Attributes:
        db_path (str): Path to the SQLite database file
        conn (sqlite3.Connection): This thread's database connection
        source (str): The acquisition route this instance was given.
        offline_properties (dict): The part of :attr:`OFFLINE_PROPERTIES` the
            open database has columns for --- what :meth:`properties` answers
            without the network, and retrieves when no properties are named.

    Examples:
        >>> from provesid import PubChemID
        >>> with PubChemID() as db:
        ...     db.get_by_cas("50-78-2")["cmpdname"]
        ...     db.cas_to_inchikey("50-78-2")
        ...     db.inchikey_to_cid("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
        ...     db.batch_cas_to_cid(["50-78-2", "50-00-0"])
        'Aspirin'
        'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
        2244
        {'50-78-2': 2244, '50-00-0': 712}
    """

    DEFAULT_DB_NAME = "pubchem_id.db"
    DEFAULT_DB_URL = "https://zenodo.org/records/18173204/files/pubchem_id.db"

    #: Where a missing database comes from. ``"ftp"`` builds it from PubChem's
    #: FTP site (:mod:`provesid.pubchem_ftp`); ``"zenodo"`` downloads a
    #: prebuilt copy.
    SOURCES = ("ftp", "zenodo")

    #: PubChem property names the local database can answer, mapped to their
    #: column in the ``compounds`` table. These are the properties that are
    #: *data* about a compound --- its identifiers, names, formula and masses.
    #: The computed descriptors (``XLogP``, ``TPSA``, ``Complexity`` and the
    #: counts) are not served from disk even by a Zenodo database that still
    #: holds them: they are PubChem's model outputs, and a user who asks for
    #: them gets PubChem's current values, labelled ``Source='online'``, or
    #: RDKit's from :meth:`descriptors`, labelled ``Source='rdkit'``.
    #: Note that ``smiles`` holds the isomeric SMILES, which is what PubChem now
    #: calls ``SMILES``; the stereochemistry-free ``ConnectivitySMILES`` is not
    #: stored locally. ``MolecularWeight`` is computed from the formula when the
    #: database is built from FTP --- PubChem's files do not carry it --- and
    #: agrees with PubChem's to the second decimal for most compounds.
    OFFLINE_PROPERTIES = {
        'MolecularFormula': 'mf',
        'MolecularWeight': 'mw',
        'SMILES': 'smiles',
        'InChI': 'inchi',
        'InChIKey': 'inchikey',
        'IUPACName': 'iupacname',
        'Title': 'cmpdname',
        'ExactMass': 'exactmass',
        'MonoisotopicMass': 'monoisotopicmass',
    }

    #: What :meth:`properties` retrieves when the caller names no properties,
    #: for a database that has every column. An open database uses
    #: :attr:`offline_properties` instead, so that a Zenodo copy without
    #: ``monoisotopicmass`` does not send every default lookup online.
    DEFAULT_PROPERTIES = tuple(OFFLINE_PROPERTIES)

    #: Type each property is normalised to, so that a table assembled from both
    #: sources is usable as one table. PUG-REST reports ``MolecularWeight`` and
    #: ``ExactMass`` as strings while the local database holds floats.
    _PROPERTY_CASTS = {
        'MolecularWeight': float,
        'ExactMass': float,
        'MonoisotopicMass': float,
        'XLogP': float,
        'TPSA': float,
        'Complexity': float,
        'Charge': int,
        'HBondDonorCount': int,
        'HBondAcceptorCount': int,
        'RotatableBondCount': int,
        'HeavyAtomCount': int,
    }

    #: Bound parameters per ``IN`` clause. SQLite's own default ceiling is 999.
    _SQL_PARAMETER_LIMIT = 500


    def __init__(
        self,
        db_path: Optional[str] = None,
        auto_download: bool = True,
        data_dir: Optional[str] = None,
        db_url: Optional[str] = None,
        redownload: bool = False,
        api: Optional['PubChemAPI'] = None,
        source: str = "ftp",
    ):
        """
        Initialize PubChemID database connection.

        Args:
            db_path (str, optional): Path to SQLite database. If None, uses default
                                    location in the persistent user dataset directory.
            auto_download (bool): If True, acquire the database from ``source``
                when it is not on disk. Default is True.
            data_dir (str, optional): Directory to store the database when
                ``db_path`` is not provided.
            db_url (str, optional): Download URL for ``source="zenodo"``. If
                None, uses the package default URL.
            redownload (bool): If True, acquire the database again even though
                one is on disk, when ``auto_download`` is enabled. With
                ``source="ftp"`` that is a rebuild from the newest snapshot.
            api (PubChemAPI, optional): Online client used by
                :meth:`properties` when the local database cannot answer a
                request. One is created on first use if none is given, so
                passing this is only needed to share a client or to configure
                its pause time.
            source (str): How a missing database is acquired, one of
                :attr:`SOURCES`. ``"ftp"`` (default) builds it from the newest
                monthly snapshot of PubChem's FTP site: 15.4 GB transferred,
                a 2.5 GB database plus 7.4 GB of free disk at peak, and about
                12 minutes of processing on top of the download time.
                ``"zenodo"`` downloads a 2.2 GB
                prebuilt copy. It describes an acquisition, not a file: a
                database already on disk is opened whichever way it was made.

        Raises:
            ValueError: If ``source`` is not one of :attr:`SOURCES`. Checked
                before anything is fetched.
            FileNotFoundError: If database file doesn't exist and auto_download is False

        Examples:
            >>> db = PubChemID()                       # the default location
            >>> db.source
            'ftp'
            >>> PubChemID(db_path="/no/such/pubchem_id.db", auto_download=False)
            Traceback (most recent call last):
            ...
            FileNotFoundError: PubChem ID database not found at /no/such/pubchem_id.db. ...
            >>> PubChemID(source="ncbi")
            Traceback (most recent call last):
            ...
            ValueError: PubChemID(source='ncbi') is not a download route. Use one of 'ftp', 'zenodo'.
        """
        self.logger = logging.getLogger(__name__)
        self.source = self._validate_source(source)
        self.db_url = db_url or self.DEFAULT_DB_URL
        self._api = api

        if db_path is None:
            base_dir = data_dir or user_dataset_path()
            db_path = os.path.join(base_dir, self.DEFAULT_DB_NAME)

        self.db_path = os.path.abspath(os.path.expanduser(db_path))

        needs_download = redownload or not os.path.exists(self.db_path)

        if needs_download:
            if auto_download:
                if redownload and os.path.exists(self.db_path):
                    self.logger.info(
                        "Forced PubChemID redownload requested for: %s", self.db_path
                    )
                else:
                    self.logger.info("Database not found at %s", self.db_path)
                self._acquire(force=redownload)
            else:
                raise FileNotFoundError(
                    f"PubChem ID database not found at {self.db_path}. "
                    "Set auto_download=True, run "
                    "provesid.pubchem_ftp.build_pubchem_id_db(), or run "
                    "PubChemID.download_database()."
                )

        # One connection per thread, released by close() or by leaving a
        # ``with`` block --- see
        # :class:`~provesid.sqlite_client.SQLiteClient`.
        self._open_database(self.db_path)
        self.offline_properties = self._available_offline_properties()

    @classmethod
    def _validate_source(cls, source: str) -> str:
        """
        Check an acquisition route name before anything is fetched.

        Args:
            source: The ``source`` argument as given.

        Returns:
            The same value, once it is known to be a route.

        Raises:
            ValueError: If it is not one of :attr:`SOURCES`.
        """
        if source in cls.SOURCES:
            return source
        options = ", ".join(repr(name) for name in cls.SOURCES)
        raise ValueError(f"PubChemID(source={source!r}) is not a download route. "
                         f"Use one of {options}.")

    def _acquire(self, *, force: bool) -> None:
        """Put a database at :attr:`db_path` by the route :attr:`source` names."""
        if self.source == "ftp":
            from .pubchem_ftp import build_pubchem_id_db

            self.logger.info("Building the PubChem ID database from PubChem's FTP site")
            build_pubchem_id_db(self.db_path, force=force)
        else:
            self.logger.info("Downloading the PubChem ID database from Zenodo")
            self.download_database(db_path=self.db_path, zenodo_url=self.db_url,
                                   force=force)

    def _available_offline_properties(self) -> Dict[str, str]:
        """
        The part of :attr:`OFFLINE_PROPERTIES` this database has columns for.

        A database built from FTP and one downloaded from Zenodo differ in one
        column --- ``monoisotopicmass`` exists only in the first --- so which
        properties can be answered from disk is a fact about the file, not the
        class.
        """
        columns = {row[1] for row in self.conn.execute("PRAGMA table_info(compounds)")}
        return {name: column for name, column in self.OFFLINE_PROPERTIES.items()
                if column in columns}

    @staticmethod
    def download_database(
        db_path: Optional[str] = None,
        zenodo_url: Optional[str] = None,
        force: bool = False,
    ) -> str:
        """
        Download PubChem ID database from Zenodo --- the ``source="zenodo"`` route.

        To build it from PubChem's own files instead, which is the default
        route, see :func:`provesid.pubchem_ftp.build_pubchem_id_db`.

        The transfer is resumable: an interrupted download leaves a ``.part``
        file beside the destination and the next call continues from it rather
        than fetching the 2.2 GB again. The file is opened and queried before
        it is moved into place, so a damaged download never replaces a working
        database.

        Args:
            db_path (str, optional): Path where to save the database. If None, uses default
                                    location in the persistent user dataset directory.
            zenodo_url (str, optional): URL to download from. If None, uses default Zenodo URL.
                                       Format: https://zenodo.org/record/XXXXXX/files/pubchem_id.db
            force (bool): If True, overwrite an existing local database file.

        Returns:
            str: Path to the downloaded database file

        Raises:
            FileExistsError: If the database exists and ``force`` is False.
            provesid.datasets.DownloadError: If the download could not be
                completed.
            RuntimeError: If the file that arrived is not the PubChem ID
                database.

        Examples:
            >>> from provesid import PubChemID
            >>> PubChemID.download_database(force=True)                  # doctest: +SKIP
            '/home/me/.local/share/provesid/pubchem_id.db'
            >>> PubChemID.download_database(db_path='/tmp/pubchem_id.db')  # doctest: +SKIP
            '/tmp/pubchem_id.db'

        Note:
            The database file is ~2.2 GB, so download may take several minutes.
        """
        logger = logging.getLogger(__name__)

        if db_path is None:
            db_path = os.path.join(
                user_dataset_path(),
                PubChemID.DEFAULT_DB_NAME,
            )
        else:
            db_path = os.path.abspath(os.path.expanduser(db_path))

        if os.path.exists(db_path) and not force:
            raise FileExistsError(
                f"Database already exists at: {db_path}. Use force=True to overwrite."
            )

        def must_be_the_compounds_database(path):
            """Reject a download that cannot answer the query this class asks."""
            connection = sqlite3.connect(path)
            try:
                count = connection.execute(
                    "SELECT COUNT(*) FROM compounds"
                ).fetchone()[0]
            except Exception as exc:
                raise RuntimeError(
                    f"Downloaded file is not a valid database: {exc}"
                ) from exc
            finally:
                connection.close()
            logger.info("Database verified: %s compounds", f"{count:,}")

        logger.info("This is a large file (~2.2 GB), please be patient.")
        return download_file(
            zenodo_url or PubChemID.DEFAULT_DB_URL,
            db_path,
            verify=must_be_the_compounds_database,
            description="PubChem ID database",
            log=logger,
        )

    def get_by_cid(self, cid: int) -> Optional[Dict[str, Any]]:
        """
        Get compound information by PubChem CID.

        Every other ``get_by_*`` method finds a CID and then returns this
        record for it.

        Args:
            cid (int): PubChem Compound ID

        Returns:
            dict: Every column of the ``compounds`` row (``cid``, ``cmpdname``,
            ``mf``, ``inchi``, ``smiles``, ``inchikey``, ``iupacname``, ``mw``,
            ``exactmass``, ``cidcdate`` and whichever others the database has;
            see :meth:`get_by_cas_batch`), plus ``cas_numbers``, the compound's
            distinct CAS numbers in database order, and ``synonyms``, at most
            100 of its names. None if the CID is not in the database.

        Examples:
            >>> db = PubChemID()
            >>> result = db.get_by_cid(2244)  # Aspirin
            >>> result['cmpdname'], result['mf'], result['cas_numbers']
            ('Aspirin', 'C9H8O4', ['50-78-2'])
            >>> result['synonyms'][:2]
            ['aspirin', 'ACETYLSALICYLIC ACID']
            >>> db.get_by_cid(999999999) is None
            True
        """
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT * FROM compounds WHERE cid = ?
        """, (cid,))

        row = cursor.fetchone()
        if not row:
            return None

        result = dict(row)

        # Add CAS numbers
        # GROUP BY because the Zenodo copy repeats some (cid, cas) pairs ---
        # aspirin's CAS is listed twice. A database built from FTP does not.
        cursor.execute("""
            SELECT cas FROM cas_numbers WHERE cid = ?
            GROUP BY cas ORDER BY MIN(id)
        """, (cid,))
        result['cas_numbers'] = [r[0] for r in cursor.fetchall()]

        # Add synonyms
        cursor.execute("""
            SELECT synonym FROM synonyms WHERE cid = ? LIMIT 100
        """, (cid,))
        result['synonyms'] = [r[0] for r in cursor.fetchall()]

        return result

    def get_by_cas(self, cas: str) -> Optional[Dict[str, Any]]:
        """
        Get compound information by CAS Registry Number.

        Args:
            cas (str): CAS Registry Number (e.g., "50-78-2")

        Returns:
            dict: The :meth:`get_by_cid` record, or None if not found. A CAS
            number PubChem gives to several compounds returns the first one.

        Examples:
            >>> db = PubChemID()
            >>> result = db.get_by_cas("50-78-2")  # Aspirin
            >>> print(result['inchi'])
            InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)
            >>> db.get_by_cas("50782") is None       # the hyphens are required
            True
        """
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT cid FROM cas_numbers WHERE cas = ? LIMIT 1
        """, (cas,))

        row = cursor.fetchone()
        if not row:
            return None

        return self.get_by_cid(row[0])

    def get_by_inchikey(self, inchikey: str) -> Optional[Dict[str, Any]]:
        """
        Get compound information by InChIKey.

        Args:
            inchikey (str): Standard InChIKey (27 characters)

        Returns:
            dict: The :meth:`get_by_cid` record, or None if not found.

        Examples:
            >>> db = PubChemID()
            >>> result = db.get_by_inchikey("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
            >>> print(result['cmpdname'])
            Aspirin
        """
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT * FROM compounds WHERE inchikey = ?
        """, (inchikey,))

        row = cursor.fetchone()
        if not row:
            return None

        cid = row['cid']
        return self.get_by_cid(cid)

    def get_by_inchi(self, inchi: str) -> Optional[Dict[str, Any]]:
        """
        Get compound information by InChI string.

        Args:
            inchi (str): Standard InChI string

        Returns:
            dict: The :meth:`get_by_cid` record, or None if not found. The
            match is exact: a truncated or non-standard InChI finds nothing.

        Examples:
            >>> db = PubChemID()
            >>> inchi = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
            >>> print(db.get_by_inchi(inchi)['cmpdname'])
            Aspirin
            >>> db.get_by_inchi("InChI=1S/C9H8O4/c1-6(10)") is None
            True
        """
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT * FROM compounds WHERE inchi = ?
        """, (inchi,))

        row = cursor.fetchone()
        if not row:
            return None

        cid = row['cid']
        return self.get_by_cid(cid)

    def get_by_smiles(self, smiles: str) -> Optional[Dict[str, Any]]:
        """
        Get compound information by SMILES string.

        The match is on the stored string, not the structure, so only
        PubChem's own SMILES for a compound finds it. :meth:`smiles_to_cas`
        compares structures instead.

        Args:
            smiles (str): SMILES string

        Returns:
            dict: The :meth:`get_by_cid` record, or None if not found.

        Examples:
            >>> db = PubChemID()
            >>> result = db.get_by_smiles("CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin
            >>> print(result['cmpdname'])
            Aspirin
            >>> db.get_by_smiles("CCO")['cid'], db.get_by_smiles("OCC")
            (702, None)
        """
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT * FROM compounds WHERE smiles = ?
        """, (smiles,))

        row = cursor.fetchone()
        if not row:
            return None

        cid = row['cid']
        return self.get_by_cid(cid)

    def search_by_name(self, name: str, exact: bool = False, limit: int = 10) -> List[Dict[str, Any]]:
        """
        Search compounds by name or synonym.

        Compound titles are searched first, then synonyms, until ``limit`` is
        reached. An exact match is case-sensitive: ``"Aspirin"`` is the title
        and ``"aspirin"`` a synonym, and both find CID 2244. A partial match
        is SQL ``LIKE``, case-insensitive for ASCII letters, and returns
        compounds in database order, not by closeness.

        Args:
            name (str): Compound name or synonym to search for
            exact (bool): If True, exact match only. If False, partial match (case-insensitive)
            limit (int): Maximum number of results to return

        Returns:
            list: :meth:`get_by_cid` records, each compound once. Empty when
            nothing matches.

        Examples:
            >>> db = PubChemID()
            >>> for r in db.search_by_name("aspirin", limit=3):
            ...     print(r['cid'], r['cmpdname'])
            2244 Aspirin
            6247 Calcium aspirin
            21975 Carbaspirin Calcium
            >>> [r['cid'] for r in db.search_by_name("aspirin", exact=True)]
            [2244]
        """
        cursor = self.conn.cursor()

        results = []

        if exact:
            # Search in main compound name
            cursor.execute("""
                SELECT cid FROM compounds WHERE cmpdname = ? LIMIT ?
            """, (name, limit))

            cids = [r[0] for r in cursor.fetchall()]

            # Also search in synonyms
            if len(cids) < limit:
                cursor.execute("""
                    SELECT DISTINCT cid FROM synonyms WHERE synonym = ? LIMIT ?
                """, (name, limit - len(cids)))
                cids.extend([r[0] for r in cursor.fetchall()])
        else:
            # Partial match with LIKE
            search_term = f"%{name}%"

            # Search in main compound name
            cursor.execute("""
                SELECT cid FROM compounds WHERE cmpdname LIKE ? LIMIT ?
            """, (search_term, limit))

            cids = [r[0] for r in cursor.fetchall()]

            # Also search in synonyms
            if len(cids) < limit:
                cursor.execute("""
                    SELECT DISTINCT cid FROM synonyms WHERE synonym LIKE ? LIMIT ?
                """, (search_term, limit - len(cids)))
                cids.extend([r[0] for r in cursor.fetchall()])

        # A compound can match both its title and a synonym.
        cids = list(dict.fromkeys(cids))

        # Get full compound info for each CID
        for cid in cids[:limit]:
            compound = self.get_by_cid(cid)
            if compound:
                results.append(compound)

        return results

    def search_by_formula(self, formula: str, limit: int = 100) -> List[Dict[str, Any]]:
        """
        Search compounds by molecular formula.

        Args:
            formula (str): Molecular formula (e.g., "C9H8O4")
            limit (int): Maximum number of results to return

        Returns:
            list: :meth:`get_by_cid` records for compounds whose formula is
            exactly ``formula``, in database order. The formula must be
            written as PubChem writes it (Hill order).

        Examples:
            >>> db = PubChemID()
            >>> results = db.search_by_formula("C9H8O4", limit=5)
            >>> len(results), all(r['mf'] == 'C9H8O4' for r in results)
            (5, True)
            >>> 'Aspirin' in [r['cmpdname'] for r in db.search_by_formula("C9H8O4")]
            True
        """
        cursor = self.conn.cursor()
        cursor.execute("""
            SELECT cid FROM compounds WHERE mf = ? LIMIT ?
        """, (formula, limit))

        results = []
        for row in cursor.fetchall():
            compound = self.get_by_cid(row[0])
            if compound:
                results.append(compound)

        return results

    # Conversion methods

    def cas_to_cid(self, cas: str) -> Optional[int]:
        """
        Convert CAS number to PubChem CID.

        Args:
            cas: CAS Registry Number, with hyphens.

        Returns:
            The CID, or None if the CAS number is not in the database.

        Examples:
            >>> db = PubChemID()
            >>> db.cas_to_cid("50-78-2")
            2244
        """
        result = self.get_by_cas(cas)
        return result['cid'] if result else None

    def cas_to_inchi(self, cas: str) -> Optional[str]:
        """
        Convert CAS number to InChI.

        Args:
            cas: CAS Registry Number, with hyphens.

        Returns:
            The standard InChI, or None if not found.

        Examples:
            >>> db = PubChemID()
            >>> db.cas_to_inchi("50-78-2")
            'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'
        """
        result = self.get_by_cas(cas)
        return result['inchi'] if result else None

    def cas_to_inchikey(self, cas: str) -> Optional[str]:
        """
        Convert CAS number to InChIKey.

        Args:
            cas: CAS Registry Number, with hyphens.

        Returns:
            The standard InChIKey, or None if not found.

        Examples:
            >>> db = PubChemID()
            >>> db.cas_to_inchikey("50-78-2")
            'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
        """
        result = self.get_by_cas(cas)
        return result['inchikey'] if result else None

    def cas_to_smiles(self, cas: str) -> Optional[str]:
        """
        Convert CAS number to SMILES.

        Args:
            cas: CAS Registry Number, with hyphens.

        Returns:
            PubChem's isomeric SMILES, or None if not found.

        Examples:
            >>> db = PubChemID()
            >>> db.cas_to_smiles("50-78-2")
            'CC(=O)OC1=CC=CC=C1C(=O)O'
        """
        result = self.get_by_cas(cas)
        return result['smiles'] if result else None

    def inchikey_to_cid(self, inchikey: str) -> Optional[int]:
        """
        Convert InChIKey to PubChem CID.

        Args:
            inchikey: Standard InChIKey.

        Returns:
            The CID, or None if not found.

        Examples:
            >>> db = PubChemID()
            >>> db.inchikey_to_cid("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
            2244
        """
        result = self.get_by_inchikey(inchikey)
        return result['cid'] if result else None

    def inchikey_to_cas(self, inchikey: str) -> Optional[List[str]]:
        """
        Convert InChIKey to CAS number(s).

        Args:
            inchikey: Standard InChIKey.

        Returns:
            The compound's CAS numbers (possibly empty), or None if the InChIKey is not found.

        Examples:
            >>> db = PubChemID()
            >>> db.inchikey_to_cas("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")
            ['50-78-2']
        """
        result = self.get_by_inchikey(inchikey)
        return result['cas_numbers'] if result else None

    def inchi_to_cid(self, inchi: str) -> Optional[int]:
        """
        Convert InChI to PubChem CID.

        Args:
            inchi: Standard InChI, matched exactly.

        Returns:
            The CID, or None if not found.

        Examples:
            >>> db = PubChemID()
            >>> db.inchi_to_cid("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
            702
        """
        result = self.get_by_inchi(inchi)
        return result['cid'] if result else None

    def inchi_to_cas(self, inchi: str) -> Optional[List[str]]:
        """
        Convert InChI to CAS number(s).

        Args:
            inchi: Standard InChI, matched exactly.

        Returns:
            The compound's CAS numbers, or None if the InChI is not found.

        Examples:
            >>> db = PubChemID()
            >>> db.inchi_to_cas("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
            ['64-17-5']
        """
        result = self.get_by_inchi(inchi)
        return result['cas_numbers'] if result else None

    def cid_to_cas(self, cid: int) -> Optional[List[str]]:
        """
        Convert PubChem CID to CAS number(s).

        Args:
            cid: PubChem Compound ID.

        Returns:
            The compound's distinct CAS numbers, or None if the CID is not found. A compound can have several: retired numbers, and numbers for mixtures PubChem maps to it.

        Examples:
            >>> db = PubChemID()
            >>> db.cid_to_cas(712)
            ['50-00-0', '30525-89-4', '53026-80-5', '8013-13-6', '12795-06-1']
        """
        result = self.get_by_cid(cid)
        return result['cas_numbers'] if result else None

    def cid_to_inchikey(self, cid: int) -> Optional[str]:
        """
        Convert PubChem CID to InChIKey.

        Args:
            cid: PubChem Compound ID.

        Returns:
            The standard InChIKey, or None if not found.

        Examples:
            >>> db = PubChemID()
            >>> db.cid_to_inchikey(2244)
            'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
        """
        result = self.get_by_cid(cid)
        return result['inchikey'] if result else None

    def cid_to_inchi(self, cid: int) -> Optional[str]:
        """
        Convert PubChem CID to InChI.

        Args:
            cid: PubChem Compound ID.

        Returns:
            The standard InChI, or None if not found.

        Examples:
            >>> db = PubChemID()
            >>> db.cid_to_inchi(702)
            'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
        """
        result = self.get_by_cid(cid)
        return result['inchi'] if result else None

    def cid_to_smiles(self, cid: int) -> Optional[str]:
        """
        Convert PubChem CID to SMILES.

        Args:
            cid: PubChem Compound ID.

        Returns:
            PubChem's isomeric SMILES, or None if not found.

        Examples:
            >>> db = PubChemID()
            >>> db.cid_to_smiles(2244)
            'CC(=O)OC1=CC=CC=C1C(=O)O'
        """
        result = self.get_by_cid(cid)
        return result['smiles'] if result else None

    def smiles_to_cid(self, smiles: str) -> Optional[int]:
        """
        Convert SMILES string to PubChem CID.

        Args:
            smiles: SMILES, matched as a string against PubChem's; see :meth:`get_by_smiles`.

        Returns:
            The CID, or None if not found.

        Examples:
            >>> db = PubChemID()
            >>> db.smiles_to_cid("CCO")
            702
        """
        result = self.get_by_smiles(smiles)
        return result['cid'] if result else None

    # Batch conversion methods

    def batch_cas_to_cid(self, cas_list: List[str]) -> Dict[str, Optional[int]]:
        """
        Convert multiple CAS numbers to CIDs.

        Args:
            cas_list (list): List of CAS numbers

        Returns:
            dict: Mapping of CAS -> CID (None if not found), in input order.

        Examples:
            >>> db = PubChemID()
            >>> results = db.batch_cas_to_cid(["50-78-2", "50-00-0"])
            >>> print(results)
            {'50-78-2': 2244, '50-00-0': 712}
        """
        results = {}
        for cas in cas_list:
            results[cas] = self.cas_to_cid(cas)
        return results

    def batch_cas_to_inchikey(self, cas_list: List[str]) -> Dict[str, Optional[str]]:
        """
        Convert multiple CAS numbers to InChIKeys.

        Args:
            cas_list (list): List of CAS numbers

        Returns:
            dict: Mapping of CAS -> InChIKey (None if not found)

        Examples:
            >>> db = PubChemID()
            >>> db.batch_cas_to_inchikey(["50-78-2", "0-00-0"])
            {'50-78-2': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N', '0-00-0': None}
        """
        results = {}
        for cas in cas_list:
            results[cas] = self.cas_to_inchikey(cas)
        return results

    def batch_cid_to_cas(self, cid_list: List[int]) -> Dict[int, Optional[List[str]]]:
        """
        Convert multiple CIDs to CAS numbers.

        Args:
            cid_list (list): List of PubChem CIDs

        Returns:
            dict: Mapping of CID -> list of CAS numbers (None if not found)

        Examples:
            >>> db = PubChemID()
            >>> db.batch_cid_to_cas([2244, 702])
            {2244: ['50-78-2'], 702: ['64-17-5']}
        """
        results = {}
        for cid in cid_list:
            results[cid] = self.cid_to_cas(cid)
        return results

    def batch_smiles_to_cid(self, smiles_list: List[str]) -> Dict[str, Optional[int]]:
        """
        Convert multiple SMILES strings to CIDs.

        Args:
            smiles_list (list): List of SMILES strings

        Returns:
            dict: Mapping of SMILES -> CID (None if not found)

        Examples:
            >>> db = PubChemID()
            >>> results = db.batch_smiles_to_cid(["CC(=O)OC1=CC=CC=C1C(=O)O", "C"])
            >>> print(results)
            {'CC(=O)OC1=CC=CC=C1C(=O)O': 2244, 'C': 297}
        """
        results = {}
        for smiles in smiles_list:
            results[smiles] = self.smiles_to_cid(smiles)
        return results

    def get_by_cas_batch(self, cas_list: List[str]) -> 'pd.DataFrame':
        """
        Get complete compound information for multiple CAS numbers as a DataFrame.

        One row per CAS number found, carrying every column of the
        ``compounds`` table. Which columns those are depends on how the
        database was made: one built from PubChem's FTP site has
        ``monoisotopicmass``, a Zenodo copy has the eight descriptor columns
        instead (``xlogp``, ``polararea`` and the like).

        Args:
            cas_list (list): List of CAS Registry Numbers

        Returns:
            pandas.DataFrame: ``cid``, ``cas`` and then the ``compounds``
            columns --- ``cmpdname``, ``mf``, ``inchi``, ``smiles``,
            ``inchikey``, ``iupacname``, ``mw``, ``exactmass``, ``cidcdate``
            and whichever others the database has. Empty, with those columns,
            when nothing is found.

        Examples:
            >>> db = PubChemID()
            >>> cas_list = ["50-78-2", "50-00-0", "64-17-5"]
            >>> df = db.get_by_cas_batch(cas_list)
            >>> print(df[['cas', 'cmpdname', 'mf', 'mw']])
                   cas      cmpdname      mf       mw
            0  50-78-2       Aspirin  C9H8O4  180.160
            1  50-00-0  Formaldehyde    CH2O   30.026
            2  64-17-5       Ethanol   C2H6O   46.070
        """
        rows = []
        for cas in cas_list:
            result = self.get_by_cas(cas)
            if result:
                rows.append({'cas': cas, **self._compound_columns(result)})
        return pd.DataFrame(rows, columns=['cid', 'cas'] + self._compound_column_names()[1:])

    def _compound_column_names(self) -> List[str]:
        """The ``compounds`` table's columns, in table order, ``cid`` first."""
        return [row[1] for row in self.conn.execute("PRAGMA table_info(compounds)")]

    def _compound_columns(self, record: Dict[str, Any]) -> Dict[str, Any]:
        """A :meth:`get_by_cid` record without its list-valued extras."""
        return {name: record.get(name) for name in self._compound_column_names()}

    def get_id_table_from_cas(self, cas: str) -> Optional['pd.DataFrame']:
        """
        Get identifier table for a CAS number (similar to ZeroPM format).

        Args:
            cas (str): CAS Registry Number

        Returns:
            pandas.DataFrame: Table with columns [cid, cas, inchi, inchikey, smiles,
                             cmpdname, mf, mw] or None if not found

        Examples:
            >>> db = PubChemID()
            >>> df = db.get_id_table_from_cas("50-78-2")
            >>> df[['cid', 'cas', 'cmpdname', 'mf', 'mw']].to_dict('records')
            [{'cid': 2244, 'cas': '50-78-2', 'cmpdname': 'Aspirin', 'mf': 'C9H8O4', 'mw': 180.16}]
            >>> db.get_id_table_from_cas("0-00-0") is None
            True
        """
        import pandas as pd

        result = self.get_by_cas(cas)
        if not result:
            return None

        # Create DataFrame with main identifiers and properties
        df = pd.DataFrame([{
            'cid': result['cid'],
            'cas': cas,
            'inchi': result.get('inchi', ''),
            'inchikey': result.get('inchikey', ''),
            'smiles': result.get('smiles', ''),
            'cmpdname': result.get('cmpdname', ''),
            'mf': result.get('mf', ''),
            'mw': result.get('mw', None)
        }])

        return df

    def batch_get_id_table_from_cas(self, cas_list: List[str]) -> 'pd.DataFrame':
        """
        Get identifier tables for multiple CAS numbers.

        Args:
            cas_list (list): List of CAS Registry Numbers

        Returns:
            pandas.DataFrame: One :meth:`get_id_table_from_cas` row per CAS
            number found; CAS numbers not found are left out. Empty, with the
            same columns, when none is found.

        Examples:
            >>> db = PubChemID()
            >>> df = db.batch_get_id_table_from_cas(["50-78-2", "0-00-0", "64-17-5"])
            >>> print(df[['cid', 'cas', 'cmpdname', 'mf']])
                cid      cas cmpdname      mf
            0  2244  50-78-2  Aspirin  C9H8O4
            1   702  64-17-5  Ethanol   C2H6O
        """
        import pandas as pd

        tables = []
        for cas in cas_list:
            df = self.get_id_table_from_cas(cas)
            if df is not None:
                tables.append(df)

        if not tables:
            # Return empty DataFrame with correct columns
            return pd.DataFrame(columns=['cid', 'cas', 'inchi', 'inchikey',
                                        'smiles', 'cmpdname', 'mf', 'mw'])

        return pd.concat(tables, ignore_index=True)

    def get_by_smiles_batch(self, smiles_list: List[str]) -> 'pd.DataFrame':
        """
        Get complete compound information for multiple SMILES strings as a DataFrame.

        One row per SMILES found, carrying the compound's first CAS number and
        every column of the ``compounds`` table; see :meth:`get_by_cas_batch`
        for how those columns depend on where the database came from.

        Args:
            smiles_list (list): List of SMILES strings

        Returns:
            pandas.DataFrame: ``cid``, ``cas`` and then the ``compounds``
            columns. Empty, with those columns, when nothing is found.

        Examples:
            >>> db = PubChemID()
            >>> smiles_list = ["CC(=O)OC1=CC=CC=C1C(=O)O", "C", "CCO"]
            >>> df = db.get_by_smiles_batch(smiles_list)
            >>> print(df[['smiles', 'cmpdname', 'mf', 'mw']])
                                 smiles cmpdname      mf       mw
            0  CC(=O)OC1=CC=CC=C1C(=O)O  Aspirin  C9H8O4  180.160
            1                         C  Methane     CH4   16.043
            2                       CCO  Ethanol   C2H6O   46.070
        """
        rows = []
        for smiles in smiles_list:
            result = self.get_by_smiles(smiles)
            if result:
                cas_numbers = result.get('cas_numbers') or [None]
                rows.append({'cas': cas_numbers[0], **self._compound_columns(result)})
        return pd.DataFrame(rows, columns=['cid', 'cas'] + self._compound_column_names()[1:])

    def smiles_to_cas(self, smiles: str) -> Optional[List[str]]:
        """
        Convert SMILES string to CAS number(s).

        Unlike :meth:`smiles_to_cid`, this compares structures: the SMILES is
        converted to a standard InChI with RDKit and looked up by that, so any
        valid SMILES for the compound finds it.

        Args:
            smiles (str): SMILES string

        Returns:
            list: List of CAS numbers, or None if not found, if RDKit cannot
            parse the SMILES, or if RDKit is not installed.

        Examples:
            >>> db = PubChemID()
            >>> db.smiles_to_cas("CC(=O)OC1=CC=CC=C1C(=O)O")  # Aspirin
            ['50-78-2']
            >>> db.smiles_to_cas("OCC"), db.smiles_to_cid("OCC")
            (['64-17-5'], None)
        """
        # First convert SMILES to InChI using RDKit
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            inchi = Chem.MolToInchi(mol)
        except Exception:
            return None

        # Then look up by InChI
        return self.inchi_to_cas(inchi)

    def name_to_cas(self, name: str, exact: bool = True) -> Optional[List[str]]:
        """
        Convert chemical name to CAS number(s).

        Args:
            name (str): Chemical name or synonym
            exact (bool): If True, exact match only. If False, returns first match from search.

        Returns:
            list: The first matching compound's CAS numbers, or None if no
            compound matches. See :meth:`search_by_name` for how names match.

        Examples:
            >>> db = PubChemID()
            >>> db.name_to_cas("aspirin")
            ['50-78-2']
            >>> db.name_to_cas("no such compound") is None
            True

        Note:
            For exact=False, only the first match from the search is returned.
            Use search_by_name() for more control over multiple matches.
        """
        results = self.search_by_name(name, exact=exact, limit=1)
        if not results:
            return None
        return results[0].get('cas_numbers')

    def formula_to_cas(self, formula: str, limit: int = 100) -> Optional[List[str]]:
        """
        Convert molecular formula to CAS numbers.

        Note: Molecular formulas are not unique - many isomers can share the same formula.
        This method returns CAS numbers for all compounds matching the formula.

        Args:
            formula (str): Molecular formula (e.g., "C9H8O4", "CH2O")
            limit (int): Maximum number of compounds to retrieve

        Returns:
            list: The distinct CAS numbers of the first ``limit`` compounds
            with this formula, sorted as strings, or None if none is found

        Examples:
            >>> db = PubChemID()
            >>> cas_list = db.formula_to_cas("C9H8O4")
            >>> "50-78-2" in cas_list, cas_list == sorted(cas_list)
            (True, True)

        Warning:
            Can return many results for common formulas. Use limit parameter to control.
        """
        results = self.search_by_formula(formula, limit=limit)
        if not results:
            return None

        # Collect all unique CAS numbers from all matching compounds
        all_cas = []
        for compound in results:
            cas_numbers = compound.get('cas_numbers', [])
            if cas_numbers:
                all_cas.extend(cas_numbers)

        # Remove duplicates and sort
        unique_cas = sorted(set(all_cas))
        return unique_cas if unique_cas else None

    def batch_smiles_to_cas(self, smiles_list: List[str]) -> Dict[str, Optional[List[str]]]:
        """
        Convert multiple SMILES strings to CAS numbers.

        Args:
            smiles_list (list): List of SMILES strings

        Returns:
            dict: Mapping of SMILES -> list of CAS numbers (None if not found)

        Examples:
            >>> db = PubChemID()
            >>> db.batch_smiles_to_cas(["OCC", "not a smiles"])
            {'OCC': ['64-17-5'], 'not a smiles': None}
        """
        return {smiles: self.smiles_to_cas(smiles) for smiles in smiles_list}

    def batch_name_to_cas(self, name_list: List[str], exact: bool = True) -> Dict[str, Optional[List[str]]]:
        """
        Convert multiple chemical names to CAS numbers.

        Args:
            name_list (list): List of chemical names
            exact (bool): If True, exact match only

        Returns:
            dict: Mapping of name -> list of CAS numbers (None if not found)

        Examples:
            >>> db = PubChemID()
            >>> db.batch_name_to_cas(["aspirin", "ethanol", "xyzzy"])
            {'aspirin': ['50-78-2'], 'ethanol': ['64-17-5'], 'xyzzy': None}
        """
        return {name: self.name_to_cas(name, exact=exact) for name in name_list}

    def batch_formula_to_cas(self, formula_list: List[str], limit: int = 100) -> Dict[str, Optional[List[str]]]:
        """
        Convert multiple molecular formulas to CAS numbers.

        Args:
            formula_list (list): List of molecular formulas
            limit (int): Maximum number of compounds per formula

        Returns:
            dict: Mapping of formula -> list of CAS numbers (None if not found)

        Examples:
            >>> db = PubChemID()
            >>> results = db.batch_formula_to_cas(["H2O", "CH4", "XeF9"])
            >>> "7732-18-5" in results["H2O"], "74-82-8" in results["CH4"], results["XeF9"]
            (True, True, None)
        """
        return {formula: self.formula_to_cas(formula, limit=limit) for formula in formula_list}

    @property
    def api(self) -> 'PubChemAPI':
        """
        The online client used to answer what the local database cannot.

        Created on first use rather than in ``__init__``, so a strictly offline
        session never builds one.

        Returns:
            The :class:`PubChemAPI` instance passed to ``__init__``, or one
            created with default settings.

        Examples:
            >>> from provesid import PubChemAPI
            >>> api = PubChemAPI()
            >>> PubChemID(api=api).api is api
            True
        """
        if self._api is None:
            self._api = PubChemAPI()
        return self._api

    def properties(self, cid: Union[int, str],
                   properties: Optional[List[str]] = None,
                   use_online_fallback: bool = True) -> Optional[Dict[str, Any]]:
        """
        Look up computed properties for one compound, offline first.

        The local database answers from disk in microseconds; the online API is
        consulted only when the local database cannot serve the request, either
        because it holds no row for this CID or because a requested property is
        not one of the columns it carries (see :attr:`offline_properties`).

        Args:
            cid: PubChem Compound ID.
            properties: Property names to retrieve, e.g.
                ``['MolecularWeight', 'XLogP']``. Defaults to every property the
                local database can answer, :attr:`offline_properties`.
            use_online_fallback: When True (default), fall back to PUG-REST for
                anything the local database cannot answer. When False, the
                lookup is strictly offline, and a request the local database
                cannot answer in full --- an unknown CID, or any property
                outside :attr:`offline_properties` --- returns None.

        Returns:
            A dict carrying ``CID``, a ``Source`` of ``'offline'`` or
            ``'online'``, and one key per property that has a value. A property
            the compound has no value for is omitted rather than set to None,
            which is how PubChem itself reports it — so ``'XLogP' not in
            result`` means PubChem computes no logP for this compound, not that
            the lookup fell short. Returns None when neither source knows the
            CID, or when ``use_online_fallback`` is False and the request needs
            the network.

        Raises:
            ValueError: If ``cid`` is not an integer, or ``properties`` is an
                empty list.
            PubChemError: If the online fallback was needed and its request
                could not be completed. An incomplete answer is never passed off
                as a complete one.

        Examples:
            >>> db = PubChemID()
            >>> db.properties(2244, ['MolecularFormula', 'MolecularWeight'])
            {'CID': 2244, 'Source': 'offline', 'MolecularFormula': 'C9H8O4', 'MolecularWeight': 180.16}
            >>> # XLogP is PubChem's model output, never served from disk
            >>> db.properties(2244, ['XLogP'])['Source']            # doctest: +SKIP
            'online'
            >>> db.properties(2244, ['XLogP'], use_online_fallback=False) is None
            True
        """
        rows = self.properties_for_cids([cid], properties,
                                        use_online_fallback=use_online_fallback)
        return rows[0] if rows else None

    def properties_for_cids(self, cids: List[Union[int, str]],
                            properties: Optional[List[str]] = None,
                            use_online_fallback: bool = True,
                            chunk_size: int = PROPERTY_CHUNK_SIZE) -> List[Dict[str, Any]]:
        """
        Look up computed properties for many compounds, offline first.

        Everything the local database can answer is read in a handful of SQL
        statements; only the remainder is requested from PubChem, in bulk, a few
        hundred compounds per request. A list of ten thousand CIDs that the
        local database covers therefore costs no network traffic at all.

        Args:
            cids: PubChem Compound IDs. Duplicates are collapsed and the order
                of first appearance is preserved.
            properties: Property names to retrieve. Defaults to
                :attr:`offline_properties`.
            use_online_fallback: When True (default), CIDs the local database
                does not cover are requested from PUG-REST.
            chunk_size: How many CIDs to put in one online request.

        Returns:
            One dict per CID that could be answered, in the order requested,
            each carrying ``CID``, a ``Source`` of ``'offline'`` or
            ``'online'``, and one key per property that has a value. CIDs
            neither source knows are omitted; use :meth:`properties_table` to
            get a row for every CID asked about.

        Raises:
            ValueError: If a CID is not an integer, ``properties`` is an empty
                list, or ``chunk_size`` is not positive.
            PubChemError: If an online request could not be completed.

        Note:
            If *any* requested property lies outside
            :attr:`offline_properties`, the whole request goes online: the
            missing property would need a request per compound anyway, so
            splitting the property list between the two sources would cost the
            same traffic and return rows assembled from two different PubChem
            snapshots.

        Examples:
            >>> db = PubChemID()
            >>> rows = db.properties_for_cids([2244, 702], ['MolecularFormula'])
            >>> for row in rows:
            ...     print(row['CID'], row['MolecularFormula'], row['Source'])
            2244 C9H8O4 offline
            702 C2H6O offline
        """
        if properties is not None and not properties:
            raise ValueError("properties must name at least one property, or be None")
        if chunk_size <= 0:
            raise ValueError(f"chunk_size must be positive, got {chunk_size}")

        requested_properties = list(properties) if properties else list(self.offline_properties)
        wanted_cids = [self._coerce_cid(cid) for cid in cids]
        wanted_cids = list(dict.fromkeys(wanted_cids))
        if not wanted_cids:
            return []

        online_only = [name for name in requested_properties
                       if name not in self.offline_properties]

        found: Dict[int, Dict[str, Any]] = {}
        if online_only:
            self.logger.debug(
                "Going straight online for %d CIDs: %s not in the local database",
                len(wanted_cids), ', '.join(online_only))
            missing = wanted_cids
        else:
            found = self._offline_properties(wanted_cids, requested_properties)
            missing = [cid for cid in wanted_cids if cid not in found]
            self.logger.debug("Served %d/%d CIDs offline", len(found), len(wanted_cids))

        if missing and use_online_fallback:
            self.logger.debug("Falling back online for %d CIDs", len(missing))
            found.update(self._online_properties(missing, requested_properties, chunk_size))

        return [found[cid] for cid in wanted_cids if cid in found]

    def properties_table(self, cids: List[Union[int, str]],
                         properties: Optional[List[str]] = None,
                         use_online_fallback: bool = True,
                         chunk_size: int = PROPERTY_CHUNK_SIZE) -> 'pd.DataFrame':
        """
        Offline-first property lookup for many compounds, as a DataFrame.

        Same lookup as :meth:`properties_for_cids`, reshaped so that every CID
        asked about has a row whether or not it could be answered. That makes
        the frame safe to concatenate or join against the caller's own table.

        Args:
            cids: PubChem Compound IDs. Duplicates are collapsed.
            properties: Property names to retrieve. Defaults to
                :attr:`offline_properties`.
            use_online_fallback: When True (default), consult PUG-REST for CIDs
                the local database does not cover.
            chunk_size: How many CIDs to put in one online request.

        Returns:
            A DataFrame with one row per distinct CID in the order requested.
            Columns are ``CID``, ``Source`` and the requested properties.
            ``Source`` reads ``'offline'``, ``'online'``, or ``'missing'`` for a
            CID neither source knows; a property with no value is NaN/None.

        Raises:
            ValueError: If a CID is not an integer, ``properties`` is an empty
                list, or ``chunk_size`` is not positive.
            PubChemError: If an online request could not be completed.

        Examples:
            >>> db = PubChemID()
            >>> table = db.properties_table([2244, 702], ['MolecularWeight'])
            >>> table[['CID', 'MolecularWeight', 'Source']].to_dict('records')
            [{'CID': 2244, 'MolecularWeight': 180.16, 'Source': 'offline'},
             {'CID': 702, 'MolecularWeight': 46.07, 'Source': 'offline'}]
        """
        requested_properties = list(properties) if properties else list(self.offline_properties)
        rows = self.properties_for_cids(cids, requested_properties,
                                        use_online_fallback=use_online_fallback,
                                        chunk_size=chunk_size)
        return self._table(cids, rows, requested_properties)

    def descriptors(self, cid: Union[int, str],
                    descriptors: Optional[List[str]] = None,
                    source: str = 'rdkit',
                    use_online_fallback: bool = True) -> Optional[Dict[str, Any]]:
        """
        Computed molecular descriptors for one compound, from RDKit or PubChem.

        The local database stores identifiers and structures, not descriptors:
        XLogP, TPSA and the counts are the output of a model run over the
        structure, and there is more than one model. This method runs one,
        and says which:

        * ``source="rdkit"`` (default) computes them with RDKit from the
          compound's stored SMILES --- no network, milliseconds. The record
          says ``Source='rdkit'``. RDKit and PubChem count some things
          differently, and the logP is a different model altogether, named
          ``MolLogP`` rather than ``XLogP``; :func:`rdkit_descriptors`
          measures how far apart they are. ``Complexity`` is not available.
        * ``source="pubchem"`` fetches PubChem's own values from PUG-REST,
          through the same path as :meth:`properties`, labelled
          ``Source='online'``. This is the only way to PubChem's ``XLogP``
          and ``Complexity``.

        Args:
            cid: PubChem Compound ID.
            descriptors: Names to compute. Defaults to every descriptor the
                source has: :data:`RDKIT_DESCRIPTORS` or
                :data:`PUBCHEM_DESCRIPTORS`.
            source: ``'rdkit'`` or ``'pubchem'``.
            use_online_fallback: For ``source="rdkit"``, whether a compound the
                local database does not hold may have its SMILES fetched from
                PubChem to compute from. When False, such a compound returns
                None. ``source="pubchem"`` is online by definition and does not
                accept False.

        Returns:
            A dict carrying ``CID``, ``Source`` and one key per descriptor that
            has a value, or None when the compound is unknown. A compound whose
            SMILES RDKit cannot parse --- a handful in a million --- or that
            has no structure comes back with ``CID`` and ``Source`` only.

        Raises:
            ValueError: If ``cid`` is not an integer, ``source`` is unknown,
                a name is not one ``source`` provides (asking RDKit for
                ``XLogP`` says to ask for ``MolLogP``), or ``source="pubchem"``
                is combined with ``use_online_fallback=False``.
            PubChemError: If an online request could not be completed.

        Examples:
            >>> db = PubChemID()
            >>> db.descriptors(2244, ['MolLogP', 'TPSA'])
            {'CID': 2244, 'Source': 'rdkit', 'MolLogP': 1.3101, 'TPSA': 63.6}
            >>> db.descriptors(2244, ['XLogP'], source='pubchem')  # doctest: +SKIP
            {'CID': 2244, 'Source': 'online', 'XLogP': 1.2}
        """
        rows = self.descriptors_for_cids([cid], descriptors, source=source,
                                         use_online_fallback=use_online_fallback)
        return rows[0] if rows else None

    def descriptors_for_cids(self, cids: List[Union[int, str]],
                             descriptors: Optional[List[str]] = None,
                             source: str = 'rdkit',
                             use_online_fallback: bool = True,
                             chunk_size: int = PROPERTY_CHUNK_SIZE) -> List[Dict[str, Any]]:
        """
        Computed molecular descriptors for many compounds, from RDKit or PubChem.

        The list form of :meth:`descriptors`. With ``source="rdkit"`` the
        SMILES of every compound in the local database are read in a handful
        of statements, and only those it lacks are fetched from PubChem, in
        bulk; with ``source="pubchem"`` the whole list goes to PUG-REST a few
        hundred compounds per request.

        Args:
            cids: PubChem Compound IDs. Duplicates are collapsed and the order
                of first appearance is preserved.
            descriptors: Names to compute; defaults to every descriptor the
                source has.
            source: ``'rdkit'`` or ``'pubchem'``.
            use_online_fallback: For ``source="rdkit"``, whether SMILES missing
                from the local database may be fetched from PubChem.
            chunk_size: How many CIDs to put in one online request.

        Returns:
            One dict per CID that could be answered, in the order requested,
            shaped as :meth:`descriptors` describes. CIDs no source knows are
            omitted; :meth:`descriptors_table` gives a row for every CID.

        Raises:
            ValueError: As for :meth:`descriptors`, or if ``chunk_size`` is not
                positive.
            PubChemError: If an online request could not be completed.

        Examples:
            >>> db = PubChemID()
            >>> for row in db.descriptors_for_cids([2244, 702], ['HeavyAtomCount']):
            ...     print(row)
            {'CID': 2244, 'Source': 'rdkit', 'HeavyAtomCount': 13}
            {'CID': 702, 'Source': 'rdkit', 'HeavyAtomCount': 3}
        """
        names = _check_descriptor_names(descriptors, source)

        if source == 'pubchem':
            if not use_online_fallback:
                raise ValueError("source='pubchem' fetches PubChem's values online; "
                                 "use_online_fallback=False contradicts it. For "
                                 "descriptors without the network use source='rdkit'.")
            return self.properties_for_cids(cids, names, chunk_size=chunk_size)

        structures = self.properties_for_cids(cids, ['SMILES'],
                                              use_online_fallback=use_online_fallback,
                                              chunk_size=chunk_size)
        rows = []
        for structure in structures:
            values = rdkit_descriptors(structure.get('SMILES'), names)
            if values is None:
                self.logger.debug("RDKit could not read the SMILES of CID %d: %r",
                                  structure['CID'], structure.get('SMILES'))
            rows.append({'CID': structure['CID'], 'Source': 'rdkit', **(values or {})})
        return rows

    def descriptors_table(self, cids: List[Union[int, str]],
                          descriptors: Optional[List[str]] = None,
                          source: str = 'rdkit',
                          use_online_fallback: bool = True,
                          chunk_size: int = PROPERTY_CHUNK_SIZE) -> 'pd.DataFrame':
        """
        Computed molecular descriptors for many compounds, as a DataFrame.

        Same lookup as :meth:`descriptors_for_cids`, with a row for every CID
        asked about, so the frame joins safely against the caller's own table.
        Because the RDKit and PubChem columns share names wherever the quantity
        is the same, two tables built with each source line up column for
        column, apart from ``MolLogP`` / ``XLogP`` and ``Complexity``.

        Args:
            cids: PubChem Compound IDs. Duplicates are collapsed.
            descriptors: Names to compute; defaults to every descriptor the
                source has.
            source: ``'rdkit'`` or ``'pubchem'``.
            use_online_fallback: For ``source="rdkit"``, whether SMILES missing
                from the local database may be fetched from PubChem.
            chunk_size: How many CIDs to put in one online request.

        Returns:
            A DataFrame with one row per distinct CID in the order requested.
            Columns are ``CID``, ``Source`` and the descriptors. ``Source``
            reads ``'rdkit'``, ``'online'``, or ``'missing'`` for a CID no
            source knows; a descriptor with no value is NaN/None.

        Raises:
            ValueError: As for :meth:`descriptors_for_cids`.
            PubChemError: If an online request could not be completed.

        Examples:
            >>> db = PubChemID()
            >>> db.descriptors_table([2244, 702], ['TPSA'])
                CID Source   TPSA
            0  2244  rdkit  63.60
            1   702  rdkit  20.23
        """
        names = _check_descriptor_names(descriptors, source)
        rows = self.descriptors_for_cids(cids, names, source=source,
                                         use_online_fallback=use_online_fallback,
                                         chunk_size=chunk_size)
        return self._table(cids, rows, names)

    def _table(self, cids: List[Union[int, str]], rows: List[Dict[str, Any]],
               names: List[str]) -> 'pd.DataFrame':
        """
        Reshape looked-up records into a frame with a row for every CID asked about.

        Args:
            cids: The CIDs as the caller gave them; duplicates are collapsed.
            rows: Records carrying ``CID``, ``Source`` and values, as
                :meth:`properties_for_cids` returns them.
            names: The value columns, in order.

        Returns:
            A DataFrame with columns ``CID``, ``Source`` and ``names``, where a
            CID with no record has ``Source='missing'``.
        """
        by_cid = {row['CID']: row for row in rows}
        records = []
        for cid in dict.fromkeys(self._coerce_cid(cid) for cid in cids):
            row = by_cid.get(cid, {'CID': cid, 'Source': 'missing'})
            records.append({'CID': cid, 'Source': row['Source'],
                            **{name: row.get(name) for name in names}})

        return pd.DataFrame(records, columns=['CID', 'Source'] + names)

    @staticmethod
    def _coerce_cid(cid: Union[int, str]) -> int:
        """
        Normalise a CID to an int so that ``2244`` and ``"2244"`` share a row.

        Args:
            cid: CID as an int or a string of digits.

        Returns:
            The CID as an int.

        Raises:
            ValueError: If ``cid`` is not an integer. PubChem answers a
                malformed CID with a blanket ``PUGREST.BadRequest`` that fails
                the whole batch, so it is worth catching here, where the
                offending value can be named.
        """
        try:
            return int(cid)
        except (TypeError, ValueError):
            raise ValueError(f"CID must be an integer, got {cid!r}")

    def _offline_properties(self, cids: List[int],
                            properties: List[str]) -> Dict[int, Dict[str, Any]]:
        """
        Read properties for the given CIDs from the local database.

        Args:
            cids: CIDs to look up, already coerced to int.
            properties: Property names, all of which must be keys of
                :attr:`offline_properties`.

        Returns:
            A dict keyed by CID, holding one record per CID present in the
            database. A record carries ``CID``, ``Source='offline'`` and the
            properties that have a value; a NULL column is left out, matching
            PubChem, which omits a property rather than reporting it as null.
        """
        columns = [self.offline_properties[name] for name in properties]
        cursor = self.conn.cursor()
        found: Dict[int, Dict[str, Any]] = {}

        # SQLite allows a limited number of bound parameters per statement
        # (999 by default), so the IN list is filled in batches.
        for start in range(0, len(cids), self._SQL_PARAMETER_LIMIT):
            batch = cids[start:start + self._SQL_PARAMETER_LIMIT]
            placeholders = ','.join('?' * len(batch))
            cursor.execute(
                f"SELECT cid, {', '.join(columns)} FROM compounds "
                f"WHERE cid IN ({placeholders})", batch)
            for row in cursor.fetchall():
                record = {'CID': row['cid'], 'Source': 'offline'}
                for name, column in zip(properties, columns):
                    value = row[column]
                    if value is not None and value != '':
                        record[name] = self._cast_property(name, value)
                found[row['cid']] = record

        return found

    def _online_properties(self, cids: List[int], properties: List[str],
                           chunk_size: int) -> Dict[int, Dict[str, Any]]:
        """
        Fetch properties for the given CIDs from PUG-REST.

        Args:
            cids: CIDs the local database could not answer.
            properties: Property names to request.
            chunk_size: How many CIDs to put in one request.

        Returns:
            A dict keyed by CID, holding one record per CID PubChem answered
            for, shaped like the offline records but with
            ``Source='online'``. PubChem returns a bare CID for a compound it
            has no record of; such a row is dropped, so the CID is reported as
            unknown rather than as a compound with no properties.

        Raises:
            PubChemError: If a request could not be completed.
        """
        rows = self.api.get_properties_for_cids(cids, properties, chunk_size=chunk_size)

        found: Dict[int, Dict[str, Any]] = {}
        for row in rows:
            values = {name: self._cast_property(name, row[name])
                      for name in properties
                      if row.get(name) is not None and row.get(name) != ''}
            if not values:
                continue
            cid = self._coerce_cid(row['CID'])
            found[cid] = {'CID': cid, 'Source': 'online', **values}

        return found

    @classmethod
    def _cast_property(cls, name: str, value: Any) -> Any:
        """
        Coerce a property value to one consistent type across both sources.

        The two sources disagree on types for the same property: PUG-REST
        returns ``MolecularWeight`` as the string ``"180.16"`` while the local
        database holds it as a float, and a table assembled from both sources
        has to be usable as one table.

        Args:
            name: Property name.
            value: Raw value from either source.

        Returns:
            The value cast to the type recorded in ``_PROPERTY_CASTS``, or
            unchanged when no cast is recorded or the cast does not apply. An
            uncastable value is returned as-is rather than discarded: a
            surprising value is more useful to the caller than a silent hole.
        """
        cast = cls._PROPERTY_CASTS.get(name)
        if cast is None:
            return value
        try:
            return cast(value)
        except (TypeError, ValueError):
            logging.debug("Could not cast %s=%r with %s", name, value, cast.__name__)
            return value

    def provenance(self) -> Dict[str, Any]:
        """
        Where this database came from and how it was built.

        A database built by :func:`provesid.pubchem_ftp.build_pubchem_id_db`
        records its PubChem release, the snapshot's timestamp, the URL and MD5
        of every source file, the row counts and the build time. That is what
        makes a lookup against it citable: the release pins down exactly which
        state of PubChem answered.

        Returns:
            A dict of the ``provenance`` table's entries, plus ``files``: one
            dict per source file with ``file``, ``url``, ``md5``, ``bytes``,
            ``lines_read`` and ``rows_kept``. Empty for a database made before
            provenance was recorded --- every Zenodo copy so far.

        Examples:
            >>> db = PubChemID()                                  # doctest: +SKIP
            >>> db.provenance()["release"]                        # doctest: +SKIP
            '2026-09-01'
        """
        tables = {row[0] for row in self.conn.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'")}
        if "provenance" not in tables:
            return {}
        record: Dict[str, Any] = dict(
            self.conn.execute("SELECT key, value FROM provenance").fetchall())
        record["files"] = [dict(row) for row in self.conn.execute(
            "SELECT * FROM provenance_files ORDER BY rowid")]
        return record

    def xrefs(self, cid: Union[int, str]) -> Dict[str, List[str]]:
        """
        Identifiers other databases give this compound, as PubChem links them.

        PubChem publishes these links itself, in the same file the CAS
        numbers come from, so they cost nothing to keep: DSSTox substance IDs
        (``dtxsid``), ChEBI IDs, ChEMBL IDs, EC numbers and UNIIs --- see
        :data:`provesid.pubchem_ftp.XREF_TYPES`.

        Args:
            cid: PubChem Compound ID.

        Returns:
            A dict from source (``"dtxsid"``, ``"chebi"``, ``"chembl"``,
            ``"ec"``, ``"unii"``) to that source's identifiers for the
            compound, sorted. Sources with none are left out, so a compound
            with no links returns ``{}``.

        Raises:
            ValueError: If ``cid`` is not an integer.
            RuntimeError: If the database has no ``xrefs`` table --- a Zenodo
                copy. The message says how to build one that does.

        Examples:
            >>> db = PubChemID()                                  # doctest: +SKIP
            >>> db.xrefs(2244)                                    # doctest: +SKIP
            {'chebi': ['CHEBI:15365'], 'chembl': ['CHEMBL25'],
             'dtxsid': ['DTXSID5020108'], 'ec': ['200-064-1'],
             'unii': ['R16CO5Y76E']}
        """
        cid = self._coerce_cid(cid)
        tables = {row[0] for row in self.conn.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'")}
        if "xrefs" not in tables:
            raise RuntimeError(
                f"{self.db_path} has no cross-references. They exist only in a "
                "database built from PubChem's FTP site: "
                "provesid.pubchem_ftp.build_pubchem_id_db(force=True)."
            )
        found: Dict[str, List[str]] = {}
        for source, identifier in self.conn.execute(
                "SELECT source, identifier FROM xrefs WHERE cid = ? "
                "ORDER BY source, identifier", (cid,)):
            found.setdefault(source, []).append(identifier)
        return found

    def get_stats(self) -> Dict[str, int]:
        """
        Get database statistics.

        Returns:
            dict: ``total_compounds``, ``total_cas_numbers`` (rows in the CAS
            table), ``compounds_with_cas``, ``total_synonyms``,
            ``compounds_with_inchikey``, ``database_path`` and
            ``database_size_mb``. The counts depend on the release.

        Examples:
            >>> db = PubChemID()
            >>> stats = db.get_stats()
            >>> print(f"Total compounds: {stats['total_compounds']:,}")  # doctest: +SKIP
            Total compounds: 1,589,910
            >>> stats['compounds_with_cas'] <= stats['total_compounds']
            True
        """
        cursor = self.conn.cursor()

        cursor.execute("SELECT COUNT(*) FROM compounds")
        total_compounds = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM cas_numbers")
        total_cas = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(DISTINCT cid) FROM cas_numbers")
        compounds_with_cas = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM synonyms")
        total_synonyms = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM compounds WHERE inchikey IS NOT NULL AND inchikey != ''")
        compounds_with_inchikey = cursor.fetchone()[0]

        return {
            'total_compounds': total_compounds,
            'total_cas_numbers': total_cas,
            'compounds_with_cas': compounds_with_cas,
            'total_synonyms': total_synonyms,
            'compounds_with_inchikey': compounds_with_inchikey,
            'database_path': self.db_path,
            'database_size_mb': os.path.getsize(self.db_path) / (1024**2)
        }
