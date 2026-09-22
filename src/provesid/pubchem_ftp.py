"""
Build the PubChem identifier database from PubChem's own FTP files.

:class:`~provesid.PubChemID` answers CAS, name, InChIKey and formula lookups
from one SQLite file, ``pubchem_id.db``. That file used to come from a manual
pipeline: someone downloaded a CSV from PubChem's classification browser,
pulled CAS numbers out of its free-text synonym column with a regular
expression, and uploaded the result to Zenodo. Nothing recorded which PubChem
release it came from, and 17 379 of its CAS numbers --- 1.25% --- fail the CAS
check digit, because ``\\d{2,7}-\\d{2}-\\d`` matches plenty of things that are
not registry numbers.

This module builds the same database from ``Compound/Extras/`` on PubChem's
FTP site, in one call, from a dated monthly snapshot:

* **Scope** comes from ``CID-Identifiers.tsv.gz``, PubChem's curated mapping of
  compounds to third-party identifiers. Its ``CAS`` rows, each checked with
  :func:`provesid.utils.check_CASRN`, decide which compounds go in --- about
  1.43 M, with 123 CAS rows rejected where the regex let 17 379 through.
* **Columns** come from one file each: title, formula and masses, isomeric
  SMILES, IUPAC name, InChI and InChIKey, creation date, and the filtered
  synonym list.
* **Cross-references** --- DTXSID, ChEBI, ChEMBL, EC and UNII --- come from the
  same identifier file at no extra cost, and go into an ``xrefs`` table.
* **Provenance** is written into the database itself: the release, every
  source file's URL and MD5, the row counts and the build time. A database on
  disk can always say where it came from.

The files are processed one at a time --- downloaded, streamed, filtered to the
compounds in scope and deleted --- so the disk needed at any moment is the
database plus the largest single file (7.4 GB, ``CID-InChI-Key.gz``), not the
15.4 GB total.

The eight computed descriptors of the old database (XLogP, TPSA, complexity,
charge and the four counts) are not in any of these files and are not stored:
they are properties of the structure rather than data about the substance, and
belong to an on-demand calculation instead.

Example:
    >>> from provesid.pubchem_ftp import build_pubchem_id_db, list_releases
    >>> list_releases()                                    # doctest: +SKIP
    ['2026-09-01', '2026-08-01', '2026-07-01', '2026-06-01', 'current']
    >>> build_pubchem_id_db()                              # doctest: +SKIP
    '/home/me/.local/share/provesid/pubchem_id.db'
"""

import datetime
import gzip
import logging
import os
import re
import sqlite3
from dataclasses import dataclass
from decimal import ROUND_HALF_UP, Decimal
from typing import Callable, Dict, Iterator, List, Optional, Set, Tuple

import requests
from tqdm import tqdm

from .datasets import DownloadError, download_file, md5_of_file, read_checksum
from .utils import check_CASRN, user_dataset_path

logger = logging.getLogger(__name__)

__all__ = [
    "FTP_ROOT",
    "LATEST",
    "CURRENT",
    "XREF_TYPES",
    "SOURCE_FILES",
    "build_pubchem_id_db",
    "list_releases",
    "resolve_release",
    "extras_url",
    "molecular_weight",
]

#: Root of PubChem's compound files. The HTTPS mirror of the FTP site answers
#: ``Range`` requests and publishes an ``.md5`` beside every file, which is
#: what :func:`~provesid.datasets.download_file` needs to resume and verify.
FTP_ROOT = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound"

#: Release name meaning "the newest monthly snapshot", the default.
LATEST = "latest"

#: Release name meaning PubChem's rolling ``Compound/Extras/``, regenerated with
#: every dump. Fresher than any snapshot, but not reproducible: the files can
#: change between two builds, or during one.
CURRENT = "current"

#: Version of the schema and the procedure this module writes, recorded in
#: every database it builds. Bump it when either changes.
BUILDER_VERSION = 1

#: Default name of the directory, beside the database, that source files are
#: downloaded into, one subdirectory per release.
DOWNLOAD_DIRNAME = "pubchem_ftp"

#: Identifier types from ``CID-Identifiers.tsv.gz`` stored in ``xrefs``, mapped
#: to the short name the table uses. These are the identifiers
#: :class:`~provesid.Search` otherwise reconciles across sources by matching
#: structures; PubChem publishes the links directly.
XREF_TYPES: Dict[str, str] = {
    "DSSTox Substance ID": "dtxsid",
    "ChEBI ID": "chebi",
    "ChEMBL ID": "chembl",
    "European Community (EC) Number": "ec",
    "UNII": "unii",
}

_CAS_TYPE = b"CAS"

#: Rows written per ``executemany``. Large enough that the statement overhead
#: vanishes, small enough that a batch of synonyms stays a few megabytes.
_BATCH_ROWS = 50_000

#: Lines between progress-bar updates while a file is streamed.
_PROGRESS_EVERY = 1_000_000

_RELEASE_PATTERN = re.compile(r"^\d{4}-\d{2}-\d{2}$")
_LISTED_RELEASE = re.compile(r'href="(\d{4}-\d{2}-\d{2})/"')


# ──────────────────────────────────────────────────────────────────────────────
# Molecular weight
#
# ``CID-Mass.gz`` carries the formula, the monoisotopic mass and the exact
# mass, but not the molecular weight, and no file on the FTP site does. So the
# weight is computed from the formula. PubChem's own weights could not be
# reproduced exactly: they are rounded to between zero and three decimals by a
# rule that no function of the formula's elements and counts fits, and 7% of
# them match no rounding of any single weight table. What *can* be matched is
# the atomic weights. The table below is IUPAC 2005 --- the values PubChem's C,
# H, N, O and P weights agree with --- with twelve elements set instead to the
# values that reproduce PubChem's weights best, fitted against the 1.58 M
# molecular weights in PubChem's 2026-01 CAS export.
# ──────────────────────────────────────────────────────────────────────────────

#: Standard atomic weights used for ``mw``, by element symbol. Elements missing
#: here (the radioactive ones without a standard weight) fall back to RDKit's
#: periodic table.
ATOMIC_WEIGHTS: Dict[str, float] = {
    "H": 1.00794, "He": 4.002602, "Li": 6.941, "Be": 9.012182, "B": 10.812,
    "C": 12.0107, "N": 14.0067, "O": 15.9994, "F": 18.9984032, "Ne": 20.1797,
    "Na": 22.98976928, "Mg": 24.3050, "Al": 26.9815386, "Si": 28.085,
    "P": 30.973762, "S": 32.067, "Cl": 35.45, "Ar": 39.948, "K": 39.0983,
    "Ca": 40.078, "Sc": 44.955912, "Ti": 47.867, "V": 50.9415, "Cr": 51.9961,
    "Mn": 54.938045, "Fe": 55.845, "Co": 58.933195, "Ni": 58.6934,
    "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.630, "As": 74.92160,
    "Se": 78.971, "Br": 79.904, "Kr": 83.798, "Rb": 85.4678, "Sr": 87.62,
    "Y": 88.90585, "Zr": 91.224, "Nb": 92.90638, "Mo": 95.95, "Ru": 101.07,
    "Rh": 102.90550, "Pd": 106.42, "Ag": 107.8682, "Cd": 112.414,
    "In": 114.818, "Sn": 118.710, "Sb": 121.760, "Te": 127.60, "I": 126.904,
    "Xe": 131.293, "Cs": 132.9054519, "Ba": 137.327, "La": 138.90547,
    "Ce": 140.116, "Pr": 140.90765, "Nd": 144.242, "Sm": 150.36,
    "Eu": 151.964, "Gd": 157.25, "Tb": 158.92535, "Dy": 162.500,
    "Ho": 164.93032, "Er": 167.259, "Tm": 168.93421, "Yb": 173.045,
    "Lu": 174.9668, "Hf": 178.486, "Ta": 180.94788, "W": 183.84,
    "Re": 186.207, "Os": 190.23, "Ir": 192.217, "Pt": 195.084,
    "Au": 196.966569, "Hg": 200.592, "Tl": 204.3833, "Pb": 207.2,
    "Bi": 208.98040, "Th": 232.03806, "Pa": 231.03588, "U": 238.02891,
}

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")
_FORMULA_SHAPE = re.compile(r"^(?:[A-Z][a-z]?\d*)+(?:[+-]\d*)?$")


def _atomic_weight(symbol: str) -> Optional[float]:
    """Weight of one element from :data:`ATOMIC_WEIGHTS`, or RDKit's, or None."""
    weight = ATOMIC_WEIGHTS.get(symbol)
    if weight is not None:
        return weight
    from rdkit import Chem  # a hard dependency, but slow to import

    try:
        weight = Chem.GetPeriodicTable().GetAtomicWeight(symbol)
    except RuntimeError:
        return None
    return weight or None


def _round_weight(value: float) -> float:
    """Round half up to two decimals, the precision PubChem uses most often."""
    return float(Decimal(repr(value)).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))


def molecular_weight(formula: str) -> Optional[float]:
    """
    Molecular weight of a PubChem molecular formula, in g/mol.

    Computed with :data:`ATOMIC_WEIGHTS` and rounded to two decimals. A charge
    suffix (``+``, ``-2``) is ignored, as PubChem ignores it: the weight of an
    ion is the weight of its atoms.

    The result agrees with PubChem's own ``MolecularWeight`` to the second
    decimal for most compounds, but PubChem rounds some weights more coarsely
    --- to one decimal, or to a whole number for compounds of lead --- and
    where it does, this is the more precise of the two.

    A formula cannot describe isotopic labelling: PubChem writes
    chloroform-*d* as ``CHCl3``. The builder corrects those compounds from
    their SMILES; this function cannot.

    Args:
        formula: A formula as PubChem writes it, e.g. ``"C9H8O4"`` or
            ``"C9H18NO4+"``.

    Returns:
        The weight, or None when the formula is empty, malformed or names an
        element with no known weight.

    Example:
        >>> molecular_weight("C9H8O4")
        180.16
        >>> molecular_weight("C9H18NO4+")
        204.24
        >>> molecular_weight("not a formula") is None
        True
    """
    if not formula or not _FORMULA_SHAPE.match(formula):
        return None
    body = re.sub(r"[+-]\d*$", "", formula)
    total = 0.0
    for symbol, count in _FORMULA_TOKEN.findall(body):
        weight = _atomic_weight(symbol)
        if weight is None:
            return None
        total += weight * (int(count) if count else 1)
    return _round_weight(total)


def _isotopic_molecular_weight(smiles: str) -> Optional[float]:
    """
    Molecular weight of an isotopically labelled structure, from its SMILES.

    Labelled atoms weigh their isotope's mass; every other atom, implicit
    hydrogens included, weighs its entry in :data:`ATOMIC_WEIGHTS`.

    Args:
        smiles: A SMILES carrying at least one isotope label.

    Returns:
        The weight rounded to two decimals, or None if RDKit cannot parse it.
    """
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    table = Chem.GetPeriodicTable()
    total = 0.0
    for atom in molecule.GetAtoms():
        symbol = atom.GetSymbol()
        isotope = atom.GetIsotope()
        weight = table.GetMassForIsotope(symbol, isotope) if isotope else _atomic_weight(symbol)
        if weight is None:
            return None
        total += weight + atom.GetTotalNumHs() * ATOMIC_WEIGHTS["H"]
    return _round_weight(total)


# ──────────────────────────────────────────────────────────────────────────────
# Releases
# ──────────────────────────────────────────────────────────────────────────────

def _root(base_url: Optional[str]) -> str:
    """The FTP root to use, read at call time so a test or mirror can move it."""
    return (base_url or FTP_ROOT).rstrip("/")


def list_releases(*, base_url: Optional[str] = None,
                  session: Optional[requests.Session] = None,
                  timeout: float = 30) -> List[str]:
    """
    The PubChem releases a database can be built from, newest first.

    Monthly snapshots live under ``Compound/Monthly/YYYY-MM-01/`` and are
    frozen once published, so a database built from one is reproducible and
    can be cited. PubChem keeps the last few months. ``"current"`` --- the
    rolling ``Compound/Extras/`` --- is always listed last.

    Args:
        base_url: PubChem compound root. Defaults to :data:`FTP_ROOT`.
        session: ``requests.Session`` to fetch through.
        timeout: Seconds to wait for the listing.

    Returns:
        Snapshot dates as ``YYYY-MM-DD`` strings, newest first, followed by
        ``"current"``.

    Raises:
        DownloadError: If the listing cannot be fetched.

    Example:
        >>> list_releases()                                  # doctest: +SKIP
        ['2026-09-01', '2026-08-01', '2026-07-01', '2026-06-01', 'current']
    """
    url = f"{_root(base_url)}/Monthly/"
    getter = session.get if session is not None else requests.get
    try:
        response = getter(url, timeout=timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise DownloadError(f"Could not list PubChem releases at {url}: {exc}",
                            url=url) from exc
    dates = sorted(set(_LISTED_RELEASE.findall(response.text)), reverse=True)
    return dates + [CURRENT]


def resolve_release(release: str = LATEST, *, base_url: Optional[str] = None,
                    session: Optional[requests.Session] = None) -> str:
    """
    Turn a release argument into a concrete release name.

    Args:
        release: ``"latest"`` for the newest monthly snapshot, ``"current"``
            for the rolling dump, or a snapshot date such as ``"2026-09-01"``.
        base_url: PubChem compound root. Defaults to :data:`FTP_ROOT`.
        session: ``requests.Session`` to fetch the listing through.

    Returns:
        ``"current"`` or a ``YYYY-MM-DD`` date. ``"latest"`` costs one request
        for the listing; the other two cost none, so a date that PubChem no
        longer keeps is only discovered when its first file is requested.

    Raises:
        ValueError: If ``release`` is none of the three forms.
        DownloadError: If ``"latest"`` was asked for and the listing cannot be
            fetched, or lists no snapshot.

    Example:
        >>> resolve_release("2026-09-01"), resolve_release("current")
        ('2026-09-01', 'current')
        >>> resolve_release()                                  # doctest: +SKIP
        '2026-09-01'
        >>> resolve_release("yesterday")
        Traceback (most recent call last):
        ...
        ValueError: release='yesterday' is not a PubChem release. ...
    """
    if release == CURRENT or _RELEASE_PATTERN.match(release or ""):
        return release
    if release != LATEST:
        raise ValueError(
            f"release={release!r} is not a PubChem release. Use {LATEST!r}, "
            f"{CURRENT!r} or a snapshot date such as '2026-09-01' "
            f"(see provesid.pubchem_ftp.list_releases())."
        )
    snapshots = [name for name in list_releases(base_url=base_url, session=session)
                 if name != CURRENT]
    if not snapshots:
        raise DownloadError(f"No monthly snapshot is listed under {_root(base_url)}/Monthly/")
    return snapshots[0]


def release_url(release: str, *, base_url: Optional[str] = None) -> str:
    """
    URL of a concrete release's directory: ``Monthly/<date>`` or the root.

    Args:
        release: ``"current"`` or a snapshot date, as :func:`resolve_release`
            returns.
        base_url: PubChem compound root. Defaults to :data:`FTP_ROOT`.

    Returns:
        The directory URL, without a trailing slash.

    Example:
        >>> release_url("2026-09-01")
        'https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Monthly/2026-09-01'
        >>> release_url("current")
        'https://ftp.ncbi.nlm.nih.gov/pubchem/Compound'
    """
    if release == CURRENT:
        return _root(base_url)
    return f"{_root(base_url)}/Monthly/{release}"


def extras_url(release: str, *, base_url: Optional[str] = None) -> str:
    """
    URL of the ``Extras/`` directory of a concrete release.

    Args:
        release: ``"current"`` or a snapshot date, as :func:`resolve_release`
            returns.
        base_url: PubChem compound root. Defaults to :data:`FTP_ROOT`.

    Returns:
        The directory URL, without a trailing slash.

    Example:
        >>> extras_url("2026-09-01")
        'https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Monthly/2026-09-01/Extras'
        >>> extras_url("current")
        'https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras'
    """
    return f"{release_url(release, base_url=base_url)}/Extras"


# ──────────────────────────────────────────────────────────────────────────────
# Source files
# ──────────────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class SourceFile:
    """
    One file of ``Compound/Extras/`` and where its fields go.

    Attributes:
        key: Short name, used in logs and in the ``provenance_files`` table.
        filename: Name of the file under ``Extras/``.
        columns: ``compounds`` columns this file fills, in the order its
            converter returns them. Empty for the files that fill other
            tables.
        fields: Tab-separated fields after the CID on each line. Usually one
            per column; ``CID-Mass.gz`` has three and fills four, because the
            molecular weight is computed from its formula.
        approx_bytes: Compressed size in the 2026-09-01 snapshot, for the
            estimate a user sees before the build starts. Later snapshots are
            a little larger.

    Example:
        >>> [source.key for source in SOURCE_FILES]
        ['identifiers', 'date', 'mass', 'smiles', 'title', 'iupac', 'inchi', 'synonyms']
        >>> SOURCE_FILES[0].filename
        'CID-Identifiers.tsv.gz'
    """

    key: str
    filename: str
    columns: Tuple[str, ...]
    fields: int
    approx_bytes: int


#: Every file the builder can read, in the order it reads them. The identifier
#: file comes first because it decides which compounds are in scope; the rest
#: are filtered by that decision.
SOURCE_FILES: Tuple[SourceFile, ...] = (
    SourceFile("identifiers", "CID-Identifiers.tsv.gz", (), 2, 98_179_955),
    SourceFile("date", "CID-Date.gz", ("cidcdate",), 1, 330_091_504),
    SourceFile("mass", "CID-Mass.gz", ("mf", "mw", "monoisotopicmass", "exactmass"), 3,
               1_390_299_686),
    SourceFile("smiles", "CID-SMILES.gz", ("smiles",), 1, 1_485_327_323),
    SourceFile("title", "CID-Title.gz", ("cmpdname",), 1, 1_888_458_761),
    SourceFile("iupac", "CID-IUPAC.gz", ("iupacname",), 1, 1_851_942_537),
    SourceFile("inchi", "CID-InChI-Key.gz", ("inchi", "inchikey"), 2, 7_361_682_757),
    SourceFile("synonyms", "CID-Synonym-filtered.gz", (), 1, 968_608_795),
)

_FILES = {source.key: source for source in SOURCE_FILES}


def selected_files(*, include_inchi: bool = True,
                   include_synonyms: bool = True) -> List[SourceFile]:
    """
    The source files a build with these options reads, in reading order.

    Args:
        include_inchi: Whether InChI and InChIKey come from PubChem's own file.
        include_synonyms: Whether the synonym list is included.

    Returns:
        The :class:`SourceFile` entries, a subset of :data:`SOURCE_FILES`.

    Example:
        >>> " ".join(f.key for f in selected_files(include_inchi=False))
        'identifiers date mass smiles title iupac synonyms'
    """
    skip = set()
    if not include_inchi:
        skip.add("inchi")
    if not include_synonyms:
        skip.add("synonyms")
    return [source for source in SOURCE_FILES if source.key not in skip]


# ──────────────────────────────────────────────────────────────────────────────
# Schema
# ──────────────────────────────────────────────────────────────────────────────

#: The ``compounds``, ``cas_numbers`` and ``synonyms`` tables are the ones the
#: Zenodo database has always had, so every query in
#: :class:`~provesid.PubChemID` runs unchanged; the eight descriptor columns
#: are gone and ``monoisotopicmass`` is new.
_SCHEMA = """
CREATE TABLE compounds (
    cid INTEGER PRIMARY KEY,
    cmpdname TEXT,
    mf TEXT,
    inchi TEXT,
    smiles TEXT,
    inchikey TEXT,
    iupacname TEXT,
    mw REAL,
    exactmass REAL,
    monoisotopicmass REAL,
    cidcdate TEXT
);
CREATE TABLE cas_numbers (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    cid INTEGER,
    cas TEXT,
    FOREIGN KEY (cid) REFERENCES compounds(cid)
);
CREATE TABLE synonyms (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    cid INTEGER,
    synonym TEXT,
    FOREIGN KEY (cid) REFERENCES compounds(cid)
);
CREATE TABLE xrefs (
    cid INTEGER,
    source TEXT,
    identifier TEXT,
    FOREIGN KEY (cid) REFERENCES compounds(cid)
);
CREATE TABLE provenance (
    key TEXT PRIMARY KEY,
    value TEXT
);
CREATE TABLE provenance_files (
    file TEXT PRIMARY KEY,
    url TEXT,
    md5 TEXT,
    bytes INTEGER,
    lines_read INTEGER,
    rows_kept INTEGER
);
"""

#: Built after the tables are filled, which is several times faster than
#: maintaining them row by row. The names are the Zenodo database's.
_INDEXES = """
CREATE INDEX idx_compounds_inchikey ON compounds(inchikey);
CREATE INDEX idx_compounds_inchi ON compounds(inchi);
CREATE INDEX idx_compounds_mf ON compounds(mf);
CREATE INDEX idx_cas_cas ON cas_numbers(cas);
CREATE INDEX idx_cas_cid ON cas_numbers(cid);
CREATE INDEX idx_synonyms_synonym ON synonyms(synonym);
CREATE INDEX idx_synonyms_cid ON synonyms(cid);
CREATE INDEX idx_xrefs_cid ON xrefs(cid);
CREATE INDEX idx_xrefs_identifier ON xrefs(identifier);
"""


# ──────────────────────────────────────────────────────────────────────────────
# The build
# ──────────────────────────────────────────────────────────────────────────────

def build_pubchem_id_db(
    db_path: Optional[str] = None,
    *,
    release: str = LATEST,
    include_inchi: bool = True,
    include_synonyms: bool = True,
    keep_downloads: bool = False,
    download_dir: Optional[str] = None,
    force: bool = False,
    progress: bool = True,
    base_url: Optional[str] = None,
    session: Optional[requests.Session] = None,
) -> str:
    """
    Build ``pubchem_id.db`` from PubChem's FTP files.

    Every compound carrying a valid CAS number in PubChem's curated identifier
    mapping is included, with its title, formula, molecular weight, exact and
    monoisotopic mass, isomeric SMILES, IUPAC name, InChI, InChIKey, creation
    date, synonyms, CAS numbers and cross-references to DSSTox, ChEBI, ChEMBL,
    EC and UNII.

    The source files are handled one at a time: each is downloaded (resumably,
    and checked against the MD5 PubChem publishes beside it), streamed through
    once to keep only the compounds in scope, and deleted. The database is
    built at ``db_path + '.tmp'`` and moved into place only when it is
    complete, so a failed or interrupted build never touches an existing
    database. A rerun starts the database over, but resumes an interrupted
    download and reuses any file ``keep_downloads=True`` left behind once its
    MD5 checks out.

    Measured costs, 2026-09-01 snapshot: 15.4 GB transferred (8.0 GB with
    ``include_inchi=False``), a 2.5 GB database, and at most the database plus
    the 7.4 GB InChI file on disk at once. Reading and writing take about 12
    minutes; the rest is the download, which is about 20 minutes at 12 MB/s
    and was five hours on a day PubChem served 0.85 MB/s. Memory stays under
    300 MB.

    Args:
        db_path: Where the database goes. Defaults to ``pubchem_id.db`` in the
            per-user dataset directory, where :class:`~provesid.PubChemID`
            looks for it.
        release: ``"latest"`` (default) for the newest monthly snapshot, a
            snapshot date such as ``"2026-09-01"``, or ``"current"`` for the
            rolling dump. A snapshot is frozen and makes the build
            reproducible; ``"current"`` is a day or so fresher but its files
            are regenerated in place and can change mid-build.
        include_inchi: Take InChI and InChIKey from PubChem's
            ``CID-InChI-Key.gz`` (default). False skips that 7.4 GB file and
            computes both from the SMILES with RDKit instead --- 7.4 GB less to
            download, but RDKit's InChI, not PubChem's, and a much longer
            build, since every structure is parsed.
        include_synonyms: Include the synonym table, which
            :meth:`~provesid.PubChemID.search_by_name` searches. False saves a
            1 GB download and most of the database's size.
        keep_downloads: Keep each source file after it has been read, in
            ``download_dir``. A later build of the same release then reuses
            them without downloading again.
        download_dir: Where source files are downloaded. Defaults to
            ``pubchem_ftp/<release>/`` beside the database.
        force: Replace a database already at ``db_path``.
        progress: Show progress bars.
        base_url: PubChem compound root. Defaults to :data:`FTP_ROOT`; a
            mirror, or a local server in a test, goes here.
        session: ``requests.Session`` to download through.

    Returns:
        The path of the finished database.

    Raises:
        FileExistsError: If a database is already at ``db_path`` and ``force``
            is False.
        ValueError: If ``release`` is not a release name.
        DownloadError: If a file cannot be downloaded or fails its checksum.

    Example:
        >>> from provesid.pubchem_ftp import build_pubchem_id_db
        >>> build_pubchem_id_db(release="2026-09-01")          # doctest: +SKIP
        '/home/me/.local/share/provesid/pubchem_id.db'
        >>> build_pubchem_id_db("/data/ids.db", include_synonyms=False,
        ...                     keep_downloads=True)            # doctest: +SKIP
        '/data/ids.db'
    """
    if db_path is None:
        db_path = os.path.join(user_dataset_path(), "pubchem_id.db")
    db_path = os.path.abspath(os.path.expanduser(db_path))
    if os.path.exists(db_path) and not force:
        raise FileExistsError(
            f"A database already exists at {db_path}. Pass force=True to rebuild it."
        )

    release = resolve_release(release, base_url=base_url, session=session)
    source_url = extras_url(release, base_url=base_url)
    if download_dir is None:
        download_dir = os.path.join(os.path.dirname(db_path), DOWNLOAD_DIRNAME, release)
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    os.makedirs(download_dir, exist_ok=True)
    files = selected_files(include_inchi=include_inchi, include_synonyms=include_synonyms)

    logger.info(
        "Building %s from PubChem release %s: %d files, about %.1f GB to download",
        db_path, release, len(files), sum(f.approx_bytes for f in files) / 1e9,
    )

    build = _Build(
        path=db_path + ".tmp",
        release=release,
        source_url=source_url,
        download_dir=download_dir,
        keep_downloads=keep_downloads,
        progress=progress,
        session=session,
    )
    try:
        build.run(files, derive_inchi=not include_inchi,
                  timestamp=_release_timestamp(release, base_url, session))
    except BaseException:
        build.discard()
        raise

    os.replace(build.path, db_path)
    if not keep_downloads:
        _remove_empty_directories(download_dir, stop_at=os.path.dirname(db_path))
    logger.info("PubChem ID database ready at %s (%.2f GB)",
                db_path, os.path.getsize(db_path) / 1e9)
    return db_path


def _remove_empty_directories(directory: str, *, stop_at: str) -> None:
    """
    Remove ``directory`` and its parents while they are empty, up to ``stop_at``.

    What a build without ``keep_downloads`` leaves is the database and nothing
    else, so the per-release download directory goes once its last file has
    been read --- and ``pubchem_ftp/`` with it, unless an earlier build kept
    files there. A directory the caller named that is not empty stays.
    """
    directory = os.path.abspath(directory)
    stop_at = os.path.abspath(stop_at)
    if os.path.commonpath([directory, stop_at]) != stop_at:
        # Somewhere else entirely: its parents are none of this build's business.
        if os.path.isdir(directory) and not os.listdir(directory):
            os.rmdir(directory)
        return
    while directory != stop_at and os.path.isdir(directory) and not os.listdir(directory):
        os.rmdir(directory)
        directory = os.path.dirname(directory)


def _release_timestamp(release: str, base_url: Optional[str],
                       session: Optional[requests.Session]) -> str:
    """
    The ``TIMESTAMP`` a monthly snapshot publishes, or ``""``.

    It records when PubChem took the snapshot, which is what a citation of the
    database wants. The rolling dump has no equivalent, and a snapshot without
    one is still usable, so a failure here is logged rather than raised.
    """
    if release == CURRENT:
        return ""
    url = f"{release_url(release, base_url=base_url)}/TIMESTAMP"
    getter = session.get if session is not None else requests.get
    try:
        response = getter(url, timeout=30)
        response.raise_for_status()
    except requests.RequestException as exc:
        logger.debug("No TIMESTAMP for release %s: %s", release, exc)
        return ""
    return response.text.strip()


class _Build:
    """
    One run of the builder: the temporary database and the files feeding it.

    A class rather than a function because the steps share a connection, the
    set of compounds in scope and the provenance being accumulated, and
    passing those through eight functions obscured more than it showed.
    """

    def __init__(self, *, path: str, release: str, source_url: str,
                 download_dir: str, keep_downloads: bool, progress: bool,
                 session: Optional[requests.Session]):
        self.path = path
        self.release = release
        self.source_url = source_url
        self.download_dir = download_dir
        self.keep_downloads = keep_downloads
        self.progress = progress
        self.session = session
        self.scope: Set[int] = set()
        self.stats: Dict[str, int] = {}
        self.connection: Optional[sqlite3.Connection] = None

    # ── lifecycle ──────────────────────────────────────────────────────────

    def run(self, files: List[SourceFile], *, derive_inchi: bool, timestamp: str) -> None:
        """Build the whole database at :attr:`path`."""
        self.discard()
        self.connection = sqlite3.connect(self.path)
        # Nothing reads this file until it is renamed into place, and a crash
        # discards it, so there is nothing for a journal to protect.
        self.connection.execute("PRAGMA journal_mode = OFF")
        self.connection.execute("PRAGMA synchronous = OFF")
        self.connection.executescript(_SCHEMA)

        for source in files:
            path, md5 = self._fetch(source)
            reader = {
                "identifiers": self._read_identifiers,
                "synonyms": self._read_synonyms,
            }.get(source.key, self._read_columns)
            lines, kept = reader(source, path)
            self.connection.execute(
                "INSERT INTO provenance_files VALUES (?, ?, ?, ?, ?, ?)",
                (source.filename, f"{self.source_url}/{source.filename}", md5,
                 os.path.getsize(path), lines, kept),
            )
            self.connection.commit()
            if not self.keep_downloads:
                os.remove(path)

        self._correct_isotopic_weights()
        if derive_inchi:
            self._derive_inchi()
        logger.info("Building indexes")
        self.connection.executescript(_INDEXES)
        self._write_provenance(files, derive_inchi=derive_inchi, timestamp=timestamp)
        self.connection.commit()
        self.connection.execute("ANALYZE")
        self.connection.close()
        self.connection = None

    def discard(self) -> None:
        """Close and delete the temporary database, if there is one."""
        if self.connection is not None:
            self.connection.close()
            self.connection = None
        if os.path.exists(self.path):
            os.remove(self.path)

    # ── downloading ────────────────────────────────────────────────────────

    def _fetch(self, source: SourceFile) -> Tuple[str, str]:
        """
        Put one source file on disk, verified, and return its path and MD5.

        A file already in the download directory is reused when its MD5 matches
        the one published beside it --- that is what ``keep_downloads=True`` is
        for --- and downloaded again otherwise.
        """
        url = f"{self.source_url}/{source.filename}"
        dest = os.path.join(self.download_dir, source.filename)
        expected = read_checksum(url + ".md5", session=self.session)
        if os.path.exists(dest):
            if md5_of_file(dest) == expected:
                logger.info("Reusing %s, already downloaded and verified", source.filename)
                return dest, expected
            logger.info("%s on disk does not match the published MD5; downloading again",
                        source.filename)
        download_file(url, dest, expected_md5=expected, progress=self.progress,
                      session=self.session, log=logger)
        return dest, expected

    # ── reading ────────────────────────────────────────────────────────────

    def _lines(self, source: SourceFile, path: str) -> Iterator[bytes]:
        """
        Every line of a gzipped source file, with a progress bar on its bytes.

        Lines are yielded as bytes: decoding only the lines that survive the
        scope filter is most of what keeps a pass over ``CID-InChI-Key.gz``
        --- 120 M lines, 1.4 M of them kept --- to minutes.
        """
        size = os.path.getsize(path)
        with open(path, "rb") as raw, gzip.GzipFile(fileobj=raw) as lines, \
                tqdm(total=size, unit="B", unit_scale=True, disable=not self.progress,
                     desc=f"Reading {source.filename}") as bar:
            for count, line in enumerate(lines, 1):
                yield line
                if count % _PROGRESS_EVERY == 0:
                    bar.update(raw.tell() - bar.n)
            bar.update(raw.tell() - bar.n)

    def _in_scope(self, line: bytes) -> Optional[int]:
        """The line's CID if that compound is in scope, else None."""
        tab = line.find(b"\t")
        if tab <= 0:
            return None
        head = line[:tab]
        if not head.isdigit():
            return None
        cid = int(head)
        return cid if cid in self.scope else None

    def _read_identifiers(self, source: SourceFile, path: str) -> Tuple[int, int]:
        """
        Decide the scope from the CAS rows, then keep the cross-references.

        Two passes over a 98 MB file, because the second needs the answer of
        the first: which compounds carry a valid CAS number is not known until
        the whole file has been read, and nothing guarantees it is sorted.
        Rows go straight to the database as they are read, and duplicates are
        removed there afterwards, rather than collected in Python first ---
        1.46 M (CID, CAS) pairs are a quarter of a gigabyte as tuples.
        """
        rejected = 0
        lines = 0
        batch: List[Tuple[int, str]] = []
        for line in self._lines(source, path):
            lines += 1
            fields = line.rstrip(b"\r\n").split(b"\t")
            if len(fields) != 3 or fields[2] != _CAS_TYPE or not fields[0].isdigit():
                continue
            cas = fields[1].decode("ascii", "replace").strip()
            if not check_CASRN(cas):
                rejected += 1
                continue
            cid = int(fields[0])
            self.scope.add(cid)
            batch.append((cid, cas))
            if len(batch) >= _BATCH_ROWS:
                self._insert("INSERT INTO cas_numbers (cid, cas) VALUES (?, ?)", batch)
        self._insert("INSERT INTO cas_numbers (cid, cas) VALUES (?, ?)", batch)
        self._deduplicate("cas_numbers", "cid, cas")

        self.stats["cas_rows_rejected"] = rejected
        logger.info("%s compounds carry a valid CAS number (%s CAS rows rejected by "
                    "the check digit)", f"{len(self.scope):,}", f"{rejected:,}")
        self.connection.executemany("INSERT INTO compounds (cid) VALUES (?)",
                                    ((cid,) for cid in sorted(self.scope)))

        wanted = {name.encode(): short for name, short in XREF_TYPES.items()}
        xref_batch: List[Tuple[int, str, str]] = []
        for line in self._lines(source, path):
            fields = line.rstrip(b"\r\n").split(b"\t")
            if len(fields) != 3 or fields[2] not in wanted:
                continue
            cid = self._in_scope(line)
            if cid is None:
                continue
            xref_batch.append((cid, wanted[fields[2]], fields[1].decode("utf-8", "replace")))
            if len(xref_batch) >= _BATCH_ROWS:
                self._insert("INSERT INTO xrefs VALUES (?, ?, ?)", xref_batch)
        self._insert("INSERT INTO xrefs VALUES (?, ?, ?)", xref_batch)
        self._deduplicate("xrefs", "cid, source, identifier")

        kept = sum(self.connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
                   for table in ("cas_numbers", "xrefs"))
        return lines, kept

    def _insert(self, statement: str, batch: List[Tuple]) -> None:
        """Write a batch of rows and empty it, so the caller can keep filling it."""
        if batch:
            self.connection.executemany(statement, batch)
            batch.clear()

    def _deduplicate(self, table: str, columns: str) -> None:
        """Keep the first of any rows that repeat ``columns``."""
        self.connection.execute(
            f"DELETE FROM {table} WHERE rowid NOT IN "
            f"(SELECT MIN(rowid) FROM {table} GROUP BY {columns})"
        )

    def _read_columns(self, source: SourceFile, path: str) -> Tuple[int, int]:
        """Fill one file's ``compounds`` columns for the compounds in scope."""
        width = source.fields
        assignments = ", ".join(f"{column} = ?" for column in source.columns)
        statement = f"UPDATE compounds SET {assignments} WHERE cid = ?"
        convert = _CONVERTERS.get(source.key, _as_text)

        batch: List[Tuple] = []
        lines = kept = 0
        for line in self._lines(source, path):
            lines += 1
            cid = self._in_scope(line)
            if cid is None:
                continue
            fields = line.rstrip(b"\r\n").split(b"\t", width)[1:]
            if len(fields) != width:
                continue
            batch.append(tuple(convert(fields)) + (cid,))
            kept += 1
            if len(batch) >= _BATCH_ROWS:
                self._insert(statement, batch)
        self._insert(statement, batch)
        return lines, kept

    def _read_synonyms(self, source: SourceFile, path: str) -> Tuple[int, int]:
        """Insert the synonyms of the compounds in scope, as PubChem orders them."""
        statement = "INSERT INTO synonyms (cid, synonym) VALUES (?, ?)"
        batch: List[Tuple[int, str]] = []
        lines = kept = 0
        for line in self._lines(source, path):
            lines += 1
            cid = self._in_scope(line)
            if cid is None:
                continue
            name = line.rstrip(b"\r\n").split(b"\t", 1)[1].decode("utf-8", "replace")
            if not name:
                continue
            batch.append((cid, name))
            kept += 1
            if len(batch) >= _BATCH_ROWS:
                self._insert(statement, batch)
        self._insert(statement, batch)
        return lines, kept

    # ── derived columns ────────────────────────────────────────────────────

    def _rows_in_pages(self, sql_filter: str, page: int = _BATCH_ROWS
                       ) -> Iterator[List[Tuple[int, str]]]:
        """
        ``(cid, smiles)`` rows matching a filter, a page at a time, by CID.

        Paging by key rather than holding one cursor open lets the caller
        update the same table between pages, and keeps memory to one page
        rather than 1.4 M rows.
        """
        last = -1
        while True:
            rows = self.connection.execute(
                f"SELECT cid, smiles FROM compounds WHERE cid > ? AND {sql_filter} "
                f"ORDER BY cid LIMIT ?", (last, page)).fetchall()
            if not rows:
                return
            yield rows
            last = rows[-1][0]

    def _correct_isotopic_weights(self) -> None:
        """
        Recompute ``mw`` from the SMILES for isotopically labelled compounds.

        ``mw`` was computed from the formula as ``CID-Mass.gz`` was read, and a
        formula cannot say that a hydrogen is deuterium: PubChem writes
        chloroform-*d* as ``CHCl3``. For the few thousand labelled compounds
        the SMILES is the only record of what the compound weighs.
        """
        corrected = 0
        for rows in self._rows_in_pages("smiles GLOB '*[[][0-9]*'"):
            updates = [(_isotopic_molecular_weight(smiles), cid) for cid, smiles in rows]
            self.connection.executemany("UPDATE compounds SET mw = ? WHERE cid = ?", updates)
            corrected += len(updates)
        self.stats["mw_from_isotopes"] = corrected

    def _derive_inchi(self) -> None:
        """Compute InChI and InChIKey from the stored SMILES with RDKit."""
        from rdkit import Chem, RDLogger

        RDLogger.DisableLog("rdApp.*")
        logger.info("Computing InChI and InChIKey from SMILES with RDKit")
        total = self.connection.execute(
            "SELECT COUNT(*) FROM compounds WHERE smiles IS NOT NULL").fetchone()[0]
        failed = 0
        with tqdm(total=total, desc="InChI from SMILES", disable=not self.progress) as bar:
            for rows in self._rows_in_pages("smiles IS NOT NULL"):
                updates = []
                for cid, smiles in rows:
                    molecule = Chem.MolFromSmiles(smiles)
                    inchi = Chem.MolToInchi(molecule) if molecule is not None else ""
                    if inchi:
                        updates.append((inchi, Chem.InchiToInchiKey(inchi), cid))
                    else:
                        failed += 1
                self.connection.executemany(
                    "UPDATE compounds SET inchi = ?, inchikey = ? WHERE cid = ?", updates)
                bar.update(len(rows))
        self.stats["inchi_failed"] = failed

    # ── provenance ─────────────────────────────────────────────────────────

    def _write_provenance(self, files: List[SourceFile], *, derive_inchi: bool,
                          timestamp: str) -> None:
        """Record what this database is, where it came from and how it was made."""
        from . import __version__

        def count(sql: str) -> int:
            return self.connection.execute(sql).fetchone()[0]

        entries = {
            "builder": "provesid.pubchem_ftp.build_pubchem_id_db",
            "builder_version": BUILDER_VERSION,
            "provesid_version": __version__,
            "built_at": datetime.datetime.now(datetime.timezone.utc)
                                .replace(microsecond=0).isoformat(),
            "release": self.release,
            "release_timestamp": timestamp,
            "source_url": self.source_url,
            "scope": "compounds with a valid CAS number in CID-Identifiers.tsv.gz",
            "files": ",".join(source.key for source in files),
            "inchi_source": "rdkit" if derive_inchi else "pubchem",
            "mw_source": "computed from mf with provesid.pubchem_ftp.ATOMIC_WEIGHTS; "
                         "from the SMILES for isotopically labelled compounds",
            "compounds": count("SELECT COUNT(*) FROM compounds"),
            "cas_numbers": count("SELECT COUNT(*) FROM cas_numbers"),
            "synonyms": count("SELECT COUNT(*) FROM synonyms"),
            "xrefs": count("SELECT COUNT(*) FROM xrefs"),
            "compounds_without_smiles": count(
                "SELECT COUNT(*) FROM compounds WHERE smiles IS NULL"),
            **self.stats,
        }
        self.connection.executemany("INSERT INTO provenance VALUES (?, ?)",
                                                    [(key, str(value)) for key, value in entries.items()])


def _as_text(fields: List[bytes]) -> List[str]:
    """Decode fields as UTF-8; IUPAC names and titles carry non-ASCII text."""
    return [field.decode("utf-8", "replace") for field in fields]


def _as_mass(fields: List[bytes]) -> List[object]:
    """
    ``CID-Mass.gz``: formula, monoisotopic mass, exact mass --- and, computed
    from the formula, the molecular weight no FTP file carries.
    """
    formula, monoisotopic, exact = fields
    formula = formula.decode("ascii", "replace")
    return [formula, molecular_weight(formula), _as_float(monoisotopic), _as_float(exact)]


def _as_float(field: bytes) -> Optional[float]:
    """A float, or None for an empty or malformed field."""
    try:
        return float(field)
    except ValueError:
        return None


_CONVERTERS: Dict[str, Callable[[List[bytes]], List[object]]] = {
    "mass": _as_mass,
}
