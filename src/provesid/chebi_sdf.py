"""
Offline access to ChEBI through its SDF release: :class:`ChebiSDF`.

EBI publishes the whole of ChEBI as one gzipped SDF file (``chebi.sdf.gz``).
:class:`ChebiSDF` downloads it once into the per-user dataset directory,
indexes it, pickles the index beside the file, and answers lookups by ChEBI
ID, name, synonym, InChIKey, formula and cross-reference with no network.

The online ChEBI 2.0 REST client lives in :mod:`provesid.chebi`.

Examples:
    >>> from provesid import ChebiSDF
    >>> chebi = ChebiSDF()                              # doctest: +SKIP
    >>> chebi.get_compound_by_id("CHEBI:15377")['ChEBI NAME']   # doctest: +SKIP
    'water'
"""

import gzip
import logging
import os
import pickle
from typing import Any, Dict, List, Optional

import pandas as pd
from tqdm import tqdm

from .datasets import CHUNK_SIZE, download_file
from .utils import user_dataset_path


class ChebiSDF:
    """
    Parser for ChEBI SDF (Structure-Data File) for offline access to ChEBI data.

    This class provides efficient querying of the ChEBI SDF file containing
    ~190,000 compounds with structure data, chemical properties, synonyms,
    and cross-references to 80+ external databases.

    The class builds an index on first use for fast lookups. The index is
    cached to disk for faster subsequent initializations.

    If the SDF file is not found, it can be automatically downloaded from the
    ChEBI FTP server.

    Attributes:
        sdf_path (str): Path to ChEBI SDF file
        index (dict): In-memory index for fast lookups

    Examples:
        >>> chebi_sdf = ChebiSDF()
        >>> compound = chebi_sdf.get_compound_by_id("CHEBI:15377")
        >>> print(compound['ChEBI NAME'])
        water
    """

    # Default download URL for ChEBI SDF file
    DEFAULT_SDF_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/chebi.sdf.gz"

    def __init__(
        self,
        sdf_path: Optional[str] = None,
        rebuild_index: bool = False,
        auto_download: bool = True,
        sdf_url: Optional[str] = None,
        data_dir: Optional[str] = None,
        redownload: bool = False,
    ):
        """
        Initialize ChebiSDF parser.

        Args:
            sdf_path (str, optional): Path to ChEBI SDF file. If None, uses
                the persistent user dataset directory.
            rebuild_index (bool): If True, rebuild index even if cache exists (default: False)
            auto_download (bool): If True, automatically download SDF file if not found (default: True)
            sdf_url (str, optional): Custom URL to download the SDF from. If None, uses default.
            data_dir (str, optional): Directory to store the SDF when
                ``sdf_path`` is not provided.
            redownload (bool): If True, force re-download when
                ``auto_download`` is enabled.
        """

        if sdf_path is None:
            base_dir = data_dir or user_dataset_path()
            sdf_path = os.path.join(base_dir, 'chebi.sdf')

        self.sdf_path = os.path.abspath(os.path.expanduser(sdf_path))
        self.index_path = self.sdf_path + '.index.pkl'
        self.sdf_url = sdf_url or self.DEFAULT_SDF_URL
        self.logger = logging.getLogger(__name__)

        needs_download = redownload or not os.path.exists(self.sdf_path)

        # Check if SDF file exists, download if needed
        if needs_download:
            if auto_download:
                if redownload and os.path.exists(self.sdf_path):
                    self.logger.info(
                        "Forced ChEBI SDF redownload requested for: %s", self.sdf_path
                    )
                else:
                    self.logger.info(f"ChEBI SDF file not found at: {self.sdf_path}")
                self.logger.info("Downloading ChEBI SDF file automatically...")
                self.download_sdf(url=self.sdf_url, force=redownload)
            else:
                raise FileNotFoundError(
                    f"ChEBI SDF file not found at: {self.sdf_path}\n"
                    f"Please run ChebiSDF.download_sdf() or set auto_download=True\n"
                    f"Or download manually from: https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/"
                )

        if redownload:
            rebuild_index = True

        # Load or build index
        if rebuild_index or not os.path.exists(self.index_path):
            self.logger.info("Building index from SDF file...")
            self.index = self._build_index()
            self._save_index()
        else:
            self.logger.info("Loading cached index...")
            self.index = self._load_index()

    def download_sdf(self, url: Optional[str] = None, force: bool = False) -> str:
        """
        Download the ChEBI SDF file from the ChEBI FTP server.

        The file is downloaded as a gzip archive (~250 MB) and automatically
        extracted to the data directory (~868 MB uncompressed).

        The transfer is resumable: an interrupted download leaves a ``.part``
        file beside the archive and the next call continues from it. The gzip
        archive is kept until the SDF has been extracted and checked, so a
        failure during extraction does not cost another download.

        Args:
            url (str, optional): URL to download from. If None, uses default ChEBI FTP URL.
            force (bool): If True, download even if file already exists (default: False)

        Returns:
            str: Path to the downloaded and extracted SDF file

        Raises:
            FileExistsError: If the file already exists and force=False
            provesid.datasets.DownloadError: If the download could not be
                completed
            RuntimeError: If the extracted file is empty
            gzip.BadGzipFile: If the archive is damaged -- gzip's CRC and
                length trailer catch a truncated transfer for free

        Examples:
            >>> chebi_sdf = ChebiSDF(auto_download=False)   # doctest: +SKIP
            >>> chebi_sdf.download_sdf()                    # doctest: +SKIP
        """
        download_url = url or self.sdf_url

        # Check if file already exists
        if os.path.exists(self.sdf_path) and not force:
            raise FileExistsError(
                f"ChEBI SDF file already exists at: {self.sdf_path}\n"
                f"Use force=True to overwrite"
            )

        # The gzip is kept beside the SDF rather than in a temporary directory,
        # so an interrupted download's .part file is found and resumed next
        # time instead of fetching the archive again.
        gz_path = self.sdf_path + ".gz"
        temp_path = self.sdf_path + ".tmp"

        try:
            download_file(
                download_url,
                gz_path,
                description="ChEBI SDF",
                log=self.logger,
            )

            self.logger.info("Download complete. Extracting gzip archive...")

            # gzip's CRC and length trailer make this the truncation check the
            # plain downloads never had: a short archive fails here rather than
            # producing a half-written SDF.
            with gzip.open(gz_path, 'rb') as f_in:
                with open(temp_path, 'wb') as f_out:
                    with tqdm(total=os.path.getsize(gz_path), unit='B',
                              unit_scale=True, desc="Extracting ChEBI SDF") as pbar:
                        while True:
                            chunk = f_in.read(CHUNK_SIZE)
                            if not chunk:
                                break
                            f_out.write(chunk)
                            pbar.update(len(chunk))

            # Verify the file is valid before it replaces anything.
            with open(temp_path, 'r', encoding='utf-8', errors='ignore') as f:
                if not any(f.readline() for _ in range(5)):
                    raise RuntimeError("Downloaded file appears to be empty")

            os.replace(temp_path, self.sdf_path)
            os.remove(gz_path)

            self.logger.info("✓ ChEBI SDF file downloaded and extracted successfully")
            self.logger.info(f"✓ File location: {self.sdf_path}")
            return self.sdf_path

        except Exception:
            # The .gz is left where it is when it downloaded cleanly, so a
            # failure in extraction does not cost another 300 MB transfer.
            if os.path.exists(temp_path):
                os.remove(temp_path)
            raise

    def _build_index(self) -> Dict:
        """
        Build index from SDF file for fast lookups.

        The file is read in **binary** mode and each line decoded individually,
        so the recorded offsets are exact byte positions.  Text mode must not be
        used here: universal-newline translation collapses ``\\r\\n`` to ``\\n``,
        which under-counts one byte per CRLF line.  The ChEBI SDF mixes line
        endings (~59 000 CRLF lines in the 2026 release), so a text-mode offset
        drifts steadily and ``get_compound_by_id`` ends up seeking into a
        *neighbouring* record — silently returning another compound's data.

        Returns:
            dict: Index containing mappings for various query types, plus a
            ``_meta`` entry recording the SDF the index was built from (see
            :meth:`_index_meta`).
        """
        index = {
            'id_to_offset': {},           # ChEBI ID -> file offset
            'name_to_ids': {},            # lowercase name -> list of ChEBI IDs
            'inchikey_to_id': {},         # InChIKey -> ChEBI ID
            'inchi_to_id': {},            # InChI -> ChEBI ID
            'formula_to_ids': {},         # formula -> list of ChEBI IDs
            'cas_to_ids': {},             # CAS number -> list of ChEBI IDs
            'synonym_to_ids': {},         # lowercase synonym -> list of ChEBI IDs
            '_meta': self._index_meta(),  # provenance/validity of this index
        }

        # Parse SDF file and build index
        with open(self.sdf_path, 'rb') as f:
            file_offset = 0
            current_mol_offset = 0
            in_mol = False
            current_data = {}

            # Use tqdm for progress bar
            file_size = os.path.getsize(self.sdf_path)
            pbar = tqdm(total=file_size, unit='B', unit_scale=True, desc="Indexing ChEBI SDF")

            for raw_line in f:
                line = raw_line.decode('utf-8', errors='ignore')

                if not in_mol and line.strip():
                    # Start of a new molecule
                    in_mol = True
                    current_mol_offset = file_offset
                    current_data = {}

                # Check for property tags
                if line.startswith('> <'):
                    field_name = line.strip()[3:-1]  # Extract field name
                    raw_value = next(f, b'')
                    file_offset += len(raw_line) + len(raw_value)
                    current_data[field_name] = raw_value.decode('utf-8', errors='ignore').strip()
                    pbar.update(len(raw_line) + len(raw_value))
                    continue

                # Check for end of molecule
                if line.startswith('$$$$'):
                    in_mol = False

                    # Index this molecule
                    if 'ChEBI ID' in current_data:
                        chebi_id = current_data['ChEBI ID']
                        index['id_to_offset'][chebi_id] = current_mol_offset

                        # Index by name
                        if 'ChEBI NAME' in current_data:
                            name_lower = current_data['ChEBI NAME'].lower()
                            if name_lower not in index['name_to_ids']:
                                index['name_to_ids'][name_lower] = []
                            index['name_to_ids'][name_lower].append(chebi_id)

                        # Index by InChIKey
                        if 'INCHIKEY' in current_data:
                            index['inchikey_to_id'][current_data['INCHIKEY']] = chebi_id

                        # Index by InChI
                        if 'INCHI' in current_data:
                            index['inchi_to_id'][current_data['INCHI']] = chebi_id

                        # Index by formula
                        if 'FORMULA' in current_data:
                            formula = current_data['FORMULA']
                            if formula not in index['formula_to_ids']:
                                index['formula_to_ids'][formula] = []
                            index['formula_to_ids'][formula].append(chebi_id)

                        # Index by CAS
                        if 'CAS Registry Numbers' in current_data:
                            cas_numbers = current_data['CAS Registry Numbers'].split(';')
                            for cas in cas_numbers:
                                cas = cas.strip()
                                if cas:
                                    if cas not in index['cas_to_ids']:
                                        index['cas_to_ids'][cas] = []
                                    index['cas_to_ids'][cas].append(chebi_id)

                        # Index by synonyms
                        if 'SYNONYM' in current_data:
                            synonyms = current_data['SYNONYM'].split(';')
                            for syn in synonyms:
                                syn_lower = syn.strip().lower()
                                if syn_lower:
                                    if syn_lower not in index['synonym_to_ids']:
                                        index['synonym_to_ids'][syn_lower] = []
                                    index['synonym_to_ids'][syn_lower].append(chebi_id)

                file_offset += len(raw_line)
                pbar.update(len(raw_line))

            pbar.close()

        self.logger.info(f"Index built: {len(index['id_to_offset'])} compounds indexed")
        return index

    # Bumped whenever the offset convention or index layout changes, so a stale
    # cached index is rebuilt rather than silently mis-read.
    INDEX_FORMAT_VERSION = 2

    def _index_meta(self) -> Dict[str, Any]:
        """Describe the SDF this index belongs to.

        The index stores raw byte offsets into the SDF, so it is only valid for
        the exact file it was built from.  Recording the file's size and the
        index format version lets :meth:`_load_index` detect a stale cache — for
        instance one written by an older release whose offsets were computed in
        text mode and are silently wrong.

        Returns:
            dict: ``format_version`` and ``sdf_size`` (bytes) of the SDF file.
        """
        return {
            'format_version': self.INDEX_FORMAT_VERSION,
            'sdf_size': os.path.getsize(self.sdf_path),
        }

    def _save_index(self):
        """Save index to disk for faster subsequent loads."""
        try:
            with open(self.index_path, 'wb') as f:
                pickle.dump(self.index, f, protocol=pickle.HIGHEST_PROTOCOL)
            self.logger.info(f"Index saved to {self.index_path}")
        except Exception as e:
            self.logger.warning(f"Failed to save index: {e}")

    def _load_index(self) -> Dict:
        """Load index from disk, rebuilding it when it does not match the SDF.

        Returns:
            dict: A validated index for the current SDF file.  When the cached
            index is missing, unreadable, written by an older format version, or
            was built from a differently sized SDF, it is rebuilt and re-saved.
        """
        try:
            with open(self.index_path, 'rb') as f:
                index = pickle.load(f)
        except Exception as e:
            self.logger.warning(f"Failed to load index: {e}. Rebuilding...")
            index = self._build_index()
            self.index = index
            self._save_index()
            return index

        expected = self._index_meta()
        meta = index.get('_meta') if isinstance(index, dict) else None
        if meta != expected:
            self.logger.warning(
                "Cached ChEBI index does not match %s (cached meta: %s, expected: %s). "
                "Rebuilding — the stale index would return the wrong compound.",
                self.sdf_path, meta, expected,
            )
            index = self._build_index()
            self.index = index
            self._save_index()
            return index

        self.logger.info(f"Index loaded: {len(index['id_to_offset'])} compounds")
        return index

    def _read_mol_at_offset(self, offset: int) -> Dict[str, str]:
        """
        Read a molecule entry from the SDF file at a specific offset.

        Reads in binary and decodes per line so that ``offset`` is interpreted
        as the exact byte position recorded by :meth:`_build_index`.

        Args:
            offset (int): File offset where molecule starts

        Returns:
            dict: Molecule data including molfile and properties
        """
        data = {'molfile': ''}

        with open(self.sdf_path, 'rb') as f:
            f.seek(offset)
            in_molfile = True

            for raw_line in f:
                line = raw_line.decode('utf-8', errors='ignore')

                if in_molfile:
                    data['molfile'] += line
                    if line.startswith('M  END'):
                        in_molfile = False
                    continue

                # Parse property tags
                if line.startswith('> <'):
                    field_name = line.strip()[3:-1]
                    value_line = next(f, b'').decode('utf-8', errors='ignore').strip()
                    data[field_name] = value_line
                    continue

                # End of molecule
                if line.startswith('$$$$'):
                    break

        return data

    def get_compound_by_id(self, chebi_id: str) -> Optional[Dict[str, str]]:
        """
        Get compound data by ChEBI ID.

        Args:
            chebi_id (str): ChEBI ID (e.g., "CHEBI:15377" or "15377")

        Returns:
            dict: Compound data, or None if not found

        Examples:
            >>> chebi_sdf = ChebiSDF()
            >>> water = chebi_sdf.get_compound_by_id("CHEBI:15377")
            >>> print(water['ChEBI NAME'])
            water
        """
        # Normalize ID format
        if not chebi_id.startswith('CHEBI:'):
            chebi_id = f'CHEBI:{chebi_id}'

        offset = self.index['id_to_offset'].get(chebi_id)
        if offset is None:
            return None

        return self._read_mol_at_offset(offset)

    def search_by_name(self, name: str, exact: bool = True) -> List[Dict[str, str]]:
        """
        Search compounds by name.

        Args:
            name (str): Compound name to search for
            exact (bool): If True, exact match; if False, partial match (default: True)

        Returns:
            list: List of matching compound data

        Examples:
            >>> chebi_sdf = ChebiSDF()
            >>> results = chebi_sdf.search_by_name("water")
            >>> print(len(results))
            1
        """
        name_lower = name.lower()
        results = []

        if exact:
            chebi_ids = self.index['name_to_ids'].get(name_lower, [])
        else:
            # Partial match
            chebi_ids = []
            for indexed_name, ids in self.index['name_to_ids'].items():
                if name_lower in indexed_name:
                    chebi_ids.extend(ids)

        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)

        return results

    def search_by_synonym(self, synonym: str, exact: bool = True) -> List[Dict[str, str]]:
        """
        Search compounds by synonym.

        Matching ignores case.

        Args:
            synonym (str): Synonym to search for
            exact (bool): If True, exact match; if False, partial match (default: True)

        Returns:
            list: List of matching compound data. A partial match returns each
            compound once, in no particular order.

        Examples:
            >>> sdf = ChebiSDF()
            >>> [c["ChEBI ID"] for c in sdf.search_by_synonym("aspirin")]
            ['CHEBI:15365']
            >>> len(sdf.search_by_synonym("aspirin", exact=False)) > 1
            True
        """
        synonym_lower = synonym.lower()
        results = []

        if exact:
            chebi_ids = self.index['synonym_to_ids'].get(synonym_lower, [])
        else:
            # Partial match
            chebi_ids = []
            for indexed_syn, ids in self.index['synonym_to_ids'].items():
                if synonym_lower in indexed_syn:
                    chebi_ids.extend(ids)
            # Remove duplicates
            chebi_ids = list(set(chebi_ids))

        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)

        return results

    def search_by_inchikey(self, inchikey: str) -> Optional[Dict[str, str]]:
        """
        Search compound by InChIKey.

        Args:
            inchikey (str): InChIKey to search for

        Returns:
            dict: Compound data, or None if not found

        Examples:
            >>> ChebiSDF().search_by_inchikey("BSYNRYMUTXBXSQ-UHFFFAOYSA-N")["ChEBI NAME"]
            'acetylsalicylic acid'
        """
        chebi_id = self.index['inchikey_to_id'].get(inchikey)
        if chebi_id:
            return self.get_compound_by_id(chebi_id)
        return None

    def search_by_inchi(self, inchi: str) -> Optional[Dict[str, str]]:
        """
        Search compound by InChI.

        Args:
            inchi (str): InChI string to search for

        Returns:
            dict: Compound data, or None if not found. The match is on the
            exact string.

        Examples:
            >>> ChebiSDF().search_by_inchi("InChI=1S/H2O/h1H2")["ChEBI ID"]
            'CHEBI:15377'
        """
        chebi_id = self.index['inchi_to_id'].get(inchi)
        if chebi_id:
            return self.get_compound_by_id(chebi_id)
        return None

    def search_by_cas(self, cas: str) -> List[Dict[str, str]]:
        """
        Search compounds by CAS Registry Number.

        Args:
            cas (str): CAS Registry Number

        Returns:
            list: List of matching compound data. Only about 29 000 of
            ChEBI's ~192 000 entries carry a CAS number.

        Examples:
            >>> [(c["ChEBI ID"], c["ChEBI NAME"]) for c in ChebiSDF().search_by_cas("7732-18-5")]
            [('CHEBI:15377', 'water'), ('CHEBI:29375', 'diprotium oxide')]
        """
        chebi_ids = self.index['cas_to_ids'].get(cas, [])
        results = []

        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)

        return results

    def search_by_formula(self, formula: str) -> List[Dict[str, str]]:
        """
        Search compounds by molecular formula.

        Args:
            formula (str): Molecular formula (e.g., "H2O")

        Returns:
            list: List of matching compound data. The formula must be written
            as ChEBI writes it.

        Examples:
            >>> [c["ChEBI NAME"] for c in ChebiSDF().search_by_formula("H2O")]
            ['water', 'diprotium oxide']
        """
        chebi_ids = self.index['formula_to_ids'].get(formula, [])
        results = []

        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)

        return results

    def filter_by_star_rating(self, min_stars: int = 3) -> List[str]:
        """
        Get ChEBI IDs of compounds with minimum star rating.

        Args:
            min_stars (int): Minimum star rating (1-3, default: 3)

        Returns:
            list: List of ChEBI IDs matching criteria

        Note:
            Reads every record in the file, which takes several seconds.

        Examples:
            >>> three_star = ChebiSDF().filter_by_star_rating(3)  # doctest: +SKIP
            >>> len(three_star), three_star[:2]                   # doctest: +SKIP
            (52903, ['CHEBI:7', 'CHEBI:8'])
        """
        matching_ids = []

        for chebi_id in tqdm(self.index['id_to_offset'].keys(), desc="Filtering by star rating"):
            compound = self.get_compound_by_id(chebi_id)
            if compound and 'STAR' in compound:
                try:
                    star = int(compound['STAR'])
                    if star >= min_stars:
                        matching_ids.append(chebi_id)
                except ValueError:
                    continue

        return matching_ids

    def get_compounds_by_ids(self, chebi_ids: List[str]) -> List[Dict[str, str]]:
        """
        Get multiple compounds by ChEBI IDs.

        Args:
            chebi_ids (list): List of ChEBI IDs

        Returns:
            list: List of compound data dictionaries, in the order asked;
            IDs not in the file are left out

        Examples:
            >>> found = ChebiSDF().get_compounds_by_ids(["CHEBI:15377", "CHEBI:0", "CHEBI:15365"])
            >>> [c["ChEBI NAME"] for c in found]
            ['water', 'acetylsalicylic acid']
        """
        results = []
        for chebi_id in chebi_ids:
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                results.append(compound)
        return results

    def export_to_dataframe(self, chebi_ids: Optional[List[str]] = None,
                           fields: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Export compounds to pandas DataFrame.

        Args:
            chebi_ids (list, optional): List of ChEBI IDs. If None, exports all compounds.
            fields (list, optional): List of field names to include. If None, includes common fields.

        Returns:
            pd.DataFrame: DataFrame with compound data

        Examples:
            >>> chebi_sdf = ChebiSDF()
            >>> df = chebi_sdf.export_to_dataframe(["CHEBI:15377", "CHEBI:16236"])
            >>> print(df[['ChEBI ID', 'ChEBI NAME', 'FORMULA']])
                  ChEBI ID ChEBI NAME FORMULA
            0  CHEBI:15377      water     H2O
            1  CHEBI:16236    ethanol   C2H6O
        """
        if fields is None:
            fields = ['ChEBI ID', 'ChEBI NAME', 'STAR', 'FORMULA', 'MASS',
                     'SMILES', 'INCHI', 'INCHIKEY', 'CAS Registry Numbers']

        if chebi_ids is None:
            chebi_ids = list(self.index['id_to_offset'].keys())

        data = []
        for chebi_id in tqdm(chebi_ids, desc="Exporting to DataFrame"):
            compound = self.get_compound_by_id(chebi_id)
            if compound:
                row = {field: compound.get(field, None) for field in fields}
                data.append(row)

        return pd.DataFrame(data)

    def get_database_stats(self) -> Dict[str, int]:
        """
        Get statistics about the ChEBI SDF database.

        Returns:
            dict: Statistics including counts of various indexed fields:
            ``total_compounds``, and the number of distinct keys in each index
            --- ``compounds_with_inchikey``, ``compounds_with_inchi`` and
            ``compounds_with_cas`` count distinct InChIKeys, InChIs and CAS
            numbers, not compounds --- plus ``unique_formulas``,
            ``indexed_names`` and ``indexed_synonyms``.

        Examples:
            >>> stats = ChebiSDF().get_database_stats()
            >>> stats["total_compounds"] > stats["compounds_with_cas"]
            True
        """
        return {
            'total_compounds': len(self.index['id_to_offset']),
            'compounds_with_inchikey': len(self.index['inchikey_to_id']),
            'compounds_with_inchi': len(self.index['inchi_to_id']),
            'compounds_with_cas': len(self.index['cas_to_ids']),
            'unique_formulas': len(self.index['formula_to_ids']),
            'indexed_names': len(self.index['name_to_ids']),
            'indexed_synonyms': len(self.index['synonym_to_ids']),
        }
