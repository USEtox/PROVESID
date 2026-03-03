"""
REACH dossier identifier lookup and conversion utilities.

This module provides the `REACHDossierID` class for reading the REACH dossier
study results Excel sheet and performing fast lookups/conversions between:

- Dossier UUID
- Substance name
- CAS number
- EC inventory number
- IUPAC name
"""

import os
import re
import zipfile
import xml.etree.ElementTree as ET

from typing import Any, Dict, List, Optional

import pandas as pd

from .utils import data_path


class REACHDossierID:
    """
    Interface for REACH dossier identifier lookup and conversion.

    The class reads `reach_study_results-dossier_info_23-05-2023.xlsx` from the
    package data directory by default and exposes methods to search and convert
    identifiers across key columns in the dataset.
    """

    DEFAULT_FILE_NAME = "reach_study_results-dossier_info_23-05-2023.xlsx"
    DEFAULT_SHEET_NAME = "Data"

    COL_DOSSIER_UUID = "DOSSIER UUID"
    COL_NAME_SUBSTANCE = "NAME_SUBSTANCE"
    COL_CAS = "CAS_NUMBER_ref_sub"
    COL_EC = "NUMBER_IN_INVENTORY_ref_sub"
    COL_IUPAC = "IUPAC_NAME_ref_sub"

    REQUIRED_COLUMNS = {
        COL_DOSSIER_UUID,
        COL_NAME_SUBSTANCE,
        COL_CAS,
        COL_EC,
        COL_IUPAC,
    }

    def __init__(
        self,
        excel_path: Optional[str] = None,
        sheet_name: str = DEFAULT_SHEET_NAME,
    ):
        """
        Initialize the REACH dossier dataset.

        Args:
            excel_path (str, optional): Path to the REACH Excel file. If None,
                uses the default file in the package data folder.
            sheet_name (str, optional): Excel sheet to read. Defaults to `Data`.

        Raises:
            FileNotFoundError: If the Excel file does not exist.
            RuntimeError: If the workbook cannot be parsed or required columns are missing.
        """
        if excel_path is None:
            excel_path = os.path.join(data_path(), self.DEFAULT_FILE_NAME)

        if not os.path.exists(excel_path):
            raise FileNotFoundError(
                f"REACH dataset not found at: {excel_path}. "
                "Please ensure the Excel file is present in the data directory."
            )

        self.excel_path = excel_path
        self.sheet_name = sheet_name
        self.df = self._load_dataframe(excel_path=excel_path, sheet_name=sheet_name)
        self._verify_columns()
        self._normalize_dataframe()

    @staticmethod
    def _normalize_text(value: Any) -> str:
        """Convert any value to stripped string, returning empty string for NA-like values."""
        if value is None:
            return ""
        text = str(value).strip()
        if text.lower() in {"nan", "none", "na"}:
            return ""
        return text

    @staticmethod
    def _normalize_name(value: Any) -> str:
        """Normalize a name-like string for case-insensitive matching."""
        text = REACHDossierID._normalize_text(value).lower()
        return re.sub(r"\s+", " ", text).strip()

    def _load_dataframe(self, excel_path: str, sheet_name: str) -> pd.DataFrame:
        """
        Load the REACH Excel sheet into a pandas DataFrame.

        This method first attempts `pandas.read_excel`. If the runtime lacks the
        Excel engine dependency (e.g., `openpyxl`), it falls back to a built-in
        XLSX parser based on `zipfile` + XML.

        Args:
            excel_path (str): Path to the Excel file.
            sheet_name (str): Worksheet name.

        Returns:
            pd.DataFrame: Loaded data.

        Raises:
            RuntimeError: If the file cannot be parsed.
        """
        try:
            return pd.read_excel(excel_path, sheet_name=sheet_name)
        except ImportError:
            return self._read_xlsx_with_stdlib(excel_path=excel_path, sheet_name=sheet_name)
        except ValueError as exc:
            # Pandas may raise ValueError for missing sheet or engine issues.
            # Try stdlib parser first; if sheet is truly missing, it will raise clearly.
            try:
                return self._read_xlsx_with_stdlib(
                    excel_path=excel_path,
                    sheet_name=sheet_name,
                )
            except Exception as fallback_exc:
                raise RuntimeError(
                    f"Failed to parse REACH workbook with pandas and stdlib fallback: {exc}"
                ) from fallback_exc
        except Exception as exc:
            raise RuntimeError(f"Failed to load REACH dataset: {exc}") from exc

    def _read_xlsx_with_stdlib(self, excel_path: str, sheet_name: str) -> pd.DataFrame:
        """
        Read a simple XLSX worksheet using standard library XML parsing.

        Args:
            excel_path (str): Path to XLSX file.
            sheet_name (str): Target worksheet name.

        Returns:
            pd.DataFrame: Parsed worksheet data.

        Raises:
            RuntimeError: If parsing fails or sheet is not found.
        """
        ns = {
            "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
            "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
        }

        def column_index(cell_ref: str) -> int:
            match = re.match(r"([A-Z]+)", cell_ref or "")
            if not match:
                return 0
            index = 0
            for char in match.group(1):
                index = index * 26 + (ord(char) - 64)
            return index

        with zipfile.ZipFile(excel_path) as archive:
            workbook = ET.fromstring(archive.read("xl/workbook.xml"))
            rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
            rel_map = {item.attrib["Id"]: item.attrib["Target"] for item in rels}

            shared_strings: List[str] = []
            if "xl/sharedStrings.xml" in archive.namelist():
                shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
                for si in shared_root.findall("a:si", ns):
                    text = "".join(token.text or "" for token in si.findall(".//a:t", ns))
                    shared_strings.append(text)

            def read_cell(cell: ET.Element) -> str:
                cell_type = cell.attrib.get("t")
                value_node = cell.find("a:v", ns)
                if value_node is None:
                    inline_node = cell.find("a:is", ns)
                    if inline_node is not None:
                        return "".join(
                            token.text or "" for token in inline_node.findall(".//a:t", ns)
                        )
                    return ""

                raw = value_node.text or ""
                if cell_type == "s" and raw.isdigit():
                    idx = int(raw)
                    return shared_strings[idx] if 0 <= idx < len(shared_strings) else ""
                return raw

            sheets = workbook.find("a:sheets", ns)
            if sheets is None:
                raise RuntimeError("Invalid XLSX file: workbook has no sheets")

            target_sheet = None
            for sheet in sheets.findall("a:sheet", ns):
                if sheet.attrib.get("name") == sheet_name:
                    target_sheet = sheet
                    break

            if target_sheet is None:
                raise RuntimeError(f"Sheet '{sheet_name}' not found in workbook")

            rel_id = target_sheet.attrib.get(
                "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
            )
            if not rel_id or rel_id not in rel_map:
                raise RuntimeError(f"Cannot resolve XML target for sheet '{sheet_name}'")

            sheet_xml = rel_map[rel_id]
            if not sheet_xml.startswith("xl/"):
                sheet_xml = f"xl/{sheet_xml}"

            sheet_root = ET.fromstring(archive.read(sheet_xml))
            row_nodes = sheet_root.findall(".//a:sheetData/a:row", ns)
            if not row_nodes:
                return pd.DataFrame()

            parsed_rows: List[List[str]] = []
            max_cols = 0

            for row in row_nodes:
                values_by_col: Dict[int, str] = {}
                for cell in row.findall("a:c", ns):
                    col = column_index(cell.attrib.get("r", ""))
                    values_by_col[col] = read_cell(cell)

                max_cols = max(max_cols, max(values_by_col.keys(), default=0))
                parsed_rows.append(values_by_col)

            materialized: List[List[str]] = []
            for values_by_col in parsed_rows:
                row_values = [values_by_col.get(idx, "") for idx in range(1, max_cols + 1)]
                materialized.append(row_values)

            headers = materialized[0] if materialized else []
            records = materialized[1:] if len(materialized) > 1 else []
            return pd.DataFrame(records, columns=headers)

    def _verify_columns(self):
        """
        Validate that required identifier columns exist.

        Raises:
            RuntimeError: If any required columns are missing.
        """
        missing = self.REQUIRED_COLUMNS - set(self.df.columns)
        if missing:
            raise RuntimeError(f"REACH dataset missing required columns: {sorted(missing)}")

    def _normalize_dataframe(self):
        """Normalize core identifier columns to clean string values."""
        for column in self.REQUIRED_COLUMNS:
            self.df[column] = self.df[column].map(self._normalize_text)

    def _records_from_dataframe(self, frame: pd.DataFrame) -> List[Dict[str, str]]:
        """Convert a DataFrame slice to list of dict records with core columns only."""
        if frame.empty:
            return []
        subset = frame[
            [
                self.COL_DOSSIER_UUID,
                self.COL_NAME_SUBSTANCE,
                self.COL_CAS,
                self.COL_EC,
                self.COL_IUPAC,
            ]
        ]
        return subset.to_dict(orient="records")

    def _unique_nonempty(self, values: List[str]) -> List[str]:
        """Return deduplicated, non-empty strings preserving order."""
        output: List[str] = []
        seen = set()
        for value in values:
            clean = self._normalize_text(value)
            if clean and clean not in seen:
                output.append(clean)
                seen.add(clean)
        return output

    def get_stats(self) -> Dict[str, Any]:
        """
        Get summary statistics for the loaded REACH dataset.

        Returns:
            dict: Summary fields including row count and non-empty ID counts.
        """
        return {
            "total_rows": int(len(self.df)),
            "rows_with_dossier_uuid": int((self.df[self.COL_DOSSIER_UUID] != "").sum()),
            "rows_with_cas": int((self.df[self.COL_CAS] != "").sum()),
            "rows_with_inventory_number": int((self.df[self.COL_EC] != "").sum()),
            "rows_with_substance_name": int((self.df[self.COL_NAME_SUBSTANCE] != "").sum()),
            "rows_with_iupac_name": int((self.df[self.COL_IUPAC] != "").sum()),
        }

    def get_by_dossier_uuid(self, dossier_uuid: str) -> Optional[Dict[str, str]]:
        """
        Get one REACH record by dossier UUID.

        Args:
            dossier_uuid (str): Dossier UUID.

        Returns:
            dict | None: Matching record or None if not found.
        """
        key = self._normalize_text(dossier_uuid)
        if not key:
            return None
        frame = self.df[self.df[self.COL_DOSSIER_UUID] == key]
        records = self._records_from_dataframe(frame)
        return records[0] if records else None

    def get_by_cas(self, cas_number: str) -> List[Dict[str, str]]:
        """
        Get all REACH records matching a CAS number.

        Args:
            cas_number (str): CAS Registry Number.

        Returns:
            list[dict]: Matching records.
        """
        key = self._normalize_text(cas_number)
        if not key:
            return []
        frame = self.df[self.df[self.COL_CAS] == key]
        return self._records_from_dataframe(frame)

    def get_by_inventory_number(self, inventory_number: str) -> List[Dict[str, str]]:
        """
        Get all REACH records matching an EC inventory number.

        Args:
            inventory_number (str): EC inventory number.

        Returns:
            list[dict]: Matching records.
        """
        key = self._normalize_text(inventory_number)
        if not key:
            return []
        frame = self.df[self.df[self.COL_EC] == key]
        return self._records_from_dataframe(frame)

    def get_by_name(
        self,
        name: str,
        exact: bool = False,
        limit: int = 20,
    ) -> List[Dict[str, str]]:
        """
        Search records by substance name.

        Args:
            name (str): Substance name text to match.
            exact (bool, optional): If True, exact case-insensitive match.
                If False, partial contains match.
            limit (int, optional): Maximum number of results.

        Returns:
            list[dict]: Matching records.
        """
        key = self._normalize_name(name)
        if not key:
            return []

        normalized_column = self.df[self.COL_NAME_SUBSTANCE].map(self._normalize_name)
        if exact:
            frame = self.df[normalized_column == key]
        else:
            frame = self.df[normalized_column.str.contains(re.escape(key), regex=True)]

        if limit > 0:
            frame = frame.head(limit)
        return self._records_from_dataframe(frame)

    def get_by_iupac_name(
        self,
        iupac_name: str,
        exact: bool = False,
        limit: int = 20,
    ) -> List[Dict[str, str]]:
        """
        Search records by IUPAC name.

        Args:
            iupac_name (str): IUPAC name text to match.
            exact (bool, optional): If True, exact case-insensitive match.
                If False, partial contains match.
            limit (int, optional): Maximum number of results.

        Returns:
            list[dict]: Matching records.
        """
        key = self._normalize_name(iupac_name)
        if not key:
            return []

        normalized_column = self.df[self.COL_IUPAC].map(self._normalize_name)
        if exact:
            frame = self.df[normalized_column == key]
        else:
            frame = self.df[normalized_column.str.contains(re.escape(key), regex=True)]

        if limit > 0:
            frame = frame.head(limit)
        return self._records_from_dataframe(frame)

    def dossier_uuid_to_cas(self, dossier_uuid: str) -> Optional[str]:
        """
        Convert dossier UUID to CAS number.

        Args:
            dossier_uuid (str): Dossier UUID.

        Returns:
            str | None: CAS number if found and non-empty.
        """
        row = self.get_by_dossier_uuid(dossier_uuid)
        if not row:
            return None
        value = self._normalize_text(row.get(self.COL_CAS))
        return value if value else None

    def dossier_uuid_to_inventory_number(self, dossier_uuid: str) -> Optional[str]:
        """
        Convert dossier UUID to EC inventory number.

        Args:
            dossier_uuid (str): Dossier UUID.

        Returns:
            str | None: Inventory number if found and non-empty.
        """
        row = self.get_by_dossier_uuid(dossier_uuid)
        if not row:
            return None
        value = self._normalize_text(row.get(self.COL_EC))
        return value if value else None

    def dossier_uuid_to_name(self, dossier_uuid: str) -> Optional[str]:
        """
        Convert dossier UUID to substance name.

        Args:
            dossier_uuid (str): Dossier UUID.

        Returns:
            str | None: Substance name if found and non-empty.
        """
        row = self.get_by_dossier_uuid(dossier_uuid)
        if not row:
            return None
        value = self._normalize_text(row.get(self.COL_NAME_SUBSTANCE))
        return value if value else None

    def cas_to_dossier_uuid(self, cas_number: str) -> List[str]:
        """
        Convert CAS number to dossier UUID values.

        Args:
            cas_number (str): CAS number.

        Returns:
            list[str]: Dossier UUID values.
        """
        rows = self.get_by_cas(cas_number)
        return self._unique_nonempty([row.get(self.COL_DOSSIER_UUID, "") for row in rows])

    def cas_to_inventory_number(self, cas_number: str) -> List[str]:
        """
        Convert CAS number to EC inventory number values.

        Args:
            cas_number (str): CAS number.

        Returns:
            list[str]: Inventory number values.
        """
        rows = self.get_by_cas(cas_number)
        return self._unique_nonempty([row.get(self.COL_EC, "") for row in rows])

    def cas_to_name(self, cas_number: str) -> List[str]:
        """
        Convert CAS number to substance names.

        Args:
            cas_number (str): CAS number.

        Returns:
            list[str]: Substance names.
        """
        rows = self.get_by_cas(cas_number)
        return self._unique_nonempty([row.get(self.COL_NAME_SUBSTANCE, "") for row in rows])

    def inventory_number_to_cas(self, inventory_number: str) -> List[str]:
        """
        Convert EC inventory number to CAS number values.

        Args:
            inventory_number (str): EC inventory number.

        Returns:
            list[str]: CAS values.
        """
        rows = self.get_by_inventory_number(inventory_number)
        return self._unique_nonempty([row.get(self.COL_CAS, "") for row in rows])

    def inventory_number_to_dossier_uuid(self, inventory_number: str) -> List[str]:
        """
        Convert EC inventory number to dossier UUID values.

        Args:
            inventory_number (str): EC inventory number.

        Returns:
            list[str]: Dossier UUID values.
        """
        rows = self.get_by_inventory_number(inventory_number)
        return self._unique_nonempty([row.get(self.COL_DOSSIER_UUID, "") for row in rows])

    def name_to_cas(self, name: str, exact: bool = False, limit: int = 20) -> List[str]:
        """
        Convert substance name to CAS number values.

        Args:
            name (str): Substance name query.
            exact (bool, optional): Exact or partial matching behavior.
            limit (int, optional): Maximum records to inspect.

        Returns:
            list[str]: CAS values.
        """
        rows = self.get_by_name(name=name, exact=exact, limit=limit)
        return self._unique_nonempty([row.get(self.COL_CAS, "") for row in rows])

    def name_to_dossier_uuid(
        self,
        name: str,
        exact: bool = False,
        limit: int = 20,
    ) -> List[str]:
        """
        Convert substance name to dossier UUID values.

        Args:
            name (str): Substance name query.
            exact (bool, optional): Exact or partial matching behavior.
            limit (int, optional): Maximum records to inspect.

        Returns:
            list[str]: Dossier UUID values.
        """
        rows = self.get_by_name(name=name, exact=exact, limit=limit)
        return self._unique_nonempty([row.get(self.COL_DOSSIER_UUID, "") for row in rows])
