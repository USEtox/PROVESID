"""Tests for REACHDossierID Excel-based identifier lookup and conversion."""

import os

import pytest

from provesid import REACHDossierID
from provesid.utils import data_path


REACH_PATH = os.path.join(data_path(), "reach_study_results-dossier_info_23-05-2023.xlsx")
REACH_EXISTS = os.path.exists(REACH_PATH)

pytestmark = pytest.mark.skipif(
    not REACH_EXISTS,
    reason="REACH dataset not found in data folder.",
)


@pytest.fixture(scope="module")
def reach_db():
    """Create one REACHDossierID instance for all tests in this module."""
    return REACHDossierID()


class TestREACHInitialization:
    """Test class initialization and basic dataset integrity."""

    def test_init_and_required_columns(self, reach_db):
        """The dataset loads and contains required identifier columns."""
        assert len(reach_db.df) > 0
        for col in reach_db.REQUIRED_COLUMNS:
            assert col in reach_db.df.columns

    def test_stats(self, reach_db):
        """Dataset statistics are returned with expected keys."""
        stats = reach_db.get_stats()
        assert "total_rows" in stats
        assert "rows_with_dossier_uuid" in stats
        assert stats["total_rows"] > 0


class TestREACHLookupsAndConversions:
    """Test lookup and conversion methods using live values from the dataset."""

    def test_get_by_dossier_uuid(self, reach_db):
        """Lookup by dossier UUID returns the corresponding record."""
        sample_uuid = (
            reach_db.df[reach_db.df[reach_db.COL_DOSSIER_UUID] != ""]
            .iloc[0][reach_db.COL_DOSSIER_UUID]
        )
        row = reach_db.get_by_dossier_uuid(sample_uuid)
        assert row is not None
        assert row[reach_db.COL_DOSSIER_UUID] == sample_uuid

    def test_cas_lookup_and_conversion(self, reach_db):
        """CAS lookup and CAS->UUID conversion are internally consistent."""
        cas_rows = reach_db.df[reach_db.df[reach_db.COL_CAS] != ""]
        if cas_rows.empty:
            pytest.skip("No CAS values available in the dataset")

        sample_cas = cas_rows.iloc[0][reach_db.COL_CAS]
        records = reach_db.get_by_cas(sample_cas)
        assert len(records) > 0
        assert all(record[reach_db.COL_CAS] == sample_cas for record in records)

        uuids = reach_db.cas_to_dossier_uuid(sample_cas)
        assert isinstance(uuids, list)
        assert len(uuids) > 0

    def test_inventory_lookup_and_conversion(self, reach_db):
        """Inventory lookup and EC->CAS conversion work for non-empty values."""
        ec_rows = reach_db.df[reach_db.df[reach_db.COL_EC] != ""]
        if ec_rows.empty:
            pytest.skip("No inventory number values available in the dataset")

        sample_ec = ec_rows.iloc[0][reach_db.COL_EC]
        records = reach_db.get_by_inventory_number(sample_ec)
        assert len(records) > 0
        assert all(record[reach_db.COL_EC] == sample_ec for record in records)

        cas_values = reach_db.inventory_number_to_cas(sample_ec)
        assert isinstance(cas_values, list)

    def test_name_search_and_conversion(self, reach_db):
        """Name search and name->UUID conversion return list outputs."""
        name_rows = reach_db.df[reach_db.df[reach_db.COL_NAME_SUBSTANCE] != ""]
        if name_rows.empty:
            pytest.skip("No NAME_SUBSTANCE values available in the dataset")

        sample_name = name_rows.iloc[0][reach_db.COL_NAME_SUBSTANCE]
        partial = sample_name[:20] if len(sample_name) > 20 else sample_name

        records = reach_db.get_by_name(partial, exact=False, limit=5)
        assert isinstance(records, list)
        assert len(records) >= 1

        uuid_values = reach_db.name_to_dossier_uuid(partial, exact=False, limit=10)
        assert isinstance(uuid_values, list)

    def test_uuid_to_fields(self, reach_db):
        """UUID conversion methods return expected scalar outputs."""
        sample_uuid = (
            reach_db.df[reach_db.df[reach_db.COL_DOSSIER_UUID] != ""]
            .iloc[0][reach_db.COL_DOSSIER_UUID]
        )

        cas_value = reach_db.dossier_uuid_to_cas(sample_uuid)
        ec_value = reach_db.dossier_uuid_to_inventory_number(sample_uuid)
        name_value = reach_db.dossier_uuid_to_name(sample_uuid)

        assert cas_value is None or isinstance(cas_value, str)
        assert ec_value is None or isinstance(ec_value, str)
        assert name_value is None or isinstance(name_value, str)
