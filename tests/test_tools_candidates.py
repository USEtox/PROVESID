"""
Tests for the candidate adapters in :mod:`provesid.tools`.

Everything here is offline: ZeroPM is a stub answering the two calls the
adapter makes, with the rows the real v0.0.4 database holds for ethanol.
"""

import pandas as pd
import pytest

from provesid.tools import candidate_from_zeropm_smiles

ETHANOL = ("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N")
ETHANOL_13C = ("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3/i2+1", "LFQSCWFLJHTTHZ-VQEHIDDOSA-N")


class _ZeroPMStub:
    """What ZeroPM answers for "CCO": 13C-ethanol's CAS sorts first."""

    TABLES = {
        "14742-23-5": [ETHANOL_13C, ETHANOL],
        "64-17-5": [ETHANOL],
    }

    def get_cas_from_smiles(self, smiles):
        return list(self.TABLES)

    def get_id_table_from_cas(self, cas):
        return pd.DataFrame(
            [{"cas": cas, "rank": rank, "inchi": inchi, "inchikey": key, "synonyms": ""}
             for rank, (inchi, key) in enumerate(self.TABLES[cas], start=1)]
        )


@pytest.mark.unit
def test_zeropm_structure_candidate_takes_the_query_structure_not_a_relative():
    candidate = candidate_from_zeropm_smiles("CCO", _ZeroPMStub())
    assert candidate["InChIKey"] == ETHANOL[1]
    assert candidate["InChI"] == ETHANOL[0]
    assert candidate["CAS_candidates"] == ["14742-23-5", "64-17-5"]
