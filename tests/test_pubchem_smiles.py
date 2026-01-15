"""
Test module for PubChemID SMILES-related methods.

Tests the following methods:
- get_by_smiles
- smiles_to_cid
- batch_smiles_to_cid
- get_by_smiles_batch
"""

import pytest
from provesid import PubChemID


@pytest.fixture
def db():
    """Create a PubChemID instance for testing."""
    return PubChemID()


def test_get_by_smiles(db):
    """Test get_by_smiles with known compound (Aspirin)."""
    # Aspirin SMILES: CC(=O)OC1=CC=CC=C1C(=O)O
    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    result = db.get_by_smiles(smiles)
    
    assert result is not None, "Should find compound for Aspirin SMILES"
    assert result['cid'] == 2244, "CID should be 2244 for Aspirin"
    assert 'cmpdname' in result
    assert 'inchikey' in result
    assert 'mf' in result
    print(f"✓ get_by_smiles found: {result['cmpdname']} (CID: {result['cid']})")


def test_get_by_smiles_not_found(db):
    """Test get_by_smiles with non-existent SMILES."""
    # Use an invalid/non-existent SMILES
    result = db.get_by_smiles("INVALID_SMILES_XYZ123")
    
    assert result is None, "Should return None for non-existent SMILES"
    print("✓ get_by_smiles correctly returns None for invalid SMILES")


def test_smiles_to_cid(db):
    """Test smiles_to_cid conversion."""
    # Aspirin SMILES
    smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    cid = db.smiles_to_cid(smiles)
    
    assert cid is not None, "Should find CID for Aspirin SMILES"
    assert cid == 2244, "CID should be 2244 for Aspirin"
    print(f"✓ smiles_to_cid: {smiles} -> CID {cid}")


def test_smiles_to_cid_not_found(db):
    """Test smiles_to_cid with non-existent SMILES."""
    cid = db.smiles_to_cid("INVALID_SMILES_XYZ123")
    
    assert cid is None, "Should return None for non-existent SMILES"
    print("✓ smiles_to_cid correctly returns None for invalid SMILES")


def test_batch_smiles_to_cid(db):
    """Test batch_smiles_to_cid with multiple SMILES."""
    smiles_list = [
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin, CID 2244
        "C",                          # Methane, CID 297
        "CCO"                         # Ethanol, CID 702
    ]
    
    results = db.batch_smiles_to_cid(smiles_list)
    
    assert len(results) == 3, "Should return 3 results"
    assert results["CC(=O)OC1=CC=CC=C1C(=O)O"] == 2244, "Aspirin CID should be 2244"
    assert results["C"] == 297, "Methane CID should be 297"
    assert results["CCO"] == 702, "Ethanol CID should be 702"
    
    print("✓ batch_smiles_to_cid results:")
    for smiles, cid in results.items():
        print(f"  {smiles[:30]:30s} -> CID {cid}")


def test_batch_smiles_to_cid_with_invalid(db):
    """Test batch_smiles_to_cid with mix of valid and invalid SMILES."""
    smiles_list = [
        "C",                          # Valid: Methane
        "INVALID_SMILES",             # Invalid
        "CCO"                         # Valid: Ethanol
    ]
    
    results = db.batch_smiles_to_cid(smiles_list)
    
    assert len(results) == 3, "Should return 3 results"
    assert results["C"] == 297, "Methane CID should be 297"
    assert results["INVALID_SMILES"] is None, "Invalid SMILES should return None"
    assert results["CCO"] == 702, "Ethanol CID should be 702"
    
    print("✓ batch_smiles_to_cid handles invalid SMILES correctly")


def test_get_by_smiles_batch(db):
    """Test get_by_smiles_batch with multiple SMILES."""
    smiles_list = [
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
        "C",                          # Methane
        "CCO"                         # Ethanol
    ]
    
    df = db.get_by_smiles_batch(smiles_list)
    
    assert len(df) == 3, "DataFrame should have 3 rows"
    assert 'cid' in df.columns, "Should have cid column"
    assert 'smiles' in df.columns, "Should have smiles column"
    assert 'cmpdname' in df.columns, "Should have cmpdname column"
    assert 'mf' in df.columns, "Should have molecular formula column"
    
    # Check specific values
    aspirin_row = df[df['cid'] == 2244]
    assert len(aspirin_row) == 1, "Should find Aspirin"
    
    methane_row = df[df['cid'] == 297]
    assert len(methane_row) == 1, "Should find Methane"
    
    ethanol_row = df[df['cid'] == 702]
    assert len(ethanol_row) == 1, "Should find Ethanol"
    
    print("✓ get_by_smiles_batch DataFrame:")
    print(df[['cid', 'smiles', 'cmpdname', 'mf', 'mw']].to_string(index=False))


def test_get_by_smiles_batch_empty(db):
    """Test get_by_smiles_batch with empty list."""
    df = db.get_by_smiles_batch([])
    
    assert len(df) == 0, "DataFrame should be empty"
    assert 'cid' in df.columns, "Should have correct columns even when empty"
    print("✓ get_by_smiles_batch handles empty list correctly")


def test_get_by_smiles_batch_all_invalid(db):
    """Test get_by_smiles_batch with all invalid SMILES."""
    smiles_list = ["INVALID1", "INVALID2", "INVALID3"]
    
    df = db.get_by_smiles_batch(smiles_list)
    
    assert len(df) == 0, "DataFrame should be empty for all invalid SMILES"
    assert 'cid' in df.columns, "Should have correct columns even when empty"
    print("✓ get_by_smiles_batch handles all invalid SMILES correctly")


def test_roundtrip_conversion(db):
    """Test roundtrip: CID -> SMILES -> CID."""
    original_cid = 2244  # Aspirin
    
    # Get SMILES from CID
    smiles = db.cid_to_smiles(original_cid)
    assert smiles is not None, "Should get SMILES for Aspirin"
    
    # Convert SMILES back to CID
    converted_cid = db.smiles_to_cid(smiles)
    assert converted_cid == original_cid, "Roundtrip conversion should yield same CID"
    
    print(f"✓ Roundtrip conversion: CID {original_cid} -> SMILES '{smiles}' -> CID {converted_cid}")


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v", "-s"])
