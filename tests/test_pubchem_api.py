#!/usr/bin/env python3
"""
Test script for the PubChem API implementation
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from provesid.pubchem import PubChemAPI, CompoundProperties, PubChemNotFoundError

def test_basic_functionality():
    """Test basic API functionality"""
    print("Testing PubChem API implementation...")
    
    # Initialize API
    api = PubChemAPI()
    
    # Test 1: Get compound by CID (aspirin)
    print("\n1. Testing get_compound_by_cid(2244) - Aspirin")
    try:
        result = api.get_compound_by_cid(2244)
        print(f"   Success: Found compound data")
        if 'PC_Compounds' in result:
            print(f"   Compound ID: {result['PC_Compounds'][0]['id']['id']['cid']}")
    except Exception as e:
        print(f"   Error: {e}")
    
    # Test 2: Get compound properties
    print("\n2. Testing get_compound_properties() - Multiple properties")
    try:
        properties = [
            CompoundProperties.MOLECULAR_FORMULA,
            CompoundProperties.MOLECULAR_WEIGHT,
            CompoundProperties.SMILES,
            CompoundProperties.INCHIKEY
        ]
        result = api.get_compound_properties(2244, properties)
        print(f"   Success: Retrieved properties")
        if 'PropertyTable' in result:
            props = result['PropertyTable']['Properties'][0]
            print(f"   Molecular Formula: {props.get('MolecularFormula', 'N/A')}")
            print(f"   Molecular Weight: {props.get('MolecularWeight', 'N/A')}")
    except Exception as e:
        print(f"   Error: {e}")
    
    # Test 3: Search by name
    print("\n3. Testing get_cids_by_name('aspirin')")
    try:
        result = api.get_cids_by_name('aspirin')
        print(f"   Success: Found CIDs")
        if 'IdentifierList' in result:
            cids = result['IdentifierList']['CID'][:5]  # Show first 5
            print(f"   First 5 CIDs: {cids}")
    except Exception as e:
        print(f"   Error: {e}")
    
    # Test 4: Similarity search
    print("\n4. Testing similarity_search() - Ethanol")
    try:
        result = api.similarity_search('CCO', threshold=95, MaxRecords=10)
        print(f"   Success: Found similar compounds")
        if 'IdentifierList' in result:
            cids = result['IdentifierList']['CID']
            print(f"   Found {len(cids)} similar compounds")
    except Exception as e:
        print(f"   Error: {e}")
    
    # Test 5: Convenience method
    print("\n5. Testing search_compound() convenience method")
    try:
        result = api.search_compound('caffeine', 'name')
        print(f"   Success: {result['success']}")
        if result['success'] and 'PC_Compounds' in result['data']:
            cid = result['data']['PC_Compounds'][0]['id']['id']['cid']
            print(f"   Caffeine CID: {cid}")
    except Exception as e:
        print(f"   Error: {e}")
    
    # Test 6: Error handling
    print("\n6. Testing error handling with invalid CID")
    try:
        result = api.get_compound_by_cid(99999999)
        print(f"   Unexpected success: {result}")
    except PubChemNotFoundError:
        print(f"   Success: Correctly handled not found error")
    except Exception as e:
        print(f"   Error: {e}")
    
    print("\nAll tests completed!")

if __name__ == "__main__":
    test_basic_functionality()
