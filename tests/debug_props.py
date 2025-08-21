#!/usr/bin/env python3
"""
Debug test to check property names
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from provesid.pubchem import PubChemAPI, CompoundProperties

def debug_properties():
    """Debug property names"""
    api = PubChemAPI()
    
    properties = [
        CompoundProperties.MOLECULAR_FORMULA,
        CompoundProperties.MOLECULAR_WEIGHT,
        CompoundProperties.SMILES,
        CompoundProperties.INCHIKEY
    ]
    
    result = api.get_compound_properties(2244, properties)
    print("Full result:", result)
    
    if 'PropertyTable' in result:
        props = result['PropertyTable']['Properties'][0]
        print("Available properties:")
        for key in props.keys():
            print(f"  {key}: {props[key]}")

if __name__ == "__main__":
    debug_properties()
