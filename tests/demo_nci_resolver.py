#!/usr/bin/env python3
"""
Demonstration script for the NCI Chemical Identifier Resolver
Shows practical usage examples for chemical identifier conversion
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from provesid.resolver import (
    NCIChemicalIdentifierResolver,
    nci_cas_to_mol,
    nci_name_to_smiles,
    nci_get_molecular_weight,
    nci_get_formula
)

def demo_nci_resolver():
    """Demonstrate practical usage of the NCI resolver"""
    print("=== NCI Chemical Identifier Resolver Demo ===\n")
    
    # Initialize the resolver
    resolver = NCIChemicalIdentifierResolver()
    
    # Demo 1: Convert common drug names to molecular structures
    print("1. Converting drug names to SMILES:")
    drugs = ['aspirin', 'caffeine', 'ibuprofen', 'acetaminophen']
    for drug in drugs:
        try:
            smiles = nci_name_to_smiles(drug)
            print(f"   {drug:15} → {smiles}")
        except Exception as e:
            print(f"   {drug:15} → Error: {e}")
    
    # Demo 2: Get molecular properties for common chemicals
    print("\n2. Molecular properties of common chemicals:")
    chemicals = ['water', 'ethanol', 'benzene', 'glucose']
    for chemical in chemicals:
        try:
            mw = nci_get_molecular_weight(chemical)
            formula = nci_get_formula(chemical)
            print(f"   {chemical:10} → Formula: {formula:10} MW: {mw}")
        except Exception as e:
            print(f"   {chemical:10} → Error: {e}")
    
    # Demo 3: Convert CAS numbers to molecular data
    print("\n3. Converting CAS numbers to molecular data:")
    cas_numbers = {
        '64-17-5': 'ethanol',
        '50-00-0': 'formaldehyde', 
        '71-43-2': 'benzene',
        '67-56-1': 'methanol'
    }
    
    for cas, name in cas_numbers.items():
        try:
            result = nci_cas_to_mol(cas)
            print(f"   {cas} ({name}):")
            print(f"      SMILES: {result.get('smiles', 'N/A')}")
            print(f"      Formula: {result.get('formula', 'N/A')}")
            print(f"      MW: {result.get('mw', 'N/A')}")
        except Exception as e:
            print(f"   {cas} ({name}) → Error: {e}")
    
    # Demo 4: Batch conversion of identifiers
    print("\n4. Batch conversion example:")
    compounds = ['aspirin', 'caffeine', 'morphine', 'nicotine']
    try:
        batch_results = resolver.batch_resolve(compounds, 'formula')
        print("   Chemical formulas:")
        for compound, formula in batch_results.items():
            if formula:
                print(f"      {compound:12} → {formula}")
            else:
                print(f"      {compound:12} → Not found")
    except Exception as e:
        print(f"   Batch conversion error: {e}")
    
    # Demo 5: Get comprehensive molecular data
    print("\n5. Comprehensive data for vitamin C:")
    try:
        mol_data = resolver.get_molecular_data('ascorbic acid')
        print(f"   SMILES: {mol_data.get('smiles', 'N/A')}")
        print(f"   InChI: {mol_data.get('stdinchi', 'N/A')[:60]}...")
        print(f"   Formula: {mol_data.get('formula', 'N/A')}")
        print(f"   MW: {mol_data.get('mw', 'N/A')}")
        names = mol_data.get('names', [])
        print(f"   Alternative names ({len(names)}): {', '.join(names[:3])}...")
    except Exception as e:
        print(f"   Error: {e}")
    
    # Demo 6: Chemical structure images
    print("\n6. Generating structure image URLs:")
    molecules = ['aspirin', 'caffeine', 'benzene']
    for molecule in molecules:
        try:
            image_url = resolver.get_image_url(molecule)
            print(f"   {molecule:10} → {image_url}")
        except Exception as e:
            print(f"   {molecule:10} → Error: {e}")
    
    print("\n=== Demo completed! ===")
    print("\nThe NCI Chemical Identifier Resolver provides:")
    print("• Conversion between different chemical identifiers")
    print("• Molecular properties (formula, molecular weight)")
    print("• Chemical names and synonyms")
    print("• Structure images")
    print("• Batch processing capabilities")
    print("• Support for CAS numbers, SMILES, InChI, and common names")

if __name__ == "__main__":
    demo_nci_resolver()
