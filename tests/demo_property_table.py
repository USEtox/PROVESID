#!/usr/bin/env python3
"""
Demonstration of the new get_property_table functionality
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from provesid import get_property_table, PubChemView

def demo_property_table():
    """Demonstrate the property table functionality"""
    print("=== PubChem PUG View Property Table Demo ===\n")
    
    # Demo 1: Aspirin melting point table
    print("1. Aspirin Melting Point Data Table:")
    mp_table = get_property_table(2244, "Melting Point")
    print(f"   Found {len(mp_table)} experimental entries\n")
    
    # Display first few rows with key columns
    for idx, (i, row) in enumerate(mp_table.head(3).iterrows()):
        print(f"   Entry {idx+1}:")
        print(f"     Original Text: {row['StringWithMarkup']}")
        print(f"     Parsed Value: {row['ExperimentalValue']} {row['Unit'] or ''}")
        print(f"     Reference: {row['FullReference'][:80]}...")
        print()
    
    # Demo 2: DMSO viscosity table
    print("2. DMSO Viscosity Data Table:")
    visc_table = get_property_table(679, "Viscosity")
    print(f"   Found {len(visc_table)} experimental entries\n")
    
    for idx, (i, row) in enumerate(visc_table.iterrows()):
        print(f"   Entry {idx+1}:")
        print(f"     Original: {row['StringWithMarkup']}")
        print(f"     Value: {row['ExperimentalValue']}")
        print(f"     Unit: {row['Unit']}")
        print(f"     Source: {row['FullReference'].split('|')[0] if '|' in row['FullReference'] else row['FullReference'][:50]}")
        print()
    
    # Demo 3: Export to CSV
    print("3. Exporting to CSV format:")
    pugview = PubChemView()
    density_table = pugview.get_property_table(702, "Density")  # ethanol
    
    if not density_table.empty:
        print("   Ethanol Density Data (CSV format):")
        print(density_table.to_csv(index=False)[:500] + "...")
    
    print("\n=== Demo completed! ===")
    print("\nThe get_property_table function provides:")
    print("• CID: PubChem Compound ID")
    print("• StringWithMarkup: Original experimental text")
    print("• ExperimentalValue: Parsed numerical value")
    print("• Unit: Extracted unit of measurement")
    print("• FullReference: Complete citation information")

if __name__ == "__main__":
    demo_property_table()
