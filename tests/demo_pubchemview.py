#!/usr/bin/env python3
"""
Demonstration script for the PubChem PUG View implementation
Shows practical usage examples for extracting experimental properties
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from provesid.pubchemview import (
    PubChemView, 
    PropertyData,
    get_experimental_property,
    get_property_values_only
)

def demo_pubchem_view():
    """Demonstrate practical usage of PubChem PUG View"""
    print("=== PubChem PUG View Experimental Properties Demo ===\n")
    
    # Initialize the PUG View client
    pugview = PubChemView()
    
    # Demo 1: Extract melting points for common compounds
    print("1. Melting points of common compounds:")
    compounds = {"aspirin": 2244, "caffeine": 2519, "ethanol": 702, "benzene": 241}
    for name, cid in compounds.items():
        try:
            mp_data = pugview.get_melting_point(cid)
            if mp_data:
                # Show the first melting point value
                print(f"   {name:10} → {mp_data[0].value}")
                if mp_data[0].unit:
                    print(f"              Unit: {mp_data[0].unit}")
            else:
                print(f"   {name:10} → No melting point data found")
        except Exception as e:
            print(f"   {name:10} → Error: {e}")
    
    # Demo 2: Detailed analysis of aspirin properties
    print("\n2. Detailed analysis of aspirin experimental properties:")
    try:
        all_props = pugview.extract_all_experimental_properties(2244)
        print(f"   Found {len(all_props)} experimental property types:")
        
        # Focus on key physical properties
        key_props = ["Melting Point", "Boiling Point", "Density", "Solubility", "Flash Point"]
        for prop in key_props:
            if prop in all_props:
                data_list = all_props[prop]
                print(f"\n   {prop}:")
                for i, data in enumerate(data_list[:2]):  # Show first 2 entries
                    print(f"     [{i+1}] {data.value}")
                    if data.reference_number:
                        print(f"         Reference: #{data.reference_number}")
                    if data.unit:
                        print(f"         Unit: {data.unit}")
                    if data.conditions:
                        print(f"         Conditions: {data.conditions}")
    except Exception as e:
        print(f"   Error: {e}")
    
    # Demo 3: Compare properties across compounds
    print("\n3. Density comparison across solvents:")
    solvents = {"water": 962, "ethanol": 702, "acetone": 180, "DMSO": 679}
    for name, cid in solvents.items():
        try:
            density_data = pugview.get_density(cid)
            if density_data:
                print(f"   {name:10} → {density_data[0].value}")
            else:
                print(f"   {name:10} → No density data")
        except Exception as e:
            print(f"   {name:10} → Error: {e}")
    
    # Demo 4: Using convenience functions
    print("\n4. Using convenience functions for quick access:")
    try:
        # Get just the values for boiling point of ethanol
        bp_values = get_property_values_only(702, "Boiling Point")
        print(f"   Ethanol boiling points ({len(bp_values)} values):")
        for value in bp_values[:3]:
            print(f"     - {value}")
        
        # Get single property with full details
        flash_point = get_experimental_property(2244, "Flash Point")
        if flash_point:
            print(f"\n   Aspirin flash point: {flash_point[0].value}")
            if flash_point[0].reference:
                print(f"   Reference: {flash_point[0].reference[:80]}...")
    except Exception as e:
        print(f"   Error: {e}")
    
    # Demo 5: Property summary and statistics
    print("\n5. Property summary example (vapor pressure of aspirin):")
    try:
        vp_summary = pugview.get_property_summary(2244, "Vapor Pressure")
        print(f"   Property: {vp_summary['property']}")
        print(f"   Number of values: {vp_summary['count']}")
        print(f"   Units found: {vp_summary['units']}")
        print(f"   Sample values:")
        for value in vp_summary['values'][:3]:
            print(f"     - {value}")
    except Exception as e:
        print(f"   Error: {e}")
    
    # Demo 6: Batch property extraction
    print("\n6. Batch extraction of multiple properties:")
    try:
        properties_of_interest = ["Melting Point", "Boiling Point", "Density", "Solubility"]
        batch_results = pugview.batch_extract_properties(2244, properties_of_interest)
        
        print(f"   Results for aspirin (CID 2244):")
        for prop_name, prop_data in batch_results.items():
            if prop_data:
                print(f"     {prop_name}: {len(prop_data)} entries")
                print(f"       First value: {prop_data[0].value}")
            else:
                print(f"     {prop_name}: No data available")
    except Exception as e:
        print(f"   Error: {e}")
    
    # Demo 7: Exploring available properties
    print("\n7. Available experimental properties for caffeine:")
    try:
        available = pugview.get_available_properties(2519)  # caffeine
        print(f"   Found {len(available)} property types:")
        for prop in available:
            print(f"     - {prop}")
    except Exception as e:
        print(f"   Error: {e}")
    
    print("\n=== Demo completed! ===")
    print("\nThe PubChem PUG View interface provides:")
    print("• Access to experimental properties not in standard PubChem API")
    print("• Structured extraction of values, units, conditions, and references")
    print("• Support for all major experimental property types")
    print("• Batch processing and convenience functions")
    print("• Robust error handling and rate limiting")
    print("• Easy export and summarization capabilities")

if __name__ == "__main__":
    demo_pubchem_view()
