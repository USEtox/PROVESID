import pandas as pd
import os
from tqdm import tqdm

from provesid import PubChemView, get_property_table

df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'PubChem_compound_dissociation_constant.csv'))
cid_list = df["Compound_CID"].to_list()

# use PubChemView to extract Dissociation Constants values for the first 20 cids
print("Extracting Dissociation Constants values for the first 20 compounds...")

# Initialize PubChemView
pugview = PubChemView()

# Create a list to store all results
all_dc_data = []

# Process first 20 CIDs
for cid in tqdm(cid_list[:20]):
    
    try:
        # Get Dissociation Constants property table for this CID
        dc_table = get_property_table(cid, "Dissociation Constants")
        
        if not dc_table.empty:
            # Add this data to our collection
            all_dc_data.append(dc_table)
            
    except Exception as e:
        print(f"  Error processing CID {cid}: {str(e)}")

# Combine all results into a single DataFrame
if all_dc_data:
    combined_dc_df = pd.concat(all_dc_data, ignore_index=True)

    print(f"\n=== SUMMARY ===")
    print(f"Total Dissociation Constants entries found: {len(combined_dc_df)}")
    print(f"Unique compounds with Dissociation Constants data: {combined_dc_df['CID'].nunique()}")
    
    # Save results to CSV
    output_file = os.path.join(os.path.dirname(__file__), 'extracted_DC_data.csv')
    combined_dc_df.to_csv(output_file, index=False)
    print(f"Results saved to: {output_file}")
    
    # Display summary statistics
    print("\n=== Dissociation Constants Value Statistics ===")
    # Convert ExperimentalValue to numeric where possible
    numeric_dc = pd.to_numeric(combined_dc_df['ExperimentalValue'], errors='coerce')
    valid_dc = numeric_dc.dropna()
    
    if len(valid_dc) > 0:
        print(f"Valid numeric Dissociation Constants values: {len(valid_dc)}")
        print(f"Dissociation Constants range: {valid_dc.min():.2e} to {valid_dc.max():.2e}")
        print(f"Mean Dissociation Constants: {valid_dc.mean():.2e}")
        print(f"Median Dissociation Constants: {valid_dc.median():.2e}")
    
    # Display first few results
    print("\n=== Sample Results ===")
    print(combined_dc_df[['CID', 'ExperimentalValue', 'Unit', 'Temperature', 'Conditions']].head(10))
    
else:
    print("No Dissociation Constants data was successfully extracted.")
