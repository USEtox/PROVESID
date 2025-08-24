import pandas as pd
import os
from tqdm import tqdm

from provesid import PubChemView, get_property_table

df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'PubChem_compound_VP.csv'))
cid_list = df["Compound_CID"].to_list()

# use PubChemView to extract Vapor Pressure values for the first 20 cids
print("Extracting Vapor Pressure values for the first 20 compounds...")

# Initialize PubChemView
pugview = PubChemView()

# Create a list to store all results
all_vp_data = []

# Process first 20 CIDs
for cid in tqdm(cid_list[:40]):
    
    try:
        # Get Vapor Pressure property table for this CID
        vp_table = get_property_table(cid, "Vapor Pressure")
        
        if not vp_table.empty:
            # Add this data to our collection
            all_vp_data.append(vp_table)
            
    except Exception as e:
        print(f"  Error processing CID {cid}: {str(e)}")

# Combine all results into a single DataFrame
if all_vp_data:
    combined_vp_df = pd.concat(all_vp_data, ignore_index=True)

    print(f"\n=== SUMMARY ===")
    print(f"Total Vapor Pressure entries found: {len(combined_vp_df)}")
    print(f"Unique compounds with Vapor Pressure data: {combined_vp_df['CID'].nunique()}")
    
    # Save results to CSV
    output_file = os.path.join(os.path.dirname(__file__), 'extracted_VP_data.csv')
    combined_vp_df.to_csv(output_file, index=False)
    print(f"Results saved to: {output_file}")
    
    # Display summary statistics
    print("\n=== Vapor Pressure Value Statistics ===")
    # Convert ExperimentalValue to numeric where possible
    numeric_vp = pd.to_numeric(combined_vp_df['ExperimentalValue'], errors='coerce')
    valid_vp = numeric_vp.dropna()
    
    if len(valid_vp) > 0:
        print(f"Valid numeric Vapor Pressure values: {len(valid_vp)}")
        print(f"Vapor Pressure range: {valid_vp.min():.2e} to {valid_vp.max():.2e}")
        print(f"Mean Vapor Pressure: {valid_vp.mean():.2e}")
        print(f"Median Vapor Pressure: {valid_vp.median():.2e}")
    
    # Display first few results
    print("\n=== Sample Results ===")
    print(combined_vp_df[['CID', 'ExperimentalValue', 'Unit']].head(10))
    
else:
    print("No Vapor Pressure data was successfully extracted.")
