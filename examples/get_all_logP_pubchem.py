import pandas as pd
import os
from tqdm import tqdm

from provesid import PubChemView, get_property_table

df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'PubChem_compound_logP.csv'))
cid_list = df["Compound_CID"].to_list()

# use PubChemView to extract logP values for the first 20 cids
print("Extracting LogP values for the first 20 compounds...")

# Initialize PubChemView
pugview = PubChemView()

# Create a list to store all results
all_logp_data = []

# Process first 20 CIDs
for cid in tqdm(cid_list):
    
    try:
        # Get LogP property table for this CID
        logp_table = get_property_table(cid, "LogP")
        
        if not logp_table.empty:
            # Add this data to our collection
            all_logp_data.append(logp_table)
            
    except Exception as e:
        print(f"  Error processing CID {cid}: {str(e)}")

# Combine all results into a single DataFrame
if all_logp_data:
    combined_logp_df = pd.concat(all_logp_data, ignore_index=True)

    print(f"\n=== SUMMARY ===")
    print(f"Total LogP entries found: {len(combined_logp_df)}")
    print(f"Unique compounds with LogP data: {combined_logp_df['CID'].nunique()}")
    
    # Save results to CSV
    output_file = os.path.join(os.path.dirname(__file__), 'extracted_logP_data.csv')
    combined_logp_df.to_csv(output_file, index=False)
    print(f"Results saved to: {output_file}")
    
    # Display summary statistics
    print("\n=== LogP Value Statistics ===")
    # Convert ExperimentalValue to numeric where possible
    numeric_logp = pd.to_numeric(combined_logp_df['ExperimentalValue'], errors='coerce')
    valid_logp = numeric_logp.dropna()
    
    if len(valid_logp) > 0:
        print(f"Valid numeric LogP values: {len(valid_logp)}")
        print(f"LogP range: {valid_logp.min():.2f} to {valid_logp.max():.2f}")
        print(f"Mean LogP: {valid_logp.mean():.2f}")
        print(f"Median LogP: {valid_logp.median():.2f}")
    
    # Display first few results
    print("\n=== Sample Results ===")
    print(combined_logp_df[['CID', 'ExperimentalValue', 'Unit']].head(10))
    
else:
    print("No LogP data was successfully extracted.")