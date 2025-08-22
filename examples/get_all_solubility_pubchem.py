import pandas as pd
import os
from tqdm import tqdm

from provesid import PubChemView, get_property_table

# Read the light CSV file (or create it if it doesn't exist)
light_csv_path = os.path.join(os.path.dirname(__file__), 'PubChem_compound_solubility_light.csv')

if not os.path.exists(light_csv_path):
    # Create the light CSV from the full CSV
    df = pd.read_csv(os.path.join(os.path.dirname(__file__), 'PubChem_compound_solubility.csv'))
    # only keep column Compound_CID Molecular_Weight Molecular_Formula XLogP
    df = df[["Compound_CID", "Molecular_Weight", "Molecular_Formula", "XLogP"]]
    # save it to the same directory
    df.to_csv(light_csv_path, index=False)
else:
    # Read the existing light CSV
    df = pd.read_csv(light_csv_path)

cid_list = df["Compound_CID"].to_list()

# use PubChemView to extract Solubility values for the first 20 cids
print("Extracting Solubility values for the first 20 compounds...")
print("Note: Solubility is reported both at a specific temperature and in a specific solvent.")
print("When no solvent is mentioned, the solubility is reported in water.")

# Initialize PubChemView
pugview = PubChemView()

# Create a list to store all results
all_solubility_data = []

# Process first 40 CIDs
for cid in tqdm(cid_list[:40]):
    
    try:
        # Get Solubility property table for this CID
        solubility_table = get_property_table(cid, "Solubility")
        
        if not solubility_table.empty:
            # Add this data to our collection
            all_solubility_data.append(solubility_table)
            
    except Exception as e:
        print(f"  Error processing CID {cid}: {str(e)}")

# Combine all results into a single DataFrame
if all_solubility_data:
    combined_solubility_df = pd.concat(all_solubility_data, ignore_index=True)

    print(f"\n=== SUMMARY ===")
    print(f"Total Solubility entries found: {len(combined_solubility_df)}")
    print(f"Unique compounds with Solubility data: {combined_solubility_df['CID'].nunique()}")
    
    # Save results to CSV
    output_file = os.path.join(os.path.dirname(__file__), 'extracted_Solubility_data.csv')
    combined_solubility_df.to_csv(output_file, index=False)
    print(f"Results saved to: {output_file}")
    
    # Display summary statistics
    print("\n=== Solubility Value Statistics ===")
    # Convert ExperimentalValue to numeric where possible
    numeric_solubility = pd.to_numeric(combined_solubility_df['ExperimentalValue'], errors='coerce')
    valid_solubility = numeric_solubility.dropna()
    
    if len(valid_solubility) > 0:
        print(f"Valid numeric Solubility values: {len(valid_solubility)}")
        print(f"Solubility range: {valid_solubility.min():.2e} to {valid_solubility.max():.2e}")
        print(f"Mean Solubility: {valid_solubility.mean():.2e}")
        print(f"Median Solubility: {valid_solubility.median():.2e}")
    
    # Display temperature and solvent information
    print("\n=== Temperature and Solvent Information ===")
    temp_count = combined_solubility_df['Temperature'].notna().sum()
    print(f"Entries with temperature information: {temp_count}")
    
    # Look for solvent information in conditions or string markup
    solvent_indicators = combined_solubility_df['StringWithMarkup'].str.contains(
        'water|ethanol|methanol|acetone|benzene|chloroform|ether|dmso|alcohol', 
        case=False, na=False
    ).sum()
    print(f"Entries with solvent indicators in text: {solvent_indicators}")
    
    # Display first few results
    print("\n=== Sample Results ===")
    display_cols = ['CID', 'ExperimentalValue', 'Unit', 'Temperature', 'StringWithMarkup']
    available_cols = [col for col in display_cols if col in combined_solubility_df.columns]
    print(combined_solubility_df[available_cols].head(10))
    
else:
    print("No Solubility data was successfully extracted.")
