import pandas as pd
import os
current_dir = os.path.dirname(os.path.abspath(__file__))

# read only the column cmpdsynonym into a df
df = pd.read_csv(os.path.join(current_dir, "PubChem_CAS_202601.csv"), usecols=["cmpdsynonym", "inchikey"])
df_10 = pd.read_csv(os.path.join(current_dir, "PubChem_CAS_202601.csv"), nrows=10)

# look for UNJJBGNPUUVVFQ-ZJUUUORDSA-N in the inchikey column and print the corresponding row
row = df[df["inchikey"].str.contains("UNJJBGNPUUVVFQ-ZJUUUORDSA-N", na=False)]
print(row)