"""
Build PubChem ID SQLite database from PubChem CAS CSV file.

This script processes the PubChem_CAS_202601.csv file to:
1. Extract CAS numbers, InChI, and InChIKey from the cmpdsynonym column
2. Create a SQLite database with identifiers and chemical properties
3. Build indexes for fast lookups

The resulting database is much smaller than the original CSV (~1.6M compounds).
"""

import sqlite3
import csv
import os
import sys
import re
from tqdm import tqdm

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from provesid.utils import data_path


# CAS number pattern: XXX-XX-X or longer variations
CAS_PATTERN = re.compile(r'\b\d{2,7}-\d{2}-\d\b')

# InChI pattern: starts with InChI=
INCHI_PATTERN = re.compile(r'\bInChI=1S?/[^\|]+')

# InChIKey pattern: 27 characters, format XXXXXXXXXXXXXX-YYYYYYYYYYYY-Z
INCHIKEY_PATTERN = re.compile(r'\b[A-Z]{14}-[A-Z]{10}-[A-Z]\b')


def extract_identifiers(synonym_string):
    """
    Extract CAS numbers, InChI, and InChIKey from synonym string.
    
    Args:
        synonym_string (str): Pipe-separated synonyms
        
    Returns:
        tuple: (cas_list, inchi_list, inchikey_list, remaining_synonyms)
    """
    if not synonym_string:
        return [], [], [], []
    
    synonyms = synonym_string.split('|')
    
    cas_numbers = []
    inchis = []
    inchikeys = []
    remaining = []
    
    for syn in synonyms:
        syn = syn.strip()
        if not syn:
            continue
        
        # Check for CAS
        cas_match = CAS_PATTERN.search(syn)
        if cas_match:
            cas_numbers.append(cas_match.group())
            continue
        
        # Check for InChIKey (check before InChI as InChI contains more text)
        inchikey_match = INCHIKEY_PATTERN.search(syn)
        if inchikey_match:
            inchikeys.append(inchikey_match.group())
            continue
        
        # Check for InChI
        inchi_match = INCHI_PATTERN.search(syn)
        if inchi_match:
            inchis.append(inchi_match.group())
            continue
        
        # Keep as regular synonym
        remaining.append(syn)
    
    return cas_numbers, inchis, inchikeys, remaining


def create_database(csv_path, db_path):
    """
    Create SQLite database from PubChem CAS CSV file.
    
    Args:
        csv_path (str): Path to PubChem_CAS_202601.csv
        db_path (str): Path to output SQLite database
    """
    print(f"Processing: {csv_path}")
    print(f"Output database: {db_path}")
    
    # Remove existing database
    if os.path.exists(db_path):
        os.remove(db_path)
    
    # Create database and tables
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Main compounds table with identifiers and properties
    cursor.execute("""
        CREATE TABLE compounds (
            cid INTEGER PRIMARY KEY,
            cmpdname TEXT,
            mf TEXT,
            inchi TEXT,
            smiles TEXT,
            inchikey TEXT,
            iupacname TEXT,
            mw REAL,
            polararea REAL,
            complexity REAL,
            xlogp REAL,
            heavycnt INTEGER,
            hbonddonor INTEGER,
            hbondacc INTEGER,
            rotbonds INTEGER,
            exactmass REAL,
            charge INTEGER,
            cidcdate TEXT
        )
    """)
    
    # CAS numbers table (one-to-many: one compound can have multiple CAS)
    cursor.execute("""
        CREATE TABLE cas_numbers (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            cid INTEGER,
            cas TEXT,
            FOREIGN KEY (cid) REFERENCES compounds(cid)
        )
    """)
    
    # Synonyms table (one-to-many)
    cursor.execute("""
        CREATE TABLE synonyms (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            cid INTEGER,
            synonym TEXT,
            FOREIGN KEY (cid) REFERENCES compounds(cid)
        )
    """)
    
    # Process CSV file
    print("Reading CSV file...")
    with open(csv_path, 'r', encoding='utf-8', errors='ignore') as f:
        reader = csv.DictReader(f)
        
        # Count total rows for progress bar
        print("Counting rows...")
        row_count = sum(1 for _ in open(csv_path, 'r', encoding='utf-8', errors='ignore')) - 1
        f.seek(0)
        next(reader)  # Skip header
        
        print(f"Processing {row_count:,} compounds...")
        
        batch_size = 10000
        compounds_batch = []
        cas_batch = []
        synonyms_batch = []
        
        for row in tqdm(reader, total=row_count, desc="Processing compounds"):
            try:
                # Extract identifiers from cmpdsynonym
                synonym_string = row.get('cmpdsynonym', '')
                cas_list, inchi_list, inchikey_list, remaining_synonyms = extract_identifiers(synonym_string)
                
                # Use InChI from synonym if not in main inchi field
                inchi = row.get('inchi', '').strip()
                if not inchi and inchi_list:
                    inchi = inchi_list[0]
                
                # Use InChIKey from synonym if not in main inchikey field
                inchikey = row.get('inchikey', '').strip()
                if not inchikey and inchikey_list:
                    inchikey = inchikey_list[0]
                
                # Add to compounds batch
                cid = int(row['cid'])
                compounds_batch.append((
                    cid,
                    row.get('cmpdname', '').strip(),
                    row.get('mf', '').strip(),
                    inchi,
                    row.get('smiles', '').strip(),
                    inchikey,
                    row.get('iupacname', '').strip(),
                    float(row['mw']) if row.get('mw') else None,
                    float(row['polararea']) if row.get('polararea') else None,
                    float(row['complexity']) if row.get('complexity') else None,
                    float(row['xlogp']) if row.get('xlogp') else None,
                    int(row['heavycnt']) if row.get('heavycnt') else None,
                    int(row['hbonddonor']) if row.get('hbonddonor') else None,
                    int(row['hbondacc']) if row.get('hbondacc') else None,
                    int(row['rotbonds']) if row.get('rotbonds') else None,
                    float(row['exactmass']) if row.get('exactmass') else None,
                    int(row['charge']) if row.get('charge') else None,
                    row.get('cidcdate', '').strip()
                ))
                
                # Add CAS numbers
                for cas in cas_list:
                    cas_batch.append((cid, cas))
                
                # Add remaining synonyms
                for syn in remaining_synonyms:
                    if syn:
                        synonyms_batch.append((cid, syn))
                
                # Insert batches when they reach batch_size
                if len(compounds_batch) >= batch_size:
                    cursor.executemany("""
                        INSERT INTO compounds VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
                    """, compounds_batch)
                    
                    if cas_batch:
                        cursor.executemany("""
                            INSERT INTO cas_numbers (cid, cas) VALUES (?,?)
                        """, cas_batch)
                    
                    if synonyms_batch:
                        cursor.executemany("""
                            INSERT INTO synonyms (cid, synonym) VALUES (?,?)
                        """, synonyms_batch)
                    
                    conn.commit()
                    compounds_batch = []
                    cas_batch = []
                    synonyms_batch = []
                
            except Exception as e:
                print(f"Error processing row for CID {row.get('cid', 'unknown')}: {e}")
                continue
        
        # Insert remaining batches
        if compounds_batch:
            cursor.executemany("""
                INSERT INTO compounds VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            """, compounds_batch)
        
        if cas_batch:
            cursor.executemany("""
                INSERT INTO cas_numbers (cid, cas) VALUES (?,?)
            """, cas_batch)
        
        if synonyms_batch:
            cursor.executemany("""
                INSERT INTO synonyms (cid, synonym) VALUES (?,?)
            """, synonyms_batch)
        
        conn.commit()
    
    # Create indexes for fast lookups
    print("Creating indexes...")
    cursor.execute("CREATE INDEX idx_compounds_inchikey ON compounds(inchikey)")
    cursor.execute("CREATE INDEX idx_compounds_inchi ON compounds(inchi)")
    cursor.execute("CREATE INDEX idx_compounds_mf ON compounds(mf)")
    cursor.execute("CREATE INDEX idx_cas_cas ON cas_numbers(cas)")
    cursor.execute("CREATE INDEX idx_cas_cid ON cas_numbers(cid)")
    cursor.execute("CREATE INDEX idx_synonyms_synonym ON synonyms(synonym)")
    cursor.execute("CREATE INDEX idx_synonyms_cid ON synonyms(cid)")
    conn.commit()
    
    # Get statistics
    cursor.execute("SELECT COUNT(*) FROM compounds")
    compound_count = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM cas_numbers")
    cas_count = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM synonyms")
    synonym_count = cursor.fetchone()[0]
    
    conn.close()
    
    print(f"\n✓ Database created successfully!")
    print(f"  - {compound_count:,} compounds")
    print(f"  - {cas_count:,} CAS numbers")
    print(f"  - {synonym_count:,} synonyms")
    print(f"  - Database size: {os.path.getsize(db_path) / (1024**3):.2f} GB")


if __name__ == '__main__':
    # Paths
    data_dir = data_path()
    csv_path = os.path.join(data_dir, 'PubChem_CAS_202601.csv')
    db_path = os.path.join(data_dir, 'pubchem_id.db')
    
    if not os.path.exists(csv_path):
        print(f"Error: CSV file not found at {csv_path}")
        sys.exit(1)
    
    create_database(csv_path, db_path)
    print("\nDone!")
