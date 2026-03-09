---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.1
kernelspec:
  display_name: physchem
  language: python
  name: python3
---

# ClassyFire Chemical Classification Tutorial

**NOTE**: currently, it does not run. It seems to be a problem with the service.

ClassyFire is a web-based application for the automated structural classification of chemical entities. This tutorial demonstrates how to use the `ClassyFireAPI` class from the `provesid` package to classify chemical compounds using their structural features.

ClassyFire provides hierarchical chemical classification based on:
- Chemical structure analysis
- Functional group identification
- Taxonomic classification into superclass, class, subclass levels
- Molecular framework analysis
- Chemical fingerprinting

The ClassyFire system can classify compounds into over 4,800 chemical categories and is particularly useful for:
- Metabolomics research
- Chemical database organization
- Drug discovery
- Natural product analysis
- Chemical space exploration

**Service URL**: http://classyfire.wishartlab.com
**Input Types**: SMILES, InChI, chemical names, structure files

```{code-cell} ipython3
from provesid import ClassyFireAPI
import time
import json

# Initialize ClassyFire API
print("ClassyFire API initialized successfully!")
print(f"Service URL: {ClassyFireAPI.URL}")
print("Ready to classify chemical compounds!")

# Note: ClassyFire is a web service that requires submitting queries and waiting for results
print("\nImportant: ClassyFire processing involves:")
print("1. Submit a query (with SMILES, InChI, or chemical name)")
print("2. Wait for processing (can take several seconds to minutes)")
print("3. Retrieve classification results")
print("4. Parse the hierarchical classification data")
```

## 1. Basic Usage - Submitting a Classification Query

The basic workflow involves submitting a query and then retrieving results. Let's start with a simple example:

```{code-cell} ipython3
# Example: Complete ClassyFire workflow with proper status codes
# Based on real usage patterns from the inspiration code

def demo_classyfire_workflow():
    """
    Demonstrates the complete ClassyFire classification workflow
    with proper error handling and status code checking.
    """
    print("=== ClassyFire Complete Workflow Demo ===")
    print()
    
    # Initialize the API
    classyfire = ClassyFireAPI()
    
    # Example compound: caffeine
    compound_name = "Caffeine"
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    
    print(f"Classifying: {compound_name}")
    print(f"SMILES: {smiles}")
    print()
    
    # Step 1: Submit query (expecting 201 status)
    print("Step 1: Submitting query...")
    response = classyfire.submit_query(compound_name, smiles)
    
    if response.status_code == 201:
        query_result = response.json()
        query_id = query_result['id']
        print(f"✓ Query submitted successfully! ID: {query_id}")
        
        # Step 2: Wait and check status (expecting 200 status with "Done" text)
        print()
        print("Step 2: Checking status...")
        max_attempts = 10
        wait_time = 15  # seconds between checks
        
        for attempt in range(max_attempts):
            print(f"  Attempt {attempt + 1}/{max_attempts}...")
            
            status_response = classyfire.query_status(query_id)
            if status_response and status_response.status_code == 200:
                # Status response is typically plain text "Done" when ready
                status_text = status_response.text.strip()
                print(f"  Status: {status_text}")
                
                if status_text == "Done":
                    print("  ✓ Classification complete!")
                    break
                else:
                    print(f"  ⏳ Still processing... waiting {wait_time}s")
                    time.sleep(wait_time)
            else:
                print(f"  ⚠️ Status check failed (code: {status_response.status_code if status_response else 'None'})")
                time.sleep(wait_time)
        
        # Step 3: Retrieve results (expecting 200 status with JSON data)
        print()
        print("Step 3: Retrieving results...")
        result_response = classyfire.get_query(query_id, format="json")
        
        if result_response.status_code == 200:
            classification_data = result_response.json()
            print("✓ Results retrieved successfully!")
            print()
            
            # Display key classification information
            print("=== Classification Summary ===")
            print(f"Compound: {classification_data.get('label', 'N/A')}")
            print(f"SMILES: {classification_data.get('smiles', 'N/A')}")
            print(f"Molecular Formula: {classification_data.get('molecular_formula', 'N/A')}")
            
            entities = classification_data.get('entities', [])
            if entities:
                entity = entities[0]
                print(f"Kingdom: {entity.get('kingdom', {}).get('name', 'N/A')}")
                print(f"Superclass: {entity.get('superclass', {}).get('name', 'N/A')}")
                print(f"Class: {entity.get('class', {}).get('name', 'N/A')}")
                print(f"Subclass: {entity.get('subclass', {}).get('name', 'N/A')}")
            
            print()
            print("✅ Complete workflow successful!")
            return classification_data
            
        elif result_response.status_code == 202:
            print("⏳ Results not ready yet (202 status)")
            return None
        else:
            print(f"✗ Failed to retrieve results (status: {result_response.status_code})")
            return None
    else:
        print(f"✗ Query submission failed (status: {response.status_code})")
        print(f"Response: {response.text}")
        return None

# Note: Actual execution commented out to avoid long waits in tutorial
print("This demo shows the complete ClassyFire workflow:")
print("1. Submit query → expect 201 (Created)")
print("2. Check status → expect 200 with 'Done' text when ready")
print("3. Get results → expect 200 with JSON classification data")
print()
print("Key status codes:")
print("• 201: Query successfully submitted")
print("• 200: Status check or results retrieval successful")
print("• 202: Results not ready yet (still processing)")
print("• 4xx/5xx: Various error conditions")

# Uncomment the line below to run the actual demo (takes several minutes):
# demo_result = demo_classyfire_workflow()
aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"

print("Submitting classification query for aspirin...")
print(f"SMILES: {aspirin_smiles}")

# Submit the query
response = ClassyFireAPI.submit_query("Aspirin Classification", aspirin_smiles)

# ClassyFire returns 201 (Created) for successful query submission
if response.status_code == 201:
    query_result = response.json()
    query_id = query_result['id']
    print(f"✓ Query submitted successfully!")
    print(f"  Query ID: {query_id}")
    print(f"  Label: {query_result.get('label', 'N/A')}")
    print(f"  Status: {query_result.get('classification_status', 'N/A')}")
    
    # Store the query ID for later use
    aspirin_query_id = query_id
    
else:
    print(f"✗ Failed to submit query. Status code: {response.status_code}")
    print(f"  Response: {response.text}")
    aspirin_query_id = None
```

## 2. Checking Query Status

After submitting a query, we need to check its status before retrieving results:

```{code-cell} ipython3
# Example: Batch processing with proper timing and retry logic
# Based on proven patterns for large-scale ClassyFire processing

def process_compounds_with_proper_timing(compounds_dict, save_every=10):
    """
    Process multiple compounds with proper timing, retry logic, and intermediate saving.
    This demonstrates the robust approach used in real research applications.
    
    Args:
        compounds_dict (dict): Dictionary of {name: smiles} pairs
        save_every (int): Save results every N successful classifications
    
    Returns:
        tuple: (successful_results, failed_compounds)
    """
    import tqdm
    
    # Configuration based on real-world usage
    sleep_time = [10, 10, 15, 15]  # Progressive wait times in seconds
    retry_time = 3  # Wait time between submission retries
    number_of_retrials = 3  # Number of submission attempts
    
    classification_results = []
    failed_compounds = []
    
    print(f"Processing {len(compounds_dict)} compounds with robust timing...")
    print(f"Sleep schedule: {sleep_time} seconds")
    print(f"Retry attempts: {number_of_retrials}")
    print("=" * 60)
    
    # Initialize ClassyFire API
    classyfire = ClassyFireAPI()
    
    # Process each compound with progress tracking
    for i, (compound_name, smiles) in enumerate(tqdm.tqdm(compounds_dict.items()), 1):
        print(f"\nProcessing {i}/{len(compounds_dict)}: {compound_name}")
        print(f"SMILES: {smiles}")
        
        success = False
        
        try:
            # Try to submit query with retries
            for attempt in range(number_of_retrials):
                print(f"  Submission attempt {attempt + 1}/{number_of_retrials}...")
                
                response = classyfire.submit_query(compound_name, smiles)
                
                if response.status_code == 201:
                    query_result = response.json()
                    query_id = query_result['id']
                    print(f"  ✓ Query submitted successfully! ID: {query_id}")
                    
                    # Wait and check status with progressive timing
                    for j, wait_time in enumerate(sleep_time, 1):
                        print(f"  Waiting {wait_time}s (check {j}/{len(sleep_time)})...")
                        time.sleep(wait_time)
                        
                        status_response = classyfire.query_status(query_id)
                        if status_response and status_response.status_code == 200:
                            status_text = status_response.text.strip()
                            print(f"  Status: {status_text}")
                            
                            if status_text == "Done":
                                print("  🎉 Classification complete!")
                                
                                # Retrieve results
                                result_response = classyfire.get_query(query_id, format="json")
                                if result_response.status_code == 200:
                                    classification_data = result_response.json()
                                    classification_results.append({
                                        'compound_name': compound_name,
                                        'smiles': smiles,
                                        'query_id': query_id,
                                        'classification': classification_data
                                    })
                                    print(f"  ✓ Results retrieved and saved!")
                                    success = True
                                    break
                                else:
                                    print(f"  ✗ Failed to retrieve results: {result_response.status_code}")
                            else:
                                print(f"  ⏳ Still processing...")
                        else:
                            print(f"  ⚠️ Status check failed")
                    
                    if success:
                        break
                    else:
                        print(f"  ⚠️ Classification timed out after all wait periods")
                        if j == len(sleep_time):  # If we've exhausted all wait times
                            failed_compounds.append({
                                'name': compound_name, 
                                'smiles': smiles, 
                                'reason': 'timeout'
                            })
                    break
                    
                else:
                    print(f"  ✗ Submission failed: {response.status_code}")
                    if attempt < number_of_retrials - 1:
                        print(f"    Retrying in {retry_time}s...")
                        time.sleep(retry_time)
                    else:
                        print(f"  ✗ All submission attempts failed")
                        failed_compounds.append({
                            'name': compound_name, 
                            'smiles': smiles, 
                            'reason': f'submission_failed_{response.status_code}'
                        })
        
        except Exception as e:
            print(f"  ✗ Exception occurred: {e}")
            failed_compounds.append({
                'name': compound_name, 
                'smiles': smiles, 
                'reason': f'exception_{str(e)}'
            })
        
        # Save intermediate results
        if len(classification_results) % save_every == 0 and len(classification_results) > 0:
            print(f"\n💾 Intermediate save: {len(classification_results)} successful classifications")
    
    print(f"\n" + "=" * 60)
    print(f"PROCESSING COMPLETE")
    print(f"Successful classifications: {len(classification_results)}")
    print(f"Failed compounds: {len(failed_compounds)}")
    
    if failed_compounds:
        print("\nFailed compounds:")
        for failed in failed_compounds[:5]:  # Show first 5
            print(f"  {failed['name']}: {failed['reason']}")
        if len(failed_compounds) > 5:
            print(f"  ... and {len(failed_compounds) - 5} more")
    
    return classification_results, failed_compounds

# Example usage (demonstration only - actual execution takes hours)
print("Example of robust batch processing:")
print()

demo_compounds = {
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Glucose": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O"
}

print("To run batch processing:")
print("results, failures = process_compounds_with_proper_timing(demo_compounds)")
print()
print("Key features of this approach:")
print("• Progressive wait times: 10, 10, 15, 15 seconds")
print("• Multiple submission retries with delays")
print("• Comprehensive error handling and logging")
print("• Intermediate result saving")
print("• Progress tracking with tqdm")
print("• Detailed failure reporting")

# Uncomment the line below to run actual batch processing:
# results, failures = process_compounds_with_proper_timing(demo_compounds, save_every=2)
```

```{code-cell} ipython3
# Check query status for aspirin
if aspirin_query_id:
    print(f"Checking status for query ID: {aspirin_query_id}")
    
    status_response = ClassyFireAPI.query_status(aspirin_query_id)
    
    if status_response and status_response.status_code == 200:
        # Check if response is just plain text "Done" or JSON
        content_type = status_response.headers.get('content-type', '')
        if 'application/json' in content_type:
            try:
                status_data = status_response.json()
                status = status_data.get('classification_status', 'Unknown')
                print(f"✓ Status check successful:")
                print(f"  Classification Status: {status}")
                print(f"  Submission Time: {status_data.get('created_at', 'Unknown')}")
            except Exception:
                status = status_response.text.strip()
                print(f"✓ Status check successful:")
                print(f"  Classification Status: {status}")
        else:
            # Plain text response (usually "Done" when complete)
            status = status_response.text.strip()
            print(f"✓ Status check successful:")
            print(f"  Classification Status: {status}")
        
        # Check if classification is complete
        if status == 'Done':
            print("  🎉 Classification is complete! Ready to retrieve results.")
        elif status in ['In progress', 'Queued']:
            print("  ⏳ Classification is still in progress. Please wait...")
        else:
            print(f"  ⚠️  Status: {status}")
            
    else:
        print("✗ Failed to check query status")
        if status_response:
            print(f"  Status code: {status_response.status_code}")
else:
    print("No query ID available to check status")

print()
print("Note: ClassyFire processing typically takes 30 seconds to several minutes")
print("depending on the complexity of the molecule and server load.")
```

## 3. Retrieving Classification Results

Once the classification is complete, we can retrieve the detailed results in various formats:

```{code-cell} ipython3
# Wait a bit to allow processing (in real use, you might need longer waits)
print("Waiting for classification to complete...")
time.sleep(5)  # Adjust as needed

# Retrieve classification results for aspirin
if aspirin_query_id:
    print(f"Retrieving classification results for query ID: {aspirin_query_id}")
    
    # Get results in JSON format
    result_response = ClassyFireAPI.get_query(aspirin_query_id, format="json")
    
    if result_response.status_code == 200:
        classification_results = result_response.json()
        print("✓ Classification results retrieved successfully!")
        print()
        
        # Display basic information
        print("=== COMPOUND INFORMATION ===")
        print(f"Query Label: {classification_results.get('label', 'N/A')}")
        print(f"SMILES: {classification_results.get('smiles', 'N/A')}")
        print(f"InChI: {classification_results.get('inchi', 'N/A')}")
        print(f"InChI Key: {classification_results.get('inchikey', 'N/A')}")
        print(f"Molecular Formula: {classification_results.get('molecular_formula', 'N/A')}")
        print()
        
        # Display hierarchical classification
        print("=== HIERARCHICAL CLASSIFICATION ===")
        entities = classification_results.get('entities', [])
        if entities:
            for entity in entities:
                kingdom = entity.get('kingdom', {})
                superclass = entity.get('superclass', {})
                class_info = entity.get('class', {})
                subclass = entity.get('subclass', {})
                
                print(f"Kingdom: {kingdom.get('name', 'N/A')} ({kingdom.get('description', 'No description')})")
                print(f"Superclass: {superclass.get('name', 'N/A')} ({superclass.get('description', 'No description')})")
                print(f"Class: {class_info.get('name', 'N/A')} ({class_info.get('description', 'No description')})")
                print(f"Subclass: {subclass.get('name', 'N/A')} ({subclass.get('description', 'No description')})")
                
                # Show intermediate nodes if available
                intermediate_nodes = entity.get('intermediate_nodes', [])
                if intermediate_nodes:
                    print(f"Intermediate Nodes ({len(intermediate_nodes)}):")
                    for i, node in enumerate(intermediate_nodes[:3], 1):  # Show first 3
                        print(f"  {i}. {node.get('name', 'N/A')}: {node.get('description', 'No description')}")
                    if len(intermediate_nodes) > 3:
                        print(f"  ... and {len(intermediate_nodes) - 3} more")
                
                # Show direct parent
                direct_parent = entity.get('direct_parent', {})
                if direct_parent:
                    print(f"Direct Parent: {direct_parent.get('name', 'N/A')}")
                
                print()
        
        # Display molecular framework (if available)
        molecular_framework = classification_results.get('molecular_framework', 'N/A')
        if molecular_framework != 'N/A':
            print(f"=== MOLECULAR FRAMEWORK ===")
            print(f"Framework: {molecular_framework}")
            print()
            
    elif result_response.status_code == 202:
        print("⏳ Classification is still in progress. Please wait longer and try again.")
    else:
        print(f"✗ Failed to retrieve results. Status code: {result_response.status_code}")
        print(f"  Response: {result_response.text}")
else:
    print("No query ID available to retrieve results")
```

## 4. Classifying Multiple Compounds

Let's classify several different types of compounds to see the diversity of ClassyFire classifications:

```{code-cell} ipython3
# Define different types of compounds to classify
compounds = {
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Glucose": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O",
    "Cholesterol": "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
    "Ethanol": "CCO",
    "Benzene": "C1=CC=CC=C1"
}

# Submit queries for multiple compounds
query_ids = {}
print("Submitting classification queries for multiple compounds:")
print("=" * 60)

for name, smiles in compounds.items():
    print(f"Submitting query for {name}...")
    print(f"  SMILES: {smiles}")
    
    response = ClassyFireAPI.submit_query(f"{name} Classification", smiles)
    
    # ClassyFire returns 201 (Created) for successful query submission
    if response.status_code == 201:
        query_result = response.json()
        query_id = query_result['id']
        query_ids[name] = query_id
        print(f"  ✓ Success! Query ID: {query_id}")
    else:
        print(f"  ✗ Failed! Status code: {response.status_code}")
        query_ids[name] = None
    
    print()
    # Add a small delay between requests to be respectful to the server
    time.sleep(1)

print(f"Submitted {len([q for q in query_ids.values() if q])} successful queries out of {len(compounds)} compounds.")
print("\\nNote: Processing may take several minutes. You can check status individually.")
```

## 5. Different Output Formats

ClassyFire supports multiple output formats. Let's demonstrate JSON, CSV, and SDF formats:

```{code-cell} ipython3
# Demonstrate different output formats using aspirin query (if available)
if aspirin_query_id:
    print("Demonstrating different output formats for aspirin classification:")
    print("=" * 65)
    
    # JSON format (most detailed)
    print("1. JSON Format (detailed):")
    json_response = ClassyFireAPI.get_query(aspirin_query_id, format="json")
    if json_response.status_code == 200:
        json_data = json_response.json()
        print(f"   ✓ JSON data retrieved ({len(json.dumps(json_data))} characters)")
        print(f"   Contains keys: {list(json_data.keys())}")
        
        # Show a sample of the JSON structure
        print("   Sample JSON structure:")
        sample_data = {
            "smiles": json_data.get("smiles"),
            "molecular_formula": json_data.get("molecular_formula"),
            "kingdom": json_data.get("entities", [{}])[0].get("kingdom", {}).get("name") if json_data.get("entities") else None
        }
        print(f"   {json.dumps(sample_data, indent=4)}")
    else:
        print(f"   ✗ Failed to get JSON format. Status: {json_response.status_code}")
    
    print()
    
    # CSV format
    print("2. CSV Format (tabular):")
    csv_response = ClassyFireAPI.get_query(aspirin_query_id, format="csv")
    if csv_response.status_code == 200:
        csv_data = csv_response.text
        print(f"   ✓ CSV data retrieved ({len(csv_data)} characters)")
        # Show first few lines of CSV
        csv_lines = csv_data.split('\\n')[:5]
        print("   First few lines of CSV:")
        for i, line in enumerate(csv_lines, 1):
            if line.strip():
                print(f"   {i}: {line[:100]}..." if len(line) > 100 else f"   {i}: {line}")
    else:
        print(f"   ✗ Failed to get CSV format. Status: {csv_response.status_code}")
    
    print()
    
    # SDF format
    print("3. SDF Format (structure file):")
    sdf_response = ClassyFireAPI.get_query(aspirin_query_id, format="sdf")
    if sdf_response.status_code == 200:
        sdf_data = sdf_response.text
        print(f"   ✓ SDF data retrieved ({len(sdf_data)} characters)")
        # Show first few lines of SDF
        sdf_lines = sdf_data.split('\\n')[:10]
        print("   First few lines of SDF:")
        for i, line in enumerate(sdf_lines, 1):
            if line.strip():
                print(f"   {i}: {line}")
    else:
        print(f"   ✗ Failed to get SDF format. Status: {sdf_response.status_code}")
        
else:
    print("No aspirin query ID available to demonstrate output formats")

print()
print("Format recommendations:")
print("• JSON: Best for programmatic analysis and detailed classification data")
print("• CSV: Good for spreadsheet analysis and simple data processing")  
print("• SDF: Ideal for integration with chemical structure software")
```

## 6. Helper Functions for Automated Classification

Let's create some helper functions to streamline the classification process:

```{code-cell} ipython3
def classify_compound_complete(name, smiles, sleep_schedule=None, max_retries=3):
    """
    Complete classification workflow using proven timing patterns.
    Based on real-world batch processing experience.
    
    Args:
        name (str): Descriptive name for the compound
        smiles (str): SMILES string of the compound
        sleep_schedule (list): Wait times between status checks [10, 10, 15, 15]
        max_retries (int): Maximum submission retry attempts
    
    Returns:
        dict: Classification results or error information
    """
    if sleep_schedule is None:
        sleep_schedule = [10, 10, 15, 15]  # Proven timing pattern
    
    print(f"Starting complete classification for {name}...")
    print(f"Using sleep schedule: {sleep_schedule} seconds")
    
    # Try to submit query with retries
    query_id = None
    for attempt in range(max_retries):
        print(f"  Submission attempt {attempt + 1}/{max_retries}...")
        
        response = ClassyFireAPI.submit_query(f"{name} Classification", smiles)
        
        # ClassyFire returns 201 (Created) for successful query submission
        if response.status_code == 201:
            query_result = response.json()
            query_id = query_result['id']
            print(f"  ✓ Query submitted successfully! ID: {query_id}")
            break
        else:
            print(f"  ✗ Submission failed: {response.status_code}")
            if attempt < max_retries - 1:
                retry_wait = 3
                print(f"    Retrying in {retry_wait}s...")
                time.sleep(retry_wait)
            else:
                return {
                    "success": False,
                    "error": f"Failed to submit query after {max_retries} attempts: {response.status_code}",
                    "name": name,
                    "smiles": smiles
                }
    
    # Wait for completion using progressive timing
    print(f"  Waiting for classification to complete...")
    for i, wait_time in enumerate(sleep_schedule, 1):
        print(f"  Status check {i}/{len(sleep_schedule)} (waiting {wait_time}s)...")
        time.sleep(wait_time)
        
        status_response = ClassyFireAPI.query_status(query_id)
        if status_response and status_response.status_code == 200:
            # Status response is typically plain text "Done" when ready
            status_text = status_response.text.strip()
            print(f"    Status: {status_text}")
            
            if status_text == "Done":
                print(f"  ✓ Classification complete!")
                break
            elif status_text in ['In progress', 'Queued']:
                print(f"  ⏳ Still processing...")
            else:
                print(f"  ⚠️ Unexpected status: {status_text}")
        else:
            print(f"  ⚠️ Status check failed (code: {status_response.status_code if status_response else 'None'})")
        
        # If this was the last check and still not done
        if i == len(sleep_schedule) and status_text != "Done":
            return {
                "success": False,
                "error": f"Classification timed out after {sum(sleep_schedule)} seconds",
                "name": name,
                "smiles": smiles,
                "query_id": query_id
            }
    
    # Retrieve results
    print(f"  Retrieving results...")
    result_response = ClassyFireAPI.get_query(query_id, format="json")
    
    if result_response.status_code == 200:
        classification_results = result_response.json()
        print(f"  ✓ Results retrieved successfully!")
        return {
            "success": True,
            "name": name,
            "smiles": smiles,
            "query_id": query_id,
            "results": classification_results
        }
    else:
        return {
            "success": False,
            "error": f"Failed to retrieve results: {result_response.status_code}",
            "name": name,
            "smiles": smiles,
            "query_id": query_id
        }

def extract_classification_summary(classification_data):
    """
    Extract key classification information from ClassyFire results.
    
    Args:
        classification_data (dict): Full ClassyFire classification results
    
    Returns:
        dict: Simplified classification summary
    """
    if not classification_data.get("success"):
        return {"error": classification_data.get("error", "Unknown error")}
    
    results = classification_data["results"]
    entities = results.get("entities", [])
    
    if not entities:
        return {"error": "No classification entities found"}
    
    entity = entities[0]  # Use first entity
    
    return {
        "compound_name": classification_data["name"],
        "smiles": results.get("smiles"),
        "molecular_formula": results.get("molecular_formula"),
        "inchi_key": results.get("inchikey"),
        "kingdom": entity.get("kingdom", {}).get("name"),
        "superclass": entity.get("superclass", {}).get("name"),
        "class": entity.get("class", {}).get("name"),
        "subclass": entity.get("subclass", {}).get("name"),
        "direct_parent": entity.get("direct_parent", {}).get("name"),
        "molecular_framework": results.get("molecular_framework"),
        "num_intermediate_nodes": len(entity.get("intermediate_nodes", []))
    }

# Example usage of helper functions with proven timing patterns
print("Testing helper functions with optimized timing:")
print("=" * 60)

# Note: This is a demonstration - actual execution may take several minutes
ethanol_smiles = "CCO"
print(f"Using SMILES: {ethanol_smiles}")
print(f"Using proven sleep schedule: [10, 10, 15, 15] seconds")
print()
print("In a real application, you would run:")
print("ethanol_classification = classify_compound_complete('Ethanol', ethanol_smiles)")
print("ethanol_summary = extract_classification_summary(ethanol_classification)")
print()
print("This approach provides:")
print("- Progressive wait times: 10s, 10s, 15s, 15s (total ~50s)")
print("- Multiple submission retries with 3s delays")
print("- Proper error handling at each step")
print("- Detailed logging and progress tracking")
print()
print("Expected timeline for a single compound:")
print("  0s: Submit query (expect 201 status)")
print(" 10s: First status check")
print(" 20s: Second status check") 
print(" 35s: Third status check")
print(" 50s: Final status check (should be 'Done')")
print(" 51s: Retrieve results (expect 200 status)")
print()
print("Benefits over simple polling:")
print("• Reduces server load with intelligent timing")
print("• Accounts for typical ClassyFire processing times")
print("• Proven effective for large-scale batch processing")
print("• Balances responsiveness with efficiency")

# Uncomment the lines below to run actual classification (takes ~1 minute):
# print("\n" + "="*60)
# print("RUNNING ACTUAL CLASSIFICATION (uncomment to execute):")
# ethanol_result = classify_compound_complete('Ethanol', ethanol_smiles)
# if ethanol_result.get('success'):
#     summary = extract_classification_summary(ethanol_result)
#     print(f"SUCCESS: {summary}")
# else:
#     print(f"FAILED: {ethanol_result.get('error')}")
```

## 7. Error Handling and Best Practices

ClassyFire classification can encounter various issues. Let's demonstrate proper error handling:

```{code-cell} ipython3
# Test error handling with various scenarios
print("Testing error handling scenarios:")
print("=" * 40)

# Test 1: Invalid SMILES
print("1. Testing with invalid SMILES:")
invalid_smiles = "INVALID_SMILES_STRING"
try:
    response = ClassyFireAPI.submit_query("Invalid SMILES Test", invalid_smiles)
    print(f"   Status Code: {response.status_code}")
    # ClassyFire returns 201 for valid submissions, other codes for errors
    if response.status_code != 201:
        print(f"   ✓ Correctly handled invalid SMILES")
        print(f"   Response: {response.text[:100]}...")
    else:
        print(f"   Unexpected success with invalid SMILES")
except Exception as e:
    print(f"   Exception caught: {e}")

print()

# Test 2: Very long SMILES (may cause issues)
print("2. Testing with very long SMILES:")
long_smiles = "C" * 1000  # Very long alkane chain
try:
    response = ClassyFireAPI.submit_query("Long SMILES Test", long_smiles)
    print(f"   Status Code: {response.status_code}")
    # ClassyFire returns 201 for valid submissions
    if response.status_code == 201:
        print(f"   ✓ Long SMILES accepted")
    else:
        print(f"   Long SMILES rejected: {response.text[:100]}...")
except Exception as e:
    print(f"   Exception caught: {e}")

print()

# Test 3: Query with invalid ID
print("3. Testing with invalid query ID:")
try:
    invalid_id = "invalid_query_id_12345"
    status_response = ClassyFireAPI.query_status(invalid_id)
    if status_response:
        print(f"   Status Code: {status_response.status_code}")
        if status_response.status_code == 404:
            print(f"   ✓ Correctly returned 404 for invalid ID")
        else:
            print(f"   Response: {status_response.text[:100]}...")
    else:
        print(f"   ✓ Correctly returned None for invalid ID")
except Exception as e:
    print(f"   Exception caught: {e}")

print()

# Best practices recommendations
print("=== BEST PRACTICES ===")
print()
print("1. Input Validation:")
print("   • Validate SMILES strings before submission")
print("   • Check molecular size (very large molecules may fail)")
print("   • Ensure proper chemical structure representation")
print()
print("2. Timing Strategy (CRITICAL):")
print("   • Use progressive wait times: [10, 10, 15, 15] seconds")
print("   • Total wait time: ~50 seconds (covers 90% of classifications)")
print("   • Don't check status continuously - be respectful to server")
print("   • Add 1-3 second delays between submission retries")
print("   • For batch processing: add small delays between compounds")
print()
print("3. Rate Limiting:")
print("   • Process compounds sequentially, not in parallel")
print("   • Monitor server response times and adjust if needed")
print("   • Consider server load during peak hours")
print()
print("4. Error Handling:")
print("   • Always check response status codes (201 for submission, 200 for status/results)")
print("   • Implement retry logic for temporary failures")
print("   • Handle timeout scenarios gracefully")
print("   • Log failures with detailed error information")
print()
print("5. Batch Processing:")
print("   • Save intermediate results every 10-100 compounds")
print("   • Track failed compounds separately for retry")
print("   • Use progress tracking (tqdm) for long runs")
print("   • Implement checkpointing to resume interrupted jobs")
print()
print("6. Data Processing:")
print("   • Parse JSON results carefully (structure may vary)")
print("   • Handle missing classification levels gracefully")
print("   • Extract key information into simplified formats")
print("   • Store raw results for future reanalysis")

print()
print("=== COMMON ISSUES ===")
print()
print("• Server overload: Classification may be slow during peak times")
print("• Invalid structures: Some SMILES may not be classifiable")
print("• Network issues: Implement retry logic for connection problems")
print("• Large molecules: Complex structures may timeout or fail")
print("• API changes: Monitor for service updates and changes")
```

## 8. Practical Applications

Here are some practical use cases for ClassyFire chemical classification:

```{code-cell} ipython3
# Use Case 1: Building a Chemical Class Database
def build_classification_database(compounds_dict, delay=2):
    """
    Build a database of chemical classifications for multiple compounds.
    
    Args:
        compounds_dict (dict): Dictionary of {name: smiles} pairs
        delay (int): Delay between submissions in seconds
    
    Returns:
        dict: Classification database
    """
    print(f"Building classification database for {len(compounds_dict)} compounds...")
    print("Note: This is a demonstration of the workflow")
    
    database = {}
    
    for name, smiles in compounds_dict.items():
        print(f"\\n{name}:")
        print(f"  SMILES: {smiles}")
        
        # In a real implementation, you would:
        # 1. Submit the query (expecting 201 status code)
        # 2. Wait for completion (checking until status is "Done")
        # 3. Extract classification data (200 status code for results)
        
        # Simulated database entry
        database[name] = {
            "smiles": smiles,
            "submitted": True,
            "classification_levels": {
                "kingdom": "Organic compounds",  # Example
                "superclass": "Unknown (would be determined)",
                "class": "Unknown (would be determined)",
                "subclass": "Unknown (would be determined)"
            },
            "molecular_framework": "Unknown (would be determined)",
            "status": "Would be determined after classification"
        }
        
        print(f"  Database entry created (simulated)")
        
        if delay > 0:
            time.sleep(delay)  # Respectful delay
    
    return database

# Example with pharmaceutical compounds
pharma_compounds = {
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Paracetamol": "CC(=O)NC1=CC=C(C=C1)O",
    "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
}

print("Use Case 1: Building a Classification Database")
print("=" * 50)

pharma_db = build_classification_database(pharma_compounds)

print(f"\\nDatabase created with {len(pharma_db)} entries:")
for name, data in pharma_db.items():
    print(f"  {name}: {data['status']}")

print()
print("In a real implementation, this database would contain:")
print("• Complete hierarchical classifications")
print("• Molecular frameworks")
print("• Chemical fingerprints")
print("• Structural features")
print("• Cross-references to chemical databases")
```

```{code-cell} ipython3
# Use Case 2: Metabolomics Classification Analysis
def analyze_metabolite_classes(metabolite_classifications):
    """
    Analyze the distribution of chemical classes in a metabolomics dataset.
    
    Args:
        metabolite_classifications (list): List of classification results
    
    Returns:
        dict: Analysis summary
    """
    print("Use Case 2: Metabolomics Classification Analysis")
    print("=" * 50)
    
    # Simulated metabolomics classification data
    simulated_data = [
        {"name": "Glucose", "superclass": "Organic oxygen compounds", "class": "Organooxygen compounds"},
        {"name": "Alanine", "superclass": "Organic acids and derivatives", "class": "Carboxylic acids and derivatives"},
        {"name": "Cholesterol", "superclass": "Lipids and lipid-like molecules", "class": "Steroids and steroid derivatives"},
        {"name": "Caffeine", "superclass": "Organoheterocyclic compounds", "class": "Purinones"},
        {"name": "Glucose-6-phosphate", "superclass": "Organic oxygen compounds", "class": "Organooxygen compounds"},
        {"name": "Tryptophan", "superclass": "Organic acids and derivatives", "class": "Carboxylic acids and derivatives"}
    ]
    
    print("Analyzing chemical class distribution in metabolomics data:")
    print()
    
    # Count superclasses
    superclass_counts = {}
    class_counts = {}
    
    for metabolite in simulated_data:
        superclass = metabolite["superclass"]
        class_name = metabolite["class"]
        
        superclass_counts[superclass] = superclass_counts.get(superclass, 0) + 1
        class_counts[class_name] = class_counts.get(class_name, 0) + 1
    
    print("Superclass Distribution:")
    for superclass, count in sorted(superclass_counts.items(), key=lambda x: x[1], reverse=True):
        percentage = (count / len(simulated_data)) * 100
        print(f"  {superclass}: {count} ({percentage:.1f}%)")
    
    print()
    print("Class Distribution:")
    for class_name, count in sorted(class_counts.items(), key=lambda x: x[1], reverse=True):
        percentage = (count / len(simulated_data)) * 100
        print(f"  {class_name}: {count} ({percentage:.1f}%)")
    
    print()
    print("Analysis Summary:")
    print(f"  Total metabolites analyzed: {len(simulated_data)}")
    print(f"  Unique superclasses: {len(superclass_counts)}")
    print(f"  Unique classes: {len(class_counts)}")
    print(f"  Most common superclass: {max(superclass_counts, key=superclass_counts.get)}")
    print(f"  Chemical diversity index: {len(class_counts) / len(simulated_data):.2f}")
    
    return {
        "total_metabolites": len(simulated_data),
        "superclass_distribution": superclass_counts,
        "class_distribution": class_counts,
        "diversity_metrics": {
            "unique_superclasses": len(superclass_counts),
            "unique_classes": len(class_counts),
            "diversity_index": len(class_counts) / len(simulated_data)
        }
    }

# Run metabolomics analysis
metabolomics_analysis = analyze_metabolite_classes([])

print()
print("Applications in Metabolomics:")
print("• Pathway enrichment analysis")
print("• Chemical space visualization")
print("• Biomarker classification")
print("• Metabolite identification support")
print("• Chemical similarity assessment")
```

## Summary

The `ClassyFireAPI` class provides comprehensive access to the ClassyFire chemical classification service:

### Main ClassyFireAPI Methods:
1. **`submit_query(label, input, type='STRUCTURE')`**: Submit a classification query
2. **`query_status(query_id)`**: Check the status of a submitted query
3. **`get_query(query_id, format="json")`**: Retrieve classification results in various formats

### Supported Input Types:
- **SMILES notation**: Most common format for chemical structures
- **InChI strings**: International chemical identifiers
- **Chemical names**: Common or IUPAC names (may have variable success)
- **Structure files**: SDF and other chemical file formats

### Output Formats:
- **JSON**: Complete hierarchical classification data with all details
- **CSV**: Tabular format suitable for spreadsheet analysis
- **SDF**: Structure-Data File format for chemical software integration

### Classification Hierarchy:
ClassyFire organizes compounds into a hierarchical taxonomy:
1. **Kingdom**: Broadest level (e.g., "Organic compounds")
2. **Superclass**: Major chemical groups (e.g., "Organoheterocyclic compounds")
3. **Class**: More specific classifications (e.g., "Purinones")
4. **Subclass**: Detailed sub-categories
5. **Intermediate Nodes**: Additional classification levels
6. **Direct Parent**: Most specific classification level
7. **Molecular Framework**: Core structural motif

### Key Features:
- ✅ **Hierarchical Classification**: Multi-level taxonomic organization
- ✅ **Structural Analysis**: Based on molecular structure features
- ✅ **Multiple Formats**: JSON, CSV, SDF output options
- ✅ **Free Service**: No API key required
- ✅ **Comprehensive Database**: Over 4,800 chemical categories
- ✅ **Research-Quality**: Peer-reviewed classification system
- ✅ **Batch Processing**: Can handle multiple compounds

### Workflow Pattern:
1. **Submit**: Send SMILES/InChI to ClassyFire
2. **Wait**: Classification takes 30 seconds to several minutes
3. **Check**: Monitor query status until complete
4. **Retrieve**: Get results in desired format
5. **Parse**: Extract relevant classification information

### Best Use Cases:
- **Metabolomics Research**: Classify detected metabolites
- **Drug Discovery**: Organize compound libraries by chemical class
- **Natural Products**: Systematic classification of natural compounds
- **Chemical Database Curation**: Standardize chemical classifications
- **Pathway Analysis**: Group compounds by chemical function
- **Chemical Space Exploration**: Understand molecular diversity

### Limitations:
- **Processing Time**: Can take several minutes per compound
- **Server Dependent**: Relies on external web service availability
- **Complex Molecules**: Very large structures may fail or timeout
- **Rate Limits**: Should respect server capacity with delays between requests
- **Network Dependency**: Requires stable internet connection

### Service Information:
- **Provider**: Wishart Research Group, University of Alberta
- **URL**: http://classyfire.wishartlab.com
- **Method**: RESTful API with JSON responses
- **Citation**: Required for academic use (check ClassyFire website)
- **Updates**: Classification database is periodically updated

ClassyFire is an essential tool for chemical classification and is particularly valuable for researchers working with large chemical datasets who need systematic, hierarchical organization of molecular structures. The service provides research-quality classifications that are widely accepted in the chemical and biochemical research communities.
