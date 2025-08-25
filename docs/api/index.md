# API Reference

This section provides comprehensive documentation for all PROVESID API modules.

## Overview

PROVESID provides Python interfaces to several major chemical databases and web services:

- **[PubChem API](pubchem.md)** - Standard PubChem REST API with enhanced data access ✨ *Enhanced*
- **[PubChem PUG View](pubchemview.md)** - Advanced property extraction from PubChem ✨ *Enhanced*
- **[NCI Resolver](nci_resolver.md)** - Chemical identifier resolution
- **[CAS Common Chemistry](cascommonchem.md)** - CAS Registry data access
- **[OPSIN](opsin.md)** - IUPAC name to structure conversion
- **[ClassyFire](classyfire.md)** - Chemical structure classification

## Quick Reference

### Basic Usage Pattern

All API classes follow a similar pattern:

```python
# Import the desired API
from provesid import PubChemAPI  # Enhanced standard API
from provesid.pubchemview import PubChemPUGViewAPI  # Advanced properties
from provesid.cascommonchem import CASCommonChem
from provesid.opsin import OPSIN
from provesid.classyfire import ClassyFireAPI

# Initialize APIs
pc = PubChemAPI()  # Standard PubChem API
pug_view = PubChemPUGViewAPI()  # Advanced PubChem properties
cas_api = CASCommonChem()
opsin = OPSIN()

# ClassyFireAPI uses static methods
result = ClassyFireAPI.submit_query("label", "CCO")
```

### Common Data Flow

```python
# Start with a compound name
compound_name = "caffeine"

# 1. Convert name to structure (if needed)
opsin = OPSIN()
structure = opsin.get_id(compound_name)

if structure['status'] == 'SUCCESS':
    smiles = structure['smiles']
    
    # 2. Get basic compound info (Enhanced PubChem API)
    pc = PubChemAPI()
    cids = pc.get_cids_by_name(compound_name)
    if cids:
        basic_info = pc.get_basic_compound_info(cids[0])
        # Direct access - no wrapper needed!
        compound = pc.get_compound_by_cid(cids[0])
    
    # 3. Get detailed experimental properties
    pug_view = PubChemPUGViewAPI() 
    properties = pug_view.get_compound_properties_by_smiles(smiles)
```

## Recent API Enhancements ✨

### PubChem API Improvements

**Enhanced Data Access**: Removed redundant wrapper structures for cleaner code:
```python
# Before: compound["PC_Compounds"][0]
# Now: compound (direct access!)
compound = pc.get_compound_by_cid(cid)
```

**New Search Methods**:
- `get_cids_by_inchikey()` - Search by InChI Key
- Improved `get_cids_by_smiles()` - Returns clean lists
- Multi-domain search capabilities

**Enhanced Pattern Recognition in PubChem View**:
- Support for "log Kow = 1.19" format
- Comprehensive solubility pattern matching
- Improved vapor pressure extraction
    smiles = structure['smiles']
    
    # 2. Get detailed properties
    pubchem = PubChemPUGViewAPI()
    properties = pubchem.get_compound_properties_by_smiles(smiles)
    
    # 3. Get CAS information
    cas_api = CASCommonChem()
    cas_data = cas_api.smiles_to_detail(smiles)
    
    # 4. Get classification
    response = ClassyFireAPI.submit_query("caffeine_classification", smiles)
    
    # 5. Get ChEBI information
    chebi = ChEBI()
    chebi_results = chebi.search_by_name(compound_name)
    if chebi_results:
        chebi_id = chebi_results[0]['chebiId']
        chebi_entity = chebi.get_complete_entity(chebi_id)
```

## Module Comparison

| Feature | PubChem API | PubChem View | CAS Common | OPSIN | ClassyFire | NCI Resolver |
|---------|-------------|--------------|------------|-------|------------|--------------|
| **Primary Use** | Standard data | Properties | Registry data | Name→Structure | Classification | ID conversion |
| **Input Types** | CID, Name, SMILES, InChI Key | CID, Name, SMILES | CAS, Name, SMILES | IUPAC names | SMILES, InChI | Various IDs |
| **Output Format** | JSON ✨ *Clean* | JSON, DataFrame | JSON | JSON | JSON, SDF, CSV | JSON |
| **Rate Limits** | Yes | Yes | Unofficial | Unofficial | Unofficial | Yes |
| **Batch Support** | Yes ✨ *Enhanced* | Yes | Manual | Yes | Manual | Yes |
| **Recent Updates** | ✨ *Major* | ✨ *Enhanced* | - | - | - | - |

## Authentication Requirements

| Service | Authentication | Notes |
|---------|---------------|-------|
| PubChem API | None | Rate limits apply, enhanced data access ✨ |
| PubChem PUG View | None | Rate limits apply, advanced properties |
| CAS Common Chemistry | None | Free tier available |
| OPSIN | None | Cambridge University service |
| ClassyFire | None | Long processing times |
| NCI Resolver | None | Rate limits apply |

## Error Handling Best Practices

All modules implement consistent error handling:

```python
try:
    result = api.method(parameter)
    
    if result:  # Check if result exists
        # Process result
        pass
    else:
        print("No data found")
        
except requests.exceptions.Timeout:
    print("Request timed out")
except requests.exceptions.ConnectionError:
    print("Connection error")
except Exception as e:
    print(f"Unexpected error: {e}")
```

## Performance Considerations

### Rate Limiting
- Implement delays between requests
- Use caching for repeated queries
- Batch process when supported

### Memory Management
- Process large datasets in chunks
- Clear unnecessary variables
- Use generators for large result sets

### Network Optimization
- Set appropriate timeouts
- Implement retry logic
- Handle network failures gracefully

## Integration Examples

### Cross-Platform Validation

```python
def validate_compound_across_platforms(identifier, id_type="name"):
    """Validate compound information across multiple platforms"""
    
    results = {}
    
    # Get structure
    if id_type == "name":
        opsin = OPSIN()
        structure = opsin.get_id(identifier)
        if structure['status'] == 'SUCCESS':
            smiles = structure['smiles']
            results['opsin'] = structure
        else:
            return None
    else:
        smiles = identifier
    
    # Get data from each platform
    platforms = {
        'pubchem': lambda: PubChemPUGViewAPI().get_compound_properties_by_smiles(smiles),
        'cas': lambda: CASCommonChem().smiles_to_detail(smiles),
        'nci': lambda: NCIResolverAPI().smiles_to_names(smiles),
        'chebi': lambda: ChEBI().search_by_name(identifier) if id_type == "name" else None
    }
    
    for platform, method in platforms.items():
        try:
            data = method()
            if data:
                results[platform] = data
        except Exception as e:
            print(f"Error with {platform}: {e}")
            results[platform] = None
    
    return results
```

### Comprehensive Compound Report

```python
def generate_compound_report(compound_name):
    """Generate comprehensive compound report"""
    
    report = {
        'compound_name': compound_name,
        'timestamp': datetime.now().isoformat(),
        'structure': None,
        'properties': None,
        'classification': None,
        'registry_data': None,
        'chebi_data': None
    }
    
    # Structure information
    opsin = OPSIN()
    structure = opsin.get_id(compound_name)
    if structure['status'] == 'SUCCESS':
        report['structure'] = structure
        smiles = structure['smiles']
        
        # Properties
        pubchem = PubChemPUGViewAPI()
        properties = pubchem.get_compound_properties_by_smiles(smiles)
        report['properties'] = properties
        
        # Registry data
        cas_api = CASCommonChem()
        cas_data = cas_api.smiles_to_detail(smiles)
        report['registry_data'] = cas_data
        
        # Classification (submit query - results may take time)
        response = ClassyFireAPI.submit_query(f"{compound_name}_report", smiles)
        if response and response.status_code == 200:
            query_id = response.json()['id']
            report['classification_query_id'] = query_id
        
        # ChEBI data
        chebi = ChEBI()
        chebi_results = chebi.search_by_name(compound_name)
        if chebi_results:
            report['chebi_data'] = chebi_results[0]
            # Get detailed entity information
            chebi_id = chebi_results[0]['chebiId']
            entity = chebi.get_complete_entity(chebi_id)
            if entity:
                report['chebi_entity'] = entity
    
    return report
```

## Troubleshooting

### Common Issues

1. **Timeout Errors**
   - Increase timeout parameters
   - Check network connectivity
   - Verify service availability

2. **Rate Limiting**
   - Implement request delays
   - Use batch operations where available
   - Monitor response headers

3. **Data Not Found**
   - Verify input format
   - Try alternative identifiers
   - Check service coverage

4. **Invalid Responses**
   - Validate input data
   - Check API documentation for changes
   - Implement response validation

### Debug Mode

Enable debug logging for detailed troubleshooting:

```python
import logging

# Enable debug logging
logging.basicConfig(level=logging.DEBUG)

# Use APIs normally - detailed logs will be shown
api = PubChemPUGViewAPI()
result = api.get_compound_properties_by_name("aspirin")
```

## Support and Resources

- **GitHub Issues**: Report bugs and request features
- **Documentation**: Comprehensive guides and examples
- **Community**: Discussions and community support
- **API Updates**: Monitor upstream API changes

For specific API documentation, see the individual module pages linked above.
