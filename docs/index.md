# PROVESID

A comprehensive Python package for chemical identifier resolution and experimental property extraction from multiple chemical databases and APIs.

## Overview

PROVESID (PROVenance and Experimental Structuring of Identifier Data) provides unified interfaces to access chemical information from various sources including:

- **PubChem**: Comprehensive chemical database with compound, substance, and bioassay information
- **PubChem PUG View**: Experimental property data extraction with full reference information
- **NCI Chemical Identifier Resolver**: Chemical structure identifier conversion service
- **CAS Common Chemistry**: Open chemistry database from the Chemical Abstracts Service
- **ClassyFire**: Chemical taxonomy and classification
- **OPSIN**: Chemical name to structure conversion

## Key Features

### 🔍 **Multi-Source Chemical Data Access**
- Unified interfaces to major chemical databases
- Consistent error handling and response formats
- Automatic rate limiting and retry mechanisms

### 🧪 **Experimental Property Extraction**
- Extract experimental properties from PubChem PUG View
- Parse values, units, and full reference information
- Generate structured DataFrames for analysis

### 🔄 **Chemical Identifier Conversion**
- Convert between SMILES, InChI, CAS numbers, and names
- Resolve chemical identifiers across different formats
- Validate and standardize chemical structures

### 📊 **Data Processing & Analysis**
- Batch processing capabilities for large datasets
- DataFrame output for easy integration with pandas
- Comprehensive error handling and logging

## Quick Start

### Installation

```bash
pip install provesid
```

### Basic Usage

```python
from provesid import PubChemAPI, NCIChemicalIdentifierResolver, PubChemView

# Get compound information from PubChem
api = PubChemAPI()
compound = api.get_compound_by_cid(2244)  # Aspirin
properties = api.get_compound_properties([2244], ['MolecularWeight', 'MolecularFormula'])

# Convert chemical identifiers
resolver = NCIChemicalIdentifierResolver()
smiles = resolver.resolve('aspirin', 'smiles')
inchi = resolver.resolve('CCO', 'stdinchi')  # Ethanol SMILES to InChI

# Extract experimental properties
view = PubChemView()
melting_points = view.get_experimental_properties(2244, 'Melting Point')
df = view.experimental_properties_to_dataframe(2244, 'Melting Point')
```

### Property Extraction Example

```python
from provesid.pubchemview import get_experimental_properties_table

# Get a comprehensive table of experimental properties
table = get_experimental_properties_table(2244, 'Boiling Point')
print(table)
#   CID StringWithMarkup          Value Unit                    Reference
# 0  2244     139 °C at 760 mmHg  139    °C    J. Chem. Eng. Data 1996, 41, 1190-1193
# 1  2244     140 °C              140    °C    Lange's Handbook of Chemistry, 1985
```

## API Documentation

Comprehensive API documentation is available for all modules:

- [PubChem API](api/pubchem.md) - Access to PubChem compound, substance, and bioassay data
- [PubChem View](api/pubchemview.md) - Experimental property extraction and reference parsing
- [NCI Resolver](api/nci_resolver.md) - Chemical identifier conversion and validation
- [Common Chemistry](api/cascommonchem.md) - CAS Common Chemistry database access
- [ClassyFire](api/classyfire.md) - Chemical classification and taxonomy
- [OPSIN](api/opsin.md) - Chemical name to structure conversion

## Examples

Explore comprehensive tutorials for each API:

- [PubChem Tutorial](examples/pubchem/pubchem_tutorial.ipynb) - Complete PubChem API guide with enhanced features
- [CAS Common Chemistry](examples/CCC/CAS_Common_Chemistry_tutorial.ipynb) - Working with CAS Registry data
- [ChEBI Tutorial](examples/ChEBI/ChEBI_tutorial.ipynb) - Chemical Entities of Biological Interest database
- [ClassyFire Tutorial](examples/ClassyFire/classyfire_tutorial.ipynb) - Chemical structure classification
- [OPSIN Tutorial](examples/OPSIN/opsin_tutorial.ipynb) - IUPAC name to structure conversion
- [Chemical ID Resolver](examples/resolver/chem_id_resolver_tutorial.ipynb) - NCI chemical identifier resolution

## Development

PROVESID is actively developed and welcomes contributions:

- [Contributing Guidelines](development/contributing.md)
- [Testing Documentation](development/testing.md)
- [API Design Principles](development/api_design.md)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use PROVESID in your research, please cite:

```bibtex
@software{provesid,
  title={PROVESID: A Python Package for Chemical Identifier Resolution and Property Extraction},
  author={PROVESID Team},
  year={2024},
  url={https://github.com/provesid/provesid}
}
```
