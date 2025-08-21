# Installation

## Requirements

PROVESID requires Python 3.8 or higher and the following dependencies:

- `requests` >= 2.25.0
- `pandas` >= 1.3.0
- `numpy` >= 1.21.0

## Install from PyPI

The easiest way to install PROVESID is using pip:

```bash
pip install provesid
```

## Install from Source

To install the latest development version from GitHub:

```bash
git clone https://github.com/provesid/provesid.git
cd provesid
pip install -e .
```

## Development Installation

For development, install with additional dependencies:

```bash
git clone https://github.com/provesid/provesid.git
cd provesid
pip install -e ".[dev]"
```

This includes:
- `pytest` for testing
- `black` for code formatting
- `flake8` for linting
- `mkdocs` for documentation

## Verify Installation

Test your installation by importing the main modules:

```python
import provesid
from provesid import PubChemAPI, NCIChemicalIdentifierResolver, PubChemView

# Check version
print(provesid.__version__)

# Quick test
api = PubChemAPI()
result = api.get_compound_by_cid(2244)  # Aspirin
print(f"Successfully retrieved compound with CID 2244")
```

## Optional Dependencies

Some features require additional packages:

### For Jupyter Notebook Support
```bash
pip install jupyter ipywidgets
```

### For Advanced Data Analysis
```bash
pip install scipy scikit-learn matplotlib seaborn
```

### For Chemical Structure Visualization
```bash
pip install rdkit-pypi
```

## Troubleshooting

### Common Issues

**Import Error**: If you get import errors, make sure you have the correct Python version and all dependencies installed:

```bash
python --version  # Should be 3.8+
pip list | grep provesid
```

**Network Issues**: PROVESID makes API calls to external services. If you're behind a firewall or proxy, you may need to configure your network settings:

```python
import requests
# Set proxy if needed
proxies = {
    'http': 'http://proxy.company.com:8080',
    'https': 'https://proxy.company.com:8080'
}
# Pass proxies to API calls as needed
```

**Rate Limiting**: If you encounter rate limiting errors, the package includes automatic retry mechanisms. You can adjust the delay between requests:

```python
from provesid import PubChemAPI

# Increase delay between requests
api = PubChemAPI(pause_time=1.0)  # 1 second between requests
```

### Getting Help

If you encounter issues:

1. Check the [GitHub Issues](https://github.com/provesid/provesid/issues) for known problems
2. Review the [API documentation](../api/pubchem.md) for usage examples
3. Create a new issue with detailed error information

## Next Steps

Once installed, head to the [Quick Start Guide](quickstart.md) to begin using PROVESID!
