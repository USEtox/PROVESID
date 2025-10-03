# Testing Guide for PROVESID

This guide explains how to run tests locally and in CI/CD environments.

## Local Testing

### Running All Tests (Including CAS API)

To run all tests including CAS Common Chemistry API tests, you need to provide a CAS API key:

#### Method 1: Environment Variable (GitHub Actions compatible)
```bash
export CCC_API_KEY="your-cas-api-key-here"
python -m pytest tests/ -v
```

#### Method 1b: Environment Variable (Local testing)
```bash
export CAS_API_KEY="your-cas-api-key-here"
python -m pytest tests/ -v
```

#### Method 2: API Key File
```bash
export CAS_API_KEY_FILE="/path/to/your/api-key-file.txt"
python -m pytest tests/ -v
```

#### Method 3: PowerShell (Windows)
```powershell
# GitHub Actions compatible
$env:CCC_API_KEY="your-cas-api-key-here"
python -m pytest tests/ -v

# Or local testing
$env:CAS_API_KEY="your-cas-api-key-here"
python -m pytest tests/ -v
```

### Running Tests Without CAS API

To skip CAS API tests (no API key required):

```bash
python -m pytest tests/ -v --ignore=tests/test_cascommonchem.py
```

## Continuous Integration (CI)

### GitHub Actions Behavior

1. **Regular CI Tests** (`test.yml`):
   - Now includes CAS API tests with repository secret
   - Run all tests (PubChem, NCI, OPSIN, CAS, etc.)
   - Uses `CCC_API_KEY` repository secret

2. **Full Tests with API Keys** (`test-with-api-keys.yml`):
   - Manual workflow dispatch only
   - Can include CAS tests if API key is configured
   - Requires repository secret `CAS_API_KEY`

### Setting Up API Key for GitHub Actions (Optional)

If you want to enable CAS tests in GitHub Actions:

1. Go to your repository settings
2. Navigate to "Secrets and variables" → "Actions"
3. Add a new repository secret:
   - Name: `CCC_API_KEY`
   - Value: Your CAS API key
4. Both workflows can now run CAS tests automatically

## Test Categories

### Core API Tests
- ✅ **PubChem API** - No authentication required
- ✅ **NCI Resolver** - No authentication required  
- ✅ **OPSIN** - No authentication required
- ✅ **ClassyFire** - No authentication required

### Authentication-Required Tests
- 🔐 **CAS Common Chemistry** - Requires API key
  - Skipped in CI if no key available
  - Can run locally with proper API key setup

## Cache Tests

All APIs include cache functionality tests:
- Service-specific caching
- Cache clearing
- Cache information retrieval
- `use_cache` property functionality

## Coverage Reports

Generate coverage reports locally:

```bash
python -m pytest tests/ --cov=src/provesid --cov-report=html --cov-report=xml
```

The HTML report will be available in `htmlcov/index.html`.

## Troubleshooting

### CAS API Tests Failing
- Verify your API key is valid
- Check that the environment variable is set correctly
- Ensure your API key has not expired

### Import Errors
- Make sure you've installed the package: `pip install -e .`
- Check that all dependencies are installed: `pip install -e ".[test]"`

### Network Timeouts
- Some tests may fail due to network issues
- Re-run the tests if you suspect network problems
- Tests include timeout handling and retry logic