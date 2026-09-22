# Testing Guide for PROVESID

Tests are run locally. There is no CI test workflow — the only GitHub Actions
workflow in this repository deploys the documentation.

## Running the Tests

```bash
python -m pytest tests/ -v
```

### Including the CAS Common Chemistry Tests

The CAS tests need an API key and are skipped without one. Provide it in any of
these ways:

```bash
# Environment variable
export CAS_API_KEY="your-cas-api-key-here"

# Or a file holding the key
export CAS_API_KEY_FILE="/path/to/your/api-key-file.txt"
```

PowerShell (Windows):

```powershell
$env:CAS_API_KEY="your-cas-api-key-here"
```

`CCC_API_KEY` is accepted as an alternative name for the same key.

### Skipping the CAS Tests Explicitly

```bash
python -m pytest tests/ -v --ignore=tests/test_cascommonchem.py
```

## Running the Docstring Examples

Every public object's docstring carries an example, and they can be run as
doctests:

```bash
python -m pytest --doctest-modules src/provesid
```

The ordinary run above never collects them. Examples that need the network
are marked `# doctest: +SKIP`, and the run enforces it: `src/conftest.py`
refuses every outbound connection, so an example that forgets the marker
fails instead of quietly asking a live service. To re-record an online
example's output, lift the block for that run:

```bash
PROVESID_DOCTEST_ALLOW_NETWORK=1 python -m pytest --doctest-modules src/provesid
```

Examples that read an offline database run against the installed copy; a
module whose database is not installed is skipped, never downloaded
(`src/conftest.py` again). The outputs shown were taken from the databases as
installed on the day they were written, so a later release can change a count
or an ordering without anything being wrong.

## Building the Documentation

```bash
mkdocs build --strict
```

The API reference is rendered from the same docstrings, so a malformed one
fails the build, and so does a link to a page or heading that does not exist.
The tutorials are read from `examples/` (`scripts/mkdocs_hooks.py`); there is
no copy under `docs/`. `scripts/validate_docs_local.sh --execute` also runs
the MyST tutorials, which needs the network.

## Test Categories

### Core API Tests
- ✅ **PubChem API** - No authentication required
- ✅ **NCI Resolver** - No authentication required  
- ✅ **OPSIN** - No authentication required
- ✅ **ClassyFire** - No authentication required

### Authentication-Required Tests
- 🔐 **CAS Common Chemistry** - Requires API key
  - Skipped automatically when no key is available

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