# API Key Configuration Guide

PROVESID supports persistent API key storage, making it easy to configure your API keys once and use them automatically in future sessions.

## Quick Start

### Setting Your CAS API Key (One-time setup)

```python
from provesid import set_cas_api_key

# Set your API key once
set_cas_api_key("your-cas-api-key-here")
```

### Using CAS API After Configuration

```python
from provesid import CASCommonChem

# No need to provide API key - it's automatically loaded!
cas_api = CASCommonChem()

# Use normally
result = cas_api.cas_to_detail("64-17-5")  # Ethanol
print(result['name'])  # "Ethanol"
```

## Configuration Management

### View Current Configuration

```python
from provesid import show_config

show_config()
```

Output:
```
PROVESID Configuration:
  Config directory: C:\Users\username\AppData\Roaming\PROVESID
  Config file: C:\Users\username\AppData\Roaming\PROVESID\config.json
  Config exists: True
  Configured services: cas
```

### Remove Stored API Key

```python
from provesid import remove_cas_api_key

remove_cas_api_key()
```

### Check if API Key is Configured

```python
from provesid import get_cas_api_key

api_key = get_cas_api_key()
if api_key:
    print("API key is configured")
else:
    print("No API key configured")
```

## API Key Priority Order

When initializing `CASCommonChem()`, the system checks for API keys in this order:

1. **Direct parameter**: `CASCommonChem(api_key="your-key")`
2. **File parameter**: `CASCommonChem(api_key_file="path/to/key.txt")`
3. **Persistent config**: Set via `set_cas_api_key("your-key")`
4. **Environment variables**: `CCC_API_KEY` or `CAS_API_KEY`

## Configuration Storage

### Windows
- Config directory: `%APPDATA%\PROVESID\`
- Config file: `%APPDATA%\PROVESID\config.json`

### Linux/macOS
- Config directory: `~/.config/provesid/`
- Config file: `~/.config/provesid/config.json`

## Security Considerations

- The API key is stored in plain text in the configuration file
- The configuration file is stored in your user directory (not system-wide)
- File permissions are set to be readable only by your user account
- For enhanced security, consider using environment variables in production

## Migration from Manual API Key Management

If you were previously providing the API key manually:

**Before:**
```python
cas_api = CASCommonChem(api_key="your-key")
```

**After (one-time setup):**
```python
# Run once to configure
from provesid import set_cas_api_key
set_cas_api_key("your-key")

# Then use anywhere without specifying the key
from provesid import CASCommonChem
cas_api = CASCommonChem()  # Automatically uses stored key
```

## Troubleshooting

### "API key is required" Error

If you get this error, it means no API key was found. The error message shows all available options:

```
API key is required for CAS Common Chemistry API v2.0.
Options:
1. Provide api_key parameter: CASCommonChem(api_key='your-key')
2. Provide api_key_file parameter: CASCommonChem(api_key_file='path/to/key.txt')
3. Set persistent API key: from provesid.config import set_cas_api_key; set_cas_api_key('your-key')
4. Set environment variable: CCC_API_KEY or CAS_API_KEY
```

### Check Configuration

```python
from provesid import show_config, get_cas_api_key

show_config()  # Shows configuration status
key = get_cas_api_key()  # Returns None if not configured
print(f"API key configured: {bool(key)}")
```

### Reset Configuration

To completely reset your configuration:

```python
from provesid import remove_cas_api_key
remove_cas_api_key()
```

Or manually delete the config file shown by `show_config()`.

## Examples

### Complete Setup Example

```python
# One-time setup
from provesid import set_cas_api_key, show_config

# Configure your API key
set_cas_api_key("your-cas-api-key-from-cas")

# Verify configuration
show_config()

# Now use normally
from provesid import CASCommonChem
cas = CASCommonChem()

# Test with a simple compound
water = cas.name_to_detail("water")
if water['found']:
    print(f"Water CAS RN: {water['rn']}")
    print(f"Molecular formula: {water['molecularFormula']}")
```

### Using with Different API Keys

```python
# Use stored API key
cas1 = CASCommonChem()

# Override with different API key for specific instance
cas2 = CASCommonChem(api_key="different-api-key")

# Use file-based key for another instance
cas3 = CASCommonChem(api_key_file="/path/to/another/key.txt")
```