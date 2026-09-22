# API keys

Only one service PROVESID talks to needs a key: **CAS Common Chemistry**.
Request one from CAS, then store it once:

```python
from provesid import set_cas_api_key

set_cas_api_key("your-cas-api-key")
```

and every later `CASCommonChem()` picks it up:

```python
from provesid import CASCommonChem

cas = CASCommonChem()
water = cas.cas_to_detail("7732-18-5")
water["found"], water["name"], water["molecularFormula"]
```

`cas_to_detail`, `name_to_detail` and `smiles_to_detail` do not raise: a
failure is reported in `found` and `status`, and a rejected key reads
`status == "Unauthorized - Check API Key"`.

## Where the key comes from

`CASCommonChem` takes the first of these that yields a key:

1. `CASCommonChem(api_key="...")`
2. `CASCommonChem(api_key_file="path/to/key.txt")`
3. the stored key, from `set_cas_api_key`
4. the `CCC_API_KEY` or `CAS_API_KEY` environment variable

A stored key outranks the environment. If you set `CAS_API_KEY` and it is being
ignored, a stored key — perhaps a placeholder — is the likely reason:
`remove_cas_api_key()` deletes it. With no key from any of them, the
constructor raises an error listing these four options.

## Where it is stored

```python
from provesid import show_config, get_cas_api_key, remove_cas_api_key

show_config()           # prints the config directory, file, and which services have keys
get_cas_api_key()       # the stored key, or None
remove_cas_api_key()
```

The key is stored in plain text in `config.json`, in
`$XDG_CONFIG_HOME/provesid` (by default `~/.config/provesid`) on Linux and
macOS and `%APPDATA%\PROVESID` on Windows. The file gets your user account's
default permissions; nothing restricts it further. On a shared machine, prefer
the environment variable or a key file only you can read.
