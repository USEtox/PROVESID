# Data sources

## Zero PM

The ZeroPM data comes from the ZeroPM global inventory repository.  

## Chebi

The ChEBI data comes from [here](https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/).

**Note:** Large offline datasets (including ChEBI SDF) are downloaded on first
use into a shared per-user PROVESID dataset directory (platform-specific via
`platformdirs`). Set `PROVESID_DATA_DIR` to override the default location.
The ChEBI index file (`chebi.sdf.index.pkl`) is also automatically created on
first use for fast queries.

Both files are excluded from version control via .gitignore.  