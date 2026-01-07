# Data sources

## Zero PM

The ZeroPM data comes from the ZeroPM global inventory repository.  

## Chebi

The ChEBI data comes from [here](https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/).

**Note:** The ChEBI SDF file (`chebi.sdf`, ~868 MB) is automatically downloaded on first use of the `ChebiSDF` class. The file is downloaded as a gzip archive (~250 MB) and automatically extracted. An index file (`chebi.sdf.index.pkl`) is also automatically created on first use for fast queries.

Both files are excluded from version control via .gitignore.  