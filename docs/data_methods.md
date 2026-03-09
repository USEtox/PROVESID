# Online and Offline Data Methods

This page explains how PROVESID separates online API access from local/offline dataset access, and how to choose the right method for each workflow.

## Why this distinction matters

PROVESID supports two execution styles:

- Online methods call remote services and return current data from external APIs.
- Offline methods query local datasets and databases for fast, reproducible lookup.

In many workflows, a practical pattern is:

1. Use local/offline lookup first for speed and reproducibility.
2. Use online services for missing records, richer metadata, or live updates.

## Online methods

These classes primarily use network APIs:

- `PubChemAPI`
- `PubChemView`
- `NCIChemicalIdentifierResolver`
- `CASCommonChem`
- `ChEBI`
- `OPSIN` and `PYOPSIN`
- `ClassyFireAPI`

### Online method characteristics

- Requires internet connectivity.
- Subject to remote service availability and response-time variance.
- May be rate-limited depending on provider.
- Usually provides the latest upstream data.

## Offline and local dataset methods

These classes read local data files/databases:

- `CheMBL` (local SQLite)
- `PubChemID` (local SQLite identifier database)
- `CompToxID` (local SQLite)
- `ZeroPM` (local SQLite)
- `REACHDossierID` (local REACH dataset)
- `ChebiSDF` (local SDF parsing)

### Offline method characteristics

- Fast and stable lookup once data is available locally.
- Better reproducibility for pipelines and batch processing.
- Large datasets may require significant storage.
- Some classes can auto-download missing datasets.

## Storage and environment strategy

Because some local datasets are large, `uv` is the recommended installation workflow:

```bash
uv pip install provesid
```

For development:

```bash
git clone https://github.com/USEtox/PROVESID.git
cd PROVESID
uv pip install -e .
```

`uv` helps avoid repeated copies across multiple virtual environments and keeps dependency management efficient while you work with large local assets.

## Practical selection guide

Use online methods when:

- you need live upstream updates,
- you need fields not present in local snapshots,
- or a local dataset for that source is not available.

Use offline methods when:

- you process many records,
- you need predictable/reproducible runs,
- you work in constrained or intermittent network environments,
- or you want lower latency per lookup.

## Related pages

- [Quick Start](quickstart.md)
- [Advanced Caching](advanced_caching.md)
- [API Overview](api/index.md)
