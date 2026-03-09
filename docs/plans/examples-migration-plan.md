# Examples Migration Plan (Notebook -> MyST/Jupytext)

This plan makes `examples/` the single source of truth for tutorials and replaces committed `.ipynb` files with git-friendly Markdown sources.

## Goals

- Keep tutorial content executable and testable locally.
- Stop committing notebook JSON as primary source.
- Remove duplication between `examples/` and `docs/examples/`.
- Keep CI docs build fast (non-executed) while enabling local execution checks.

## Scope

### Current tutorial notebooks to migrate

- `examples/CCC/CAS_Common_Chemistry_tutorial.ipynb`
- `examples/ChEBI/ChEBI_tutorial.ipynb`
- `examples/ChEBI/chebi_sdf_tutorial.ipynb`
- `examples/ClassyFire/classyfire_tutorial.ipynb`
- `examples/OPSIN/opsin_tutorial.ipynb`
- `examples/pubchem/pubchem_tutorial.ipynb`
- `examples/resolver/chem_id_resolver_tutorial.ipynb`
- `examples/chembl/chembl_tutorial.ipynb`
- `examples/pubchemview/pubchem_view_tutorial.ipynb`
- `examples/zeropm/zeropm-example.ipynb`

### Duplicate docs copies to retire

- `docs/examples/CCC/CAS_Common_Chemistry_tutorial.ipynb`
- `docs/examples/ChEBI/ChEBI_tutorial.ipynb`
- `docs/examples/ClassyFire/classyfire_tutorial.ipynb`
- `docs/examples/OPSIN/opsin_tutorial.ipynb`
- `docs/examples/pubchem/pubchem_tutorial.ipynb`
- `docs/examples/resolver/chem_id_resolver_tutorial.ipynb`

## Migration Steps

1. Pilot convert two notebooks (`CCC` and `pubchem`) to MyST Markdown.
2. Validate round-trip export (`.md` -> `.ipynb`) for those pilots.
3. Update MkDocs nav for pilot tutorials from `.ipynb` to `.md`.
4. Convert the remaining eight tutorial notebooks.
5. Update all tutorial links in docs to `.md` targets.
6. Remove duplicated `docs/examples/` notebook copies.
7. Remove committed `.ipynb` tutorial sources from `examples/` after validation.
8. Generate `.ipynb` on demand (local command or CI artifact), not as source.

## Validation Checklist

- `mkdocs build --strict` passes.
- Every tutorial renders from Markdown in docs.
- Round-trip conversion succeeds for all migrated tutorials.
- Local execution check runs at least one happy path per tutorial.
- No tutorial notebook JSON files remain committed as source.

## Tooling

Use scripts added in `scripts/`:

- `scripts/convert_notebooks_to_myst.sh`
- `scripts/generate_notebooks_from_myst.sh`
