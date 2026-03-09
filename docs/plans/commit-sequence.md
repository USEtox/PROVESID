# Suggested Commit Sequence

This sequence keeps review clean and separates CI changes, docs fixes, migration tooling, and full tutorial conversion.

## Commit 1 - Disable test workflows

```bash
git add .github/workflows/test.yml .github/workflows/test-with-api-keys.yml
git commit -m "ci: disable GitHub test workflows during active development"
```

## Commit 2 - Stabilize strict docs build

```bash
git add docs/api/index.md docs/api/pubchem.md docs/api/pubchemview.md docs/api/nci_resolver.md docs/index.md docs/quickstart.md mkdocs.yml src/provesid/chembl.py docs/api/API_Documentation.md docs/api/PUBCHEM_API_DOCUMENTATION.md
git commit -m "docs: fix strict build issues and remove stale API pages"
```

## Commit 3 - Add migration/validation tooling

```bash
git add pyproject.toml uv.lock DEVELOPMENT.md scripts/convert_notebooks_to_myst.sh scripts/generate_notebooks_from_myst.sh scripts/validate_docs_local.sh docs/plans/examples-migration-plan.md docs/plans/docs-refresh-plan.md docs/plans/commit-sequence.md
git commit -m "docs: add markdown tutorial migration and local validation tooling"
```

## Commit 4 - Convert examples to Markdown sources and remove notebooks

```bash
git add examples/CCC/CAS_Common_Chemistry_tutorial.md examples/ChEBI/ChEBI_tutorial.md examples/ChEBI/chebi_sdf_tutorial.md examples/ClassyFire/classyfire_tutorial.md examples/OPSIN/opsin_tutorial.md examples/chembl/chembl_tutorial.md examples/pubchem/pubchem_tutorial.md examples/pubchemview/pubchem_view_tutorial.md examples/resolver/chem_id_resolver_tutorial.md examples/zeropm/zeropm-example.md examples/CCC/CAS_Common_Chemistry_tutorial.ipynb examples/ChEBI/ChEBI_tutorial.ipynb examples/ChEBI/chebi_sdf_tutorial.ipynb examples/ClassyFire/classyfire_tutorial.ipynb examples/OPSIN/opsin_tutorial.ipynb examples/chembl/chembl_tutorial.ipynb examples/pubchem/pubchem_tutorial.ipynb examples/pubchemview/pubchem_view_tutorial.ipynb examples/resolver/chem_id_resolver_tutorial.ipynb examples/zeropm/zeropm-example.ipynb
git commit -m "examples: migrate tutorials to MyST markdown and remove committed notebooks"
```

## Commit 5 - Deduplicate docs examples using symlinks

```bash
git add docs/examples/CCC docs/examples/ChEBI docs/examples/ClassyFire docs/examples/OPSIN docs/examples/chembl docs/examples/pubchem docs/examples/pubchemview docs/examples/resolver docs/examples/zeropm docs/examples/CCC/CAS_Common_Chemistry_tutorial.ipynb docs/examples/ChEBI/ChEBI_tutorial.ipynb docs/examples/ClassyFire/classyfire_tutorial.ipynb docs/examples/OPSIN/opsin_tutorial.ipynb docs/examples/pubchem/pubchem_tutorial.ipynb docs/examples/resolver/chem_id_resolver_tutorial.ipynb
git commit -m "docs: deduplicate tutorials by linking docs/examples to examples sources"
```

## Optional Commit - Existing README changes

You already had local `README.md` changes before this migration work. Keep them separate:

```bash
git add README.md
git commit -m "docs: update README project description and install notes"
```

## Final Verification

```bash
uv run --extra docs mkdocs build --strict
./scripts/validate_docs_local.sh
```
