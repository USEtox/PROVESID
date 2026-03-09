# Documentation Refresh Plan

This plan updates outdated docs to match current PROVESID functionality and aligns docs authoring with the Markdown-first examples workflow.

## Objectives

- Reflect current features across all API pages and guides.
- Remove stale references and duplicated legacy pages.
- Keep docs maintainable with executable local validation.

## Workstreams

## 1. Build Stability

- Keep CI docs build non-executed and strict.
- Fix malformed Markdown and nav/content drift issues.
- Ensure all files under `docs/` are intentional and maintained.

## 2. API Reference Consistency

- Review and update:
  - `docs/api/pubchem.md`
  - `docs/api/pubchemview.md`
  - `docs/api/chembl.md`
  - `docs/api/chebi.md`
  - `docs/api/cascommonchem.md`
  - `docs/api/classyfire.md`
  - `docs/api/opsin.md`
  - `docs/api/nci_resolver.md`
- Ensure method names and return shapes match `src/provesid/`.
- Use one consistent style across pages (manual narrative + auto-doc where reliable).

## 3. Guides Refresh

- Update:
  - `docs/index.md`
  - `docs/installation.md`
  - `docs/quickstart.md`
  - `docs/advanced_caching.md`
- Reflect offline-first behavior where applicable.
- Link guides to migrated tutorial Markdown pages.

## 4. Contributor Workflow

- Update `DEVELOPMENT.md` and `README.md` with:
  - Markdown-first tutorial authoring policy.
  - Local docs validation commands.
  - Optional notebook generation command.

## Execution Milestones

1. Stabilize strict docs build and remove stale files.
2. Complete tutorial format migration pilot.
3. Refresh API pages in two batches.
4. Refresh guides and contributor docs.
5. Final strict-build verification and content spot-check.

## Acceptance Criteria

- Docs build in strict mode without warnings/errors.
- Tutorials are maintained in text-first format.
- API docs align with current code and examples.
- Contributor workflow is documented and reproducible.
