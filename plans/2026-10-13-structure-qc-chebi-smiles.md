# Structure quality control in `Search`: unparsable source SMILES, echoed support, silent conflicts

**Date:** 2026-10-13
**Author:** Ali A Eftekhari + Claude Code
**Reported from:** the structure audit of the 630 Sanofi CAS in the SANOFI USEtox project
(`SANOFI/structure_audit.py`, `SANOFI/data/structure_audit/README.md`)
**Scope:** `search.py` (`Search`, `_resolve_cas`, result finalisation, consensus), `tools.py`
(`make_candidate`, `candidate_from_chebi_row`, `apply_candidate_to_result`,
`compute_consensus`), `chebi_sdf.py`; downstream `proves.service.Service.resolve`
**Status:** Implemented on branch `structure-qc` (2026-10-13), except as noted:

- §3: done. The scan also found 19,228 ChEBI SMILES that parse but are not the record's
  compound (lost stereo, or lost ring H and still parse: guanine `Nc1nc(=O)c2ncnc2n1`), so a
  molfile with the stated key is preferred over those too. The `ChebiSDF` molfile was
  unreadable for 55% of records (blank title line skipped by the index) and is fixed on read.
  §3 step 4 (`_first_smiles_from_candidates`) is moot: ChEMBL is now asked by InChIKey.
- §4: done for CAS, including the per-cluster option (ChEMBL is asked for every InChIKey
  found). Not applied to InChI queries, where ChEMBL is asked for the query's own structure,
  which is not an echo. DTXSID queries had the same echo pattern (every source asked by
  CompTox's InChIKey; cefquinome's DTXSID came back wrong with support 2 and no conflict).
  Decided for correctness: the other sources are asked for CompTox's CAS number and vote
  from their own records; those without it, and ChEMBL, are asked by InChIKey and marked
  derived. A 1-vs-1 split goes to CompTox, flagged; a DTXSID whose CAS no other source
  holds has support 1 and `strict` refuses it.
- §5: done; `on_conflict="refuse"` in `strict` also refuses 2-vs-1 queries the majority gets
  right (tirapazamine, dalbavancin).
- §6: done, except the `substance_type` hint: no offline source states one.
- §7: done. `identity_consistent` compares InChIKey skeleton and protonation, not stereo.
- §8: `validate=` and `dataset_versions` done. The PROVES change lives in another repository
  and is not made here.
- §9: `scripts/qc_sources.py` and `tests/test_structure_qc.py` done.

---

## 1. Summary

**Why it matters.** A `Search` hit is used as a structure by everything downstream (PROVES,
USEtoxInput, EXProves), and nothing in the result says whether that structure is usable,
repaired, contested or echoed.

**How it was found.** We audited the 630 Sanofi CAS end to end. Six defects in provesid
surfaced; the first is the one that started the audit.

| # | Defect | Effect today (0.8.0) |
|---|---|---|
| A | ChEBI's own SMILES field is unparsable for 6.4 % of ChEBI (12 227 records, 1 825 with a CAS), almost all because the aromatic ring N–H is omitted. `Search` passes it through. | The result has `SMILES` set to an invalid string, `canonical_smiles = None`, and a normal-looking confidence (valsartan 0.82). 30 of 630 Sanofi CAS get no usable structure; `PROVES.resolve()` hands them on with `smiles = None`. |
| B | For CAS (and InChI) queries, ChEMBL is searched with the *first* candidate's SMILES and then counted as an independent vote and as source support. | Any wrong first hit gets a ChEMBL "confirmation". With the pre-0.6.0 ChEBI index bug this turned 8 wrong ChEBI hits into 2-vs-2 ties that ChEBI won on reputation. |
| C | Disagreement between sources is resolved silently: a losing cluster leaves no trace in the result. | Cefquinome sulfate (2 CAS) resolves to a wrong CompTox record although PubChem has the right one; the output shows `n_source_support = 1`, confidence 0.765 and no sign of a conflict. |
| D | A structure-less record (CompTox mixture / UVCB entries) cannot cluster with a structured record of the same substance, and wins the tie on reputation. | 9 Sanofi CAS (eprinomectin, gentamicin, neomycin, milbemycin oxime, …) return no structure at confidence 0.765, although PubChem has one. |
| E | Consistency between a source's stated InChIKey / formula / mass and its structure is checked only at `log.debug`, and the source InChIKey wins. | A result can carry an InChIKey that does not belong to its SMILES, or (A) an InChIKey with no SMILES at all. |
| F | `Search` returns unvalidated claims; validation (`resolve_cascade`, `mw_within`) is opt-in elsewhere, and outputs carry no dataset version. | PROVES uses `Search` directly. Lookups made before 0.6.0 (the 2026-05-04 Sanofi file, the 2026-06-27 ecotox QSAR inputs) cannot be traced to the buggy build from the files alone. |

**Proposal in one sentence.** Repair structures at candidate level (molfile → InChI →
N–H repair, accepted only if the InChIKey matches); never let an invalid SMILES occupy the
result; never let a structure-less record beat a structured one; count only independent sources; report conflicts and QC status as columns; and offer
validation inside `Search`.

---

## 2. Context: what the audit found downstream

The release structures of USEtoxInput were right: it resolves through the validated cascade
(`resolve_table` → `resolve_cascade` with an MW/name check). The damage came from places that
used `Search` hits as they came:

- **A wrong-structure lookup reached QSAR runs.**
  - The structures behind the ecotoxicity QSAR predictions (T.E.S.T. / VEGA, 2026-06-27) came from a lookup made before the 0.6.0 fixes (2026-08-04: ChEBI index offsets counted in text mode, and a lone ChEBI hit outranking a three-source consensus).
  - For 83 of the 630 Sanofi CAS, the effect QSAR was run on another compound: paracetamol as C9H9NO3, trimethoprim as C13H26O2, 5-fluorouracil as undecane, atomoxetine as azetidin-2-one, cefodizime as TCEP.
  - 39 of the shipped HC20s contain such a prediction.
  - That table had no column that would have flagged the rows.
- **Names from a pre-0.6.0 `Search` were mistaken for the data owner's names.** The file was `SANOFI/data/sanofi_extracted_ids.csv` (2026-05-04). It gave another compound for 12 CAS. The 8 of them sourced from ChEBI all showed the same pattern:
  - `source = ChEBI`, `source_match_scores = {'ChEBI': 1.0, 'ChEMBL': 1.0, 'CompTox': 0.0, 'PubChemID': 0.0}`, `confidence 0.675`;
  - two independent registries disagreed completely, and the result carried no flag.
- **provesid 0.8.0 still fails today, in different ways:**
  - 30 Sanofi CAS have an unusable structure (§3);
  - 2 resolve to a wrong CompTox record (§5);
  - 33 return no structure at all. 21 of these have a structure in the release: 12 found online through PubChem, and 9 that an earlier provesid resolved offline but 0.8.0's default no longer does (§6).

---

## 3. Defect A — unparsable source SMILES

### 3.1 Evidence

Valsartan, CAS 137862-53-4, CHEBI:9927. The `SMILES` field of the local `chebi.sdf` is:

    CCCCC(=O)N(Cc1ccc(-c2ccccc2-c2nnnn2)cc1)[C@H](C(=O)O)C(C)C

The tetrazole `c2nnnn2` has no `[nH]`, so RDKit cannot kekulise it. The SDF's own molfile and
InChI describe the right molecule.

    Search("cas").search(["137862-53-4"])
    # name valsartan | SMILES <the string above> | canonical_smiles None | InChIKey ACWBQPMHZXGDFX-QFIPXVFZSA-N
    # confidence 0.8175 | source ChEBI | source_match_scores {'ChEBI': 1.0, 'CompTox': 1.0, 'PubChemID': 0.45}

CompTox, which agrees, has a parsable structure: `Search(..., sources=["comptox"])` returns
`canonical_smiles = ...c2nnn[nH]2...`. It is not used, because ChEBI filled the field first.

### 3.2 Scope (local `chebi.sdf`, 2026-05-04 build, 192 479 records)

| | records |
|---|---|
| with a `SMILES` field | 192 445 |
| SMILES fails RDKit sanitisation | **12 227 (6.4 %)** |
| — kekulisation (aromatic N–H / ring H omitted) | 12 118 |
| — valence (radicals, hypervalent inorganics, e.g. `O=[Cl]=O`, `C#N=S`) | 108 |
| — of them with a CAS number | 1 825 |
| molfile parses | 12 128; InChIKey = ChEBI's for 11 727 |
| InChI parses | 11 797; InChIKey = ChEBI's for 11 782 |
| CAS-bearing kekulisation failures recovered from the molfile with ChEBI's InChIKey | 1 781 of 1 797 |

**Sanofi set: 30 of 630 CAS are affected, and the molfile recovers all 30 with the matching
InChIKey.** They are ordinary drugs:

- valsartan, irbesartan, lansoprazole, dexlansoprazole, omeprazole
- sumatriptan, zolmitriptan, naratriptan
- allopurinol, azathioprine, 5-fluorouracil, entecavir, aciclovir, methotrexate
- mebendazole, albendazole, oxfendazole
- theophylline, sildenafil, carvedilol, riboflavin phosphate
- daptomycin, insulin glargine, buserelin, dihydroergotamine, …

Mostly tetrazoles, imidazoles, purines and pyrimidinones.

The scan is reproducible with the script in the appendix. Re-run it on a current ChEBI
download before deciding how much repair is needed: a newer export may differ.

### 3.3 Mechanism

**1. The molfile is read but dropped.** `ChebiSDF._read_mol_at_offset` returns
`molfile` + properties. `candidate_from_chebi_row` passes only `row["SMILES"]`, `INCHI` and
`INCHIKEY` to `make_candidate`.

**2. The invalid SMILES wins the field.** `_finalise_result` applies candidates in priority
order with `apply_candidate_to_result`, which **never overwrites**. The consensus source's
invalid SMILES occupies `result["SMILES"]`, and a compatible CompTox / PubChem SMILES can no
longer fill it.

**3. Nothing marks the failure.** `normalize_structure(result["SMILES"])` returns all `None`.
`canonical_smiles`, `kekulized_smiles` and the RDKit InChIKey are empty, but `InChIKey` and
`name` come from the source, so the row looks resolved. The confidence is unaffected.

**4. The ChEMBL lookup fails as well.** For CAS queries, `_first_smiles_from_candidates`
returns that same invalid string, which is then used to search ChEMBL, so ChEMBL finds
nothing.

**5. PROVES keeps the row with no structure.** `proves.service.Service.resolve` builds its
anchor from `canonical_smiles` when the column exists. It drops only rows without an
InChIKey, so the chemical enters the pipeline with `smiles = None`.

### 3.4 Proposed fix (layered, each step useful alone)

1. **Repair at candidate level** in `make_candidate`, so every source benefits. Try, in order:
   1. the SMILES as given;
   2. the source molfile (ChEBI: pass `row["molfile"]` through `candidate_from_chebi_row`);
   3. the source InChI;
   4. an RDKit repair of aromatic N–H: parse with `sanitize=False`, then try adding one H to each aromatic N until `SanitizeMol` passes — the RDKit-cookbook "fix aromatic nitrogens" approach.

   **Accept a repaired structure only when its InChIKey equals the source's stated InChIKey**
   (full key, or the skeleton when the source states no stereo). Otherwise mark the candidate
   `structure_status = "invalid"` and give it no structure.

   For ChEBI, consider using the molfile *first*. The molfile is the curated record; the
   SMILES is a derived export.
2. **Record where the structure came from.** Store `structure_origin` (`smiles` / `molfile` / `inchi` / `nh_repair`) and the untouched source string (`source_smiles_raw`) on the candidate.
3. **Never place an invalid SMILES in the result.** `apply_candidate_to_result` should treat an unparsable SMILES as missing, so the next compatible source fills the field.
4. **`_first_smiles_from_candidates` returns the first *parsable* SMILES.**
5. **Report the status per row.** Add `structure_status`: `valid`, `repaired:<origin>`, `invalid` (source structure unusable and no fallback) or `none`. Rows with `invalid` should log at WARNING. `canonical_smiles = None` with a filled `InChIKey` should be impossible without such a flag.

---

## 4. Defect B — ChEMBL support is an echo, not a vote

### 4.1 Mechanism

`_resolve_cas` (and `_resolve_inchi`) works in two steps:

1. collects CAS hits from ChEBI, CompTox, PubChemID and ZeroPM;
2. then runs `self._collect("smiles", _first_smiles_from_candidates(hits), sources=["chembl"])`.

ChEMBL holds no CAS numbers, so its candidate is by construction a copy of the first source's
structure, ChEBI's whenever ChEBI answered. That candidate then:

- votes in `compute_consensus`;
- raises `n_source_support`, and with it `_SUPPORT_FACTOR` in `_compute_confidence`;
- fills fields in `_finalise_result`.

### 4.2 Evidence

All 8 wrong ChEBI-sourced rows of the pre-0.6.0 Sanofi lookup had ChEBI + ChEMBL = 1.0 against
CompTox + PubChem = 0.0. In the 2-vs-2 tie, the `_PRIORITY_TOLERANCE` reputation tie-break
picked ChEBI.

| CAS | Correct (CompTox, PubChem, release) | Returned (ChEBI + echoed ChEMBL) |
|---|---|---|
| 171500-79-1 | dalbavancin | empagliflozin |
| 864070-44-0 | empagliflozin | efinaconazole |
| 366789-02-8 | rivaroxaban | 5-aminonicotinic acid |
| 82248-59-7 | atomoxetine HCl | azetidin-2-one |
| 153832-46-3 | ertapenem | N-acetyl-Tyr-Gly-Gly |
| 69739-16-8 | cefodizime | TCEP |
| 168273-06-1 | rimonabant | "S9947" (another structure) |
| 343306-79-6 | sugammadex sodium | lumacaftor |

The ChEBI index bug that caused these is fixed. The echo is not: any future wrong first hit,
from any source, gets the same free corroboration.

### 4.3 Proposed fix

**Mark derived candidates and keep them out of the vote.** Give them `derived_from = "<source>"`
and exclude them from `compute_consensus` and `n_source_support`. Let them fill only fields the
independent sources lack (ChEMBL ID, synonyms, ATC).

**Optionally, query ChEMBL per cluster rather than per query.** Run the ChEMBL lookup with
every distinct candidate structure. ChEMBL then contributes cross-references to each cluster
without tipping the vote.

---

## 5. Defect C — conflicts are resolved silently

### 5.1 Evidence (provesid 0.8.0)

Cefquinome sulfate, CAS 118443-89-3 and 123766-80-3.

| source | structure |
|---|---|
| CompTox | `...C1=NSC(N)=N1...` — a 1,2,4-thiadiazole, parent C22H24N7O5S2⁺ |
| PubChem | cefquinome, with an aminothiazole, parent C23H24N6O5S2 |
| ChEBI | — |

The two land in different clusters. CompTox's single-source cluster outranks PubChem's on
base reputation.

    query 118443-89-3 | Cefquinome sulfate | source CompTox | confidence 0.765
    source_match_scores {'CompTox': 1.0} | n_source_support 1

Nothing in the row says that PubChem returned a different compound.

**Wrong registry records seen in the same audit.** Each is detectable only by cross-source
comparison:

| CAS | registry | problem |
|---|---|---|
| 27314-97-2 tirapazamine | CompTox | C8H7N3O2; correct is C7H6N4O2 (PubChem, ChEBI) |
| 164579-32-2 pantoprazole sodium | CompTox | a 2,3-*dihydro*benzimidazolide, C16H17F2N3O4S; correct is C16H15F2N3O4S |
| 68880-55-7 spiramycin adipate | CompTox | parent C44H76N2O14; spiramycin I is C43H74N2O14 (PubChem) |
| 171500-79-1 dalbavancin | CompTox | C87 component; PubChem and ChEBI give C88H100Cl2N10O28 |

### 5.2 Proposed fix

**Keep the evidence of a conflict in the result**, even with `n_hits=1`:

- `n_clusters`;
- `runner_up_inchikey`, `runner_up_sources` and `runner_up_confidence` for the best losing cluster;
- `conflict = True` when a losing cluster from an *independent* source has a different InChIKey skeleton (not a salt/stereo variant). "Independent" is in the sense of §4.

**Penalise the confidence of a contested result.** The current `_SUPPORT_FACTOR` penalises a
lone source but not a contested one. Multiply the confidence by a factor that falls with the
runner-up's support, so that 1-vs-1 and 2-vs-2 splits cannot reach the "trusted" range. When
the pick is made by the reputation tie-break, also set `conflict`.

**Add a "strict" option that refuses to choose.** Under `preset="strict"`, a conflicted query
should return no structure (or every cluster with `hit_rank`) rather than a winner.

---

## 6. Defect D — a structure-less record outranks a structure

### 6.1 Evidence (provesid 0.8.0, default sources)

    Search("cas").search(["123997-26-2", "1403-66-3"])
    # Eprinomectin | SMILES None | InChIKey None | DTXSID8037542 | CompTox | confidence 0.765 | n_source_support 1
    # Gentamicin   | SMILES None | InChIKey None | DTXSID5034642 | CompTox | confidence 0.765 | n_source_support 1

    Search("cas", n_hits="all").search(["123997-26-2"])
    # hit_rank 0  CompTox    0.765  no structure
    # hit_rank 1  PubChemID  0.765  C[C@H]1/C=C\C=C/2\CO[C@H]3...   (eprinomectin B1a)

**Why it happens.**

- CompTox represents these as mixtures / UVCBs, with a DTXSID but no structure.
- `_candidate_cluster_keys` keys a candidate by InChIKey, then SMILES, then name. A structured
  candidate gets InChIKey keys; a structure-less one gets only a name key. The two never
  intersect, so the same substance forms two clusters.
- They tie at 0.765, and the reputation order puts CompTox first.
- With `sources="all"`, ZeroPM corroborates PubChem and the structured cluster wins. The
  USEtoxInput release, resolved with an earlier provesid, has structures for all 9 Sanofi CAS
  in this position. Eprinomectin is the top-ranked API by freshwater ecotoxicity CF in that
  project.

**Affected Sanofi CAS (9):** eprinomectin, milbemycin oxime, gentamicin, neomycin, neomycin
sulfate, benzododecinium chloride, simethicone, streptokinase, expanded perlite. For the last
three, "no structure" is arguably the right answer; that is a decision for an explicit rule,
not for a tie-break.

### 6.2 Proposed fix

- **Merge across the structure gap.** Also emit the normalised-name key (and the DTXSID / CAS
  key) for structured candidates, so a structure-less record joins the structured cluster of
  the same substance. It then adds its DTXSID and support instead of competing.
- **Prefer structure in ties.** Among clusters within the priority tolerance, prefer one that
  carries a valid structure.
- **Make the mixture case explicit.** When the only answer is structure-less, return it with
  `structure_status = "none"` and a `substance_type` hint (mixture / UVCB / polymer / protein)
  where the source states one. Lower the confidence, because a result without a structure
  should not read as 0.765.

---

## 7. Defect E — identity consistency is debug-only

In `_finalise_result`, a source InChIKey that differs from the InChIKey RDKit derives from the
chosen SMILES only produces `log.debug("InChIKey mismatch ...")`. `pick_first` keeps the source
key. Formula and mass have the same gap: `molecular_mass` comes from the first source that
states one, whatever structure it describes.

**Proposal:**

- **Make the structure authoritative.** When a valid structure exists, derive `InChIKey`, `molecular_formula` and `molecular_mass` from it. Keep the source values only in `source_inchikey` / `source_formula`.
- **Add `identity_consistent`.** It checks structure ↔ stated InChIKey ↔ stated formula. A mismatch logs at WARNING and becomes a QC flag.

---

## 8. Defect F — validation and provenance live outside `Search`

**Validation is opt-in elsewhere.** `resolve_cascade` and `mw_within` validate hits against a
reference MW or name, but only callers that know about them use them. `PROVES.resolve()` and
ad-hoc scripts (such as the one that produced the Sanofi file) call `Search` directly.

**Proposals:**

- **A `validate` option in `Search`.** For example `validate={"mw": col, "formula": col, "name": col}` with the `mw_within` semantics: reject a hit that fails and fall through to the next cluster. Record `validated_by` as `resolve_cascade` does.
- **PROVES.** Use `resolve_cascade` (or `validate=`) and drop or flag rows whose `structure_status` is not `valid`. Today it keeps any row with an InChIKey.
- **Dataset provenance in `df.attrs`.** Add `dataset_versions`: per source, the file date/size or hash (for ChEBI, the `_index_meta` values and the SDF date) plus the provesid version. The 2026-05-04 Sanofi lookup could then be traced to the pre-fix index without guessing from file dates.

---

## 9. Quality control: a source-level QC report and a regression set

### 9.1 Source-level QC script

Add a `scripts/qc_sources.py`, run after every dataset fetch, that writes one table per source:

- **Structure parsing:** structures that do not parse, and which fallback recovers them (the appendix scan, generalised).
- **Self-consistency:** stated InChIKey / formula / mass vs the structure.
- **Ambiguous CAS:** CAS numbers mapping to more than one structure skeleton within a source.
- **Cross-source conflicts:** for a sample of CAS, how often independent sources disagree.

Publishing the counts in the CHANGELOG with each dataset bump would show regressions such as
the byte-offset bug immediately.

### 9.2 Regression tests

Each case is from this audit. Expected values are verified against PubChem and the
USEtoxInput release.

| CAS | expected parent | guards against |
|---|---|---|
| 137862-53-4 valsartan | C24H29N5O3, skeleton ACWBQPMHZXGDFX | A (tetrazole N–H) |
| 315-30-0 allopurinol | C5H4N4O | A |
| 446-86-6 azathioprine | C9H7N7O2S | A |
| 51-21-8 5-fluorouracil | C4H3FN2O2 | A |
| 103577-45-3 lansoprazole | C16H14F3N3O2S | A |
| 118443-89-3 cefquinome sulfate | C23H24N6O5S2 (or `conflict = True`) | C (wrong CompTox record) |
| 27314-97-2 tirapazamine | C7H6N4O2 | C (wrong CompTox formula) |
| 171500-79-1 dalbavancin | C88H100Cl2N10O28 | B, index regressions |
| 864070-44-0 empagliflozin | C23H27ClO7 | B, index regressions |
| 69739-16-8 cefodizime | C20H20N6O7S4 | B, index regressions |
| 82248-59-7 atomoxetine HCl | C17H21NO | B, index regressions |
| 343306-79-6 sugammadex sodium | the γ-cyclodextrin derivative, not lumacaftor | B, index regressions |
| 123997-26-2 eprinomectin | eprinomectin B1a, C50H75NO14 (not a structure-less hit) | D |
| 1403-66-3 gentamicin | a gentamicin C component (not a structure-less hit) | D |

**Assertions:**

- the parent formula / skeleton;
- `structure_status == "valid"` or `"repaired:molfile"`;
- `canonical_smiles` is not None whenever `InChIKey` is set;
- `conflict` is set for cefquinome until CompTox fixes its record.

---

## 10. Suggested order

1. **§3 A-fix steps 1–4:** the molfile fallback, and no invalid SMILES in the result. Smallest change, largest gain: about 1 800 CAS-bearing ChEBI records get a structure.
2. **§3 step 5 and §7:** the `structure_status` and `identity_consistent` columns. Cheap, and makes every later change visible.
3. **§4:** independent-source counting for ChEMBL.
4. **§5 and §6:** conflict reporting, the confidence penalty and structure-preferring clustering.
5. **§8:** the `validate=` option, PROVES switching to validated resolution, and `dataset_versions`.
6. **§9:** the QC script and the regression set (start the regression tests with step 1).

---

## Appendix — reproduce the ChEBI scan

```python
"""Which ChEBI SMILES fail RDKit, and does the molfile or the InChI recover them?"""
import csv, sys
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

def records(path):
    buf = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("$$$$"):
                yield "".join(buf); buf = []
            else:
                buf.append(line)

def prop(block, name):
    i = block.find(f"> <{name}>")
    if i < 0:
        return None
    j = block.find("\n", i); k = block.find("\n\n", j)
    return block[j + 1:k if k > 0 else None].strip() or None

rows = []
for block in records(sys.argv[1]):              # ~/.local/share/provesid/chebi.sdf
    smi = prop(block, "SMILES")
    if not smi or Chem.MolFromSmiles(smi) is not None:
        continue
    ik = prop(block, "INCHIKEY")
    m1 = Chem.MolFromMolBlock(block[: block.find("M  END") + 6]) if "M  END" in block else None
    m2 = Chem.MolFromInchi(prop(block, "INCHI") or "")
    rows.append([prop(block, "ChEBI ID"), prop(block, "ChEBI NAME"),
                 prop(block, "CAS Registry Numbers"), smi, ik,
                 m1 is not None and Chem.MolToInchiKey(m1) == ik,
                 m2 is not None and Chem.MolToInchiKey(m2) == ik])
with open(sys.argv[2], "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["chebi_id", "name", "cas", "smiles", "inchikey",
                "molfile_recovers", "inchi_recovers"])
    w.writerows(rows)
print(len(rows), "unparsable SMILES")
```
