"""Quality-control report for the offline sources behind ``Search``.

Run after every dataset fetch.  A source can change under the resolver without
any code changing --- a new ChEBI export, a CompTox refresh --- and the 0.6.0
ChEBI index bug showed that such a change is invisible until someone audits
results by hand.  This writes one table per check and source, and a summary
whose counts can go into the CHANGELOG with the dataset bump:

``<source>_structure_parsing.csv``
    Records whose SMILES RDKit cannot read, why (``kekulize``, ``valence``,
    ``other``), and which fallback rebuilds it with the source's own InChIKey
    (``molfile``, ``inchi``, ``nh_repair``, or none) --- what
    ``provesid.tools.repair_structure`` does at query time.
``<source>_self_consistency.csv``
    Records whose structure does not have the InChIKey skeleton, the formula,
    or (within 0.5 Da) the mass the same record states.
``<source>_ambiguous_cas.csv``
    CAS numbers the source gives to more than one InChIKey skeleton.
``cross_source_conflicts.csv``
    For a random sample of CompTox CAS numbers, the ``Search`` rows where
    independent sources return different compounds (``conflict``).
``summary.json``
    The counts, per source and check.

Sources are read directly --- ChEBI's SDF record by record, the others' SQLite
tables --- so the report covers every record, not just the ones a query
happens to reach.  ChEBI takes a few minutes and CompTox about half an hour;
``--limit`` caps the records read per source for a quick look.  PubChem is
read only for its CAS-bearing compounds, and is the slowest.

Usage::

    uv run python scripts/qc_sources.py --out qc/2026-10-13
    uv run python scripts/qc_sources.py --sources chebi --out qc/chebi
    uv run python scripts/qc_sources.py --sources comptox --limit 50000 --cas-sample 0
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import random
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

from provesid.datasets import data_directory, dataset_files
from provesid.tools import (
    extract_cas_values,
    formula_element_counts,
    is_missing,
    repair_structure,
    structure_formula,
    to_float,
)

RDLogger.DisableLog("rdApp.*")
log = logging.getLogger("qc_sources")

SOURCES = ("chebi", "comptox", "pubchem", "zeropm")


# ── Readers: each yields one dict per record ─────────────────────────────────
#
# Keys: id, cas (list), smiles, molfile, inchi, inchikey, formula, mass.  A
# reader leaves out what its source does not state.


def _sdf_property(block: str, name: str) -> Optional[str]:
    """The first line of an SDF data item, or None."""
    start = block.find(f"> <{name}>")
    if start < 0:
        return None
    value_start = block.find("\n", start) + 1
    value_end = block.find("\n", value_start)
    value = block[value_start:value_end if value_end > 0 else None].strip()
    return value or None


def read_chebi(path: str) -> Iterator[Dict[str, Any]]:
    """Every record of ChEBI's SDF, molfile included."""
    lines: List[str] = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("$$$$"):
                lines.append(line)
                continue
            block = "".join(lines)
            lines = []
            end = block.find("M  END")
            yield {
                "id": _sdf_property(block, "ChEBI ID"),
                "name": _sdf_property(block, "ChEBI NAME"),
                "cas": extract_cas_values(_sdf_property(block, "CAS Registry Numbers")),
                "smiles": _sdf_property(block, "SMILES"),
                "molfile": block[: end + 6] if end >= 0 else None,
                "inchi": _sdf_property(block, "INCHI"),
                "inchikey": _sdf_property(block, "INCHIKEY"),
                "formula": _sdf_property(block, "FORMULA"),
            }


def read_comptox(path: str) -> Iterator[Dict[str, Any]]:
    """Every substance in the CompTox table that has a structure."""
    db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    rows = db.execute(
        "SELECT DTXSID, PREFERRED_NAME, CASRN, SMILES, INCHIKEY, MOLECULAR_FORMULA, "
        "AVERAGE_MASS FROM chemicals WHERE SMILES IS NOT NULL"
    )
    for dtxsid, name, casrn, smiles, inchikey, formula, mass in rows:
        yield {
            "id": dtxsid, "name": name, "cas": extract_cas_values(casrn),
            "smiles": smiles, "inchikey": inchikey, "formula": formula, "mass": mass,
        }


def read_pubchem(path: str) -> Iterator[Dict[str, Any]]:
    """The PubChem compounds that carry a CAS number, with their numbers."""
    db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    rows = db.execute(
        "SELECT c.cid, c.cmpdname, group_concat(n.cas, ' '), c.smiles, c.inchi, "
        "c.inchikey, c.mf, c.mw FROM compounds c JOIN cas_numbers n ON n.cid = c.cid "
        "GROUP BY c.cid"
    )
    for cid, name, cas, smiles, inchi, inchikey, formula, mass in rows:
        yield {
            "id": cid, "name": name, "cas": extract_cas_values(cas), "smiles": smiles,
            "inchi": inchi, "inchikey": inchikey, "formula": formula, "mass": mass,
        }


def read_zeropm(path: str) -> Iterator[Dict[str, Any]]:
    """ZeroPM's substances: an InChI and the InChIKey it states for it."""
    db = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    for inchi_id, inchi, inchikey in db.execute("SELECT inchi_id, inchi, inchikey FROM substances"):
        mol = Chem.MolFromInchi(inchi) if inchi else None
        yield {
            "id": inchi_id, "cas": [], "inchi": inchi, "inchikey": inchikey,
            # The structure is the InChI: one that does not parse is reported
            # as an unparsable structure, like a SMILES.
            "smiles": Chem.MolToSmiles(mol) if mol is not None else (inchi or None),
        }


READERS = {"chebi": read_chebi, "comptox": read_comptox, "pubchem": read_pubchem, "zeropm": read_zeropm}


# ── Checks ───────────────────────────────────────────────────────────────────


def failure_kind(smiles: str) -> str:
    """Why RDKit rejects a SMILES: ``kekulize``, ``valence``, ``syntax`` or ``other``."""
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return "syntax"
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:  # RDKit raises several sanitisation exception types
        text = str(exc).lower()
        if "kekulize" in text:
            return "kekulize"
        if "valence" in text:
            return "valence"
    return "other"


def check_source(records: Iterator[Dict[str, Any]], limit: Optional[int]) -> Dict[str, Any]:
    """Run the per-record checks over one source.

    Returns:
        ``{"parsing": rows, "consistency": rows, "ambiguous_cas": rows,
        "counts": {...}}``.  ``parsing`` also lists SMILES that parse but
        were replaced by the molfile because they are not the compound the
        record's InChIKey names (``failure`` = ``wrong_structure``).
    """
    counts: Counter = Counter()
    parsing: List[Dict[str, Any]] = []
    consistency: List[Dict[str, Any]] = []
    skeletons_by_cas: Dict[str, set] = defaultdict(set)

    for n, rec in enumerate(records):
        if limit is not None and n >= limit:
            break
        counts["records"] += 1
        if rec.get("cas"):
            counts["with_cas"] += 1
        raw = rec.get("smiles")
        if is_missing(raw):
            counts["no_structure"] += 1
            continue
        counts["with_structure"] += 1

        smiles, origin, status = repair_structure(
            raw, molfile=rec.get("molfile"), inchi=rec.get("inchi"), inchikey=rec.get("inchikey")
        )
        if status != "valid":
            kind = "wrong_structure" if Chem.MolFromSmiles(raw) is not None else failure_kind(raw)
            counts["unparsable"] += 1
            counts[f"unparsable_{kind}"] += 1
            counts[f"recovered_by_{origin}" if origin else "unrecovered"] += 1
            if rec.get("cas"):
                counts["unparsable_with_cas"] += 1
                if origin is None:
                    counts["unrecovered_with_cas"] += 1
            parsing.append({
                "id": rec["id"], "name": rec.get("name"), "cas": "; ".join(rec.get("cas") or []),
                "smiles": raw, "inchikey": rec.get("inchikey"), "failure": kind,
                "recovered_by": origin or "", "recovered_smiles": smiles or "",
            })

        mol = Chem.MolFromSmiles(smiles) if smiles else None
        structure_key = Chem.MolToInchiKey(mol) if mol is not None else None
        stated_key = rec.get("inchikey")
        skeleton = (structure_key or stated_key or "")[:14]
        for cas in rec.get("cas") or []:
            if skeleton:
                skeletons_by_cas[cas].add(skeleton)

        if mol is None:
            continue
        problems = []
        if not is_missing(stated_key) and (
            structure_key[:14] != str(stated_key)[:14] or structure_key[-1] != str(stated_key)[-1]
        ):
            problems.append("inchikey")
        stated_formula = formula_element_counts(rec.get("formula"))
        computed_formula = structure_formula(smiles)
        computed_counts = formula_element_counts(computed_formula)
        if stated_formula is not None and computed_counts is not None and stated_formula != computed_counts:
            problems.append("formula")
        stated_mass = to_float(rec.get("mass"))
        computed_mass = Descriptors.MolWt(mol) if stated_mass is not None else None
        if stated_mass is not None and abs(stated_mass - computed_mass) > 0.5:
            problems.append("mass")
        for problem in problems:
            counts[f"inconsistent_{problem}"] += 1
        if problems:
            counts["inconsistent"] += 1
            consistency.append({
                "id": rec["id"], "name": rec.get("name"), "cas": "; ".join(rec.get("cas") or []),
                "problems": "+".join(problems), "smiles": smiles,
                "stated_inchikey": stated_key, "structure_inchikey": structure_key,
                "stated_formula": rec.get("formula"), "structure_formula": computed_formula,
                "stated_mass": stated_mass,
                "structure_mass": round(computed_mass, 3) if computed_mass is not None else None,
            })

    ambiguous = [
        {"cas": cas, "n_skeletons": len(skeletons), "skeletons": "; ".join(sorted(skeletons))}
        for cas, skeletons in sorted(skeletons_by_cas.items()) if len(skeletons) > 1
    ]
    counts["ambiguous_cas"] = len(ambiguous)
    return {"parsing": parsing, "consistency": consistency, "ambiguous_cas": ambiguous,
            "counts": dict(counts)}


def cross_source_conflicts(n: int, seed: int) -> Dict[str, Any]:
    """Resolve a random sample of CompTox CAS numbers and keep the conflicts."""
    from provesid import Search
    from provesid.comptox import CompToxID

    with CompToxID() as comptox:
        cas_numbers = [
            row[0] for row in comptox.conn.execute(
                "SELECT CASRN FROM chemicals WHERE CASRN IS NOT NULL"
            )
        ]
    sample = random.Random(seed).sample(cas_numbers, min(n, len(cas_numbers)))
    with Search("cas", show_progress=True) as search:
        frame = search.search(sample)
    conflicts = frame[frame["conflict"].fillna(False).astype(bool)]
    columns = ["query", "name", "InChIKey", "source", "n_source_support", "confidence",
               "runner_up_inchikey", "runner_up_sources", "runner_up_confidence",
               "structure_status", "identity_consistent"]
    counts = {
        "sampled": len(frame),
        "conflict": int(len(conflicts)),
        "structure_status": frame["structure_status"].value_counts().to_dict(),
        "identity_inconsistent": int((frame["identity_consistent"] == False).sum()),  # noqa: E712
    }
    return {"rows": conflicts[columns].to_dict("records"), "counts": counts}


# ── Output ───────────────────────────────────────────────────────────────────


def write_rows(path: Path, rows: List[Dict[str, Any]]) -> None:
    """Write dict rows as CSV; an empty table still gets its file."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        if not rows:
            return
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sources", default="chebi,comptox",
                        help=f"comma-separated, of {', '.join(SOURCES)} (default: chebi,comptox)")
    parser.add_argument("--out", default="qc", help="output directory (default: ./qc)")
    parser.add_argument("--limit", type=int, default=None, help="records to read per source")
    parser.add_argument("--cas-sample", type=int, default=500,
                        help="CompTox CAS numbers to resolve for the cross-source check (0 skips it)")
    parser.add_argument("--seed", type=int, default=13, help="seed for the CAS sample")
    parser.add_argument("--data-dir", default=None, help="dataset directory (default: provesid's)")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    directory = data_directory(args.data_dir)
    summary: Dict[str, Any] = {"data_dir": directory}

    for source in [s.strip() for s in args.sources.split(",") if s.strip()]:
        if source not in READERS:
            parser.error(f"unknown source {source!r}; choose from {', '.join(SOURCES)}")
        files = dataset_files(source, directory, include_extras=False)
        if not files:
            log.warning("%s is not installed in %s; skipped.", source, directory)
            continue
        path = files[0]
        log.info("%s: reading %s", source, path)
        result = check_source(READERS[source](path), args.limit)
        write_rows(out / f"{source}_structure_parsing.csv", result["parsing"])
        write_rows(out / f"{source}_self_consistency.csv", result["consistency"])
        write_rows(out / f"{source}_ambiguous_cas.csv", result["ambiguous_cas"])
        summary[source] = {"file": path, **result["counts"]}
        log.info("%s: %s", source, json.dumps(result["counts"], sort_keys=True))

    if args.cas_sample > 0:
        log.info("Resolving %d CompTox CAS numbers for the cross-source check", args.cas_sample)
        result = cross_source_conflicts(args.cas_sample, args.seed)
        write_rows(out / "cross_source_conflicts.csv", result["rows"])
        summary["cross_source"] = result["counts"]
        log.info("cross-source: %s", json.dumps(result["counts"], sort_keys=True, default=str))

    with open(out / "summary.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True, default=str)
    log.info("Wrote %s", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
