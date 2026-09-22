"""The sources :class:`~provesid.search.Search` queries, as a table.

Each lookup takes one source client and one :class:`Query` and returns that
source's candidates for it, best first, already adapted by the ``tools``
``candidate_from_*`` helpers.  A lookup that finds nothing returns an empty
list.  Lookups do not catch exceptions, bar one: an online service's "not
found" is a miss rather than a failure, so the online lookups turn
:class:`~provesid.http.NotFoundError` into an empty list.  The one driver
that calls them, ``Search._collect``, logs a failing source and carries on
with the rest, so that policy is written once instead of once per rung.

The table is keyed by *lookup kind* first and source second.  The first five
sources are offline databases; the last two are web services, asked only
when ``Search(online_fallback=True)`` and no offline source answered:

======================  ======  =======  =======  ======  ======  ==============  ======
kind                    chebi   comptox  pubchem  zeropm  chembl  pubchem_online  cactus
======================  ======  =======  =======  ======  ======  ==============  ======
``cas``                 yes     yes      yes      yes     --      yes             yes
``inchikey``            yes     yes      yes      yes     yes     yes             yes
``inchikey_skeleton``   yes     yes      yes      --      --      --              --
``inchi``               yes     --       yes      yes     --      yes             yes
``smiles``              --      yes      yes      yes     yes     yes             yes
``dtxsid``              --      yes      --       --      --      yes             --
``name``                yes     yes      yes      yes     yes     yes             yes
``fuzzy_name``          yes     yes      yes      yes     yes     --              --
``formula``             yes     yes      yes      --      --      --              --
======================  ======  =======  =======  ======  ======  ==============  ======

A gap means the source has no index for that identifier.  ``Search`` reaches
it through an identifier it does have instead: ChEMBL through the SMILES
another source found for a CAS number, ChEBI through the InChIKey of a SMILES
query, and so on.  Those routes are the resolver's business, since they
depend on what the other sources answered, so they live in ``search.py``.

The online gaps are deliberate.  Neither service has a fuzzy or prefix search
worth a network round trip, and a formula names thousands of PubChem
compounds.  PubChem reaches a DTXSID through its synonyms, where EPA's DSSTox
deposit puts it; CACTUS does not know DTXSIDs at all.

Before this module existed, ``search.py`` held each row of this table as a
hand-written ``if client is not None: try: ... except: log`` block, 47 of
them in nine methods.  Adding a source meant editing every one.
It now means adding a column here.

Example::

    from provesid.sources import LOOKUPS, Query

    lookup = LOOKUPS["cas"]["pubchem"]
    candidates = lookup(pubchem_id, Query("50-78-2"))
    candidates[0]["InChIKey"]   # 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N'
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional

from .http import NotFoundError
from .tools import (
    candidate_from_cactus,
    candidate_from_chebi_row,
    candidate_from_chembl_row,
    candidate_from_comptox_row,
    candidate_from_pubchem_online,
    candidate_from_pubchem_row,
    candidate_from_zeropm_name_table,
    candidate_from_zeropm_smiles,
    is_missing,
)

log = logging.getLogger(__name__)

Candidate = Dict[str, Any]

#: Every offline source key, in the order results are pooled and reported.
SOURCE_KEYS: List[str] = ["chebi", "comptox", "pubchem", "zeropm", "chembl"]

#: The web services ``Search(online_fallback=True)`` asks when every offline
#: source missed, pooled and reported after the offline ones.
ONLINE_SOURCE_KEYS: List[str] = ["pubchem_online", "cactus"]

#: Display names, as they appear in ``source_details`` and in log lines.
SOURCE_DISPLAY: Dict[str, str] = {
    "chebi": "ChEBI",
    "comptox": "CompTox",
    "pubchem": "PubChemID",
    "zeropm": "ZeroPM",
    "chembl": "ChEMBL",
    "pubchem_online": "PubChem (online)",
    "cactus": "CACTUS",
}


@dataclass(frozen=True)
class Query:
    """One identifier to look up, plus what the lookups need around it.

    Attributes:
        value: The identifier as the source is asked for it.
        label: The name a ZeroPM candidate is given.  ZeroPM answers with a
            table rather than a row, and the candidate built from it is named
            after the user's query, which is not always ``value``: a DTXSID
            query reaches ZeroPM by the InChIKey CompTox found, and the
            candidate still carries the DTXSID.  Defaults to ``value``.
        k: How many candidates to take from each source.  ``1`` for the
            identifier lookups, ``top_k_per_source`` for names and formulas.
        fuzzy_cutoff: Score cut-off in [0, 100] that ZeroPM's fuzzy name
            retrieval applies; only ``fuzzy_name`` reads it.

    Example:
        >>> Query("50-78-2").label
        '50-78-2'
        >>> Query("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", label="DTXSID5020108").label
        'DTXSID5020108'
    """

    value: str
    label: Optional[str] = None
    k: int = 1
    fuzzy_cutoff: float = 80.0

    def __post_init__(self) -> None:
        """Default ``label`` to ``value``."""
        if self.label is None:
            object.__setattr__(self, "label", self.value)


Lookup = Callable[[Any, Query], List[Candidate]]


# ── Adapting what a client returns ────────────────────────────────────────────


def _one(adapt: Callable[[Dict[str, Any]], Candidate], row: Optional[Dict[str, Any]]) -> List[Candidate]:
    """Adapt a single-row answer, which is ``None`` or empty on a miss."""
    return [adapt(row)] if row else []


def _top(
    adapt: Callable[[Dict[str, Any]], Candidate],
    rows: Optional[List[Dict[str, Any]]],
    k: int,
) -> List[Candidate]:
    """Adapt the first ``k`` rows of a ranked answer."""
    return [adapt(row) for row in (rows or [])[:k]]


def _chembl(client: Any) -> Callable[[Dict[str, Any]], Candidate]:
    """ChEMBL's adapter, which needs the client to fetch synonyms."""
    return lambda row: candidate_from_chembl_row(row, client)


def _zeropm_table(label: str, table: Any) -> List[Candidate]:
    """Adapt a ZeroPM lookup table into at most one candidate."""
    cand = candidate_from_zeropm_name_table(label, table)
    return [cand] if cand is not None else []


def rank_rows_by_completeness(rows: Optional[List[Dict[str, Any]]]) -> List[Dict[str, Any]]:
    """Sort source rows by number of non-null fields, most complete first.

    A formula matches many compounds and the sources return them in no useful
    order, so the formula lookups rank by how much each row says before
    taking the top ``k``.

    Args:
        rows: Raw source rows, or None.

    Returns:
        A new list ordered by descending completeness; stable for ties.

    Example:
        >>> rank_rows_by_completeness([{"a": 1, "b": None}, {"a": 1, "b": 2}])
        [{'a': 1, 'b': 2}, {'a': 1, 'b': None}]
    """
    return sorted(
        rows or [],
        key=lambda row: sum(1 for v in row.values() if not is_missing(v)),
        reverse=True,
    )


# ── Source-specific searches with no public method ───────────────────────────


def comptox_skeleton_search(comptox: Any, skeleton: str) -> List[Dict[str, Any]]:
    """Search CompTox for InChIKeys sharing a 14-character skeleton.

    Runs a ``LIKE 'skeleton%'`` query on CompTox's SQLite table, since the
    client has no public prefix search.

    Args:
        comptox: An open :class:`~provesid.CompToxID` client.
        skeleton: The 14-character InChIKey connectivity block.

    Returns:
        Up to 20 matching rows; empty on a miss or a failed query.
    """
    try:
        cur = comptox._conn.execute(
            "SELECT * FROM chemicals WHERE INCHIKEY LIKE ? LIMIT 20",
            (f"{skeleton}%",),
        )
        cols = [d[0] for d in cur.description]
        return [dict(zip(cols, row)) for row in cur.fetchall()]
    except Exception as exc:
        log.warning("CompTox skeleton search (SQL) failed: %s", exc)
        return []


def pubchem_skeleton_search(pubchem: Any, skeleton: str) -> List[Dict[str, Any]]:
    """Search PubChemID for InChIKeys sharing a 14-character skeleton.

    Args:
        pubchem: An open :class:`~provesid.PubChemID` client.
        skeleton: The 14-character InChIKey connectivity block.

    Returns:
        Up to 20 matching rows; empty on a miss or a failed query.
    """
    try:
        cur = pubchem._conn.execute(
            "SELECT * FROM compounds WHERE inchikey LIKE ? LIMIT 20",
            (f"{skeleton}%",),
        )
        cols = [d[0] for d in cur.description]
        return [dict(zip(cols, row)) for row in cur.fetchall()]
    except Exception as exc:
        log.warning("PubChemID skeleton search (SQL) failed: %s", exc)
        return []


def chebi_skeleton_search(chebi: Any, skeleton: str) -> List[Dict[str, Any]]:
    """Search ChEBI's in-memory InChIKey index for a 14-character skeleton.

    Args:
        chebi: A loaded :class:`~provesid.ChebiSDF` client.
        skeleton: The 14-character InChIKey connectivity block.

    Returns:
        Up to 20 matching compounds; empty on a miss or a failed scan.
    """
    try:
        results = []
        for inchikey, chebi_id in chebi.index.get("inchikey_to_id", {}).items():
            if inchikey.startswith(skeleton):
                compound = chebi.get_compound_by_id(chebi_id)
                if compound:
                    results.append(compound)
                if len(results) >= 20:
                    break
        return results
    except Exception as exc:
        log.warning("ChEBI skeleton search failed: %s", exc)
        return []


def _exact_or_skeleton(
    client: Any,
    inchikey: str,
    search: Callable[[Any, str], List[Dict[str, Any]]],
) -> Optional[Dict[str, Any]]:
    """The exact InChIKey row if there is one, else the first skeleton match."""
    row = client.get_by_inchikey(inchikey)
    if row is None:
        rows = search(client, inchikey[:14])
        row = rows[0] if rows else None
    return row


# ── ZeroPM's fuzzy retrieval ─────────────────────────────────────────────────


def _zeropm_fuzzy(client: Any, q: Query) -> List[Candidate]:
    """ZeroPM's fuzzy name retrieval, labelled with what it actually matched.

    ZeroPM is the only source that does true fuzzy *retrieval*; the others
    are substring-matched with ``exact=False``.  So it is the one that can
    reach a typo like "asprin" -> "aspirin".  The candidate is named after
    the matched name rather than the query, and carries ZeroPM's own
    similarity as ``query_match_score`` so that it is not re-derived from a
    name the candidate was just given.
    """
    table = client.get_id_table_from_similar_name(
        q.value, number_of_results=q.k, score_cutoff=q.fuzzy_cutoff
    )
    if table is None or table.empty:
        return []
    matched_name = str(table["matched_name"].iloc[0])
    candidates = _zeropm_table(matched_name, table)
    for cand in candidates:
        cand["query_match_score"] = float(table["match_score"].iloc[0]) / 100.0
    return candidates


# ── The online services ──────────────────────────────────────────────────────

#: The PUG-REST properties an online PubChem candidate is built from.
PUBCHEM_ONLINE_PROPERTIES: List[str] = [
    "Title", "IUPACName", "MolecularFormula", "SMILES", "InChI", "InChIKey",
    "MolecularWeight",
]


def _pubchem_online(cids_for: Callable[[Any, str], Any]) -> Lookup:
    """A PUG-REST lookup: identifier -> CIDs -> one candidate per CID.

    Costs two requests plus one per CID taken: the CID lookup, one property
    table for every CID at once, and each CID's synonyms, which are where
    PubChem keeps CAS numbers.  A CID PubChem answers with no properties is
    dropped, as :meth:`~provesid.PubChemID.properties_for_cids` drops it.

    Args:
        cids_for: ``(api, value) -> CIDs``, the one call that differs between
            identifier kinds.

    Returns:
        A lookup for the table.
    """

    def lookup(api: Any, q: Query) -> List[Candidate]:
        try:
            cids = cids_for(api, q.value)
        except NotFoundError:
            return []
        # PubChem answers some unknown identifiers with CID 0 instead of a 404.
        cids = [cid for cid in (cids if isinstance(cids, list) else []) if cid][: q.k]
        if not cids:
            return []
        rows = api.get_properties_for_cids(cids, PUBCHEM_ONLINE_PROPERTIES)
        return [
            candidate_from_pubchem_online(row, api.get_compound_synonyms(row["CID"]))
            for row in rows
            if len(row) > 1
        ]

    return lookup


def _cactus(resolver: Any, q: Query) -> List[Candidate]:
    """Ask the NCI/CADD resolver for a structure and its names.

    CACTUS takes any identifier it recognises as the same URL segment, so one
    function serves every kind it has a row for.  Two requests: ``smiles``,
    then ``names``.  When an identifier is ambiguous CACTUS lists one SMILES
    per line; the first is taken, as the first CID is from PubChem.  A name
    list that fails to come back costs the names, not the structure.
    """
    try:
        answer = resolver.resolve(q.value, "smiles")
    except NotFoundError:
        return []
    lines = [line.strip() for line in str(answer or "").splitlines() if line.strip()]
    if not lines:
        return []
    try:
        names = resolver.resolve(q.value, "names").splitlines()
    except NotFoundError:
        names = []
    return [candidate_from_cactus(lines[0], names)]


def _cas_or_name_cids(api: Any, value: str) -> Any:
    """CIDs whose synonyms include ``value`` exactly: PubChem's name index."""
    return api.get_cids_by_name(value, name_type="complete")


# ── The table ────────────────────────────────────────────────────────────────

#: ``LOOKUPS[kind][source](client, query)`` returns that source's candidates.
LOOKUPS: Dict[str, Dict[str, Lookup]] = {
    "cas": {
        "chebi": lambda c, q: _top(candidate_from_chebi_row, c.search_by_cas(q.value), q.k),
        # A retired or alternate number, when it is no chemical's own CASRN.
        "comptox": lambda c, q: _one(
            candidate_from_comptox_row,
            c.get_by_casrn(q.value) or c.get_by_alternate_casrn(q.value),
        ),
        "pubchem": lambda c, q: _one(candidate_from_pubchem_row, c.get_by_cas(q.value)),
        "zeropm": lambda c, q: _zeropm_table(q.label, c.get_id_table_from_cas(q.value)),
        "pubchem_online": _pubchem_online(_cas_or_name_cids),
        "cactus": _cactus,
    },
    "inchikey": {
        "chebi": lambda c, q: _one(candidate_from_chebi_row, c.search_by_inchikey(q.value)),
        "comptox": lambda c, q: _one(candidate_from_comptox_row, c.get_by_inchikey(q.value)),
        "pubchem": lambda c, q: _one(candidate_from_pubchem_row, c.get_by_inchikey(q.value)),
        "zeropm": lambda c, q: _zeropm_table(q.label, c.get_id_table_from_inchikey(q.value)),
        "chembl": lambda c, q: _one(_chembl(c), c.search_by_inchikey(q.value)),
        "pubchem_online": _pubchem_online(lambda api, v: api.get_cids_by_inchikey(v)),
        "cactus": _cactus,
    },
    # ``value`` is the full InChIKey; the exact key is preferred where the
    # source can be asked for it, and the 14-character skeleton otherwise.
    "inchikey_skeleton": {
        "chebi": lambda c, q: _top(
            candidate_from_chebi_row, chebi_skeleton_search(c, q.value[:14]), 1
        ),
        "comptox": lambda c, q: _one(
            candidate_from_comptox_row, _exact_or_skeleton(c, q.value, comptox_skeleton_search)
        ),
        "pubchem": lambda c, q: _one(
            candidate_from_pubchem_row, _exact_or_skeleton(c, q.value, pubchem_skeleton_search)
        ),
    },
    "inchi": {
        "chebi": lambda c, q: _one(candidate_from_chebi_row, c.search_by_inchi(q.value)),
        "pubchem": lambda c, q: _one(candidate_from_pubchem_row, c.get_by_inchi(q.value)),
        "zeropm": lambda c, q: _zeropm_table(q.label, c.get_id_table_from_inchi(q.value)),
        "pubchem_online": _pubchem_online(lambda api, v: api.get_cids_by_inchi(v)),
        "cactus": _cactus,
    },
    "smiles": {
        "comptox": lambda c, q: _one(candidate_from_comptox_row, c.get_by_smiles(q.value)),
        "pubchem": lambda c, q: _one(candidate_from_pubchem_row, c.get_by_smiles(q.value)),
        "zeropm": lambda c, q: [
            cand for cand in [candidate_from_zeropm_smiles(q.value, c)] if cand is not None
        ],
        "chembl": lambda c, q: _one(_chembl(c), c.search_by_smiles(q.value)),
        "pubchem_online": _pubchem_online(lambda api, v: api.get_cids_by_smiles(v)),
        "cactus": _cactus,
    },
    "dtxsid": {
        "comptox": lambda c, q: _one(candidate_from_comptox_row, c.get_by_dtxsid(q.value)),
        "pubchem_online": _pubchem_online(_cas_or_name_cids),
    },
    # Exact name, falling back to the source's synonym or preferred-name
    # index where it keeps them apart.
    "name": {
        "chebi": lambda c, q: _top(
            candidate_from_chebi_row,
            c.search_by_name(q.value, exact=True) or c.search_by_synonym(q.value, exact=True),
            q.k,
        ),
        # Every name CompTox holds, synonyms included, via its name index.
        "comptox": lambda c, q: _top(
            candidate_from_comptox_row, c.search_by_name(q.value, exact=True, limit=q.k), q.k
        ),
        "pubchem": lambda c, q: _top(
            candidate_from_pubchem_row, c.search_by_name(q.value, exact=True, limit=q.k), q.k
        ),
        "zeropm": lambda c, q: _zeropm_table(q.label, c.get_id_table_from_name(q.value)),
        "chembl": lambda c, q: _top(
            _chembl(c), c.search_by_name(q.value, limit=q.k, exact=True), q.k
        ),
        "pubchem_online": _pubchem_online(_cas_or_name_cids),
        "cactus": _cactus,
    },
    # Substring matching everywhere but ZeroPM, which retrieves by similarity.
    "fuzzy_name": {
        "chebi": lambda c, q: _top(
            candidate_from_chebi_row,
            c.search_by_name(q.value, exact=False) or c.search_by_synonym(q.value, exact=False),
            q.k,
        ),
        "comptox": lambda c, q: _top(
            candidate_from_comptox_row, c.search_by_name(q.value, exact=False, limit=q.k), q.k
        ),
        "pubchem": lambda c, q: _top(
            candidate_from_pubchem_row, c.search_by_name(q.value, exact=False, limit=q.k), q.k
        ),
        "zeropm": _zeropm_fuzzy,
        "chembl": lambda c, q: _top(
            _chembl(c), c.search_by_name(q.value, limit=q.k, exact=False), q.k
        ),
    },
    "formula": {
        "chebi": lambda c, q: _top(
            candidate_from_chebi_row, rank_rows_by_completeness(c.search_by_formula(q.value)), q.k
        ),
        "comptox": lambda c, q: _top(
            candidate_from_comptox_row, rank_rows_by_completeness(c.search_by_formula(q.value)), q.k
        ),
        "pubchem": lambda c, q: _top(
            candidate_from_pubchem_row, rank_rows_by_completeness(c.search_by_formula(q.value)), q.k
        ),
    },
}
