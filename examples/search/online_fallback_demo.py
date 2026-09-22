"""Online fallback — ask PubChem and CACTUS what the offline databases miss.

``Search`` is offline by default: it opens no socket, and a query none of the
installed databases holds comes back as an empty row.  With
``online_fallback=True`` such a query — and only such a query — is retried
against PubChem's PUG-REST service and the NCI/CADD Chemical Identifier
Resolver (CACTUS).  Rows the network supplied say so in ``source`` and
``source_details``, and ``df.attrs`` counts how many queries went online.

Run with::

    uv run python examples/search/online_fallback_demo.py
"""

import logging

from provesid import Search

# Each fallback is logged at DEBUG; switch it on to watch the network being used.
logging.basicConfig(level=logging.WARNING)
logging.getLogger("provesid.search").setLevel(logging.DEBUG)

# Old datasets often carry CAS numbers that CAS has since deleted or merged.
# CompTox lists most of them beside the current number, so they are answered
# offline too; the network is asked only about what no database holds.
cas_list = [
    "50-78-2",     # Aspirin: in every offline database, so never asked online
    "39400-72-1",  # A retired CAS number of atrazine: CompTox lists it, so offline
    "0000-00-0",   # Not a compound: every source misses, online included
]

with Search("cas", online_fallback=True, show_progress=False) as s:
    df = s.search(cas_list)

print(
    df[["query", "name", "InChIKey", "source", "n_source_support", "confidence"]]
    .to_string(index=False)
)

# How much network did the batch use?
print("\nOffline sources:  ", df.attrs["sources_available"])
print("Queries sent online:", df.attrs["online_fallbacks"])
print("...answered online: ", df.attrs["online_resolved"])

# Which source said what, for one row.
for source, detail in df.iloc[1]["source_details"].items():
    print(f"  {source:18s} found={detail['found']!s:5s} fields={detail['fields']}")
