"""Presets — name a whole matching policy in one word.

``Search`` has a dozen arguments that decide what counts as a match and what
is returned.  ``Search.PRESETS`` bundles them into three named policies:

- ``"balanced"`` (the default): exact matching, one row per query.
- ``"strict"``: as balanced, but only structures two independent databases
  agree on (``min_source_support=2``).
- ``"recall"``: fuzzy names, InChIKey-skeleton and Tanimoto widening, ZeroPM
  queried, and every plausible compound returned (``n_hits="all"``).

Any argument passed explicitly overrides the preset, and every result frame
records the preset and the settings it ran with in ``df.attrs``, so a saved
table says how it was made.

Run with::

    uv run python examples/search/presets_demo.py
"""

import pandas as pd

from provesid import Search

# The presets are plain dicts of constructor arguments: inspect them side by side.
print(pd.DataFrame(Search.PRESETS).to_string())

names = ["atrazine", "paracetamol", "atrazin"]  # the last one is a typo

for preset in ["balanced", "strict", "recall"]:
    # n_hits=1 overrides recall's "all", so the three runs line up row for row.
    with Search("name", preset=preset, n_hits=1, show_progress=False) as s:
        df = s.search(names)
    print(f"\n=== preset={preset!r} ===")
    print(
        df[["query", "name", "InChIKey", "n_source_support", "confidence"]]
        .to_string(index=False)
    )

# The frame remembers how it was produced.
print("\npreset:  ", df.attrs["preset"])
print("settings:", df.attrs["settings"])
