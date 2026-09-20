"""
Decide what PROVESID downloads, instead of finding out afterwards.

``Search("cas").search("50-00-0")`` on a clean machine used to fetch about
32 GB -- ChEBI, CompTox, PubChem and ChEMBL, each of whose clients defaults to
``auto_download=True`` -- without announcing the total, asking, or offering to
run on the sources already present.  On a laptop that is often the whole free
disk.

So the datasets are now something you install deliberately:

* ``datasets.status()`` -- what is on disk, and what it costs.
* ``datasets.plan(...)`` -- what a download would transfer, before it starts.
* ``datasets.fetch(...)`` -- install by name, resumable and verified.
* ``datasets.remove(...)`` -- reclaim the space, by name.

and ``Search`` no longer downloads behind your back: ``datasets="present"`` is
the default.

This script downloads nothing.  The sections that would (``fetch``) print the
call instead, and the ``Search`` sections point at an empty temporary directory
to show what a clean machine now does.

Run with::

    python examples/datasets/dataset_manager_demo.py
"""

import logging
import tempfile

import pandas as pd

from provesid import Search, datasets
from provesid.datasets import MissingDatasetError

logging.basicConfig(level=logging.WARNING, format="%(levelname)-7s %(message)s")
pd.set_option("display.width", 200)
pd.set_option("display.max_colwidth", 60)


def rule(title):
    print(f"\n{'─' * 78}\n{title}\n")


def main():
    # ── 1. What is on disk ───────────────────────────────────────────────────
    rule("1. datasets.status() -- built from filenames and stat calls, so it is "
         "instant\n   even with 30 GB of ChEMBL in the directory.")

    status = datasets.status()
    print(status[["dataset", "present", "files", "size", "release"]].to_string(index=False))
    print(f"\n   data directory: {status.attrs['data_dir']}")
    print(f"   occupied:       {datasets.human_bytes(status.attrs['total_bytes'])}")

    # A leftover `.part` from an interrupted download is counted here but does
    # not make a dataset "present" -- which is what explains a full disk after
    # a download that was cancelled.

    # ── 2. What a clean machine would pay ────────────────────────────────────
    rule("2. datasets.plan() -- the number nobody used to be told.")

    with tempfile.TemporaryDirectory() as empty:
        todo = datasets.plan(datasets.DEFAULT_DATASETS, data_dir=empty)
        print(todo[["dataset", "action", "download", "installed", "role"]].to_string(index=False))
        print(f"\n   transfer:      {datasets.human_bytes(todo.attrs['total_download_bytes'])}")
        print(f"   on disk after: {datasets.human_bytes(todo.attrs['total_resident_bytes'])}")
        print(f"   needed at peak:{datasets.human_bytes(todo.attrs['peak_bytes'])}"
              "   <- ChEMBL's archive sits beside the database it extracts into")

        print("\n   ChEMBL is 87% of that, and in a search it only *enriches* a "
              "structure\n   the other sources already found.  Three sources are a "
              "reasonable install:")
        three = datasets.plan(["pubchem", "comptox", "chebi"], data_dir=empty)
        print(f"   {datasets.human_bytes(three.attrs['total_download_bytes'])} to "
              f"download, {datasets.human_bytes(three.attrs['total_resident_bytes'])} on disk.")

        # ── 3. Search on a machine with nothing installed ────────────────────
        rule("3. Search(datasets='present') -- the default.  Nothing is downloaded;\n"
             "   the search runs on what is there and says what is missing.")

        search = Search("cas", data_dir=empty, show_progress=False)
        frame = search.search("50-00-0")
        print(f"   sources available:   {search.sources_available}")
        print(f"   sources unavailable: {search.sources_unavailable}")
        print(f"   result:              name={frame.iloc[0]['name']!r}, "
              f"confidence={frame.iloc[0]['confidence']}")
        print("\n   (the WARNING lines above name each missing dataset, its size, "
              "and the\n   fetch call that installs it)")

        # ── 4. When a partial answer is worse than none ──────────────────────
        rule("4. Search(datasets='required') -- refuse rather than answer from "
             "fewer sources.")

        try:
            Search("cas", datasets="required", data_dir=empty)
        except MissingDatasetError as exc:
            print(exc)

        print("\n   It raises in the constructor, before the first query, and "
              "downloads nothing.")

    # ── 5. Installing and reclaiming ─────────────────────────────────────────
    rule("5. fetch and remove -- printed rather than run, since these move "
         "gigabytes.")

    print("""    from provesid import datasets

    datasets.fetch("pubchem")                  # one dataset, by name
    datasets.fetch(["pubchem", "comptox", "chebi"])   # skips what is present
    datasets.fetch("chembl", force=True)       # re-download; picks up the
                                               # current release

    datasets.remove("chembl", dry_run=True)    # list what would go
    datasets.remove("chembl")                  # 27.7 GiB back, and every
                                               # release and extract with it

    Search("cas", datasets="auto")             # the old behaviour, if you
                                               # want it""")

    rule("6. ChEMBL is worth a second look.")
    print("""   27.7 GiB for a source that only adds ChEMBL IDs to structures the
   others already found.  The package reads 8 of ChEMBL's 74 tables, and an
   extract of just those answers identically:

       from provesid import CheMBL
       CheMBL.compact(remove_source=True)      # 27.7 GiB -> 2.4 GiB

   which turns a full install from ~32 GB into ~6.5 GB.""")


if __name__ == "__main__":
    main()
