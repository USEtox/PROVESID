"""
Build pubchem_id.db from PubChem's FTP site, from the command line.

A thin wrapper over :func:`provesid.pubchem_ftp.build_pubchem_id_db`, for the
one job that is easier from a shell than from Python: refreshing the copy on
Zenodo that ``PubChemID(source="zenodo")`` downloads. Build it, check
``PubChemID(db_path=...).provenance()``, upload the file.

Usage:
    python scripts/build_pubchem_id_db.py                      # newest snapshot
    python scripts/build_pubchem_id_db.py --release 2026-09-01 --out ./pubchem_id.db
    python scripts/build_pubchem_id_db.py --list               # available releases
"""

import argparse
import logging
import sys

from provesid.pubchem_ftp import LATEST, build_pubchem_id_db, list_releases


def main(argv=None) -> int:
    """
    Parse the arguments and run the build.

    Args:
        argv: Command-line arguments, without the program name. Defaults to
            ``sys.argv[1:]``.

    Returns:
        The process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", default=None,
                        help="database path (default: the per-user dataset directory)")
    parser.add_argument("--release", default=LATEST,
                        help="'latest' (default), 'current', or a snapshot date YYYY-MM-DD")
    parser.add_argument("--list", action="store_true",
                        help="list the available releases and exit")
    parser.add_argument("--no-inchi", action="store_true",
                        help="skip CID-InChI-Key.gz (7.4 GB) and compute InChI with RDKit")
    parser.add_argument("--no-synonyms", action="store_true",
                        help="leave out the synonym table")
    parser.add_argument("--keep-downloads", action="store_true",
                        help="keep the source files, so a rebuild does not fetch them again")
    parser.add_argument("--force", action="store_true",
                        help="replace an existing database")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    if args.list:
        for release in list_releases():
            print(release)
        return 0

    build_pubchem_id_db(
        args.out,
        release=args.release,
        include_inchi=not args.no_inchi,
        include_synonyms=not args.no_synonyms,
        keep_downloads=args.keep_downloads,
        force=args.force,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
