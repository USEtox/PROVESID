"""
MkDocs hooks for the PROVESID documentation. Wired up in ``mkdocs.yml`` under
``hooks:``.

**The tutorials are served from ``examples/``.** They live in
``examples/<area>/`` next to the scripts they go with, and that is the only
copy. Each is added to the site at the path it would have had inside
``docs/`` --- ``examples/pubchem/pubchem_tutorial.ipynb`` --- so the nav and the
links in ``docs/`` name them as if they were there. A tutorial is any ``.md``
or ``.ipynb`` file under ``examples/`` except a folder's ``README.md``, which
describes the scripts beside it for someone reading the repository. The
``.py`` demos are not pages: mkdocs-jupyter would render them as notebooks,
which they are not.

**One upstream warning is dropped.** ``mkdocs.yml`` asks mkdocstrings for
``docstring_style: auto`` with per-style options. mkdocstrings-python 2.0
fills in defaults for the Sphinx style too, then warns that griffe's Sphinx
parser does not accept one of them (``warn_missing_types``) --- once per
module, which ``--strict`` turns into a failed build. The package has no
Sphinx-style docstrings, so that one message is filtered and nothing else.
"""

import logging
from pathlib import Path

from mkdocs.plugins import event_priority
from mkdocs.structure.files import File

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples"
TUTORIAL_SUFFIXES = {".md", ".ipynb"}


class _DropSphinxOptionWarning(logging.Filter):
    def filter(self, record):
        return "unsupported sphinx docstring parser option" not in record.getMessage()


logging.getLogger(
    "mkdocs.plugins.mkdocstrings_handlers.python._internal.handler"
).addFilter(_DropSphinxOptionWarning())


def _tutorials():
    for path in sorted(EXAMPLES_DIR.rglob("*")):
        if (path.suffix in TUTORIAL_SUFFIXES
                and path.name != "README.md"
                and ".ipynb_checkpoints" not in path.parts):
            yield path


# Ahead of mkdocs-jupyter's own on_files (priority 0), which is what turns a
# jupytext .md or an .ipynb into a notebook page: it has to see these files.
@event_priority(50)
def on_files(files, config):
    """Add every tutorial under ``examples/`` to the site's files."""
    project_root = EXAMPLES_DIR.parent
    for path in _tutorials():
        files.append(File(
            path.relative_to(project_root).as_posix(),
            src_dir=str(project_root),
            dest_dir=config["site_dir"],
            use_directory_urls=config["use_directory_urls"],
        ))
    return files


def on_serve(server, config, builder):
    """Rebuild on ``mkdocs serve`` when a tutorial changes."""
    server.watch(str(EXAMPLES_DIR))
    return server
