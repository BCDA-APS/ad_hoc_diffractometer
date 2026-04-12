# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import pathlib
import re
import sys
import tomllib
from importlib.metadata import version as _version

# -- Path setup --------------------------------------------------------------

root_path = pathlib.Path(__file__).parent.parent.parent
sys.path.insert(0, str(root_path / "src"))

with open(root_path / "pyproject.toml", "rb") as _f:
    _toml = tomllib.load(_f)
_metadata = _toml["project"]

# -- Project information -----------------------------------------------------

project = "ad_hoc_diffractometer"
author = _metadata["authors"][0]["name"]
copyright = f"2026, {author}"
description = _metadata["description"]

release = _version("ad_hoc_diffractometer")
version = ".".join(release.split(".")[:2])

# version_match: used by the pydata-sphinx-theme version switcher.
# Set DOC_VERSION in CI to "latest" (main branch) or the tag (e.g. "0.3.0").
# Falls back to release so local builds still work.
switcher_version_match = os.environ.get("DOC_VERSION", release)

# -- General configuration ---------------------------------------------------

extensions = [
    "autoapi.extension",
    "myst_nb",
    "sphinx.ext.graphviz",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_tabs.tabs",
]

exclude_patterns = [
    "**.ipynb_checkpoints",
    "_build",
    "Thumbs.db",
    ".DS_Store",
    # roadmap.md lives in the repo root (not docs/source); guard against
    # accidental inclusion if the file is ever copied back here.
    "roadmap.md",
]

myst_enable_extensions = ["colon_fence"]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "restructuredtext",
    ".ipynb": "myst-nb",
}

templates_path = ["_templates"]

# -- Notebook execution (myst-nb) --------------------------------------------

# Execute notebooks during the doc build.  All notebooks in docs/source/
# use only ad_hoc_diffractometer + NumPy + matplotlib — no hardware required.
nb_execution_mode = "auto"  # execute only notebooks without stored outputs
nb_execution_timeout = 120  # seconds per cell

# -- AutoAPI -----------------------------------------------------------------

# Keep generated files inside the build tree, not in the source tree.
# Without this, sphinx-autoapi writes docs/source/autoapi/ which must not
# be committed to git.
autoapi_root = "autoapi"
autoapi_dirs = [str(root_path / "src")]
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    # "imported-members" is intentionally omitted: including it causes
    # "more than one target" cross-reference warnings for every symbol
    # re-exported via __init__.py.
]
autoapi_member_order = "alphabetical"
autoapi_python_class_content = "both"
autoapi_template_dir = "_templates/autoapi"
autoapi_add_toctree_entry = False  # toctree entry added manually in api.rst
suppress_warnings = ["autoapi.python_import_resolution"]


def autoapi_skip_member(app, what, name, obj, skip, options):
    """Skip logger instances and private _version from the autoapi output."""
    if what == "data" and name.endswith(".logger"):
        return True
    if "_version" in name:
        return True
    return skip


def _shorten_type_str(s: str) -> str:
    """Replace dotted.module.ClassName with ClassName in type annotation strings.

    Used in literal code spans where Sphinx's domain resolver does not run.
    Bare names (``str``, ``int``, ``None``) are left unchanged.
    """
    return re.sub(
        r"(?<![.\w])([a-zA-Z_]\w*(?:\.[a-zA-Z_]\w*)+)",
        lambda m: m.group(0).split(".")[-1],
        s,
    )


def _tilde_type_str(s: str) -> str:
    """Prepend ~ to dotted.module.ClassName in type annotation strings.

    Used inside ``.. py:*::`` directives so Sphinx resolves the full cross-
    reference but displays only the short name (the ``~`` prefix convention).
    Bare names (``str``, ``int``, ``None``) are left unchanged.
    """
    return re.sub(
        r"(?<![.\w~])([a-zA-Z_]\w*(?:\.[a-zA-Z_]\w*)+)",
        lambda m: "~" + m.group(0),
        s,
    )


def _prepare_jinja_env(jinja_env) -> None:
    """Register custom Jinja2 filters for autoapi templates."""
    jinja_env.filters["shorten_type"] = _shorten_type_str
    jinja_env.filters["tilde_type"] = _tilde_type_str


def setup(app):
    """Connect autoapi events and register Jinja2 filters."""
    app.connect("autoapi-skip-member", autoapi_skip_member)
    app.config.autoapi_prepare_jinja_env = _prepare_jinja_env


# -- Intersphinx -------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
}

# -- Napoleon (NumPy/Google docstring style) ---------------------------------

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True

# -- HTML output -------------------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

html_theme_options = {
    "github_url": "https://github.com/prjemian/ad_hoc_diffractometer",
    "navbar_end": ["version-switcher", "theme-switcher", "navbar-icon-links"],
    "switcher": {
        "json_url": (
            "https://prjemian.github.io/ad_hoc_diffractometer"
            "/latest/_static/switcher.json"
        ),
        "version_match": switcher_version_match,
    },
    "show_version_warning_banner": True,
}

html_title = f"{project} {version}"

# -- Copy-button (sphinx-copybutton) -----------------------------------------

copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
