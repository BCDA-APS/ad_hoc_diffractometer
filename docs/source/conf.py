# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import pathlib
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
    "myst_parser",
    "nbsphinx",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinx_design",
]

exclude_patterns = [
    "**.ipynb_checkpoints",
    "_build",
    "Thumbs.db",
    ".DS_Store",
    # Notebooks require pandoc to convert; excluded until CI installs pandoc
    # and notebook content is reviewed for Phase 2.
    "*.ipynb",
]

myst_enable_extensions = ["colon_fence"]

source_suffix = [".rst", ".md"]

templates_path = ["_templates"]

# Do not execute notebooks during doc build — notebooks may require
# hardware or data not available in CI.
nb_execution_mode = "off"

# -- AutoAPI -----------------------------------------------------------------

autoapi_dirs = [str(root_path / "src")]
autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "imported-members",
]
autoapi_member_order = "alphabetical"
autoapi_python_class_content = "both"
autoapi_template_dir = "_templates/autoapi"
suppress_warnings = ["autoapi.python_import_resolution"]

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
