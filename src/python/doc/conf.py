import sys
import os

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
# Path to Gudhi.so from source path
sys.path.insert(0, os.path.abspath("."))

extensions = [
    "matplotlib.sphinxext.plot_directive",
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinxcontrib.bibtex",
    "sphinx_paramlinks",
]

autodoc_class_signature = "separated"

bibtex_bibfiles = ["../../biblio/bibliography.bib"]

templates_path = ["_templates"]

master_doc = "index"

import gudhi

# General information about the project.
project = gudhi.__name__
copyright = f"{gudhi.__copyright__} - MIT"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
version = gudhi.__version__

language = "en"

exclude_patterns = ["_build", "*.inc"]

pygments_style = "sphinx"

html_theme = "pydata_sphinx_theme"

html_context = {"default_mode": "light"}

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    "show_nav_level": 2,
    "navbar_start": ["navbar_start_gudhi.html"],
    "navbar_center": ["navbar_center_gudhi.html"],
    "navbar_end": ["navbar_end_gudhi.html"],
    "navbar_persistent": [],
    "secondary_sidebar_items": ["page-toc"],
    "switcher": {
        "json_url": "https://gudhi.inria.fr/switcher.json",
        "version_match": version,
    },
}

html_sidebars = {"**": ["version-switcher", "sidebar-gudhi.html", "search-gudhi-field.html"]}

html_title = f"{project} v{version} documentation"

html_static_path = ["_static"]

html_css_files = [
    "python_gudhi.css",
]
