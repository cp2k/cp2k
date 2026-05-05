import os
import sys

sys.path.insert(0, os.path.join(os.path.abspath(".."), "python"))
import dftd3


project = "s-dftd3"
author = "Sebastian Ehlert"
copyright = f"2019-2024, {author}"

version = dftd3.__version__
release = version

extensions = [
    "sphinx_design",
    "sphinx_copybutton",
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.extlinks",
    "sphinxcontrib.bibtex",
]

html_theme = "sphinx_book_theme"
html_title = "Simple DFT-D3"

html_theme_options = {
    "repository_url": "https://github.com/dftd3/simple-dftd3",
    "repository_branch": "main",
    "use_repository_button": True,
    "use_edit_page_button": True,
    "use_download_button": False,
    "path_to_docs": "doc",
}

html_css_files = []
html_static_path = ["_static"]
templates_path = ["_templates"]
locale_dirs = ["locales"]
autodoc_mock_imports = ["dftd3.library", "numpy", "ase", "qcelemental", "pyscf"]
bibtex_bibfiles = ["_static/references.bib"]

master_doc = "index"
