# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


# Disable benign warning: https://github.com/sphinx-doc/sphinx/blob/d3c91f951255c6729a53e38c895ddc0af036b5b9/sphinx/domains/python.py#L1283
import logging

logging.getLogger("sphinx.sphinx.domains.python").setLevel(logging.ERROR)

# add_module_names = False

project = "CP2K"
copyright = "2023, CP2K Developers"
author = "CP2K Developers"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["myst_parser", "sphinx_rtd_theme"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

suppress_warnings = ["app", "ref", "index"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_favicon = "_static/favicon.png"

# add_module_names = False

# https://myst-parser.readthedocs.io/en/v0.16.1/syntax/optional.html#syntax-header-anchors

# EOF
