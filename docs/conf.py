# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "CP2K"
copyright = "2000-2025, CP2K Developers"
author = "CP2K Developers"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx_rtd_theme",
    "sphinx_reredirects",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinxcontrib.youtube",
    "sphinxcontrib.mermaid",
]

myst_enable_extensions = [
    "attrs_inline",
    "dollarmath",
    "smartquotes",
    "strikethrough",
]

myst_fence_as_directive = ["mermaid"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "README.md"]

suppress_warnings = ["ref"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_favicon = "_static/favicon.png"
html_copy_source = False

add_module_names = False

# https://myst-parser.readthedocs.io/en/v0.16.1/syntax/optional.html#syntax-header-anchors

html_context = {
    "display_github": True,
    "github_user": "cp2k",
    "github_repo": "cp2k",
    "github_version": "master",
    "conf_py_path": "/docs/",
}

# https://sphinx-rtd-theme.readthedocs.io/en/stable/configuring.html

html_theme_options = {
    "collapse_navigation": False,
}

redirects = {
    "getting-started/CMake": "build-from-source.html",
    "getting-started/spack": "build-with-spack.html",
}

# EOF
