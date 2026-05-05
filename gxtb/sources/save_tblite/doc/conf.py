# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.

import os
import subprocess
import sys

_dir = os.path.abspath(os.path.dirname(__file__))

import tblite
from tblite.library import get_version

project = "tblite"
author = "Sebastian Ehlert"
copyright = f"2021-2023, {author}"

version = "{}.{}.{}".format(*get_version())
release = version


extensions = [
    "breathe",
    "chemiscope.sphinx",
    "myst_nb",
    "sphinx_design",
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.extlinks",
    "sphinxcontrib.bibtex",
]

breathe_default_project = project
breathe_use_project_refids = True
breathe_projects = {
    project: "_doxygen/xml",
}
breathe_domain_by_extension = {
    "h": "c",
}
breathe_show_include = True

bibtex_bibfiles = ["_static/references.bib"]

jupyter_execute_notebooks = "force"
myst_enable_extensions = [
    "amsmath",
    "dollarmath",
    "colon_fence",
]

autosummary_generate = True

napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True

c_id_attributes = [
    "TBLITE_API_ENTRY",
    "TBLITE_API_CALL",
]

templates_path = ["_templates"]
source_suffix = [".rst"]
master_doc = "index"
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
pygments_style = "default"


html_theme = "sphinx_book_theme"
html_title = project

_extra_navbar = """
<div class="sd-fs-4">
<a href="https://github.com/orgs/tblite/discussions" target="_blank">
    <i class="fa fa-comments"></i>
</a>
<a href="https://github.com/tblite" target="_blank">
    <i class="fab fa-github"></i>
</a>
<a href="https://twitter.com/grimmelab" target="_blank">
    <i class="fab fa-twitter"></i>
</a>
</div>
"""

html_theme_options = {
    "repository_url": "https://github.com/tblite/tblite",
    "repository_branch": "main",
    "use_repository_button": True,
    "use_edit_page_button": True,
    "use_download_button": False,
    "path_to_docs": "doc",
    "extra_navbar": _extra_navbar,
}

html_static_path = ["_static"]


def run_doxygen(folder):
    """Run the doxygen make command in the designated folder"""

    try:
        retcode = subprocess.call("set -ex; cd %s; doxygen Doxyfile" % folder, shell=True)
        if retcode < 0:
            sys.stderr.write("doxygen terminated by signal %s" % (-retcode))
    except OSError as e:
        sys.stderr.write("doxygen execution failed: %s" % e)


def generate_doxygen_xml(app):
    """Run the doxygen make commands if we're on the ReadTheDocs server"""

    read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

    if read_the_docs_build:

        run_doxygen(_dir)


def setup(app):

    # Add hook for building doxygen xml when needed
    app.connect("builder-inited", generate_doxygen_xml)
