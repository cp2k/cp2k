# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?= -W -n --keep-going --jobs 16
SPHINXBUILD   ?= sphinx-build
PAGEFIND      ?= python3 -m pagefind
PAGEFINDOPTS  ?=
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	./fix_github_links.py "$(BUILDDIR)"

.PHONY: help html pagefind Makefile

# Index only after Sphinx exits, so their memory footprints do not overlap.
html: Makefile
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@$(MAKE) --no-print-directory pagefind

# Also usable after a direct sphinx-build invocation (run from docs/).
pagefind:
	@$(PAGEFIND) --site "$(BUILDDIR)/html" $(PAGEFINDOPTS) --serve

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
