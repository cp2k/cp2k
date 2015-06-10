#!/bin/bash

# output basename
BASENAME=$(basename ${PWD})

# temporary file
TEMPLATE=$(mktemp --tmpdir=. --suffix=.tex)

# dump pandoc template for latex
pandoc -D latex > ${TEMPLATE}

# adjust the template
sed -i \
  -e 's/\(\\documentclass\[.\+\]{.\+}\)/\1\n\\pagenumbering{gobble}/' \
  -e 's/\\usepackage{listings}/\\usepackage{listings}\\lstset{basicstyle=\\footnotesize\\ttfamily}/' \
  ${TEMPLATE}

# cleanup markup and pipe into pandoc using the template
sed \
  -e 's/\[\[.\+\](.\+)\]//' \
  -e '/!\[.\+\](.\+)/{n;d}' \
  README.md | tee >( \
pandoc \
  --latex-engine=xelatex \
  --template=${TEMPLATE} --listings \
  -f markdown_github+implicit_figures \
  -V documentclass=scrartcl \
  -V classoption=DIV=45 \
  -V linkcolor=black \
  -V citecolor=black \
  -V urlcolor=black \
  -o documentation/${BASENAME}.pdf) | \
pandoc \
  -f markdown_github+implicit_figures \
  -o documentation/${BASENAME}.docx

# remove temporary file
rm ${TEMPLATE}
