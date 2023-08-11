# CP2K Documentation

These are the source of the [CP2K manual](https://manual.cp2k.org/trunk). They are published daily by [this script](../tools/docker/scripts/test_manual.sh).

To build  a local version of the manual perform the following steps:

1. Install the required Python packaged:

   `pip3 install sphinx myst-parser sphinx_rtd_theme lxml`

1. Build a CP2K binary and use it to generate the `cp2k_input.xml` and `references.html` files:

   `../exe/local/cp2k.psmp --xml`

1. Generate Markdown pages from the aforementioned files:

   `./generate_input_reference.py ./cp2k_input.xml ./references.html`

1. Run Sphinx:

   `make html`

1. Browse the HTML output in the `_build/html` directory.
