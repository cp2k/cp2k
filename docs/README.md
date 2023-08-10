# CP2K Documentation

These are the source of the [CP2K manual](https://manual.cp2k.org/trunk). They are published daily by [this script](../tools/docker/scripts/test_manual.sh).

To build  a local version of the manual perform the following steps:

1. Install the required Python packaged:

   `pip3 install sphinx myst-parser sphinx_rtd_theme lxml`

1. Build a CP2K binary and generate `cp2k_input.xml` file:

   `cp2k/exe/local/cp2k.psmp --xml`

1. Generate the Markdown files from the XML file:

   `cp2k/docs/generate_input_reference.py ./cp2k_input.xml ./references.html`

1. And finally run Sphinx:

   `sphinx-build cp2k/docs/ output_dir -W -n --keep-going --jobs 16`
