# CP2K Documentation

These are the source of the [CP2K manual](https://manual.cp2k.org/trunk). They are published daily
by [this script](../tools/docker/scripts/test_manual.sh).

To build a local version of the manual perform the following steps:

1. Install the required Python packaged:

   `pip3 install -r ./docs/requirements.txt`

1. (optional) Build a CP2K binary and use it to generate the `cp2k_input.xml` and `references.html`
   files:

   `../exe/local/cp2k.psmp --xml`

1. (optional) Generate Markdown pages from the aforementioned files:

   `./generate_input_reference.py ./cp2k_input.xml ./references.html`

1. Run Sphinx:

   `make html`

1. Browse the HTML output in the `_build/html` directory.

> \[!TIP\]
>
> While the first invocation of Sphinx can be quite slow, subsequent builds are significantly faster
> thanks to its doctree cache. Nevertheless, a build with the full input reference can take several
> minutes and requires a lot of memory. So for development it's advisable to build without the input
> reference. To check cross-references one can generate the input reference and then remove all
> pages except the relevant ones.
