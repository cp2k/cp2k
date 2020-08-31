#!/bin/sh -l
# author: Tiziano Müller
# SPDX-License-Identifier: MIT

set -o errexit
set -o nounset

PYTHON=$(command -v python3 || true)

if [ x"${PYTHON}" = "x" ] ; then
    echo "ERROR: the python3 executable could not be found, but is required for the complete build process"
    exit 1
fi

SCRIPTDIR=$(CDPATH='' cd -- "$(dirname -- "$0")" && pwd)

# here we assume that this script is located in the tools/ subfolder
ENVDIR="${SCRIPTDIR}/../.pre-commit-env"

if [ -f "${ENVDIR}/bootstrap_done" ] ; then
    echo "The environment in <CP2K-SRC-DIR>/.pre-commit-env/ is ready to use..."
    exit 0
fi

echo "Installing all required tools to run pre-commit in <CP2K-SRC-DIR>/.pre-commit-env/"

echo "Creating the virtualenv for pre-commit..."
if [ ! -d "${ENVDIR}" ] ; then
    "${PYTHON}" -m venv "${ENVDIR}"
fi

# resolve the relative path now that we have it:
ENVDIR=$(CDPATH='' cd -- "${ENVDIR}" && pwd)

# shellcheck source=/dev/null
. "${ENVDIR}/bin/activate"

# change the command to the prefixed one for the rest of the script:
PYTHON=$(command -v python)
PIP=$(command -v pip || true)

if [ x"${PIP}" = "x" ] ; then
    echo "Bootstrapping the pip command inside the virtual environment..."

    FETCHER="$(command -v curl || true)"
    FETCHER_OPTS=""

    if [ x"${FETCHER}" = "x" ] ; then
        FETCHER="$(command -v wget || true)"
        FETCHER_OPTS="-O-"
    fi

    if [ x"${FETCHER}" = "x" ] ; then
        echo "ERROR: neither wget nor curl seem to be available, please download"
        echo "    https://bootstrap.pypa.io/get-pip.py"
        echo "manually and run:"
        echo "    python3 get-pip.py"
        echo "to install the pip command, then rerun this script."
        rm -rf "${ENVDIR}"  # cleaning up to avoid picking up the unprefixed pip on the rerun
        exit 1
    fi

    "${FETCHER}" ${FETCHER_OPTS} "https://bootstrap.pypa.io/get-pip.py" | "${PYTHON}"
fi

PIP=$(command -v pip)

echo "Installing pre-commit..."
"${PIP}" install pre-commit
# this will also install things like shellcheck and a pre-built node if necessary

PRE_COMMIT=$(command -v pre-commit)

echo "Letting pre-commit setup the rest..."
"${PRE_COMMIT}" install --install-hooks

echo "successfully bootstrapped at $(date) by $0" > "${ENVDIR}/bootstrap_done"

echo ""
echo "All done: pre-commit will now run on each 'git commit' on the files you are committing!"

echo ""
echo "Should you need to run it manually, use"
echo ""
echo "    ${ENVDIR}/bin/pre-commit"
echo ""
echo "or add the directory ${ENVDIR}/bin to your \$PATH variable."
echo "Happy hacking! :)"
