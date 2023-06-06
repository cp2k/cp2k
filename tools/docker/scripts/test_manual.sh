#!/bin/bash -e

# author: Ole Schuett

if (($# != 1)); then
  echo "Usage: test_manual.sh <ADD_EDIT_LINKS>"
  exit 1
fi

ADD_EDIT_LINKS=$1

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Compiling CP2K =========="
cd /opt/cp2k
echo -n "Compiling cp2k... "
if make -j VERSION=psmp &> make.out; then
  echo "done."
else
  echo -e "failed.\n\n"
  tail -n 100 make.out
  echo -e "\nSummary: Compilation failed."
  echo -e "Status: FAILED\n"
  exit 0
fi

echo -e "\n========== Installing Dependencies =========="
apt-get update -qq
apt-get install -qq --no-install-recommends \
  default-jre-headless \
  libsaxonhe-java \
  python3 \
  python3-pip
rm -rf /var/lib/apt/lists/*

pip3 install --quiet sphinx myst-parser sphinx_rtd_theme lxml

echo -e "\n========== Generating Manual =========="

mkdir -p /workspace/artifacts/manual
cd /workspace/artifacts/manual

/opt/cp2k/exe/local/cp2k.psmp --version
/opt/cp2k/exe/local/cp2k.psmp --xml

TOOLS=/opt/cp2k/tools
cp ${TOOLS}/manual/favicon.png .
cp ${TOOLS}/manual/toggle_folding.js .

set +e # disable error trapping for remainder of script
(
  set -e                            # abort if error is encountered
  sed -i 's/\x0/?/g' cp2k_input.xml # replace null bytes which would crash saxon
  SAXON="java -jar /usr/share/java/Saxon-HE.jar"
  $SAXON -o:index.html ./cp2k_input.xml ${TOOLS}/manual/cp2k_input.xsl add_edit_links="${ADD_EDIT_LINKS}"
  $SAXON -o:cp2k.vim ./cp2k_input.xml ${TOOLS}/input_editing/vim/vim.xsl
)
EXIT_CODE=$?
if ((EXIT_CODE)); then
  echo -e "\nSummary: Saxon run failed."
  echo -e "Status: FAILED\n"
  exit
fi

echo -e "\n========== Generating New Manual =========="
(
  set -e # abort if error is encountered
  /opt/cp2k/docs/generate_input_reference.py ./cp2k_input.xml ./references.html
  sphinx-build /opt/cp2k/docs/ /workspace/artifacts/manual/new --jobs 16
)
EXIT_CODE=$?
if ((EXIT_CODE)); then
  echo -e "\nSummary: Sphinx build failed"
  echo -e "Status: FAILED\n"
  exit
fi

echo -e "\nSummary: Manual generation works fine."
echo -e "Status: OK\n"

#EOF
