#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Compiling CP2K =========="
cd /workspace/cp2k
touch src/cp2k_info.F # ensure latest REVISION is picked up.
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

echo -e "\n========== Generating Manual =========="

mkdir -p /workspace/artifacts/manual
cd /workspace/artifacts/manual

/workspace/cp2k/exe/local/cp2k.psmp --version
/workspace/cp2k/exe/local/cp2k.psmp --xml

TOOLS=/workspace/cp2k/tools
cp ${TOOLS}/manual/favicon.png .
cp ${TOOLS}/manual/toggle_folding.js .

set +e # disable error trapping for remainder of script
(
  set -e                            # abort if error is encountered
  sed -i 's/\x0/?/g' cp2k_input.xml # replace null bytes which would crash saxon
  SAXON="java -jar /usr/share/java/Saxon-HE.jar"
  $SAXON -o:index.html ./cp2k_input.xml ${TOOLS}/manual/cp2k_input.xsl add_edit_links=yes
  $SAXON -o:cp2k.vim ./cp2k_input.xml ${TOOLS}/input_editing/vim/vim.xsl
)
EXIT_CODE=$?

echo ""
if ((EXIT_CODE)); then
  echo "Summary: Something is wrong."
  echo "Status: FAILED"
else
  echo "Summary: Manual generation works fine."
  echo "Status: OK"
fi

#EOF
