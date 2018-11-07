#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Compiling CP2K =========="
cd /workspace/cp2k
touch src/cp2k_info.F  # ensure latest REVISION is picked up.
make -j VERSION="sopt" cp2k

echo -e "\n========== Generating Manual =========="

mkdir -p /workspace/artifacts/manual
cd /workspace/artifacts/manual

/workspace/cp2k/exe/local/cp2k.sopt --version
/workspace/cp2k/exe/local/cp2k.sopt --xml

TOOLS=/workspace/cp2k/tools
cp ${TOOLS}/manual/favicon.png .
cp ${TOOLS}/manual/toggle_folding.js .

set +e # disable error trapping for remainder of script
(
set -e # abort if error is encountered
SAXON="java -jar /usr/share/java/Saxon-HE.jar"
$SAXON -o:index.html ./cp2k_input.xml ${TOOLS}/manual/cp2k_input.xsl add_edit_links=yes
$SAXON -o:cp2k.vim   ./cp2k_input.xml ${TOOLS}/input_editing/vim/vim.xsl
)
EXIT_CODE=$?

echo ""
if (( EXIT_CODE )); then
    echo "Summary: Something is wrong."
    echo "Status: FAILED"
else
    echo "Summary: Manual generation works fine."
    echo "Status: OK"
fi

#EOF
