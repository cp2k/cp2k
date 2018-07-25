#!/bin/bash

# author: Ole Schuett

set -e

echo -e "\n========== Copying Changed Files =========="
rsync --exclude="*~"          \
      --exclude=".*/"         \
      --exclude="*.pyc"       \
      --exclude=/cp2k/obj/    \
      --exclude=/cp2k/lib/    \
      --exclude=/cp2k/exe/    \
      --ignore-times          \
      --update                \
      --verbose               \
      --recursive             \
      --checksum              \
      /opt/cp2k-local/  /opt/cp2k-master/

rsync --exclude="*~"          \
      --exclude=".*/"         \
      --exclude="*.pyc"       \
      --ignore-times          \
      --update                \
      --verbose               \
      --recursive             \
      --checksum              \
      /opt/cp2k-local/cp2k/tools/toolchain/  /opt/cp2k-toolchain/

echo -e "\n========== Updating Toolchain =========="
cd /opt/cp2k-toolchain/
./install_cp2k_toolchain.sh --install-all --with-make=no
source ./install/setup

echo -e "\n========== Compiling CP2K =========="
cd /opt/cp2k-master/cp2k/makefiles
echo ${CP2K_REVISION} > ../REVISION
touch ../src/cp2k_info.F  # ensure latest REVISION is picked up.
make -j VERSION="sopt" cp2k

echo -e "\n========== Generating Manual =========="

mkdir -p /opt/cp2k_test_artifacts/manual
cd /opt/cp2k_test_artifacts/manual

/opt/cp2k-master/cp2k/exe/local/cp2k.sopt --version
/opt/cp2k-master/cp2k/exe/local/cp2k.sopt --xml

TOOLS=/opt/cp2k-master/cp2k/tools
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
if (( $EXIT_CODE )); then
    echo "Summary: Something is wrong."
    echo "Status: FAILED"
else
    echo "Summary: Manual generation works fine."
    echo "Status: OK"
fi

#EOF
