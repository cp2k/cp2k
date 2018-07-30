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

echo -e "\n========== Compiling CP2K =========="
source /opt/cp2k-toolchain/install/setup
cd /opt/cp2k-master/cp2k/makefiles
make -j VERSION=pdbg cp2k

echo -e "\n========== Installing CP2K =========="
cat > /usr/bin/cp2k << EndOfMessage
#!/bin/bash -e
source /opt/cp2k-toolchain/install/setup
/opt/cp2k-master/cp2k/exe/local/cp2k.pdbg "\$@"
EndOfMessage
chmod +x /usr/bin/cp2k

echo -e "\n========== Installing AiiDA-CP2K plugin =========="
cd /opt/
git clone https://github.com/cp2k/aiida-cp2k.git
pip install ./aiida-cp2k/

echo -e "\n========== Configuring AiiDA =========="
for i in $(dirname $(which mpirun))/* ; do ln -sf $i /usr/bin/; done
SUDO="sudo -u ubuntu -H"
cd /opt/aiida-cp2k/test/
$SUDO ./configure_aiida.sh

echo -e "\n========== Running AiiDA-CP2K Tests =========="
cd  /opt/aiida-cp2k/test/

set +e # disable error trapping for remainder of script
(
set -e # abort on error
$SUDO ./run_tests.sh
)

EXIT_CODE=$?

echo ""

AIIDA_COMMIT=`git rev-parse --short HEAD`
if (( $EXIT_CODE )); then
    echo "Summary: Something is wrong with aiida-cp2k commit ${AIIDA_COMMIT}."
    echo "Status: FAILED"
else
    echo "Summary: aiida-cp2k commit ${AIIDA_COMMIT} works fine."
    echo "Status: OK"
fi

#EOF
