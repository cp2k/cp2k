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

echo -e "\n========== Installing i-Pi =========="
cd /opt/i-pi
git pull
pip install .

echo -e "\n========== Running i-Pi Tests =========="

cd  /opt/i-pi/examples/cp2k/nvt-cl
set +e # disable error trapping for remainder of script

# launch cp2k
(
  mkdir -p run_1
  cd run_1
  echo 42 > cp2k_exit_code
  sleep 2 # give i-pi some time to startup
  mpiexec -np 2 /opt/cp2k-master/cp2k/exe/local/cp2k.pdbg ../in.cp2k
  echo $? > cp2k_exit_code
) &

# launch i-pi
sed -i "s/total_steps>1000/total_steps>10/" input.xml
/usr/local/bin/i-pi input.xml
IPI_EXIT_CODE=$?

wait # for cp2k to shutdown
CP2K_EXIT_CODE=`cat ./run_1/cp2k_exit_code`

echo ""
echo "CP2K exit code: " $CP2K_EXIT_CODE
echo "i-Pi exit code: " $IPI_EXIT_CODE

IPI_REVISION=`git rev-parse --short HEAD`
if (( $IPI_EXIT_CODE || $CP2K_EXIT_CODE )); then
    echo "Summary: Something is wrong with i-Pi commit ${IPI_REVISION}."
    echo "Status: FAILED"
else
    echo "Summary: i-Pi commit ${IPI_REVISION} works fine."
    echo "Status: OK"
fi

#EOF

