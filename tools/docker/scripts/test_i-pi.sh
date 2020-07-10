#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Compiling CP2K =========="
cd /workspace/cp2k
echo -n "Compiling cp2k... "
if make -j VERSION=pdbg &> make.out ; then
    echo "done."
else
    echo -e "failed.\n\n"
    tail -n 100 make.out
    echo -e "\nSummary: Compilation failed."
    echo -e "Status: FAILED\n"
    exit 0
fi

echo -e "\n========== Installing i-Pi =========="
cd /opt/i-pi
git pull --quiet
pip install --quiet .

echo -e "\n========== Running i-Pi Tests =========="

cd  /opt/i-pi/examples/cp2k/nvt-cl
set +e # disable error trapping for remainder of script

TIMEOUT_SEC="300"
ulimit -t ${TIMEOUT_SEC}  # Limit cpu time.

# launch cp2k
(
  mkdir -p run_1
  cd run_1
  echo 42 > cp2k_exit_code
  sleep 10 # give i-pi some time to startup
  export OMP_NUM_THREADS=2
  mpiexec -np 2 /workspace/cp2k/exe/local/cp2k.pdbg ../in.cp2k
  echo $? > cp2k_exit_code
) &

# launch i-pi
sed -i "s/total_steps>1000/total_steps>10/" input.xml
# Limit walltime too, because waiting for a connection consumes no cpu time.
timeout ${TIMEOUT_SEC} /usr/local/bin/i-pi input.xml
IPI_EXIT_CODE=$?

wait # for cp2k to shutdown
CP2K_EXIT_CODE=$(cat ./run_1/cp2k_exit_code)

echo ""
echo "CP2K exit code: ${CP2K_EXIT_CODE}"
echo "i-Pi exit code: ${IPI_EXIT_CODE}"

IPI_REVISION=$(git rev-parse --short HEAD)
if (( IPI_EXIT_CODE )) || (( CP2K_EXIT_CODE )) ; then
    echo "Summary: Something is wrong with i-Pi commit ${IPI_REVISION}."
    echo "Status: FAILED"
else
    echo "Summary: i-Pi commit ${IPI_REVISION} works fine."
    echo "Status: OK"
fi

#EOF
