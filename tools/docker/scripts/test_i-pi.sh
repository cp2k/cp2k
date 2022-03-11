#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Compiling CP2K =========="
cd /opt/cp2k
echo -n "Compiling cp2k... "
if make -j VERSION=sdbg &> make.out; then
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
  python3 \
  python3-pip \
  python3-wheel \
  python3-setuptools
rm -rf /var/lib/apt/lists/*

echo -e "\n========== Installing i-Pi =========="
git clone --quiet --depth=1 --single-branch -b master https://github.com/i-pi/i-pi.git /opt/i-pi
cd /opt/i-pi
pip3 install --quiet .

echo -e "\n========== Running i-Pi Tests =========="

cd /opt/i-pi/examples/cp2k/nvt-cl
set +e # disable error trapping for remainder of script

TIMEOUT_SEC="300"
ulimit -t ${TIMEOUT_SEC} # Limit cpu time.

# launch cp2k
(
  mkdir -p run_1
  cd run_1
  echo 42 > cp2k_exit_code
  sleep 10 # give i-pi some time to startup
  export OMP_NUM_THREADS=2
  /opt/cp2k/exe/local/cp2k.sdbg ../in.cp2k
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
if ((IPI_EXIT_CODE)) || ((CP2K_EXIT_CODE)); then
  echo -e "\nnSummary: Something is wrong with i-Pi commit ${IPI_REVISION}."
  echo -e "Status: FAILED\n"
else
  echo -e "\nSummary: i-Pi commit ${IPI_REVISION} works fine."
  echo -e "Status: OK\n"
fi

#EOF
