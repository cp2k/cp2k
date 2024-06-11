#!/bin/bash -e

# author: Sebastian Seidenath
# Based on a script by Ole Schuett

# Compile CP2K.
./build_cp2k_cmake.sh "ubuntu" || exit 0

echo -e "\n========== Running i-Pi Protocol Tests =========="

export CP2K_DATA_DIR="/opt/cp2k/data"
TIMEOUT_SEC="300"
ulimit -t ${TIMEOUT_SEC} # Limit cpu time.
export OMP_NUM_THREADS=2

# launch cp2k in client mode
(
  mkdir -p run_client
  cd run_client
  echo 42 > cp2k_client_exit_code
  sleep 10 # give server some time to startup
  /opt/cp2k/exe/local/cp2k.ssmp /opt/cp2k/tests/i-PI/ipi_client.inp
  echo $? > cp2k_client_exit_code
) &

# launch cp2k in server mode
mkdir -p run_server
cd run_server
/opt/cp2k/exe/local/cp2k.ssmp /opt/cp2k/tests/i-PI/ipi_server.inp
SERVER_EXIT_CODE=$?

wait # for cp2k client to shutdown
cd /opt/cp2k
CLIENT_EXIT_CODE=$(cat ./run_client/cp2k_client_exit_code)

echo ""
echo "Client CP2K exit code: ${CLIENT_EXIT_CODE}"
echo "Server CP2K exit code: ${SERVER_EXIT_CODE}"

if ((SERVER_EXIT_CODE)) || ((CLIENT_EXIT_CODE)); then
  echo -e "\nSummary: Something is wrong with i-Pi functionality."
  echo -e "Status: FAILED\n"
else
  echo -e "\nSummary: i-Pi communication works fine."
  echo -e "Status: OK\n"
fi

#EOF
