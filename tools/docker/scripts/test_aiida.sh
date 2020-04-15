#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

echo -e "\n========== Compiling CP2K =========="
cd /workspace/cp2k
echo -n "Compiling cp2k... "
if make -j VERSION=pdbg &> /dev/null ; then
   echo "done."
else
   echo -e "failed.\n\n"
   echo "Summary: Compilation failed."
   echo "Status: FAILED"
   exit
fi

echo -e "\n========== Installing AiiDA-CP2K plugin =========="
cd /opt/aiida-cp2k/
git pull
pip3 install ./[test]

echo -e "\n========== Configuring AiiDA =========="
AS_UBUNTU_USER="sudo -u ubuntu -H"

#update reentry cache
$AS_UBUNTU_USER reentry scan

# start RabbitMQ
service rabbitmq-server start

# start and configure PostgreSQL
service postgresql start

# setup code
cat > /usr/bin/cp2k << EndOfMessage
#!/bin/bash -e
source /opt/cp2k-toolchain/install/setup
/workspace/cp2k/exe/local/cp2k.pdbg "\$@"
EndOfMessage
chmod +x /usr/bin/cp2k

echo -e "\n========== Running AiiDA-CP2K Tests =========="
cd  /opt/aiida-cp2k/

set +e # disable error trapping for remainder of script
(
set -e # abort on error
ulimit -t 1800  # abort after 30 minutes
$AS_UBUNTU_USER py.test
)

EXIT_CODE=$?

echo ""

AIIDA_COMMIT=$(git rev-parse --short HEAD)
if (( EXIT_CODE )); then
    echo "Summary: Something is wrong with aiida-cp2k commit ${AIIDA_COMMIT}."
    echo "Status: FAILED"
else
    echo "Summary: aiida-cp2k commit ${AIIDA_COMMIT} works fine."
    echo "Status: OK"
fi

#EOF
