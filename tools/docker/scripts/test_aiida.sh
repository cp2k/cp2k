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
sudo -u postgres psql -d template1 -c "CREATE USER aiida WITH PASSWORD 'aiida_db_passwd';"
sudo -u postgres psql -d template1 -c "CREATE DATABASE aiidadb OWNER aiida;"
sudo -u postgres psql -d template1 -c "GRANT ALL PRIVILEGES ON DATABASE aiidadb to aiida;"

# setup aiida user
$AS_UBUNTU_USER verdi setup                       \
      --non-interactive                           \
      --email aiida@localhost                     \
      --first-name Some                           \
      --last-name Body                            \
      --institution XYZ                           \
      --db-backend django                         \
      --db-username aiida                         \
      --db-password aiida_db_passwd               \
      --db-name aiidadb                           \
      --db-host localhost                         \
      --db-port 5432                              \
      --repository /home/ubuntu/aiida_repository  \
      --profile default

# start aiida daemon
$AS_UBUNTU_USER verdi profile setdefault default
$AS_UBUNTU_USER verdi daemon start

# setup local computer
$AS_UBUNTU_USER mkdir -p /home/ubuntu/aiida_run

$AS_UBUNTU_USER verdi computer setup    \
      --non-interactive                 \
      --label localhost                 \
      --hostname localhost              \
      --transport local                 \
      --scheduler direct                \
      --work-dir /home/ubuntu/aiida_run

$AS_UBUNTU_USER verdi computer configure local localhost --non-interactive --safe-interval 0.0
$AS_UBUNTU_USER verdi computer test localhost

# setup code
cat > /usr/bin/cp2k << EndOfMessage
#!/bin/bash -e
source /opt/cp2k-toolchain/install/setup
/workspace/cp2k/exe/local/cp2k.pdbg "\$@"
EndOfMessage
chmod +x /usr/bin/cp2k

$AS_UBUNTU_USER verdi code setup        \
      --non-interactive                 \
      --label cp2k                      \
      --computer localhost              \
      --remote-abs-path /usr/bin/cp2k   \
      --input-plugin cp2k

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
