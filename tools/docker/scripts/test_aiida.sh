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
export DEBIAN_FRONTEND=noninteractive
export DEBCONF_NONINTERACTIVE_SEEN=true
apt-get install -qq --no-install-recommends \
  python3-setuptools \
  python3-wheel \
  python3-pip \
  python3-dev \
  postgresql \
  rabbitmq-server \
  sudo \
  ssh
rm -rf /var/lib/apt/lists/*

# Some buggy Python packages open utf8 files in text mode.
# As a workaround we set locale.getpreferredencoding() to utf8.
export LANG="en_US.UTF-8" LANGUAGE="en_US:en" LC_ALL="en_US.UTF-8"

# create ubuntu user with sudo powers
adduser --disabled-password --gecos "" ubuntu
echo "ubuntu ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers

# link mpi executables into path
MPI_INSTALL_DIR=$(dirname "$(command -v mpirun)")
for i in "${MPI_INSTALL_DIR}"/*; do ln -sf "$i" /usr/bin/; done

echo -e "\n========== Installing AiiDA-CP2K plugin =========="
git clone --quiet https://github.com/aiidateam/aiida-cp2k.git /opt/aiida-cp2k/
cd /opt/aiida-cp2k/
pip3 install './[test]'

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
/opt/cp2k/exe/local/cp2k.sdbg "\$@"
EndOfMessage
chmod +x /usr/bin/cp2k

echo -e "\n========== Running AiiDA-CP2K Tests =========="
set +e # disable error trapping for remainder of script
(
  set -e         # abort on error
  ulimit -t 1800 # abort after 30 minutes
  $AS_UBUNTU_USER py.test
)

EXIT_CODE=$?

AIIDA_COMMIT=$(git rev-parse --short HEAD)
if ((EXIT_CODE)); then
  echo -e "\nSummary: Something is wrong with aiida-cp2k commit ${AIIDA_COMMIT}."
  echo -e "Status: FAILED\n"
else
  echo -e "\nSummary: aiida-cp2k commit ${AIIDA_COMMIT} works fine."
  echo -e "Status: OK\n"
fi

#EOF
