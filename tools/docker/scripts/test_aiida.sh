#!/bin/bash -e

# author: Ole Schuett

# Compile CP2K.
./build_cp2k_cmake.sh "ubuntu" "ssmp" || exit 0

echo -e "\n========== Installing Dependencies =========="
apt-get update -qq
export DEBIAN_FRONTEND=noninteractive
export DEBCONF_NONINTERACTIVE_SEEN=true
apt-get install -qq --no-install-recommends \
  python3-setuptools \
  python3-wheel \
  python3-pip \
  python3-venv \
  python3-dev \
  python3-reentry \
  postgresql \
  libpq-dev \
  rabbitmq-server \
  sudo \
  git \
  ssh
rm -rf /var/lib/apt/lists/*

# Create and activate a virtual environment for Python packages.
python3 -m venv /opt/venv
export PATH="/opt/venv/bin:$PATH"

# Some buggy Python packages open utf8 files in text mode.
# As a workaround we set locale.getpreferredencoding() to utf8.
export LANG="en_US.UTF-8" LANGUAGE="en_US:en" LC_ALL="en_US.UTF-8"
locale-gen ${LANG}

# link mpi executables into path
MPI_INSTALL_DIR=$(dirname "$(command -v mpiexec)")
for i in "${MPI_INSTALL_DIR}"/*; do ln -sf "$i" /usr/bin/; done

# Pick a compiler (needed to build some Python packages)
export CC=gcc

echo -e "\n========== Installing AiiDA-CP2K plugin =========="
git clone --quiet https://github.com/aiidateam/aiida-cp2k.git /opt/aiida-cp2k/
cd /opt/aiida-cp2k/
pip3 install './[dev]'

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
export OMP_NUM_THREADS=2
/opt/cp2k/build/bin/cp2k.ssmp "\$@"
EndOfMessage
chmod +x /usr/bin/cp2k

echo -e "\n========== Running AiiDA-CP2K Tests =========="
set +e # disable error trapping for remainder of script
(
  set -e         # abort on error
  ulimit -t 1800 # abort after 30 minutes
  $AS_UBUNTU_USER /opt/venv/bin/py.test
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
