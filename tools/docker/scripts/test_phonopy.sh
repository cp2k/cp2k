#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

# Compile CP2K.
./build_cp2k_cmake.sh "toolchain_all" "ssmp" || exit 0

# Fake installation of data files.
mkdir -p ./share/cp2k
ln -s ../../data ./share/cp2k/data

echo -e "\n========== Installing Dependencies =========="
apt-get update -qq
apt-get install -qq --no-install-recommends \
  git \
  python3 \
  python3-venv \
  python3-pip \
  python3-wheel \
  python3-setuptools
rm -rf /var/lib/apt/lists/*

# Create and activate a virtual environment for Python packages.
python3 -m venv /opt/venv
export PATH="/opt/venv/bin:$PATH"

echo -e "\n========== Installing Phonopy =========="
git clone --quiet --depth=1 --single-branch -b develop https://github.com/phonopy/phonopy.git /opt/phonopy
cd /opt/phonopy/
pip3 install ".[cp2k]"

# Workaround https://github.com/hgrecco/pint/issues/1974
pip3 install Pint==0.24.4

echo -e "\n========== Running Phonopy Test =========="
mkdir tmp
cd tmp
cp ../example/Si-CP2K/Si.inp .

set +e # disable error trapping for remainder of script
(
  set -e # abort on error
  # Following example/Si-CP2K/README.md
  phonopy --cp2k -c Si.inp -d --dim="2 2 2"
  OMP_NUM_THREADS=2 /opt/cp2k/build/bin/cp2k.ssmp Si-supercell-001.inp &> cp2k.out
  phonopy --cp2k -f Si-supercell-001-forces-1_0.xyz
)
EXIT_CODE=$?

PHONOPY_REVISION=$(git rev-parse --short HEAD)
if ((EXIT_CODE)); then
  echo -e "\nSummary: Something is wrong with Phonopy commit ${PHONOPY_REVISION}."
  echo -e "Status: FAILED\n"
else
  echo -e "\nSummary: Phonopy commit ${PHONOPY_REVISION} works fine."
  echo -e "Status: OK\n"
fi

#EOF
