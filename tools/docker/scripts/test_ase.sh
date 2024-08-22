#!/bin/bash -e

# author: Ole Schuett

# Compile CP2K.
./build_cp2k_cmake.sh "ubuntu" "ssmp" || exit 0

# Fake installation of data files.
mkdir -p ./share/cp2k
ln -s ../../data ./share/cp2k/data

cat > /usr/bin/cp2k_shell << EndOfMessage
#!/bin/bash -e
export OMP_NUM_THREADS=1
/opt/cp2k/build/bin/cp2k.ssmp --shell "\$@"
EndOfMessage
chmod +x /usr/bin/cp2k_shell

# The cp2k main binary is used by ase/test/cp2k/cp2k_dcd.py.
# https://gitlab.com/ase/ase/merge_requests/1109
cat > /usr/bin/cp2k << EndOfMessage
#!/bin/bash -e
export OMP_NUM_THREADS=1
/opt/cp2k/build/bin/cp2k.ssmp "\$@"
EndOfMessage
chmod +x /usr/bin/cp2k

mkdir -p ~/.config/ase
cat > ~/.config/ase/config.ini << EndOfMessage
[cp2k]
cp2k_shell = /usr/bin/cp2k_shell
cp2k_main = /usr/bin/cp2k
EndOfMessage

echo -e "\n========== Installing Dependencies =========="
apt-get update -qq
apt-get install -qq --no-install-recommends \
  git \
  python3 \
  python3-dev \
  python3-venv \
  python3-pip \
  python3-wheel \
  python3-setuptools \
  build-essential
rm -rf /var/lib/apt/lists/*

# Create and activate a virtual environment for Python packages.
python3 -m venv /opt/venv
export PATH="/opt/venv/bin:$PATH"

echo -e "\n========== Installing ASE =========="
git clone --quiet --depth=1 --single-branch -b master https://gitlab.com/ase/ase.git /opt/ase
cd /opt/ase/
pip3 install ".[test]"

echo -e "\n========== Running ASE Tests =========="
ASE_REVISION=$(git rev-parse --short HEAD)
echo -

# Make test temp files available as artifacts.
export PYTEST_DEBUG_TEMPROOT=/workspace/artifacts
mkdir -p ${PYTEST_DEBUG_TEMPROOT}

if ase test -j 0 -c cp2k calculator/cp2k; then
  echo -e "\nSummary: ASE commit ${ASE_REVISION} works fine."
  echo -e "Status: OK\n"
else
  echo -e "\nSummary: Something is wrong with ASE commit ${ASE_REVISION}."
  echo -e "Status: FAILED\n"
fi

#EOF
