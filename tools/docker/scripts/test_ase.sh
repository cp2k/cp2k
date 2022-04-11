#!/bin/bash -e

# author: Ole Schuett

# shellcheck disable=SC1091
source /opt/cp2k-toolchain/install/setup

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

cat > /usr/bin/cp2k_shell << EndOfMessage
#!/bin/bash -e
source /opt/cp2k-toolchain/install/setup
export OMP_NUM_THREADS=1
/opt/cp2k/exe/local/cp2k_shell.sdbg "\$@"
EndOfMessage
chmod +x /usr/bin/cp2k_shell

# The cp2k main binary is used by ase/test/cp2k/cp2k_dcd.py.
# https://gitlab.com/ase/ase/merge_requests/1109
cat > /usr/bin/cp2k << EndOfMessage
#!/bin/bash -e
source /opt/cp2k-toolchain/install/setup
export OMP_NUM_THREADS=1
/opt/cp2k/exe/local/cp2k.sdbg "\$@"
EndOfMessage
chmod +x /usr/bin/cp2k

mkdir -p ~/.config/ase
cat > ~/.config/ase/ase.conf << EndOfMessage
[executables]
cp2k = /usr/bin/cp2k_shell
cp2k_main = /usr/bin/cp2k
EndOfMessage

echo -e "\n========== Installing Dependencies =========="
apt-get update -qq
apt-get install -qq --no-install-recommends \
  python3 \
  python3-dev \
  python3-pip \
  python3-wheel \
  python3-setuptools \
  build-essential
rm -rf /var/lib/apt/lists/*

echo -e "\n========== Installing ASE =========="
git clone --quiet --depth=1 --single-branch -b master https://gitlab.com/ase/ase.git /opt/ase
cd /opt/ase/
pip3 install ".[test]"

echo -e "\n========== Running ASE Tests =========="
ASE_REVISION=$(git rev-parse --short HEAD)
echo -

if ase test -j 0 -c cp2k calculator/cp2k; then
  echo -e "\nSummary: ASE commit ${ASE_REVISION} works fine."
  echo -e "Status: OK\n"
else
  echo -e "\nSummary: Something is wrong with ASE commit ${ASE_REVISION}."
  echo -e "Status: FAILED\n"
fi

#EOF
