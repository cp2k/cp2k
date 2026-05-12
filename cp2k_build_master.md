# CP2K Build Instructions

### 1. Install system packages

```bash
sudo apt update
sudo apt upgrade
sudo apt install -y cmake ninja-build make bzip2 xz-utils patch libgsl-dev pkg-config libeigen3-dev 
sudo apt install -y libfftw3-dev libfftw3-mpi-dev libopenmpi-dev gfortran gcc g++ libscalapack-openmpi-dev libscalapack-mpi-dev
```

### 2. Checkout branch / set up git

configure:
git config --global user.name "Johann Pototschnig"
git config --global user.email "j.pototschnig@hzdr.de"
git config --local user.name "Johann Pototschnig"
git config --local user.email "j.pototschnig@hzdr.de"
export GIT_COMMITTER_NAME="Johann Pototschnig"
export GIT_COMMITTER_EMAIL="j.pototschnig@hzdr.de"
use this always as committer and author

If it is a new setup:
```bash
cd /workspace
git checkout master
create new branch with jvp and terok in the name
```
otherwise stay on branch

### 3. Toolchain Installation

Run the toolchain installation script from the correct directory:

```bash
cd /workspace/tools/toolchain
./install_cp2k_toolchain.sh --with-openmpi=system --with-gcc=system --with-cmake=system --with-ninja=system --with-openblas=install --with-tblite=install
```

Note: Must run from `/workspace/tools/toolchain` directory (not from `/workspace`).
Note: Increase timeout time.
Note: If there is timeout or error delete the latest directory in the /tools/toolchain/install folder

### 4. Build CP2K

Always use build_cp2k.sh

Run the make_pretty.sh script
Always run the toolchain script first to set up the correct environment. 
Follow build instructions from the ouput of the toolchain script. 
If there is an error look into the build.log and try to fix it before recompilation. 

NOTE: Never use make_cp2k.sh


### 5. Running Regression Tests

Source the toolchain setup to ensure libraries are found:

```bash
source /workspace/tools/toolchain/install/setup
```

Run the regression tests:

```bash
cd /workspace/tests
python3 do_regtest.py /workspace/build/bin psmp
```

### REMARK:
ALways build and run tests before pushing changes

