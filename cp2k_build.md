# Instructions to build and test CP2K

NOTE: Only use 16 MPI processes and 1 OpenMP thread for building and testing!

## install / check system packages
sudo apt update
sudo apt upgrade
sudo apt install build-essential gcc g++ gfortran cmake make ninja-build git pkg-config libspglib-f08-dev 
sudo apt install python3-dev bzip2 less nano ca-certificates openmpi-bin openmpi-common libopenmpi-dev

NOTE: Do NOT install system FFTW packages. The toolchain installs its own.


## run toolchain script
cd tools/toolchain
./install_cp2k_toolchain.sh -j 16 \
    --with-tblite=install \
    --with-gcc=system \
    --with-cmake=system \
    --with-ninja=system \
    --with-openblas=install \
    --with-libxc=install \
    --with-libint=install \
    --with-fftw=install \
    --with-libxsmm=install \
    --with-spglib=install \
    --with-gauxc=install

## build CP2K
source ./install/setup
./build_cp2k.sh -j 16

NOTE: never use make_cp2k.sh
NOTE: run always setup script before, see toolchain script output
NOTE: If it is a rebuild and the build system has not been modified use the --rebuild-only flag

## run regression tests
source ./install/cp2k_env
./tests/do_regtest.py ./install/bin psmp

NOTE: Must source cp2k_env before running tests to set library paths.

# Instructions to build ASE

## install / check system packages
sudo apt install python3-pip python3-numpy python3-scipy python3-matplotlib python3-tk python3-flask python3-pytest python3-spglib python3-pytest-xdist

##get ase
git clone --recursive https://gitlab.com/jpoto/ase.git

## install ASE
python -m venv ~/venvs/ase
source ~/venvs/ase/bin/activate
pip install -e /workspace
pip install pytest pytest-xdist

NOTE: Install ASE with path (not -e .) to avoid package discovery conflicts.
NOTE: Create venv outside workspace.

## run ASE tests
ase test

## run ASE CP2K test
source ~/venvs/ase/bin/activate
source /workspace/cp2k/install/cp2k_env

Create ~/.config/ase/config.ini:
```
[cp2k]
cp2k_shell = /workspace/cp2k/install/bin/cp2k.psmp -s
cp2k_main = /workspace/cp2k/install/bin/cp2k.psmp
```

OMP_NUM_THREADS=1 ase test -c cp2k

# Run only CP2K calculator tests (without regression suite)
# Use OMP_NUM_THREADS=1 to avoid threading issues
OMP_NUM_THREADS=1 python -m pytest ase/test/calculator/cp2k/ \
    --calculators=cp2k

NOTE: cp2k_shell is not a separate binary. Use cp2k.psmp -s for shell mode.

# Coding Conventions for ASE

Reference: https://wiki.fysik.dtu.dk/ase/development/python_codingstandard.html

## Key Rules

- Use 4 spaces per indentation level, no tabs
- Maximum line length is 78 characters
- Use "StudlyCaps" for class names
- Use "lowercase" or "lowercase_with_underscores" for function, method, and variable names
- Use `'single quotes'` for string literals, `"""triple double quotes"""` for docstrings
- No one-liner compound statements (no `if x: return`)
- Run ruff on your code before committing:
  ```bash
  ruff check --fix filename.py
  ruff format filename.py
  ```
- Use NumPy/SciPy convention for docstrings

