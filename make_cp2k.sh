#!/usr/bin/env bash

# Purpose: Build CP2K using Spack and CMake locally within the folder CP2K_ROOT
#          which defaults to the current working directory. This should be a
#          "cp2k/" folder containing the CP2K source tree which can be download with
#
#          > git clone https://github.com/cp2k/cp2k.git cp2k
#
#          and after
#
#          > cd cp2k
#
#          this script can be run in subshell using the default options with
#
#          > ./make_cp2k.sh
#
#          The flags -h or --help print the available options.
#
#          The first run will take longer as it will build all CP2K dependencies
#          with Spack. The Spack installation is kept fully local in the subfolder
#          cp2k/spack which corresponds to the tools/toolchain/install folder
#          created by the CP2K toolchain.
#
#          Subsequent runs of the script will use the software stack from that
#          cp2k/spack/ folder.
#
#          A rebuild of all CP2K dependencies can be enforced simply by removing
#          or renaming the folder cp2k/spack. The latter allows for keeping different
#          software stacks (see also -bd and -bd_only flags).
#
#          It is recommended to store already compiled packages in a cache (default).
#          This will accelerate the (re)build of the CP2K dependencies with Spack
#          significantly (see -uc flag).
#
#          After the CP2K dependencies are built with Spack, CP2K itself is built
#          and installed using CMake in the subfolders cp2k/build and cp2k/install,
#          respectively.
#
#          Subsequent runs of the script will use the CMake configuration in the
#          subfolder cp2k/build. A rebuild of CP2K from scratch can be enforced
#          by removing or renaming that subfolder.
#
#          A CP2K regression run can be launched automatically by adding the flag
#          -t "" (or --test ""). This flag expects a string with the TESTOPTS, e.g.
#
#          > ./make_cp2k.sh -t "--maxtasks 8 --restrictdir QS/regtest-gpw-1"
#
#          Alternatively, the script cp2k/install/bin/run_tests can be launched
#          after a successful CP2K build. Usage examples are printed at the end.
#
#          The actual work is performed by the modules in the folder
#          cp2k/tools/make_cp2k which are sourced by this script.

# Authors: Matthias Krack (MK)

# Version: 2.5

# The script has to be executed and must not be sourced, because it would
# otherwise modify the calling shell (shell options, variables, exit)
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
  echo "ERROR: $(basename "${BASH_SOURCE[0]}") must be executed, not sourced"
  return 1
fi

# Facilitate the deugging of this script
set -uo pipefail

# Retrieve script name (used by the modules)
# shellcheck disable=SC2034
SCRIPT_NAME="$(basename "${BASH_SOURCE[0]}")"

# Locate the modules of this script, which are sourced in the following order
MAKE_CP2K_MODULE_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/tools/make_cp2k"
MAKE_CP2K_MODULES="common options settings help checks spack cmake install tests gromacs ase"
for module_name in ${MAKE_CP2K_MODULES}; do
  if [[ ! -f "${MAKE_CP2K_MODULE_PATH}/${module_name}.sh" ]]; then
    echo "ERROR: The module file ${MAKE_CP2K_MODULE_PATH}/${module_name}.sh was not found"
    exit 1
  fi
done
unset module_name

# The modules have to be sourced at top level (not within a function),
# because they share variables including associative arrays

# Environment checks and helper functions (retry_download, canonical_build_type)
# shellcheck source=/dev/null
source "${MAKE_CP2K_MODULE_PATH}/common.sh"

# Default values and parsing of the command line flags
# shellcheck source=/dev/null
source "${MAKE_CP2K_MODULE_PATH}/options.sh"

# Derived settings: CMake feature flags and build modes (ASAN, conventions, coverage)
# shellcheck source=/dev/null
source "${MAKE_CP2K_MODULE_PATH}/settings.sh"

# Usage information (-h | --help)
# shellcheck source=/dev/null
source "${MAKE_CP2K_MODULE_PATH}/help.sh"

# Settings printout, consistency checks, CUDA settings, and available resources
# shellcheck source=/dev/null
source "${MAKE_CP2K_MODULE_PATH}/checks.sh"

# Build of the CP2K dependencies with Spack
# shellcheck source=/dev/null
source "${MAKE_CP2K_MODULE_PATH}/spack.sh"

# CMake configuration, build, and installation of CP2K
# shellcheck source=/dev/null
source "${MAKE_CP2K_MODULE_PATH}/cmake.sh"

# Runtime setup: library paths, launch and run_tests scripts, benchmark script
# shellcheck source=/dev/null
source "${MAKE_CP2K_MODULE_PATH}/install.sh"

# Optional regression tests, coverage report, and performance run
# shellcheck source=/dev/null
source "${MAKE_CP2K_MODULE_PATH}/tests.sh"

# Optional GROMACS/CP2K build and test
# shellcheck source=/dev/null
source "${MAKE_CP2K_MODULE_PATH}/gromacs.sh"

# Optional ASE build and test
# shellcheck source=/dev/null
source "${MAKE_CP2K_MODULE_PATH}/ase.sh"
