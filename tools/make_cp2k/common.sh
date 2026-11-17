# shellcheck shell=bash disable=SC2034,SC2153,SC2154
#
# make_cp2k.sh module: Environment checks and helper functions
#
# This file is sourced by make_cp2k.sh at top level and shares all variables
# with the other modules. It is not meant to be executed on its own.

# Check platform
if [[ "$(uname -s)" == "Darwin" ]]; then
  echo "ERROR: A build for Apple macOS (Darwin, Homebrew) is not supported yet,"
  echo "       but CP2K can be built within a container using podman."
  echo "       Check the Dockerfiles in the folder cp2k/tools/docker like"
  echo "       Dockerfile.test_spack_openmpi-psmp and the usage description in their header"
  exit 1
fi

# Help finding ldconfig
if ! command -v ldconfig &> /dev/null; then
  PATH="${PATH}:/sbin"
fi

# Check if all mandatory packages are installed in the environment
for package in awk bzip2 find g++ gcc gfortran git gzip ldconfig make patch python3 tar wget xz; do
  if ! command -v "${package}" &> /dev/null; then
    echo "ERROR: The package \"${package}\" is mandatory to build CP2K with Spack/CMake"
    echo "       Install the missing package and re-run the script"
    exit 1
  fi
done

# Disable python bytecode generation to suppress __pycache__ folder creation
export PYTHONDONTWRITEBYTECODE=1

# Abort git transfers slower than 1 KB/s for 60 s, since git has no timeout,
# so that a stalled clone (e.g. py-torch with submodules) fails and is retried
export GIT_HTTP_LOW_SPEED_LIMIT=1000
export GIT_HTTP_LOW_SPEED_TIME=60

# Run a command that downloads data with up to 3 attempts and a backoff
# of 60 s and 120 s to bridge transient network issues
retry_download() {
  local attempt delay=60
  for attempt in 1 2 3; do
    "$@" && return 0
    if ((attempt < 3)); then
      echo "WARNING: Attempt ${attempt} of 3 failed, retrying in ${delay} seconds: $*"
      sleep "${delay}"
      delay=$((delay * 2))
    fi
  done
  return 1
}

# Print OS info if available
if [[ -f /etc/os-release ]]; then
  echo ""
  cat /etc/os-release
fi

# Check bash version
if ((BASH_VERSINFO < 4)); then
  echo ""
  echo "ERROR: The employed bash version ${BASH_VERSION} is too old (from 2004)"
  echo "       Install a newer bash version"
  if [[ "$(uname -s)" == "Darwin" ]]; then
    echo "HINT:  Use the bash command provided by Homebrew with Apple"
    echo "       /opt/homebrew/bin should precede /bin/bash in the PATH"
  fi
  exit 1
elif ((BASH_VERSINFO < 5)); then
  echo -e "\nWARNING: The employed bash version ${BASH_VERSION} is quite old (from 2009)"
else
  echo -e "\nINFO: Using bash version ${BASH_VERSION}"
fi

# Check if the python3 version is new enough for spack
if ! python3 -c "import sys; sys.exit(not(sys.version_info >= (3, 10)))"; then
  echo ""
  echo "ERROR: Python version is NOT >= 3.10 (needed for Spack)"
  echo "       Found only $(python3 -V)"
  exit 1
else
  echo -e "\nINFO: Found $(python3 -V)"
fi

# Convert a CMake build type to its canonical spelling (case-insensitive input),
# unknown build types are returned unchanged
canonical_build_type() {
  case "${1,,}" in
    debug) echo "Debug" ;;
    release) echo "Release" ;;
    relwithdebinfo) echo "RelWithDebInfo" ;;
    *) echo "${1}" ;;
  esac
}
