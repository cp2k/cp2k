#!/usr/bin/env bash
#
# build_deps.sh: Build the CP2K dependencies with Spack (step 2 of make_cp2k.sh)
#
# Usage: build_deps.sh CONFIG_FILE...
#
# The settings are read from the configuration file options.conf written by
# make_cp2k.sh (see also make_cp2k.sh --print_config). The dependencies are
# built in the folder SPACK_BUILD_PATH only if they are not complete yet. The
# settings provided by the dependency build are written to CONFIG_PATH/deps.conf.

# shellcheck disable=SC2154 # the variables are defined by the configuration files

set -uo pipefail

# shellcheck source=/dev/null
source "$(dirname "$(realpath "${BASH_SOURCE[0]}")")/common.sh"
load_config "$@"

### Build CP2K dependencies with Spack if needed or requested ###

# An incomplete build already contains the Spack installation, source cache,
# environment, and every successfully installed dependency. Re-enter the
# setup/install path until the completion marker is written.
if [[ ! -f "${SPACK_BUILD_PATH}/BUILD_DEPENDENCIES_COMPLETED" ]]; then

  if [[ -d "${SPACK_BUILD_PATH}" ]]; then
    echo "INFO: Resuming incomplete CP2K dependency build in ${SPACK_BUILD_PATH}"
  else
    echo "INFO: Creating local Spack build directory ${SPACK_BUILD_PATH}"
  fi
  mkdir -p "${SPACK_BUILD_PATH}"
  cd "${SPACK_BUILD_PATH}" || exit 1

  # Reset the spack environment
  unset SPACK_ENV

  # Install Spack only when the local bootstrap is absent. wget --continue
  # also retains a partially downloaded release archive across interruptions.
  if [[ ! -d "${SPACK_ROOT}" ]]; then
    echo "Installing Spack ${SPACK_VERSION}"
    if ! retry_download wget -nv -c --timeout=60 --tries=2 "https://github.com/spack/spack/releases/download/v${SPACK_VERSION}/spack-${SPACK_VERSION}.tar.gz"; then
      echo "ERROR: Downloading Spack ${SPACK_VERSION} failed"
      exit 1
    fi
    if ! tar -xzf "spack-${SPACK_VERSION}.tar.gz"; then
      echo "ERROR: Extracting Spack ${SPACK_VERSION} failed"
      exit 1
    fi
    if ! mv -f "${SPACK_BUILD_PATH}/spack-${SPACK_VERSION}" "${SPACK_ROOT}"; then
      echo "ERROR: Installing Spack ${SPACK_VERSION} failed"
      exit 1
    fi
    rm -f "spack-${SPACK_VERSION}.tar.gz"
  fi
  export PATH="${SPACK_ROOT}/bin:${PATH}"

  # Create cache folder for user configuration
  mkdir -p "${SPACK_USER_CACHE_PATH}"

  # Set Spack cache (folder is the default)
  export SPACK_CACHE="${SPACK_CACHE:-file://${CP2K_ROOT}/spack_cache}"

  # Prepare for package caching
  if [[ "${USE_CACHE}" == @("folder"|"minio") ]]; then
    # Create the venv only once, then reuse it on a resumed dependency build
    if [[ ! -x "${SPACK_BUILD_PATH}/venv/bin/python3" ]]; then
      if command -v python3 -m venv --help &> /dev/null; then
        echo "Installing virtual environment for Python packages"
        if ! python3 -m venv "${SPACK_BUILD_PATH}/venv"; then
          echo "ERROR: The creation of a virtual environment (venv) for Python packages failed"
          exit 1
        fi
      else
        echo "ERROR: python3 -m venv was not found"
        exit 1
      fi
    fi
    export PATH="${SPACK_BUILD_PATH}/venv/bin:${PATH}"

    # Avoid contacting PyPI again when the previously prepared venv is intact
    if ! python3 -c 'from importlib.metadata import version; assert version("boto3") == "1.38.11"; assert version("google-cloud-storage") == "3.1.0"' &> /dev/null; then
      if ! python3 -m pip --version &> /dev/null; then
        echo "ERROR: python3 -m pip was not found"
        exit 1
      fi
      if ! retry_download python3 -m pip install "${VERBOSE_FLAG}" --upgrade pip \
        boto3==1.38.11 google-cloud-storage==3.1.0; then
        echo "ERROR: The installation of Python packages for the Spack cache failed"
        exit 1
      fi
    fi
  fi

  # Prepare for using an MinIO object storage as cache
  if [[ "${USE_CACHE}" == "minio" ]]; then
    if [[ "${IN_CONTAINER}" == "yes" ]]; then
      echo "ERROR: An external MinIO object storage container cannot be used as cache within a container"
      echo "       Use \"folder\" or \"no\" caching instead of \"MinIO\""
      exit 1
    fi
    # MinIO requires podman
    if command -v podman &> /dev/null; then
      # Check podman version
      ((VERBOSE > 0)) && podman version
      PODMAN_VERSION=$(podman version --format '{{.Client.Version}}' | awk -F'[.-]' '{print $1}')
      if ((PODMAN_VERSION < 4)); then
        echo "WARNING: Outdated podman version $(podman version --format '{{.Client.Version}}') found"
      fi
      export SPACK_CACHE="${SPACK_CACHE:-s3://spack-cache --s3-endpoint-url=http://localhost:9000}"
      if ! "${CP2K_ROOT}"/tools/docker/spack_cache_start.sh; then
        echo "ERROR: Could not start (new) spack cache"
        echo ""
        echo "An error message starting with \"Error: initial journal cursor: ...\" indicates that the"
        echo "journald logging does not work. Try to switch podman (3.x) to file-based logs by creating"
        echo "the file \"~/.config/containers/containers.conf\" with the following two lines:"
        echo "[containers]"
        echo "log_driver = \"k8s-file\""
        exit 1
      fi
    else
      echo "INFO: podman was not found"
      echo "      Install the package podman to use an MinIO object storage as cache"
      echo "      This accelerates a rebuild of the CP2K dependencies, significantly"
      echo "      The folder ${SPACK_CACHE} will be used as cache"
      USE_CACHE="folder"
    fi
  fi

  # Initialize Spack shell hooks
  # shellcheck source=/dev/null
  source "${SPACK_ROOT}/share/spack/setup-env.sh"

  # Add spack cache as mirror if possible
  if [[ "${USE_CACHE}" == @("folder"|"minio") ]]; then
    echo "SPACK_CACHE = \"${SPACK_CACHE}\""
    spack mirror list
    if ! spack mirror list | grep -q "local-cache"; then
      echo "Setting up spack cache"
      if ((VERBOSE > 0)); then
        "${CP2K_ROOT}"/tools/docker/scripts/setup_spack_cache.sh
      else
        "${CP2K_ROOT}"/tools/docker/scripts/setup_spack_cache.sh &> /dev/null
      fi
    fi
  else
    echo "INFO: A spack cache is NOT used"
  fi

  # Prepare the CP2K spack configuration file
  CP2K_DEPS_FILE="${CP2K_ROOT}/tools/spack/cp2k_deps_${CP2K_VERSION:0:1}${CP2K_VERSION:4}.yaml"
  if [[ -f "${CP2K_DEPS_FILE}" ]]; then
    sed -E \
      -e "s|root: /opt/spack|root: ${SPACK_ROOT}/opt/spack|" \
      -e "/\"build_type=/s|build_type=[^\"]*|build_type=${DEPS_BUILD_TYPE}|" \
      "${CP2K_DEPS_FILE}" > "${CP2K_CONFIG_FILE}"
  else
    echo -e "\nERROR: The spack configuration file ${CP2K_DEPS_FILE} was not found\n"
    exit 1
  fi

  # Apply selected MPI type if needed
  # MPICH is selected by default for psmp and no change needed for ssmp
  if [[ "${MPI_MODE}" == "openmpi" ]]; then
    sed -E \
      -e '/\s*-\s+"mpich@/ s/^ /#/' \
      -e '/\s*#\s*-\s+"openmpi@/ s/#/ /' \
      -e '/\s*-\s+mpich/ s/mpich$/openmpi/' \
      -i "${CP2K_CONFIG_FILE}"
  fi

  # Activate CUDA in the spack configuration file if requested
  if ((CUDA_SM_CODE > 0)); then
    # The generated MPI stack is not CUDA-aware, so COSMA has to stage MPI
    # transfers through host memory instead of passing device pointers.
    sed -E \
      -e "0,/~cuda/s//+cuda cuda_arch=${CUDA_SM_CODE}/" \
      -e 's/"~cuda\s+~gpu_direct"/"\+cuda ~gpu_direct"/' \
      -e '/\s*#\s*-\s+"fabrics=efa,ucx"/ s/#/ /' \
      -e "/^[[:space:]]+libxc:/{n; n; s/- \"~cuda\"/- \"+cuda cuda_arch=${CUDA_SM_CODE}\"/}" \
      -i "${CP2K_CONFIG_FILE}"
    # Building libfabric with CUDA causes problems
    # sed -E -e 's/"~cuda\s+~gdrcopy"/"\+cuda \+gdrcopy"/' -i "${CP2K_CONFIG_FILE}"
    sed -E -e 's/"~cuda\s+~gdrcopy"/"\~cuda"/' -i "${CP2K_CONFIG_FILE}"
    echo -e "\nLibxc will be built with CUDA support (cuda_arch=${CUDA_SM_CODE})"
    if [[ -n "${CUDA_VERSION:-}" ]]; then
      # Set CUDA SM code
      sed -E -e "s/spec:\s+cuda@[.0-9]*/spec: cuda@${CUDA_VERSION}/" -i "${CP2K_CONFIG_FILE}"
    fi
    if [[ -n "${CUDA_HOME:-}" ]]; then
      sed -E -e "s|prefix: /usr/local/cuda|prefix: ${CUDA_HOME}|" -i "${CP2K_CONFIG_FILE}"
    fi
  else
    # CUDA is required for LibGint (make_cp2k.sh sets -DCP2K_USE_LIBGINT=OFF accordingly)
    sed -E \
      -e 's/"~cuda\s+~gdrcopy"/"\~cuda"/' \
      -e '/\s*-\s+"libgint@/ s/^ /#/' \
      -i "${CP2K_CONFIG_FILE}"
  fi

  # Activate OpenCL support if requested
  if [[ "${USE_OPENCL}" == "yes" ]]; then
    sed -E \
      -e 's/"~opencl"/"+opencl"/' \
      -e '/\s*#\s*-\s+"libxstream@/ s/#/ /' \
      -i "${CP2K_CONFIG_FILE}"
  fi

  # Apply Cray specific adaptation of the spack configuration if requested (CSCS)
  if [[ "${CRAY}" == "yes" ]]; then
    sed -E \
      -e '/\s*#\s*-\s+"netmod=ofi"/ s/#/ /' \
      -e 's/~xpmem/+xpmem/' \
      -e 's/"libfabric@[.0-9]*"/"libfabric@1.22.0"/' \
      -i "${CP2K_CONFIG_FILE}"
  fi

  # Apply feature selection to spack configuration file
  if [[ -n "${SED_PATTERN_LIST}" ]]; then
    eval sed -E "${SED_PATTERN_LIST}" -i "${CP2K_CONFIG_FILE}"
  fi

  # Find all compilers
  echo "Searching for GCC compilers"
  if ! spack compiler find &> /dev/null; then
    echo "ERROR: The compiler detection of spack failed"
    exit 1
  fi
  spack compiler list

  # Retrieve the newest compiler version found by spack
  GCC_VERSION_NEWEST="$(spack compilers | sed -nE 's/.*gcc@([0-9]+(\.[0-9]+)*).*/\1/p' | sort -V | tail -n 1)"
  echo "The newest GCC compiler version found by spack is ${GCC_VERSION_NEWEST}"
  GCC_VERSION_NEWEST="$(echo "${GCC_VERSION_NEWEST}" | sed -E -e 's/.*@([0-9]+).*/\1/' | cut -d. -f1)"

  # Check if the newest compiler version found on the host system is new enough
  GCC_VERSION_MINIMUM="10"
  if ((GCC_VERSION_NEWEST < GCC_VERSION_MINIMUM)); then
    echo "INFO: The newest GCC compiler version ${GCC_VERSION_NEWEST} is too old,"
    echo "      because at least version ${GCC_VERSION_MINIMUM} is required"
    GCC_VERSION="14"
    echo "INFO: The lastest stable GCC version ${GCC_VERSION} will be built and used by spack"
  fi

  # Add a specific GCC version to the spack configuration if requested
  if [[ "${GCC_VERSION}" == "auto" ]]; then
    echo "Spack will automatically select the GCC version which is not necessarily the newest one found"
  else
    sed -E -e "s/gcc@10:/gcc@${GCC_VERSION}/" -i "${CP2K_CONFIG_FILE}"
  fi

  # Create CP2K environment if needed
  if spack env list | grep -q "${CP2K_ENV}"; then
    echo "The Spack environment \"${CP2K_ENV}\" exists already"
  else
    cat "${CP2K_CONFIG_FILE}"
    echo "The Spack environment \"${CP2K_ENV}\" does NOT exist"
    if spack env create "${CP2K_ENV}" "${CP2K_CONFIG_FILE}" &> /dev/null; then
      echo "The Spack environment \"${CP2K_ENV}\" has been created"
    else
      echo "ERROR: The creation of the Spack environment \"${CP2K_ENV}\" failed"
      exit 1
    fi
  fi

  # Activate CP2K environment if needed
  if spack env activate "${CP2K_ENV}" &> /dev/null; then
    echo "The Spack environment \"${CP2K_ENV}\" has been activated"
  else
    echo "ERROR: The activation of the Spack environment \"${CP2K_ENV}\" failed"
    exit 1
  fi

  # Update the repo builtin
  if ! retry_download spack repo update; then
    echo "ERROR: The update of the repo builtin failed"
    exit 1
  fi

  # CUDA-enabled libxc needs a C++ compiler, which the builtin recipe does not declare yet
  if ((CUDA_SM_CODE > 0)); then
    LIBXC_PACKAGE_FILE="$(find -L "${SPACK_USER_CACHE_PATH}/package_repos" -path "*/builtin/packages/libxc/package.py" -print -quit)"
    if [[ -f "${LIBXC_PACKAGE_FILE}" ]] && ! grep -q "type=\"build\", when=\"+cuda\"" "${LIBXC_PACKAGE_FILE}"; then
      sed -i \
        -e 's/^\(    depends_on("c", type="build")\)$/\1\n    depends_on("cxx", type="build", when="+cuda")/' \
        "${LIBXC_PACKAGE_FILE}"
      echo "The builtin spack recipe of libxc has been patched to add the cxx build dependency for CUDA builds"
    fi
  fi

  # Add the local CP2K development Spack repository when missing
  export CP2K_REPO="cp2k_dev"
  if ! spack repo list | grep -q "${CP2K_REPO}"; then
    spack repo add --scope "env:${CP2K_ENV}" "${CP2K_ROOT}/tools/spack/spack_repo/${CP2K_REPO}"
  else
    echo "The CP2K development repo ${CP2K_REPO} is already present ... skipping add"
  fi

  spack -e "${CP2K_ENV}" repo list

  # Find all external packages
  if [[ "${USE_EXTERNALS}" == "yes" ]]; then
    if ! spack -C "${SPACK_USER_CONFIG_PATH}" external find; then
      echo -e "\nERROR: The detection of externals by spack failed"
      exit 1
    fi
    # Concretize CP2K dependencies
    if ! spack -e "${CP2K_ENV}" concretize --fresh --jobs "${NUM_PROCS}"; then
      echo -e "\nERROR: The spack concretize for environment \"${CP2K_ENV}\" failed"
      echo ""
      echo "HINT: The (-ue | --use_externals) flags can cause conflicts with outdated"
      echo "      packages on the host system, e.g. old python or gcc versions"
      echo ""
      exit 1
    fi
  else
    # If the installation of a specific GCC version for building CP2K is requested,
    # then an installation of the same GCC version on the host system has to be removed
    # from the configuration for a static build of CP2K
    if [[ "${CP2K_VERSION}" == *"-static" ]]; then
      if spack compilers | grep -q "gcc@${GCC_VERSION}"; then
        echo "INFO: A static build of CP2K with GCC ${GCC_VERSION} is requested without using"
        echo "      any externals which requires to build that GCC version with spack"
        echo "      while removing a gcc@${GCC_VERSION} installed on the host from the spack"
        echo "      configuration at the same time"
        if ! spack compiler remove "gcc@${GCC_VERSION}"; then
          echo -e "\nERROR: Removing the gcc@${GCC_VERSION} compiler failed"
          exit 1
        else
          echo -e "\nINFO: gcc@${GCC_VERSION} compiler removed"
        fi
        spack compiler list
      fi
    fi
    if ! spack -e "${CP2K_ENV}" concretize --fresh --jobs "${NUM_PROCS}"; then
      echo -e "\nERROR: The spack concretize for environment \"${CP2K_ENV}\" failed"
      exit 1
    fi
  fi

  ((VERBOSE > 0)) && spack find -c

  # Install CP2K dependencies via Spack (a retry skips already installed packages)
  if ! retry_download spack -e "${CP2K_ENV}" install -j "${NUM_PROCS}" -p "${NUM_PACKAGES}" "${VERBOSE_SPACK}"; then
    echo "ERROR: Building the CP2K dependencies with spack failed"
    if [[ "${USE_EXTERNALS}" == "yes" ]]; then
      echo "HINT:  Try to re-run the build without the (-ue | --use_externals) flag which avoids"
      echo "       errors or conflicts caused by externals from the host system"
    fi
    exit 1
  fi

  # Find all compilers
  if ! spack compiler find; then
    echo "ERROR: The compiler detection of spack failed"
    exit 1
  fi

  # Fix Libs list in elpa pkg-config file when CUDA is used
  if ((CUDA_SM_CODE > 0)); then
    ELPA_PKG_CONFIG_FILE="$(find -L "${SPACK_ROOT}/opt/spack/view" -name "elpa*.pc")"
    if [[ -f "${ELPA_PKG_CONFIG_FILE}" ]]; then
      sed -E -e 's|Libs: |Libs: -L/usr/local/cuda/lib64 |' -i "${ELPA_PKG_CONFIG_FILE}"
    fi
  fi

  # Return from spack folder after all installations are done
  cd "${CP2K_ROOT}" || exit 1

  # Make a note of the successful build
  touch "${SPACK_BUILD_PATH}/BUILD_DEPENDENCIES_COMPLETED"

  echo -e "\n*** Installation of CP2K dependencies completed ***\n"

fi

# PyTorch's TorchConfig.cmake is buried in the Python site-packages directory
# shellcheck disable=SC2034 # written to deps.conf
Torch_DIR="$(dirname "$(find "${SPACK_ROOT}" ! -type l -name TorchConfig.cmake | tail -n 1)")"

# Write the settings provided by the dependency build
write_config Torch_DIR > "${CONFIG_PATH}/deps.conf"

### End of build CP2K dependencies ###
