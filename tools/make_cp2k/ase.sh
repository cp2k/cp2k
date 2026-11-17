# shellcheck shell=bash disable=SC2034,SC2153,SC2154
#
# make_cp2k.sh module: Optional ASE build and test
#
# This file is sourced by make_cp2k.sh at top level and shares all variables
# with the other modules. It is not meant to be executed on its own.

# Optionally, build CP2K with ASE support
if [[ -n "${ASE_VERSION}" ]]; then

  # Download ASE
  ASE_ROOT="${CMAKE_BUILD_PATH}"/ase
  [[ -d "${ASE_ROOT}" ]] && rm -rf "${ASE_ROOT}"
  echo -e "\n*** Downloading ASE ${ASE_VERSION} ***\n"
  # Spack's OpenSSL libs on LD_LIBRARY_PATH can break the system git-remote-https
  # helper (mismatched libldap), so clone with a clean library path
  if ! LD_LIBRARY_PATH="" retry_download git clone -b "${ASE_VERSION}" -c advice.detachedHead=false --depth=1 "${VERBOSE_FLAG}" --single-branch \
    https://gitlab.com/ase/ase.git "${ASE_ROOT}"; then
    echo -e "\nERROR: Downloading ASE ${ASE_VERSION} failed"
    exit 1
  fi
  cd "${ASE_ROOT}" || exit 1
  ASE_REVISION="$(git rev-parse --short HEAD)"

  # Install ASE
  echo -e "\n*** Installing ASE ${ASE_VERSION} ***\n"
  if ! python3 -m venv "${INSTALL_PREFIX}"/ase; then
    echo -e "\nERROR: The creation of the virtual environment for ASE failed"
    exit 1
  fi
  export PATH="${INSTALL_PREFIX}/ase/bin:${PATH}"
  if ! retry_download "${INSTALL_PREFIX}"/ase/bin/python3 -m pip install --ignore-installed "${VERBOSE_FLAG}" ".[test]"; then
    echo -e "\nERROR: The ASE installation (venv) failed"
    exit 1
  fi
  cat << *** > "${INSTALL_PREFIX}"/ase/config.ini
[cp2k]
cp2k_shell = ${INSTALL_PREFIX}/bin/cp2k_shell
cp2k_main =  ${INSTALL_PREFIX}/bin/cp2k
***
  # Install additional packages for ASE
  if ! retry_download "${INSTALL_PREFIX}"/ase/bin/python3 -m pip install --ignore-installed "${VERBOSE_FLAG}" matplotlib numpy packaging six spglib; then
    echo -e "\nERROR: The installation of additional packages for ASE failed"
    exit 1
  fi
  echo ""
  echo "*** The ASE/CP2K installation can be tested with"
  if [[ "${IN_CONTAINER}" == "yes" ]]; then
    echo "    podman run -it --rm ${IMAGE_TAG} ${LAUNCH_SCRIPT} test_ase"
  else
    echo "    ${LAUNCH_SCRIPT} test_ase"
  fi
  echo ""

  # Test the ASE installation
  echo "${LAUNCH_SCRIPT} ase test -j 0 -c cp2k calculator/cp2k" > "${INSTALL_PREFIX}"/ase/bin/test_ase
  chmod 750 "${INSTALL_PREFIX}"/ase/bin/test_ase
  if [[ "${TEST_ASE}" == "yes" ]]; then
    echo -e "\n*** Running ASE tests ***\n"
    if [[ "${IN_CONTAINER}" == "yes" ]]; then
      export PYTEST_DEBUG_TEMPROOT="/workspace/artifacts"
      mkdir -p ${PYTEST_DEBUG_TEMPROOT}
    fi
    if "${INSTALL_PREFIX}"/ase/bin/test_ase; then
      echo -e "\nSummary: ASE commit ${ASE_REVISION} works fine"
      echo -e "Status: OK\n"
      exit 0
    else
      echo -e "\nSummary: Something is wrong with ASE commit ${ASE_REVISION}"
      echo -e "Status: FAILED\n"
      exit 0
    fi
  fi

fi
