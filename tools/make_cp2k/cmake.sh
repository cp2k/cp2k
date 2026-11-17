# shellcheck shell=bash disable=SC2034,SC2153,SC2154
#
# make_cp2k.sh module: CMake configuration, build, and installation of CP2K
#
# This file is sourced by make_cp2k.sh at top level and shares all variables
# with the other modules. It is not meant to be executed on its own.

### Build CP2K ###

if spack env list | grep -q "${CP2K_ENV}"; then
  spack env list
else
  echo "ERROR: No Spack environment \"${CP2K_ENV}\" found"
  exit 1
fi

# Activate spack environment
eval "$(spack env activate --sh "${CP2K_ENV}")"

spack env status

# CMake configuration step
export CMAKE_BUILD_PATH="${BUILD_PATH}/build"

# PyTorch's TorchConfig.cmake is buried in the Python site-packages directory
Torch_DIR="$(dirname "$(find "${SPACK_ROOT}" ! -type l -name TorchConfig.cmake | tail -n 1)")"
export Torch_DIR

if [[ ! -d "${CMAKE_BUILD_PATH}" ]]; then
  mkdir -p "${CMAKE_BUILD_PATH}"
  case "${CP2K_VERSION}" in
    pdbg | psmp)
      # shellcheck disable=SC2086
      cmake -S "${CP2K_ROOT}" -B "${CMAKE_BUILD_PATH}" --preset "${CMAKE_PRESET}" \
        -GNinja \
        -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} \
        -DCMAKE_BUILD_TYPE="${CP2K_BUILD_TYPE}" \
        -DCMAKE_Fortran_COMPILER_LAUNCHER="${Fortran_COMPILER_LAUNCHER}" \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_INSTALL_MESSAGE="${INSTALL_MESSAGE}" \
        -DCMAKE_SKIP_RPATH="ON" \
        -DCMAKE_VERBOSE_MAKEFILE="${VERBOSE_MAKEFILE}" \
        ${CMAKE_FEATURE_FLAGS} \
        ${CMAKE_CUDA_FLAGS} \
        -Wno-error=dev |&
        tee "${CMAKE_BUILD_PATH}/cmake.log"
      EXIT_CODE=$?
      ;;
    sdbg | ssmp)
      # shellcheck disable=SC2086
      cmake -S "${CP2K_ROOT}" -B "${CMAKE_BUILD_PATH}" --preset "${CMAKE_PRESET}" \
        -GNinja \
        -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} \
        -DCMAKE_BUILD_TYPE="${CP2K_BUILD_TYPE}" \
        -DCMAKE_Fortran_COMPILER_LAUNCHER="${Fortran_COMPILER_LAUNCHER}" \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_INSTALL_MESSAGE="${INSTALL_MESSAGE}" \
        -DCMAKE_SKIP_RPATH="ON" \
        -DCMAKE_VERBOSE_MAKEFILE="${VERBOSE_MAKEFILE}" \
        ${CMAKE_FEATURE_FLAGS} \
        ${CMAKE_CUDA_FLAGS} \
        -Wno-error=dev |&
        tee "${CMAKE_BUILD_PATH}/cmake.log"
      EXIT_CODE=$?
      ;;
    ssmp-static)
      # Find some static libraries in advance
      LIBOPENBLAS=$(find -L "${SPACK_ROOT}"/opt/spack/view -name libopenblas.a)
      OPENBLAS_INCLUDE_DIR="$(
        dirname "$(find -L "${SPACK_ROOT}"/opt/spack/view -name cblas.h -print -quit)"
      )"
      LIBM="$(find /usr -name libm.a 2> /dev/null)"
      # shellcheck disable=SC2086
      cmake -S "${CP2K_ROOT}" -B "${CMAKE_BUILD_PATH}" --preset "${CMAKE_PRESET}" \
        -GNinja \
        -DBUILD_SHARED_LIBS="OFF" \
        -DCMAKE_BUILD_TYPE="${CP2K_BUILD_TYPE}" \
        -DCMAKE_EXE_LINKER_FLAGS="-static" \
        -DCMAKE_FIND_LIBRARY_SUFFIXES=".a" \
        -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" \
        -DCMAKE_INSTALL_LIBDIR="lib" \
        -DCMAKE_INSTALL_MESSAGE="${INSTALL_MESSAGE}" \
        -DCMAKE_SKIP_RPATH="ON" \
        -DCMAKE_VERBOSE_MAKEFILE="${VERBOSE_MAKEFILE}" \
        -DCP2K_BLAS_VENDOR="OpenBLAS" \
        -DCP2K_BLAS_INCLUDE_DIRS="${OPENBLAS_INCLUDE_DIR}" \
        -DCP2K_BLAS_LINK_LIBRARIES="${LIBOPENBLAS};${LIBM}" \
        -DCP2K_LAPACK_LINK_LIBRARIES="${LIBOPENBLAS};${LIBM}" \
        ${CMAKE_FEATURE_FLAGS} \
        -Wno-error=dev |&
        tee "${CMAKE_BUILD_PATH}/cmake.log"
      EXIT_CODE=$?
      # It is almost impossible to avoid that shared libraries are pulled in
      # when the CP2K dependencies are built and thus just change all shared
      # to static libraries hoping that these are also available
      sed -E -e "s/\.so/.a/g" -i build/build.ninja
      ;;
    *)
      echo "ERROR: Unknown CP2K version ${CP2K_VERSION} specified"
      exit 1
      ;;
  esac
  if ((EXIT_CODE != 0)); then
    echo "ERROR: The CMake configuration step failed with the error code ${EXIT_CODE}"
    echo "       You can try to remove the build folder with 'rm -rf build' and re-run"
    echo "       or even start the whole CP2K installation from scratch with the \"-bd\" flag"
    exit "${EXIT_CODE}"
  fi
fi

# CMake build step
echo -e "\n*** Compiling CP2K ***\n"
cmake --build "${CMAKE_BUILD_PATH}" --parallel "${NUM_PROCS}" -- "${VERBOSE_FLAG}" |& tee "${CMAKE_BUILD_PATH}"/ninja.log
EXIT_CODE=${PIPESTATUS[0]}
if ((EXIT_CODE != 0)); then
  echo "ERROR: The CMake build step failed with the error code ${EXIT_CODE}"
  exit "${EXIT_CODE}"
fi

# CMake install step
echo -e "\n*** Installing CP2K ***\n"
cmake --install "${CMAKE_BUILD_PATH}" |& tee "${CMAKE_BUILD_PATH}"/install.log
EXIT_CODE=${PIPESTATUS[0]}
if ((EXIT_CODE != 0)); then
  echo -e "\nERROR: The CMake installation step failed with the error code ${EXIT_CODE}"
  exit "${EXIT_CODE}"
fi

# Install Skala resource files if needed
export SKALA_MODEL="(not available)"
if spack location -i gauxc &> /dev/null; then
  GAUXC_PATH="share/gauxc/onedft_models"
  GAUXC_PREFIX="$(spack location -i gauxc)"
  GAUXC_MODEL_SOURCE_PATH="${GAUXC_PREFIX%/}/${GAUXC_PATH}"
  if [[ -d "${GAUXC_MODEL_SOURCE_PATH}" ]]; then
    GAUXC_MODEL_TARGET_PATH="${INSTALL_PREFIX%/}/${GAUXC_PATH}"
    mkdir -p "${GAUXC_MODEL_TARGET_PATH}"
    if compgen -G "${GAUXC_MODEL_SOURCE_PATH}/*.fun" > /dev/null; then
      cp "${GAUXC_MODEL_SOURCE_PATH}"/*.fun "${GAUXC_MODEL_TARGET_PATH}/"
    else
      echo -e "\nERROR: No GauXC model files found in source folder ${GAUXC_MODEL_SOURCE_PATH}"
      exit 1
    fi
    shopt -s nullglob
    MATCHES=("${GAUXC_MODEL_TARGET_PATH}"/skala*)
    shopt -u nullglob
    if ((${#MATCHES[@]} > 0)); then
      export SKALA_MODEL="${MATCHES[${#MATCHES[@]} - 1]}"
    else
      echo -e "\nERROR: Failed to resolve SKALA_MODEL in target path ${GAUXC_MODEL_TARGET_PATH}"
      exit 1
    fi
  fi
fi
echo -e "\nSKALA_MODEL = ${SKALA_MODEL}"

# Collect and compress all log files when building within a container
if [[ "${IN_CONTAINER}" == "yes" ]]; then
  if ! cat "${CMAKE_BUILD_PATH}"/cmake.log \
    "${CMAKE_BUILD_PATH}"/ninja.log \
    "${CMAKE_BUILD_PATH}"/install.log |
    gzip > "${CP2K_ROOT}"/build_cp2k.log.gz; then
    echo -e "\nERROR: The compressed log file generation failed"
    exit 1
  fi
fi
