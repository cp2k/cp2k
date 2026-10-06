#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.

# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

# libtorch is built from the PyTorch sources with only the features needed by
# CP2K/FTorch/SKALA. 2.6.0 is the last PyTorch release officially paired with
# CUDA 12.4, which matches the system CUDA toolkit used by this toolchain.
libtorch_ver="2.6.0"
libtorch_rev="v${libtorch_ver}"

# shellcheck source=/dev/null
source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_libtorch" ] && rm "${BUILDDIR}/setup_libtorch"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

case "${with_libtorch}" in
  __INSTALL__)
    echo "==================== Building libtorch from source ===================="
    pkg_install_dir="${INSTALLDIR}/libtorch-${libtorch_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"

    if [ -z "${ARCH_NUM:-}" ] || [ "${ARCH_NUM}" = "no" ]; then
      report_error ${LINENO} "ARCH_NUM is not set; building libtorch needs a target GPU architecture (see --gpu-ver)."
    fi
    # TORCH_CUDA_ARCH_LIST wants a decimal form (e.g. 8.6, 8.0, 9.0).
    if [[ "${ARCH_NUM}" == *.* ]]; then
      torch_cuda_arch="${ARCH_NUM}"
    else
      torch_cuda_arch="${ARCH_NUM:0:1}.${ARCH_NUM:1}"
    fi
    torch_build_jobs="${NPROCS_OVERWRITE:-$(get_nprocs)}"

    if verify_checksums "${install_lock_file}"; then
      echo "libtorch-${libtorch_ver} is already installed, skipping it."
    else
      # Python build dependencies required by the PyTorch build system are
      # installed into a dedicated virtual environment to avoid clashing with
      # the (possibly externally-managed) system Python.
      torch_venv="${BUILDDIR}/pytorch-venv"
      if [ ! -x "${torch_venv}/bin/python3" ]; then
        python3 -m venv "${torch_venv}" > /dev/null 2>&1 ||
          report_error ${LINENO} "Failed to create a Python virtual environment at ${torch_venv}."
      fi
      "${torch_venv}/bin/python3" -m pip install --quiet --disable-pip-version-check \
        packaging pyyaml setuptools typing_extensions filelock sympy jinja2 \
        networkx > /dev/null 2>&1 ||
        report_error ${LINENO} "Failed to install the Python build dependencies of PyTorch."
      torch_python="${torch_venv}/bin/python3"

      # Fetch PyTorch sources with submodules. A full (non-shallow) submodule
      # checkout is used because some submodules are nested and would otherwise
      # be left empty, breaking the CMake configure step.
      src_dir="${BUILDDIR}/pytorch-${libtorch_ver}"
      if [ ! -d "${src_dir}/.git" ]; then
        rm -rf "${src_dir}"
        git clone --depth 1 --branch "${libtorch_rev}" \
          https://github.com/pytorch/pytorch.git "${src_dir}" > git.log 2>&1 ||
          report_error ${LINENO} "Failed to clone PyTorch ${libtorch_rev}; see ${BUILDDIR}/git.log"
      fi
      git -C "${src_dir}" submodule sync --recursive >> git.log 2>&1 || true
      git -C "${src_dir}" submodule update --init --recursive --force >> git.log 2>&1 ||
        report_error ${LINENO} "Failed to update PyTorch submodules; see ${BUILDDIR}/git.log"

      # Reduce the build to the features CP2K/FTorch/SKALA actually use:
      # no Python bindings, no distributed/NCCL/GLOO/MPI, no quantized kernels,
      # no MKLDNN/NNPACK, no tests, no cuDNN (the model uses no convolution),
      # single CUDA architecture. The operator set is kept complete so the JIT
      # can load arbitrary TorchScript models without "unknown builtin op".
      build_dir="${BUILDDIR}/pytorch-build-${libtorch_ver}"
      rm -rf "${build_dir}"
      mkdir -p "${build_dir}"
      cd "${build_dir}"

      echo "Building libtorch ${libtorch_ver} from source (CUDA arch ${torch_cuda_arch})"
      # Configure with CMake directly (BUILD_PYTHON=OFF), then build the
      # 'install' target. This installs into CMAKE_INSTALL_PREFIX so the
      # headers, libraries and CMake package files land next to the other
      # toolchain packages, which is what CP2K and FTorch consume.
      # Put the build venv first on PATH so CMake finds its Python, and pass
      # the interpreter explicitly. Environment variables are NOT forwarded to
      # CMake by the 2.6.0 build system, so every toggle is a -D option.
      export PATH="${torch_venv}/bin:${PATH}"
      export MAX_JOBS="${torch_build_jobs}"
      export TORCH_CUDA_ARCH_LIST="${torch_cuda_arch}"
      export CFLAGS="" CXXFLAGS=""
      cmake -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DPYTHON_EXECUTABLE="${torch_python}" \
        -DPython_EXECUTABLE="${torch_python}" \
        -DBUILD_PYTHON=OFF \
        -DBUILD_TEST=OFF \
        -DUSE_CUDA=ON \
        -DUSE_CUDNN=OFF \
        -DUSE_CUSPARSELT=OFF \
        -DUSE_CUDSS=OFF \
        -DUSE_CUFILE=OFF \
        -DUSE_DISTRIBUTED=OFF \
        -DUSE_NCCL=OFF \
        -DUSE_GLOO=OFF \
        -DUSE_MPI=OFF \
        -DUSE_TENSORPIPE=OFF \
        -DUSE_MKLDNN=OFF \
        -DUSE_NNPACK=OFF \
        -DUSE_QNNPACK=OFF \
        -DUSE_FBGEMM=OFF \
        -DUSE_KINETO=OFF \
        -DUSE_NUMA=OFF \
        -DUSE_ITT=OFF \
        -DUSE_OPENMP=ON \
        -DUSE_CUDA_STATIC_LINK=OFF \
        -DREL_WITH_DEB_INFO=OFF \
        "${src_dir}" > configure.log 2>&1 || {
        tail_excerpt configure.log
        report_error ${LINENO} "Configuring libtorch from source failed; see ${build_dir}/configure.log"
      }
      cmake --build . --target install -j "${torch_build_jobs}" > build.log 2>&1 || {
        tail_excerpt build.log
        report_error ${LINENO} "Building libtorch from source failed; see ${build_dir}/build.log"
      }

      if [ ! -f "${pkg_install_dir}/share/cmake/Torch/TorchConfig.cmake" ]; then
        report_error ${LINENO} "From-source libtorch install is missing ${pkg_install_dir}/share/cmake/Torch/TorchConfig.cmake"
      fi

      write_checksums "${install_lock_file}" "${SCRIPT_DIR}/stage6/$(basename "${SCRIPT_NAME}")"

      # Remove the PyTorch sources and build tree now that libtorch is
      # installed; they are large and no longer needed.
      echo "Deleting PyTorch build directory and sources ..."
      cd "${BUILDDIR}"
      rm -rf "${build_dir}" "${src_dir}"
      echo "PyTorch build directory and sources deleted."
    fi
    ;;

  __SYSTEM__)
    echo "==================== Finding libtorch from system paths ===================="
    check_lib -ltorch "libtorch"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libtorch.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking libtorch to user paths ===================="
    pkg_install_dir="${with_libtorch}"
    # use the lib64 directory if present (multi-abi distros may link lib/ to lib32/ instead)
    LIBTORCH_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && LIBTORCH_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${LIBTORCH_LIBDIR}"
    ;;
esac

if [ "$with_libtorch" != "__DONTUSE__" ]; then
  cat << EOF > "${BUILDDIR}/setup_libtorch"
export LIBTORCH_VER="${libtorch_ver}"
EOF
  if [ "$with_libtorch" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_libtorch"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path PKG_CONFIG_PATH "${pkg_install_dir}/lib/pkgconfig"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_libtorch" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_libtorch"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "libtorch"
