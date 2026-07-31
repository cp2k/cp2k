#!/bin/bash -e

# TODO: Review and if possible fix shellcheck errors.
# shellcheck disable=all

[ "${BASH_SOURCE[0]}" ] && SCRIPT_NAME="${BASH_SOURCE[0]}" || SCRIPT_NAME=$0
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_NAME")/.." && pwd -P)"

ftorch_ver="1.0.0"
ftorch_sha256="c4b6741e582623b7ecaecd59d02f779e8a6f6017f8068c85da8a034f468df375"
ftorch_pkg="FTorch-${ftorch_ver}.tar.gz"
ftorch_urlpath="https://github.com/Cambridge-ICCS/FTorch/archive/refs/tags"

skala_model_ver="1.1"
skala_model_pkg="skala-${skala_model_ver}.fun"
skala_model_sha256="0c8432ac3f03c8f1276372df9aca5b7ee7f8939d47a8789eb158976e89aa0606"
skala_model_urlpath="https://huggingface.co/microsoft/skala-${skala_model_ver}/resolve/main"

source "${SCRIPT_DIR}"/common_vars.sh
source "${SCRIPT_DIR}"/tool_kit.sh
source "${SCRIPT_DIR}"/signal_trap.sh
source "${INSTALLDIR}"/toolchain.conf
source "${INSTALLDIR}"/toolchain.env

[ -f "${BUILDDIR}/setup_ftorch" ] && rm "${BUILDDIR}/setup_ftorch"

! [ -d "${BUILDDIR}" ] && mkdir -p "${BUILDDIR}"
cd "${BUILDDIR}"

# Skala model is installed by shared script in main toolchain

retrieve_github_archive() {
  local __sha256="$1"
  local __filename="$2"
  local __urlpath="$3"
  local __outfile="$4"
  local __attempt
  if ! [ -f "${__outfile}" ]; then
    for __attempt in 1 2 3 4 5; do
      download_pkg_from_urlpath "${__sha256}" "${__filename}" "${__urlpath}" "${__outfile}" && return
      echo "Download attempt ${__attempt} for ${__filename} failed."
      sleep $((__attempt * 5))
    done
    return 1
  elif ! checksum "${__sha256}" "${__outfile}"; then
    echo "${__outfile} is found but checksum is wrong; delete and re-download"
    rm -vf "${__outfile}"
    for __attempt in 1 2 3 4 5; do
      download_pkg_from_urlpath "${__sha256}" "${__filename}" "${__urlpath}" "${__outfile}" && return
      echo "Download attempt ${__attempt} for ${__filename} failed."
      sleep $((__attempt * 5))
    done
    return 1
  else
    echo "${__outfile} is found and checksum is right"
  fi
}

cat > skala_ftorch.cpp << 'SKALA_CXX_EOF'
#include <torch/torch.h>
#include <torch/script.h>
#include <vector>
#include <string>
#include <cstdint>

namespace {

const char* FEATURE_KEYS[] = {
    "density", "grad", "kin", "grid_coords", "grid_weights",
    "coarse_0_atomic_coords", "atomic_grid_weights",
    "atomic_grid_sizes", "atomic_grid_size_bound_shape"
};
const int NUM_FEATURES = 9;

torch::jit::script::Module load_skala_model(
    const std::string& model_path,
    const std::vector<std::string>& serialized_features) {
  torch::jit::script::Module module;
  try {
    torch::jit::ExtraFilesMap extras;
    for (int i = 0; i < NUM_FEATURES; ++i) {
      extras[FEATURE_KEYS[i]] = serialized_features[i];
    }
    module = torch::jit::load(model_path, torch::kCPU, false, extras);
    module.eval();
  } catch (const c10::Error& e) {
    std::cerr << "[skala_ftorch ERROR] loading model: " << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }
  return module;
}

} // anonymous namespace

extern "C" {

void* skala_load_model(const char* model_path,
                       const char* const* serialized_features,
                       int n_features) {
  std::vector<std::string> feats;
  feats.reserve(n_features);
  for (int i = 0; i < n_features; ++i) {
    feats.emplace_back(serialized_features[i]);
  }
  torch::jit::script::Module* mod =
      new torch::jit::script::Module(load_skala_model(model_path, feats));
  return static_cast<void*>(mod);
}

void skala_unload_model(void* model_ptr) {
  delete static_cast<torch::jit::script::Module*>(model_ptr);
}

void skala_eval_exc(
    void* model_ptr,
    const double* density_ptr, const int64_t* density_shape, int density_ndim,
    const double* grid_weights_ptr, const int64_t* grid_weights_shape, int grid_weights_ndim,
    const double* grid_coords_ptr, const int64_t* grid_coords_shape, int grid_coords_ndim,
    const double* atomic_coords_ptr, const int64_t* atomic_coords_shape, int atomic_coords_ndim,
    const double* atomic_grid_weights_ptr, const int64_t* agw_shape, int agw_ndim,
    const int* atomic_grid_sizes_ptr, const int* atomic_grid_size_bound_ptr, int n_grid_types,
    double* out_exc, double* out_vxc) {

  auto* mod = static_cast<torch::jit::script::Module*>(model_ptr);

  auto dtype_opts = torch::TensorOptions().dtype(torch::kFloat64).layout(torch::kStrided);
  torch::Device device(torch::kCPU);

  std::vector<int64_t> density_shape_vec(density_shape, density_shape + density_ndim);
  std::vector<int64_t> grid_weights_shape_vec(grid_weights_shape, grid_weights_shape + grid_weights_ndim);
  std::vector<int64_t> grid_coords_shape_vec(grid_coords_shape, grid_coords_shape + grid_coords_ndim);
  std::vector<int64_t> atomic_coords_shape_vec(atomic_coords_shape, atomic_coords_shape + atomic_coords_ndim);

  torch::Tensor density = torch::from_blob(const_cast<double*>(density_ptr), density_shape_vec, dtype_opts).to(device);
  torch::Tensor grid_weights = torch::from_blob(const_cast<double*>(grid_weights_ptr), grid_weights_shape_vec, dtype_opts).to(device);
  torch::Tensor grid_coords = torch::from_blob(const_cast<double*>(grid_coords_ptr), grid_coords_shape_vec, dtype_opts).to(device);
  torch::Tensor atomic_coords = torch::from_blob(const_cast<double*>(atomic_coords_ptr), atomic_coords_shape_vec, dtype_opts).to(device);
  torch::Tensor atomic_grid_weights = torch::from_blob(const_cast<double*>(atomic_grid_weights_ptr), torch::ArrayRef<int64_t>(agw_shape, agw_ndim), dtype_opts).to(device);

  torch::Tensor atomic_grid_sizes = torch::from_blob(
      const_cast<int*>(atomic_grid_sizes_ptr), {n_grid_types}, torch::kInt32).to(device);
  torch::Tensor atomic_grid_size_bound = torch::from_blob(
      const_cast<int*>(atomic_grid_size_bound_ptr), {n_grid_types}, torch::kInt32).to(device);

  c10::Dict<torch::Tensor, torch::Tensor> dict;
  dict.insert(torch::tensor(1, torch::kInt32), density);
  dict.insert(torch::tensor(2, torch::kInt32), grid_weights);
  dict.insert(torch::tensor(3, torch::kInt32), grid_coords);
  dict.insert(torch::tensor(4, torch::kInt32), atomic_coords);
  dict.insert(torch::tensor(5, torch::kInt32), atomic_grid_weights);
  dict.insert(torch::tensor(6, torch::kInt32), atomic_grid_sizes);
  dict.insert(torch::tensor(7, torch::kInt32), atomic_grid_size_bound);

  c10::IValue input_val(dict);
  torch::NoGradGuard no_grad;

  c10::IValue output = mod->forward({input_val});

  torch::Tensor exc;
  torch::Tensor vxc;

  if (output.isDict()) {
    auto out_dict = output.toDict();
    if (out_dict.contains(torch::tensor(0, torch::kInt32))) {
      exc = out_dict.at(torch::tensor(0, torch::kInt32)).toTensor();
    }
    if (out_dict.contains(torch::tensor(1, torch::kInt32))) {
      vxc = out_dict.at(torch::tensor(1, torch::kInt32)).toTensor();
    }
  }

  if (exc.defined()) {
    torch::Tensor E_xc = (exc * grid_weights).sum();
    *out_exc = E_xc.item<double>();
  } else {
    *out_exc = 0.0;
  }

  if (vxc.defined()) {
    torch::Tensor weighted_vxc = vxc * grid_weights;
    weighted_vxc.copy_(weighted_vxc.sum(torch::kNone));
    *out_vxc = weighted_vxc.item<double>();
  } else {
    *out_vxc = 0.0;
  }
}

} // extern "C"
SKALA_CXX_EOF

cat > skala_ftorch.f90 << 'SKALA_F90_EOF'
module skala_ftorch
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int64_t, c_int, c_char, c_associated, c_null_ptr, c_loc
  implicit none

  type :: skala_model
    type(c_ptr) :: p = c_null_ptr
  end type skala_model

  interface
    function c_skala_load_model(model_path, serialized_features, n_features) &
        bind(c, name='skala_load_model') result(m)
      use, intrinsic :: iso_c_binding, only: c_char, c_int, c_ptr
      character(c_char), intent(in) :: model_path(*)
      type(c_ptr), intent(in) :: serialized_features(*)
      integer(c_int), value, intent(in) :: n_features
      type(c_ptr) :: m
    end function c_skala_load_model

    subroutine c_skala_unload_model(model_ptr) bind(c, name='skala_unload_model')
      use, intrinsic :: iso_c_binding, only: c_ptr
      type(c_ptr), value, intent(in) :: model_ptr
    end subroutine c_skala_unload_model

    subroutine c_skala_eval_exc(model_ptr, &
        density, density_shape, density_ndim, &
        grid_weights, grid_weights_shape, grid_weights_ndim, &
        grid_coords, grid_coords_shape, grid_coords_ndim, &
        atomic_coords, atomic_coords_shape, atomic_coords_ndim, &
        atomic_grid_weights, agw_shape, agw_ndim, &
        atomic_grid_sizes, atomic_grid_size_bound, n_grid_types, &
        out_exc, out_vxc) bind(c, name='skala_eval_exc')
      use, intrinsic :: iso_c_binding, only: c_double, c_int64_t, c_int, c_ptr
      type(c_ptr), value, intent(in) :: model_ptr
      type(c_ptr), value, intent(in) :: density
      integer(c_int64_t), intent(in) :: density_shape(*)
      integer(c_int), value, intent(in) :: density_ndim
      type(c_ptr), value, intent(in) :: grid_weights
      integer(c_int64_t), intent(in) :: grid_weights_shape(*)
      integer(c_int), value, intent(in) :: grid_weights_ndim
      type(c_ptr), value, intent(in) :: grid_coords
      integer(c_int64_t), intent(in) :: grid_coords_shape(*)
      integer(c_int), value, intent(in) :: grid_coords_ndim
      type(c_ptr), value, intent(in) :: atomic_coords
      integer(c_int64_t), intent(in) :: atomic_coords_shape(*)
      integer(c_int), value, intent(in) :: atomic_coords_ndim
      type(c_ptr), value, intent(in) :: atomic_grid_weights
      integer(c_int64_t), intent(in) :: agw_shape(*)
      integer(c_int), value, intent(in) :: agw_ndim
      integer(c_int), intent(in) :: atomic_grid_sizes(*)
      integer(c_int), intent(in) :: atomic_grid_size_bound(*)
      integer(c_int), value, intent(in) :: n_grid_types
      real(c_double), intent(out) :: out_exc
      real(c_double), intent(out) :: out_vxc
    end subroutine c_skala_eval_exc
  end interface

  private
  public :: skala_model, &
            skala_load, skala_unload, &
            skala_eval_exc

contains

  subroutine skala_load(path, n_features, model)
    character(len=*), intent(in) :: path
    integer, intent(in) :: n_features
    type(skala_model), intent(out) :: model
    model%p = c_skala_load_model(trim(path)//c_null_char, &
                                 c_null_ptr, int(n_features, c_int))
  end subroutine skala_load

  subroutine skala_unload(model)
    type(skala_model), intent(inout) :: model
    if (c_associated(model%p)) then
      call c_skala_unload_model(model%p)
      model%p = c_null_ptr
    end if
  end subroutine skala_unload

  subroutine skala_eval_exc(model, &
      density, grid_weights, grid_coords, &
      atomic_coords, atomic_grid_weights, &
      atomic_grid_sizes, atomic_grid_size_bound, &
      E_xc, V_xc)
    type(skala_model), intent(in) :: model
    real(8), target, intent(in) :: density(*)
    real(8), target, intent(in) :: grid_weights(*)
    real(8), target, intent(in) :: grid_coords(*)
    real(8), target, intent(in) :: atomic_coords(*)
    real(8), target, intent(in) :: atomic_grid_weights(*)
    integer, intent(in) :: atomic_grid_sizes(*)
    integer, intent(in) :: atomic_grid_size_bound(*)
    real(8), intent(out) :: E_xc
    real(8), intent(out) :: V_xc

    integer(c_int64_t) :: density_shape(4), gw_shape(1), gc_shape(2)
    integer(c_int64_t) :: ac_shape(2), agw_shape(1)
    integer(c_int) :: density_ndim, gw_ndim, gc_ndim, ac_ndim, agw_ndim
    integer(c_int) :: n_grid_types
    real(8) :: out_exc, out_vxc

    density_ndim = 4
    density_shape = [1_8, 1_8, 1_8, 1_8]
    gw_ndim = 1
    gw_shape = [1_8]
    gc_ndim = 2
    gc_shape = [1_8, 3_8]
    ac_ndim = 2
    ac_shape = [1_8, 3_8]
    agw_ndim = 1
    agw_shape = [1_8]
    n_grid_types = size(atomic_grid_sizes)

    call c_skala_eval_exc(model%p, &
        c_loc(density), density_shape, density_ndim, &
        c_loc(grid_weights), gw_shape, gw_ndim, &
        c_loc(grid_coords), gc_shape, gc_ndim, &
        c_loc(atomic_coords), ac_shape, ac_ndim, &
        c_loc(atomic_grid_weights), agw_shape, agw_ndim, &
        atomic_grid_sizes, atomic_grid_size_bound, n_grid_types, &
        out_exc, out_vxc)
    E_xc = out_exc
    V_xc = out_vxc
  end subroutine skala_eval_exc
end module skala_ftorch
SKALA_F90_EOF

case "${with_skala_ftorch}" in
  __INSTALL__)
    echo "==================== Installing FTorch + skala bindings ===================="
    pkg_install_dir="${INSTALLDIR}/ftorch-${ftorch_ver}"
    install_lock_file="${pkg_install_dir}/install_successful"

    if verify_checksums "${install_lock_file}"; then
      echo "ftorch-${ftorch_ver} is already installed, skipping it."
    else
      retrieve_github_archive "${ftorch_sha256}" "v${ftorch_ver}.tar.gz" \
        "${ftorch_urlpath}" "${ftorch_pkg}"
      # Skala model is now installed by shared install_skala.sh script

      echo "Installing from scratch into ${pkg_install_dir}"
      rm -rf "FTorch-${ftorch_ver}" "${pkg_install_dir}"
      tar -xzf "${ftorch_pkg}"

      mkdir -p build_ftorch
      cd build_ftorch

      cmake "${BUILDDIR}/FTorch-${ftorch_ver}" \
        -DCMAKE_INSTALL_PREFIX="${pkg_install_dir}" \
        -DCMAKE_INSTALL_LIBDIR=lib \
        -DCMAKE_C_COMPILER="${CC}" \
        -DCMAKE_CXX_COMPILER="${CXX}" \
        -DCMAKE_Fortran_COMPILER="${FC}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_SHARED_LIBS=ON \
        > configure.log 2>&1 || tail_excerpt configure.log
      make -j $(get_nprocs) > make.log 2>&1 || tail_excerpt make.log
      make install > install.log 2>&1 || tail_excerpt install.log

      cd "${BUILDDIR}"

       export SKALA_BINDINGS_SRC_DIR="${BUILDDIR}"

       echo "Building scala C++ bindings..."
       if ! "${CXX}" -shared -fPIC -std=c++17 \
         -I"${INSTALLDIR}/libtorch-${libtorch_ver}/include/torch/csrc/api/include" \
         -I"${INSTALLDIR}/libtorch-${libtorch_ver}/include" \
         -I"${INSTALLDIR}/libtorch-${libtorch_ver}/include/torch" \
         -L"${INSTALLDIR}/libtorch-${libtorch_ver}/lib" \
         -o libskala_ftorch_bindings.so \
         "${BUILDDIR}/skala_ftorch.cpp" \
         -ltorch -ltorch_cpu -lc10 \
         -lpthread -ldl \
         -Wl,-rpath,"${INSTALLDIR}/libtorch-${libtorch_ver}/lib" \
         > build_bindings.log 2>&1; then
         echo "WARNING: Failed to build scala C++ bindings. See build_bindings.log for details."
         echo "Skala-ftorch functionality will be limited without the C++ bindings."
         echo "The Fortran interface file has been installed, but runtime will fail without the C++ library."
       fi

       mkdir -p "${pkg_install_dir}/lib"
       if [ -f libskala_ftorch_bindings.so ]; then
         cp libskala_ftorch_bindings.so "${pkg_install_dir}/lib/"
       else
         echo "WARNING: libskala_ftorch_bindings.so not found, skipping copy"
       fi

       mkdir -p "${pkg_install_dir}/include"
       cp "${BUILDDIR}/skala_ftorch.f90" "${pkg_install_dir}/include/"

      # Copy skala model from shared location (like other dependencies)
      mkdir -p "${pkg_install_dir}/share/skala/onedft_models"
      if [ -f "${INSTALLDIR}/skala-${skala_model_ver}/share/skala/onedft_models/${skala_model_pkg}" ]; then
        cp "${INSTALLDIR}/skala-${skala_model_ver}/share/skala/onedft_models/${skala_model_pkg}" \
          "${pkg_install_dir}/share/skala/onedft_models/${skala_model_pkg}"
      fi

      write_checksums "${install_lock_file}" \
        "${SCRIPT_DIR}/stage6/$(basename ${SCRIPT_NAME})" \
         "${BUILDDIR}/${ftorch_pkg}"
    fi
    ;;
  __SYSTEM__)
    echo "==================== Finding FTorch from system paths ===================="
    check_lib -lftorch "FTorch"
    pkg_install_dir="$(dirname $(dirname $(find_in_paths "libftorch.*" $LIB_PATHS)))"
    ;;
  __DONTUSE__) ;;

  *)
    echo "==================== Linking skala-ftorch to user paths ===================="
    pkg_install_dir="${with_skala_ftorch}"
    FTORCH_LIBDIR="${pkg_install_dir}/lib"
    [ -d "${pkg_install_dir}/lib64" ] && FTORCH_LIBDIR="${pkg_install_dir}/lib64"
    check_dir "${FTORCH_LIBDIR}"
    check_dir "${pkg_install_dir}/include"
    ;;
esac

if [ "${with_skala_ftorch}" != "__DONTUSE__" ]; then
  if [ -n "${pkg_install_dir:-}" ]; then
    # Use shared skala model location
    skala_model="${INSTALLDIR}/skala-${skala_model_ver}/share/skala/onedft_models/${skala_model_pkg}"
  else
    skala_model=""
  fi
  cat << EOF > "${BUILDDIR}/setup_ftorch"
export FTORCH_VER="${ftorch_ver}"
export FTORCH_ROOT="${pkg_install_dir}"
export SKALA_MODEL="${skala_model}"
EOF
  if [ "${with_skala_ftorch}" != "__SYSTEM__" ]; then
    cat << EOF >> "${BUILDDIR}/setup_ftorch"
prepend_path LD_LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path LD_RUN_PATH "${pkg_install_dir}/lib"
prepend_path LIBRARY_PATH "${pkg_install_dir}/lib"
prepend_path CMAKE_PREFIX_PATH "${pkg_install_dir}"
EOF
  fi
  filter_setup "${BUILDDIR}/setup_ftorch" "${SETUPFILE}"
fi

load "${BUILDDIR}/setup_ftorch"
write_toolchain_env "${INSTALLDIR}"

cd "${ROOTDIR}"
report_timing "ftorch"