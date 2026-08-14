#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

# Locate a symmetrix checkout + build tree and derive the include directories
# and libraries the CP2K symmetrix backend needs.
#
# libsymmetrix currently has no install()/package-config rules, so it cannot be
# consumed with a plain find_package() against an install prefix. Instead the
# user points at the checkout and its build directory:
#
# -DCP2K_SYMMETRIX_ROOT=<symmetrix checkout> -DCP2K_SYMMETRIX_BUILD_DIR=<build
# dir>   (default: ${CP2K_SYMMETRIX_ROOT}/build)
#
# built beforehand with e.g.
#
# cmake -S <checkout>/libsymmetrix -B <build dir> \ -DCMAKE_BUILD_TYPE=Release
# -DSYMMETRIX_KOKKOS=OFF     # host build cmake --build <build dir> -j
#
# (add -DSYMMETRIX_KOKKOS=ON and the Kokkos/CUDA options for a Kokkos build, and
# configure CP2K with -DCP2K_SYMMETRIX_KOKKOS=ON to match).
#
# Explicitly set CP2K_SYMMETRIX_INCLUDE_DIR + CP2K_SYMMETRIX_LIBRARIES always
# win and bypass this module entirely (see the top-level CMakeLists.txt); they
# remain the escape hatch for unusual layouts.
#
# Results (normal variables, consumed by src/CMakeLists.txt):
# CP2K_SYMMETRIX_INCLUDE_DIR   list of include directories
# CP2K_SYMMETRIX_LIBRARIES     list of libraries / imported targets to link

include(FindPackageHandleStandardArgs)

set(CP2K_SYMMETRIX_ROOT
    ""
    CACHE PATH "Path to a symmetrix source checkout")
set(CP2K_SYMMETRIX_BUILD_DIR
    ""
    CACHE PATH
          "libsymmetrix build directory (default: CP2K_SYMMETRIX_ROOT/build)")

# Fall back to the environment, which is how the CP2K toolchain
# (install_symmetrix.sh) advertises the checkout it built.
if(NOT CP2K_SYMMETRIX_ROOT AND DEFINED ENV{CP2K_SYMMETRIX_ROOT})
  set(CP2K_SYMMETRIX_ROOT "$ENV{CP2K_SYMMETRIX_ROOT}")
endif()
if(NOT CP2K_SYMMETRIX_BUILD_DIR AND DEFINED ENV{CP2K_SYMMETRIX_BUILD_DIR})
  set(CP2K_SYMMETRIX_BUILD_DIR "$ENV{CP2K_SYMMETRIX_BUILD_DIR}")
endif()

set(_sym_root "${CP2K_SYMMETRIX_ROOT}")
set(_sym_build "${CP2K_SYMMETRIX_BUILD_DIR}")

# Last, look for an installed prefix that keeps the checkout+build layout, which
# is what the spack package stages (libsymmetrix has no install rules of its
# own, so there is nothing for find_package to consume). Searching
# CMAKE_PREFIX_PATH and the usual roots means a spack- or module-provided
# symmetrix needs no manual paths.
if(NOT _sym_root)
  find_path(
    _sym_prefix_hint
    NAMES libsymmetrix/source/mace.hpp
    HINTS ${CMAKE_PREFIX_PATH} ENV CMAKE_PREFIX_PATH ENV SYMMETRIX_ROOT
    PATH_SUFFIXES symmetrix
    DOC "Prefix holding a symmetrix checkout plus its libsymmetrix build")
  if(_sym_prefix_hint)
    set(_sym_root "${_sym_prefix_hint}")
  endif()
endif()

if(_sym_root AND NOT _sym_build)
  set(_sym_build "${_sym_root}/build")
endif()

if(NOT _sym_root)
  message(
    FATAL_ERROR
      "CP2K_USE_SYMMETRIX=ON needs either CP2K_SYMMETRIX_ROOT (a symmetrix "
      "checkout with a completed libsymmetrix build, see "
      "https://github.com/wcwitt/symmetrix), a symmetrix prefix on "
      "CMAKE_PREFIX_PATH (as installed by the spack package), or explicit "
      "CP2K_SYMMETRIX_INCLUDE_DIR + CP2K_SYMMETRIX_LIBRARIES.")
endif()

if(NOT EXISTS "${_sym_root}/libsymmetrix/source/mace.hpp")
  message(
    FATAL_ERROR
      "CP2K_SYMMETRIX_ROOT=${_sym_root} does not look like a symmetrix "
      "checkout (missing libsymmetrix/source/mace.hpp).")
endif()

find_library(CP2K_SYMMETRIX_CORE_LIBRARY symmetrix PATHS "${_sym_build}")
if(NOT CP2K_SYMMETRIX_CORE_LIBRARY)
  message(
    FATAL_ERROR
      "libsymmetrix.a not found in CP2K_SYMMETRIX_BUILD_DIR=${_sym_build}. "
      "Build it first (SYMMETRIX_KOKKOS matching CP2K_SYMMETRIX_KOKKOS):\n"
      "  cmake -S ${_sym_root}/libsymmetrix -B ${_sym_build} "
      "-DCMAKE_BUILD_TYPE=Release -DSYMMETRIX_KOKKOS=${CP2K_SYMMETRIX_KOKKOS}\n"
      "  cmake --build ${_sym_build} -j")
endif()

find_library(
  CP2K_SYMMETRIX_SPHERICART_LIBRARY sphericart
  PATHS "${_sym_build}/external/sphericart/sphericart"
        "${_sym_build}/external/sphericart")
if(NOT CP2K_SYMMETRIX_SPHERICART_LIBRARY)
  message(
    FATAL_ERROR
      "libsphericart.a not found under ${_sym_build}/external/sphericart "
      "(it is built as part of libsymmetrix).")
endif()

# Include directories: the checkout's headers plus build-generated headers.
# Entries that a given symmetrix version/configuration does not have are skipped
# (e.g. the Kokkos trees in a host-only build).
set(_sym_incs_rel
    # host API (mace.hpp) + bundled deps
    "libsymmetrix/source"
    "libsymmetrix/external/json/single_include"
    "libsymmetrix/external/sphericart/sphericart/include"
    # vendored Kokkos + kokkos-kernels source trees (Kokkos builds)
    "libsymmetrix/external/kokkos/core/src"
    "libsymmetrix/external/kokkos/containers/src"
    "libsymmetrix/external/kokkos/algorithms/src"
    "libsymmetrix/external/kokkos/simd/src"
    "libsymmetrix/external/kokkos/tpls/desul/include"
    "libsymmetrix/external/kokkos/tpls/mdspan/include"
    "libsymmetrix/external/kokkos-kernels/common/src"
    "libsymmetrix/external/kokkos-kernels/common/impl"
    "libsymmetrix/external/kokkos-kernels/batched"
    "libsymmetrix/external/kokkos-kernels/batched/dense/src"
    "libsymmetrix/external/kokkos-kernels/batched/dense/impl"
    "libsymmetrix/external/kokkos-kernels/batched/sparse/src"
    "libsymmetrix/external/kokkos-kernels/batched/sparse/impl"
    "libsymmetrix/external/kokkos-kernels/blas/src"
    "libsymmetrix/external/kokkos-kernels/blas/impl"
    "libsymmetrix/external/kokkos-kernels/blas/eti"
    "libsymmetrix/external/kokkos-kernels/blas/tpls"
    "libsymmetrix/external/kokkos-kernels/lapack/src"
    "libsymmetrix/external/kokkos-kernels/lapack/impl"
    "libsymmetrix/external/kokkos-kernels/lapack/eti"
    "libsymmetrix/external/kokkos-kernels/lapack/tpls"
    "libsymmetrix/external/kokkos-kernels/graph/src"
    "libsymmetrix/external/kokkos-kernels/graph/impl"
    "libsymmetrix/external/kokkos-kernels/graph/eti"
    "libsymmetrix/external/kokkos-kernels/sparse/src"
    "libsymmetrix/external/kokkos-kernels/sparse/impl"
    "libsymmetrix/external/kokkos-kernels/sparse/eti"
    "libsymmetrix/external/kokkos-kernels/sparse/tpls"
    "libsymmetrix/external/kokkos-kernels/ode/src"
    "libsymmetrix/external/kokkos-kernels/ode/impl")
set(_sym_build_incs_rel
    # build-generated headers (configuration and ETI stubs)
    "external/sphericart/sphericart/include"
    "external/kokkos"
    "external/kokkos/core/src"
    "external/kokkos/containers/src"
    "external/kokkos/algorithms/src"
    "external/kokkos/simd/src"
    "external/kokkos-kernels"
    "external/kokkos-kernels/batched/eti"
    "external/kokkos-kernels/blas/eti"
    "external/kokkos-kernels/lapack/eti"
    "external/kokkos-kernels/graph/eti"
    "external/kokkos-kernels/sparse/eti")

set(CP2K_SYMMETRIX_INCLUDE_DIR "")
foreach(_d IN LISTS _sym_incs_rel)
  if(IS_DIRECTORY "${_sym_root}/${_d}")
    list(APPEND CP2K_SYMMETRIX_INCLUDE_DIR "${_sym_root}/${_d}")
  endif()
endforeach()
foreach(_d IN LISTS _sym_build_incs_rel)
  if(IS_DIRECTORY "${_sym_build}/${_d}")
    list(APPEND CP2K_SYMMETRIX_INCLUDE_DIR "${_sym_build}/${_d}")
  endif()
endforeach()

set(CP2K_SYMMETRIX_LIBRARIES "${CP2K_SYMMETRIX_CORE_LIBRARY}")

if(CP2K_SYMMETRIX_KOKKOS)
  if(NOT EXISTS "${_sym_root}/libsymmetrix/source/mace_kokkos.hpp")
    message(
      FATAL_ERROR
        "CP2K_SYMMETRIX_KOKKOS=ON but ${_sym_root}/libsymmetrix/source has no "
        "mace_kokkos.hpp; update the symmetrix checkout.")
  endif()
  foreach(_kklib kokkoskernels kokkoscontainers kokkoscore kokkossimd)
    find_library(
      CP2K_SYMMETRIX_${_kklib}_LIBRARY ${_kklib}
      PATHS "${_sym_build}/external/kokkos-kernels"
            "${_sym_build}/external/kokkos/containers/src"
            "${_sym_build}/external/kokkos/core/src"
            "${_sym_build}/external/kokkos/simd/src")
    if(NOT CP2K_SYMMETRIX_${_kklib}_LIBRARY)
      message(
        FATAL_ERROR
          "lib${_kklib}.a not found under ${_sym_build}/external; rebuild "
          "libsymmetrix with -DSYMMETRIX_KOKKOS=ON to match "
          "CP2K_SYMMETRIX_KOKKOS=ON.")
    endif()
    list(APPEND CP2K_SYMMETRIX_LIBRARIES "${CP2K_SYMMETRIX_${_kklib}_LIBRARY}")
  endforeach()
endif()

list(APPEND CP2K_SYMMETRIX_LIBRARIES "${CP2K_SYMMETRIX_SPHERICART_LIBRARY}")

if(CP2K_SYMMETRIX_KOKKOS)
  # A CUDA-enabled Kokkos additionally needs the CUDA runtime/driver and the
  # TPLs kokkos-kernels was built against.
  set(_sym_kokkos_cfg "${_sym_build}/external/kokkos/KokkosCore_config.h")
  set(_sym_kokkos_cuda FALSE)
  if(EXISTS "${_sym_kokkos_cfg}")
    file(STRINGS "${_sym_kokkos_cfg}" _sym_cuda_line
         REGEX "#define KOKKOS_ENABLE_CUDA")
    if(_sym_cuda_line)
      set(_sym_kokkos_cuda TRUE)
    endif()
  endif()
  if(_sym_kokkos_cuda)
    find_package(CUDAToolkit REQUIRED)
    list(APPEND CP2K_SYMMETRIX_LIBRARIES CUDA::cudart CUDA::cublas
         CUDA::cusparse CUDA::cuda_driver)
    if(TARGET CUDA::nvJitLink)
      list(APPEND CP2K_SYMMETRIX_LIBRARIES CUDA::nvJitLink)
    endif()
  endif()
endif()

find_package_handle_standard_args(
  Symmetrix DEFAULT_MSG CP2K_SYMMETRIX_CORE_LIBRARY
  CP2K_SYMMETRIX_SPHERICART_LIBRARY)

mark_as_advanced(
  CP2K_SYMMETRIX_CORE_LIBRARY CP2K_SYMMETRIX_SPHERICART_LIBRARY
  CP2K_SYMMETRIX_kokkoskernels_LIBRARY CP2K_SYMMETRIX_kokkoscontainers_LIBRARY
  CP2K_SYMMETRIX_kokkoscore_LIBRARY CP2K_SYMMETRIX_kokkossimd_LIBRARY)
