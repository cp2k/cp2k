#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

#
# CMake recipes https://github.com/eth-cscs/cmake-recipes
#
# Copyright (c) 2018-2019, ETH Zurich BSD 3-Clause License. All rights reserved.
#
# Author: Teodor Nikolov (teodor.nikolov22@gmail.com)
#

#[=======================================================================[.rst:
FindMKL
-------

Find the oneMKL configuration requested by CP2K and provide these targets:

``cp2k::BLAS::MKL::MKL``
  BLAS/LAPACK/FFT link closure for the selected integer, threading and link mode.

``cp2k::BLAS::MKL::blas`` and ``cp2k::BLAS::MKL::lapack``
  Aliases of ``cp2k::BLAS::MKL::MKL``.

``cp2k::BLAS::MKL::scalapack_link``
  ScaLAPACK/BLACS/MPI link closure, available for MPI builds.

Search hints
^^^^^^^^^^^^

``MKL_ROOT``, ``MKLROOT`` and the ``MKLROOT`` environment variable are searched.

CP2K options
^^^^^^^^^^^^

``CP2K_MKL_LINK``
  ``dynamic`` (default) or ``static``.

``CP2K_MKL_MPI``
  ``auto`` (default), ``openmpi``, ``intelmpi`` or ``mpich``. Relevant only for
  MPI builds. Set this explicitly when cross compiling.

The selected oneMKL interface and threading layer follow ``CP2K_BLAS_INTERFACE``
and ``CP2K_BLAS_THREADING``.

#]=======================================================================]

# Copyright (c) 2022- ETH Zurich
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

include(FindPackageHandleStandardArgs)

if(NOT CMAKE_Fortran_COMPILER_LOADED)
  message(FATAL_ERROR "FindMKL requires Fortran to be enabled.")
endif()

# CP2K configuration
set(CP2K_MKL_LINK
    "dynamic"
    CACHE STRING "oneMKL link mode")
set_property(CACHE CP2K_MKL_LINK PROPERTY STRINGS dynamic static)
if(NOT CP2K_MKL_LINK MATCHES "^(dynamic|static)$")
  message(FATAL_ERROR "CP2K_MKL_LINK must be dynamic or static")
endif()

if(NOT CP2K_BLAS_INTERFACE MATCHES "^(32bits|64bits)$")
  message(FATAL_ERROR "CP2K_BLAS_INTERFACE must be 32bits or 64bits")
endif()

if(CP2K_BLAS_INTERFACE STREQUAL "64bits")
  set(_mkl_int_interface ilp64)
else()
  set(_mkl_int_interface lp64)
endif()

# GNU Fortran needs the GNU Fortran interface. Other compilers use Intel's
# generic Fortran interface. The GNU interface is not available on macOS.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND NOT APPLE)
  set(_mkl_interface_name mkl_gf_${_mkl_int_interface})
else()
  set(_mkl_interface_name mkl_intel_${_mkl_int_interface})
endif()

# Resolve only the requested threading layer. CP2K itself uses OpenMP, so do not
# mix GNU and Intel OpenMP runtimes in the same process.
set(_mkl_threading_supported TRUE)
set(_mkl_thread_dependency_found TRUE)
set(_mkl_thread_libraries)
set(_mkl_use_openmp FALSE)

if(CP2K_BLAS_THREADING STREQUAL "sequential")
  set(_mkl_thread_name mkl_sequential)
elseif(CP2K_BLAS_THREADING STREQUAL "tbb-thread")
  set(_mkl_thread_name mkl_tbb_thread)
  find_package(TBB QUIET CONFIG HINTS "${TBB_ROOT}" "$ENV{TBBROOT}")
  if(TARGET TBB::tbb)
    list(APPEND _mkl_thread_libraries TBB::tbb)
  else()
    set(_mkl_thread_dependency_found FALSE)
  endif()
elseif(CP2K_BLAS_THREADING STREQUAL "gnu-thread")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(_mkl_thread_name mkl_gnu_thread)
    set(_mkl_use_openmp TRUE)
  else()
    set(_mkl_threading_supported FALSE)
  endif()
elseif(CP2K_BLAS_THREADING STREQUAL "intel-thread")
  if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
    set(_mkl_thread_name mkl_intel_thread)
    set(_mkl_use_openmp TRUE)
  else()
    set(_mkl_threading_supported FALSE)
  endif()
elseif(CP2K_BLAS_THREADING MATCHES "^(thread|openmp)$")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(_mkl_thread_name mkl_gnu_thread)
    set(_mkl_use_openmp TRUE)
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
    set(_mkl_thread_name mkl_intel_thread)
    set(_mkl_use_openmp TRUE)
  else()
    set(_mkl_threading_supported FALSE)
  endif()
else()
  set(_mkl_threading_supported FALSE)
endif()

if(_mkl_use_openmp)
  find_package(OpenMP QUIET COMPONENTS Fortran)
  if(TARGET OpenMP::OpenMP_Fortran)
    list(APPEND _mkl_thread_libraries OpenMP::OpenMP_Fortran)
  else()
    set(_mkl_thread_dependency_found FALSE)
  endif()
endif()

# Search roots. Keep both historical MKL_ROOT and Intel's MKLROOT spelling.
set(_mkl_root_hints "${MKL_ROOT}" "${MKLROOT}" "$ENV{MKLROOT}")
list(REMOVE_ITEM _mkl_root_hints "")
list(REMOVE_DUPLICATES _mkl_root_hints)

# Determine architecture and exact library file names. Exact names let static
# and dynamic selection remain local to this finder rather than modifying
# CMAKE_FIND_LIBRARY_SUFFIXES globally.
set(_mkl_libpath_suffix intel64)
if(CMAKE_SIZEOF_VOID_P EQUAL 4)
  set(_mkl_libpath_suffix ia32)
endif()
set(_mkl_libpath_suffixes ${_mkl_libpath_suffix})

if(WIN32)
  list(APPEND _mkl_libpath_suffixes ${_mkl_libpath_suffix}_win)
  set(_mkl_libname_prefix "")
  if(CP2K_MKL_LINK STREQUAL "static")
    set(_mkl_library_suffix ".lib")
  else()
    set(_mkl_library_suffix "_dll.lib")
  endif()
elseif(APPLE)
  list(APPEND _mkl_libpath_suffixes ${_mkl_libpath_suffix}_mac)
  set(_mkl_libname_prefix "lib")
  if(CP2K_MKL_LINK STREQUAL "static")
    set(_mkl_library_suffix ".a")
  else()
    set(_mkl_library_suffix ".dylib")
  endif()
else()
  list(APPEND _mkl_libpath_suffixes ${_mkl_libpath_suffix}_lin)
  set(_mkl_libname_prefix "lib")
  if(CP2K_MKL_LINK STREQUAL "static")
    set(_mkl_library_suffix ".a")
  else()
    set(_mkl_library_suffix ".so")
  endif()
endif()

set(_mkl_search_paths)
foreach(_mkl_root IN LISTS _mkl_root_hints)
  list(APPEND _mkl_search_paths "${_mkl_root}" "${_mkl_root}/lib"
       "${_mkl_root}/mkl/lib" "${_mkl_root}/compiler/lib")
endforeach()
list(REMOVE_DUPLICATES _mkl_search_paths)

function(_cp2k_mkl_find_library _var _name)
  # These are finder results, not user inputs. Re-search when MKL_ROOT or the
  # requested link mode changes in an existing build directory.
  unset(${_var} CACHE)
  unset(${_var})
  find_library(
    ${_var}
    NAMES "${_mkl_libname_prefix}${_name}${_mkl_library_suffix}"
    HINTS ${_mkl_search_paths}
    PATH_SUFFIXES ${_mkl_libpath_suffixes})
  mark_as_advanced(${_var})
  set(${_var}
      "${${_var}}"
      PARENT_SCOPE)
endfunction()

unset(_mkl_include_dir CACHE)
find_path(
  _mkl_include_dir mkl.h
  HINTS ${_mkl_root_hints}
  PATH_SUFFIXES include mkl/include)
mark_as_advanced(_mkl_include_dir)
set(CP2K_MKL_INCLUDE_DIRS "${_mkl_include_dir}")

_cp2k_mkl_find_library(_mkl_interface_library ${_mkl_interface_name})
if(_mkl_threading_supported)
  _cp2k_mkl_find_library(_mkl_thread_library ${_mkl_thread_name})
endif()
_cp2k_mkl_find_library(_mkl_core_library mkl_core)

find_package(Threads QUIET)

# Resolve the BLACS implementation only when CP2K actually needs ScaLAPACK.
set(_mkl_mpi_interface_found TRUE)
if(CP2K_USE_MPI)
  set(CP2K_MKL_MPI
      "auto"
      CACHE STRING "oneMKL BLACS MPI interface")
  set_property(CACHE CP2K_MKL_MPI PROPERTY STRINGS auto openmpi intelmpi mpich)
  if(NOT CP2K_MKL_MPI MATCHES "^(auto|openmpi|intelmpi|mpich)$")
    message(
      FATAL_ERROR "CP2K_MKL_MPI must be auto, openmpi, intelmpi, or mpich")
  endif()

  set(_mkl_mpi_interface ${CP2K_MKL_MPI})
  if(_mkl_mpi_interface STREQUAL "auto")
    if(CMAKE_CROSSCOMPILING)
      set(_mkl_mpi_interface_found FALSE)
    elseif(MPI_Fortran_LIBRARY_VERSION_STRING MATCHES "Open MPI")
      set(_mkl_mpi_interface openmpi)
    elseif(MPI_Fortran_LIBRARY_VERSION_STRING MATCHES
           "Intel\\(R\\) MPI|Intel MPI|MPICH|HYDRA")
      set(_mkl_mpi_interface intelmpi)
    else()
      set(_mkl_mpi_interface_found FALSE)
    endif()
  elseif(_mkl_mpi_interface STREQUAL "mpich")
    set(_mkl_mpi_interface intelmpi)
  endif()

  if(_mkl_mpi_interface_found)
    if(_mkl_mpi_interface STREQUAL "openmpi")
      set(_mkl_blacs_name mkl_blacs_openmpi_${_mkl_int_interface})
    elseif(APPLE)
      set(_mkl_blacs_name mkl_blacs_mpich_${_mkl_int_interface})
    else()
      set(_mkl_blacs_name mkl_blacs_intelmpi_${_mkl_int_interface})
    endif()

    _cp2k_mkl_find_library(_mkl_scalapack_library
                           mkl_scalapack_${_mkl_int_interface})
    _cp2k_mkl_find_library(_mkl_blacs_library ${_mkl_blacs_name})
  endif()
endif()

set(_mkl_required_vars
    CP2K_MKL_INCLUDE_DIRS
    _mkl_interface_library
    _mkl_thread_library
    _mkl_core_library
    Threads_FOUND
    _mkl_threading_supported
    _mkl_thread_dependency_found)
if(CP2K_USE_MPI)
  list(APPEND _mkl_required_vars _mkl_mpi_interface_found
       _mkl_scalapack_library _mkl_blacs_library MPI_Fortran_FOUND)
endif()

set(_mkl_failure_reason)
if(NOT _mkl_threading_supported)
  set(_mkl_failure_reason
      "CP2K_BLAS_THREADING=${CP2K_BLAS_THREADING} is incompatible with ${CMAKE_Fortran_COMPILER_ID} Fortran and oneMKL."
  )
elseif(CP2K_USE_MPI AND NOT _mkl_mpi_interface_found)
  set(_mkl_failure_reason
      "Cannot determine the oneMKL BLACS MPI interface. Set CP2K_MKL_MPI explicitly, especially when cross compiling."
  )
elseif(NOT _mkl_thread_dependency_found)
  set(_mkl_failure_reason
      "The runtime dependency for CP2K_BLAS_THREADING=${CP2K_BLAS_THREADING} was not found."
  )
endif()

find_package_handle_standard_args(
  MKL REQUIRED_VARS ${_mkl_required_vars} REASON_FAILURE_MESSAGE
                    "${_mkl_failure_reason}")

function(_cp2k_mkl_link_group _output)
  # oneMKL static archives contain circular references. CP2K already requires
  # CMake 3.24, which provides LINK_GROUP:RESCAN for this purpose.
  if(CP2K_MKL_LINK STREQUAL "static"
     AND (CMAKE_Fortran_LINK_GROUP_USING_RESCAN_SUPPORTED
          OR CMAKE_LINK_GROUP_USING_RESCAN_SUPPORTED))
    list(JOIN ARGN "," _mkl_group_items)
    set(${_output}
        "$<LINK_GROUP:RESCAN,${_mkl_group_items}>"
        PARENT_SCOPE)
  else()
    set(${_output}
        "${ARGN}"
        PARENT_SCOPE)
  endif()
endfunction()

if(MKL_FOUND)
  set(_mkl_base_libraries ${_mkl_interface_library} ${_mkl_thread_library}
                          ${_mkl_core_library})
  _cp2k_mkl_link_group(_mkl_base_link ${_mkl_base_libraries})

  set(_mkl_system_libraries Threads::Threads)
  if(UNIX)
    list(APPEND _mkl_system_libraries m)
    if(CP2K_MKL_LINK STREQUAL "static" AND CMAKE_DL_LIBS)
      list(APPEND _mkl_system_libraries ${CMAKE_DL_LIBS})
    endif()
  endif()

  if(NOT TARGET cp2k::BLAS::MKL::MKL)
    add_library(cp2k::BLAS::MKL::MKL INTERFACE IMPORTED)
    add_library(cp2k::BLAS::MKL::blas ALIAS cp2k::BLAS::MKL::MKL)
    add_library(cp2k::BLAS::MKL::lapack ALIAS cp2k::BLAS::MKL::MKL)
  endif()
  set_target_properties(
    cp2k::BLAS::MKL::MKL
    PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${CP2K_MKL_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES
      "${_mkl_base_link};${_mkl_thread_libraries};${_mkl_system_libraries}")

  if(CP2K_USE_MPI)
    set(_mkl_cluster_libraries
        ${_mkl_scalapack_library} ${_mkl_interface_library}
        ${_mkl_thread_library} ${_mkl_core_library} ${_mkl_blacs_library})
    _cp2k_mkl_link_group(_mkl_cluster_link ${_mkl_cluster_libraries})

    if(NOT TARGET cp2k::BLAS::MKL::scalapack_link)
      add_library(cp2k::BLAS::MKL::scalapack_link INTERFACE IMPORTED)
    endif()
    set_target_properties(
      cp2k::BLAS::MKL::scalapack_link
      PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${CP2K_MKL_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES
        "${_mkl_cluster_link};MPI::MPI_Fortran;${_mkl_thread_libraries};${_mkl_system_libraries}"
    )
  endif()

  set(CP2K_BLAS_VENDOR "MKL")
endif()

mark_as_advanced(CP2K_BLAS_VENDOR)
