#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
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

The following conventions are used:

intel / INTEL  - Bindings for everything except GNU Fortran
gf / GF        - GNU Fortran bindings
seq / SEQ      - sequential MKL
omp / OMP      - threaded MKL with OpenMP back end
tbb / TBB      - threaded MKL with TBB back end
32bit / 32BIT  - MKL 32 bit integer interface (used most often)
64bit / 64BIT  - MKL 64 bit integer interface
mpich / MPICH  - MPICH / IntelMPI BLACS back end
ompi / OMPI    - OpenMPI BLACS back end
st / ST        - static libraries
dyn / DYN      - dynamic libraries

The module attempts to define a target for each MKL configuration. The
configuration will not be available if there are missing library files or a
missing dependency.

MKL Link line advisor:
https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

Note: Mixing GCC and Intel OpenMP backends is a bad idea.

Search variables
^^^^^^^^^^^^^^^^

``MKLROOT``
Environment variable set to MKL's root directory

``MKL_ROOT``
CMake variable set to MKL's root directory

Example usage
^^^^^^^^^^^^^

To Find MKL:

find_package(MKL REQUIRED)

To check if target is available:

if (TARGET MKL::scalapack_mpich_intel_32bit_omp_dyn)
  ...
endif()

To link to an available target (see list below):

target_link_libraries(... MKL::scalapack_mpich_intel_32bit_omp_dyn)

Note: dependencies are handled for you (MPI, OpenMP, ...)

the target MKL::blas, MKL::MKL, MKL::lapack also include all necessary libraries
for linking.
MKL::MKL is also used by the cmake module provided by intel.

MKL::scalapack_link gives all libraries needed for scalapack.

Imported targets
^^^^^^^^^^^^^^^^

MKL (BLAS, LAPACK, FFT) targets:

MKL::[gf|intel]_[32bit|64bit]_[seq|omp|tbb]_[st|dyn] e.g.

MKL::mkl_intel_32bit_omp_dyn

BLACS targets:

MKL::blacs_[mpich|ompi]_[gf|intel]_[32bit|64bit]_[seq|omp|tbb]_[st|dyn] e.g.

MKL::blacs_intel_mpich_32bit_seq_st


ScaLAPACK targets:

MKL::scalapack_[mpich|ompi]_[gf|intel]_[32bit|64bit]_[seq|omp|tbb]_[st|dyn] e.g.

MKL::scalapack_mpich_intel_64bit_omp_dyn

Result variables
^^^^^^^^^^^^^^^^

MKL_FOUND

Not supported
^^^^^^^^^^^^^

- F95 interfaces

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

if(NOT
   (CMAKE_C_COMPILER_LOADED
    OR CMAKE_CXX_COMPILER_LOADED
    OR CMAKE_Fortran_COMPILER_LOADED))
  message(FATAL_ERROR "FindMKL requires Fortran, C, or C++ to be enabled.")
endif()

# Dependencies
#
find_package(Threads)
find_package(MPI COMPONENTS CXX C Fortran)
find_package(OpenMP COMPONENTS CXX C Fortran)

# If MKL_ROOT is not set, set it via the env variable MKLROOT.
#
if(NOT DEFINED MKL_ROOT)
  set(MKL_ROOT
      $ENV{MKLROOT}
      CACHE PATH "MKL's root directory.")
endif()

# Determine MKL's library folder
#
set(_mkl_libpath_suffix "intel64")
if(CMAKE_SIZEOF_VOID_P EQUAL 4) # 32 bit
  set(_mkl_libpath_suffix "ia32")
endif()

if(WIN32)
  list(APPEND _mkl_libpath_suffix_list ${_mkl_libpath_suffix})
  string(APPEND _mkl_libpath_suffix "_win")
  list(APPEND _mkl_libpath_suffix_list ${_mkl_libpath_suffix})
  set(_mkl_libname_prefix "")
  set(_mkl_shared_lib "_dll.lib")
  set(_mkl_static_lib ".lib")
elseif(APPLE)
  list(APPEND _mkl_libpath_suffix_list ${_mkl_libpath_suffix})
  string(APPEND _mkl_libpath_suffix "_mac")
  list(APPEND _mkl_libpath_suffix_list ${_mkl_libpath_suffix})
  set(_mkl_libname_prefix "lib")
  set(_mkl_shared_lib ".dylib")
  set(_mkl_static_lib ".a")
else() # LINUX
  list(APPEND _mkl_libpath_suffix_list ${_mkl_libpath_suffix})
  string(APPEND _mkl_libpath_suffix "_lin")
  list(APPEND _mkl_libpath_suffix_list ${_mkl_libpath_suffix})
  set(_mkl_libname_prefix "lib")
  set(_mkl_shared_lib ".so")
  set(_mkl_static_lib ".a")
endif()
set(_mkl_search_paths "${MKL_ROOT}" "${MKL_ROOT}/lib" "${MKL_ROOT}/mkl/lib"
                      "${MKL_ROOT}/compiler/lib")

# Functions: finds both static and shared MKL libraries
#
function(__mkl_find_library _varname _libname)
  find_library(
    ${_varname}_DYN
    NAMES ${_mkl_libname_prefix}${_libname}${_mkl_shared_lib}
    HINTS ${_mkl_search_paths}
    PATH_SUFFIXES ${_mkl_libpath_suffix_list})
  mark_as_advanced(${_varname}_DYN)
  find_library(
    ${_varname}_ST
    NAMES ${_mkl_libname_prefix}${_libname}${_mkl_static_lib}
    HINTS ${_mkl_search_paths}
    PATH_SUFFIXES ${_mkl_libpath_suffix_list})
  mark_as_advanced(${_varname}_ST)
endfunction()

# Find MKL headers
#
find_path(CP2K_MKL_INCLUDE_DIRS mkl.h HINTS ${MKL_ROOT}/include
                                            ${MKL_ROOT}/mkl/include)
mark_as_advanced(CP2K_MKL_INCLUDE_DIRS)

# Group flags for static libraries on Linux (GNU, PGI, ICC -> same linker)
#
if(UNIX AND NOT APPLE)
  set(_mkl_linker_pre_flags_ST "-Wl,--start-group")
  set(_mkl_linker_post_flags_ST "-Wl,--end-group")
endif()

# Core MKL
#
__mkl_find_library(MKL_CORE_LIB mkl_core)

# Interface
#
__mkl_find_library(MKL_INTERFACE_INTEL_32BIT_LIB mkl_intel_lp64)
__mkl_find_library(MKL_INTERFACE_INTEL_64BIT_LIB mkl_intel_ilp64)
if(NOT APPLE
   AND CMAKE_Fortran_COMPILER_LOADED
   AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  __mkl_find_library(MKL_INTERFACE_GF_32BIT_LIB mkl_gf_lp64)
  __mkl_find_library(MKL_INTERFACE_GF_64BIT_LIB mkl_gf_ilp64)
endif()

# Threading
#
__mkl_find_library(MKL_SEQ_LIB mkl_sequential)
if(NOT APPLE
   AND (CMAKE_C_COMPILER_ID STREQUAL "GNU"
        OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
        OR CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"))
  __mkl_find_library(MKL_OMP_LIB mkl_gnu_thread)
else()
  __mkl_find_library(MKL_OMP_LIB mkl_intel_thread)
endif()
__mkl_find_library(MKL_TBB_LIB mkl_tbb_thread)

# BLACS
#
if(APPLE)
  __mkl_find_library(MKL_BLACS_MPICH_32BIT_LIB mkl_blacs_mpich_lp64)
  __mkl_find_library(MKL_BLACS_MPICH_64BIT_LIB mkl_blacs_mpich_ilp64)
else()
  __mkl_find_library(MKL_BLACS_MPICH_32BIT_LIB mkl_blacs_intelmpi_lp64)
  __mkl_find_library(MKL_BLACS_MPICH_64BIT_LIB mkl_blacs_intelmpi_ilp64)
endif()
__mkl_find_library(MKL_BLACS_OMPI_32BIT_LIB mkl_blacs_openmpi_lp64)
__mkl_find_library(MKL_BLACS_OMPI_64BIT_LIB mkl_blacs_openmpi_ilp64)

# ScaLAPACK
#
__mkl_find_library(MKL_SCALAPACK_32BIT_LIB mkl_scalapack_lp64)
__mkl_find_library(MKL_SCALAPACK_64BIT_LIB mkl_scalapack_ilp64)

# Check if core libs were found
#
find_package_handle_standard_args(MKL REQUIRED_VARS CP2K_MKL_INCLUDE_DIRS
                                                    Threads_FOUND)

# Sequential has no threading dependency. There is currently no TBB module
# shipped with CMake. The dependency is not accounted for. (FIXME)
#
set(_mkl_dep_found_SEQ TRUE)
set(_mkl_dep_found_TBB TRUE)
if(TARGET OpenMP::OpenMP_CXX)
  set(_mkl_dep_OMP ${OpenMP_CXX_LIBRARIES})
  set(_mkl_dep_found_OMP TRUE)
endif()

# Define all blas, blacs and scalapack
#
foreach(_libtype "ST" "DYN")
  set(_mkl_core_lib ${MKL_CORE_LIB_${_libtype}})
  foreach(_bits "32BIT" "64BIT")
    set(_mkl_scalapack_lib ${MKL_SCALAPACK_${_bits}_LIB_${_libtype}})
    foreach(_iface "INTEL" "GF")
      set(_mkl_interface_lib
          ${MKL_INTERFACE_${_iface}_${_bits}_LIB_${_libtype}})
      foreach(_threading "SEQ" "OMP" "TBB")
        set(_mkl_threading_lib ${MKL_${_threading}_LIB_${_libtype}})

        string(TOLOWER "${_iface}_${_bits}_${_threading}_${_libtype}"
                       _tgt_config)
        set(_mkl_tgt cp2k::BLAS::MKL::${_tgt_config})

        if(MKL_FOUND
           AND _mkl_interface_lib
           AND _mkl_threading_lib
           AND _mkl_core_lib
           AND _mkl_dep_found_${_threading}
           AND NOT TARGET ${_mkl_tgt})
          set(_mkl_libs
              "${_mkl_linker_pre_flags_${_threading}}"
              "${_mkl_interface_lib}"
              "${_mkl_threading_lib}"
              "${_mkl_core_lib}"
              "${_mkl_linker_post_flags_${_threading}}"
              "${_mkl_dep_${_threading}}"
              "Threads::Threads")
          add_library(${_mkl_tgt} INTERFACE IMPORTED)
          set_target_properties(
            ${_mkl_tgt}
            PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_MKL_INCLUDE_DIRS}"
                       INTERFACE_LINK_LIBRARIES "${_mkl_libs}")
        endif()

        foreach(_mpi_impl "MPICH" "OMPI")
          set(_mkl_blacs_lib ${MKL_BLACS_${_mpi_impl}_${_bits}_LIB_${_libtype}})

          string(
            TOLOWER "${_mpi_impl}_${_iface}_${_bits}_${_threading}_${_libtype}"
                    _tgt_config)

          set(_scalapack_tgt cp2k::BLAS::MKL::scalapack_${_tgt_config})

          if(_mkl_blacs_lib
             AND TARGET ${_mkl_tgt}
             AND TARGET MPI::MPI_CXX
             AND NOT TARGET cp2k::BLAS::MKL::blacs_${_tgt_config})
            set(_blacs_libs
                "${_mkl_linker_pre_flags_${_libtype}}"
                "${_mkl_interface_lib}"
                "${_mkl_threading_lib}"
                "${_mkl_core_lib}"
                "${_mkl_blacs_lib}"
                "${_mkl_linker_post_flags_${_libtype}}"
                "MPI::MPI_CXX"
                "${_mkl_dep_${_threading}}"
                "Threads::Threads")
            add_library(cp2k::BLAS::MKL::blacs_${_tgt_config} INTERFACE
                        IMPORTED)
            set_target_properties(
              cp2k::BLAS::MKL::blacs_${_tgt_config}
              PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                         "${CP2K_MKL_INCLUDE_DIRS}" INTERFACE_LINK_LIBRARIES
                                                    "${_mkl_blacs_lib}")
          endif()

          if(_mkl_scalapack_lib AND NOT TARGET
                                    cp2k::BLAS::MKL::scalapack_${_tgt_config})
            set(_scalapack_libs "${_mkl_scalapack_lib}" "${_blacs_tgt}")
            add_library(cp2k::BLAS::MKL::scalapack_${_tgt_config} INTERFACE
                        IMPORTED)
            set_target_properties(
              cp2k::BLAS::MKL::scalapack_${_tgt_config}
              PROPERTIES INTERFACE_LINK_LIBRARIES "${_scalapack_libs}")
          endif()
        endforeach()
      endforeach()
    endforeach()
  endforeach()
endforeach()

if(MKL_FOUND)
  # BLAS in the Intel MKL 10+ library?

  # the findMKL package finds all possible combination and define target for
  # each of them we just need to find which compiler we use, mpi etc...

  if(CMAKE_Fortran_COMPILER_LOADED
     AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"
     AND NOT APPLE)
    set(BLAS_mkl_INTFACE "gf")
  else()
    set(BLAS_mkl_INTFACE "intel")
  endif()

  if(CP2K_BLAS_THREADING MATCHES "thread|gnu-thread")
    set(BLAS_mkl_thread__ "omp")
  endif()

  if(CP2K_BLAS_THREADING MATCHES "sequential")
    set(BLAS_mkl_thread__ "seq")
  endif()

  if(CP2K_BLAS_THREADING MATCHES "intel-thread")
    set(BLAS_mkl_thread__ "intel")
  endif()

  if(CP2K_BLAS_THREADING MATCHES "tbb")
    set(BLAS_mkl_thread__ "tbb")
  endif()

  if(CP2K_BLAS_INTERFACE MATCHES "64bits")
    set(BLAS_mkl_ILP_MODE "64bit")
  else()
    set(BLAS_mkl_ILP_MODE "32bit")
  endif()

  get_target_property(
    MKL_BLAS_INCLUDE_DIRS
    cp2k::BLAS::MKL::${BLAS_mkl_INTFACE}_${BLAS_mkl_ILP_MODE}_${BLAS_mkl_thread__}_dyn
    INTERFACE_INCLUDE_DIRECTORIES)
  get_target_property(
    MKL_BLAS_LIBRARIES
    cp2k::BLAS::MKL::${BLAS_mkl_INTFACE}_${BLAS_mkl_ILP_MODE}_${BLAS_mkl_thread__}_dyn
    INTERFACE_LINK_LIBRARIES)
  if(NOT TARGET cp2k::BLAS::MKL::blas)
    add_library(cp2k::BLAS::MKL::MKL INTERFACE IMPORTED)
    add_library(cp2k::BLAS::MKL::blas ALIAS cp2k::BLAS::MKL::MKL)
    # create a empty lapack
    add_library(cp2k::BLAS::MKL::lapack INTERFACE IMPORTED)
  endif()
  set_target_properties(
    cp2k::BLAS::MKL::MKL
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_MKL_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${MKL_BLAS_LIBRARIES}")

  if("${MPI_Fortran_LIBRARY_VERSION_STRING}" MATCHES "Open MPI")
    set(__mkl_mpi_ver_ "ompi")
  else()
    set(__mkl_mpi_ver_ "mpich")
  endif()

  get_target_property(
    __mkl_scalapack_inc
    cp2k::BLAS::MKL::scalapack_${__mkl_mpi_ver_}_${BLAS_mkl_INTFACE}_${BLAS_mkl_ILP_MODE}_${BLAS_mkl_thread__}_dyn
    INTERFACE_INCLUDE_DIRECTORIES)
  get_target_property(
    __mkl_scalapack_lib
    cp2k::BLAS::MKL::scalapack_${__mkl_mpi_ver_}_${BLAS_mkl_INTFACE}_${BLAS_mkl_ILP_MODE}_${BLAS_mkl_thread__}_dyn
    INTERFACE_LINK_LIBRARIES)
  get_target_property(
    __mkl_blacs_inc
    cp2k::BLAS::MKL::blacs_${__mkl_mpi_ver_}_${BLAS_mkl_INTFACE}_${BLAS_mkl_ILP_MODE}_${BLAS_mkl_thread__}_dyn
    INTERFACE_INCLUDE_DIRECTORIES)
  get_target_property(
    __mkl_blacs_lib
    cp2k::BLAS::MKL::blacs_${__mkl_mpi_ver_}_${BLAS_mkl_INTFACE}_${BLAS_mkl_ILP_MODE}_${BLAS_mkl_thread__}_dyn
    INTERFACE_LINK_LIBRARIES)
  if(NOT TARGET cp2k::BLAS::MKL::scalapack_link)
    add_library(cp2k::BLAS::MKL::scalapack_link INTERFACE IMPORTED)
    set_target_properties(
      cp2k::BLAS::MKL::scalapack_link
      PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${__mkl_scalapack_inc}"
                 INTERFACE_LINK_LIBRARIES
                 "${__mkl_scalapack_lib};${__mkl_blacs_lib}")
  endif()
  unset(BLAS_mkl_ILP_MODE)
  unset(BLAS_mkl_INTFACE)
  unset(BLAS_mkl_thread__)
  unset(BLAS_mkl_OMP)
  unset(BLAS_mkl_OS_NAME)
  unset(__mkl_blacs_lib)
  unset(__mkl_blacs_inc)
  unset(__mkl_scalapack_lib)
  unset(__mkl_scalapack_inc)
  set(CP2K_BLAS_VENDOR "MKL")
  mark_as_advanced(CP2K_BLAS_VENDOR)
  mark_as_advanced(CP2K_MKL_FOUND)
endif()
