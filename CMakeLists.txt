#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

cmake_minimum_required(VERSION 3.22)

# include our cmake snippets
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# =================================================================================================
# REQUIRE OUT-OF-SOURCE BUILDS
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
  message(
    FATAL_ERROR
      "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory."
  )
endif()

# =================================================================================================
# PROJECT AND VERSION
include(CMakeDependentOption)
include(GitSubmodule)
include(CustomTargets)

cmake_policy(SET CMP0048 NEW)

if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()

set(VERSION_MAJOR 23)
# anything above 12 is to indicate that it is devel-branch. The git hash commit
# is the only important information to give. It is retrieved few lines below.
set(VERSION_MINOR 1)
set(VERSION_PATCH 00)

project(
  cp2k
  DESCRIPTION "CP2K"
  HOMEPAGE_URL "https://www.cp2k.org"
  VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
  LANGUAGES Fortran C CXX)

set(cp2k_APIVERSION ${cp2k_VERSION_MAJOR}.${cp2k_VERSION_MINOR})

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

if(NOT DEFINED CMAKE_CUDA_STANDARD)
  set(CMAKE_CUDA_STANDARD 11)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

# set language and standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_STANDARD 11)

find_package(PkgConfig)

# ##############################################################################
# Define the paths for static libraries and executables
# ##############################################################################
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${cp2k_BINARY_DIR}/lib
    CACHE PATH "Single output directory for building all libraries.")

# Search for rocm in common locations
foreach(__var ROCM_ROOT CRAY_ROCM_ROOT ORNL_ROCM_ROOT CRAY_ROCM_PREFIX
              ROCM_PREFIX CRAY_ROCM_DIR)
  if($ENV{${__var}})
    list(APPEND CMAKE_PREFIX_PATH $ENV{__var})
    set(ROCM_PATH
        $ENV{__var}
        CACHE PATH "Path to ROCm installation")
  endif()
endforeach()

# =================================================================================================
# OPTIONS

# option(CP2K_BUILD_DBCSR "Build dbcsr before building cp2k" OFF)
option(CP2K_ENABLE_REGTESTS
       "Enable installation of the binaries for running regtests afterwards"
       OFF)
option(CMAKE_POSITION_INDEPENDENT_CODE "Enable position independent code" ON)
option(CP2K_ENABLE_CONSISTENCY_CHECKS
       "Check that the list of compiled files and files contained in src match"
       OFF)
option(CP2K_DEBUG_MODE "Enable several additional options for debugging cp2k."
       OFF)
option(CP2K_USE_SIRIUS "Enable plane wave dft calculations with sirius" OFF)
option(CP2K_USE_FFTW3 "Use fftw3 for the calculating fast fourier transforms"
       ON)
option(CP2K_USE_MPI "Enable MPI support" ON)
option(CP2K_USE_ELPA "Enable elpa support" OFF)
option(CP2K_USE_PEXSI "Enable pexsi support" OFF)
option(CP2K_USE_SUPERLU "Enable superlu support" OFF)
option(CP2K_USE_COSMA "COSMA is a drop in replacement of scalapack dgemm" OFF)
option(CP2K_USE_LIBINT2 "Enable libint2 support" OFF)
option(CP2K_USE_PLUMED "Enable plumed2 support" OFF)
option(CP2K_USE_VORI "Enable libvori support" OFF)
option(CP2K_USE_PEXSI "Enable pexsi support" OFF)
option(CP2K_USE_QUIP "Enable quip support" OFF)
option(CP2K_USE_SPGLIB "Enable spglib support" OFF)
option(CP2K_USE_LIBXC "Enable libxc support" OFF)
option(CP2K_USE_LIBTORCH "Enable libtorch support" OFF)
option(CP2K_USE_STATIC_BLAS "Link against static version of BLAS/LAPACK" OFF)
option(CP2K_USE_SPLA
       "Use SPLA offloading gemm feature to the GPU if it is beneficial. " OFF)
option(CP2K_USE_METIS "enable metis library support" OFF)
option(CP2K_USE_LIBXSMM "Use libxsmm for small gemms (supports x86 platforms)"
       OFF)
option(BUILD_SHARED_LIBS "Build cp2k shared library" ON)
option(
  CP2K_USE_FFTW3_WITH_MKL
  "MKL has its own compatible implementation of the fftw. This option when ON will use the original implementation of the fftw library"
  OFF)
cmake_dependent_option(CP2K_ENABLE_ELPA_OPENMP_SUPPORT
                       "Enable elpa openmp support" ON "CP2K_USE_ELPA" OFF)
cmake_dependent_option(CP2K_ENABLE_FFTW3_OPENMP_SUPPORT
                       "Enable FFTW openmp support" ON "CP2K_USE_FFTW3" OFF)
cmake_dependent_option(CP2K_ENABLE_FFTW3_THREADS_SUPPORT
                       "Enable FFTW THREADS support" OFF "CP2K_USE_FFTW3" OFF)
cmake_dependent_option(CP2K_USE_MPI_F08 "Enable MPI Fortran 2008 interface" OFF
                       "CP2K_USE_MPI" OFF)

cmake_dependent_option(
  CP2K_USE_CUSOLVER_MP
  "Use Nvidia gpu accelerated eigensolver. Only active when CUDA is ON" OFF
  "CP2K_USE_ACCEL MATCHES \"CUDA\"" OFF)

set(CP2K_BLAS_VENDOR
    "auto"
    CACHE STRING "Blas library for computations on host")

set(CP2K_SCALAPACK_VENDOR_LIST "MKL" "SCI" "GENERIC")
set(CP2K_SCALAPACK_VENDOR
    "GENERIC"
    CACHE STRING "scalapack vendor/generic backend")
set_property(CACHE CP2K_SCALAPACK_VENDOR PROPERTY STRINGS
                                                  ${CP2K_SCALAPACK_VENDOR_LIST})

if(NOT ${CP2K_SCALAPACK_VENDOR} IN_LIST CP2K_SCALAPACK_VENDOR_LIST)
  message(FATAL_ERROR "Invalid scalapack vendor backend")
endif()

set(CP2K_BUILD_OPTIONS_LIST "CUSTOM" "DEFAULT" "MINIMAL" "FULL" "SERIAL")

set(CP2K_BUILD_OPTIONS
    "CUSTOM"
    CACHE
      STRING
      "Build cp2k with a predefined set of dependencies. The default setting is full user control"
)
set_property(CACHE CP2K_BUILD_OPTIONS PROPERTY STRINGS
                                               ${CP2K_BUILD_OPTIONS_LIST})

string(TOUPPER ${CP2K_BUILD_OPTIONS} cp2k_build_options_up)

if(NOT ${cp2k_build_options_up} IN_LIST CP2K_BUILD_OPTIONS_LIST)
  message(FATAL_ERROR "Invalid value for cp2k build options")
endif()

if(${cp2k_build_options_up} STREQUAL "MINIMAL")
  set(CP2K_USE_FFTW3 ON)
  set(CP2K_USE_MPI ON)
endif()

if(${cp2k_build_options_up} STREQUAL "FULL")
  set(CP2K_USE_FFTW3 ON)
  set(CP2K_USE_MPI ON)
  set(CP2K_USE_SIRIUS ON)
  set(CP2K_USE_LIBXSMM ON)
  set(CP2K_USE_VORI ON)
  set(CP2K_USE_COSMA ON)
  set(CP2K_USE_SPLA ON)
  set(CP2K_USE_LIBXC ON)
  set(CP2K_USE_LIBINT2 ON)
  set(CP2K_USE_ELPA ON)
  set(CP2K_USE_SPGLIB ON)
  set(CP2K_USE_LIBTORCH ON)
  set(CP2K_USE_METIS ON)
  set(CP2K_USE_QUIP ON)
  set(CP2K_USE_PEXSI ON)
  set(CP2K_USE_SUPERLU ON)
endif()

if(${cp2k_build_options_up} STREQUAL "SERIAL")
  set(CP2K_USE_FFTW3 ON)
  set(CP2K_USE_MPI OFF)
  set(CP2K_USE_SIRIUS ON)
  set(CP2K_USE_LIBXSMM ON)
  set(CP2K_USE_VORI ON)
  set(CP2K_USE_COSMA ON)
  set(CP2K_USE_SPLA ON)
  set(CP2K_USE_LIBXC ON)
  set(CP2K_USE_LIBINT2 ON)
  set(CP2K_USE_ELPA ON)
  set(CP2K_USE_SPGLIB ON)
  set(CP2K_USE_LIBTORCH ON)
  set(CP2K_USE_METIS ON)
  set(CP2K_USE_QUIP ON)
  set(CP2K_USE_PEXSI ON)
  set(CP2K_USE_SUPERLU ON)
endif()

if(${cp2k_build_options_up} STREQUAL "DEFAULT")
  set(CP2K_USE_FFTW3 ON)
  set(CP2K_USE_MPI ON)
  set(CP2K_USE_COSMA ON)
  set(CP2K_USE_LIBXSMM ON)
  set(CP2K_USE_LIBXC ON)
  set(CP2K_USE_LIBINT2 ON)
  set(CP2K_USE_SPGLIB ON)
endif()

# ##############################################################################
# # gpu related options                                                    # #
# ##############################################################################

set(CP2K_SUPPORTED_ACCELERATION_TARGETS CUDA HIP NONE)
set(CP2K_SUPPORTED_CUDA_ARCHITECTURES K20X K40 K80 P100 V100 A100)
set(CP2K_SUPPORTED_HIP_ARCHITECTURES
    Mi50
    Mi100
    Mi210
    Mi250
    K20X
    K40
    K80
    P100
    V100
    A100)

set(CP2K_WITH_GPU
    "P100"
    CACHE STRING
          "Set the CUDA GPU architecture if HIP is enabled (default: P100)")

set_property(
  CACHE CP2K_WITH_GPU PROPERTY STRINGS ${CP2K_SUPPORTED_CUDA_ARCHITECTURES}
                               ${CP2K_SUPPORTED_HIP_ARCHITECTURES})

set(CP2K_USE_ACCEL
    "NONE"
    CACHE STRING "Set hardware acceleartion support: CUDA, HIP")

set_property(CACHE CP2K_USE_ACCEL
             PROPERTY STRINGS ${CP2K_SUPPORTED_ACCELERATION_TARGETS})

cmake_dependent_option(CP2K_USE_NVHPC OFF "Enable Nvidia NVHPC kit"
                       "(NOT CP2K_USE_ACCEL MATCHES \"CUDA\")" OFF)

cmake_dependent_option(
  CP2K_USE_SPLA_GEMM_OFFLOADING ON
  "Enable SPLA dgemm offloading (only valid with gpu support on)"
  "(NOT CP2K_USE_ACCEL MATCHES \"NONE\") AND (CP2K_USE_SPLA)" OFF)

# ##############################################################################
#
# GPU debug options
#
# ##############################################################################

cmake_dependent_option(
  CP2K_DISABLE_GRID_GPU
  OFF
  "disable the hardware accelerated backend for grid related functions. It is only effective when general gpu support is enabled."
  "CP2K_DEBUG_MODE"
  OFF)

cmake_dependent_option(
  CP2K_DISABLE_PW_GPU
  OFF
  "disable the ffts accelerated (mostly GPU) backend. It is only effective when general gpu support is enabled."
  "CP2K_DEBUG_MODE"
  OFF)
cmake_dependent_option(
  CP2K_DISABLE_DBM_GPU
  OFF
  "disable the dbm accelerated (mostly GPU) backend. It is only effective when general gpu support is enabled."
  "CP2K_DEBUG_MODE"
  OFF)
cmake_dependent_option(
  CP2K_DBCSR_CPU_ONLY "Use DBCSR compiled without GPU support." OFF
  "(NOT CP2K_USE_ACCEL MATCHES \"NONE\") AND (CP2K_DEBUG_MODE)" OFF)
cmake_dependent_option(DBCSR_USE_OPENMP "Enable openmp support in DBCSR" ON
                       "CP2K_BUILD_DBCSR" OFF)
# ##############################################################################
# specific variables for the regtests. Binaries will be created with an
# extension               #
# ##############################################################################

set(__cp2k_ext "")
set(__cp2k_cmake_name "")

if(CP2K_USE_MPI)
  set(__cp2k_ext "psmp")
else()
  set(__cp2k_ext "ssmp")
endif()

set(__cp2k_cmake_name "local")

if(CP2K_USE_ACCEL MATCHES "CUDA")
  set(__cp2k_cmake_name "local_cuda")
endif()

if(CP2K_USE_ACCEL MATCHES "HIP")
  set(__cp2k_cmake_name "local_hip")
endif()

if(CP2K_CMAKE_SUFFIX)
  string(APPEND __cp2k_cmake_name "${CP2K_CMAKE_SUFFIX}")
endif()

# we can run the src consistency checks without actually searching for any
# dependencies.

if(CP2K_ENABLE_CONSISTENCY_CHECKS)
  add_subdirectory(src)
  # it is better to simply rm -Rf build but if someone wants to do something
  # like
  #
  # cmake -DCP2K_ENABLE_CONSISTENCY_CHECKS=ON .. cmake ..
  #
  # he/she can

  set(CP2K_ENABLE_CONSISTENCY_CHECKS
      OFF
      CACHE BOOL "" FORCE)
  return()
endif()

# for the time being I change the bin directory to point to
# ../exe/build-cmake-{cuda,hip,cpu} when we want to run the regtests afterwards.
# This solution does not require modifying the regtests scripts at all and cmake
# does not mind this change when we want to install cp2k with a command like
# `make install`
if(CP2K_ENABLE_REGTESTS)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
      "${cp2k_SOURCE_DIR}/exe/${__cp2k_cmake_name}"
      CACHE PATH "Single output directory for building all executables.")
else()
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
      ${cp2k_BINARY_DIR}/bin
      CACHE PATH "Single output directory for building all executables.")

endif()

# Python
#
# this module looks preferably for version 3 of Python. If not found, version 2
# is searched. In CMake 3.15, if a python virtual environment is activated, it
# will search the virtual environment for a python interpreter before searching
# elsewhere in the system. In CMake <3.15, the system is searched before the
# virtual environment.

if(NOT Python_EXECUTABLE)
  # If the python interpreter isn't specified as a command line option, look for
  # it:
  find_package(
    Python
    COMPONENTS Interpreter
    REQUIRED)
endif()

# get the git hash Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
  OUTPUT_VARIABLE CP2K_GIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
  COMMAND hostnamectl --transient
  WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
  OUTPUT_VARIABLE CP2K_HOST_NAME
  OUTPUT_STRIP_TRAILING_WHITESPACE)

add_custom_target(
  AlwaysCheckGit
  COMMAND
    ${CMAKE_COMMAND} -DRUN_CHECK_GIT_VERSION=1
    -Dpre_configure_dir=${pre_configure_dir}
    -Dpost_configure_file=${post_configure_dir}
    -DGIT_HASH_CACHE=${GIT_HASH_CACHE} -P ${CURRENT_LIST_DIR}/CheckGit.cmake
  BYPRODUCTS ${post_configure_file})

# MPI

if(CP2K_USE_MPI)
  get_property(REQUIRED_MPI_COMPONENTS GLOBAL PROPERTY ENABLED_LANGUAGES)
  list(REMOVE_ITEM REQUIRED_MPI_COMPONENTS CUDA) # CUDA does not have an MPI
  # component
  if(NOT CMAKE_CROSSCOMPILING) # when cross compiling, assume the users know
    # what they are doing
    set(MPI_DETERMINE_LIBRARY_VERSION TRUE)
  endif()
  find_package(
    MPI
    COMPONENTS ${REQUIRED_MPI_COMPONENTS}
    REQUIRED)

  if(NOT MPI_Fortran_HAVE_F90_MODULE)
    message(
      FATAL_ERROR
        "\
The listed MPI implementation does not provide the required mpi.mod interface. \
When using the GNU compiler in combination with Intel MPI, please use the \
Intel MPI compiler wrappers. Check the INSTALL.md for more information.")
  endif()
  if("${MPI_Fortran_LIBRARY_VERSION_STRING}" MATCHES "Open MPI v2.1"
     OR "${MPI_Fortran_LIBRARY_VERSION_STRING}" MATCHES "Open MPI v3.1")
    message(
      WARNING
        "RMA with ${MPI_Fortran_LIBRARY_VERSION_STRING} is not supported due to issues with its implementation."
        " Please use a newer version of OpenMPI or switch to MPICH if you plan on using MPI-RMA."
    )
  endif()
else()
  if(CP2K_USE_SIRIUS
     OR CP2K_USE_COSMA
     OR CP2K_USE_SPLA
     OR CP2K_USE_ELPA
     OR CP2K_USE_METIS
     OR CP2K_USE_PLUMED
     OR CP2K_USE_PEXSI
     OR CP2K_USE_QUIP)
    message(
      WARNING
        "SIRIUS, COSMA, SPLA, ELPA, METIS, PLUMED, PEXSI and QUIP require MPI support.\n"
        "Either set -DCP2K_USE_MPI to ON or turn the following dependencies off\n\n"
        "List of dependencies with mpi support\n")
    foreach(
      __libs
      sirius
      COSMA
      SPLA
      ELPA
      METIS
      PLUMED
      PEXSI
      QUIP)
      if(CP2K_USE_${__libs})
        message(WARNING " - ${__libs}")
      endif()
    endforeach()
    message(FATAL_ERROR " -- ")
  endif()
endif()

# BLAS & LAPACK, PkgConfig
find_package(Lapack REQUIRED) # also calls find_package(BLAS)

# SMM (Small Matrix-Matrix multiplication)
if(CP2K_USE_LIBXSMM)
  find_package(LibXSMM REQUIRED)
  message(STATUS "Using libxsmm for Small Matrix Multiplication")
endif()

# in practice it is always for any decent configuration. But I add a flags to
# turn it off
if(CP2K_USE_MPI)
  find_package(SCALAPACK REQUIRED)
endif()

# CUDA / ROCM easy for cuda a moving target for hip

if((CP2K_USE_ACCEL MATCHES CUDA) OR (CP2K_USE_ACCEL MATCHES HIP))
  set(CP2K_GPU_ARCH_NUMBER_K20X 35)
  set(CP2K_GPU_ARCH_NUMBER_K40 35)
  set(CP2K_GPU_ARCH_NUMBER_K80 37)
  set(CP2K_GPU_ARCH_NUMBER_P100 60)
  set(CP2K_GPU_ARCH_NUMBER_V100 70)
  set(CP2K_GPU_ARCH_NUMBER_A100 80)
  set(CP2K_GPU_ARCH_NUMBER_Mi50 gfx906)
  set(CP2K_GPU_ARCH_NUMBER_Mi100 gfx908)
  set(CP2K_GPU_ARCH_NUMBER_Mi200 gfx90a)
  set(CP2K_GPU_ARCH_NUMBER_Mi250 gfx90a)
endif()

set(CP2K_USE_HIP OFF)
set(CP2K_USE_CUDA OFF)

if(CP2K_USE_ACCEL MATCHES "CUDA")
  # P100 is the default target.
  set(CMAKE_CUDA_ARCHITECTURES 60)

  # allow for unsupported compilers (gcc/cuda version mismatch)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -allow-unsupported-compiler")

  enable_language(CUDA)
  if(CP2K_USE_NVHPC)
    find_package(NVHPC REQUIRED COMPONENTS CUDA MATH HOSTUTILS NCCL)
  else()
    find_package(CUDAToolkit REQUIRED)
  endif()

  list(FIND CP2K_SUPPORTED_CUDA_ARCHITECTURES ${CP2K_WITH_GPU}
       CP2K_GPU_SUPPORTED)

  if(CP2K_GPU_SUPPORTED EQUAL -1)
    message(
      FATAL_ERROR
        "GPU architecture (${CP2K_WITH_GPU}) is not supported. Please choose from: ${CP2K_SUPPORTED_CUDA_ARCHITECTURES}"
    )
  endif()

  set(CMAKE_CUDA_ARCHITECTURES ${CP2K_GPU_ARCH_NUMBER_${CP2K_WITH_GPU}})

  message(STATUS "GPU target architecture: ${CP2K_WITH_GPU}\n"
                 "GPU architecture number: ${CMAKE_CUDA_ARCHITECTURES}\n"
                 "GPU profiling enabled: ${CP2K_WITH_CUDA_PROFILING}\n")

  if(WITH_CUDA_PROFILING)
    find_library(
      CUDA_NVTOOLSEXT nvToolsExt
      PATHS ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES}
      DOC "Building with CUDA profiling requires the nvToolsExt CUDA library"
          REQUIRED)
    message(STATUS "Found nvToolsExt: ${CUDA_NVTOOLSEXT}")
  endif()

  set(CP2K_USE_CUDA ON)

  if(CP2K_USE_CUSOLVER_MP)
    find_package(CuSolverMP REQUIRED)
  endif()
  message(STATUS ``"-- CUDA compiler and libraries found")
elseif(CP2K_USE_ACCEL MATCHES "HIP")
  enable_language(HIP)
  # Find hip
  find_package(hipfft REQUIRED IMPORTED CONFIG)
  find_package(hipblas REQUIRED IMPORTED CONFIG)

  set(CMAKE_HIP_ARCHITECTURES gfx801 gfx900 gfx90a)
  if(NOT CMAKE_BUILD_TYPE)
    set(HIP_RELEASE_OPTIONS "-O3 -g -DNDEBUG --std=c++11")
  elseif(${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
    set(HIP_RELEASE_OPTIONS "-O3 -g -DNDEBUG --std=c++11")
  elseif(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    set(HIP_RELEASE_OPTIONS "-O3 -DNDEBUG --std=c++11")
  elseif(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(HIP_RELEASE_OPTIONS "-O0 -g --std=c++11")
  endif()

  set(CMAKE_HIP_ARCGITECTURES "${CP2K_GPU_ARCH_NUMBER_${CP2K_WITH_GPU}}")
  set(CP2K_USE_HIP ON)
endif()

# PACKAGE DISCOVERY (compiler configuration can impact package discovery)
find_package(OpenMP REQUIRED COMPONENTS Fortran C CXX)

find_package(DBCSR 2.6 REQUIRED)

# ==================================
if(CP2K_USE_ELPA)
  find_package(Elpa REQUIRED)
endif()

if(CP2K_USE_LIBXC)
  find_package(LibXC 6 REQUIRED EXACT)
endif()

if(CP2K_USE_COSMA)
  find_package(cosma REQUIRED)

  # check that cosma::cosma_pxgemm and cosma::cosma_prefixed_pxgemm exist
  if(NOT TARGET cosma::cosma_pxgemm OR NOT TARGET cosma::cosma_prefixed_pxgemm)
    message(
      FATAL_ERROR
        " COSMA needs to be build with scalapack offloading support. COSTA_SCALAPACK and COSMA_SCALAPACK should probably be set properly"
    )
  endif()
endif()

if(CP2K_USE_VORI)
  find_package(LibVORI REQUIRED)
endif()

# FFTW3

# we set this variable to ON when we want fftw3 support (with or without MKL).
set(CP2K_USE_FFTW3_ OFF)
if(CP2K_USE_FFTW3)
  if(NOT CP2K_BLAS_VENDOR MATCHES "MKL" OR CP2K_USE_FFTW3_WITH_MKL)
    find_package(Fftw REQUIRED)
    if(CP2K_ENABLE_FFTW3_THREADS_SUPPORT AND CP2K_ENABLE_FFTW3_OPENMP_SUPPORT)
      message(
        FATAL_ERROR
          "Fftw3 threads and openmp supports can not be used at the same time")
    endif()

    if((CP2K_ENABLE_FFTW3_THREADS_SUPPORT) AND (NOT TARGET
                                                cp2k::FFTW3::fftw3_threads))
      message(
        FATAL_ERROR
          "fftw3 was compiled without multithreading support (--enable-threads option in the fftw build system)."
      )
    endif()

    if((CP2K_ENABLE_FFTW3_OPENMP_SUPPORT) AND (NOT TARGET cp2k::FFTW3::fftw3_omp
                                              ))
      message(
        FATAL_ERROR
          "fftw3 was compiled without openmp support  (--enable-openmp option in the fftw build system)."
      )
    endif()
    set(CP2K_USE_FFTW3_ ON)
  else()
    message("-- Using the MKL implementation of FFTW3.")
  endif()
endif()

# QUIP
if(CP2K_USE_QUIP)
  find_package(Quip REQUIRED)
endif()

# libint

if(CP2K_USE_LIBINT2)
  find_package(Libint2 REQUIRED)
endif()

# spglib
if(CP2K_USE_SPGLIB)
  find_package(LibSPG REQUIRED)
endif()

if(CP2K_USE_SPLA)
  find_package(SPLA REQUIRED)
  get_target_property(SPLA_INCLUDE_DIRS SPLA::spla
                      INTERFACE_INCLUDE_DIRECTORIES)
  if(NOT SPLA_INCLUDE_DIRS)
    set(SPLA_INCLUDE_DIRS "/usr/include;/usr/include/spla")
  endif()

  if(NOT SPLA_GPU_BACKEND AND CP2K_USE_GEMM_OFFLOADING)
    set(CP2K_USE_GEMM_OFFLOADING OFF)
    message(
      FATAL_ERROR
        "SPLA should be compiled with GPU support if the gemm offloading is requested. Use -DCP2K_USE_GEMM_OFFLOADING=OFF otherwise"
    )
  endif()
endif()
# SIRIUS

if(CP2K_USE_SIRIUS)
  find_package(sirius REQUIRED)
endif()

if(CP2K_USE_SUPERLU)
  find_package(SuperLU REQUIRED)
endif()

if(CP2K_USE_METIS)
  find_package(Metis)
endif()

if(CP2K_USE_PTSCOTCH)
  find_package(Ptscotch REQUIRED)
endif()

if(CP2K_USE_PEXSI)
  # PEXSI 1.2 uses cmake as build system
  find_package(PEXSI REQUIRED)
endif()

if(CP2K_USE_PLUMED)
  find_package(Plumed REQUIRED)
endif()

if(CP2K_USE_LIBTORCH)
  find_package(Torch REQUIRED)
endif()

if(CP2K_USE_MPI_F08 AND NOT MPI_Fortran_HAVE_F08_MODULE)
  message(
    FATAL_ERROR
      "The Fortran 2008 interface is not supported by the MPI implementation found by cmake."
  )
endif()

# OPTION HANDLING

# make sure that the default build type is RELEASE
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(
    STATUS
      "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "${default_build_type}"
      CACHE STRING
            "Choose the type of build, options are: Debug Release Coverage."
            FORCE)
  # set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "Coverage")
endif()

# compiler configuration could have impacted package discovery (above)
include(CompilerConfiguration)
include(CheckCompilerSupport)

if(CP2K_BUILD_DBCSR)
  add_subdirectory(exts/dbcsr/src)
  add_library(DBCSR::dbcsr ALIAS dbcsr)
endif()

# subdirectories
add_subdirectory(src)

include(GNUInstallDirs)

get_target_property(CP2K_LIBS cp2k_link_libs INTERFACE_LINK_LIBRARIES)
configure_file(cmake/libcp2k.pc.in libcp2k.pc @ONLY)

message(
  ""
  "--------------------------------------------------------------------\n"
  "-                                                                  -\n"
  "-               Summary of enabled dependencies                    -\n"
  "-                                                                  -\n"
  "--------------------------------------------------------------------\n\n")

message(
  "  - BLAS AND LAPACK\n\n"
  "   - vendor :              ${CP2K_BLAS_VENDOR}\n"
  "   - include directories : ${CP2K_BLAS_INCLUDE_DIR} ${LAPACK_INCLUDE_DIR}\n"
  "   - libraries :           ${CP2K_BLAS_LINK_LIBRARIES} ${CP2K_LAPACK_LINK_LIBRARIES}\n\n"
)

if(CP2K_USE_MPI)
  message("  - MPI\n" "   - include directories :  ${MPI_INCLUDE_DIRS}\n"
          "   - libraries :           ${MPI_LIBRARIES}\n\n")

  if(CP2K_USE_MPI_F08)
    message("   - MPI_08 :              ON\n")
  endif()

  if(MPI_Fortran_HAVE_F08_MODULE AND NOT CP2K_USE_MPI_F08)
    message(
      "  - MPI_08 is supposed by MPI but turned off by default.\n"
      "    To use it add -DCP2K_USE_MPI_F08=ON to the cmake command line\n\n")
  endif()

  message("  - SCALAPACK:\n"
          "    - libraries : ${CP2K_SCALAPACK_LINK_LIBRARIES}\n\n")
endif()

if((CP2K_USE_ACCEL MATCHES "CUDA") OR (CP2K_USE_ACCEL MATCHES "HIP"))

  message("  - Hardware Acceleration:\n")
  if(CP2K_USE_ACCEL MATCHES "CUDA")
    message(
      "   - CUDA:\n" "     - GPU target architecture : ${CP2K_WITH_GPU}\n"
      "     - GPU architecture number : ${CP2K_ACC_ARCH_NUMBER}"
      "     - GPU profiling enabled :   ${CP2K_WITH_CUDA_PROFILING}\n\n")
  endif()

  if(CP2K_USE_ACCEL MATCHES "HIP")
    message("   - HIP:\n" "    - GPU target architecture : ${CP2K_WITH_GPU}\n"
            "    - GPU architecture number : ${CP2K_ACC_ARCH_NUMBER}\n\n")
  endif()

  message(
    "Note : Enabling hardware acceleration enable acceleration of the grid, pw, and dbm modules by default\n"
    "    - GRID module : ${CP2K_USE_GRID_GPU}\n"
    "    - PW   module : ${CP2K_USE_PW_GPU}\n"
    "    - DBM  module : ${CP2K_USE_DBM_GPU}\n\n")
endif()

if(CP2K_USE_LIBXC)
  message(
    "  - LIBXC (note to package managers : libxc can be build with cmake as well)"
    "   - include directories : ${CP2K_LIBXC_INCLUDE_DIRS}\n"
    "   - libraries : ${CP2K_LIBXC_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_LIBTORCH)
  message("  - LIBTORCH\n" "    - libraries : ${CP2K_LIBTORCH_LIBRARIES}\n\n")
endif()

if(CP2K_USE_FFTW3)
  message("  - FFTW3\n"
          "    - include directories : ${CP2K_FFTW3_INCLUDE_DIRS}\n"
          "    - libraries : ${CP2K_FFTW3_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_LIBXSMM)
  message(
    "  - libxsmm\n"
    "    - include directories : ${CP2K_LIBXSMM_INCLUDE_DIRS}\n"
    "    - libraries :           ${CP2K_LIBXSMMEXT_LINK_LIBRARIES};${CP2K_LIBXSMMF_LINK_LIBRARIES}\n\n"
  )
endif()

if(CP2K_USE_SPLA)
  message(" - SPLA\n" "   - include directories : ${SPLA_INCLUDE_DIRS}\n"
          "   - lbraries : ${SPLA_LIBRARIES}\n\n")
endif()

if(CP2K_USE_SIRIUS)
  message(
    " - SIRIUS :\n"
    "   - include directories :  ${SIRIUS_INCLUDE_DIRS}\n"
    "   - libraries           :  ${SIRIUS_LIBRARIES}\n"
    "   - dependencies :\n"
    "       - spla\n"
    "       - SpFFT\n"
    "       - SPGLIB\n"
    "       - LibXC\n"
    "       - fftw3\n"
    "       - hdf5\n"
    "       - GSL\n\n")

  if(CP2K_USE_VDWXC)
    message("     - VDWXC\n")
  endif()
endif()

if(CP2K_USE_COSMA)
  message(" - COSMA\n" "   - include directories : ${COSMA_INCLUDE_DIRS}\n"
          "   - libraries           : ${COSMA_LIBRARIES}\n\n")
endif()

if(CP2K_USE_QUIP)
  message(" - QUIP\n"
          "   - include directories : ${CP2K_LIBQUIP_INCLUDE_DIRS}\n"
          "   - libraries :           ${CP2K_LIBQUIP_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_PEXSI)
  message(" - PEXSI\n\n")
endif()

if(CP2K_USE_LIBINT2)
  message(" - libint2\n"
          "   - include directories : ${CP2K_LIBINT2_INCLUDE_DIRS}\n"
          "   - libraries :           ${CP2K_LIBINT2_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_VORI)
  message(" - libvori\n"
          "   - include directories : ${CP2K_LIBVORI_INCLUDE_DIRS}\n"
          "   - libraries :           ${CP2K_LIBVORI_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_SPGLIB)
  message(" - spglib\n"
          "   - include directories : ${CP2K_LIBSPG_INCLUDE_DIRS}\n"
          "   - libraries           : ${CP2K_LIBSPG_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_ELPA)
  message(" - ELPA\n" "   - include directories : ${CP2K_ELPA_INCLUDE_DIRS}\n"
          "   - libraries           :  ${CP2K_ELPA_LINK_LIBRARIES})\n\n")
endif()

if(CP2K_USE_SUPERLU)
  message(" - superlu\n"
          "   - include directories : ${CP2K_SUPERLU_INCLUDE_DIRS}\n"
          "   - libraries           : ${CP2K_SUPERLU_LINK_LIBRARIES}\n\n")
endif()

message(
  "--------------------------------------------------------------------\n"
  "-                                                                  -\n"
  "-        List of dependencies not included in this build           -\n"
  "-                                                                  -\n"
  "--------------------------------------------------------------------\n")

if(NOT CP2K_USE_MPI)
  message("   - MPI")
endif()

if(NOT CP2K_USE_SIRIUS)
  message("   - SIRIUS")
endif()

if(NOT CP2K_USE_SPGLIB)
  message("   - SPGLIB")
endif()

if(NOT CP2K_USE_COSMA)
  message("   - COSMA")
endif()

if(NOT CP2K_USE_SPLA)
  message("   - SPLA")
endif()

if(${CP2K_USE_ACCEL} MATCHES "NONE")
  message("   - GPU acceleration is disabled")
endif()

if(NOT CP2K_USE_ELPA)
  message("   - ELPA")
endif()

if(NOT CP2K_USE_PLUMMED)
  message("   - PLUMED")
endif()

if(NOT CP2K_USE_QUIP)
  message("   - QUIP")
endif()

if(NOT CP2K_USE_LIBXSMM)
  message("   - LIBXSMM")
endif()

if(NOT CP2K_USE_LIBINT2)
  message("   - LIBINT2")
endif()

if(NOT CP2K_USE_LIBXC)
  message("   - LIBXC")
endif()

if(NOT CP2K_USE_VORI)
  message("   - LIBVORI")
endif()

if(NOT CP2K_USE_FFTW3)
  message("   - FFTW3")
endif()

if(NOT CP2K_USE_PEXSI)
  message("   - PEXSI")
endif()

if(NOT CP2K_USE_SUPERLU)
  message("   - SUPERLU")
endif()

if(NOT CP2K_USE_LIBTORCH)
  message("   - libtorch")
endif()

if(CP2K_ENABLE_REGTESTS OR cp2k_TESTS)
  message("\n\n" "To run the regtests you need to run the following commands\n"
          "\n\n cd ..\n" " export CP2K_DATA_DIR=${CMAKE_SOURCE_DIR}/data/\n"
          " ./tests/do_regtest.py ${__cp2k_cmake_name} ${__cp2k_ext}\n\n")
endif()

# files needed for cmake

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/cp2kConfigVersion.cmake"
  VERSION "${CP2K_VERSION}"
  COMPATIBILITY SameMajorVersion)

configure_file("${PROJECT_SOURCE_DIR}/cmake/cp2kConfig.cmake.in"
               "${PROJECT_BINARY_DIR}/cp2kConfig.cmake" @ONLY)

install(FILES "${PROJECT_BINARY_DIR}/cp2kConfig.cmake"
              "${PROJECT_BINARY_DIR}/cp2kConfigVersion.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/cp2k")

install(FILES "${PROJECT_BINARY_DIR}/libcp2k.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

install(
  DIRECTORY "${PROJECT_SOURCE_DIR}/cmake/modules"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/cp2k"
  FILES_MATCHING
  PATTERN "*.cmake")
