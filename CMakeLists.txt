#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

cmake_minimum_required(VERSION 3.24)

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

# !!! Keep version in sync with cp2k_info.F !!!
project(
  cp2k
  DESCRIPTION "CP2K"
  HOMEPAGE_URL "https://www.cp2k.org"
  VERSION "2025.2"
  LANGUAGES Fortran C CXX)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")

# set language and standard.
#
# cmake does not provide any mechanism to set the fortran standard. Adding the
# `std` option to compiler flags is the only way to control it. So leave them
# be.
#

if(NOT DEFINED CMAKE_CUDA_STANDARD)
  set(CMAKE_CUDA_STANDARD 14)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

if(NOT DEFINED CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 11)
  set(CMAKE_C_STANDARD_REQUIRED ON)
endif()

if(NOT DEFINED CMAKE_HIP_STANDARD)
  set(CMAKE_HIP_STANDARD 14)
  set(CMAKE_HIP_STANDARD_REQUIRED ON)
endif()

# remove NDEBUG flag
string(REPLACE "-DNDEBUG" "" CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")
string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELEASE
               "${CMAKE_CXX_FLAGS_RELEASE}")
string(REPLACE "-DNDEBUG" "" CMAKE_Fortran_FLAGS_RELEASE
               "${CMAKE_Fortran_FLAGS_RELEASE}")
string(REPLACE "-DNDEBUG" "" CMAKE_C_FLAGS_RELWITHDEBINFO
               "${CMAKE_C_FLAGS_RELWITHDEBINFO}")
string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_RELWITHDEBINFO
               "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
string(REPLACE "-DNDEBUG" "" CMAKE_Fortran_FLAGS_RELWITHDEBINFO
               "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
string(REPLACE "-DNDEBUG" "" CMAKE_C_FLAGS_MINSIZEREL
               "${CMAKE_C_FLAGS_MINSIZEREL}")
string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_MINSIZEREL
               "${CMAKE_CXX_FLAGS_MINSIZEREL}")
string(REPLACE "-DNDEBUG" "" CMAKE_Fortran_FLAGS_MINSIZEREL
               "${CMAKE_Fortran_FLAGS_MINSIZEREL}")

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

set(CMAKE_INSTALL_LIBDIR
    "lib"
    CACHE PATH "Default installation directory for libraries")

# =================================================================================================
# OPTIONS
option(CMAKE_POSITION_INDEPENDENT_CODE "Enable position independent code" ON)

option(CP2K_ENABLE_CONSISTENCY_CHECKS
       "Check that the list of compiled files and files contained in src match"
       OFF)

option(CP2K_USE_EVERYTHING
       "Enable all dependencies. They can be individually turned off again."
       OFF)

option(CP2K_USE_ACE "Enable ACE support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_DEEPMD "Enable DeePMD support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_DFTD4 "Enable DFTD4 support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_TBLITE "Enable TBLITE support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_FFTW3 "Enable FFTW3 support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_GREENX "Enable GreenX support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_GRPP "Enable GRPP support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_HDF5 "Enable HDF5 support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_LIBINT2 "Enable Libint2 support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_LIBTORCH "Enable LibTorch support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_LIBXC "Enable LibXC support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_MPI "Enable MPI support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_PEXSI "Enable PEXSI support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_SPGLIB "Enable Spglib support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_TREXIO "Enable TREXIO support" ${CP2K_USE_EVERYTHING})
option(CP2K_USE_VORI "Enable Libvori support" ${CP2K_USE_EVERYTHING})

option(CP2K_USE_STATIC_BLAS "Link against static version of BLAS/LAPACK" OFF)
option(BUILD_SHARED_LIBS "Build CP2K shared library" ON)
option(
  CP2K_USE_FFTW3_WITH_MKL
  "MKL has its own compatible implementation of the FFTW3 library. This option, when ON, will use the separate and original FFTW3 library."
  OFF)

# MPI-enabled options
cmake_dependent_option(CP2K_USE_COSMA "Enable COSMA support"
                       ${CP2K_USE_EVERYTHING} "CP2K_USE_MPI" OFF)
cmake_dependent_option(CP2K_USE_DLAF "Enable DLA-Future support"
                       ${CP2K_USE_EVERYTHING} "CP2K_USE_MPI" OFF)
cmake_dependent_option(CP2K_USE_ELPA "Enable ELPA support"
                       ${CP2K_USE_EVERYTHING} "CP2K_USE_MPI" OFF)
cmake_dependent_option(CP2K_USE_LIBSMEAGOL "Enable LibSMEAGOL support"
                       ${CP2K_USE_EVERYTHING} "CP2K_USE_MPI" OFF)
cmake_dependent_option(CP2K_USE_MPI_F08 "Enable MPI Fortran 2008 interface"
                       ${CP2K_USE_EVERYTHING} "CP2K_USE_MPI" OFF)
cmake_dependent_option(CP2K_USE_PEXSI "Enable PEXSI support"
                       ${CP2K_USE_EVERYTHING} "CP2K_USE_MPI" OFF)
cmake_dependent_option(CP2K_USE_PLUMED "Enable PLUMED2 support"
                       ${CP2K_USE_EVERYTHING} "CP2K_USE_MPI" OFF)
cmake_dependent_option(CP2K_USE_SIRIUS "Enable SIRIUS support"
                       ${CP2K_USE_EVERYTHING} "CP2K_USE_MPI" OFF)
cmake_dependent_option(CP2K_USE_SPLA "Enable SpLA support"
                       ${CP2K_USE_EVERYTHING} "CP2K_USE_MPI" OFF)
cmake_dependent_option(
  CP2K_USE_LIBXSMM "Enable libxsmm support" ${CP2K_USE_EVERYTHING}
  "NOT CP2K_USE_ACCEL MATCHES \"OPENCL\"" ON)
cmake_dependent_option(CP2K_USE_LIBVDWXC "Enable libvdwxc support with SIRIUS"
                       ${CP2K_USE_EVERYTHING} "CP2K_USE_SIRIUS" OFF)

# do not turn it on if sirius is built with the toolchain
cmake_dependent_option(
  CP2K_USE_SIRIUS_DFTD4 "Enable dftd3 and dftd4 sirius support" OFF
  "CP2K_USE_SIRIUS" OFF)
cmake_dependent_option(
  CP2K_USE_SIRIUS_VCSQNM "Enable variable cell shape relaxation in SIRIUS" OFF
  "CP2K_USE_SIRIUS" OFF)

cmake_dependent_option(CP2K_USE_SIRIUS_NLCG "Enable nlcg method in SIRIUS" OFF
                       "CP2K_USE_SIRIUS" OFF)

cmake_dependent_option(
  CP2K_DBCSR_USE_CPU_ONLY "Disable the DBCSR accelerated backend" OFF
  "NOT CP2K_USE_ACCEL MATCHES \"OPENCL\"" OFF)

cmake_dependent_option(
  CP2K_ENABLE_DBM_GPU "Disable the dbm accelerated backend (mostly GPU)." ON
  "CP2K_USE_ACCEL" OFF)

cmake_dependent_option(
  CP2K_ENABLE_GRID_GPU "Disable acceleration for grid related functions." ON
  "CP2K_USE_ACCEL MATCHES \"HIP|CUDA\"" OFF)

cmake_dependent_option(
  CP2K_ENABLE_PW_GPU "Disable the ffts accelerated backend (mostly GPU)." ON
  "CP2K_USE_ACCEL MATCHES \"HIP|CUDA\"" OFF)

cmake_dependent_option(
  CP2K_USE_UNIFIED_MEMORY "Use CPU/GPU unified memory (Mi250x onwards)" OFF
  "CP2K_USE_ACCEL MATCHES \"HIP\"" OFF)

cmake_dependent_option(CP2K_ENABLE_ELPA_OPENMP_SUPPORT
                       "Enable ELPA OpenMP support" ON "CP2K_USE_ELPA" OFF)
cmake_dependent_option(CP2K_ENABLE_FFTW3_OPENMP_SUPPORT
                       "Enable FFTW3 OpenMP support" ON "CP2K_USE_FFTW3" OFF)
cmake_dependent_option(CP2K_ENABLE_FFTW3_THREADS_SUPPORT
                       "Enable FFTW3 threads support" OFF "CP2K_USE_FFTW3" OFF)

cmake_dependent_option(
  CP2K_USE_CUSOLVER_MP
  "Use Nvidia GPU accelerated eigensolver. Only active when CUDA is ON" OFF
  "CP2K_USE_ACCEL MATCHES \"CUDA\"" OFF)

cmake_dependent_option(CP2K_USE_NVHPC OFF "Enable Nvidia NVHPC kit"
                       "(NOT CP2K_USE_ACCEL MATCHES \"CUDA\")" OFF)

cmake_dependent_option(
  CP2K_USE_SPLA_GEMM_OFFLOADING ON
  "Enable SpLA dgemm offloading (only valid with GPU support on)"
  "(NOT CP2K_USE_ACCEL MATCHES \"NONE\") AND (CP2K_USE_SPLA)" OFF)

cmake_dependent_option(
  CP2K_USE_LIBVDWXC OFF
  "Compile CP2K with libvdwxc support when SIRIUS is compiled with it"
  "CP2K_USE_SIRIUS" OFF)

set(CP2K_BLAS_VENDOR
    "auto"
    CACHE STRING "BLAS library for computations on host")

set(CP2K_SCALAPACK_VENDOR_LIST "MKL" "SCI" "GENERIC" "auto")
set(CP2K_SCALAPACK_VENDOR
    "auto"
    CACHE STRING "ScaLAPACK vendor/generic backend")
set_property(CACHE CP2K_SCALAPACK_VENDOR PROPERTY STRINGS
                                                  ${CP2K_SCALAPACK_VENDOR_LIST})

if(DEFINED CP2K_SCALAPACK_VENDOR)
  if(NOT ${CP2K_SCALAPACK_VENDOR} IN_LIST CP2K_SCALAPACK_VENDOR_LIST)
    message(FATAL_ERROR "An invalid ScaLAPACK vendor backend was specified")
  endif()
endif()

set(CP2K_DATA_DIR
    "default"
    CACHE STRING "Set the location for CP2K data")

# ##############################################################################
# GPU related options
# ##############################################################################

set(CP2K_SUPPORTED_ACCELERATION_TARGETS CUDA HIP OPENCL NONE)
set(CP2K_SUPPORTED_CUDA_ARCHITECTURES
    K20X
    K40
    K80
    P100
    V100
    A100
    H100
    A40)
set(CP2K_SUPPORTED_HIP_ARCHITECTURES
    Mi50
    Mi100
    Mi210
    Mi250
    Mi300
    K20X
    K40
    K80
    P100
    V100
    A100
    H100
    A40)

set(CP2K_WITH_GPU
    "NONE"
    CACHE STRING
          "Set the CUDA GPU architecture if HIP is enabled (default: NONE)")

set_property(
  CACHE CP2K_WITH_GPU PROPERTY STRINGS ${CP2K_SUPPORTED_CUDA_ARCHITECTURES}
                               ${CP2K_SUPPORTED_HIP_ARCHITECTURES})

set(CP2K_USE_ACCEL
    "NONE"
    CACHE STRING "Set hardware acceleration support: CUDA, HIP, OPENCL")

set_property(CACHE CP2K_USE_ACCEL
             PROPERTY STRINGS ${CP2K_SUPPORTED_ACCELERATION_TARGETS})

# ##############################################################################
# specific variables for the regtests. Binaries will be created with an
# extension
# ##############################################################################

set(__cp2k_ext "")

if(CP2K_USE_MPI)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(__cp2k_ext "pdbg")
  else()
    set(__cp2k_ext "psmp")
  endif()
else()
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(__cp2k_ext "sdbg")
  else()
    set(__cp2k_ext "ssmp")
  endif()
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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${cp2k_BINARY_DIR}/bin
    CACHE PATH "Single output directory for building all executables.")

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
  get_target_property(MPI_INCLUDE_DIRS MPI::MPI_Fortran
                      INTERFACE_INCLUDE_DIRECTORIES)
  if(NOT MPI_Fortran_HAVE_F90_MODULE)
    message(
      FATAL_ERROR
        "The listed MPI implementation does not provide the required mpi.mod interface."
        " When using the GNU compiler in combination with Intel MPI, please use the"
        " Intel MPI compiler wrappers. Check the INSTALL.md for more information."
    )
  endif()
  if("${MPI_Fortran_LIBRARY_VERSION_STRING}" MATCHES "Open MPI v2.1"
     OR "${MPI_Fortran_LIBRARY_VERSION_STRING}" MATCHES "Open MPI v3.1")
    message(
      WARNING
        "RMA with ${MPI_Fortran_LIBRARY_VERSION_STRING} is not supported due to issues with its implementation."
        " Please use a newer version of OpenMPI or switch to MPICH if you plan on using MPI-RMA."
    )
  endif()
endif()

# BLAS & LAPACK, PkgConfig
find_package(Lapack REQUIRED) # also calls find_package(BLAS)

# SMM (Small Matrix-Matrix multiplication)
if(CP2K_USE_LIBXSMM)
  find_package(LibXSMM REQUIRED)
  message(STATUS "Using LIBXSMM for Small Matrix Multiplication")
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
  set(CP2K_GPU_ARCH_NUMBER_H100 90)
  set(CP2K_GPU_ARCH_NUMBER_A40 86)
  set(CP2K_GPU_ARCH_NUMBER_Mi50 gfx906)
  set(CP2K_GPU_ARCH_NUMBER_Mi100 gfx908)
  set(CP2K_GPU_ARCH_NUMBER_Mi200 gfx90a)
  set(CP2K_GPU_ARCH_NUMBER_Mi250 gfx90a)
  set(CP2K_GPU_ARCH_NUMBER_Mi300 gfx942)

  # CMAKE_HIP_ARCHITECTURES and CMAKE_CUDA_ARCHITECTURES are the prefered
  # mechanism to set the GPU architecture. We still offer the CP2K_WITH_GPU
  # option to avoid breaking the ci/cd or any other scripts based on this
  # option.

  # check that CMAKE_{HIP|CUDA}_ARCHITECTURES or CP2K_WITH_GPU are given
  if((NOT DEFINED CMAKE_HIP_ARCHITECTURES)
     AND (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
     AND (CP2K_WITH_GPU MATCHES NONE))
    message(
      FATAL_ERROR
        "----------------------------------------------------------------\n"
        "                                                                \n"
        "CMAKE_{HIP|CUDA}_ARCHITECTURES or CP2K_WITH_GPU should be given.\n"
        "                                                                \n"
        "----------------------------------------------------------------\n")
  endif()

  set(_ignore_with_gpu_option false)
  if((DEFINED CMAKE_HIP_ARCHITECTURES) OR (DEFINED CMAKE_CUDA_ARCHITECTURES))
    set(_ignore_with_gpu_option true)
  endif()

  if(NOT _ignore_with_gpu_option)
    message(
      STATUS
        "\nCP2K_WITH_GPU is deprecated in favor of CMAKE_HIP_ARCHITECTURES or CMAKE_CUDA_ARCHITECTURES\n"
    )
    if(CP2K_USE_ACCEL MATCHES CUDA)
      list(FIND CP2K_SUPPORTED_CUDA_ARCHITECTURES ${CP2K_WITH_GPU}
           CP2K_GPU_SUPPORTED)

      if(CP2K_GPU_SUPPORTED EQUAL -1)
        message(
          FATAL_ERROR
            "GPU architecture (${CP2K_WITH_GPU}) is not supported. Please choose from: ${CP2K_SUPPORTED_CUDA_ARCHITECTURES}"
        )
      endif()
      set(CMAKE_CUDA_ARCHITECTURES ${CP2K_GPU_ARCH_NUMBER_${CP2K_WITH_GPU}})
    else(CP2K_USE_ACCEL MATCHES HIP)
      list(FIND CP2K_SUPPORTED_HIP_ARCHITECTURES ${CP2K_WITH_GPU}
           CP2K_GPU_SUPPORTED)
      if(CP2K_GPU_SUPPORTED EQUAL -1)
        message(
          FATAL_ERROR
            "GPU architecture (${CP2K_WITH_GPU}) is not supported. Please choose from: ${CP2K_SUPPORTED_HIP_ARCHITECTURES}"
        )
      endif()
      set(CMAKE_HIP_ARCHITECTURES "${CP2K_GPU_ARCH_NUMBER_${CP2K_WITH_GPU}}")
    endif()
  endif()
endif()

set(CP2K_USE_HIP OFF)
set(CP2K_USE_CUDA OFF)
set(CP2K_USE_OPENCL OFF)

if(CP2K_USE_ACCEL MATCHES "CUDA")
  option(CP2K_WITH_CUDA_PROFILING "Enable CUDA profiling" OFF)
  # P100 is the default target.

  # allow for unsupported compilers (gcc/cuda version mismatch)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -allow-unsupported-compiler")

  enable_language(CUDA)
  if(CP2K_USE_NVHPC)
    find_package(NVHPC REQUIRED COMPONENTS CUDA MATH HOSTUTILS NCCL)
  else()
    find_package(CUDAToolkit REQUIRED)
  endif()

  message("\n-----------------------------------------------------------")
  message("-                          CUDA                           -")
  message("-----------------------------------------------------------\n")

  message(STATUS "GPU architecture number: ${CMAKE_CUDA_ARCHITECTURES}")
  message(STATUS "GPU profiling enabled: ${CP2K_WITH_CUDA_PROFILING}")

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

  message(STATUS "CUDA compiler and libraries found\n")

elseif(CP2K_USE_ACCEL MATCHES "HIP")
  message("\n------------------------------------------------------------")
  message("-                          HIP                             -")
  message("------------------------------------------------------------\n")
  message(INFO "${CMAKE_HIP_ARCHITECTURES}")
  enable_language(HIP)

  if(CMAKE_HIP_PLATFORM MATCHES "nvidia")
    find_package(CUDAToolkit)
  endif()

  if(NOT CMAKE_BUILD_TYPE AND (CMAKE_HIP_PLATFORM MATCHES "amd"))
    set(CMAKE_HIP_FLAGS "-O3")
  elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_HIP_FLAGS "-O2 -g")
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_HIP_FLAGS "-O3")
  elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_HIP_FLAGS "-O0 -g")
  endif()

  # Find HIP
  find_package(hipfft REQUIRED IMPORTED CONFIG)
  find_package(hipblas REQUIRED IMPORTED CONFIG)

  set(CP2K_USE_HIP ON)

  # use hardware atomic operations on Mi250X.
  if(NOT CMAKE_HIP_PLATFORM OR (CMAKE_HIP_PLATFORM MATCHES "amd"))
    set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -munsafe-fp-atomics")
  endif()

  # add the Mi300A parameters when available
  if(CP2K_USE_UNIFIED_MEMORY)
    if(CMAKE_HIP_ARCHITECTURES MATCHES "gfx90a")
      set(CMAKE_HIP_ARCHITECTURES "gfx90a:xnack+")
    endif()

    if(CMAKE_HIP_ARCHITECTURES MATCHES "gfx942")
      set(CMAKE_HIP_ARCHITECTURES "gfx942:xnack+")
    endif()
  endif()
elseif(CP2K_USE_ACCEL MATCHES "OPENCL")
  find_package(OpenCL REQUIRED)
  set(CP2K_USE_OPENCL ON)
endif()

message("\n------------------------------------------------------------")
message("-                          OPENMP                          -")
message("------------------------------------------------------------\n")

# PACKAGE DISCOVERY (compiler configuration can impact package discovery)
find_package(OpenMP REQUIRED COMPONENTS Fortran C CXX)

find_package(DBCSR 2.6 REQUIRED)

# ==================================
if(CP2K_USE_ELPA)
  find_package(Elpa REQUIRED)
endif()

if(CP2K_USE_LIBXC)
  find_package(Libxc 7 REQUIRED CONFIG)
endif()

if(CP2K_USE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS C Fortran)
endif()

if(CP2K_USE_COSMA)
  find_package(cosma REQUIRED)

  get_target_property(CP2K_COSMA_INCLUDE_DIRS cosma::cosma
                      INTERFACE_INCLUDE_DIRECTORIES)
  get_target_property(CP2K_COSMA_LINK_LIBRARIES cosma::cosma
                      INTERFACE_LINK_LIBRARIES)

  # check that cosma::cosma_pxgemm_cpp and cosma::cosma_prefixed_pxgemm exist
  if(NOT TARGET cosma::cosma_pxgemm_cpp OR NOT TARGET
                                           cosma::cosma_prefixed_pxgemm)
    message(
      FATAL_ERROR
        " COSMA needs to be build with scalapack offloading support. COSTA_SCALAPACK and COSMA_SCALAPACK should probably be set properly"
    )
  endif()
endif()

if(CP2K_USE_VORI)
  find_package(LibVORI REQUIRED)
endif()

if(CP2K_USE_DLAF)
  find_package(DLAFFortran 0.4.0 REQUIRED)

  get_target_property(CP2K_DLAF_INCLUDE_DIRS DLAF::DLAF
                      INTERFACE_INCLUDE_DIRECTORIES)
  get_target_property(CP2K_DLAF_LINK_LIBRARIES DLAF::dlaf.prop
                      INTERFACE_LINK_LIBRARIES)
  message("${CP2K_DLAF_INCLUDE_DIRS} ${CP2K_DLAF_LINK_LIBRARIES}")
endif()

# FFTW3

# we set this variable to ON when we want FFTW3 support (with or without MKL)
set(CP2K_USE_FFTW3_ OFF)
if(CP2K_USE_FFTW3)
  if(NOT CP2K_BLAS_VENDOR MATCHES "MKL" OR CP2K_USE_FFTW3_WITH_MKL)
    find_package(Fftw REQUIRED)
    if(CP2K_ENABLE_FFTW3_THREADS_SUPPORT AND CP2K_ENABLE_FFTW3_OPENMP_SUPPORT)
      message(
        FATAL_ERROR
          "FFTW3 threads and openmp supports can not be used at the same time")
    endif()

    if((CP2K_ENABLE_FFTW3_THREADS_SUPPORT) AND (NOT TARGET
                                                cp2k::FFTW3::fftw3_threads))
      message(
        FATAL_ERROR
          "FFTW3 was compiled without multithreading support (--enable-threads option in the FFTW3 build system)."
      )
    endif()

    if((CP2K_ENABLE_FFTW3_OPENMP_SUPPORT) AND (NOT TARGET cp2k::FFTW3::fftw3_omp
                                              ))
      message(
        FATAL_ERROR
          "FFTW3 was compiled without openmp support (--enable-openmp option in the FFTW3 build system)."
      )
    endif()
    set(CP2K_USE_FFTW3_ ON)
  else()
    message("-- Using the MKL implementation of FFTW3.")
    foreach(DIR ${CP2K_BLAS_INCLUDE_DIRS})
      if(EXISTS ${DIR}/fftw)
        list(APPEND CP2K_FFTW3_INCLUDE_DIRS ${DIR}/fftw)
        include_directories(${DIR}/fftw)
      endif()
    endforeach()
    set(CP2K_USE_FFTW3_MKL_ ON)
  endif()
endif()

# Libint
if(CP2K_USE_LIBINT2)
  find_package(Libint2 REQUIRED)
endif()

# Spglib
if(CP2K_USE_SPGLIB)
  find_package(Spglib REQUIRED CONFIG)
  get_target_property(CP2K_SPGLIB_INCLUDE_DIRS Spglib::symspg
                      INTERFACE_INCLUDE_DIRECTORIES)
endif()

if(CP2K_USE_LIBSMEAGOL)
  find_package(libsmeagol REQUIRED)
endif()

if(CP2K_USE_SPLA)
  find_package(SPLA REQUIRED CONFIG)
  get_target_property(CP2K_SPLA_INCLUDE_DIRS SPLA::spla
                      INTERFACE_INCLUDE_DIRECTORIES)
  get_target_property(CP2K_SPLA_LINK_LIBRARIES SPLA::spla
                      INTERFACE_LINK_LIBRARIES)
  if(NOT SPLA_GPU_BACKEND AND CP2K_USE_GEMM_OFFLOADING)
    set(CP2K_USE_GEMM_OFFLOADING OFF)
    message(
      FATAL_ERROR
        "SpLA should be compiled with GPU support if the gemm offloading is requested. Use -DCP2K_USE_GEMM_OFFLOADING=OFF otherwise"
    )
  endif()
endif()

if(CP2K_USE_DFTD4)
  find_package(dftd4 REQUIRED)
endif()

if(CP2K_USE_DEEPMD)
  find_package(DeePMD REQUIRED CONFIG)
endif()

if(CP2K_USE_PEXSI)
  find_package(PEXSI REQUIRED CONFIG)
  get_target_property(CP2K_PEXSI_INCLUDE_DIRS PEXSI::PEXSI
                      INTERFACE_INCLUDE_DIRECTORIES)
  get_target_property(CP2K_PEXSI_LINK_LIBRARIES PEXSI::PEXSI
                      INTERFACE_LINK_LIBRARIES)
endif()

if(CP2K_USE_ACE)
  find_package(ACE REQUIRED)
  get_target_property(CP2K_ACE_INCLUDE_DIRS ACE::pace
                      INTERFACE_INCLUDE_DIRECTORIES)
  get_target_property(CP2K_ACE_LINK_LIBRARIES ACE::pace
                      INTERFACE_LINK_LIBRARIES)
  get_target_property(LIB_PATH ACE::yaml-cpp-pace INTERFACE_LINK_LIBRARIES)
  list(APPEND CP2K_ACE_LINK_LIBRARIES ${LIB_PATH})
  get_target_property(LIB_PATH ACE::cnpy INTERFACE_LINK_LIBRARIES)
  list(APPEND CP2K_ACE_LINK_LIBRARIES ${LIB_PATH})
endif()

if(CP2K_USE_TBLITE)
  find_package(mctc-lib REQUIRED)
  find_package(toml-f REQUIRED)
  find_package(s-dftd3 REQUIRED)
  find_package(tblite REQUIRED)
  target_link_libraries(
    tblite::tblite INTERFACE mctc-lib::mctc-lib dftd4::dftd4 toml-f::toml-f
                             s-dftd3::s-dftd3)
endif()

# SIRIUS
if(CP2K_USE_SIRIUS)
  find_package(sirius 7.7.0 REQUIRED)
endif()

if(CP2K_USE_PLUMED)
  find_package(Plumed REQUIRED)
endif()

if(CP2K_USE_LIBTORCH)
  find_package(Torch REQUIRED)
  get_target_property(CP2K_TORCH_CXX_STANDARD torch CXX_STANDARD)
  if(${CP2K_TORCH_CXX_STANDARD} VERSION_GREATER ${CMAKE_CXX_STANDARD})
    set(CMAKE_CXX_STANDARD ${CP2K_TORCH_CXX_STANDARD})
    message("-- Torch setting CMAKE_CXX_STANDARD to ${CP2K_TORCH_CXX_STANDARD}")
  endif()
endif()

if(CP2K_USE_MPI_F08 AND NOT MPI_Fortran_HAVE_F08_MODULE)
  message(
    FATAL_ERROR
      "The Fortran 2008 interface is not supported by the MPI implementation found by cmake."
  )
endif()

if(CP2K_USE_TREXIO)
  find_package(TrexIO REQUIRED)
endif()

if(CP2K_USE_GREENX)
  find_package(greenX REQUIRED CONFIG)
  get_target_property(CP2K_GREENX_INCLUDE_DIRS greenX::GXCommon
                      INTERFACE_INCLUDE_DIRECTORIES)
  get_target_property(CP2K_GREENX_LINK_LIBRARIES greenX::GXCommon LOCATION)
  get_target_property(LIB_GXMiniMax greenX::LibGXMiniMax LOCATION)
  list(APPEND CP2K_GREENX_LINK_LIBRARIES ${LIB_GXMiniMax})
  get_target_property(LIB_GXAC greenX::LibGXAC LOCATION)
  list(APPEND CP2K_GREENX_LINK_LIBRARIES ${LIB_GXAC})
  # GreenX with AC requires GMP library
  find_library(GMP_LIBRARY NAMES gmp REQUIRED)
  get_filename_component(GMP_LIBRARY_PATH ${GMP_LIBRARY} DIRECTORY)
  set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-L${GMP_LIBRARY_PATH}")
  set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS}
                                "-L${GMP_LIBRARY_PATH}")
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

include(GNUInstallDirs)

# subdirectories
add_subdirectory(src)

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
  "  - BLAS\n" # keep line break
  "    - vendor: ${CP2K_BLAS_VENDOR}\n"
  "    - include directories: ${CP2K_BLAS_INCLUDE_DIRS}\n"
  "    - libraries: ${CP2K_BLAS_LINK_LIBRARIES}\n\n")

message("  - LAPACK\n" # keep line break
        "    - include directories: ${CP2K_LAPACK_INCLUDE_DIRS}\n"
        "    - libraries: ${CP2K_LAPACK_LINK_LIBRARIES}\n\n")

if(CP2K_USE_MPI)
  message("  - MPI\n" # keep line break
          "    - include directories: ${MPI_INCLUDE_DIRS}\n"
          "    - libraries: ${MPI_LIBRARIES}\n\n")
  if(CP2K_USE_MPI_F08)
    message("  - MPI_F08: ON\n")
  endif()
  if(MPI_Fortran_HAVE_F08_MODULE AND NOT CP2K_USE_MPI_F08)
    message(
      "  - MPI_F08 is supposed by MPI but turned off by default.\n"
      "    To use it add -DCP2K_USE_MPI_F08=ON to the cmake command line\n\n")
  endif()
  message(
    "  - ScaLAPACK\n" # keep line break
    "    - vendor: ${CP2K_SCALAPACK_VENDOR}\n"
    "    - include directories: ${CP2K_SCALAPACK_INCLUDE_DIRS}\n"
    "    - libraries: ${CP2K_SCALAPACK_LINK_LIBRARIES}\n\n")
endif()

if((CP2K_USE_ACCEL MATCHES "CUDA") OR (CP2K_USE_ACCEL MATCHES "HIP"))

  message("  - Hardware Acceleration:\n")
  if(CP2K_USE_ACCEL MATCHES "CUDA")
    message("   - CUDA:\n" # let below line separate
            "     - GPU architecture number: ${CMAKE_CUDA_ARCHITECTURES}\n"
            "     - GPU profiling enabled: ${CP2K_WITH_CUDA_PROFILING}\n\n")
  endif()

  if(CP2K_USE_ACCEL MATCHES "HIP")
    message(
      "   - HIP:\n" # let below line separate
      "     - GPU target architecture: ${CP2K_WITH_GPU}\n"
      "     - GPU architecture number: ${CP2K_ACC_ARCH_NUMBER}\n"
      "     - FLAGS: ${CMAKE_HIP_FLAGS}")
  endif()

  message(
    "      - GPU accelerated modules\n"
    "        - PW     module: ${CP2K_ENABLE_PW_GPU}\n"
    "        - GRID   module: ${CP2K_ENABLE_GRID_GPU}\n"
    "        - DBM    module: ${CP2K_ENABLE_DBM_GPU}\n\n")
endif()

if(CP2K_USE_CUSOLVER_MP)
  message(
    "   - CUSolverMP: \n"
    "     - library: ${CP2K_CUSOLVER_MP_LINK_LIBRARIES} \n"
    "     - include: ${CP2K_CUSOLVER_MP_INCLUDE_DIRS} \n"
    "     - CAL library: ${CP2K_CAL_LINK_LIBRARIES} \n"
    "     - CAL include: ${CP2K_CAL_INCLUDE_DIRS} \n"
    "     - ucc library: ${CP2K_UCC_LINK_LIBRARIES} \n"
    "     - ucx library: ${CP2K_UCX_LINK_LIBRARIES} \n"
    "     - ucc include: ${CP2K_UCC_INCLUDE_DIRS} \n")
endif()

if(CP2K_USE_LIBXC)
  message(
    "  - LibXC\n" # let below line separate
    "    - version:  ${Libxc_VERSION}\n"
    "    - include directories: ${Libxc_INCLUDE_DIRS}\n"
    "    - libraries: ${Libxc_LIBRARIES}\n\n")
endif()

if(CP2K_USE_SPGLIB)
  message("  - Spglib\n"
          "    - include directories: ${CP2K_SPGLIB_INCLUDE_DIRS}\n\n")
endif()

if(CP2K_USE_LIBTORCH)
  message("  - LibTorch\n" "    - extra CXX flags: ${TORCH_CXX_FLAGS}\n"
          "    - include directories: ${TORCH_INCLUDE_DIRS}\n"
          "    - libraries: ${TORCH_LIBRARY}\n")
endif()

if(CP2K_USE_HDF5)
  message(
    "  - HDF5\n" # let below line separate
    "    - version: ${HDF5_VERSION}\n"
    "    - include directories: ${HDF5_INCLUDE_DIRS}\n"
    "    - libraries: ${HDF5_LIBRARIES}\n\n")
endif()

if(CP2K_USE_FFTW3)
  message("  - FFTW3\n"
          "    - include directories: ${CP2K_FFTW3_INCLUDE_DIRS}\n"
          "    - libraries: ${CP2K_FFTW3_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_PLUMED)
  message("  - PLUMED\n"
          "    - include directories: ${CP2K_PLUMED_INCLUDE_DIRS}\n"
          "    - libraries: ${CP2K_PLUMED_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_LIBXSMM)
  message(
    "  - LIBXSMM\n"
    "    - include directories: ${CP2K_LIBXSMM_INCLUDE_DIRS}\n"
    "    - libraries: ${CP2K_LIBXSMMEXT_LINK_LIBRARIES};${CP2K_LIBXSMMF_LINK_LIBRARIES}\n\n"
  )
endif()

if(CP2K_USE_SPLA)
  message(" - SpLA\n" # keep line break
          "   - include directories: ${CP2K_SPLA_INCLUDE_DIRS}\n"
          "   - libraries: ${CP2K_SPLA_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_DFTD4)
  message(" - DFTD4\n" "   - include directories :  ${dftd4_INCLUDE_DIRS}\n"
          "   - libraries           :  ${dftd4_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_DEEPMD)
  message(" - DeePMD\n\n")
endif()

if(CP2K_USE_PEXSI)
  message(" - PEXSI\n" # keep line break
          "   - include directories: ${CP2K_PEXSI_INCLUDE_DIRS}\n"
          "   - libraries: ${CP2K_PEXSI_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_ACE)
  message(" - ACE\n" # keep line break
          "   - include directories: ${CP2K_ACE_INCLUDE_DIRS}\n"
          "   - libraries: ${CP2K_ACE_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_LIBSMEAGOL)
  message(" - LibSMEAGOL\n"
          "   - include directories: ${CP2K_LIBSMEAGOL_INCLUDE_DIRS}\n"
          "   - libraries: ${CP2K_LIBSMEAGOL_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_TBLITE)
  message(" - TBLITE :\n"
          "   - include directories :  ${tblite_INCLUDE_DIRS}\n"
          "   - tblite libraries    :  ${tblite_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_SIRIUS)
  message(" - SIRIUS\n"
          "   - include directories: ${CP2K_SIRIUS_INCLUDE_DIRS}\n"
          "   - libraries: ${CP2K_SIRIUS_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_COSMA)
  message(" - COSMA\n" # keep line break
          "   - include directories: ${CP2K_COSMA_INCLUDE_DIRS}\n"
          "   - libraries: ${CP2K_COSMA_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_LIBINT2)
  message(" - Libint2\n"
          "   - include directories: ${CP2K_LIBINT2_INCLUDE_DIRS}\n"
          "   - libraries: ${CP2K_LIBINT2_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_VORI)
  message(" - Libvori\n" "   - libraries: ${CP2K_LIBVORI_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_ELPA)
  message(" - ELPA\n" # keep line break
          "   - include directories: ${CP2K_ELPA_INCLUDE_DIRS}\n"
          "   - libraries: ${CP2K_ELPA_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_DLAF)
  message(" - DLA-Future\n"
          "   - include directories: ${CP2K_DLAF_INCLUDE_DIRS}\n"
          "   - libraries: ${CP2K_DLAF_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_GRPP)
  message(" - GRPP\n")
endif()

if(CP2K_USE_TREXIO)
  message(" - TREXIO\n"
          "   - include directories: ${CP2K_TREXIO_INCLUDE_DIRS}\n"
          "   - libraries: ${CP2K_TREXIO_LINK_LIBRARIES}\n\n")
endif()

if(CP2K_USE_GREENX)
  message(" - GreenX\n"
          "   - include directories: ${CP2K_GREENX_INCLUDE_DIRS}\n"
          "   - libraries: ${CP2K_GREENX_LINK_LIBRARIES}\n\n")
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

if(NOT CP2K_USE_DFTD4)
  message("   - DFTD4")
endif()

if(NOT CP2K_USE_DEEPMD)
  message("   - DeePMD")
endif()

if(NOT CP2K_USE_PEXSI)
  message("   - PEXSI")
endif()

if(NOT CP2K_USE_ACE)
  message("   - ACE (libpace)")
endif()

if(NOT CP2K_USE_TBLITE)
  message("   - TBLITE")
endif()

if(NOT CP2K_USE_SIRIUS)
  message("   - SIRIUS")
endif()

if(NOT CP2K_USE_SPGLIB)
  message("   - Spglib")
endif()

if(NOT CP2K_USE_LIBSMEAGOL)
  message("   - LibSMEAGOL")
endif()

if(NOT CP2K_USE_COSMA)
  message("   - COSMA")
endif()

if(NOT CP2K_USE_SPLA)
  message("   - SpLA")
endif()

if(NOT CP2K_USE_HDF5)
  message("   - HDF5")
endif()

if(${CP2K_USE_ACCEL} MATCHES "NONE")
  message("   - GPU acceleration is disabled")
endif()

if(NOT CP2K_USE_ELPA)
  message("   - ELPA")
endif()

if(NOT CP2K_USE_DLAF)
  message("   - DLA-Future")
endif()

if(NOT CP2K_USE_PLUMED)
  message("   - PLUMED")
endif()

if(NOT CP2K_USE_LIBXSMM)
  message("   - LIBXSMM")
endif()

if(NOT CP2K_USE_LIBINT2)
  message("   - Libint2")
endif()

if(NOT CP2K_USE_LIBXC)
  message("   - LibXC")
endif()

if(NOT CP2K_USE_VORI)
  message("   - Libvori")
endif()

if(NOT CP2K_USE_FFTW3)
  message("   - FFTW3")
endif()

if(NOT CP2K_USE_LIBTORCH)
  message("   - LibTorch")
endif()

if(NOT CP2K_USE_TREXIO)
  message("   - TREXIO")
endif()

if(NOT CP2K_USE_GREENX)
  message("   - GreenX")
endif()

if(NOT CP2K_USE_GRPP)
  message("   - GRPP")
endif()

message(
  "\n\n" # let below line separate
  "To run the regtests you need to run the following commands\n"
  "\n\n cd ..\n" # let below line separate
  " export CP2K_DATA_DIR=${CMAKE_SOURCE_DIR}/data/\n"
  " ./tests/do_regtest.py ${cp2k_BINARY_DIR}/bin ${__cp2k_ext}\n\n")

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
