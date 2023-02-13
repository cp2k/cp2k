#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(fypp-sources)
include(GNUInstallDirs) # required to get a proper LIBDIR variable
include(CMakePackageConfigHelpers)

add_fypp_sources(CP2K_SRCS ${CP2K_SRCS_F})

# set the __SHORT_FILE__ per file for CP2K sources
foreach(cp2k_src ${CP2K_SRCS})
  # add_fypp_sources returns a path in the current binary dir
  get_filename_component(short_file "${cp2k_src}" NAME)
  set_source_files_properties(
    ${cp2k_src} PROPERTIES COMPILE_DEFINITIONS __SHORT_FILE__="${short_file}")
endforeach()

add_library(cp2k "${CP2K_SRCS};${CP2K_SRCS_C};${CP2K_SRCS_GPU}")

target_compile_features(cp2k PUBLIC cxx_std_14)
target_compile_features(cp2k PUBLIC c_std_99)
target_compile_features(cp2k PUBLIC cuda_std_11)

set_target_properties(cp2k PROPERTIES CXX_STANDARD 14)
set_target_properties(cp2k PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "CXX")

# =================================================================================================
# main CP2K OBJECT LIBRARY
if(CP2K_USE_ACCEL MATCHES "CUDA")
  set_property(TARGET cp2k PROPERTY CUDA_ARCHITECTURES ${CP2K_ACC_ARCH_NUMBER})
  target_compile_definitions(
    cp2k PUBLIC __OFFLOAD_CUDA $<$<COMPILE_LANGUAGE:CUDA>:__OFFLOAD_CUDA>)
  include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
elseif(CP2K_USE_ACCEL MATCHES "HIP")
  set_property(TARGET cp2k PROPERTY HIP_ARCHITECTURES ${CP2K_ACC_ARCH_NUMBER})
endif()

add_library(cp2k_link_libs INTERFACE)
target_link_libraries(
  cp2k_link_libs
  INTERFACE
    $<$<BOOL:${CP2K_USE_SIRIUS}>:sirius::sirius>
    $<$<BOOL:${CP2K_USE_VORI}>:CP2K::VORI::vori>
    $<$<BOOL:${CP2K_USE_PEXSI}>:CP2K::PEXSI::pexsi>
    $<$<BOOL:${CP2K_USE_PEXSI}>:CP2K::ptscotch::ptscotch>
    $<$<BOOL:${CP2K_USE_SPGLIB}>:CP2K::LIBSPG::libspg>
    $<$<BOOL:${CP2K_USE_LIBXC}>:CP2K::Libxc::xc>
    $<$<BOOL:${CP2K_USE_ELPA}>:CP2K::ELPA::elpa>
    $<$<BOOL:${CP2K_USE_FFTW3}>:CP2K::FFTW3::fftw3>
    $<$<BOOL:${CP2K_ENABLE_FFTW3_THREADS_SUPPORT}>:CP2K::FFTW3::fftw3_threads>
    $<$<BOOL:${CP2K_ENABLE_FFTW3_OPENMP_SUPPORT}>:CP2K::FFTW3::fftw3_omp>
    $<$<BOOL:${CP2K_USE_SPLA}>:SPLA::spla>
    $<$<BOOL:${CP2K_USE_LIBINT2}>:CP2K::Libint2::int2>
    $<$<BOOL:${CP2K_USE_COSMA}>:cosma::cosma_prefixed_pxgemm>
    $<$<BOOL:${CP2K_USE_COSMA}>:cosma::cosma>
    $<$<BOOL:${CP2K_USE_COSMA}>:costa::costa>
    DBCSR::dbcsr
    $<$<BOOL:${CP2K_USE_TORCH}>:${TORCH_LIBRARIES}>
    $<$<BOOL:${CP2K_USE_CUDA}>:CUDA::cufft>
    $<$<BOOL:${CP2K_USE_CUDA}>:CUDA::cufftw>
    $<$<BOOL:${CP2K_USE_CUDA}>:CUDA::cublas>
    $<$<BOOL:${CP2K_USE_CUDA}>:CUDA::cudart>
    $<$<BOOL:${CP2K_USE_HIP}>:hip::hipfft>
    $<$<BOOL:${CP2K_USE_HIP}>:roc::hipblas>
    $<$<BOOL:${CP2K_USE_HIP}>:hip::device>
    $<$<BOOL:${CP2K_USE_MPI}>:CP2K::SCALAPACK::scalapack>
    $<$<BOOL:${CP2K_USE_LIBXSMM}>:CP2K::LibXSMM::libxsmmf>
    $<$<BOOL:${CP2K_USE_LIBXSMM}>:CP2K::LibXSMM::libxsmmext>
    CP2K::LAPACK::lapack
    CP2K::BLAS::blas
    $<$<BOOL:${CP2K_USE_MPI}>:MPI::MPI_Fortran>
    $<$<BOOL:${CP2K_USE_MPI}>:MPI::MPI_C>
    $<$<BOOL:${CP2K_USE_MPI}>:MPI::MPI_CXX>
    OpenMP::OpenMP_Fortran
    OpenMP::OpenMP_C
    OpenMP::OpenMP_CXX)

# mix the target and variables pointing to the include directories.

include_directories(
  ${CMAKE_MPI_INCLUDE_DIRECTORIES}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/base>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/common>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/dbcsrx>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/motion>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/grid>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/dbm>)

if(DEFINED ENV{CPATH})
  string(REPLACE ":" ";" INCLUDE_LIST $ENV{CPATH})
  include_directories(${INCLUDE_LIST})
endif()

string(TIMESTAMP CP2K_TIMESTAMP "%Y-%m-%d %H:%M:%S")

target_link_libraries(cp2k PUBLIC cp2k_link_libs)
target_compile_definitions(
  cp2k
  PUBLIC $<$<BOOL:${CP2K_USE_MPI}>:__parallel>
         $<$<BOOL:${CP2K_USE_MPI}>:__SCALAPACK>
         $<$<BOOL:${CP2K_ENABLE_F08_MPI}>:__MPI_08>
         __COMPILE_DATE=\"${CP2K_TIMESTAMP}\"
         __COMPILE_HOST=\"${CP2K_HOST_NAME}\"
         __COMPILE_REVISION=\"${CP2K_GIT_HASH}\"
         __DATA_DIR=\"${CMAKE_INSTALL_FULL_DATAROOTDIR}/cp2k/data\"
         __COMPILE_ARCH=\"${CMAKE_SYSTEM_PROCESSOR}\"
         $<$<CONFIG:Release>:NDEBUG>
         $<$<CONFIG:COVERAGE>:__NO_ABORT>
         $<$<CONFIG:DEBUG>:__HAS_IEEE_EXCEPTIONS>
         $<$<CONFIG:DEBUG>:__CHECK_DIAG>
         $<$<BOOL:${CP2K_USE_PEXSI}>:__PEXSI>
         $<$<BOOL:${CP2K_USE_PLUMED2}>:__PLUMED2>
         $<$<BOOL:${CP2K_USE_QUIP}>:__QUIP>
         $<$<BOOL:${CP2K_USE_MAXWELL}>:__LIBMAXWELL>
         $<$<BOOL:${CP2K_USE_VORI}>:__LIBVORI>
         $<$<BOOL:${CP2K_USE_SPGLIB}>:__SPGLIB>
         $<$<BOOL:${CP2K_USE_SIRIUS}>:__SIRIUS>
         $<$<BOOL:${CP2K_USE_SPLA}>:__SPLA>
         $<$<BOOL:${CP2K_USE_SpFFT}>:__SPFFT>
         $<$<BOOL:${CP2K_USE_SPLA_GEMM_OFFLOADING}>:__OFFLOAD_GEMM>
         $<$<BOOL:${CP2K_USE_ELPA}>:__ELPA>
         $<$<BOOL:${CP2K_USE_LIBXC}>:__LIBXC>
         $<$<BOOL:${CP2K_USE_FFTW3}>:__FFTW3>
         $<$<BOOL:${CP2K_USE_LIBINT2}>:__LIBINT>
         $<$<BOOL:${CP2K_USE_PEXSI}>:__LIBPEXSI>
         $<$<BOOL:${CP2K_USE_LIBTORCH}>:__LIBTORCH>
         $<$<BOOL:${CP2K_USE_COSMA}>:__COSMA>
         $<$<BOOL:${CP2K_USE_LIBXSMM}>:__LIBXSMM>
         $<$<STREQUAL:${CP2K_BLAS_VENDOR},MKL>:__MKL>
         $<$<STREQUAL:${CP2K_BLAS_VENDOR},Apple>:__ACCELERATE>
         $<$<BOOL:${CP2K_USE_CUDA}>:__OFFLOAD_CUDA>
         $<$<COMPILE_LANGUAGE:CUDA>:__OFFLOAD_CUDA>
         $<$<BOOL:${CP2K_USE_HIP}>:__OFFLOAD_HIP>
         $<$<BOOL:${CP2K_USE_HIP}>:__HIP_PLATFORM_AMD__>
         $<$<COMPILE_LANGUAGE:HIP>:__OFFLOAD_HIP>
         $<$<COMPILE_LANGUAGE:HIP>:__HIP_PLATFORM_AMD__>)

if(CP2K_USE_CUDA OR CP2K_USE_HIP)
  # these checks are only relevant when the debug mode is on

  if(NOT CP2K_DBCSR_USE_CPU_ONLY)
    target_compile_definitions(cp2k PUBLIC __DBCSR_ACC)
  endif()

  # the next three are off by default
  if(CP2K_DISABLE_GRID_GPU)
    target_compile_definitions(cp2k PUBLIC __NO_OFFLOAD_GRID)
  endif()

  if(CP2K_DISABLE_PW_GPU)
    target_compile_definitions(cp2k PUBLIC __NO_OFFLOAD_PW)
  endif()

  if(CP2K_DISABLE_DBM_GPU)
    target_compile_definitions(cp2k PUBLIC __NO_OFFLOAD_DBM)
  endif()
endif()

if(CP2K_USE_HIP)
  foreach(__src ${CP2K_SRCS_GPU})
    set_source_files_properties(${__src} PROPERTIES LANGUAGE HIP)
  endforeach()
  set_target_properties(cp2k PROPERTIES HIP_ARCHITECTURES
                                        "${CP2K_ACC_ARCH_NUMBER}")
endif()

list(
  APPEND
  __CP2K_APPS
  "memory_utilities_unittest"
  "parallel_rng_types_unittest"
  "graph"
  "dumpdcd"
  "xyz2dcd"
  "libcp2k_unittest"
  "nequip_unittest")

add_executable(cp2k-bin start/cp2k.F)
set_target_properties(
  cp2k-bin PROPERTIES LINKER_LANGUAGE Fortran) # always use the Fortran compiler
                                               # for
# linking

add_executable(memory_utilities_unittest common/memory_utilities_unittest.F)
add_executable(parallel_rng_types_unittest common/parallel_rng_types_unittest.F)
add_executable(graph metadyn_tools/graph.F)
add_executable(dumpdcd motion/dumpdcd.F)
add_executable(xyz2dcd motion/xyz2dcd.F)
add_executable(libcp2k_unittest start/libcp2k_unittest.c)
add_executable(nequip_unittest nequip_unittest.F)

set_target_properties(cp2k-bin PROPERTIES CXX_STANDARD 14)
set_target_properties(cp2k-bin PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES
                                          "CXX")
set_target_properties(cp2k-bin PROPERTIES LINKER_LANGUAGE "CXX")
set_target_properties(cp2k-bin PROPERTIES OUTPUT_NAME cp2k)
target_link_libraries(cp2k-bin PUBLIC cp2k)

foreach(_app ${__CP2K_APPS})
  set_target_properties(${_app} PROPERTIES LINKER_LANGUAGE "CXX")
  set_target_properties(${_app} PROPERTIES CXX_STANDARD 14)
  set_target_properties(${_app} PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES
                                           "CXX")
  target_compile_features(${_app} PUBLIC c_std_99)
  if(CP2K_USE_ACCEL MATCHES HIP)
    set_target_properties(${_app} PROPERTIES HIP_ARCHITECTURES
                                             ${CP2K_ACC_ARCH_NUMBER})
  endif()
  target_link_libraries(${_app} PUBLIC cp2k)
endforeach()

#
# apps included in the grid module
#

list(APPEND __GRID_APPS "grid_miniapp" "grid_unittest")

add_executable(grid_miniapp grid/grid_miniapp.c)
target_link_libraries(grid_miniapp PUBLIC cp2k)

add_executable(grid_unittest grid/grid_unittest.c)
target_link_libraries(grid_unittest PUBLIC cp2k)

foreach(_app ${__GRID_APPS})
  set_target_properties(${_app} PROPERTIES LINKER_LANGUAGE C)
  target_compile_features(${_app} PUBLIC c_std_99)

  target_link_libraries(
    ${_app}
    PUBLIC $<$<BOOL:${CP2K_USE_CUDA}>:CUDA::cudart>
           $<$<BOOL:${CP2K_USE_HIP}>:hip::device>
           $<$<BOOL:${CP2K_USE_LIBXSMM}>:CP2K::LibXSMM::libxsmmf>
           $<$<BOOL:${CP2K_USE_LIBXSMM}>:CP2K::LibXSMM::libxsmm>
           CP2K::BLAS::blas
           OpenMP::OpenMP_C
           m)

  target_include_directories(
    ${_app} PUBLIC $<$<BOOL:${CP2K_USE_CUDA}>:${CMAKE_CUDA_INCLUDE_DIRECTORIES}>
                   $<$<BOOL:${CP2K_USE_HIP}>:${CMAKE_HIP_INCLUDE_DIRECTORIES}>)

  target_compile_definitions(
    ${_app}
    PUBLIC $<$<BOOL:${CP2K_USE_CUDA}>:__OFFLOAD_CUDA>
           $<$<BOOL:${CP2K_USE_HIP}>:__OFFLOAD_HIP>
           $<$<STREQUAL:CP2K_BLAS_VENDOR,"MKL">:__MKL>
           $<$<BOOL:${CP2K_USE_LIBXSMM}>:__LIBXSMM>)

  target_compile_features(${_app} PUBLIC $<$<BOOL:CP2K_USE_CUDA>:cuda_std_11>)

  if(CP2K_USE_CUDA)
    set_target_properties(${_app} PROPERTIES CUDA_ARCHITECTURES
                                             ${CP2K_ACC_ARCH_NUMBER})
  endif()

  # cmake 3.21 is needed for this
  if(CP2K_USE_HIP)
    set_target_properties(${_app} PROPERTIES HIP_ARCHITECTURES
                                             ${CP2K_ACC_ARCH_NUMBER})
  endif()
endforeach()

#
# apps included in dbm
#

add_executable(dbm_miniapp dbm/dbm_miniapp.c)
target_link_libraries(dbm_miniapp PUBLIC cp2k)
set_target_properties(dbm_miniapp PROPERTIES LINKER_LANGUAGE C)
target_compile_features(dbm_miniapp PUBLIC c_std_99)

# TODO target_link_libraries( dbm_miniapp PUBLIC
# $<$<BOOL:${CP2K_USE_CUDA}>:CUDA::cudart>
# $<$<BOOL:${CP2K_USE_HIP}>:hip::device>
# $<$<BOOL:${CP2K_USE_LIBXSMM}>:CP2K::LibXSMM::libxsmm>
# $<$<BOOL:${CP2K_USE_LIBXSMM}>:CP2K::LibXSMM::libxsmmext>
# $<$<BOOL:${CP2K_USE_MPI}>:MPI::MPI_C> CP2K::BLAS::blas OpenMP::OpenMP_C m)

target_compile_definitions(
  dbm_miniapp
  PUBLIC $<$<BOOL:${CP2K_USE_CUDA}>:__OFFLOAD_CUDA>
         $<$<BOOL:${CP2K_USE_HIP}>:__OFFLOAD_HIP>
         $<$<BOOL:${CP2K_USE_LIBXSMM}>:__LIBXSMM>)

if(CP2K_DISABLE_DBM_GPU)
  if(CP2K_USE_CUDA)
    target_compile_features(dbm_miniapp PUBLIC cuda_std_11)
    set_target_properties(dbm_miniapp PROPERTIES CUDA_ARCHITECTURES
                                                 ${CP2K_ACC_ARCH_NUMBER})
  endif()

  # cmake 3.21 is needed for this
  if(CP2K_USE_HIP)
    set_target_properties(dbm_miniapp PROPERTIES HIP_ARCHITECTURES
                                                 ${CP2K_ACC_ARCH_NUMBER})
  endif()
endif()

# installation

foreach(_app ${__GRID_APPS} ${__CP2K_APPS} cp2k-bin)
  install(
    TARGETS ${_app}
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    INCLUDES
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cp2k")
endforeach()

install(DIRECTORY "${CMAKE_SOURCE_DIR}/data"
        DESTINATION "${CMAKE_INSTALL_FULL_DATAROOTDIR}/cp2k")
install(TARGETS cp2k LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
