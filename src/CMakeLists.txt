#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(fypp-sources)
include(GNUInstallDirs) # required to get a proper LIBDIR variable
include(CMakePackageConfigHelpers)
configure_file(base/base_uses.f90 base/base_uses.f90 @ONLY)

# take all fortran files we can find.
file(GLOB_RECURSE CP2K_SRCS_F "*.F")

# search for c files that are not included in the dependencies

set(CP2K_SRCS_C "sockets.c")
set(TMPFILES "")
foreach(__dir "pw" "pw/fpga" "pw/gpu" "start" "base")
  file(GLOB_RECURSE TMPFILES "${__dir}/*.c")
  list(APPEND CP2K_SRCS_C "${TMPFILES}")
  set(TMPFILES "")
endforeach()

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

# remove those files from the cp2k library as they contain a main function

foreach(__app __${__CP2K_APPS})
  list(REMOVE_ITEM CP2K_SRCS_F "${__app}.F")
endforeach()
set(CP2K_SRCS_GPU "pw/gpu/pw_gpu_kernels.cu")

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

set_target_properties(cp2k PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "CXX")
set_target_properties(cp2k PROPERTIES Fortran_MODULE_DIRECTORY mod_files)

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

add_subdirectory(offload)
add_subdirectory(grid)
add_subdirectory(dbm)

# mix the target and variables pointing to the include directories.

include_directories(
  ${CMAKE_MPI_INCLUDE_DIRECTORIES}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/base>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/common>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/dbcsrx>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/motion>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/dbm>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files>)

if(DEFINED ENV{CPATH})
  string(REPLACE ":" ";" INCLUDE_LIST $ENV{CPATH})
  include_directories(${INCLUDE_LIST})
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
    $<$<BOOL:${CP2K_USE_TORCH}>:${TORCH_LIBRARIES}>
    $<$<BOOL:${CP2K_USE_COSMA}>:cosma::cosma_prefixed_pxgemm>
    $<$<BOOL:${CP2K_USE_COSMA}>:cosma::cosma>
    DBCSR::dbcsr
    cp2k_grid
    cp2k_dbm
    $<$<BOOL:${CP2K_USE_CUDA}>:CUDA::cufft>
    $<$<BOOL:${CP2K_USE_CUDA}>:CUDA::cufftw>
    $<$<BOOL:${CP2K_USE_CUDA}>:CUDA::cublas>
    $<$<BOOL:${CP2K_USE_CUDA}>:CUDA::cudart>
    $<$<BOOL:${CP2K_USE_HIP}>:hip::hipfft>
    $<$<BOOL:${CP2K_USE_HIP}>:roc::hipblas>
    $<$<BOOL:${CP2K_USE_HIP}>:hip::host>
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

  target_compile_definitions(
    cp2k
    PUBLIC $<$<NOT:$<BOOL:${CP2K_DBCSR_USE_CPU_ONLY}>>:__DBCSR_ACC>
           $<$<BOOL:${CP2K_DISABLE_GRID_GPU}>:__NO_OFFLOAD_GRID>
           $<$<BOOL:${CP2K_DISABLE_PW_GPU}>:__NO_OFFLOAD_PW>
           $<$<BOOL:${CP2K_DISABLE_DBM_GPU}>:__NO_OFFLOAD_DBM>)
endif()

if(CP2K_USE_HIP)
  set_source_files_properties(${CP2K_SRCS_GPU} PROPERTIES LANGUAGE HIP)
  set_target_properties(cp2k PROPERTIES HIP_ARCHITECTURES
                                        "${CP2K_ACC_ARCH_NUMBER}")
endif()

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

# installation

foreach(_app ${__CP2K_APPS} cp2k-bin)
  install(
    TARGETS ${__CP2K_APPS}
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    INCLUDES
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cp2k")
endforeach()

install(DIRECTORY "${CMAKE_SOURCE_DIR}/data"
        DESTINATION "${CMAKE_INSTALL_FULL_DATAROOTDIR}/cp2k")

install(
  TARGETS cp2k cp2k_link_libs cp2k_grid cp2k_offload cp2k_dbm
  EXPORT cp2k_targets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(
  EXPORT cp2k_targets
  FILE cp2kTargets.cmake
  NAMESPACE cp2k::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/cp2k")
install(
  DIRECTORY "${PROJECT_BINARY_DIR}/src/mod_files"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cp2k"
  FILES_MATCHING
  PATTERN "*.mod")
