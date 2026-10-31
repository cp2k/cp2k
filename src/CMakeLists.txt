#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: GPL-2.0-or-later                                                     !
#!-------------------------------------------------------------------------------------------------!

include(fypp-sources)
include(cp2k_utils)
include(GNUInstallDirs) # required to get a proper LIBDIR variable
include(CMakePackageConfigHelpers)
configure_file(base/base_uses.f90 base/base_uses.f90 @ONLY)

# Source lists are kept separately from source selection and build logic.
include("${CMAKE_CURRENT_LIST_DIR}/CMakeSources.cmake")

if(CP2K_USE_PW_GPU)
  list(APPEND CP2K_SRCS_F pw/pw_gpu.F)
endif()

list(APPEND CP2K_SRCS_C "${CP2K_OFFLOAD_SRCS_C}" "${CP2K_MPIWRAP_SRCS_C}")

# ##############################################################################
# gpu code
# ##############################################################################

set(CP2K_SRCS_GPU "")

if(CP2K_USE_MPI)
  list(APPEND CP2K_PROGS_F mpiwrap/message_passing_thread_unittest.F
       fm/cp_fm_diag_utils_unittest.F)
else()
  list(APPEND CP2K_EXCLUDED_SRCS_F mpiwrap/message_passing_thread_unittest.F
       fm/cp_fm_diag_utils_unittest.F)
endif()

if(CP2K_USE_CUDA
   OR CP2K_USE_HIP
   OR CP2K_USE_OPENCL)
  if(CP2K_ENABLE_PW_GPU)
    list(APPEND CP2K_SRCS_GPU "pw/gpu/pw_gpu_kernels.cu")
    list(APPEND CP2K_SRCS_C ${CP2K_PW_SRCS_C})
  endif()

  if(CP2K_ENABLE_DBM_GPU)
    list(APPEND CP2K_DBM_SRCS_C dbm/dbm_multiply_gpu.c)
    if(CP2K_USE_OPENCL)
      set(CP2K_DBM_OPENCL_KERNEL dbm/dbm_multiply_opencl.cl)
      set(CP2K_DBM_OPENCL_GENHDR
          "${CMAKE_CURRENT_BINARY_DIR}/${CP2K_DBM_OPENCL_KERNEL}.h")
      # cmake-lint: disable=E1126
      add_custom_command(
        COMMAND
          ${CP2K_OPENCL_SCRIPT} -b 6 -p \"\" -I "${CP2K_OPENCL_BACKEND}"
          ${CMAKE_CURRENT_SOURCE_DIR}/${CP2K_DBM_OPENCL_KERNEL}
          "${CP2K_DBM_OPENCL_GENHDR}"
        DEPENDS ${CP2K_OPENCL_SCRIPT} ${CP2K_OPENCL_COMMON}
                ${CMAKE_CURRENT_SOURCE_DIR}/${CP2K_DBM_OPENCL_KERNEL}
        OUTPUT "${CP2K_DBM_OPENCL_GENHDR}"
        COMMENT "Generating OpenCL kernel for DBM")
      add_custom_target(cp2k_dbm_opencl_kernel
                        DEPENDS "${CP2K_DBM_OPENCL_GENHDR}")
      include_directories(${CMAKE_CURRENT_BINARY_DIR}/dbm)
      include_directories(${CP2K_OPENCL_BACKEND})
      if(CP2K_OPENCL_SAMPLES)
        include_directories(${CP2K_OPENCL_SAMPLES})
      endif()
      list(APPEND CP2K_DBM_SRCS_C dbm/dbm_multiply_opencl.c)
    else()
      list(APPEND CP2K_SRCS_GPU ${CP2K_DBM_SRCS_GPU})
    endif()
  endif()

  if(CP2K_ENABLE_GRID_GPU)
    list(APPEND CP2K_SRCS_GPU ${CP2K_GRID_SRCS_GPU})
  endif()
  if(CP2K_USE_HIP)
    set_source_files_properties(${CP2K_SRCS_GPU} PROPERTIES LANGUAGE HIP)
  endif()
endif()

list(APPEND CP2K_SRCS_C ${CP2K_DBM_SRCS_C} ${CP2K_GRID_SRCS_C})

# ##############################################################################
# Check that the registered source files match the source tree. When enabled,
# CMake runs only these consistency checks.
# ##############################################################################

if(CP2K_ENABLE_CONSISTENCY_CHECKS)
  cp2k_compare_src_with_list(
    "${CP2K_FPGA_SRC_C};${CP2K_SRCS_C};${CP2K_MPIWRAP_SRCS_C};${CP2K_OFFLOAD_SRCS_C};${CP2K_DBM_SRCS_C};${CP2K_GRID_SRCS_C};${CP2K_PW_SRCS_C};${CP2K_DBM_SRCS_GPU_C};${CP2K_PROGS_C}"
    "*.c"
    "${CP2K_EXCLUDED_SRCS_C}")
  cp2k_compare_src_with_list(
    "${CP2K_SRCS_F};${CP2K_OFFLOAD_SRCS_F};${CP2K_PROGS_F}" "*.F"
    "${CP2K_EXCLUDED_SRCS_F}")
  cp2k_compare_src_with_list(
    "${CP2K_GRID_SRCS_GPU};${CP2K_GRID_SRCS_HIP};${CP2K_DBM_SRCS_GPU};${CP2K_PW_SRCS_GPU}"
    "*.cu"
    "")
  cp2k_compare_src_with_list("${CP2K_SRCS_CPP}" "*.cpp" "")
  return()
endif()

# ##############################################################################
# alias targets encapsulating cuda and hip libraries
# ##############################################################################

# create some short version of these libraries to ease things afterwards
if(CP2K_USE_CUDA)
  add_library(cp2k_cuda_libs INTERFACE)
  if(CP2K_USE_NVHPC)
    target_link_libraries(cp2k_cuda_libs INTERFACE NVHPC::MATH NVHPC::CUDA)
  else()
    target_link_libraries(
      cp2k_cuda_libs
      INTERFACE
        $<$<BOOL:${CP2K_USE_CUSOLVER_MP}>:cp2k::CUSOLVER_MP::cusolver_mp>
        CUDA::cufft CUDA::cufftw CUDA::cublas CUDA::cudart CUDA::cuda_driver)
  endif()
  if(CP2K_WITH_GPU_PROFILING)
    target_link_libraries(cp2k_cuda_libs INTERFACE CUDA::nvtx3)
  endif()
endif()

if(CP2K_USE_HIP)
  add_library(cp2k_hip_libs INTERFACE)
  target_link_libraries(
    cp2k_hip_libs
    INTERFACE hip::hipfft roc::hipblas hip::host
              $<$<STREQUAL:"${CMAKE_HIP_PLATFORM}","nvidia">:CUDA::toolkit
              CUDA::cudart CUDA::cuda_driver>)
endif()
add_library(cp2k_linalg_libs INTERFACE)
target_link_libraries(
  cp2k_linalg_libs
  INTERFACE cp2k::LAPACK::lapack
            $<$<BOOL:${CP2K_USE_OPENCL}>:OpenCL::OpenCL>
            $<$<BOOL:${CP2K_USE_OPENCL}>:libxstream::libxstream>
            $<$<BOOL:${CP2K_USE_LIBXS}>:libxs::libxs>
            $<$<BOOL:${CP2K_USE_LIBXSMM}>:libxsmm::libxsmm>
            cp2k::BLAS::blas
            $<$<BOOL:${CP2K_USE_CUDA}>:cp2k_cuda_libs>
            $<$<BOOL:${CP2K_USE_HIP}>:cp2k_hip_libs>
            $<$<BOOL:${CP2K_USE_MPI}>:MPI::MPI_Fortran>
            $<$<BOOL:${CP2K_USE_MPI}>:MPI::MPI_C>
            $<$<BOOL:${CP2K_USE_MPI}>:MPI::MPI_CXX>
            OpenMP::OpenMP_Fortran
            OpenMP::OpenMP_C
            OpenMP::OpenMP_CXX)

if(CMAKE_HIP_PLATFORM MATCHES "amd")
  set(CP2K_GPU_PLATFORM_DFLAGS "__HIP_PLATFORM_AMD__")
else()
  set(CP2K_HIP_INCLUDES "${CP2K_HIP_PATH}/include")
  set(CP2K_GPU_PLATFORM_DFLAGS "__HIP_PLATFORM_NVIDIA__")
endif()
set(CP2K_GPU_DFLAGS
    $<$<BOOL:${CP2K_USE_UNIFIED_MEMORY}>:__OFFLOAD_UNIFIED_MEMORY>
    $<$<BOOL:${CP2K_USE_OPENCL}>:__OFFLOAD_OPENCL>
    $<$<BOOL:${CP2K_USE_CUDA}>:__OFFLOAD_CUDA>
    $<$<COMPILE_LANGUAGE:CUDA>:__OFFLOAD_CUDA>
    $<$<BOOL:${CP2K_USE_HIP}>:__OFFLOAD_HIP
    ${CP2K_GPU_PLATFORM_DFLAGS}>
    $<$<COMPILE_LANGUAGE:HIP>:__OFFLOAD_HIP
    ${CP2K_GPU_PLATFORM_DFLAGS}>)

# ##############################################################################
# set properties for the cp2k target.
# ##############################################################################

add_fypp_sources(CP2K_SRCS ${CP2K_SRCS_F})
add_fypp_sources(CP2K_PROGS_F_FYPP ${CP2K_PROGS_F})

# Intel oneAPI workarounds for known compiler optimization bugs
if(CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
  set_source_files_properties(
    ${CMAKE_CURRENT_BINARY_DIR}/hfxbase/hfx_contraction_methods.F
    PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -O1")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set_source_files_properties(
      ${CMAKE_CURRENT_BINARY_DIR}/aobasis/paw_basis_types.F
      ${CMAKE_CURRENT_BINARY_DIR}/mp2_optimize_ri_basis.F
      PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -O1")
    set_source_files_properties(grid/dgemm/grid_dgemm_prepare_pab.c
                                PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS} -O2")
  endif()
endif()

# set the __SHORT_FILE__ per file for CP2K sources
foreach(cp2k_src ${CP2K_SRCS} ${CP2K_PROGS_F_FYPP})
  # add_fypp_sources returns a path in the current binary dir
  file(RELATIVE_PATH short_file "${CMAKE_BINARY_DIR}/src" "${cp2k_src}")
  set_source_files_properties(
    ${cp2k_src} PROPERTIES COMPILE_DEFINITIONS __SHORT_FILE__="${short_file}")
endforeach()
foreach(cp2k_prog ${CP2K_PROGS_F})
  set_source_files_properties(
    ${cp2k_prog} PROPERTIES COMPILE_DEFINITIONS __SHORT_FILE__="${cp2k_prog}")
endforeach()

if(BUILD_SHARED_LIBS)
  add_library(cp2k SHARED
              "${CP2K_SRCS};${CP2K_SRCS_C};${CP2K_SRCS_GPU};${CP2K_SRCS_CPP}")
else()
  add_library(cp2k STATIC
              "${CP2K_SRCS};${CP2K_SRCS_C};${CP2K_SRCS_GPU};${CP2K_SRCS_CPP}")
endif()

target_compile_features(cp2k PUBLIC cuda_std_11)
target_include_directories(cp2k PUBLIC $<INSTALL_INTERFACE:include/cp2k>)

if(TARGET cp2k_dbm_opencl_kernel)
  add_dependencies(cp2k cp2k_dbm_opencl_kernel)
endif()

set_target_properties(
  cp2k
  PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
             Fortran_MODULE_DIRECTORY mod_files
             POSITION_INDEPENDENT_CODE ON
             VERSION ${cp2k_VERSION})

if(CP2K_USE_ACCEL MATCHES "CUDA")
  set_property(TARGET cp2k PROPERTY CUDA_ARCHITECTURES
                                    "${CMAKE_CUDA_ARCHITECTURES}")
  include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
elseif(CP2K_USE_ACCEL MATCHES "HIP")
  set_property(TARGET cp2k PROPERTY HIP_ARCHITECTURES
                                    "${CMAKE_HIP_ARCHITECTURES}")
elseif(CP2K_USE_ACCEL MATCHES "OPENCL")
  include_directories(${CP2K_OPENCL_BACKEND})
  if(CP2K_OPENCL_SAMPLES)
    include_directories(${CP2K_OPENCL_SAMPLES})
  endif()
  include_directories(${OpenCL_INCLUDE_DIR})
endif()

# mix the target and variables pointing to the include directories.

if(CP2K_USE_DLAF)
  include_directories(${CP2K_DLAF_INCLUDE_DIRS})
endif()

include_directories(
  ${CMAKE_MPI_INCLUDE_DIRECTORIES}
  $<$<BOOL:${CUDAToolkit_FOUND}>:${CUDAToolkit_INCLUDE_DIRS}>
  $<$<BOOL:${CP2K_USE_HIP}>:${CMAKE_HIP_INCLUDE_DIRECTORIES}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/base>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/common>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/motion>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/dbm>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod_files>)

add_library(cp2k_link_libs INTERFACE)

# the order is important so keep it that way
target_link_libraries(
  cp2k_link_libs
  INTERFACE
    $<$<BOOL:${CP2K_USE_DFTD4}>:mctc-lib::mctc-lib>
    $<$<BOOL:${CP2K_USE_DFTD4}>:dftd4::dftd4>
    $<$<BOOL:${CP2K_USE_DEEPMD}>:DeePMD::deepmd
    DeePMD::deepmd_c>
    $<$<BOOL:${CP2K_USE_PEXSI}>:PEXSI::PEXSI>
    $<$<BOOL:${CP2K_USE_ACE}>:pace::pace>
    $<$<BOOL:${CP2K_USE_TBLITE}>:cp2k::tblite>
    $<$<BOOL:${CP2K_USE_TBLITE}>:s-dftd3::s-dftd3>
    $<$<BOOL:${CP2K_USE_SIRIUS}>:sirius::sirius>
    $<$<BOOL:${CP2K_USE_VORI}>:cp2k::VORI::vori>
    $<$<BOOL:${CP2K_USE_SPGLIB}>:Spglib::symspg>
    $<$<BOOL:${CP2K_USE_LIBSMEAGOL}>:cp2k::libsmeagol::smeagol>
    $<$<BOOL:${CP2K_USE_LIBXC}>:Libxc::xcf03>
    $<$<BOOL:${CP2K_USE_GAUXC}>:gauxc::gauxc>
    $<$<BOOL:${CP2K_USE_ELPA}>:cp2k::ELPA::elpa>
    $<$<BOOL:${CP2K_USE_PLUMED}>:cp2k::plumed::plumed>
    $<$<BOOL:${CP2K_USE_FFTW3_}>:cp2k::FFTW3::fftw3>
    $<$<AND:$<BOOL:${CP2K_USE_FFTW3_}>,$<BOOL:${CP2K_ENABLE_FFTW3_THREADS_SUPPORT}>>:cp2k::FFTW3::fftw3_threads>
    $<$<AND:$<BOOL:${CP2K_USE_FFTW3_}>,$<BOOL:${CP2K_ENABLE_FFTW3_OPENMP_SUPPORT}>>:cp2k::FFTW3::fftw3_omp>
    $<$<BOOL:${CP2K_USE_SPLA}>:SPLA::spla>
    $<$<BOOL:${CP2K_USE_MIMIC}>:CP2K::MIMIC::mclf>
    $<$<BOOL:${CP2K_USE_MIMIC}>:CP2K::MIMIC::mcl>
    $<$<BOOL:${CP2K_USE_LIBINT2}>:Libint2::int2>
    $<$<BOOL:${CP2K_USE_LIBFCI}>:libfci::fci>
    $<$<BOOL:${CP2K_USE_LIBGINT}>:cp2k::LibGint::libGint>
    $<$<BOOL:${CP2K_USE_COSMA}>:cp2k::cosma>
    $<$<BOOL:${CP2K_USE_DLAF}>:${CP2K_DLAF_FORTRAN_TARGET}>
    $<$<BOOL:${CP2K_USE_GREENX}>:cp2k::greenx>
    $<$<BOOL:${CP2K_USE_TREXIO}>:cp2k::trexio::trexio>
    $<$<BOOL:${CP2K_USE_HDF5}>:HDF5::HDF5
    hdf5::hdf5_fortran>
    wignernj::wignernj
    DBCSR::dbcsr
    $<$<BOOL:${CP2K_USE_MPI}>:cp2k::SCALAPACK::scalapack>
    cp2k_linalg_libs
    # Torch includes some LAPACK routines, but with floating-point exceptions.
    $<$<BOOL:${CP2K_USE_LIBTORCH}>:${TORCH_LIBRARIES}>
    $<$<BOOL:${CP2K_USE_OPENPMD}>:openPMD::openPMD>)

set(CP2K_BUILD_INFO_FILE "${CMAKE_CURRENT_BINARY_DIR}/cp2k_build_info.h")
file(
  WRITE "${CP2K_BUILD_INFO_FILE}"
  "#define __COMPILE_DATE \"unknown\"\n#define __COMPILE_REVISION \"unknown\"\n"
)
add_custom_target(
  cp2k_build_info
  COMMAND
    ${CMAKE_COMMAND} -DOUTPUT_FILE=${CP2K_BUILD_INFO_FILE}
    -DSOURCE_DIR=${PROJECT_SOURCE_DIR} -P
    ${PROJECT_SOURCE_DIR}/cmake/GenerateCp2kBuildInfo.cmake
  BYPRODUCTS ${CP2K_BUILD_INFO_FILE}
  COMMENT "Updating CP2K build timestamp and git revision"
  VERBATIM)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/cp2k_info.F
                            PROPERTIES OBJECT_DEPENDS "${CP2K_BUILD_INFO_FILE}")

if(CP2K_DATA_DIR STREQUAL "default")
  set(CP2K_DATA_DIR "${CMAKE_INSTALL_FULL_DATAROOTDIR}/${PROJECT_NAME}/data")
endif()

target_link_libraries(cp2k PUBLIC cp2k_link_libs)

if(CP2K_USE_ACCEL MATCHES "CUDA"
   AND CP2K_LIBXC_CUDA_SUPPORT
   AND NOT CP2K_ENABLE_LIBXC_GPU)
  target_compile_definitions(cp2k PUBLIC __NO_OFFLOAD_LIBXC)
endif()
target_compile_definitions(
  cp2k
  PUBLIC
    $<$<BOOL:${CP2K_USE_MPI}>:__parallel>
    $<$<BOOL:${CP2K_USE_MPI_F08}>:__MPI_F08>
    __COMPILE_HOST=\"${CP2K_HOST_NAME}\"
    __DATA_DIR=\"${CP2K_DATA_DIR}\"
    __COMPILE_DIR=\"${CMAKE_BINARY_DIR}\"
    __COMPILE_ARCH=\"${CMAKE_SYSTEM_PROCESSOR}\"
    $<$<BOOL:${CP2K_USE_PLUMED}>:__PLUMED2>
    $<$<BOOL:${CP2K_USE_MAXWELL}>:__LIBMAXWELL>
    $<$<BOOL:${CP2K_USE_VORI}>:__LIBVORI>
    $<$<BOOL:${CP2K_USE_SPGLIB}>:__SPGLIB>
    $<$<BOOL:${CP2K_USE_LIBSMEAGOL}>:__SMEAGOL>
    $<$<BOOL:${CP2K_USE_TBLITE}>:__TBLITE>
    $<$<BOOL:${CP2K_USE_TBLITE}>:__S_DFTD3>
    $<$<BOOL:${CP2K_USE_TBLITE}>:__DFTD4>
    $<$<BOOL:${CP2K_USE_DFTD4}>:__DFTD4>
    $<$<BOOL:${CP2K_USE_DEEPMD}>:__DEEPMD>
    $<$<BOOL:${CP2K_USE_PEXSI}>:__PEXSI>
    $<$<BOOL:${CP2K_USE_ACE}>:__ACE>
    $<$<BOOL:${CP2K_USE_SIRIUS}>:__SIRIUS>
    $<$<BOOL:${CP2K_USE_LIBVDWXC}>:__LIBVDWXC>
    $<$<BOOL:${CP2K_USE_GAUXC}>:__GAUXC>
    $<$<AND:$<BOOL:${CP2K_USE_GAUXC}>,$<BOOL:${CP2K_GAUXC_HAS_ONEDFT_GRADIENT_FIX}>>:__GAUXC_ONEDFT_GRADIENT_FIX>
    # GauXC's CMake package exports CUTLASS even when its Fortran config does
    # not.
    $<$<AND:$<BOOL:${CP2K_USE_GAUXC}>,$<BOOL:${GAUXC_HAS_CUTLASS}>>:GAUXC_HAS_CUTLASS>
    # GauXC feature flags GAUXC_HAS_MPI, GAUXC_HAS_ONEDFT, GAUXC_HAS_HOST,
    # GAUXC_HAS_DEVICE, GAUXC_HAS_CUDA, and GAUXC_HAS_HDF5 are provided by
    # gauxc/gauxc_config.f in the Fortran sources.
    $<$<BOOL:${CP2K_USE_SKALA_FTORCH}>:__SKALA_FTORCH>
    $<$<OR:$<BOOL:${CP2K_USE_GAUXC}>,$<BOOL:${CP2K_USE_SKALA_FTORCH}>>:__SKALA>
    $<$<BOOL:${CP2K_USE_SPLA}>:__SPLA>
    $<$<BOOL:${CP2K_USE_SPLA_GEMM_OFFLOADING}>:__OFFLOAD_GEMM>
    $<$<BOOL:${CP2K_USE_ELPA}>:__ELPA>
    $<$<BOOL:${CP2K_USE_DLAF}>:__DLAF>
    $<$<BOOL:${CP2K_USE_LIBXC}>:__LIBXC>
    $<$<BOOL:${CP2K_USE_HDF5}>:__HDF5>
    $<$<BOOL:${CP2K_USE_OPENPMD}>:__OPENPMD>
    $<$<BOOL:${CP2K_USE_TREXIO}>:__TREXIO>
    $<$<BOOL:${CP2K_USE_GREENX}>:__GREENX>
    $<$<OR:$<BOOL:${CP2K_USE_FFTW3_}>,$<BOOL:${CP2K_USE_FFTW3_MKL_}>>:__FFTW3>
    $<$<BOOL:${CP2K_USE_LIBINT2}>:__LIBINT>
    $<$<BOOL:${CP2K_USE_LIBGINT}>:__LIBGINT>
    $<$<BOOL:${CP2K_USE_LIBFCI}>:__LIBFCI>
    $<$<BOOL:${CP2K_USE_LIBTORCH}>:__LIBTORCH>
    $<$<BOOL:${CP2K_USE_LIBTORCH_CUDA}>:__LIBTORCH_CUDA>
    $<$<BOOL:${CP2K_USE_COSMA}>:__COSMA>
    $<$<BOOL:${CP2K_USE_LIBXS}>:__LIBXS>
    $<$<BOOL:${CP2K_USE_LIBXS}>:__BLAS>
    $<$<BOOL:${CP2K_USE_LIBXSMM}>:__LIBXSMM>
    $<$<BOOL:${CP2K_USE_OPENCL}>:__LIBXSTREAM>
    $<$<BOOL:${CP2K_USE_SIRIUS_DFTD3}>:__SIRIUS_DFTD3>
    $<$<BOOL:${CP2K_USE_SIRIUS_DFTD4}>:__SIRIUS_DFTD4>
    $<$<BOOL:${CP2K_USE_SIRIUS_NLCG}>:__SIRIUS_NLCG>
    $<$<BOOL:${CP2K_USE_SIRIUS_VCSQNM}>:__SIRIUS_VCSQNM>
    $<$<BOOL:${CP2K_USE_MIMIC}>:__MIMIC>
    $<$<STREQUAL:"${CP2K_BLAS_VENDOR}","MKL">:__MKL>
    $<$<STREQUAL:"${CP2K_BLAS_VENDOR}","Apple">:__ACCELERATE>
    $<$<BOOL:${CP2K_USE_CUSOLVER_MP}>:__CUSOLVERMP>
    $<$<BOOL:${CP2K_CUSOLVERMP_USE_NCCL}>:__CUSOLVERMP_NCCL>
    $<$<AND:$<STREQUAL:"${CP2K_BLAS_VENDOR}","OpenBLAS">,$<STREQUAL:"${CP2K_BLAS_INTERFACE}","32bits">>:__OPENBLAS>
    $<$<BOOL:${CP2K_USE_CRAY_PM_ENERGY}>:__CRAY_PM_ENERGY>
    $<$<BOOL:${CP2K_USE_CRAY_PM_ACCEL_ENERGY}>:__CRAY_PM_ACCEL_ENERGY>
    ${CP2K_GPU_DFLAGS})

if(CP2K_USE_CUDA
   OR CP2K_USE_HIP
   OR CP2K_USE_OPENCL)
  target_compile_definitions(
    cp2k
    PUBLIC $<$<BOOL:${CP2K_WITH_GPU_PROFILING}>:__OFFLOAD_PROFILING>
           $<$<NOT:$<BOOL:${CP2K_DBCSR_USE_CPU_ONLY}>>:__DBCSR_ACC>
           $<$<NOT:$<BOOL:${CP2K_ENABLE_ELPA_GPU}>>:__NO_OFFLOAD_ELPA>
           $<$<NOT:$<BOOL:${CP2K_ENABLE_GRID_GPU}>>:__NO_OFFLOAD_GRID>
           $<$<NOT:$<BOOL:${CP2K_ENABLE_DBM_GPU}>>:__NO_OFFLOAD_DBM>
           $<$<NOT:$<BOOL:${CP2K_ENABLE_PW_GPU}>>:__NO_OFFLOAD_PW>)
endif()

add_dependencies(cp2k cp2k_build_info)

# add an alias for the FetchContent. Note that FindPackage will automatically
# have this alias available.

add_library(cp2k::cp2k ALIAS cp2k)

if(CP2K_USE_LIBXS)
  message(STATUS "Adding libxs_jit.F from dependency libxs for compilation")
  get_target_property(_libxs_incdirs libxs::libxs INTERFACE_INCLUDE_DIRECTORIES)
  find_path(CP2K_LIBXS_JIT_DIR "libxs_jit.F" HINTS ${_libxs_incdirs} REQUIRED)
  target_sources(cp2k PRIVATE "${CP2K_LIBXS_JIT_DIR}/libxs_jit.F")
  unset(_libxs_incdirs)
endif()

# ##############################################################################
# applications
# ##############################################################################

list(
  APPEND
  __CP2K_APPS
  "orbital_transformation_matrices_unittest"
  "gemm_square_unittest"
  "memory_utilities_unittest"
  "cell_opt_utils_unittest"
  "parallel_rng_types_unittest"
  "upf_inverse_unittest"
  "atom_upf_unittest"
  "xc_semilocal_unittest"
  "upf_projector_integrals_unittest"
  "upf_local_integrals_unittest"
  "gapw_unittest"
  "kpoint_lattice_fft_unittest"
  "kpoint_smearing_unittest"
  "topology_wilson_unittest"
  "kubo_unittest"
  "spectral_localizer_unittest"
  "quadratic_pseudospectrum_unittest"
  "pw_unittest"
  "topology_symmetry_unittest"
  "topology_curvature_unittest"
  "topology_snapshot_unittest"
  "topology_phasons"
  "kpsym_unittest"
  "qs_ot_minimizer_unittest"
  "optbas_frontier_orbitals_utils_unittest"
  "low_rank_preconditioner_unittest"
  "ot_covariant_preconditioner_unittest"
  "qs_charge_mixing_unittest"
  "qs_fod_unittest"
  "qs_ot_complex_ref_unittest"
  "wannier90_projections_unittest"
  "cell_types_unittest"
  "hfx_pair_list_methods_unittest"
  "graph"
  "dumpdcd"
  "xyz2dcd"
  "libcp2k_unittest"
  "nequip_unittest"
  "dbt_unittest"
  "dbt_tas_unittest"
  "gx_ac_unittest")

if(CP2K_USE_MPI)
  list(APPEND __CP2K_APPS "message_passing_thread_unittest"
       "cp_fm_diag_utils_unittest")
endif()

add_executable(cp2k-bin start/cp2k.F)
# for linking

add_executable(orbital_transformation_matrices_unittest
               aobasis/orbital_transformation_matrices_unittest.F)
add_executable(gemm_square_unittest common/gemm_square_unittest.F)
add_executable(memory_utilities_unittest common/memory_utilities_unittest.F)
add_executable(parallel_rng_types_unittest common/parallel_rng_types_unittest.F)
add_executable(upf_inverse_unittest pseudopotentials/upf_inverse_unittest.F)
add_executable(atom_upf_unittest pseudopotentials/atom_upf_unittest.F)
add_executable(xc_semilocal_unittest xc/xc_semilocal_unittest.F)
add_executable(upf_projector_integrals_unittest
               pseudopotentials/upf_projector_integrals_unittest.F)
add_executable(upf_local_integrals_unittest
               pseudopotentials/upf_local_integrals_unittest.F)
add_executable(gapw_unittest gapw_unittest.F)
add_executable(kpoint_lattice_fft_unittest kpoint_lattice_fft_unittest.F)
add_executable(kpoint_smearing_unittest kpoint_smearing_unittest.F)
add_executable(topology_wilson_unittest topology_wilson_unittest.F)
add_executable(kubo_unittest kubo_unittest.F)
add_executable(spectral_localizer_unittest spectral_localizer_unittest.F)
add_executable(quadratic_pseudospectrum_unittest
               quadratic_pseudospectrum_unittest.F)
add_executable(pw_unittest pw/pw_unittest.F)
add_executable(topology_symmetry_unittest topology_symmetry_unittest.F)
add_executable(topology_curvature_unittest topology_curvature_unittest.F)
add_executable(topology_snapshot_unittest topology_snapshot_unittest.F)
add_executable(topology_phasons topology_phasons.F)
add_executable(kpsym_unittest kpsym_unittest.F)
add_executable(qs_ot_minimizer_unittest qs_ot_minimizer_unittest.F)
add_executable(optbas_frontier_orbitals_utils_unittest
               optbas_frontier_orbitals_utils_unittest.F)
add_executable(low_rank_preconditioner_unittest
               low_rank_preconditioner_unittest.F)
add_executable(ot_covariant_preconditioner_unittest
               ot_covariant_preconditioner_unittest.F)
add_executable(qs_charge_mixing_unittest qs_charge_mixing_unittest.F)
add_executable(qs_fod_unittest qs_fod_unittest.F)
add_executable(qs_ot_complex_ref_unittest qs_ot_complex_ref_unittest.F)
add_executable(wannier90_projections_unittest wannier90_projections_unittest.F)
add_executable(cell_types_unittest subsys/cell_types_unittest.F)
add_executable(hfx_pair_list_methods_unittest hfx_pair_list_methods_unittest.F)
add_executable(cell_opt_utils_unittest motion/cell_opt_utils_unittest.F)
add_executable(graph metadyn_tools/graph.F)
add_executable(dumpdcd motion/dumpdcd.F)
add_executable(xyz2dcd motion/xyz2dcd.F)
add_executable(libcp2k_unittest start/libcp2k_unittest.c)
add_executable(nequip_unittest nequip_unittest.F)
add_executable(dbt_unittest dbt/dbt_unittest.F)
add_executable(dbt_tas_unittest dbt/tas/dbt_tas_unittest.F)
add_executable(gx_ac_unittest gx_ac_unittest.F)
if(CP2K_USE_MPI)
  add_executable(message_passing_thread_unittest
                 mpiwrap/message_passing_thread_unittest.F)
  add_executable(cp_fm_diag_utils_unittest fm/cp_fm_diag_utils_unittest.F)
endif()

set_target_properties(
  cp2k-bin
  PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
             POSITION_INDEPENDENT_CODE ON
             LINKER_LANGUAGE "Fortran"
             OUTPUT_NAME "cp2k.${__cp2k_ext}")

target_link_libraries(cp2k-bin PRIVATE stdc++)
target_link_libraries(cp2k-bin PUBLIC cp2k)

foreach(_app ${__CP2K_APPS})
  # libcp2k_unittest is a C main program; keep linking it with CXX
  if(_app STREQUAL "libcp2k_unittest")
    set_target_properties(
      ${_app}
      PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
                 POSITION_INDEPENDENT_CODE ON
                 LINKER_LANGUAGE "CXX"
                 OUTPUT_NAME "${_app}.${__cp2k_ext}")
  else()
    set_target_properties(
      ${_app}
      PROPERTIES IMPORTED_LINK_INTERFACE_LANGUAGES "CXX"
                 POSITION_INDEPENDENT_CODE ON
                 LINKER_LANGUAGE "Fortran"
                 OUTPUT_NAME "${_app}.${__cp2k_ext}")
    target_link_libraries(${_app} PRIVATE stdc++)
  endif()

  if(CP2K_USE_ACCEL MATCHES HIP)
    set_property(TARGET ${_app} PROPERTY HIP_ARCHITECTURES
                                         "${CMAKE_HIP_ARCHITECTURES}")
  endif()

  if(CP2K_USE_ACCEL MATCHES "CUDA"
     AND CP2K_LIBXC_CUDA_SUPPORT
     AND NOT CP2K_ENABLE_LIBXC_GPU)
    target_compile_definitions(${_app} PUBLIC __NO_OFFLOAD_LIBXC)
  endif()

  target_link_libraries(${_app} PUBLIC cp2k)
endforeach()

# ##############################################################################
# different miniapps
# ##############################################################################

add_executable(
  dbm_miniapp dbm/dbm_miniapp.c ${CP2K_DBM_SRCS_C} ${CP2K_MPIWRAP_SRCS_C}
              ${CP2K_DBM_SRCS_GPU} ${CP2K_OFFLOAD_SRCS_C})
if(TARGET cp2k_dbm_opencl_kernel)
  add_dependencies(dbm_miniapp cp2k_dbm_opencl_kernel)
endif()
target_compile_definitions(
  dbm_miniapp
  PRIVATE $<$<BOOL:${CP2K_USE_LIBXS}>:__LIBXS>
          $<$<BOOL:${CP2K_USE_LIBXS}>:__BLAS>
          $<$<BOOL:${CP2K_USE_LIBXSMM}>:__LIBXSMM>
          $<$<STREQUAL:"${CP2K_BLAS_VENDOR}","MKL">:__MKL>)

if(CP2K_DISABLE_DBM_GPU)
  target_compile_definitions(dbm_miniapp PRIVATE __NO_OFFLOAD_DBM)
else()
  target_compile_definitions(dbm_miniapp PRIVATE ${CP2K_GPU_DFLAGS})
endif()

foreach(__app grid_miniapp grid_unittest)
  add_executable(
    ${__app}
    "grid/${__app}.c"
    ${CP2K_GRID_SRCS_C}
    $<$<OR:$<BOOL:${CP2K_USE_CUDA}>,$<BOOL:${CP2K_USE_HIP}>>:${CP2K_GRID_SRCS_GPU}>
    $<$<BOOL:${CP2K_USE_HIP}>:${CP2K_GRID_SRCS_HIP}>
    ${CP2K_MPIWRAP_SRCS_C}
    ${CP2K_OFFLOAD_SRCS_C})
  target_compile_definitions(
    ${__app}
    PRIVATE $<$<BOOL:${CP2K_USE_LIBXS}>:__LIBXS>
            $<$<BOOL:${CP2K_USE_LIBXS}>:__BLAS>
            $<$<BOOL:${CP2K_USE_LIBXSMM}>:__LIBXSMM>
            $<$<STREQUAL:"${CP2K_BLAS_VENDOR}","MKL">:__MKL>)

  if(CP2K_DISABLE_GRID_GPU)
    target_compile_definitions(${__app} PUBLIC __NO_OFFLOAD_GRID)
  else()
    target_compile_definitions(${__app} PUBLIC ${CP2K_GPU_DFLAGS})
  endif()
endforeach()

# ##############################################################################
# common to all c executables
# ##############################################################################

# list of all targets that share the same properties. simplify code as most
# cmake function can work on list of targets

foreach(__app grid_miniapp grid_unittest dbm_miniapp)
  set_target_properties(
    ${__app}
    PROPERTIES POSITION_INDEPENDENT_CODE ON
               LINKER_LANGUAGE CXX
               POSITION_INDEPENDENT_CODE ON
               OUTPUT_NAME "${__app}.${__cp2k_ext}")

  if(CP2K_USE_HIP)
    set_property(TARGET ${__app} PROPERTY HIP_ARCHITECTURES
                                          "${CMAKE_HIP_ARCHITECTURES}")
  endif()

  if(CP2K_USE_CUDA)
    set_target_properties(${__app} PROPERTIES CUDA_ARCHITECTURES
                                              "${CMAKE_CUDA_ARCHITECTURES}")
  endif()

  target_compile_features(${__app} PUBLIC cuda_std_14 c_std_99 cxx_std_14)

  target_include_directories(
    ${__app}
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
           $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/base>
           $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/common>
           $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
           $<$<BOOL:${CP2K_USE_CUDA}>:${CMAKE_CUDA_INCLUDE_DIRECTORIES}>
           $<$<BOOL:${CP2K_USE_HIP}>:${CMAKE_HIP_INCLUDE_DIRECTORIES}>)

  target_link_libraries(${__app} PUBLIC cp2k_linalg_libs m)
endforeach()

# ##############################################################################
# installation
# ##############################################################################

set(__CP2K_INSTALL_APPS ${__CP2K_APPS} cp2k-bin grid_miniapp grid_unittest
                        dbm_miniapp)
if(NOT CP2K_INSTALL_UNITTESTS)
  list(FILTER __CP2K_INSTALL_APPS EXCLUDE REGEX "_unittest$")
endif()
install(TARGETS ${__CP2K_INSTALL_APPS})

install(DIRECTORY "${CMAKE_SOURCE_DIR}/data/" DESTINATION "${CP2K_DATA_DIR}")

# ##############################################################################
# cp2kConfig.cmake, etc...
# ##############################################################################

set(cp2k_export_libs cp2k cp2k_link_libs cp2k_linalg_libs)
if(CP2K_USE_CUDA)
  list(APPEND cp2k_export_libs cp2k_cuda_libs)
endif()

if(CP2K_USE_HIP)
  list(APPEND cp2k_export_libs cp2k_hip_libs)
endif()

install(
  TARGETS ${cp2k_export_libs}
  EXPORT cp2k_targets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")

install(
  EXPORT cp2k_targets
  FILE cp2kTargets.cmake
  NAMESPACE cp2k::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

install(FILES start/libcp2k.h
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}")

if(NOT CMAKE_INSTALL_Fortran_MODULES)
  set(CMAKE_INSTALL_Fortran_MODULES
      "${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}/mod_files"
  )
endif()

install(
  DIRECTORY "${PROJECT_BINARY_DIR}/src/mod_files/"
  DESTINATION "${CMAKE_INSTALL_Fortran_MODULES}"
  FILES_MATCHING
  PATTERN "*.mod")
