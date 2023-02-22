#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: BSD-3-Clause                                                         !
#!-------------------------------------------------------------------------------------------------!

set(GRID_SRCS_C
    common/grid_basis_set.c
    common/grid_library.c
    common/grid_sphere_cache.c
    cpu/coefficients.c
    cpu/collocation_integration.c
    cpu/grid_collocate_dgemm.c
    cpu/grid_context_cpu.c
    cpu/grid_integrate_dgemm.c
    cpu/grid_prepare_pab_dgemm.c
    cpu/non_orthorombic_corrections.c
    cpu/tensor_local.c
    cpu/utils.c
    ref/grid_ref_collocate.c
    ref/grid_ref_integrate.c
    ref/grid_ref_prepare_pab.c
    ref/grid_ref_task_list.c
    grid_replay.c
    grid_task_list.c)

if(CP2K_USE_ACCEL MATCHES "HIP|CUDA")
  set(GRID_SRCS_GPU gpu/grid_gpu_collocate.cu gpu/grid_gpu_integrate.cu
                    gpu/grid_gpu_task_list.cu)

  if(CP2K_USE_ACCEL MATCHES "HIP")
    list(
      APPEND
      GRID_SRCS_GPU
      "hip/grid_hip_collocate.cu;hip/grid_hip_context.cu;hip/grid_hip_integrate.cu"
    )
    set_source_files_properties(${GRID_SRCS_GPU} PROPERTIES LANGUAGE HIP)
  endif()
endif()

if(CP2K_DISABLE_GRID_GPU)
  target_compile_definitions(cp2k_grid PUBLIC __NO_OFFLOAD_GRID)
endif()

add_library(cp2k_grid "${GRID_SRCS_C};${GRID_SRCS_GPU}")

set_property(TARGET cp2k_grid PROPERTY CUDA_ARCHITECTURES
                                       ${CP2K_ACC_ARCH_NUMBER})
set_property(TARGET cp2k_grid PROPERTY HIP_ARCHITECTURES
                                       ${CP2K_ACC_ARCH_NUMBER})
set_target_properties(cp2k_grid PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(cp2k_grid PUBLIC cuda_std_14)
target_compile_features(cp2k_grid PUBLIC c_std_99)
target_compile_features(cp2k_grid PUBLIC cxx_std_14)

target_link_libraries(
  cp2k_grid
  PUBLIC cp2k_offload
         $<$<BOOL:${CP2K_USE_MPI}>:MPI::MPI_C>
         $<$<BOOL:${CP2K_USE_LIBXSMM}>:CP2K::LibXSMM::libxsmmf>
         $<$<BOOL:${CP2K_USE_LIBXSMM}>:CP2K::LibXSMM::libxsmmext>
         $<$<BOOL:${CP2K_USE_HIP}>:hip::host>
         OpenMP::OpenMP_C
         CP2K::BLAS::blas
         m)

target_compile_definitions(
  cp2k_grid
  PRIVATE $<$<CONFIG:Release>:NDEBUG>
          $<$<BOOL:${CP2K_USE_HIP}>:__HIP_PLATFORM_AMD__
          __OFFLOAD_HIP>
          $<$<COMPILE_LANGUAGE:HIP>:__HIP_PLATFORM_AMD__
          __OFFLOAD_HIP>
          $<$<BOOL:${CP2K_USE_CUDA}>:__OFFLOAD_CUDA>
          $<$<COMPILE_LANGUAGE:CUDA>:__OFFLOAD_CUDA>
          $<$<STREQUAL:${CP2K_BLAS_VENDOR},"MKL">:__MKL>
          $<$<BOOL:${CP2K_USE_LIBXSMM}>:__LIBXSMM>
          $<$<BOOL:${CP2K_DISABLE_GRID_GPU}>:__NO_OFFLOAD_GRID>)
target_include_directories(
  cp2k_grid
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
         $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/base>
         $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/common>
         $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
         $<$<BOOL:${CP2K_USE_CUDA}>:${CMAKE_CUDA_INCLUDE_DIRECTORIES}>
         $<$<BOOL:${CP2K_USE_HIP}>:${CMAKE_HIP_INCLUDE_DIRECTORIES}>)

foreach(__app grid_miniapp grid_unittest)
  add_executable(${__app} ${__app}.c)
  set_target_properties(${__app} PROPERTIES LINKER_LANGUAGE C)
  target_link_libraries(${__app} PUBLIC cp2k_grid)
endforeach()
install(
  TARGETS grid_miniapp grid_unittest
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  INCLUDES
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/cp2k")
install(TARGETS cp2k_grid LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
