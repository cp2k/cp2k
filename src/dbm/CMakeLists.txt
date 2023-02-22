#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: BSD-3-Clause                                                         !
#!-------------------------------------------------------------------------------------------------!

set(DBM_SRCS_C
    dbm_distribution.c
    dbm_library.c
    dbm_matrix.c
    dbm_mempool.c
    dbm_miniapp.c
    dbm_mpi.c
    dbm_multiply.c
    dbm_multiply_comm.c
    dbm_multiply_cpu.c
    dbm_multiply_gpu.c
    dbm_shard.c)

set(DBM_SRCS_GPU "")

if(CP2K_USE_ACCEL MATCHES "CUDA|HIP")
  if(NOT CP2K_DISABLE_DBM_GPU)
    set(DBM_SRCS_GPU dbm_multiply_gpu_kernel.cu)
  endif()
endif()

add_library(cp2k_dbm ${DBM_SRCS_C};${DBM_SRCS_GPU})

if(CP2K_DISABLE_DBM_GPU)
  target_compile_definitions(cp2k_dbm __NO_OFFLOAD_DBM)
endif()

target_compile_features(cp2k_dbm PUBLIC cxx_std_14)
target_compile_features(cp2k_dbm PUBLIC c_std_99)
target_compile_features(cp2k_dbm PUBLIC cuda_std_11)

if(CP2K_USE_HIP)
  set_target_properties(cp2k_dbm PROPERTIES HIP_ARCHITECTURES
                                            ${CP2K_ACC_ARCH_NUMBER})
  set_source_files_properties(${DBM_SRCS_GPU} PROPERTIES LANGUAGE HIP)
endif()

if(CP2K_USE_CUDA)
  set_target_properties(cp2k_dbm PROPERTIES CUDA_ARCHITECTURES
                                            ${CP2K_ACC_ARCH_NUMBER})
endif()

add_executable(dbm_miniapp dbm_miniapp.c)
set_target_properties(dbm_miniapp PROPERTIES LINKER_LANGUAGE C)
target_compile_features(dbm_miniapp PUBLIC c_std_99)

target_compile_definitions(
  cp2k_dbm
  PUBLIC
    $<$<BOOL:${CP2K_DISABLE_DBM_GPU}>:
    __NO_OFFLOAD_DBM>
    $<$<AND:$<NOT:$<BOOL:${CP2K_DISABLE_DBM_GPU}>>,$<BOOL:${CP2K_USE_CUDA}>>:__OFFLOAD_CUDA>
    $<$<AND:$<NOT:$<BOOL:${CP2K_DISABLE_DBM_GPU}>>,$<BOOL:${CP2K_USE_HIP}>>:__OFFLOAD_HIP>
    $<$<BOOL:${CP2K_USE_LIBXSMM}>:__LIBXSMM>)

target_link_libraries(
  cp2k_dbm
  PUBLIC
    cp2k_offload
    $<$<AND:$<NOT:$<BOOL:${CP2K_DISABLE_DBM_GPU}>>,$<BOOL:${CP2K_USE_CUDA}>>:CUDA::cudart>
    $<$<AND:$<NOT:$<BOOL:${CP2K_DISABLE_DBM_GPU}>>,$<BOOL:${CP2K_USE_HIP}>>:hip::host>
    $<$<BOOL:${CP2K_USE_LIBXSMM}>:CP2K::LibXSMM::libxsmm>
    $<$<BOOL:${CP2K_USE_LIBXSMM}>:CP2K::LibXSMM::libxsmmext>
    $<$<BOOL:${CP2K_USE_MPI}>:MPI::MPI_C>
    CP2K::BLAS::blas
    OpenMP::OpenMP_C
    m)

target_link_libraries(dbm_miniapp PUBLIC cp2k_dbm)

target_compile_definitions(
  dbm_miniapp
  PUBLIC
    $<$<BOOL:${CP2K_DISABLE_DBM_GPU}>:
    __NO_OFFLOAD_DBM>
    $<$<AND:$<NOT:$<BOOL:${CP2K_DISABLE_DBM_GPU}>>,$<BOOL:${CP2K_USE_CUDA}>>:__OFFLOAD_CUDA>
    $<$<AND:$<NOT:$<BOOL:${CP2K_DISABLE_DBM_GPU}>>,$<BOOL:${CP2K_USE_HIP}>>:__OFFLOAD_HIP>
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
