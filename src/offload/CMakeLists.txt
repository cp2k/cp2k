#!-------------------------------------------------------------------------------------------------!
#!   CP2K: A general program to perform molecular dynamics simulations                             !
#!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                  !
#!                                                                                                 !
#!   SPDX-License-Identifier: BSD-3-Clause                                                         !
#!-------------------------------------------------------------------------------------------------!

add_library(cp2k_offload offload_buffer.c offload_library.c)

target_link_libraries(
  cp2k_offload PUBLIC $<$<BOOL:${CP2K_USE_CUDA}>:CUDA::cudart>
                      $<$<BOOL:${CP2K_USE_HIP}>:hip::host>)

target_compile_definitions(
  cp2k_offload PUBLIC $<$<BOOL:${CP2K_USE_CUDA}>:__OFFLOAD_CUDA>
                      $<$<BOOL:${CP2K_USE_HIP}>:__OFFLOAD_HIP>)

install(TARGETS cp2k_offload LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")
