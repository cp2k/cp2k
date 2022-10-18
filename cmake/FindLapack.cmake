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

# check for blas first. Most of the vendor libraries bundle lapack and blas in
# the same library. (MKL, and OPENBLAS do)

find_package(PkgConfig)
find_package(Blas REQUIRED)

if(CP2K_BLAS_FOUND)
  # LAPACK in the Intel MKL 10+ library?
  if(CP2K_BLAS_VENDOR MATCHES "MKL"
     OR CP2K_BLAS_VENDOR MATCHES "OpenBLAS"
     OR CP2K_BLAS_VENDOR MATCHES "Arm"
     OR CP2K_BLAS_VENDOR MATCHES "SCI"
     OR CP2K_BLAS_VENDOR MATCHES "FlexiBLAS")
    # we just need to create the interface that's all
    set(CP2K_LAPACK_FOUND TRUE)
    get_target_property(CP2K_LAPACK_INCLUDE_DIRS CP2K_BLAS::blas
                        INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(CP2K_LAPACK_LIBRARIES CP2K_BLAS::blas
                        INTERFACE_LINK_LIBRARIES)
  else()

    # we might get lucky to find a pkgconfig package for lapack (fedora provides
    # one for instance)
    if(PKG_CONFIG_FOUND)
      pkg_check_modules(CP2K_LAPACK lapack)
    endif()

    if(NOT CP2K_LAPACK_FOUND)
      find_library(
        CP2K_LAPACK_LIBRARIES
        NAMES "lapack" "lapack64"
        PATH_SUFFIXES "openblas" "openblas64" "openblas-pthread"
                      "openblas-openmp" "lib" "lib64"
        NO_DEFAULT_PATH)
    endif()
  endif()
endif()

# check if found
find_package_handle_standard_args(Lapack REQUIRED_VARS CP2K_LAPACK_LIBRARIES)

if(CP2K_LAPACK_FOUND)
  if(NOT TARGET CP2K_LAPACK::lapack)
    add_library(CP2K_LAPACK::lapack INTERFACE IMPORTED)
    set_property(TARGET CP2K_LAPACK::lapack PROPERTY INTERFACE_LINK_LIBRARIES
                                                     ${CP2K_LAPACK_LIBRARIES})
    if(CP2K_LAPACK_INCLUDE_DIRS)
      set_property(
        TARGET CP2K_LAPACK::lapack PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                            ${CP2K_LAPACK_INCLUDE_DIRS})
    endif()
    add_library(CP2K_LAPACK::LAPACK ALIAS CP2K_LAPACK::lapack)
  endif()
endif()

# prevent clutter in cache
mark_as_advanced(CP2K_LAPACK_FOUND CP2K_LAPACK_LIBRARIES
                 CP2K_LAPACK_INCLUDE_DIRS)
