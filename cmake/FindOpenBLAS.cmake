# Copyright (c) 2019 ETH Zurich
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

# .rst: FindOPENBLAS
# -----------
#
# This module tries to find the OPENBLAS library.
#
# The following variables are set
#
# ::
#
# OPENBLAS_FOUND           - True if openblas is found OPENBLAS_LIBRARIES - The
# required libraries OPENBLAS_INCLUDE_DIRS    - The required include directory
#
# The following import target is created
#
# ::
#
# OpenBLAS::openblas

# set paths to look for library from ROOT variables.If new policy is set,
# find_library() automatically uses them.

include(cp2k_utils)
include(FindPackageHandleStandardArgs)

find_package(PkgConfig)

cp2k_set_default_paths(OPENBLAS)

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_OPENBLAS openblas)
endif()

# try the openblas module of openblas library Maybe we are lucky it is installed
# find_package(OPENBLAS QUIET)

if(NOT CP2K_OPENBLAS_FOUND)
  cp2k_find_libraries(OPENBLAS "openblas")
  cp2k_find_libraries(OPENBLAS64 "openblas64")
  cp2k_find_libraries(OPENBLAS_THREADS "openblas_threads")
  cp2k_find_libraries(OPENBLAS_THREADS "openblas64_threads")
  cp2k_find_libraries(OPENBLAS_THREADS "openblas64_omp")
  cp2k_find_libraries(OPENBLAS_THREADS "openblas_omp")
endif()

cp2k_include_dirs(OPENBLAS "cblas.h")

# check if found
find_package_handle_standard_args(
  OpenBLAS REQUIRED_VARS CP2K_OPENBLAS_INCLUDE_DIRS
                         CP2K_OPENBLAS_LINK_LIBRARIES)

# add target to link against
if(CP2K_OPENBLAS_FOUND AND NOT TARGET CP2K_OpenBLAS::openblas)
  add_library(CP2K_OpenBLAS::openblas INTERFACE IMPORTED)
  set_property(
    TARGET CP2K_OpenBLAS::openblas PROPERTY INTERFACE_LINK_LIBRARIES
                                            ${CP2K_OPENBLAS_LINK_LIBRARIES})
  if(CP2K_OPENBLAS_INCLUDE_DIRS)
    set_property(
      TARGET CP2K_OpenBLAS::openblas PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                              ${CP2K_OPENBLAS_INCLUDE_DIRS})
  endif()
  add_library(CP2K_OpenBLAS::blas ALIAS CP2K_OpenBLAS::openblas)
  set(CP2K_BLAS_VENDOR "OpenBLAS")
endif()

# prevent clutter in cache
mark_as_advanced(CP2K_BLAS_VENDOR CP2K_OPENBLAS_FOUND
                 CP2K_OPENBLAS_LINK_LIBRARIES CP2K_OPENBLAS_INCLUDE_DIRS)
