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

# .rst: FindGenericBLAS
# -----------
#
# This module tries to find the GenericBLAS library.
#
# The following variables are set
#
# ::
#
# GenericBLAS_FOUND           - True if blas is found GenericBLAS_LIBRARIES -
# The required libraries GenericBLAS_INCLUDE_DIRS    - The required include
# directory
#
# The following import target is created
#
# ::
#
# GenericBLAS::blas

# set paths to look for library from ROOT variables.If new policy is set,
# find_library() automatically uses them.
if(NOT POLICY CMP0074)
  set(_GenericBLAS_PATHS ${GenericBLAS_ROOT} $ENV{GenericBLAS_ROOT})
endif()

find_library(
  GenericBLAS_LIBRARIES
  NAMES "blas"
  HINTS ${_GenericBLAS_PATHS})
find_library(
  # optinally look for cblas library - not required
  GenericBLAS_CBLAS_LIBRARIES
  NAMES "cblas"
  HINTS ${_GenericBLAS_PATHS})
find_path(
  GenericBLAS_INCLUDE_DIRS
  NAMES "cblas.h"
  HINTS ${_GenericBLAS_PATHS})

# check if found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  GenericBLAS REQUIRED_VARS GenericBLAS_INCLUDE_DIRS GenericBLAS_LIBRARIES)
if(GenericBLAS_CBLAS_LIBRARIES)
  list(APPEND GenericBLAS_LIBRARIES ${GenericBLAS_CBLAS_LIBRARIES})
endif()

# add target to link against
if(GenericBLAS_FOUND)
  if(NOT TARGET GenericBLAS::blas)
    add_library(GenericBLAS::blas INTERFACE IMPORTED)
  endif()
  set_property(TARGET GenericBLAS::blas PROPERTY INTERFACE_LINK_LIBRARIES
                                                 ${GenericBLAS_LIBRARIES})
  set_property(TARGET GenericBLAS::blas PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                                 ${GenericBLAS_INCLUDE_DIRS})
endif()

# prevent clutter in cache
mark_as_advanced(GenericBLAS_FOUND GenericBLAS_LIBRARIES
                 GenericBLAS_INCLUDE_DIRS GenericBLAS_CBLAS_LIBRARIES)
