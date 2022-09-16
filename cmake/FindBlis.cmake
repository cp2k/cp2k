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

# .rst: FindBLIS
# -----------
#
# This module tries to find the BLIS library.
#
# The following variables are set
#
# ::
#
# BLIS_FOUND           - True if blis is found BLIS_LIBRARIES       - The
# required libraries BLIS_INCLUDE_DIRS    - The required include directory
#
# The following import target is created
#
# ::
#
# BLIS::blis

# set paths to look for library from ROOT variables.If new policy is set,
# find_library() automatically uses them.

find_package(PkgConfig)

if(DEFINED BLIS_ROOT)
  list(APPEND _BLIS_PATHS ${BLIS_ROOT} $ENV{BLIS_ROOT})
endif()

if(DEFINED AOCL_ROOT)
  list(APPEND _BLIS_PATHS ${AOCL_ROOT} $ENV{AOCL_ROOT})
endif()

# one day blis will have a pkg-config file
if(PKG_CONFIG_FOUND)
  pkg_check_modules(BLIS blis)
endif()

if(NOT BLIS_FOUND)
  find_library(
    BLIS_LIBRARIES
    NAMES "blis-mt" "blis"
    HINTS ${_BLIS_PATHS}
    PATH_SUFFIXES "blis/lib" "blis/lib64" "blis")
endif()

find_path(
  BLIS_INCLUDE_DIRS
  NAMES "blis.h"
  HINTS ${_BLIS_PATHS}
  PATH_SUFFIXES "blis" "blis/include" "include/blis")

# check if found
include(FindPackageHandleStandardArgs)
if(BLIS_INCLUDE_DIRS)
  find_package_handle_standard_args(BLIS REQUIRED_VARS BLIS_INCLUDE_DIRS
                                                       BLIS_LIBRARIES)
else()
  find_package_handle_standard_args(BLIS REQUIRED_VARS BLIS_LIBRARIES)
endif()

if(NOT BLIS_FOUND)
  set(BLIS_FOUND ON)
endif()

# add target to link against
if(BLIS_FOUND AND NOT TARGET BLIS::blis)
  add_library(BLIS::blis INTERFACE IMPORTED)
endif()
set_property(TARGET BLIS::blis PROPERTY INTERFACE_LINK_LIBRARIES
                                        ${BLIS_LIBRARIES})
if(BLIS_INCLUDE_DIRS)
  set_property(TARGET BLIS::blis PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                          ${BLIS_INCLUDE_DIRS})
endif()
endif()
# prevent clutter in cache
mark_as_advanced(BLIS_FOUND BLIS_LIBRARIES BLIS_INCLUDE_DIRS)
