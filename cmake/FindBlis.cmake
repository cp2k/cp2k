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

find_package(PkgConfig)
include(cp2k_utils)
include(FindPackageHandleStandardArgs)

cp2k_set_default_paths(BLIS "BLIS")

if(DEFINED AOCL_ROOT)
  list(CP2K_BLIS_PREFIX "${AOCL_ROOT}" "$ENV{AOCL_ROOT}")
endif()

# one day blis will have a pkg-config file
if(PKG_CONFIG_FOUND)
  pkg_check_modules(BLIS IMPORTED_TARGET GLOBAL blis)
endif()

if(NOT CP2K_BLIS_FOUND)
  cp2k_find_libraries(BLIS "blis")
endif()

if(NOT CP2K_BLIS_INCLUDE_DIRS)
  cp2k_include_dirs(BLIS "blis.h")
endif()

# check if found
if(CP2K_BLIS_INCLUDE_DIRS)
  find_package_handle_standard_args(
    BLIS REQUIRED_VARS CP2K_BLIS_FOUND CP2K_BLIS_INCLUDE_DIRS
                       CP2K_BLIS_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(BLIS REQUIRED_VARS CP2K_BLIS_FOUND
                                                       CP2K_BLIS_LINK_LIBRARIES)
endif()

# add target to link against
if(CP2K_BLIS_FOUND AND NOT TARGET CP2K_BLIS::blis)
  add_library(CP2K_BLIS::blis INTERFACE IMPORTED)
endif()

set_property(TARGET CP2K_BLIS::blis PROPERTY INTERFACE_LINK_LIBRARIES
                                             ${CP2K_BLIS_LINK_LIBRARIES})

if(BLIS_INCLUDE_DIRS)
  set_property(TARGET CP2K_BLIS::blis PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                               ${CP2K_BLIS_INCLUDE_DIRS})
endif()

# prevent clutter in cache
mark_as_advanced(CP2K_BLIS_FOUND CP2K_BLIS_LINK_LIBRARIES
                 CP2K_BLIS_INCLUDE_DIRS)
