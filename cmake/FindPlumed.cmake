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
include(cp2k_utils)
find_package(PkgConfig)

cp2k_set_default_paths(PLUMED "Plumed")

# First try with pkg
if(PKG_CONFIG_FOUND)
  # plumed has a pkg-config module
  pkg_check_module(CP2K_PLUMED IMPORTED_TARGET GLOBAL "plumed")
endif()

if(NOT ${CP2K_PLUMED_FOUND})
  cp2k_find_libraries(PLUMED "plumed")
endif()

cp2k_include_dirs(PLUMED "plumed.h plumed/plumed.h")

if(CP2K_PLUMED_INCLUDE_DIRS)
  find_package_handle_standard_args(Plumed DEFAULT_MSG CP2K_PLUMED_INCLUDE_DIRS
                                    CP2K_PLUMED_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(Plumed DEFAULT_MSG
                                    CP2K_PLUMED_LINK_LIBRARIES)
endif()

if(CP2K_PLUMED_FOUND AND NOT TARGET CP2K_plumed::plumed)
  add_library(CP2K_plumed::plumed INTERFACE IMPORTED)
  set_target_properties(
    CP2K_plumed::plumed
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_PLUMED_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${CP2K_PLUMED_LINK_LIBRARIES}")
endif()

mark_as_advanced(CP2K_PLUMED_LINK_LIBRARIES CP2K_PLUMED_INCLUDE_DIRS
                 CP2K_PLUMED_FOUND)
