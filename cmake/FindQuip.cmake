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
find_package(PkgConfig)
cp2k_set_default_paths(LIBQUIP)

pkg_search_module(CP2K_LIBQUIP quip)

if(CP2K_LIBQUIP_FOUND)
  cp2k_find_libraries(LIBQUIP "quip")
  cp2k_include_dirs(LIBQUIP "quip.h")
endif()
find_package_handle_standard_args(Quip DEFAULT_MSG CP2K_LIBQUIP_INCLUDE_DIRS
                                  CP2K_QUIP_LINK_LIBRARIES)

if(CP2K_LIBQUIP_FOUND AND NOT TARGET CP2K_quip::quip)
  add_library(CP2K_quip::quip INTERFACE IMPORTED)
  set_target_properties(
    CP2K_quip::quip
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${CP2K_LIBQUIP_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${CP2K_LIBQUIP_LINK_LIBRARIES}")
endif()
