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
find_package(PkgConfig REQUIRED)

cp2k_set_default_paths(LIBXSMM)
set(CP2K_LIBXSMMEXT_PREFIX "${CP2K_LIBXSMM_PREFIX}")
set(CP2K_LIBXSMMF_PREFIX "${CP2K_LIBXSMM_PREFIX}")
set(CP2K_LIBXSMMNOBLAS_PREFIX "${CP2K_LIBXSMM_PREFIX}")

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CP2K_LIBXSMM IMPORTED_TARGET GLOBAL libxsmm)
  pkg_check_modules(CP2K_LIBXSMMEXT IMPORTED_TARGET GLOBAL libxsmmext)
  pkg_check_modules(CP2K_LIBXSMMF IMPORTED_TARGET GLOBAL libxsmmf)
  pkg_check_modules(CP2K_LIBXSMMNOBLAS IMPORTED_TARGET GLOBAL libxsmmnoblas)
endif()

foreach(__lib libxsmm libxsmmf libxsmmext libxsmmnoblas)
  string(TOUPPER "${__lib}" __lib_search_up)
  if(NOT CP2K_${__lib_search_up}_FOUND)
    cp2k_find_libraries(${__lib_search_up} ${__lib})
  endif()
endforeach()

if(CP2K_LIBXSMM_INCLUDE_DIRS)
  find_package_handle_standard_args(
    LibXSMM
    DEFAULT_MSG
    CP2K_LIBXSMM_INCLUDE_DIRS
    CP2K_LIBXSMMNOBLAS_LINK_LIBRARIES
    CP2K_LIBXSMMEXT_LINK_LIBRARIES
    CP2K_LIBXSMMF_LINK_LIBRARIES
    CP2K_LIBXSMM_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(
    LibXSMM DEFAULT_MSG CP2K_LIBXSMMNOBLAS_LINK_LIBRARIES
    CP2K_LIBXSMMEXT_LINK_LIBRARIES CP2K_LIBXSMMF_LINK_LIBRARIES
    CP2K_LIBXSMM_LINK_LIBRARIES)
  set(CP2K_LIBXSMM_INCLUDE_DIRS "/usr/include;/usr/include/libxsmm")
endif()

if(NOT TARGET CP2K_LibXSMM::libxsmm)
  foreach(__lib libxsmm libxsmmf libxsmmext libxsmmnoblas)
    string(TOUPPER "CP2K_${__lib}" __lib_search_up)
    if(${__lib_search_up}_FOUND AND NOT TARGET CP2K_LibXSMM::${__lib})
      add_library(CP2K_LibXSMM::${__lib} INTERFACE IMPORTED)
      set_target_properties(
        CP2K_LibXSMM::${__lib}
        PROPERTIES INTERFACE_LINK_LIBRARIES
                   "${${__lib_search_up}_LINK_LIBRARIES}")
    endif()
  endforeach()
  message("LibXSMM : ${CP2K_LIBXSMM_INCLUDE_DIRS}")
  if(CP2K_LIBXSMM_INCLUDE_DIRS)
    set_target_properties(
      CP2K_LibXSMM::libxsmm PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                       "${CP2K_LIBXSMM_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(
  CP2K_LIBXSMM_INCLUDE_DIRS CP2K_LIBXSMMNOBLAS_LINK_LIBRARIES
  CP2K_LIBXSMMEXT_LINK_LIBRARIES CP2K_LIBXSMMF_LINK_LIBRARIES
  CP2K_LIBXSMM_LINK_LIBRARIES)
