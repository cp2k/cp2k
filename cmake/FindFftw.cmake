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

# * Find the FFTW library
#
# Usage: find_package(FFTW [REQUIRED] [QUIET] )

include(FindPackageHandleStandardArgs)
include(cp2k_utils)

cp2k_set_default_paths(FFTW3 "Fftw")

# Check if we can use PkgConfig
find_package(PkgConfig)

# First try with pkg
if(PKG_CONFIG_FOUND)
  pkg_search_module(CP2K_FFTW3 fftw3)
  pkg_search_module(CP2K_FFTW3F fftw3f)
  pkg_search_module(CP2K_FFTW3L fftw3l)
  pkg_search_module(CP2K_FFTW3Q fftw3q)
endif()

foreach(_lib fftw3 fftw3f fftw3l fftw3q)
  if(NOT CP2K_${__lib_up}_FOUND)
    cp2k_find_libraries("${__lib_up}" ${_lib})
  endif()

  # OMP variant
  foreach(_subtype "mpi" "omp" "threads")
    string(TOUPPER "${_lib}_${_subtype}" _sub_lib)
    cp2k_find_libraries("${_sub_lib}" ${_lib}_${_subtype})
  endforeach()
endforeach()

if(NOT CP2K_FFTW3_INCLUDE_DIRS)
  cp2k_include_dirs(FFTW3 "fftw3.h;fftw3/fftw3.h")
endif()

if(CP2K_FFTW3_INCLUDE_DIRS)
  find_package_handle_standard_args(Fftw DEFAULT_MSG CP2K_FFTW3_INCLUDE_DIRS
                                    CP2K_FFTW3_LINK_LIBRARIES)
else()
  find_package_handle_standard_args(Fftw DEFAULT_MSG CP2K_FFTW3_LINK_LIBRARIES)
endif()

foreach(lib_name "fftw3" "fftw3l" "fftw3q" "fftw3f")
  string(TOUPPER "${lib_name}" __lib_name_up)

  if(CP2K_${__lib_name_up}_FOUND AND NOT TARGET CP2K_FFTW3::${lib_name})
    add_library(CP2K_FFTW3::${lib_name} INTERFACE IMPORTED)
    # we do not recheck if the libraries are found when pkg_config is
    # successful.
    set_target_properties(
      CP2K_FFTW3::${lib_name}
      PROPERTIES INTERFACE_LINK_LIBRARIES
                 "${CP2K_${__lib_name_up}_LINK_LIBRARIES}")

    if(CP2K_FFTW3_INCLUDE_DIRS)
      set_target_properties(
        CP2K_FFTW3::${lib_name} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                           "${CP2K_FFTW3_INCLUDE_DIRS}")
    endif()

    foreach(sub_type "threads" "mpi" "omp")
      string(TOUPPER "${lib_name}_${sub_type}" __libs)
      if(CP2K_${__libs}_FOUND AND NOT TARGET
                                  CP2K_FFTW3::${lib_name}_${sub_type})
        add_library(CP2K_FFTW3::${lib_name}_${sub_type} INTERFACE IMPORTED)
        set_target_properties(
          CP2K_FFTW3::${lib_name}_${sub_type}
          PROPERTIES INTERFACE_LINK_LIBRARIES
                     "${CP2K_${__libs}_LINK_LIBRARIES}")
      endif()
    endforeach()
  endif()
endforeach()

set(CP2K_FFTW3_FOUND ON)
mark_as_advanced(
  CP2K_FFTW3_FOUND
  CP2K_FFTW3_PREFIX
  CP2K_FFTW3_INCLUDE_DIRS
  CP2K_FFTW3_MPI
  FFTW3_OMP
  CP2K_FFTW3_THREADS
  CP2K_FFTW3Q_OMP
  CP2K_FFTW3Q_THREADS
  CP2K_FFTW3F_MPI
  CP2K_FFTW3_OMP
  CP2K_FFTW3F_THREADS
  CP2K_FFTW3L_MPI
  CP2K_FFTW3L_OMP
  CP2K_FFTW3L_THREADS)
