# Copyright (c) 2019 ETH Zurich, Simon Frasch
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

cp2k_set_default_paths(SCALAPACK "SCALAPACK")

# check if we have mkl as blas library or not and pick the scalapack from mkl
# distro if found

if(CP2K_SCALAPACK_VENDOR STREQUAL "GENERIC")
  if(TARGET CP2K_MKL::scalapack_link)
    message("-----------------------------------------------------------------")
    message("-                  FindScalapack warning                        -")
    message("-----------------------------------------------------------------")
    message("                                                                 ")
    message(
      WARNING
        "You may want to use mkl implementation of scalapack. To do this add -DSCALAPACK_VENDOR=MKL to the cmake command line"
    )
  endif()

  if(TARGET CP2K_SCI::scalapack_link)
    message("-----------------------------------------------------------------")
    message("-                  FindScalapack warning                        -")
    message("-----------------------------------------------------------------")
    message("                                                                 ")
    message(
      WARNING
        "You may want to use Cray implementation of scalapack. To do this add -DSCALAPACK_VENDOR=SCI to the cmake command line"
    )
    message("                                                                 ")
    message("                                                                 ")
  endif()

  # try to detect location with pkgconfig
  find_package(PkgConfig QUIET)
  if(PKG_CONFIG_FOUND)
    pkg_check_modules(CP2K_SCALAPACK "scalapack")
  endif()

  # this should be enough for detecting scalapack compiled by hand. If scalapack
  # is vendor specific then we sahould have a target blas::scalapack available.
  # it removes the problem of modifying too many files when we add a vendor
  # specific blas/lapack/scalapack implementation

  if(NOT CP2K_SCALAPACK_FOUND)
    cp2k_find_libraries(SCALAPACK "scalapack")
  endif()
elseif(TARGET CP2K_MKL::scalapack_link)
  # we have mkl check for the different mkl target
  get_target_property(CP2K_SCALAPACK_LINK_LIBRARIES CP2K_MKL::scalapack_link
                      INTERFACE_LINK_LIBRARIES)
  set(CP2K_SCALAPACK_FOUND yes)
elseif(TARGET CP2K_SCI::scalapack_link)
  # we have mkl check for the different mkl target
  get_target_property(CP2K_SCALAPACK_LINK_LIBRARIES CP2K_SCI::scalapack_link
                      INTERFACE_LINK_LIBRARIES)
  set(CP2K_SCALAPACK_FOUND yes)
endif()

# check if found
find_package_handle_standard_args(SCALAPACK
                                  REQUIRED_VARS CP2K_SCALAPACK_LINK_LIBRARIES)
# prevent clutter in cache

# add target to link against
if(CP2K_SCALAPACK_FOUND)

  if(NOT TARGET CP2K_SCALAPACK::scalapack)
    add_library(CP2K_SCALAPACK::scalapack INTERFACE IMPORTED)
  endif()

  set_property(
    TARGET CP2K_SCALAPACK::scalapack PROPERTY INTERFACE_LINK_LIBRARIES
                                              ${CP2K_SCALAPACK_LINK_LIBRARIES})
endif()
mark_as_advanced(CP2K_SCALAPACK_FOUND CP2K_SCALAPACK_LINK_LIBRARIES)
