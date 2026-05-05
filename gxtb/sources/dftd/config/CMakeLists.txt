# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

option(BUILD_SHARED_LIBS "Whether the libraries built should be shared" FALSE)

option(WITH_API "Enable support for C API via iso_c_binding module" TRUE)
option(WITH_OpenMP "Enable support for shared memory parallelisation with OpenMP" TRUE)
option(WITH_ILP64 "Enable support for ILP64 BLAS/LAPACK calls" FALSE)
if(NOT DEFINED "${PROJECT_NAME}-dependeny-method")
  set(
    "${PROJECT_NAME}-dependency-method"
    "subproject" "cmake" "pkgconf" "fetch"
  )
endif()

set(
  "${PROJECT_NAME}-module-dir"
  "${PROJECT_NAME}/${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}"
  CACHE STRING
  "Subdirectory to install generated module files to"
)
set(
  module-dir
  "${${PROJECT_NAME}-module-dir}"
)
set(module-dir "${module-dir}" PARENT_SCOPE)

# Compiler-specific configurations

# gfortran "Interface mismatch" bug (119928)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "15.0"
     AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "15.2")
    message(FATAL_ERROR
      "Unsupported GNU Fortran version ${CMAKE_Fortran_COMPILER_VERSION}. "
      "This release has a known issue (fixed in 15.2). "
      "Please use GCC >= 15.2 or < 15.0, or a different Fortran compiler.")
  endif()
endif()

if(
  "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI"
  OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NVHPC"
  OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Flang"
)
  set(
    CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -Mbackslash -Mallocatable=03"
    PARENT_SCOPE
  )
endif()

# Set build type as CMake does not provide defaults
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(
    CMAKE_BUILD_TYPE "RelWithDebInfo"
    CACHE STRING "Build type to be used."
    FORCE
  )
  message(
    STATUS
    "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified."
  )
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" PARENT_SCOPE)
install(
  DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/template.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMinorVersion
)
install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/template.pc"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
  @ONLY
)
install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
)
