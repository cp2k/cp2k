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

cmake_minimum_required(VERSION 3.14)

project(
  "dftd4"
  LANGUAGES "Fortran"
  VERSION "4.0.2"
  DESCRIPTION "Generally Applicable Atomic-Charge Dependent London Dispersion Correction"
)

# Follow GNU conventions for installing directories
include(GNUInstallDirs)

# General configuration information
add_subdirectory("config")

if(NOT TARGET "OpenMP::OpenMP_Fortran" AND WITH_OpenMP)
  find_package("OpenMP" REQUIRED)
endif()

if(WITH_ILP64 AND BLAS_LIBRARIES)
  message(STATUS "Using LAPACK/BLAS ILP64 interface")
elseif(WITH_ILP64)
  message(FATAL_ERROR "ILP64 support needs BLAS_LIBRARIES")
endif()

if(NOT TARGET "BLAS::BLAS")
  find_package("custom-blas" REQUIRED)
endif()

# Collect subprojects
if(NOT TARGET "mctc-lib::mctc-lib")
  find_package("mctc-lib")
endif()
if(NOT TARGET "mstore::mstore")
  find_package("mstore")
endif()
if(NOT TARGET "multicharge::multicharge")
  find_package("multicharge")
endif()
set(
  lib-deps
  "mctc-lib::mctc-lib"
  "multicharge::multicharge"
)

# Collect source of the project
set(srcs)
add_subdirectory("src")

# DFT-D4 library target
add_library(
  "${PROJECT_NAME}-lib"
  "${srcs}"
)
set_target_properties(
  "${PROJECT_NAME}-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)
target_link_libraries(
  "${PROJECT_NAME}-lib"
  PUBLIC
  "${lib-deps}"
  "$<$<BOOL:${WITH_OpenMP}>:OpenMP::OpenMP_Fortran>"
)
target_include_directories(
  "${PROJECT_NAME}-lib"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)
target_include_directories(
  "${PROJECT_NAME}-lib"
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dftd4
)

if(WITH_ILP64)
  target_compile_definitions(
    "${PROJECT_NAME}-lib"
    PUBLIC -DIK=i8)
else()
  target_compile_definitions(
    "${PROJECT_NAME}-lib"
    PUBLIC -DIK=i4)
endif()

# Add example application
add_subdirectory("app")

# Export targets for other projects
add_library("${PROJECT_NAME}" INTERFACE)
target_link_libraries("${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}-lib")
install(
  TARGETS
  "${PROJECT_NAME}"
  "${PROJECT_NAME}-lib"
  EXPORT
  "${PROJECT_NAME}-targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)
if(WITH_API)
  enable_language("C")
  install(
    DIRECTORY
    "${CMAKE_CURRENT_SOURCE_DIR}/include/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )
endif()
# Package license files
install(
  FILES
  "COPYING"
  "COPYING.LESSER"
  DESTINATION "${CMAKE_INSTALL_DATADIR}/licenses/${PROJECT_NAME}"
)

# add the testsuite
enable_testing()
add_subdirectory("test")
