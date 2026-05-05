# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.

# Handling of subproject dependencies
macro(
  "tblite_find_package"
  package
  methods
  url
  revision
)
  string(TOLOWER "${package}" _pkg_lc)
  string(TOUPPER "${package}" _pkg_uc)

  foreach(method in ITEMS ${methods})

    if(TARGET "${package}::${package}")
      break()
    endif()

    if("${method}" STREQUAL "cmake")
      message(STATUS "${package}: Find installed package")
      find_package("${package}" CONFIG QUIET)
      if("${package}_FOUND")
        message(STATUS "${package}: Found installed package")
        break()
      endif()
    endif()

    if("${method}" STREQUAL "pkgconf")
      find_package(PkgConfig QUIET)
      pkg_check_modules("${_pkg_uc}" QUIET "${package}")
      if("${_pkg_uc}_FOUND")
        message(STATUS "Found ${package} via pkg-config")

        add_library("${package}::${package}" INTERFACE IMPORTED)
        target_link_libraries(
          "${package}::${package}"
          INTERFACE
          "${${_pkg_uc}_LINK_LIBRARIES}"
        )
        target_include_directories(
          "${package}::${package}"
          INTERFACE
          "${${_pkg_uc}_INCLUDE_DIRS}"
        )
        break()
      endif()
    endif()

    if("${method}" STREQUAL "subproject")
      set("${_pkg_uc}_SOURCE_DIR" "${PROJECT_SOURCE_DIR}/subprojects/${package}")
      set("${_pkg_uc}_BINARY_DIR" "${PROJECT_BINARY_DIR}/subprojects/${package}")
      if(EXISTS "${${_pkg_uc}_SOURCE_DIR}/CMakeLists.txt")
        message(STATUS "Include ${package} from subprojects")
        add_subdirectory(
          "${${_pkg_uc}_SOURCE_DIR}"
          "${${_pkg_uc}_BINARY_DIR}"
        )

        add_library("${package}::${package}" INTERFACE IMPORTED)
        target_link_libraries("${package}::${package}" INTERFACE "${package}")

        # We need the module directory in the subproject before we finish the configure stage
        if(NOT EXISTS "${${_pkg_uc}_BINARY_DIR}/include")
          file(MAKE_DIRECTORY "${${_pkg_uc}_BINARY_DIR}/include")
        endif()

        break()
      endif()
    endif()

    if("${method}" STREQUAL "fetch")
      message(STATUS "Retrieving ${package} from ${url}")
      include(FetchContent)
      FetchContent_Declare(
        "${_pkg_lc}"
        GIT_REPOSITORY "${url}"
        GIT_TAG "${revision}"
      )
      FetchContent_MakeAvailable("${_pkg_lc}")

      add_library("${package}::${package}" INTERFACE IMPORTED)
      target_link_libraries("${package}::${package}" INTERFACE "${package}")

      # We need the module directory in the subproject before we finish the configure stage
      FetchContent_GetProperties("${_pkg_lc}" BINARY_DIR "${_pkg_uc}_BINARY_DIR")
      if(NOT EXISTS "${${_pkg_uc}_BINARY_DIR}/include")
        file(MAKE_DIRECTORY "${${_pkg_uc}_BINARY_DIR}/include")
      endif()

      break()
    endif()

  endforeach()

  unset(_pkg_lc)
  unset(_pkg_uc)

  if(NOT TARGET "${package}::${package}")
    message(FATAL_ERROR "Could not find dependency ${package}")
  endif()

endmacro()
