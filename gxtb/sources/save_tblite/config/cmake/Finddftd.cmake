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

# The dftd repository uses project("dftd4") in its CMakeLists.txt,
# so CMake targets are named dftd4/dftd4-lib. We fetch using the internal
# name and then create a dftd::dftd alias for tblite to consume.
set(_lib "dftd4")
set(_pkg "DFTD")
set(_url "git@github.com:grimme-lab/dftd.git")
set(_rev "HEAD")

if(NOT DEFINED "${_pkg}_FIND_METHOD")
  if(DEFINED "${PROJECT_NAME}-dependency-method")
    set("${_pkg}_FIND_METHOD" "${${PROJECT_NAME}-dependency-method}")
  else()
    set("${_pkg}_FIND_METHOD" "cmake" "pkgconf" "subproject" "fetch")
  endif()
  set("_${_pkg}_FIND_METHOD")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/tblite-utils.cmake")

tblite_find_package("${_lib}" "${${_pkg}_FIND_METHOD}" "${_url}" "${_rev}")

# Create dftd::dftd alias pointing to dftd4::dftd4
if(TARGET "dftd4::dftd4" AND NOT TARGET "dftd::dftd")
  add_library("dftd::dftd" INTERFACE IMPORTED)
  target_link_libraries("dftd::dftd" INTERFACE "dftd4::dftd4")
endif()
if(TARGET "dftd::dftd")
  set(dftd_FOUND TRUE)
endif()

if(DEFINED "_${_pkg}_FIND_METHOD")
  unset("${_pkg}_FIND_METHOD")
  unset("_${_pkg}_FIND_METHOD")
endif()
unset(_lib)
unset(_pkg)
unset(_url)
unset(_rev)
