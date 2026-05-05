# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

if(NOT BLAS_FOUND)
  if("${BLA_VENDOR}" MATCHES "^Intel" OR DEFINED ENV{MKLROOT})
    # C must be enabled to use MKL
    # https://cmake.org/cmake/help/v3.14/module/FindBLAS.html#result-variables
    enable_language("C")
  endif()
  find_package("BLAS" REQUIRED)
  if(NOT TARGET "BLAS::BLAS")
    add_library("BLAS::BLAS" INTERFACE IMPORTED)
    target_link_libraries("BLAS::BLAS" INTERFACE "${BLAS_LIBRARIES}")
  endif()
endif()
