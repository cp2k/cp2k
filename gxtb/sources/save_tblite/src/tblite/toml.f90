! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/toml.f90
!> Provides reexports of the [TOML Fortran](https://toml-f.readthedocs.io) modules

!> Proxy module for TOML library implementation
module tblite_toml
   use tomlf, only : toml_table, toml_array, toml_error, toml_dump, toml_parse, &
      & toml_key, get_value, set_value, add_table, add_array, len
   use tomlf_build, only : merge_table
   use tomlf_type, only : toml_value
   implicit none
   private

   public :: toml_table, toml_array, toml_error, toml_dump, toml_parse, toml_key, &
      & toml_value, len
   public :: get_value, set_value, add_table, add_array, merge_table
end module tblite_toml
