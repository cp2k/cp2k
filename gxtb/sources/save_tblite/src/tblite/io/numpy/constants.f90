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

!> @file src/tblite/io/numpy/constants.f90
!> Constants for npy and npz file format

!> Headers for numpy file formats and zip file format
module tblite_io_numpy_constants
   use mctc_env, only : sp, dp, i2, i4
   implicit none
   private

   public :: &
       magic_number, magic_string, &
       type_i1, type_i2, type_i4, type_i8, type_rsp, type_rdp, type_csp, type_cdp, &
       zip_global_sig, zip_local_sig, zip_footer_sig, zip_min_version

   character(len=*), parameter :: &
       type_i1 = "<i1", type_i2 = "<i2", type_i4 = "<i4", type_i8 = "<i8", &
       type_rsp = "<f4", type_rdp = "<f8", type_csp = "<c8", type_cdp = "<c16"

   character(len=*), parameter :: &
       & magic_number = char(int(z"93")), &
       & magic_string = "NUMPY"

   integer(i4), parameter :: &
       & zip_global_sig = int(z'02014b50', i4), &
       & zip_local_sig = int(z'04034b50', i4), &
       & zip_footer_sig = int(z'06054b50', i4)

   integer(i2), parameter :: zip_min_version = 20_i2

end module tblite_io_numpy_constants