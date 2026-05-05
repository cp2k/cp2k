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

!> @dir tblite/io/numpy
!> Contains npy input and output routines

!> @file tblite/io/numpy.f90
!> Provides npy input and output routines

!> Implementation of npy input and output routines
module tblite_io_numpy
   use tblite_io_numpy_load, only : load_npy
   use tblite_io_numpy_loadz, only : load_npz
   use tblite_io_numpy_save, only : save_npy
   use tblite_io_numpy_savez, only : save_npz
   implicit none
   public
end module tblite_io_numpy