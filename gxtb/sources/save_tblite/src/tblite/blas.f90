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

!> @dir tblite/blas
!> Contains wrappers for the basic linear algebra subprograms

!> @file tblite/blas.f90
!> Provides high-level interfaces to the basic linear algebra subprograms

!> Proxy module to reexport high-level basic linear algebra subprogram wrappers
module tblite_blas
   use tblite_blas_level1, only : dot => wrap_dot, axpy => wrap_axpy
   use tblite_blas_level2, only : gemv => wrap_gemv, symv => wrap_symv
   use tblite_blas_level3, only : gemm => wrap_gemm, symm => wrap_symm
   implicit none
   public
end module tblite_blas
