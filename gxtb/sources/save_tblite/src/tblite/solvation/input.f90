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

!> @file tblite/solvation/input.f90
!> Provides a collection of all input configurations

!> Collection of the configuration types for all available implicit solvation models
module tblite_solvation_input
   use tblite_solvation_alpb, only : alpb_input
   use tblite_solvation_cpcm, only : cpcm_input
   use tblite_solvation_cds,  only : cds_input
   use tblite_solvation_shift,  only : shift_input
   implicit none
   private


   !> Collection of possible solvation models
   type, public :: solvation_input
      !> Input for CPCM solvation model
      type(cpcm_input), allocatable :: cpcm
      !> Input for ALPB solvation model
      type(alpb_input), allocatable :: alpb
      !> Input for CDS model 
      type(cds_input), allocatable :: cds
      !> Input for solvation shift 
      type(shift_input), allocatable :: shift
   end type solvation_input

end module tblite_solvation_input
