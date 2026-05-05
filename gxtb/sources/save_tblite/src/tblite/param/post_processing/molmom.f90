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

!> @file tblite/param/post_processing/molmom.f90
!> Provides model for the calculation of molecular multipole moments for the converged singlepoint

module tblite_param_molecular_moments
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_symbols, only : to_number
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table, toml_array
   implicit none
   private

   public :: count

   character(len=*), parameter :: k_dipm = "dipole", k_qp = "quadrupole", &
      & k_key = "molecular-multipole"

   !> Parametrization record specifying the dispersion model
   type, public, extends(serde_record) :: molecular_multipole_record
      !> Compute density-based xtbml features
      logical :: moldipm = .false.
      !> Return vectorial information additional to norm of the corresponding multipole moments
      logical :: molqp = .false.


   contains
      generic :: load => load_from_array
      generic :: dump => dump_to_array
      !> Read parametrization data from TOML data structure
      procedure :: load_from_toml
      !> Write parametrization data to TOML data structure
      procedure :: dump_to_toml
      !> Read parametrization data from parameter array
      procedure, private :: load_from_array
      !> Write parametrization data to parameter array
      procedure, private :: dump_to_array
      procedure :: populate_default_param
   end type


   !> Masking for the dispersion model
   type, public :: molmom_mask
   end type molmom_mask


   interface count
      module procedure :: count_mask
   end interface count


contains

subroutine populate_default_param(param)
   class(molecular_multipole_record), intent(inout) :: param

   param%moldipm = .true.
   param%molqp= .true.


end subroutine

!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(molecular_multipole_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   type(toml_table), pointer :: child
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   integer :: stat
   call get_value(table, k_key, child, requested=.false.)
   if (.not.associated(child)) return
   call get_value(child, k_dipm, self%moldipm, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read entry for molecular dipole, boolean expected")
      return
   end if

   call get_value(child, k_qp, self%molqp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read entry for molecular quadrupole, boolean expected")
      return
   end if

end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(molecular_multipole_record), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   call add_table(table, k_key, child)

   call set_value(child, k_dipm, self%moldipm)
   call set_value(child, k_qp, self%molqp)

end subroutine dump_to_toml


!> Read parametrization data from parameter array
subroutine load_from_array(self, array, offset, base, mask, error)
   class(molecular_multipole_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: offset
   type(molecular_multipole_record), intent(in) :: base
   type(molmom_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   select type(self)
   type is (molecular_multipole_record)
      self = base
   end select

end subroutine load_from_array

!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, offset, mask, error)
   class(molecular_multipole_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: offset
   type(molmom_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

end subroutine dump_to_array

elemental function count_mask(mask) result(ncount)
   type(molmom_mask), intent(in) :: mask
   integer :: ncount
   ncount = 0
end function count_mask


end module tblite_param_molecular_moments
