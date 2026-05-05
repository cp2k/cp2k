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

!> @file tblite/param/post_processing/xtbml.f90
!> Provides paramerts to address the xtbml features

!> Definition of the dispersion corrections
module tblite_param_xtbml_features
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_symbols, only : to_number
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table, toml_array, add_array
   implicit none
   private

   public :: count

   character(len=*), parameter :: k_xtbmlgeometry = "geometry", k_xtbmldensity = "density", &
   & k_xtbmlorbital = "orbital", k_xtbmlenergy = "energy", k_xtbmlconvolution = "convolution", k_xtbmla = "a", &
   & k_tensor = "tensorial-output", k_xtbml = "xtbml"

   !> Parametrization record specifying the dispersion model
   type, public, extends(serde_record) :: xtbml_features_record
      !> Compute geometry-based xtbml features
      logical :: xtbml_geometry = .false.
      !> Compute density-based xtbml features
      logical :: xtbml_density = .false.
      !> Return vectorial information additional to norm of the corresponding multipole moments
      logical :: xtbml_tensor = .false.
      !> Compute orbital energy based xtbml features
      logical :: xtbml_orbital_energy = .false.
      !> Compute energy based features, necessary for partitioning weights
      logical :: xtbml_energy = .false.
      !> Compute extended feature i.e. do CN weigthed real space convolution
      logical :: xtbml_convolution = .false.
      !> Scaling for logistic function, convolution over an array of values is supported
      real(wp), allocatable :: xtbml_a(:)

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
   type, public :: ml_features_mask
   end type ml_features_mask


   interface count
      module procedure :: count_mask
   end interface count


contains

   subroutine populate_default_param(param, tensor)
      !> Instance of the parametrization data
      class(xtbml_features_record), intent(inout) :: param
      !> Return vectorial information additional and norm of the corresponding multipole moments
      logical, optional :: tensor
      ! Compute geometry-based xtbml features
      param%xtbml_geometry = .true.
      ! Compute density-based xtbml features
      param%xtbml_density = .true.
      ! Return vectorial information additional to norm of the corresponding multipole moments
      if (present(tensor)) then
         param%xtbml_tensor = tensor
      else
         param%xtbml_tensor = .false.
      end if
      ! Compute orbital energy based xtbml features
      param%xtbml_orbital_energy = .true.
      ! Compute energy based features, necessary for partitioning weights
      param%xtbml_energy = .true.
      ! Compute extended feature i.e. do CN weigthed real space convolution
      param%xtbml_convolution = .true.
      ! Scaling for logistic function, convolution over an array of values is supported
      param%xtbml_a = [1.0_wp]

   end subroutine

!> Read parametrization data from TOML data structure
   subroutine load_from_toml(self, table, error)
      use tblite_toml, only : len
      !> Instance of the parametrization data
      class(xtbml_features_record), intent(inout) :: self
      !> Data structure
      type(toml_table), intent(inout) :: table
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(toml_array), pointer :: array
      type(toml_table), pointer :: child
      integer :: stat, i, stat2
      character(len=:), allocatable :: method

      call get_value(table, k_xtbml, child, requested=.false., stat=stat)
      if (stat /= 0 .or. .not.associated(child)) then
         call fatal_error(error, "Cannot read entry for xtbml table")
         return
      end if

      call get_value(child, k_xtbmlgeometry, self%xtbml_geometry, .false., stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot read entry for xtbml geometry based features, boolean expected")
         return
      end if

      call get_value(child, k_xtbmldensity, self%xtbml_density, .false., stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot read entry for xtbml density based features, boolean expected")
         return
      end if
      if (self%xtbml_density) then
         call get_value(child, k_tensor, self%xtbml_tensor, .false., stat=stat)
         if (stat /= 0) then
            call fatal_error(error, "Cannot read entry for xtbml tensorial-output, boolean expected")
            return
         end if
      end if

      call get_value(child, k_xtbmlorbital, self%xtbml_orbital_energy, .false., stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot read entry for xtbml orbital energy based features, boolean expected")
         return
      end if

      call get_value(child, k_xtbmlenergy, self%xtbml_energy, .false., stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot read entry for xtbml geometry based features, boolean expected")
         return
      end if

      call get_value(child, k_xtbmlconvolution, self%xtbml_convolution, .false., stat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Cannot read entry for xtbml convolution, boolean expected")
         return
      end if
      if (self%xtbml_convolution) then
         if (child%has_key(k_xtbmla)) then
            call get_value(child, k_xtbmla, array, stat=stat)
            if (.not. associated(array)) then
               if (allocated(self%xtbml_a)) deallocate(self%xtbml_a)
               allocate(self%xtbml_a(1))
               call get_value(child, k_xtbmla, self%xtbml_a(1), stat=stat)
               if (stat /= 0) then
                  call fatal_error(error, "Cannot read entry for xtbml a, float expected")
                  return
               end if
            else
               if (allocated(self%xtbml_a)) deallocate(self%xtbml_a)
               allocate(self%xtbml_a(len(array)))
               do i = 1, size(self%xtbml_a)
                  call get_value(array, i, self%xtbml_a(i))
               end do
            end if
         else
            self%xtbml_a = [1.0_wp]
         end if

      end if

   end subroutine load_from_toml


!> Write parametrization data to TOML datastructure
   subroutine dump_to_toml(self, table, error)
      !> Instance of the parametrization data
      class(xtbml_features_record), intent(in) :: self
      !> Data structure
      type(toml_table), intent(inout) :: table
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(toml_table), pointer :: child

      type(toml_array), pointer :: array

      integer :: i

      call add_table(table, k_xtbml, child)

      call set_value(child, k_xtbmlgeometry, self%xtbml_geometry)
      call set_value(child, k_xtbmldensity, self%xtbml_density)
      call set_value(child, k_tensor, self%xtbml_tensor)
      call set_value(child, k_xtbmlorbital, self%xtbml_orbital_energy)
      call set_value(child, k_xtbmlenergy, self%xtbml_energy)
      call set_value(child, k_xtbmlconvolution, self%xtbml_convolution)
      call add_array(child, k_xtbmla, array)
      do i = 1, size(self%xtbml_a)
         call set_value(array, i, self%xtbml_a(i))
      end do


   end subroutine dump_to_toml


!> Read parametrization data from parameter array
   subroutine load_from_array(self, array, offset, base, mask, error)
      class(xtbml_features_record), intent(inout) :: self
      real(wp), intent(in) :: array(:)
      integer, intent(inout) :: offset
      type(xtbml_features_record), intent(in) :: base
      type(ml_features_mask), intent(in) :: mask
      type(error_type), allocatable, intent(out) :: error
      !same content as for other param files
      select type(self)
       type is (xtbml_features_record)
         self = base
      end select

   end subroutine load_from_array

!> Write parametrization data to parameter array
   subroutine dump_to_array(self, array, offset, mask, error)
      class(xtbml_features_record), intent(in) :: self
      real(wp), intent(inout) :: array(:)
      integer, intent(inout) :: offset
      type(ml_features_mask), intent(in) :: mask
      type(error_type), allocatable, intent(out) :: error
      !same content as for other param files
   end subroutine dump_to_array

   elemental function count_mask(mask) result(ncount)
      type(ml_features_mask), intent(in) :: mask
      integer :: ncount
      ncount = 0
      !same content as for other param files
   end function count_mask


end module tblite_param_xtbml_features
