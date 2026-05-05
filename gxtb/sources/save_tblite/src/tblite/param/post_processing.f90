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

!> @file tblite/param/post_processing.f90
!> Provides a list-type container to store various post processing param records

!> Definition of the dispersion corrections
module tblite_param_post_processing
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_symbols, only : to_number
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, get_value, set_value, add_table, toml_array, toml_key
   use tblite_param_molecular_moments, only : molecular_multipole_record
   use tblite_param_xtbml_features, only : xtbml_features_record
   implicit none
   private

   public :: molecular_multipole_record, post_processing_param_list
   !> Parametrization record specifying the dispersion model
   type :: post_processing_record
      class(serde_record), allocatable :: record
   end type

   type :: post_processing_param_list
      type(post_processing_record), allocatable :: list(:)
      integer :: n = 0
   contains
      private
      generic, public :: load => load_from_toml
      generic, public :: dump => dump_to_toml
      procedure, public :: push
      procedure, public :: get_n_records
      procedure :: load_from_toml
      procedure :: dump_to_toml
   end type

contains

function get_n_records(self) result(n)
   class(post_processing_param_list), intent(in) :: self
   integer n
   n = self%n
end function

subroutine push(self, record)
   class(post_processing_param_list), intent(inout) :: self
   class(serde_record), allocatable, intent(inout) :: record

   if (.not.allocated(self%list)) call resize(self%list)
   if (self%n >= size(self%list)) then
      call resize(self%list)
   end if

   self%n = self%n + 1
   call move_alloc(record, self%list(self%n)%record)
end subroutine push

subroutine resize(list, n)
   !> Instance of the array to be resized
   type(post_processing_record), allocatable, intent(inout) :: list(:)
   !> Dimension of the final array size
   integer, intent(in), optional :: n

   type(post_processing_record), allocatable :: tmp(:)
   integer :: this_size, new_size, item
   integer, parameter :: initial_size = 20

   if (allocated(list)) then
      this_size = size(list, 1)
      call move_alloc(list, tmp)
   else
      this_size = initial_size
   end if

   if (present(n)) then
      new_size = n
   else
      new_size = this_size + this_size/2 + 1
   end if

   allocate(list(new_size))

   if (allocated(tmp)) then
      this_size = min(size(tmp, 1), size(list, 1))
      do item = 1, this_size
         call move_alloc(tmp(item)%record, list(item)%record)
      end do
      deallocate(tmp)
   end if

end subroutine resize

subroutine dump_to_toml(self, table, error)
   !> List of all element records
   class(post_processing_param_list), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii

   do ii = 1, self%get_n_records()
      associate(rec => self%list(ii)%record)
         call rec%dump(table, error)
      end associate
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine dump_to_toml

   !> Deserialize records from a table by iterating over all entries
subroutine load_from_toml(self, table, error)
   !> List of all element records
   class(post_processing_param_list), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii
   type(toml_key), allocatable :: list(:)
   character(len=:), allocatable :: key
   call table%get_keys(list)
   do ii = 1, size(list)
      key = trim(adjustl(list(ii)%key))
      select case(key)
      case("molecular-multipole")
         block
            type(molecular_multipole_record), allocatable :: tmp_record
            class(serde_record), allocatable :: cont
            allocate(tmp_record)
            call tmp_record%load(table, error) 
            call move_alloc(tmp_record, cont)
            call self%push(cont)
         end block
      case("xtbml")
         block
            type(xtbml_features_record), allocatable :: tmp_record
            class(serde_record), allocatable :: cont
            allocate(tmp_record)
            call tmp_record%load(table, error) 
            call move_alloc(tmp_record, cont)
            call self%push(cont)
         end block
      case default
         return
      end select
   end do
   
   if (allocated(error)) return
end subroutine load_from_toml

end module tblite_param_post_processing
