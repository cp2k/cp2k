! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

!> Implementation of a parameter database for damping parameters
module dftd3_app_toml
   use mctc_env, only : error_type, fatal_error
   use dftd3, only : d3_param, damping_param, rational_damping_param, new_rational_damping, &
      & zero_damping_param, new_zero_damping, mzero_damping_param, new_mzero_damping, &
      & optimizedpower_damping_param, new_optimizedpower_damping, &
      & cso_damping_param, new_cso_damping
   use tomlf, only : toml_table, toml_array, toml_key, toml_error, toml_parse, &
      & get_value, len
   implicit none
   private

   public :: param_database

   !> Individual damping parameter record
   type :: param_record
      !> Functional name identifying this record
      character(len=:), allocatable :: key
      !> Damping function identifier for this record
      character(len=:), allocatable :: id
      !> Actual damping parameters
      type(d3_param) :: param
      !> Name of the damping function
      character(len=:), allocatable :: damping
      !> Reference to publication
      character(len=:), allocatable :: doi
   end type param_record

   !> Damping parameter database
   type :: param_database
      !> List of supported damping functions
      type(param_record), allocatable :: defaults(:)
      !> List of damping parameter records
      type(param_record), allocatable :: records(:)
      !> Mask for default damping functions in queries
      logical, allocatable :: mask(:)
   contains
      !> Reading of damping parameter data
      generic :: load => load_from_file, load_from_unit, load_from_toml
      !> Read damping parameter data from file
      procedure, private :: load_from_file
      !> Read damping parameter data from formatted unit
      procedure, private :: load_from_unit
      !> Read damping parameter data from TOML data structure
      procedure, private :: load_from_toml
      !> Get parameters
      procedure :: get
   end type param_database

contains

!> Read damping parameter data from file
subroutine load_from_file(self, file, error)
   !> Instance of the damping parameter data
   class(param_database), intent(inout) :: self
   !> File name
   character(len=*), intent(in) :: file
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: unit
   logical :: exist

   inquire(file=file, exist=exist)
   if (.not.exist) then
     call fatal_error(error, "Could not find parameter file '"//file//"'")
     return
   end if

   open(file=file, newunit=unit)
   call self%load(unit, error)
   close(unit)
end subroutine load_from_file


!> Read damping parameter data from file
subroutine load_from_unit(self, unit, error)
   !> Instance of the damping parameter data
   class(param_database), intent(inout) :: self
   !> File name
   integer, intent(in) :: unit
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_error), allocatable :: parse_error
   type(toml_table), allocatable :: table

   call toml_parse(table, unit, parse_error)

   if (allocated(parse_error)) then
      allocate(error)
      call move_alloc(parse_error%message, error%message)
      return
   end if

   call self%load(table, error)
   if (allocated(error)) return

end subroutine load_from_unit

!> Read damping parameter data from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the damping parameter data
   class(param_database), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   call get_value(table, "default", child)
   call load_default(self, child, error)
   if (allocated(error)) return

   call get_value(table, "parameter", child)
   call load_parameter(self, child, error)
   if (allocated(error)) return
end subroutine load_from_toml

!> Read the defaults from the TOML table
subroutine load_default(self, table, error)
   !> Instance of the damping parameter data
   type(param_database), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child, child2
   type(toml_array), pointer :: children
   type(toml_key), allocatable :: keys(:)
   integer :: ik
   type(param_record) :: stub
   character(len=:), allocatable :: val

   call get_value(table, "d3", children)
   call get_value(table, "parameter", child)
   call get_value(child, "d3", child2)

   call child2%get_keys(keys)
   call resize(self%defaults, size(keys))
   do ik = 1, size(keys)
      call get_value(child2, keys(ik)%key, child)
      call load_record(self%defaults(ik), child, stub, error)
      self%defaults(ik)%key = ""
      if (allocated(error)) exit
   end do

   allocate(self%mask(size(keys)), source=.false.)
   do ik = 1, len(children)
      call get_value(children, ik, val)
      associate(id => get_record(self%defaults, "", val))
         if (id > 0) self%mask(id) = .true.
      end associate
   end do
end subroutine load_default

!> Deserialize the parameter subtable into a list of records
subroutine load_parameter(self, table, error)
   !> Instance of the damping parameter data
   type(param_database), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(param_record) :: stub
   type(toml_key), allocatable :: keys(:), list(:)
   type(toml_table), pointer :: child, child2
   integer :: nr, ik, iv, id

   nr = 0
   call table%get_keys(keys)
   call resize(self%records, size(keys))

   records: do ik = 1, size(keys)
      call get_value(table, keys(ik)%key, child)
      call get_value(child, "d3", child2)

      call child2%get_keys(list)
      if (nr + size(list) > size(self%records)) call resize(self%records)
      do iv = 1, size(list)
         id = get_record(self%defaults, "", list(iv)%key)
         call get_value(child2, list(iv)%key, child)
         if (id > 0) then
            call load_record(self%records(iv+nr), child, self%defaults(id), error)
         else
            call load_record(self%records(iv+nr), child, stub, error)
         end if
         self%records(iv+nr)%key = keys(ik)%key
         if (allocated(error)) exit records
      end do
      nr = nr + size(list)

   end do records
   if (allocated(error)) return
   call resize(self%records, nr)
end subroutine load_parameter

!> Deserialize a record from a TOML table
subroutine load_record(record, table, default, error)
   !> Instance of the damping parameter data
   type(param_record), intent(inout) :: record
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Default values
   type(param_record), intent(in) :: default
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call table%get_key(record%id)
   call get_value(table, "damping", record%damping, default%damping)
   call get_value(table, "doi", record%doi, default%doi)
   call get_value(table, "s6", record%param%s6, default%param%s6)
   call get_value(table, "s8", record%param%s8, default%param%s8)
   call get_value(table, "s9", record%param%s9, default%param%s9)
   call get_value(table, "a1", record%param%a1, default%param%a1)
   call get_value(table, "a2", record%param%a2, default%param%a2)
   call get_value(table, "rs6", record%param%rs6, default%param%rs6)
   call get_value(table, "rs8", record%param%rs8, default%param%rs8)
   call get_value(table, "alp", record%param%alp, default%param%alp)
   call get_value(table, "bet", record%param%bet, default%param%bet)
end subroutine load_record

!> Load damping parameters from data base
subroutine get(self, param, method, damping)
   !> Instance of the damping parameter data
   class(param_database), intent(inout) :: self
   !> Damping parameters
   class(damping_param), allocatable, intent(out) :: param
   !> Method identifier
   character(len=*), intent(in) :: method
   !> Damping function identifier
   character(len=*), intent(in), optional :: damping

   integer :: ir, id

   if (present(damping)) then
      ir = get_record(self%records, method, damping)
   else
      ir = 0
      do id = 1, size(self%defaults)
         if (self%mask(id)) then
            ir = get_record(self%records, method, self%defaults(id)%id)
            if (ir > 0) exit
         end if
      end do
   end if
   if (ir == 0) return

   associate(record => self%records(ir))
      select case(record%damping)
      case("rational")
         block
            type(rational_damping_param), allocatable :: tmp
            allocate(tmp)
            call new_rational_damping(tmp, record%param)
            call move_alloc(tmp, param)
         end block
      case("zero")
         block
            type(zero_damping_param), allocatable :: tmp
            allocate(tmp)
            call new_zero_damping(tmp, record%param)
            call move_alloc(tmp, param)
         end block
      case("mzero")
         block
            type(mzero_damping_param), allocatable :: tmp
            allocate(tmp)
            call new_mzero_damping(tmp, record%param)
            call move_alloc(tmp, param)
         end block
      case("optimizedpower")
         block
            type(optimizedpower_damping_param), allocatable :: tmp
            allocate(tmp)
            call new_optimizedpower_damping(tmp, record%param)
            call move_alloc(tmp, param)
         end block
      case("cso")
         block
            type(cso_damping_param), allocatable :: tmp
            allocate(tmp)
            call new_cso_damping(tmp, record%param)
            call move_alloc(tmp, param)
         end block
      end select
   end associate
end subroutine get


!> Find a record in the record list
pure function get_record(record, key, id) result(pos)
   !> Instance of the parameters
   type(param_record), intent(in) :: record(:)
   !> Key to find
   character(len=*), intent(in) :: key
   !> Identifier to find
   character(len=*), intent(in) :: id
   !> Position in records
   integer :: pos

   integer :: ii

   pos = 0
   do ii = 1, size(record)
      if (record(ii)%key == key .and. record(ii)%id == id) then
         pos = ii
         exit
      end if
   end do
end function get_record


!> Reallocate list of records
pure subroutine resize(var, n)
   !> Instance of the array to be resized
   type(param_record), allocatable, intent(inout) :: var(:)
   !> Dimension of the final array size
   integer, intent(in), optional :: n

   type(param_record), allocatable :: tmp(:)
   integer :: this_size, new_size
   integer, parameter :: initial_size = 16

   if (allocated(var)) then
      this_size = size(var, 1)
      call move_alloc(var, tmp)
   else
      this_size = initial_size
   end if

   if (present(n)) then
      new_size = n
   else
      new_size = this_size + this_size/2 + 1
   end if

   allocate(var(new_size))

   if (allocated(tmp)) then
      this_size = min(size(tmp, 1), size(var, 1))
      var(:this_size) = tmp(:this_size)
      deallocate(tmp)
   end if
end subroutine resize

end module dftd3_app_toml
