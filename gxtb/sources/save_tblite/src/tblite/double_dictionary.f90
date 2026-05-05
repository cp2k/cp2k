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

!> @dir tblite/double_dictionary
!> Contains implemenation of a dictionary with strings as keys and arrays of doubles as entries.

!> @file tblite/double_dictionery.f90
!> Implements double dictionary type
module tblite_double_dictionary
   use mctc_env_accuracy, only : wp, i8
   use mctc_env, only : error_type, fatal_error
   use tblite_toml, only : toml_array, toml_table, toml_key, add_table, set_value, toml_error
   use tblite_toml, only : toml_dump, add_array,  get_value, toml_parse
   use tblite_io_numpy, only : save_npz, load_npz
   use tblite_io_numpy_loadz, only : get_npz_descriptor
   use tblite_io_numpy_zip, only : zip_file, list_zip_file
   implicit none
   private

   public :: double_dictionary_type

   type :: double_record
      character(len=:), allocatable :: label
      real(wp), allocatable :: array1(:)
      real(wp), allocatable :: array2(:, :)
      real(wp), allocatable :: array3(:, :, :)
   contains
      generic :: assignment(=) => copy_record
      generic :: operator(==) => equal_record
      procedure :: copy_record
      procedure :: equal_record
   end type double_record

   type :: double_dictionary_type
      integer :: n = 0
      type(double_record), allocatable :: record(:)
   contains
      private
      generic, public :: initialize_entry => ini_label, ini_1d, ini_2d, ini_3d
      procedure :: ini_label
      procedure :: ini_1d
      procedure :: ini_2d
      procedure :: ini_3d
      generic, public :: add_entry =>  add_1d, add_2d, add_3d
      procedure :: add_1d
      procedure :: add_2d
      procedure :: add_3d
      !procedure :: update_entry ! label and array pair
      generic, public :: get_entry =>  get_1d_index, get_2d_index, get_3d_index, get_1d_label, get_2d_label, get_3d_label !check
      procedure :: get_1d_index
      procedure :: get_2d_index
      procedure :: get_3d_index
      procedure :: get_1d_label
      procedure :: get_2d_label
      procedure :: get_3d_label
      generic, public :: update => update_1d, update_2d, update_3d
      procedure :: update_1d
      procedure :: update_2d
      procedure :: update_3d
      procedure, public :: get_label ! return label label
      procedure :: push
      procedure, public :: get_n_entries
      generic, public :: concatenate => concatenate_overwrite
      procedure :: concatenate_overwrite
      generic, public :: assignment(=) => copy
      procedure :: copy
      generic, public :: operator(+) => combine_dict
      procedure :: combine_dict
      generic, public :: remove_entry => remove_entry_label, remove_entry_index
      procedure :: remove_entry_label
      procedure :: remove_entry_index
      generic, public :: get_index => return_label_index
      procedure :: return_label_index
      procedure :: dump_to_file
      generic, public :: dump => dump_to_file
      generic, public :: operator(==) => equal_dict
      procedure :: equal_dict
      generic, public :: load => load_from_file
      procedure :: load_from_file

   end type double_dictionary_type

contains

!> Check if two double records are equal
pure function equal_record(lhs, rhs) result(equal)
   !> First record to compare
   class(double_record), intent(in) :: lhs
   !> Second record to compare
   class(double_record), intent(in) :: rhs
   !> Result of the comparison
   logical :: equal
   equal = .false.

   if (lhs%label /= rhs%label) then
      return
   end if

   if (allocated(lhs%array1) .and. allocated(rhs%array1)) then
      equal = all(lhs%array1 == rhs%array1)
      return
   end if
      
   if (allocated(lhs%array2) .and. allocated(rhs%array2)) then
      equal = all(lhs%array2 == rhs%array2)
      return
   end if

   if (allocated(lhs%array3) .and. allocated(rhs%array3)) then
      equal = all(lhs%array3 == rhs%array3)
      return
   end if

end function equal_record


!> Check if two double dictionaries are equal
pure function equal_dict(lhs, rhs) result(equal)
   !> First dictionary to compare
   class(double_dictionary_type), intent(in) :: lhs
   !> Second dictionary to compare
   class(double_dictionary_type), intent(in) :: rhs
   !> Result of the comparison
   logical :: equal

   integer :: i
   equal = .false.

   if (lhs%n /= rhs%n) then
      return
   end if

   do i = 1, lhs%n
      if (.not.(lhs%record(i) == rhs%record(i))) return
   end do

   equal = .true.
end function equal_dict

!> Read double dictionary data from file
subroutine load_from_file(self, filename, error)
   !> Instance of the parametrization data
   class(double_dictionary_type), intent(inout) :: self
   !> File name
   character(len=*), intent(in) :: filename
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: io, stat, irec, it
   integer, allocatable :: vshape(:)
   logical :: exist
   type(zip_file) :: zip
   character(len=:), allocatable :: msg, label, vtype
   character(len=*), parameter :: prefix = "tblite0_"

   inquire(file=filename, exist=exist)
   if (.not. exist) then
      call fatal_error(error, "File '"//filename//"' does not exist")
      return
   end if

   open(newunit=io, file=filename, form="unformatted", access="stream", iostat=stat)
   call list_zip_file(io, filename, zip, stat, msg)
   close(io)
   if (stat /= 0) then
      if (.not.allocated(msg)) msg = "Failed to read zip file '"//filename//"'"
      call fatal_error(error, msg)
      return
   end if

   if (.not.allocated(zip%records)) then
      call fatal_error(error, "No records found in file '"//filename//"'")
      return
   end if

   do irec = 1, size(zip%records)
      associate(path => zip%records(irec)%path)
         if (len(path) < len(prefix) + 4) then
            call fatal_error(error, "File '"//filename//"::"//path// &
               & "' does not match expected format")
            exit
         end if
         label = path(len(prefix) + 1:len(path) - 4)
         if (path(:len(prefix)) /= prefix) then
            call fatal_error(error, "Unknown file type '"//filename//"::"//path//"'")
            exit
         end if
         call get_npz_descriptor(filename, prefix//label, vtype, vshape, stat, msg)
         if (stat == 0 .and. allocated(vshape)) then
            select case(size(vshape))
            case default
               call fatal_error(error, "Unknown file type '"//filename//"::"//path//"'")
               exit
            case (1)
               call self%push(label, it)
               call load_npz(filename, prefix//label, self%record(it)%array1, stat, msg)
            case (2)
               call self%push(label, it)
               call load_npz(filename, prefix//label, self%record(it)%array2, stat, msg)
            case (3)
               call self%push(label, it)
               call load_npz(filename, prefix//label, self%record(it)%array3, stat, msg)
            end select
         end if
         if (stat /= 0) then
            if (.not.allocated(msg)) then
               msg = "Failed to load file '"//filename//"::"//path//"'"
            end if
            call fatal_error(error, msg)
            exit
         end if
      end associate
   end do

end subroutine load_from_file

!> Write double dictionary data to file
subroutine dump_to_file(self, file, error)
   !> Instance of the parametrization data
   class(double_dictionary_type), intent(in) :: self
   !> File name
   character(len=*), intent(in) :: file
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: prefix = "tblite0_"

   integer :: it, stat, io
   logical :: exist
   real(kind=wp), allocatable :: array1(:), array2(:, :), array3(:, :, :)
   character(len=:), allocatable :: msg


   inquire(file=file, exist=exist)
   if (exist) then
      open(file=file, newunit=io)
      close(io, status='delete')
   end if

   do it = 1, self%n
      associate(record => self%record(it))
         if (allocated(record%array1)) then
            call save_npz(file, prefix//record%label, record%array1, stat, msg)
            if (stat /= 0) exit
            cycle
         end if

         if (allocated(record%array2)) then
            call save_npz(file, prefix//record%label, record%array2, stat, msg)
            if (stat /= 0) exit
            cycle
         end if

         if (allocated(record%array3)) then
            call save_npz(file, prefix//record%label, record%array3, stat, msg)
            if (stat /= 0) exit
            cycle
         end if
      end associate
   end do
   if (stat /= 0) then
      if (.not.allocated(msg)) then
         msg = "Failed to write file '"//file//"'"
      end if
      call fatal_error(error, msg)
   end if
end subroutine dump_to_file


pure subroutine remove_entry_index(self, index)
   class(double_dictionary_type), intent(inout) :: self
   integer, intent(in) :: index
   integer :: old_n, i, it
   type(double_record), allocatable :: tmp(:)

   if (index > self%n) return
   call move_alloc(self%record, tmp)
   old_n = self%n
   self%n = self%n - 1
   
   allocate(self%record(self%n))
   it = 1
   do i = 1, old_n
      if (i == index) cycle
      call move_record(tmp(i), self%record(it))
      it = it + 1
   end do
 
end subroutine remove_entry_index

pure subroutine remove_entry_label(self, label)
   class(double_dictionary_type), intent(inout) :: self
   character(len=*), intent(in) :: label
   integer :: it 
   it = return_label_index(self, label)
   if (it /= 0) then
      call self%remove_entry_index(it) 
   end if
end subroutine remove_entry_label

pure function get_n_entries(self) result(n)
   class(double_dictionary_type), intent(in) :: self
   integer :: n
   n = self%n
end function get_n_entries

pure subroutine copy(to, from)
   class(double_dictionary_type), intent(inout) :: to
   type(double_dictionary_type), intent(in) :: from
   integer :: n_entries, i
   to%n = from%n
   if (allocated(to%record)) deallocate(to%record)
   if (from%n > 0) then 
      allocate(to%record(size(from%record)))
      n_entries = from%n

      do i = 1, n_entries
         to%record(i) = from%record(i)
      end do
   end if

end subroutine

pure subroutine move_record(from, to)
   class(double_record), intent(inout) :: to
   type(double_record), intent(inout) :: from

   call move_alloc(from%label, to%label)
   if (allocated(to%array1)) deallocate(to%array1)
   if (allocated(from%array1)) then
      call move_alloc(from%array1, to%array1)
   end if
   if (allocated(to%array2)) deallocate(to%array2)
   if (allocated(from%array2)) then
      call move_alloc(from%array2, to%array2)
   end if
   if (allocated(to%array3)) deallocate(to%array3)
   if (allocated(from%array3)) then
      call move_alloc(from%array3, to%array3)
   end if
end subroutine move_record

pure subroutine copy_record(to, from)
   class(double_record), intent(inout) :: to
   type(double_record), intent(in) :: from

   if (allocated(to%label)) deallocate(to%label)
   to%label = from%label
   if (allocated(from%array1)) to%array1 = from%array1
   if (allocated(from%array2)) to%array2 = from%array2
   if (allocated(from%array3)) to%array3 = from%array3
end subroutine

pure subroutine concatenate_overwrite(self, dict2)
   class(double_dictionary_type), intent(inout) :: self
   type(double_dictionary_type), intent(in) :: dict2
   type(double_dictionary_type) :: tmp_dict
   tmp_dict = self + dict2
   self = tmp_dict

end subroutine

pure function combine_dict(self, dict2) result(new_dict)
   class(double_dictionary_type), intent(in) :: self
   type(double_dictionary_type), intent(in) :: dict2
   type(double_dictionary_type) :: new_dict

   integer :: it, i, n_entries
   new_dict = self
   associate(dict => dict2)
      n_entries = dict%n
      do i = 1, n_entries
         call new_dict%push(dict%record(i)%label, it)
         new_dict%record(it) = dict%record(i)
      end do
   end associate

end function

pure subroutine update_1d(self, label, array)
   class(double_dictionary_type), intent(inout) :: self
   character(len=*), intent(in) :: label
   real(wp), intent(in) :: array(:)

   integer :: it

   it = return_label_index(self, label)
   if (it /= 0) then
      associate(record => self%record(it))
         if (allocated(record%array1)) deallocate(record%array1)
         if (allocated(record%array2)) deallocate(record%array2)
         if (allocated(record%array3)) deallocate(record%array3)
         record%array1 = array
      end associate
   else
      return
   end if
end subroutine

pure subroutine update_2d(self, label, array)
   class(double_dictionary_type), intent(inout) :: self
   character(len=*), intent(in) :: label
   real(wp), intent(in) :: array(:, :)
   integer :: it

   it = return_label_index(self, label)
   if (it /= 0) then
      associate(record => self%record(it))
         if (allocated(record%array1)) deallocate(record%array1)
         if (allocated(record%array2)) deallocate(record%array2)
         if (allocated(record%array3)) deallocate(record%array3)
         record%array2 = array
      end associate
   else
      return
   end if
end subroutine

pure subroutine update_3d(self, label, array)
   class(double_dictionary_type), intent(inout) :: self
   character(len=*), intent(in) :: label
   real(wp), intent(in) :: array(:, :, :)
   integer :: it

   it = return_label_index(self, label)
   if (it /= 0) then
      associate(record => self%record(it))
         if (allocated(record%array1)) deallocate(record%array1)
         if (allocated(record%array2)) deallocate(record%array2)
         if (allocated(record%array3)) deallocate(record%array3)
         record%array3 = array
      end associate
   end if
end subroutine update_3d

!> Get the label of a given index
pure subroutine get_label(self, index, label)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(in) :: self
   !> Index of the label in the dictionary
   integer, intent(in) :: index
   !> Label to be returned
   character(len=:), allocatable, intent(out) :: label

   if (index > self%n) return

   label = self%record(index)%label
end subroutine get_label

!> Get a 1D array from the dictionary
pure subroutine get_1d_label(self, label, array)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(in) :: self
   !> Label for the array
   character(len=*), intent(in) :: label
   !> 1D array to be returned
   real(wp), allocatable, intent(out) :: array(:)
   integer :: it
   it = return_label_index(self, label)
   if (it == 0) return
   call self%get_entry(it, array)

end subroutine get_1d_label

!> Get a 2D array from the dictionary
pure subroutine get_2d_label(self, label, array)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(in) :: self
   !> Label for the array
   character(len=*), intent(in) :: label
   !> 2D array to be returned
   real(wp), allocatable, intent(out) :: array(:,:)
   integer :: it
   it = return_label_index(self, label)
   if (it == 0) return
   call self%get_entry(it, array)
end subroutine get_2d_label

!> Get a 3D array from the dictionary
pure subroutine get_3d_label(self, label, array)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(in) :: self
   !> Label for the array
   character(len=*), intent(in) :: label
   !> 3D array to be returned
   real(wp), allocatable, intent(out) :: array(:, :, :)
   integer :: it
   it = return_label_index(self, label)
   if (it == 0) return
   call self%get_entry(it, array)
end subroutine get_3d_label

!> Get a 1D array from the dictionary
pure subroutine get_1d_index(self, index, array)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(in) :: self
   !> Index of the array in the dictionary
   integer, intent(in) :: index
   !> 1D array to be returned
   real(wp), allocatable, intent(out) :: array(:)

   if (index > self%n) return
   if (index <= 0) return
   associate(rec => self%record(index))
      if (allocated(rec%array1)) then
         array = rec%array1
      end if
   end associate

end subroutine get_1d_index

!> Get a 2D array from the dictionary
pure subroutine get_2d_index(self, index, array)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(in) :: self
   !> Index of the array in the dictionary
   integer, intent(in) :: index
   !> 2D array to be returned
   real(wp), allocatable, intent(out) :: array(:,:)

   if (index > self%n) return
   if (index <= 0) return
   associate(rec => self%record(index))
      if (allocated(rec%array2)) then
         array = rec%array2
      end if
   end associate

end subroutine get_2d_index

!> Get a 3D array from the dictionary
pure subroutine get_3d_index(self, index, array)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(in) :: self
   !> Index of the array in the dictionary
   integer, intent(in) :: index
   !> 3D array to be returned
   real(wp), allocatable, intent(out) :: array(:, :, :)

   if (index > self%n) return
   if (index <= 0) return
   associate(rec => self%record(index))
      if (allocated(rec%array3)) then
         array = rec%array3
      end if
   end associate

end subroutine get_3d_index

!> Add a new entry to the dictionary
pure subroutine push(self, label, it)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(inout) :: self
   !> Label for the array
   character(len=*), intent(in) :: label
   !> Index of the label in the dictionary
   integer, intent(out) :: it

   if (.not.allocated(self%record)) call resize(self%record)
   it = find(self%record(:self%n), label)

   if (it == 0) then
      if (self%n >= size(self%record)) then
         call resize(self%record)
      end if

      self%n = self%n + 1
      it = self%n
      self%record(it)%label = label
   end if
end subroutine push

!> Initialize a label in the dictionary
pure subroutine ini_label(self, label)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(inout) :: self
   !> Label for the array
   character(len=*), intent(in) :: label
   integer :: it

   call self%push(label, it)
end subroutine ini_label

!> Initialize a 1D array in the dictionary
pure subroutine ini_1d(self, label, ndim1)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(inout) :: self
   !> Label for the array
   character(len=*), intent(in) :: label
   !> Dimensions of the 1D array
   integer, intent(in) :: ndim1
   integer :: it

   call self%push(label, it)

   associate(record => self%record(it))
      allocate(record%array1(ndim1), source=0.0_wp)
   end associate
end subroutine ini_1d

!> Initialize a 2D array in the dictionary
pure subroutine ini_2d(self, label, ndim1, ndim2)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(inout) :: self
   !> Label for the array
   character(len=*), intent(in) :: label
   !> Dimensions of the 2D array
   integer, intent(in) :: ndim1, ndim2

   integer :: it

   call self%push(label, it)

   associate(record => self%record(it))
      allocate(record%array2(ndim1, ndim2), source=0.0_wp)
   end associate
end subroutine ini_2d

!> Initialize a 3D array in the dictionary
pure subroutine ini_3d(self, label, ndim1, ndim2, ndim3)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(inout) :: self
   !> Label for the array
   character(len=*), intent(in) :: label
   !> Dimensions of the 3D array
   integer, intent(in) :: ndim1, ndim2, ndim3
   integer :: it

   call self%push(label, it)

   associate(record => self%record(it))
      allocate(record%array3(ndim1, ndim2, ndim3), source=0.0_wp)
   end associate
end subroutine ini_3d

!> Add a 1D array to the dictionary
pure subroutine add_1d(self, label, array)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(inout) :: self
   !> Label for the array
   character(len=*), intent(in) :: label
   !> 1D array to be added
   real(wp), intent(in) :: array(:)
   integer :: it

   call self%push(label, it)

   associate(record => self%record(it))
      record%array1 = array
      if (allocated(record%array2)) deallocate(record%array2)
      if (allocated(record%array3)) deallocate(record%array3)
   end associate
end subroutine add_1d

!> Add a 2D array to the dictionary
pure subroutine add_2d(self, label, array)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(inout) :: self
   !> Label for the array
   character(len=*), intent(in) :: label
   !> 2D array to be added
   real(wp), intent(in) :: array(:,:)
   integer :: it

   call self%push(label, it)

   associate(record => self%record(it))
      record%array2 = array
      if (allocated(record%array1)) deallocate(record%array1)
      if (allocated(record%array3)) deallocate(record%array3)
   end associate

end subroutine add_2d

!> Add a 3D array to the dictionary
pure subroutine add_3d(self, label, array)
   !> Instance of the double dictionary
   class(double_dictionary_type), intent(inout) :: self
   !> Label for the array
   character(len=*), intent(in) :: label
   !> 3D array to be added
   real(wp), intent(in) :: array(:, :, :)
   integer :: it

   call self%push(label, it)

   associate(record => self%record(it))
      record%array3 = array
      if (allocated(record%array1)) deallocate(record%array1)
      if (allocated(record%array2)) deallocate(record%array2)
   end associate
end subroutine add_3d

!> Find the index of a label in the dictionary
pure function return_label_index(self, label) result(it)
   class(double_dictionary_type), intent(in) :: self
   character(len=*), intent(in) :: label
   integer :: it
   it = 0
   if (self%n <= 0) return
   it = find(self%record(:self%n), label)
end function return_label_index

!> Find the index of a label in the dictionary
pure function find(record, label) result(pos)
   !> List of double arrays
   type(double_record), intent(in) :: record(:)
   !> Label to search for
   character(len=*), intent(in) :: label
   !> Position of the label in the list
   integer :: pos

   integer :: i

   pos = 0

   do i = size(record), 1, -1
      if (allocated(record(i)%label)) then
         if (label == record(i)%label) then
            pos = i
            exit
         end if
      end if
   end do

end function find

!> Reallocate list of double arrays
pure subroutine resize(var, n)
   !> Instance of the array to be resized
   type(double_record), allocatable, intent(inout) :: var(:)
   !> Dimension of the final array size
   integer, intent(in), optional :: n

   type(double_record), allocatable :: tmp(:)
   integer :: this_size, new_size, it
   integer, parameter :: initial_size = 20

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
      do it = 1, this_size
         call move_record(tmp(it), var(it))
      end do
      deallocate(tmp)
   end if

end subroutine resize

end module tblite_double_dictionary