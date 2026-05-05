! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Implementation of the argument list processor.
module dftd4_argument
   implicit none
   private

   public :: argument_list, len
   public :: argument_count_interface, get_argument_interface


   !> Internal representation of the command line arguments
   type :: argument
      private
      !> Actual payload of the argument
      character(len=:), allocatable :: raw
   end type argument

   !> Constructor for the argument representation
   interface argument
      module procedure :: new_argument
   end interface argument


   !> Argument list class
   type :: argument_list
      !> Name of the invoked executable, if available
      character(len=:), allocatable :: prog
      !> Number of arguments
      integer :: nargs = 0
      !> Array of arguments in this list
      type(argument), allocatable :: argv(:)
   contains
      !> Append a command line argument
      procedure :: push_back
      !> Display debug information on this instance
      procedure :: info
      !> Get command line argument
      procedure :: get
   end type argument_list

   !> Constructor for the argument list
   interface argument_list
      module procedure :: new_argument_list
   end interface argument_list

   interface len
      module procedure :: get_length
   end interface len


   abstract interface
      !> Interface of the argument counter
      function argument_count_interface() result(argument_count)
         !> Number of available arguments
         integer :: argument_count
      end function argument_count_interface

      !> Interface of the argument getter
      subroutine get_argument_interface(idx, arg)
         !> Index of the argument to retrieve, range 0 to argument_counter()
         integer, intent(in) :: idx
         !> Returned argument payload, allocation status is used to signal errors
         character(len=:), allocatable, intent(out) :: arg
      end subroutine get_argument_interface
   end interface

   !> Token identifyin response files
   character(len=*), parameter :: response_token = "@"

   !> Initial size for the resizer of the argument list
   integer, parameter :: initial_size = 20


contains


!> Create a new argument from a raw payload
pure function new_argument(raw) result(new)
   !> Raw argument value
   character(len=*), intent(in) :: raw
   !> Representation of the argument
   type(argument) :: new

   new%raw = raw
end function new_argument

!> Constructor of the argument list
function new_argument_list(argument_counter, argument_getter) result(new)
   !> Argument counter interface
   procedure(argument_count_interface), optional :: argument_counter
   !> Argument getter interface
   procedure(get_argument_interface), optional :: argument_getter
   !> Newly created argument list
   type(argument_list) :: new

   intrinsic :: present

   if (present(argument_getter) .and. present(argument_counter)) then
      call make_argument_list(new, argument_counter, argument_getter)
   else
      call make_argument_list(new, default_argument_count, get_default_argument)
   end if
end function new_argument_list

!> Internal constructor of the argument list
subroutine make_argument_list(self, argument_counter, argument_getter)
   !> Instance of the argument list to be created
   type(argument_list), intent(out) :: self
   !> Argument counter interface
   procedure(argument_count_interface) :: argument_counter
   !> Argument getter interface
   procedure(get_argument_interface) :: argument_getter

   integer :: iarg, narg, info
   character(len=:), allocatable :: arg
   intrinsic :: allocated

   info = 0
   narg = argument_counter()
   self%nargs = 0
   call resize(self%argv, narg)
   call argument_getter(0, self%prog)
   do iarg = 1, narg
      call argument_getter(iarg, arg)
      if (.not.allocated(arg)) return
      if (is_response_file(arg)) then
         call get_response_file(self, arg(2:), info)
         if (info == 0) cycle
      end if
      call push_back(self, arg)
   end do

end subroutine make_argument_list


!> Check if an argument represents a response file
pure function is_response_file(arg) result(is_resp)
   !> Argument of interest
   character(len=*), intent(in) :: arg
   !> Whether the argument could be a response file or not
   logical :: is_resp
   intrinsic :: len

   if (len(arg) > 1) then
      is_resp = arg(1:1) == response_token
   else
      is_resp = .false.
   end if
end function is_response_file

!> Recursively consume a response file and append it to the argument list
recursive subroutine get_response_file(self, resp, stat)
   !> Instance of the argument list
   class(argument_list), intent(inout) :: self
   !> Name of the response file to be appended
   character(len=*), intent(in) :: resp
   !> Status of reading the reponse file
   integer, intent(out) :: stat

   integer :: unit, info, istat
   logical :: opened
   character(len=:), allocatable :: arg

   inquire(file=resp, opened=opened)
   if (opened) then
      stat = 1
      return
   end if

   open(file=resp, unit=unit, iostat=info, status='old', action='read')
   do while(info == 0)
      call getline(unit, arg, info)
      if (info /= 0) exit
      if (is_response_file(arg)) then
         call get_response_file(self, arg(2:), istat)
         if (istat == 0) cycle
      end if
      call push_back(self, arg)
   end do
   close(unit, iostat=stat)
   if (info /= 0) then
      stat = merge(0, info, is_iostat_end(info))
   end if
end subroutine get_response_file

!> Consume a whole line from a formatted unit
subroutine getline(unit, line, iostat, iomsg)
   !> Formatted IO unit
   integer, intent(in) :: unit
   !> Line to read
   character(len=:), allocatable, intent(out) :: line
   !> Status of operation
   integer, intent(out) :: iostat
   !> Error message
   character(len=:), allocatable, optional :: iomsg

   integer, parameter :: bufsize = 512
   character(len=bufsize) :: buffer, msg
   integer :: size, stat
   intrinsic :: is_iostat_eor, present, trim

   allocate(character(len=0) :: line)
   do
      read(unit, '(a)', advance='no', iostat=stat, iomsg=msg, size=size) &
         & buffer
      if (stat > 0) exit
      line = line // buffer(:size)
      if (stat < 0) exit
   end do

   if (is_iostat_eor(stat)) stat = 0

   if (stat /= 0) then
      if (present(iomsg)) iomsg = trim(msg)
   end if
   iostat = stat

end subroutine getline


!> Append a string to the argument list
subroutine push_back(self, string)
   !> Instance of the argument list
   class(argument_list), intent(inout) :: self
   !> String representing the argument
   character(len=*), intent(in) :: string
   intrinsic :: size

   self%nargs = self%nargs + 1
   if (self%nargs > size(self%argv)) call resize(self%argv)
   self%argv(self%nargs) = argument(string)

end subroutine push_back

!> Reallocate list of arguments
pure subroutine resize(list, n)
   !> Instance of the array to be resized
   type(argument), allocatable, intent(inout) :: list(:)
   !> Dimension of the final array size
   integer, intent(in), optional :: n

   type(argument), allocatable :: tmp(:)
   integer :: this_size, new_size, iv
   intrinsic :: allocated, size, move_alloc, present, min

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
      do iv = 1, this_size
         call move_alloc(tmp(iv)%raw, list(iv)%raw)
      end do
      deallocate(tmp)
   end if
end subroutine resize


!> Display debug information on an argument list instance
subroutine info(self, unit)
   !> Instance of the argument list
   class(argument_list), intent(in) :: self
   !> Formatted unit for output
   integer, intent(in) :: unit

   character(len=*), parameter :: fmt = '("#", *(1x, g0))'
   integer :: iarg
   intrinsic :: allocated

   if (allocated(self%prog)) then
      write(unit, fmt) self%prog
   end if

   if (allocated(self%argv)) then
      write(unit, fmt) self%nargs, "arguments provided"
      do iarg = 1, self%nargs
         write(unit, fmt) iarg, "/", self%nargs, "->", self%argv(iarg)%raw
      end do
   end if
end subroutine info


!> Default argument counter using the intrinsic command_argument_count procedure
function default_argument_count() result(argument_count)
   !> Number of available arguments
   integer :: argument_count

   intrinsic :: command_argument_count

   argument_count = command_argument_count()
end function default_argument_count

!> Default argument getter using the intrinsic get_command_argument procedure
subroutine get_default_argument(idx, arg)
   !> Index of the argument to retrieve, range 0 to argument_counter()
   integer, intent(in) :: idx
   !> Returned argument payload, allocation status is used to signal errors
   character(len=:), allocatable, intent(out) :: arg

   integer :: length, stat
   intrinsic :: get_command_argument

   call get_command_argument(idx, length=length, status=stat)
   if (stat /= 0) then
      return
   endif

   allocate(character(len=length) :: arg, stat=stat)
   if (stat /= 0) then
      return
   endif

   if (length > 0) then
      call get_command_argument(idx, arg, status=stat)
      if (stat /= 0) then
         deallocate(arg)
         return
      end if
   end if
end subroutine get_default_argument


pure subroutine get(self, idx, arg)
   class(argument_list), intent(in) :: self
   character(len=:), allocatable, intent(out) :: arg
   integer, intent(in) :: idx

   if (idx > 0 .and. idx <= self%nargs) arg = self%argv(idx)%raw
end subroutine get

pure function get_length(self) result(length)
   class(argument_list), intent(in) :: self
   integer :: length
   length = self%nargs
end function get_length

end module dftd4_argument
