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

!> @dir tblite/io
!> Contains interfaces to IO formats

!> @file tblite/io/tag.f90
!> Provides tagged input and output routines

!> Implementation of tagged input and output files
module tblite_io_tag
   use mctc_io_utils, only : getline
   use mctc_env, only : sp, dp
   implicit none
   private

   public :: tagged_data, tagged_entry
   public :: write_tagged, read_tagged

   interface write_tagged
      module procedure :: write_tagged
      module procedure :: write_tagged_dp0
      module procedure :: write_tagged_dp1
      module procedure :: write_tagged_dp2
   end interface

   interface read_tagged
      module procedure :: read_tagged
   end interface

   character(len=*), parameter :: tag_header = &
      & '(a,t20,":",a,":",i0,":",*(i0:,","))'

   type :: tagged_entry
      character(len=:), allocatable :: tag
      class(*), allocatable :: raw(:)
      integer, allocatable :: shape(:)
   end type tagged_entry

   type :: tagged_data
      type(tagged_entry), allocatable :: val(:)
   contains
      procedure :: get
      !> Reading of tagged data
      generic :: load => load_from_file, load_from_unit
      !> Read tagged data from file
      procedure, private :: load_from_file
      !> Read tagged data from formatted unit
      procedure, private :: load_from_unit
      !> Writing of tagged data
      generic :: dump => dump_to_file, dump_to_unit
      !> Write tagged data to file
      procedure, private :: dump_to_file
      !> Write tagged data to formatted unit
      procedure, private :: dump_to_unit
   end type tagged_data

contains

subroutine write_tagged_dp0(unit, tag, val, stat)
   integer, intent(in) :: unit
   character(len=*), intent(in) :: tag
   real(dp), intent(in) :: val
   integer, intent(out) :: stat

   type(tagged_entry) :: tmp

   tmp%tag = tag
   tmp%raw = [val]
   allocate(tmp%shape(0))
   call write_tagged(unit, tmp, stat)
end subroutine write_tagged_dp0

subroutine write_tagged_dp1(unit, tag, val, stat)
   integer, intent(in) :: unit
   character(len=*), intent(in) :: tag
   real(dp), intent(in) :: val(:)
   integer, intent(out) :: stat

   type(tagged_entry) :: tmp

   tmp%tag = tag
   tmp%raw = val
   tmp%shape = shape(val)
   call write_tagged(unit, tmp, stat)
end subroutine write_tagged_dp1

subroutine write_tagged_dp2(unit, tag, val, stat)
   integer, intent(in) :: unit
   character(len=*), intent(in) :: tag
   real(dp), intent(in) :: val(:, :)
   integer, intent(out) :: stat

   type(tagged_entry) :: tmp

   tmp%tag = tag
   tmp%raw = reshape(val, [size(val)])
   tmp%shape = shape(val)
   call write_tagged(unit, tmp, stat)
end subroutine write_tagged_dp2

subroutine write_tagged(unit, val, stat)
   integer, intent(in) :: unit
   type(tagged_entry), intent(in) :: val
   integer, intent(out) :: stat

   select type(raw => val%raw)
   type is(real(dp))
      write(unit, tag_header, iostat=stat) val%tag, "real", size(val%shape), val%shape
      write(unit, '(3es24.16)') raw
   end select
end subroutine write_tagged

subroutine read_tagged(unit, val, stat)
   integer, intent(in) :: unit
   type(tagged_entry), intent(out) :: val
   integer, intent(out) :: stat

   integer :: idel, jdel, kdel, vrank, i
   character(len=:), allocatable :: line

   call getline(unit, line, stat)
   if (is_iostat_end(stat)) return
   idel = scan(line, ":")
   jdel = scan(line(idel+1:), ":") + idel

   val%tag = trim(line(:idel-1))

   kdel = scan(line(jdel+1:), ":") + jdel
   read(line(jdel+1:kdel-1), *) vrank
   allocate(val%shape(vrank))
   if (vrank > 0) then
      do i = kdel+1, len(line)
         if (line(i:i) == ",") line(i:i) = " "
      end do
      read(line(kdel+1:), *) val%shape
   end if

   select case(line(idel+1:jdel-1))
   case("real")
      block
         integer :: i, n
         real(dp), allocatable :: tmp(:)
         n = max(product(val%shape), 1)
         allocate(tmp(n))
         read(unit, *) (tmp(i), i = 1, n)
         call move_alloc(tmp, val%raw)
      end block
   end select
end subroutine read_tagged


subroutine load_from_unit(self, unit, stat)
   class(tagged_data), intent(out) :: self
   integer, intent(in) :: unit
   integer, intent(out) :: stat

   integer :: it

   call resize(self%val)

   it = 0
   stat = 0
   do while(stat == 0)
      if (it > size(self%val)) call resize(self%val)
      it = it + 1
      call read_tagged(unit, self%val(it), stat)
      if (stat /= 0) it = it - 1
   end do
   call resize(self%val, it)

   if (is_iostat_end(stat)) stat = 0
end subroutine load_from_unit

subroutine dump_to_unit(self, unit, stat)
   class(tagged_data), intent(in) :: self
   integer, intent(in) :: unit
   integer, intent(out) :: stat

   integer :: it

   do it = 1, size(self%val)
      call write_tagged(unit, self%val(it), stat)
      if (stat /= 0) exit
   end do
end subroutine dump_to_unit

subroutine get(self, tag, ptr)
   class(tagged_data), intent(in), target :: self
   character(len=*), intent(in) :: tag
   type(tagged_entry), pointer, intent(out) :: ptr

   integer :: i

   nullify(ptr)
   do i = 1, size(self%val)
      if (self%val(i)%tag == tag) then
         ptr => self%val(i)
         exit
      end if
   end do
end subroutine get


!> Reallocate list of tagged entries
subroutine resize(var, n)
   !> Instance of the array to be resized
   type(tagged_entry), allocatable, intent(inout) :: var(:)
   !> Dimension of the final array size
   integer, intent(in), optional :: n

   type(tagged_entry), allocatable :: tmp(:)
   integer :: this_size, new_size, i
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
      do i = 1, this_size
         call move_alloc(tmp(i)%tag, var(i)%tag)
         call move_alloc(tmp(i)%raw, var(i)%raw)
         call move_alloc(tmp(i)%shape, var(i)%shape)
      end do
      deallocate(tmp)
   end if
end subroutine resize


!> Read parametrization data from file
subroutine load_from_file(self, file, stat)
   !> Instance of the parametrization data
   class(tagged_data), intent(inout) :: self
   !> File name
   character(len=*), intent(in) :: file
   !> Error handling
   integer, intent(out) :: stat

   integer :: unit
   logical :: exist

   inquire(file=file, exist=exist)
   if (.not.exist) then
     stat = 1
     return
   end if

   open(file=file, newunit=unit)
   call self%load(unit, stat)
   close(unit)
end subroutine load_from_file


!> Write parametrization data to file
subroutine dump_to_file(self, file, stat)
   !> Instance of the parametrization data
   class(tagged_data), intent(in) :: self
   !> File name
   character(len=*), intent(in) :: file
   !> Error handling
   integer, intent(out) :: stat

   integer :: unit

   open(file=file, newunit=unit)
   call self%dump(unit, stat)
   close(unit)
   if (stat /= 0) return

end subroutine dump_to_file


end module tblite_io_tag
