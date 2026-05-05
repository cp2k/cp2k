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

!> @file tblite/container/cache.f90
!> Provides a cache for use with interaction containers

!> Definition of restart data cache
module tblite_container_cache
   implicit none
   private

   public :: resize

   !> Restart data for an interaction container
   type, public :: container_cache
      !> Label identifying this contribution
      character(len=:), allocatable :: label
      !> Actual restart data
      class(*), allocatable :: raw
   end type container_cache


   !> Reallocate list of containers
   interface resize
      module procedure :: resize_cache
   end interface


contains


!> Reallocate list of containers
subroutine resize_cache(list, n)
   !> Instance of the array to be resized
   type(container_cache), allocatable, intent(inout) :: list(:)
   !> Dimension of the final array size
   integer, intent(in), optional :: n

   type(container_cache), allocatable :: tmp(:)
   integer, parameter :: initial_size = 20
   integer :: this_size, new_size, item

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
         call move_alloc(tmp(item)%raw, list(item)%raw)
         call move_alloc(tmp(item)%label, list(item)%label)
      end do
      deallocate(tmp)
   end if

end subroutine resize_cache


end module tblite_container_cache
