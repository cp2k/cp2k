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

!> @file tblite/output/property.f90
!> Provides helper routines to print properties together with units

!> Simple derived type with formatted IO procedure to provide simple debug printing
module tblite_output_property
   use mctc_env, only : wp
   implicit none
   private

   public :: write(formatted)

   type, public :: property
      character(len=:), allocatable :: label
      real(wp) :: value
      character(len=:), allocatable :: unit
   end type property

   interface property
      module procedure :: new_property
   end interface property

    interface write(formatted)
        module procedure :: write_formatted
    end interface


contains


pure function new_property(label, value, unit) result(new)
   character(len=*), intent(in) :: label
   real(wp), intent(in) :: value
   character(len=*), intent(in) :: unit
   type(property) :: new

   new%label = label
   new%value = value
   new%unit = unit
end function new_property


subroutine write_formatted(prop, unit, iotype, v_list, iostat, iomsg)
   class(property), intent(in) :: prop
   integer, intent(in) :: unit
   character(len=*), intent(in) :: iotype
   integer, intent(in) :: v_list(:)
   integer, intent(out) :: iostat
   character(len=*), intent(inout) :: iomsg

   write(unit, '(a, t25, es20.13, 1x, a)') prop%label, prop%value, prop%unit
end subroutine write_formatted


end module tblite_output_property
