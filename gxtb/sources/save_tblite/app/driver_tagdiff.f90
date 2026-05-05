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

!> Implementation of the driver entry points for handling data outputs
module tblite_driver_tagdiff
   use, intrinsic :: iso_fortran_env, only : output_unit
   use mctc_env, only : error_type, fatal_error, wp
   use tblite_cli, only : tagdiff_config
   use tblite_io_tag, only : tagged_data, tagged_entry
   use tblite_output_format, only : format_string
   implicit none
   private

   public :: main

   interface main
      module procedure :: tagdiff_main
   end interface



contains


subroutine tagdiff_main(config, error)
   type(tagdiff_config), intent(in) :: config
   type(error_type), allocatable, intent(out) :: error

   type(tagged_data), target :: actual, reference
   type(tagged_entry), pointer :: ptr
   integer :: it, iv, stat

   call actual%load(config%actual, stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read data from '"//config%actual//"'")
      return
   end if
   call reference%load(config%reference, stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read data from '"//config%reference//"'")
      return
   end if

   stat = 0
   do it = 1, size(reference%val)
      call actual%get(reference%val(it)%tag, ptr)
      if (.not.associated(ptr)) then
         if (.not.config%fit) then
            write(output_unit, '(a)') &
               & "Data field '"//reference%val(it)%tag//"' not found in actual data"
         end if
         stat = stat + 1
         cycle
      end if
      select type(aval => ptr%raw)
      type is(real(wp))
         select type(rval => reference%val(it)%raw)
         type is(real(wp))
            if (config%fit) then
               do iv = 1, size(rval)
                  write(output_unit, '(2es24.16)') rval(iv), aval(iv)
               end do
            else
               write(output_unit, '(a)') &
                  & "Difference in '"//reference%val(it)%tag//"' is "//&
                  & format_string(norm2(rval - aval), '(f10.6)')
            end if
         end select
      end select
   end do
   if (stat /= 0) then
      call fatal_error(error, "Data fields from reference are missing in actual data")
      return
   end if

end subroutine tagdiff_main


end module tblite_driver_tagdiff
