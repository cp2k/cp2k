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

!> @file tblite/utils/string.f90
!> Contains string conversion utilities.

module tblite_utils_string
   implicit none
   private

   public :: compact, lowercase

contains

!> Convert input string to lowercase
pure function compact(str) result(lcstr)

   !> Input string
   character(len=*), intent(in) :: str

   !> Lowercase version of string
   character(len=:), allocatable :: lcstr

   integer :: ilen, i, j, iav
   integer, parameter :: offset = iachar('A') - iachar('a')

   ilen = len(str)
   lcstr = ''

   j = 0
   do i = 1, ilen
      iav = iachar(str(i:i))
      if (any(iav == [iachar('-'), iachar(','), iachar(' ')])) cycle
      j = j + 1
      if (iav >= iachar('A') .and. iav <= iachar('Z')) then
         lcstr = lcstr // achar(iav - offset)
      else
         lcstr = lcstr // str(i:i)
      end if
   end do
   lcstr = trim(lcstr)

end function compact


!> Convert string to lower case
pure function lowercase(str) result(lcstr)
   character(len=*), intent(in)  :: str
   character(len=len_trim(str)) :: lcstr
   integer :: ilen, ioffset, iquote, i, iav, iqc

   ilen=len_trim(str)
   ioffset=iachar('A')-iachar('a')
   iquote=0
   lcstr=str
   do i=1, ilen
      iav=iachar(str(i:i))
      if(iquote==0 .and. (iav==34 .or.iav==39)) then
         iquote=1
         iqc=iav
        cycle
      endif
      if(iquote==1 .and. iav==iqc) then
         iquote=0
         cycle
      endif
      if (iquote==1) cycle
      if(iav >= iachar('A') .and. iav <= iachar('Z')) then
         lcstr(i:i)=achar(iav-ioffset)
      else
         lcstr(i:i)=str(i:i)
      endif
   enddo

end function lowercase

end module tblite_utils_string