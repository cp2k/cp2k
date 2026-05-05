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

module dftd4_utils
   use mctc_env, only : wp
   use mctc_io_math, only : matinv_3x3
   implicit none

   public :: lowercase, wrap_to_central_cell


contains


subroutine wrap_to_central_cell(xyz, lattice, periodic)
   !DEC$ ATTRIBUTES DLLEXPORT :: wrap_to_central_cell
   real(wp), intent(inout) :: xyz(:, :)
   real(wp), intent(in) :: lattice(:, :)
   logical, intent(in) :: periodic(:)
   real(wp) :: invlat(3, 3), vec(3)
   integer :: iat

   if (.not.any(periodic)) return

   invlat = matinv_3x3(lattice)
   do iat = 1, size(xyz, 2)
      vec(:) = matmul(invlat, xyz(:, iat))
      vec(:) = shift_back_abc(vec)
      xyz(:, iat) = matmul(lattice, vec)
   end do

end subroutine wrap_to_central_cell


elemental function shift_back_abc(in) result(out)
   !> fractional coordinate in (-∞,+∞)
   real(wp),intent(in) :: in
   !> fractional coordinate in [0,1)
   real(wp) :: out
   real(wp),parameter :: p_pbc_eps = 1.0e-14_wp
   out = in
   if(in < (0.0_wp - p_pbc_eps)) &
      out = in + real(ceiling(-in),wp)
   if(in > (1.0_wp + p_pbc_eps)) &
      out = in - real(floor  ( in),wp)
   if (abs(in - 1.0_wp) < p_pbc_eps) &
      out = in - 1.0_wp
end function shift_back_abc


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


end module dftd4_utils
