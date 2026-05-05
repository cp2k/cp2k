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

!> @file tblite/wavefunction/fermi.f90
!> Provides a Fermi based electronic filling function

!> Implementation of a Fermi distribution for filling the electronic levels
module tblite_wavefunction_fermi
   use mctc_env, only : wp
   implicit none
   private

   public :: get_fermi_filling


contains

subroutine get_fermi_filling(nel, kt, emo, homo, focc, e_fermi)
   real(wp), intent(in) :: nel
   real(wp), intent(in) :: emo(:)
   real(wp), intent(in) :: kt
   integer, intent(out) :: homo
   real(wp), intent(out) :: focc(:)
   real(wp), intent(out) :: e_fermi

   real(wp) :: etmp

   e_fermi = 0.0_wp

   call get_aufbau_filling(nel, homo, focc)

   ! Optimize the Fermilevel if there is a finite temperature
   if (homo > 0 .and. kt .gt. 0) then
      call get_fermi_filling_(homo, kt, emo, focc, etmp)
      e_fermi = 0.5_wp * etmp
   end if

end subroutine get_fermi_filling

subroutine get_aufbau_filling(nel, homo, occ)
   real(wp), intent(in) :: nel
   integer, intent(out) :: homo
   real(wp), intent(out) :: occ(:)

   occ(:) = 0.0_wp
   homo = floor(nel)
   occ(:min(homo, size(occ))) = 1.0_wp
   if (homo < size(occ)) occ(homo+1) = mod(nel, 1.0_wp)
   homo = merge(homo+1, homo, mod(nel, 1.0_wp) > 0.5_wp)
end subroutine get_aufbau_filling

subroutine get_fermi_filling_(homo, kt, emo, occ, e_fermi)
   integer, intent(in) :: homo
   real(wp), intent(in) :: emo(:)
   real(wp), intent(in) :: kt
   real(wp), intent(out) :: occ(:)
   real(wp), intent(out) :: e_fermi

   real(wp) :: occt, total_number
   real(wp) :: total_dfermi, dfermifunct, fermifunct, change_fermi
   integer, parameter :: max_cycle = 200
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: sqrttiny = sqrt(tiny(1.0_wp))
   integer :: ncycle, iao

   e_fermi = 0.5*(emo(max(homo, 1))+emo(min(homo+1, size(emo))))
   occt = homo

   do ncycle = 1, max_cycle
      total_number = 0.0
      total_dfermi = 0.0
      do iao = 1, size(emo)
         fermifunct = 0.0
         if((emo(iao)-e_fermi)/kt.lt.50) then
            fermifunct = 1.0/(exp((emo(iao)-e_fermi)/kt)+1.0)
            dfermifunct = exp((emo(iao)-e_fermi)/kt) / &
               & (kt*(exp((emo(iao)-e_fermi)/kt)+1.0)**2)
         else
            dfermifunct = 0.0
         end if
         occ(iao) = fermifunct
         total_number = total_number + fermifunct
         total_dfermi = total_dfermi + dfermifunct
      end do
      if (total_dfermi > sqrttiny) then
          change_fermi = (occt-total_number)/total_dfermi
      else
          change_fermi = 0.0_wp
      end if
      e_fermi = e_fermi+change_fermi
      if (abs(occt-total_number) <= thr) exit
   end do

end subroutine get_fermi_filling_

end module tblite_wavefunction_fermi
