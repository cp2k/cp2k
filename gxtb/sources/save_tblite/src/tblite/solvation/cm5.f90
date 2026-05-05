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

!> @file tblite/data/cm5.f90
!> Provides CM5 charges and their derivatives
!> for GFN1-xTB implementation of ALPB/GBSA

!> Implements CM5 charges and derivatives.
!> See: J. Chem. Theory Comput. 2012, 8, 2, 527–541
!> https://doi.org/10.1021/ct200866d
module tblite_solvation_cm5
   use mctc_env, only: error_type, fatal_error
   use mctc_env, only: wp
   use mctc_io,only: structure_type
   use mctc_data_atomicrad, only: get_atomic_rad
   use mctc_io_convert, only: aatoau, kcaltoau
   implicit none
   private

   public :: get_cm5_charges

   integer,parameter :: max_elements = 118

   !> CM5 model atomic parameters
   real(wp), parameter :: cm5_a0(max_elements) = (/ &
      & 0.0056_wp,-0.1543_wp, 0.0000_wp, 0.0333_wp,-0.1030_wp,-0.0446_wp, &
      &-0.1072_wp,-0.0802_wp,-0.0629_wp,-0.1088_wp, 0.0184_wp, 0.0000_wp, &
      &-0.0726_wp,-0.0790_wp,-0.0756_wp,-0.0565_wp,-0.0444_wp,-0.0767_wp, &
      & 0.0130_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      &-0.0512_wp,-0.0557_wp,-0.0533_wp,-0.0399_wp,-0.0313_wp,-0.0541_wp, &
      & 0.0092_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      &-0.0361_wp,-0.0393_wp,-0.0376_wp,-0.0281_wp,-0.0220_wp,-0.0381_wp, &
      & 0.0065_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp,-0.0255_wp,-0.0277_wp,-0.0265_wp,-0.0198_wp, &
      &-0.0155_wp,-0.0269_wp, 0.0046_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp,-0.0179_wp,-0.0195_wp, &
      &-0.0187_wp,-0.0140_wp,-0.0110_wp,-0.0189_wp /)

  !> CM5 alpha parameter
  real(wp), parameter :: cm5_alpha = 2.4740_wp/aatoau

contains

!> Get CM5 charges and derivatives
subroutine get_cm5_charges(mol, cm5, dcm5dr)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> CM5 charges
   real(wp), intent(out) :: cm5(mol%nat)
   !> CM5 charge derivatives
   real(wp), intent(out) :: dcm5dr(3, mol%nat, mol%nat)

   ! Pair parameters set up from cm5_a0
   real(wp), allocatable :: pairpar(:,:)
   integer :: iat, jat, nati, natj
   real(wp) :: pij, pji, rab(3), dist, bkk, bkkda, radi, radj

   cm5(:) = 0.0_wp
   dcm5dr(:,:,:) = 0.0_wp

   ! Set up pair parameters
   allocate (pairpar(max_elements,max_elements),source=0.0_wp)
   do iat = 1,max_elements
      do jat = iat+1,max_elements
         pairpar(iat,jat) = cm5_a0(iat)-cm5_a0(jat)
      end do
   end do
   pairpar(1,6) = 0.0502_wp
   pairpar(1,7) = 0.1747_wp
   pairpar(1,8) = 0.1671_wp
   pairpar(6,7) = 0.0556_wp
   pairpar(6,8) = 0.0234_wp
   pairpar(7,8) = -0.0346_wp
   do iat = 1,max_elements
      do jat = iat+1,max_elements
         pairpar(jat,iat) = -pairpar(iat,jat)
      end do
   end do

   ! Calculate the actual charges
   do iat = 1, mol%nat
      nati = mol%num(mol%id(iat))
      radi = get_atomic_rad(nati)
      do jat = 1,iat-1
         natj = mol%num(mol%id(jat))
         if (nati == natj) cycle
         radj = get_atomic_rad(natj)
         pij = pairpar(nati,natj)
         pji = pairpar(natj,nati)
         rab = mol%xyz(:,iat)-mol%xyz(:,jat)
         dist = norm2(rab)
         bkk = exp(-cm5_alpha*(dist-radi-radj))
         bkkda = bkk*cm5_alpha/dist
         cm5(iat) = cm5(iat)+bkk*pij
         cm5(jat) = cm5(jat)+bkk*pji
         dcm5dr(:,iat,iat) = dcm5dr(:,iat,iat)-bkkda*rab*pij
         dcm5dr(:,jat,jat) = dcm5dr(:,jat,jat)+bkkda*rab*pji
         dcm5dr(:,iat,jat) = dcm5dr(:,iat,jat)-bkkda*rab*pji
         dcm5dr(:,jat,iat) = dcm5dr(:,jat,iat)+bkkda*rab*pij
      end do
   end do

   deallocate (pairpar)

end subroutine get_cm5_charges

end module tblite_solvation_cm5
