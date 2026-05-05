! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

module dftd3_model
   use ieee_arithmetic, only : ieee_is_nan
   use dftd3_data, only : get_r4r2_val, get_vdw_rad
   use dftd3_reference
   use mctc_data, only : get_covalent_rad
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none
   private

   public :: d3_model, new_d3_model


   !> Base D3 dispersion model to evaluate C6 coefficients
   type :: d3_model

      !> Covalent radii for coordination number
      real(wp), allocatable :: rcov(:)

      !> Van-der-Waals radii for damping function
      real(wp), allocatable :: rvdw(:, :)

      !> Expectation values for C8 extrapolation
      real(wp), allocatable :: r4r2(:)

      !> Weighting factor for CN interpolation
      real(wp), allocatable :: wf

      !> Number of reference systems
      integer, allocatable :: ref(:)

      !> Reference coordination numbers
      real(wp), allocatable :: cn(:, :)

      !> Reference C6 coefficients
      real(wp), allocatable :: c6(:, :, :, :)

   contains

      !> Generate weights for all reference systems
      procedure :: weight_references

      !> Evaluate C6 coefficient
      procedure :: get_atomic_c6

   end type d3_model


   !> Default weighting factor for coordination number interpolation
   real(wp), parameter :: wf_default = 4.0_wp


contains


!> Create new dispersion model from molecular structure input
subroutine new_d3_model(self, mol, wf)

   !> Instance of the dispersion model
   type(d3_model), intent(out) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Weighting factor for coordination number interpolation
   real(wp), intent(in), optional :: wf

   integer :: isp, izp, iref, jsp, jzp, jref
   integer :: mref

   call init_reference_c6

   if (present(wf)) then
      self%wf = wf
   else
      self%wf = wf_default
   end if

   allocate(self%ref(mol%nid))
   self%ref(:) = number_of_references(mol%num)
   mref = maxval(self%ref)

   allocate(self%rcov(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      self%rcov(isp) = get_covalent_rad(izp)
   end do

   allocate(self%r4r2(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      self%r4r2(isp) = get_r4r2_val(izp)
   end do

   allocate(self%rvdw(mol%nid, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, isp
         jzp = mol%num(jsp)
         self%rvdw(jsp, isp) = get_vdw_rad(jzp, izp)
         self%rvdw(isp, jsp) = self%rvdw(jsp, isp)
      end do
   end do

   allocate(self%cn(mref, mol%nid), source=0.0_wp)
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do iref = 1, self%ref(isp)
         self%cn(iref, isp) = reference_cn(iref, izp)
      end do
   end do

   allocate(self%c6(mref, mref, mol%nid, mol%nid), source=0.0_wp)
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, isp
         jzp = mol%num(jsp)
         do iref = 1, self%ref(isp)
            do jref = 1, self%ref(jsp)
               self%c6(jref, iref, jsp, isp) = get_c6(jref, iref, jzp, izp)
               self%c6(iref, jref, isp, jsp) = self%c6(jref, iref, jsp, isp)
            end do
         end do
      end do
   end do

end subroutine new_d3_model


!> Calculate the weights of the reference system and the derivatives w.r.t.
!> coordination number for later use.
subroutine weight_references(self, mol, cn, gwvec, gwdcn)

   !> Instance of the dispersion model
   class(d3_model), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Coordination number of every atom.
   real(wp), intent(in) :: cn(:)

   !> weighting for the atomic reference systems
   real(wp), intent(out) :: gwvec(:, :)

   !> derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(out), optional :: gwdcn(:, :)

   integer :: iat, izp, iref, mref
   real(wp) :: norm, dnorm, gw, expw, expd, gwk, dgwk
   real(wp) :: wf2, dcn
   real(wp), allocatable :: gwt(:)

   mref = maxval(self%ref)

   if (present(gwdcn)) then
      gwvec(:, :) = 0.0_wp
      gwdcn(:, :) = 0.0_wp
      wf2 = 2 * self%wf

      !$omp parallel default(none) &
      !$omp shared(gwvec, gwdcn, mol, self, cn, wf2, mref) &
      !$omp private(iat, izp, iref, norm, dnorm, gw, expw, expd, gwk, dgwk, dcn, gwt)
      allocate(gwt(mref))
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         norm = 0.0_wp
         dnorm = 0.0_wp
         do iref = 1, self%ref(izp)
            gw = weight_cn(self%wf, cn(iat), self%cn(iref, izp))
            gwt(iref) = gw
            norm = norm + gw
            dnorm = dnorm + wf2 * (self%cn(iref, izp) - cn(iat)) * gw
         end do
         norm = 1.0_wp / norm
         do iref = 1, self%ref(izp)
            expw = gwt(iref)
            dcn = self%cn(iref, izp) - cn(iat)
            expd = wf2 * dcn * expw
            gwk = expw * norm
            if (is_exceptional(gwk)) then
               if (maxval(self%cn(:self%ref(izp), izp)) == self%cn(iref, izp)) then
                  gwk = 1.0_wp
               else
                  gwk = 0.0_wp
               end if
            end if
            gwvec(iref, iat) = gwk

            dgwk = expd * norm - expw * dnorm * norm**2
            if (is_exceptional(dgwk)) then
               dgwk = 0.0_wp
            end if
            gwdcn(iref, iat) = dgwk
         end do
      end do
      !$omp end do
      deallocate(gwt)
      !$omp end parallel

   else

      gwvec(:, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(gwvec, mol, self, cn, mref) &
      !$omp private(iat, izp, iref, norm, gw, expw, gwk, gwt)
      allocate(gwt(mref))
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         norm = 0.0_wp
         do iref = 1, self%ref(izp)
            gw = weight_cn(self%wf, cn(iat), self%cn(iref, izp))
            gwt(iref) = gw
            norm = norm + gw
         end do
         norm = 1.0_wp / norm
         do iref = 1, self%ref(izp)
            gwk = gwt(iref) * norm
            if (is_exceptional(gwk)) then
               if (maxval(self%cn(:self%ref(izp), izp)) == self%cn(iref, izp)) then
                  gwk = 1.0_wp
               else
                  gwk = 0.0_wp
               end if
            end if
            gwvec(iref, iat) = gwk
         end do
      end do
      !$omp end do
      deallocate(gwt)
      !$omp end parallel
   end if

end subroutine weight_references


!> Check whether we are dealing with an exceptional value, NaN or Inf
elemental function is_exceptional(val)
   real(wp), intent(in) :: val
   logical :: is_exceptional
   is_exceptional = ieee_is_nan(val) .or. abs(val) > huge(val)
end function is_exceptional


!> Calculate atomic dispersion coefficients and their derivatives w.r.t.
!> the coordination number.
subroutine get_atomic_c6(self, mol, gwvec, gwdcn, c6, dc6dcn)

   !> Instance of the dispersion model
   class(d3_model), intent(in) :: self

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :)

   !> Derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in), optional :: gwdcn(:, :)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(out) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(out), optional :: dc6dcn(:, :)

   integer :: iat, jat, izp, jzp, iref, jref
   real(wp) :: refc6, dc6, dc6dcni, dc6dcnj

   if (present(gwdcn).and.present(dc6dcn)) then
      c6(:, :) = 0.0_wp
      dc6dcn(:, :) = 0.0_wp

      !$omp parallel do schedule(runtime) default(none) &
      !$omp shared(c6, dc6dcn, mol, self, gwvec, gwdcn) &
      !$omp private(iat, jat, izp, jzp, iref, jref, refc6, dc6, dc6dcni, dc6dcnj)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            dc6 = 0.0_wp
            dc6dcni = 0.0_wp
            dc6dcnj = 0.0_wp
            do iref = 1, self%ref(izp)
               do jref = 1, self%ref(jzp)
                  refc6 = self%c6(iref, jref, izp, jzp)
                  dc6 = dc6 + gwvec(iref, iat) * gwvec(jref, jat) * refc6
                  dc6dcni = dc6dcni + gwdcn(iref, iat) * gwvec(jref, jat) * refc6
                  dc6dcnj = dc6dcnj + gwvec(iref, iat) * gwdcn(jref, jat) * refc6
               end do
            end do
            c6(iat, jat) = dc6
            c6(jat, iat) = dc6
            dc6dcn(iat, jat) = dc6dcni
            dc6dcn(jat, iat) = dc6dcnj
         end do
      end do

   else

      c6(:, :) = 0.0_wp

      !$omp parallel do schedule(runtime) default(none) &
      !$omp shared(c6, mol, self, gwvec) &
      !$omp private(iat, jat, izp, jzp, iref, jref, refc6, dc6)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            dc6 = 0.0_wp
            do iref = 1, self%ref(izp)
               do jref = 1, self%ref(jzp)
                  refc6 = self%c6(iref, jref, izp, jzp)
                  dc6 = dc6 + gwvec(iref, iat) * gwvec(jref, jat) * refc6
               end do
            end do
            c6(iat, jat) = dc6
            c6(jat, iat) = dc6
         end do
      end do
   end if

end subroutine get_atomic_c6


elemental function weight_cn(wf,cn,cnref) result(cngw)
   real(wp),intent(in) :: wf, cn, cnref
   real(wp) :: cngw
   intrinsic :: exp
   cngw = exp ( -wf * ( cn - cnref )**2 )
end function weight_cn


end module dftd3_model
