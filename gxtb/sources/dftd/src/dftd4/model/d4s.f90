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

!> Definition of the D4S dispersion model for the evaluation of C6 coefficients.
module dftd4_model_d4s
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
   use ieee_arithmetic, only : ieee_is_nan
   use dftd4_model_type, only : dispersion_model, d4_qmod
   use dftd4_cache, only : dispersion_cache
   use dftd4_damping_type, only : twobody_damping_function, threebody_damping_function
   use dftd4_data, only : get_covalent_rad, get_r4r2_val, get_wfpair_val, &
      & get_effective_charge, get_electronegativity, get_hardness
   use dftd4_model_reference_d4
   use dftd4_model_utils, only : is_exceptional, weight_cn, zeta, dzeta
   use dftd4_integrator_type, only : integrator_type
   use dftd4_integrator_trapezoid, only : trapezoid_integrator, new_trapezoid_integrator
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use multicharge, only : new_eeq2019_model, new_eeqbc2025_model
   implicit none
   private

   public :: d4s_model, new_d4s_model


   !> D4S dispersion model to evaluate C6 coefficients
   type, extends(dispersion_model) :: d4s_model

      !> Weighting factors for CN interpolation
      real(wp), allocatable :: wf(:, :)

   contains

      !> Update cache with dispersion coefficients and properties
      procedure :: update

      !> Evaluate atomic polarizabilities from cache
      procedure :: get_polarizabilities

      !> Get two-body dispersion coefficients for an atom pair
      procedure :: get_2b_coeffs

      !> Get two-body dispersion coefficients and derivatives for an atom pair
      procedure :: get_2b_derivs

      !> Get three-body dispersion coefficients for an atom triple
      procedure :: get_3b_coeffs

      !> Get three-body dispersion coefficients and derivatives for an atom triple
      procedure :: get_3b_derivs

      !> Calculate damping radius for two-body interactions
      procedure :: get_2b_rdamp

      !> Calculate damping radius for three-body interactions
      procedure :: get_3b_rdamp

   end type d4s_model


   !> Default maximum charge scaling height for partial charge extrapolation
   real(wp), parameter :: ga_default = 3.0_wp

   !> Default charge scaling steepness for partial charge extrapolation
   real(wp), parameter :: gc_default = 2.0_wp

   !> Number of imaginary frequency integration points
   integer, parameter :: ngrid = 23

   !> Default two-body damping function for D4S
   integer, parameter :: default_damping_2b = twobody_damping_function%rational

   !> Default three-body damping function for D4S
   integer, parameter :: default_damping_3b = threebody_damping_function%zero_avg

   !> Imaginary frequencies for integration
   real(wp), parameter :: freq(ngrid) = [ &
      & 0.000001_wp, 0.050000_wp, 0.100000_wp, &
      & 0.200000_wp, 0.300000_wp, 0.400000_wp, &
      & 0.500000_wp, 0.600000_wp, 0.700000_wp, &
      & 0.800000_wp, 0.900000_wp, 1.000000_wp, &
      & 1.200000_wp, 1.400000_wp, 1.600000_wp, &
      & 1.800000_wp, 2.000000_wp, 2.500000_wp, &
      & 3.000000_wp, 4.000000_wp, 5.000000_wp, &
      & 7.500000_wp, 10.00000_wp]

   character(len=*), parameter :: d4s_label = "D4S"

contains


!> Create new D4S dispersion model from molecular structure input
subroutine new_d4s_model(error, d4, mol, ga, gc, qmod)
   !DEC$ ATTRIBUTES DLLEXPORT :: new_d4_model

   !> Instance of the dispersion model
   type(d4s_model), intent(out) :: d4

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Charge scaling height
   real(wp), intent(in), optional :: ga

   !> Charge scaling steepness
   real(wp), intent(in), optional :: gc

   !> Charge model selection
   integer, intent(in), optional :: qmod

   integer :: isp, izp, iref, jsp, jzp, jref
   integer :: mref, tmp_qmod
   real(wp) :: aiw(ngrid), c6
   real(wp), parameter :: thopi = 3.0_wp/pi
   type(trapezoid_integrator), allocatable :: integrator

   ! check for unsupported elements (104 (Rf) - 111 (Rg))
   do isp = 1, mol%nid
      if (mol%num(isp) > 103 .and. mol%num(isp) < 112) then
         call fatal_error(error, "Structure contains unsupported element '"//trim(mol%sym(isp))//"'")
         return
      end if
   end do

   d4%label = d4s_label

   d4%default_damping_2b = default_damping_2b
   d4%default_damping_3b = default_damping_3b

   d4%ncoup = mol%nat
   d4%ngrid = ngrid
   allocate(integrator)
   call new_trapezoid_integrator(integrator, d4%ngrid, freq)
   call move_alloc(integrator, d4%integrator)

   if (present(ga)) then
      d4%ga = ga
   else
      d4%ga = ga_default
   end if

   if (present(gc)) then
      d4%gc = gc
   else
      d4%gc = gc_default
   end if

   allocate(d4%wf(mol%nid, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, mol%nid
         jzp = mol%num(jsp)
         d4%wf(isp, jsp) = get_wfpair_val(izp, jzp)
      end do 
   end do

   allocate(d4%rcov(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      d4%rcov(isp) = get_covalent_rad(izp)
   end do

   allocate(d4%en(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      d4%en(isp) = get_electronegativity(izp)
   end do

   allocate(d4%zeff(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      d4%zeff(isp) = get_effective_charge(izp)
   end do

   allocate(d4%eta(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      d4%eta(isp) = get_hardness(izp)
   end do

   allocate(d4%r4r2(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      d4%r4r2(isp) = get_r4r2_val(izp)
   end do

   allocate(d4%ref(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      d4%ref(isp) = get_nref(izp)
   end do

   mref = maxval(d4%ref)
   allocate(d4%cn(mref, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      call set_refcn(d4%cn(:, isp), izp)
   end do

   if (present(qmod)) then
      tmp_qmod = qmod
   else
      tmp_qmod = d4_qmod%eeq
   end if

   allocate(d4%q(mref, mol%nid))
   allocate(d4%aiw(ngrid, mref, mol%nid))
   select case(tmp_qmod)
   case default
      call fatal_error(error, "Unsupported option for charge model.")
      return
   case(d4_qmod%eeq)
      do isp = 1, mol%nid
         izp = mol%num(isp)
         call set_refq_eeq(d4%q(:, isp), izp)
         call set_refalpha_eeq(d4%aiw(:, :, isp), d4%ga, d4%gc, izp)
      end do
      ! Setup EEQ model
      call new_eeq2019_model(mol, d4%mchrg, error)
      if(allocated(error)) return
   case(d4_qmod%eeqbc)
      do isp = 1, mol%nid
         izp = mol%num(isp)
         call set_refq_eeqbc(d4%q(:, isp), izp)
         call set_refalpha_eeqbc(d4%aiw(:, :, isp), d4%ga, d4%gc, izp)
      end do
      ! Setup EEQBC model
      call new_eeqbc2025_model(mol, d4%mchrg, error)  
      if(allocated(error)) return
   case(d4_qmod%gfn2)
      do isp = 1, mol%nid
         izp = mol%num(isp)
         call set_refq_gfn2(d4%q(:, isp), izp)
         call set_refalpha_gfn2(d4%aiw(:, :, isp), d4%ga, d4%gc, izp)
      end do
   end select

   allocate(d4%ngw(mref, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      call set_refgw(d4%ngw(:, isp), izp)
   end do

   allocate(d4%c6(mref, mref, mol%nid, mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, isp
         jzp = mol%num(jsp)
         do iref = 1, d4%ref(isp)
            do jref = 1, d4%ref(jsp)
               aiw(:) = d4%aiw(:, iref, isp) * d4%aiw(:, jref, jsp)
               c6 = thopi * d4%integrator%integrate(aiw)
               d4%c6(jref, iref, jsp, isp) = c6
               d4%c6(iref, jref, isp, jsp) = c6
            end do
         end do
      end do
   end do

end subroutine new_d4s_model


!> Update dispersion cache with precomputed coefficients and properties
subroutine update(self, mol, cache, cn, q, grad, only_c6)
   !DEC$ ATTRIBUTES DLLEXPORT :: update

   !> Instance of the dispersion model
   class(d4s_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dispersion cache to populate
   type(dispersion_cache), intent(inout) :: cache

   !> Coordination number of every atom
   real(wp), intent(in) :: cn(:)

   !> Partial charge of every atom
   real(wp), intent(in) :: q(:)

   !> Whether to compute derivatives
   logical, intent(in), optional :: grad

   !> Whether to compute only C6 coefficients
   logical, intent(in), optional :: only_c6

   logical :: do_grad, c6_only
   integer :: mref
   real(wp), allocatable :: gwvec(:, :, :), gwdcn(:, :, :), &
      & gwdq(:, :, :)

   mref = maxval(self%ref)
   do_grad = .false.
   if (present(grad)) do_grad = grad
   c6_only = .false.
   if (present(only_c6)) c6_only = only_c6

   if (.not. c6_only) then
      if (.not. allocated(cache%adiw)) allocate(cache%adiw(self%ngrid, mol%nat))
      if (.not. allocated(cache%aqiw)) allocate(cache%aqiw(self%ngrid, mol%nat))
      if (do_grad) then
         if (.not. allocated(cache%dadiwdcn)) &
            & allocate(cache%dadiwdcn(self%ngrid, mol%nat))
         if (.not. allocated(cache%dadiwdq)) &
            & allocate(cache%dadiwdq(self%ngrid, mol%nat))
         if (.not. allocated(cache%daqiwdcn)) &
            & allocate(cache%daqiwdcn(self%ngrid, mol%nat))
         if (.not. allocated(cache%daqiwdq)) &
            & allocate(cache%daqiwdq(self%ngrid, mol%nat))
      end if
   end if

   if (.not. allocated(cache%c6)) allocate(cache%c6(mol%nat, mol%nat))
   if (do_grad) then
      if (.not. allocated(cache%dc6dcn)) allocate(cache%dc6dcn(mol%nat, mol%nat))
      if (.not. allocated(cache%dc6dq)) allocate(cache%dc6dq(mol%nat, mol%nat))
   end if

   allocate(gwvec(mref, mol%nat, self%ncoup))
   if (do_grad) then
      allocate(gwdcn(mref, mol%nat, self%ncoup), gwdq(mref, mol%nat, self%ncoup))
      call weight_references(self, mol, cn, q, gwvec, gwdcn, gwdq)
   else
      call weight_references(self, mol, cn, q, gwvec)
   end if

   if (.not. c6_only) then
      if (do_grad) then
         call get_atomic_pol(self, mol, gwvec, gwdcn, gwdq, cache%adiw, cache%aqiw, &
            & cache%dadiwdcn, cache%dadiwdq, cache%daqiwdcn, cache%daqiwdq)
      else
         call get_atomic_pol(self, mol, gwvec, adiw=cache%adiw, aqiw=cache%aqiw)
      end if
   end if

   if (do_grad) then
      call get_atomic_c6(self, mol, gwvec, gwdcn, gwdq, &
         & cache%c6, cache%dc6dcn, cache%dc6dq)
   else
      call get_atomic_c6(self, mol, gwvec, c6=cache%c6)
   end if

end subroutine update


!> Calculate the weights of the reference system and the derivatives w.r.t.
!> coordination number for later use.
subroutine weight_references(self, mol, cn, q, gwvec, gwdcn, gwdq)

   !> Instance of the dispersion model
   type(d4s_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Coordination number of every atom
   real(wp), intent(in) :: cn(:)

   !> Partial charge of every atom
   real(wp), intent(in) :: q(:)

   !> Pairwise weighting for the atomic reference systems
   real(wp), intent(out) :: gwvec(:, :, :)

   !> derivative of the pairwise weighting function w.r.t. the coordination number
   real(wp), intent(out), optional :: gwdcn(:, :, :)

   !> derivative of the pairwise weighting function w.r.t. the charge scaling
   real(wp), intent(out), optional :: gwdq(:, :, :)

   integer :: iat, izp, iref, igw, jat, jzp
   real(wp) :: norm, dnorm, gw, expw, expd, gwk, dgwk, wf, zi, gi, maxcn
   real(wp), parameter :: eps_norm = tiny(1._wp)**0.5_wp

   if (present(gwdcn) .and. present(gwdq)) then
      gwvec(:, :, :) = 0.0_wp
      gwdcn(:, :, :) = 0.0_wp
      gwdq(:, :, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(gwvec, gwdcn, gwdq, mol, self, cn, q) &
      !$omp private(iat, izp, iref, igw, zi, gi, jat, jzp) &
      !$omp private(norm, dnorm, gw, expw, expd, gwk, dgwk, wf, maxcn)      
      do iat = 1, mol%nat
         izp = mol%id(iat)
         zi = self%zeff(izp)
         gi = self%eta(izp) * self%gc

         do jat = 1, mol%nat
            jzp = mol%id(jat)

            norm = 0.0_wp
            dnorm = 0.0_wp
            do iref = 1, self%ref(izp)
               do igw = 1, self%ngw(iref, izp)
                  wf = igw * self%wf(izp, jzp)
                  gw = weight_cn(wf, cn(iat), self%cn(iref, izp))
                  norm = norm + gw
                  dnorm = dnorm + 2*wf * (self%cn(iref, izp) - cn(iat)) * gw
               end do
            end do

            if (abs(norm) > eps_norm) then
               norm = 1.0_wp / norm
            else
               norm = 0.0_wp
            end if

            do iref = 1, self%ref(izp)
               expw = 0.0_wp
               expd = 0.0_wp
               do igw = 1, self%ngw(iref, izp)
                  wf = igw * self%wf(izp, jzp)
                  gw = weight_cn(wf, cn(iat), self%cn(iref, izp))
                  expw = expw + gw
                  expd = expd + 2*wf * (self%cn(iref, izp) - cn(iat)) * gw
               end do

               gwk = expw * norm
               if (is_exceptional(gwk) .or. norm == 0.0_wp) then
                  maxcn = maxval(self%cn(:self%ref(izp), izp))
                  if (abs(maxcn - self%cn(iref, izp)) < 1e-12_wp) then
                     gwk = 1.0_wp
                  else
                     gwk = 0.0_wp
                  end if
               end if

               gwvec(iref, iat, jat) = gwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
               gwdq(iref, iat, jat) = gwk * dzeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
               
               ! This expression behaves differently for -O0 and -O3 optimization
               ! levels for tiny norm values. It yields incorrect values for -O0 
               ! due to fp underflow. These tiny values can occur if the CN is very
               ! far from all reference CNs. To ensure consistent behavior, we set 
               ! the norm to zero above.
               dgwk = norm * (expd - expw * dnorm * norm)
               if (is_exceptional(dgwk) .or. norm == 0.0_wp) then
                  dgwk = 0.0_wp
               end if
               gwdcn(iref, iat, jat) = dgwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
            end do

         end do
      end do

   else

      gwvec(:, :, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(gwvec, mol, self, cn, q) &
      !$omp private(iat, izp, iref, igw, zi, gi, jat, jzp) &
      !$omp private(norm, gw, expw, gwk, wf, maxcn)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         zi = self%zeff(izp)
         gi = self%eta(izp) * self%gc

         do jat = 1, mol%nat
            jzp = mol%id(jat)

            norm = 0.0_wp
            do iref = 1, self%ref(izp)
               do igw = 1, self%ngw(iref, izp)
                  wf = igw * self%wf(izp, jzp)
                  norm = norm + weight_cn(wf, cn(iat), self%cn(iref, izp))
               end do
            end do

            if (abs(norm) > eps_norm) then
               norm = 1.0_wp / norm
            else
               norm = 0.0_wp
            end if

            do iref = 1, self%ref(izp)
               expw = 0.0_wp
               do igw = 1, self%ngw(iref, izp)
                  wf = igw * self%wf(izp, jzp)
                  expw = expw + weight_cn(wf, cn(iat), self%cn(iref, izp))
               end do

               gwk = expw * norm
               if (is_exceptional(gwk) .or. norm == 0.0_wp) then
                  maxcn = maxval(self%cn(:self%ref(izp), izp))
                  if (abs(maxcn - self%cn(iref, izp)) < 1e-12_wp) then
                     gwk = 1.0_wp
                  else
                     gwk = 0.0_wp
                  end if
               end if

               gwvec(iref, iat, jat) = gwk * zeta(self%ga, gi, self%q(iref, izp)+zi, q(iat)+zi)
            end do
         end do
      end do
   end if

end subroutine weight_references


!> Calculate atomic dispersion coefficients and their derivatives w.r.t.
!> the coordination numbers and atomic partial charges.
subroutine get_atomic_c6(self, mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   !> Instance of the dispersion model
   type(d4s_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Pairwise weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :, :)

   !> Derivative of the pairwise weighting function w.r.t. the coordination number
   real(wp), intent(in), optional :: gwdcn(:, :, :)

   !> Derivative of the pairwise weighting function w.r.t. the partial charge
   real(wp), intent(in), optional :: gwdq(:, :, :)

   !> C6 coefficients for all atom pairs.
   real(wp), intent(out) :: c6(:, :)

   !> Derivative of the C6 w.r.t. the coordination number
   real(wp), intent(out), optional :: dc6dcn(:, :)

   !> Derivative of the C6 w.r.t. the partial charge
   real(wp), intent(out), optional :: dc6dq(:, :)

   integer :: iat, jat, izp, jzp, iref, jref
   real(wp) :: refc6, dc6, dc6dcni, dc6dcnj, dc6dqi, dc6dqj

   if (present(gwdcn).and.present(dc6dcn) &
      & .and.present(gwdq).and.present(dc6dq)) then
      c6(:, :) = 0.0_wp
      dc6dcn(:, :) = 0.0_wp
      dc6dq(:, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(c6, dc6dcn, dc6dq, mol, self, gwvec, gwdcn, gwdq) &
      !$omp private(iat, jat, izp, jzp, iref, jref, refc6, dc6, dc6dqi, dc6dqj, &
      !$omp& dc6dcni, dc6dcnj)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            dc6 = 0.0_wp
            dc6dcni = 0.0_wp
            dc6dcnj = 0.0_wp
            dc6dqi = 0.0_wp
            dc6dqj = 0.0_wp
            do iref = 1, self%ref(izp)
               do jref = 1, self%ref(jzp)
                  refc6 = self%c6(iref, jref, izp, jzp)
                  dc6 = dc6 + gwvec(iref, iat, jat) * gwvec(jref, jat, iat) * refc6
                  dc6dcni = dc6dcni + gwdcn(iref, iat, jat) * gwvec(jref, jat, iat) * refc6
                  dc6dcnj = dc6dcnj + gwvec(iref, iat, jat) * gwdcn(jref, jat, iat) * refc6
                  dc6dqi = dc6dqi + gwdq(iref, iat, jat) * gwvec(jref, jat, iat) * refc6
                  dc6dqj = dc6dqj + gwvec(iref, iat, jat) * gwdq(jref, jat, iat) * refc6
               end do
            end do
            c6(iat, jat) = dc6
            c6(jat, iat) = dc6
            dc6dcn(iat, jat) = dc6dcni
            dc6dcn(jat, iat) = dc6dcnj
            dc6dq(iat, jat) = dc6dqi
            dc6dq(jat, iat) = dc6dqj
         end do
      end do

   else

      c6(:, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
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
                  dc6 = dc6 + gwvec(iref, iat, jat) * gwvec(jref, jat, iat) * refc6
               end do
            end do
            c6(iat, jat) = dc6
            c6(jat, iat) = dc6
         end do
      end do
   end if

end subroutine get_atomic_c6


!> Compute atomic polarizabilities from weighted references
subroutine get_atomic_pol(self, mol, gwvec, gwdcn, gwdq, adiw, aqiw, &
      & dadiwdcn, dadiwdq, daqiwdcn, daqiwdq)

   !> Instance of the dispersion model
   type(d4s_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Weighting function for the atomic reference systems
   real(wp), intent(in) :: gwvec(:, :, :)

   !> Derivative of the weighting function w.r.t. the coordination number
   real(wp), intent(in), optional :: gwdcn(:, :, :)

   !> Derivative of the weighting function w.r.t. the partial charge
   real(wp), intent(in), optional :: gwdq(:, :, :)

   !> Dipole-dipole dynamic polarizabilities
   real(wp), intent(out), optional :: adiw(:, :)

   !> Quadrupole-quadrupole dynamic polarizabilities
   real(wp), intent(out), optional :: aqiw(:, :)

   !> Derivative of DD polarizabilities w.r.t. coordination number
   real(wp), intent(out), optional :: dadiwdcn(:, :)

   !> Derivative of DD polarizabilities w.r.t. partial charge
   real(wp), intent(out), optional :: dadiwdq(:, :)

   !> Derivative of QQ polarizabilities w.r.t. coordination number
   real(wp), intent(out), optional :: daqiwdcn(:, :)

   !> Derivative of QQ polarizabilities w.r.t. partial charge
   real(wp), intent(out), optional :: daqiwdq(:, :)

   integer :: iat, izp, iref

   if (present(gwdcn).and.present(dadiwdcn).and.present(daqiwdcn) &
      & .and.present(gwdq).and.present(dadiwdq).and.present(daqiwdq)) then
      adiw(:, :) = 0.0_wp
      aqiw(:, :) = 0.0_wp
      dadiwdcn(:, :) = 0.0_wp
      dadiwdq(:, :) = 0.0_wp
      daqiwdcn(:, :) = 0.0_wp
      daqiwdq(:, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) shared(self, adiw, aqiw) &
      !$omp shared(mol, gwvec, gwdcn, gwdq, dadiwdcn, dadiwdq, daqiwdcn, daqiwdq) &
      !$omp private(iat, izp, iref)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do iref = 1, self%ref(izp)
            adiw(:, iat) = adiw(:, iat) + gwvec(iref, iat, iat) * self%aiw(:, iref, izp)
            dadiwdcn(:, iat) = dadiwdcn(:, iat) + gwdcn(iref, iat, iat) * self%aiw(:, iref, izp)
            dadiwdq(:, iat) = dadiwdq(:, iat) + gwdq(iref, iat, iat) * self%aiw(:, iref, izp)
         end do
         aqiw(:, iat) = adiw(:, iat) * self%r4r2(izp)
         daqiwdcn(:, iat) = dadiwdcn(:, iat) * self%r4r2(izp)
         daqiwdq(:, iat) = dadiwdq(:, iat) * self%r4r2(izp)
      end do
   else

      adiw(:, :) = 0.0_wp
      aqiw(:, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(adiw, aqiw, mol, self, gwvec) &
      !$omp private(iat, izp, iref)

      do iat = 1, mol%nat
         izp = mol%id(iat)
         do iref = 1, self%ref(izp)
            adiw(:, iat) = adiw(:, iat) + gwvec(iref, iat, iat) * self%aiw(:, iref, izp)
         end do
         aqiw(:, iat) = adiw(:, iat) * self%r4r2(izp)
      end do
   end if

end subroutine get_atomic_pol


!> Extract static polarizabilities from cache
subroutine get_polarizabilities(self, cache, alpha, alphaqq, &
   & dadcn, dadq, daqqdcn, daqqdq)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_polarizabilities

   !> Instance of the dispersion model
   class(d4s_model), intent(in) :: self

   !> Dispersion cache containing polarizabilities
   type(dispersion_cache), intent(in) :: cache

   !> Static dipole-dipole polarizabilities for all atoms
   real(wp), intent(out) :: alpha(:)

   !> Static quadrupole-quadrupole polarizabilities for all atoms
   real(wp), intent(out) :: alphaqq(:)

   !> Derivative of dipole polarizibility w.r.t. coordination number
   real(wp), intent(out), optional :: dadcn(:)

   !> Derivative of dipole polarizibility w.r.t. partial charge
   real(wp), intent(out), optional :: dadq(:)

   !> Derivative of quadrupole polarizibility w.r.t. coordination number
   real(wp), intent(out), optional :: daqqdcn(:)

   !> Derivative of quadrupole polarizibility w.r.t. partial charge
   real(wp), intent(out), optional :: daqqdq(:)

   integer :: iat

   do iat = 1, size(alpha)
      alpha(iat) = cache%adiw(1, iat)
   end do

   do iat = 1, size(alphaqq)
      alphaqq(iat) = cache%aqiw(1, iat)
   end do

   if (present(dadcn)) then
      do iat = 1, size(dadcn)
         dadcn(iat) = cache%dadiwdcn(1, iat)
      end do
   end if

   if (present(dadq)) then
      do iat = 1, size(dadq)
         dadq(iat) = cache%dadiwdq(1, iat)
      end do
   end if

   if (present(daqqdcn)) then
      do iat = 1, size(daqqdcn)
         daqqdcn(iat) = cache%daqiwdcn(1, iat)
      end do
   end if

   if (present(daqqdq)) then
      do iat = 1, size(daqqdq)
         daqqdq(iat) = cache%daqiwdq(1, iat)
      end do
   end if

end subroutine get_polarizabilities


!> Get two-body dispersion coefficients for atom pair ij
subroutine get_2b_coeffs(self, cache, iat, jat, izp, jzp, c6, c8)
   !> Dispersion model
   class(d4s_model), intent(in) :: self
   !> Dispersion cache
   type(dispersion_cache), intent(in) :: cache
   !> First atom index
   integer, intent(in) :: iat
   !> Second atom index
   integer, intent(in) :: jat
   !> Atomic number of first atom
   integer, intent(in) :: izp
   !> Atomic number of second atom
   integer, intent(in) :: jzp
   !> C6 dispersion coefficient
   real(wp), intent(out) :: c6
   !> C8 dispersion coefficient
   real(wp), intent(out) :: c8
   
   c6 = cache%c6(iat, jat)
   c8 = 3.0_wp * c6 * self%r4r2(izp) * self%r4r2(jzp)

end subroutine get_2b_coeffs



!> Get two-body dispersion coefficients for atom pair ij with derivatives
subroutine get_2b_derivs(self, cache, iat, jat, izp, jzp, c6, c8, &
   & dc6dcni, dc6dqi, dc6dcnj, dc6dqj, dc8dcni, dc8dqi, dc8dcnj, dc8dqj)
   !> Dispersion model
   class(d4s_model), intent(in) :: self
   !> Dispersion cache
   type(dispersion_cache), intent(in) :: cache
   !> First atom index
   integer, intent(in) :: iat
   !> Second atom index
   integer, intent(in) :: jat
   !> Atomic number of first atom
   integer, intent(in) :: izp
   !> Atomic number of second atom
   integer, intent(in) :: jzp
   !> C6 dispersion coefficient
   real(wp), intent(out) :: c6
   !> C8 dispersion coefficient
   real(wp), intent(out) :: c8
   !> Derivative of C6 w.r.t the coordination number of atom i
   real(wp), intent(out) :: dc6dcni
   !> Derivative of C6 w.r.t the partial charge of atom i
   real(wp), intent(out) :: dc6dqi
   !> Derivative of C6 w.r.t the coordination number of atom i
   real(wp), intent(out) :: dc6dcnj
   !> Derivative of C6 w.r.t the partial charge of atom i
   real(wp), intent(out) :: dc6dqj
   !> Derivative of C8 w.r.t the coordination number of atom i
   real(wp), intent(out) :: dc8dcni
   !> Derivative of C8 w.r.t the partial charge of atom i
   real(wp), intent(out) :: dc8dqi
   !> Derivative of C8 w.r.t the coordination number of atom i
   real(wp), intent(out) :: dc8dcnj
   !> Derivative of C8 w.r.t the partial charge of atom i
   real(wp), intent(out) :: dc8dqj
   
   real(wp) :: scale
   
   scale = 3.0_wp * self%r4r2(izp) * self%r4r2(jzp)

   c6 = cache%c6(iat, jat)
   c8 = scale * c6
   
   dc6dcni = cache%dc6dcn(iat, jat)
   dc8dcni = scale * dc6dcni
   dc6dcnj = cache%dc6dcn(jat, iat)
   dc8dcnj = scale * dc6dcnj
   
   dc6dqi = cache%dc6dq(iat, jat)
   dc8dqi = scale * dc6dqi
   dc6dqj = cache%dc6dq(jat, iat)
   dc8dqj = scale * dc6dqj

end subroutine get_2b_derivs


!> Get three-body dispersion coefficient for atom triple ijk
subroutine get_3b_coeffs(self, cache, iat, jat, kat, c9)
   !> Dispersion model
   class(d4s_model), intent(in) :: self
   !> Dispersion cache
   type(dispersion_cache), intent(in) :: cache
   !> First atom index
   integer, intent(in) :: iat
   !> Second atom index
   integer, intent(in) :: jat
   !> Third atom index
   integer, intent(in) :: kat
   !> C9 dispersion coefficient
   real(wp), intent(out) :: c9

   c9 = -sqrt(abs(cache%c6(iat, jat) * cache%c6(iat, kat) * cache%c6(jat, kat)))

end subroutine get_3b_coeffs


!> Get three-body dispersion coefficient for atom triple ijk with derivatives
subroutine get_3b_derivs(self, cache, iat, jat, kat, c9, &
   & dc9dcni, dc9dqi, dc9dcnj, dc9dqj, dc9dcnk, dc9dqk)
   !> Dispersion model
   class(d4s_model), intent(in) :: self
   !> Dispersion cache
   type(dispersion_cache), intent(in) :: cache
   !> First atom index
   integer, intent(in) :: iat
   !> Second atom index
   integer, intent(in) :: jat
   !> Third atom index
   integer, intent(in) :: kat
   !> C9 dispersion coefficient
   real(wp), intent(out) :: c9
   !> Derivative of C9 w.r.t the coordination number of atom i
   real(wp), intent(out) :: dc9dcni
   !> Derivative of C9 w.r.t the partial charge of atom i
   real(wp), intent(out) :: dc9dqi
   !> Derivative of C9 w.r.t the coordination number of atom j
   real(wp), intent(out) :: dc9dcnj
   !> Derivative of C9 w.r.t the partial charge of atom j
   real(wp), intent(out) :: dc9dqj
   !> Derivative of C9 w.r.t the coordination number of atom k
   real(wp), intent(out) :: dc9dcnk
   !> Derivative of C9 w.r.t the partial charge of atom k
   real(wp), intent(out) :: dc9dqk

   real(wp) :: c6ij, c6ik, c6jk

   c6ij = cache%c6(iat, jat)
   c6ik = cache%c6(iat, kat)
   c6jk = cache%c6(jat, kat)

   c9 = -sqrt(abs(c6ij * c6ik * c6jk))

   c6ij = 0.5_wp * c9 / c6ij
   c6ik = 0.5_wp * c9 / c6ik
   c6jk = 0.5_wp * c9 / c6jk

   dc9dcni = cache%dc6dcn(iat, jat) * c6ij + cache%dc6dcn(iat, kat) * c6ik
   dc9dqi = cache%dc6dq(iat, jat) * c6ij + cache%dc6dq(iat, kat) * c6ik
      
   dc9dcnj = cache%dc6dcn(jat, iat) * c6ij + cache%dc6dcn(jat, kat) * c6jk
   dc9dqj = cache%dc6dq(jat, iat) * c6ij + cache%dc6dq(jat, kat) * c6jk
      
   dc9dcnk = cache%dc6dcn(kat, iat) * c6ik + cache%dc6dcn(kat, jat) * c6jk
   dc9dqk = cache%dc6dq(kat, iat) * c6ik + cache%dc6dq(kat, jat) * c6jk

end subroutine get_3b_derivs


!> Calculate damping radius for two-body interactions
pure function get_2b_rdamp(self, izp, jzp) result(rdamp)
   !> Dispersion model
   class(d4s_model), intent(in) :: self
   !> Atomic number of first atom
   integer, intent(in) :: izp
   !> Atomic number of second atom
   integer, intent(in) :: jzp
   !> Damping radius
   real(wp) :: rdamp

   rdamp = sqrt(3.0_wp * self%r4r2(izp) * self%r4r2(jzp))

end function get_2b_rdamp


!> Calculate damping radius for three-body interactions
pure function get_3b_rdamp(self, izp, jzp) result(rdamp)
   !> Dispersion model
   class(d4s_model), intent(in) :: self
   !> Atomic number of first atom
   integer, intent(in) :: izp
   !> Atomic number of second atom
   integer, intent(in) :: jzp
   !> Three-body damping radius
   real(wp) :: rdamp

   rdamp = sqrt(3.0_wp * self%r4r2(izp) * self%r4r2(jzp))

end function get_3b_rdamp

end module dftd4_model_d4s
