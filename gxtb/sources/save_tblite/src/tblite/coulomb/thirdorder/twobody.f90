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

!> @file tblite/coulomb/thirdorder.f90
!> Provides an two-body third-order tight-binding interaction

!> Isotropic third-order two-body correction
module tblite_coulomb_thirdorder_twobody
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use tblite_blas, only : gemv, symv
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_thirdorder_type, only : thirdorder_type
   use tblite_scf_potential, only : potential_type
   use tblite_utils_average, only : average_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_wignerseitz, only : wignerseitz_cell

   implicit none
   private

   public :: new_twobody_thirdorder

   !> Two-body correction for third-order charge expansion
   type, public, extends(thirdorder_type) :: twobody_thirdorder
      !> Averager for the hubbard parameters
      type(average_type), allocatable :: average
      !> Hubbard parameter for each shell and species
      real(wp), allocatable :: hubbard(:, :)
      !> CN-dependence of Hubbard parameter for each species
      real(wp), allocatable :: hubbard_cn(:)
      !> Exponent of radius dependent hubbard scaling
      real(wp) :: hubbard_exp
      !> Exponent for two-body third-order interaction
      real(wp) :: texp
      !> Onsite scaling factor for two-body interaction
      real(wp) :: onsite_scale
      !> Offsite scaling factor for two-body interaction
      real(wp) :: offsite_scale
      !> Long-range cutoff
      real(wp) :: rcut
   contains
      !> Update container cache
      procedure :: update
      !> Evaluate selfconsistent energy of the interaction
      procedure :: get_energy
      !> Evaluate charge dependent potential shift from the interaction
      procedure :: get_potential
      !> Evaluate gradient of charge dependent potential shift from the interaction
      procedure :: get_potential_gradient
      !> Evaluate gradient contributions from the selfconsistent interaction
      procedure :: get_gradient
      !> Evaluate tau matrix
      procedure :: get_tau_matrix
      !> Evaluate derivatives of tau matrix
      procedure :: get_tau_derivs
   end type twobody_thirdorder

   character(len=*), parameter :: label = "two-body third-order tight-binding"
   real(wp), parameter :: onethird = 1.0_wp / 3.0_wp
   real(wp), parameter :: cn_reg = 1e-6_wp

contains


!> Create new two-body third-order contribution
subroutine new_twobody_thirdorder(self, mol, texp, onsite_scale, offsite_scale, &
   & hubbard, average, hubbard_derivs, hubbard_cn, nshell)
   !> Instance of the third-oder tight-binding container
   type(twobody_thirdorder), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Exponent for two-body interaction
   real(wp), intent(in) :: texp
   !> Onsite scaling factor for two-body interaction
   real(wp), intent(in) :: onsite_scale
   !> Offsite scaling factor for two-body interaction
   real(wp), intent(in) :: offsite_scale
   !> Hubbard parameters
   real(wp), intent(in) :: hubbard(:, :)
   !> Averaging function for Hubbard parameter of a shell-pair
   type(average_type), intent(in) :: average
   !> Hubbard derivatives
   real(wp), intent(in) :: hubbard_derivs(:, :)
   !> CN-dependence of Hubbard parameters
   real(wp), intent(in), optional :: hubbard_cn(:)
   !> Number of shells for each species
   integer, intent(in), optional :: nshell(:)

   integer :: ind, iat

   self%label = label

   self%hubbard = hubbard
   self%hubbard_derivs = hubbard_derivs

   self%texp = texp
   self%onsite_scale = onsite_scale
   self%offsite_scale = offsite_scale

   self%shell_resolved = present(nshell)
   if (present(nshell)) then
      self%nsh_at = nshell(mol%id)
   else
      self%nsh_at = spread(1, 1, mol%nat)
   end if
   allocate(self%ish_at(mol%nat))
   ind = 0
   do iat = 1, mol%nat
      self%ish_at(iat) = ind
      ind = ind + self%nsh_at(iat)
   end do

   if (present(hubbard_cn)) then
      self%hubbard_cn = hubbard_cn
      self%cn_dep = .true.
   else
      allocate(self%hubbard_cn(mol%nid))
      self%hubbard_cn = 0.0_wp
      self%cn_dep = .false.
   end if

   self%average = average

end subroutine new_twobody_thirdorder


!> Update container cache
subroutine update(self, mol, cache)
   !> Instance of the third-order tight-binding container
   class(twobody_thirdorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   
   type(coulomb_cache), pointer :: ptr

   call taint(cache, ptr)

   if (.not.allocated(ptr%taumat)) then
      allocate(ptr%taumat(sum(self%nsh_at), sum(self%nsh_at)))
   end if
   call self%get_tau_matrix(mol, ptr, ptr%taumat)

   if (.not.allocated(ptr%vvec)) then
      allocate(ptr%vvec(sum(self%nsh_at)))
   end if

end subroutine update


!> Evaluate tau matrix
subroutine get_tau_matrix(self, mol, cache, taumat)
   !> Instance of the third-order tight-binding container
   class(twobody_thirdorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Tau matrix
   real(wp), contiguous, intent(out) :: taumat(:, :)

   taumat(:, :) = 0.0_wp

   if (any(mol%periodic)) then
      call get_taumat_3d(mol, self%nsh_at, self%ish_at, self%hubbard, &
         & self%hubbard_cn, self%average, cache%cn, self%texp, &
         & self%onsite_scale, self%offsite_scale, cache%wsc, taumat)
   else
      call get_taumat_0d(mol, self%nsh_at, self%ish_at, self%hubbard, &
         & self%hubbard_cn, self%average, cache%cn, self%texp, &
         & self%onsite_scale, self%offsite_scale, taumat)
   end if

end subroutine get_tau_matrix


!> Evaluate tau matrix for finite systems
subroutine get_taumat_0d(mol, nsh_at, ish_at, hubbard, hubbard_cn, &
   & average, cn, texp, onsite_scale, offsite_scale, taumat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nsh_at(:)
   !> Index offset for each shell
   integer, intent(in) :: ish_at(:)
   !> Hubbard parameter for each shell
   real(wp), intent(in) :: hubbard(:, :)
   !> CN-dependence of Hubbard parameters
   real(wp), intent(in) :: hubbard_cn(:)
   !> Averaging function for Hubbard parameter of a shell-pair
   type(average_type), intent(in) :: average
   !> Coordination number
   real(wp), intent(in) :: cn(:)
   !> Exponent for two-body interaction
   real(wp), intent(in) :: texp
   !> Onsite scaling factor for two-body interaction
   real(wp), intent(in) :: onsite_scale
   !> Offsite scaling factor for two-body interaction
   real(wp), intent(in) :: offsite_scale
   !> Coulomb matrix
   real(wp), intent(inout) :: taumat(:, :)

   integer :: iat, jat, izp, jzp, ii, jj, ish, jsh
   real(wp) :: sqrtcni, sqrtcnj, scalei, scalej
   real(wp) :: vec(3), r1, gam, hubi, hubj, poly, tmp 

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(taumat, mol, nsh_at, ish_at, hubbard, hubbard_cn, cn) &
   !$omp shared(average, texp, onsite_scale, offsite_scale) &
   !$omp private(iat, izp, ii, ish, jat, jzp, jj, jsh, sqrtcni, sqrtcnj) & 
   !$omp private(scalei, scalej, hubi, hubj, gam, vec, r1, poly, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = ish_at(iat)
      ! Coordination number scaling
      sqrtcni = sqrt(cn(iat) + cn_reg**2) - cn_reg
      scalei = 1.0_wp + hubbard_cn(izp) * sqrtcni
      do jat = 1, iat-1
         jzp = mol%id(jat)
         jj = ish_at(jat)
         ! Coordination number scaling
         sqrtcnj = sqrt(cn(jat) + cn_reg**2) - cn_reg
         scalej = 1.0_wp + hubbard_cn(jzp) * sqrtcnj
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         do ish = 1, nsh_at(iat)
            hubi = hubbard(ish, izp) * scalei
            do jsh = 1, nsh_at(jat)
               hubj = hubbard(jsh, jzp) * scalej
               gam = 1.0_wp / average%value(hubi, hubj)
               ! Polynomial expansion of short-range tau function
               poly = offsite_scale * gam * r1 &
                  & * (1.0_wp - 0.5_wp * texp * gam * r1)
               tmp = poly * exp(-texp * gam * r1)
               !$omp atomic
               taumat(jj+jsh, ii+ish) = taumat(jj+jsh, ii+ish) + tmp
               !$omp atomic
               taumat(ii+ish, jj+jsh) = taumat(ii+ish, jj+jsh) + tmp
            end do
         end do
      end do
      do ish = 1, nsh_at(iat)
         hubi = hubbard(ish, izp) * scalei
         do jsh = 1, ish-1
            hubj = hubbard(jsh, izp) * scalei
            gam = 1.0_wp/average%value(hubi, hubj)
            tmp = - onsite_scale / gam**2.0_wp
            !$omp atomic
            taumat(ii+jsh, ii+ish) = taumat(ii+jsh, ii+ish) + tmp
            !$omp atomic
            taumat(ii+ish, ii+jsh) = taumat(ii+ish, ii+jsh) + tmp
         end do
         hubi = hubbard(ish, izp) * scalei
         tmp = - onsite_scale * hubi**2.0_wp
         !$omp atomic
         taumat(ii+ish, ii+ish) = taumat(ii+ish, ii+ish) + tmp
      end do
   end do

end subroutine get_taumat_0d


!> Evaluate tau matrix for periodic systems
subroutine get_taumat_3d(mol, nsh_at, ish_at, hubbard, hubbard_cn, &
   & average, cn, texp, onsite_scale, offsite_scale, wsc, taumat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nsh_at(:)
   !> Index offset for each shell
   integer, intent(in) :: ish_at(:)
   !> Hubbard parameter for each shell
   real(wp), intent(in) :: hubbard(:, :)
   !> CN-dependence of Hubbard parameters
   real(wp), intent(in) :: hubbard_cn(:)
   !> Averaging function for Hubbard parameter of an shell-pair
   type(average_type), intent(in) :: average
   !> Coordination number
   real(wp), intent(in) :: cn(:)
   !> Exponent for two-body interaction
   real(wp), intent(in) :: texp
   !> Onsite scaling factor for two-body interaction
   real(wp), intent(in) :: onsite_scale
   !> Offsite scaling factor for two-body interaction
   real(wp), intent(in) :: offsite_scale
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Coulomb matrix
   real(wp), intent(inout) :: taumat(:, :)

   taumat(:, :) = 0.0_wp

end subroutine get_taumat_3d


!> Evaluate selfconsistent energy of the interaction
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the third-order tight-binding container
   class(twobody_thirdorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic and tight-binding energy
   real(wp), intent(inout) :: energies(:)

   integer :: iat, izp, ii, ish, jat, jzp, jj, jsh
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   if (self%shell_resolved) then
      call symv(ptr%taumat, wfn%qsh(:, 1), ptr%vvec, alpha=onethird)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = self%ish_at(iat)
         do ish = 1, self%nsh_at(iat)
            energies(iat) = energies(iat) + ptr%vvec(ii+ish) &
               & * wfn%qsh(ii+ish, 1)**2 * self%hubbard_derivs(ish, izp)
         end do
      end do
   else
      call symv(ptr%taumat, wfn%qat(:, 1), ptr%vvec, alpha=onethird)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         energies(iat) = energies(iat) + ptr%vvec(iat) &
            & * wfn%qat(iat, 1)**2 * self%hubbard_derivs(1, izp)
      end do
   end if
end subroutine get_energy


!> Evaluate charge dependent potential shift from the interaction
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the third-order tight-binding container
   class(twobody_thirdorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   integer :: iat, izp, ii, ish
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   if (self%shell_resolved) then
      call symv(ptr%taumat, wfn%qsh(:, 1), ptr%vvec, alpha=onethird)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = self%ish_at(iat)
         do ish = 1, self%nsh_at(iat)
            pot%vsh(ii+ish, 1) = pot%vsh(ii+ish, 1) + 2.0_wp * ptr%vvec(ii+ish) &
               & * wfn%qsh(ii+ish, 1) * self%hubbard_derivs(ish, izp)
            ! Overwrite temp. vector intermediate square charge term
            ptr%vvec(ii+ish) = wfn%qsh(ii+ish, 1)**2 &
               & * self%hubbard_derivs(ish, izp)
         end do
      end do
      call symv(ptr%taumat, ptr%vvec, pot%vsh(:, 1), beta=1.0_wp, alpha=onethird)
   else
      call symv(ptr%taumat, wfn%qat(:, 1), ptr%vvec, alpha=onethird)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         pot%vat(iat, 1) = pot%vat(iat, 1) + 2.0_wp * ptr%vvec(iat) &
            & * wfn%qat(iat, 1) * self%hubbard_derivs(1, izp)
         ptr%vvec(iat) = wfn%qat(iat, 1)**2 * self%hubbard_derivs(1, izp)
      end do
      call symv(ptr%taumat, ptr%vvec, pot%vat(:, 1), beta=1.0_wp, alpha=onethird)
   end if
end subroutine get_potential


!> Evaluate derivatives of tau matrix
subroutine get_tau_derivs(self, mol, cache, qat, qsh, dtaudr, dtaudL, dtaudcn)
   !> Instance of the third-order tight-binding container
   class(twobody_thirdorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)
   !> Derivative of interactions with respect to cartesian displacements
   real(wp), contiguous, intent(out) :: dtaudr(:, :, :)
   !> Derivative of interactions with respect to strain deformations
   real(wp), contiguous, intent(out) :: dtaudL(:, :, :)
   !> Coulomb matrix derivative w.r.t. the coordination numbers
   real(wp), contiguous, intent(out) :: dtaudcn(:, :)

   if(self%shell_resolved) then
      if (any(mol%periodic)) then
         call get_dtaumat_3d(mol, self%nsh_at, self%ish_at, self%hubbard, & 
            & self%hubbard_cn, self%average, cache%cn, self%hubbard_derivs, &
            & self%texp, self%onsite_scale, self%offsite_scale, cache%wsc, &
            & qsh, dtaudr, dtaudL, dtaudcn)
      else
         call get_dtaumat_0d(mol, self%nsh_at, self%ish_at, self%hubbard, & 
            & self%hubbard_cn, self%average, cache%cn, self%hubbard_derivs, & 
            & self%texp, self%onsite_scale, self%offsite_scale, qsh, &
            & dtaudr, dtaudL, dtaudcn)
      end if
   else
      if (any(mol%periodic)) then
         call get_dtaumat_3d(mol, self%nsh_at, self%ish_at, self%hubbard, & 
            & self%hubbard_cn, self%average, cache%cn, self%hubbard_derivs, &
            & self%texp, self%onsite_scale, self%offsite_scale, cache%wsc, &
            & qat, dtaudr, dtaudL, dtaudcn)
      else
         call get_dtaumat_0d(mol, self%nsh_at, self%ish_at, self%hubbard, & 
            & self%hubbard_cn, self%average, cache%cn, self%hubbard_derivs, &
            & self%texp, self%onsite_scale, self%offsite_scale, qat, &
            & dtaudr, dtaudL, dtaudcn)
      end if
   end if 

end subroutine get_tau_derivs


!> Evaluate uncontracted derivatives of tau matrix for finite system
subroutine get_dtaumat_0d(mol, nshell, offset, hubbard, hubbard_cn, &
   & average, cn, hubbard_derivs, texp, onsite_scale, offsite_scale, &
   & qvec, dtaudr, dtaudL, dtaudcn)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Index offset for each shell
   integer, intent(in) :: offset(:)
   !> Hubbard parameter for each shell and species
   real(wp), intent(in) :: hubbard(:, :)
   !> CN-dependence of Hubbard parameters
   real(wp), intent(in) :: hubbard_cn(:)
   !> Averaging function for Hubbard parameter of a shell-pair
   type(average_type), intent(in) :: average
   !> Coordination number
   real(wp), intent(in) :: cn(:)
   !> Hubbard parameter derivatives
   real(wp), intent(in) :: hubbard_derivs(:, :)
   !> Exponent for two-body interaction
   real(wp), intent(in) :: texp
   !> Onsite scaling factor for two-body interaction
   real(wp), intent(in) :: onsite_scale
   !> Offsite scaling factor for two-body interaction
   real(wp), intent(in) :: offsite_scale
   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)
   !> Derivative of interactions with respect to cartesian displacements
   real(wp), intent(out) :: dtaudr(:, :, :)
   !> Derivative of interactions with respect to strain deformations
   real(wp), intent(out) :: dtaudL(:, :, :)
   !> Coulomb matrix derivative w.r.t. the coordination numbers
   real(wp), intent(out) :: dtaudcn(:, :)

   integer :: iat, jat, izp, jzp, ii, jj, ish, jsh
   real(wp) :: vec(3), r1, r2, damp, hubi, hubj, dhubidcn, dhubjdcn
   real(wp) :: sqrtcni, sqrtcnj, scalei, scalej, dscaleidcn, dscalejdcn 
   real(wp) :: gam, gam2, dgami, dgamj, dpoly, dtmp, dtmpdgam, dG(3), dS(3, 3)
   real(wp), allocatable :: didr(:, :, :), didL(:, :, :), didcn(:, :)

   dtaudr(:, :, :) = 0.0_wp
   dtaudL(:, :, :) = 0.0_wp
   dtaudcn(:, :) = 0.0_wp

   !$omp parallel default(none) shared(dtaudr, dtaudL, dtaudcn)&
   !$omp shared(mol, qvec, nshell, offset, hubbard, hubbard_cn, average) &
   !$omp shared(cn, hubbard_derivs, texp, onsite_scale, offsite_scale) &
   !$omp private(iat, izp, ii, ish, jat, jzp, jj, jsh, r1, r2, vec, damp) &
   !$omp private(sqrtcni, sqrtcnj, scalei, scalej, dscaleidcn, dscalejdcn) &
   !$omp private(hubi, hubj, dhubidcn, dhubjdcn, gam, gam2, dgami, dgamj) &
   !$omp private(dpoly, dtmp, dtmpdgam, dG, dS, didr, didL, didcn)
   allocate(didr, source=dtaudr)
   allocate(didL, source=dtaudL)
   allocate(didcn, source=dtaudcn)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      ! Coordination number scaling
      sqrtcni = sqrt(cn(iat) + cn_reg**2) - cn_reg
      scalei = 1.0_wp + hubbard_cn(izp) * sqrtcni
      dscaleidcn = 0.5_wp * hubbard_cn(izp) / (sqrtcni + cn_reg)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         jj = offset(jat)
         ! Coordination number scaling
         sqrtcnj = sqrt(cn(jat) + cn_reg**2) - cn_reg
         scalej = 1.0_wp + hubbard_cn(jzp) * sqrtcnj
         dscalejdcn = 0.5_wp * hubbard_cn(jzp) / (sqrtcnj + cn_reg)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = dot_product(vec, vec)
         r1 = sqrt(r2)
         do ish = 1, nshell(iat)
            ! CN dependent hubbard parameter and its CN derivative
            hubi = hubbard(ish, izp) * scalei
            dhubidcn = hubbard(ish, izp) * dscaleidcn
            do jsh = 1, nshell(jat)
               ! CN dependent hubbard parameter and its CN derivative
               hubj = hubbard(jsh, jzp) * scalej
               dhubjdcn = hubbard(jsh, jzp) * dscalejdcn
               ! Inverse average hubbard parameter and its CN derivative
               gam = 1.0_wp / average%value(hubi, hubj)
               gam2 = gam * gam
               dgami = -gam2 * average%deriv(hubi, hubj) * dhubidcn
               dgamj = -gam2 * average%deriv(hubj, hubi) * dhubjdcn
               ! Derivative of short-range tau function
               dpoly = offsite_scale * (1.0_wp - 2.0_wp * texp * gam * r1 &
                  & + 0.5_wp * texp*texp * gam2 * r2)
               ! Direct position derivative of short-range tau function
               dtmp = dpoly * gam * exp(-texp * gam * r1)
               ! Implicit inverse averaged hubbard parameter derivative
               dtmpdgam = dpoly * r1 * exp(-texp * gam * r1)
               dG = dtmp*vec/r1
               dS = spread(dG, 1, 3) * spread(vec, 2, 3)
               ! Build once-contracted derivatives ∂τ · q,
               ! accounting for opposite changes in ij and ji elements
               didr(:, iat, jj+jsh) = +dG*qvec(ii+ish) + didr(:, iat, jj+jsh)
               didr(:, jat, ii+ish) = -dG*qvec(jj+jsh) + didr(:, jat, ii+ish)
               didr(:, jat, jj+jsh) = -dG*qvec(ii+ish) + didr(:, jat, jj+jsh)
               didr(:, iat, ii+ish) = +dG*qvec(jj+jsh) + didr(:, iat, ii+ish)
               didL(:, :, jj+jsh) = +dS*qvec(ii+ish) + didL(:, :, jj+jsh)
               didL(:, :, ii+ish) = +dS*qvec(jj+jsh) + didL(:, :, ii+ish)
               didcn(jat, ii+ish) = +dtmpdgam*dgamj*qvec(jj+jsh) + didcn(jat, ii+ish)
               didcn(iat, jj+jsh) = +dtmpdgam*dgami*qvec(ii+ish) + didcn(iat, jj+jsh)
               didcn(iat, ii+ish) = +dtmpdgam*dgami*qvec(jj+jsh) + didcn(iat, ii+ish)
               didcn(jat, jj+jsh) = +dtmpdgam*dgamj*qvec(ii+ish) + didcn(jat, jj+jsh)
            end do
         end do
      end do
      do ish = 1, nshell(iat)
         ! CN dependent hubbard parameter and its CN derivative
         hubi = hubbard(ish, izp) * scalei
         dhubidcn = hubbard(ish, izp) * dscaleidcn
         do jsh = 1, ish-1
            ! CN dependent hubbard parameter and its CN derivative
            hubj = hubbard(jsh, izp) * scalei
            dhubjdcn = hubbard(jsh, izp) * dscaleidcn
            ! Inverse average onsite hubbard parameter and its CN derivative
            gam = average%value(hubi, hubj)
            dgami = average%deriv(hubi, hubj) * dhubidcn
            dgamj = average%deriv(hubj, hubi) * dhubjdcn
            ! Implicit inverse averaged hubbard parameter derivative
            dtmpdgam = - onsite_scale * 2.0_wp * gam 
            didcn(iat, ii+ish) = +dtmpdgam*dgamj*qvec(ii+jsh) + didcn(iat, ii+ish)
            didcn(iat, ii+jsh) = +dtmpdgam*dgami*qvec(ii+ish) + didcn(iat, ii+jsh)
            didcn(iat, ii+ish) = +dtmpdgam*dgami*qvec(ii+jsh) + didcn(iat, ii+ish)
            didcn(iat, ii+jsh) = +dtmpdgam*dgamj*qvec(ii+ish) + didcn(iat, ii+jsh)
         end do
         dtmpdgam = - onsite_scale * 2.0_wp * hubi
         didcn(iat, ii+ish) = +dtmpdgam*dhubidcn*qvec(ii+ish) + didcn(iat, ii+ish)
      end do
   end do
   !$omp critical (get_dtaumat_0d_)
   dtaudr(:, :, :) = dtaudr + didr
   dtaudL(:, :, :) = dtaudL + didL
   dtaudcn(:, :) = dtaudcn + didcn
   !$omp end critical (get_dtaumat_0d_)
   deallocate(didL, didr, didcn)
   !$omp end parallel

end subroutine get_dtaumat_0d


!> Evaluate uncontracted derivatives of tau matrix for finite system
subroutine get_dtaumat_3d(mol, nshell, offset, hubbard, hubbard_cn, &
   & average, cn, hubbard_derivs, texp, onsite_scale, offsite_scale, &
   & wsc, qvec, dtaudr, dtaudL, dtaudcn)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Index offset for each shell
   integer, intent(in) :: offset(:)
   !> Hubbard parameter for each shell and species
   real(wp), intent(in) :: hubbard(:, :)
   !> CN-dependence of Hubbard parameters
   real(wp), intent(in) :: hubbard_cn(:)
   !> Averaging function for Hubbard parameter of a shell-pair
   type(average_type), intent(in) :: average
   !> Coordination number
   real(wp), intent(in) :: cn(:)
   !> Hubbard parameter derivatives
   real(wp), intent(in) :: hubbard_derivs(:, :)
   !> Exponent for two-body interaction
   real(wp), intent(in) :: texp
   !> Onsite scaling factor for two-body interaction
   real(wp), intent(in) :: onsite_scale
   !> Offsite scaling factor for two-body interaction
   real(wp), intent(in) :: offsite_scale
   !> Wigner-Seitz image information
   type(wignerseitz_cell), intent(in) :: wsc
   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)
   !> Derivative of interactions with respect to cartesian displacements
   real(wp), intent(out) :: dtaudr(:, :, :)
   !> Derivative of interactions with respect to strain deformations
   real(wp), intent(out) :: dtaudL(:, :, :)
   !> Coulomb matrix derivative w.r.t. the coordination numbers
   real(wp), intent(out) :: dtaudcn(:, :)

   dtaudr(:, :, :) = 0.0_wp
   dtaudL(:, :, :) = 0.0_wp
   dtaudcn(:, :) = 0.0_wp

end subroutine get_dtaumat_3d


!> Evaluate gradient of charge dependent potential shift from the interaction
subroutine get_potential_gradient(self, mol, cache, wfn, pot)
   !> Instance of the third-order tight-binding container
   class(twobody_thirdorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

end subroutine get_potential_gradient


!> Evaluate gradient contributions from the selfconsistent interaction
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the third-order tight-binding container
   class(twobody_thirdorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the third-order tight-binding energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the third-order tight-binding energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: ndim, iat, izp, ii, ish
   real(wp), allocatable :: dtaudr(:, :, :), dtaudL(:, :, :)
   real(wp), allocatable :: dtaudcn(:, :), dedcn(:)
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   ndim = sum(self%nsh_at)
   allocate(dtaudr(3, mol%nat, ndim), dtaudL(3, 3, ndim), dtaudcn(mol%nat, ndim))
   
   call self%get_tau_derivs(mol, ptr, wfn%qat(:, 1), wfn%qsh(:, 1), &
      & dtaudr, dtaudL, dtaudcn)

   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = self%ish_at(iat)
      do ish = 1, self%nsh_at(iat)
         ptr%vvec(ii+ish) = wfn%qsh(ii+ish, 1)**2 * self%hubbard_derivs(ish, izp)
      end do
   end do

   call gemv(dtaudr, ptr%vvec, gradient, beta=1.0_wp, alpha=onethird)
   call gemv(dtaudL, ptr%vvec, sigma, beta=1.0_wp, alpha=onethird)

   ! Derivatives due to CN dependent Hubbard parameters
   if (self%cn_dep) then
      allocate(dedcn(mol%nat))
      call gemv(dtaudcn, ptr%vvec, dedcn)
      call gemv(ptr%dcndr, dedcn, gradient, beta=1.0_wp, alpha=onethird)
      call gemv(ptr%dcndL, dedcn, sigma, beta=1.0_wp, alpha=onethird)
   end if

end subroutine get_gradient

!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the coulomb cache
   type(coulomb_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(coulomb_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

!> Return reference to container cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(coulomb_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(coulomb_cache)
      ptr => target
   end select
end subroutine view

end module tblite_coulomb_thirdorder_twobody
