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

!> @file tblite/repulsion/gxtb.f90
!> Provides a screened Coulomb + Pauli repulsion interaction as used in g-xTB

!> Classical Coulomb + Pauli repulsion interaction as used in the g-xTB Hamiltonian
module tblite_repulsion_gxtb
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count
   use tblite_blas, only : symv
   use tblite_container, only : container_cache
   use tblite_cutoff, only : get_lattice_points
   use tblite_repulsion_type, only : repulsion_type
   use tblite_repulsion_cache, only : repulsion_cache
   use tblite_utils_average, only : average_type, new_average, average_id
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_repulsion_gxtb

   !> Container to evaluate classical repulsion interactions for the g-xTB Hamiltonian
   type, public, extends(repulsion_type) :: gxtb_repulsion
      !> Exponent for the damping of the repulsion
      real(wp), allocatable :: alpha(:)
      !> Effective nuclear charge
      real(wp), allocatable :: zeff(:)
      !> Linear CN dependence of damping exponent
      real(wp), allocatable :: kcn(:)
      !> Linear atomic charge dependence of effective nuclear charge
      real(wp), allocatable :: kq(:)
      !> Geometric mean of fitted atomic radii offset for the repulsion
      real(wp), allocatable :: roffset(:, :)
      !> Distance exponent for repulsion damping
      real(wp) :: kexp
      !> Arithmetic average of first-order repulsion expansion coefficient
      real(wp), allocatable :: k1(:, :)
      !> Arithmetic average of second-order repulsion expansion coefficient
      real(wp), allocatable :: k2(:, :)
      !> Third-order repulsion expansion coefficient
      real(wp) :: k3
      !> Fourth-order repulsion expansion coefficient
      real(wp) :: k4
      !> Short-range repulsion correction prefactor
      real(wp) :: kshort
      !> Short-range repulsion correction damping exponent scaling
      real(wp) :: kshort_alpha
      !> Short-range repulsion correction distance damping exponent 
      real(wp) :: kshort_exp
      !> Pairwise scaled van-der-Waals radii
      real(wp), allocatable :: rvdw(:, :)
      !> Coordination number for modifying the damping exponent
      class(ncoord_type), allocatable :: ncoord
   contains
      !> Update repulsion cache
      procedure :: update
      !> Evaluate selfconsistent energy of the interaction
      procedure :: get_energy
      !> Evaluate charge dependent potential shift from the interaction
      procedure :: get_potential
      !> Evaluate gradient contributions from the selfconsistent interaction
      procedure :: get_gradient
   end type gxtb_repulsion

   character(len=*), parameter :: label = "self-consistent screened Coulomb + Pauli repulsion"

   real(wp), parameter :: cn_reg = 1e-6_wp

contains


!> Create a new instance of the g-xTB repulsion container
subroutine new_repulsion_gxtb(self, mol, alpha, zeff, kcn, kq, roffset, &
   & kexp, k1, k2, k2_light, k3, k4, kshort, kshort_alpha, kshort_exp, &
   & rvdw, cn_rcov, cn_exp, error, cutoff)
   !> Instance of the repulsion container
   type(gxtb_repulsion), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Exponent for the damping of the repulsion
   real(wp), intent(in) :: alpha(:)
   !> Effective nuclear charge
   real(wp), intent(in) :: zeff(:)
   !> Linear CN dependence of effective nuclear charge
   real(wp), intent(in) :: kcn(:)
   !> Linear atomic charge dependence of damping exponent
   real(wp), intent(in) :: kq(:)
   !> Fitted atomic radii offset for the repulsion
   real(wp), intent(in) :: roffset(:)
   !> Distance exponent for repulsion damping
   real(wp), intent(in) :: kexp
   !> First-order repulsion expansion coefficients
   real(wp), intent(in) :: k1(:)
   !> Second-order repulsion expansion coefficient
   !> for all elements except H and He
   real(wp), intent(in) :: k2
   !> Second-order repulsion expansion coefficient
   !> for light elements (H or He)
   real(wp), intent(in) :: k2_light
   !> Third-order repulsion expansion coefficient
   real(wp), intent(in) :: k3
   !> Fourth-order repulsion expansion coefficient
   real(wp), intent(in) :: k4
   !> Short-range repulsion correction prefactor
   real(wp), intent(in) :: kshort
   !> Short-range repulsion correction damping exponent scaling
   real(wp), intent(in) :: kshort_alpha
   !> Short-range repulsion correction distance damping exponent 
   real(wp), intent(in) :: kshort_exp
   !> Pairwise scaled van-der-Waals radii
   real(wp), intent(in) :: rvdw(:, :)
   !> Fitted covalent atomic radii for repulsion CN
   real(wp), intent(in) :: cn_rcov(:)
   !> Exponent of the repulsion CN
   real(wp), intent(in) :: cn_exp
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Real-space cutoff
   real(wp), intent(in), optional :: cutoff

   integer :: isp, izp, jsp, jzp
   real(wp) :: k2i, k2j
   type(average_type), allocatable :: average

   self%label = label

   self%alpha = alpha
   self%zeff = zeff
   self%kcn = kcn
   self%kq = kq
   
   allocate(average)
   call new_average(average, average_id%geometric)
   allocate(self%roffset(mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         self%roffset(jsp, isp) = average%value(roffset(isp),roffset(jsp))
      end do
   end do

   self%kexp = kexp

   ! Expansion coefficients for the repulsion
   allocate(self%k1(mol%nid, mol%nid), self%k2(mol%nid, mol%nid))
   call new_average(average, average_id%arithmetic)
   do isp = 1, mol%nid
      izp = mol%num(isp)
      k2i = merge(k2, k2_light, izp > 2)
      do jsp = 1, mol%nid
         jzp = mol%num(jsp)
         k2j = merge(k2, k2_light, jzp > 2)
         self%k1(jsp, isp) = average%value(k1(isp), k1(jsp))
         self%k2(jsp, isp) = average%value(k2i, k2j)
      end do
   end do
   self%k3 = k3
   self%k4 = k4

   self%kshort = kshort
   self%kshort_alpha = kshort_alpha
   self%kshort_exp = kshort_exp

   self%rvdw = rvdw

   if (present(cutoff)) self%cutoff = cutoff

   ! Create the coordination number 
   call new_ncoord(self%ncoord, mol, cn_count_type=cn_count%erf, &
      & error=error, kcn=cn_exp, rcov=cn_rcov)

end subroutine new_repulsion_gxtb


!> Update repulsion cache
subroutine update(self, mol, cache)
   !> Instance of the classical repulsion
   class(gxtb_repulsion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache

   real(wp), allocatable :: trans(:, :)
   type(repulsion_cache), pointer :: ptr

   call taint(cache, ptr)

   ! Allocate temporary arrays for energy and potential evaluation
   if (.not. allocated(ptr%scaled_zeff)) then
      allocate(ptr%scaled_zeff(mol%nat))
   end if
   if (.not. allocated(ptr%vvec)) then
      allocate(ptr%vvec(mol%nat))
   end if

   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff, trans)

   ! Calculate the scaled damping exponents
   if (.not. allocated(ptr%scaled_alpha)) then
      allocate(ptr%scaled_alpha(mol%nat))
   end if
   if (.not. allocated(ptr%scaled_dalphadcn)) then
      allocate(ptr%scaled_dalphadcn(mol%nat))
   end if
   call self%ncoord%get_coordination_number(mol, trans, ptr%vvec)
   call get_scaled_alpha(mol, self%alpha, self%kcn, ptr%vvec, &
      & ptr%scaled_alpha, ptr%scaled_dalphadcn)

   ! Setup the repulsion matrix
   if (.not. allocated(ptr%rmat)) then
      allocate(ptr%rmat(mol%nat, mol%nat))
   end if
   ptr%rmat(:, :) = 0.0_wp
   call get_repulsion_matrix(mol, ptr%scaled_alpha, self%roffset, &
      & self%kexp, self%k1, self%k2, self%k3, self%k4, self%kshort, &
      & self%kshort_alpha, self%kshort_exp, self%rvdw, trans, &
      & self%cutoff, ptr%rmat)

end subroutine update


!> Evaluate the scaled damping exponent
subroutine get_scaled_alpha(mol, alpha, kcn, cn, scaled_alpha, &
   & scaled_dalphadcn)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Exponent for the damping of the repulsion
   real(wp), intent(in) :: alpha(:)
   !> Linear CN dependence of damping exponent
   real(wp), intent(in) :: kcn(:)
   !> Coordination number
   real(wp), intent(in) :: cn(:)
   !> Scaled damping exponent
   real(wp), intent(out) :: scaled_alpha(:)
   !> Derivative of the scaled damping exponent w.r.t. the coordination number
   real(wp), intent(out), optional :: scaled_dalphadcn(:)

   integer :: iat, izp
   real(wp) :: sqrtcni

   scaled_alpha(:) = 0.0_wp
   if (present(scaled_dalphadcn)) then
      scaled_dalphadcn(:) = 0.0_wp
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Coordination number scaling with regularization
         sqrtcni = sqrt(cn(iat) + cn_reg**2) - cn_reg
         scaled_alpha(iat) = alpha(izp) * (1.0_wp + kcn(izp) * sqrtcni)
         scaled_dalphadcn(iat) = alpha(izp) * kcn(izp) / (2.0_wp * (sqrtcni + cn_reg))
      end do
   else
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Coordination number scaling with regularization
         sqrtcni = sqrt(cn(iat) + cn_reg**2) - cn_reg
         scaled_alpha(iat) = alpha(izp) * (1.0_wp + kcn(izp) * sqrtcni)
      end do 
   end if

end subroutine get_scaled_alpha


!> Evaluate the scaled effective nuclear charges
subroutine get_scaled_zeff(mol, zeff, kq, qat, scaled_zeff)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Base effective nuclear charge
   real(wp), intent(in) :: zeff(:)
   !> Linear atomic charge dependence of the effective nuclear charge
   real(wp), intent(in) :: kq(:)
   !> Atomic charges
   real(wp), intent(in) :: qat(:)
   !> Scaled effective nuclear charges
   real(wp), intent(out) :: scaled_zeff(:)

   integer :: iat, izp
   real(wp) :: scale

   scaled_zeff(:) = 0.0_wp
   do iat = 1, mol%nat
      izp = mol%id(iat)
      scale = 1.0_wp - kq(izp) * qat(iat)
      scaled_zeff(iat) = zeff(izp) * scale
   end do 

end subroutine get_scaled_zeff


!> Evaluates the repulsion matrix
subroutine get_repulsion_matrix(mol, scaled_alpha, roffset, kexp, k1, k2, k3, k4, &
   & kshort, kshort_alpha, kshort_exp, rvdw, trans, cutoff, rmat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> CN-scaled exponent for the damping of the repulsion
   real(wp), intent(in) :: scaled_alpha(:)
   !> Geometric mean of fitted atomic radii offset for the repulsion
   real(wp), intent(in) :: roffset(:, :) 
   !> Distance exponent for repulsion damping
   real(wp), intent(in) :: kexp
   !> Arithmetic average of first-order repulsion expansion coefficient
   real(wp), intent(in) :: k1(:, :)
   !> Second-order repulsion expansion coefficient
   real(wp), intent(in) :: k2(:, :)
   !> Third-order repulsion expansion coefficient
   real(wp), intent(in) :: k3
   !> Fourth-order repulsion expansion coefficient
   real(wp), intent(in) :: k4
   !> Short-range repulsion correction prefactor
   real(wp), intent(in) :: kshort
   !> Short-range repulsion correction damping exponent scaling
   real(wp), intent(in) :: kshort_alpha
   !> Short-range repulsion correction distance damping exponent 
   real(wp), intent(in) :: kshort_exp
   !> Pairwise scaled van-der-Waals radii
   real(wp), intent(in) :: rvdw(:, :)
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Repulsion matrix
   real(wp), intent(inout) :: rmat(:, :)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r1, r2, inv_r1, rel_r1, rel_r2, rel_r3, rel_r4, vec(3), inv_poly
   real(wp) :: rdamp, exa, cutoff2, alphai, alphaj, gpt_alpha

   cutoff2 = cutoff**2

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(rmat, mol, trans, cutoff2, scaled_alpha, roffset, kexp) &
   !$omp shared(k1, k2, k3, k4, kshort, kshort_alpha, kshort_exp, rvdw) &
   !$omp private(iat, jat, izp, jzp, itr, alphai, alphaj, gpt_alpha, r1, r2) &
   !$omp private(vec, inv_r1, rel_r1, rel_r2, rel_r3, rel_r4, inv_poly, rdamp, exa)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      alphai = scaled_alpha(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         alphaj = scaled_alpha(jat)
         ! Gaussian product exponent
         gpt_alpha = alphai * alphaj / (alphai + alphaj)
         do itr = 1, size(trans, dim=2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = sum(vec**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)
            inv_r1 = 1.0_wp / r1

            ! Interatomic distances relativ to the pairwise van-der-Waals radii
            rel_r1 = rvdw(jzp, izp) * inv_r1
            rel_r2 = rel_r1 * rel_r1
            rel_r3 = rel_r2 * rel_r1
            rel_r4 = rel_r2 * rel_r2

            ! Inverse power series expansion for the Coulomb and Pauli repulsion
            inv_poly = 1.0_wp + k1(jzp, izp) * inv_r1 + k2(jzp, izp) * rel_r2 &
               & + k3 * rel_r3 + k4 * rel_r4

            ! Exponential damping with short-range correction
            rdamp = (r1 + roffset(jzp, izp))
            exa = exp(-gpt_alpha * rdamp**kexp) &
               & + kshort * exp(-kshort_alpha * gpt_alpha * rdamp**kshort_exp)

            !$omp atomic
            rmat(jat, iat) = rmat(jat, iat) + exa * inv_poly
            !$omp atomic
            rmat(iat, jat) = rmat(iat, jat) + exa * inv_poly
         end do
      end do
   end do

end subroutine get_repulsion_matrix


!> Evaluate selfconsistent repulsion interaction energy
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the repulsion container
   class(gxtb_repulsion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Repulsion energy
   real(wp), intent(inout) :: energies(:)

   integer :: iat, izp
   real(wp) :: scale
   type(repulsion_cache), pointer :: ptr

   call view(cache, ptr)

   ! Calculate current effective nuclear charges
   call get_scaled_zeff(mol, self%zeff, self%kq, wfn%qat(:, 1), ptr%scaled_zeff)

   call symv(ptr%rmat, ptr%scaled_zeff, ptr%vvec, alpha=0.5_wp)
   do iat = 1, mol%nat
      energies(iat) = energies(iat) + ptr%scaled_zeff(iat) * ptr%vvec(iat)
   end do

end subroutine get_energy


!> Evaluate charge dependent potential shift for the selfconsistent repulsion
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the repulsion container
   class(gxtb_repulsion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   integer :: iat, izp
   real(wp) :: scale
   type(repulsion_cache), pointer :: ptr

   call view(cache, ptr)

   ! Calculate current effective nuclear charges
   call get_scaled_zeff(mol, self%zeff, self%kq, wfn%qat(:, 1), ptr%scaled_zeff)

   call symv(ptr%rmat, ptr%scaled_zeff, ptr%vvec)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      pot%vat(iat, 1) = pot%vat(iat, 1) - self%zeff(izp) * self%kq(izp) &
         & * ptr%vvec(iat)
   end do
end subroutine get_potential


!> Evaluate gradient contributions from the selfconsistent repulsion
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the repulsion container
   class(gxtb_repulsion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the repulsion energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   real(wp), allocatable :: dEdcn(:), trans(:, :)
   type(repulsion_cache), pointer :: ptr

   call view(cache, ptr)

   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff, trans)

   ! Calculate current effective nuclear charges
   call get_scaled_zeff(mol, self%zeff, self%kq, wfn%qat(:, 1), ptr%scaled_zeff)

   allocate(dEdcn(mol%nat), source=0.0_wp)
   call get_repulsion_derivs(mol, ptr%scaled_alpha, ptr%scaled_dalphadcn, &
      & ptr%scaled_zeff, self%roffset, self%kexp, self%k1, self%k2, self%k3,&
      & self%k4, self%kshort, self%kshort_alpha, self%kshort_exp, self%rvdw, &
      & trans, self%cutoff, gradient, sigma, dEdcn)

   call self%ncoord%add_coordination_number_derivs(mol, trans, dEdcn, gradient, sigma)

end subroutine get_gradient


!> Evaluates the gradient of the repulsion energy
subroutine get_repulsion_derivs(mol, scaled_alpha, scaled_dalphadcn, scaled_zeff, &
   & roffset, kexp, k1, k2, k3, k4, kshort, kshort_alpha, kshort_exp, rvdw, &
   & trans, cutoff, gradient, sigma, dEdcn)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Scaled damping exponent
   real(wp), intent(in) :: scaled_alpha(:)
   !> Derivative of the scaled damping exponent w.r.t. the coordination number
   real(wp), intent(in) :: scaled_dalphadcn(:)
   !> Scaled effective nuclear charges
   real(wp), intent(in) :: scaled_zeff(:)
   !> Geometric mean of fitted atomic radii offset for the repulsion
   real(wp), intent(in) :: roffset(:, :)
   !> Distance exponent for repulsion damping
   real(wp), intent(in) :: kexp
   !> Arithmetic average of first-order repulsion expansion coefficient
   real(wp), intent(in) :: k1(:, :)
   !> Second-order repulsion expansion coefficient
   real(wp), intent(in) :: k2(:, :)
   !> Third-order repulsion expansion coefficient
   real(wp), intent(in) :: k3
   !> Fourth-order repulsion expansion coefficient
   real(wp), intent(in) :: k4
   !> Short-range repulsion correction prefactor
   real(wp), intent(in) :: kshort
   !> Short-range repulsion correction damping exponent scaling
   real(wp), intent(in) :: kshort_alpha
   !> Short-range repulsion correction distance damping exponent 
   real(wp), intent(in) :: kshort_exp
   !> Pairwise scaled van-der-Waals radii
   real(wp), intent(in) :: rvdw(:, :)
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Molecular gradient of the repulsion energy
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivatives of the repulsion energy
   real(wp), intent(inout) :: sigma(:, :)
   !> Repulsion energy derivative w.r.t. the coordination numbers
   real(wp), contiguous, intent(inout) :: dEdcn(:)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r1, r2, vec(3), inv_r1, rel_r1, rel_r2, rel_r3, rel_r4, cutoff2

   real(wp) :: alphai, alphaj, dalphaidcn, dalphajdcn
   real(wp) :: gpt_alpha, dgpti, dgptj, prod_zeff, zeffi, zeffj
   real(wp) :: inv_poly, dinv_poly, exa, exashort, dexadrdamp, dexadgpt, rdamp
   real(wp) :: rdamp_kexp_m1, rdamp_kexp, rdamp_kshort_exp_m1, rdamp_kshort_exp
   real(wp) :: dEdgpt, dG(3), dS(3, 3)

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable ::  gradient_local(:, :), sigma_local(:, :), dEdcn_local(:)

   cutoff2 = cutoff**2

   !$omp parallel default(none) &
   !$omp shared(gradient, sigma, dEdcn, mol, trans, cutoff2, scaled_alpha) &
   !$omp shared(scaled_dalphadcn, scaled_zeff, roffset, kexp) &
   !$omp shared(k1, k2, k3, k4, kshort, kshort_alpha, kshort_exp, rvdw) &
   !$omp private(iat, jat, izp, jzp, itr, zeffi, zeffj, alphai, alphaj) &
   !$omp private(dalphaidcn, dalphajdcn, gpt_alpha, dgpti, dgptj, prod_zeff) &
   !$omp private(vec, r1, r2, inv_r1, rel_r1, rel_r2, rel_r3, rel_r4) &
   !$omp private(inv_poly, dinv_poly, rdamp, rdamp_kexp_m1, rdamp_kexp) &
   !$omp private(rdamp_kshort_exp_m1, rdamp_kshort_exp, exa, exashort) &
   !$omp private(dexadrdamp, dexadgpt, dEdgpt, dG, dS) &
   !$omp private(gradient_local, sigma_local, dEdcn_local)
   allocate(gradient_local(size(gradient, 1), size(gradient, 2)), source=0.0_wp)
   allocate(sigma_local(size(sigma, 1), size(sigma, 2)), source=0.0_wp)
   allocate(dEdcn_local(size(dEdcn)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      zeffi = scaled_zeff(iat)
      alphai = scaled_alpha(iat)
      dalphaidcn = scaled_dalphadcn(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         zeffj = scaled_zeff(jat)
         alphaj = scaled_alpha(jat)
         dalphajdcn = scaled_dalphadcn(jat)
         ! Gaussian product type exponent and derivative
         gpt_alpha = alphai * alphaj / (alphai + alphaj)
         dgpti = (alphaj**2 / (alphai + alphaj)**2)
         dgptj = (alphai**2 / (alphai + alphaj)**2)
         ! Product of the effective nuclear charges
         prod_zeff = zeffi * zeffj
         do itr = 1, size(trans, dim=2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = sum(vec**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)
            inv_r1 = 1.0_wp / r1

            ! Interatomic distances relativ to the pairwise van-der-Waals radii
            rel_r1 = rvdw(jzp, izp) * inv_r1
            rel_r2 = rel_r1 * rel_r1
            rel_r3 = rel_r2 * rel_r1
            rel_r4 = rel_r2 * rel_r2

            ! Inverse power series expansion for the Coulomb and Pauli repulsion
            inv_poly = 1.0_wp + k1(jzp, izp) * inv_r1 + k2(jzp, izp) * rel_r2 &
               & + k3 * rel_r3 + k4 * rel_r4
            dinv_poly = - k1(jzp, izp) * inv_r1 - 2.0_wp * k2(jzp, izp) * rel_r2 &
               & - 3.0_wp * k3 * rel_r3 - 4.0_wp * k4 * rel_r4
            dinv_poly = dinv_poly * inv_r1

            ! Damped interatomic distance for exponential damping
            rdamp = (r1 + roffset(jzp, izp))
            rdamp_kexp_m1 = rdamp**(kexp - 1.0_wp)
            rdamp_kexp = rdamp_kexp_m1 * rdamp
            rdamp_kshort_exp_m1 = rdamp**(kshort_exp - 1.0_wp)
            rdamp_kshort_exp = rdamp_kshort_exp_m1 * rdamp

            ! Exponential damping and short-range correction
            exa = exp(-gpt_alpha * rdamp_kexp) 
            exashort = kshort * exp(-kshort_alpha * gpt_alpha * rdamp_kshort_exp)

            ! Derivative of exponential damping and short-range correction 
            ! w.r.t. the damped interatomic distance
            dexadrdamp = -gpt_alpha * (kexp * rdamp_kexp_m1 * exa &
               & + kshort_alpha * kshort_exp * rdamp_kshort_exp_m1 * exashort)

            ! Derivative of exponential damping and short-range correction 
            ! w.r.t. the Gaussian product exponent
            dexadgpt = - rdamp_kexp * exa &
               & - kshort_alpha * rdamp_kshort_exp * exashort

            ! Direct derivative due to te expanded repulsion and exponential damping
            dG = (dinv_poly * (exa + exashort) + inv_poly * dexadrdamp) &
               & * prod_zeff * vec * inv_r1
            dS = spread(dG, 1, 3) * spread(vec, 2, 3)
            if (iat /= jat) then
               gradient_local(:, iat) = gradient_local(:, iat) + dG
               gradient_local(:, jat) = gradient_local(:, jat) - dG
               sigma_local(:, :) = sigma_local + dS
            else
               sigma_local(:, :) = sigma_local + 0.5_wp * dS
            end if

            ! Implicit derivative w.r.t. the exponent CN dependence
            dEdgpt = prod_zeff * inv_poly * dexadgpt
            dEdcn_local(iat) = dEdcn_local(iat) + dEdgpt * dgpti * dalphaidcn
            if (iat /= jat) then
               dEdcn_local(jat) = dEdcn_local(jat) + dEdgpt * dgptj * dalphajdcn
            end if

         end do
      end do
   end do
   !$omp end do
   !$omp critical (get_repulsion_derivs_)
   gradient(:, :) = gradient + gradient_local
   sigma(:, :) = sigma + sigma_local
   dEdcn(:) = dEdcn + dEdcn_local
   !$omp end critical (get_repulsion_derivs_)
   deallocate(gradient_local, sigma_local, dEdcn_local)
   !$omp end parallel

end subroutine get_repulsion_derivs


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(repulsion_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(repulsion_cache), allocatable :: tmp
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
   type(repulsion_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(repulsion_cache)
      ptr => target
   end select
end subroutine view

end module tblite_repulsion_gxtb