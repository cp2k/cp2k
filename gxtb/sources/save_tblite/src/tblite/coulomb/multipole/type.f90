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

!> @file tblite/coulomb/multipole/type.f90
!> Provides an implemenation of a damped multipole based second-order electrostatic

!> Anisotropic second-order electrostatics using a damped multipole expansion
module tblite_coulomb_multipole_type
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use mctc_io_constants, only : pi
   use tblite_blas, only : dot, gemv, symv, gemm
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_ewald, only : get_dir_cutoff, get_rec_cutoff
   use tblite_coulomb_type, only : coulomb_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_wignerseitz, only : wignerseitz_cell
   implicit none
   private


   !> Container to handle multipole electrostatics
   type, public, extends(coulomb_type), abstract :: damped_multipole
      !> Damping prefactor for inverse quadratic interactions
      real(wp) :: kdmp3 = 0.0_wp
      !> Damping prefactor for inverse cubic interactions
      real(wp) :: kdmp5 = 0.0_wp
      !> Damping prefactor for inverse quartic interactions
      real(wp) :: kdmp7 = 0.0_wp
      !> Damping prefactor for inverse quintic interactions
      real(wp) :: kdmp9 = 0.0_wp
      !> Damping function exponent for inverse quadratic interactions
      real(wp) :: kexp3 = 0.0_wp
      !> Damping function exponent for inverse cubic interactions
      real(wp) :: kexp5 = 0.0_wp
      !> Damping function exponent for inverse quartic interactions
      real(wp) :: kexp7 = 0.0_wp 
      !> Damping function exponent for inverse quintic interactions
      real(wp) :: kexp9 = 0.0_wp
      !> Kernel for on-site dipole exchange-correlation
      real(wp), allocatable :: dkernel(:)
      !> Kernel for on-site quadrupolar exchange-correlation
      real(wp), allocatable :: qkernel(:)
      !> Scaling factor for atomic dipole moments
      real(wp), allocatable :: dip_scale(:)
   contains
      !> Return dependency on density
      procedure :: variable_info
      !> Get anisotropic electrostatic energy
      procedure :: get_energy
      !> Get anisotropic electrostatic potential
      procedure :: get_potential
      !> Get derivatives of anisotropic electrostatics
      procedure :: get_gradient
      ! These additional functions are necessary as the atomic contributions to AES are not 
      ! the same in this implementation and in the GFN2 paper,
      ! for details see https://github.com/tblite/tblite/pull/224/files#r1970341792
      !> Get only AXC part of the anisotropic energy
      procedure :: get_energy_axc
      !> Get AES energy of the anisotropic electrostatics
      procedure :: get_energy_aes
      !> Get interaction matrix for all multipole moments up to quadrupoles
      procedure :: get_multipole_matrix
      !> Evaluate multipole damping radius
      procedure(get_mrad_pair), deferred :: get_mrad_pair
      !> Evaluate multipole damping radius derivatives
      procedure(get_mrad_derivs), deferred :: get_mrad_derivs
      !> Evaluate pairwise damping factors
      procedure(get_damping_pair), deferred :: get_damping_pair
      !> Evaluate derivatives of damping factors
      procedure(get_damping_derivs), deferred :: get_damping_derivs
   end type damped_multipole


   abstract interface
      !> Evaluate multipole damping radius for a diatomic pair
      pure subroutine get_mrad_pair(self, cache, iat, jat, isp, jsp, mrad)
         import :: damped_multipole, coulomb_cache, wp
         !> Instance of the multipole containeris
         class(damped_multipole), intent(in) :: self
         !> Reusable data container
         type(coulomb_cache), intent(in) :: cache
         !> Index of atom i
         integer, intent(in) :: iat
         !> Index of atom j
         integer, intent(in) :: jat
         !> Species of atom i
         integer, intent(in) :: isp
         !> Species of atom j
         integer, intent(in) :: jsp
         !> Multipole damping radius for the pair
         real(wp), intent(out) :: mrad
      end subroutine get_mrad_pair

      !> Evaluate multipole damping radius derivatives for a diatomic pair
      pure subroutine get_mrad_derivs(self, cache, iat, jat, isp, jsp, mrad, &
         & dmraddcni, dmraddcnj)
         import :: damped_multipole, coulomb_cache, wp
         !> Instance of the multipole containeris
         class(damped_multipole), intent(in) :: self
         !> Reusable data container
         type(coulomb_cache), intent(in) :: cache
         !> Index of atom i
         integer, intent(in) :: iat
         !> Index of atom j
         integer, intent(in) :: jat
         !> Species of atom i
         integer, intent(in) :: isp
         !> Species of atom j
         integer, intent(in) :: jsp
         !> Multipole damping radius for the pair
         real(wp), intent(out) :: mrad
         !> Derivative of the multipole damping radius w.r.t. the CN of atom i
         real(wp), intent(out) :: dmraddcni
         !> Derivative of the multipole damping radius w.r.t. the CN of atom j
         real(wp), intent(out) :: dmraddcnj
      end subroutine get_mrad_derivs

      !> Evaluate damping function for a diatomic pair
      pure subroutine get_damping_pair(self, r, rinv, mrad, fdmp3, fdmp5, &
         & fdmp7, fdmp9)
         import :: damped_multipole, wp
         !> Instance of the multipole containeris
         class(damped_multipole), intent(in) :: self
         !> Interatomic distance
         real(wp), intent(in) :: r
         !> Inverse interatomic distance
         real(wp), intent(in) :: rinv
         !> Pairwise reference radius
         real(wp), intent(in) :: mrad
         !> Damping factor for inverse quadratic interactions
         real(wp), intent(out) :: fdmp3
         !> Damping factor for inverse cubic interactions
         real(wp), intent(out) :: fdmp5
         !> Damping factor for inverse quartic interactions
         real(wp), intent(out) :: fdmp7
         !> Damping factor for inverse quintic interactions
         real(wp), intent(out) :: fdmp9
      end subroutine get_damping_pair

      !> Evaluate derivatives of the damping function for a diatomic pair
      pure subroutine get_damping_derivs(self, r, rinv, mrad, fdmp3, fdmp5, &
         & fdmp7, fdmp9, dfdmp3dr, dfdmp5dr, dfdmp7dr, dfdmp9dr, &
         & dfdmp3dmrad, dfdmp5dmrad, dfdmp7dmrad, dfdmp9dmrad)
         import :: damped_multipole, wp
         !> Instance of the multipole container
         class(damped_multipole), intent(in) :: self
         !> Interatomic distance
         real(wp), intent(in) :: r
         !> Inverse interatomic distance
         real(wp), intent(in) :: rinv
         !> Averaged pairwise reference radius
         real(wp), intent(in) :: mrad
         !> Damping factor for inverse quadratic interactions
         real(wp), intent(out) :: fdmp3
         !> Damping factor for inverse cubic interactions
         real(wp), intent(out) :: fdmp5
         !> Damping factor for inverse quartic interactions
         real(wp), intent(out) :: fdmp7
         !> Damping factor for inverse quintic interactions
         real(wp), intent(out) :: fdmp9
         !> Derivative of damping factor for inverse quadratic interactions w.r.t. the distance
         real(wp), intent(out) :: dfdmp3dr
         !> Derivative of damping factor for inverse cubic interactions w.r.t. the distance
         real(wp), intent(out) :: dfdmp5dr
         !> Derivative of damping factor for inverse quartic interactions w.r.t. the distance
         real(wp), intent(out) :: dfdmp7dr
         !> Derivative of damping factor for inverse quintic interactions w.r.t. the distance
         real(wp), intent(out) :: dfdmp9dr
         !> Derivative of damping factor for inverse quadratic interactions w.r.t. the radius
         real(wp), intent(out) :: dfdmp3dmrad
         ! Derivative of damping factor for inverse cubic interactions w.r.t. the radius
         real(wp), intent(out) :: dfdmp5dmrad
         !> Derivative of damping factor for inverse quartic interactions w.r.t. the radius
         real(wp), intent(out) :: dfdmp7dmrad
         !> Derivative of damping factor for inverse quintic interactions w.r.t. the radius
         real(wp), intent(out) :: dfdmp9dmrad
      end subroutine get_damping_derivs
   end interface

   real(wp), parameter :: unity(3, 3) = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))
   real(wp), parameter :: conv = 100*eps
   real(wp), parameter :: cn_reg = 1e-6_wp

contains


!> Get anisotropic electrostatic energy
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   real(wp), allocatable :: vs(:), vd(:, :), vq(:, :)
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   allocate(vs(mol%nat), vd(3, mol%nat), vq(6, mol%nat))

   ! Damped multipole electrostatics
   call gemv(ptr%amat_sd, wfn%qat(:, 1), vd)
   call gemv(ptr%amat_dd, wfn%dpat(:, :, 1), vd, beta=1.0_wp, alpha=0.5_wp)
   call gemv(ptr%amat_sq, wfn%qat(:, 1), vq)
   call gemv(ptr%amat_dq, wfn%qpat(:, :, 1), vd, beta=1.0_wp)
   call gemv(ptr%amat_qq, wfn%qpat(:, :, 1), vq, beta=1.0_wp, alpha=0.5_wp/3.0_wp)

   energies(:) = energies + sum(wfn%dpat(:, :, 1) * vd, 1) &
      & + sum(wfn%qpat(:, :, 1) * vq, 1)

   ! Onsite multipole exchange-correlation kernel
   if (allocated(self%dkernel)) then
      call get_kernel_energy(mol, self%dkernel, wfn%dpat(:, :, 1), energies)
   end if
   if (allocated(self%qkernel)) then
      call get_kernel_energy(mol, self%qkernel, wfn%qpat(:, :, 1), energies)
   end if
end subroutine get_energy


!> Get anisotropic electrostatic energy
subroutine get_energy_aes(self, mol, cache, wfn, energies)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   real(wp), allocatable :: vd(:, :), vq(:, :)
   real(wp), allocatable :: mur(:), t1(:)
   real(wp), allocatable :: e01(:), e11(:), e02(:), e12(:), e22(:)
   type(coulomb_cache), pointer :: ptr
   integer :: i

   call view(cache, ptr)

   allocate(vd(3, mol%nat), vq(6, mol%nat))
   allocate(mur(mol%nat), t1(mol%nat), source=0.0_wp)
   allocate(e01(mol%nat), e11(mol%nat), e02(mol%nat), e12(mol%nat), e22(mol%nat), &
      source=0.0_wp)

   ! charge-dipole: e01 = q^T (A_sd^T d) + d^T (A_sd q)
   vd = 0.0_wp
   mur = 0.0_wp
   call gemv(ptr%amat_sd, wfn%qat(:, 1), vd)
   do i = 1, 3
      call gemv(ptr%amat_sd(i,:,:), wfn%dpat(i,:,1), mur, beta=1.0_wp, trans="T")
   end do
   e01 = mur*wfn%qat(:, 1) + sum(wfn%dpat(:,:,1) * vd, 1)

   ! dipole-dipole: e11 = d^T (A_dd d)
   vd = 0.0_wp
   call gemv(ptr%amat_dd, wfn%dpat(:, :, 1), vd, beta=1.0_wp)
   e11 = sum(wfn%dpat(:, :, 1) * vd, 1)

   ! charge-quadrupole: e02 = q^T (A_sq^T Q) + Q^T (A_sq q)
   vq = 0.0_wp
   t1 = 0.0_wp
   call gemv(ptr%amat_sq, wfn%qat(:, 1), vq)
   do i = 1, 6
      call gemv(ptr%amat_sq(i,:,:), wfn%qpat(i,:,1), t1, beta=1.0_wp, trans="T")
   end do
   e02 = t1*wfn%qat(:,1) + sum(wfn%qpat(:,:,1) * vq, 1)

   ! dipole-quadrupole: e12 = d^T (A_dq Q) + Q^T (A_dq^T d)
   vd = 0.0_wp
   vq = 0.0_wp
   call gemv(ptr%amat_dq, wfn%qpat(:, :, 1), vd, beta=1.0_wp)
   call gemv(ptr%amat_dq, wfn%dpat(:, :, 1), vq, beta=1.0_wp, trans="T")
   e12 = sum(wfn%dpat(:,:,1) * vd, 1) + sum(wfn%qpat(:,:,1) * vq, 1)

   ! quadrupole-quadrupole: e22 = Q^T (A_qq Q)
   vq = 0.0_wp
   call gemv(ptr%amat_qq, wfn%qpat(:, :, 1), vq, beta=1.0_wp)
   e22 = sum(wfn%qpat(:, :, 1) * vq, 1) / 3.0_wp

   energies(:) = energies(:) + 0.5_wp * (e01 + e11 + e02 + e12 + e22)
end subroutine get_energy_aes


!> Get only AXC part of the anisotropic energy
subroutine get_energy_axc(self, mol, wfn, energies)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)

   call get_kernel_energy(mol, self%dkernel, wfn%dpat(:, :, 1), energies)
   call get_kernel_energy(mol, self%qkernel, wfn%qpat(:, :, 1), energies)

end subroutine get_energy_axc


!> Get multipolar anisotropic exchange-correlation kernel energy
subroutine get_kernel_energy(mol, kernel, mpat, energies)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole kernel
   real(wp), intent(in) :: kernel(:)
   !> Atomic multipole momemnt
   real(wp), intent(in) :: mpat(:, :)
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)

   integer :: iat, izp
   real(wp) :: mpt(size(mpat, 1)), mpscale(size(mpat, 1))

   mpscale(:) = 1
   if (size(mpat, 1) == 6) mpscale([2, 4, 5]) = 2

   do iat = 1, mol%nat
      izp = mol%id(iat)
      mpt(:) = mpat(:, iat) * mpscale
      energies(iat) = energies(iat) + kernel(izp) * dot_product(mpt, mpat(:, iat))
   end do
end subroutine get_kernel_energy


!> Get anisotropic electrostatic potential
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   ! Damped multipole electrostatics
   call gemv(ptr%amat_sd, wfn%qat(:, 1), pot%vdp(:, :, 1), beta=1.0_wp)
   call gemv(ptr%amat_sd, wfn%dpat(:, :, 1), pot%vat(:, 1), beta=1.0_wp, trans="T")

   call gemv(ptr%amat_dd, wfn%dpat(:, :, 1), pot%vdp(:, :, 1), beta=1.0_wp)

   call gemv(ptr%amat_sq, wfn%qat(:, 1), pot%vqp(:, :, 1), beta=1.0_wp)
   call gemv(ptr%amat_sq, wfn%qpat(:, :, 1), pot%vat(:, 1), beta=1.0_wp, trans="T")

   call gemv(ptr%amat_dq, wfn%qpat(:, :, 1), pot%vdp(:, :, 1), beta=1.0_wp)
   call gemv(ptr%amat_dq, wfn%dpat(:, :, 1), pot%vqp(:, :, 1), beta=1.0_wp, trans="T")

   call gemv(ptr%amat_qq, wfn%qpat(:, :, 1), pot%vqp(:, :, 1), beta=1.0_wp, alpha=1.0_wp/3.0_wp)

   ! Onsite multipole exchange-correlation kernel
   if (allocated(self%dkernel)) then
      call get_kernel_potential(mol, self%dkernel, wfn%dpat(:, :, 1), pot%vdp(:, :, 1))
   end if
   if (allocated(self%qkernel)) then
      call get_kernel_potential(mol, self%qkernel, wfn%qpat(:, :, 1), pot%vqp(:, :, 1))
   end if
end subroutine get_potential


!> Get multipolar anisotropic potential contribution
subroutine get_kernel_potential(mol, kernel, mpat, vm)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Multipole kernel
   real(wp), intent(in) :: kernel(:)
   !> Atomic multipole momemnt
   real(wp), intent(in) :: mpat(:, :)
   !> Potential shoft on atomic multipole moment
   real(wp), intent(inout) :: vm(:, :)

   integer :: iat, izp
   real(wp) :: mpscale(size(mpat, 1))

   mpscale(:) = 1
   if (size(mpat, 1) == 6) mpscale([2, 4, 5]) = 2

   do iat = 1, mol%nat
      izp = mol%id(iat)
      vm(:, iat) = vm(:, iat) + 2*kernel(izp) * mpat(:, iat) * mpscale
   end do
end subroutine get_kernel_potential


!> Get derivatives of anisotropic electrostatics
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the multipole electrostatics energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the multipole electrostatics energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   real(wp), allocatable :: dEdcn(:)
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   allocate(dEdcn(mol%nat))
   dEdcn = 0.0_wp

   call get_multipole_gradient(self, mol, ptr, &
      & wfn%qat(:, 1), wfn%dpat(:, :, 1), wfn%qpat(:, :, 1), &
      & dEdcn, gradient, sigma)

   call gemv(ptr%dcndr, dEdcn, gradient, beta=1.0_wp)
   call gemv(ptr%dcndL, dEdcn, sigma, beta=1.0_wp)
end subroutine get_gradient


!> Get real lattice vectors
subroutine get_dir_trans(lattice, alpha, conv, trans)
   !> Lattice parameters
   real(wp), intent(in) :: lattice(:, :)
   !> Parameter for Ewald summation
   real(wp), intent(in) :: alpha
   !> Tolerance for Ewald summation
   real(wp), intent(in) :: conv
   !> Translation vectors
   real(wp), allocatable, intent(out) :: trans(:, :)

   call get_lattice_points([.true.], lattice, get_dir_cutoff(alpha, conv), trans)

end subroutine get_dir_trans

!> Get reciprocal lattice translations
subroutine get_rec_trans(lattice, alpha, volume, conv, trans)
   !> Lattice parameters
   real(wp), intent(in) :: lattice(:, :)
   !> Parameter for Ewald summation
   real(wp), intent(in) :: alpha
   !> Cell volume
   real(wp), intent(in) :: volume
   !> Tolerance for Ewald summation
   real(wp), intent(in) :: conv
   !> Translation vectors
   real(wp), allocatable, intent(out) :: trans(:, :)

   real(wp) :: rec_lat(3, 3)

   rec_lat = twopi*transpose(matinv_3x3(lattice))
   call get_lattice_points([.true.], rec_lat, get_rec_cutoff(alpha, volume, conv), trans)
   trans = trans(:, 2:)

end subroutine get_rec_trans


!> Get interaction matrix for all multipole moments up to quadrupoles
subroutine get_multipole_matrix(self, mol, cache, amat_sd, amat_dd, &
   & amat_sq, amat_dq, amat_qq)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(in) :: cache
   !> Interaction matrix for charges and dipoles
   real(wp), contiguous, intent(out) :: amat_sd(:, :, :)
   !> Interaction matrix for dipoles and dipoles
   real(wp), contiguous, intent(out) :: amat_dd(:, :, :, :)
   !> Interaction matrix for charges and quadrupoles
   real(wp), contiguous, intent(out) :: amat_sq(:, :, :)
   !> Interaction matrix for dipoles and quadrupoles
   real(wp), contiguous, intent(out) :: amat_dq(:, :, :, :)
   !> Interaction matrix for quadrupoles and quadrupoles
   real(wp), contiguous, intent(out) :: amat_qq(:, :, :, :)

   amat_sd(:, :, :) = 0.0_wp
   amat_dd(:, :, :, :) = 0.0_wp
   amat_sq(:, :, :) = 0.0_wp
   amat_dq(:, :, :, :) = 0.0_wp
   amat_qq(:, :, :, :) = 0.0_wp
   if (any(mol%periodic)) then
      call get_multipole_matrix_3d(self, mol, cache, cache%wsc, cache%alpha, &
         & amat_sd, amat_dd, amat_sq, amat_dq, amat_qq)
   else
      call get_multipole_matrix_0d(self, mol, cache, &
         & amat_sd, amat_dd, amat_sq, amat_dq, amat_qq)
   end if
end subroutine get_multipole_matrix

!> Calculate the multipole interaction matrix for finite systems
subroutine get_multipole_matrix_0d(self, mol, cache, &
   & amat_sd, amat_dd, amat_sq, amat_dq, amat_qq)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(in) :: cache
   !> Interaction matrix for charges and dipoles
   real(wp), intent(inout) :: amat_sd(:, :, :)
   !> Interaction matrix for dipoles and dipoles
   real(wp), intent(inout) :: amat_dd(:, :, :, :)
   !> Interaction matrix for charges and quadrupoles
   real(wp), intent(inout) :: amat_sq(:, :, :)
   !> Interaction matrix for dipoles and quadrupoles
   real(wp), intent(inout) :: amat_dq(:, :, :, :)
   !> Interaction matrix for quadrupoles and quadrupoles
   real(wp), intent(inout) :: amat_qq(:, :, :, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: r1, vec(3), tdd(3,3), tsq(6), tdq(3,6), tqq(6,6)
   real(wp) :: x2, y2, z2, xy, xz, yz, g1, g2, g3, g5, g7, g9
   real(wp) :: scale_di, scale_dj, sqrtcni, sqrtcnj
   real(wp) :: dg3, dg5, dg7, dg9, mrad, fdmp3, fdmp5, fdmp7, fdmp9
   
   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(self, cache, mol, amat_sd, amat_dd, amat_sq, amat_dq, amat_qq) &
   !$omp private(iat, jat, izp, jzp, r1, vec, g1, g2, g3, g5, g7, g9, sqrtcni) & 
   !$omp private(sqrtcnj, scale_di, scale_dj, mrad, fdmp3, fdmp5, fdmp7, fdmp9) & 
   !$omp private(x2, y2, z2, xy, xz, yz, dg3, dg5, dg7, dg9, tdd, tsq, tdq, tqq)
   do iat = 1, mol%nat
      do jat = 1, mol%nat
         if (iat == jat) cycle
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         izp = mol%id(iat)
         jzp = mol%id(jat)
         r1 = norm2(vec)
         g1 = 1.0_wp / r1
         g2 = g1 * g1
         g3 = g1 * g2
         g5 = g3 * g2
         g7 = g5 * g2
         g9 = g7 * g2

         ! Get the atomic dipole scaling factors
         sqrtcni = sqrt(cache%cn(iat) + cn_reg**2) - cn_reg
         sqrtcnj = sqrt(cache%cn(jat) + cn_reg**2) - cn_reg
         scale_di = 1.0_wp + self%dip_scale(izp) * sqrtcni
         scale_dj = 1.0_wp + self%dip_scale(jzp) * sqrtcnj

         ! Get the multipole damping radii
         call self%get_mrad_pair(cache, iat, jat, izp, jzp, mrad)

         ! Get the pairwise damping functions
         call self%get_damping_pair(r1, g1, mrad, fdmp3, fdmp5, &
            & fdmp7, fdmp9)

         ! Charge-dipole interaction
         dg3 = g3 * fdmp3 * scale_dj
         amat_sd(:, jat, iat) = amat_sd(:, jat, iat) + vec * dg3

         ! Dipole-dipole interaction
         x2 = vec(1)*vec(1)
         y2 = vec(2)*vec(2)
         z2 = vec(3)*vec(3)
         xy = vec(1)*vec(2)
         xz = vec(1)*vec(3)
         yz = vec(2)*vec(3)
         dg3 = g3 * fdmp5 * scale_dj * scale_di
         dg5 = g5 * fdmp5 * scale_dj * scale_di
         ! x-component of the dipole
         tdd(1, 1) = -3.0_wp*dg5*x2 + dg3
         tdd(2, 1) = -3.0_wp*dg5*xy
         tdd(3, 1) = -3.0_wp*dg5*xz
         ! y-component of the dipole
         tdd(1, 2) = -3.0_wp*dg5*xy
         tdd(2, 2) = -3.0_wp*dg5*y2 + dg3
         tdd(3, 2) = -3.0_wp*dg5*yz
         ! z-component of the dipole
         tdd(1, 3) = -3.0_wp*dg5*xz
         tdd(2, 3) = -3.0_wp*dg5*yz
         tdd(3, 3) = -3.0_wp*dg5*z2 + dg3

         amat_dd(:, jat, :, iat) = amat_dd(:, jat, :, iat) + tdd

         ! Charge-quadrupole interaction packed with doubled off-diagonal elements
         dg5 = g5 * fdmp5
         tsq(1) = dg5*x2
         tsq(2) = 2.0_wp*dg5*xy
         tsq(3) = dg5*y2
         tsq(4) = 2.0_wp*dg5*xz
         tsq(5) = 2.0_wp*dg5*yz
         tsq(6) = dg5*z2

         amat_sq(:, jat, iat) = amat_sq(:, jat, iat) + tsq 

         ! Dipole-quadrupole interaction packed with doubled off-diagonal elements
         dg7 = g7 * fdmp7 * scale_dj
         dg5 = g5 * fdmp7 * scale_dj
         ! xx-component of the quadrupole
         tdq(1, 1) = -5.0_wp*dg7*x2*vec(1) + 3.0_wp*dg5*vec(1)
         tdq(2, 1) = -5.0_wp*dg7*xy*vec(1) + dg5*vec(2)
         tdq(3, 1) = -5.0_wp*dg7*xz*vec(1) + dg5*vec(3)
         ! xy-component of the quadrupole
         tdq(1, 2) = -10.0_wp*dg7*x2*vec(2) + 2.0_wp*dg5*vec(2)
         tdq(2, 2) = -10.0_wp*dg7*y2*vec(1) + 2.0_wp*dg5*vec(1)
         tdq(3, 2) = -10.0_wp*dg7*xz*vec(2)
         ! yy-component of the quadrupole
         tdq(1, 3) = -5.0_wp*dg7*xy*vec(2) + dg5*vec(1)
         tdq(2, 3) = -5.0_wp*dg7*y2*vec(2) + 3.0_wp*dg5*vec(2)
         tdq(3, 3) = -5.0_wp*dg7*yz*vec(2) + dg5*vec(3)
         ! xz-component of the quadrupole
         tdq(1, 4) = -10.0_wp*dg7*x2*vec(3) + 2.0_wp*dg5*vec(3)
         tdq(2, 4) = -10.0_wp*dg7*xy*vec(3)
         tdq(3, 4) = -10.0_wp*dg7*z2*vec(1) + 2.0_wp*dg5*vec(1)
         ! yz-component of the quadrupole
         tdq(1, 5) = -10.0_wp*dg7*xy*vec(3)
         tdq(2, 5) = -10.0_wp*dg7*y2*vec(3) + 2.0_wp*dg5*vec(3)
         tdq(3, 5) = -10.0_wp*dg7*z2*vec(2) + 2.0_wp*dg5*vec(2)
         ! zz-component of the quadrupole
         tdq(1, 6) = -5.0_wp*dg7*xz*vec(3) + dg5*vec(1)
         tdq(2, 6) = -5.0_wp*dg7*yz*vec(3) + dg5*vec(2)
         tdq(3, 6) = -5.0_wp*dg7*z2*vec(3) + 3.0_wp*dg5*vec(3)

         ! Negative sign consistent with the definition of the interatomic vector j to i
         amat_dq(:, jat, :, iat) = amat_dq(:, jat, :, iat) - tdq

         ! Quadrupole-quadrupole interaction packed with doubled off-diagonal elements
         dg9 = g9 * fdmp9
         dg7 = g7 * fdmp9
         dg5 = g5 * fdmp9
         ! xx-component of the quadrupole
         tqq(1, 1) = 35.0_wp*dg9*x2*x2 - 30.0_wp*dg7*x2 + 4.5_wp*dg5
         tqq(2, 1) = 70.0_wp*dg9*x2*xy - 30.0_wp*dg7*xy
         tqq(3, 1) = 35.0_wp*dg9*x2*y2 - 5.0_wp*dg7*(x2+y2) + 1.5_wp*dg5
         tqq(4, 1) = 70.0_wp*dg9*x2*xz - 30.0_wp*dg7*xz
         tqq(5, 1) = 70.0_wp*dg9*x2*yz - 10.0_wp*dg7*yz
         tqq(6, 1) = 35.0_wp*dg9*x2*z2 - 5.0_wp*dg7*(x2+z2) + 1.5_wp*dg5
         ! xy-component of the quadrupole
         tqq(1, 2) = 70.0_wp*dg9*x2*xy - 30.0_wp*dg7*xy
         tqq(2, 2) = 140.0_wp*dg9*x2*y2 - 20.0_wp*dg7*(x2+y2) + 6.0_wp*dg5
         tqq(3, 2) = 70.0_wp*dg9*y2*xy - 30.0_wp*dg7*xy
         tqq(4, 2) = 140.0_wp*dg9*x2*yz - 20.0_wp*dg7*yz
         tqq(5, 2) = 140.0_wp*dg9*y2*xz - 20.0_wp*dg7*xz
         tqq(6, 2) = 70.0_wp*dg9*z2*xy - 10.0_wp*dg7*xy
         ! yy-component of the quadrupole
         tqq(1, 3) = 35.0_wp*dg9*x2*y2 - 5.0_wp*dg7*(x2+y2) + 1.5_wp*dg5
         tqq(2, 3) = 70.0_wp*dg9*y2*xy - 30.0_wp*dg7*xy
         tqq(3, 3) = 35.0_wp*dg9*y2*y2 - 30.0_wp*dg7*y2 + 4.5_wp*dg5
         tqq(4, 3) = 70.0_wp*dg9*y2*xz - 10.0_wp*dg7*xz
         tqq(5, 3) = 70.0_wp*dg9*y2*yz - 30.0_wp*dg7*yz
         tqq(6, 3) = 35.0_wp*dg9*y2*z2 - 5.0_wp*dg7*(y2+z2) + 1.5_wp*dg5
         ! xz-component of the quadrupole
         tqq(1, 4) = 70.0_wp*dg9*x2*xz - 30.0_wp*dg7*xz
         tqq(2, 4) = 140.0_wp*dg9*x2*yz - 20.0_wp*dg7*yz
         tqq(3, 4) = 70.0_wp*dg9*y2*xz - 10.0_wp*dg7*xz
         tqq(4, 4) = 140.0_wp*dg9*x2*z2 - 20.0_wp*dg7*(x2+z2) + 6.0_wp*dg5
         tqq(5, 4) = 140.0_wp*dg9*z2*xy - 20.0_wp*dg7*xy
         tqq(6, 4) = 70.0_wp*dg9*z2*xz - 30.0_wp*dg7*xz
         ! yz-component of the quadrupole
         tqq(1, 5) = 70.0_wp*dg9*x2*yz - 10.0_wp*dg7*yz
         tqq(2, 5) = 140.0_wp*dg9*y2*xz - 20.0_wp*dg7*xz
         tqq(3, 5) = 70.0_wp*dg9*y2*yz - 30.0_wp*dg7*yz
         tqq(4, 5) = 140.0_wp*dg9*z2*xy - 20.0_wp*dg7*xy
         tqq(5, 5) = 140.0_wp*dg9*y2*z2 - 20.0_wp*dg7*(y2+z2) + 6.0_wp*dg5
         tqq(6, 5) = 70.0_wp*dg9*z2*yz - 30.0_wp*dg7*yz
         ! zz-component of the quadrupole
         tqq(1, 6) = 35.0_wp*dg9*x2*z2 - 5.0_wp*dg7*(x2+z2) + 1.5_wp*dg5
         tqq(2, 6) = 70.0_wp*dg9*z2*xy - 10.0_wp*dg7*xy
         tqq(3, 6) = 35.0_wp*dg9*y2*z2 - 5.0_wp*dg7*(y2+z2) + 1.5_wp*dg5
         tqq(4, 6) = 70.0_wp*dg9*z2*xz - 30.0_wp*dg7*xz
         tqq(5, 6) = 70.0_wp*dg9*z2*yz - 30.0_wp*dg7*yz
         tqq(6, 6) = 35.0_wp*dg9*z2*z2 - 30.0_wp*dg7*z2 + 4.5_wp*dg5

         amat_qq(:, jat, :, iat) = amat_qq(:, jat, :, iat) + tqq

      end do
   end do

end subroutine get_multipole_matrix_0d

!> Evaluate multipole interaction matrix under 3D periodic boundary conditions
subroutine get_multipole_matrix_3d(self, mol, cache, wsc, alpha, &
   & amat_sd, amat_dd, amat_sq, amat_dq, amat_qq)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(in) :: cache
   !> Wigner-Seitz cell images
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence parameter for Ewald sum
   real(wp), intent(in) :: alpha
   !> Interaction matrix for charges and dipoles
   real(wp), intent(inout) :: amat_sd(:, :, :)
   !> Interaction matrix for dipoles and dipoles
   real(wp), intent(inout) :: amat_dd(:, :, :, :)
   !> Interaction matrix for charges and quadrupoles
   real(wp), intent(inout) :: amat_sq(:, :, :)
   !> Interaction matrix for dipoles and quadrupoles
   real(wp), intent(inout) :: amat_dq(:, :, :, :)
   !> Interaction matrix for quadrupoles and quadrupoles
   real(wp), intent(inout) :: amat_qq(:, :, :, :)

   integer :: iat, jat, izp, jzp, img, k
   real(wp) :: vec(3), rr, mrad, wsw, vol, sqrtcni, sqrtcnj, scale_di, scale_dj
   real(wp) :: d_sd(3), d_dd(3, 3), d_sq(6), r_sd(3), r_dd(3, 3), r_sq(6)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, alpha, conv, dtrans)
   call get_rec_trans(mol%lattice, alpha, vol, conv, rtrans)

   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(self, cache, amat_sd, amat_dd, amat_sq, amat_dq, amat_qq) &
   !$omp shared(mol, wsc, vol, alpha, rtrans, dtrans) &
   !$omp private(iat, jat, izp, jzp, img, vec, mrad, wsw, sqrtcni, sqrtcnj) &
   !$omp private(scale_di, scale_dj, d_sd, d_dd, d_sq, r_sd, r_dd, r_sq)
   do iat = 1, mol%nat
      do jat = 1, mol%nat
         izp = mol%id(iat)
         jzp = mol%id(jat)

         ! Get the atomic dipole scaling factors
         sqrtcni = sqrt(cache%cn(iat) + cn_reg**2) - cn_reg
         sqrtcnj = sqrt(cache%cn(jat) + cn_reg**2) - cn_reg
         scale_di = 1.0_wp + self%dip_scale(izp) * sqrtcni
         scale_dj = 1.0_wp + self%dip_scale(jzp) * sqrtcnj

         ! Get the multipole damping radii
         call self%get_mrad_pair(cache, iat, jat, izp, jzp, mrad)

         ! Iterator over all images of jat
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            
            call get_amat_sdq_rec_3d(vec, vol, alpha, rtrans, r_sd, r_dd, r_sq)
            call get_amat_sdq_dir_3d(self, vec, mrad, alpha, dtrans, d_sd, d_dd, d_sq)

            amat_sd(:, jat, iat) = amat_sd(:, jat, iat) + wsw * (d_sd + r_sd) * scale_dj
            amat_dd(:, jat, :, iat) = amat_dd(:, jat, :, iat) &
               & + wsw * (r_dd + d_dd) * scale_dj * scale_di
            amat_sq(:, jat, iat) = amat_sq(:, jat, iat) + wsw * (r_sq + d_sq)
         end do
      end do
   end do

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(self, cache, amat_sd, amat_dd, amat_sq, mol, vol, alpha) &
   !$omp private(iat, izp, rr, k, sqrtcni, scale_di)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ! Get the atomic dipole scaling factor
      sqrtcni = sqrt(cache%cn(iat) + cn_reg**2) - cn_reg
      scale_di = 1.0_wp + self%dip_scale(izp) * sqrtcni

      ! dipole-dipole selfenergy: -2/3·α³/sqrt(π) Σ(i) μ²(i)
      rr = -2.0_wp/3.0_wp * alpha**3 / sqrtpi
      do k = 1, 3
         amat_dd(k, iat, k, iat) = amat_dd(k, iat, k, iat) + 2*rr * scale_di**2
      end do

      ! charge-quadrupole selfenergy: 4/9·α³/sqrt(π) Σ(i) q(i)Tr(θi)
      ! (no actual contribution since quadrupoles are traceless)
      rr = 4.0_wp/9.0_wp * alpha**3 / sqrtpi
      amat_sq([1, 3, 6], iat, iat) = amat_sq([1, 3, 6], iat, iat) + rr
   end do
end subroutine get_multipole_matrix_3d

pure subroutine get_amat_sdq_rec_3d(rij, vol, alp, trans, amat_sd, amat_dd, amat_sq)
   !> Interatomic distance vector
   real(wp), intent(in) :: rij(3)
   !> Cell volume
   real(wp), intent(in) :: vol
   !> Ewald summation separation parameter
   real(wp), intent(in) :: alp
   !> Reciprocal space translation vectors
   real(wp), intent(in) :: trans(:, :)
   !> Interaction matrix for charges and dipoles
   real(wp), intent(out) :: amat_sd(:)
   !> Interaction matrix for dipoles and dipoles
   real(wp), intent(out) :: amat_dd(:, :)
   !> Interaction matrix for charges and quadrupoles
   real(wp), intent(out) :: amat_sq(:)

   integer :: itr
   real(wp) :: fac, vec(3), g2, expk, sink, cosk, gv

   amat_sd = 0.0_wp
   amat_dd = 0.0_wp
   amat_sq = 0.0_wp
   fac = 4*pi/vol

   amat_dd(1, 1) = fac/6.0_wp
   amat_dd(2, 2) = fac/6.0_wp
   amat_dd(3, 3) = fac/6.0_wp
   amat_sq([1, 3, 6]) = -fac/9.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      expk = fac*exp(-0.25_wp*g2/(alp*alp))/g2
      sink = sin(gv)*expk
      cosk = cos(gv)*expk

      amat_sd(:) = amat_sd + 2*vec*sink
      amat_dd(:, :) = amat_dd + spread(vec, 1, 3) * spread(vec, 2, 3) * cosk
      amat_sq(1) = amat_sq(1) +   vec(1)*vec(1)*cosk
      amat_sq(2) = amat_sq(2) + 2*vec(1)*vec(2)*cosk
      amat_sq(3) = amat_sq(3) +   vec(2)*vec(2)*cosk
      amat_sq(4) = amat_sq(4) + 2*vec(1)*vec(3)*cosk
      amat_sq(5) = amat_sq(5) + 2*vec(2)*vec(3)*cosk
      amat_sq(6) = amat_sq(6) +   vec(3)*vec(3)*cosk
   end do

end subroutine get_amat_sdq_rec_3d

pure subroutine get_amat_sdq_dir_3d(self, rij, mrad, alp, trans, &
   & amat_sd, amat_dd, amat_sq)
   class(damped_multipole), intent(in) :: self
   !> Interatomic distance vector
   real(wp), intent(in) :: rij(3)
   !> Multipole damping radius for the pair of atoms
   real(wp), intent(in) :: mrad
   !> Ewald summation separation parameter
   real(wp), intent(in) :: alp
   !> Real-space translation vectors
   real(wp), intent(in) :: trans(:, :)
   !> Interaction matrix for charges and dipoles
   real(wp), intent(out) :: amat_sd(:)
   !> Interaction matrix for dipoles and dipoles
   real(wp), intent(out) :: amat_dd(:, :)
   !> Interaction matrix for charges and quadrupoles
   real(wp), intent(out) :: amat_sq(:)

   integer :: itr
   real(wp) :: vec(3), r1, tmp, fdmp3, fdmp5, fdmp7, fdmp9
   real(wp) :: g1, g2, g3, g5, arg, arg2, alp2, e1, e2, erft, expt

   amat_sd = 0.0_wp
   amat_dd = 0.0_wp
   amat_sq = 0.0_wp
   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      g1 = 1.0_wp/r1
      g2 = g1 * g1
      g3 = g2 * g1
      g5 = g3 * g2

      ! Get the pairwise damping functions
      call self%get_damping_pair(r1, g1, mrad, fdmp3, fdmp5, &
         & fdmp7, fdmp9)

      arg = r1*alp
      arg2 = arg*arg
      expt = exp(-arg2)/sqrtpi
      erft = -erf(arg)*g1
      e1 = g1*g1 * (erft + 2*expt*alp)
      e2 = g1*g1 * (e1 + 4*expt*alp2*alp/3)

      tmp = fdmp3 * g3 + e1
      amat_sd = amat_sd + vec * tmp
      amat_dd(:, :) = amat_dd(:, :) + unity * (fdmp5*g3 + e1) &
         & - spread(vec, 1, 3) * spread(vec, 2, 3) * (3 * (g5*fdmp5 + e2))
      amat_sq(1) = amat_sq(1) +   vec(1)*vec(1)*(g5*fdmp5 + e2) - (fdmp5*g3 + e1)/3.0_wp
      amat_sq(2) = amat_sq(2) + 2*vec(1)*vec(2)*(g5*fdmp5 + e2)
      amat_sq(3) = amat_sq(3) +   vec(2)*vec(2)*(g5*fdmp5 + e2) - (fdmp5*g3 + e1)/3.0_wp
      amat_sq(4) = amat_sq(4) + 2*vec(1)*vec(3)*(g5*fdmp5 + e2)
      amat_sq(5) = amat_sq(5) + 2*vec(2)*vec(3)*(g5*fdmp5 + e2)
      amat_sq(6) = amat_sq(6) +   vec(3)*vec(3)*(g5*fdmp5 + e2) - (fdmp5*g3 + e1)/3.0_wp
   end do

end subroutine get_amat_sdq_dir_3d


!> Calculate derivatives of multipole interactions
subroutine get_multipole_gradient(self, mol, cache, qat, dpat, qpat, dEdcn, gradient, sigma)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(in) :: cache
   !> Atomic partial charges
   real(wp), contiguous, intent(in) :: qat(:)
   !> Atomic dipole moments
   real(wp), contiguous, intent(in) :: dpat(:, :)
   !> Atomic quadrupole moments
   real(wp), contiguous, intent(in) :: qpat(:, :)
   !> Derivative of the energy w.r.t. the coordination number
   real(wp), contiguous, intent(inout) :: dEdcn(:)
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   if (any(mol%periodic)) then
      call get_multipole_gradient_3d(self, mol, cache, qat, dpat, qpat, &
         & cache%wsc, cache%alpha, dEdcn, gradient, sigma)
   else
      call get_multipole_gradient_0d(self, mol, cache, qat, dpat, qpat, &
         & dEdcn, gradient, sigma)
   end if
end subroutine get_multipole_gradient

!> Evaluate multipole derivatives for finite systems
subroutine get_multipole_gradient_0d(self, mol, cache, qat, dpat, qpat, &
   & dEdcn, gradient, sigma)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(in) :: cache
   !> Atomic partial charges
   real(wp), contiguous, intent(in) :: qat(:)
   !> Atomic dipole moments
   real(wp), contiguous, intent(in) :: dpat(:, :)
   !> Atomic quadrupole moments
   real(wp), contiguous, intent(in) :: qpat(:, :)
   !> Derivative of the energy w.r.t. the coordination number
   real(wp), contiguous, intent(inout) :: dEdcn(:)
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r1, r2, g1, g2, g3, g5, g7, g9, g11
   real(wp) :: dG(3), dS(3, 3), dEdmrad, dEdsci, dEdscj, mrad, dmraddcni, dmraddcnj
   real(wp) :: fdmp3, fdmp5, fdmp7, fdmp9, dfdmp3dr, dfdmp5dr, dfdmp7dr, dfdmp9dr
   real(wp) :: dfdmp3dmrad, dfdmp5dmrad, dfdmp7dmrad, dfdmp9dmrad
   real(wp) :: esd, edd, esq, edq_g7, edq_g5, eqq_g9, eqq_g7, eqq_g5
   real(wp) :: sqrtcni, sqrtcnj, scale_di, scale_dj, dscale_didcn, dscale_djdcn
   real(wp) :: scale_dpati(3), scale_dpatj(3), dpiv, dpjv, dpi0v, dpj0v
   real(wp) :: dpidpj, dpi0dpj, dpj0dpi, qpiv(3), qpjv(3), vqpiv, vqpjv
   real(wp) :: qpidpj(3), qpjdpi(3), dpi0qpjv, dpj0qpiv, qpiqpjv(3), qpjqpiv(3)
   real(wp) :: tr_qpiqpj, dkernel3, dkernel5, dkernel7, dkernel9

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: dEdcn_local(:), gradient_local(:, :), sigma_local(:, :)

   !$omp parallel default(none) &
   !$omp shared(self, mol, cache, qat, dpat, qpat, dEdcn, gradient, sigma) &
   !$omp private(iat, jat, izp, jzp, r1, r2, vec, g1, g2, g3, g5, g7, g9, g11) &
   !$omp private(dG, dS, dEdmrad, dEdsci, dEdscj, mrad, dmraddcni, dmraddcnj) &
   !$omp private(fdmp3, fdmp5, fdmp7, fdmp9, dfdmp3dr, dfdmp5dr, dfdmp7dr, dfdmp9dr) &
   !$omp private(dfdmp3dmrad, dfdmp5dmrad, dfdmp7dmrad, dfdmp9dmrad) &
   !$omp private(esd, edd, esq, edq_g7, edq_g5, eqq_g9, eqq_g7, eqq_g5) &
   !$omp private(sqrtcni, sqrtcnj, scale_di, scale_dj, dscale_didcn, dscale_djdcn) &
   !$omp private(scale_dpati, scale_dpatj, dpiv, dpjv, dpi0v, dpj0v) &
   !$omp private(dpidpj, dpi0dpj, dpj0dpi, qpiv, qpjv, vqpiv, vqpjv) &
   !$omp private(qpidpj, qpjdpi, dpi0qpjv, dpj0qpiv, qpiqpjv, qpjqpiv) &
   !$omp private(tr_qpiqpj, dkernel3, dkernel5, dkernel7, dkernel9) &
   !$omp private(dEdcn_local, gradient_local, sigma_local)
   allocate(dEdcn_local(size(dEdcn)), source=0.0_wp)
   allocate(gradient_local(size(gradient, 1), size(gradient, 2)), source=0.0_wp)
   allocate(sigma_local(size(sigma, 1), size(sigma, 2)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)

      ! Get the scaled atomic dipole moments for atom i
      sqrtcni = sqrt(cache%cn(iat) + cn_reg**2) - cn_reg
      scale_di = 1.0_wp + self%dip_scale(izp) * sqrtcni
      dscale_didcn = self%dip_scale(izp) * 0.5_wp / (sqrtcni + cn_reg)
      scale_dpati(:) = dpat(:, iat) * scale_di

      do jat = 1, iat - 1
         jzp = mol%id(jat)

         ! Get the scaled atomic dipole moments for atom j
         sqrtcnj = sqrt(cache%cn(jat) + cn_reg**2) - cn_reg
         scale_dj = 1.0_wp + self%dip_scale(jzp) * sqrtcnj
         dscale_djdcn = self%dip_scale(jzp) * 0.5_wp / (sqrtcnj + cn_reg)
         scale_dpatj(:) = dpat(:, jat) * (1.0_wp + self%dip_scale(jzp) * sqrtcnj)

         ! Reset derivative accumulators
         dG(:) = 0.0_wp
         dS(:, :) = 0.0_wp
         dEdmrad = 0.0_wp
         dEdsci = 0.0_wp
         dEdscj = 0.0_wp

         ! Precompute the interatomic vector and distance
         vec(:) = mol%xyz(:, jat)-mol%xyz(:, iat)
         r1 = norm2(vec)
         r2 = r1 * r1
         ! Precompute powers of the inverse distance
         g1 = 1.0_wp/r1
         g2 = g1 * g1
         g3 = g2 * g1
         g5 = g3 * g2
         g7 = g5 * g2
         g9 = g7 * g2
         g11 = g9 * g2

         ! Get the multipole damping radii and their derivatives
         call self%get_mrad_derivs(cache, iat, jat, izp, jzp, mrad, &
            & dmraddcni, dmraddcnj)

         ! Get the pairwise damping functions and their derivatives
         call self%get_damping_derivs(r1, g1, mrad, fdmp3, fdmp5, &
            & fdmp7, fdmp9, dfdmp3dr, dfdmp5dr, dfdmp7dr, dfdmp9dr, &
            & dfdmp3dmrad, dfdmp5dmrad, dfdmp7dmrad, dfdmp9dmrad)

         ! Skip atom pairs with negligible multipole interactions
         if (abs(fdmp3) < eps) cycle

         ! Charge-dipole interaction
         ! Scalar products of (scaled) dipole moment and the interatomic vector
         dpiv = dot_product(scale_dpati, vec)
         dpi0v = dot_product(dpat(:, iat), vec)
         dpjv = dot_product(scale_dpatj, vec)
         dpj0v = dot_product(dpat(:, jat), vec)
         ! Scalar energy numerator for charge-dipole interaction
         esd = (qat(jat)*dpiv - qat(iat)*dpjv)
         ! Radial derivative of the charge-dipole interaction kernel
         dkernel3 = - 3.0_wp * g5 * fdmp3 + dfdmp3dr * g3 * g1
         ! Charge-dipole derivative w.r.t. the positions
         dG(:) = dG - dkernel3 * vec * esd &
            & + fdmp3 * g3 * (qat(iat)*scale_dpatj(:) - qat(jat)*scale_dpati(:))
         ! Charge-dipole derivative w.r.t. the multipole damping radius
         dEdmrad = dEdmrad + esd * dfdmp3dmrad * g3
         ! Charge-dipole derivative w.r.t. the dipole scaling factors
         dEdsci = dEdsci + qat(jat) * fdmp3 * g3 * dpi0v
         dEdscj = dEdscj - qat(iat) * fdmp3 * g3 * dpj0v

         ! Skip atom pairs with negligible multipole interactions
         if (abs(fdmp5) < eps) then 
            ! Calculate stress tensor contribution from the gradient for all terms
            dS(:, :) = dS - 0.5_wp * (spread(vec, 1, 3) * spread(dG, 2, 3) &
               & + spread(dG, 1, 3) * spread(vec, 2, 3))

            ! Accumulate contributions
            dEdcn_local(iat) = dEdcn_local(iat) + dEdmrad * dmraddcni &
               & + dEdsci * dscale_didcn
            dEdcn_local(jat) = dEdcn_local(jat) + dEdmrad * dmraddcnj &
               & + dEdscj * dscale_djdcn
            gradient_local(:, iat) = gradient_local(:, iat) + dG
            gradient_local(:, jat) = gradient_local(:, jat) - dG
            sigma_local(:, :) = sigma_local + dS
            cycle
         end if

         ! Dipole-dipole interaction
         ! Scalar product of (scaled) dipole moments
         dpidpj = dot_product(scale_dpatj, scale_dpati)
         dpi0dpj = dot_product(dpat(:, iat), scale_dpatj)
         dpj0dpi = dot_product(dpat(:, jat), scale_dpati)
         ! Scalar energy numerator for dipole-dipole interaction
         edd = dpidpj * r2 - 3.0_wp * dpjv * dpiv
         ! Radial derivative of the dipole-dipole interaction kernel
         dkernel5 = - 5.0_wp * g7 * fdmp5 + dfdmp5dr * g5 * g1
         ! Dipole-dipole derivative w.r.t. the positions
         dG(:) = dG - dkernel5 * vec * edd & 
            & - 2.0_wp * fdmp5 * g5 * dpidpj * vec & 
            & + 3.0_wp * fdmp5 * g5 * (dpiv * scale_dpatj(:) + dpjv * scale_dpati(:))
         ! Dipole-dipole derivative w.r.t. the multipole damping radius
         dEdmrad = dEdmrad + edd * g5 * dfdmp5dmrad
         ! Dipole-dipole derivative w.r.t. the dipole scaling factors
         dEdsci = dEdsci + fdmp5 * g5 * (r2 * dpi0dpj - 3.0_wp * dpjv * dpi0v)
         dEdscj = dEdscj + fdmp5 * g5 * (r2 * dpj0dpi - 3.0_wp * dpiv * dpj0v)


         ! Charge-quadrupole interaction
         ! Tensor products of quadrupole moments and interatomic vector
         qpiv(1) = qpat(1,iat)*vec(1) + qpat(2,iat)*vec(2) + qpat(4,iat)*vec(3)
         qpiv(2) = qpat(2,iat)*vec(1) + qpat(3,iat)*vec(2) + qpat(5,iat)*vec(3)
         qpiv(3) = qpat(4,iat)*vec(1) + qpat(5,iat)*vec(2) + qpat(6,iat)*vec(3)
         qpjv(1) = qpat(1,jat)*vec(1) + qpat(2,jat)*vec(2) + qpat(4,jat)*vec(3)
         qpjv(2) = qpat(2,jat)*vec(1) + qpat(3,jat)*vec(2) + qpat(5,jat)*vec(3)
         qpjv(3) = qpat(4,jat)*vec(1) + qpat(5,jat)*vec(2) + qpat(6,jat)*vec(3)
         ! Scalar vector-quadrupole-vector products
         vqpiv = dot_product(vec, qpiv)
         vqpjv = dot_product(vec, qpjv)
         ! Scalar energy numerator for the charge-quadrupole interaction
         esq = qat(jat) * vqpiv + qat(iat) * vqpjv
         ! Charge-quadrupole derivative w.r.t. the positions
         dG(:) = dG - dkernel5 * vec * esq &
            & - 2.0_wp * fdmp5 * g5 * (qat(iat) * qpjv(:) + qat(jat) * qpiv(:))
         ! Charge-quadrupole derivative w.r.t. the multipole damping radius
         dEdmrad = dEdmrad + esq * g5 * dfdmp5dmrad

         ! Skip atom pairs with negligible multipole interactions
         if (abs(fdmp7) < eps) then 
            ! Calculate stress tensor contribution from the gradient for all terms
            dS(:, :) = dS - 0.5_wp * (spread(vec, 1, 3) * spread(dG, 2, 3) &
               & + spread(dG, 1, 3) * spread(vec, 2, 3))

            ! Accumulate contributions
            dEdcn_local(iat) = dEdcn_local(iat) + dEdmrad * dmraddcni &
               & + dEdsci * dscale_didcn
            dEdcn_local(jat) = dEdcn_local(jat) + dEdmrad * dmraddcnj &
               & + dEdscj * dscale_djdcn
            gradient_local(:, iat) = gradient_local(:, iat) + dG
            gradient_local(:, jat) = gradient_local(:, jat) - dG
            sigma_local(:, :) = sigma_local + dS
            cycle
         end if
         
         ! Dipole-quadrupole interaction
         ! Tensor product of quadrupole and dipole moments
         qpidpj(1) = qpat(1,jat)*scale_dpati(1) + qpat(2,jat)*scale_dpati(2) &
            & + qpat(4,jat)*scale_dpati(3)
         qpidpj(2) = qpat(2,jat)*scale_dpati(1) + qpat(3,jat)*scale_dpati(2) &
            & + qpat(5,jat)*scale_dpati(3)
         qpidpj(3) = qpat(4,jat)*scale_dpati(1) + qpat(5,jat)*scale_dpati(2) &
            & + qpat(6,jat)*scale_dpati(3)
         qpjdpi(1) = qpat(1,iat)*scale_dpatj(1) + qpat(2,iat)*scale_dpatj(2) &
            & + qpat(4,iat)*scale_dpatj(3)
         qpjdpi(2) = qpat(2,iat)*scale_dpatj(1) + qpat(3,iat)*scale_dpatj(2) &
            & + qpat(5,iat)*scale_dpatj(3)
         qpjdpi(3) = qpat(4,iat)*scale_dpatj(1) + qpat(5,iat)*scale_dpatj(2) &
            & + qpat(6,iat)*scale_dpatj(3)
         ! Scalar product of unscaled dipole moment and quadrupole vector product
         dpi0qpjv = dot_product(dpat(:, iat), qpjv)
         dpj0qpiv = dot_product(dpat(:, jat), qpiv)
         ! Scalar energy numerator contributions for dipole-quadrupole interaction
         edq_g7 = - 5.0_wp * (dpiv * vqpjv - dpjv * vqpiv)
         edq_g5 = 2.0_wp * (dot_product(scale_dpati, qpjv) &
            & - dot_product(scale_dpatj, qpiv))
         ! Radial derivative of the dipole-quadrupole interaction kernel
         dkernel7 = - 7.0_wp * g9 * fdmp7 + dfdmp7dr * g7 * g1
         dkernel5 = - 5.0_wp * g7 * fdmp7 + dfdmp7dr * g5 * g1
         ! Dipole-quadrupole derivative w.r.t. the positions
         dG(:) = dG + vec * (edq_g7 * dkernel7 + edq_g5 * dkernel5) &
            & - 5.0_wp * fdmp7 * g7 * (scale_dpati(:)*vqpjv + 2.0_wp * dpiv * qpjv) &
            & + 5.0_wp * fdmp7 * g7 * (scale_dpatj(:)*vqpiv + 2.0_wp * dpjv * qpiv) &
            & + 2.0_wp * fdmp7 * g5 * (qpidpj(:) - qpjdpi(:))
         ! Dipole-quadrupole derivative w.r.t. the multipole damping radius
         dEdmrad = dEdmrad - (edq_g7 * g7 + edq_g5 * g5) * dfdmp7dmrad
         ! Dipole-quadrupole derivative w.r.t. the dipole scaling factors
         dEdsci = dEdsci - fdmp7 * (-5.0_wp*g7 * vqpjv * dpi0v + 2.0_wp*g5 * dpi0qpjv)
         dEdscj = dEdscj - fdmp7 * (+5.0_wp*g7 * vqpiv * dpj0v - 2.0_wp*g5 * dpj0qpiv)

         ! Skip atom pairs with negligible multipole interactions
         if (abs(fdmp9) < eps) then 
            ! Calculate stress tensor contribution from the gradient for all terms
            dS(:, :) = dS - 0.5_wp * (spread(vec, 1, 3) * spread(dG, 2, 3) &
               & + spread(dG, 1, 3) * spread(vec, 2, 3))

            ! Accumulate contributions
            dEdcn_local(iat) = dEdcn_local(iat) + dEdmrad * dmraddcni &
               & + dEdsci * dscale_didcn
            dEdcn_local(jat) = dEdcn_local(jat) + dEdmrad * dmraddcnj &
               & + dEdscj * dscale_djdcn
            gradient_local(:, iat) = gradient_local(:, iat) + dG
            gradient_local(:, jat) = gradient_local(:, jat) - dG
            sigma_local(:, :) = sigma_local + dS
            cycle
         end if

         ! Quadrupole-quadrupole interaction
         ! Tensor products of quadrupole moments and the quadrupole-vector products
         qpiqpjv(1) = qpat(1,iat)*qpjv(1) + qpat(2,iat)*qpjv(2) + qpat(4,iat)*qpjv(3)
         qpiqpjv(2) = qpat(2,iat)*qpjv(1) + qpat(3,iat)*qpjv(2) + qpat(5,iat)*qpjv(3)
         qpiqpjv(3) = qpat(4,iat)*qpjv(1) + qpat(5,iat)*qpjv(2) + qpat(6,iat)*qpjv(3)
         qpjqpiv(1) = qpat(1,jat)*qpiv(1) + qpat(2,jat)*qpiv(2) + qpat(4,jat)*qpiv(3)
         qpjqpiv(2) = qpat(2,jat)*qpiv(1) + qpat(3,jat)*qpiv(2) + qpat(5,jat)*qpiv(3)
         qpjqpiv(3) = qpat(4,jat)*qpiv(1) + qpat(5,jat)*qpiv(2) + qpat(6,jat)*qpiv(3)
         ! Trace of quadrupole-quadrupole product
         tr_qpiqpj = qpat(1,iat)*qpat(1,jat) + 2.0_wp*qpat(2,iat)*qpat(2,jat) & 
            & + qpat(3,iat)*qpat(3,jat) + 2.0_wp*qpat(4,iat)*qpat(4,jat) &
            & + 2.0_wp*qpat(5,iat)*qpat(5,jat) + qpat(6,iat)*qpat(6,jat)
         ! Scalar energy numerator contributions for quadrupole-quadrupole interaction
         eqq_g9 = (35.0_wp/3.0_wp) * vqpiv * vqpjv
         eqq_g7 = -(20.0_wp/3.0_wp) * dot_product(qpiv, qpjv)
         eqq_g5 = 1.0_wp * tr_qpiqpj
         ! Radial derivative of the quadrupole-quadrupole interaction kernel
         dkernel9 = - 9.0_wp * g11 * fdmp9 + dfdmp9dr * g9 * g1
         dkernel7 = - 7.0_wp * g9 * fdmp9 + dfdmp9dr * g7 * g1
         dkernel5 = - 5.0_wp * g7 * fdmp9 + dfdmp9dr * g5 * g1
         ! Quadrupole-quadrupole derivative w.r.t. the positions
         dG(:) = dG - vec * (eqq_g9 * dkernel9 + eqq_g7 * dkernel7 + eqq_g5 * dkernel5) &
            & - (70.0_wp/3.0_wp) * g9 * fdmp9 * (vqpjv * qpiv + vqpiv * qpjv) & 
            & + (20.0_wp/3.0_wp) * g7 * fdmp9 * (qpiqpjv + qpjqpiv)
         ! Quadrupole-quadrupole derivative w.r.t. the multipole damping radius
         dEdmrad = dEdmrad + (eqq_g9 * g9 + eqq_g7 * g7 + eqq_g5 * g5) * dfdmp9dmrad

         ! Calculate stress tensor contribution from the gradient for all terms
         dS(:, :) = dS - 0.5_wp * (spread(vec, 1, 3) * spread(dG, 2, 3) &
            & + spread(dG, 1, 3) * spread(vec, 2, 3))

         ! Accumulate contributions
         dEdcn_local(iat) = dEdcn_local(iat) + dEdmrad * dmraddcni &
            & + dEdsci * dscale_didcn
         dEdcn_local(jat) = dEdcn_local(jat) + dEdmrad * dmraddcnj &
            & + dEdscj * dscale_djdcn
         gradient_local(:, iat) = gradient_local(:, iat) + dG
         gradient_local(:, jat) = gradient_local(:, jat) - dG
         sigma_local(:, :) = sigma_local + dS

      end do
   end do
   !$omp end do
   !$omp critical (get_multipole_gradient_0d_)
   dEdcn(:) = dEdcn + dEdcn_local
   gradient(:, :) = gradient + gradient_local
   sigma(:, :) = sigma + sigma_local
   !$omp end critical (get_multipole_gradient_0d_)
   deallocate(dEdcn_local, gradient_local, sigma_local)
   !$omp end parallel

end subroutine get_multipole_gradient_0d

!> Evaluate multipole derivatives under 3D periodic boundary conditions
subroutine get_multipole_gradient_3d(self, mol, cache, qat, dpat, qpat, &
   & wsc, alpha, dEdcn, gradient, sigma)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(in) :: cache
   !> Atomic partial charges
   real(wp), contiguous, intent(in) :: qat(:)
   !> Atomic dipole moments
   real(wp), contiguous, intent(in) :: dpat(:, :)
   !> Atomic quadrupole moments
   real(wp), contiguous, intent(in) :: qpat(:, :)
   !> Wigner-Seitz cell images
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence parameter for Ewald sum
   real(wp), intent(in) :: alpha
   !> Derivative of the energy w.r.t. the coordination number
   real(wp), contiguous, intent(inout) :: dEdcn(:)
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: vec(3), dG(3), dGr(3), dGd(3), dS(3, 3), dSr(3, 3), dSd(3, 3)
   real(wp) :: wsw, vol, dEdmrad, dEdmradd, dEdsci, dEdscj, dEdscid, dEdscjd, rr
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)
   real(wp) :: sqrtcni, sqrtcnj, scale_di, scale_dj, dscale_didcn, dscale_djdcn
   real(wp) :: scale_dpati(3), scale_dpatj(3), mrad, dmraddcni, dmraddcnj

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: dEdcn_local(:), gradient_local(:, :), sigma_local(:, :)

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, alpha, conv, dtrans)
   call get_rec_trans(mol%lattice, alpha, vol, conv, rtrans)

   !$omp parallel default(none) shared(dEdcn, gradient, sigma) &
   !$omp shared(self, mol, cache, wsc, vol, alpha, rtrans, dtrans, qat, dpat, qpat) &
   !$omp private(iat, jat, izp, jzp, img, wsw, vec, scale_dpati, scale_dpatj) &
   !$omp private(sqrtcni, sqrtcnj, scale_di, scale_dj, dscale_didcn, dscale_djdcn) &
   !$omp private(mrad, dmraddcni, dmraddcnj, rr, dG, dGd, dGr, dS, dSd, dSr) &
   !$omp private(dEdmrad, dEdmradd, dEdsci, dEdscj, dEdscid, dEdscjd) & 
   !$omp private(dEdcn_local, gradient_local, sigma_local)
   allocate(dEdcn_local(size(dEdcn)), source=0.0_wp)
   allocate(gradient_local(size(gradient, 1), size(gradient, 2)), source=0.0_wp)
   allocate(sigma_local(size(sigma, 1), size(sigma, 2)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)

      ! Get the scaled atomic dipole moments for atom i
      sqrtcni = sqrt(cache%cn(iat) + cn_reg**2) - cn_reg
      scale_di = 1.0_wp + self%dip_scale(izp) * sqrtcni
      dscale_didcn = self%dip_scale(izp) * 0.5_wp / (sqrtcni + cn_reg)
      scale_dpati(:) = dpat(:, iat) * scale_di

      do jat = 1, iat - 1
         jzp = mol%id(jat)

         ! Get the scaled atomic dipole moments for atom j
         sqrtcnj = sqrt(cache%cn(jat) + cn_reg**2) - cn_reg
         scale_dj = 1.0_wp + self%dip_scale(jzp) * sqrtcnj
         dscale_djdcn = self%dip_scale(jzp) * 0.5_wp / (sqrtcnj + cn_reg)
         scale_dpatj(:) = dpat(:, jat) * (1.0_wp + self%dip_scale(jzp) * sqrtcnj)

         ! Reset derivative accumulators
         dG(:) = 0.0_wp
         dS(:, :) = 0.0_wp
         dEdmrad = 0.0_wp
         dEdsci = 0.0_wp
         dEdscj = 0.0_wp

         ! Get the multipole damping radii and their derivatives
         call self%get_mrad_derivs(cache, iat, jat, izp, jzp, mrad, &
            & dmraddcni, dmraddcnj)

         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            ! Precompute the interatomic vector and distance
            vec(:) = mol%xyz(:, jat)-mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))

            ! Compute reciprocal and direct space pairwise contributions
            call get_damat_sdq_rec_3d(vec, qat(iat), qat(jat), scale_dpati, &
               & scale_dpatj, qpat(:, iat), qpat(:, jat), vol, alpha, rtrans, &
               & dGr, dSr)
            call get_damat_sdq_dir_3d(self, vec, qat(iat), qat(jat), scale_dpati, &
               & scale_dpatj, dpat(:, iat), dpat(:, jat), qpat(:, iat), qpat(:, jat), &
               & mrad, alpha, dtrans, dEdmradd, dEdscid, dEdscjd, dGd, dSd)
            
            ! Scale derivative contributions with their number in the Wigner-Seitz cell
            dG = dG + (dGd + dGr) * wsw
            dS = dS + (dSd + dSr) * wsw
            dEdmrad = dEdmrad + dEdmradd * wsw
            dEdsci = dEdsci + dEdscid * wsw
            dEdscj = dEdscj + dEdscjd * wsw
         end do

         ! Accumulate contributions
         dEdcn_local(iat) = dEdcn_local(iat) + dEdmrad * dmraddcni &
            & + dEdsci * dscale_didcn
         dEdcn_local(jat) = dEdcn_local(jat) + dEdmrad * dmraddcnj &
            & + dEdscj * dscale_djdcn
         gradient_local(:, iat) = gradient_local(:, iat) + dG
         gradient_local(:, jat) = gradient_local(:, jat) - dG
         sigma_local(:, :) = sigma_local + dS
      end do

      ! Onsite correction for self interactions of atom i with itself
      ! Get the multipole damping radii and their derivatives
      call self%get_mrad_derivs(cache, iat, iat, izp, izp, mrad, dmraddcni, dmraddcnj)

      ! Reset derivative accumulators
      dS(:, :) = 0.0_wp
      dEdmrad = 0.0_wp
      dEdsci = 0.0_wp
      dEdscj = 0.0_wp

      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         ! Precompute the interatomic vector and distance
         vec(:) = wsc%trans(:, wsc%tridx(img, iat, iat))

         ! Compute reciprocal and direct space corrections
         call get_damat_sdq_rec_3d(vec, qat(iat), qat(iat), scale_dpati, &
            & scale_dpati, qpat(:, iat), qpat(:, iat), vol, alpha, rtrans, &
            & dGr, dSr)
         call get_damat_sdq_dir_3d(self, vec, qat(iat), qat(iat), scale_dpati, &
            & scale_dpati, dpat(:, iat), dpat(:, iat), qpat(:, iat), qpat(:, iat), &
            & mrad, alpha, dtrans, dEdmradd, dEdscid, dEdscjd, dGd, dSd)

         ! Scale derivative contributions with their number in the Wigner-Seitz cell
         dS = dS + (dSd + dSr) * wsw
         dEdmrad = dEdmrad + dEdmradd * wsw
         dEdsci = dEdsci + dEdscid * wsw
      end do

      ! Dipole-dipole self-interaction correction
      rr = -2.0_wp/3.0_wp * alpha**3 / sqrtpi
      dEdsci = dEdsci + 2.0_wp * rr * scale_di * dot_product(dpat(:, iat), dpat(:, iat))

      !> Accumulate contributions (gradient is zero due to symmetry)
      dEdcn_local(iat) = dEdcn_local(iat) + dEdmrad * dmraddcni &
         & + dEdsci * dscale_didcn
      sigma_local = sigma_local + 0.5_wp * dS
   end do
   !$omp end do
   !$omp critical (get_multipole_gradient_3d_)
   dEdcn(:) = dEdcn + dEdcn_local
   gradient(:, :) = gradient + gradient_local
   sigma(:, :) = sigma + sigma_local
   !$omp end critical (get_multipole_gradient_3d_)
   deallocate(dEdcn_local, gradient_local, sigma_local)
   !$omp end parallel

end subroutine get_multipole_gradient_3d

pure subroutine get_damat_sdq_rec_3d(rij, qi, qj, mi, mj, ti, tj, vol, alp, trans, dG, dS)
   !> Interatomic distance vector
   real(wp), intent(in) :: rij(3)
   !> Charge of atom i
   real(wp), intent(in) :: qi
   !> Charge of atom j
   real(wp), intent(in) :: qj
   !> Dipole moment of atom i
   real(wp), intent(in) :: mi(3)
   !> Dipole moment of atom j
   real(wp), intent(in) :: mj(3)
   !> Quadrupole moment of atom i
   real(wp), intent(in) :: ti(6)
   !> Quadrupole moment of atom j
   real(wp), intent(in) :: tj(6)
   !> Cell volume
   real(wp), intent(in) :: vol
   !> Ewald summation separation parameter
   real(wp), intent(in) :: alp
   !> Reciprocal space translation vectors
   real(wp), intent(in) :: trans(:, :)
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), intent(out) :: dG(3)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), intent(out) :: dS(3, 3)

   integer :: itr
   real(wp) :: fac, vec(3), g2, gv, expk, alp2, sink, cosk, dpiqj, qidpj
   real(wp) :: qpiqj, qiqpj, dpiv, dpjv

   dG(:) = 0.0_wp
   dS(:, :) = 0.0_wp
   fac = 4*pi/vol
   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      expk = fac*exp(-0.25_wp*g2/alp2)/g2
      cosk = cos(gv)*expk
      sink = sin(gv)*expk

      dpiqj = dot_product(vec, mi)*qj
      qidpj = dot_product(vec, mj)*qi

      dG(:) = dG - 2*vec*cosk * (dpiqj - qidpj)
      dS(:, :) = dS + 2 * sink * (dpiqj - qidpj) &
         & * ((2.0_wp/g2 + 0.5_wp/alp2) * spread(vec, 1, 3)*spread(vec, 2, 3) - unity)

      dpiv = dot_product(mi, vec)
      dpjv = dot_product(mj, vec)

      dG(:) = dG + vec*sink*dpiv*dpjv
      dS(:, :) = dS + cosk * dpiv*dpjv &
         & * ((2.0_wp/g2 + 0.5_wp/alp2) * spread(vec, 1, 3)*spread(vec, 2, 3) - unity)

      qiqpj = qi*(tj(1)*vec(1)*vec(1) + tj(3)*vec(2)*vec(2) + tj(6)*vec(3)*vec(3) &
         & + 2*tj(2)*vec(1)*vec(2) + 2*tj(4)*vec(1)*vec(3) + 2*tj(5)*vec(2)*vec(3))
      qpiqj = qj*(ti(1)*vec(1)*vec(1) + ti(3)*vec(2)*vec(2) + ti(6)*vec(3)*vec(3) &
         & + 2*ti(2)*vec(1)*vec(2) + 2*ti(4)*vec(1)*vec(3) + 2*ti(5)*vec(2)*vec(3))

      dG(:) = dG + vec * sink * (qiqpj + qpiqj)
      dS(:, :) = dS + cosk * (qiqpj + qpiqj) &
         & * ((2.0_wp/g2 + 0.5_wp/alp2) * spread(vec, 1, 3)*spread(vec, 2, 3) - unity)
   end do

end subroutine get_damat_sdq_rec_3d

pure subroutine get_damat_sdq_dir_3d(self, rij, qi, qj, scmi, scmj, mi, mj, ti, tj, &
   & mrad, alp, trans, dEdmrad, dEdsci, dEdscj, dG, dS)
   !> Instance of the multipole container
   class(damped_multipole), intent(in) :: self
   !> Interatomic distance vector
   real(wp), intent(in) :: rij(3)
   !> Charge of atom i
   real(wp), intent(in) :: qi
   !> Charge of atom j
   real(wp), intent(in) :: qj
   !> Scaled dipole moment of atom i
   real(wp), intent(in) :: scmi(3)
   !> Scaled dipole moment of atom j
   real(wp), intent(in) :: scmj(3)
   !> Unscaled dipole moment of atom i
   real(wp), intent(in) :: mi(3)
   !> Unscaled dipole moment of atom j
   real(wp), intent(in) :: mj(3)
   !> Quadrupole moment of atom i
   real(wp), intent(in) :: ti(6)
   !> Quadrupole moment of atom j
   real(wp), intent(in) :: tj(6)
   !> Multipole damping radius
   real(wp), intent(in) :: mrad
   !> Ewald summation separation parameter
   real(wp), intent(in) :: alp
   !> Direct space translation vectors
   real(wp), intent(in) :: trans(:, :)
   !> Derivative of the energy w.r.t. the multipole damping radius
   real(wp), intent(out) :: dEdmrad
   !> Derivative of the energy w.r.t. the dipole scaling factor of atom i
   real(wp), intent(out) :: dEdsci
   !> Derivative of the energy w.r.t. the dipole scaling factor of atom j
   real(wp), intent(out) :: dEdscj
   !> Derivative of the energy w.r.t. atomic displacements
   real(wp), intent(out) :: dG(3)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), intent(out) :: dS(3, 3)

   integer :: itr, k
   real(wp) :: vec(3), r1, r2, g1, g3, g5, g7
   real(wp) :: fdmp3, fdmp5, fdmp7, fdmp9, dfdmp3dr, dfdmp5dr, dfdmp7dr, dfdmp9dr
   real(wp) :: dfdmp3dmrad, dfdmp5dmrad, dfdmp7dmrad, dfdmp9dmrad
   real(wp) :: esd, edd, esq, g_sd(3), g_dd(3), g_sq(3)
   real(wp) :: dpiv, dpjv, dpi0v, dpj0v, dpidpj, dpi0dpj, dpj0dpi
   real(wp) :: qpiv(3), qpjv(3), vqpiv, vqpjv
   real(wp) :: dkernel3, dkernel5
   real(wp) :: alp2, arg, arg2, erft, expt, e1, e2, e3, tabc(3, 3, 3), tab(3, 3)

   dEdmrad = 0.0_wp
   dEdsci = 0.0_wp
   dEdscj = 0.0_wp
   dG(:) = 0.0_wp
   dS(:, :) = 0.0_wp

   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      ! Precompute the interatomic vector and distance
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      r2 = r1 * r1
      ! Precompute powers of the inverse distance
      g1 = 1.0_wp/r1
      g3 = g1 * g1 * g1
      g5 = g3 * g1 * g1
      g7 = g5 * g1 * g1

      ! Get the pairwise damping functions and their derivatives
      call self%get_damping_derivs(r1, g1, mrad, fdmp3, fdmp5, &
         & fdmp7, fdmp9, dfdmp3dr, dfdmp5dr, dfdmp7dr, dfdmp9dr, &
         & dfdmp3dmrad, dfdmp5dmrad, dfdmp7dmrad, dfdmp9dmrad)

      ! Ewald damping function
      arg = r1*alp
      arg2 = arg*arg
      erft = -erf(arg)*g1
      expt = exp(-arg2)/sqrtpi
      e1 = g1*g1 * (erft + expt*(2*alp2)/alp)
      e2 = g1*g1 * (e1 + expt*(2*alp2)**2/(3*alp))
      e3 = g1*g1 * (e2 + expt*(2*alp2)**3/(15*alp))

      ! Charge-dipole interaction
      ! Scalar products of (scaled) dipole moment and the interatomic vector
      dpiv = dot_product(scmi, vec)
      dpi0v = dot_product(mi, vec)
      dpjv = dot_product(scmj, vec)
      dpj0v = dot_product(mj, vec)
      ! Scalar energy numerator for charge-dipole interaction
      esd = dpiv*qj - dpjv*qi
      ! Radial derivative of the charge-dipole interaction kernel
      dkernel3 = - 3.0_wp * g5 * fdmp3 + dfdmp3dr * g3 * g1
      ! Charge-dipole derivative w.r.t. the positions
      g_sd(:)  = - dkernel3 * esd * vec + fdmp3 * g3 * (qi * scmj - qj * scmi)
      ! Ewald Correction Force for charge-dipole interaction
      block
         integer :: a, b
         do b = 1, 3
            do a = 1, 3
               tab(a, b) = 3*vec(a)*vec(b)*e2
            end do
            tab(b, b) = tab(b, b) - e1
         end do
      end block
      do k = 1, 3
         g_sd(k) = g_sd(k) &
            & + qj * (tab(1, k) * scmi(1) &
            &       + tab(2, k) * scmi(2) &
            &       + tab(3, k) * scmi(3))&
            & - qi * (tab(1, k) * scmj(1) &
            &       + tab(2, k) * scmj(2) &
            &       + tab(3, k) * scmj(3))
      end do
      dG(:) = dG + g_sd
      ! Charge-dipole stress tensor contribution
      dS(:, :) = dS - 0.5_wp * (spread(vec, 1, 3) * spread(g_sd, 2, 3) &
         & + spread(g_sd, 1, 3) * spread(vec, 2, 3))
      ! Charge-dipole derivative w.r.t. the multipole damping radius
      dEdmrad = dEdmrad + esd * dfdmp3dmrad * g3
      ! Charge-dipole derivative w.r.t. the dipole scaling factors
      dEdsci = dEdsci + qj * (fdmp3 * g3 + e1) * dpi0v
      dEdscj = dEdscj - qi * (fdmp3 * g3 + e1) * dpj0v


      ! Dipole-dipole interaction
      ! Scalar product of (scaled) dipole moments
      dpidpj = dot_product(scmj, scmi)
      dpi0dpj = dot_product(mi, scmj)
      dpj0dpi = dot_product(mj, scmi)
      ! Scalar energy numerator for dipole-dipole interaction
      edd = dpidpj * r2 - 3.0_wp * dpjv * dpiv
      ! Radial derivative of the dipole-dipole interaction kernel
      dkernel5 = - 5.0_wp * g7 * fdmp5 + dfdmp5dr * g5 * g1
      ! Dipole-dipole derivative w.r.t. the positions
      g_dd(:) = - dkernel5 * vec * edd &
         & - 2.0_wp * fdmp5 * g5 * dpidpj * vec &
         & + 3.0_wp * fdmp5 * g5 * (dpiv * scmj + dpjv * scmi)
      ! Ewald Correction Force for dipole-dipole interaction
      block
         integer :: a, b, c
         do c = 1, 3
            do b = 1, 3
               do a = 1, 3
                  tabc(a, b, c) = - 15*vec(a)*vec(b)*vec(c)*e3
               end do
            end do
            do a = 1, 3
               tabc(a, a, c) = tabc(a, a, c) + 3*e2*vec(c)
               tabc(c, a, c) = tabc(c, a, c) + 3*e2*vec(a)
               tabc(a, c, c) = tabc(a, c, c) + 3*e2*vec(a)
            end do
         end do
      end block
      do k = 1, 3
         g_dd(k) = g_dd(k) &
            & + mi(1)*(tabc(1, 1, k) * scmj(1) &
            &        + tabc(2, 1, k) * scmj(2) &
            &        + tabc(3, 1, k) * scmj(3))&
            & + mi(2)*(tabc(1, 2, k) * scmj(1) &
            &        + tabc(2, 2, k) * scmj(2) &
            &        + tabc(3, 2, k) * scmj(3))&
            & + mi(3)*(tabc(1, 3, k) * scmj(1) &
            &        + tabc(2, 3, k) * scmj(2) &
            &        + tabc(3, 3, k) * scmj(3))
      end do
      dG(:) = dG + g_dd
      ! Dipole-dipole stress tensor contribution
      dS(:, :) = dS - 0.5_wp * (spread(vec, 1, 3) * spread(g_dd, 2, 3) &
         & + spread(g_dd, 1, 3) * spread(vec, 2, 3))
      ! Dipole-dipole derivative w.r.t. the multipole damping radius
      dEdmrad = dEdmrad + edd * dfdmp5dmrad * g5
      ! Dipole-dipole derivative w.r.t. the dipole scaling factors
      dedsci = dedsci + (fdmp5 * g3 + e1) * dpi0dpj &
         & - 3.0_wp * (fdmp5 * g5 + e2) * dpjv * dpi0v
      dedscj = dedscj + (fdmp5 * g3 + e1) * dpj0dpi &
         & - 3.0_wp * (fdmp5 * g5 + e2) * dpiv * dpj0v


      ! Charge-quadrupole interaction
      ! Tensor products of quadrupole moments and interatomic vector
      qpiv(1) = ti(1)*vec(1) + ti(2)*vec(2) + ti(4)*vec(3)
      qpiv(2) = ti(2)*vec(1) + ti(3)*vec(2) + ti(5)*vec(3)
      qpiv(3) = ti(4)*vec(1) + ti(5)*vec(2) + ti(6)*vec(3)
      qpjv(1) = tj(1)*vec(1) + tj(2)*vec(2) + tj(4)*vec(3)
      qpjv(2) = tj(2)*vec(1) + tj(3)*vec(2) + tj(5)*vec(3)
      qpjv(3) = tj(4)*vec(1) + tj(5)*vec(2) + tj(6)*vec(3)
      ! Scalar vector-quadrupole-vector products
      vqpiv = dot_product(vec, qpiv)
      vqpjv = dot_product(vec, qpjv)
      ! Scalar energy numerator for the charge-quadrupole interaction
      esq = qj * vqpiv + qi * vqpjv
      ! Charge-quadrupole derivative w.r.t. the positions
      g_sq(:) = - dkernel5 * vec * esq &
         & - 2.0_wp * fdmp5 * g5 * (qi * qpjv(:) + qj * qpiv(:))
      ! Ewald Correction Force for charge-quadrupole interaction
      do k = 1, 3
         g_sq(k) = g_sq(k) + (&
            & - qi * (tabc(1, 1, k) * tj(1) &
            &     + 2*tabc(2, 1, k) * tj(2) &
            &     + 2*tabc(3, 1, k) * tj(4) &
            &     +   tabc(2, 2, k) * tj(3) &
            &     + 2*tabc(3, 2, k) * tj(5) &
            &     +   tabc(3, 3, k) * tj(6))&
            & - qj * (tabc(1, 1, k) * ti(1) &
            &     + 2*tabc(2, 1, k) * ti(2) &
            &     + 2*tabc(3, 1, k) * ti(4) &
            &     +   tabc(2, 2, k) * ti(3) &
            &     + 2*tabc(3, 2, k) * ti(5) &
            &     +   tabc(3, 3, k) * ti(6)))/3.0_wp
      end do
      dG(:) = dG + g_sq
      ! Charge-quadrupole stress tensor contribution
      dS(:, :) = dS - 0.5_wp * (spread(vec, 1, 3) * spread(g_sq, 2, 3) &
         & + spread(g_sq, 1, 3) * spread(vec, 2, 3))
      ! Charge-quadrupole derivative w.r.t. the multipole damping radius
      dEdmrad = dEdmrad + esq * dfdmp5dmrad * g5
   end do

end subroutine get_damat_sdq_dir_3d


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, atom_resolved
   !> Instance of the electrostatic container
   class(damped_multipole), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(dipole=atom_resolved, quadrupole=atom_resolved)
end function variable_info


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

end module tblite_coulomb_multipole_type
