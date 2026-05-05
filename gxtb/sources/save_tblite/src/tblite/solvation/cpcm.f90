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

!> @file tblite/solvation/cpcm.f90
!> Provides a polarizable continuum model

!> Implicit solvation model based on a polarizable dielectric continuum
module tblite_solvation_cpcm
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_blas, only : dot, gemv
   use tblite_container_cache, only : container_cache
   use tblite_mesh_lebedev, only : grid_size, get_angular_grid, list_bisection
   use tblite_scf_info, only : scf_info, atom_resolved
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_solvation_cpcm_dd
   use tblite_solvation_data, only : get_vdw_rad_cosmo
   use tblite_solvation_type, only : solvation_type
   implicit none
   private

   public :: cpcm_solvation, new_cpcm, cpcm_input


   !> Input for CPCM solvation
   type :: cpcm_input
      !> Dielectric constant
      real(wp) :: dielectric_const
      !> Number of grid points for each atom
      integer :: nang = grid_size(6)
      !> Scaling of van-der-Waals radii
      real(wp) :: rscale = 1.0_wp
      !> Accuracy for iterative solver
      real(wp) :: conv = 1.0e-8_wp
      !> Regularization parameter
      real(wp) :: eta = 2.0_wp
      !> Maximum angular momentum of basis functions
      integer :: lmax = 6
      !> Van-der-Waals radii for all species
      real(wp), allocatable :: rvdw(:)
   end type cpcm_input


   !> Definition of polarizable continuum model
   type, extends(solvation_type) :: cpcm_solvation
      !> Actual domain decomposition calculator
      type(domain_decomposition) :: dd
      !> Dielectric function
      real(wp) :: keps
      !> Van-der-Waal radii for all atoms
      real(wp), allocatable :: rvdw(:)
   contains
      !> Update cache from container
      procedure :: update
      !> Return dependency on density
      procedure :: variable_info
      !> Get electric field energy
      procedure :: get_energy
      !> Get electric field potential
      procedure :: get_potential
      !> Get electric field gradient
      procedure :: get_gradient
   end type cpcm_solvation

   !> Provide constructor for CPCM solvation
   interface cpcm_solvation
      module procedure :: create_cpcm
   end interface


   !> Restart data for CPCM calculation
   type :: cpcm_cache
      !> Actual domain decomposition calculator
      type(domain_decomposition) :: dd
      !> Electrostatic potential phi(ncav)
      real(wp), allocatable :: phi(:)
      !> Psi vector psi(nylm, n)
      real(wp), allocatable :: psi(:, :)
      !> CPCM solution sigma(nylm, n)
      real(wp), allocatable :: sigma(:, :)
      !> CPCM adjoint solution s(nylm, n)
      real(wp), allocatable :: s(:, :)
      !> Interaction matrix with surface charges jmat(ncav, n)
      real(wp), allocatable :: jmat(:, :)
   end type cpcm_cache


   !> Identifier for container
   character(len=*), parameter :: label = "polarizable continuum model"

   real(wp), parameter :: alpha_alpb = 0.571412_wp
   integer, parameter :: ndiis = 25


contains


!> Create new electric field container
subroutine new_cpcm(self, mol, input, error)
   !> Instance of the solvation model
   type(cpcm_solvation), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for CPCM solvation
   type(cpcm_input), intent(in) :: input
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: igrid, stat, iat, izp
   real(wp), allocatable :: ang_grid(:, :), ang_weight(:)
   type(domain_decomposition_input) :: ddinput

   self%label = label

   allocate(self%rvdw(mol%nat))
   if (allocated(input%rvdw)) then
      self%rvdw(:) = input%rscale*input%rvdw(mol%id)
   else
      do iat = 1, mol%nat
         izp = mol%num(mol%id(iat))
         self%rvdw(iat) = input%rscale*get_vdw_rad_cosmo(izp)
      end do
   end if
   self%keps = -0.5_wp * (1.0_wp/input%dielectric_const - 1.0_wp) / (1.0_wp + alpha_alpb)
   ! choose the lebedev grid with number of points closest to nAng:
   igrid = list_bisection(grid_size, input%nang)
   allocate(ang_grid(3, grid_size(igrid)))
   allocate(ang_weight(grid_size(igrid)))
   call get_angular_grid(igrid, ang_grid, ang_weight, stat)
   if (stat /= 0) then
      call fatal_error(error, "Could not initialize angular grid for CPCM model")
   end if

   ddinput = domain_decomposition_input(lmax=input%lmax, conv=input%conv, eta=input%eta)

   call new_domain_decomposition(self%dd, ddinput, self%rvdw, ang_weight, ang_grid)

end subroutine new_cpcm


!> Type constructor for CPCM splvation
function create_cpcm(mol, input) result(self)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for CPCM solvation
   type(cpcm_input), intent(in) :: input
   !> Instance of the solvation model
   type(cpcm_solvation) :: self

   type(error_type), allocatable :: error

   call new_cpcm(self, mol, input, error)
end function create_cpcm


!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the solvation model
   class(cpcm_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(cpcm_cache), pointer :: ptr

   call taint(cache, ptr)

   if (allocated(ptr%phi)) deallocate(ptr%phi)
   if (allocated(ptr%psi)) deallocate(ptr%psi)
   if (allocated(ptr%sigma)) deallocate(ptr%sigma)
   if (allocated(ptr%s)) deallocate(ptr%s)
   if (allocated(ptr%jmat)) deallocate(ptr%jmat)

   ptr%dd = self%dd

   call ddupdate(ptr%dd, mol%xyz)

   allocate(ptr%phi(ptr%dd%ncav), ptr%psi(ptr%dd%nylm, ptr%dd%nat))
   allocate(ptr%jmat(ptr%dd%ncav, ptr%dd%nat))

   call get_coulomb_matrix(mol%xyz, ptr%dd%ccav, ptr%jmat)
end subroutine update


!> Get electric field energy
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the solvation model
   class(cpcm_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Solvation free energy
   real(wp), intent(inout) :: energies(:)

   type(cpcm_cache), pointer :: ptr

   call view(cache, ptr)

   energies(:) = energies + self%keps * sum(ptr%sigma * ptr%psi, 1)
end subroutine get_energy


!> Get electric field potential
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the solvation model
   class(cpcm_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   real(wp) :: xx(1, 1)
   logical :: restart
   type(cpcm_cache), pointer :: ptr

   call view(cache, ptr)

   restart = allocated(ptr%sigma)
   if (.not.allocated(ptr%sigma)) then
      allocate(ptr%sigma(ptr%dd%nylm, ptr%dd%nat))
   end if

   call get_phi(wfn%qat(:, 1), ptr%jmat, ptr%phi)

   call solve_cosmo_direct(ptr%dd, .true., ptr%phi, xx, ptr%sigma, restart)

   restart = allocated(ptr%s)
   if (.not.allocated(ptr%s)) then
      allocate(ptr%s(ptr%dd%nylm, ptr%dd%nat))
   end if

   call get_psi(wfn%qat(:, 1), ptr%psi)

   ! solve adjoint ddCOSMO equation to get full potential contributions
   call solve_cosmo_adjoint(ptr%dd, ptr%psi, ptr%s, restart)

   ! we abuse Phi to store the unpacked and scaled value of s
   call get_zeta(ptr%dd, self%keps, ptr%s, ptr%phi)
   ! and contract with the Coulomb matrix
   call gemv(ptr%jmat, ptr%phi, pot%vat(:, 1), alpha=-1.0_wp, &
      & beta=1.0_wp, trans='t')

   pot%vat(:, 1) = pot%vat(:, 1) + (self%keps * sqrt(4*pi)) * ptr%sigma(1, :)

end subroutine get_potential


!> Get electric field gradient
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the solvation model
   class(cpcm_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the solvation free energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the solvation free energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: ii, iat, ig
   real(wp), allocatable :: gx(:, :), zeta(:), ef(:, :)
   type(cpcm_cache), pointer :: ptr

   call view(cache, ptr)

   allocate(gx(3, mol%nat), zeta(ptr%dd%ncav), ef(3, max(mol%nat, ptr%dd%ncav)))

   call solve_cosmo_adjoint(ptr%dd, ptr%psi, ptr%s, .true., &
      & accuracy=ptr%dd%conv*1e-3_wp)

   ! reset Phi
   call get_phi(wfn%qat(:, 1), ptr%jmat, ptr%phi)

   ! now call the routine that computes the ddcosmo specific contributions to the forces.
   call get_deriv(ptr%dd, self%keps, ptr%phi, ptr%sigma, ptr%s, gx)

   ! form the "zeta" intermediate
   call get_zeta(ptr%dd, self%keps, ptr%s, zeta)

   ! 1. solute's electric field at the cav points times zeta:
   !    compute the electric field
   call efld(mol%nat, wfn%qat(:, 1), ptr%dd%xyz, ptr%dd%ncav, ptr%dd%ccav, ef)

   ! contract it with the zeta intermediate
   ii = 0
   do iat = 1, ptr%dd%nat
      do ig = 1, ptr%dd%ngrid
         if (ptr%dd%ui(ig, iat) > 0.0_wp) then
            ii = ii + 1
            gx(:, iat) = gx(:, iat) + zeta(ii)*ef(:, ii)
         end if
      end do
   end do

   ! 2. "zeta's" electric field at the nuclei times the charges.
   !    compute the "electric field"
   call efld(ptr%dd%ncav, zeta, ptr%dd%ccav, mol%nat, ptr%dd%xyz, ef)

   ! contract it with the solute's charges.
   do iat = 1, ptr%dd%nat
      gx(:, iat) = gx(:, iat) + ef(:, iat)*wfn%qat(iat, 1)
   end do

   gradient(:, :) = gradient(:, :) + gx
end subroutine get_gradient


!> Return dependency on density
pure function variable_info(self) result(info)
   !> Instance of the solvation model
   class(cpcm_solvation), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=atom_resolved)
end function variable_info


subroutine taint(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(cpcm_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(cpcm_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

subroutine view(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(cpcm_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(cpcm_cache)
      ptr => target
   end select
end subroutine view


!> Evaluate the Coulomb interactions between the atomic sides (xyz) and the
!> surface elements of the cavity (ccav).
subroutine get_coulomb_matrix(xyz, ccav, jmat)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: ccav(:, :)
   real(wp), intent(inout) :: jmat(:, :)

   integer :: ic, j
   real(wp) :: vec(3), d2, d

   jmat(:, :) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(ccav, xyz, jmat) private(ic, j, vec, d2, d)
   do ic = 1, size(ccav, 2)
      do j = 1, size(xyz, 2)
         vec(:) = ccav(:, ic) - xyz(:, j)
         d2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         d = sqrt(d2)
         jmat(ic, j) = 1.0_wp / d
      end do
   end do

end subroutine get_coulomb_matrix


!> Routine to compute the psi vector
subroutine get_psi(charge, psi)
   real(wp), intent(in) :: charge(:)
   real(wp), intent(out) :: psi(:, :)

   integer :: iat
   real(wp), parameter :: fac = sqrt(4*pi)

   psi(:,:) = 0.0_wp

   do iat = 1, size(charge)
      psi(1, iat) = fac*charge(iat)
   end do

end subroutine get_psi


!> Routine to compute the potential vector
subroutine get_phi(charge, jmat, phi)
   real(wp), intent(in) :: charge(:)
   real(wp), intent(in) :: jmat(:, :)
   real(wp), intent(out) :: phi(:)

   phi(:) = 0.0_wp

   call gemv(jmat, charge, phi)

end subroutine get_phi


!> Computes the electric field produced by the sources src (nsrc point charges
!  with coordinates csrc) at the ntrg target points ctrg:
subroutine efld(nsrc, src, csrc, ntrg, ctrg, ef)
   integer, intent(in) :: nsrc, ntrg
   real(wp), intent(in) :: src(:)
   real(wp), intent(in) :: csrc(:, :)
   real(wp), intent(in) :: ctrg(:, :)
   real(wp), intent(inout) :: ef(:, :)

   integer :: i, j
   real(wp) :: vec(3), r2, rr, r3, f
   real(wp), parameter :: zero=0.0_wp

   ef(:, :) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp reduction(+:ef) shared(ntrg, nsrc, ctrg, csrc, src) &
   !$omp private(j, i, f, vec, r2, rr, r3)
   do j = 1, ntrg
      do i = 1, nsrc
         vec(:) = ctrg(:, j) - csrc(:, i)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         rr = sqrt(r2)
         r3 = r2*rr
         f = src(i)/r3
         ef(:, j) = ef(:, j) + f*vec
      end do
   end do

end subroutine efld


end module tblite_solvation_cpcm
