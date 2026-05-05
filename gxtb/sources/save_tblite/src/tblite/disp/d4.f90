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

!> @file tblite/disp/d4.f90
!> Provides a proxy for the [DFT-D4 dispersion correction](https://dftd4.readthedocs.io)

!> Generally applicable charge-dependent London-dispersion correction, DFT-D4.
module tblite_disp_d4
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use dftd4, only : dispersion_model, d4_model, d4s_model, d4srev_model, &
      & dftd_models, new_d4_model, new_d4s_model, new_d4srev_model, &
      & damping_type, new_damping, twobody_damping_function, &
      & threebody_damping_function, param_type, get_damping_params, &
      & get_damping_function_id, realspace_cutoff
   use dftd4_model, only : d4_qmod
   use dftd4_cache, only : dftd4_cache_type => dispersion_cache
   use tblite_blas, only : dot, gemv
   use tblite_container_cache, only : container_cache
   use tblite_disp_cache, only : dispersion_cache
   use tblite_disp_type, only : dispersion_type
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_cutoff, only : get_lattice_points
   implicit none
   private

   public :: new_d4_dispersion, new_d4s_dispersion, new_d4srev_dispersion
   public :: twobody_damping_function, threebody_damping_function
   public :: get_damping_function_id


   !> Container for self-consistent D4 dispersion interactions
   type, public, extends(dispersion_type) :: d4_dispersion
      !> Instance of the actual D4 dispersion model
      class(dispersion_model), allocatable :: model
      !> Damping function
      type(damping_type) :: damp
      !> Damping parameters
      type(param_type) :: param
      !> Selected real space cutoffs for this instance
      type(realspace_cutoff) :: cutoff
      !> Coordination number instance
      class(ncoord_type), allocatable :: ncoord
   contains
      !> Update dispersion cache
      procedure :: update
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
      !> Evaluate non-selfconsistent part of the dispersion correction
      procedure :: get_engrad
      !> Evaluate selfconsistent energy of the dispersion correction
      procedure :: get_energy
      !> Evaluate charge dependent potential shift from the dispersion correction
      procedure :: get_potential
      !> Evaluate gradient contributions from the selfconsistent dispersion correction
      procedure :: get_gradient
   end type d4_dispersion

   character(len=*), parameter :: label_d4 = "self-consistent DFT-D4 dispersion"
   character(len=*), parameter :: label_d4s = "self-consistent DFT-D4S dispersion"
   character(len=*), parameter :: label_d4srev = "self-consistent DFT-D4Srev dispersion"


contains


!> Create a new instance of a self-consistent D4 dispersion correction
subroutine new_d4_dispersion(self, mol, damping_2b, damping_3b, s6, s8, &
   & a1, a2, a3, a4, s9, error)
   !> Instance of the dispersion correction
   type(d4_dispersion), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Two-body damping function
   integer, intent(in), optional :: damping_2b
   !> Three-body damping function
   integer, intent(in), optional :: damping_3b
   !> Damping parameters
   real(wp), intent(in), optional :: s6, s8, a1, a2, a3, a4, s9
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(d4_model), allocatable :: tmp
   integer :: d2b, d3b

   self%label = label_d4

   ! Create a new instance of the D4 model
   allocate(tmp)
   call new_d4_model(error, tmp, mol, qmod=d4_qmod%gfn2)
   if(allocated(error)) return
   call move_alloc(tmp, self%model)

   ! Setup the damping function for the D4 model
   if (present(damping_2b)) then
      d2b = damping_2b
   else
      d2b = twobody_damping_function%rational
   end if
   if (present(damping_3b)) then
      d3b = damping_3b
   else
      d3b = threebody_damping_function%zero_avg
   end if
   call new_damping(error, self%damp, d2b, d3b)
   if (allocated(error)) return

   ! Setup the specific damping parameters and fill up with defaults
   self%param = param_type(s6=s6, s8=s8, s9=s9, a1=a1, a2=a2, a3=a3, a4=a4, &
      & alp=16.0_wp)

   self%cutoff = realspace_cutoff(disp3=25.0_wp, disp2=50.0_wp)

   ! Setup the coordination number for the dispersion correction
   call new_ncoord(self%ncoord, mol, cn_count%dftd4, error, &
      & cutoff=self%cutoff%cn, rcov=self%model%rcov, en=self%model%en)
end subroutine new_d4_dispersion


!> Create a new instance of a self-consistent D4S dispersion correction
subroutine new_d4s_dispersion(self, mol, damping_2b, damping_3b, s6, s8, &
   & a1, a2, a3, a4, s9, error)
   !> Instance of the dispersion correction
   type(d4_dispersion), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Two-body damping function
   integer, intent(in), optional :: damping_2b
   !> Three-body damping function
   integer, intent(in), optional :: damping_3b
   !> Damping parameters
   real(wp), intent(in), optional :: s6, s8, a1, a2, a3, a4, s9
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(d4s_model), allocatable :: tmp
   integer :: d2b, d3b

   self%label = label_d4s

   ! Create a new instance of the D4S model
   allocate(tmp)
   call new_d4s_model(error, tmp, mol, qmod=d4_qmod%gfn2)
   if(allocated(error)) return
   call move_alloc(tmp, self%model)

   ! Setup the damping function for the D4S model
   if (present(damping_2b)) then
      d2b = damping_2b
   else
      d2b = twobody_damping_function%rational
   end if
   if (present(damping_3b)) then
      d3b = damping_3b
   else
      d3b = threebody_damping_function%zero_avg
   end if
   call new_damping(error, self%damp, d2b, d3b)
   if (allocated(error)) return

   ! Setup the specific damping parameters
   self%param = param_type(s6=s6, s8=s8, s9=s9, a1=a1, a2=a2, a3=a3, a4=a4, &
      & alp=16.0_wp)

   self%cutoff = realspace_cutoff(disp3=25.0_wp, disp2=50.0_wp)

   ! Setup the coordination number for the dispersion correction
   call new_ncoord(self%ncoord, mol, cn_count%dftd4, error, &
      & cutoff=self%cutoff%cn, rcov=self%model%rcov, en=self%model%en)
end subroutine new_d4s_dispersion


!> Create a new instance of a revised self-consistent D4Srev dispersion correction
subroutine new_d4srev_dispersion(self, mol, damping_2b, damping_3b, s6, s8, &
   & a1, a2, a3, a4, s9, kcn, error)
   !> Instance of the dispersion correction
   type(d4_dispersion), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Two-body damping function
   integer, intent(in), optional :: damping_2b
   !> Three-body damping function
   integer, intent(in), optional :: damping_3b
   !> Damping parameters
   real(wp), intent(in), optional :: s6, s8, a1, a2, a3, a4, s9
   !> Coordination number exponent
   real(wp), intent(in), optional :: kcn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(d4srev_model), allocatable :: tmp
   integer :: d2b, d3b

   self%label = label_d4srev

   ! Create a new instance of the revised D4Srev model
   allocate(tmp)
   call new_d4srev_model(error, tmp, mol, qmod=d4_qmod%gxtb)
   if(allocated(error)) return
   call move_alloc(tmp, self%model)

   ! Setup the damping function for the revised D4Srev model
   if (present(damping_2b)) then
      d2b = damping_2b
   else
      d2b = twobody_damping_function%screened
   end if
   if (present(damping_3b)) then
      d3b = damping_3b
   else
      d3b = threebody_damping_function%screened
   end if
   call new_damping(error, self%damp, d2b, d3b)
   if (allocated(error)) return

   ! Setup the specific damping parameters and fill up with defaults
   self%param = param_type(s6=s6, s8=s8, s9=s9, a1=a1, a2=a2, a3=a3, a4=a4, &
      & alp=16.0_wp)

   self%cutoff = realspace_cutoff(disp3=25.0_wp, disp2=50.0_wp)

   ! Setup the coordination number for the dispersion correction
   call new_ncoord(self%ncoord, mol, cn_count%erf, error, &
      & kcn=kcn, cutoff=self%cutoff%cn, rcov=self%model%rcov)
end subroutine new_d4srev_dispersion


!> Update dispersion cache
subroutine update(self, mol, cache)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   
   real(wp), allocatable :: lattr(:, :)
   type(dispersion_cache), pointer :: ptr

   call taint(cache, ptr)

   ! Allocate and compute coordination numbers
   if (.not.allocated(ptr%cn)) allocate(ptr%cn(mol%nat))
   if (.not.allocated(ptr%dcndr)) allocate(ptr%dcndr(3, mol%nat, mol%nat))
   if (.not.allocated(ptr%dcndL)) allocate(ptr%dcndL(3, 3, mol%nat))

   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%cn, lattr)
   call self%ncoord%get_coordination_number(mol, lattr, ptr%cn, ptr%dcndr, ptr%dcndL)

   ! Allocate and build damped dispersion matrices
   if (.not.allocated(ptr%dispmat6)) allocate(ptr%dispmat6(mol%nat, mol%nat))
   if (.not.allocated(ptr%dispmat8)) allocate(ptr%dispmat8(mol%nat, mol%nat))

   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%disp2, lattr)
   call get_dispersion_matrix(mol, self%model, self%damp, self%param, lattr, &
      & self%cutoff%disp2, ptr%dispmat6, ptr%dispmat8)
end subroutine update


!> Evaluate non-selfconsistent part of the dispersion correction
subroutine get_engrad(self, mol, cache, energies, gradient, sigma)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Dispersion energy
   real(wp), intent(inout) :: energies(:)
   !> Dispersion gradient
   real(wp), contiguous, intent(inout), optional :: gradient(:, :)
   !> Dispersion virial
   real(wp), contiguous, intent(inout), optional :: sigma(:, :)

   type(dispersion_cache), pointer :: ptr

   call view(cache, ptr)

   call get_dispersion_nonsc(mol, self%model, self%damp, self%param, self%cutoff, &
      & ptr, energies, gradient, sigma)

end subroutine get_engrad


!> Evaluate selfconsistent energy of the dispersion correction
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Dispersion energy
   real(wp), intent(inout) :: energies(:)

   type(dispersion_cache), pointer :: ptr
   integer :: iat, jat, izp, jzp
   real(wp) :: c6, c8

   call view(cache, ptr)

   call self%model%update(mol, ptr%disp_cache, ptr%cn, wfn%qat(:, 1))

   !$omp parallel do schedule(runtime) default(none) &
   !$omp reduction(+:energies) shared(self, mol, ptr) &
   !$omp private(iat, jat, izp, jzp, c6, c8)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         call self%model%get_2b_coeffs(ptr%disp_cache, iat, jat, izp, jzp, c6, c8)
         energies(iat) = energies(iat) + 0.5_wp * &
            & (c6 * ptr%dispmat6(iat, jat) + c8 * ptr%dispmat8(iat, jat))
      end do
   end do

end subroutine get_energy


!> Evaluate charge dependent potential shift from the dispersion correction
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   type(dispersion_cache), pointer :: ptr
   integer :: iat, jat, izp, jzp
   real(wp) :: c6, c8, dc6dcni, dc6dqi, dc6dcnj, dc6dqj
   real(wp) :: dc8dcni, dc8dqi, dc8dcnj, dc8dqj
   real(wp), allocatable :: tmp_vat(:)

   call view(cache, ptr)
   call self%model%update(mol, ptr%disp_cache, ptr%cn, wfn%qat(:, 1), grad=.true.)

   allocate(tmp_vat(mol%nat), source=0.0_wp)
   !$omp parallel do schedule(runtime) default(none) &
   !$omp reduction(+:tmp_vat) shared(self, mol, ptr) &
   !$omp private(iat, jat, izp, jzp, c6, c8, dc6dcni, dc6dqi, dc6dcnj, dc6dqj, &
   !$omp& dc8dcni, dc8dqi, dc8dcnj, dc8dqj)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         call self%model%get_2b_derivs(ptr%disp_cache, iat, jat, izp, jzp, c6, c8, &
            & dc6dcni, dc6dqi, dc6dcnj, dc6dqj, dc8dcni, dc8dqi, dc8dcnj, dc8dqj)
         tmp_vat(iat) = tmp_vat(iat) &
            & + dc6dqi * ptr%dispmat6(iat, jat) + dc8dqi * ptr%dispmat8(iat, jat)
      end do
   end do
   pot%vat(:, 1) = pot%vat(:, 1) + tmp_vat

end subroutine get_potential


!> Evaluate gradient contributions from the selfconsistent dispersion correction
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Dispersion gradient
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Dispersion virial
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   real(wp), allocatable :: dEdcn(:), dEdq(:), energies(:)
   real(wp), allocatable :: lattr(:, :)
   type(dispersion_cache), pointer :: ptr

   call view(cache, ptr)

   call self%model%update(mol, ptr%disp_cache, ptr%cn, wfn%qat(:, 1), grad=.true.)

   allocate(energies(mol%nat), dEdcn(mol%nat), dEdq(mol%nat))
   energies(:) = 0.0_wp
   dEdcn(:) = 0.0_wp
   dEdq(:) = 0.0_wp

   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%disp2, lattr)
   call self%model%get_dispersion2(mol, ptr%disp_cache, self%damp, self%param, &
      & lattr, self%cutoff%disp2, energies, dEdcn, dEdq, gradient, sigma)

   call gemv(ptr%dcndr, dEdcn, gradient, beta=1.0_wp)
   call gemv(ptr%dcndL, dEdcn, sigma, beta=1.0_wp)
end subroutine get_gradient


!> Build the damped dispersion matrices for C6/R^6 and C8/R^8 terms.
subroutine get_dispersion_matrix(mol, disp, damp, param, trans, cutoff, &
   & dispmat6, dispmat8)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Instance of the dispersion model
   class(dispersion_model), intent(in) :: disp
   !> Damping function
   type(damping_type), intent(in) :: damp
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Damped R^-6 dispersion matrix
   real(wp), intent(out) :: dispmat6(:, :)
   !> Damped R^-8 dispersion matrix
   real(wp), intent(out) :: dispmat8(:, :)

   integer :: iat, jat, izp, jzp, jtr
   real(wp) :: vec(3), r2, cutoff2, rdamp, d6, d8, dE6, dE8

   dispmat6(:, :) = 0.0_wp
   dispmat8(:, :) = 0.0_wp
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, param, disp, damp, trans, cutoff2, dispmat6, dispmat8) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, rdamp, d6, d8, dE6, dE8)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)
         rdamp = disp%get_2b_rdamp(izp, jzp)
         dE6 = 0.0_wp
         dE8 = 0.0_wp
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle

            call damp%damping_2b%get_2b_damp(param, r2, rdamp, d6, d8)

            dE6 = dE6 - d6
            dE8 = dE8 - d8
         end do

         dispmat6(iat, jat) = dE6
         dispmat6(jat, iat) = dE6
         dispmat8(iat, jat) = dE8
         dispmat8(jat, iat) = dE8
      end do
   end do
end subroutine get_dispersion_matrix


!> Wrapper to handle the evaluation of non-self-consistent dispersion energy and derivatives
subroutine get_dispersion_nonsc(mol, disp, damp, param, cutoff, cache, energies, &
   & gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Dispersion model
   class(dispersion_model), intent(in) :: disp
   !> Damping function
   type(damping_type), intent(in) :: damp
   !> Damping parameters
   type(param_type), intent(in) :: param
   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff
   !> Cached data between different dispersion runs
   type(dispersion_cache), intent(inout) :: cache
   !> Dispersion energy
   real(wp), intent(inout) :: energies(:)
   !> Dispersion gradient
   real(wp), intent(inout), contiguous, optional :: gradient(:, :)
   !> Dispersion virial
   real(wp), intent(inout), contiguous, optional :: sigma(:, :)

   logical :: grad
   real(wp), allocatable :: qat(:)
   real(wp), allocatable :: dEdcn(:), dEdq(:)
   real(wp), allocatable :: lattr(:, :)
   type(dftd4_cache_type) :: nonsc_cache

   grad = present(gradient) .and. present(sigma)

   ! Update model cache with zero charges for the non-SC three-body term
   allocate(qat(mol%nat), source=0.0_wp)
   call disp%update(mol, nonsc_cache, cache%cn, qat, grad=grad)

   if (grad) then
      allocate(dEdcn(mol%nat), dEdq(mol%nat))
      dEdcn(:) = 0.0_wp
      dEdq(:) = 0.0_wp
   end if

   call get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
   call disp%get_dispersion3(mol, nonsc_cache, damp, param, &
      & lattr, cutoff%disp3, energies, dEdcn, dEdq, gradient, sigma)

   if (grad) then
      call gemv(cache%dcndr, dEdcn, gradient, beta=1.0_wp)
      call gemv(cache%dcndL, dEdcn, sigma, beta=1.0_wp)
   end if

end subroutine get_dispersion_nonsc


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, atom_resolved
   !> Instance of the dispersion container
   class(d4_dispersion), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=atom_resolved)
end function variable_info


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(dispersion_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(dispersion_cache), allocatable :: tmp
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
   type(dispersion_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(dispersion_cache)
      ptr => target
   end select
end subroutine view

end module tblite_disp_d4
