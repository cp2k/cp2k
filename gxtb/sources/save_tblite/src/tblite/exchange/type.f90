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

!> @file tblite/exchange/type.f90
!> Provides a base class for calculating approximated Fock exchange. 

!> Base-class for approximated Fock exchange interaction 
module tblite_exchange_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_blas, only : gemm, dot
   use tblite_container_cache, only : container_cache
   use tblite_container_type, only : container_type
   use tblite_exchange_cache, only : exchange_cache
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private


   !> General Mulliken approximated Fock exchange
   type, public, extends(container_type), abstract :: exchange_type
      !> Full-range admixture of exchange
      real(wp) :: frscale
      !> Range separation parameter 
      real(wp) :: omega
      !> Scaling of the range separated exchange
      real(wp) :: lrscale
      !> Number of AOS
      integer :: nao
      !> Number of shells
      integer :: nsh
      !> Maximum number of shells per species
      integer :: maxsh
      !> Number of shells for each species
      integer, allocatable :: nsh_id(:)
      !> Number of spherical atomic orbitals for each shell
      integer, allocatable :: nao_sh(:)
      !> Index offset for each atom in the shell space
      integer, allocatable :: ish_at(:)
      !> Index offset for each shell in the atomic orbital space
      integer, allocatable :: iao_sh(:)
   contains
      !> Update container cache
      procedure :: update
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
      !> Evaluate selfconsistent and overlap-dependent energy of the interaction
      procedure :: get_energy_w_overlap
      !> Evaluate density dependent potential
      procedure :: get_potential_w_overlap
      !> Evaluate gradient contributions from the selfconsistent interaction
      procedure :: get_gradient_w_overlap
      !> Information on container
      procedure :: info
      !> Evaluate Mulliken Fock exchange gamma matrix
      procedure(get_mulliken_Kmatrix), deferred :: get_mulliken_Kmatrix
      !> Evaluate onsite Fock exchange gamma matrix
      procedure(get_onsite_Kmatrix), deferred :: get_onsite_Kmatrix
      !> Evaluate bond-order correlation correction gamma matrix
      procedure(get_bocorr_Kmatrix), deferred :: get_bocorr_Kmatrix
      !> Evaluate the gradient of the Mulliken exchange energy
      procedure(get_mulliken_derivs), deferred :: get_mulliken_derivs
      !> Evaluate the gradient of the bond-order correlation correction energy
      procedure(get_bocorr_derivs), deferred :: get_bocorr_derivs
      !> Calculate exchange contribution to the Fock matrix
      procedure(get_KFock), deferred :: get_KFock
      !> Calculate exchange contribution to the gradient
      procedure(get_KGrad), deferred :: get_KGrad
   end type exchange_type

   abstract interface
   
      !> Evaluate Mulliken exchange gamma matrix
      subroutine get_mulliken_Kmatrix(self, mol, cache)
         import :: exchange_type, structure_type, exchange_cache
         !> Instance of the exchange container
         class(exchange_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container with intermediates and the final exchange matrix
         type(exchange_cache), intent(inout) :: cache
      end subroutine get_mulliken_Kmatrix

      !> Evaluate onsite exchange lambda matrix
      subroutine get_onsite_Kmatrix(self, mol, wfn, cache)
         import :: exchange_type, structure_type, wavefunction_type, exchange_cache
         !> Instance of the exchange container
         class(exchange_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Wavefunction data
         type(wavefunction_type), intent(in) :: wfn
         !> Reusable data container with intermediates and the final exchange matrices
         type(exchange_cache), intent(inout) :: cache
      end subroutine get_onsite_Kmatrix

      !> Evaluate bond-order correlation correction xi matrix
      subroutine get_bocorr_Kmatrix(self, mol, cache)
         import :: exchange_type, structure_type, exchange_cache, wp
         !> Instance of the exchange container
         class(exchange_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container with intermediates and the final exchange matrix
         type(exchange_cache), intent(inout) :: cache
      end subroutine get_bocorr_Kmatrix

      !> Evaluate the gradient of the Mulliken exchange energy
      subroutine get_mulliken_derivs(self, mol, cache, mulliken_grad, gradient, sigma)
         import :: exchange_type, structure_type, exchange_cache, wp
         !> Instance of the exchange container
         class(exchange_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container
         type(exchange_cache), intent(inout) :: cache
         !> Operator gradient w.r.t. the Mulliken gamma matrix
         real(wp), contiguous, intent(in) :: mulliken_grad(:, :)
         !> Molecular gradient of the exchange energy
         real(wp), contiguous, intent(inout) :: gradient(:, :)
         !> Strain derivatives of the exchange energy
         real(wp), contiguous, intent(inout) :: sigma(:, :)
      end subroutine get_mulliken_derivs

      !> Evaluate the gradient of the bond-order correlation correction energy
      subroutine get_bocorr_derivs(self, mol, cache, bocorr_grad, gradient, sigma)
         import :: exchange_type, structure_type, exchange_cache, wp
         !> Instance of the exchange container
         class(exchange_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container
         type(exchange_cache), intent(inout) :: cache
         !> Operator gradient w.r.t. the bond-order correlation matrix
         real(wp), contiguous, intent(in) :: bocorr_grad(:, :)
         !> Molecular gradient of the exchange energy
         real(wp), contiguous, intent(inout) :: gradient(:, :)
         !> Strain derivatives of the exchange energy
         real(wp), contiguous, intent(inout) :: sigma(:, :)
      end subroutine get_bocorr_derivs

      !> Calculate exchange contribution to the Fock matrix and atomic potential
      subroutine get_KFock(self, mol, cache, density, overlap)
         import :: exchange_type, structure_type, exchange_cache, wp
         !> Instance of the exchange container
         class(exchange_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container with intermediates and the final Fock matrix/potential
         type(exchange_cache), intent(inout) :: cache
         !> Density matrix
         real(wp), intent(in) :: density(:, :, :)
         !> Overlap matrix
         real(wp), intent(in) :: overlap(:, :)
      end subroutine get_KFock

      !> Calculate exchange contribution to the gradient
      subroutine get_KGrad(self, mol, cache, density, overlap, mulliken_grad, &
         & bocorr_grad, ao_grad)
         import :: exchange_type, structure_type, exchange_cache, wp
         !> Instance of the exchange container
         class(exchange_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container with intermediates
         type(exchange_cache), intent(inout) :: cache
         !> Density matrix
         real(wp), intent(in) :: density(:, :, :)
         !> Overlap matrix
         real(wp), intent(in) :: overlap(:, :)
         !> Operator gradient w.r.t. the Mulliken gamma matrix
         real(wp), contiguous, intent(out) :: mulliken_grad(:, :)
         !> Operator gradient w.r.t. the bond-order correlation matrix
         real(wp), contiguous, intent(out) :: bocorr_grad(:, :)
         !> Orbital gradient contribution to the energy-weighted density matrix
         real(wp), contiguous, intent(inout) :: ao_grad(:, :)
      end subroutine get_KGrad

   end interface

contains


!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the exchange container
   class(exchange_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(exchange_cache), pointer :: ptr

   call taint(cache, ptr)
   call ptr%update(mol)

   ! Mulliken exchange matrix
   if (.not.allocated(ptr%g_mulliken)) then
      allocate(ptr%g_mulliken(self%nsh, self%nsh))
   end if
   call self%get_mulliken_Kmatrix(mol, ptr)

   ! Onsite and rotational invariance correction exchange matrices
   if (.not.allocated(ptr%g_onsfx)) then
      allocate(ptr%g_onsfx(self%maxsh, self%maxsh, mol%nat))
   end if
   if (.not.allocated(ptr%dgdq_onsfx)) then
      allocate(ptr%dgdq_onsfx(self%maxsh, self%maxsh, self%nsh))
   end if
   if (.not.allocated(ptr%g_onsri)) then
      allocate(ptr%g_onsri(self%maxsh, mol%nat))
   end if
   if (.not.allocated(ptr%dgdq_onsri)) then
      allocate(ptr%dgdq_onsri(self%maxsh, self%nsh))
   end if
 
   ! Bond-order correlation correction matrix
   if (.not.allocated(ptr%g_bocorr)) then
      allocate(ptr%g_bocorr(mol%nat, mol%nat))
   end if 
   call self%get_bocorr_Kmatrix(mol, ptr)

end subroutine update


!> Evaluate selfconsistent energy of the interaction
subroutine get_energy_w_overlap(self, mol, cache, wfn, overlap, energies)
   !> Instance of the exchange container
   class(exchange_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Overlap matrix
   real(wp), contiguous, intent(in) :: overlap(:, :)
   !> Exchange energy
   real(wp), intent(inout) :: energies(:)
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(exchange_cache), pointer :: ptr
   integer :: spin, iat, izp, is, ish, ii, iao
   real(wp) :: tmp_energy

   call view(cache, ptr)

   if (.not.allocated(ptr%prev_F) .or. .not.allocated(ptr%prev_vsh)) then
      ! Allocate temporary Fock matrix/potential storage of for the energy evaluation
      allocate(ptr%prev_F(self%nao, self%nao, wfn%nspin))
      allocate(ptr%prev_vsh(self%nsh, wfn%nspin))
   end if

   ! Build the charge-dependent onsite exchange matrices with current charges
   call self%get_onsite_Kmatrix(mol, wfn, ptr)
   ! Calculate the Fock matrix contribution for the initial current density 
   call self%get_KFock(mol, ptr, wfn%density, overlap)

   !$omp parallel do collapse(2) schedule(runtime) default(none) &
   !$omp shared(self, mol, ptr, wfn, energies) &
   !$omp private(spin, iat, izp, is, ish, ii, iao, tmp_energy)
   do spin = 1, size(wfn%density, 3)
      do iat = 1, mol%nat
         tmp_energy = 0.0_wp
         izp = mol%id(iat)
         is = self%ish_at(iat)
         do ish = 1, self%nsh_id(izp)
            ii = self%iao_sh(is+ish)
            do iao = 1, self%nao_sh(is + ish)
               tmp_energy = tmp_energy + 0.5_wp &
                  & * dot(ptr%prev_F(:, ii+iao, spin), wfn%density(:, ii+iao, spin))
            end do
         end do
         !$omp atomic
         energies(iat) = energies(iat) + tmp_energy
      end do
   end do
   
end subroutine get_energy_w_overlap


!> Evaluate density dependent potential 
subroutine get_potential_w_overlap(self, mol, cache, wfn, overlap, pot)
   !> Instance of the exchange container
   class(exchange_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Tight-binding wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential shifts
   type(potential_type), intent(inout) :: pot
   !> Overlap matrix
   real(wp), contiguous, intent(in) :: overlap(:, :)
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(exchange_cache), pointer :: ptr

   call view(cache, ptr)

   ! First iteration, the temporary Fock matrix/potential storage is not allocated
   if (.not.allocated(ptr%prev_F) .or. .not.allocated(ptr%prev_vsh)) then
      ! Allocate temporary Fock matrix/potential storage of for the energy evaluation
      allocate(ptr%prev_F(self%nao, self%nao, wfn%nspin))
      allocate(ptr%prev_vsh(self%nsh, wfn%nspin))
      ! Build the charge-dependent onsite exchange matrices with guess charges
      call self%get_onsite_Kmatrix(mol, wfn, ptr)
      ! Calculate the Fock matrix contribution for the initial guess density
      call self%get_KFock(mol, ptr, wfn%density, overlap)
   end if

   ! Add the exchange contribution from the previous iteration energy evaluation
   pot%kao(:, :, :) = pot%kao + ptr%prev_F
   pot%vsh(:, :) = pot%vsh + ptr%prev_vsh

end subroutine get_potential_w_overlap


!> Evaluate gradient contributions from the selfconsistent interaction
subroutine get_gradient_w_overlap(self, mol, cache, wfn, overlap, &
   & ao_grad, gradient, sigma)
   !> Instance of the exchange container
   class(exchange_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Overlap matrix
   real(wp), contiguous, intent(in) :: overlap(:, :)
   !> Orbital gradient contribution to the energy-weighted density matrix
   real(wp), contiguous, intent(inout) :: ao_grad(:, :)
   !> Molecular gradient of the exchange energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   type(exchange_cache), pointer :: ptr
   real(wp), allocatable :: mulliken_grad(:, :), bocorr_grad(:, :)

   call view(cache, ptr)

   allocate(mulliken_grad(self%nao, self%nao), bocorr_grad(self%nao, self%nao))

   ! Calculate direct and indirect overlap matrix derivatives of the exchange energy
   call self%get_KGrad(mol, ptr, wfn%density, overlap, mulliken_grad, bocorr_grad, &
      & ao_grad)

   ! Calculate gradient and strain contributions of the Mulliken approximate exchange
   call self%get_mulliken_derivs(mol, ptr, mulliken_grad, gradient, sigma)

   ! Calculate gradient and strain contributions of the bond-order correlation correction
   call self%get_bocorr_derivs(mol, ptr, bocorr_grad, gradient, sigma)

end subroutine get_gradient_w_overlap


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, shell_resolved, orbital_resolved, not_used
   !> Instance of the exchange container
   class(exchange_type), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=shell_resolved, dipole=not_used, quadrupole=not_used, &
      & density=orbital_resolved)
end function variable_info


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(exchange_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(exchange_cache), allocatable :: tmp
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
   type(exchange_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(exchange_cache)
      ptr => target
   end select
end subroutine view


!> Information on container
pure function info(self, verbosity, indent) result(str)
   use tblite_output_format, only : format_string
   !> Instance of the exchange container
   class(exchange_type), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   character(len=*), parameter :: nl = new_line('a'), marker = " * "

   if (allocated(self%label)) then
      str = self%label
   else
      str = "Approximated Fock exchange"
   end if

   if (verbosity > 1) then
      str = str // nl // indent // "Using the shell resolved chemical hardness"
      if (self%omega.ne.0.0_wp .and. self%lrscale.ne.0.0_wp) then
         str= str // nl // indent // "Range separted exchange is used:" // &
            & nl // indent // marker // "Full-range scale: " // format_string(self%frscale,'(f5.2)') // &
            & nl // indent // marker // "Omega: " // format_string(self%omega,'(f5.2)') // &
            & nl // indent // marker // "Long-range scale: " // format_string(self%lrscale,'(f5.2)')
      else
         str = str // nl // indent // "Full range exchange is used" // &
            & nl // indent // marker // "Full-range scale: " // format_string(self%frscale,'(f5.2)')
      end if
   end if 
end function info

end module tblite_exchange_type
