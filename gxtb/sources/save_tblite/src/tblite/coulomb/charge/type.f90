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

!> @file tblite/coulomb/charge/type.f90
!> Provides a general base class for electrostatic interactions

!> General isotropic second-order electrostatics model
module tblite_coulomb_charge_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_blas, only : gemv, symv
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_type, only : coulomb_type
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private


   !> General second-order electrostatics
   type, public, extends(coulomb_type), abstract :: coulomb_charge_type
      !> Chemical hardness for each shell and species
      real(wp), allocatable :: hubbard(:, :)
      !> Whether the second order contribution is shell-dependent
      logical :: shell_resolved
      !> Number of shells for each atom
      integer, allocatable :: nshell(:)
      !> Index offset for each shell
      integer, allocatable :: offset(:)
      !> Shell-distributed effective valence reference charges
      !> for first/pseudo-second-order offsite tight-binding
      real(wp), allocatable :: refqsh(:)
      !> Logical indicating whether there is a CN dependence 
      logical :: cn_dep
   contains
      !> Update container cache
      procedure :: update
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
      !> Evaluate selfconsistent energy of the interaction
      procedure :: get_energy
      !> Evaluate charge dependent potential shift from the interaction
      procedure :: get_potential
      !> Evaluate gradient of the charge dependent potential shift from the interaction
      procedure :: get_potential_gradient
      !> Evaluate gradient contributions from the selfconsistent interaction
      procedure :: get_gradient
      !> Evaluate Coulomb matrix
      procedure(get_coulomb_matrix), deferred :: get_coulomb_matrix
      !> Evaluate uncontracted derivatives of Coulomb matrix
      procedure(get_coulomb_derivs), deferred :: get_coulomb_derivs
   end type coulomb_charge_type

   abstract interface
      !> Evaluate coulomb matrix
      subroutine get_coulomb_matrix(self, mol, cache, amat)
         import :: coulomb_charge_type, structure_type, coulomb_cache, wp
         !> Instance of the electrostatic container
         class(coulomb_charge_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container
         type(coulomb_cache), intent(inout) :: cache
         !> Coulomb matrix
         real(wp), contiguous, intent(out) :: amat(:, :)
      end subroutine get_coulomb_matrix

      !> Evaluate uncontracted derivatives of Coulomb matrix
      subroutine get_coulomb_derivs(self, mol, cache, qat, qsh, dadr, dadL, &
         & atrace, dadcn, datracedcn)
         import :: coulomb_charge_type, structure_type, coulomb_cache, wp
         !> Instance of the electrostatic container
         class(coulomb_charge_type), intent(in) :: self
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Reusable data container
         type(coulomb_cache), intent(inout) :: cache
         !> Atomic partial charges
         real(wp), intent(in) :: qat(:)
         !> Shell-resolved partial charges
         real(wp), intent(in) :: qsh(:)
         !> Derivative of interactions with respect to cartesian displacements
         real(wp), contiguous, intent(out) :: dadr(:, :, :)
         !> Derivative of interactions with respect to strain deformations
         real(wp), contiguous, intent(out) :: dadL(:, :, :)
         !> On-site derivatives with respect to cartesian displacements
         real(wp), contiguous, intent(out) :: atrace(:, :)
         !> Coulomb matrix derivative w.r.t. the coordination numbers
         real(wp), contiguous, intent(out) :: dadcn(:, :)
         !> On-site derivatives with respect to coordination numbers
         real(wp), contiguous, intent(out) :: datracedcn(:)
      end subroutine get_coulomb_derivs
   end interface

contains


!> Update container cache
subroutine update(self, mol, cache)
   !> Instance of the electrostatic container
   class(coulomb_charge_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(coulomb_cache), pointer :: ptr

   call taint(cache, ptr)
   call ptr%update(mol)

   if (.not.allocated(ptr%amat)) then
      allocate(ptr%amat(sum(self%nshell), sum(self%nshell)))
   end if
   call self%get_coulomb_matrix(mol, ptr, ptr%amat)

   if (.not.allocated(ptr%vvec)) then
      allocate(ptr%vvec(sum(self%nshell)))
   end if

end subroutine update


!> Evaluate selfconsistent energy of the interaction
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the electrostatic container
   class(coulomb_charge_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatic energy
   real(wp), intent(inout) :: energies(:)

   integer :: iat, jat, ii, jj, ish, jsh
   type(coulomb_cache), pointer :: ptr
   real(wp) :: tmp

   call view(cache, ptr)

   if(self%shell_resolved) then
      call symv(ptr%amat, wfn%qsh(:, 1), ptr%vvec, alpha=0.5_wp)
      do iat = 1, mol%nat
         ii = self%offset(iat)
         do ish = 1, self%nshell(iat)
            energies(iat) = energies(iat) + ptr%vvec(ii+ish) * wfn%qsh(ii+ish, 1)
         end do
      end do
      ! Offsite first/pseudo-second-order reference shell charge correction#
      if (allocated(self%refqsh)) then
         !$omp parallel do default(none) schedule(runtime) &
         !$omp shared(energies, self, mol, ptr, wfn) &
         !$omp private(iat, ii, ish, jat, jj, jsh, tmp) 
         do iat = 1, mol%nat
            ii = self%offset(iat)
            tmp = 0.0_wp 
            do ish = 1, self%nshell(iat)
               do jat = 1, iat-1
                  jj = self%offset(jat)
                  do jsh = 1, self%nshell(jat)
                     tmp = tmp + ptr%amat(jj+jsh, ii+ish) &
                        & * (self%refqsh(ii+ish) * wfn%qsh(jj+jsh, 1) &
                        & + self%refqsh(jj+jsh) * wfn%qsh(ii+ish, 1))
                  end do
               end do
            end do
            energies(iat) = energies(iat) + tmp
         end do
      end if
   else
      call symv(ptr%amat, wfn%qat(:, 1), ptr%vvec, alpha=0.5_wp)
      do iat = 1, mol%nat
         ii = self%offset(iat)
         energies(iat) = energies(iat) + ptr%vvec(iat) * wfn%qat(iat, 1)
      end do
   end if
end subroutine get_energy


!> Evaluate charge dependent potential shift from the interaction
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the electrostatic container
   class(coulomb_charge_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   integer :: iat, ii, ish, jsh
   type(coulomb_cache), pointer :: ptr
   real(wp) :: tmp

   call view(cache, ptr)

   if(self%shell_resolved) then
      call symv(ptr%amat, wfn%qsh(:, 1), pot%vsh(:, 1), beta=1.0_wp)
      ! Offsite first/pseudo-second-order reference shell charge correction
      if(allocated(self%refqsh)) then
         call symv(ptr%amat, self%refqsh, pot%vsh(:, 1), beta=1.0_wp)
         ! Remove onsite contributions
         do iat = 1, mol%nat
            ii = self%offset(iat)
            do ish = 1, self%nshell(iat)
               tmp = 0.0_wp
               do jsh = 1, self%nshell(iat)
                  tmp = tmp + ptr%amat(ii + ish, ii + jsh) * self%refqsh(ii + jsh)
               end do
               pot%vsh(ii + ish, 1) = pot%vsh(ii+ish, 1) - tmp
            end do
         end do
      end if
   else
      call symv(ptr%amat, wfn%qat(:, 1), pot%vat(:, 1), beta=1.0_wp)
   end if

end subroutine get_potential


!> Evaluate gradient of the charge dependent potential shift 
subroutine get_potential_gradient(self, mol, cache, wfn, pot)
   !> Instance of the electrostatic container
   class(coulomb_charge_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   integer :: ic, jc, iat, ndim, ii, ish
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), datr(:, :), &
      & dadcn(:, :), datrdcn(:), tmpdq(:)
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   ndim = sum(self%nshell)
   allocate(dadr(3, mol%nat, ndim), dadL(3, 3, ndim), datr(3, ndim), &
      & dadcn(ndim, ndim), datrdcn(ndim), tmpdq(ndim))

   ! Get Coulomb matrix derivative contracted with the atom-resolved charges
   call self%get_coulomb_derivs(mol, ptr, wfn%qat(:, 1), wfn%qsh(:, 1), &
      & dadr, dadL, datr, dadcn, datrdcn)

   if(self%shell_resolved) then
      ! Off-diagonal Coulomb matrix derivative
      pot%dvshdr(:, :, :, 1) = dadr
      
      do ic = 1, 3 
         do iat = 1, mol%nat
            ii = self%offset(iat)
            do ish = 1, self%nshell(iat)
               ! Diagonal Coulomb matrix derivative
               pot%dvshdr(ic, iat, ii+ish, 1) = - sum(dadr(ic, :, ii+ish))
            end do 
            ! Charge derivative
            tmpdq = wfn%dqshdr(ic, iat, :, 1)
            call symv(ptr%amat, tmpdq, pot%dvshdr(ic, iat, :, 1), beta=1.0_wp)
         end do
         
         ! Coulomb matrix derivative 
         pot%dvshdL(ic, :, :, 1) = pot%dvshdL(ic, :, :, 1) + dadL(ic, :, :)
         do jc = 1, 3
            ! Charge derivative
            tmpdq = wfn%dqshdL(ic, jc, :, 1)
            call symv(ptr%amat, tmpdq, pot%dvshdL(ic, jc, :, 1), beta=1.0_wp)
         end do
      end do

   else
      ! Off-diagonal Coulomb matrix derivative
      pot%dvatdr(:, :, :, 1) = dadr

      do ic = 1, 3 
         do iat = 1, mol%nat
            ! Diagonal Coulomb matrix derivative
            pot%dvatdr(ic, iat, iat, 1) = - sum(dadr(ic, :, iat))

            ! Charge derivative
            tmpdq = wfn%dqatdr(ic, iat, :, 1)
            call symv(ptr%amat, tmpdq, pot%dvatdr(ic, iat, :, 1), beta=1.0_wp)
         end do
         
         ! Coulomb matrix derivative 
         pot%dvatdL(ic, :, :, 1) = pot%dvatdL(ic, :, :, 1) + dadL(ic, :, :)
         do jc = 1, 3
            ! Charge derivative
            tmpdq = wfn%dqatdL(ic, jc, :, 1)
            call symv(ptr%amat, tmpdq, pot%dvatdL(ic, jc, :, 1), beta=1.0_wp)
         end do
      end do
   end if

end subroutine get_potential_gradient


!> Evaluate gradient contributions from the selfconsistent interaction
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the electrostatic container
   class(coulomb_charge_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the electrostatic energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the electrostatic energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: ndim, iat, ii, ish
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), datr(:, :)
   real(wp), allocatable :: dadcn(:, :), datrdcn(:), dedcn(:)
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   ndim = sum(self%nshell)
   allocate(dadr(3, mol%nat, ndim), dadL(3, 3, ndim), datr(3, ndim), &
      & dadcn(mol%nat, ndim), datrdcn(ndim))
   
   call self%get_coulomb_derivs(mol, ptr, wfn%qat(:, 1), wfn%qsh(:, 1), &
      & dadr, dadL, datr, dadcn, datrdcn)

   call gemv(dadr, wfn%qsh(:, 1), gradient, beta=1.0_wp)
   call gemv(dadL, wfn%qsh(:, 1), sigma, beta=1.0_wp, alpha=0.5_wp)

   if (allocated(self%refqsh)) then
      ! offsite first/pseudo-second-order reference shell charge correction
      call gemv(dadr, self%refqsh, gradient, beta=1.0_wp)
      ! Remove onsite contributions 
      do iat = 1, mol%nat
         ii = self%offset(iat)
         do ish = 1, self%nshell(iat)
            gradient(:, iat) = gradient(:, iat) + datr(:, ii+ish) * self%refqsh(ii+ish)
         end do  
      end do
      call gemv(dadL, self%refqsh, sigma, beta=1.0_wp, alpha=0.5_wp)
   end if

   ! possible derivatives due to CN dependent Hubbard parameters
   if (self%cn_dep) then
      allocate(dedcn(mol%nat))
      call gemv(dadcn, wfn%qsh(:, 1), dedcn)
      if (allocated(self%refqsh)) then
         ! offsite first/pseudo-second-order reference shell charge correction
         call gemv(dadcn, self%refqsh, dedcn, beta=1.0_wp)
         ! Remove onsite contributions 
         do iat = 1, mol%nat
            ii = self%offset(iat)
            do ish = 1, self%nshell(iat)
               dedcn(iat) = dedcn(iat) + datrdcn(ii+ish) * self%refqsh(ii+ish) 
            end do  
         end do
      end if
      call gemv(ptr%dcndr, dedcn, gradient, beta=1.0_wp)
      call gemv(ptr%dcndL, dedcn, sigma, beta=1.0_wp)
   end if

end subroutine get_gradient


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, atom_resolved, shell_resolved
   !> Instance of the electrostatic container
   class(coulomb_charge_type), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=merge(shell_resolved, atom_resolved, self%shell_resolved))
end function variable_info


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
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


end module tblite_coulomb_charge_type
