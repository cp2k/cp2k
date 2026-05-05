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

!> @file tblite/acp/type.f90
!> Provides the atomic correction potential for xTB.

!> Implementation of the atomic correction potential
module tblite_acp_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_acp_cache, only : acp_cache
   use tblite_adjlist, only : adjacency_list
   use tblite_basis_cache, only : basis_cache
   use tblite_basis_type, only : basis_type, new_basis, cgto_container
   use tblite_blas, only : gemm
   use tblite_integral_overlap, only : overlap_cgto, overlap_grad_cgto, &
      & maxl, msao, smap
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_acp, get_acp, get_acp_gradient

   type, public :: acp_type
      !> Auxiliary basis set for ACP projector
      class(basis_type), allocatable :: auxbas
      !> Energy levels of the auxiliary basis functions
      real(wp), allocatable :: levels(:, :)
   contains
      !> Update the ACP cache
      procedure :: update
   end type acp_type

contains


!> Factory for a new atomic correction potential object
subroutine new_acp(self, mol, nproj, cgtp, levels, accuracy)
   !> ACP object
   type(acp_type), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of projectors per species
   integer, intent(in) :: nproj(:)
   !> Contracte gaussian type projector functions for each projector and species
   type(cgto_container), allocatable, intent(in) :: cgtp(:, :)
   !> Energy levels of the ACP auxiliary basis functions
   real(wp), intent(in) :: levels(:, :)
   !> Optional accuracy specification
   real(wp), intent(in), optional :: accuracy

   allocate(self%auxbas)
   call new_basis(self%auxbas, mol, nproj, cgtp, accuracy=accuracy)

   self%levels = levels

end subroutine new_acp


! Update ACP cache
subroutine update(self, mol, cache)
   !> Instance of the basis type
   class(acp_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(acp_cache), intent(inout) :: cache
   call self%auxbas%update(mol, cache%auxbas, .false.)

end subroutine update


subroutine get_acp(mol, trans, list, bas, bcache, acp, acache, hamiltonian)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(inout) :: acache
   !> Effective Hamiltonian
   real(wp), intent(inout) :: hamiltonian(:, :)

   real(wp), allocatable :: pv_overlap(:, :)

   if (.not. allocated(acache%scaled_pv_overlap)) then
      allocate(acache%scaled_pv_overlap(acp%auxbas%nao, bas%nao), source=0.0_wp)
   end if
   allocate(pv_overlap(acp%auxbas%nao, bas%nao), source=0.0_wp)

   ! Obtain the (scaled) projector-valence overlap matrix
   if (any(mol%periodic)) then
      call get_pv_overlap_3d(mol, trans, list, bas, bcache, acp, acache, &
         & pv_overlap)
   else
      call get_pv_overlap_0d(mol, bas, bcache, acp, acache, pv_overlap)
   end if

   ! Contract over the auxiliary projectors and add the ACP to the Hamiltonian
   call gemm(amat=pv_overlap, bmat=acache%scaled_pv_overlap, cmat=hamiltonian, &
      & transa='T', beta=1.0_wp)

end subroutine get_acp


subroutine get_pv_overlap_0d(mol, bas, bcache, acp, acache, pv_overlap)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(inout) :: acache
   !> Projector-valence overlap matrix
   real(wp), intent(out) :: pv_overlap(:, :)

   integer :: iat, jat, isp, jsp, itr, k, img, inl
   integer :: ish, jproj, is, js, ii, jj, iao, jao, nao, ij, iaosh, jaoproj
   real(wp) :: r2, vec(3), level
   real(wp), allocatable :: stmp(:)

   allocate(stmp(msao(bas%maxl)*msao(acp%auxbas%maxl)))
   
   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, bcache, acp, acache, pv_overlap) &
   !$omp private(iat, jat, isp, jsp, itr, inl, img, is, js, ish, jproj, ii, jj) &
   !$omp private(iao, jao, iaosh, jaoproj, nao, ij, r2, vec, stmp, level)
   do iat = 1, mol%nat
      isp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         jsp = mol%id(jat)
         js = acp%auxbas%ish_at(jat)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         ! Loop over valence shells
         do ish = 1, bas%nsh_id(isp)
            ii = bas%iao_sh(is+ish)
            ! Loop over projectors
            do jproj = 1, acp%auxbas%nsh_id(jsp)
               jj = acp%auxbas%iao_sh(js+jproj)

               call overlap_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                  & bas%cgto(ish, isp)%raw, acache%auxbas%cgto(jproj, jat), &
                  & bcache%cgto(ish, iat), r2, vec, acp%auxbas%intcut, stmp)

               level = acp%levels(jproj, jsp)
               
               nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
               do iao = 1, msao(bas%cgto(ish, isp)%raw%ang)
                  do jao = 1, nao
                     ij = jao + nao*(iao-1)

                     ! Store scaled overlap for later gradient calculation
                     acache%scaled_pv_overlap(jj+jao, ii+iao) = level * stmp(ij)

                     pv_overlap(jj+jao, ii+iao) = stmp(ij)
                  end do
               end do

            end do
         end do

      end do
   end do

end subroutine get_pv_overlap_0d


subroutine get_pv_overlap_3d(mol, trans, list, bas, bcache, acp, acache, &
   & pv_overlap)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(inout) :: acache
   !> Projector-valence overlap matrix
   real(wp), intent(out) :: pv_overlap(:, :)

   pv_overlap = 0.0_wp
   acache%scaled_pv_overlap = 0.0_wp 

end subroutine get_pv_overlap_3d


subroutine get_acp_gradient(mol, trans, list, bas, bcache, acp, acache, wfn, &
   & dEdcnbas, dEdqbas, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(in) :: acache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Derivative of the electronic energy w.r.t. basis set coordination numbers
   real(wp), intent(inout) :: dEdcnbas(:)
   !> Derivative of the electronic energy w.r.t. basis set atomic charges
   real(wp), intent(inout) :: dEdqbas(:)
   !> Derivative of the electronic energy w.r.t. coordinate displacements
   real(wp), intent(inout) :: gradient(:, :)
   !> Derivative of the electronic energy w.r.t. strain deformations
   real(wp), intent(inout) :: sigma(:, :)

   if (any(mol%periodic)) then
      call get_pv_overlap_deriv_3d(mol, trans, list, bas, bcache, acp, acache, &
         & wfn%density, dEdcnbas, dEdqbas, gradient, sigma)
   else
      call get_pv_overlap_deriv_0d(mol, bas, bcache, acp, acache, wfn%density, &
         & dEdcnbas, dEdqbas, gradient, sigma)
   end if

end subroutine get_acp_gradient


subroutine get_pv_overlap_deriv_0d(mol, bas, bcache, acp, acache, pmat, &
   & dEdcnbas, dEdqbas, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(in) :: acache
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :, :)
   !> Derivative of the electronic energy w.r.t. basis set coordination numbers
   real(wp), intent(inout) :: dEdcnbas(:)
   !> Derivative of the electronic energy w.r.t. basis set atomic charges
   real(wp), intent(inout) :: dEdqbas(:)
   !> Derivative of the electronic energy w.r.t. coordinate displacements
   real(wp), intent(inout) :: gradient(:, :)
   !> Derivative of the electronic energy w.r.t. strain deformations
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, isp, jsp, spin, nspin
   integer :: ish, jproj, is, js, ii, jj, iao, jao, nao, ij
   real(wp) :: r2, vec(3), dG(3), dcnbasi, dqbasi, tmp
   real(wp), allocatable :: stmp(:), dstmp(:, :), spmat(:, :)
   real(wp), allocatable :: dstmpdqeffi(:), dstmpdqeffj(:)
   logical :: compute_qeff_grad

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: dEdcnbas_local(:), dEdqbas_local(:), &
      & gradient_local(:, :), sigma_local(:, :)

   ! Determine before the loop if effective charge derivatives are needed
   compute_qeff_grad = bas%charge_dependent .or. acp%auxbas%charge_dependent

   nspin = size(pmat, 3)

   allocate(spmat(acp%auxbas%nao, bas%nao), stmp(msao(bas%maxl)*msao(acp%auxbas%maxl)), &
      & dstmp(3, msao(bas%maxl)*msao(acp%auxbas%maxl)))
   
   if (compute_qeff_grad) then
      allocate(dstmpdqeffi(msao(bas%maxl)*msao(acp%auxbas%maxl)), &
         & dstmpdqeffj(msao(bas%maxl)*msao(acp%auxbas%maxl)))
   end if

   do spin = 1, nspin
      ! Precalculate scaled projector-density matrix product
      call gemm(amat=acache%scaled_pv_overlap, bmat=pmat(:, :, spin), &
         & cmat=spmat(:, :), beta=0.0_wp)

      !$omp parallel default(none) shared(dEdcnbas, dEdqbas, gradient, sigma) &
      !$omp shared(mol, bas, bcache, acp, acache, spmat, compute_qeff_grad) &
      !$omp private(iat, jat, isp, jsp, is, js, ish, jproj, ii, jj, iao, jao, nao, ij) &
      !$omp private(r2, vec, stmp, dstmp, dstmpdqeffi, dstmpdqeffj, dcnbasi, dqbasi) &
      !$omp private(tmp, dG, dEdcnbas_local, dEdqbas_local, gradient_local, sigma_local)
      if (compute_qeff_grad) then
         allocate(dEdcnbas_local(size(dEdcnbas)), source=0.0_wp)
         allocate(dEdqbas_local(size(dEdqbas)), source=0.0_wp)
      end if
      allocate(gradient_local(size(gradient,1), size(gradient,2)), source=0.0_wp)
      allocate(sigma_local(size(sigma,1), size(sigma,2)), source=0.0_wp)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         isp = mol%id(iat)
         is = bas%ish_at(iat)
         dcnbasi = 0.0_wp
         dqbasi = 0.0_wp
         do jat = 1, mol%nat
            !if (iat == jat) cycle
            jsp = mol%id(jat)
            js = acp%auxbas%ish_at(jat)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            ! Loop over valence shells
            dG(:) = 0.0_wp
            do ish = 1, bas%nsh_id(isp)
               ii = bas%iao_sh(is+ish)
               ! Loop over projectors
               do jproj = 1, acp%auxbas%nsh_id(jsp)
                  jj = acp%auxbas%iao_sh(js+jproj)

                  ! Calculate overlap integral derivatives
                  call overlap_grad_cgto(acp%auxbas%cgto(jproj, jsp)%raw, &
                     & bas%cgto(ish, isp)%raw, acache%auxbas%cgto(jproj, jat), &
                     & bcache%cgto(ish, iat), r2, vec, acp%auxbas%intcut, &
                     & stmp, dstmp, dstmpdqeffj, dstmpdqeffi)

                  ! Contract derivatives with the SxP intermediate
                  nao = msao(acp%auxbas%cgto(jproj, jsp)%raw%ang)
                  do iao = 1, msao(bas%cgto(ish, isp)%raw%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)
                        dG(:) = dG + 2.0_wp * spmat(jj+jao, ii+iao) * dstmp(:, ij)
                        if (compute_qeff_grad) then
                           tmp = 2.0_wp * spmat(jj+jao, ii+iao) * dstmpdqeffi(ij)
                           dcnbasi = dcnbasi + bcache%cgto(ish, iat)%dqeffdcn * tmp
                           dqbasi = dqbasi + bcache%cgto(ish, iat)%dqeffdq * tmp
                        end if
                     end do
                  end do

               end do
            end do
            ! Collect gradient and sigma once per atom pair
            gradient_local(:, iat) = gradient_local(:, iat) + dG
            gradient_local(:, jat) = gradient_local(:, jat) - dG
            sigma_local(:, :) = sigma_local + 0.5_wp * (spread(vec, 1, 3) &
               & * spread(dG, 2, 3) + spread(dG, 1, 3) * spread(vec, 2, 3))
         end do
         ! Collect effective charge derivatives once per atom
         if (compute_qeff_grad) then
            dEdcnbas_local(iat) = dEdcnbas_local(iat) + dcnbasi
            dEdqbas_local(iat) = dEdqbas_local(iat) + dqbasi
         end if
      end do
      !$omp end do
      !$omp critical
      if (compute_qeff_grad) then
         dEdcnbas(:) = dEdcnbas + dEdcnbas_local
         dEdqbas(:) = dEdqbas + dEdqbas_local
      end if
      gradient(:, :) = gradient + gradient_local
      sigma(:, :) = sigma + sigma_local
      !$omp end critical
      if (compute_qeff_grad) then
         deallocate(dEdcnbas_local, dEdqbas_local)
      end if
      deallocate(gradient_local, sigma_local)
      !$omp end parallel
   end do

end subroutine get_pv_overlap_deriv_0d


subroutine get_pv_overlap_deriv_3d(mol, trans, list, bas, bcache, acp, acache, &
   & pmat, dEdcnbas, dEdqbas, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Atomic correction potential data
   type(acp_type), intent(in) :: acp
   !> Atomic correction potential cache
   type(acp_cache), intent(in) :: acache
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :, :)
   !> Derivative of the electronic energy w.r.t. basis set coordination numbers
   real(wp), intent(inout) :: dEdcnbas(:)
   !> Derivative of the electronic energy w.r.t. basis set atomic charges
   real(wp), intent(inout) :: dEdqbas(:)
   !> Derivative of the ACP energy w.r.t. coordinate displacements
   real(wp), intent(inout) :: gradient(:, :)
   !> Derivative of the ACP energy w.r.t. strain deformations
   real(wp), intent(inout) :: sigma(:, :)


end subroutine get_pv_overlap_deriv_3d


end module tblite_acp_type
