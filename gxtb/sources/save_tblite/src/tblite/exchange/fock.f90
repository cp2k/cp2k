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

!> @file tblite/exchange/fock.f90
!> Provides matrix multiplication algorithm for approximated Fock exchange

!> Approximated Fock exchange based on matrix multiplication algorithm
module tblite_exchange_fock
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_blas, only: gemm, symm, axpy
   use tblite_basis_type, only : basis_type
   use tblite_container_cache, only : container_cache
   use tblite_exchange_cache, only : exchange_cache
   use tblite_exchange_type, only : exchange_type
   use tblite_utils_average, only : average_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_wignerseitz, only : wignerseitz_cell
   implicit none
   private

   public :: new_exchange_fock

   type, public, extends(exchange_type) :: exchange_fock
      !> Averaged Hubbard parameter for each shell and species
      real(wp), allocatable :: hubbard(:, :, :, :)
      !> One center exchange integrals
      real(wp), allocatable :: onecxints(:, :, :)
      !> Diagonal scaling of the Fock exchange
      real(wp) :: ondiag_scale
      !> Off-diagonal scaling of the Fock exchange
      real(wp), allocatable :: offdiag_scale(:, :, :, :)
      !> Exponent of radius dependent hubbard scaling
      real(wp) :: hubbard_exp
      !> Radius prefactor of radius dependent hubbard scaling
      real(wp) :: hubbard_exp_r0
      !> Smoothening exponent (1 = Mataga-Nishimoto, 2 = Klopman-Ohno)
      real(wp) :: gexp
      !> Pairwise radii for approximate exchange integrals
      real(wp), allocatable :: rad(:, :)
      !> Charge-dependence of the onsite Fock exchange
      real(wp), allocatable :: kq(:, :)
      !> Bond-order correlation scaling factor for each atom pair
      real(wp), allocatable :: corr_scale(:, :)
      !> Bond-order correlation damping exponent
      real(wp) :: corr_exp
      !> Bond-order correlation radius for each atom pair
      real(wp), allocatable :: corr_rad(:, :)
      !> Long-range cutoff
      real(wp) :: rcut
   contains
      !> Evaluate Mulliken Fock exchange gamma matrix
      procedure :: get_mulliken_Kmatrix
      !> Evaluate onsite Fock exchange gamma matrix
      procedure :: get_onsite_Kmatrix
      !> Evaluate bond-order correlation correction gamma matrix
      procedure :: get_bocorr_Kmatrix
      !> Evaluate the gradient of the Mulliken exchange energy
      procedure :: get_mulliken_derivs
      !> Evaluate the gradient of the bond-order correlation correction energy
      procedure :: get_bocorr_derivs
      !> Calculate exchange contribution to the Fock matrix
      procedure :: get_KFock
      !> Calculate exchange contribution to the gradient
      procedure :: get_KGrad
   end type exchange_fock

   real(wp), parameter :: sqrtpi = sqrt(pi)
   character(len=*), parameter :: label = "Mulliken and onsite Fock exchange + bond-order correlation"

contains


!> Create a new approximate exchange container
subroutine new_exchange_fock(self, mol, bas, hubbard, hubbard_average, &
   & avg_exponents, ondiag_scale, offdiag_scale, offdiag_average, hubbard_exp, &
   & hubbard_exp_r0, rad, gexp, onecxints, kq, corr_scale, corr_scale_average, &
   & corr_exp, corr_rad, corr_rad_average, frscale, omega, lrscale)
   !> Instance of the Fock exchange container
   type(exchange_fock), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Description of the basis set
   type(basis_type), intent(in) :: bas
   !> Hubbard parameter for all shells and species
   real(wp), intent(in) :: hubbard(:, :)
   !> Averaging function for Hubbard parameter of a shell-pair
   type(average_type), intent(in) :: hubbard_average
   !> Averaging exponents for all shells and species
   real(wp), intent(in) :: avg_exponents(:, :)
   !> Diagonal scaling of the Fock exchange
   real(wp), intent(in) :: ondiag_scale
   !> Off-diagonal scaling of the Fock exchange
   real(wp), intent(in) :: offdiag_scale(:, :)
   !> Averaging function for the off-diagonal scaling of a shell-pair
   type(average_type), intent(in) :: offdiag_average
   !> Exponent of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp
   !> Radius prefactor of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp_r0
   !> Radius for hubbard scaling
   real(wp), intent(in) :: rad(:, :)
   !> Smoothening exponent; 1 = Mataga-Nishimoto, 2 = Klopman-Ohno
   real(wp), intent(in) :: gexp
   !> One center exchange integrals
   real(wp), intent(in) :: onecxints(:, :, :)
   !> Charge dependence of the onsite Fock exchange 
   real(wp), intent(in) :: kq(:, :)
   !> Bond-order correlation scaling factor for each atoa
   real(wp), intent(in) :: corr_scale(:)
   !> Averaging function for correlation scaling of an atom-pair
   type(average_type), intent(in) :: corr_scale_average
   !> Bond-order correlation damping exponent
   real(wp), intent(in) :: corr_exp
   !> Bond-order correlation radius for each atom
   real(wp), intent(in) :: corr_rad(:)
   !> Averaging function for correlation radius of an atom-pair
   type(average_type), intent(in) :: corr_rad_average
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Optional range separation parameter
   real(wp), intent(in), optional :: omega
   !> Optional long range scaling for range seperated exchange
   real(wp), intent(in), optional :: lrscale

   integer :: isp, jsp, ish, jsh
   real(wp) :: inter_avg

   self%label = label

   allocate(self%nsh_id(mol%nid), self%nao_sh(bas%nsh), self%ish_at(mol%nat), &
      & self%iao_sh(bas%nsh))
   self%nao = bas%nao
   self%nsh = bas%nsh
   self%nsh_id = bas%nsh_id
   self%nao_sh = bas%nao_sh
   self%ish_at = bas%ish_at
   self%iao_sh = bas%iao_sh
   self%maxsh = maxval(bas%nsh_id)

   self%frscale = frscale
   if (present(omega).and.present(lrscale)) then
      self%omega = omega
      self%lrscale = lrscale
   end if

   allocate(self%offdiag_scale(self%maxsh, self%maxsh, mol%nid, mol%nid))
   self%offdiag_scale = 0.0_wp
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         do ish = 1, bas%nsh_id(isp)
            do jsh = 1, bas%nsh_id(jsp)
               self%offdiag_scale(jsh, ish, jsp, isp) = offdiag_average%value( &
                  & offdiag_scale(ish, isp), offdiag_scale(jsh, jsp))
            end do
         end do
      end do
   end do
   self%ondiag_scale  = ondiag_scale

   self%gexp = gexp
   self%rcut = 10.0_wp

   self%hubbard_exp = hubbard_exp
   self%hubbard_exp_r0 = hubbard_exp_r0
   self%rad = rad

   allocate(self%hubbard(self%maxsh, self%maxsh, mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         self%hubbard(:, :, jsp, isp) = 0.0_wp
         do ish = 1, bas%nsh_id(isp)
            do jsh = 1, bas%nsh_id(jsp)
               inter_avg = max(avg_exponents(ish, isp), avg_exponents(jsh, jsp))
               self%hubbard(jsh, ish, jsp, isp) = hubbard_average%value(&
                  & hubbard(ish, isp), hubbard(jsh, jsp), inter_avg)
            end do
         end do
      end do
   end do

   self%kq = kq
   self%onecxints = onecxints

   allocate(self%corr_scale(mol%nid, mol%nid), self%corr_rad(mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         self%corr_scale(isp, jsp) = &
            & corr_scale_average%value(corr_scale(isp), corr_scale(jsp))
         self%corr_rad(isp, jsp) = &
            & corr_rad_average%value(corr_rad(isp), corr_rad(jsp))
         end do 
   end do 
   self%corr_exp = corr_exp

end subroutine new_exchange_fock


!> Calculate exchange contribution to the Fock matrix and atomic potential
subroutine get_KFock(self, mol, cache, density, overlap)
   !> Instance of the Fock exchange container
   class(exchange_fock), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container with intermediates and the final Fock matrix/potential
   type(exchange_cache), intent(inout) :: cache
   !> Density matrix
   real(wp), intent(in) :: density(:, :, :)
   !> Overlap matrix
   real(wp), intent(in) :: overlap(:, :)

   integer :: spin, iat, izp, ish, is, ii, iao, jsh, jj, jao
   real(wp) :: gfx, gri, tmp
   real(wp), allocatable :: tmpA(:,:), tmpB(:,:)
   real(wp), allocatable :: diagP(:), gdiagP(:)
   real(wp), allocatable :: tmpSP(:), gdiagSP(:)
   real(wp), allocatable :: tmpSPS(:), gdiagSPS(:)

   real(wp) :: spin_factor

   allocate(tmpA(self%nao, self%nao), tmpB(self%nao, self%nao), diagP(self%nao), &
      & gdiagP(self%nao), tmpSP(self%nao), gdiagSP(self%nao), tmpSPS(self%nao), &
      & gdiagSPS(self%nao), source = 0.0_wp)

   ! Select spin factor to cancel the quadratic dependence of the exchange energy
   ! on the occupation numbers (0.5 for restricted, and 1.0 for unrestricted)
   spin_factor = 0.5_wp
   if(size(density, 3) .gt. 1) then
      spin_factor = 1.0_wp
   end if

   cache%prev_F(:, :, :) = 0.0_wp
   cache%prev_vsh(:, :) = 0.0_wp

   ! Evaluate the Fock matrix contribution for Mulliken
   ! and onsite approximated Fock exchange for a symmetric density matrix
   do spin = 1, size(density, 3)

      ! Intermediate A = S x P
      call symm(amat=overlap, bmat=density(:, :, spin), cmat=tmpA)

      ! Collect P diagonal onsite correction
      do iao = 1, self%nao
         diagP(iao) = density(iao, iao, spin)
      end do
      call onsite_fx_symv(mol%nat, mol%id, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, cache%g_onsfx, diagP, gdiagP)

      ! Collect S x P diagonal onsite correction
      do iao = 1, self%nao
         tmpSP(iao) = tmpA(iao, iao)
      end do
      call onsite_fx_symv(mol%nat, mol%id, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, cache%g_onsfx, tmpSP, gdiagSP)

      ! Collect S * P summation
      tmpSP = 0.0_wp
      tmpSP(:) = sum(density(:, :, spin) * overlap, dim=2)

      ! Collect S x P x S diagonal onsite correction
      tmpB(:, :) = tmpA * overlap
      tmpSPS(:) = sum(tmpB, dim=2)
      call onsite_fx_symv(mol%nat, mol%id, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, cache%g_onsfx, tmpSPS, gdiagSPS)

      ! Intermediate F = 1/2 * (S x P) x S
      call gemm(amat=tmpA, bmat=overlap, alpha = 0.5_wp, &
         & cmat=cache%prev_F(:, :, spin))

      ! Shell potential from self-consistent onsite term
      !$omp parallel do default(none) schedule(runtime) shared(mol, self, cache) &
      !$omp shared(density, tmpA, diagP, tmpSP, tmpSPS, spin, spin_factor) &
      !$omp private(iat, izp, is, ish, jsh, ii, jj, iao, jao, gfx, gri)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is  = self%ish_at(iat)
         do ish = 1, self%nsh_id(izp)
            ii = self%iao_sh(is+ish)
            do jsh = 1, self%nsh_id(izp)
               jj = self%iao_sh(is+jsh)

               gfx = 0.0_wp
               gri = 0.0_wp
               do iao = 1, self%nao_sh(is+ish)
                  do jao = 1, self%nao_sh(is+jsh)
                     gfx = gfx - spin_factor * 0.25_wp * ( &
                        & + (density(jj+jao, ii+iao, spin) * cache%prev_F(jj+jao, ii+iao, spin)) &
                        & + 0.5_wp * (tmpA(ii+iao, jj+jao) * tmpA(ii+iao, jj+jao)) &
                        & + 0.25_wp * (tmpSPS(jj+jao) * diagP(ii+iao)) &
                        & + 0.5_wp * (tmpSP(jj+jao) * tmpSP(ii+iao)) &
                        & + 0.25_wp * (diagP(jj+jao) * tmpSPS(ii+iao)) )

                     if (ish == jsh) then
                        gri = gri + spin_factor * ( &
                           & + (density(jj+jao, ii+iao, spin) * cache%prev_F(jj+jao, ii+iao, spin)) &
                           & + 0.5_wp * (tmpA(ii+iao, jj+jao) * tmpA(jj+jao, ii+iao)) )
                     end if
                  end do
               end do

               if (ish == jsh) then
                  cache%prev_vsh(is+ish, 1) = cache%prev_vsh(is+ish, 1) &
                     & + cache%dgdq_onsri(ish, is+ish) * gri & 
                     & + cache%dgdq_onsfx(jsh, ish, is+ish) * gfx
               else
                  cache%prev_vsh(is+ish, 1) = cache%prev_vsh(is+ish, 1) + &
                     & cache%dgdq_onsfx(jsh, ish, is+ish) * gfx
                  cache%prev_vsh(is+jsh, 1) = cache%prev_vsh(is+jsh, 1) + &
                     & cache%dgdq_onsfx(jsh, ish, is+jsh) * gfx
               end if
            end do
         end do
      end do

      ! Apply Mulliken, onsite, and bond-order correlation matrices as g * (S x P)
      tmpB = 0.0_wp
      call shell_hadamard_add(self%nsh, self%nao_sh, self%iao_sh, cache%g_mulliken, &
         & tmpA, 1.0_wp, tmpB)
      call atom_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, cache%g_bocorr, tmpA, -4.0_wp, tmpB)
      call onsite_fx_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsfx, tmpA, 0.5_wp, tmpB, trans_src=.true.)
      call onsite_ri_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsri, tmpA, -2.0_wp, tmpB)
      tmpA = tmpB

      ! Apply Mulliken, onsite, and bond-order correlation matrices as g * (S x P x S)
      tmpB = 0.0_wp
      call onsite_fx_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsfx, cache%prev_F(:, :, spin), &
         & 0.5_wp, tmpB)
      call onsite_ri_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsri, cache%prev_F(:, :, spin), &
         & -2.0_wp, tmpB)
      call shell_hadamard_add(self%nsh, self%nao_sh, self%iao_sh, cache%g_mulliken, &
         & cache%prev_F(:, :, spin), 1.0_wp, tmpB)
      cache%prev_F(:, :, spin) = tmpB

      ! Apply Mulliken, onsite, and bond-order correlation matrices as g * P
      tmpB = 0.0_wp
      call shell_hadamard_add(self%nsh, self%nao_sh, self%iao_sh, cache%g_mulliken, &
         & density(:, :, spin), 1.0_wp, tmpB)
      call onsite_fx_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsfx, density(:, :, spin), &
         & 0.5_wp, tmpB)
      call onsite_ri_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsri, density(:, :, spin), &
         & -2.0_wp, tmpB)

      ! Add P diagonal onsite correction
      do iao = 1, self%nao
         call axpy(xvec=overlap(:, iao), yvec=tmpA(:, iao), &
            & alpha=0.25_wp*gdiagP(iao))
      end do 

      ! Intermediate A += 1/2 * S x (g * P)
      call symm(amat=overlap, bmat=tmpB, alpha=0.5_wp, &
         & cmat=tmpA, beta=1.0_wp)

      ! Add intermediate F += (S x X) x S
      call gemm(amat=tmpA, bmat=overlap, cmat=cache%prev_F(:, :, spin), &
         & beta = 1.0_wp)

      ! Add S x P and S x P x S diagonal onsite corrections
      do iao = 1, self%nao
         call axpy(xvec=overlap(:, iao), yvec=cache%prev_F(:, iao, spin), &
            & alpha=0.5_wp*gdiagSP(iao))
         cache%prev_F(iao, iao, spin) = cache%prev_F(iao, iao, spin) &
            & + 0.25_wp * gdiagSPS(iao)
      end do

      ! Save symmetrized Fock matrix for energy evaluation
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(self, cache, spin_factor, spin) &
      !$omp private(ii, jj, tmp)
      do ii = 1, self%nao
         cache%prev_F(ii, ii, spin) = -0.5_wp * spin_factor &
            & * cache%prev_F(ii, ii, spin)
         do jj = 1, ii-1
            tmp = -0.25_wp * spin_factor * &
               & (cache%prev_F(jj, ii, spin) + cache%prev_F(ii, jj, spin))
            cache%prev_F(jj, ii, spin) = tmp
            cache%prev_F(ii, jj, spin) = tmp
         end do
      end do
   end do

end subroutine get_KFock


!> Calculate exchange contribution to the gradient
subroutine get_KGrad(self, mol, cache, density, overlap, mulliken_grad, &
   & bocorr_grad, ao_grad)
   !> Instance of the exchange container
   class(exchange_fock), intent(in) :: self
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

   integer :: spin, iao
   real(wp), allocatable :: tmpA(:, :), tmpB(:, :), tmpC(:, :), tmpD(:, :)
   real(wp), allocatable :: tmpVec(:), diagP(:), diagSP(:), tmpSPSvec(:)
   real(wp) :: spin_factor

   mulliken_grad = 0.0_wp
   bocorr_grad = 0.0_wp

   allocate(tmpA(self%nao, self%nao), tmpB(self%nao, self%nao), &
      & tmpC(self%nao, self%nao), tmpD(self%nao, self%nao), &
      & tmpVec(self%nao), tmpSPSvec(self%nao), diagP(self%nao), &
      & diagSP(self%nao), source = 0.0_wp)

   ! Select spin factor to cancel the quadratic dependence of the exchange energy
   ! on the occupation numbers (0.5 for restricted, and 1.0 for unrestricted)
   spin_factor = 0.5_wp
   if(size(density, 3) .gt. 1) then
      spin_factor = 1.0_wp
   end if

   ! Evaluate the operator and overlap derivative matrices for Mulliken
   ! and onsite approximated Fock exchange for a symmetric density matrix
   do spin = 1, size(density, 3)
      
      ! Intermediate A = S x P
      call symm(amat=overlap, bmat=density(:, :, spin), cmat=tmpA)

      ! Collect S x P diagonal onsite correction
      tmpVec(:) = 0.0_wp
      do iao = 1, self%nao
         tmpVec(iao) = tmpA(iao, iao)
      end do
      call onsite_fx_symv(mol%nat, mol%id, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, cache%g_onsfx, tmpVec, diagSP)

      ! Collect P diagonal onsite correction
      tmpVec(:) = 0.0_wp
      do iao = 1, self%nao
         tmpVec(iao) = density(iao, iao, spin)
      end do
      call onsite_fx_symv(mol%nat, mol%id, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, cache%g_onsfx, tmpVec, diagP)

      ! Out-of-place transpose of A for onsite correction
      tmpC = transpose(tmpA)

      ! Operator derivatives (P x S) * (S x P)
      mulliken_grad = mulliken_grad + 0.5_wp * spin_factor * tmpA * tmpC
      bocorr_grad = bocorr_grad + 2.0_wp * spin_factor * tmpA * tmpC

      ! Apply Mulliken, onsite, bond-order correlation matrices
      ! as g * (S x P) and g * (P x S)
      tmpB = 0.0_wp
      call shell_hadamard_add(self%nsh, self%nao_sh, self%iao_sh, cache%g_mulliken, &
         & tmpA, 1.0_wp, tmpB)
      call atom_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, cache%g_bocorr, tmpA, -4.0_wp, tmpB)
      call onsite_fx_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsfx, tmpC, 0.5_wp, tmpB)
      call onsite_ri_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsri, tmpA, -2.0_wp, tmpB)

      ! Add diagP * S x P and (P x S) * diagP onsite corrections
      tmpD = 0.0_wp
      do iao = 1, self%nao
         call axpy(xvec=density(:, iao, spin), yvec=tmpD(:, iao), &
            & alpha=0.5_wp * diagSP(iao))
         call axpy(xvec=tmpC(:, iao), yvec=tmpD(:, iao), &
            & alpha=0.5_wp * diagP(iao))
      end do

      ! Apply Mulliken and onsite matrices as g * P
      tmpC = 0.0_wp
      call shell_hadamard_add(self%nsh, self%nao_sh, self%iao_sh, cache%g_mulliken, &
         & density(:, :, spin), 1.0_wp, tmpC)
      call onsite_fx_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsfx, density(:, :, spin), 0.5_wp, tmpC)
      call onsite_ri_hadamard_add(mol%nat, mol%id, self%nsh_id, self%nao_sh, &
         & self%ish_at, self%iao_sh, cache%g_onsri, density(:, :, spin), -2.0_wp, tmpC)

      ! Add intermediate B += S x (g * P)
      call symm(amat=overlap, bmat=tmpC, cmat=tmpB, beta=1.0_wp)

      ! Add intermediate D += P x (S x X)
      call symm(amat=density(:, :, spin), bmat=tmpB, cmat=tmpD, beta=1.0_wp)

      ! Out-of-place transpose of D for gradient contribution
      tmpB = transpose(tmpD)

      ! Symmetrized overlap derivative contribution
      ! 0.25 arises since we need to treat the upper and lower triangular parts
      ! as indpendent variables to get the correct energy weighed density matrix
      ao_grad = ao_grad + 0.25_wp * spin_factor * (tmpD + tmpB)

      ! Intermediate B = (S x P) x S
      call gemm(amat=tmpA, bmat=overlap, cmat=tmpB)

      ! Operator derivatives P * (S x P) x S
      mulliken_grad = mulliken_grad + 0.5_wp * spin_factor * density(:, :, spin) * tmpB

      ! Collect S x P x S diagonal vector
      tmpSPSvec(:) = 0.0_wp
      do iao = 1, self%nao
        tmpSPSvec(iao) = tmpB(iao, iao)
      end do
   end do

end subroutine get_KGrad


!> Evaluate Mulliken exchange gamma matrix
subroutine get_mulliken_Kmatrix(self, mol, cache)
   !> Instance of the exchange container
   class(exchange_fock), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container with intermediates and the final exchange matrix
   type(exchange_cache), intent(inout) :: cache

   if (any(mol%periodic)) then
      call get_gmulliken_3d(mol, self%nsh_id, self%ish_at, self%hubbard, &
         & self%ondiag_scale, self%offdiag_scale, self%hubbard_exp, &
         & self%hubbard_exp_r0, self%rad, self%gexp, self%frscale, self%omega, &
         & self%lrscale, self%rcut, cache%wsc, cache%alpha, cache%g_mulliken)
   else
      call get_gmulliken_0d(mol, self%nsh_id, self%ish_at, self%hubbard, &
         & self%ondiag_scale, self%offdiag_scale, self%hubbard_exp, &
         & self%hubbard_exp_r0, self%rad, self%gexp, self%frscale, self%omega, &
         & self%lrscale, cache%g_mulliken)
   end if

end subroutine get_mulliken_Kmatrix


!> Evaluate range separated exchange matrix for finite systems
subroutine get_gmulliken_0d(mol, nsh_id, ish_at, hubbard, ondiag_scale, &
   & offdiag_scale, hubbard_exp, hubbard_exp_r0, rad, gexp, frscale, &
   & omega, lrscale, g_mulliken)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Hubbard parameter parameter for each shell
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Diagonal scaling of the Fock exchange
   real(wp), intent(in) :: ondiag_scale
   !> Off-diagonal scaling of the Fock exchange
   real(wp), intent(in) :: offdiag_scale(:, :, :, :)
   !> Exponent of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp
   !> Radius prefactor of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp_r0
   !> Radius for hubbard scaling
   real(wp), intent(in) :: rad(:, :)
   !> Exponent of exchange kernel
   real(wp), intent(in) :: gexp
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Range separation parameter
   real(wp), intent(in) :: omega
   !> Long-range scaling factor
   real(wp), intent(in) :: lrscale
   !> Mulliken exchange matrix
   real(wp), intent(out) :: g_mulliken(:, :)

   integer :: iat, jat, izp, jzp, is, js, ish, jsh
   real(wp) :: vec(3), r1, r1g, gam, denom, rsh, scale, damp

   g_mulliken(:, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) shared(g_mulliken) &
   !$omp shared(mol, nsh_id, ish_at, hubbard, gexp, hubbard_exp, hubbard_exp_r0) &
   !$omp shared(rad, offdiag_scale, ondiag_scale, frscale, omega, lrscale) &
   !$omp private(iat, izp, is, ish, jat, jzp, js, jsh) &
   !$omp private(gam, vec, r1, r1g, denom, rsh, scale, damp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = ish_at(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         js = ish_at(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         r1g = r1**gexp
         damp = exp (-(hubbard_exp + hubbard_exp_r0 * rad(izp, jzp)) * r1)
         do ish = 1, nsh_id(izp)
            do jsh = 1, nsh_id(jzp)
               scale = offdiag_scale(jsh, ish, jzp, izp) / damp
               gam = hubbard(jsh, ish, jzp, izp) * scale
               denom = (r1g + gam**(-gexp))**(1.0_wp/gexp)
               rsh = (frscale+lrscale*erf(omega*r1)) / denom

               g_mulliken(js+jsh, is+ish) = rsh
               g_mulliken(is+ish, js+jsh) = rsh
            end do
         end do
      end do
      ! Onsite terms
      do ish = 1, nsh_id(izp)
         do jsh = 1, ish-1
            scale = offdiag_scale(jsh, ish, izp, izp)
            gam = hubbard(jsh, ish, izp, izp) * scale * frscale

            g_mulliken(is+jsh, is+ish) = gam
            g_mulliken(is+ish, is+jsh) = gam
         end do
         ! Diagonal elements
         gam = hubbard(ish, ish, izp, izp) * ondiag_scale * frscale

         g_mulliken(is+ish, is+ish) = gam
      end do
   end do 

end subroutine get_gmulliken_0d


!> Evaluate range separated exchange matrix for finite systems
subroutine get_gmulliken_3d(mol, nsh_id, ish_at, hubbard, ondiag_scale, &
   & offdiag_scale, hubbard_exp, hubbard_exp_r0, rad, gexp, frscale, &
   & omega, lrscale, rcut, wsc, alpha, g_mulliken)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Hubbard parameter parameter for each shell
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Diagonal scaling of the Fock exchange
   real(wp), intent(in) :: ondiag_scale
   !> Off-diagonal scaling of the Fock exchange
   real(wp), intent(in) :: offdiag_scale(:, :, :, :)
   !> Exponent of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp
   !> Radius prefactor of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp_r0
   !> Radius for hubbard scaling
   real(wp), intent(in) :: rad(:, :)
   !> Exponent of exchange kernel
   real(wp), intent(in) :: gexp
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Range separation parameter
   real(wp), intent(in) :: omega
   !> Long-range scaling factor
   real(wp), intent(in) :: lrscale
   !> Long-range cutoff
   real(wp), intent(in) :: rcut
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence factor
   real(wp), intent(in) :: alpha
   !> Mulliken exchange matrix
   real(wp), intent(out) :: g_mulliken(:, :)

   g_mulliken(:, :) = 0.0_wp

end subroutine get_gmulliken_3d


!> Evaluate the gradient of the Mulliken exchange energy
subroutine get_mulliken_derivs(self, mol, cache, mulliken_grad, gradient, sigma)
   !> Instance of the exchange container
   class(exchange_fock), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(exchange_cache), intent(inout) :: cache
   !> Operator gradient w.r.t. the gamma Mulliken matrix
   real(wp), contiguous, intent(in) :: mulliken_grad(:, :)
   !> Molecular gradient of the exchange energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   if (any(mol%periodic)) then
      call get_gmulliken_derivs_3d(mol, self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, &
         & self%hubbard, self%ondiag_scale, self%offdiag_scale, self%hubbard_exp, &
         & self%hubbard_exp_r0, self%rad, self%gexp, self%frscale, self%omega, &
         & self%lrscale, cache%wsc, cache%alpha, mulliken_grad, gradient, sigma)
   else
      call get_gmulliken_derivs_0d(mol, self%nsh_id, self%nao_sh, self%ish_at, self%iao_sh, &
         & self%hubbard, self%ondiag_scale, self%offdiag_scale, self%hubbard_exp, &
         & self%hubbard_exp_r0, self%rad, self%gexp, self%frscale, self%omega, &
         & self%lrscale, mulliken_grad, gradient, sigma)
   end if

end subroutine get_mulliken_derivs


!> Evaluate derivatives of Mulliken exchange matrix for finite systems (0D)
subroutine get_gmulliken_derivs_0d(mol, nsh_id, nao_sh, ish_at, iao_sh, hubbard, &
   & ondiag_scale, offdiag_scale, hubbard_exp, hubbard_exp_r0, rad, &
   & gexp, frscale, omega, lrscale, mulliken_grad, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Hubbard parameter parameter for each shell
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Diagonal scaling of the Fock exchange (unused here)
   real(wp), intent(in) :: ondiag_scale
   !> Off-diagonal scaling of the Fock exchange
   real(wp), intent(in) :: offdiag_scale(:, :, :, :)
   !> Exponent of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp
   !> Radius prefactor of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp_r0
   !> Radius for hubbard scaling
   real(wp), intent(in) :: rad(:, :)
   !> Exponent of exchange kernel
   real(wp), intent(in) :: gexp
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Range separation parameter
   real(wp), intent(in) :: omega
   !> Long-range scaling factor
   real(wp), intent(in) :: lrscale
   !> Operator gradient w.r.t. the Mulliken gamma matrix
   real(wp), intent(in) :: mulliken_grad(:, :)
   !> Molecular gradient of the exchange energy
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, izp, jzp, is, js, ii, jj, ish, jsh, iao, jao
   real(wp) :: vec(3), r1, r1g, scale, damp, exparg, rsh, drsh, gam
   real(wp) :: denom, denom_pow, denom_deriv, tmp, shell_grad, dG(3)

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: gradient_local(:, :)

   !$omp parallel default(none) &
   !$omp shared(mol, nsh_id, nao_sh, ish_at, iao_sh, hubbard) &
   !$omp shared(gexp, hubbard_exp, hubbard_exp_r0, rad, offdiag_scale) &
   !$omp shared(frscale, omega, lrscale, mulliken_grad, gradient, sigma) &
   !$omp private(iat, izp, is, ii, ish, iao, jat, jzp, js, jj, jsh, jao) &
   !$omp private(vec, r1, r1g, scale, damp, exparg, rsh, drsh, gam, denom) &
   !$omp private(denom_pow, denom_deriv, tmp, shell_grad, dG, gradient_local)
   allocate(gradient_local(3, mol%nat), source = 0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = ish_at(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         js = ish_at(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         r1g = r1**gexp

         ! Radius dependent hardness damping factor
         exparg = (hubbard_exp + hubbard_exp_r0 * rad(izp, jzp))
         damp = exp(-exparg * r1)
         
         ! Range-separation factor and derivative
         rsh = (frscale+lrscale*erf(omega*r1)) 
         drsh = lrscale * 2.0_wp*omega / sqrtpi * exp(-(omega*r1)**2)

         do ish = 1, nsh_id(izp)
            ii = iao_sh(is+ish)
            do jsh = 1, nsh_id(jzp)
               jj = iao_sh(js+jsh)

               ! Diagonal/off-diagonal scaling factor for hubbard parameter
               scale = offdiag_scale(jsh, ish, jzp, izp) / damp

               ! Scaled hubbard parameter
               gam = hubbard(jsh, ish, jzp, izp) * scale

               ! Damped coulomb interaction denominator
               denom = r1g + gam**(-gexp)
               denom_pow = denom**(1.0_wp/gexp)
               denom_deriv = denom_pow * denom
               
               ! Derivative of range-separation factor
               tmp = drsh / denom_pow

               ! Derivative of damped coulomb interaction denominator
               tmp = tmp + rsh * ( -r1g/r1 + exparg * gam**(-gexp)) / denom_deriv

               ! Collect all operator gradient contributions per shell pair
               shell_grad = 0.0_wp
               do iao = 1, nao_sh(is + ish)
                  do jao = 1, nao_sh(js + jsh)
                     shell_grad = shell_grad + mulliken_grad(ii+iao, jj+jao)
                  end do
               end do
               
               ! Add operator contribution to the gradient
               dG = shell_grad * tmp * vec(:)/r1
               gradient_local(:, iat) = gradient_local(:, iat) + dG
               gradient_local(:, jat) = gradient_local(:, jat) - dG
            end do
         end do
      end do
   end do
   !$omp critical (get_gmulliken_derivs_0d_)
   gradient(:, :) = gradient + gradient_local
   !$omp end critical (get_gmulliken_derivs_0d_)
   deallocate(gradient_local)
   !$omp end parallel

end subroutine get_gmulliken_derivs_0d


!> Evaluate derivatives of Mulliken exchange matrix for periodic systems (3D)
subroutine get_gmulliken_derivs_3d(mol, nsh_id, nao_sh, ish_at, iao_sh, hubbard, &
   & ondiag_scale, offdiag_scale, hubbard_exp, hubbard_exp_r0, rad, &
   & gexp, frscale, omega, lrscale, wsc, alpha, mulliken_grad, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Hubbard parameter parameter for each shell
   real(wp), intent(in) :: hubbard(:, :, :, :)
   !> Diagonal scaling of the Fock exchange
   real(wp), intent(in) :: ondiag_scale
   !> Off-diagonal scaling of the Fock exchange
   real(wp), intent(in) :: offdiag_scale(:, :, :, :)
   !> Exponent of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp
   !> Radius prefactor of radius dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp_r0
   !> Radius for hubbard scaling
   real(wp), intent(in) :: rad(:, :)
   !> Exponent of exchange kernel
   real(wp), intent(in) :: gexp
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Range separation parameter
   real(wp), intent(in) :: omega
   !> Long-range scaling factor
   real(wp), intent(in) :: lrscale
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence factor
   real(wp), intent(in) :: alpha
   !> Operator gradient w.r.t. the Mulliken gamma matrix
   real(wp), intent(in) :: mulliken_grad(:, :)
   !> Molecular gradient of the exchange energy
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), intent(inout) :: sigma(:, :)

   gradient(:, :) = gradient + 0.0_wp
   sigma(:, :) = gradient + 0.0_wp
end subroutine get_gmulliken_derivs_3d


!> Evaluate onsite exchange gamma matrix
subroutine get_onsite_Kmatrix(self, mol, wfn, cache)
   !> Instance of the exchange container
   class(exchange_fock), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Reusable data container with intermediates and the final exchange matrices
   type(exchange_cache), intent(inout) :: cache

   call get_gons(mol, self%nsh_id, self%ish_at, self%onecxints, self%frscale, &
      & self%kq, wfn%qsh(:, 1), cache%g_onsfx, cache%g_onsri, cache%dgdq_onsfx, &
      & cache%dgdq_onsri)

end subroutine get_onsite_Kmatrix


!> Evaluate onsite exchange and rotational invariance correction matrices
subroutine get_gons(mol, nsh_id, ish_at, onecxints, frscale, kq, qsh, &
   & g_onsfx, g_onsri, dgdq_onsfx, dgdq_onsri)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> One center exchange integrals
   real(wp), intent(in) :: onecxints(:, :, :)
   !> Full-range scale for K
   real(wp), intent(in) :: frscale
   !> Charge-dependence of the effective Fock exchange 
   real(wp), intent(in) :: kq(:, :)
   !> Shell-resolved charges
   real(wp), intent(in) :: qsh(:)
   !> Onsite exchange matrix
   real(wp), intent(out) :: g_onsfx(:, :, :)
   !> Onsite rotational invariance correction matrix
   real(wp), intent(out) :: g_onsri(:, :)
   !> Charge-derivative of the onsite exchange matrix
   real(wp), intent(out) :: dgdq_onsfx(:, :, :)
   !> Charge-derivative of the onsite rotational invariance correction matrix
   real(wp), intent(out) :: dgdq_onsri(:, :)

   integer :: iat, izp, is, ish, jsh
   real(wp) :: gam, dgami, dgamj, denom

   g_onsfx(:, :, :) = 0.0_wp
   g_onsri(:, :) = 0.0_wp
   dgdq_onsfx(:, :, :) = 0.0_wp
   dgdq_onsri(:, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(g_onsfx, g_onsri, dgdq_onsfx, dgdq_onsri) &
   !$omp shared(onecxints, mol, frscale, kq, qsh, nsh_id, ish_at) &
   !$omp private(iat, izp, is, ish, jsh, gam, dgami, dgamj, denom)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = ish_at(iat)
      do ish = 1, nsh_id(izp)
         do jsh = 1, nsh_id(izp)

            ! Scaled onsite exchange integral
            gam = frscale * onecxints(jsh, ish, izp) * (1.0_wp - 0.5_wp * &
               & (kq(ish, izp) * qsh(is+ish) + kq(jsh, izp) * qsh(is+jsh)))
            dgami = -0.5_wp * frscale * kq(ish, izp) * onecxints(jsh, ish, izp)
            dgamj = -0.5_wp * frscale * kq(jsh, izp) * onecxints(jsh, ish, izp)

            ! Onsite correction exchange matrix
            g_onsfx(jsh, ish, iat) = gam
            if (ish == jsh) then
               dgdq_onsfx(jsh, ish, is+ish) = dgami + dgamj
            else
               dgdq_onsfx(jsh, ish, is+ish) = dgami
               dgdq_onsfx(jsh, ish, is+jsh) = dgamj
            end if

            ! Rotational invariance correction matrix
            if (ish == jsh .and. ish > 1) then
               denom = 1.0_wp / real(4*ish - 2, wp)
               g_onsri(ish, iat) = gam * denom
               dgdq_onsri(ish, is+ish) = (dgami + dgamj) * denom
            end if
         end do
      end do
   end do 

end subroutine get_gons


!> Evaluate bond-order correlation correction matrix
subroutine get_bocorr_Kmatrix(self, mol, cache)
   !> Instance of the exchange container
   class(exchange_fock), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(exchange_cache), intent(inout) :: cache

   if (any(mol%periodic)) then
      call get_gbocorr_3d(mol, self%corr_scale, self%corr_exp, self%corr_rad, &
         & self%rad, self%rcut, cache%wsc, cache%alpha, cache%g_bocorr)
   else
      call get_gbocorr_0d(mol, self%corr_scale, self%corr_exp, self%corr_rad, &
         & self%rad, cache%g_bocorr)
   end if

end subroutine get_bocorr_Kmatrix


!> Evaluate bond-order correlation correction matrix for finite systems
subroutine get_gbocorr_0d(mol, corr_scale, corr_exp, corr_rad, &
   & rad, g_bocorr)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Bond-order correlation scaling factor for each atom pair
   real(wp), intent(in) :: corr_scale(:, :)
   !> Bond-order correlation damping exponent
   real(wp), intent(in) :: corr_exp
   !> Bond-order correlation radius for each atom pair
   real(wp), intent(in) :: corr_rad(:, :)
   !> Reference van-der-Waals radius
   real(wp), intent(in) :: rad(:, :)
   !> bond_order correlation correction matrix
   real(wp), intent(out) :: g_bocorr(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r1, arg, damp, corr

   g_bocorr(:, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) shared(g_bocorr) &
   !$omp shared(mol, corr_rad, corr_exp, corr_scale, rad) &
   !$omp private(iat, izp, jat, jzp, vec, r1, arg, damp, corr)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         arg = corr_exp * (r1 - corr_rad(izp, jzp)) / rad(izp, jzp)
         damp = 0.5_wp * (1.0_wp + erf(-arg))
         corr = corr_scale(izp, jzp) * damp

         g_bocorr(jat, iat) = corr
         g_bocorr(iat, jat) = corr
      end do
   end do

end subroutine get_gbocorr_0d


!> Evaluate bond-order correlation correction matrix for periodic systems (3D)
subroutine get_gbocorr_3d(mol, corr_scale, corr_exp, corr_rad, &
   & rad, rcut, wsc, alpha, g_bocorr)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Bond-order correlation scaling factor for each atom pair
   real(wp), intent(in) :: corr_scale(:, :)
   !> Bond-order correlation damping exponent
   real(wp), intent(in) :: corr_exp
   !> Bond-order correlation radius for each atom pair
   real(wp), intent(in) :: corr_rad(:, :)
   !> Reference van-der-Waals radius
   real(wp), intent(in) :: rad(:, :)
   !> Long-range cutoff
   real(wp), intent(in) :: rcut
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence factor
   real(wp), intent(in) :: alpha
   !> Bond-order correlation correction matrix
   real(wp), intent(out) :: g_bocorr(:, :)

   g_bocorr(:, :) = 0.0_wp

end subroutine get_gbocorr_3d


!> Evaluate the gradient of the bond-order correlation energy
subroutine get_bocorr_derivs(self, mol, cache, bocorr_grad, gradient, sigma)
   !> Instance of the exchange container
   class(exchange_fock), intent(in) :: self
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

   if (any(mol%periodic)) then
      call get_gbocorr_derivs_3d(mol, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, self%corr_scale, self%corr_exp, self%corr_rad, &
         & self%rad, cache%wsc, cache%alpha, bocorr_grad, gradient, sigma)
   else
      call get_gbocorr_derivs_0d(mol, self%nsh_id, self%nao_sh, self%ish_at, &
         & self%iao_sh, self%corr_scale, self%corr_exp, self%corr_rad, &
         & self%rad, bocorr_grad, gradient, sigma)
   end if

end subroutine get_bocorr_derivs

!> Evaluate derivatives of Mulliken exchange matrix for finite systems (0D)
subroutine get_gbocorr_derivs_0d(mol, nsh_id, nao_sh, ish_at, iao_sh, corr_scale, &
   & corr_exp, corr_rad, rad, bocorr_grad, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Bond-order correlation scaling factor for each atom pair
   real(wp), intent(in) :: corr_scale(:, :)
   !> Bond-order correlation damping exponent
   real(wp), intent(in) :: corr_exp
   !> Bond-order correlation radius for each atom pair
   real(wp), intent(in) :: corr_rad(:, :)
   !> Reference van-der-Waals radius
   real(wp), intent(in) :: rad(:, :)
   !> Operator gradient w.r.t. the bond-order correlation matrix
   real(wp), intent(in) :: bocorr_grad(:, :)
   !> Molecular gradient of the exchange energy
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), intent(inout) :: sigma(:, :)
   
   integer :: iat, jat, izp, jzp, is, js, ii, jj, ish, jsh, iao, jao
   real(wp) :: vec(3), r1, r1g, arg, ddamp, dcorr, atom_grad, dG(3)

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: gradient_local(:, :)

   !$omp parallel default(none) &
   !$omp shared(mol, nsh_id, nao_sh, ish_at, iao_sh, corr_rad, corr_exp) &
   !$omp shared(corr_scale, rad, bocorr_grad, gradient, sigma) &
   !$omp private(iat, izp, is, ii, ish, iao, jat, jzp, js, jj, jsh, jao) &
   !$omp private(vec, r1, arg, ddamp, dcorr, atom_grad, dG, gradient_local)
   allocate(gradient_local(3, mol%nat), source = 0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = ish_at(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         js = ish_at(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         arg = corr_exp * (r1 - corr_rad(izp, jzp)) / rad(izp, jzp)
         ddamp = corr_exp * exp(-arg**2) / (sqrtpi * rad(izp, jzp))
         dcorr = corr_scale(izp, jzp) * ddamp

         ! Collect all operator gradient contributions per atom pair
         atom_grad = 0.0_wp
         do ish = 1, nsh_id(izp)
            ii = iao_sh(is+ish)
            do jsh = 1, nsh_id(jzp)
               jj = iao_sh(js+jsh)
               do iao = 1, nao_sh(is + ish)
                  do jao = 1, nao_sh(js + jsh)
                     atom_grad = atom_grad + bocorr_grad(ii+iao, jj+jao)
                  end do
               end do
            end do
         end do
         ! Add operator contribution to the gradient
         dG(:) = atom_grad * dcorr * vec/r1
         gradient_local(:, iat) = gradient_local(:, iat) + dG
         gradient_local(:, jat) = gradient_local(:, jat) - dG
      end do
   end do
   !$omp critical (get_gmulliken_derivs_0d_)
   gradient(:, :) = gradient + gradient_local
   !$omp end critical (get_gmulliken_derivs_0d_)
   deallocate(gradient_local)
   !$omp end parallel

end subroutine get_gbocorr_derivs_0d


!> Evaluate derivatives of bond-order correlation matrix for periodic systems (3D)
subroutine get_gbocorr_derivs_3d(mol, nsh_id, nao_sh, ish_at, iao_sh, corr_scale, &
   & corr_exp, corr_rad, rad, wsc, alpha, bocorr_grad, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Bond-order correlation scaling factor for each atom pair
   real(wp), intent(in) :: corr_scale(:, :)
   !> Bond-order correlation damping exponent
   real(wp), intent(in) :: corr_exp
   !> Bond-order correlation radius for each atom pair
   real(wp), intent(in) :: corr_rad(:, :)
   !> Reference van-der-Waals radius
   real(wp), intent(in) :: rad(:, :)
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence factor
   real(wp), intent(in) :: alpha
   !> Operator gradient w.r.t. the bond-order correlation matrix
   real(wp), intent(in) :: bocorr_grad(:, :)

   !> Molecular gradient of the exchange energy
   real(wp), intent(inout) :: gradient(:, :)
   !> Strain derivatives of the exchange energy
   real(wp), intent(inout) :: sigma(:, :)

   gradient(:, :) = gradient + 0.0_wp
   sigma(:, :) = gradient + 0.0_wp
end subroutine get_gbocorr_derivs_3d


subroutine shell_hadamard_add(nsh, nao_sh, iao_sh, g_sh, src, alpha, dst, trans_src)
   !> Number of shells in the molecule
   integer, intent(in) :: nsh
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Shell-resolved exchange matrix: [nsh, nsh]
   real(wp), intent(in) :: g_sh(:, :)
   !> Source AO matrix
   real(wp), intent(in) :: src(:, :)
   !> Prefactor
   real(wp), intent(in) :: alpha
   !> Destination AO matrix
   real(wp), intent(inout) :: dst(:, :)
   !> Whether src should be read as transposed
   logical, intent(in), optional :: trans_src

   integer :: ish, jsh, ii, jj, ni, nj, iao, jao
   real(wp) :: scale
   logical :: trans

   trans = .false.
   if (present(trans_src)) trans = trans_src

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(nsh, nao_sh, iao_sh, g_sh, src, alpha, dst, trans) &
   !$omp private(ish, jsh, ii, jj, ni, nj, iao, jao, scale)
   do ish = 1, nsh
      ii = iao_sh(ish)
      ni = nao_sh(ish)

      do jsh = 1, nsh
         jj = iao_sh(jsh)
         nj = nao_sh(jsh)

         scale = alpha * g_sh(jsh, ish)
         if (abs(scale) < epsilon(1.0_wp)) cycle

         if (.not.trans) then
            dst(jj+1:jj+nj, ii+1:ii+ni) = dst(jj+1:jj+nj, ii+1:ii+ni) + &
               scale * src(jj+1:jj+nj, ii+1:ii+ni)
         else
            do iao = 1, ni
               do jao = 1, nj
                  dst(jj+jao, ii+iao) = dst(jj+jao, ii+iao) + &
                     scale * src(ii+iao, jj+jao)
               end do
            end do
         end if
      end do
   end do

end subroutine shell_hadamard_add


subroutine atom_hadamard_add(nat, id, nsh_id, nao_sh, ish_at, iao_sh, &
   & g_at, src, alpha, dst, trans_src)
   !> Number of atoms in the molecule
   integer, intent(in) :: nat
   !> Species identifier for each atom
   integer, intent(in) :: id(:)
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Atom-resolved matrix: [nat, nat]
   real(wp), intent(in) :: g_at(:, :)
   !> Source AO matrix
   real(wp), intent(in) :: src(:, :)
   !> Prefactor
   real(wp), intent(in) :: alpha
   !> Destination AO matrix
   real(wp), intent(inout) :: dst(:, :)
   !> Whether src should be read as transposed
   logical, intent(in), optional :: trans_src

   integer :: iat, jat, izp, jzp, is, js, ish, jsh, ii, jj, ni, nj, iao, jao
   real(wp) :: scale
   logical :: trans

   trans = .false.
   if (present(trans_src)) trans = trans_src

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(nat, id, nsh_id, nao_sh, ish_at, iao_sh, g_at, src, alpha, dst, trans) &
   !$omp private(iat, jat, izp, jzp, is, js, ish, jsh, ii, jj, ni, nj, iao, jao, scale)
   do iat = 1, nat
      izp = id(iat)
      is  = ish_at(iat)

      do jat = 1, nat
         jzp = id(jat)
         js  = ish_at(jat)

         scale = alpha * g_at(jat, iat)
         if (abs(scale) < epsilon(1.0_wp)) cycle

         do ish = 1, nsh_id(izp)
            ii = iao_sh(is + ish)
            ni = nao_sh(is + ish)

            do jsh = 1, nsh_id(jzp)
               jj = iao_sh(js + jsh)
               nj = nao_sh(js + jsh)

               if (.not.trans) then
                  dst(jj+1:jj+nj, ii+1:ii+ni) = dst(jj+1:jj+nj, ii+1:ii+ni) + &
                     & scale * src(jj+1:jj+nj, ii+1:ii+ni)
               else
                  do iao = 1, ni
                     do jao = 1, nj
                        dst(jj+jao, ii+iao) = dst(jj+jao, ii+iao) + &
                           & scale * src(ii+iao, jj+jao)
                     end do
                  end do
               end if
            end do
         end do
      end do
   end do

end subroutine atom_hadamard_add


subroutine onsite_fx_hadamard_add(nat, id, nsh_id, nao_sh, ish_at, iao_sh, &
   & g_onsfx, src, alpha, dst, trans_src)
   !> Number of atoms in the molecule
   integer, intent(in) :: nat
   !> Species identifier for each atom
   integer, intent(in) :: id(:)
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Compact shell-block onsite exchange matrix: [maxsh, maxsh, nat]
   real(wp), intent(in) :: g_onsfx(:, :, :)
   !> Source AO matrix
   real(wp), intent(in) :: src(:, :)
   !> Prefactor
   real(wp), intent(in) :: alpha
   !> Destination AO matrix
   real(wp), intent(inout) :: dst(:, :)
   !> Whether src should be read as transposed
   logical, intent(in), optional :: trans_src

   integer :: iat, izp, is, ish, jsh, ii, jj, ni, nj, iao, jao
   real(wp) :: scale
   logical :: trans

   trans = .false.
   if (present(trans_src)) trans = trans_src

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(nat, id, ish_at, nsh_id, iao_sh, nao_sh, alpha, g_onsfx, src, dst, trans) &
   !$omp private(iat, izp, is, ish, jsh, ii, jj, ni, nj, iao, jao, scale)
   do iat = 1, nat
      izp = id(iat)
      is  = ish_at(iat)

      do ish = 1, nsh_id(izp)
         ii = iao_sh(is + ish)
         ni = nao_sh(is + ish)

         do jsh = 1, nsh_id(izp)
            jj = iao_sh(is + jsh)
            nj = nao_sh(is + jsh)

            scale = alpha * g_onsfx(jsh, ish, iat)
            if (abs(scale) < epsilon(1.0_wp)) cycle

            if (.not.trans) then
               dst(jj+1:jj+nj, ii+1:ii+ni) = dst(jj+1:jj+nj, ii+1:ii+ni) + &
                  & scale * src(jj+1:jj+nj, ii+1:ii+ni)
            else
               do iao = 1, ni
                  do jao = 1, nj
                     dst(jj+jao, ii+iao) = dst(jj+jao, ii+iao) + &
                        & scale * src(ii+iao, jj+jao)
                  end do
               end do
            end if
         end do
      end do
   end do

end subroutine onsite_fx_hadamard_add


subroutine onsite_fx_symv(nat, id, nsh_id, nao_sh, ish_at, iao_sh, &
   & g_onsfx, xvec, yvec)
   !> Number of atoms in the molecule
   integer, intent(in) :: nat
   !> Species identifier for each atom
   integer, intent(in) :: id(:)
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Compact diagonal storage of the onsite exchange matrix: [maxsh, maxsh, nat]
   real(wp), intent(in) :: g_onsfx(:, :, :)
   !> Input vector to be contracted: [nao]
   real(wp), intent(in) :: xvec(:)
   !> Output vector: [nao]
   real(wp), intent(out) :: yvec(:)

   integer :: iat, izp, is, ish, jsh, ii, jj, ni, nj
   real(wp) :: shsum

   yvec = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(nat, id, ish_at, nsh_id, iao_sh, nao_sh, g_onsfx, xvec, yvec) &
   !$omp private(iat, izp, is, ish, jsh, ii, jj, ni, nj, shsum)
   do iat = 1, nat
      izp = id(iat)
      is  = ish_at(iat)

      do ish = 1, nsh_id(izp)
         ii = iao_sh(is+ish)
         ni = nao_sh(is+ish)

         shsum = sum(xvec(ii+1:ii+ni))
         if (abs(shsum) < epsilon(1.0_wp)) cycle

         do jsh = 1, nsh_id(izp)
            jj = iao_sh(is+jsh)
            nj = nao_sh(is+jsh)

            yvec(jj+1:jj+nj) = yvec(jj+1:jj+nj) + g_onsfx(jsh, ish, iat) * shsum
         end do
      end do
   end do

end subroutine onsite_fx_symv


subroutine onsite_ri_hadamard_add(nat, id, nsh_id, nao_sh, ish_at, iao_sh, &
   & g_onsri, src, alpha, dst)
   !> Number of atoms in the molecule
   integer, intent(in) :: nat
   !> Species identifier for each atom
   integer, intent(in) :: id(:)
   !> Number of shells for each species
   integer, intent(in) :: nsh_id(:)
   !> Number of spherical atomic orbitals for each shell
   integer, intent(in) :: nao_sh(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Index offset for each shell in the atomic orbital space
   integer, intent(in) :: iao_sh(:)
   !> Compact shell-diagonal onsite rotational invariance matrix: [maxsh, nat]
   real(wp), intent(in) :: g_onsri(:, :)
   !> Source AO matrix
   real(wp), intent(in) :: src(:, :)
   !> Prefactor
   real(wp), intent(in) :: alpha
   !> Destination AO matrix
   real(wp), intent(inout) :: dst(:, :)

   integer :: iat, izp, is, ish, ii, ni
   real(wp) :: scale

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(nat, id, ish_at, nsh_id, iao_sh, nao_sh, alpha, g_onsri, src, dst) &
   !$omp private(iat, izp, is, ish, ii, ni, scale)
   do iat = 1, nat
      izp = id(iat)
      is  = ish_at(iat)

      do ish = 1, nsh_id(izp)
         ii = iao_sh(is + ish)
         ni = nao_sh(is + ish)

         scale = alpha * g_onsri(ish, iat)
         if (abs(scale) < epsilon(1.0_wp)) cycle

         dst(ii+1:ii+ni, ii+1:ii+ni) = dst(ii+1:ii+ni, ii+1:ii+ni) + &
            & scale * src(ii+1:ii+ni, ii+1:ii+ni)
      end do
   end do

end subroutine onsite_ri_hadamard_add


end module tblite_exchange_fock
