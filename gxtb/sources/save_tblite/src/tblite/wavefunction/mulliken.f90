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

!> @file tblite/wavefunction/mulliken.f90
!> Provides Mulliken population analysis

!> Wavefunction analysis via Mulliken populations
module tblite_wavefunction_mulliken
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_blas, only: gemm
   use tblite_wavefunction_spin, only : updown_to_magnet
   implicit none
   private

   public :: get_mulliken_shell_charges, get_mulliken_atomic_multipoles
   public :: get_molecular_dipole_moment, get_molecular_quadrupole_moment
   public :: get_mayer_bond_orders, get_mulliken_shell_multipoles

contains


subroutine get_mulliken_shell_charges(bas, smat, pmat, n0sh, qsh)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: smat(:, :)
   real(wp), intent(in) :: pmat(:, :, :)
   real(wp), intent(in) :: n0sh(:)
   real(wp), intent(out) :: qsh(:, :)

   integer :: iao, jao, spin
   real(wp) :: pao

   qsh(:, :) = 0.0_wp
   !$omp parallel do default(none) collapse(2) schedule(runtime) reduction(+:qsh) &
   !$omp shared(bas, pmat, smat) private(spin, iao, jao, pao)
   do spin = 1, size(pmat, 3)
      do iao = 1, bas%nao
         pao = 0.0_wp
         do jao = 1, bas%nao
            pao = pao + pmat(jao, iao, spin) * smat(jao, iao)
         end do
         qsh(bas%ao2sh(iao), spin) = qsh(bas%ao2sh(iao), spin) - pao
      end do
   end do

   call updown_to_magnet(qsh)
   qsh(:, 1) = qsh(:, 1) + n0sh

end subroutine get_mulliken_shell_charges


subroutine get_mulliken_atomic_multipoles(bas, mpmat, pmat, mpat)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: mpmat(:, :, :)
   real(wp), intent(in) :: pmat(:, :, :)
   real(wp), intent(out) :: mpat(:, :, :)

   integer :: iao, jao, spin
   real(wp) :: pao(size(mpmat, 1))

   mpat(:, :, :) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) reduction(+:mpat) &
   !$omp shared(bas, pmat, mpmat) private(spin, iao, jao, pao)
   do spin = 1, size(pmat, 3)
      do iao = 1, bas%nao
         pao(:) = 0.0_wp
         do jao = 1, bas%nao
            pao(:) = pao + pmat(jao, iao, spin) * mpmat(:, jao, iao)
         end do
         mpat(:, bas%ao2at(iao), spin) = mpat(:, bas%ao2at(iao), spin) - pao
      end do
   end do

   call updown_to_magnet(mpat)

end subroutine get_mulliken_atomic_multipoles

subroutine get_mulliken_shell_multipoles(bas, mpmat, pmat, mpat)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: mpmat(:, :, :)
   real(wp), intent(in) :: pmat(:, :, :)
   real(wp), intent(out) :: mpat(:, :, :)

   integer :: iao, jao, spin
   real(wp) :: pao(size(mpmat, 1))
   mpat(:, :, :) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) reduction(+:mpat) &
   !$omp shared(bas, pmat, mpmat) private(spin, iao, jao, pao)
   do spin = 1, size(pmat, 3)
      do iao = 1, bas%nao
         pao(:) = 0.0_wp
         do jao = 1, bas%nao
            pao(:) = pao + pmat(jao, iao, spin) * mpmat(:, jao, iao)
         end do
         mpat(:, bas%ao2sh(iao), spin) = mpat(:, bas%ao2sh(iao), spin) - pao
      end do
   end do

   call updown_to_magnet(mpat)

end subroutine get_mulliken_shell_multipoles

subroutine get_molecular_dipole_moment(mol, qat, dpat, dpmom)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: dpat(:, :)
   real(wp), intent(out) :: dpmom(:)

   integer :: iat

   dpmom(:) = 0.0_wp
   do iat = 1, mol%nat
      dpmom(:) = dpmom + mol%xyz(:, iat) * qat(iat) + dpat(:, iat)
   end do
end subroutine get_molecular_dipole_moment

subroutine get_molecular_quadrupole_moment(mol, qat, dpat, qpat, qpmom)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: qat(:)
   real(wp), intent(in) :: dpat(:, :)
   real(wp), intent(in) :: qpat(:, :)
   real(wp), intent(out) :: qpmom(:)

   integer :: iat
   real(wp) :: vec(3), cart(6), tr

   qpmom(:) = 0.0_wp
   do iat = 1, mol%nat
      vec(:) = mol%xyz(:, iat)*qat(iat)
      cart([1, 3, 6]) = mol%xyz(:, iat) * (vec + 2*dpat(:, iat))
      cart(2) = mol%xyz(1, iat) * (vec(2) + dpat(2, iat)) + dpat(1, iat)*mol%xyz(2, iat)
      cart(4) = mol%xyz(1, iat) * (vec(3) + dpat(3, iat)) + dpat(1, iat)*mol%xyz(3, iat)
      cart(5) = mol%xyz(2, iat) * (vec(3) + dpat(3, iat)) + dpat(2, iat)*mol%xyz(3, iat)
      tr = 0.5_wp * (cart(1) + cart(3) + cart(6))
      cart(1) = 1.5_wp * cart(1) - tr
      cart(2) = 3.0_wp * cart(2)
      cart(3) = 1.5_wp * cart(3) - tr
      cart(4) = 3.0_wp * cart(4)
      cart(5) = 3.0_wp * cart(5)
      cart(6) = 1.5_wp * cart(6) - tr
      qpmom(:) = qpmom(:) + qpat(:, iat) + cart
   end do
end subroutine get_molecular_quadrupole_moment


!> Evaluate Wiberg/Mayer bond orders including factor 2 scaling for open-shell cases
subroutine get_mayer_bond_orders(mol, bas, smat, pmat, mbo)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap matrix
   real(wp), intent(in) :: smat(:, :)
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :, :)
   !> Wiberg/Mayer bond orders
   real(wp), intent(out) :: mbo(:, :, :)

   integer :: iao, jao, iat, jat, ish, jsh, ii, jj, is, js, spin, nspin
   real(wp) :: pao, uhf_scale
   real(wp), allocatable :: psmat(:, :)

   allocate(psmat(bas%nao, bas%nao))

   nspin = size(pmat, 3)

   ! Scaling factor of 2 for unrestricted case
   uhf_scale = 1.0_wp
   if (nspin == 2) then
      uhf_scale = 2.0_wp
   end if

   mbo(:, :, :) = 0.0_wp
   do spin = 1, nspin
      call gemm(pmat(:, :, spin), smat, psmat)
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, bas, psmat, mbo, spin, uhf_scale) &
      !$omp private(iat, jat, is, js, ish, jsh, ii, jj, iao, jao, pao)
      do iat = 1, mol%nat
         is = bas%ish_at(iat)
         do jat = 1, iat - 1
            pao = 0.0_wp
            js = bas%ish_at(jat)
            do ish = 1, bas%nsh_at(iat)
            ii = bas%iao_sh(is + ish)
               do jsh = 1, bas%nsh_at(jat)
                  jj = bas%iao_sh(js + jsh)
                  do iao = 1, bas%nao_sh(is + ish)
                     do jao = 1, bas%nao_sh(js + jsh)
                        pao = pao + psmat(jj+jao, ii+iao) * psmat(ii+iao, jj+jao)
                     end do
                  end do
               end do
            end do 
            mbo(jat, iat, spin) = pao * uhf_scale
            mbo(iat, jat, spin) = pao * uhf_scale
         end do
      end do
   end do

   call updown_to_magnet(mbo)
end subroutine get_mayer_bond_orders


end module tblite_wavefunction_mulliken
