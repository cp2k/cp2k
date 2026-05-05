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

!> @file tblite/wavefunction/type.f90
!> Provides a wavefunction type for holding density-related information

!> Declaration of a tight-binding wavefunction
module tblite_wavefunction_type
   use mctc_env, only : wp
   use tblite_blas, only : gemm
   implicit none
   private

   public :: new_wavefunction
   public :: get_density_matrix, get_alpha_beta_occupation

   !> Tight-binding wavefunction
   type, public :: wavefunction_type
      !> Electronic temperature
      real(wp) :: kt = 0.0_wp
      !> Number of electrons in this wavefunction
      real(wp) :: nocc = 0.0_wp
      !> Number of unpaired electrons in this wavefunction
      real(wp) :: nuhf = 0.0_wp
      !> Number of spin channels
      integer :: nspin = 1
      !> Number of electrons
      real(wp), allocatable :: nel(:)
      !> Reference occupation number for each atom, shape: [nat]
      real(wp), allocatable :: n0at(:)
      !> Reference occupation number for each shell, shape: [nsh]
      real(wp), allocatable :: n0sh(:)

      !> Density matrix, shape: [nao, nao, spin]
      real(wp), allocatable :: density(:, :, :)
      !> Orbital coefficients, shape: [nao, nao, spin]
      real(wp), allocatable :: coeff(:, :, :)
      !> Orbital energies, eigenvalues, shape: [nao, spin]
      real(wp), allocatable :: emo(:, :)
      !> Occupation numbers, shape: [nao, spin]
      real(wp), allocatable :: focc(:, :)

      !> Number of electrons for each atom, shape: [nat, spin]
      real(wp), allocatable :: qat(:, :)
      !> Number of electrons for each shell, shape: [nsh, spin]
      real(wp), allocatable :: qsh(:, :)

      !> Atomic dipole moments for each atom, shape: [3, nat, spin]
      real(wp), allocatable :: dpat(:, :, :)
      !> Atomic quadrupole moments for each atom, shape: [6, nat, spin]
      real(wp), allocatable :: qpat(:, :, :)

      !> Derivative of atomic charges w.r.t. the positions: [3, nat, nat, spin]
      real(wp), allocatable :: dqatdr(:, :, :, :)
      !> Derivative of atomic charges w.r.t. the lattice vectors: [3, 3, nat, spin]
      real(wp), allocatable :: dqatdL(:, :, :, :)

      !> Derivative of shell charges w.r.t. the positions: [3, nat, nsh, spin]
      real(wp), allocatable :: dqshdr(:, :, :, :)
      !> Derivative of shell charges w.r.t. the lattice vectors: [3, 3, nsh, spin]
      real(wp), allocatable :: dqshdL(:, :, :, :)
   end type wavefunction_type

contains


subroutine new_wavefunction(self, nat, nsh, nao, nspin, kt, grad)
   type(wavefunction_type), intent(out) :: self
   integer, intent(in) :: nat
   integer, intent(in) :: nsh
   integer, intent(in) :: nao
   integer, intent(in) :: nspin
   real(wp), intent(in) :: kt
   !> Flag to indicate if a wavefunction gradient is requested
   logical, intent(in), optional :: grad

   self%nspin = nspin
   self%kt = kt

   allocate(self%nel(max(2, nspin)))

   allocate(self%n0at(nat))
   allocate(self%n0sh(nsh))

   allocate(self%density(nao, nao, nspin))
   allocate(self%coeff(nao, nao, nspin))
   allocate(self%emo(nao, nspin))
   allocate(self%focc(nao, max(2, nspin)))

   allocate(self%qat(nat, nspin))
   allocate(self%qsh(nsh, nspin))

   allocate(self%dpat(3, nat, nspin))
   allocate(self%qpat(6, nat, nspin))

   ! Check if a wavefunction gradient is requested
   if(present(grad)) then
      if(grad) then
         allocate(self%dqatdr(3, nat, nat, nspin))
         allocate(self%dqatdL(3, 3, nat, nspin))

         allocate(self%dqshdr(3, nat, nsh, nspin))
         allocate(self%dqshdL(3, 3, nsh, nspin))

         self%dqatdr(:, :, :, :) = 0.0_wp
         self%dqatdL(:, :, :, :) = 0.0_wp
         self%dqshdr(:, :, :, :) = 0.0_wp
         self%dqshdL(:, :, :, :) = 0.0_wp
      end if
   end if

   self%density(:, :, :) = 0.0_wp
   self%coeff(:, :, :) = 0.0_wp
   self%qat(:, :) = 0.0_wp
   self%qsh(:, :) = 0.0_wp
   self%dpat(:, :, :) = 0.0_wp
   self%qpat(:, :, :) = 0.0_wp
end subroutine new_wavefunction

subroutine get_density_matrix(focc, coeff, pmat)
   real(wp), intent(in) :: focc(:)
   real(wp), contiguous, intent(in) :: coeff(:, :)
   real(wp), contiguous, intent(out) :: pmat(:, :)

   real(wp), allocatable :: scratch(:, :)
   integer :: iao, jao

   allocate(scratch(size(pmat, 1), size(pmat, 2)))
   !$omp parallel do collapse(2) default(none) schedule(runtime) &
   !$omp shared(scratch, coeff, focc, pmat) private(iao, jao)
   do iao = 1, size(pmat, 1)
      do jao = 1, size(pmat, 2)
         scratch(jao, iao) = coeff(jao, iao) * focc(iao)
      end do
   end do
   call gemm(scratch, coeff, pmat, transb='t')
end subroutine get_density_matrix


!> Split an real occupation number into alpha and beta space.
!>
!> This routine does not perform any checks on the condition
!> ``mod(nocc, 2) == 0 .eqv. mod(nuhf, 2) == 0`` and will yield fractional
!> occupations in case those condtions are not fullfilled.
!> However, it will avoid creating negative occupation numbers.
subroutine get_alpha_beta_occupation(nocc, nuhf, nalp, nbet)
   real(wp), intent(in) :: nocc
   real(wp), intent(in) :: nuhf
   real(wp), intent(out) :: nalp
   real(wp), intent(out) :: nbet

   real(wp) :: ntmp, diff

   ! make sure we cannot get a negative occupation here
   diff = min(nuhf, nocc)
   ntmp = nocc - diff

   nalp = ntmp / 2 + diff
   nbet = ntmp / 2
end subroutine get_alpha_beta_occupation

end module tblite_wavefunction_type
