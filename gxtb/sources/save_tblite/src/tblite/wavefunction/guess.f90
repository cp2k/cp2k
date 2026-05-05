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

!> @file tblite/wavefunction/guess.f90
!> Provides guesses for the wavefunction

!> Implementation of the guess wavefunctions
module tblite_wavefunction_guess
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use multicharge, only : mchrg_model_type, get_eeq_charges, get_eeqbc_charges
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_h0, only : get_occupation
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: sad_guess, eeq_guess, eeqbc_guess, shell_partition
   interface sad_guess
      module procedure sad_guess_qat
      module procedure sad_guess_qsh
   end interface sad_guess

   interface eeq_guess
      module procedure eeq_guess_qat
      module procedure eeq_guess_qsh
   end interface eeq_guess

   interface eeqbc_guess
      module procedure eeqbc_guess_qat
      module procedure eeqbc_guess_qsh
   end interface eeqbc_guess

contains

subroutine sad_guess_qat(mol, charges, dpat)
   type(structure_type), intent(in) :: mol
   real(wp), intent(inout) :: charges(:), dpat(:, :)

   dpat(:, :) = 0.0_wp
   charges = mol%charge / mol%nat
end subroutine sad_guess_qat

subroutine sad_guess_qsh(mol, calc, wfn)
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(inout) :: wfn

   call sad_guess_qat(mol, wfn%qat(:, 1), wfn%dpat(:, :, 1))
   call shell_partition(mol, calc, wfn)
end subroutine sad_guess_qsh

subroutine eeq_guess_qat(mol, charges, dpat, error)
   type(structure_type), intent(in) :: mol
   real(wp), intent(inout) :: charges(:), dpat(:, :)
   type(error_type), intent(out), allocatable :: error

   dpat(:, :) = 0.0_wp
   call get_eeq_charges(mol, error, charges)
end subroutine eeq_guess_qat

subroutine eeq_guess_qsh(mol, calc, wfn, error)
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(inout) :: wfn
   type(error_type), intent(out), allocatable :: error

   call eeq_guess_qat(mol, wfn%qat(:, 1), wfn%dpat(:, :, 1), error)
   call shell_partition(mol, calc, wfn)
end subroutine eeq_guess_qsh

subroutine eeqbc_guess_qat(mol, charges, dpat, error)
   type(structure_type), intent(in) :: mol
   real(wp), intent(inout) :: charges(:), dpat(:, :)
   type(error_type), intent(out), allocatable :: error

   dpat(:, :) = 0.0_wp
   call get_eeqbc_charges(mol, error, charges)
end subroutine eeqbc_guess_qat

subroutine eeqbc_guess_qsh(mol, calc, wfn, error)
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(inout) :: wfn
   type(error_type), intent(out), allocatable :: error

   call eeqbc_guess_qat(mol, wfn%qat(:, 1), wfn%dpat(:, :, 1), error)
   call shell_partition(mol, calc, wfn)
end subroutine eeqbc_guess_qsh

subroutine shell_partition(mol, calc, wfn)
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(inout) :: wfn

   integer :: iat, ii, ish, spin

   call get_occupation(mol, calc%bas, calc%h0, wfn%nocc, wfn%n0at, wfn%n0sh)
   do spin = 1, size(wfn%qat, 2)
      do iat = 1, size(wfn%qat, 1)
         ii = calc%bas%ish_at(iat)
         do ish = 1, calc%bas%nsh_at(iat)
            wfn%qsh(ii+ish, spin) = (wfn%n0sh(ii+ish) / wfn%n0at(iat)) * wfn%qat(iat, spin)
         end do
      end do
   end do
end subroutine shell_partition

end module tblite_wavefunction_guess
