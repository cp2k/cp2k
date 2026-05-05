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

module test_mixer
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_context_type, only : context_type
   use tblite_scf_iterator, only : new_iterator
   use tblite_scf_mixer, only : new_mixer
   use tblite_scf_mixer_broyden, only : broyden_input
   use tblite_scf_mixer_diis, only : diis_input
   use tblite_scf_mixer_input, only : mixer_mode, mixer_input_container
   use tblite_scf_mixer_simple, only : simple_input
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, eeq_guess, sad_guess, eeqbc_guess
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_ipea1, only : new_ipea1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_mixer

   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = 10*sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 3.166808578545117e-06_wp

contains


!> Collect all exported unit tests
subroutine collect_mixer(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("convergence-gfn1-default", test_gfn1_m01_default), &
      new_unittest("convergence-gfn1-diis-density", test_gfn1_m01_diis_density), &
      new_unittest("convergence-gfn1-diis-potential", test_gfn1_m01_diis_potential), &
      new_unittest("convergence-ipea1-default", test_ipea1_m01_default), &
      new_unittest("convergence-ipea1-diis-density", test_ipea1_m01_diis_density), &
      new_unittest("convergence-ipea1-diis-potential", test_ipea1_m01_diis_potential), &
      new_unittest("convergence-gfn2-default", test_gfn2_m01_default), &
      new_unittest("convergence-gfn2-diis-density", test_gfn2_m01_diis_density), &
      new_unittest("convergence-gfn2-diis-potential", test_gfn2_m01_diis_potential) &
      ]

end subroutine collect_mixer

subroutine test_convergence(error, mol, calc, ref, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> xTB Calculator
   type(xtb_calculator), intent(in) :: calc

   !> Reference energy
   real(wp), intent(in) :: ref

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   type(context_type) :: ctx
   type(wavefunction_type) :: wfn
   integer :: nspin
   real(wp) :: energy
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   if (mol%uhf > 0) then
      nspin = 2
   else
      nspin = 1
   end if

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & nspin, calc%default_etemp * kt)
   call sad_guess(mol, calc, wfn)
   if (allocated(error)) return
   energy = 0.0_wp
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

   call check(error, energy, ref, thr=thr_)

end subroutine test_convergence


subroutine test_gfn1_m01_default(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   real(wp), parameter :: ref = -33.04034511578_wp

   call get_structure(mol, "MB16-43", "01")

   call new_gfn1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   call test_convergence(error, mol, calc, ref, thr_in=thr1)

end subroutine test_gfn1_m01_default

subroutine test_gfn1_m01_diis_density(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(mixer_input_container), allocatable :: mixer_inputs(:)
   real(wp), parameter :: ref = -33.04034511578_wp

   call get_structure(mol, "MB16-43", "01")

   call new_gfn1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   call new_mixer(calc%iterator%mixers(2)%raw, diis_input(mode=mixer_mode%density, &
      & output_fraction=0.5_wp, damp=0.5_wp, precollect=2, start=3))

   call test_convergence(error, mol, calc, ref, thr_in=thr1)

end subroutine test_gfn1_m01_diis_density

subroutine test_gfn1_m01_diis_potential(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(mixer_input_container), allocatable :: mixer_inputs(:)
   real(wp), parameter :: ref = -33.04034511578_wp

   call get_structure(mol, "MB16-43", "01")

   call new_gfn1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   call new_mixer(calc%iterator%mixers(1)%raw, simple_input(mode=mixer_mode%potential, &
      & damp=0.4_wp, start=1))
   call new_mixer(calc%iterator%mixers(2)%raw, diis_input(mode=mixer_mode%potential, &
      & output_fraction=0.5_wp, damp=1.0_wp, precollect=2, start=4))

   call test_convergence(error, mol, calc, ref, thr_in=thr1)

end subroutine test_gfn1_m01_diis_potential

subroutine test_ipea1_m01_default(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   real(wp), parameter :: ref = -35.41368160407_wp

   call get_structure(mol, "MB16-43", "01")

   call new_ipea1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   call test_convergence(error, mol, calc, ref, thr_in=thr1)

end subroutine test_ipea1_m01_default

subroutine test_ipea1_m01_diis_density(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(mixer_input_container), allocatable :: mixer_inputs(:)
   real(wp), parameter :: ref = -35.41368160407_wp

   call get_structure(mol, "MB16-43", "01")

   call new_ipea1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   call new_mixer(calc%iterator%mixers(2)%raw, diis_input(mode=mixer_mode%density, &
      & output_fraction=1.0_wp, damp=0.9_wp, precollect=3, start=5))

   call test_convergence(error, mol, calc, ref, thr_in=thr1)

end subroutine test_ipea1_m01_diis_density

subroutine test_ipea1_m01_diis_potential(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(mixer_input_container), allocatable :: mixer_inputs(:)
   real(wp), parameter :: ref = -35.41368160407_wp

   call get_structure(mol, "MB16-43", "01")

   call new_ipea1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   call new_mixer(calc%iterator%mixers(1)%raw, simple_input(mode=mixer_mode%potential, &
      & damp=0.1_wp, start=1))
   call new_mixer(calc%iterator%mixers(2)%raw, diis_input(mode=mixer_mode%potential, &
      & output_fraction=0.5_wp, damp=1.0_wp, precollect=2, start=4))

   call test_convergence(error, mol, calc, ref, thr_in=thr1)

end subroutine test_ipea1_m01_diis_potential

subroutine test_gfn2_m01_default(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   real(wp), parameter :: ref = -30.34890234779_wp

   call get_structure(mol, "MB16-43", "01")

   call new_gfn2_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   call test_convergence(error, mol, calc, ref, thr_in=thr1)

end subroutine test_gfn2_m01_default

subroutine test_gfn2_m01_diis_density(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(mixer_input_container), allocatable :: mixer_inputs(:)
   real(wp), parameter :: ref = -30.34890234779_wp

   call get_structure(mol, "MB16-43", "01")

   call new_gfn2_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   call new_mixer(calc%iterator%mixers(2)%raw, diis_input(mode=mixer_mode%density, &
      & output_fraction=0.7_wp, damp=1.0_wp, precollect=2, start=4))

   call test_convergence(error, mol, calc, ref, thr_in=thr1)

end subroutine test_gfn2_m01_diis_density

subroutine test_gfn2_m01_diis_potential(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(mixer_input_container), allocatable :: mixer_inputs(:)
   real(wp), parameter :: ref = -30.34890234779_wp

   call get_structure(mol, "MB16-43", "01")

   call new_gfn2_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   call new_mixer(calc%iterator%mixers(1)%raw, simple_input(mode=mixer_mode%potential, &
      & damp=0.1_wp, start=1))
   call new_mixer(calc%iterator%mixers(2)%raw, diis_input(mode=mixer_mode%potential, &
      & output_fraction=1.0_wp, damp=1.0_wp, precollect=2, start=4))

   call test_convergence(error, mol, calc, ref, thr_in=thr1)

end subroutine test_gfn2_m01_diis_potential

end module test_mixer
