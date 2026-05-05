! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

module test_pairwise
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use dftd3
   implicit none
   private

   public :: collect_pairwise

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   type(realspace_cutoff), parameter :: cutoff = &
      & realspace_cutoff(cn=30_wp, disp2=60.0_wp, disp3=15.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_pairwise(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("PBE-D3(BJ)", test_pbed3_bj_mb01), &
      & new_unittest("B97-D3(BJ)-ATM", test_b97d3_bj_atm_mb02), &
      & new_unittest("TPSS-D3(BJ)-ATM", test_tpssd3_bj_atm_ammonia), &
      & new_unittest("HF-D3(0)", test_hfd3_zero_mb03), &
      & new_unittest("HSE06-D3(0)-ATM", test_hse06d3_zero_atm_mb04), &
      & new_unittest("BLYP-D3(0)-ATM", test_blypd3_zero_atm_urea), &
      & new_unittest("B3LYP-D3M(0)", test_b3lypd3_mzero_mb05), &
      & new_unittest("PBE0-D3M(0)-ATM", test_pbe0d3_mzero_atm_mb06), &
      & new_unittest("BP86-D3M(0)-ATM", test_bpd3_mzero_atm_co2), &
      & new_unittest("PBE0-D3(op)", test_revtpssd3_op_mb07), &
      & new_unittest("PBE0-D3(op)-ATM", test_revpbed3_op_atm_mb08) &
      & ]

end subroutine collect_pairwise


subroutine test_dftd3_pairwise(error, mol, param)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   type(d3_model) :: d3
   real(wp) :: energy
   real(wp), allocatable :: energy2(:, :), energy3(:, :)

   allocate(energy2(mol%nat, mol%nat), energy3(mol%nat, mol%nat))
   call new_d3_model(d3, mol)
   call get_dispersion(mol, d3, param, cutoff, energy)
   call get_pairwise_dispersion(mol, d3, param, cutoff, energy2, energy3)

   call check(error, energy, sum(energy2) + sum(energy3), thr=thr)
   if (allocated(error)) then
      print*,energy, sum(energy2) + sum(energy3)
   end if

end subroutine test_dftd3_pairwise


subroutine test_pbed3_bj_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.4289_wp, s8 = 0.7875_wp, a2 = 4.4407_wp)

   call get_structure(mol, "MB16-43", "01")
   call test_dftd3_pairwise(error, mol, param)

end subroutine test_pbed3_bj_mb01


subroutine test_b97d3_bj_atm_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.5545_wp, s8 = 2.2609_wp, a2 = 3.2297_wp)

   call get_structure(mol, "MB16-43", "02")
   call test_dftd3_pairwise(error, mol, param)

end subroutine test_b97d3_bj_atm_mb02


subroutine test_tpssd3_bj_atm_ammonia(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.4535_wp, s8 = 1.9435_wp, a2 = 4.4752_wp)

   call get_structure(mol, "X23", "ammonia")
   call test_dftd3_pairwise(error, mol, param)

end subroutine test_tpssd3_bj_atm_ammonia


subroutine test_hfd3_zero_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param = zero_damping_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, rs8 = 1.0_wp, alp = 14.0_wp, &
      & rs6 = 1.158_wp, s8 = 1.746_wp)

   call get_structure(mol, "MB16-43", "03")
   call test_dftd3_pairwise(error, mol, param)

end subroutine test_hfd3_zero_mb03


subroutine test_hse06d3_zero_atm_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param = zero_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, rs8 = 1.0_wp, alp = 14.0_wp, &
      & rs6 = 1.129_wp, s8 = 0.109_wp)

   call get_structure(mol, "MB16-43", "04")
   call test_dftd3_pairwise(error, mol, param)

end subroutine test_hse06d3_zero_atm_mb04


subroutine test_blypd3_zero_atm_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param = zero_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, rs8 = 1.0_wp, alp = 14.0_wp, &
      & rs6 = 1.094_wp, s8 = 1.682_wp)

   call get_structure(mol, "X23", "urea")
   call test_dftd3_pairwise(error, mol, param)

end subroutine test_blypd3_zero_atm_urea


subroutine test_b3lypd3_mzero_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(mzero_damping_param) :: param = mzero_damping_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, rs8 = 1.0_wp, alp = 14.0_wp, &
      & rs6 = 1.338153_wp, s8 = 1.532981_wp, bet = 0.013988_wp)

   call get_structure(mol, "MB16-43", "05")
   call test_dftd3_pairwise(error, mol, param)

end subroutine test_b3lypd3_mzero_mb05


subroutine test_pbe0d3_mzero_atm_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(mzero_damping_param) :: param = mzero_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, rs8 = 1.0_wp, alp = 14.0_wp, &
      & rs6 = 2.077949_wp, s8 = 0.000081_wp, bet = 0.116755_wp)

   call get_structure(mol, "MB16-43", "06")
   call test_dftd3_pairwise(error, mol, param)

end subroutine test_pbe0d3_mzero_atm_mb06


subroutine test_bpd3_mzero_atm_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(mzero_damping_param) :: param = mzero_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, rs8 = 1.0_wp, alp = 14.0_wp, &
      & rs6 = 1.233460_wp, s8 = 1.945174_wp, bet = 0.000000_wp)

   call get_structure(mol, "X23", "CO2")
   call test_dftd3_pairwise(error, mol, param)

end subroutine test_bpd3_mzero_atm_co2


subroutine test_revtpssd3_op_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(optimizedpower_damping_param) :: param = optimizedpower_damping_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & s8 = 0.27632_wp, a1 = 0.700_wp, a2 = 2.500_wp, bet = 8.0_wp)

   call get_structure(mol, "MB16-43", "07")
   call test_dftd3_pairwise(error, mol, param)

end subroutine test_revtpssd3_op_mb07


subroutine test_revpbed3_op_atm_mb08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(optimizedpower_damping_param) :: param = optimizedpower_damping_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & s8 = 1.44765_wp, a1 = 0.600_wp, a2 = 2.50_wp, bet = 0.0_wp)

   call get_structure(mol, "MB16-43", "08")
   call test_dftd3_pairwise(error, mol, param)

end subroutine test_revpbed3_op_atm_mb08


end module test_pairwise
