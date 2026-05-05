! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

module test_dftd3
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use dftd3
   implicit none
   private

   public :: collect_dftd3

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_dftd3(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("PBE-D3(BJ)", test_pbed3bj_mb01), &
      & new_unittest("B97-D3(BJ)", test_b97d3bj_mb02), &
      & new_unittest("TPSS-D3(BJ)", test_tpssd3bj_mb03), &
      & new_unittest("PWPB95-D3(BJ)", test_pwpb95d3bj_mb04), &
      & new_unittest("B2PLYP-D3(BJ)", test_b2plypd3bj_mb05), &
      & new_unittest("PW6B95-D3(BJ)", test_pw6b95d3bj_mb06), &
      & new_unittest("OLYP-D3(BJ)", test_olypd3bj_mb07), &
      & new_unittest("PBE0-D3(BJ)", test_pbe0d3bj_mb08), &
      & new_unittest("RPBE-D3(0)", test_rpbed3zero_mb09), &
      & new_unittest("B2GPPLYP-D3(0)", test_b2gpplypd3zero_mb10), &
      & new_unittest("SSB-D3(0)", test_ssbd3zero_mb11), &
      & new_unittest("B1B95-D3(0)", test_b1b95d3zero_mb12), &
      & new_unittest("M06L-D3(0)", test_m06ld3zero_mb13), &
      & new_unittest("TPSSh-D3(0)", test_tpsshd3zero_mb14), &
      & new_unittest("HF-D3(0)", test_hfd3zero_mb15), &
      & new_unittest("CAM-B3LYP-D3(0)", test_camb3lypd3zero_mb16), &
      & new_unittest("DSD-BLYP-D3(BJ)-ATM", test_dsdblypd3bjatm_mb17), &
      & new_unittest("mPWLYP-D3(BJ)-ATM", test_mpwlypd3bjatm_mb18), &
      & new_unittest("TPSS0-D3(BJ)-ATM", test_tpss0d3bjatm_mb19), &
      & new_unittest("BPBE-D3(BJ)-ATM", test_bpbed3bjatm_mb20), &
      & new_unittest("PWB6K-D3(BJ)-ATM", test_pwb6kd3bjatm_mb21), &
      & new_unittest("lc-wPBE-D3(BJ)-ATM", test_lcwpbed3bjatm_mb22), &
      & new_unittest("PW1PW-D3(BJ)-ATM", test_pw1pwd3bjatm_mb23), &
      & new_unittest("mPW1B95-D3(BJ)-ATM", test_mpw1b95d3bjatm_mb24), &
      & new_unittest("BLYP-D3(0)-ATM", test_blypd3zeroatm_mb25), &
      & new_unittest("revPBE-D3(0)-ATM", test_revpbed3zeroatm_mb26), &
      & new_unittest("HSE06-D3(0)-ATM", test_hse06d3zeroatm_mb27), &
      & new_unittest("BMK-D3(0)-ATM", test_bmkd3zeroatm_mb28), &
      & new_unittest("M06-D3(0)-ATM", test_m06d3zeroatm_mb29), &
      & new_unittest("B3LYP-D3(0)-ATM", test_b3lypd3zeroatm_mb30), &
      & new_unittest("BP-D3(0)-ATM", test_bpd3zeroatm_mb31), &
      & new_unittest("M05-D3(0)-ATM", test_m05d3zeroatm_mb32), &
      & new_unittest("PBE-D3(0M)", test_pbed3zerom_mb33), &
      & new_unittest("BLYP-D3(0M)", test_blypd3zerom_mb34), &
      & new_unittest("B97-D3(0M)", test_b97dd3zerom_mb35), &
      & new_unittest("lc-wPBE-D3(0M)", test_lcwpbed3zerom_mb36), &
      & new_unittest("B97h-D3(op)", test_b97hd3op_mb37), &
      & new_unittest("TPSSh-D3(op)", test_tpsshd3op_mb38), &
      & new_unittest("B3LYP-D3(CSO)", test_b3lypd3cso_mb01), &
      & new_unittest("PBE-D3(CSO)-ATM", test_pbed3cso_mb02), &
      & new_unittest("PBE-D3(BJ) Actinides", test_pbed3bj_actinides) &
      & ]

end subroutine collect_dftd3


subroutine test_dftd3_gen(error, mol, param, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Expected dispersion energy
   real(wp), intent(in) :: ref

   type(d3_model) :: d3
   real(wp) :: energy

   call new_d3_model(d3, mol)
   call get_dispersion(mol, d3, param, realspace_cutoff(), energy)

   call check(error, energy, ref, thr=thr)
   if (allocated(error)) then
      print*,energy
   end if

end subroutine test_dftd3_gen


subroutine test_numgrad(error, mol, param)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   integer :: iat, ic
   type(d3_model) :: d3
   real(wp) :: energy, er, el, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   call new_d3_model(d3, mol)

   do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_dispersion(mol, d3, param, realspace_cutoff(), er)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_dispersion(mol, d3, param, realspace_cutoff(), el)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(er - el)/step
      end do
   end do

   call get_dispersion(mol, d3, param, realspace_cutoff(), energy, gradient, sigma)

   if (any(abs(gradient - numgrad) > thr2)) then
      call test_failed(error, "Gradient of dispersion energy does not match")
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_numgrad


subroutine test_numsigma(error, mol, param)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   integer :: ic, jc
   type(d3_model) :: d3
   real(wp) :: energy, er, el, sigma(3, 3), eps(3, 3), numsigma(3, 3)
   real(wp), allocatable :: gradient(:, :), xyz(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), xyz(3, mol%nat))
   call new_d3_model(d3, mol)

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   do ic = 1, 3
      do jc = 1, 3
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call get_dispersion(mol, d3, param, realspace_cutoff(), er)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call get_dispersion(mol, d3, param, realspace_cutoff(), el)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         numsigma(jc, ic) = 0.5_wp*(er - el)/step
      end do
   end do

   call get_dispersion(mol, d3, param, realspace_cutoff(), energy, gradient, sigma)

   if (any(abs(sigma - numsigma) > thr2)) then
      call test_failed(error, "Strain derivatives do not match")
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_numsigma



subroutine test_pbed3bj_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.4289_wp, s8 = 0.7875_wp, a2 = 4.4407_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -1.7882220155186028E-002_wp)

end subroutine test_pbed3bj_mb01


subroutine test_b97d3bj_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.5545_wp, s8 = 2.2609_wp, a2 = 3.2297_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -5.9503371130070551E-002_wp)

end subroutine test_b97d3bj_mb02


subroutine test_tpssd3bj_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.4535_wp, s8 = 1.9435_wp, a2 = 4.4752_wp)

   call get_structure(mol, "MB16-43", "03")
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -2.8684059240803347E-002_wp)

end subroutine test_tpssd3bj_mb03


subroutine test_pwpb95d3bj_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 0.82_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.0000_wp, s8 = 0.2904_wp, a2 = 7.3141_wp)

   call get_structure(mol, "MB16-43", "04")
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -1.2270948604225146E-002_wp)

end subroutine test_pwpb95d3bj_mb04


subroutine test_b2plypd3bj_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 0.64_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.3065_wp, s8 = 0.9147_wp, a2 = 5.0570_wp)

   call get_structure(mol, "MB16-43", "05")
   call new_rational_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_b2plypd3bj_mb05


subroutine test_pw6b95d3bj_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.2076_wp, s8 = 0.7257_wp, a2 = 6.3750_wp)

   call get_structure(mol, "MB16-43", "06")
   call new_rational_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_pw6b95d3bj_mb06


subroutine test_olypd3bj_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.5299_wp, s8 = 2.6205_wp, a2 = 2.8065_wp)

   call get_structure(mol, "MB16-43", "07")
   call new_rational_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_olypd3bj_mb07


subroutine test_pbe0d3bj_mb08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.4145_wp, s8 = 1.2177_wp, a2 = 4.8593_wp)

   call get_structure(mol, "MB16-43", "08")
   call new_rational_damping(param, inp)
   call test_numsigma(error, mol, param)

end subroutine test_pbe0d3bj_mb08


subroutine test_rpbed3zero_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 0.872_wp, s8 = 0.514_wp)

   call get_structure(mol, "MB16-43", "09")
   call new_zero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -2.0178760785797962E-002_wp)

end subroutine test_rpbed3zero_mb09


subroutine test_b2gpplypd3zero_mb10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 0.56_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.586_wp, s8 = 0.760_wp)

   call get_structure(mol, "MB16-43", "10")
   call new_zero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -8.2668192717285059E-003_wp)

end subroutine test_b2gpplypd3zero_mb10


subroutine test_ssbd3zero_mb11(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.215_wp, s8 = 0.663_wp)

   call get_structure(mol, "MB16-43", "11")
   call new_zero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -8.7923023966613584E-003_wp)

end subroutine test_ssbd3zero_mb11


subroutine test_b1b95d3zero_mb12(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.613_wp, s8 = 1.868_wp)

   call get_structure(mol, "MB16-43", "12")
   call new_zero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -1.5146271685761267E-002_wp)

end subroutine test_b1b95d3zero_mb12


subroutine test_m06ld3zero_mb13(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.581_wp, s8 = 0.000_wp)

   call get_structure(mol, "MB16-43", "13")
   call new_zero_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_m06ld3zero_mb13


subroutine test_tpsshd3zero_mb14(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.223_wp, s8 = 1.219_wp)

   call get_structure(mol, "MB16-43", "14")
   call new_zero_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_tpsshd3zero_mb14


subroutine test_hfd3zero_mb15(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.158_wp, s8 = 1.746_wp)

   call get_structure(mol, "MB16-43", "15")
   call new_zero_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_hfd3zero_mb15


subroutine test_camb3lypd3zero_mb16(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.378_wp, s8 = 1.217_wp)

   call get_structure(mol, "MB16-43", "16")
   call new_zero_damping(param, inp)
   call test_numsigma(error, mol, param)

end subroutine test_camb3lypd3zero_mb16


subroutine test_dsdblypd3bjatm_mb17(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 0.5_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.0_wp, s8 = 0.2130_wp, a2 = 6.0519_wp)

   call get_structure(mol, "MB16-43", "17")
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -1.3592755832923201E-002_wp)

end subroutine test_dsdblypd3bjatm_mb17


subroutine test_mpwlypd3bjatm_mb18(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.4831_wp, s8 = 2.0077_wp, a2 = 3.5043_wp)

   call get_structure(mol, "MB16-43", "18")
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -4.4058898934399460E-002_wp)

end subroutine test_mpwlypd3bjatm_mb18


subroutine test_tpss0d3bjatm_mb19(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.3768_wp, s8 = 1.2576_wp, a2 = 4.5865_wp)

   call get_structure(mol, "MB16-43", "19")
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -4.1318505521837870E-002_wp)

end subroutine test_tpss0d3bjatm_mb19


subroutine test_bpbed3bjatm_mb20(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.4567_wp, s8 = 4.0728_wp, a2 = 4.3908_wp)

   call get_structure(mol, "MB16-43", "20")
   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -3.4953087919973355E-002_wp)

end subroutine test_bpbed3bjatm_mb20


subroutine test_pwb6kd3bjatm_mb21(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.1865_wp, s8 = 0.9383_wp, a2 = 7.7627_wp)

   call get_structure(mol, "MB16-43", "21")
   call new_rational_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_pwb6kd3bjatm_mb21


subroutine test_lcwpbed3bjatm_mb22(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.3919_wp, s8 = 1.8541_wp, a2 = 5.0897_wp)

   call get_structure(mol, "MB16-43", "22")
   call new_rational_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_lcwpbed3bjatm_mb22


subroutine test_pw1pwd3bjatm_mb23(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.3807_wp, s8 = 2.3363_wp, a2 = 5.8844_wp)

   call get_structure(mol, "MB16-43", "23")
   call new_rational_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_pw1pwd3bjatm_mb23


subroutine test_mpw1b95d3bjatm_mb24(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.1955_wp, s8 = 1.0508_wp, a2 = 6.4177_wp)

   call get_structure(mol, "MB16-43", "24")
   call new_rational_damping(param, inp)
   call test_numsigma(error, mol, param)

end subroutine test_mpw1b95d3bjatm_mb24


subroutine test_blypd3zeroatm_mb25(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.094_wp, s8 = 1.682_wp)

   call get_structure(mol, "MB16-43", "25")
   call new_zero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -1.8876259384422459E-002_wp)

end subroutine test_blypd3zeroatm_mb25


subroutine test_revpbed3zeroatm_mb26(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 0.923_wp, s8 = 1.010_wp)

   call get_structure(mol, "MB16-43", "26")
   call new_zero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -2.1820118878237430E-002_wp)

end subroutine test_revpbed3zeroatm_mb26


subroutine test_hse06d3zeroatm_mb27(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.129_wp, s8 = 0.109_wp)

   call get_structure(mol, "MB16-43", "27")
   call new_zero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -6.1849950573474831E-003_wp)

end subroutine test_hse06d3zeroatm_mb27


subroutine test_bmkd3zeroatm_mb28(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.931_wp, s8 = 2.168_wp)

   call get_structure(mol, "MB16-43", "28")
   call new_zero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -1.5629608155684784E-002_wp)

end subroutine test_bmkd3zeroatm_mb28


subroutine test_m06d3zeroatm_mb29(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.325_wp, s8 = 0.000_wp)

   call get_structure(mol, "MB16-43", "29")
   call new_zero_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_m06d3zeroatm_mb29


subroutine test_b3lypd3zeroatm_mb30(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.261_wp, s8 = 1.105_wp)

   call get_structure(mol, "MB16-43", "30")
   call new_zero_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_b3lypd3zeroatm_mb30


subroutine test_bpd3zeroatm_mb31(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.139_wp, s8 = 1.683_wp)

   call get_structure(mol, "MB16-43", "31")
   call new_zero_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_bpd3zeroatm_mb31


subroutine test_m05d3zeroatm_mb32(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.355_wp, s8 = 1.279_wp)

   call get_structure(mol, "MB16-43", "32")
   call new_zero_damping(param, inp)
   call test_numsigma(error, mol, param)

end subroutine test_m05d3zeroatm_mb32


subroutine test_pbed3zerom_mb33(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(mzero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 2.340218_wp, s8 = 0.0_wp, bet = 0.129434_wp)

   call get_structure(mol, "MB16-43", "33")
   call new_mzero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -6.4839784232922359E-002_wp)

end subroutine test_pbed3zerom_mb33


subroutine test_blypd3zerom_mb34(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(mzero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.279637_wp, s8 = 1.841686_wp, bet = 0.014370_wp)

   call get_structure(mol, "MB16-43", "34")
   call new_mzero_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -4.1835874470780833E-002_wp)

end subroutine test_blypd3zerom_mb34


subroutine test_b97dd3zerom_mb35(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(mzero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.338153_wp, s8 = 1.532981_wp, bet = 0.013988_wp)

   call get_structure(mol, "MB16-43", "35")
   call new_mzero_damping(param, inp)
   call test_numgrad(error, mol, param)

end subroutine test_b97dd3zerom_mb35


subroutine test_lcwpbed3zerom_mb36(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(mzero_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, rs8 = 1.0_wp, &
      & rs6 = 1.366361_wp, s8 = 1.280619_wp, bet = 0.003160_wp)

   call get_structure(mol, "MB16-43", "36")
   call new_mzero_damping(param, inp)
   call test_numsigma(error, mol, param)

end subroutine test_lcwpbed3zerom_mb36


subroutine test_b97hd3op_mb37(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(optimizedpower_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 0.97388_wp, s9 = 1.0_wp, alp = 14.0_wp, bet = 6.0_wp, &
      & a1 = 0.150_wp, s8 = 0.0_wp, a2 = 4.25_wp)

   call get_structure(mol, "MB16-43", "37")
   call new_optimizedpower_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -2.7952861273999385E-002_wp)
   if (allocated(error)) return
   call test_numgrad(error, mol, param)

end subroutine test_b97hd3op_mb37


subroutine test_tpsshd3op_mb38(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(optimizedpower_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, bet = 8.0_wp, &
      & a1 = 0.575_wp, s8 = 0.43185_wp, a2 = 3.00_wp)

   call get_structure(mol, "MB16-43", "38")
   call new_optimizedpower_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -9.0101761599590859E-003_wp)
   if (allocated(error)) return
   call test_numsigma(error, mol, param)

end subroutine test_tpsshd3op_mb38


subroutine test_b3lypd3cso_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(cso_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.86_wp, a2 = 2.5_wp, rs6 = 0.0_wp, rs8 = 6.25_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_cso_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -3.8002950817452329E-002_wp)
   if (allocated(error)) return
   call test_numgrad(error, mol, param)

end subroutine test_b3lypd3cso_mb01


subroutine test_pbed3cso_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(cso_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 14.0_wp, &
      & a1 = 0.24_wp, a2 = 2.5_wp, rs6 = 0.0_wp, rs8 = 6.25_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_cso_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -3.8229908827789316E-002_wp)
   if (allocated(error)) return
   call test_numsigma(error, mol, param)

end subroutine test_pbed3cso_mb02


subroutine test_pbed3bj_actinides(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp = d3_param(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 14.0_wp, &
      & a1 = 0.4289_wp, s8 = 0.7875_wp, a2 = 4.4407_wp)

   ! Molecular structure data 
   mol%nat = 17
   mol%nid = 17
   mol%id = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, &
      & 12, 13, 14, 15, 16, 17]
   mol%num = [87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, &
      & 98, 99, 100, 101, 102, 103]
   mol%xyz = reshape([ &
      & 0.98692316414074_wp, 6.12727238368797_wp,-6.67861597188102_wp, &
      & 3.63898862390869_wp, 5.12109301182962_wp, 3.01908613326278_wp, &
      & 5.14503571563551_wp,-3.97172984617710_wp, 3.82011791828867_wp, &
      & 6.71986847575494_wp, 1.71382138402812_wp, 3.92749159076307_wp, &
      & 4.13783589704826_wp,-2.10695793491818_wp, 0.19753203068899_wp, &
      & 8.97685097698326_wp,-3.08813636191844_wp,-4.45568615593938_wp, &
      & 12.5486412940776_wp,-1.77128765259458_wp, 0.59261498922861_wp, &
      & 7.82051475868325_wp,-3.97159756604558_wp,-0.53637703616916_wp, &
      &-0.43444574624893_wp,-1.69696511583960_wp,-1.65898182093050_wp, &
      &-4.71270645149099_wp,-0.11534827468942_wp, 2.84863373521297_wp, &
      &-2.52061680335614_wp, 1.82937752749537_wp,-2.10366982879172_wp, &
      & 0.13551154616576_wp, 7.99805359235043_wp,-1.55508522619903_wp, &
      & 3.91594542499717_wp,-1.72975169129597_wp,-5.07944366756113_wp, &
      &-1.03393930231679_wp, 4.69307230054046_wp, 0.02656940927472_wp, &
      & 6.20675384557240_wp, 4.24490721493632_wp,-0.71004195169885_wp, &
      & 7.04586341131562_wp, 5.20053667939076_wp,-7.51972863675876_wp, &
      & 2.01082807362334_wp, 1.34838807211157_wp,-4.70482633508447_wp],&
      & [3, 17])
   mol%periodic = [.false.]
   allocate(mol%lattice(0, 0))

   call new_rational_damping(param, inp)
   call test_dftd3_gen(error, mol, param, -1.4131143363689097E-001_wp)

end subroutine test_pbed3bj_actinides

end module test_dftd3
