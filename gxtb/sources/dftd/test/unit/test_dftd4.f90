! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module test_dftd4
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use dftd4
   implicit none
   private

   public :: collect_dftd4

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_dftd4(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("PBE-D4", test_pbed4_mb01), &
      & new_unittest("PBE-D4-EEQBC", test_pbed4_eeqbc_mb01), &
      & new_unittest("PBE-D4S", test_pbed4s_mb01), &
      & new_unittest("PBE-D4S-EEQBC", test_pbed4s_eeqbc_mb01), &
      & new_unittest("B97-D4", test_b97d4_mb02), &
      & new_unittest("B97-D4S", test_b97d4s_mb02), &
      & new_unittest("TPSS-D4", test_tpssd4_mb03), &
      & new_unittest("TPSS-D4S", test_tpssd4s_mb03), &
      & new_unittest("PWPB95-D4", test_pwpb95d4_mb04), &
      & new_unittest("PWPB95-D4S", test_pwpb95d4s_mb04), &
      & new_unittest("B2PLYP-D4", test_b2plypd4_mb05), &
      & new_unittest("B2PLYP-D4S", test_b2plypd4s_mb05), &
      & new_unittest("PW6B95-D4", test_pw6b95d4_mb06), &
      & new_unittest("PW6B95-D4S", test_pw6b95d4s_mb06), &
      & new_unittest("OLYP-D4", test_olypd4_mb07), &
      & new_unittest("OLYP-D4S", test_olypd4s_mb07), &
      & new_unittest("PBE0-D4", test_pbe0d4_mb08), &
      & new_unittest("PBE0-D4S", test_pbe0d4s_mb08), &
      & new_unittest("RPBE-D4-ATM", test_rpbed4atm_mb09), &
      & new_unittest("RPBE-D4S-ATM", test_rpbed4satm_mb09), &
      & new_unittest("B2GPPLYP-D4-ATM", test_b2gpplypd4atm_mb10), &
      & new_unittest("B2GPPLYP-D4S-ATM", test_b2gpplypd4satm_mb10), &
      & new_unittest("LH14t-calPBE-D4-ATM", test_lh14tcalpbed4atm_mb11), &
      & new_unittest("LH14t-calPBE-D4S-ATM", test_lh14tcalpbed4satm_mb11), &
      & new_unittest("B1B95-D4-ATM", test_b1b95d4atm_mb12), &
      & new_unittest("B1B95-D4S-ATM", test_b1b95d4satm_mb12), &
      & new_unittest("M06L-D4-ATM", test_m06ld4atm_mb13), &
      & new_unittest("M06L-D4S-ATM", test_m06ld4satm_mb13), &
      & new_unittest("TPSSh-D4-ATM", test_tpsshd4atm_mb14), &
      & new_unittest("TPSSh-D4S-ATM", test_tpsshd4satm_mb14), &
      & new_unittest("HF-D4-ATM", test_hfd4atm_mb15), &
      & new_unittest("HF-D4S-ATM", test_hfd4satm_mb15), &
      & new_unittest("CAM-B3LYP-D4-ATM", test_camb3lypd4atm_mb16), &
      & new_unittest("CAM-B3LYP-D4S-ATM", test_camb3lypd4satm_mb16), &
      & new_unittest("r2SCAN-3c", test_r2scan3c_mb01), &
      & new_unittest("r2SCAN-3c-D4S", test_r2scan3c_d4s_mb01), &
      & new_unittest("TPSSh-D4-ATM-AmF3", test_tpsshd4atm_amf3), &
      & new_unittest("TPSSh-D4S-ATM-AmF3", test_tpsshd4satm_amf3), &
      & new_unittest("Actinides-D4", test_actinides_d4), &
      & new_unittest("Actinides-D4S", test_actinides_d4s) &
      & ]

end subroutine collect_dftd4


subroutine test_dftd4_gen(error, mol, d4, damp, param, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Damping function
   type(damping_type), intent(in) :: damp
   
   !> Damping parameters
   type(param_type), intent(in) :: param

   !> Expected dispersion energy
   real(wp), intent(in) :: ref

   real(wp) :: energy

   call get_dispersion(mol, d4, damp, param, realspace_cutoff(), energy)

   call check(error, energy, ref, thr=thr)
   if (allocated(error)) then
      call test_failed(error, "Dispersion energy does not match")
      print*,energy
   end if

end subroutine test_dftd4_gen


subroutine test_numgrad(error, mol, d4, damp, param)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Damping function
   type(damping_type), intent(in) :: damp

   !> Damping parameters
   type(param_type), intent(in) :: param

   integer :: iat, ic
   real(wp) :: energy, er, el, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))

   do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_dispersion(mol, d4, damp, param, realspace_cutoff(), er)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_dispersion(mol, d4, damp, param, realspace_cutoff(), el)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(er - el)/step
      end do
   end do

   call get_dispersion(mol, d4, damp, param, realspace_cutoff(), energy, &
      & gradient, sigma)

   if (any(abs(gradient - numgrad) > thr2)) then
      call test_failed(error, "Gradient of dispersion energy does not match")
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_numgrad


subroutine test_numsigma(error, mol, d4, damp, param)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Dispersion model
   class(dispersion_model), intent(in) :: d4

   !> Damping function
   type(damping_type), intent(in) :: damp

   !> Damping parameters
   type(param_type), intent(in) :: param

   integer :: ic, jc
   real(wp) :: energy, er, el, sigma(3, 3), eps(3, 3), numsigma(3, 3)
   real(wp), allocatable :: gradient(:, :), xyz(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), xyz(3, mol%nat))

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   do ic = 1, 3
      do jc = 1, 3
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call get_dispersion(mol, d4, damp, param, realspace_cutoff(), er)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call get_dispersion(mol, d4, damp, param, realspace_cutoff(), el)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         numsigma(jc, ic) = 0.5_wp*(er - el)/step
      end do
   end do

   call get_dispersion(mol, d4, damp, param, realspace_cutoff(), energy, &
      & gradient, sigma)

   if (any(abs(sigma - numsigma) > thr2)) then
      call test_failed(error, "Strain derivatives do not match")
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_numsigma


subroutine test_pbed4_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4, damp, param, -1.8578752883363366E-002_wp)

end subroutine test_pbed4_mb01

subroutine test_pbed4_eeqbc_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol, qmod=d4_qmod%eeqbc)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4, damp, param, -1.8374771323104909E-002_wp)

end subroutine test_pbed4_eeqbc_mb01

subroutine test_pbed4s_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4s, damp, param, -1.9870451633183694E-002_wp)

end subroutine test_pbed4s_mb01

subroutine test_pbed4s_eeqbc_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4s_model(error, d4s, mol, qmod=d4_qmod%eeqbc)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4s, damp, param, -1.9652490628871744E-002_wp)

end subroutine test_pbed4s_eeqbc_mb01

subroutine test_b97d4_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.69460052_wp, a1 = 0.28904684_wp, a2 = 4.13407323_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4, damp, param, -8.9181168937810723E-002_wp)

end subroutine test_b97d4_mb02

subroutine test_b97d4s_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.69460052_wp, a1 = 0.28904684_wp, a2 = 4.13407323_wp)

   call get_structure(mol, "MB16-43", "02")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4s, damp, param, -0.10354699373838475_wp)

end subroutine test_b97d4s_mb02

subroutine test_tpssd4_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.91130849_wp, a1 = 0.43332851_wp, a2 = 4.56986797_wp)

   call get_structure(mol, "MB16-43", "03")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4, damp, param, -2.4695638764787930E-002_wp)

end subroutine test_tpssd4_mb03

subroutine test_tpssd4s_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.91130849_wp, a1 = 0.43332851_wp, a2 = 4.56986797_wp)

   call get_structure(mol, "MB16-43", "03")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4s, damp, param, -2.9321970435591389E-002_wp)

end subroutine test_tpssd4s_mb03

subroutine test_pwpb95d4_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 0.82_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = -0.34639127_wp, a1 = 0.41080636_wp, a2 = 3.83878274_wp)

   call get_structure(mol, "MB16-43", "04")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4, damp, param, -9.5128100471706181E-003_wp)

end subroutine test_pwpb95d4_mb04

subroutine test_pwpb95d4s_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 0.82_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = -0.34639127_wp, a1 = 0.41080636_wp, a2 = 3.83878274_wp)

   call get_structure(mol, "MB16-43", "04")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4s, damp, param, -1.0003412660121589E-002_wp)

end subroutine test_pwpb95d4s_mb04

subroutine test_b2plypd4_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 0.64_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.15117773_wp, a1 = 0.42666167_wp, a2 = 4.73635790_wp)

   call get_structure(mol, "MB16-43", "05")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_numgrad(error, mol, d4, damp, param)

end subroutine test_b2plypd4_mb05

subroutine test_b2plypd4s_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 0.64_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.15117773_wp, a1 = 0.42666167_wp, a2 = 4.73635790_wp)

   call get_structure(mol, "MB16-43", "05")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_numgrad(error, mol, d4s, damp, param)

end subroutine test_b2plypd4s_mb05

subroutine test_pw6b95d4_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = -0.31629935_wp, a1 = 0.03999357_wp, a2 = 5.83690254_wp)

   call get_structure(mol, "MB16-43", "06")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_numgrad(error, mol, d4, damp, param)

end subroutine test_pw6b95d4_mb06

subroutine test_pw6b95d4s_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = -0.31629935_wp, a1 = 0.03999357_wp, a2 = 5.83690254_wp)

   call get_structure(mol, "MB16-43", "06")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_numgrad(error, mol, d4s, damp, param)

end subroutine test_pw6b95d4s_mb06

subroutine test_olypd4_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 2.74836820_wp, a1 = 0.60184498_wp, a2 = 2.53292167_wp)

   call get_structure(mol, "MB16-43", "07")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_numgrad(error, mol, d4, damp, param)

end subroutine test_olypd4_mb07

subroutine test_olypd4s_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 2.74836820_wp, a1 = 0.60184498_wp, a2 = 2.53292167_wp)

   call get_structure(mol, "MB16-43", "07")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_numgrad(error, mol, d4s, damp, param)

end subroutine test_olypd4s_mb07


subroutine test_pbe0d4_mb08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.20065498_wp, a1 = 0.40085597_wp, a2 = 5.02928789_wp)

   call get_structure(mol, "MB16-43", "08")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_numsigma(error, mol, d4, damp, param)

end subroutine test_pbe0d4_mb08

subroutine test_pbe0d4s_mb08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.20065498_wp, a1 = 0.40085597_wp, a2 = 5.02928789_wp)

   call get_structure(mol, "MB16-43", "08")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_numsigma(error, mol, d4s, damp, param)

end subroutine test_pbe0d4s_mb08


subroutine test_rpbed4atm_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.31183787_wp, a1 = 0.46169493_wp, a2 = 3.15711757_wp)

   call get_structure(mol, "MB16-43", "09")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4, damp, param, -4.5140422485299259E-002_wp)

end subroutine test_rpbed4atm_mb09

subroutine test_rpbed4satm_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.31183787_wp, a1 = 0.46169493_wp, a2 = 3.15711757_wp)

   call get_structure(mol, "MB16-43", "09")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4s, damp, param, -4.7626123029052128E-002_wp)

end subroutine test_rpbed4satm_mb09


subroutine test_b2gpplypd4atm_mb10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 0.56_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 0.94633372_wp, a1 = 0.42907301_wp, a2 = 5.18802602_wp)

   call get_structure(mol, "MB16-43", "10")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4, damp, param, -9.6812427202205668E-003_wp)

end subroutine test_b2gpplypd4atm_mb10

subroutine test_b2gpplypd4satm_mb10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 0.56_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 0.94633372_wp, a1 = 0.42907301_wp, a2 = 5.18802602_wp)

   call get_structure(mol, "MB16-43", "10")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4s, damp, param, -1.1490632914927183E-002_wp)

end subroutine test_b2gpplypd4satm_mb10

subroutine test_lh14tcalpbed4atm_mb11(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.27677253_wp, a1 = 0.38128670_wp, a2 = 4.91698883_wp)

   call get_structure(mol, "MB16-43", "11")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4, damp, param, -1.7460015867914524E-002_wp)

end subroutine test_lh14tcalpbed4atm_mb11

subroutine test_lh14tcalpbed4satm_mb11(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.27677253_wp, a1 = 0.38128670_wp, a2 = 4.91698883_wp)

   call get_structure(mol, "MB16-43", "11")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4s, damp, param, -2.0329508047867848E-002_wp)

end subroutine test_lh14tcalpbed4satm_mb11

subroutine test_b1b95d4atm_mb12(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.27701162_wp, a1 = 0.40554715_wp, a2 = 4.63323074_wp)

   call get_structure(mol, "MB16-43", "12")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4, damp, param, -2.5712178361964221E-002_wp)

end subroutine test_b1b95d4atm_mb12

subroutine test_b1b95d4satm_mb12(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.27701162_wp, a1 = 0.40554715_wp, a2 = 4.63323074_wp)

   call get_structure(mol, "MB16-43", "12")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4s, damp, param, -2.6528352555725724E-002_wp)

end subroutine test_b1b95d4satm_mb12


subroutine test_m06ld4atm_mb13(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 0.59493760_wp, a1 = 0.71422359_wp, a2 = 6.35314182_wp)

   call get_structure(mol, "MB16-43", "13")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_numgrad(error, mol, d4, damp, param)

end subroutine test_m06ld4atm_mb13

subroutine test_m06ld4satm_mb13(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 0.59493760_wp, a1 = 0.71422359_wp, a2 = 6.35314182_wp)

   call get_structure(mol, "MB16-43", "13")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_numgrad(error, mol, d4s, damp, param)

end subroutine test_m06ld4satm_mb13

subroutine test_tpsshd4atm_mb14(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.85897750_wp, a1 = 0.44286966_wp, a2 = 4.60230534_wp)

   call get_structure(mol, "MB16-43", "14")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_numgrad(error, mol, d4, damp, param)

end subroutine test_tpsshd4atm_mb14

subroutine test_tpsshd4satm_mb14(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.85897750_wp, a1 = 0.44286966_wp, a2 = 4.60230534_wp)

   call get_structure(mol, "MB16-43", "14")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_numgrad(error, mol, d4s, damp, param)

end subroutine test_tpsshd4satm_mb14

subroutine test_hfd4atm_mb15(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.61679827_wp, a1 = 0.44959224_wp, a2 = 3.35743605_wp)

   call get_structure(mol, "MB16-43", "15")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_numgrad(error, mol, d4, damp, param)

end subroutine test_hfd4atm_mb15

subroutine test_hfd4satm_mb15(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.61679827_wp, a1 = 0.44959224_wp, a2 = 3.35743605_wp)

   call get_structure(mol, "MB16-43", "15")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_numgrad(error, mol, d4s, damp, param)

end subroutine test_hfd4satm_mb15

subroutine test_camb3lypd4atm_mb16(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.74407961_wp, a1 = 0.40137870_wp, a2 = 5.18731225_wp)

   call get_structure(mol, "MB16-43", "16")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_numsigma(error, mol, d4, damp, param)

end subroutine test_camb3lypd4atm_mb16

subroutine test_camb3lypd4satm_mb16(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.74407961_wp, a1 = 0.40137870_wp, a2 = 5.18731225_wp)

   call get_structure(mol, "MB16-43", "16")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_numsigma(error, mol, d4s, damp, param)

end subroutine test_camb3lypd4satm_mb16

subroutine test_r2scan3c_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref = -5.7825025556386862E-003_wp

   param = param_type(&
      & s6 = 1.00_wp, s9 = 2.00_wp, alp = 16.0_wp, &
      & s8 = 0.00_wp, a1 = 0.42_wp, a2 = 5.65_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4_model(error, d4, mol, ga=2.0_wp, gc=1.0_wp)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4, damp, param, ref)

end subroutine test_r2scan3c_mb01

subroutine test_r2scan3c_d4s_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp
   
   real(wp), parameter :: ref = -6.1176284425895639E-003_wp

   param = param_type(&
      & s6 = 1.00_wp, s9 = 2.00_wp, alp = 16.0_wp, &
      & s8 = 0.00_wp, a1 = 0.42_wp, a2 = 5.65_wp)

   call get_structure(mol, "MB16-43", "01")
   call new_d4s_model(error, d4s, mol, ga=2.0_wp, gc=1.0_wp)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4s, damp, param, ref)

end subroutine test_r2scan3c_d4s_mb01


subroutine test_tpsshd4atm_amf3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref = -2.4882226918209061E-003_wp

   integer, parameter :: nat = 4
   integer, parameter :: num(nat) = [95, 9, 9, 9]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & -1.13163973200000_wp, -2.17446990100000_wp, +1.10012477100000_wp, &
      & -4.66377948900000_wp, -3.12947883400000_wp, -0.36987606800000_wp, &
      & -0.19032564300000_wp, +1.36339950600000_wp, -0.36521789300000_wp, &
      & +1.46283310800000_wp, -4.75734549200000_wp, -0.36503081000000_wp],&
      & [3, nat])

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.85897750_wp, a1 = 0.44286966_wp, a2 = 4.60230534_wp)

   call new(mol, num, xyz)
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return

   call test_dftd4_gen(error, mol, d4, damp, param, ref)
   if (allocated(error)) return

   call test_numgrad(error, mol, d4, damp, param)
   if (allocated(error)) return

   call test_numsigma(error, mol, d4, damp, param)

end subroutine test_tpsshd4atm_amf3

subroutine test_tpsshd4satm_amf3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref = -2.5523952752362009E-003_wp

   integer, parameter :: nat = 4
   integer, parameter :: num(nat) = [95, 9, 9, 9]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & -1.13163973200000_wp, -2.17446990100000_wp, +1.10012477100000_wp, &
      & -4.66377948900000_wp, -3.12947883400000_wp, -0.36987606800000_wp, &
      & -0.19032564300000_wp, +1.36339950600000_wp, -0.36521789300000_wp, &
      & +1.46283310800000_wp, -4.75734549200000_wp, -0.36503081000000_wp],&
      & [3, nat])

   param = param_type(&
      & s6 = 1.0_wp, s9 = 1.0_wp, alp = 16.0_wp, &
      & s8 = 1.85897750_wp, a1 = 0.44286966_wp, a2 = 4.60230534_wp)

   call new(mol, num, xyz)
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return

   call test_dftd4_gen(error, mol, d4s, damp, param, ref)
   if (allocated(error)) return

   call test_numgrad(error, mol, d4s, damp, param)
   if (allocated(error)) return

   call test_numsigma(error, mol, d4s, damp, param)

end subroutine test_tpsshd4satm_amf3

subroutine test_actinides_d4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref = -0.17966420554540324_wp
   
   integer, parameter :: nat = 17
   integer, parameter :: num(nat) = [&
      & 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
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
      & [3, nat])

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call new(mol, num, xyz)
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4, damp, param, ref)

end subroutine test_actinides_d4

subroutine test_actinides_d4s(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   real(wp), parameter :: ref = -0.18578133252612403_wp
   
   integer, parameter :: nat = 17
   integer, parameter :: num(nat) = [&
      & 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
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
      & [3, nat])

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call new(mol, num, xyz)
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4s, damp, param, ref)

end subroutine test_actinides_d4s


end module test_dftd4
