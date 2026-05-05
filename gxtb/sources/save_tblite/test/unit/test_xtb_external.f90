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

module test_xtb_external
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_container, only : container_type, container_cache
   use tblite_context_type, only : context_type
   use tblite_external_field, only : electric_field
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_xtb_external

   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: thr2 = 1e+4_wp*sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 3.166808578545117e-06_wp

  real(wp), parameter :: aatoau = 1.0_wp / 0.529177249_wp, &
     & ctoau = 1.0_wp / 1.60217653e-19_wp, jtoau = 1.0_wp / 4.3597441775e-18_wp
  !> Convert V/Å = J/(C·Å) to atomic units
  real(wp), parameter :: vatoau = jtoau / (ctoau * aatoau)

  type, extends(container_type) :: empty_interaction
  end type empty_interaction

contains


!> Collect all exported unit tests
subroutine collect_xtb_external(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gfn1-efield", test_e_mb01), &
      new_unittest("gfn2-efield", test_e_mb02), &
      new_unittest("gfn1-dipole", test_d_mb03), &
      new_unittest("gfn2-dipole", test_d_mb04), &
      new_unittest("gfn1-empty", test_g_mb05), &
      new_unittest("gfn2-empty", test_g_mb06) &
      ]

end subroutine collect_xtb_external


subroutine test_e_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   real(wp) :: energy
   real(wp), parameter :: ref1 = -33.092759826817499_wp, ref0 = -33.040345115781605_wp

   call get_structure(mol, "MB16-43", "01")
   energy = 0.0_wp

   call new_gfn1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)

   cont = electric_field([-2.0_wp, 0.0_wp, 0.0_wp]*vatoau)
   call calc%push_back(cont)

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

   call check(error, energy, ref1, thr=thr)
   if (allocated(error)) return

   call calc%pop(cont)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

   call check(error, energy, ref0, thr=thr)

end subroutine test_e_mb01


subroutine test_e_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   real(wp), parameter :: ref1 = -24.157134025160921_wp, ref0 = -24.069929678894123_wp

   call get_structure(mol, "MB16-43", "02")
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp

   call new_gfn2_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)

   cont = electric_field([0.0_wp, sqrt(2.0_wp), -sqrt(2.0_wp)]*vatoau)
   call calc%push_back(cont)

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, verbosity=0)

   call check(error, energy, ref1, thr=thr)
   if (allocated(error)) return

   call calc%pop(cont)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, verbosity=0)

   call check(error, energy, ref0, thr=thr)

end subroutine test_e_mb02


subroutine test_d_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn, wfn0
   class(container_type), allocatable :: cont
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp) :: energy, efield(3), er, el, numdip(3), dipole(3)
   integer :: i

   call get_structure(mol, "MB16-43", "03")
   energy = 0.0_wp
   efield(:) = 0.0_wp

   call new_gfn1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call new_wavefunction(wfn0, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)

   cont = electric_field(efield)
   call calc%push_back(cont)

   call xtb_singlepoint(ctx, mol, calc, wfn0, acc, energy, verbosity=0)
   dipole(:) = matmul(mol%xyz, wfn0%qat(:, 1)) + sum(wfn0%dpat(:, :, 1), 2)

   do i = 1, 3
      wfn = wfn0
      efield(i) = efield(i) + step
      call calc%pop(cont)
      cont = electric_field(efield)
      call calc%push_back(cont)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, er, verbosity=0)

      wfn = wfn0
      efield(i) = efield(i) - 2*step
      call calc%pop(cont)
      cont = electric_field(efield)
      call calc%push_back(cont)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, el, verbosity=0)

      efield(i) = efield(i) + step
      numdip(i) = -0.5_wp * (er - el) / step
   end do

   if (any(abs(dipole - numdip) > thr2)) then
      call test_failed(error, "Numerical dipole moment does not match")
      print '(3es21.14)', dipole
      print '("---")'
      print '(3es21.14)', numdip
      print '("---")'
      print '(3es21.14)', dipole - numdip
   end if

end subroutine test_d_mb03


subroutine test_d_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn, wfn0
   class(container_type), allocatable :: cont
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp) :: energy, efield(3), er, el, numdip(3), dipole(3)
   integer :: i

   call get_structure(mol, "MB16-43", "04")
   energy = 0.0_wp
   efield(:) = 0.0_wp

   call new_gfn2_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call new_wavefunction(wfn0, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)

   cont = electric_field(efield)
   call calc%push_back(cont)

   call xtb_singlepoint(ctx, mol, calc, wfn0, acc, energy, verbosity=0)
   dipole(:) = matmul(mol%xyz, wfn0%qat(:, 1)) + sum(wfn0%dpat(:, :, 1), 2)

   do i = 1, 3
      wfn = wfn0
      efield(i) = efield(i) + step
      call calc%pop(cont)
      cont = electric_field(efield)
      call calc%push_back(cont)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, er, verbosity=0)

      wfn = wfn0
      efield(i) = efield(i) - 2*step
      call calc%pop(cont)
      cont = electric_field(efield)
      call calc%push_back(cont)
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, el, verbosity=0)

      efield(i) = efield(i) + step
      numdip(i) = -0.5_wp * (er - el) / step
   end do

   if (any(abs(dipole - numdip) > thr2)) then
      call test_failed(error, "Numerical dipole moment does not match")
      print '(3es21.14)', dipole
      print '("---")'
      print '(3es21.14)', numdip
      print '("---")'
      print '(3es21.14)', dipole - numdip
   end if

end subroutine test_d_mb04


subroutine test_g_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   real(wp), parameter :: eref = -27.735987624622005_wp, gref(3, 16) = reshape([&
      & -1.82353787492949E-3_wp,  1.31818459922823E-3_wp,  8.58226495653576E-3_wp, &
      &  1.97688739182230E-2_wp, -7.29458073763057E-3_wp, -6.04380003974638E-3_wp, &
      &  3.98810014714115E-3_wp,  6.14906826920453E-3_wp,  4.80384729197490E-3_wp, &
      &  9.11662446246199E-3_wp,  8.11888967395009E-4_wp, -2.06033497236509E-3_wp, &
      &  6.93297807104095E-3_wp, -9.99636455500038E-4_wp, -5.74433211393415E-3_wp, &
      & -4.86754354610952E-3_wp, -3.06270172653747E-3_wp, -2.01091116475248E-3_wp, &
      &  4.16440905040429E-3_wp, -2.03405181467356E-3_wp, -6.07808932105647E-4_wp, &
      & -2.20381696875807E-2_wp, -1.11474305059426E-4_wp,  7.88311942639446E-3_wp, &
      & -1.13525078678607E-3_wp, -9.03126254289776E-3_wp,  5.41440972040784E-3_wp, &
      & -4.31763236474853E-4_wp,  7.28641048680264E-3_wp,  1.46388718395281E-3_wp, &
      & -9.52395151098753E-3_wp, -3.06196332824414E-4_wp, -1.18009384848677E-2_wp, &
      & -3.17217740137865E-3_wp,  1.23412949417049E-2_wp, -6.93057638998609E-3_wp, &
      &  3.84703556856836E-3_wp, -5.30965159378415E-3_wp, -4.76649943206558E-3_wp, &
      &  4.01813275719477E-3_wp, -6.60996461366484E-3_wp, -4.94608420633252E-3_wp, &
      &  4.66628771990608E-3_wp,  2.54601905375790E-3_wp,  1.53173370986307E-3_wp, &
      & -1.35100476506938E-2_wp,  4.30665380447907E-3_wp,  1.52320234470266E-2_wp],&
      & shape(gref))

   call get_structure(mol, "MB16-43", "05")
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_gfn2_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)

   cont = empty_interaction()
   call calc%push_back(cont)

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, verbosity=0)

   call check(error, energy, eref, thr=thr)
   if (allocated(error)) return
   call check(error, all(abs(gradient - gref) < thr))

end subroutine test_g_mb05


subroutine test_g_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   real(wp), parameter :: eref = -18.559531670207463_wp, gref(3, 16) = reshape([&
      & -7.68855208823154E-3_wp,  2.10638867503445E-2_wp,  2.41215155941334E-3_wp, &
      &  2.90805965762717E-3_wp, -1.75018015056621E-3_wp,  3.85582760443651E-3_wp, &
      &  1.54219205661951E-3_wp,  1.32376377402360E-3_wp,  1.01824888645767E-3_wp, &
      & -1.91452656544107E-2_wp, -2.23384392341465E-2_wp,  5.48586109290690E-3_wp, &
      & -4.39746104836314E-3_wp, -1.10529660508857E-2_wp, -2.03434329217776E-2_wp, &
      &  9.62864183290068E-3_wp,  9.33495782979855E-3_wp,  4.28513096166196E-3_wp, &
      &  3.40301832758327E-3_wp,  1.33326298386767E-2_wp,  3.81296500108361E-4_wp, &
      & -2.32159630655584E-3_wp,  8.79100677206436E-4_wp, -2.03120821307283E-3_wp, &
      &  1.73423639495923E-2_wp,  1.52735549282184E-2_wp, -5.01999960210040E-3_wp, &
      & -3.35073564349097E-3_wp, -1.59744035347341E-2_wp,  8.31528670874199E-3_wp, &
      & -1.82578281101229E-3_wp, -3.86682679585475E-3_wp, -1.47462197318047E-3_wp, &
      &  1.33591887148489E-3_wp, -9.88764121832669E-4_wp,  6.43834617770471E-3_wp, &
      & -2.84651079600129E-3_wp,  1.62990353769770E-3_wp,  1.10448375328198E-3_wp, &
      &  4.00339318812389E-3_wp, -1.03144794325501E-3_wp,  4.05047193028005E-3_wp, &
      &  1.72567953487094E-3_wp, -6.60586563316507E-3_wp, -3.98074228816423E-3_wp, &
      & -3.13363070736831E-4_wp,  7.71096128474042E-4_wp, -4.49710017669782E-3_wp],&
      & shape(gref))

   call get_structure(mol, "MB16-43", "06")
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_gfn2_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)

   cont = empty_interaction()
   call calc%push_back(cont)

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, verbosity=0)

   call check(error, energy, eref, thr=thr)
   if (allocated(error)) return
   call check(error, all(abs(gradient - gref) < thr))

end subroutine test_g_mb06


end module test_xtb_external
