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

module test_gxtb
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use multicharge, only : get_eeqbc_charges
   use tblite_context_type, only : context_type
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gxtb, only : new_gxtb_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_gxtb

   real(wp), parameter :: acc = 0.00001_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 3.166808578545117e-06_wp

contains


!> Collect all exported unit tests
subroutine collect_gxtb(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-h2", test_e_h2), &
      new_unittest("energy-lih", test_e_lih), &
      new_unittest("energy-water", test_e_water), &
      new_unittest("energy-no", test_e_no), &
      new_unittest("energy-chex-rad", test_e_chex_rad), &
      new_unittest("energy-s2", test_e_s2), &
      new_unittest("energy-cecl3", test_e_cecl3), &
      new_unittest("energy-ce2", test_e_ce2), &
      new_unittest("energy-mb01", test_e_mb01), &
      new_unittest("energy-mb02", test_e_mb02), &
      new_unittest("gradient-h2", test_g_h2), &
      new_unittest("gradient-lih", test_g_lih), &
      new_unittest("gradient-water", test_g_water), &
      new_unittest("gradient-no", test_g_no), &
      new_unittest("gradient-s2", test_g_s2), &
      new_unittest("gradient-pcl", test_g_pcl), &
      new_unittest("gradient-sih4", test_g_sih4), &
      new_unittest("gradient-cecl3", test_g_cecl3), &
      new_unittest("gradient-ce2", test_g_ce2), &
      new_unittest("gradient-panp", test_g_panp), &
      new_unittest("gradient-mb02", test_g_mb02) &
      !new_unittest("sigma-lih", test_s_lih), &
      !new_unittest("sigma-mb03", test_s_mb03) &
      ]

end subroutine collect_gxtb


subroutine test_e_gen(error, mol, ref, thr_in, guess_qat, guess_qsh, guess_density, guess_dp, guess_qp)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reference g-xTB energy
   real(wp), intent(in) :: ref
   !> Test threshold
   real(wp), intent(in), optional :: thr_in
   !> Optional guess for atomic charges
   real(wp), intent(in), optional :: guess_qat(:, :)
   !> Optional guess for shell charges
   real(wp), intent(in), optional :: guess_qsh(:, :)
   !> Optional guess for density matrix
   real(wp), intent(in), optional :: guess_density(:, :, :)
   !> Optional guess for dipole moments
   real(wp), intent(in), optional :: guess_dp(:, :, :)
   !> Optional guess for quadrupole moments
   real(wp), intent(in), optional :: guess_qp(:, :, :)

   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn, wfn_aux
   real(wp) :: energy
   integer :: nspin
   real(wp) :: thr_

   thr_ = thr1
   if (present(thr_in)) thr_ = thr_in

   energy = 0.0_wp

   nspin = 1
   if (mol%uhf > 0) then
      nspin = 2
   end if

   ! Setup g-xTB calculator
   call new_gxtb_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   ! Obtain EEQBC charges for charge adaptation
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return

   ! Setup wavefunction 
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & nspin, calc%default_etemp * kt, .false.)

   ! Use EEQBC as a guess
   if (present(guess_qat) .and. present(guess_qsh) .and. present(guess_density)) then
      wfn%qat(:, :) = guess_qat(:, :)
      wfn%qsh(:, :) = guess_qsh(:, :)
      wfn%dpat(:, :, :) = guess_dp(:, :, :)
      wfn%qpat(:, :, :) = guess_qp(:, :, :)
      wfn%density(:, :, :) = guess_density(:, :, :)
   else
      wfn%qat(:, 1) = wfn_aux%qat(:, 1)
   end if

   ! Perform g-xTB calculation
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0, wfn_aux=wfn_aux)

   ! Check result
   call check(error, energy, ref, thr=thr_)
   if (allocated(error)) then
      call test_failed(error, "Energy does not match reference within threshold")
      print *, ref, energy, ref - energy
      return
   end if

end subroutine test_e_gen


subroutine test_numgrad(error, mol, thr_in, four_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   !> Logical to activate four point stencil
   logical, intent(in), optional :: four_in

   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn, wfn_aux
   integer :: iat, ic, nspin
   real(wp) :: energy, er, el, er2, el2, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   logical :: four_point
   real(wp) :: thr_, step

   thr_ = thr1
   if (present(thr_in)) thr_ = thr_in

   four_point = .false.
   if (present(four_in)) four_point = four_in

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   nspin = 1
   if (mod(mol%uhf, 2) == 1) then
      nspin = 2
   end if

   ! Setup g-xTB calculator
   call new_gxtb_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   ! Obtain EEQBC charges for charge adaptation
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return

   ! Setup wavefunction 
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & nspin, calc%default_etemp * kt, .false.)

   ! Use EEQBC as a guess
   wfn%qat(:, 1) = wfn_aux%qat(:, 1)

   if (.not. four_point) then
      step = 5.0e-5_wp
      do iat = 1, mol%nat
         do ic = 1, 3
            er = 0.0_wp
            el = 0.0_wp

            ! Right hand
            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
            ! Update the EEQBC charges
            call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
            if (allocated(error)) return
            wfn%qat(:, 1) = wfn_aux%qat(:, 1)
            ! Perform g-xTB single point calculation
            call xtb_singlepoint(ctx, mol, calc, wfn, acc, er, verbosity=0, &
               & wfn_aux=wfn_aux)

            ! Left hand
            mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
            ! Update the EEQBC charges
            call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
            if (allocated(error)) return
            wfn%qat(:, 1) = wfn_aux%qat(:, 1)
            ! Perform g-xTB single point calculation
            call xtb_singlepoint(ctx, mol, calc, wfn, acc, el, verbosity=0, &
               & wfn_aux=wfn_aux)

            mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
            numgrad(ic, iat) = 0.5_wp*(er - el)/step
         end do
      end do
   else
      step = 2.0e-3_wp
      do iat = 1, mol%nat
         do ic = 1, 3
            ! f(x + 2h)
            mol%xyz(ic, iat) = mol%xyz(ic, iat) + 2*step
            ! Update the EEQBC charges
            call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
            if (allocated(error)) return
            wfn%qat(:, 1) = wfn_aux%qat(:, 1)
            ! Perform g-xTB single point calculation
            call xtb_singlepoint(ctx, mol, calc, wfn, acc, er2, verbosity=0, &
               & wfn_aux=wfn_aux)

            ! f(x + h)
            mol%xyz(ic, iat) = mol%xyz(ic, iat) - step
            ! Update the EEQBC charges
            call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
            if (allocated(error)) return
            wfn%qat(:, 1) = wfn_aux%qat(:, 1)
            ! Perform g-xTB single point calculation
            call xtb_singlepoint(ctx, mol, calc, wfn, acc, er, verbosity=0, &
               & wfn_aux=wfn_aux)

            ! f(x - h)
            mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
            ! Update the EEQBC charges
            call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
            if (allocated(error)) return
            wfn%qat(:, 1) = wfn_aux%qat(:, 1)
            ! Perform g-xTB single point calculation
            call xtb_singlepoint(ctx, mol, calc, wfn, acc, el, verbosity=0, &
               & wfn_aux=wfn_aux)

            ! f(x - 2h)
            mol%xyz(ic, iat) = mol%xyz(ic, iat) - step
            ! Update the EEQBC charges
            call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
            if (allocated(error)) return
            wfn%qat(:, 1) = wfn_aux%qat(:, 1)
            ! Perform g-xTB single point calculation
            call xtb_singlepoint(ctx, mol, calc, wfn, acc, el2, verbosity=0, &
               & wfn_aux=wfn_aux)

            mol%xyz(ic, iat) = mol%xyz(ic, iat) + 2*step
            numgrad(ic, iat) = (-er2 + 8.0_wp*er - 8.0_wp*el + el2) / (12.0_wp*step)
         end do
      end do
   end if

   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .true.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1), wfn_aux%dqatdr(:, :, :, 1), &
      & wfn_aux%dqatdL(:, :, :, 1))
   if (allocated(error)) return

   ! Perform g-xTB analytical gradient calculation
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient=gradient, &
      & sigma=sigma, verbosity=0, wfn_aux=wfn_aux)

   if (any(abs(gradient - numgrad) > thr_)) then
      call test_failed(error, "Gradient of g-xTB energy does not match")
      write(*,*) "numerical gradient:"
      print'(3es21.14)', numgrad
      write(*,*) "analytical gradient:"
      print'(3es21.14)', gradient
      write(*,*) "difference:"
      print'(3es21.14)', numgrad-gradient
   end if

end subroutine test_numgrad


subroutine test_e_h2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol
   ! Reference energy calculated with g-xTB development code
   real(wp), parameter :: ref = -1.16302850380698_wp

   ! Fully converged values from the g-xTB development code
   real(wp), parameter :: qat(2, 1) = reshape([&
      &-6.66133814775094E-16_wp, 4.44089209850063E-16_wp], shape(qat))

   real(wp), parameter :: qsh(2, 1) = reshape([&
      &-6.66133814775094E-16_wp, 4.44089209850063E-16_wp], shape(qsh))
      
   real(wp), parameter :: dp(3, 2, 1) = reshape([&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,-2.85449063469693E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 2.85449063469693E-01_wp],&
      & shape(dp))

   real(wp), parameter :: qp(6, 2, 1) = reshape([&
      & 1.23307601331135E-01_wp, 0.00000000000000E+00_wp, 1.23307601331135E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,-2.46615202662270E-01_wp, &
      & 1.23307601331135E-01_wp, 0.00000000000000E+00_wp, 1.23307601331135E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,-2.46615202662270E-01_wp],&
      & shape(qp))

   real(wp), parameter :: density(2, 2, 1) = reshape([&
      & 5.93683766916993E-01_wp, 5.93683766916992E-01_wp, 5.93683766916992E-01_wp, &
      & 5.93683766916992E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "H2")
   call test_e_gen(error, mol, ref, guess_qat=qat, guess_qsh=qsh, &
      & guess_density=density, guess_dp=dp, guess_qp=qp)

end subroutine test_e_h2

subroutine test_e_lih(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol
   ! Reference energy calculated with g-xTB development code
   real(wp), parameter :: ref = -8.06914465311961_wp

   ! Fully converged values from the g-xTB development code
   real(wp), parameter :: qat(2, 1) = reshape([&
      & 3.90304067136331E-01_wp,-3.90304067184860E-01_wp], shape(qat))

   real(wp), parameter :: qsh(3, 1) = reshape([&
      & 1.88324695751579E-01_wp, 2.01979371384751E-01_wp,-3.90304067184860E-01_wp], &
      & shape(qsh))
      
   real(wp), parameter :: dp(3, 2, 1) = reshape([&
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,-1.45652592802741E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 8.24364077964541E-02_wp],&
      & shape(dp))

   real(wp), parameter :: qp(6, 2, 1) = reshape([&
      & 1.71404874585399E+00_wp, 0.00000000000000E+00_wp, 1.71404874585399E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,-3.42809749170798E+00_wp, &
      &-1.17983418798047E-01_wp, 0.00000000000000E+00_wp,-1.17983418798047E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 2.35966837596095E-01_wp],&
      & shape(qp))

   real(wp), parameter :: density(5, 5, 1) = reshape([&
      & 7.43135813030401E-02_wp, 0.00000000000000E+00_wp, 1.15037971862012E-01_wp, &
      & 0.00000000000000E+00_wp, 2.77066832001789E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.15037971862012E-01_wp, 0.00000000000000E+00_wp, &
      & 1.78079628758031E-01_wp, 0.00000000000000E+00_wp, 4.28901499091321E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 2.77066832001789E-01_wp, &
      & 0.00000000000000E+00_wp, 4.28901499091321E-01_wp, 0.00000000000000E+00_wp, &
      & 1.03300134443079E+00_wp], shape(density))

   call get_structure(mol, "MB16-43", "LiH")
   call test_e_gen(error, mol, ref, guess_qat=qat, guess_qsh=qsh, &
      & guess_density=density, guess_dp=dp, guess_qp=qp)

end subroutine test_e_lih

subroutine test_e_water(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol
   ! Reference energy calculated with g-xTB development code
   real(wp), parameter :: ref = -76.4373874264764_wp 

   ! Fully converged values from the g-xTB development code
   real(wp), parameter :: qat(3, 1) = reshape([&
      & 3.14505586823761E-01_wp, 3.14503636947414E-01_wp,-6.29009223821014E-01_wp], &
      & shape(qat))

   real(wp), parameter :: qsh(4, 1) = reshape([&
      & 3.14505586823761E-01_wp, 3.14503636947414E-01_wp,-1.33128876658769E-01_wp, &
      &-4.95880347162244E-01_wp], shape(qsh))

   real(wp), parameter :: dp(3, 3, 1) = reshape([&
      & 7.43760389653602E-03_wp, 1.58888911895518E-16_wp,-3.64734942988183E-02_wp, &
      & 3.68734220821483E-02_wp,-2.58872250415853E-16_wp, 5.13562031559125E-03_wp, &
      & 1.03871485779780E-01_wp, 4.73065563308248E-20_wp,-7.34826380849913E-02_wp],&
      & shape(dp))

   real(wp), parameter :: qp(6, 3, 1) = reshape([&
      &-7.23750258344499E-02_wp,-1.07897772988824E-17_wp,-9.76720678380104E-02_wp, &
      & 6.47260769724112E-03_wp, 2.48224251346313E-16_wp, 1.70047093672460E-01_wp, &
      & 1.39095605920762E-01_wp, 3.87206514992289E-16_wp,-9.76753929921613E-02_wp, &
      & 8.11703026504840E-02_wp, 1.18123346164607E-16_wp,-4.14202129286005E-02_wp, &
      & 2.55501303889414E-02_wp, 2.24766134170308E-16_wp,-7.50517241196755E-02_wp, &
      & 3.38609130034574E-02_wp,-1.08229979047532E-16_wp, 4.95015937307342E-02_wp],&
      & shape(qp))
      
   real(wp), parameter :: density(6, 6, 1) = reshape([&
      & 4.00153978367091E-01_wp,-9.25481137319279E-02_wp, 3.73691113449640E-02_wp, &
      &-5.84123241887087E-16_wp,-6.25268345664895E-01_wp, 5.96674492392101E-02_wp, &
      &-9.25481137319279E-02_wp, 4.00158958706178E-01_wp, 3.73744967919657E-02_wp, &
      & 3.19295858595147E-16_wp, 1.51987133225997E-01_wp, 6.09447264285695E-01_wp, &
      & 3.73691113449640E-02_wp, 3.73744967919657E-02_wp, 1.77094623264867E+00_wp, &
      & 2.12432427866977E-16_wp, 3.16552287956886E-01_wp,-4.47523974626552E-01_wp, &
      &-5.84123241887087E-16_wp, 3.19295858595147E-16_wp, 2.12432427866977E-16_wp, &
      & 2.00000000000000E+00_wp,-3.94781182860990E-17_wp,-3.83146515704489E-16_wp, &
      &-6.25268345664895E-01_wp, 1.51987133225997E-01_wp, 3.16552287956886E-01_wp, &
      &-3.94781182860990E-17_wp, 1.05658003470220E+00_wp,-1.93369385497883E-01_wp, &
      & 5.96674492392101E-02_wp, 6.09447264285695E-01_wp,-4.47523974626552E-01_wp, &
      &-3.83146515704489E-16_wp,-1.93369385497883E-01_wp, 1.19318510023331E+00_wp],&
      & shape(density))

   call get_structure(mol, "ICE10", "gas")
   call test_e_gen(error, mol, ref, guess_qat=qat, guess_qsh=qsh, &
      & guess_density=density, guess_dp=dp, guess_qp=qp)

end subroutine test_e_water


subroutine test_e_no(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol
   ! Reference energy calculated with g-xTB development code
   real(wp), parameter :: ref = -129.893702527621_wp

   ! Fully converged from the g-xTB development code (tight convergence)
   real(wp), parameter :: qat(2, 2) = reshape([&
      & 1.12022036457523E-01_wp,-1.12022036557734E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp], shape(qat))

   real(wp), parameter :: qsh(4, 2) = reshape([&
      &-3.73400121557711E-01_wp, 4.85422158015234E-01_wp,-2.06420661593868E-01_wp, &
      & 9.43986250361339E-02_wp,-4.03682735064992E-03_wp, 6.72069840990768E-01_wp, &
      &-5.83211004786133E-03_wp, 3.37799096407743E-01_wp], shape(qsh))

   real(wp), parameter :: dp(3, 2, 2) = reshape([&
      & 3.98923026791198E-16_wp,-3.65511014840819E-16_wp, 7.13720793034213E-01_wp, &
      &-2.03524943525522E-16_wp, 2.11974094790780E-17_wp,-3.71745293836024E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp],&
      & shape(dp))

   real(wp), parameter :: qp(6, 2, 2) = reshape([&
      & 4.95203611294404E-03_wp,-3.41537083793069E-01_wp,-4.10411576464851E-01_wp, &
      &-1.10502615194510E-16_wp, 3.51306894422144E-16_wp, 4.05459540351908E-01_wp, &
      &-1.58143221148122E-01_wp,-9.91669659519196E-02_wp,-2.78746081513047E-01_wp, &
      & 1.03696212467558E-16_wp,-1.06720912738827E-16_wp, 4.36889302661170E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp],&
      & shape(qp))
      
   real(wp), parameter :: density(8, 8, 2) = reshape([&
      & 9.42019518198553E-01_wp,-6.70410354606576E-17_wp,-3.13829009362565E-01_wp, &
      & 2.12934722021574E-16_wp,-1.91142083958341E-01_wp,-4.81260324790394E-17_wp, &
      &-2.54599896675051E-02_wp,-2.79792496841097E-16_wp,-6.70410354606576E-17_wp, &
      & 8.87344818286031E-01_wp, 8.35353299149472E-17_wp, 3.56202052037748E-01_wp, &
      & 7.45157465640904E-17_wp,-1.45904439038172E-01_wp,-1.23063006538454E-17_wp, &
      &-2.89571919393296E-01_wp,-3.13829009362565E-01_wp, 8.35353299149472E-17_wp, &
      & 4.49237020657981E-01_wp,-1.49073003510408E-17_wp,-1.32690121489473E-03_wp, &
      &-2.92468370377978E-16_wp,-4.06976408233437E-01_wp,-1.67376946335903E-17_wp, &
      & 2.12934722021574E-16_wp, 3.56202052037748E-01_wp,-1.49073003510408E-17_wp, &
      & 4.54146263391490E-01_wp,-2.16312011776505E-17_wp,-2.89571919393295E-01_wp, &
      & 1.17820930332047E-17_wp, 2.06261238564159E-01_wp,-1.91142083958341E-01_wp, &
      & 7.45157465640904E-17_wp,-1.32690121489473E-03_wp,-2.16312011776505E-17_wp, &
      & 1.01793801587268E+00_wp, 6.85109393296945E-17_wp, 2.62429579626541E-01_wp, &
      &-9.80060047866476E-18_wp,-4.81260324790394E-17_wp,-1.45904439038172E-01_wp, &
      &-2.92468370377978E-16_wp,-2.89571919393295E-01_wp, 6.85109393296945E-17_wp, &
      & 9.55267322853733E-01_wp, 2.04809613232682E-18_wp, 2.35405427962641E-01_wp, &
      &-2.54599896675051E-02_wp,-1.23063006538454E-17_wp,-4.06976408233437E-01_wp, &
      & 1.17820930332047E-17_wp, 2.62429579626541E-01_wp, 2.04809613232682E-18_wp, &
      & 5.34554559386748E-01_wp,-7.55440264848595E-18_wp,-2.79792496841097E-16_wp, &
      &-2.89571919393296E-01_wp,-1.67376946335903E-17_wp, 2.06261238564159E-01_wp, &
      &-9.80060047866476E-18_wp, 2.35405427962641E-01_wp,-7.55440264848595E-18_wp, &
      & 6.68976745288669E-01_wp, 9.54493642978193E-01_wp, 1.60726885403408E-16_wp, &
      &-3.10425436698734E-01_wp, 3.16456449124610E-16_wp,-2.03716266983524E-01_wp, &
      &-2.25467194570894E-16_wp,-1.78381990679257E-02_wp,-2.45109398946379E-16_wp, &
      & 1.60726885403408E-16_wp, 2.44195332069810E-01_wp,-2.14807774282083E-16_wp, &
      &-3.05522867574033E-05_wp, 4.69399662581865E-17_wp, 3.65849585508874E-01_wp, &
      & 7.04266904153635E-18_wp,-1.06669780669403E-05_wp,-3.10425436698734E-01_wp, &
      &-2.14807774282083E-16_wp, 4.40523282932497E-01_wp,-1.38915448849928E-16_wp, &
      &-4.71240562857907E-03_wp,-1.27421506841548E-16_wp,-4.09806783985608E-01_wp, &
      & 1.36996215798866E-16_wp, 3.16456449124610E-16_wp,-3.05522867574033E-05_wp, &
      &-1.38915448849928E-16_wp, 2.44232488527490E-01_wp,-4.17545553420317E-16_wp, &
      &-1.06669780670643E-05_wp,-8.93591097772286E-18_wp, 3.65862558256945E-01_wp, &
      &-2.03716266983524E-01_wp, 4.69399662581865E-17_wp,-4.71240562857907E-03_wp, &
      &-4.17545553420317E-16_wp, 1.03061284626693E+00_wp, 1.00334623483500E-16_wp, &
      & 2.54769754979054E-01_wp, 3.30064951269781E-16_wp,-2.25467194570894E-16_wp, &
      & 3.65849585508874E-01_wp,-1.27421506841548E-16_wp,-1.06669780670643E-05_wp, &
      & 1.00334623483500E-16_wp, 5.48110069569196E-01_wp,-4.32899001663326E-17_wp, &
      & 3.66079823576893E-05_wp,-1.78381990679257E-02_wp, 7.04266904153635E-18_wp, &
      &-4.09806783985608E-01_wp,-8.93591097772286E-18_wp, 2.54769754979054E-01_wp, &
      &-4.32899001663326E-17_wp, 5.36711332400866E-01_wp, 3.21389716405877E-18_wp, &
      &-2.45109398946379E-16_wp,-1.06669780669403E-05_wp, 1.36996215798866E-16_wp, &
      & 3.65862558256945E-01_wp, 3.30064951269781E-16_wp, 3.66079823576893E-05_wp, &
      & 3.21389716405877E-18_wp, 5.48065548419011E-01_wp],&
      & shape(density))

   call get_structure(mol, "MB16-43", "NO")
   call test_e_gen(error, mol, ref, thr_in=thr1*10, guess_qat=qat, guess_qsh=qsh, &
      & guess_density=density, guess_dp=dp, guess_qp=qp)

end subroutine test_e_no


subroutine test_e_chex_rad(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol
   ! Reference energy calculated with g-xTB development code
   real(wp), parameter :: ref = -234.293753601443_wp 

   call get_structure(mol, "RC21", "8e")
   call test_e_gen(error, mol, ref, thr_in=thr1*10)

end subroutine test_e_chex_rad


subroutine test_e_s2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol
   ! Reference energy calculated with g-xTB development code
   real(wp), parameter :: ref = -796.322753612399_wp

   ! Fully converged values from the g-xTB development code
   real(wp), parameter :: qat(2, 1) = reshape([&
      &-4.23824725315214E-11_wp,-4.23819451755847E-11_wp], shape(qat))

   real(wp), parameter :: qsh(6, 1) = reshape([&
      &-1.65829471166743E-01_wp, 8.50037149362413E-02_wp, 8.08257561881188E-02_wp, &
      &-1.65829471166743E-01_wp, 8.50037149362421E-02_wp, 8.08257561881186E-02_wp],&
      & shape(qsh))
      
   real(wp), parameter :: dp(3, 2, 1) = reshape([&
      & 7.20390933021236E-16_wp,-1.08429022192345E-16_wp, 1.33456044924941E-01_wp, &
      &-1.99692331381165E-16_wp, 4.15387300445189E-16_wp,-1.33456044924939E-01_wp],&
      & shape(dp))

   real(wp), parameter :: qp(6, 2, 1) = reshape([&
      & 7.78513750067228E-01_wp, 2.80109502236468E-01_wp,-1.71390962399832E+00_wp, &
      & 4.15465547416762E-16_wp,-4.58804982640748E-16_wp, 9.35395873931098E-01_wp, &
      & 7.78513750067223E-01_wp, 2.80109502236468E-01_wp,-1.71390962399832E+00_wp, &
      & 7.54473207287948E-16_wp,-6.71563036244977E-16_wp, 9.35395873931100E-01_wp],&
      & shape(qp))

   real(wp), parameter :: density(18, 18, 1) = reshape([&
      & 2.01512275725746E+00_wp,-1.07774714878549E-17_wp,-3.99929207996199E-01_wp, &
      &-1.07496635420861E-16_wp,-3.67096317766298E-05_wp, 3.69391318119056E-18_wp, &
      &-3.51352576193404E-02_wp, 9.93512287536696E-18_wp,-1.63321743037866E-04_wp, &
      &-2.61461514520067E-01_wp, 1.28530688970932E-17_wp, 8.46120373097255E-02_wp, &
      &-5.21844262454934E-18_wp, 6.52391207996257E-06_wp, 1.06118614489111E-18_wp, &
      & 1.18677041290426E-02_wp,-2.50242368836971E-17_wp, 2.90249899222578E-05_wp, &
      &-1.07774714878549E-17_wp, 1.99493935742289E+00_wp,-2.19088386634370E-17_wp, &
      &-1.31464071393950E-01_wp,-3.64173280263401E-18_wp,-4.24834987435869E-02_wp, &
      & 1.86358830907554E-17_wp, 9.70517739656199E-03_wp,-5.61682172442967E-18_wp, &
      & 4.89202686255866E-18_wp,-3.83022248598425E-01_wp,-6.18616087262823E-17_wp, &
      & 1.32488995300194E-01_wp, 6.14828476880349E-18_wp,-1.51958552995893E-01_wp, &
      &-7.39258219650094E-17_wp, 1.18778354329892E-02_wp, 5.10782472349948E-18_wp, &
      &-3.99929207996199E-01_wp,-2.19088386634370E-17_wp, 8.16076839534090E-01_wp, &
      &-5.52547822483447E-18_wp, 6.02142834131760E-06_wp, 6.24158884175913E-19_wp, &
      & 5.39351517895339E-02_wp,-5.65731270603935E-17_wp, 2.67894316750492E-05_wp, &
      &-8.46120373097259E-02_wp, 5.12171500684398E-17_wp,-7.72403991855101E-01_wp, &
      &-2.99328951662798E-16_wp, 3.33878825672404E-08_wp, 3.64964487414511E-18_wp, &
      & 4.74250299828801E-02_wp, 8.36200608354454E-17_wp, 1.48543226707051E-07_wp, &
      &-1.07496635420861E-16_wp,-1.31464071393950E-01_wp,-5.52547822483447E-18_wp, &
      & 8.25167815308294E-01_wp, 1.52700642525368E-18_wp, 9.70517739656189E-03_wp, &
      & 2.52725824056651E-17_wp, 4.38734823638329E-02_wp,-6.16880272532067E-18_wp, &
      &-2.90756569844945E-17_wp, 1.32488995300194E-01_wp,-2.41054922702186E-16_wp, &
      & 7.95869099449793E-01_wp, 8.20376359559986E-18_wp, 1.18778354329889E-02_wp, &
      & 3.65818932892587E-18_wp,-4.62691906925265E-02_wp, 7.12070034072669E-17_wp, &
      &-3.67096317766298E-05_wp,-3.64173280263401E-18_wp, 6.02142834131760E-06_wp, &
      & 1.52700642525368E-18_wp, 6.72102565545150E-10_wp, 2.76110047099558E-19_wp, &
      & 5.71975865091167E-07_wp, 3.14157124137708E-20_wp, 2.99019513932858E-09_wp, &
      & 6.52391208026026E-06_wp, 3.37378592651583E-18_wp,-3.33878823276038E-08_wp, &
      & 7.46960457287212E-19_wp,-1.48925383029250E-10_wp, 2.97486784648590E-19_wp, &
      &-3.20635330261507E-07_wp,-9.49097733847636E-20_wp,-6.62571427782813E-10_wp, &
      & 3.69391318119056E-18_wp,-4.24834987435869E-02_wp, 6.24158884175913E-19_wp, &
      & 9.70517739656189E-03_wp, 2.76110047099558E-19_wp, 1.16687129118781E-02_wp, &
      &-3.66039855748310E-18_wp,-1.01852810953662E-03_wp, 1.01662808695120E-19_wp, &
      &-1.13357499977961E-18_wp, 1.51958552995893E-01_wp, 9.76554020426011E-18_wp, &
      &-1.18778354329889E-02_wp,-5.60779408961709E-20_wp, 4.23049731461357E-03_wp, &
      & 7.72673376676950E-18_wp,-7.46279735413830E-04_wp, 1.26845085436653E-18_wp, &
      &-3.51352576193404E-02_wp, 1.86358830907554E-17_wp, 5.39351517895339E-02_wp, &
      & 2.52725824056651E-17_wp, 5.71975865091167E-07_wp,-3.66039855748310E-18_wp, &
      & 3.73731480108770E-03_wp,-1.65005415356533E-18_wp, 2.54472983634185E-06_wp, &
      & 1.18677041290421E-02_wp,-4.35755627238332E-17_wp,-4.74250299828799E-02_wp, &
      & 1.48400353402036E-17_wp,-3.20635330245983E-07_wp,-1.56199455323153E-18_wp, &
      & 2.76687920574955E-03_wp, 3.55065882562968E-18_wp,-1.42651174861375E-06_wp, &
      & 9.93512287536696E-18_wp, 9.70517739656199E-03_wp,-5.65731270603935E-17_wp, &
      & 4.38734823638329E-02_wp, 3.14157124137708E-20_wp,-1.01852810953662E-03_wp, &
      &-1.65005415356533E-18_wp, 2.60581698338230E-03_wp,-3.88842350006346E-19_wp, &
      & 1.59208575921780E-17_wp,-1.18778354329890E-02_wp, 4.13691518309280E-17_wp, &
      & 4.62691906925265E-02_wp, 4.98794187554439E-19_wp,-7.46279735413839E-04_wp, &
      &-4.70316368471004E-18_wp,-2.40992381100263E-03_wp, 3.83847936960920E-18_wp, &
      &-1.63321743037866E-04_wp,-5.61682172442967E-18_wp, 2.67894316750492E-05_wp, &
      &-6.16880272532067E-18_wp, 2.99019513932858E-09_wp, 1.01662808695120E-19_wp, &
      & 2.54472983634185E-06_wp,-3.88842350006346E-19_wp, 1.33034263364430E-08_wp, &
      & 2.90249899220495E-05_wp, 7.20517863867866E-19_wp,-1.48543226693076E-07_wp, &
      &-6.87814222377941E-18_wp,-6.62571427756790E-10_wp, 4.16305480672137E-19_wp, &
      &-1.42651174861486E-06_wp, 3.32126706815553E-19_wp,-2.94779095388117E-09_wp, &
      &-2.61461514520067E-01_wp, 4.89202686255866E-18_wp,-8.46120373097259E-02_wp, &
      &-2.90756569844945E-17_wp, 6.52391208026026E-06_wp,-1.13357499977961E-18_wp, &
      & 1.18677041290421E-02_wp, 1.59208575921780E-17_wp, 2.90249899220495E-05_wp, &
      & 2.01512275725746E+00_wp,-1.27481751540811E-17_wp, 3.99929207996199E-01_wp, &
      & 1.76360118407346E-16_wp,-3.67096317767693E-05_wp,-2.07958009321102E-18_wp, &
      &-3.51352576193406E-02_wp,-7.97027096320685E-17_wp,-1.63321743037760E-04_wp, &
      & 1.28530688970932E-17_wp,-3.83022248598425E-01_wp, 5.12171500684398E-17_wp, &
      & 1.32488995300194E-01_wp, 3.37378592651583E-18_wp, 1.51958552995893E-01_wp, &
      &-4.35755627238332E-17_wp,-1.18778354329890E-02_wp, 7.20517863867866E-19_wp, &
      &-1.27481751540811E-17_wp, 1.99493935742289E+00_wp, 8.00541663295006E-17_wp, &
      &-1.31464071393950E-01_wp,-2.25179226325293E-20_wp, 4.24834987435870E-02_wp, &
      & 9.86239829642528E-17_wp,-9.70517739656185E-03_wp, 1.87223706676972E-17_wp, &
      & 8.46120373097255E-02_wp,-6.18616087262823E-17_wp,-7.72403991855101E-01_wp, &
      &-2.41054922702186E-16_wp,-3.33878823276038E-08_wp, 9.76554020426011E-18_wp, &
      &-4.74250299828799E-02_wp, 4.13691518309280E-17_wp,-1.48543226693076E-07_wp, &
      & 3.99929207996199E-01_wp, 8.00541663295006E-17_wp, 8.16076839534090E-01_wp, &
      & 4.25322705176787E-17_wp,-6.02142834153533E-06_wp, 2.90803512893118E-18_wp, &
      &-5.39351517895342E-02_wp,-7.71639069575584E-17_wp,-2.67894316750773E-05_wp, &
      &-5.21844262454934E-18_wp, 1.32488995300194E-01_wp,-2.99328951662798E-16_wp, &
      & 7.95869099449793E-01_wp, 7.46960457287212E-19_wp,-1.18778354329889E-02_wp, &
      & 1.48400353402036E-17_wp, 4.62691906925265E-02_wp,-6.87814222377941E-18_wp, &
      & 1.76360118407346E-16_wp,-1.31464071393950E-01_wp, 4.25322705176787E-17_wp, &
      & 8.25167815308295E-01_wp, 8.88552356043173E-18_wp,-9.70517739656187E-03_wp, &
      &-3.48680276562384E-17_wp,-4.38734823638328E-02_wp, 6.96815674075468E-17_wp, &
      & 6.52391207996257E-06_wp, 6.14828476880349E-18_wp, 3.33878825672404E-08_wp, &
      & 8.20376359559986E-18_wp,-1.48925383029250E-10_wp,-5.60779408961709E-20_wp, &
      &-3.20635330245983E-07_wp, 4.98794187554439E-19_wp,-6.62571427756790E-10_wp, &
      &-3.67096317767693E-05_wp,-2.25179226325293E-20_wp,-6.02142834153533E-06_wp, &
      & 8.88552356043173E-18_wp, 6.72102565549143E-10_wp,-4.48462041503179E-19_wp, &
      & 5.71975865103421E-07_wp,-4.41673869063451E-19_wp, 2.99019513933563E-09_wp, &
      & 1.06118614489111E-18_wp,-1.51958552995893E-01_wp, 3.64964487414511E-18_wp, &
      & 1.18778354329889E-02_wp, 2.97486784648590E-19_wp, 4.23049731461357E-03_wp, &
      &-1.56199455323153E-18_wp,-7.46279735413839E-04_wp, 4.16305480672137E-19_wp, &
      &-2.07958009321102E-18_wp, 4.24834987435870E-02_wp, 2.90803512893118E-18_wp, &
      &-9.70517739656187E-03_wp,-4.48462041503179E-19_wp, 1.16687129118781E-02_wp, &
      & 6.34175960452835E-18_wp,-1.01852810953664E-03_wp,-1.55049130225332E-19_wp, &
      & 1.18677041290426E-02_wp,-7.39258219650094E-17_wp, 4.74250299828801E-02_wp, &
      & 3.65818932892587E-18_wp,-3.20635330261507E-07_wp, 7.72673376676950E-18_wp, &
      & 2.76687920574955E-03_wp,-4.70316368471004E-18_wp,-1.42651174861486E-06_wp, &
      &-3.51352576193406E-02_wp, 9.86239829642528E-17_wp,-5.39351517895342E-02_wp, &
      &-3.48680276562384E-17_wp, 5.71975865103421E-07_wp, 6.34175960452835E-18_wp, &
      & 3.73731480108775E-03_wp, 5.48606965350553E-18_wp, 2.54472983634505E-06_wp, &
      &-2.50242368836971E-17_wp, 1.18778354329892E-02_wp, 8.36200608354454E-17_wp, &
      &-4.62691906925265E-02_wp,-9.49097733847636E-20_wp,-7.46279735413830E-04_wp, &
      & 3.55065882562968E-18_wp,-2.40992381100263E-03_wp, 3.32126706815553E-19_wp, &
      &-7.97027096320685E-17_wp,-9.70517739656185E-03_wp,-7.71639069575584E-17_wp, &
      &-4.38734823638328E-02_wp,-4.41673869063451E-19_wp,-1.01852810953664E-03_wp, &
      & 5.48606965350553E-18_wp, 2.60581698338229E-03_wp,-3.95709921802641E-18_wp, &
      & 2.90249899222578E-05_wp, 5.10782472349948E-18_wp, 1.48543226707051E-07_wp, &
      & 7.12070034072669E-17_wp,-6.62571427782813E-10_wp, 1.26845085436653E-18_wp, &
      &-1.42651174861375E-06_wp, 3.83847936960920E-18_wp,-2.94779095388117E-09_wp, &
      &-1.63321743037760E-04_wp, 1.87223706676972E-17_wp,-2.67894316750773E-05_wp, &
      & 6.96815674075468E-17_wp, 2.99019513933563E-09_wp,-1.55049130225332E-19_wp, &
      & 2.54472983634505E-06_wp,-3.95709921802641E-18_wp, 1.33034263364267E-08_wp],&
      & shape(density))

   call get_structure(mol, "MB16-43", "S2")
   call test_e_gen(error, mol, ref, thr_in=thr1*10, guess_qat=qat, guess_qsh=qsh, &
      & guess_density=density, guess_dp=dp, guess_qp=qp)

end subroutine test_e_s2


subroutine test_e_cecl3(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol
   ! Reference energy calculated with g-xTB development code (stricter convergence)
   real(wp), parameter :: ref = -1855.84992878162_wp 

   ! Fully converged values from the g-xTB development code
   real(wp), parameter :: qat(4, 2) = reshape([&
      & 1.76149086526932E+00_wp,-5.87174062893260E-01_wp,-5.87964325848794E-01_wp, &
      &-5.86352476679365E-01_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp], shape(qat))

   real(wp), parameter :: qsh(13, 2) = reshape([&
      & 8.56111616941588E-01_wp, 1.08157069403099E-01_wp, 8.39880071173497E-01_wp, &
      &-4.26578922488661E-02_wp,-9.50257492879911E-02_wp,-5.37611554235073E-01_wp, &
      & 4.54632406298039E-02_wp,-9.62011513600360E-02_wp,-5.37310669549760E-01_wp, &
      & 4.55474950610019E-02_wp,-9.54523802276321E-02_wp,-5.36431942194229E-01_wp, &
      & 4.55318457424962E-02_wp,-4.80198892142368E-04_wp, 2.02908435901956E-03_wp, &
      & 8.82622777976660E-03_wp, 1.04031437192929E+00_wp,-4.87663452155784E-04_wp, &
      &-1.35348821728360E-02_wp, 8.79790216689616E-04_wp,-5.51731744499628E-04_wp, &
      &-2.08559823405556E-02_wp, 4.76435941974605E-04_wp,-5.19725766091605E-04_wp, &
      &-1.68053375010775E-02_wp, 7.09611642610962E-04_wp], shape(qsh))

   real(wp), parameter :: dp(3, 4, 2) = reshape([&
      & 9.93296371203635E-02_wp, 1.68223053140515E-01_wp,-2.59177340439747E-01_wp, &
      &-1.87670116724596E-01_wp,-2.73518574388596E-01_wp,-1.49288398132009E-01_wp, &
      &-2.04014954986284E-01_wp, 2.75361082437328E-01_wp, 1.61555539492173E-01_wp, &
      & 3.19829776865433E-01_wp,-9.79364767646937E-02_wp, 1.56313884843267E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp],&
      & shape(dp))

   real(wp), parameter :: qp(6, 4, 2) = reshape([&
      &-7.79025899224358E-01_wp, 2.93628532591290E-01_wp,-1.87211885973218E-01_wp, &
      &-9.66216061189907E-01_wp,-1.26193462425494E+00_wp, 9.66237785197578E-01_wp, &
      & 1.66472202112189E-03_wp,-2.19517452595636E-02_wp,-8.66870189910252E-03_wp, &
      &-2.20927749725244E-02_wp,-3.13751546430074E-02_wp, 7.00397987797885E-03_wp, &
      &-2.71675086777812E-03_wp, 3.33171880750437E-02_wp,-2.09838581426460E-02_wp, &
      & 5.34229826270407E-03_wp,-3.31619231860499E-02_wp, 2.37006090104224E-02_wp, &
      &-4.05514795491584E-02_wp, 1.50748714007049E-02_wp, 2.04623838256381E-02_wp, &
      &-2.82623824959281E-02_wp,-3.97145101398323E-03_wp, 2.00890957235274E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp],&
      & shape(qp))
      
   real(wp), parameter :: density(43, 43, 2) = reshape([&
      & 6.91508514884115E-03_wp,-4.64128428381965E-04_wp, 7.39851890191156E-04_wp, &
      &-3.08180901288227E-04_wp,-2.70834565861992E-03_wp, 6.73888500345714E-03_wp, &
      &-5.28983463453815E-03_wp, 4.82025506336992E-03_wp, 9.20947878209081E-04_wp, &
      &-2.62666746814417E-03_wp,-1.03344549755776E-02_wp,-2.13949303343975E-03_wp, &
      &-2.86177079831766E-03_wp,-2.19360855296546E-03_wp, 3.65331330489966E-03_wp, &
      & 1.10745293250135E-02_wp,-8.43269606675483E-04_wp, 3.10523011533081E-02_wp, &
      & 2.30387204910565E-02_wp, 2.04828269566997E-02_wp, 1.17529497110506E-03_wp, &
      & 1.12189815116146E-03_wp,-3.37676390958414E-04_wp, 7.56850361168490E-04_wp, &
      &-4.90010996783169E-04_wp,-1.05403814114163E-03_wp,-3.44887534829157E-02_wp, &
      &-1.34351071293633E-02_wp, 2.17894300125732E-02_wp,-1.32317738878282E-03_wp, &
      & 9.23295714966318E-04_wp,-5.90066311108052E-04_wp,-6.04014519638357E-04_wp, &
      &-5.29941099561523E-04_wp,-9.36835419320956E-04_wp, 9.73171784325228E-03_wp, &
      &-1.32361102143554E-02_wp,-4.00527620687556E-02_wp,-7.26493843775791E-04_wp, &
      &-2.84864087459862E-04_wp,-5.84338260819475E-04_wp, 1.06457926932973E-03_wp, &
      & 1.20863288757621E-03_wp,-4.64128428381965E-04_wp, 2.35817284818712E-03_wp, &
      & 3.92587857259299E-04_wp,-1.97055165549110E-04_wp,-2.70852366440473E-03_wp, &
      & 7.57822519784544E-04_wp, 3.58502225940096E-05_wp,-1.71837573182039E-03_wp, &
      &-1.11630435040151E-04_wp, 5.43646927839607E-04_wp, 1.70307470883770E-04_wp, &
      & 3.32526426234747E-03_wp,-3.25191253821528E-03_wp, 6.16979032414891E-03_wp, &
      &-2.47482306012547E-03_wp,-2.24705243703616E-04_wp, 6.95417057774431E-03_wp, &
      &-1.08633255493783E-02_wp,-1.76598581111579E-02_wp,-1.93609631896134E-02_wp, &
      &-8.72287154795862E-04_wp,-5.65631496294248E-04_wp, 1.93375486670924E-04_wp, &
      &-6.13369523145670E-04_wp,-7.09403501203854E-05_wp,-7.79253334557071E-03_wp, &
      &-1.13057986258420E-02_wp,-1.77696502981528E-02_wp, 2.23419867269140E-02_wp, &
      &-8.58013966216079E-04_wp, 6.57139908879595E-04_wp,-1.13105836830082E-04_wp, &
      &-7.25066400813867E-04_wp, 1.72200596151693E-04_wp, 2.47719300655217E-03_wp, &
      & 1.50974473410191E-02_wp, 4.19682133691405E-03_wp, 1.19489386836249E-02_wp, &
      &-4.06772542652344E-04_wp,-1.49188497758960E-04_wp,-1.57316635327317E-05_wp, &
      &-4.42293296849883E-04_wp,-5.24583626459510E-04_wp, 7.39851890191156E-04_wp, &
      & 3.92587857259299E-04_wp, 1.91536792757105E-03_wp, 2.47696727692812E-04_wp, &
      &-1.35579823647709E-03_wp,-1.33623661961980E-03_wp, 2.74709544437153E-04_wp, &
      &-5.62967318479205E-04_wp, 4.27220362358588E-04_wp, 1.26153513658508E-03_wp, &
      & 2.11075105964037E-03_wp,-7.92507173129347E-03_wp,-3.74937418369353E-03_wp, &
      &-9.22900502394054E-03_wp,-1.87857625789110E-03_wp, 2.91417129341749E-03_wp, &
      & 4.48629524106340E-03_wp,-1.94279743453715E-02_wp, 1.08495811731455E-02_wp, &
      &-1.31613564875938E-02_wp,-6.31970934436649E-04_wp,-1.51927600393151E-04_wp, &
      & 4.76566544317181E-04_wp,-6.47001035862719E-05_wp, 2.63809182457078E-04_wp, &
      &-3.70243694837620E-03_wp,-1.60655815988062E-02_wp, 1.10901045966389E-02_wp, &
      & 1.06050890565114E-02_wp,-6.86780652659270E-04_wp,-1.07487287567163E-05_wp, &
      &-5.80145263261248E-04_wp,-3.80501932542685E-05_wp,-1.88723964487023E-04_wp, &
      &-3.49727589459609E-03_wp, 4.62054256448873E-03_wp, 1.19028870012948E-02_wp, &
      &-1.89335924402475E-02_wp,-4.11420701809817E-04_wp,-8.15775787828162E-05_wp, &
      &-5.09350857470944E-04_wp, 4.03614321206780E-05_wp, 4.69680812585431E-04_wp, &
      &-3.08180901288227E-04_wp,-1.97055165549110E-04_wp, 2.47696727692812E-04_wp, &
      & 2.47716044945277E-03_wp,-1.15072008734581E-03_wp,-1.49843085926392E-03_wp, &
      &-1.80460860339231E-06_wp, 1.76900836181693E-03_wp, 2.09637914900878E-03_wp, &
      &-2.83985103032821E-03_wp, 1.88829189961650E-03_wp, 4.80175313084219E-03_wp, &
      &-2.23538245636728E-03_wp, 1.98252167593746E-03_wp, 3.75685726704417E-03_wp, &
      &-1.16225091398365E-03_wp, 4.60018757159525E-03_wp,-1.97342688404856E-02_wp, &
      &-1.21032106402909E-02_wp, 5.11250827437741E-03_wp,-2.83717679476789E-04_wp, &
      &-5.88102839630281E-04_wp, 1.40865590757356E-04_wp,-9.26778827042833E-05_wp, &
      & 7.36601869053522E-04_wp, 5.33024933890096E-03_wp, 2.09425872447759E-02_wp, &
      & 1.10099056243336E-02_wp, 3.15743420378640E-03_wp, 3.35355290917966E-04_wp, &
      &-7.16937718373765E-04_wp, 5.04918708188073E-05_wp, 2.25742271534700E-04_wp, &
      & 6.87322165500211E-04_wp,-8.83698060891256E-03_wp, 1.35558312928542E-02_wp, &
      &-2.03523547681287E-02_wp,-2.27896294371807E-02_wp,-6.80990386544041E-04_wp, &
      &-4.08274695309232E-04_wp,-1.56609377614967E-04_wp, 1.00472117204186E-03_wp, &
      & 6.42112571717576E-04_wp,-2.70834565861992E-03_wp,-2.70852366440473E-03_wp, &
      &-1.35579823647709E-03_wp,-1.15072008734581E-03_wp, 4.40698408457315E-02_wp, &
      & 5.53684835982152E-03_wp, 2.57660969737744E-03_wp, 7.74565445322245E-03_wp, &
      &-4.29172335721274E-04_wp, 8.87308179509967E-04_wp, 7.59819987168787E-03_wp, &
      &-2.90143366303692E-03_wp, 1.65065497842907E-02_wp,-8.56341761190990E-03_wp, &
      &-9.61754361712255E-04_wp, 2.03346391796122E-03_wp, 2.24529473315599E-02_wp, &
      & 6.64065788186632E-02_wp, 8.88286037294378E-02_wp, 6.09775054127693E-04_wp, &
      & 1.00856343231612E-03_wp, 2.43202024557563E-03_wp, 2.80107501851078E-04_wp, &
      & 1.28797142871154E-03_wp,-1.29393640101282E-03_wp,-2.40535809207919E-02_wp, &
      & 7.93344271687818E-02_wp, 7.02584853006239E-02_wp,-1.26549252593272E-02_wp, &
      & 1.71977098435945E-03_wp,-2.33627951322494E-03_wp, 1.31736619868509E-04_wp, &
      & 1.23881950039292E-03_wp, 1.31126852060688E-03_wp,-1.35887084557653E-02_wp, &
      & 7.65772988251802E-02_wp, 3.30174815943171E-02_wp, 6.43422108828071E-02_wp, &
      &-1.00280882990665E-03_wp,-7.03928758350215E-04_wp, 3.69302881743411E-04_wp, &
      &-1.38995877590136E-03_wp,-2.37514750915699E-03_wp, 6.73888500345714E-03_wp, &
      & 7.57822519784544E-04_wp,-1.33623661961980E-03_wp,-1.49843085926392E-03_wp, &
      & 5.53684835982152E-03_wp, 3.15062482341216E-02_wp, 5.47450138896558E-03_wp, &
      &-5.00796520472959E-03_wp,-7.28454099059860E-03_wp, 2.79752158035175E-03_wp, &
      &-1.71486204161579E-02_wp, 7.94977673651953E-03_wp, 1.65032263005495E-03_wp, &
      & 3.93584481005497E-03_wp, 8.76221073515787E-04_wp, 9.74767674057578E-03_wp, &
      & 2.07568900075376E-02_wp, 6.76261685007141E-02_wp,-2.17931392833033E-02_wp, &
      & 7.90366353904093E-02_wp, 2.79402996272633E-03_wp, 6.30314688624698E-04_wp, &
      &-1.78080649698502E-03_wp, 6.08345016822846E-04_wp,-4.94440977938366E-04_wp, &
      & 1.84723470423765E-02_wp,-6.44273965614169E-02_wp, 3.79307105674980E-02_wp, &
      & 6.86738378968837E-02_wp,-2.66785075883609E-03_wp, 3.39078332490191E-04_wp, &
      &-2.13068196013103E-03_wp,-3.19055500196838E-04_wp,-6.05876862607398E-04_wp, &
      &-5.22575415051881E-03_wp, 3.76217881067743E-02_wp,-1.71192654078014E-02_wp, &
      & 4.37862245341372E-02_wp,-7.45924910929455E-04_wp,-5.93745257585652E-04_wp, &
      & 5.99043014817094E-04_wp, 9.36639684742331E-05_wp,-1.11598042499230E-03_wp, &
      &-5.28983463453815E-03_wp, 3.58502225940096E-05_wp, 2.74709544437153E-04_wp, &
      &-1.80460860339231E-06_wp, 2.57660969737744E-03_wp, 5.47450138896558E-03_wp, &
      & 2.93992940596931E-02_wp, 3.63929394233870E-03_wp,-7.16353689179876E-04_wp, &
      & 3.39231449448403E-03_wp,-1.74717335334318E-04_wp, 3.77016184389070E-03_wp, &
      &-1.41840963660445E-05_wp, 5.13900921974701E-03_wp, 6.46729385593851E-04_wp, &
      &-1.71066903016120E-02_wp,-7.55036597130801E-03_wp, 1.56196290712313E-02_wp, &
      &-9.66954896066762E-02_wp, 1.11349868418996E-02_wp, 1.22601538765852E-04_wp, &
      &-2.11297496649442E-03_wp,-1.37121718553120E-03_wp,-1.46955469231190E-03_wp, &
      &-3.20354448730695E-05_wp,-9.72011422056808E-03_wp, 3.77941388679876E-04_wp, &
      & 8.63683575762467E-02_wp,-5.52544715963979E-03_wp, 2.37974792536134E-04_wp, &
      &-1.95351964708941E-03_wp,-1.29213692990610E-03_wp, 1.69539315415403E-03_wp, &
      &-2.70533130730944E-04_wp,-9.91706774102228E-03_wp,-5.51719125082554E-03_wp, &
      & 8.65389931089010E-02_wp, 3.02004039852569E-03_wp, 3.82254772656266E-04_wp, &
      & 9.48497348959330E-04_wp,-1.24471946550559E-03_wp,-2.38179325663893E-03_wp, &
      & 4.90483620803742E-05_wp, 4.82025506336992E-03_wp,-1.71837573182039E-03_wp, &
      &-5.62967318479205E-04_wp, 1.76900836181693E-03_wp, 7.74565445322245E-03_wp, &
      &-5.00796520472959E-03_wp, 3.63929394233870E-03_wp, 3.44559711845734E-02_wp, &
      & 4.73877632173727E-03_wp,-1.08018769154115E-02_wp,-4.09594036775609E-03_wp, &
      &-1.24377496915671E-02_wp, 8.23484487858084E-04_wp,-3.00725027685771E-03_wp, &
      & 1.05005114623486E-02_wp, 5.68409823837784E-03_wp, 1.41430409071168E-02_wp, &
      & 7.87648700134890E-02_wp,-1.31132113589655E-02_wp, 5.05354123682068E-03_wp, &
      & 1.56860476651673E-03_wp, 4.95051949598724E-04_wp,-1.40540813271128E-03_wp, &
      &-4.36472660545697E-06_wp,-1.55186598747247E-03_wp,-1.19396609932784E-02_wp, &
      & 7.87414454388724E-02_wp,-2.99966905662299E-02_wp,-2.09673880851737E-03_wp, &
      & 1.14415784568630E-03_wp,-1.63860849716878E-04_wp, 1.53063976201816E-03_wp, &
      &-6.84468500357535E-04_wp, 1.76369974612941E-03_wp, 2.14468930226159E-02_wp, &
      & 3.28265297953579E-02_wp, 4.80550508488050E-02_wp,-1.04816251636682E-01_wp, &
      &-1.57983304804792E-03_wp,-8.84168815739565E-05_wp,-2.38428556173981E-03_wp, &
      & 4.96039250557894E-04_wp, 2.58030054106487E-03_wp, 9.20947878209081E-04_wp, &
      &-1.11630435040151E-04_wp, 4.27220362358588E-04_wp, 2.09637914900878E-03_wp, &
      &-4.29172335721274E-04_wp,-7.28454099059860E-03_wp,-7.16353689179876E-04_wp, &
      & 4.73877632173727E-03_wp, 4.33658763762548E-02_wp,-3.56787563595711E-03_wp, &
      &-6.98093184673465E-04_wp, 7.94226083545157E-03_wp,-6.12346751361982E-03_wp, &
      &-1.95149951492134E-03_wp, 7.85741728972451E-03_wp,-1.75116019159196E-03_wp, &
      &-8.86026341157948E-03_wp, 3.47893285568553E-02_wp,-3.47127043956547E-02_wp, &
      &-8.92145888526571E-02_wp,-1.31598382335613E-03_wp,-3.50664103863002E-04_wp, &
      &-9.20913975699812E-05_wp,-1.31550407059796E-03_wp,-1.78591244473431E-03_wp, &
      &-8.62500575397342E-03_wp,-3.79154983079839E-02_wp, 3.52606725581625E-02_wp, &
      &-9.70310288993085E-02_wp, 1.09295280591295E-03_wp,-3.65715279297370E-04_wp, &
      &-3.71310843920043E-04_wp, 1.86764638864081E-03_wp,-2.30237562464351E-03_wp, &
      & 2.20223974235199E-02_wp, 7.33725876257370E-02_wp,-7.21093547466876E-02_wp, &
      &-4.68711684094018E-02_wp,-2.11281902812233E-03_wp,-1.44674360192421E-03_wp, &
      & 1.06346755131592E-04_wp, 2.30301478670150E-03_wp, 5.27801624728434E-04_wp, &
      &-2.62666746814417E-03_wp, 5.43646927839607E-04_wp, 1.26153513658508E-03_wp, &
      &-2.83985103032821E-03_wp, 8.87308179509967E-04_wp, 2.79752158035175E-03_wp, &
      & 3.39231449448403E-03_wp,-1.08018769154115E-02_wp,-3.56787563595711E-03_wp, &
      & 2.93612618474126E-02_wp, 2.31391722613461E-02_wp,-6.60986227282958E-02_wp, &
      & 1.24986207340805E-03_wp,-6.24507721166043E-02_wp,-1.03214399155631E-02_wp, &
      & 1.32639693946498E-02_wp, 4.71236001984045E-04_wp,-4.68879492302323E-02_wp, &
      &-2.47260223304392E-03_wp, 4.41570070864646E-02_wp, 1.02315953426925E-03_wp, &
      &-1.82803803735718E-03_wp,-8.23219426039424E-04_wp, 9.92990884166748E-05_wp, &
      & 2.01374669314284E-03_wp, 3.07721358115098E-04_wp,-5.84884265387732E-02_wp, &
      &-1.64182438222803E-02_wp,-5.32622591594028E-02_wp, 5.34721367813905E-04_wp, &
      & 1.82873388025052E-03_wp, 5.62879688324057E-04_wp, 1.47738399320061E-05_wp, &
      &-2.63393191871144E-03_wp,-1.29176122752557E-04_wp, 3.91140896726203E-02_wp, &
      & 5.34659737495685E-02_wp, 5.84897036304593E-02_wp,-1.07874392638474E-03_wp, &
      &-7.53049287002383E-04_wp, 3.71126983773621E-04_wp,-2.11151989500910E-03_wp, &
      &-2.94549789451623E-03_wp,-1.03344549755776E-02_wp, 1.70307470883770E-04_wp, &
      & 2.11075105964037E-03_wp, 1.88829189961650E-03_wp, 7.59819987168787E-03_wp, &
      &-1.71486204161579E-02_wp,-1.74717335334318E-04_wp,-4.09594036775609E-03_wp, &
      &-6.98093184673465E-04_wp, 2.31391722613461E-02_wp, 5.91406735819692E-02_wp, &
      &-1.22440708946439E-01_wp, 9.13015853525735E-04_wp,-1.30859206658707E-01_wp, &
      &-2.79627346033901E-02_wp, 1.48749197875323E-02_wp, 3.28851241135510E-04_wp, &
      &-9.20016419634780E-02_wp,-9.28743574717102E-03_wp,-3.36525729953606E-02_wp, &
      &-1.36607297660148E-03_wp,-4.04468477623580E-03_wp,-5.61533155460872E-04_wp, &
      &-1.79865060685337E-03_wp, 1.74946817518995E-03_wp, 1.11503351299903E-04_wp, &
      & 7.63323219735722E-02_wp, 4.93070855422708E-03_wp,-3.14965189527595E-02_wp, &
      & 2.32956349384313E-03_wp,-9.71795793234602E-04_wp, 2.00788514025407E-03_wp, &
      &-3.47553523950925E-04_wp, 1.90670496921082E-03_wp, 8.63729250654217E-04_wp, &
      & 2.78574653236564E-02_wp, 8.43210392843945E-03_wp, 4.10238241669585E-02_wp, &
      &-7.73909621631388E-04_wp,-1.64696146897236E-03_wp, 1.82835861706770E-03_wp, &
      & 3.86716929434750E-04_wp,-3.09145824765568E-03_wp,-2.13949303343975E-03_wp, &
      & 3.32526426234747E-03_wp,-7.92507173129347E-03_wp, 4.80175313084219E-03_wp, &
      &-2.90143366303692E-03_wp, 7.94977673651953E-03_wp, 3.77016184389070E-03_wp, &
      &-1.24377496915671E-02_wp, 7.94226083545157E-03_wp,-6.60986227282958E-02_wp, &
      &-1.22440708946439E-01_wp, 4.33566898170077E-01_wp, 1.95813898813488E-03_wp, &
      & 4.51517201212699E-01_wp, 8.78248098592704E-02_wp,-9.51337187011603E-02_wp, &
      &-1.25511734610843E-03_wp,-7.63119679638036E-03_wp,-1.11538571281217E-02_wp, &
      &-2.66649662958480E-04_wp,-4.41948382247958E-03_wp, 6.85887018088822E-03_wp, &
      & 7.11138177727583E-03_wp, 3.37362565469666E-03_wp, 1.44852243036690E-03_wp, &
      & 5.29637457862276E-04_wp, 8.80439773935814E-03_wp, 3.91331273267315E-02_wp, &
      & 9.30637689348436E-03_wp, 4.90262202098411E-05_wp,-2.82389446940278E-03_wp, &
      &-2.61620860312370E-03_wp, 3.78836342506481E-03_wp,-7.48580181756997E-04_wp, &
      &-1.70277108432316E-04_wp, 1.13608539646023E-02_wp,-6.50953860124206E-02_wp, &
      & 2.26686761566592E-02_wp, 6.83488977352320E-04_wp, 3.49268135061221E-03_wp, &
      &-3.31406784596862E-03_wp,-2.97731898645141E-03_wp, 4.21572546220010E-03_wp, &
      &-2.86177079831766E-03_wp,-3.25191253821528E-03_wp,-3.74937418369353E-03_wp, &
      &-2.23538245636728E-03_wp, 1.65065497842907E-02_wp, 1.65032263005495E-03_wp, &
      &-1.41840963660445E-05_wp, 8.23484487858084E-04_wp,-6.12346751361982E-03_wp, &
      & 1.24986207340805E-03_wp, 9.13015853525735E-04_wp, 1.95813898813488E-03_wp, &
      & 1.54119321480556E-02_wp, 2.26611010331695E-03_wp, 2.67954178420509E-04_wp, &
      &-4.53822937519479E-03_wp, 7.02110491915086E-04_wp, 5.25286901117441E-02_wp, &
      & 3.34564803756770E-02_wp, 3.58398520224896E-02_wp, 1.95120372343408E-03_wp, &
      & 1.52226885827680E-03_wp,-7.20172591204559E-04_wp, 1.03692381433069E-03_wp, &
      &-7.60687881040046E-04_wp,-3.00205239528501E-04_wp, 4.94122247974750E-02_wp, &
      & 6.09409119288440E-03_wp,-4.30859858617094E-02_wp, 2.32658889807607E-03_wp, &
      &-7.53514102241452E-04_wp, 1.24050141092100E-03_wp, 8.55035288845360E-04_wp, &
      & 1.67231401206569E-04_wp,-2.42764062226375E-04_wp,-2.43253126656368E-02_wp, &
      & 5.40302200976721E-03_wp, 6.11538184986990E-02_wp, 1.67995511737816E-03_wp, &
      & 5.55405529262713E-04_wp, 1.20746182337580E-03_wp,-9.68246628988539E-04_wp, &
      &-1.52223779069168E-03_wp,-2.19360855296546E-03_wp, 6.16979032414891E-03_wp, &
      &-9.22900502394054E-03_wp, 1.98252167593746E-03_wp,-8.56341761190990E-03_wp, &
      & 3.93584481005497E-03_wp, 5.13900921974701E-03_wp,-3.00725027685771E-03_wp, &
      &-1.95149951492134E-03_wp,-6.24507721166043E-02_wp,-1.30859206658707E-01_wp, &
      & 4.51517201212699E-01_wp, 2.26611010331695E-03_wp, 4.94207676645429E-01_wp, &
      & 9.42895852759774E-02_wp,-9.92737254972440E-02_wp,-1.38059090374086E-03_wp, &
      & 2.33082309751693E-03_wp,-3.34569171741958E-02_wp, 7.72602495057111E-03_wp, &
      &-4.36466559079192E-03_wp, 6.92263473084404E-03_wp, 7.01850986280783E-03_wp, &
      & 3.33674783782453E-03_wp, 1.34879934634476E-03_wp, 1.35486725625938E-03_wp, &
      & 8.94652619621598E-03_wp,-6.96972599091764E-02_wp, 1.83033351486987E-03_wp, &
      & 2.45042748729643E-04_wp, 1.87685383965128E-04_wp,-6.56150423365870E-04_wp, &
      & 1.97823885362361E-03_wp,-1.35160374576819E-03_wp,-1.19390525889912E-03_wp, &
      & 2.31843308381987E-02_wp, 2.60927043452444E-02_wp, 2.77810746244067E-02_wp, &
      & 2.20034973314928E-04_wp, 4.37709869636345E-03_wp,-5.12575192802934E-03_wp, &
      &-6.39253233006904E-03_wp, 3.81725763353449E-03_wp, 3.65331330489966E-03_wp, &
      &-2.47482306012547E-03_wp,-1.87857625789110E-03_wp, 3.75685726704417E-03_wp, &
      &-9.61754361712255E-04_wp, 8.76221073515787E-04_wp, 6.46729385593851E-04_wp, &
      & 1.05005114623486E-02_wp, 7.85741728972451E-03_wp,-1.03214399155631E-02_wp, &
      &-2.79627346033901E-02_wp, 8.78248098592704E-02_wp, 2.67954178420509E-04_wp, &
      & 9.42895852759774E-02_wp, 3.69230914926252E-02_wp,-1.38056289411748E-02_wp, &
      &-1.76980300375920E-04_wp,-3.51744465800494E-03_wp,-2.38840898506098E-03_wp, &
      & 7.69775758691067E-02_wp, 1.36644796891834E-03_wp, 1.37586040259405E-03_wp, &
      & 7.48766171935598E-04_wp, 2.13490442713015E-03_wp, 2.25243325147397E-03_wp, &
      &-5.92615522027911E-04_wp, 1.22050037972373E-02_wp,-5.67480870018552E-04_wp, &
      &-5.95872033137470E-02_wp, 1.86239368261565E-03_wp,-4.55773812572026E-04_wp, &
      & 4.50816048701901E-04_wp, 1.66782105440837E-03_wp,-1.27145991723743E-03_wp, &
      & 1.03577075031939E-04_wp, 4.75823268913222E-02_wp,-1.10489742679580E-02_wp, &
      &-6.98074449911554E-02_wp,-2.04606981318891E-03_wp,-4.80369404248456E-06_wp, &
      &-2.27506888247999E-03_wp, 4.87712090937551E-04_wp, 2.67023870340793E-03_wp, &
      & 1.10745293250135E-02_wp,-2.24705243703616E-04_wp, 2.91417129341749E-03_wp, &
      &-1.16225091398365E-03_wp, 2.03346391796122E-03_wp, 9.74767674057578E-03_wp, &
      &-1.71066903016120E-02_wp, 5.68409823837784E-03_wp,-1.75116019159196E-03_wp, &
      & 1.32639693946498E-02_wp, 1.48749197875323E-02_wp,-9.51337187011603E-02_wp, &
      &-4.53822937519479E-03_wp,-9.92737254972440E-02_wp,-1.38056289411748E-02_wp, &
      & 4.32218247474175E-02_wp,-4.67504351396502E-04_wp, 2.48296300979312E-02_wp, &
      & 7.10728100204719E-02_wp, 3.34202618958169E-02_wp, 2.32679657233186E-03_wp, &
      & 5.56849442015036E-04_wp,-1.25012981617534E-03_wp, 9.41209570149681E-04_wp, &
      &-3.02768365096295E-04_wp,-3.55593063923384E-04_wp,-4.93732028737240E-02_wp, &
      &-4.86975181966102E-02_wp, 4.52877356366966E-02_wp,-2.27875367667909E-03_wp, &
      & 2.50120395418680E-03_wp,-5.97905730924909E-05_wp,-2.37569126804638E-03_wp, &
      &-2.62248026531996E-04_wp,-2.86511742670158E-04_wp, 5.67693761700763E-02_wp, &
      &-2.51531663908916E-02_wp,-3.76022596630043E-02_wp,-2.46662270361714E-03_wp, &
      &-2.04534482380589E-03_wp, 2.56487472738511E-04_wp, 2.33410259062064E-03_wp, &
      &-2.89952769295361E-04_wp,-8.43269606675483E-04_wp, 6.95417057774431E-03_wp, &
      & 4.48629524106340E-03_wp, 4.60018757159525E-03_wp, 2.24529473315599E-02_wp, &
      & 2.07568900075376E-02_wp,-7.55036597130801E-03_wp, 1.41430409071168E-02_wp, &
      &-8.86026341157948E-03_wp, 4.71236001984045E-04_wp, 3.28851241135510E-04_wp, &
      &-1.25511734610843E-03_wp, 7.02110491915086E-04_wp,-1.38059090374086E-03_wp, &
      &-1.76980300375920E-04_wp,-4.67504351396502E-04_wp, 9.79337804800354E-01_wp, &
      &-6.74219468612087E-02_wp,-4.24289558786069E-02_wp,-4.55035171770853E-02_wp, &
      &-3.44483283793832E-03_wp,-3.19796187466709E-03_wp, 1.08310917393414E-03_wp, &
      &-2.17887581858298E-03_wp, 1.36316065785928E-03_wp,-1.74330998715765E-03_wp, &
      &-4.58677875448408E-03_wp,-4.64473277852141E-03_wp,-2.29331504156114E-02_wp, &
      & 1.08820023278782E-03_wp, 8.19971279780699E-04_wp, 8.60563407195275E-04_wp, &
      & 3.74866221979298E-04_wp,-2.27285687292247E-03_wp,-4.66589175119681E-04_wp, &
      &-2.18088387623203E-02_wp,-4.84923080911854E-03_wp, 2.45274067402458E-03_wp, &
      & 2.10826591425014E-03_wp, 5.15305049073426E-04_wp, 8.24138929388585E-04_wp, &
      & 7.13070191665346E-04_wp, 1.14992353955420E-03_wp, 3.10523011533081E-02_wp, &
      &-1.08633255493783E-02_wp,-1.94279743453715E-02_wp,-1.97342688404856E-02_wp, &
      & 6.64065788186632E-02_wp, 6.76261685007141E-02_wp, 1.56196290712313E-02_wp, &
      & 7.87648700134890E-02_wp, 3.47893285568553E-02_wp,-4.68879492302323E-02_wp, &
      &-9.20016419634780E-02_wp,-7.63119679638036E-03_wp, 5.25286901117441E-02_wp, &
      & 2.33082309751693E-03_wp,-3.51744465800494E-03_wp, 2.48296300979312E-02_wp, &
      &-6.74219468612087E-02_wp, 8.89616589730809E-01_wp,-1.30912516202371E-02_wp, &
      &-1.38437828918995E-02_wp, 1.53689321441061E-02_wp, 1.30549601532400E-02_wp, &
      &-1.30573212534305E-02_wp,-1.98870685257982E-05_wp,-2.41087721216451E-02_wp, &
      & 5.83917498577022E-03_wp, 4.54851475235891E-03_wp,-1.29884316231089E-02_wp, &
      &-1.20861867873460E-02_wp,-6.93167404011562E-04_wp, 2.31083585867816E-03_wp, &
      & 3.77275766987838E-04_wp,-9.37624650851280E-04_wp,-9.82782656149650E-04_wp, &
      &-1.79306697564405E-02_wp,-4.45044026214113E-02_wp,-4.51857000032706E-03_wp, &
      & 3.23050764298899E-02_wp, 5.49801573247503E-03_wp, 2.23738676373525E-03_wp, &
      & 2.14802250375719E-03_wp, 2.15130914664620E-04_wp,-2.63692688715765E-04_wp, &
      & 2.30387204910565E-02_wp,-1.76598581111579E-02_wp, 1.08495811731455E-02_wp, &
      &-1.21032106402909E-02_wp, 8.88286037294378E-02_wp,-2.17931392833033E-02_wp, &
      &-9.66954896066762E-02_wp,-1.31132113589655E-02_wp,-3.47127043956547E-02_wp, &
      &-2.47260223304392E-03_wp,-9.28743574717102E-03_wp,-1.11538571281217E-02_wp, &
      & 3.34564803756770E-02_wp,-3.34569171741958E-02_wp,-2.38840898506098E-03_wp, &
      & 7.10728100204719E-02_wp,-4.24289558786069E-02_wp,-1.30912516202371E-02_wp, &
      & 9.01707082609818E-01_wp,-9.56241147634736E-03_wp, 1.42821778264060E-03_wp, &
      & 2.13457060026934E-02_wp, 1.22148139653361E-02_wp, 1.49853782223083E-02_wp, &
      &-4.43860671944334E-04_wp, 9.69085929894596E-04_wp,-1.16033523110767E-02_wp, &
      & 2.04795466204749E-02_wp,-1.36569680115221E-02_wp,-7.21219809087939E-04_wp, &
      & 7.26272817679995E-05_wp,-1.52973074976070E-03_wp, 5.08502235285356E-04_wp, &
      &-2.09194699062579E-03_wp, 1.70651438123184E-03_wp,-1.64247478153722E-02_wp, &
      & 2.26152122557670E-02_wp,-9.39648238256889E-03_wp, 5.30552877341724E-04_wp, &
      & 1.21492917179406E-04_wp,-1.14204591564980E-03_wp, 3.16310816330874E-04_wp, &
      & 1.51717062160637E-03_wp, 2.04828269566997E-02_wp,-1.93609631896134E-02_wp, &
      &-1.31613564875938E-02_wp, 5.11250827437741E-03_wp, 6.09775054127693E-04_wp, &
      & 7.90366353904093E-02_wp, 1.11349868418996E-02_wp, 5.05354123682068E-03_wp, &
      &-8.92145888526571E-02_wp, 4.41570070864646E-02_wp,-3.36525729953606E-02_wp, &
      &-2.66649662958480E-04_wp, 3.58398520224896E-02_wp, 7.72602495057111E-03_wp, &
      & 7.69775758691067E-02_wp, 3.34202618958169E-02_wp,-4.55035171770853E-02_wp, &
      &-1.38437828918995E-02_wp,-9.56241147634736E-03_wp, 9.01322110906132E-01_wp, &
      & 2.36536811545839E-02_wp, 3.21991095436961E-04_wp,-8.65332379118361E-03_wp, &
      & 1.29451244742844E-02_wp, 1.70363842050354E-02_wp,-2.30363036376635E-02_wp, &
      & 1.65447353542449E-02_wp, 2.35204047336634E-03_wp,-5.05728278786472E-02_wp, &
      & 5.43932631484292E-03_wp,-9.04894894323018E-05_wp, 2.37278451351742E-03_wp, &
      & 2.62727660060946E-03_wp,-2.77813812297788E-03_wp, 1.34190846756862E-02_wp, &
      & 9.68094225852710E-03_wp,-1.32881870951947E-02_wp,-3.36511000867629E-03_wp, &
      &-2.43229230059397E-03_wp,-1.01164160126711E-03_wp,-7.48265888202462E-04_wp, &
      & 2.35043460435739E-03_wp, 1.60863721335192E-03_wp, 1.17529497110506E-03_wp, &
      &-8.72287154795862E-04_wp,-6.31970934436649E-04_wp,-2.83717679476789E-04_wp, &
      & 1.00856343231612E-03_wp, 2.79402996272633E-03_wp, 1.22601538765852E-04_wp, &
      & 1.56860476651673E-03_wp,-1.31598382335613E-03_wp, 1.02315953426925E-03_wp, &
      &-1.36607297660148E-03_wp,-4.41948382247958E-03_wp, 1.95120372343408E-03_wp, &
      &-4.36466559079192E-03_wp, 1.36644796891834E-03_wp, 2.32679657233186E-03_wp, &
      &-3.44483283793832E-03_wp, 1.53689321441061E-02_wp, 1.42821778264060E-03_wp, &
      & 2.36536811545839E-02_wp, 9.68828816917606E-04_wp, 2.22946953098977E-04_wp, &
      &-5.13301913686200E-04_wp, 3.54712425095362E-04_wp, 1.78410932887196E-05_wp, &
      &-9.90869530174335E-04_wp,-4.65549049501994E-04_wp,-1.08192661475894E-03_wp, &
      &-5.13523422536957E-03_wp, 2.00184140256020E-04_wp, 8.28467001722436E-05_wp, &
      & 1.24196822502666E-04_wp, 6.25153681108038E-05_wp,-1.66773420133365E-04_wp, &
      & 9.05681680054237E-04_wp,-1.85522769600568E-03_wp,-1.54510379312106E-03_wp, &
      &-8.07845307739065E-04_wp, 5.56544110078614E-05_wp,-1.75760754114177E-05_wp, &
      & 6.87759158596898E-05_wp, 1.51322538300238E-04_wp, 3.34778995454745E-05_wp, &
      & 1.12189815116146E-03_wp,-5.65631496294248E-04_wp,-1.51927600393151E-04_wp, &
      &-5.88102839630281E-04_wp, 2.43202024557563E-03_wp, 6.30314688624698E-04_wp, &
      &-2.11297496649442E-03_wp, 4.95051949598724E-04_wp,-3.50664103863002E-04_wp, &
      &-1.82803803735718E-03_wp,-4.04468477623580E-03_wp, 6.85887018088822E-03_wp, &
      & 1.52226885827680E-03_wp, 6.92263473084404E-03_wp, 1.37586040259405E-03_wp, &
      & 5.56849442015036E-04_wp,-3.19796187466709E-03_wp, 1.30549601532400E-02_wp, &
      & 2.13457060026934E-02_wp, 3.21991095436961E-04_wp, 2.22946953098977E-04_wp, &
      & 8.45466612613266E-04_wp, 2.05801176749084E-04_wp, 4.33613614775327E-04_wp, &
      &-3.53468261639624E-04_wp, 6.90995475941657E-04_wp,-2.30536758883944E-03_wp, &
      &-2.21987307540992E-04_wp,-6.65357590846186E-04_wp,-5.73476407233950E-05_wp, &
      & 4.56304342717178E-05_wp,-8.45751742463379E-05_wp, 5.05323434555636E-05_wp, &
      &-1.40434659510226E-04_wp,-3.66882568279996E-04_wp,-3.35446115555736E-03_wp, &
      &-3.70549396201395E-04_wp,-2.78558396867084E-04_wp, 1.62696676108879E-04_wp, &
      & 1.32878211231700E-04_wp,-6.24456776946116E-05_wp,-4.09465025176957E-05_wp, &
      & 1.51524435468932E-04_wp,-3.37676390958414E-04_wp, 1.93375486670924E-04_wp, &
      & 4.76566544317181E-04_wp, 1.40865590757356E-04_wp, 2.80107501851078E-04_wp, &
      &-1.78080649698502E-03_wp,-1.37121718553120E-03_wp,-1.40540813271128E-03_wp, &
      &-9.20913975699812E-05_wp,-8.23219426039424E-04_wp,-5.61533155460872E-04_wp, &
      & 7.11138177727583E-03_wp,-7.20172591204559E-04_wp, 7.01850986280783E-03_wp, &
      & 7.48766171935598E-04_wp,-1.25012981617534E-03_wp, 1.08310917393414E-03_wp, &
      &-1.30573212534305E-02_wp, 1.22148139653361E-02_wp,-8.65332379118361E-03_wp, &
      &-5.13301913686200E-04_wp, 2.05801176749084E-04_wp, 5.61141640965953E-04_wp, &
      & 1.30790843839943E-04_wp, 2.05141072904138E-04_wp, 7.03407855049382E-04_wp, &
      &-5.10132597923450E-04_wp, 1.58086708677862E-03_wp, 1.55347294906608E-03_wp, &
      &-8.38605640641629E-05_wp,-8.01750894411576E-05_wp,-1.15990126177669E-04_wp, &
      & 5.02969544916235E-05_wp, 4.77941967414101E-06_wp, 6.82450253750889E-04_wp, &
      & 1.42540445979515E-03_wp, 1.25335498180113E-03_wp,-6.61458893577938E-04_wp, &
      &-7.14060994174372E-05_wp, 3.86875972970652E-05_wp,-1.39324488358490E-04_wp, &
      &-1.15338852384393E-04_wp, 8.54210495205125E-05_wp, 7.56850361168490E-04_wp, &
      &-6.13369523145670E-04_wp,-6.47001035862719E-05_wp,-9.26778827042833E-05_wp, &
      & 1.28797142871154E-03_wp, 6.08345016822846E-04_wp,-1.46955469231190E-03_wp, &
      &-4.36472660545697E-06_wp,-1.31550407059796E-03_wp, 9.92990884166748E-05_wp, &
      &-1.79865060685337E-03_wp, 3.37362565469666E-03_wp, 1.03692381433069E-03_wp, &
      & 3.33674783782453E-03_wp, 2.13490442713015E-03_wp, 9.41209570149681E-04_wp, &
      &-2.17887581858298E-03_wp,-1.98870685257982E-05_wp, 1.49853782223083E-02_wp, &
      & 1.29451244742844E-02_wp, 3.54712425095362E-04_wp, 4.33613614775327E-04_wp, &
      & 1.30790843839943E-04_wp, 4.82943053337423E-04_wp, 2.51166086052836E-04_wp, &
      &-7.79526292526259E-04_wp,-7.94734488889291E-04_wp, 2.61362627131244E-04_wp, &
      &-3.32328556211764E-03_wp, 1.04443160535004E-04_wp, 6.48879501415757E-07_wp, &
      & 3.48543604806478E-06_wp, 1.06422065011750E-04_wp,-1.50279663061911E-04_wp, &
      & 9.18878709482770E-04_wp, 2.34666691620089E-04_wp,-3.00287949888014E-04_wp, &
      &-1.87675831672112E-03_wp,-4.37473660543290E-05_wp, 1.47390311484202E-05_wp, &
      &-9.08774681400845E-05_wp, 3.16501856547077E-05_wp, 1.29899358262604E-04_wp, &
      &-4.90010996783169E-04_wp,-7.09403501203854E-05_wp, 2.63809182457078E-04_wp, &
      & 7.36601869053522E-04_wp,-1.29393640101282E-03_wp,-4.94440977938366E-04_wp, &
      &-3.20354448730695E-05_wp,-1.55186598747247E-03_wp,-1.78591244473431E-03_wp, &
      & 2.01374669314284E-03_wp, 1.74946817518995E-03_wp, 1.44852243036690E-03_wp, &
      &-7.60687881040046E-04_wp, 1.34879934634476E-03_wp, 2.25243325147397E-03_wp, &
      &-3.02768365096295E-04_wp, 1.36316065785928E-03_wp,-2.41087721216451E-02_wp, &
      &-4.43860671944334E-04_wp, 1.70363842050354E-02_wp, 1.78410932887196E-05_wp, &
      &-3.53468261639624E-04_wp, 2.05141072904138E-04_wp, 2.51166086052836E-04_wp, &
      & 1.00086837618209E-03_wp,-2.26253345482604E-03_wp, 1.21437349393744E-03_wp, &
      & 1.43726633423116E-03_wp,-2.97126073135876E-03_wp, 1.86692789766110E-04_wp, &
      &-1.00760302788538E-04_wp, 4.75140035698790E-05_wp, 1.30570945945150E-04_wp, &
      &-4.85382282645348E-05_wp, 2.19789828956273E-03_wp, 5.36001223564970E-03_wp, &
      &-1.07860478362868E-04_wp,-2.11500688326205E-03_wp,-3.10082485911504E-04_wp, &
      &-1.22377900628337E-04_wp,-1.21743448062207E-04_wp, 4.07538919273474E-05_wp, &
      & 4.44903889252788E-05_wp,-1.05403814114163E-03_wp,-7.79253334557071E-03_wp, &
      &-3.70243694837620E-03_wp, 5.33024933890096E-03_wp,-2.40535809207919E-02_wp, &
      & 1.84723470423765E-02_wp,-9.72011422056808E-03_wp,-1.19396609932784E-02_wp, &
      &-8.62500575397342E-03_wp, 3.07721358115098E-04_wp, 1.11503351299903E-04_wp, &
      & 5.29637457862276E-04_wp,-3.00205239528501E-04_wp, 1.35486725625938E-03_wp, &
      &-5.92615522027911E-04_wp,-3.55593063923384E-04_wp,-1.74330998715765E-03_wp, &
      & 5.83917498577022E-03_wp, 9.69085929894596E-04_wp,-2.30363036376635E-02_wp, &
      &-9.90869530174335E-04_wp, 6.90995475941657E-04_wp, 7.03407855049382E-04_wp, &
      &-7.79526292526259E-04_wp,-2.26253345482604E-03_wp, 9.80374960721459E-01_wp, &
      & 6.95673556711425E-02_wp, 3.39764120112855E-02_wp,-4.74420206878530E-02_wp, &
      & 3.77754294309058E-03_wp,-2.82730124351925E-03_wp, 1.53500484114958E-03_wp, &
      & 1.89164342083146E-03_wp, 1.39975897057004E-03_wp,-3.03459364770215E-03_wp, &
      & 1.59267193057405E-02_wp, 1.87802567572018E-02_wp, 5.08789290648573E-03_wp, &
      &-1.68949061274845E-03_wp,-6.53393729519532E-04_wp,-1.24410297693552E-03_wp, &
      &-1.84667925671337E-03_wp,-3.94588332180406E-04_wp,-3.44887534829157E-02_wp, &
      &-1.13057986258420E-02_wp,-1.60655815988062E-02_wp, 2.09425872447759E-02_wp, &
      & 7.93344271687818E-02_wp,-6.44273965614169E-02_wp, 3.77941388679876E-04_wp, &
      & 7.87414454388724E-02_wp,-3.79154983079839E-02_wp,-5.84884265387732E-02_wp, &
      & 7.63323219735722E-02_wp, 8.80439773935814E-03_wp, 4.94122247974750E-02_wp, &
      & 8.94652619621598E-03_wp, 1.22050037972373E-02_wp,-4.93732028737240E-02_wp, &
      &-4.58677875448408E-03_wp, 4.54851475235891E-03_wp,-1.16033523110767E-02_wp, &
      & 1.65447353542449E-02_wp,-4.65549049501994E-04_wp,-2.30536758883944E-03_wp, &
      &-5.10132597923450E-04_wp,-7.94734488889291E-04_wp, 1.21437349393744E-03_wp, &
      & 6.95673556711425E-02_wp, 8.84180392749141E-01_wp,-5.73002967894826E-03_wp, &
      & 1.37764859760830E-02_wp, 1.70763017686579E-02_wp,-1.40683480660463E-02_wp, &
      & 1.28751370446314E-02_wp, 2.60640487898893E-04_wp, 2.31072626539267E-02_wp, &
      & 9.93166183841885E-03_wp,-2.98032130752656E-03_wp,-3.93047927589246E-02_wp, &
      &-1.30198398411196E-02_wp, 1.41350626127810E-03_wp,-5.48195230523568E-04_wp, &
      & 2.00987442039049E-03_wp, 3.42911184709526E-03_wp, 1.23007486722984E-04_wp, &
      &-1.34351071293633E-02_wp,-1.77696502981528E-02_wp, 1.10901045966389E-02_wp, &
      & 1.10099056243336E-02_wp, 7.02584853006239E-02_wp, 3.79307105674980E-02_wp, &
      & 8.63683575762467E-02_wp,-2.99966905662299E-02_wp, 3.52606725581625E-02_wp, &
      &-1.64182438222803E-02_wp, 4.93070855422708E-03_wp, 3.91331273267315E-02_wp, &
      & 6.09409119288440E-03_wp,-6.96972599091764E-02_wp,-5.67480870018552E-04_wp, &
      &-4.86975181966102E-02_wp,-4.64473277852141E-03_wp,-1.29884316231089E-02_wp, &
      & 2.04795466204749E-02_wp, 2.35204047336634E-03_wp,-1.08192661475894E-03_wp, &
      &-2.21987307540992E-04_wp, 1.58086708677862E-03_wp, 2.61362627131244E-04_wp, &
      & 1.43726633423116E-03_wp, 3.39764120112855E-02_wp,-5.73002967894826E-03_wp, &
      & 9.01059729280098E-01_wp, 8.62886738828925E-03_wp,-3.56781117512672E-04_wp, &
      &-2.30967637120700E-02_wp,-1.61244951739643E-02_wp, 1.68999501118351E-02_wp, &
      & 4.38643500293381E-04_wp, 1.88978502103454E-02_wp,-2.14362405207621E-02_wp, &
      &-1.69515855882754E-02_wp,-3.46245748403255E-02_wp, 1.95194307175564E-03_wp, &
      & 3.15789256203551E-04_wp, 1.90274135264279E-04_wp, 3.95026807626506E-03_wp, &
      & 3.86843136045921E-03_wp, 2.17894300125732E-02_wp, 2.23419867269140E-02_wp, &
      & 1.06050890565114E-02_wp, 3.15743420378640E-03_wp,-1.26549252593272E-02_wp, &
      & 6.86738378968837E-02_wp,-5.52544715963979E-03_wp,-2.09673880851737E-03_wp, &
      &-9.70310288993085E-02_wp,-5.32622591594028E-02_wp,-3.14965189527595E-02_wp, &
      & 9.30637689348436E-03_wp,-4.30859858617094E-02_wp, 1.83033351486987E-03_wp, &
      &-5.95872033137470E-02_wp, 4.52877356366966E-02_wp,-2.29331504156114E-02_wp, &
      &-1.20861867873460E-02_wp,-1.36569680115221E-02_wp,-5.05728278786472E-02_wp, &
      &-5.13523422536957E-03_wp,-6.65357590846186E-04_wp, 1.55347294906608E-03_wp, &
      &-3.32328556211764E-03_wp,-2.97126073135876E-03_wp,-4.74420206878530E-02_wp, &
      & 1.37764859760830E-02_wp, 8.62886738828925E-03_wp, 9.00665213971733E-01_wp, &
      &-2.32504572549556E-02_wp,-5.61086445597062E-04_wp,-1.03203103337986E-02_wp, &
      &-1.38491527921108E-02_wp, 1.79764926201737E-02_wp, 1.33708457306665E-02_wp, &
      & 1.11144849655078E-02_wp,-9.19866444507808E-03_wp,-2.07177508980419E-03_wp, &
      &-1.79485494310802E-03_wp,-9.30706264008207E-04_wp,-4.94260832249692E-04_wp, &
      & 2.70908666784429E-03_wp, 1.56234233918496E-03_wp,-1.32317738878282E-03_wp, &
      &-8.58013966216079E-04_wp,-6.86780652659270E-04_wp, 3.35355290917966E-04_wp, &
      & 1.71977098435945E-03_wp,-2.66785075883609E-03_wp, 2.37974792536134E-04_wp, &
      & 1.14415784568630E-03_wp, 1.09295280591295E-03_wp, 5.34721367813905E-04_wp, &
      & 2.32956349384313E-03_wp, 4.90262202098411E-05_wp, 2.32658889807607E-03_wp, &
      & 2.45042748729643E-04_wp, 1.86239368261565E-03_wp,-2.27875367667909E-03_wp, &
      & 1.08820023278782E-03_wp,-6.93167404011562E-04_wp,-7.21219809087939E-04_wp, &
      & 5.43932631484292E-03_wp, 2.00184140256020E-04_wp,-5.73476407233950E-05_wp, &
      &-8.38605640641629E-05_wp, 1.04443160535004E-04_wp, 1.86692789766110E-04_wp, &
      & 3.77754294309058E-03_wp, 1.70763017686579E-02_wp,-3.56781117512672E-04_wp, &
      &-2.32504572549556E-02_wp, 9.73649840453167E-04_wp,-2.66784471781253E-04_wp, &
      & 5.37753156590185E-04_wp, 3.79172244412516E-04_wp,-2.14001521425141E-05_wp, &
      &-1.08428096201032E-03_wp,-1.81560152724172E-03_wp,-1.02967867264894E-03_wp, &
      & 1.85080213781225E-03_wp, 1.23988018064649E-04_wp, 2.72412354977432E-05_wp, &
      & 9.71754964880210E-05_wp,-1.46386688796714E-05_wp,-7.69755308471498E-05_wp, &
      & 9.23295714966318E-04_wp, 6.57139908879595E-04_wp,-1.07487287567163E-05_wp, &
      &-7.16937718373765E-04_wp,-2.33627951322494E-03_wp, 3.39078332490191E-04_wp, &
      &-1.95351964708941E-03_wp,-1.63860849716878E-04_wp,-3.65715279297370E-04_wp, &
      & 1.82873388025052E-03_wp,-9.71795793234602E-04_wp,-2.82389446940278E-03_wp, &
      &-7.53514102241452E-04_wp, 1.87685383965128E-04_wp,-4.55773812572026E-04_wp, &
      & 2.50120395418680E-03_wp, 8.19971279780699E-04_wp, 2.31083585867816E-03_wp, &
      & 7.26272817679995E-05_wp,-9.04894894323018E-05_wp, 8.28467001722436E-05_wp, &
      & 4.56304342717178E-05_wp,-8.01750894411576E-05_wp, 6.48879501415757E-07_wp, &
      &-1.00760302788538E-04_wp,-2.82730124351925E-03_wp,-1.40683480660463E-02_wp, &
      &-2.30967637120700E-02_wp,-5.61086445597062E-04_wp,-2.66784471781253E-04_wp, &
      & 8.52283031176346E-04_wp, 2.20685847996954E-04_wp,-4.57974385017563E-04_wp, &
      &-3.99086844991340E-04_wp,-1.96684355217297E-03_wp, 3.16068630072855E-03_wp, &
      & 3.82500426394768E-03_wp, 2.13899757107132E-03_wp,-1.32999227993988E-04_wp, &
      &-2.33025107255889E-05_wp,-6.08938528767931E-05_wp,-2.32989247203373E-04_wp, &
      &-1.72466828137293E-04_wp,-5.90066311108052E-04_wp,-1.13105836830082E-04_wp, &
      &-5.80145263261248E-04_wp, 5.04918708188073E-05_wp, 1.31736619868509E-04_wp, &
      &-2.13068196013103E-03_wp,-1.29213692990610E-03_wp, 1.53063976201816E-03_wp, &
      &-3.71310843920043E-04_wp, 5.62879688324057E-04_wp, 2.00788514025407E-03_wp, &
      &-2.61620860312370E-03_wp, 1.24050141092100E-03_wp,-6.56150423365870E-04_wp, &
      & 4.50816048701901E-04_wp,-5.97905730924909E-05_wp, 8.60563407195275E-04_wp, &
      & 3.77275766987838E-04_wp,-1.52973074976070E-03_wp, 2.37278451351742E-03_wp, &
      & 1.24196822502666E-04_wp,-8.45751742463379E-05_wp,-1.15990126177669E-04_wp, &
      & 3.48543604806478E-06_wp, 4.75140035698790E-05_wp, 1.53500484114958E-03_wp, &
      & 1.28751370446314E-02_wp,-1.61244951739643E-02_wp,-1.03203103337986E-02_wp, &
      & 5.37753156590185E-04_wp, 2.20685847996954E-04_wp, 6.11130765891932E-04_wp, &
      &-1.47297018815758E-04_wp, 1.23701157200696E-04_wp,-1.25333154224695E-03_wp, &
      & 2.66387587541066E-04_wp, 1.26877729245354E-04_wp, 2.10248906709898E-03_wp, &
      & 1.31163094446301E-05_wp,-2.06466519835329E-05_wp, 7.42984219798043E-05_wp, &
      &-6.15057484205189E-05_wp,-1.53448287155134E-04_wp,-6.04014519638357E-04_wp, &
      &-7.25066400813867E-04_wp,-3.80501932542685E-05_wp, 2.25742271534700E-04_wp, &
      & 1.23881950039292E-03_wp,-3.19055500196838E-04_wp, 1.69539315415403E-03_wp, &
      &-6.84468500357535E-04_wp, 1.86764638864081E-03_wp, 1.47738399320061E-05_wp, &
      &-3.47553523950925E-04_wp, 3.78836342506481E-03_wp, 8.55035288845360E-04_wp, &
      & 1.97823885362361E-03_wp, 1.66782105440837E-03_wp,-2.37569126804638E-03_wp, &
      & 3.74866221979298E-04_wp,-9.37624650851280E-04_wp, 5.08502235285356E-04_wp, &
      & 2.62727660060946E-03_wp, 6.25153681108038E-05_wp, 5.05323434555636E-05_wp, &
      & 5.02969544916235E-05_wp, 1.06422065011750E-04_wp, 1.30570945945150E-04_wp, &
      & 1.89164342083146E-03_wp, 2.60640487898893E-04_wp, 1.68999501118351E-02_wp, &
      &-1.38491527921108E-02_wp, 3.79172244412516E-04_wp,-4.57974385017563E-04_wp, &
      &-1.47297018815758E-04_wp, 5.67497382271228E-04_wp,-2.66471163188258E-04_wp, &
      &-2.05536209231146E-05_wp,-2.24672310279949E-03_wp,-9.83364562929856E-04_wp, &
      &-9.36373281820405E-04_wp, 1.09877909822793E-04_wp, 6.51748658577424E-05_wp, &
      &-7.84130470331478E-06_wp, 3.18695207242608E-05_wp, 1.13231975167471E-04_wp, &
      &-5.29941099561523E-04_wp, 1.72200596151693E-04_wp,-1.88723964487023E-04_wp, &
      & 6.87322165500211E-04_wp, 1.31126852060688E-03_wp,-6.05876862607398E-04_wp, &
      &-2.70533130730944E-04_wp, 1.76369974612941E-03_wp,-2.30237562464351E-03_wp, &
      &-2.63393191871144E-03_wp, 1.90670496921082E-03_wp,-7.48580181756997E-04_wp, &
      & 1.67231401206569E-04_wp,-1.35160374576819E-03_wp,-1.27145991723743E-03_wp, &
      &-2.62248026531996E-04_wp,-2.27285687292247E-03_wp,-9.82782656149650E-04_wp, &
      &-2.09194699062579E-03_wp,-2.77813812297788E-03_wp,-1.66773420133365E-04_wp, &
      &-1.40434659510226E-04_wp, 4.77941967414101E-06_wp,-1.50279663061911E-04_wp, &
      &-4.85382282645348E-05_wp, 1.39975897057004E-03_wp, 2.31072626539267E-02_wp, &
      & 4.38643500293381E-04_wp, 1.79764926201737E-02_wp,-2.14001521425141E-05_wp, &
      &-3.99086844991340E-04_wp, 1.23701157200696E-04_wp,-2.66471163188258E-04_wp, &
      & 9.79650159782266E-04_wp, 1.37106988591821E-03_wp,-3.44483198091035E-04_wp, &
      &-4.36588489382340E-03_wp,-8.48619259743212E-04_wp, 8.58217515795934E-06_wp, &
      &-6.62858192445359E-05_wp, 1.05720918721642E-04_wp, 2.47636898695574E-04_wp, &
      & 3.18217233255463E-05_wp,-9.36835419320956E-04_wp, 2.47719300655217E-03_wp, &
      &-3.49727589459609E-03_wp,-8.83698060891256E-03_wp,-1.35887084557653E-02_wp, &
      &-5.22575415051881E-03_wp,-9.91706774102228E-03_wp, 2.14468930226159E-02_wp, &
      & 2.20223974235199E-02_wp,-1.29176122752557E-04_wp, 8.63729250654217E-04_wp, &
      &-1.70277108432316E-04_wp,-2.42764062226375E-04_wp,-1.19390525889912E-03_wp, &
      & 1.03577075031939E-04_wp,-2.86511742670158E-04_wp,-4.66589175119681E-04_wp, &
      &-1.79306697564405E-02_wp, 1.70651438123184E-03_wp, 1.34190846756862E-02_wp, &
      & 9.05681680054237E-04_wp,-3.66882568279996E-04_wp, 6.82450253750889E-04_wp, &
      & 9.18878709482770E-04_wp, 2.19789828956273E-03_wp,-3.03459364770215E-03_wp, &
      & 9.93166183841885E-03_wp, 1.88978502103454E-02_wp, 1.33708457306665E-02_wp, &
      &-1.08428096201032E-03_wp,-1.96684355217297E-03_wp,-1.25333154224695E-03_wp, &
      &-2.05536209231146E-05_wp, 1.37106988591821E-03_wp, 9.79735154589452E-01_wp, &
      &-2.24094419081517E-02_wp, 3.40734106064249E-02_wp, 8.19089172143707E-02_wp, &
      & 2.08377776931181E-03_wp, 8.38445537400908E-04_wp, 1.55818848538031E-03_wp, &
      &-3.26712171644324E-03_wp,-3.45181082893583E-03_wp, 9.73171784325228E-03_wp, &
      & 1.50974473410191E-02_wp, 4.62054256448873E-03_wp, 1.35558312928542E-02_wp, &
      & 7.65772988251802E-02_wp, 3.76217881067743E-02_wp,-5.51719125082554E-03_wp, &
      & 3.28265297953579E-02_wp, 7.33725876257370E-02_wp, 3.91140896726203E-02_wp, &
      & 2.78574653236564E-02_wp, 1.13608539646023E-02_wp,-2.43253126656368E-02_wp, &
      & 2.31843308381987E-02_wp, 4.75823268913222E-02_wp, 5.67693761700763E-02_wp, &
      &-2.18088387623203E-02_wp,-4.45044026214113E-02_wp,-1.64247478153722E-02_wp, &
      & 9.68094225852710E-03_wp,-1.85522769600568E-03_wp,-3.35446115555736E-03_wp, &
      & 1.42540445979515E-03_wp, 2.34666691620089E-04_wp, 5.36001223564970E-03_wp, &
      & 1.59267193057405E-02_wp,-2.98032130752656E-03_wp,-2.14362405207621E-02_wp, &
      & 1.11144849655078E-02_wp,-1.81560152724172E-03_wp, 3.16068630072855E-03_wp, &
      & 2.66387587541066E-04_wp,-2.24672310279949E-03_wp,-3.44483198091035E-04_wp, &
      &-2.24094419081517E-02_wp, 9.06985946300356E-01_wp, 4.93086775727770E-03_wp, &
      & 6.81988568603333E-03_wp,-2.76733692321061E-02_wp,-1.32975683918741E-02_wp, &
      &-5.75745265195767E-03_wp,-8.97405825882184E-04_wp,-8.63422682880298E-03_wp, &
      &-1.32361102143554E-02_wp, 4.19682133691405E-03_wp, 1.19028870012948E-02_wp, &
      &-2.03523547681287E-02_wp, 3.30174815943171E-02_wp,-1.71192654078014E-02_wp, &
      & 8.65389931089010E-02_wp, 4.80550508488050E-02_wp,-7.21093547466876E-02_wp, &
      & 5.34659737495685E-02_wp, 8.43210392843945E-03_wp,-6.50953860124206E-02_wp, &
      & 5.40302200976721E-03_wp, 2.60927043452444E-02_wp,-1.10489742679580E-02_wp, &
      &-2.51531663908916E-02_wp,-4.84923080911854E-03_wp,-4.51857000032706E-03_wp, &
      & 2.26152122557670E-02_wp,-1.32881870951947E-02_wp,-1.54510379312106E-03_wp, &
      &-3.70549396201395E-04_wp, 1.25335498180113E-03_wp,-3.00287949888014E-04_wp, &
      &-1.07860478362868E-04_wp, 1.87802567572018E-02_wp,-3.93047927589246E-02_wp, &
      &-1.69515855882754E-02_wp,-9.19866444507808E-03_wp,-1.02967867264894E-03_wp, &
      & 3.82500426394768E-03_wp, 1.26877729245354E-04_wp,-9.83364562929856E-04_wp, &
      &-4.36588489382340E-03_wp, 3.40734106064249E-02_wp, 4.93086775727770E-03_wp, &
      & 9.03728560633972E-01_wp,-1.07443586780546E-02_wp,-5.83120430092661E-04_wp, &
      & 7.58962005017377E-03_wp,-1.49660914596034E-02_wp,-2.62733890723684E-02_wp, &
      &-1.06415646433602E-03_wp,-4.00527620687556E-02_wp, 1.19489386836249E-02_wp, &
      &-1.89335924402475E-02_wp,-2.27896294371807E-02_wp, 6.43422108828071E-02_wp, &
      & 4.37862245341372E-02_wp, 3.02004039852569E-03_wp,-1.04816251636682E-01_wp, &
      &-4.68711684094018E-02_wp, 5.84897036304593E-02_wp, 4.10238241669585E-02_wp, &
      & 2.26686761566592E-02_wp, 6.11538184986990E-02_wp, 2.77810746244067E-02_wp, &
      &-6.98074449911554E-02_wp,-3.76022596630043E-02_wp, 2.45274067402458E-03_wp, &
      & 3.23050764298899E-02_wp,-9.39648238256889E-03_wp,-3.36511000867629E-03_wp, &
      &-8.07845307739065E-04_wp,-2.78558396867084E-04_wp,-6.61458893577938E-04_wp, &
      &-1.87675831672112E-03_wp,-2.11500688326205E-03_wp, 5.08789290648573E-03_wp, &
      &-1.30198398411196E-02_wp,-3.46245748403255E-02_wp,-2.07177508980419E-03_wp, &
      & 1.85080213781225E-03_wp, 2.13899757107132E-03_wp, 2.10248906709898E-03_wp, &
      &-9.36373281820405E-04_wp,-8.48619259743212E-04_wp, 8.19089172143707E-02_wp, &
      & 6.81988568603333E-03_wp,-1.07443586780546E-02_wp, 8.78160059707328E-01_wp, &
      & 8.67164406929543E-03_wp, 7.25402962375287E-04_wp, 1.48330126334421E-02_wp, &
      &-1.41349407234096E-02_wp,-2.62925614731553E-02_wp,-7.26493843775791E-04_wp, &
      &-4.06772542652344E-04_wp,-4.11420701809817E-04_wp,-6.80990386544041E-04_wp, &
      &-1.00280882990665E-03_wp,-7.45924910929455E-04_wp, 3.82254772656266E-04_wp, &
      &-1.57983304804792E-03_wp,-2.11281902812233E-03_wp,-1.07874392638474E-03_wp, &
      &-7.73909621631388E-04_wp, 6.83488977352320E-04_wp, 1.67995511737816E-03_wp, &
      & 2.20034973314928E-04_wp,-2.04606981318891E-03_wp,-2.46662270361714E-03_wp, &
      & 2.10826591425014E-03_wp, 5.49801573247503E-03_wp, 5.30552877341724E-04_wp, &
      &-2.43229230059397E-03_wp, 5.56544110078614E-05_wp, 1.62696676108879E-04_wp, &
      &-7.14060994174372E-05_wp,-4.37473660543290E-05_wp,-3.10082485911504E-04_wp, &
      &-1.68949061274845E-03_wp, 1.41350626127810E-03_wp, 1.95194307175564E-03_wp, &
      &-1.79485494310802E-03_wp, 1.23988018064649E-04_wp,-1.32999227993988E-04_wp, &
      & 1.31163094446301E-05_wp, 1.09877909822793E-04_wp, 8.58217515795934E-06_wp, &
      & 2.08377776931181E-03_wp,-2.76733692321061E-02_wp,-5.83120430092661E-04_wp, &
      & 8.67164406929543E-03_wp, 9.70600913019732E-04_wp, 4.29953384215367E-04_wp, &
      & 3.37914579921590E-04_wp,-1.08472699745849E-04_wp, 7.46780711376324E-06_wp, &
      &-2.84864087459862E-04_wp,-1.49188497758960E-04_wp,-8.15775787828162E-05_wp, &
      &-4.08274695309232E-04_wp,-7.03928758350215E-04_wp,-5.93745257585652E-04_wp, &
      & 9.48497348959330E-04_wp,-8.84168815739565E-05_wp,-1.44674360192421E-03_wp, &
      &-7.53049287002383E-04_wp,-1.64696146897236E-03_wp, 3.49268135061221E-03_wp, &
      & 5.55405529262713E-04_wp, 4.37709869636345E-03_wp,-4.80369404248456E-06_wp, &
      &-2.04534482380589E-03_wp, 5.15305049073426E-04_wp, 2.23738676373525E-03_wp, &
      & 1.21492917179406E-04_wp,-1.01164160126711E-03_wp,-1.75760754114177E-05_wp, &
      & 1.32878211231700E-04_wp, 3.86875972970652E-05_wp, 1.47390311484202E-05_wp, &
      &-1.22377900628337E-04_wp,-6.53393729519532E-04_wp,-5.48195230523568E-04_wp, &
      & 3.15789256203551E-04_wp,-9.30706264008207E-04_wp, 2.72412354977432E-05_wp, &
      &-2.33025107255889E-05_wp,-2.06466519835329E-05_wp, 6.51748658577424E-05_wp, &
      &-6.62858192445359E-05_wp, 8.38445537400908E-04_wp,-1.32975683918741E-02_wp, &
      & 7.58962005017377E-03_wp, 7.25402962375287E-04_wp, 4.29953384215367E-04_wp, &
      & 3.06994281233379E-04_wp,-6.88833076815888E-05_wp,-2.72138600854253E-04_wp, &
      & 1.36166182734544E-04_wp,-5.84338260819475E-04_wp,-1.57316635327317E-05_wp, &
      &-5.09350857470944E-04_wp,-1.56609377614967E-04_wp, 3.69302881743411E-04_wp, &
      & 5.99043014817094E-04_wp,-1.24471946550559E-03_wp,-2.38428556173981E-03_wp, &
      & 1.06346755131592E-04_wp, 3.71126983773621E-04_wp, 1.82835861706770E-03_wp, &
      &-3.31406784596862E-03_wp, 1.20746182337580E-03_wp,-5.12575192802934E-03_wp, &
      &-2.27506888247999E-03_wp, 2.56487472738511E-04_wp, 8.24138929388585E-04_wp, &
      & 2.14802250375719E-03_wp,-1.14204591564980E-03_wp,-7.48265888202462E-04_wp, &
      & 6.87759158596898E-05_wp,-6.24456776946116E-05_wp,-1.39324488358490E-04_wp, &
      &-9.08774681400845E-05_wp,-1.21743448062207E-04_wp,-1.24410297693552E-03_wp, &
      & 2.00987442039049E-03_wp, 1.90274135264279E-04_wp,-4.94260832249692E-04_wp, &
      & 9.71754964880210E-05_wp,-6.08938528767931E-05_wp, 7.42984219798043E-05_wp, &
      &-7.84130470331478E-06_wp, 1.05720918721642E-04_wp, 1.55818848538031E-03_wp, &
      &-5.75745265195767E-03_wp,-1.49660914596034E-02_wp, 1.48330126334421E-02_wp, &
      & 3.37914579921590E-04_wp,-6.88833076815888E-05_wp, 5.89819010053817E-04_wp, &
      & 2.52240499515333E-04_wp,-4.23234562990438E-04_wp, 1.06457926932973E-03_wp, &
      &-4.42293296849883E-04_wp, 4.03614321206780E-05_wp, 1.00472117204186E-03_wp, &
      &-1.38995877590136E-03_wp, 9.36639684742331E-05_wp,-2.38179325663893E-03_wp, &
      & 4.96039250557894E-04_wp, 2.30301478670150E-03_wp,-2.11151989500910E-03_wp, &
      & 3.86716929434750E-04_wp,-2.97731898645141E-03_wp,-9.68246628988539E-04_wp, &
      &-6.39253233006904E-03_wp, 4.87712090937551E-04_wp, 2.33410259062064E-03_wp, &
      & 7.13070191665346E-04_wp, 2.15130914664620E-04_wp, 3.16310816330874E-04_wp, &
      & 2.35043460435739E-03_wp, 1.51322538300238E-04_wp,-4.09465025176957E-05_wp, &
      &-1.15338852384393E-04_wp, 3.16501856547077E-05_wp, 4.07538919273474E-05_wp, &
      &-1.84667925671337E-03_wp, 3.42911184709526E-03_wp, 3.95026807626506E-03_wp, &
      & 2.70908666784429E-03_wp,-1.46386688796714E-05_wp,-2.32989247203373E-04_wp, &
      &-6.15057484205189E-05_wp, 3.18695207242608E-05_wp, 2.47636898695574E-04_wp, &
      &-3.26712171644324E-03_wp,-8.97405825882184E-04_wp,-2.62733890723684E-02_wp, &
      &-1.41349407234096E-02_wp,-1.08472699745849E-04_wp,-2.72138600854253E-04_wp, &
      & 2.52240499515333E-04_wp, 1.08240225555745E-03_wp, 4.36693908059331E-04_wp, &
      & 1.20863288757621E-03_wp,-5.24583626459510E-04_wp, 4.69680812585431E-04_wp, &
      & 6.42112571717576E-04_wp,-2.37514750915699E-03_wp,-1.11598042499230E-03_wp, &
      & 4.90483620803742E-05_wp, 2.58030054106487E-03_wp, 5.27801624728434E-04_wp, &
      &-2.94549789451623E-03_wp,-3.09145824765568E-03_wp, 4.21572546220010E-03_wp, &
      &-1.52223779069168E-03_wp, 3.81725763353449E-03_wp, 2.67023870340793E-03_wp, &
      &-2.89952769295361E-04_wp, 1.14992353955420E-03_wp,-2.63692688715765E-04_wp, &
      & 1.51717062160637E-03_wp, 1.60863721335192E-03_wp, 3.34778995454745E-05_wp, &
      & 1.51524435468932E-04_wp, 8.54210495205125E-05_wp, 1.29899358262604E-04_wp, &
      & 4.44903889252788E-05_wp,-3.94588332180406E-04_wp, 1.23007486722984E-04_wp, &
      & 3.86843136045921E-03_wp, 1.56234233918496E-03_wp,-7.69755308471498E-05_wp, &
      &-1.72466828137293E-04_wp,-1.53448287155134E-04_wp, 1.13231975167471E-04_wp, &
      & 3.18217233255463E-05_wp,-3.45181082893583E-03_wp,-8.63422682880298E-03_wp, &
      &-1.06415646433602E-03_wp,-2.62925614731553E-02_wp, 7.46780711376324E-06_wp, &
      & 1.36166182734544E-04_wp,-4.23234562990438E-04_wp, 4.36693908059331E-04_wp, &
      & 9.39787253869348E-04_wp, 7.12005636844234E-03_wp,-4.51310995476668E-04_wp, &
      & 7.03700924044002E-04_wp,-2.94299682722306E-04_wp,-2.72611929639262E-03_wp, &
      & 6.83963834567585E-03_wp,-5.31303890635983E-03_wp, 4.86750935672929E-03_wp, &
      & 9.51752248608211E-04_wp,-2.43673258096656E-03_wp,-9.06588468766871E-03_wp, &
      &-3.79091435809209E-04_wp,-2.36936459767103E-03_wp,-3.07632771594462E-04_wp, &
      & 3.36923184494226E-03_wp, 8.99317619769735E-03_wp,-1.09159706780144E-03_wp, &
      & 3.20293364885497E-02_wp, 2.32797921511728E-02_wp, 2.11137416010826E-02_wp, &
      & 1.14014058698088E-03_wp, 1.14151258746189E-03_wp,-3.01493559574927E-04_wp, &
      & 7.61823290430089E-04_wp,-4.78022645283907E-04_wp,-1.33275714243650E-03_wp, &
      &-3.51386849120022E-02_wp,-1.39271250257825E-02_wp, 2.23819405753520E-02_wp, &
      &-1.29940451794042E-03_wp, 9.04158960138589E-04_wp,-5.82863081578004E-04_wp, &
      &-5.86305027750962E-04_wp,-5.20787359796830E-04_wp,-1.20082444862482E-03_wp, &
      & 1.01678825709110E-02_wp,-1.38192687806438E-02_wp,-4.08823673808927E-02_wp, &
      &-7.15273229046639E-04_wp,-2.62747591641126E-04_wp,-5.93120966859311E-04_wp, &
      & 1.02934074705867E-03_wp, 1.20857782799951E-03_wp,-4.51310995476668E-04_wp, &
      & 2.28710717524929E-03_wp, 5.41555843395809E-04_wp,-2.60454779703653E-04_wp, &
      &-2.64761188452365E-03_wp, 5.74936806185157E-04_wp,-4.73552549544657E-05_wp, &
      &-1.70452741731383E-03_wp,-2.03891627313829E-04_wp, 1.06857204336503E-03_wp, &
      & 1.23576920153875E-03_wp,-1.23906835103024E-03_wp,-2.62937535044638E-03_wp, &
      & 6.80225565224807E-04_wp,-2.77935804798826E-03_wp, 6.09598203037294E-04_wp, &
      & 7.15238672648235E-03_wp,-1.15189563998192E-02_wp,-1.67088324005519E-02_wp, &
      &-1.96392475312164E-02_wp,-8.08703676888246E-04_wp,-6.27557604748140E-04_wp, &
      & 1.20819061332534E-04_wp,-6.29338098715479E-04_wp,-7.69573481465763E-05_wp, &
      &-8.03102508565910E-03_wp,-1.18519097351434E-02_wp,-1.74731513337701E-02_wp, &
      & 2.21934553111561E-02_wp,-8.42367890453742E-04_wp, 6.54572917365500E-04_wp, &
      &-9.79001513299533E-05_wp,-7.38113854803026E-04_wp, 1.66609406392804E-04_wp, &
      & 2.55368882386110E-03_wp, 1.44149977964325E-02_wp, 4.82113496013389E-03_wp, &
      & 1.15441397373419E-02_wp,-3.91702493559700E-04_wp,-1.85560136617751E-04_wp, &
      & 2.73984936091776E-05_wp,-3.83899683641624E-04_wp,-5.51201130588527E-04_wp, &
      & 7.03700924044002E-04_wp, 5.41555843395809E-04_wp, 1.62759376797712E-03_wp, &
      & 3.60831805144608E-04_wp,-1.63343542947862E-03_wp,-1.08136358505995E-03_wp, &
      & 3.26048543245310E-04_wp,-6.81691710601174E-04_wp, 5.74633761645855E-04_wp, &
      &-1.94158932807316E-04_wp,-5.97682056710771E-04_wp, 1.49160479473176E-03_wp, &
      &-2.93069131723140E-03_wp, 1.04587547201306E-03_wp, 1.49439124494011E-04_wp, &
      & 6.08564565890287E-04_wp, 4.56133621083857E-03_wp,-1.88851874367886E-02_wp, &
      & 8.77992563214724E-03_wp,-1.26538770093818E-02_wp,-7.09842496552069E-04_wp, &
      & 2.78133818721783E-06_wp, 6.06385583335837E-04_wp, 3.06221366687789E-06_wp, &
      & 2.86568552326527E-04_wp,-3.81506806848197E-03_wp,-1.56903914048125E-02_wp, &
      & 9.92777410612344E-03_wp, 1.08044878880641E-02_wp,-6.58362610172832E-04_wp, &
      &-2.66475310005962E-05_wp,-5.84081941461458E-04_wp, 1.60898571495431E-05_wp, &
      &-1.95975356661681E-04_wp,-3.62978868139369E-03_wp, 5.13925269875120E-03_wp, &
      & 1.02205811651839E-02_wp,-1.81593535427986E-02_wp,-3.90294960024259E-04_wp, &
      & 9.37999421072550E-06_wp,-5.75461842467355E-04_wp,-4.74477134909322E-05_wp, &
      & 5.45960645253513E-04_wp,-2.94299682722306E-04_wp,-2.60454779703653E-04_wp, &
      & 3.60831805144608E-04_wp, 2.43634077876406E-03_wp,-1.06882707972060E-03_wp, &
      &-1.65294406463820E-03_wp,-7.33698398795506E-05_wp, 1.79929225703692E-03_wp, &
      & 2.13659141476864E-03_wp,-1.80821386649232E-03_wp, 2.37533963974315E-03_wp, &
      & 6.75003092400962E-04_wp,-1.81229609926642E-03_wp,-1.71562917162580E-03_wp, &
      & 2.35017040539554E-03_wp,-2.96433667057058E-04_wp, 4.73945317634163E-03_wp, &
      &-1.99644809784157E-02_wp,-1.13633087400760E-02_wp, 4.51839061161885E-03_wp, &
      &-2.45726954688940E-04_wp,-6.32157381091241E-04_wp, 8.19701431661231E-05_wp, &
      &-1.19768205556007E-04_wp, 7.02087727366023E-04_wp, 5.48748252491780E-03_wp, &
      & 2.07344238174008E-02_wp, 1.12782212772816E-02_wp, 2.59863950050453E-03_wp, &
      & 3.30597200248198E-04_wp,-6.86166560341841E-04_wp, 6.16125477932459E-05_wp, &
      & 2.02646773182076E-04_wp, 6.65298187726416E-04_wp,-9.07158739473559E-03_wp, &
      & 1.33594259106603E-02_wp,-1.99274174073717E-02_wp,-2.34198210219147E-02_wp, &
      &-6.68613238440462E-04_wp,-4.37632419523053E-04_wp,-1.19585665770434E-04_wp, &
      & 1.01952206317582E-03_wp, 5.93888771707889E-04_wp,-2.72611929639262E-03_wp, &
      &-2.64761188452365E-03_wp,-1.63343542947862E-03_wp,-1.06882707972060E-03_wp, &
      & 4.26628823325425E-02_wp, 5.49747786037204E-03_wp, 2.52276056753839E-03_wp, &
      & 7.36515555688035E-03_wp,-4.06093748220504E-04_wp,-3.22070270286136E-04_wp, &
      & 3.89590843248208E-03_wp, 4.64033408061550E-03_wp, 1.33984178956232E-02_wp, &
      & 7.59055322688277E-04_wp, 8.14808504820114E-04_wp,-2.97425133542682E-04_wp, &
      & 2.23748997638881E-02_wp, 6.59864384170869E-02_wp, 8.69171133301889E-02_wp, &
      & 1.24058683561526E-03_wp, 8.43408209350482E-04_wp, 2.57087685798228E-03_wp, &
      & 4.72280375210331E-04_wp, 1.34942381809520E-03_wp,-1.20043153488693E-03_wp, &
      &-2.39746385993084E-02_wp, 7.91327539671594E-02_wp, 6.88040159830834E-02_wp, &
      &-1.27823662557746E-02_wp, 1.67224821871277E-03_wp,-2.26628110062760E-03_wp, &
      & 1.02410080994947E-04_wp, 1.25422685863442E-03_wp, 1.23413770534980E-03_wp, &
      &-1.35647375904467E-02_wp, 7.56592855070334E-02_wp, 3.15195470481471E-02_wp, &
      & 6.45164556377193E-02_wp,-9.39640106534622E-04_wp,-6.34568272551357E-04_wp, &
      & 3.16393662007085E-04_wp,-1.37698019287476E-03_wp,-2.23650399659197E-03_wp, &
      & 6.83963834567585E-03_wp, 5.74936806185157E-04_wp,-1.08136358505995E-03_wp, &
      &-1.65294406463820E-03_wp, 5.49747786037204E-03_wp, 3.04747475756913E-02_wp, &
      & 5.04578105860948E-03_wp,-4.64579944366578E-03_wp,-7.19057126374986E-03_wp, &
      & 3.26703944092933E-03_wp,-1.16593175914037E-02_wp,-1.76046268027161E-03_wp, &
      & 1.32231539977624E-03_wp,-5.03958017389164E-03_wp,-8.83040300561631E-04_wp, &
      & 9.84460157129005E-03_wp, 2.07361060924392E-02_wp, 6.71002401640998E-02_wp, &
      &-1.96682582060302E-02_wp, 7.82218736578697E-02_wp, 2.84314414197796E-03_wp, &
      & 3.82184720899220E-04_wp,-1.93416453481232E-03_wp, 4.83644145237919E-04_wp, &
      &-5.30311714097208E-04_wp, 1.83841031236885E-02_wp,-6.46154022227629E-02_wp, &
      & 3.76013827348195E-02_wp, 6.76244364285764E-02_wp,-2.56243047957179E-03_wp, &
      & 3.58655742369494E-04_wp,-2.02294049636478E-03_wp,-3.48741562764879E-04_wp, &
      &-5.80912494235757E-04_wp,-5.21088322546849E-03_wp, 3.66481612496406E-02_wp, &
      &-1.62019398806568E-02_wp, 4.22691225405726E-02_wp,-7.23873827592798E-04_wp, &
      &-6.94871436601380E-04_wp, 6.92791431294828E-04_wp, 2.14712336521353E-04_wp, &
      &-1.19295873522822E-03_wp,-5.31303890635983E-03_wp,-4.73552549544657E-05_wp, &
      & 3.26048543245310E-04_wp,-7.33698398795506E-05_wp, 2.52276056753839E-03_wp, &
      & 5.04578105860948E-03_wp, 2.84678442250546E-02_wp, 3.55720765445210E-03_wp, &
      &-7.76386327641950E-04_wp, 3.40221326056398E-03_wp, 1.80491786167642E-03_wp, &
      &-1.74606241999143E-03_wp, 2.05720567011699E-04_wp,-9.55019540171823E-04_wp, &
      &-6.73908728045511E-04_wp,-1.24968339010139E-02_wp,-7.48356414874469E-03_wp, &
      & 1.50312072896683E-02_wp,-9.52279454445640E-02_wp, 1.07162323238982E-02_wp, &
      & 2.74271146447398E-04_wp,-2.28230451023414E-03_wp,-1.54891427594768E-03_wp, &
      &-1.52565406983092E-03_wp,-8.33157635163222E-05_wp,-9.69794575898898E-03_wp, &
      & 5.97329063417440E-04_wp, 8.53965235416169E-02_wp,-6.23157387188532E-03_wp, &
      & 2.39631392231160E-04_wp,-1.86189175940878E-03_wp,-1.22978136847514E-03_wp, &
      & 1.62253133956297E-03_wp,-2.57580865019132E-04_wp,-9.88315017273462E-03_wp, &
      &-6.24835765122572E-03_wp, 8.60672855726232E-02_wp, 3.14814661832583E-03_wp, &
      & 3.77982171071824E-04_wp, 9.34039516698034E-04_wp,-1.19873916027241E-03_wp, &
      &-2.29602069996173E-03_wp, 4.64419548681912E-05_wp, 4.86750935672929E-03_wp, &
      &-1.70452741731383E-03_wp,-6.81691710601174E-04_wp, 1.79929225703692E-03_wp, &
      & 7.36515555688035E-03_wp,-4.64579944366578E-03_wp, 3.55720765445210E-03_wp, &
      & 3.34375762557645E-02_wp, 4.74783621825158E-03_wp,-9.19753494895621E-03_wp, &
      &-4.80342888989154E-03_wp,-5.05146595247601E-03_wp, 6.81729766686798E-04_wp, &
      & 2.02662065988195E-03_wp, 9.45853097923228E-03_wp, 3.79403643541078E-03_wp, &
      & 1.40905798811335E-02_wp, 7.86252725507647E-02_wp,-1.34685683635130E-02_wp, &
      & 5.53966602554030E-03_wp, 1.51709969704829E-03_wp, 4.93883626873902E-04_wp, &
      &-1.33728432187232E-03_wp, 9.51928118476924E-06_wp,-1.49346423200351E-03_wp, &
      &-1.19096387956701E-02_wp, 7.76490796065693E-02_wp,-3.00304843081706E-02_wp, &
      &-2.01301060251737E-03_wp, 1.10201275650136E-03_wp,-1.69176114470135E-04_wp, &
      & 1.46303965116822E-03_wp,-6.46390044540534E-04_wp, 1.68271220553567E-03_wp, &
      & 2.13613158631452E-02_wp, 3.25299720293667E-02_wp, 4.67414968162166E-02_wp, &
      &-1.03864190087322E-01_wp,-1.50087725727626E-03_wp, 3.62419676815490E-05_wp, &
      &-2.39989501736821E-03_wp, 3.71610548046474E-04_wp, 2.60086011058214E-03_wp, &
      & 9.51752248608211E-04_wp,-2.03891627313829E-04_wp, 5.74633761645855E-04_wp, &
      & 2.13659141476864E-03_wp,-4.06093748220504E-04_wp,-7.19057126374986E-03_wp, &
      &-7.76386327641950E-04_wp, 4.74783621825158E-03_wp, 4.18714250311281E-02_wp, &
      &-2.52341070566793E-03_wp, 7.67210121001825E-04_wp, 2.11987419218247E-03_wp, &
      &-4.97769776555575E-03_wp,-5.83490623273465E-03_wp, 5.49862790914475E-03_wp, &
      &-4.23717011058369E-04_wp,-8.82914933642690E-03_wp, 3.37935036734026E-02_wp, &
      &-3.37085595218908E-02_wp,-8.79319171026842E-02_wp,-1.20607136069377E-03_wp, &
      &-4.47455992873397E-04_wp,-1.80831386425085E-04_wp,-1.30457397783449E-03_wp, &
      &-1.72606095327237E-03_wp,-8.62164655307693E-03_wp,-3.71397556081792E-02_wp, &
      & 3.53408959400390E-02_wp,-9.54603800516629E-02_wp, 1.04148539958071E-03_wp, &
      &-3.58853556247105E-04_wp,-3.58903980971413E-04_wp, 1.79457627302079E-03_wp, &
      &-2.19066150559248E-03_wp, 2.19877634531918E-02_wp, 7.20879088530275E-02_wp, &
      &-7.10369517108405E-02_wp,-4.72439825881444E-02_wp,-2.04362393783481E-03_wp, &
      &-1.49241291128878E-03_wp, 1.88746337135601E-04_wp, 2.32194941901059E-03_wp, &
      & 4.24534276701425E-04_wp,-2.43673258096656E-03_wp, 1.06857204336503E-03_wp, &
      &-1.94158932807316E-04_wp,-1.80821386649232E-03_wp,-3.22070270286136E-04_wp, &
      & 3.26703944092933E-03_wp, 3.40221326056398E-03_wp,-9.19753494895621E-03_wp, &
      &-2.52341070566793E-03_wp, 1.18754497209645E-02_wp, 3.12284761791560E-03_wp, &
      &-1.79483956434786E-03_wp, 1.08426244847884E-03_wp, 2.95149027174660E-03_wp, &
      & 9.39530660209130E-04_wp,-8.20551951959828E-04_wp, 2.51220067000251E-05_wp, &
      &-3.43256234107025E-02_wp,-1.00847266379413E-02_wp, 3.53479628929656E-02_wp, &
      & 3.47939388044115E-04_wp,-7.43623514605931E-04_wp, 3.68443102295418E-06_wp, &
      & 3.27464019601924E-04_wp, 1.58744013078058E-03_wp, 5.59681254123894E-04_wp, &
      &-4.32983265074521E-02_wp,-1.56787994955002E-02_wp,-3.91999035934349E-02_wp, &
      & 4.04261706357615E-04_wp, 1.24003860524734E-03_wp, 2.73469045608191E-04_wp, &
      & 2.83112214438342E-04_wp,-2.02065478669260E-03_wp,-8.80007601248442E-04_wp, &
      & 3.13091576817181E-02_wp, 3.67758744262592E-02_wp, 5.10602001473590E-02_wp, &
      &-6.95430188414474E-04_wp,-1.53919367966631E-04_wp,-6.53034622285015E-05_wp, &
      &-2.04327901205466E-03_wp,-1.86505173435443E-03_wp,-9.06588468766871E-03_wp, &
      & 1.23576920153875E-03_wp,-5.97682056710771E-04_wp, 2.37533963974315E-03_wp, &
      & 3.89590843248208E-03_wp,-1.16593175914037E-02_wp, 1.80491786167642E-03_wp, &
      &-4.80342888989154E-03_wp, 7.67210121001825E-04_wp, 3.12284761791560E-03_wp, &
      & 1.53464378090052E-02_wp, 4.96985215764110E-04_wp, 1.15257059393314E-03_wp, &
      & 3.14378857274412E-04_wp,-1.87096617631923E-03_wp,-8.95788759613569E-03_wp, &
      &-9.71281399247800E-04_wp,-6.97162391500141E-02_wp,-2.42710173825222E-02_wp, &
      &-2.43258124412676E-02_wp,-1.88935364209727E-03_wp,-1.84650957848923E-03_wp, &
      & 8.39626386203248E-04_wp,-8.32108813509586E-04_wp, 1.57482645531049E-03_wp, &
      &-6.34506035159502E-04_wp, 6.51239650221619E-02_wp, 1.21437671288196E-04_wp, &
      &-2.51829912413059E-02_wp, 1.91223161721948E-03_wp,-1.05936663656557E-03_wp, &
      & 1.27335881729851E-03_wp, 3.48810273263677E-04_wp, 1.31361034587710E-03_wp, &
      & 5.34558322462221E-05_wp, 2.45085826562595E-02_wp,-2.56042002348986E-03_wp, &
      & 4.06759272671501E-02_wp,-4.26658307535336E-04_wp,-4.73464627516159E-04_wp, &
      & 6.75174714704217E-04_wp,-5.61493981444976E-04_wp,-1.60711005567859E-03_wp, &
      &-3.79091435809209E-04_wp,-1.23906835103024E-03_wp, 1.49160479473176E-03_wp, &
      & 6.75003092400962E-04_wp, 4.64033408061550E-03_wp,-1.76046268027161E-03_wp, &
      &-1.74606241999143E-03_wp,-5.05146595247601E-03_wp, 2.11987419218247E-03_wp, &
      &-1.79483956434786E-03_wp, 4.96985215764110E-04_wp, 5.57587375664243E-03_wp, &
      &-1.29796423437543E-04_wp,-7.71055759065487E-04_wp,-1.85519945370246E-04_wp, &
      &-9.53869678302510E-04_wp, 3.34456557574598E-05_wp,-2.14763141592654E-02_wp, &
      & 4.22704175886859E-02_wp,-1.51686485770944E-02_wp,-7.82696633809294E-04_wp, &
      & 8.26114803513319E-04_wp, 1.18017442302747E-03_wp, 5.55319395078278E-04_wp, &
      & 3.30610124385097E-04_wp,-2.20043853657731E-04_wp,-3.86345196809904E-03_wp, &
      & 4.91636332105934E-02_wp, 3.71602880165389E-04_wp,-1.49962559622102E-04_wp, &
      &-1.27713933502452E-03_wp,-1.02227528479917E-03_wp, 9.54381580813492E-04_wp, &
      &-4.76500881652060E-05_wp, 5.47852591736715E-04_wp,-6.37022248510328E-03_wp, &
      &-1.46041445381265E-02_wp,-4.95678503503976E-03_wp, 1.86282736000768E-04_wp, &
      &-2.24409080145977E-05_wp, 1.12001295917723E-04_wp, 5.82371418669850E-04_wp, &
      & 4.18386118275000E-04_wp,-2.36936459767103E-03_wp,-2.62937535044638E-03_wp, &
      &-2.93069131723140E-03_wp,-1.81229609926642E-03_wp, 1.33984178956232E-02_wp, &
      & 1.32231539977624E-03_wp, 2.05720567011699E-04_wp, 6.81729766686798E-04_wp, &
      &-4.97769776555575E-03_wp, 1.08426244847884E-03_wp, 1.15257059393314E-03_wp, &
      &-1.29796423437543E-04_wp, 9.94185681796658E-03_wp,-7.18121111829918E-05_wp, &
      &-1.96697334064887E-04_wp,-2.73593266925918E-03_wp, 1.26668468221258E-03_wp, &
      & 4.19015258851335E-02_wp, 2.73109370928081E-02_wp, 2.86043328592255E-02_wp, &
      & 1.52385712756030E-03_wp, 1.17176473508348E-03_wp,-5.69308150559547E-04_wp, &
      & 8.01134071944117E-04_wp,-5.92916919148144E-04_wp,-8.92849680640423E-04_wp, &
      & 4.00952178344970E-02_wp, 6.11598531351291E-03_wp,-3.43458403737120E-02_wp, &
      & 1.81169955224588E-03_wp,-6.03449935835499E-04_wp, 9.59524078160904E-04_wp, &
      & 6.67148018086848E-04_wp, 1.39319355064774E-04_wp,-8.42490345787972E-04_wp, &
      &-1.91489694891317E-02_wp, 5.62405201973460E-03_wp, 4.94699911982638E-02_wp, &
      & 1.30383854676361E-03_wp, 4.09728762394107E-04_wp, 9.57203891809881E-04_wp, &
      &-7.51198337433286E-04_wp,-1.22239364856491E-03_wp,-3.07632771594462E-04_wp, &
      & 6.80225565224807E-04_wp, 1.04587547201306E-03_wp,-1.71562917162580E-03_wp, &
      & 7.59055322688277E-04_wp,-5.03958017389164E-03_wp,-9.55019540171823E-04_wp, &
      & 2.02662065988195E-03_wp,-5.83490623273465E-03_wp, 2.95149027174660E-03_wp, &
      & 3.14378857274412E-04_wp,-7.71055759065487E-04_wp,-7.18121111829918E-05_wp, &
      & 6.02070086335599E-03_wp,-1.97516410409185E-04_wp, 3.97735945558793E-04_wp, &
      & 2.34648509825606E-05_wp,-1.51973938081867E-02_wp, 2.87609703820858E-02_wp, &
      &-9.50155906323720E-03_wp,-5.48957432010783E-04_wp, 5.53339286205348E-04_wp, &
      & 8.02603107385276E-04_wp, 3.77924222566519E-04_wp, 1.95211507047451E-04_wp, &
      & 6.69838181723184E-04_wp,-5.69855677424794E-03_wp,-3.24906186581363E-02_wp, &
      &-4.84839911255815E-03_wp,-6.17444284658637E-05_wp, 1.02522733289301E-03_wp, &
      & 4.82528207292135E-04_wp,-5.52940806333149E-04_wp,-4.61445881759826E-04_wp, &
      &-4.40961160825250E-04_wp, 1.09352791492342E-03_wp, 5.78051213871221E-02_wp, &
      &-1.68963645227524E-03_wp,-1.50061542669221E-04_wp, 5.19503852174865E-04_wp, &
      &-1.09043230139903E-03_wp,-1.81994295729175E-03_wp,-3.20961296343504E-05_wp, &
      & 3.36923184494226E-03_wp,-2.77935804798826E-03_wp, 1.49439124494011E-04_wp, &
      & 2.35017040539554E-03_wp, 8.14808504820114E-04_wp,-8.83040300561631E-04_wp, &
      &-6.73908728045511E-04_wp, 9.45853097923228E-03_wp, 5.49862790914475E-03_wp, &
      & 9.39530660209130E-04_wp,-1.87096617631923E-03_wp,-1.85519945370246E-04_wp, &
      &-1.96697334064887E-04_wp,-1.97516410409185E-04_wp, 1.14266891748387E-02_wp, &
      & 3.82848567040289E-03_wp, 4.43421419479474E-04_wp,-4.12246471602370E-03_wp, &
      & 9.67883281895291E-03_wp, 5.65375190106829E-02_wp, 1.60363985580091E-03_wp, &
      & 1.52137011719090E-04_wp,-3.37196482824384E-04_wp, 1.19561502927330E-03_wp, &
      & 1.49002901462518E-03_wp,-9.35607155517462E-04_wp, 7.87103418955764E-03_wp, &
      & 3.89420518191703E-03_wp,-4.71640208254581E-02_wp, 1.35122588661134E-03_wp, &
      &-2.00224744625050E-04_wp, 5.23077500267552E-04_wp, 8.65873402966272E-04_wp, &
      &-8.10924344558873E-04_wp, 1.04269735966288E-03_wp, 3.42456107696032E-02_wp, &
      &-2.45166757199998E-03_wp,-6.15903243688644E-02_wp,-1.62313671252700E-03_wp, &
      &-5.52953181551917E-04_wp,-1.20152015219740E-03_wp, 1.02978296745931E-03_wp, &
      & 1.53238668144887E-03_wp, 8.99317619769735E-03_wp, 6.09598203037294E-04_wp, &
      & 6.08564565890287E-04_wp,-2.96433667057058E-04_wp,-2.97425133542682E-04_wp, &
      & 9.84460157129005E-03_wp,-1.24968339010139E-02_wp, 3.79403643541078E-03_wp, &
      &-4.23717011058369E-04_wp,-8.20551951959828E-04_wp,-8.95788759613569E-03_wp, &
      &-9.53869678302510E-04_wp,-2.73593266925918E-03_wp, 3.97735945558793E-04_wp, &
      & 3.82848567040289E-03_wp, 1.48313494807019E-02_wp, 2.13909072727687E-04_wp, &
      & 2.60340640603924E-02_wp, 4.55835040340949E-02_wp, 3.11676865563128E-02_wp, &
      & 1.24738784694517E-03_wp, 1.53260469644755E-03_wp, 1.46419590702231E-06_wp, &
      & 1.22724071642238E-03_wp,-1.01556221394515E-04_wp, 6.87853764950857E-04_wp, &
      &-3.92830274661842E-02_wp,-4.33396652146917E-02_wp, 3.84108882966145E-02_wp, &
      &-1.77145693527340E-03_wp, 1.76952445247401E-03_wp,-2.92072966954730E-04_wp, &
      &-1.42087387638032E-03_wp,-3.57640518895915E-04_wp, 2.43800444355134E-04_wp, &
      & 4.74547939918725E-02_wp,-2.86052824976296E-02_wp,-2.74882135918198E-02_wp, &
      &-1.80961917041824E-03_wp,-9.88501127912692E-04_wp,-3.96125716686334E-04_wp, &
      & 1.22747256066821E-03_wp, 4.59125102539591E-04_wp,-1.09159706780144E-03_wp, &
      & 7.15238672648235E-03_wp, 4.56133621083857E-03_wp, 4.73945317634163E-03_wp, &
      & 2.23748997638881E-02_wp, 2.07361060924392E-02_wp,-7.48356414874469E-03_wp, &
      & 1.40905798811335E-02_wp,-8.82914933642690E-03_wp, 2.51220067000251E-05_wp, &
      &-9.71281399247800E-04_wp, 3.34456557574598E-05_wp, 1.26668468221258E-03_wp, &
      & 2.34648509825606E-05_wp, 4.43421419479474E-04_wp, 2.13909072727687E-04_wp, &
      & 9.79896480854331E-01_wp,-6.57322957596467E-02_wp,-4.15963753309213E-02_wp, &
      &-4.43149361693328E-02_wp,-3.26277967081380E-03_wp,-3.02728686916666E-03_wp, &
      & 1.02915003219649E-03_wp,-2.06316441181171E-03_wp, 1.29122305576298E-03_wp, &
      &-1.73472826985636E-03_wp,-4.78997789547741E-03_wp,-4.81001243016759E-03_wp, &
      &-2.26440560771941E-02_wp, 1.03307475400793E-03_wp, 8.28947276525380E-04_wp, &
      & 8.36491920736309E-04_wp, 3.56449278824868E-04_wp,-2.26774951942169E-03_wp, &
      &-4.67636429378859E-04_wp,-2.15271159950398E-02_wp,-5.04314540727691E-03_wp, &
      & 2.25472983416135E-03_wp, 2.05817627414404E-03_wp, 5.05903385749172E-04_wp, &
      & 7.91175763385960E-04_wp, 7.16658955566675E-04_wp, 1.18813026352344E-03_wp, &
      & 3.20293364885497E-02_wp,-1.15189563998192E-02_wp,-1.88851874367886E-02_wp, &
      &-1.99644809784157E-02_wp, 6.59864384170869E-02_wp, 6.71002401640998E-02_wp, &
      & 1.50312072896683E-02_wp, 7.86252725507647E-02_wp, 3.37935036734026E-02_wp, &
      &-3.43256234107025E-02_wp,-6.97162391500141E-02_wp,-2.14763141592654E-02_wp, &
      & 4.19015258851335E-02_wp,-1.51973938081867E-02_wp,-4.12246471602370E-03_wp, &
      & 2.60340640603924E-02_wp,-6.57322957596467E-02_wp, 8.97759129210232E-01_wp, &
      &-9.01117356106460E-03_wp,-1.45174747420640E-02_wp, 1.54984336912709E-02_wp, &
      & 1.20456184521357E-02_wp,-1.34358656069674E-02_wp,-3.45683269750024E-04_wp, &
      &-2.38133130559744E-02_wp, 5.97235898784159E-03_wp, 3.58761398700636E-03_wp, &
      &-1.13532557461847E-02_wp,-9.58740505511901E-03_wp,-6.72831985359694E-04_wp, &
      & 2.30563240863010E-03_wp, 4.32903388983687E-04_wp,-1.03486087879590E-03_wp, &
      &-9.31983071527857E-04_wp,-1.75697167675842E-02_wp,-4.80545714771638E-02_wp, &
      &-4.04120386990838E-03_wp, 2.91141295701072E-02_wp, 5.40940849482020E-03_wp, &
      & 2.05121472162809E-03_wp, 2.29185700862415E-03_wp, 4.08333274592717E-04_wp, &
      &-4.29139345849666E-04_wp, 2.32797921511728E-02_wp,-1.67088324005519E-02_wp, &
      & 8.77992563214724E-03_wp,-1.13633087400760E-02_wp, 8.69171133301889E-02_wp, &
      &-1.96682582060302E-02_wp,-9.52279454445640E-02_wp,-1.34685683635130E-02_wp, &
      &-3.37085595218908E-02_wp,-1.00847266379413E-02_wp,-2.42710173825222E-02_wp, &
      & 4.22704175886859E-02_wp, 2.73109370928081E-02_wp, 2.87609703820858E-02_wp, &
      & 9.67883281895291E-03_wp, 4.55835040340949E-02_wp,-4.15963753309213E-02_wp, &
      &-9.01117356106460E-03_wp, 8.98867435484964E-01_wp,-5.85200876311569E-03_wp, &
      &-4.20573766553356E-05_wp, 2.32344402256585E-02_wp, 1.40802879004261E-02_wp, &
      & 1.57179505560624E-02_wp, 4.80633106127915E-05_wp, 1.29675551059358E-03_wp, &
      &-1.01166046687542E-02_wp, 1.64479228406389E-02_wp,-1.13616600785127E-02_wp, &
      &-6.58277179844695E-04_wp,-9.90104964157021E-05_wp,-1.69338120880881E-03_wp, &
      & 8.97862798560837E-04_wp,-2.17183687326237E-03_wp, 1.74301293583753E-03_wp, &
      &-1.23079993061857E-02_wp, 1.53478202076274E-02_wp,-5.76126756443682E-03_wp, &
      & 6.00834384066742E-04_wp, 6.74812958013988E-04_wp,-1.66126666824316E-03_wp, &
      &-2.78697676847693E-04_wp, 2.05969982294557E-03_wp, 2.11137416010826E-02_wp, &
      &-1.96392475312164E-02_wp,-1.26538770093818E-02_wp, 4.51839061161885E-03_wp, &
      & 1.24058683561526E-03_wp, 7.82218736578697E-02_wp, 1.07162323238982E-02_wp, &
      & 5.53966602554030E-03_wp,-8.79319171026842E-02_wp, 3.53479628929656E-02_wp, &
      &-2.43258124412676E-02_wp,-1.51686485770944E-02_wp, 2.86043328592255E-02_wp, &
      &-9.50155906323720E-03_wp, 5.65375190106829E-02_wp, 3.11676865563128E-02_wp, &
      &-4.43149361693328E-02_wp,-1.45174747420640E-02_wp,-5.85200876311569E-03_wp, &
      & 9.09495449859558E-01_wp, 2.36136986030882E-02_wp,-3.50464496910378E-04_wp, &
      &-9.07663625892514E-03_wp, 1.23559429245915E-02_wp, 1.65996348740507E-02_wp, &
      &-2.28367325786917E-02_wp, 1.33600671759115E-02_wp, 1.41674515334049E-03_wp, &
      &-5.45952997012663E-02_wp, 5.34526356976743E-03_wp,-4.98434975453687E-06_wp, &
      & 2.40242802702843E-03_wp, 2.48340309494439E-03_wp,-2.72169849948320E-03_wp, &
      & 1.34726136560787E-02_wp, 1.29098424984618E-02_wp,-9.80559410487820E-03_wp, &
      &-5.05472700198872E-03_wp,-2.40189133071662E-03_wp,-1.16055916082677E-03_wp, &
      &-5.81589521031109E-04_wp, 2.45136518477340E-03_wp, 1.40102242660422E-03_wp, &
      & 1.14014058698088E-03_wp,-8.08703676888246E-04_wp,-7.09842496552069E-04_wp, &
      &-2.45726954688940E-04_wp, 8.43408209350482E-04_wp, 2.84314414197796E-03_wp, &
      & 2.74271146447398E-04_wp, 1.51709969704829E-03_wp,-1.20607136069377E-03_wp, &
      & 3.47939388044115E-04_wp,-1.88935364209727E-03_wp,-7.82696633809294E-04_wp, &
      & 1.52385712756030E-03_wp,-5.48957432010783E-04_wp, 1.60363985580091E-03_wp, &
      & 1.24738784694517E-03_wp,-3.26277967081380E-03_wp, 1.54984336912709E-02_wp, &
      &-4.20573766553356E-05_wp, 2.36136986030882E-02_wp, 9.16664335908820E-04_wp, &
      & 2.19194016483456E-04_wp,-4.79576529056885E-04_wp, 3.36656008799403E-04_wp, &
      & 1.96148834151238E-05_wp,-9.41200990317532E-04_wp,-2.88523880443020E-04_wp, &
      &-1.27997465783394E-03_wp,-5.00411339353546E-03_wp, 1.98057688893439E-04_wp, &
      & 6.86374428440874E-05_wp, 1.12832721024584E-04_wp, 8.04650240365074E-05_wp, &
      &-1.60848310539632E-04_wp, 9.23735830986981E-04_wp,-1.64728596302177E-03_wp, &
      &-1.88915026257937E-03_wp,-5.38651668362386E-04_wp, 6.00370689453202E-05_wp, &
      & 1.06799866488623E-05_wp, 4.33426688768389E-05_wp, 1.15971638594665E-04_wp, &
      & 5.65401762358019E-05_wp, 1.14151258746189E-03_wp,-6.27557604748140E-04_wp, &
      & 2.78133818721783E-06_wp,-6.32157381091241E-04_wp, 2.57087685798228E-03_wp, &
      & 3.82184720899220E-04_wp,-2.28230451023414E-03_wp, 4.93883626873902E-04_wp, &
      &-4.47455992873397E-04_wp,-7.43623514605931E-04_wp,-1.84650957848923E-03_wp, &
      & 8.26114803513319E-04_wp, 1.17176473508348E-03_wp, 5.53339286205348E-04_wp, &
      & 1.52137011719090E-04_wp, 1.53260469644755E-03_wp,-3.02728686916666E-03_wp, &
      & 1.20456184521357E-02_wp, 2.32344402256585E-02_wp,-3.50464496910378E-04_wp, &
      & 2.19194016483456E-04_wp, 7.84936468283224E-04_wp, 1.80915755471877E-04_wp, &
      & 4.08650652107623E-04_wp,-3.38274598597018E-04_wp, 7.13200633425615E-04_wp, &
      &-2.49071515250449E-03_wp, 1.68275851590183E-04_wp,-7.40984829435843E-04_wp, &
      &-5.87717634459168E-05_wp, 5.60579840864022E-05_wp,-6.75347477155226E-05_wp, &
      & 1.75895781302209E-05_wp,-1.27743776572833E-04_wp,-3.13016321904110E-04_wp, &
      &-3.56800269797281E-03_wp, 2.99814266890537E-04_wp,-6.76380286992065E-04_wp, &
      & 1.50521404520208E-04_wp, 8.35092699986451E-05_wp,-1.78734093583186E-05_wp, &
      & 5.07786354518432E-06_wp, 9.85537457449779E-05_wp,-3.01493559574927E-04_wp, &
      & 1.20819061332534E-04_wp, 6.06385583335837E-04_wp, 8.19701431661231E-05_wp, &
      & 4.72280375210331E-04_wp,-1.93416453481232E-03_wp,-1.54891427594768E-03_wp, &
      &-1.33728432187232E-03_wp,-1.80831386425085E-04_wp, 3.68443102295418E-06_wp, &
      & 8.39626386203248E-04_wp, 1.18017442302747E-03_wp,-5.69308150559547E-04_wp, &
      & 8.02603107385276E-04_wp,-3.37196482824384E-04_wp, 1.46419590702231E-06_wp, &
      & 1.02915003219649E-03_wp,-1.34358656069674E-02_wp, 1.40802879004261E-02_wp, &
      &-9.07663625892514E-03_wp,-4.79576529056885E-04_wp, 1.80915755471877E-04_wp, &
      & 5.17997801064808E-04_wp, 1.22195967332093E-04_wp, 1.92582788745741E-04_wp, &
      & 6.94244426429266E-04_wp,-7.43888299800962E-04_wp, 1.89961114903928E-03_wp, &
      & 1.41680506335413E-03_wp,-8.73158452291938E-05_wp,-6.27090831700251E-05_wp, &
      &-9.95063274235428E-05_wp, 1.69168895582038E-05_wp, 1.01932389773469E-05_wp, &
      & 6.86048755790714E-04_wp, 1.12236019091068E-03_wp, 1.84140089469675E-03_wp, &
      &-1.05661020220776E-03_wp,-7.39580793665369E-05_wp,-5.70885224633915E-06_wp, &
      &-9.38004629219951E-05_wp,-6.43822060302723E-05_wp, 4.02430151216053E-05_wp, &
      & 7.61823290430089E-04_wp,-6.29338098715479E-04_wp, 3.06221366687789E-06_wp, &
      &-1.19768205556007E-04_wp, 1.34942381809520E-03_wp, 4.83644145237919E-04_wp, &
      &-1.52565406983092E-03_wp, 9.51928118476924E-06_wp,-1.30457397783449E-03_wp, &
      & 3.27464019601924E-04_wp,-8.32108813509586E-04_wp, 5.55319395078278E-04_wp, &
      & 8.01134071944117E-04_wp, 3.77924222566519E-04_wp, 1.19561502927330E-03_wp, &
      & 1.22724071642238E-03_wp,-2.06316441181171E-03_wp,-3.45683269750024E-04_wp, &
      & 1.57179505560624E-02_wp, 1.23559429245915E-02_wp, 3.36656008799403E-04_wp, &
      & 4.08650652107623E-04_wp, 1.22195967332093E-04_wp, 4.58450888737387E-04_wp, &
      & 2.38862168100573E-04_wp,-7.44842411579835E-04_wp,-8.81279754289441E-04_wp, &
      & 4.46531638551081E-04_wp,-3.31795032496775E-03_wp, 9.74089459436360E-05_wp, &
      & 4.58191515709555E-06_wp, 8.56172572143386E-06_wp, 8.84178196668619E-05_wp, &
      &-1.38692150236970E-04_wp, 9.28641711068648E-04_wp, 1.20019581993013E-04_wp, &
      & 5.49071489883944E-06_wp,-2.03354621399353E-03_wp,-4.17039948838420E-05_wp, &
      &-5.12570162749364E-06_wp,-6.69133655858396E-05_wp, 5.12877509169561E-05_wp, &
      & 1.03773290535138E-04_wp,-4.78022645283907E-04_wp,-7.69573481465763E-05_wp, &
      & 2.86568552326527E-04_wp, 7.02087727366023E-04_wp,-1.20043153488693E-03_wp, &
      &-5.30311714097208E-04_wp,-8.33157635163222E-05_wp,-1.49346423200351E-03_wp, &
      &-1.72606095327237E-03_wp, 1.58744013078058E-03_wp, 1.57482645531049E-03_wp, &
      & 3.30610124385097E-04_wp,-5.92916919148144E-04_wp, 1.95211507047451E-04_wp, &
      & 1.49002901462518E-03_wp,-1.01556221394515E-04_wp, 1.29122305576298E-03_wp, &
      &-2.38133130559744E-02_wp, 4.80633106127915E-05_wp, 1.65996348740507E-02_wp, &
      & 1.96148834151238E-05_wp,-3.38274598597018E-04_wp, 1.92582788745741E-04_wp, &
      & 2.38862168100573E-04_wp, 9.54162694066734E-04_wp,-2.25215401884667E-03_wp, &
      & 1.15317863417743E-03_wp, 1.46641094619476E-03_wp,-2.97128161929302E-03_wp, &
      & 1.74461165929101E-04_wp,-9.53916547354364E-05_wp, 4.56136607798075E-05_wp, &
      & 1.19018583406337E-04_wp,-4.22196797139872E-05_wp, 2.15520665339914E-03_wp, &
      & 5.25039772955339E-03_wp, 1.28691457765932E-05_wp,-2.17948090039576E-03_wp, &
      &-2.95522304533291E-04_wp,-1.25578343830325E-04_wp,-1.08474701764776E-04_wp, &
      & 5.08840121951604E-05_wp, 3.72930282297249E-05_wp,-1.33275714243650E-03_wp, &
      &-8.03102508565910E-03_wp,-3.81506806848197E-03_wp, 5.48748252491780E-03_wp, &
      &-2.39746385993084E-02_wp, 1.83841031236885E-02_wp,-9.69794575898898E-03_wp, &
      &-1.19096387956701E-02_wp,-8.62164655307693E-03_wp, 5.59681254123894E-04_wp, &
      &-6.34506035159502E-04_wp,-2.20043853657731E-04_wp,-8.92849680640423E-04_wp, &
      & 6.69838181723184E-04_wp,-9.35607155517462E-04_wp, 6.87853764950857E-04_wp, &
      &-1.73472826985636E-03_wp, 5.97235898784159E-03_wp, 1.29675551059358E-03_wp, &
      &-2.28367325786917E-02_wp,-9.41200990317532E-04_wp, 7.13200633425615E-04_wp, &
      & 6.94244426429266E-04_wp,-7.44842411579835E-04_wp,-2.25215401884667E-03_wp, &
      & 9.81042080440019E-01_wp, 6.76631480778566E-02_wp, 3.31263507473601E-02_wp, &
      &-4.62192368725004E-02_wp, 3.56676285631852E-03_wp,-2.66831117200006E-03_wp, &
      & 1.45049683061849E-03_wp, 1.78512184150084E-03_wp, 1.32158515813909E-03_wp, &
      &-3.02443053170287E-03_wp, 1.59234232521478E-02_wp, 1.87596405543643E-02_wp, &
      & 4.69726399080004E-03_wp,-1.69035008417164E-03_wp,-6.62367944133660E-04_wp, &
      &-1.25089500655171E-03_wp,-1.79892475752792E-03_wp,-3.50141084196518E-04_wp, &
      &-3.51386849120022E-02_wp,-1.18519097351434E-02_wp,-1.56903914048125E-02_wp, &
      & 2.07344238174008E-02_wp, 7.91327539671594E-02_wp,-6.46154022227629E-02_wp, &
      & 5.97329063417440E-04_wp, 7.76490796065693E-02_wp,-3.71397556081792E-02_wp, &
      &-4.32983265074521E-02_wp, 6.51239650221619E-02_wp,-3.86345196809904E-03_wp, &
      & 4.00952178344970E-02_wp,-5.69855677424794E-03_wp, 7.87103418955764E-03_wp, &
      &-3.92830274661842E-02_wp,-4.78997789547741E-03_wp, 3.58761398700636E-03_wp, &
      &-1.01166046687542E-02_wp, 1.33600671759115E-02_wp,-2.88523880443020E-04_wp, &
      &-2.49071515250449E-03_wp,-7.43888299800962E-04_wp,-8.81279754289441E-04_wp, &
      & 1.15317863417743E-03_wp, 6.76631480778566E-02_wp, 8.93154842519609E-01_wp, &
      &-4.10700162646009E-03_wp, 1.29570327653019E-02_wp, 1.67088054396957E-02_wp, &
      &-1.36859171542851E-02_wp, 1.27380822222116E-02_wp,-3.28692789808616E-05_wp, &
      & 2.27063683650581E-02_wp, 9.65825016803663E-03_wp,-5.19049786298460E-03_wp, &
      &-3.84519947027201E-02_wp,-1.24584223110816E-02_wp, 1.38574000819706E-03_wp, &
      &-6.85719533754288E-04_wp, 2.13681379765512E-03_wp, 3.54363522635215E-03_wp, &
      &-6.38061625562823E-06_wp,-1.39271250257825E-02_wp,-1.74731513337701E-02_wp, &
      & 9.92777410612344E-03_wp, 1.12782212772816E-02_wp, 6.88040159830834E-02_wp, &
      & 3.76013827348195E-02_wp, 8.53965235416169E-02_wp,-3.00304843081706E-02_wp, &
      & 3.53408959400390E-02_wp,-1.56787994955002E-02_wp, 1.21437671288196E-04_wp, &
      & 4.91636332105934E-02_wp, 6.11598531351291E-03_wp,-3.24906186581363E-02_wp, &
      & 3.89420518191703E-03_wp,-4.33396652146917E-02_wp,-4.81001243016759E-03_wp, &
      &-1.13532557461847E-02_wp, 1.64479228406389E-02_wp, 1.41674515334049E-03_wp, &
      &-1.27997465783394E-03_wp, 1.68275851590183E-04_wp, 1.89961114903928E-03_wp, &
      & 4.46531638551081E-04_wp, 1.46641094619476E-03_wp, 3.31263507473601E-02_wp, &
      &-4.10700162646009E-03_wp, 9.08224732889478E-01_wp, 9.12433476671870E-03_wp, &
      &-3.41613247939058E-04_wp,-2.27460176280414E-02_wp,-1.60912589585535E-02_wp, &
      & 1.70166983094865E-02_wp, 3.01226291286036E-04_wp, 1.87142472813712E-02_wp, &
      &-2.32116276965120E-02_wp,-2.18422956719651E-02_wp,-3.26436278753523E-02_wp, &
      & 1.96946150443605E-03_wp, 5.26688252429663E-04_wp,-1.74729268447544E-05_wp, &
      & 3.65712487530910E-03_wp, 4.03104249035507E-03_wp, 2.23819405753520E-02_wp, &
      & 2.21934553111561E-02_wp, 1.08044878880641E-02_wp, 2.59863950050453E-03_wp, &
      &-1.27823662557746E-02_wp, 6.76244364285764E-02_wp,-6.23157387188532E-03_wp, &
      &-2.01301060251737E-03_wp,-9.54603800516629E-02_wp,-3.91999035934349E-02_wp, &
      &-2.51829912413059E-02_wp, 3.71602880165389E-04_wp,-3.43458403737120E-02_wp, &
      &-4.84839911255815E-03_wp,-4.71640208254581E-02_wp, 3.84108882966145E-02_wp, &
      &-2.26440560771941E-02_wp,-9.58740505511901E-03_wp,-1.13616600785127E-02_wp, &
      &-5.45952997012663E-02_wp,-5.00411339353546E-03_wp,-7.40984829435843E-04_wp, &
      & 1.41680506335413E-03_wp,-3.31795032496775E-03_wp,-2.97128161929302E-03_wp, &
      &-4.62192368725004E-02_wp, 1.29570327653019E-02_wp, 9.12433476671870E-03_wp, &
      & 9.09747668217021E-01_wp,-2.27488865780810E-02_wp,-5.32740463510600E-04_wp, &
      &-1.00490947362311E-02_wp,-1.37415971532467E-02_wp, 1.76756468802079E-02_wp, &
      & 1.35400690416689E-02_wp, 8.14153009473268E-03_wp,-1.12403003698801E-02_wp, &
      &-3.83955909200850E-03_wp,-1.74863768724042E-03_wp,-9.80560541161186E-04_wp, &
      &-3.97380267628569E-04_wp, 2.70852057675908E-03_wp, 1.45358953653843E-03_wp, &
      &-1.29940451794042E-03_wp,-8.42367890453742E-04_wp,-6.58362610172832E-04_wp, &
      & 3.30597200248198E-04_wp, 1.67224821871277E-03_wp,-2.56243047957179E-03_wp, &
      & 2.39631392231160E-04_wp, 1.10201275650136E-03_wp, 1.04148539958071E-03_wp, &
      & 4.04261706357615E-04_wp, 1.91223161721948E-03_wp,-1.49962559622102E-04_wp, &
      & 1.81169955224588E-03_wp,-6.17444284658637E-05_wp, 1.35122588661134E-03_wp, &
      &-1.77145693527340E-03_wp, 1.03307475400793E-03_wp,-6.72831985359694E-04_wp, &
      &-6.58277179844695E-04_wp, 5.34526356976743E-03_wp, 1.98057688893439E-04_wp, &
      &-5.87717634459168E-05_wp,-8.73158452291938E-05_wp, 9.74089459436360E-05_wp, &
      & 1.74461165929101E-04_wp, 3.56676285631852E-03_wp, 1.67088054396957E-02_wp, &
      &-3.41613247939058E-04_wp,-2.27488865780810E-02_wp, 9.21471327475897E-04_wp, &
      &-2.50223062018044E-04_wp, 5.11478452729843E-04_wp, 3.56651216258034E-04_wp, &
      &-2.01067903244126E-05_wp,-1.11114144302852E-03_wp,-1.77564006933414E-03_wp, &
      &-9.69235274563376E-04_wp, 1.79530473016562E-03_wp, 1.18044156805790E-04_wp, &
      & 2.24225330283302E-05_wp, 9.58261702773522E-05_wp,-6.29010792725501E-06_wp, &
      &-7.34189784171004E-05_wp, 9.04158960138589E-04_wp, 6.54572917365500E-04_wp, &
      &-2.66475310005962E-05_wp,-6.86166560341841E-04_wp,-2.26628110062760E-03_wp, &
      & 3.58655742369494E-04_wp,-1.86189175940878E-03_wp,-1.69176114470135E-04_wp, &
      &-3.58853556247105E-04_wp, 1.24003860524734E-03_wp,-1.05936663656557E-03_wp, &
      &-1.27713933502452E-03_wp,-6.03449935835499E-04_wp, 1.02522733289301E-03_wp, &
      &-2.00224744625050E-04_wp, 1.76952445247401E-03_wp, 8.28947276525380E-04_wp, &
      & 2.30563240863010E-03_wp,-9.90104964157021E-05_wp,-4.98434975453687E-06_wp, &
      & 6.86374428440874E-05_wp, 5.60579840864022E-05_wp,-6.27090831700251E-05_wp, &
      & 4.58191515709555E-06_wp,-9.53916547354364E-05_wp,-2.66831117200006E-03_wp, &
      &-1.36859171542851E-02_wp,-2.27460176280414E-02_wp,-5.32740463510600E-04_wp, &
      &-2.50223062018044E-04_wp, 8.06895572993571E-04_wp, 2.13209637497826E-04_wp, &
      &-4.34071069253727E-04_wp,-3.77193211186177E-04_wp,-1.92654296930593E-03_wp, &
      & 3.14504756099470E-03_wp, 3.64690872990471E-03_wp, 2.20085105687747E-03_wp, &
      &-1.22517267268654E-04_wp,-1.16285810102388E-05_wp,-6.52491764156857E-05_wp, &
      &-2.32417372341244E-04_wp,-1.56685832977589E-04_wp,-5.82863081578004E-04_wp, &
      &-9.79001513299533E-05_wp,-5.84081941461458E-04_wp, 6.16125477932459E-05_wp, &
      & 1.02410080994947E-04_wp,-2.02294049636478E-03_wp,-1.22978136847514E-03_wp, &
      & 1.46303965116822E-03_wp,-3.58903980971413E-04_wp, 2.73469045608191E-04_wp, &
      & 1.27335881729851E-03_wp,-1.02227528479917E-03_wp, 9.59524078160904E-04_wp, &
      & 4.82528207292135E-04_wp, 5.23077500267552E-04_wp,-2.92072966954730E-04_wp, &
      & 8.36491920736309E-04_wp, 4.32903388983687E-04_wp,-1.69338120880881E-03_wp, &
      & 2.40242802702843E-03_wp, 1.12832721024584E-04_wp,-6.75347477155226E-05_wp, &
      &-9.95063274235428E-05_wp, 8.56172572143386E-06_wp, 4.56136607798075E-05_wp, &
      & 1.45049683061849E-03_wp, 1.27380822222116E-02_wp,-1.60912589585535E-02_wp, &
      &-1.00490947362311E-02_wp, 5.11478452729843E-04_wp, 2.13209637497826E-04_wp, &
      & 5.86625165555958E-04_wp,-1.45137738036140E-04_wp, 1.20314703232990E-04_wp, &
      &-1.26165346840580E-03_wp, 3.43704475815414E-04_wp, 4.41151041086971E-06_wp, &
      & 2.16340104418119E-03_wp, 1.38249008810169E-05_wp,-1.14143326968819E-05_wp, &
      & 6.32898306173221E-05_wp,-6.57957605354162E-05_wp,-1.36875169534557E-04_wp, &
      &-5.86305027750962E-04_wp,-7.38113854803026E-04_wp, 1.60898571495431E-05_wp, &
      & 2.02646773182076E-04_wp, 1.25422685863442E-03_wp,-3.48741562764879E-04_wp, &
      & 1.62253133956297E-03_wp,-6.46390044540534E-04_wp, 1.79457627302079E-03_wp, &
      & 2.83112214438342E-04_wp, 3.48810273263677E-04_wp, 9.54381580813492E-04_wp, &
      & 6.67148018086848E-04_wp,-5.52940806333149E-04_wp, 8.65873402966272E-04_wp, &
      &-1.42087387638032E-03_wp, 3.56449278824868E-04_wp,-1.03486087879590E-03_wp, &
      & 8.97862798560837E-04_wp, 2.48340309494439E-03_wp, 8.04650240365074E-05_wp, &
      & 1.75895781302209E-05_wp, 1.69168895582038E-05_wp, 8.84178196668619E-05_wp, &
      & 1.19018583406337E-04_wp, 1.78512184150084E-03_wp,-3.28692789808616E-05_wp, &
      & 1.70166983094865E-02_wp,-1.37415971532467E-02_wp, 3.56651216258034E-04_wp, &
      &-4.34071069253727E-04_wp,-1.45137738036140E-04_wp, 5.40221204066931E-04_wp, &
      &-2.59624443494387E-04_wp,-4.07412132404515E-05_wp,-2.31771909002759E-03_wp, &
      &-6.97059294014967E-04_wp,-1.10499621300400E-03_wp, 1.02396762885453E-04_wp, &
      & 4.30213740407625E-05_wp, 1.05577544367327E-05_wp, 5.27478245342355E-05_wp, &
      & 9.04761248419200E-05_wp,-5.20787359796830E-04_wp, 1.66609406392804E-04_wp, &
      &-1.95975356661681E-04_wp, 6.65298187726416E-04_wp, 1.23413770534980E-03_wp, &
      &-5.80912494235757E-04_wp,-2.57580865019132E-04_wp, 1.68271220553567E-03_wp, &
      &-2.19066150559248E-03_wp,-2.02065478669260E-03_wp, 1.31361034587710E-03_wp, &
      &-4.76500881652060E-05_wp, 1.39319355064774E-04_wp,-4.61445881759826E-04_wp, &
      &-8.10924344558873E-04_wp,-3.57640518895915E-04_wp,-2.26774951942169E-03_wp, &
      &-9.31983071527857E-04_wp,-2.17183687326237E-03_wp,-2.72169849948320E-03_wp, &
      &-1.60848310539632E-04_wp,-1.27743776572833E-04_wp, 1.01932389773469E-05_wp, &
      &-1.38692150236970E-04_wp,-4.22196797139872E-05_wp, 1.32158515813909E-03_wp, &
      & 2.27063683650581E-02_wp, 3.01226291286036E-04_wp, 1.76756468802079E-02_wp, &
      &-2.01067903244126E-05_wp,-3.77193211186177E-04_wp, 1.20314703232990E-04_wp, &
      &-2.59624443494387E-04_wp, 9.33862836396317E-04_wp, 1.33668771560144E-03_wp, &
      &-3.22448173345013E-04_wp,-4.38999821387560E-03_wp,-8.07257441219475E-04_wp, &
      & 5.60364474238915E-06_wp,-6.18831095843887E-05_wp, 9.70229104584636E-05_wp, &
      & 2.32508777451199E-04_wp, 3.39849143563926E-05_wp,-1.20082444862482E-03_wp, &
      & 2.55368882386110E-03_wp,-3.62978868139369E-03_wp,-9.07158739473559E-03_wp, &
      &-1.35647375904467E-02_wp,-5.21088322546849E-03_wp,-9.88315017273462E-03_wp, &
      & 2.13613158631452E-02_wp, 2.19877634531918E-02_wp,-8.80007601248442E-04_wp, &
      & 5.34558322462221E-05_wp, 5.47852591736715E-04_wp,-8.42490345787972E-04_wp, &
      &-4.40961160825250E-04_wp, 1.04269735966288E-03_wp, 2.43800444355134E-04_wp, &
      &-4.67636429378859E-04_wp,-1.75697167675842E-02_wp, 1.74301293583753E-03_wp, &
      & 1.34726136560787E-02_wp, 9.23735830986981E-04_wp,-3.13016321904110E-04_wp, &
      & 6.86048755790714E-04_wp, 9.28641711068648E-04_wp, 2.15520665339914E-03_wp, &
      &-3.02443053170287E-03_wp, 9.65825016803663E-03_wp, 1.87142472813712E-02_wp, &
      & 1.35400690416689E-02_wp,-1.11114144302852E-03_wp,-1.92654296930593E-03_wp, &
      &-1.26165346840580E-03_wp,-4.07412132404515E-05_wp, 1.33668771560144E-03_wp, &
      & 9.80347238799697E-01_wp,-2.18568952276218E-02_wp, 3.31155785434019E-02_wp, &
      & 7.98627529129182E-02_wp, 1.97136489453409E-03_wp, 7.92990830610551E-04_wp, &
      & 1.47395028414642E-03_wp,-3.09029625894063E-03_wp,-3.26271359901284E-03_wp, &
      & 1.01678825709110E-02_wp, 1.44149977964325E-02_wp, 5.13925269875120E-03_wp, &
      & 1.33594259106603E-02_wp, 7.56592855070334E-02_wp, 3.66481612496406E-02_wp, &
      &-6.24835765122572E-03_wp, 3.25299720293667E-02_wp, 7.20879088530275E-02_wp, &
      & 3.13091576817181E-02_wp, 2.45085826562595E-02_wp,-6.37022248510328E-03_wp, &
      &-1.91489694891317E-02_wp, 1.09352791492342E-03_wp, 3.42456107696032E-02_wp, &
      & 4.74547939918725E-02_wp,-2.15271159950398E-02_wp,-4.80545714771638E-02_wp, &
      &-1.23079993061857E-02_wp, 1.29098424984618E-02_wp,-1.64728596302177E-03_wp, &
      &-3.56800269797281E-03_wp, 1.12236019091068E-03_wp, 1.20019581993013E-04_wp, &
      & 5.25039772955339E-03_wp, 1.59234232521478E-02_wp,-5.19049786298460E-03_wp, &
      &-2.32116276965120E-02_wp, 8.14153009473268E-03_wp,-1.77564006933414E-03_wp, &
      & 3.14504756099470E-03_wp, 3.43704475815414E-04_wp,-2.31771909002759E-03_wp, &
      &-3.22448173345013E-04_wp,-2.18568952276218E-02_wp, 9.14808367283382E-01_wp, &
      & 7.48068598583395E-03_wp, 5.09970306011787E-03_wp,-2.71571345830003E-02_wp, &
      &-1.35284870829017E-02_wp,-5.25812007711735E-03_wp,-4.54834981829022E-04_wp, &
      &-8.92063226570178E-03_wp,-1.38192687806438E-02_wp, 4.82113496013389E-03_wp, &
      & 1.02205811651839E-02_wp,-1.99274174073717E-02_wp, 3.15195470481471E-02_wp, &
      &-1.62019398806568E-02_wp, 8.60672855726232E-02_wp, 4.67414968162166E-02_wp, &
      &-7.10369517108405E-02_wp, 3.67758744262592E-02_wp,-2.56042002348986E-03_wp, &
      &-1.46041445381265E-02_wp, 5.62405201973460E-03_wp, 5.78051213871221E-02_wp, &
      &-2.45166757199998E-03_wp,-2.86052824976296E-02_wp,-5.04314540727691E-03_wp, &
      &-4.04120386990838E-03_wp, 1.53478202076274E-02_wp,-9.80559410487820E-03_wp, &
      &-1.88915026257937E-03_wp, 2.99814266890537E-04_wp, 1.84140089469675E-03_wp, &
      & 5.49071489883944E-06_wp, 1.28691457765932E-05_wp, 1.87596405543643E-02_wp, &
      &-3.84519947027201E-02_wp,-2.18422956719651E-02_wp,-1.12403003698801E-02_wp, &
      &-9.69235274563376E-04_wp, 3.64690872990471E-03_wp, 4.41151041086971E-06_wp, &
      &-6.97059294014967E-04_wp,-4.38999821387560E-03_wp, 3.31155785434019E-02_wp, &
      & 7.48068598583395E-03_wp, 9.07039961834103E-01_wp,-6.65625977227704E-03_wp, &
      &-4.68411115220014E-04_wp, 8.51223110027637E-03_wp,-1.56557287161269E-02_wp, &
      &-2.67310457965848E-02_wp,-5.13497748132490E-05_wp,-4.08823673808927E-02_wp, &
      & 1.15441397373419E-02_wp,-1.81593535427986E-02_wp,-2.34198210219147E-02_wp, &
      & 6.45164556377193E-02_wp, 4.22691225405726E-02_wp, 3.14814661832583E-03_wp, &
      &-1.03864190087322E-01_wp,-4.72439825881444E-02_wp, 5.10602001473590E-02_wp, &
      & 4.06759272671501E-02_wp,-4.95678503503976E-03_wp, 4.94699911982638E-02_wp, &
      &-1.68963645227524E-03_wp,-6.15903243688644E-02_wp,-2.74882135918198E-02_wp, &
      & 2.25472983416135E-03_wp, 2.91141295701072E-02_wp,-5.76126756443682E-03_wp, &
      &-5.05472700198872E-03_wp,-5.38651668362386E-04_wp,-6.76380286992065E-04_wp, &
      &-1.05661020220776E-03_wp,-2.03354621399353E-03_wp,-2.17948090039576E-03_wp, &
      & 4.69726399080004E-03_wp,-1.24584223110816E-02_wp,-3.26436278753523E-02_wp, &
      &-3.83955909200850E-03_wp, 1.79530473016562E-03_wp, 2.20085105687747E-03_wp, &
      & 2.16340104418119E-03_wp,-1.10499621300400E-03_wp,-8.07257441219475E-04_wp, &
      & 7.98627529129182E-02_wp, 5.09970306011787E-03_wp,-6.65625977227704E-03_wp, &
      & 8.85714404731904E-01_wp, 8.42605643174469E-03_wp, 4.42603414028923E-05_wp, &
      & 1.51188725605635E-02_wp,-1.31893721956638E-02_wp,-2.63691966696696E-02_wp, &
      &-7.15273229046639E-04_wp,-3.91702493559700E-04_wp,-3.90294960024259E-04_wp, &
      &-6.68613238440462E-04_wp,-9.39640106534622E-04_wp,-7.23873827592798E-04_wp, &
      & 3.77982171071824E-04_wp,-1.50087725727626E-03_wp,-2.04362393783481E-03_wp, &
      &-6.95430188414474E-04_wp,-4.26658307535336E-04_wp, 1.86282736000768E-04_wp, &
      & 1.30383854676361E-03_wp,-1.50061542669221E-04_wp,-1.62313671252700E-03_wp, &
      &-1.80961917041824E-03_wp, 2.05817627414404E-03_wp, 5.40940849482020E-03_wp, &
      & 6.00834384066742E-04_wp,-2.40189133071662E-03_wp, 6.00370689453202E-05_wp, &
      & 1.50521404520208E-04_wp,-7.39580793665369E-05_wp,-4.17039948838420E-05_wp, &
      &-2.95522304533291E-04_wp,-1.69035008417164E-03_wp, 1.38574000819706E-03_wp, &
      & 1.96946150443605E-03_wp,-1.74863768724042E-03_wp, 1.18044156805790E-04_wp, &
      &-1.22517267268654E-04_wp, 1.38249008810169E-05_wp, 1.02396762885453E-04_wp, &
      & 5.60364474238915E-06_wp, 1.97136489453409E-03_wp,-2.71571345830003E-02_wp, &
      &-4.68411115220014E-04_wp, 8.42605643174469E-03_wp, 9.22078584954116E-04_wp, &
      & 4.11819260580167E-04_wp, 3.19311018631667E-04_wp,-1.03752914123972E-04_wp, &
      & 1.06588815353226E-05_wp,-2.62747591641126E-04_wp,-1.85560136617751E-04_wp, &
      & 9.37999421072550E-06_wp,-4.37632419523053E-04_wp,-6.34568272551357E-04_wp, &
      &-6.94871436601380E-04_wp, 9.34039516698034E-04_wp, 3.62419676815490E-05_wp, &
      &-1.49241291128878E-03_wp,-1.53919367966631E-04_wp,-4.73464627516159E-04_wp, &
      &-2.24409080145977E-05_wp, 4.09728762394107E-04_wp, 5.19503852174865E-04_wp, &
      &-5.52953181551917E-04_wp,-9.88501127912692E-04_wp, 5.05903385749172E-04_wp, &
      & 2.05121472162809E-03_wp, 6.74812958013988E-04_wp,-1.16055916082677E-03_wp, &
      & 1.06799866488623E-05_wp, 8.35092699986451E-05_wp,-5.70885224633915E-06_wp, &
      &-5.12570162749364E-06_wp,-1.25578343830325E-04_wp,-6.62367944133660E-04_wp, &
      &-6.85719533754288E-04_wp, 5.26688252429663E-04_wp,-9.80560541161186E-04_wp, &
      & 2.24225330283302E-05_wp,-1.16285810102388E-05_wp,-1.14143326968819E-05_wp, &
      & 4.30213740407625E-05_wp,-6.18831095843887E-05_wp, 7.92990830610551E-04_wp, &
      &-1.35284870829017E-02_wp, 8.51223110027637E-03_wp, 4.42603414028923E-05_wp, &
      & 4.11819260580167E-04_wp, 2.86722146733572E-04_wp,-6.45163188426943E-05_wp, &
      &-2.51940997896695E-04_wp, 1.26351255067754E-04_wp,-5.93120966859311E-04_wp, &
      & 2.73984936091776E-05_wp,-5.75461842467355E-04_wp,-1.19585665770434E-04_wp, &
      & 3.16393662007085E-04_wp, 6.92791431294828E-04_wp,-1.19873916027241E-03_wp, &
      &-2.39989501736821E-03_wp, 1.88746337135601E-04_wp,-6.53034622285015E-05_wp, &
      & 6.75174714704217E-04_wp, 1.12001295917723E-04_wp, 9.57203891809881E-04_wp, &
      &-1.09043230139903E-03_wp,-1.20152015219740E-03_wp,-3.96125716686334E-04_wp, &
      & 7.91175763385960E-04_wp, 2.29185700862415E-03_wp,-1.66126666824316E-03_wp, &
      &-5.81589521031109E-04_wp, 4.33426688768389E-05_wp,-1.78734093583186E-05_wp, &
      &-9.38004629219951E-05_wp,-6.69133655858396E-05_wp,-1.08474701764776E-04_wp, &
      &-1.25089500655171E-03_wp, 2.13681379765512E-03_wp,-1.74729268447544E-05_wp, &
      &-3.97380267628569E-04_wp, 9.58261702773522E-05_wp,-6.52491764156857E-05_wp, &
      & 6.32898306173221E-05_wp, 1.05577544367327E-05_wp, 9.70229104584636E-05_wp, &
      & 1.47395028414642E-03_wp,-5.25812007711735E-03_wp,-1.56557287161269E-02_wp, &
      & 1.51188725605635E-02_wp, 3.19311018631667E-04_wp,-6.45163188426943E-05_wp, &
      & 5.64230404713589E-04_wp, 2.39346126111102E-04_wp,-4.00869731963038E-04_wp, &
      & 1.02934074705867E-03_wp,-3.83899683641624E-04_wp,-4.74477134909322E-05_wp, &
      & 1.01952206317582E-03_wp,-1.37698019287476E-03_wp, 2.14712336521353E-04_wp, &
      &-2.29602069996173E-03_wp, 3.71610548046474E-04_wp, 2.32194941901059E-03_wp, &
      &-2.04327901205466E-03_wp,-5.61493981444976E-04_wp, 5.82371418669850E-04_wp, &
      &-7.51198337433286E-04_wp,-1.81994295729175E-03_wp, 1.02978296745931E-03_wp, &
      & 1.22747256066821E-03_wp, 7.16658955566675E-04_wp, 4.08333274592717E-04_wp, &
      &-2.78697676847693E-04_wp, 2.45136518477340E-03_wp, 1.15971638594665E-04_wp, &
      & 5.07786354518432E-06_wp,-6.43822060302723E-05_wp, 5.12877509169561E-05_wp, &
      & 5.08840121951604E-05_wp,-1.79892475752792E-03_wp, 3.54363522635215E-03_wp, &
      & 3.65712487530910E-03_wp, 2.70852057675908E-03_wp,-6.29010792725501E-06_wp, &
      &-2.32417372341244E-04_wp,-6.57957605354162E-05_wp, 5.27478245342355E-05_wp, &
      & 2.32508777451199E-04_wp,-3.09029625894063E-03_wp,-4.54834981829022E-04_wp, &
      &-2.67310457965848E-02_wp,-1.31893721956638E-02_wp,-1.03752914123972E-04_wp, &
      &-2.51940997896695E-04_wp, 2.39346126111102E-04_wp, 1.01787209031264E-03_wp, &
      & 4.18069424036841E-04_wp, 1.20857782799951E-03_wp,-5.51201130588527E-04_wp, &
      & 5.45960645253513E-04_wp, 5.93888771707889E-04_wp,-2.23650399659197E-03_wp, &
      &-1.19295873522822E-03_wp, 4.64419548681912E-05_wp, 2.60086011058214E-03_wp, &
      & 4.24534276701425E-04_wp,-1.86505173435443E-03_wp,-1.60711005567859E-03_wp, &
      & 4.18386118275000E-04_wp,-1.22239364856491E-03_wp,-3.20961296343504E-05_wp, &
      & 1.53238668144887E-03_wp, 4.59125102539591E-04_wp, 1.18813026352344E-03_wp, &
      &-4.29139345849666E-04_wp, 2.05969982294557E-03_wp, 1.40102242660422E-03_wp, &
      & 5.65401762358019E-05_wp, 9.85537457449779E-05_wp, 4.02430151216053E-05_wp, &
      & 1.03773290535138E-04_wp, 3.72930282297249E-05_wp,-3.50141084196518E-04_wp, &
      &-6.38061625562823E-06_wp, 4.03104249035507E-03_wp, 1.45358953653843E-03_wp, &
      &-7.34189784171004E-05_wp,-1.56685832977589E-04_wp,-1.36875169534557E-04_wp, &
      & 9.04761248419200E-05_wp, 3.39849143563926E-05_wp,-3.26271359901284E-03_wp, &
      &-8.92063226570178E-03_wp,-5.13497748132490E-05_wp,-2.63691966696696E-02_wp, &
      & 1.06588815353226E-05_wp, 1.26351255067754E-04_wp,-4.00869731963038E-04_wp, &
      & 4.18069424036841E-04_wp, 8.89205957396783E-04_wp],&
      & shape(density))
   
   call get_structure(mol, "f-block", "CeCl3")
   call test_e_gen(error, mol, ref, thr_in=thr1*10, guess_qat=qat, guess_qsh=qsh, &
      & guess_density=density, guess_dp=dp, guess_qp=qp)

end subroutine test_e_cecl3

subroutine test_e_ce2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol
   ! Reference energy calculated with g-xTB development code
   real(wp), parameter :: ref = -949.879382362896_wp 

   ! Fully converged values from the g-xTB development code
   real(wp), parameter :: qat(2, 1) = reshape([&
      & 4.68505789719131E-12_wp,-4.68625138694279E-12_wp], &
      & shape(qat))

   real(wp), parameter :: qsh(8, 1) = reshape([&
      &-8.71550788639256E-02_wp,-7.18167264352547E-02_wp, 4.08861975655191E-01_wp, &
      &-2.49890170351326E-01_wp,-8.71550788674514E-02_wp,-7.18167264347737E-02_wp, &
      & 4.08861975654052E-01_wp,-2.49890170356514E-01_wp], shape(qsh))

   real(wp), parameter :: dp(3, 2, 1) = reshape([&
      & 1.34724669685084E-16_wp,-1.30204626767757E-15_wp,-4.18312171738696E+00_wp, &
      & 1.76015419373933E-16_wp,-1.11705729020889E-15_wp, 4.18312171738804E+00_wp],&
      & shape(dp))

   real(wp), parameter :: qp(6, 2, 1) = reshape([&
      & 1.26126515740945E+00_wp, 1.32068945012275E-15_wp, 1.26126515740945E+00_wp, &
      & 2.95233242730535E-16_wp,-1.68397021221537E-15_wp,-2.52253031481891E+00_wp, &
      & 1.26126515741057E+00_wp, 1.34498779587772E-15_wp, 1.26126515741055E+00_wp, &
      &-4.77906567964199E-16_wp, 1.23213993471570E-15_wp,-2.52253031482113E+00_wp],&
      & shape(qp))
      
   real(wp), parameter :: density(32, 32, 1) = reshape([&
      & 6.36026010729421E-01_wp, 2.61536441758601E-17_wp, 6.78407335711013E-03_wp, &
      & 4.26676263633402E-18_wp,-1.98302631888483E-17_wp, 2.84822020488551E-17_wp, &
      & 5.16421185242582E-03_wp, 1.14475139662911E-18_wp,-3.83723305640967E-17_wp, &
      &-8.20891177498193E-18_wp, 5.67656570173949E-17_wp,-7.14740715342786E-17_wp, &
      &-1.03520371538161E-02_wp, 2.79545753689680E-18_wp,-2.09340112174854E-17_wp, &
      &-1.51182987629402E-29_wp, 6.36026010731052E-01_wp, 1.96757264581643E-17_wp, &
      &-6.78407335661569E-03_wp,-6.80418305457184E-18_wp, 1.49244852314233E-17_wp, &
      &-7.18530087839911E-17_wp, 5.16421185232193E-03_wp,-1.03972787565093E-17_wp, &
      &-1.03034746724869E-17_wp, 4.30343990420659E-17_wp,-2.47513382605561E-17_wp, &
      & 7.29241956472982E-17_wp, 1.03520371538272E-02_wp, 9.50478865760614E-18_wp, &
      & 3.12096455677403E-17_wp,-3.39346874026659E-31_wp, 2.61536441758601E-17_wp, &
      & 4.87261748321654E-02_wp, 1.42965026117582E-17_wp,-9.55441696196146E-17_wp, &
      & 4.67352574887384E-18_wp, 1.18968212390530E-01_wp,-1.70728188200136E-17_wp, &
      &-1.04083408558608E-16_wp, 2.78544624815793E-17_wp,-5.57954803532025E-16_wp, &
      &-2.53615849041892E-17_wp, 1.45926947360124E-01_wp,-2.67150392268047E-17_wp, &
      &-1.14430766822974E-16_wp, 7.34231138185011E-17_wp, 5.13311415360163E-20_wp, &
      &-9.53957447531848E-18_wp, 4.87261748321289E-02_wp,-1.85257725541819E-17_wp, &
      &-7.20493261191624E-17_wp, 2.67964560345503E-18_wp,-1.18968212390598E-01_wp, &
      & 1.56905261145317E-17_wp, 1.00613961606655E-16_wp, 4.53843954983374E-17_wp, &
      &-5.49434910072695E-16_wp, 1.78524787628581E-17_wp, 1.45926947360386E-01_wp, &
      & 9.66556684429014E-17_wp,-1.29963933509606E-16_wp,-4.68595854827922E-17_wp, &
      & 5.18651720592585E-20_wp, 6.78407335711013E-03_wp, 1.42965026117582E-17_wp, &
      & 7.23612722408470E-05_wp,-3.71617888654111E-19_wp, 4.33749398783186E-18_wp, &
      & 3.45285588350078E-17_wp, 5.50832693121079E-05_wp,-1.00623717339668E-18_wp, &
      & 8.46885794700693E-17_wp,-8.75590914902289E-20_wp, 5.77969127155390E-18_wp, &
      & 4.12178755635086E-17_wp,-1.10418407835988E-04_wp,-1.21941497761440E-18_wp, &
      & 9.65701108795419E-17_wp,-1.65647414197021E-31_wp, 6.78407335712753E-03_wp, &
      & 1.42274069014201E-17_wp,-7.23612722355731E-05_wp,-4.89704423532376E-19_wp, &
      & 4.70820015206975E-18_wp,-3.49911668565306E-17_wp, 5.50832693109998E-05_wp, &
      & 9.07546517214644E-19_wp, 8.49879715927887E-17_wp, 4.59019780725063E-19_wp, &
      &-5.43821547771220E-18_wp, 4.27580774117376E-17_wp, 1.10418407836107E-04_wp, &
      &-1.14785093191382E-18_wp,-9.64605074106201E-17_wp,-7.67247719260822E-33_wp, &
      & 4.26676263633402E-18_wp,-9.55441696196146E-17_wp,-3.71617888654111E-19_wp, &
      & 4.87261748321650E-02_wp,-4.55327150264907E-18_wp,-1.42579377556993E-16_wp, &
      & 2.44172800744103E-18_wp, 1.18968212390530E-01_wp,-1.55717661130493E-17_wp, &
      &-1.87855166361028E-19_wp,-7.46485666884351E-18_wp,-1.13136325496879E-16_wp, &
      & 3.74900858695076E-18_wp, 1.45926947360124E-01_wp,-2.25688448939951E-17_wp, &
      & 5.14646449135489E-16_wp, 1.09284561265826E-18_wp,-9.21938897339646E-17_wp, &
      &-1.29309845068585E-18_wp, 4.87261748321285E-02_wp,-7.57809024724549E-18_wp, &
      & 1.32914605131794E-16_wp, 6.09730216460253E-18_wp,-1.18968212390597E-01_wp, &
      &-1.56358560639567E-17_wp, 3.76988216285415E-18_wp, 1.20697971714131E-17_wp, &
      &-1.55698816381950E-16_wp,-2.32694314675766E-17_wp, 1.45926947360385E-01_wp, &
      & 2.77965952914473E-17_wp, 4.70175930933210E-16_wp,-1.98302631888483E-17_wp, &
      & 4.67352574887384E-18_wp, 4.33749398783186E-18_wp,-4.55327150264907E-18_wp, &
      & 1.07609530465231E-03_wp, 8.31262442978668E-18_wp,-1.42049417564082E-17_wp, &
      &-1.03609087795711E-17_wp, 2.01304045972900E-02_wp, 7.16533147663587E-19_wp, &
      & 1.22398975552022E-03_wp, 2.30155589166069E-18_wp,-3.03489343211416E-17_wp, &
      &-6.57092612771501E-18_wp, 2.28970509350203E-02_wp, 1.27508104680788E-22_wp, &
      &-2.15546036602755E-17_wp, 5.07766499624638E-18_wp,-4.86572067851129E-18_wp, &
      &-4.25571158609913E-18_wp, 1.07609530465305E-03_wp,-3.22699707633591E-18_wp, &
      &-1.18327968399213E-17_wp, 9.46906247451355E-18_wp, 2.01304045973042E-02_wp, &
      & 1.14644400634831E-17_wp,-1.22398975552264E-03_wp, 7.13010004780492E-18_wp, &
      & 2.97898986809078E-17_wp,-1.06737035373510E-17_wp,-2.28970509350658E-02_wp, &
      &-1.01213150778144E-20_wp, 2.84822020488551E-17_wp, 1.18968212390530E-01_wp, &
      & 3.45285588350078E-17_wp,-1.42579377556993E-16_wp, 8.31262442978668E-18_wp, &
      & 2.90468841606159E-01_wp,-4.19716446081246E-17_wp,-3.06301703403424E-17_wp, &
      & 1.00526731779722E-17_wp,-1.36228394286136E-15_wp,-6.54458968166054E-17_wp, &
      & 3.56290394779373E-01_wp,-6.46508075802931E-17_wp,-1.78663047920054E-17_wp, &
      & 1.13346347882525E-16_wp, 1.25328412697737E-19_wp,-5.86651767886133E-17_wp, &
      & 1.18968212390441E-01_wp,-4.48546040561146E-17_wp,-8.72661729585440E-17_wp, &
      & 3.44443279687886E-18_wp,-2.90468841606325E-01_wp, 3.80222517365649E-17_wp, &
      & 2.77555756156289E-17_wp, 5.28531762240424E-17_wp,-1.34148205356663E-15_wp, &
      & 4.71119123131707E-17_wp, 3.56290394780011E-01_wp, 2.35415528184118E-16_wp, &
      &-5.08621862915076E-17_wp,-4.84897181188293E-17_wp, 1.26632283911305E-19_wp, &
      & 5.16421185242582E-03_wp,-1.70728188200136E-17_wp, 5.50832693121079E-05_wp, &
      & 2.44172800744103E-18_wp,-1.42049417564082E-17_wp,-4.19716446081246E-17_wp, &
      & 4.19308072419084E-05_wp, 5.88635135711540E-18_wp,-2.63029923050218E-16_wp, &
      &-6.66522418401444E-20_wp,-1.55131657306094E-17_wp,-5.23466086776269E-17_wp, &
      &-8.40533438328713E-05_wp, 7.23152185588851E-18_wp,-2.98995345915026E-16_wp, &
      &-9.73492174866038E-32_wp, 5.16421185243906E-03_wp,-1.71254162592444E-17_wp, &
      &-5.50832693080934E-05_wp, 2.35183748725598E-18_wp,-1.39227506499343E-17_wp, &
      & 4.16194954343530E-17_wp, 4.19308072410649E-05_wp,-5.96147723349914E-18_wp, &
      &-2.62802018034202E-16_wp, 3.49417712241383E-19_wp, 1.57731059242531E-17_wp, &
      &-5.11741672519457E-17_wp, 8.40533438329617E-05_wp, 7.28599825578279E-18_wp, &
      & 2.99078778907051E-16_wp, 2.05851637011369E-32_wp, 1.14475139662911E-18_wp, &
      &-1.04083408558608E-16_wp,-1.00623717339668E-18_wp, 1.18968212390530E-01_wp, &
      &-1.03609087795711E-17_wp,-3.06301703403424E-17_wp, 5.88635135711540E-18_wp, &
      & 2.90468841606159E-01_wp,-2.38732038902015E-17_wp,-4.58660943678613E-19_wp, &
      &-1.73658082122632E-17_wp, 1.19492227582965E-16_wp, 9.30438081675999E-18_wp, &
      & 3.56290394779374E-01_wp,-3.90128251179812E-17_wp, 1.25654370115600E-15_wp, &
      &-6.60457914048679E-18_wp,-9.45951840004271E-17_wp,-3.05827902782579E-18_wp, &
      & 1.18968212390441E-01_wp,-1.77462059035397E-17_wp, 1.07616334190651E-17_wp, &
      & 1.48116799764167E-17_wp,-2.90468841606325E-01_wp,-2.40296837886836E-17_wp, &
      & 9.20441925479386E-18_wp, 2.86090782132920E-17_wp, 3.16979027818925E-19_wp, &
      &-5.69647978967860E-17_wp, 3.56290394780011E-01_wp, 5.17767271018625E-17_wp, &
      & 1.14796595884750E-15_wp,-3.83723305640967E-17_wp, 2.78544624815793E-17_wp, &
      & 8.46885794700693E-17_wp,-1.55717661130493E-17_wp, 2.01304045972900E-02_wp, &
      & 1.00526731779722E-17_wp,-2.63029923050218E-16_wp,-2.38732038902015E-17_wp, &
      & 3.76577416051013E-01_wp, 1.34041121706213E-17_wp, 2.28970509350203E-02_wp, &
      &-1.35355565563586E-16_wp,-5.73147686070294E-16_wp, 8.55364297958466E-17_wp, &
      & 4.28332785594342E-01_wp, 2.38528104931958E-21_wp,-7.06293910700988E-17_wp, &
      & 3.54146539703156E-17_wp,-9.45700620452786E-17_wp,-1.00053432449881E-17_wp, &
      & 2.01304045973039E-02_wp, 8.50836366019283E-17_wp,-2.18654451169821E-16_wp, &
      & 7.18952644834967E-18_wp, 3.76577416051279E-01_wp, 2.14464105513279E-16_wp, &
      &-2.28970509350655E-02_wp,-4.50284848129400E-17_wp, 5.62689863632493E-16_wp, &
      & 8.78619320330032E-18_wp,-4.28332785595195E-01_wp,-1.89338404034914E-19_wp, &
      &-8.20891177498193E-18_wp,-5.57954803532025E-16_wp,-8.75590914902289E-20_wp, &
      &-1.87855166361028E-19_wp, 7.16533147663587E-19_wp,-1.36228394286136E-15_wp, &
      &-6.66522418401444E-20_wp,-4.58660943678613E-19_wp, 1.34041121706213E-17_wp, &
      & 6.38962548470201E-30_wp, 8.15010741557540E-19_wp,-1.67098364533632E-15_wp, &
      & 1.33609252221827E-19_wp,-5.62595587841840E-19_wp, 1.52463224286485E-17_wp, &
      &-2.57182745703147E-33_wp,-8.20891177500258E-18_wp,-5.57954803531607E-16_wp, &
      & 8.75590914838955E-20_wp,-1.87855166361148E-19_wp, 7.16533147664106E-19_wp, &
      & 1.36228394286214E-15_wp,-6.66522418391771E-20_wp, 4.58660943678858E-19_wp, &
      & 1.34041121706306E-17_wp, 6.29854570369907E-30_wp,-8.15010741559064E-19_wp, &
      &-1.67098364533931E-15_wp,-1.33609252222772E-19_wp,-5.62595587842681E-19_wp, &
      &-1.52463224286792E-17_wp,-2.41331863526078E-33_wp, 5.67656570173949E-17_wp, &
      &-2.53615849041892E-17_wp, 5.77969127155390E-18_wp,-7.46485666884351E-18_wp, &
      & 1.22398975552022E-03_wp,-6.54458968166054E-17_wp,-1.55131657306094E-17_wp, &
      &-1.73658082122632E-17_wp, 2.28970509350203E-02_wp, 8.15010741557540E-19_wp, &
      & 1.39221025790323E-03_wp,-8.92560038611505E-17_wp,-3.58110200130173E-17_wp, &
      &-1.43196099637199E-17_wp, 2.60439346356453E-02_wp, 1.45032334194132E-22_wp, &
      & 5.48043297420948E-17_wp,-2.49019022938537E-17_wp,-6.38051543194414E-18_wp, &
      &-7.12640124771324E-18_wp, 1.22398975552107E-03_wp, 7.12304734745942E-17_wp, &
      &-1.28150022279692E-17_wp, 1.63513899389317E-17_wp, 2.28970509350364E-02_wp, &
      & 1.30400691554111E-17_wp,-1.39221025790598E-03_wp,-8.37638429387216E-17_wp, &
      & 3.51751526380207E-17_wp,-1.89862575151744E-17_wp,-2.60439346356971E-02_wp, &
      &-1.15123501739130E-20_wp,-7.14740715342786E-17_wp, 1.45926947360124E-01_wp, &
      & 4.12178755635086E-17_wp,-1.13136325496879E-16_wp, 2.30155589166069E-18_wp, &
      & 3.56290394779373E-01_wp,-5.23466086776269E-17_wp, 1.19492227582965E-16_wp, &
      &-1.35355565563586E-16_wp,-1.67098364533632E-15_wp,-8.92560038611505E-17_wp, &
      & 4.37027409584127E-01_wp,-7.75690231491895E-17_wp, 1.79448060873895E-16_wp, &
      &-2.89525100666792E-17_wp, 1.53728397822103E-19_wp,-1.78369440546588E-16_wp, &
      & 1.45926947360015E-01_wp,-5.38838493320487E-17_wp,-3.69653666937220E-17_wp, &
      &-3.66978998200006E-18_wp,-3.56290394779576E-01_wp, 4.57742663566695E-17_wp, &
      &-1.24900090270330E-16_wp,-8.28562749591696E-17_wp,-1.64546795384943E-15_wp, &
      & 6.67674555525528E-17_wp, 4.37027409584909E-01_wp, 2.87029801737001E-16_wp, &
      & 1.28365872778174E-16_wp, 1.08505945453719E-16_wp, 1.55327732149428E-19_wp, &
      &-1.03520371538161E-02_wp,-2.67150392268047E-17_wp,-1.10418407835988E-04_wp, &
      & 3.74900858695076E-18_wp,-3.03489343211416E-17_wp,-6.46508075802931E-17_wp, &
      &-8.40533438328713E-05_wp, 9.30438081675999E-18_wp,-5.73147686070294E-16_wp, &
      & 1.33609252221827E-19_wp,-3.58110200130173E-17_wp,-7.75690231491895E-17_wp, &
      & 1.68491023049652E-04_wp, 1.13901506157173E-17_wp,-6.52288591182026E-16_wp, &
      & 2.86366518214818E-31_wp,-1.03520371538426E-02_wp,-2.66096038455873E-17_wp, &
      & 1.10418407827941E-04_wp, 3.92920065328916E-18_wp,-3.09146068578970E-17_wp, &
      & 6.53567161361734E-17_wp,-8.40533438311804E-05_wp,-9.15378554962150E-18_wp, &
      &-5.73604538202771E-16_wp,-7.00433143079671E-19_wp, 3.52899510447828E-17_wp, &
      &-7.99192670097916E-17_wp,-1.68491023049833E-04_wp, 1.12809487216846E-17_wp, &
      & 6.52121343701652E-16_wp, 4.26293631624487E-32_wp, 2.79545753689680E-18_wp, &
      &-1.14430766822974E-16_wp,-1.21941497761440E-18_wp, 1.45926947360124E-01_wp, &
      &-6.57092612771501E-18_wp,-1.78663047920054E-17_wp, 7.23152185588851E-18_wp, &
      & 3.56290394779374E-01_wp, 8.55364297958466E-17_wp,-5.62595587841840E-19_wp, &
      &-1.43196099637199E-17_wp, 1.79448060873895E-16_wp, 1.13901506157173E-17_wp, &
      & 4.37027409584127E-01_wp, 8.27464454221775E-17_wp, 1.54128218664302E-15_wp, &
      &-6.70990634680478E-18_wp,-1.10795250815539E-16_wp,-3.76613897685730E-18_wp, &
      & 1.45926947360015E-01_wp,-1.56297649804834E-17_wp,-4.12090196051392E-18_wp, &
      & 1.81793703057341E-17_wp,-3.56290394779576E-01_wp, 8.53444908447583E-17_wp, &
      & 1.12901822855475E-17_wp, 2.81106557954278E-17_wp, 4.02277783928515E-17_wp, &
      &-6.98506337535823E-17_wp, 4.37027409584909E-01_wp,-6.70901855233864E-17_wp, &
      & 1.40810023687713E-15_wp,-2.09340112174854E-17_wp, 7.34231138185011E-17_wp, &
      & 9.65701108795419E-17_wp,-2.25688448939951E-17_wp, 2.28970509350203E-02_wp, &
      & 1.13346347882525E-16_wp,-2.98995345915026E-16_wp,-3.90128251179812E-17_wp, &
      & 4.28332785594342E-01_wp, 1.52463224286485E-17_wp, 2.60439346356453E-02_wp, &
      &-2.89525100666792E-17_wp,-6.52288591182026E-16_wp, 8.27464454221775E-17_wp, &
      & 4.87201216522647E-01_wp, 2.71310501563042E-21_wp,-5.76243595832713E-17_wp, &
      & 8.20223493062356E-17_wp,-1.07809666882713E-16_wp,-1.62373940123227E-17_wp, &
      & 2.28970509350361E-02_wp,-5.13486508439339E-18_wp,-2.48521077750157E-16_wp, &
      & 2.00362064254378E-17_wp, 4.28332785594644E-01_wp, 2.43939237482200E-16_wp, &
      &-2.60439346356967E-02_wp, 7.37887815336974E-17_wp, 6.40393485354896E-16_wp, &
      &-4.55205213367850E-18_wp,-4.87201216523616E-01_wp,-2.15360355038639E-19_wp, &
      &-1.51182987629402E-29_wp, 5.13311415360163E-20_wp,-1.65647414197021E-31_wp, &
      & 5.14646449135489E-16_wp, 1.27508104680788E-22_wp, 1.25328412697737E-19_wp, &
      &-9.73492174866038E-32_wp, 1.25654370115600E-15_wp, 2.38528104931958E-21_wp, &
      &-2.57182745703147E-33_wp, 1.45032334194132E-22_wp, 1.53728397822103E-19_wp, &
      & 2.86366518214818E-31_wp, 1.54128218664302E-15_wp, 2.71310501563042E-21_wp, &
      & 5.43570208734199E-30_wp,-1.51518595187965E-29_wp, 5.13311415359784E-20_wp, &
      & 1.48060142070126E-31_wp, 5.14646449135104E-16_wp, 1.27508104648926E-22_wp, &
      &-1.25328412697920E-19_wp,-5.87042029346142E-32_wp,-1.25654370115671E-15_wp, &
      & 2.38528104932060E-21_wp, 3.92400910952554E-32_wp,-1.45032334145789E-22_wp, &
      & 1.53728397821890E-19_wp,-4.92467847758371E-31_wp, 1.54128218664578E-15_wp, &
      &-2.71310501558058E-21_wp, 4.96600393445783E-30_wp, 6.36026010731052E-01_wp, &
      &-9.53957447531848E-18_wp, 6.78407335712753E-03_wp, 1.09284561265826E-18_wp, &
      &-2.15546036602755E-17_wp,-5.86651767886133E-17_wp, 5.16421185243906E-03_wp, &
      &-6.60457914048679E-18_wp,-7.06293910700988E-17_wp,-8.20891177500258E-18_wp, &
      & 5.48043297420948E-17_wp,-1.78369440546588E-16_wp,-1.03520371538426E-02_wp, &
      &-6.70990634680478E-18_wp,-5.76243595832713E-17_wp,-1.51518595187965E-29_wp, &
      & 6.36026010732683E-01_wp,-1.60174921930038E-17_wp,-6.78407335663309E-03_wp, &
      &-9.97810007827365E-18_wp, 1.32001447600840E-17_wp, 1.52943700534162E-17_wp, &
      & 5.16421185233517E-03_wp,-2.64794821941287E-18_wp,-4.25605351784397E-17_wp, &
      & 4.30343990421766E-17_wp,-2.27900109851700E-17_wp,-3.39711733648326E-17_wp, &
      & 1.03520371538538E-02_wp,-5.75226095149484E-22_wp, 6.78999939336255E-17_wp, &
      &-3.69994885462225E-31_wp, 1.96757264581643E-17_wp, 4.87261748321289E-02_wp, &
      & 1.42274069014201E-17_wp,-9.21938897339646E-17_wp, 5.07766499624638E-18_wp, &
      & 1.18968212390441E-01_wp,-1.71254162592444E-17_wp,-9.45951840004271E-17_wp, &
      & 3.54146539703156E-17_wp,-5.57954803531607E-16_wp,-2.49019022938537E-17_wp, &
      & 1.45926947360015E-01_wp,-2.66096038455873E-17_wp,-1.10795250815539E-16_wp, &
      & 8.20223493062356E-17_wp, 5.13311415359784E-20_wp,-1.60174921930038E-17_wp, &
      & 4.87261748320923E-02_wp,-1.84566768438457E-17_wp,-6.86990462394901E-17_wp, &
      & 3.08378485082934E-18_wp,-1.18968212390509E-01_wp, 1.56379286752774E-17_wp, &
      & 9.36750677027476E-17_wp, 5.29445869870659E-17_wp,-5.49434910072282E-16_wp, &
      & 1.73927961525273E-17_wp, 1.45926947360276E-01_wp, 9.65502330616313E-17_wp, &
      &-1.19389523555573E-16_wp,-5.54588209705638E-17_wp, 5.18651720592461E-20_wp, &
      &-6.78407335661569E-03_wp,-1.85257725541819E-17_wp,-7.23612722355731E-05_wp, &
      &-1.29309845068585E-18_wp,-4.86572067851129E-18_wp,-4.48546040561146E-17_wp, &
      &-5.50832693080934E-05_wp,-3.05827902782579E-18_wp,-9.45700620452786E-17_wp, &
      & 8.75590914838955E-20_wp,-6.38051543194414E-18_wp,-5.38838493320487E-17_wp, &
      & 1.10418407827941E-04_wp,-3.76613897685730E-18_wp,-1.07809666882713E-16_wp, &
      & 1.48060142070126E-31_wp,-6.78407335663309E-03_wp,-1.84566768438457E-17_wp, &
      & 7.23612722302993E-05_wp,-1.17501191581495E-18_wp,-5.23642684272253E-18_wp, &
      & 4.53172120776096E-17_wp,-5.50832693069852E-05_wp, 3.15696968400295E-18_wp, &
      &-9.48694541679830E-17_wp,-4.59019780691566E-19_wp, 6.03903963812852E-18_wp, &
      &-5.54240511801882E-17_wp,-1.10418407828060E-04_wp,-3.83770302256159E-18_wp, &
      & 1.07700063413821E-16_wp,-8.39048739596107E-33_wp,-6.80418305457184E-18_wp, &
      &-7.20493261191624E-17_wp,-4.89704423532376E-19_wp, 4.87261748321285E-02_wp, &
      &-4.25571158609913E-18_wp,-8.72661729585440E-17_wp, 2.35183748725598E-18_wp, &
      & 1.18968212390441E-01_wp,-1.00053432449881E-17_wp,-1.87855166361148E-19_wp, &
      &-7.12640124771324E-18_wp,-3.69653666937220E-17_wp, 3.92920065328916E-18_wp, &
      & 1.45926947360015E-01_wp,-1.62373940123227E-17_wp, 5.14646449135104E-16_wp, &
      &-9.97810007827365E-18_wp,-6.86990462394901E-17_wp,-1.17501191581495E-18_wp, &
      & 4.87261748320920E-02_wp,-7.28053033069308E-18_wp, 7.76014005222563E-17_wp, &
      & 6.00741164441656E-18_wp,-1.18968212390508E-01_wp,-1.00694331958914E-17_wp, &
      & 3.76988216285108E-18_wp, 1.17313417502787E-17_wp,-8.64667514399829E-17_wp, &
      &-2.34496235339006E-17_wp, 1.45926947360276E-01_wp, 2.14651444097584E-17_wp, &
      & 4.70175930932858E-16_wp, 1.49244852314233E-17_wp, 2.67964560345503E-18_wp, &
      & 4.70820015206975E-18_wp,-7.57809024724549E-18_wp, 1.07609530465305E-03_wp, &
      & 3.44443279687886E-18_wp,-1.39227506499343E-17_wp,-1.77462059035397E-17_wp, &
      & 2.01304045973039E-02_wp, 7.16533147664106E-19_wp, 1.22398975552107E-03_wp, &
      &-3.66978998200006E-18_wp,-3.09146068578970E-17_wp,-1.56297649804834E-17_wp, &
      & 2.28970509350361E-02_wp, 1.27508104648926E-22_wp, 1.32001447600840E-17_wp, &
      & 3.08378485082934E-18_wp,-5.23642684272253E-18_wp,-7.28053033069308E-18_wp, &
      & 1.07609530465380E-03_wp, 1.64119455657820E-18_wp,-1.15506057334516E-17_wp, &
      & 1.68543595984857E-17_wp, 2.01304045973181E-02_wp, 1.14644400634911E-17_wp, &
      &-1.22398975552348E-03_wp, 1.15875417413683E-18_wp, 3.03555712176634E-17_wp, &
      &-1.97325423901384E-17_wp,-2.28970509350817E-02_wp,-1.01213150778506E-20_wp, &
      &-7.18530087839911E-17_wp,-1.18968212390598E-01_wp,-3.49911668565306E-17_wp, &
      & 1.32914605131794E-16_wp,-3.22699707633591E-18_wp,-2.90468841606325E-01_wp, &
      & 4.16194954343530E-17_wp, 1.07616334190651E-17_wp, 8.50836366019283E-17_wp, &
      & 1.36228394286214E-15_wp, 7.12304734745942E-17_wp,-3.56290394779576E-01_wp, &
      & 6.53567161361734E-17_wp,-4.12090196051392E-18_wp,-5.13486508439339E-18_wp, &
      &-1.25328412697920E-19_wp, 1.52943700534162E-17_wp,-1.18968212390509E-01_wp, &
      & 4.53172120776096E-17_wp, 7.76014005222563E-17_wp, 1.64119455657820E-18_wp, &
      & 2.90468841606490E-01_wp,-3.83744009103750E-17_wp,-1.07651953126608E-17_wp, &
      & 4.22831335559007E-17_wp, 1.34148205356740E-15_wp,-5.28964889711605E-17_wp, &
      &-3.56290394780214E-01_wp,-2.36121436740097E-16_wp, 2.88749796181813E-17_wp, &
      &-5.97217646795548E-17_wp,-1.26632283911453E-19_wp, 5.16421185232193E-03_wp, &
      & 1.56905261145317E-17_wp, 5.50832693109998E-05_wp, 6.09730216460253E-18_wp, &
      &-1.18327968399213E-17_wp, 3.80222517365649E-17_wp, 4.19308072410649E-05_wp, &
      & 1.48116799764167E-17_wp,-2.18654451169821E-16_wp,-6.66522418391771E-20_wp, &
      &-1.28150022279692E-17_wp, 4.57742663566695E-17_wp,-8.40533438311804E-05_wp, &
      & 1.81793703057341E-17_wp,-2.48521077750157E-16_wp,-5.87042029346142E-32_wp, &
      & 5.16421185233517E-03_wp, 1.56379286752774E-17_wp,-5.50832693069852E-05_wp, &
      & 6.00741164441656E-18_wp,-1.15506057334516E-17_wp,-3.83744009103750E-17_wp, &
      & 4.19308072402213E-05_wp,-1.48868058528040E-17_wp,-2.18426546153779E-16_wp, &
      & 3.49417712234010E-19_wp, 1.30749424216024E-17_wp, 4.69467077825026E-17_wp, &
      & 8.40533438312708E-05_wp, 1.82338467056468E-17_wp, 2.48604510742080E-16_wp, &
      & 5.58716419482001E-32_wp,-1.03972787565093E-17_wp, 1.00613961606655E-16_wp, &
      & 9.07546517214644E-19_wp,-1.18968212390597E-01_wp, 9.46906247451355E-18_wp, &
      & 2.77555756156289E-17_wp,-5.96147723349914E-18_wp,-2.90468841606325E-01_wp, &
      & 7.18952644834967E-18_wp, 4.58660943678858E-19_wp, 1.63513899389317E-17_wp, &
      &-1.24900090270330E-16_wp,-9.15378554962150E-18_wp,-3.56290394779576E-01_wp, &
      & 2.00362064254378E-17_wp,-1.25654370115671E-15_wp,-2.64794821941287E-18_wp, &
      & 9.36750677027476E-17_wp, 3.15696968400295E-18_wp,-1.18968212390508E-01_wp, &
      & 1.68543595984857E-17_wp,-1.07651953126608E-17_wp,-1.48868058528040E-17_wp, &
      & 2.90468841606490E-01_wp, 7.34600634682007E-18_wp,-9.20441925479908E-18_wp, &
      &-2.75946599399649E-17_wp,-9.20934748341699E-18_wp, 5.68142026296744E-17_wp, &
      &-3.56290394780214E-01_wp,-3.28001084092886E-17_wp,-1.14796595884815E-15_wp, &
      &-1.03034746724869E-17_wp, 4.53843954983374E-17_wp, 8.49879715927887E-17_wp, &
      &-1.56358560639567E-17_wp, 2.01304045973042E-02_wp, 5.28531762240424E-17_wp, &
      &-2.62802018034202E-16_wp,-2.40296837886836E-17_wp, 3.76577416051279E-01_wp, &
      & 1.34041121706306E-17_wp, 2.28970509350364E-02_wp,-8.28562749591696E-17_wp, &
      &-5.73604538202771E-16_wp, 8.53444908447583E-17_wp, 4.28332785594644E-01_wp, &
      & 2.38528104932060E-21_wp,-4.25605351784397E-17_wp, 5.29445869870659E-17_wp, &
      &-9.48694541679830E-17_wp,-1.00694331958914E-17_wp, 2.01304045973181E-02_wp, &
      & 4.22831335559007E-17_wp,-2.18426546153779E-16_wp, 7.34600634682007E-18_wp, &
      & 3.76577416051544E-01_wp, 2.14464105513430E-16_wp,-2.28970509350817E-02_wp, &
      & 7.47080579163392E-18_wp, 5.63146715764963E-16_wp, 8.59425425215767E-18_wp, &
      &-4.28332785595496E-01_wp,-1.89338404035048E-19_wp, 4.30343990420659E-17_wp, &
      &-5.49434910072695E-16_wp, 4.59019780725063E-19_wp, 3.76988216285415E-18_wp, &
      & 1.14644400634831E-17_wp,-1.34148205356663E-15_wp, 3.49417712241383E-19_wp, &
      & 9.20441925479386E-18_wp, 2.14464105513279E-16_wp, 6.29854570369907E-30_wp, &
      & 1.30400691554111E-17_wp,-1.64546795384943E-15_wp,-7.00433143079671E-19_wp, &
      & 1.12901822855475E-17_wp, 2.43939237482200E-16_wp, 3.92400910952554E-32_wp, &
      & 4.30343990421766E-17_wp,-5.49434910072282E-16_wp,-4.59019780691566E-19_wp, &
      & 3.76988216285108E-18_wp, 1.14644400634911E-17_wp, 1.34148205356740E-15_wp, &
      & 3.49417712234010E-19_wp,-9.20441925479908E-18_wp, 2.14464105513430E-16_wp, &
      & 6.32075444526208E-30_wp,-1.30400691554368E-17_wp,-1.64546795385237E-15_wp, &
      & 7.00433143079628E-19_wp, 1.12901822855678E-17_wp,-2.43939237482686E-16_wp, &
      & 3.56842536578744E-32_wp,-2.47513382605561E-17_wp, 1.78524787628581E-17_wp, &
      &-5.43821547771220E-18_wp, 1.20697971714131E-17_wp,-1.22398975552264E-03_wp, &
      & 4.71119123131707E-17_wp, 1.57731059242531E-17_wp, 2.86090782132920E-17_wp, &
      &-2.28970509350655E-02_wp,-8.15010741559064E-19_wp,-1.39221025790598E-03_wp, &
      & 6.67674555525528E-17_wp, 3.52899510447828E-17_wp, 2.81106557954278E-17_wp, &
      &-2.60439346356967E-02_wp,-1.45032334145789E-22_wp,-2.27900109851700E-17_wp, &
      & 1.73927961525273E-17_wp, 6.03903963812852E-18_wp, 1.17313417502787E-17_wp, &
      &-1.22398975552348E-03_wp,-5.28964889711605E-17_wp, 1.30749424216024E-17_wp, &
      &-2.75946599399649E-17_wp,-2.28970509350817E-02_wp,-1.30400691554368E-17_wp, &
      & 1.39221025790873E-03_wp, 6.12752946300728E-17_wp,-3.46540836697844E-17_wp, &
      & 3.27773033469162E-17_wp, 2.60439346357485E-02_wp, 1.15123501739801E-20_wp, &
      & 7.29241956472982E-17_wp, 1.45926947360386E-01_wp, 4.27580774117376E-17_wp, &
      &-1.55698816381950E-16_wp, 7.13010004780492E-18_wp, 3.56290394780011E-01_wp, &
      &-5.11741672519457E-17_wp, 3.16979027818925E-19_wp,-4.50284848129400E-17_wp, &
      &-1.67098364533931E-15_wp,-8.37638429387216E-17_wp, 4.37027409584909E-01_wp, &
      &-7.99192670097916E-17_wp, 4.02277783928515E-17_wp, 7.37887815336974E-17_wp, &
      & 1.53728397821890E-19_wp,-3.39711733648326E-17_wp, 1.45926947360276E-01_wp, &
      &-5.54240511801882E-17_wp,-8.64667514399829E-17_wp, 1.15875417413683E-18_wp, &
      &-3.56290394780214E-01_wp, 4.69467077825026E-17_wp,-9.20934748341699E-18_wp, &
      & 7.47080579163392E-18_wp,-1.64546795385237E-15_wp, 6.12752946300728E-17_wp, &
      & 4.37027409585691E-01_wp, 2.89380045597980E-16_wp,-1.08544100079492E-17_wp, &
      & 5.76465385328068E-18_wp, 1.55327732149280E-19_wp, 1.03520371538272E-02_wp, &
      & 9.66556684429014E-17_wp, 1.10418407836107E-04_wp,-2.32694314675766E-17_wp, &
      & 2.97898986809078E-17_wp, 2.35415528184118E-16_wp, 8.40533438329617E-05_wp, &
      &-5.69647978967860E-17_wp, 5.62689863632493E-16_wp,-1.33609252222772E-19_wp, &
      & 3.51751526380207E-17_wp, 2.87029801737001E-16_wp,-1.68491023049833E-04_wp, &
      &-6.98506337535823E-17_wp, 6.40393485354896E-16_wp,-4.92467847758371E-31_wp, &
      & 1.03520371538538E-02_wp, 9.65502330616313E-17_wp,-1.10418407828060E-04_wp, &
      &-2.34496235339006E-17_wp, 3.03555712176634E-17_wp,-2.36121436740097E-16_wp, &
      & 8.40533438312708E-05_wp, 5.68142026296744E-17_wp, 5.63146715764963E-16_wp, &
      & 7.00433143079628E-19_wp,-3.46540836697844E-17_wp, 2.89380045597980E-16_wp, &
      & 1.68491023050014E-04_wp,-6.97414318596541E-17_wp,-6.40226237874498E-16_wp, &
      &-2.30909057673474E-31_wp, 9.50478865760614E-18_wp,-1.29963933509606E-16_wp, &
      &-1.14785093191382E-18_wp, 1.45926947360385E-01_wp,-1.06737035373510E-17_wp, &
      &-5.08621862915076E-17_wp, 7.28599825578279E-18_wp, 3.56290394780011E-01_wp, &
      & 8.78619320330032E-18_wp,-5.62595587842681E-19_wp,-1.89862575151744E-17_wp, &
      & 1.28365872778174E-16_wp, 1.12809487216846E-17_wp, 4.37027409584909E-01_wp, &
      &-4.55205213367850E-18_wp, 1.54128218664578E-15_wp,-5.75226095149484E-22_wp, &
      &-1.19389523555573E-16_wp,-3.83770302256159E-18_wp, 1.45926947360276E-01_wp, &
      &-1.97325423901384E-17_wp, 2.88749796181813E-17_wp, 1.82338467056468E-17_wp, &
      &-3.56290394780214E-01_wp, 8.59425425215767E-18_wp, 1.12901822855678E-17_wp, &
      & 3.27773033469162E-17_wp,-1.08544100079492E-17_wp,-6.97414318596541E-17_wp, &
      & 4.37027409585691E-01_wp, 2.02083120326713E-17_wp, 1.40810023687965E-15_wp, &
      & 3.12096455677403E-17_wp,-4.68595854827922E-17_wp,-9.64605074106201E-17_wp, &
      & 2.77965952914473E-17_wp,-2.28970509350658E-02_wp,-4.84897181188293E-17_wp, &
      & 2.99078778907051E-16_wp, 5.17767271018625E-17_wp,-4.28332785595195E-01_wp, &
      &-1.52463224286792E-17_wp,-2.60439346356971E-02_wp, 1.08505945453719E-16_wp, &
      & 6.52121343701652E-16_wp,-6.70901855233864E-17_wp,-4.87201216523616E-01_wp, &
      &-2.71310501558058E-21_wp, 6.78999939336255E-17_wp,-5.54588209705638E-17_wp, &
      & 1.07700063413821E-16_wp, 2.14651444097584E-17_wp,-2.28970509350817E-02_wp, &
      &-5.97217646795548E-17_wp, 2.48604510742080E-16_wp,-3.28001084092886E-17_wp, &
      &-4.28332785595496E-01_wp,-2.43939237482686E-16_wp, 2.60439346357485E-02_wp, &
      & 5.76465385328068E-18_wp,-6.40226237874498E-16_wp, 2.02083120326713E-17_wp, &
      & 4.87201216524586E-01_wp, 2.15360355039118E-19_wp,-3.39346874026659E-31_wp, &
      & 5.18651720592585E-20_wp,-7.67247719260822E-33_wp, 4.70175930933210E-16_wp, &
      &-1.01213150778144E-20_wp, 1.26632283911305E-19_wp, 2.05851637011369E-32_wp, &
      & 1.14796595884750E-15_wp,-1.89338404034914E-19_wp,-2.41331863526078E-33_wp, &
      &-1.15123501739130E-20_wp, 1.55327732149428E-19_wp, 4.26293631624487E-32_wp, &
      & 1.40810023687713E-15_wp,-2.15360355038639E-19_wp, 4.96600393445783E-30_wp, &
      &-3.69994885462225E-31_wp, 5.18651720592461E-20_wp,-8.39048739596107E-33_wp, &
      & 4.70175930932858E-16_wp,-1.01213150778506E-20_wp,-1.26632283911453E-19_wp, &
      & 5.58716419482001E-32_wp,-1.14796595884815E-15_wp,-1.89338404035048E-19_wp, &
      & 3.56842536578744E-32_wp, 1.15123501739801E-20_wp, 1.55327732149280E-19_wp, &
      &-2.30909057673474E-31_wp, 1.40810023687965E-15_wp, 2.15360355039118E-19_wp, &
      & 4.53689242217229E-30_wp], shape(density))

   call get_structure(mol, "f-block", "Ce2")
   call test_e_gen(error, mol, ref, thr_in=thr1, guess_qat=qat, guess_qsh=qsh, &
      & guess_density=density, guess_dp=dp, guess_qp=qp)

end subroutine test_e_ce2

subroutine test_e_mb01(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol
   ! Reference energy calculated with g-xTB development code
   real(wp), parameter :: ref = -1278.38829904121_wp 

   call get_structure(mol, "MB16-43", "01")
   call test_e_gen(error, mol, ref, thr_in=thr1*10)

end subroutine test_e_mb01

subroutine test_e_mb02(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol
   ! Reference energy calculated with g-xTB development code
   real(wp), parameter :: ref = -1522.66113337614_wp 

   call get_structure(mol, "MB16-43", "02")
   call test_e_gen(error, mol, ref, thr_in=thr1*10)

end subroutine test_e_mb02


subroutine test_g_h2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_numgrad(error, mol, thr_in=thr1*10, four_in=.true.)

end subroutine test_g_h2

subroutine test_g_lih(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numgrad(error, mol, thr_in=thr1, four_in=.true.)

end subroutine test_g_lih

subroutine test_g_water(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol

   call get_structure(mol, "ICE10", "gas")
   call test_numgrad(error, mol, thr_in=thr1, four_in=.true.)

end subroutine test_g_water

subroutine test_g_no(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "NO")
   call test_numgrad(error, mol, thr_in=thr1, four_in=.true.)

end subroutine test_g_no

subroutine test_g_s2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "S2")
   call test_numgrad(error, mol, thr_in=thr1*10, four_in=.true.)

end subroutine test_g_s2

subroutine test_g_pcl(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "PCl")
   call test_numgrad(error, mol, thr_in=thr1*10, four_in=.true.)

end subroutine test_g_pcl

subroutine test_g_sih4(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_numgrad(error, mol, thr_in=thr1*10, four_in=.true.)

end subroutine test_g_sih4

subroutine test_g_cecl3(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol

   call get_structure(mol, "f-block", "CeCl3")
   call test_numgrad(error, mol, thr_in=thr2, four_in=.true.)

end subroutine test_g_cecl3

subroutine test_g_ce2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol

   call get_structure(mol, "f-block", "Ce2")
   call test_numgrad(error, mol, thr_in=thr1*100, four_in=.true.)

end subroutine test_g_ce2

subroutine test_g_panp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol

   call get_structure(mol, "f-block", "PaNp")
   call test_numgrad(error, mol, thr_in=thr1*10, four_in=.true.)

end subroutine test_g_panp

subroutine test_g_mb02(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_numgrad(error, mol, thr_in=thr2)

end subroutine test_g_mb02




! subroutine test_s_gen(error, mol)
!    !> Error handling
!    type(error_type), allocatable, intent(out) :: error
!    !> Molecular structure data
!    type(structure_type), intent(inout) :: mol

!    type(context_type) :: ctx
!    type(xtb_calculator) :: calc
!    type(wavefunction_type) :: wfn, wfn_aux
!    real(wp) :: energy, er, el
!    real(wp), allocatable :: gradient(:, :), numsigma(:, :), sigma(:, :)
!    real(wp), allocatable :: eps(:, :), xyz(:, :), lat(:, :)
!    real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
!    real(wp), parameter :: unity(3, 3) = reshape(&
!    & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
!    real(wp), parameter :: step = 1.0e-6_wp
!    integer :: ic, jc

!    allocate(gradient(3, mol%nat), numsigma(3, 3), sigma(3, 3), &
!       & eps(3, 3), xyz(3, mol%nat), lat(3, 3))
!    allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
!    energy = 0.0_wp
!    gradient(:, :) = 0.0_wp
!    sigma(:, :) = 0.0_wp

!    ! Setup g-xTB calculator
!    call new_gxtb_calculator(calc, mol, error)
!    if (allocated(error)) return

!    ! Setup auxiliary wavefunction
!    call new_wavefunction(wfn_aux, mol%nat, calc%bas%nsh, 0, 1, 0.0_wp, .true.)

!    ! Setup g-xTB wavefunction
!    call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
!       & 1, calc%default_etemp * kt, .false.)

!    if(.not.derive_basis) then
!       ! Obtain EEQBC charges and derivatives
!       call get_charges(calc%charge_model, mol, error, wfn_aux%qat(:, 1))
!       if (allocated(error)) return
      
!       ! Calculate the coordination number required for basis set scaling
!       if (allocated(calc%bas%ncoord)) then
!          allocate(cn(mol%nat))
!          call calc%bas%ncoord%get_cn(mol, cn) 
!       end if

!       ! Scale the basis set only once with its charge and CN dependence
!       call calc%bas%scale_basis(mol, wfn_aux%qat(:, 1), cn, .true.)   
!    end if 

!    eps(:, :) = unity
!    xyz(:, :) = mol%xyz
!    if (any(mol%periodic)) &
!    lat(:, :) = mol%lattice
!    do ic = 1, 3
!       do jc = 1, 3
!          er = 0.0_wp
!          el = 0.0_wp

!          ! Right hand
!          eps(jc, ic) = eps(jc, ic) + step
!          mol%xyz(:, :) = matmul(eps, xyz)
!          if (any(mol%periodic)) &
!          mol%lattice(:, :) = matmul(eps, lat)

!          ! Obtain EEQBC charges and derivatives
!          call get_charges(calc%charge_model, mol, error, wfn_aux%qat(:, 1))
!          if (allocated(error)) return

!          if(derive_basis) then
!             ! Calculate the coordination number required for basis set scaling
!             call calc%bas%ncoord%get_cn(mol, cn) 
!             ! Scale the basis set with its charge and CN dependence
!             call calc%bas%scale_basis(mol, wfn_aux%qat(:, 1), cn, .true.)
!          end if

!          ! Use EEQBC charges as a guess
!          wfn%qat(:, 1) = wfn_aux%qat(:, 1)

!          ! Perform g-xTB calculation
!          call xtb_singlepoint(ctx, mol, calc, wfn, acc, er, verbosity=0, wfn_aux=wfn_aux)
         
!          ! Left hand
!          eps(jc, ic) = eps(jc, ic) - 2*step
!          mol%xyz(:, :) = matmul(eps, xyz)
!          if (any(mol%periodic)) &
!          mol%lattice(:, :) = matmul(eps, lat)

!          ! Obtain EEQBC charges and derivatives
!          call get_charges(calc%charge_model, mol, error, wfn_aux%qat(:, 1))
!          if (allocated(error)) return

!          if(derive_basis) then
!             ! Calculate the coordination number required for basis set scaling
!             call calc%bas%ncoord%get_cn(mol, cn) 
!             ! Scale the basis set with its charge and CN dependence
!             call calc%bas%scale_basis(mol, wfn_aux%qat(:, 1), cn, .true.)
!          end if 

!          ! Use EEQBC charges as a guess
!          wfn%qat(:, 1) = wfn_aux%qat(:, 1)

!          ! Perform g-xTB calculation
!          call xtb_singlepoint(ctx, mol, calc, wfn, acc, el, verbosity=0, wfn_aux=wfn_aux)
         
!          eps(jc, ic) = eps(jc, ic) + step
!          mol%xyz(:, :) = xyz
!          if (any(mol%periodic)) &
!          mol%lattice(:, :) = lat
!          numsigma(jc, ic)  = 0.5_wp*(er - el)/step
!       end do
!    end do

!    ! Obtain EEQBC charges and derivatives
!    call get_charges(calc%charge_model, mol, error, wfn_aux%qat(:, 1), &
!       & wfn_aux%dqatdr(:, :, :, 1), wfn_aux%dqatdL(:, :, :, 1))
!    if (allocated(error)) return
 
!    if(derive_basis) then
!       ! Calculate the coordination number with gradient required for basis set scaling
!       call calc%bas%ncoord%get_cn(mol, cn, dcndr, dcndL)
!       ! Scale the basis set with its charge and CN dependence
!       call calc%bas%scale_basis(mol, wfn_aux%qat(:, 1), cn, .true., &
!       & dqatdr=wfn_aux%dqatdr(:, :, :, 1), dqatdL=wfn_aux%dqatdL(:, :, :, 1), &
!       & dcndr=dcndr, dcndL=dcndL)
!    else
!       ! Calculate the coordination number required for basis set scaling
!       call calc%bas%ncoord%get_cn(mol, cn)
!       ! Scale the basis set with its charge and CN dependence
!       call calc%bas%scale_basis(mol, wfn_aux%qat(:, 1), cn, .true.)
!    end if

!    ! g-xTB calculation
!    call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, &
!       & verbosity=0, wfn_aux=wfn_aux)

!    if (any(abs(sigma - numsigma) > thr2)) then
!       call test_failed(error, "Sigma of energy does not match")
!       print'(3es21.14)', sigma
!       print'("---")'
!       print'(3es21.14)', numsigma
!       print'("---")'
!       print'(3es21.14)', sigma-numsigma
!    end if

! end subroutine test_s_gen

! subroutine test_s_lih(error)

!    !> Error handling
!    type(error_type), allocatable, intent(out) :: error

!    type(structure_type) :: mol

!    call get_structure(mol, "MB16-43", "LiH")
!    call test_s_gen(error, mol)

! end subroutine test_s_lih

! subroutine test_s_mb03(error)

!    !> Error handling
!    type(error_type), allocatable, intent(out) :: error

!    type(structure_type) :: mol

!    call get_structure(mol, "MB16-43", "03")
!    call test_s_gen(error, mol)

! end subroutine test_s_mb03

end module test_gxtb