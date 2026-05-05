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

module test_ceh
   use mctc_data_paulingen, only : get_pauling_en
   use mctc_data_covrad, only : get_covalent_rad
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, &
      & check, test_failed
   use mctc_io, only : structure_type, new
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use mstore, only : get_structure
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_cache, only : basis_cache 
   use tblite_basis_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_blas, only: gemv
   use tblite_ceh_singlepoint, only : ceh_singlepoint
   use tblite_ceh_ceh, only : ceh_h0spec, new_ceh_calculator
   use tblite_container, only : container_type
   use tblite_context_type, only : context_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_external_field, only : electric_field
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : tb_hamiltonian, new_hamiltonian, get_hamiltonian, &
      & get_selfenergy
   implicit none
   private

   public :: collect_ceh

   real(wp), parameter :: kt = 3.166808578545117e-06_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = 1e2*sqrt(epsilon(1.0_wp))

contains

!> Collect all exported unit tests
subroutine collect_ceh(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("q-mol-h2", test_q_h2), &
      new_unittest("q-mol-lih", test_q_lih), &
      new_unittest("q-mol-sih4", test_q_sih4), &
      new_unittest("q-mol-cecl3", test_q_cecl3), &
      new_unittest("q-mol-accl6", test_q_accl6), &
      new_unittest("q-mol-panp", test_q_panp), &
      new_unittest("q-mol-mb01", test_q_mb01), &
      new_unittest("q-mol-mb02", test_q_mb02), &
      new_unittest("q-mol-mb03", test_q_mb03), &
      new_unittest("q-mol-mb04", test_q_mb04), &
      new_unittest("q-chrgd-efield-mol", test_q_ef_chrg_mb01), &
      new_unittest("d-mol", test_d_mb01), &
      new_unittest("d-field-mol", test_d_field_mb04), &
      new_unittest("d-field-change-mol", test_d_hcn) &
      ]

end subroutine collect_ceh


subroutine test_q_gen(error, mol, ref, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Reference CEH charges
   real(wp), intent(in) :: ref(:)

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   integer :: i
   real(wp), allocatable :: cn(:)
   real(wp), parameter :: accuracy = 1e-8_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate(cn(mol%nat))

   call new_ceh_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)
   ctx%verbosity = 0
   call ceh_singlepoint(ctx, calc, mol, wfn, accuracy)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if

   do i = 1, mol%nat
      call check(error, wfn%qat(i,1), ref(i), thr=thr_)
      if (allocated(error)) then
         print '(3es21.13)',  wfn%qat(i,1), ref(i), &
            & wfn%qat(i,1) - ref(i)
         return
      end if
   enddo

end subroutine test_q_gen


subroutine test_q_h2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! calculated with GP3 standalone (full matrix diagonalization)
   real(wp), parameter :: charges(2) = reshape([& 
      & 0.0000000000000_wp, 0.00000000000000_wp],&
      & shape(charges))

   call get_structure(mol, "MB16-43", "H2")
   call test_q_gen(error, mol, charges, thr_in=thr2)

end subroutine test_q_h2

subroutine test_q_lih(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! calculated with GP3 standalone (full matrix diagonalization)
   real(wp), parameter :: charges(2) = reshape([&
      & 0.452928818506870_wp, -0.452928818506872_wp],&
      & shape(charges))

   call get_structure(mol, "MB16-43", "LiH")
   call test_q_gen(error, mol, charges, thr_in=thr2)

end subroutine test_q_lih

subroutine test_q_sih4(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! calculated with GP3 standalone (full matrix diagonalization)
   real(wp), parameter :: charges(5) = reshape([&
      & 0.268915553151106_wp, -0.067228888686111_wp, -0.067228888686110_wp, &
      &-0.067228888686110_wp, -0.067228888686110_wp], shape(charges))

   call get_structure(mol, "MB16-43", "SiH4")
   call test_q_gen(error, mol, charges, thr_in=thr2)

end subroutine test_q_sih4

subroutine test_q_cecl3(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! calculated with GP3 standalone (full matrix diagonalization)
   real(wp), parameter :: charges(4) = reshape([&
      & 0.941257219125013_wp, -0.312389885339237_wp, -0.316670447603892_wp, &
      &-0.312196886181877_wp], shape(charges))

   call get_structure(mol, "f-block", "CeCl3")
   call test_q_gen(error, mol, charges, thr_in=thr2)

end subroutine test_q_cecl3

subroutine test_q_accl6(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! calculated with GP3 standalone (full matrix diagonalization)
   real(wp), parameter :: charges(7) = reshape([&
      &  0.288810442976471_wp, -0.04820498331048_wp, -0.048312400345871_wp, &
      & -0.047844815375011_wp, -0.04788094096742_wp, -0.048243946234083_wp, &
      & -0.048323356740988_wp], shape(charges))

   call get_structure(mol, "f-block", "AcCl6")
   call test_q_gen(error, mol, charges, thr_in=thr2)

end subroutine test_q_accl6

subroutine test_q_panp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! calculated with GP3 standalone (full matrix diagonalization)
   real(wp), parameter :: charges(2) = reshape([&
      & -0.490095851470747_wp, 0.490095852181341_wp], shape(charges))

   call get_structure(mol, "f-block", "PaNp")
   call test_q_gen(error, mol, charges, thr_in=thr2)

end subroutine test_q_panp

subroutine test_q_mb01(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! calculated with GP3 standalone (full matrix diagonalization)
   real(wp), parameter :: charges(16) = reshape([&
      &  0.50933743182523_wp, -0.063950757122457_wp, -0.447396210062547_wp, &     
      & -0.06007626073194_wp, -0.228995347028063_wp,  0.081936572631240_wp, &
      & -0.04029707489635_wp, -0.384822906853029_wp, -0.214508333206973_wp, &     
      &  0.14648324095015_wp,  0.090840217217610_wp,  0.034875957186194_wp, &
      & -0.05930815144452_wp,  0.133798380818110_wp, -0.063944989141738_wp, &
      &  0.56602822987298_wp], shape(charges))

   call get_structure(mol, "MB16-43", "01")
   call test_q_gen(error, mol, charges, thr_in=thr2)

end subroutine test_q_mb01

subroutine test_q_mb02(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! calculated with GP3 standalone (full matrix diagonalization)
   real(wp), parameter :: charges(16) = reshape([&
      &-0.085528757393776_wp, -0.069997806147367_wp, -0.234351594245719_wp, &
      &-0.167766349660902_wp,  0.432632782697259_wp,  0.171300109658293_wp, &
      &-0.088439445405627_wp, -0.045209869368921_wp,  0.391614741682698_wp, &
      & 0.148323670914553_wp, -0.084649073627366_wp,  0.371028349885743_wp, &
      &-0.347496611788866_wp, -0.086873903658569_wp, -0.001912803785713_wp, &
      &-0.302673439755712_wp], shape(charges))

   call get_structure(mol, "MB16-43", "02")
   call test_q_gen(error, mol, charges, thr_in=thr2)

end subroutine test_q_mb02

subroutine test_q_mb03(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! calculated with GP3 standalone (full matrix diagonalization)
   real(wp), parameter :: charges(16) = reshape([&
      & 0.079458459812547_wp, -0.514984986958578_wp,  0.027119814956294_wp, &
      & 0.294936546309030_wp,  0.396051282334921_wp,  0.033859363521586_wp, &
      &-0.261951753885935_wp,  0.026314869404998_wp,  0.037212079428442_wp, &
      &-0.005505860141304_wp, -0.364181487304452_wp, -0.142392172388313_wp, &     
      &-0.286416123354516_wp,  0.100899051623118_wp,  0.558735071392839_wp, &     
      & 0.020845845249307_wp], shape(charges))

   call get_structure(mol, "MB16-43", "03")
   call test_q_gen(error, mol, charges, thr_in=thr2)

end subroutine test_q_mb03

subroutine test_q_mb04(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! calculated with GP3 standalone (full matrix diagonalization)
   real(wp), parameter :: charges(16) = reshape([&
      &-0.018039040440414_wp, -0.194798993283821_wp,  -0.076478151155481_wp, &
      &-0.167853095449230_wp,  0.304619635473503_wp,  -0.022689763513518_wp, &
      & 0.019714590105775_wp,  0.000713111871502_wp,  -0.047985458556249_wp, &
      &-0.104183443069051_wp, -0.161265309374659_wp,   0.538802940610920_wp, &
      & 0.190213172949517_wp, -0.317990923047223_wp,   0.034501918582603_wp, &
      & 0.022718808299253_wp], shape(charges))

   call get_structure(mol, "MB16-43", "04")
   call test_q_gen(error, mol, charges, thr_in=thr2)

end subroutine test_q_mb04

subroutine test_q_ef_chrg_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp), parameter :: accuracy = 1e-8_wp
   class(container_type), allocatable :: cont      
   ! calculated with GP3 standalone (full matrix diagonalization)
   real(wp), parameter :: ref(16) = reshape([&
      &-5.42237346896788_wp, -0.77304500586496_wp,   2.5895850175165_wp, &
      &-0.92233780581096_wp,  6.99267602990832_wp,   0.4742366118103_wp, &
      &-0.11849846722517_wp,  4.22307140408149_wp,   1.5873873640455_wp, &
      & 0.31672778030780_wp,  0.99906183347402_wp, -10.5405405662106_wp, &
      &-3.80217066006454_wp,  1.92138378495190_wp,   3.8481284909192_wp, &
      & 0.62670765712891_wp], shape(ref))

   real(wp) :: efield(3)
   integer :: i

   efield = 0.0_wp
   efield(3) = 0.2_wp

   call get_structure(mol, "MB16-43", "01")
   mol%charge = 2.0_wp
   call new_ceh_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)
   cont = electric_field(efield)
   call calc%push_back(cont)
   ctx%verbosity = 0
   call ceh_singlepoint(ctx, calc, mol, wfn, accuracy)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if

   do i = 1, mol%nat
      call check(error, wfn%qat(i,1), ref(i), thr=5e-6_wp, &
         & message="Calculated charge does not match reference")
      if (allocated(error)) return
   enddo

end subroutine test_q_ef_chrg_mb01

subroutine test_d_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp) :: dipole(3), tmp(3)
   real(wp), parameter :: accuracy = 1e-8_wp
   ! calculated with GP3 standalone (full matrix diagonalization)
   real(wp), parameter :: ref(3) = reshape([&
      & 0.201892497728508_wp, -1.15519893399684_wp, -1.91423938957019_wp],&
      & shape(ref))

   call get_structure(mol, "MB16-43", "01")

   call new_ceh_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)
   ctx%verbosity = 0
   call ceh_singlepoint(ctx, calc, mol, wfn, accuracy)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if
   tmp = 0.0_wp
   dipole = 0.0_wp
   call gemv(mol%xyz, wfn%qat(:, 1), tmp)
   dipole(:) = tmp + sum(wfn%dpat(:, :, 1), 2)

   if (any(abs(dipole - ref) > 1e-5_wp)) then
      call test_failed(error, "Numerical dipole moment does not match")
      print '(3es21.14)', dipole
      print '("---")'
      print '(3es21.14)', ref
      print '("---")'
      print '(3es21.14)', dipole - ref
   end if

end subroutine test_d_mb01

subroutine test_d_field_mb04(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   real(wp) :: energy, efield(3), dipole(3), tmp(3)
   real(wp), parameter :: accuracy = 1e-8_wp
   ! calculated with GP3 standalone (full matrix diagonalization)
   real(wp), parameter :: ref(3) = reshape([& 
      & -7.6402587223855_wp,  83.5065044491344_wp,  0.55047274934631_wp],&
      & shape(ref))

   call get_structure(mol, "MB16-43", "04")
   energy = 0.0_wp
   efield = 0.0_wp
   efield(2) = 0.2_wp

   call new_ceh_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)

   cont = electric_field(efield)
   call calc%push_back(cont)

   ctx%verbosity = 0
   call ceh_singlepoint(ctx, calc, mol, wfn, accuracy)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if
   tmp = 0.0_wp
   dipole = 0.0_wp
   call gemv(mol%xyz, wfn%qat(:, 1), tmp)
   dipole(:) = tmp + sum(wfn%dpat(:, :, 1), 2)

   if (any(abs(dipole - ref) > 1e-5_wp)) then
      call test_failed(error, "Numerical dipole moment does not match")
      print '(3es21.14)', dipole
      print '("---")'
      print '(3es21.14)', ref
      print '("---")'
      print '(3es21.14)', dipole - ref
   end if

end subroutine test_d_field_mb04

subroutine test_d_hcn(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol1,mol2
   type(xtb_calculator) :: calc1,calc2
   type(wavefunction_type) :: wfn1,wfn2
   class(container_type), allocatable :: cont1,cont2
   real(wp), parameter :: accuracy = 1e-8_wp
   real(wp) :: efield(3), dip1(3), dip2(3), tmp(3)
   integer, parameter :: num(3) = reshape([ &
      7, &
      6, &
      1], shape(num))
   integer, parameter :: nat = 3
   real(wp) :: xyz(3, nat) = reshape([ &
      & -0.09604091224796_wp,  0.0_wp, 0.0_wp, &
      &  2.09604091224796_wp,  0.0_wp, 0.0_wp, &
      &  4.10859879422050_wp,  0.0_wp, 0.0_wp], &
      & shape(xyz))

   ctx%verbosity = 0
   call new(mol1, num, xyz)
   efield = 0.0_wp
   efield(1) = -0.1_wp
   call new_ceh_calculator(calc1, mol1, error)
   if (allocated(error)) return
   call new_wavefunction(wfn1, mol1%nat, calc1%bas%nsh, calc1%bas%nao, &
      & 1, calc1%default_etemp * kt)
   cont1 = electric_field(efield)
   call calc1%push_back(cont1)
   call ceh_singlepoint(ctx, calc1, mol1, wfn1, accuracy)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if
   tmp = 0.0_wp
   dip1 = 0.0_wp
   call gemv(mol1%xyz, wfn1%qat(:, 1), tmp)
   dip1(:) = tmp + sum(wfn1%dpat(:, :, 1), 2)

   xyz(1, :) = xyz(1, :) - 1.0_wp
   call new(mol2, num, xyz)
   call new_ceh_calculator(calc2, mol2, error)
   if (allocated(error)) return
   call new_wavefunction(wfn2, mol2%nat, calc2%bas%nsh, calc2%bas%nao, &
      & 1, calc2%default_etemp * kt)
   cont2 = electric_field(efield)
   call calc2%push_back(cont2)
   call ceh_singlepoint(ctx, calc2, mol2, wfn2, accuracy)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if
   tmp = 0.0_wp
   dip2 = 0.0_wp
   call gemv(mol2%xyz, wfn2%qat(:, 1), tmp)
   dip2(:) = tmp + sum(wfn2%dpat(:, :, 1), 2)

   if (any(abs(dip1 - dip2) > 1e-7_wp)) then
      call test_failed(error, "Numerical dipole moment does not match")
      print '(3es21.14)', dip1
      print '("---")'
      print '(3es21.14)', dip2
      print '("---")'
      print '(3es21.14)', dip1 - dip2
   end if

end subroutine test_d_hcn

end module test_ceh
