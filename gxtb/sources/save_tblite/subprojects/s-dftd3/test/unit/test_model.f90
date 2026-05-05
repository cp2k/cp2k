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

module test_model
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io_structure, only : structure_type
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count
   use mctc_data, only : get_covalent_rad
   use mstore, only : get_structure
   use dftd3_cutoff, only : get_lattice_points
   use dftd3_data, only : get_vdw_rad, get_r4r2_val
   use dftd3_model
   implicit none
   private

   public :: collect_model

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_model(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("r4r2-val", test_r4r2_val), &
      & new_unittest("vdw-rad", test_vdw_rad), &
      & new_unittest("gw-mb01", test_gw_mb01), &
      & new_unittest("gw-mb02", test_gw_mb02), &
      & new_unittest("gw-mb03", test_gw_mb03), &
      & new_unittest("gw-amf3", test_gw_amf3), &
      & new_unittest("dgw-mb04", test_dgw_mb04), &
      & new_unittest("dgw-mb05", test_dgw_mb05) &
      & ]

end subroutine collect_model


subroutine test_gw_gen(error, mol, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Reference Gaussian weights
   real(wp), intent(in) :: ref(:, :)

   type(d3_model) :: d3
   real(wp), allocatable :: cn(:), rcov(:), gwvec(:, :)
   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :)
   class(ncoord_type), allocatable :: ncoord

   call new_d3_model(d3, mol)

   allocate(rcov(mol%nid), cn(mol%nat), gwvec(maxval(d3%ref), mol%nat))
   rcov(:) = get_covalent_rad(mol%num)

   call new_ncoord(ncoord, mol, cn_count%exp, error, cutoff=cutoff, rcov=rcov)
   if (allocated(error)) return
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call ncoord%get_coordination_number(mol, lattr, cn)

   call d3%weight_references(mol, cn, gwvec)

   if (any(abs(gwvec - ref) > thr)) then
      call test_failed(error, "Gaussian weights do not match")
      print'(3es21.14)', gwvec
   end if

end subroutine test_gw_gen


subroutine test_dgw_gen(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   integer :: iat, mref
   type(d3_model) :: d3
   real(wp), allocatable :: cn(:), rcov(:), gwvec(:, :), gwdcn(:, :)
   real(wp), allocatable :: gwr(:, :), gwl(:, :), numdcn(:, :)
   real(wp), parameter :: cutoff = 30.0_wp, lattr(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   class(ncoord_type), allocatable :: ncoord

   call new_d3_model(d3, mol)

   mref = maxval(d3%ref)
   allocate(rcov(mol%nid), cn(mol%nat), gwvec(mref, mol%nat), &
      & gwdcn(mref, mol%nat), gwr(mref, mol%nat), gwl(mref, mol%nat), &
      & numdcn(mref, mol%nat))
   rcov(:) = get_covalent_rad(mol%num)

   call new_ncoord(ncoord, mol, cn_count%exp, error, cutoff=cutoff, rcov=rcov)
   if (allocated(error)) return
   call ncoord%get_coordination_number(mol, lattr, cn)

   do iat = 1, mol%nat
      cn(iat) = cn(iat) + step
      call d3%weight_references(mol, cn, gwr)
      cn(iat) = cn(iat) - 2*step
      call d3%weight_references(mol, cn, gwl)
      cn(iat) = cn(iat) + step
      gwdcn(:, :) = 0.5_wp*(gwr - gwl)/step
      numdcn(:, iat) = gwdcn(:, iat)
      gwdcn(:, iat) = 0.0_wp
      if (any(abs(gwdcn) > thr)) then
         call test_failed(error, "Unexpected non-zero gradient element found")
         exit
      end if
   end do
   if (allocated(error)) return

   call d3%weight_references(mol, cn, gwvec, gwdcn)

   if (any(abs(gwdcn - numdcn) > thr2)) then
      call test_failed(error, "Gaussian weights derivatives do not match")
      print'(3es21.14)', gwdcn
      print'("---")'
      print'(3es21.14)', numdcn
      print'("---")'
      print'(3es21.14)', gwdcn - numdcn
   end if

end subroutine test_dgw_gen


subroutine test_gw_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(5, 16) = reshape([&
      & 4.61254014807976E-13_wp, 9.99999999999539E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 9.78431945472983E-01_wp, &
      & 2.15680545270172E-02_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 9.33252077840319E-08_wp, 1.55830681937747E-02_wp, &
      & 9.84416838481017E-01_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 9.99424904108747E-01_wp, 5.75095891252906E-04_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.35771400228363E-02_wp, &
      & 9.86422859977164E-01_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 9.82992148346892E-01_wp, 1.70078516531077E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 9.99519469064248E-01_wp, 4.80530935751615E-04_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.13181694792597E-07_wp, &
      & 1.71503960869602E-02_wp, 9.82849490731345E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.25926325849160E-25_wp, 6.73263145629432E-14_wp, &
      & 1.94165275506323E-05_wp, 9.99980583472382E-01_wp, 0.00000000000000E+00_wp, &
      & 9.86403420777318E-01_wp, 1.35965792226822E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 9.83377538259043E-01_wp, &
      & 1.66224617409573E-02_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 6.63636803493899E-06_wp, 9.99993363631965E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 4.78432084484299E-38_wp, 4.72470789879862E-24_wp, 2.64845507076682E-13_wp, &
      & 7.08386079833514E-06_wp, 9.99992916138937E-01_wp, 5.57929648633356E-26_wp, &
      & 1.48261370770972E-14_wp, 2.19715394953033E-06_wp, 1.59978977357256E-01_wp, &
      & 8.40018825488779E-01_wp, 1.11473605172390E-26_wp, 1.33471958830444E-14_wp, &
      & 8.80046323582265E-06_wp, 9.99991199536751E-01_wp, 0.00000000000000E+00_wp, &
      & 3.64404060381414E-41_wp, 1.64269207706493E-24_wp, 4.50618875164815E-11_wp, &
      & 9.99999999954938E-01_wp, 0.00000000000000E+00_wp], shape(ref))

   call get_structure(mol, "MB16-43", "01")
   call test_gw_gen(error, mol, ref)

end subroutine test_gw_mb01


subroutine test_gw_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(5, 16) = reshape([&
      & 9.76766039630376E-01_wp, 2.32339603696241E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 4.00470863539171E-21_wp, &
      & 3.29660160936946E-09_wp, 9.99999996703398E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 2.66415020809245E-25_wp, 4.94912877792885E-14_wp, &
      & 5.12174748254717E-06_wp, 2.60883387434529E-01_wp, 7.39111490817939E-01_wp, &
      & 1.58493787301861E-14_wp, 6.37917541998686E-06_wp, 9.99993620824564E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 3.61190804043511E-31_wp, &
      & 1.74570190245310E-14_wp, 9.99999999999982E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 9.88324870235887E-01_wp, 1.16751297641131E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 9.78240220424054E-01_wp, 2.17597795759463E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 9.78865958465924E-01_wp, &
      & 2.11340415340759E-02_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 7.88983927118728E-43_wp, 1.15799483190410E-27_wp, &
      & 3.54874780230275E-15_wp, 4.41039604316648E-06_wp, 9.99995589603953E-01_wp, &
      & 9.89809692831196E-01_wp, 1.01903071688044E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 5.76231094820794E-28_wp, &
      & 4.15283489804264E-16_wp, 1.67406060912827E-07_wp, 3.29932277130195E-02_wp, &
      & 9.67006604880919E-01_wp, 7.51015818491867E-13_wp, 9.99999999999249E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 6.11050341020826E-06_wp, 9.99993889496590E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 9.77573846919768E-01_wp, &
      & 2.24261530802315E-02_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 9.99856142640480E-01_wp, 1.43857359520024E-04_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 6.59504190399977E-08_wp, 1.33596297866852E-02_wp, 9.86640304262896E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp], shape(ref))

   call get_structure(mol, "MB16-43", "02")
   call test_gw_gen(error, mol, ref)

end subroutine test_gw_mb02


subroutine test_gw_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(5, 16) = reshape([&
      & 3.22985602658197E-29_wp, 5.04319495924654E-17_wp, 4.91517377051623E-08_wp, &
      & 1.20911775082148E-02_wp, 9.87908773340047E-01_wp, 1.76431281193830E-14_wp, &
      & 6.73115291498182E-06_wp, 9.99993268847067E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 9.85298197844874E-01_wp, 1.47018021551265E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.01470792750331E-15_wp, 9.99999999999999E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 8.43067420549278E-45_wp, &
      & 2.19850383316351E-21_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 3.13391714659211E-34_wp, 7.15449821607983E-20_wp, &
      & 9.85785135620685E-09_wp, 9.99999990142149E-01_wp, 0.00000000000000E+00_wp, &
      & 2.56990631870506E-41_wp, 3.99321363949352E-26_wp, 4.60982169768610E-14_wp, &
      & 1.23938288501679E-05_wp, 9.99987606171104E-01_wp, 9.95550037550049E-01_wp, &
      & 4.44996244995094E-03_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 9.83710989152179E-01_wp, 1.62890108478206E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 9.82790640769192E-01_wp, 1.72093592308080E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 6.65159426282988E-06_wp, &
      & 9.99993348405737E-01_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 8.28164794680018E-15_wp, 4.75279368702095E-06_wp, &
      & 9.99995247206305E-01_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 6.13460633571205E-41_wp, 7.68433459322833E-26_wp, 7.11247426761060E-14_wp, &
      & 1.53704140324608E-05_wp, 9.99984629585896E-01_wp, 9.91995069127435E-01_wp, &
      & 8.00493087256450E-03_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 5.69442053542165E-11_wp, 9.99999999943056E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 9.98369395301452E-01_wp, 1.63060469854839E-03_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp], shape(ref))

   call get_structure(mol, "MB16-43", "03")
   call test_gw_gen(error, mol, ref)

end subroutine test_gw_mb03

subroutine test_gw_amf3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(6, 4) = reshape([&
      & 3.01777419522501E-16_wp, 3.48560287705282E-08_wp, 6.05573875449400E-03_wp, &
      & 9.93942098041223E-01_wp, 2.12834822159420E-06_wp, 3.22312554200316E-14_wp, & 
      & 1.83164825589304E-02_wp, 9.81683517441070E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.0000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.83168129150391E-02_wp, 9.81683187084961E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.83165559699980E-02_wp, 9.81683444030002E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp], shape(ref))

   !> Molecular structure data 
   mol%nat = 4
   mol%nid = 2
   mol%id = [1, 2, 2, 2]
   mol%num = [95, 9]
   mol%xyz = reshape([ &
      & -1.13163973200000_wp, -2.17446990100000_wp, +1.10012477100000_wp, &
      & -4.66377948900000_wp, -3.12947883400000_wp, -0.36987606800000_wp, &
      & -0.19032564300000_wp, +1.36339950600000_wp, -0.36521789300000_wp, &
      & +1.46283310800000_wp, -4.75734549200000_wp, -0.36503081000000_wp], &
      & [3, 4])
   mol%periodic = [.false.]

   call test_gw_gen(error, mol, ref)

end subroutine test_gw_amf3


subroutine test_dgw_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_dgw_gen(error, mol)

end subroutine test_dgw_mb04


subroutine test_dgw_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_dgw_gen(error, mol)

end subroutine test_dgw_mb05


subroutine test_r4r2_val(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check(error, get_r4r2_val("H"), get_r4r2_val(1))
   if (allocated(error)) return
   call check(error, get_r4r2_val("Og"), get_r4r2_val(118))
   if (allocated(error)) return
   call check(error, get_r4r2_val("X"), get_r4r2_val(-1))

end subroutine test_r4r2_val


subroutine test_vdw_rad(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check(error, get_vdw_rad("He", "Rn"), get_vdw_rad(2, 86))
   if (allocated(error)) return
   call check(error, get_vdw_rad("S", "Fl"), get_vdw_rad(16, 114))
   if (allocated(error)) return
   call check(error, get_vdw_rad("Am", "U"), get_vdw_rad(95, 92))
   if (allocated(error)) return
   call check(error, get_vdw_rad("Og", "Cn"), get_vdw_rad(118, 112))
   if (allocated(error)) return
   call check(error, get_vdw_rad("X", "X"), get_vdw_rad(-1, -1))

end subroutine test_vdw_rad


end module test_model
