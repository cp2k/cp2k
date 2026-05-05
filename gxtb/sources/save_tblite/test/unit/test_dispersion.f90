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

module test_dispersion
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_container, only : container_cache
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_disp, only : dispersion_type, d3_dispersion, d4_dispersion, &
      & new_d3_dispersion, new_d4_dispersion, new_d4s_dispersion, &
      & new_d4srev_dispersion, twobody_damping_function, &
      & threebody_damping_function
   implicit none
   private

   public :: collect_dispersion


   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains


!> Collect all exported unit tests
subroutine collect_dispersion(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("e-d3-mb01", test_e_d3_mb01), &
      new_unittest("e-d4-mb01", test_e_d4_mb01), &
      new_unittest("e-d4s-mb01", test_e_d4s_mb01), &
      new_unittest("e-d4srev-mb01", test_e_d4srev_mb01), &
      new_unittest("e-d3-mb02", test_e_d3_mb02), &
      new_unittest("e-d4-mb02", test_e_d4_mb02), &
      new_unittest("e-d4s-mb02", test_e_d4s_mb02), &
      new_unittest("e-d4srev-mb02", test_e_d4srev_mb02), &
      new_unittest("p-d4-mb03", test_p_d4_mb03), &
      new_unittest("p-d4s-mb03", test_p_d4s_mb03), &
      new_unittest("p-d4srev-mb03", test_p_d4srev_mb03), &
      new_unittest("p-d4-mb04", test_p_d4_mb04), &
      new_unittest("p-d4srev-mb04", test_p_d4srev_mb04) &
      ]

end subroutine collect_dispersion



subroutine test_e(error, mol, disp, qat, ref, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Instance of dispersion type
   class(dispersion_type), intent(in) :: disp

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Reference energy
   real(wp), intent(in) :: ref

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   type(wavefunction_type) :: wfn
   type(container_cache) :: cache
   real(wp) :: energies(mol%nat)
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   energies = 0.0_wp
   wfn%qat = reshape(qat, [size(qat), 1])
   call disp%update(mol, cache)

   ! Get non-selfconsistent energy contribution
   call disp%get_engrad(mol, cache, energies)

   ! Get selfconsistent energy contribution
   call disp%get_energy(mol, cache, wfn, energies)

   if (abs(sum(energies) - ref) > thr_) then
      call test_failed(error, "Energy does not match reference")
      print *, sum(energies), ref, abs(sum(energies) - ref)
   end if

end subroutine test_e


subroutine test_p(error, mol, disp, qat, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Instance of dispersion type
   class(dispersion_type), intent(in) :: disp

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   type(container_cache) :: cache
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp), allocatable :: vat(:)
   real(wp) :: er(mol%nat), el(mol%nat)
   integer :: ii
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   wfn%qat = reshape(qat, [size(qat), 1])
   allocate(pot%vat(size(qat, 1), 1))
   call disp%update(mol, cache)

   allocate(vat(mol%nat), source=0.0_wp)
   do ii = 1, mol%nat
      er = 0.0_wp
      el = 0.0_wp
      wfn%qat(ii, 1) = wfn%qat(ii, 1) + step
      call disp%get_energy(mol, cache, wfn, er)

      wfn%qat(ii, 1) = wfn%qat(ii, 1) - 2*step
      call disp%get_energy(mol, cache, wfn, el)

      wfn%qat(ii, 1) = wfn%qat(ii, 1) + step
      vat(ii) = 0.5_wp*(sum(er) - sum(el))/step
   end do

   pot%vat(:, :) = 0.0_wp
   call disp%get_potential(mol, cache, wfn, pot)

   if (any(abs([pot%vat] - vat) > thr_)) then
      call test_failed(error, "Potential does not match")
      print '(3es20.13)', pot%vat
      print '(a)', "---"
      print '(3es20.13)', vat
      print '(a)', "---"
      print '(3es20.13)', [pot%vat] - vat
   end if
end subroutine test_p


subroutine test_e_d3_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_type), allocatable :: disp
   type(d3_dispersion), allocatable :: d3

   real(wp), parameter :: qat(16) = [&
      & 4.07331798531438E+0_wp, 2.56451650429167E+0_wp, 2.97676954448882E+0_wp, &
      & 2.88833112759592E+0_wp, 2.59127476011008E+0_wp, 2.63750279425510E+0_wp, &
      & 3.56149571036025E+0_wp, 2.89090958281373E+0_wp, 4.53815592283277E+0_wp, &
      & 2.46342847720303E+0_wp, 2.72707461251522E+0_wp, 3.52933932564532E+0_wp, &
      & 3.66934919146868E+0_wp, 3.66697827876019E+0_wp, 3.37809284756764E+0_wp, &
      & 4.57932013403544E+0_wp]

   call get_structure(mol, "MB16-43", "01")

   ! PBE-D3-ATM
   allocate(d3)
   call new_d3_dispersion(d3, mol, s6=1.0_wp, s8=0.7875_wp, &
      & a1=0.4289_wp, a2=4.4407_wp, s9=1.0_wp, error=error)
   if(allocated(error)) then 
      call test_failed(error, "D3 model could not be created")
      return
   end if
   call move_alloc(d3, disp)

   call test_e(error, mol, disp, qat, -1.7817346201130551E-02_wp, thr_in=thr)

end subroutine test_e_d3_mb01

subroutine test_e_d4_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_type), allocatable :: disp
   type(d4_dispersion), allocatable :: d4

   real(wp), parameter :: qat(16) = [&
      & 4.07331798531438E+0_wp, 2.56451650429167E+0_wp, 2.97676954448882E+0_wp, &
      & 2.88833112759592E+0_wp, 2.59127476011008E+0_wp, 2.63750279425510E+0_wp, &
      & 3.56149571036025E+0_wp, 2.89090958281373E+0_wp, 4.53815592283277E+0_wp, &
      & 2.46342847720303E+0_wp, 2.72707461251522E+0_wp, 3.52933932564532E+0_wp, &
      & 3.66934919146868E+0_wp, 3.66697827876019E+0_wp, 3.37809284756764E+0_wp, &
      & 4.57932013403544E+0_wp]

   call get_structure(mol, "MB16-43", "01")

   ! PBE-D4-ATM
   allocate(d4)
   call new_d4_dispersion(d4, mol, s6=1.0_wp, s8=0.95948085_wp, &
      & a1=0.38574991_wp, a2=4.80688534_wp, s9=1.0_wp, error=error)
   if(allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call move_alloc(d4, disp)

   call test_e(error, mol, disp, qat, -2.1345348473924487E-03_wp, thr_in=thr)

end subroutine test_e_d4_mb01

subroutine test_e_d4s_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_type), allocatable :: disp
   type(d4_dispersion), allocatable :: d4s

   real(wp), parameter :: qat(16) = [&
      & 4.07331798531438E+0_wp, 2.56451650429167E+0_wp, 2.97676954448882E+0_wp, &
      & 2.88833112759592E+0_wp, 2.59127476011008E+0_wp, 2.63750279425510E+0_wp, &
      & 3.56149571036025E+0_wp, 2.89090958281373E+0_wp, 4.53815592283277E+0_wp, &
      & 2.46342847720303E+0_wp, 2.72707461251522E+0_wp, 3.52933932564532E+0_wp, &
      & 3.66934919146868E+0_wp, 3.66697827876019E+0_wp, 3.37809284756764E+0_wp, &
      & 4.57932013403544E+0_wp]

   call get_structure(mol, "MB16-43", "01")

   ! PBE-D4S-ATM
   allocate(d4s)
   call new_d4s_dispersion(d4s, mol, s6=1.0_wp, s8=0.95948085_wp, &
      & a1=0.38574991_wp, a2=4.80688534_wp, s9=1.0_wp, error=error)
   if(allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call move_alloc(d4s, disp)

   call test_e(error, mol, disp, qat, -2.5035710126985792E-03_wp, thr_in=thr)

end subroutine test_e_d4s_mb01

subroutine test_e_d4srev_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_type), allocatable :: disp
   type(d4_dispersion), allocatable :: d4srev

   real(wp), parameter :: qat(16) = [&
      &  9.115222032717870E-1_wp, -7.820978364181719E-2_wp, -8.085043301444808E-1_wp, &
      & -1.370008894958108E-1_wp, -4.522670977383800E-1_wp,  2.366794200051805E-1_wp, &
      & -1.056480024872486E-1_wp, -7.542974658029240E-1_wp, -6.084417930861774E-1_wp, &
      &  3.252603520948056E-1_wp,  2.590209057852254E-1_wp, -3.245713596482531E-1_wp, &
      &  9.588528411034369E-2_wp,  7.104890489542336E-1_wp, -3.468331265256146E-1_wp, &
      &  1.076916633846320E+0_wp]

   call get_structure(mol, "MB16-43", "01")

   ! g-xTB D4Srev-ATM
   allocate(d4srev)
   call new_d4srev_dispersion(d4srev, mol, twobody_damping_function%screened, &
      & threebody_damping_function%screened, s6=1.0_wp, s8=0.0_wp, &
      & a1=0.0_wp, a2=8.5302496868_wp, a3=0.5395219242_wp, &
      & a4=0.65400000000_wp, s9=1.8264106910_wp, kcn=7.5_wp, error=error)
   if(allocated(error)) then 
      call test_failed(error, "D4srev model could not be created")
      return
   end if
   call move_alloc(d4srev, disp)

   call test_e(error, mol, disp, qat, -1.804914398629546E-02_wp, thr_in=thr)

end subroutine test_e_d4srev_mb01

subroutine test_e_d3_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_type), allocatable :: disp
   type(d3_dispersion), allocatable :: d3

   real(wp), parameter :: qat(16) = [&
      & 3.43501192886252E+0_wp, 5.02404916999371E+0_wp, 4.72253865039692E+0_wp, &
      & 3.52096661104217E+0_wp, 4.76023330956437E+0_wp, 3.46119195863261E+0_wp, &
      & 3.17361370475619E+0_wp, 2.90775065608382E+0_wp, 4.94595287355805E+0_wp, &
      & 3.32592657749444E+0_wp, 4.54348353409109E+0_wp, 4.30924297002105E+0_wp, &
      & 3.47420343563851E+0_wp, 2.82302349370343E+0_wp, 6.67552064739394E+0_wp, &
      & 4.23898159491675E+0_wp]

   call get_structure(mol, "MB16-43", "02")

   ! PBE-D3-ATM
   allocate(d3)
   call new_d3_dispersion(d3, mol, s6=1.0_wp, s8=0.7875_wp, &
      & a1=0.4289_wp, a2=4.4407_wp, s9=1.0_wp, error=error)
   if(allocated(error)) then 
      call test_failed(error, "D3 model could not be created")
      return
   end if
   call move_alloc(d3, disp)

   call test_e(error, mol, disp, qat, -2.4863195815287387E-02_wp, thr_in=thr)

end subroutine test_e_d3_mb02

subroutine test_e_d4_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_type), allocatable :: disp
   type(d4_dispersion), allocatable :: d4

   real(wp), parameter :: qat(16) = [&
      & 3.43501192886252E+0_wp, 5.02404916999371E+0_wp, 4.72253865039692E+0_wp, &
      & 3.52096661104217E+0_wp, 4.76023330956437E+0_wp, 3.46119195863261E+0_wp, &
      & 3.17361370475619E+0_wp, 2.90775065608382E+0_wp, 4.94595287355805E+0_wp, &
      & 3.32592657749444E+0_wp, 4.54348353409109E+0_wp, 4.30924297002105E+0_wp, &
      & 3.47420343563851E+0_wp, 2.82302349370343E+0_wp, 6.67552064739394E+0_wp, &
      & 4.23898159491675E+0_wp]

   call get_structure(mol, "MB16-43", "02")

   ! PBE-D4-ATM
   allocate(d4)
   call new_d4_dispersion(d4, mol, s6=1.0_wp, s8=0.95948085_wp, &
      & a1=0.38574991_wp, a2=4.80688534_wp, s9=1.0_wp, error=error)
   if(allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call move_alloc(d4, disp)

   call test_e(error, mol, disp, qat, -3.2658253557002496E-03_wp, thr_in=thr)

end subroutine test_e_d4_mb02

subroutine test_e_d4s_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_type), allocatable :: disp
   type(d4_dispersion), allocatable :: d4s

   real(wp), parameter :: qat(16) = [&
      & 3.43501192886252E+0_wp, 5.02404916999371E+0_wp, 4.72253865039692E+0_wp, &
      & 3.52096661104217E+0_wp, 4.76023330956437E+0_wp, 3.46119195863261E+0_wp, &
      & 3.17361370475619E+0_wp, 2.90775065608382E+0_wp, 4.94595287355805E+0_wp, &
      & 3.32592657749444E+0_wp, 4.54348353409109E+0_wp, 4.30924297002105E+0_wp, &
      & 3.47420343563851E+0_wp, 2.82302349370343E+0_wp, 6.67552064739394E+0_wp, &
      & 4.23898159491675E+0_wp]

   call get_structure(mol, "MB16-43", "02")

   ! PBE-D4S-ATM
   allocate(d4s)
   call new_d4s_dispersion(d4s, mol, s6=1.0_wp, s8=0.95948085_wp, &
      & a1=0.38574991_wp, a2=4.80688534_wp, s9=1.0_wp, error=error)
   if(allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call move_alloc(d4s, disp)

   call test_e(error, mol, disp, qat, -4.0744895782408065E-03_wp, thr_in=thr)

end subroutine test_e_d4s_mb02

subroutine test_e_d4srev_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_type), allocatable :: disp
   type(d4_dispersion), allocatable :: d4srev

   real(wp), parameter :: qat(16) = [&
      & -1.097288775394150E-1_wp, -4.735023909155970E-1_wp, -1.707389583649015E-1_wp, &
      & -6.433491759460952E-1_wp,  8.965254778844154E-1_wp,  3.351267282573084E-1_wp, &
      & -6.657210518576462E-2_wp, -5.047599279669956E-2_wp,  5.146699702222322E-1_wp, &
      &  3.455270064324552E-1_wp,  8.365308656920356E-2_wp,  7.840093894994580E-1_wp, &
      & -6.236807887791764E-1_wp, -7.429581513991601E-2_wp, -6.465871716055793E-2_wp, &
      & -6.825088374693600E-1_wp]

   call get_structure(mol, "MB16-43", "02")

   ! g-xTB D4Srev-ATM
   allocate(d4srev)
   call new_d4srev_dispersion(d4srev, mol, twobody_damping_function%screened, &
      & threebody_damping_function%screened, s6=1.0_wp, s8=0.0_wp, &
      & a1=0.0_wp, a2=8.5302496868_wp, a3=0.5395219242_wp, &
      & a4=0.65400000000_wp, s9=1.8264106910_wp, kcn=7.5_wp, error=error)
   if(allocated(error)) then 
      call test_failed(error, "D4srev model could not be created")
      return
   end if
   call move_alloc(d4srev, disp)

   call test_e(error, mol, disp, qat, -1.929815348508024E-02_wp, thr_in=thr)

end subroutine test_e_d4srev_mb02


subroutine test_p_d4_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_type), allocatable :: disp
   type(d4_dispersion), allocatable :: d4

   real(wp), parameter :: qat(16) = [&
      & 4.94764986698701E+0_wp, 3.95122747438791E+0_wp, 4.57075556245289E+0_wp, &
      & 5.46368225070994E+0_wp, 8.24269261139398E+0_wp, 5.68405471762112E+0_wp, &
      & 5.51002325309604E+0_wp, 4.75597020148093E+0_wp, 4.21190195089894E+0_wp, &
      & 4.32836770082885E+0_wp, 4.12684869499911E+0_wp, 5.15171226623248E+0_wp, &
      & 4.83223856996055E+0_wp, 3.02638025720185E+0_wp, 4.05683426506167E+0_wp, &
      & 4.63569783992096E+0_wp]

   call get_structure(mol, "MB16-43", "03")

   ! PBE-D4-ATM
   allocate(d4)
   call new_d4_dispersion(d4, mol, s6=1.0_wp, s8=0.95948085_wp, &
      & a1=0.38574991_wp, a2=4.80688534_wp, s9=1.0_wp, error=error)
   if(allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call move_alloc(d4, disp)

   call test_p(error, mol, disp, qat, thr_in=thr1)

end subroutine test_p_d4_mb03

subroutine test_p_d4s_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_type), allocatable :: disp
   type(d4_dispersion), allocatable :: d4s

   real(wp), parameter :: qat(16) = [&
      & 4.94764986698701E+0_wp, 3.95122747438791E+0_wp, 4.57075556245289E+0_wp, &
      & 5.46368225070994E+0_wp, 8.24269261139398E+0_wp, 5.68405471762112E+0_wp, &
      & 5.51002325309604E+0_wp, 4.75597020148093E+0_wp, 4.21190195089894E+0_wp, &
      & 4.32836770082885E+0_wp, 4.12684869499911E+0_wp, 5.15171226623248E+0_wp, &
      & 4.83223856996055E+0_wp, 3.02638025720185E+0_wp, 4.05683426506167E+0_wp, &
      & 4.63569783992096E+0_wp]

   call get_structure(mol, "MB16-43", "03")
   
   ! PBE-D4S-ATM
   allocate(d4s)
   call new_d4s_dispersion(d4s, mol, s6=1.0_wp, s8=0.95948085_wp, &
      & a1=0.38574991_wp, a2=4.80688534_wp, s9=1.0_wp, error=error)
   if(allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call move_alloc(d4s, disp)

   call test_p(error, mol, disp, qat, thr_in=thr1)

end subroutine test_p_d4s_mb03

subroutine test_p_d4srev_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_type), allocatable :: disp
   type(d4_dispersion), allocatable :: d4srev

   real(wp), parameter :: qat(16) = [&
      & -6.722325895098202E-2_wp, -9.343346477695722E-1_wp,  7.423184694799456E-2_wp, &
      &  6.003342781579326E-1_wp,  8.155559152718005E-1_wp,  1.119087497546364E-1_wp, &
      & -4.000789327736840E-1_wp,  8.016756865222574E-2_wp,  8.005299562801238E-2_wp, &
      &  1.006389816802352E-1_wp, -6.704006636400752E-1_wp, -4.326713538195071E-1_wp, &
      & -5.263477368917056E-1_wp,  2.692434056872415E-1_wp,  8.190366482982609E-1_wp, &
      &  7.988620327393425E-2_wp]

   call get_structure(mol, "MB16-43", "03")

   ! g-xTB D4Srev-ATM
   allocate(d4srev)
   call new_d4srev_dispersion(d4srev, mol, twobody_damping_function%screened, &
      & threebody_damping_function%screened, s6=1.0_wp, s8=0.0_wp, &
      & a1=0.0_wp, a2=8.5302496868_wp, a3=0.5395219242_wp, &
      & a4=0.65400000000_wp, s9=1.8264106910_wp, kcn=7.5_wp, error=error)
   if(allocated(error)) then 
      call test_failed(error, "D4srev model could not be created")
      return
   end if
   call move_alloc(d4srev, disp)

   call test_p(error, mol, disp, qat, thr_in=thr1)

end subroutine test_p_d4srev_mb03

subroutine test_p_d4_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_type), allocatable :: disp
   type(d4_dispersion), allocatable :: d4

   real(wp), parameter :: qat(16) = [&
      &-8.99890404486076E-2_wp, 9.42168087556583E-2_wp,-1.49387631509499E-1_wp, &
      &-2.99114121895542E-1_wp, 4.85527734875224E-1_wp,-6.83156326406137E-2_wp, &
      & 1.50011293889337E-2_wp, 2.79368544459465E-1_wp,-1.24072452878322E-1_wp, &
      &-9.36760994051244E-2_wp,-2.19062123031622E-1_wp, 2.14538817685587E-1_wp, &
      & 3.06156726072831E-1_wp,-3.86105514712244E-1_wp,-1.51265171389388E-3_wp, &
      & 3.64255069977693E-2_wp]

   call get_structure(mol, "MB16-43", "04")
   
   ! PBE-D4-ATM
   allocate(d4)
   call new_d4_dispersion(d4, mol, s6=1.0_wp, s8=0.95948085_wp, &
      & a1=0.38574991_wp, a2=4.80688534_wp, s9=1.0_wp, error=error)
   if(allocated(error)) then 
      call test_failed(error, "D4 model could not be created")
      return
   end if
   call move_alloc(d4, disp)

   call test_p(error, mol, disp, qat, thr_in=thr1)

end subroutine test_p_d4_mb04

subroutine test_p_d4s_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_type), allocatable :: disp
   type(d4_dispersion), allocatable :: d4s

   real(wp), parameter :: qat(16) = [&
      &-8.99890404486076E-2_wp, 9.42168087556583E-2_wp,-1.49387631509499E-1_wp, &
      &-2.99114121895542E-1_wp, 4.85527734875224E-1_wp,-6.83156326406137E-2_wp, &
      & 1.50011293889337E-2_wp, 2.79368544459465E-1_wp,-1.24072452878322E-1_wp, &
      &-9.36760994051244E-2_wp,-2.19062123031622E-1_wp, 2.14538817685587E-1_wp, &
      & 3.06156726072831E-1_wp,-3.86105514712244E-1_wp,-1.51265171389388E-3_wp, &
      & 3.64255069977693E-2_wp]

   call get_structure(mol, "MB16-43", "04")

   ! PBE-D4S-ATM
   allocate(d4s)
   call new_d4s_dispersion(d4s, mol, s6=1.0_wp, s8=0.95948085_wp, &
      & a1=0.38574991_wp, a2=4.80688534_wp, s9=1.0_wp, error=error)
   if(allocated(error)) then 
      call test_failed(error, "D4S model could not be created")
      return
   end if
   call move_alloc(d4s, disp)

   call test_p(error, mol, disp, qat, thr_in=thr1)

end subroutine test_p_d4s_mb04


subroutine test_p_d4srev_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(dispersion_type), allocatable :: disp
   type(d4_dispersion), allocatable :: d4srev

   real(wp), parameter :: qat(16) = [&
      & -6.417476509346409E-2_wp, -2.256661931683622E-1_wp, -3.458626714493751E-2_wp, &
      & -3.860600924041810E-1_wp,  7.238991472693167E-1_wp, -5.777084491650908E-2_wp, &
      &  1.196120411671119E-1_wp, -5.461117600438992E-3_wp, -1.794514383018986E-2_wp, &
      & -1.211514276086560E-1_wp, -3.493736689115619E-1_wp,  9.244684713156893E-1_wp, &
      &  1.394030055851094E-1_wp, -6.975412178494735E-1_wp,  1.149356570285273E-1_wp, &
      & -6.258758429064248E-2_wp]

   call get_structure(mol, "MB16-43", "04")

   ! g-xTB D4Srev-ATM
   allocate(d4srev)
   call new_d4srev_dispersion(d4srev, mol, twobody_damping_function%screened, &
      & threebody_damping_function%screened, s6=1.0_wp, s8=0.0_wp, &
      & a1=0.0_wp, a2=8.5302496868_wp, a3=0.5395219242_wp, &
      & a4=0.65400000000_wp, s9=1.8264106910_wp, kcn=7.5_wp, error=error)
   if(allocated(error)) then 
      call test_failed(error, "D4srev model could not be created")
      return
   end if
   call move_alloc(d4srev, disp)

   call test_p(error, mol, disp, qat, thr_in=thr1)

end subroutine test_p_d4srev_mb04


end module test_dispersion
