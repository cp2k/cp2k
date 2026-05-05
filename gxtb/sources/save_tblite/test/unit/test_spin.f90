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

module test_spin
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_type, container_cache
   use tblite_context_type, only : context_type
   use tblite_data_spin, only : get_spin_constant
   use tblite_spin, only : spin_polarization, new_spin_polarization
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_spin

   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 3.166808578545117e-06_wp

contains


!> Collect all exported unit tests
subroutine collect_spin(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gfn1-e-spin", test_e_crcp2), &
      new_unittest("gfn2-e-spin", test_e_p10), &
      new_unittest("gfn1-g-spin", test_g_p10), &
      new_unittest("gfn2-g-spin", test_g_crcp2) &
      ]

end subroutine collect_spin


subroutine test_e_gen(error, mol, calc, ref1, ref0)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure 
   type(structure_type), intent(inout) :: mol
   !> XTB calculator
   type(xtb_calculator), intent(inout) :: calc
   !> Reference energy for spin-polarized UHF
   real(wp), intent(in) :: ref1
   !> Reference energy without spin-polarization
   real(wp), intent(in) :: ref0

   type(context_type) :: ctx
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   integer :: tmpuhf
   real(wp), allocatable :: wll(:, :, :)
   real(wp) :: energy

   energy = 0.0_wp

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 2, calc%default_etemp * kt)

   ! Add spin polarization
   if (.not. allocated(calc%spin_polarization)) then
      call get_spin_constants(mol, calc%bas, wll)
      allocate(calc%spin_polarization)
      call new_spin_polarization(calc%spin_polarization, mol, wll, calc%bas%nsh_id)
   end if

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)
   call check(error, energy, ref1, thr=thr1)
   if (allocated(error)) return

   ! Check low-spin
   if(mod(mol%uhf, 2) == 0) then
      tmpuhf = mol%uhf
      mol%uhf = 0
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
         & 2, calc%default_etemp * kt)
   
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)
      call check(error, energy, ref0, thr=thr1)
      if (allocated(error)) return
      mol%uhf = tmpuhf
   end if
   
   ! Alternatively, add the spin-polarization to the interaction list
   deallocate(calc%spin_polarization)
   block
      type(spin_polarization), allocatable :: spin
      allocate(spin)
      call new_spin_polarization(spin, mol, wll, calc%bas%nsh_id)
      call move_alloc(spin, cont)
      call calc%push_back(cont)
   end block

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)
   call check(error, energy, ref1, thr=thr1)
   if (allocated(error)) return

   call calc%pop(cont)
   if(mod(mol%uhf, 2) == 0) then
      tmpuhf = mol%uhf
      mol%uhf = 0

      call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)
      call check(error, energy, ref0, thr=thr1)
      mol%uhf = tmpuhf
   end if


end subroutine test_e_gen


subroutine test_g_gen(error, mol, calc, eref, gref, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure 
   type(structure_type), intent(inout) :: mol
   !> XTB calculator
   type(xtb_calculator), intent(inout) :: calc
   !> Reference energy for spin-polarized UHF
   real(wp), intent(in) :: eref
   !> Reference gradient for spin-polarized UHF
   real(wp), intent(in) :: gref(:, :)
   !> Test threshold
   real(wp), intent(in), optional :: thr_in
 
   type(context_type) :: ctx
   type(wavefunction_type) :: wfn, wfni
   type(structure_type) :: moli
   class(container_type), allocatable :: cont
   integer :: iat, ic
   real(wp) :: energy, sigma(3, 3), er, el, thr_
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp), allocatable :: wll(:, :, :)
   ! Step size is square root of SCF convergence threshold
   ! (SCF energy convergence threshold is 1e-6 * acc)
   real(wp), parameter :: step = 1.0e-4_wp

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   energy = 0.0_wp
   numgrad(:, :) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 2, calc%default_etemp * kt)

   ! Add spin polarization
   call get_spin_constants(mol, calc%bas, wll)
   allocate(calc%spin_polarization)
   call new_spin_polarization(calc%spin_polarization, mol, wll, calc%bas%nsh_id)

   ! ! Numerical gradient
   ! do iat = 1, mol%nat
   !    do ic = 1, 3
   !       er = 0.0_wp
   !       el = 0.0_wp
   !       ! Right-hand side
   !       moli = mol
   !       wfni = wfn
   !       moli%xyz(ic, iat) = mol%xyz(ic, iat) + step
   !       call xtb_singlepoint(ctx, moli, calc, wfni, acc, er, verbosity=0)

   !       ! Left-hand side
   !       moli = mol
   !       wfni = wfn
   !       moli%xyz(ic, iat) = mol%xyz(ic, iat) - step
   !       call xtb_singlepoint(ctx, moli, calc, wfni, acc, el, verbosity=0)

   !       numgrad(ic, iat) = 0.5_wp*(er - el)/step
   !    end do
   ! end do

   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, verbosity=0)

   call check(error, energy, eref, thr=thr)
   if (allocated(error)) return

   if (any(abs(gradient - gref) > thr_)) then
      call test_failed(error, "Gradient of energy does not match")
      write(*,*) "reference gradient:"
      print'(3es21.14)', gref
      write(*,*) "analytical gradient:"
      print'(3es21.14)', gradient
      write(*,*) "difference:"
      print'(3es21.14)', gradient-gref
   end if

end subroutine test_g_gen


subroutine test_e_p10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   real(wp), parameter :: ref1 = -10.801224611773206_wp, ref0 = -10.7897113668574_wp

   call rse43_p10(mol)

   call new_gfn2_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   call test_e_gen(error, mol, calc, ref1, ref0)

end subroutine test_e_p10

subroutine test_e_crcp2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   real(wp), parameter :: ref1 = -28.376485853038645_wp, ref0 = -28.3496138337329_wp

   call crcp2(mol)
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp

   call new_gfn1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   call test_e_gen(error, mol, calc, ref1, ref0)

end subroutine test_e_crcp2


subroutine test_g_p10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   real(wp), parameter :: eref = -11.539671328635730_wp, gref(3, 8) = reshape([&
      &  4.6250747795231898E-03_wp,  3.0613008404354290E-03_wp, -4.1003763872489694E-17_wp, &
      & -5.8127688918931083E-03_wp,  7.0212976481872583E-03_wp,  7.4006869452289192E-17_wp, &
      &  8.4737687859435182E-03_wp, -7.8529734509620933E-03_wp,  1.4182424493755559E-17_wp, &
      &  1.6732818000371087E-04_wp, -2.6987199782906278E-03_wp,  8.6784398428552196E-18_wp, &
      & -2.4727166634656915E-03_wp,  1.1335882595829142E-03_wp,  1.1165068841281504E-17_wp, &
      & -1.2010615143533908E-03_wp, -5.3284767411932365E-04_wp,  2.1499927548673772E-03_wp, &
      & -1.2010615143534362E-03_wp, -5.3284767411941928E-04_wp, -2.1499927548674475E-03_wp, &
      & -2.5785631614047800E-03_wp,  4.0120202928582558E-04_wp, -2.7346723399901103E-18_wp], &
      & shape(gref))


   call rse43_p10(mol)
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_gfn1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   call test_g_gen(error, mol, calc, eref, gref, thr_in=thr1)

end subroutine test_g_p10


subroutine test_g_crcp2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   real(wp), parameter :: eref = -28.471439648669477_wp, gref(3, 21) = reshape([&
      & -3.3346209907941793E-14_wp,  9.4095128281949936E-15_wp,  7.3031040873146682E-04_wp, &
      &  2.0737422989057993E-14_wp, -1.3810438730716156E-03_wp,  1.4534865767653804E-02_wp, &
      &  1.2023217179885587E-02_wp,  8.5569295722914113E-05_wp,  9.6517319007189957E-04_wp, &
      &  3.7815059360380221E-05_wp,  1.8331263028967730E-03_wp, -7.4190049135553562E-03_wp, &
      & -3.7815059340987442E-05_wp,  1.8331263028886257E-03_wp, -7.4190049135963677E-03_wp, &
      & -1.2023217179908601E-02_wp,  8.5569295703461453E-05_wp,  9.6517319009760687E-04_wp, &
      & -2.1579772050526842E-15_wp, -4.7641226402683678E-03_wp,  1.7211385778858280E-03_wp, &
      &  1.8164603206477625E-03_wp, -4.7321956256350887E-03_wp, -1.2363412674730632E-04_wp, &
      &  4.6032135883547536E-04_wp, -4.5828980598753965E-03_wp, -1.7331139247164979E-03_wp, &
      & -4.6032135883693128E-04_wp, -4.5828980598751675E-03_wp, -1.7331139247122576E-03_wp, &
      & -1.8164603206444060E-03_wp, -4.7321956256338536E-03_wp, -1.2363412674987575E-04_wp, &
      &  1.2023217179875622E-02_wp, -8.5569295726056334E-05_wp,  9.6517319007691270E-04_wp, &
      &  2.1969354427773201E-14_wp,  1.3810438730645139E-03_wp,  1.4534865767642511E-02_wp, &
      &  3.7815059363820869E-05_wp, -1.8331263028954043E-03_wp, -7.4190049135521609E-03_wp, &
      &  1.8164603206502138E-03_wp,  4.7321956256355432E-03_wp, -1.2363412674731755E-04_wp, &
      & -1.2023217179899397E-02_wp, -8.5569295707092256E-05_wp,  9.6517319009984900E-04_wp, &
      & -2.2096020250407485E-15_wp,  4.7641226402693401E-03_wp,  1.7211385778886621E-03_wp, &
      & -3.7815059345164216E-05_wp, -1.8331263028875707E-03_wp, -7.4190049135908296E-03_wp, &
      &  4.6032135883590795E-04_wp,  4.5828980598754208E-03_wp, -1.7331139247177656E-03_wp, &
      & -1.8164603206467446E-03_wp,  4.7321956256343913E-03_wp, -1.2363412674968964E-04_wp, &
      & -4.6032135883749967E-04_wp,  4.5828980598752534E-03_wp, -1.7331139247138050E-03_wp], &
      & shape(gref))

   call crcp2(mol)
   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_gfn2_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return

   call test_g_gen(error, mol, calc, eref, gref, thr_in=thr1)

end subroutine test_g_crcp2


subroutine get_spin_constants(mol, bas, wll)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set data
   type(basis_type), intent(in) :: bas
   !> Spin constants
   real(wp), allocatable, intent(out) :: wll(:, :, :)

   integer :: izp, ish, jsh, il, jl

   allocate(wll(bas%nsh, bas%nsh, mol%nid), source=0.0_wp)

   do izp = 1, mol%nid
      do ish = 1, bas%nsh_id(izp)
         il = bas%cgto(ish, izp)%raw%ang
         do jsh = 1, bas%nsh_id(izp)
            jl = bas%cgto(jsh, izp)%raw%ang
            wll(jsh, ish, izp) = get_spin_constant(jl, il, mol%num(izp))
         end do
      end do
   end do
end subroutine get_spin_constants


subroutine rse43_p10(self)
   type(structure_type), intent(out) :: self
   integer, parameter :: nat = 8
   character(len=*), parameter :: sym(nat) = [character(len=4)::&
      & "C", "C", "O", "H", "H", "H", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
   & -1.97051959765227E+00_wp,   -8.65723337874754E-01_wp,    0.00000000000000E+00_wp, &     
   &  3.50984622791913E-01_wp,    6.86290619844032E-01_wp,    0.00000000000000E+00_wp, &      
   &  2.50609985217434E+00_wp,   -9.34496149122418E-01_wp,    0.00000000000000E+00_wp, &      
   & -1.83649606109455E+00_wp,   -2.90299181092583E+00_wp,    0.00000000000000E+00_wp, &      
   & -3.80466245712260E+00_wp,    3.49832428602470E-02_wp,    0.00000000000000E+00_wp, &      
   &  3.73555581511497E-01_wp,    1.94431040908594E+00_wp,   -1.66596178649581E+00_wp, &      
   &  3.73555581511497E-01_wp,    1.94431040908594E+00_wp,    1.66596178649581E+00_wp, &      
   &  4.00748247788016E+00_wp,    9.33166170468600E-02_wp,    0.00000000000000E+00_wp], &      
      & shape(xyz))
   integer, parameter :: uhf = 1
   call new(self, sym, xyz, uhf=uhf)
end subroutine rse43_p10


subroutine crcp2(self)
   type(structure_type), intent(out) :: self
   integer, parameter :: nat = 21
   character(len=*), parameter :: sym(nat) = [character(len=4)::&
      & "Cr", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "C", "C", "C", &
      & "H", "C", "H", "C", "H", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
  &  0.00000000000000E+00_wp,    0.00000000000000E+00_wp,   -6.04468452830504E-02_wp, &      
  &  0.00000000000000E+00_wp,    3.19613712523833E+00_wp,    2.30877824528580E+00_wp, &      
  &  2.18828801115897E+00_wp,    3.32943780995850E+00_wp,    7.02499485857345E-01_wp, &      
  &  1.33235791539260E+00_wp,    3.55640652898451E+00_wp,   -1.83908673090077E+00_wp, &      
  & -1.33235791539260E+00_wp,    3.55640652898451E+00_wp,   -1.83908673090077E+00_wp, &      
  & -2.18828801115897E+00_wp,    3.32943780995850E+00_wp,    7.02499485857345E-01_wp, &      
  &  0.00000000000000E+00_wp,    3.10509505378016E+00_wp,    4.34935395653655E+00_wp, &      
  &  4.13810718850644E+00_wp,    3.28428734944129E+00_wp,    1.31235006648465E+00_wp, &      
  &  2.52190264478215E+00_wp,    3.60569548880831E+00_wp,   -3.50208900904436E+00_wp, &      
  & -2.52190264478215E+00_wp,    3.60569548880831E+00_wp,   -3.50208900904436E+00_wp, &      
  & -4.13810718850644E+00_wp,    3.28428734944129E+00_wp,    1.31235006648465E+00_wp, &      
  &  2.18828801115897E+00_wp,   -3.32943780995850E+00_wp,    7.02499485857345E-01_wp, &      
  &  0.00000000000000E+00_wp,   -3.19613712523833E+00_wp,    2.30877824528580E+00_wp, &      
  &  1.33235791539260E+00_wp,   -3.55640652898451E+00_wp,   -1.83908673090077E+00_wp, &      
  &  4.13810718850644E+00_wp,   -3.28428734944129E+00_wp,    1.31235006648465E+00_wp, &      
  & -2.18828801115897E+00_wp,   -3.32943780995850E+00_wp,    7.02499485857345E-01_wp, &      
  &  0.00000000000000E+00_wp,   -3.10509505378016E+00_wp,    4.34935395653655E+00_wp, &      
  & -1.33235791539260E+00_wp,   -3.55640652898451E+00_wp,   -1.83908673090077E+00_wp, &      
  &  2.52190264478215E+00_wp,   -3.60569548880831E+00_wp,   -3.50208900904436E+00_wp, &      
  & -4.13810718850644E+00_wp,   -3.28428734944129E+00_wp,    1.31235006648465E+00_wp, &      
  & -2.52190264478215E+00_wp,   -3.60569548880831E+00_wp,   -3.50208900904436E+00_wp], &      
 & shape(xyz))
   integer, parameter :: uhf = 2
   call new(self, sym, xyz, uhf=uhf)
end subroutine crcp2


end module test_spin
