module test_halogen
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io_structure, only : structure_type, new
   use mstore, only : get_structure
   use tblite_container, only : container_cache
   use tblite_classical_halogen
   implicit none
   private

   public :: collect_halogen

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   !> Strength of the halogen bond
   real(wp), parameter :: halogen_bond(*) = 0.1_wp * [ &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.381742_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.321944_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.220000_wp, &
      & 0.000000_wp]

contains


!> Collect all exported unit tests
subroutine collect_halogen(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("Br2-NH3", test_br2nh3), &
      & new_unittest("Br2-OCN2", test_br2och2), &
      & new_unittest("FI-NCH", test_finch) &
      & ]

end subroutine collect_halogen


subroutine test_br2nh3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(halogen_correction) :: xb
   type(container_cache) :: cache
   real(wp), allocatable :: energy(:)
   real(wp), parameter :: trans(3, 1) = 0.0_wp, cutoff = sqrt(400.0_wp)

   call BrBr_NH3(mol)
   call new_halogen_correction(xb, mol, 0.44_wp, 1.3_wp, halogen_bond(mol%num))

   allocate(energy(mol%nat))

   call check(error, count(xb%halogen), 1)
   if (allocated(error)) return

   call check(error, count(xb%acceptor), 1)
   if (allocated(error)) return

   energy = 0.0_wp
   call xb%get_engrad(mol, cache, energy)

   call check(error, sum(energy), 2.4763110097465683E-3_wp, thr=thr)
   if (allocated(error)) return

end subroutine test_br2nh3


subroutine test_br2och2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(halogen_correction) :: xb
   type(container_cache) :: cache
   integer :: iat, ic
   real(wp) :: sigma(3, 3)
   real(wp), allocatable :: energy(:), er(:), el(:)
   real(wp), parameter :: trans(3, 1) = 0.0_wp, cutoff = sqrt(400.0_wp)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp), parameter :: step = 1.0e-6_wp

   call BrBr_OCH2(mol)
   call new_halogen_correction(xb, mol, 0.44_wp, 1.3_wp, halogen_bond(mol%num))

   allocate(energy(mol%nat), er(mol%nat), el(mol%nat))

   call check(error, count(xb%halogen), 1)
   if (allocated(error)) return

   call check(error, count(xb%acceptor), 1)
   if (allocated(error)) return

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call xb%get_engrad(mol, cache, energy, gradient, sigma)

   call check(error, sum(energy), -6.7587305781592112E-4_wp, thr=thr)
   if (allocated(error)) return

   do iat = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call xb%get_engrad(mol, cache, er)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call xb%get_engrad(mol, cache, el)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   if (any(abs(gradient - numgrad) > thr2)) then
      call test_failed(error, "Gradient of dispersion energy does not match")
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_br2och2


subroutine test_finch(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(halogen_correction) :: xb
   type(container_cache) :: cache
   integer :: ic, jc
   real(wp) :: sigma(3, 3), eps(3, 3), numsigma(3, 3)
   real(wp), allocatable :: energy(:), er(:), el(:)
   real(wp), parameter :: trans(3, 1) = 0.0_wp, cutoff = sqrt(400.0_wp)
   real(wp), allocatable :: gradient(:, :), xyz(:, :)
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))

   call FI_NCH(mol)
   call new_halogen_correction(xb, mol, 0.44_wp, 1.3_wp, halogen_bond(mol%num))

   allocate(energy(mol%nat), er(mol%nat), el(mol%nat))

   call check(error, count(xb%halogen), 1)
   if (allocated(error)) return

   call check(error, count(xb%acceptor), 1)
   if (allocated(error)) return

   allocate(gradient(3, mol%nat), xyz(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call xb%get_engrad(mol, cache, energy, gradient, sigma)

   call check(error, sum(energy), 1.1857937381795408E-2_wp, thr=thr)
   if (allocated(error)) return

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   do ic = 1, 3
      do jc = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call xb%get_engrad(mol, cache, er)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call xb%get_engrad(mol, cache, el)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         numsigma(jc, ic) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   if (any(abs(sigma - numsigma) > thr2)) then
      call test_failed(error, "Strain derivatives do not match")
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_finch


subroutine BrBr_NH3(self)
   type(structure_type), intent(out) :: self
   integer, parameter :: nat = 6
   character(len=*), parameter :: sym(nat) = [character(len=4)::&
      & "Br", "Br", "N", "H", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      &  0.00000000000000_wp,  0.00000000000000_wp,  3.11495251300000_wp, &
      &  0.00000000000000_wp,  0.00000000000000_wp, -1.25671880600000_wp, &
      &  0.00000000000000_wp,  0.00000000000000_wp, -6.30201130100000_wp, &
      &  0.00000000000000_wp,  1.78712709700000_wp, -6.97470840000000_wp, &
      & -1.54769692500000_wp, -0.89356260400000_wp, -6.97470840000000_wp, &
      &  1.54769692500000_wp, -0.89356260400000_wp, -6.97470840000000_wp],&
      & shape(xyz))
   call new(self, sym, xyz)
end subroutine BrBr_NH3

subroutine BrBr_OCH2(self)
   type(structure_type), intent(out) :: self
   integer, parameter :: nat = 6
   character(len=*), parameter :: sym(nat) = [character(len=4)::&
      & "Br", "Br", "O", "C", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & -1.78533374700000_wp, -3.12608299900000_wp,  0.00000000000000_wp, &
      &  0.00000000000000_wp,  0.81604226400000_wp,  0.00000000000000_wp, &
      &  2.65828699900000_wp,  5.29707580600000_wp,  0.00000000000000_wp, &
      &  4.88597158600000_wp,  4.86116137300000_wp,  0.00000000000000_wp, &
      &  5.61550975300000_wp,  2.90822215900000_wp,  0.00000000000000_wp, &
      &  6.28907612600000_wp,  6.39963643500000_wp,  0.00000000000000_wp],&
      & shape(xyz))
   call new(self, sym, xyz)
end subroutine BrBr_OCH2


subroutine FI_NCH(self)
   type(structure_type), intent(out) :: self
   integer, parameter :: nat = 5
   character(len=*), parameter :: sym(nat) = [character(len=4)::&
      & "F", "I", "N", "C", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      &  0.00000000000000_wp,  0.00000000000000_wp,  4.37637862700000_wp, &
      &  0.00000000000000_wp,  0.00000000000000_wp,  0.69981844700000_wp, &
      &  0.00000000000000_wp,  0.00000000000000_wp, -4.24181123900000_wp, &
      &  0.00000000000000_wp,  0.00000000000000_wp, -6.39520691700000_wp, &
      &  0.00000000000000_wp,  0.00000000000000_wp, -8.41387269200000_wp],&
      & shape(xyz))
   call new(self, sym, xyz)
end subroutine FI_NCH


end module test_halogen
