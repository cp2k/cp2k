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

module test_periodic
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use dftd4
   implicit none
   private

   public :: collect_periodic

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: thr3 = 100*sqrt(epsilon(1.0_wp))
   type(realspace_cutoff), parameter :: cutoff = &
      & realspace_cutoff(cn=30_wp, disp2=60.0_wp, disp3=15.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_periodic(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("PBE-D4", test_pbed4_acetic), &
      & new_unittest("PBE-D4S", test_pbed4s_acetic), &
      & new_unittest("BLYP-D4", test_blypd4_adaman), &
      & new_unittest("BLYP-D4S", test_blypd4s_adaman), &
      & new_unittest("TPSS-D4", test_tpssd4_ammonia), &
      & new_unittest("TPSS-D4S", test_tpssd4s_ammonia), &
      & new_unittest("SCAN-D4", test_scand4_anthracene), &
      & new_unittest("SCAN-D4S", test_scand4s_anthracene) &
      & ]

end subroutine collect_periodic


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

   call get_dispersion(mol, d4, damp, param, cutoff, energy)

   call check(error, energy, ref, thr=thr)
   if (allocated(error)) then
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
         call get_dispersion(mol, d4, damp, param, cutoff, er)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_dispersion(mol, d4, damp, param, cutoff, el)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(er - el)/step
      end do
   end do

   call get_dispersion(mol, d4, damp, param, cutoff, energy, gradient, sigma)

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
   real(wp) :: energy, er, el, sigma(3, 3), eps(3, 3), numsigma(3, 3), lattice(3, 3)
   real(wp), allocatable :: gradient(:, :), xyz(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), parameter :: step = 1.0e-7_wp

   allocate(gradient(3, mol%nat), xyz(3, mol%nat))

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattice(:, :) = mol%lattice
   do ic = 1, 3
      do jc = 1, 3
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         call get_dispersion(mol, d4, damp, param, cutoff, er)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         call get_dispersion(mol, d4, damp, param, cutoff, el)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         mol%lattice(:, :) = lattice
         numsigma(jc, ic) = 0.5_wp*(er - el)/step
      end do
   end do

   call get_dispersion(mol, d4, damp, param, cutoff, energy, gradient, sigma)

   if (any(abs(sigma - numsigma) > thr3)) then
      call test_failed(error, "Strain derivatives do not match")
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_numsigma


subroutine test_pbed4_acetic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "X23", "acetic")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4, damp, param, -6.6969773229895183E-002_wp)

end subroutine test_pbed4_acetic

subroutine test_pbed4s_acetic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 0.95948085_wp, a1 = 0.38574991_wp, a2 = 4.80688534_wp)

   call get_structure(mol, "X23", "acetic")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4s, damp, param, -6.8611894281210548E-002_wp)

end subroutine test_pbed4s_acetic


subroutine test_blypd4_adaman(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 2.34076671_wp, a1 = 0.44488865_wp, a2 = 4.09330090_wp)

   call get_structure(mol, "X23", "adaman")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4, damp, param, -0.23629687693703993_wp)

end subroutine test_blypd4_adaman

subroutine test_blypd4s_adaman(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 2.34076671_wp, a1 = 0.44488865_wp, a2 = 4.09330090_wp)

   call get_structure(mol, "X23", "adaman")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_dftd4_gen(error, mol, d4s, damp, param, -0.25250595410032312_wp)

end subroutine test_blypd4s_adaman


subroutine test_tpssd4_ammonia(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.76596355_wp, a1 = 0.42822303_wp, a2 = 4.54257102_wp)

   call get_structure(mol, "X23", "ammonia")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_numgrad(error, mol, d4, damp, param)

end subroutine test_tpssd4_ammonia

subroutine test_tpssd4s_ammonia(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.76596355_wp, a1 = 0.42822303_wp, a2 = 4.54257102_wp)

   call get_structure(mol, "X23", "ammonia")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_numgrad(error, mol, d4s, damp, param)

end subroutine test_tpssd4s_ammonia


subroutine test_scand4_anthracene(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4_model) :: d4
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.46126056_wp, a1 = 0.62930855_wp, a2 = 6.31284039_wp)

   call get_structure(mol, "X23", "anthracene")
   call new_d4_model(error, d4, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4%default_damping_2b, d4%default_damping_3b)
   if (allocated(error)) return
   call test_numsigma(error, mol, d4, damp, param)

end subroutine test_scand4_anthracene

subroutine test_scand4s_anthracene(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(d4s_model) :: d4s
   type(param_type) :: param
   type(damping_type) :: damp

   param = param_type(&
      & s6 = 1.0_wp, s9 = 0.0_wp, alp = 16.0_wp, &
      & s8 = 1.46126056_wp, a1 = 0.62930855_wp, a2 = 6.31284039_wp)

   call get_structure(mol, "X23", "anthracene")
   call new_d4s_model(error, d4s, mol)
   if (allocated(error)) return
   call new_damping(error, damp, d4s%default_damping_2b, d4s%default_damping_3b)
   if (allocated(error)) return
   call test_numsigma(error, mol, d4s, damp, param)

end subroutine test_scand4s_anthracene

end module test_periodic
