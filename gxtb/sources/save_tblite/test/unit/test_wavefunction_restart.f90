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

module test_wavefunction_restart
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_wavefunction_restart, only : load_wavefunction, save_wavefunction
   use tblite_io_numpy, only : save_npz
   implicit none
   private

   public :: collect_wavefunction_restart

   real(wp), parameter :: thr = epsilon(1.0_wp)

contains


!> Collect all exported unit tests
subroutine collect_wavefunction_restart(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("write-read", test_write_read), &
      new_unittest("write-twice", test_write_twice), &
      new_unittest("file-not-found", test_file_not_found, should_fail=.true.), &
      new_unittest("nat-mismatch", test_nat_mismatch, should_fail=.true.), &
      new_unittest("nsh-mismatch", test_nsh_mismatch, should_fail=.true.), &
      new_unittest("nao-mismatch", test_nao_mismatch, should_fail=.true.), &
      new_unittest("nspin-mismatch", test_nspin_mismatch, should_fail=.true.), &
      new_unittest("kt-mismatch", test_kt_mismatch, should_fail=.true.), &
      new_unittest("dipole-mismatch", test_dipole_mismatch, should_fail=.true.), &
      new_unittest("quadrupole-mismatch", test_quadrupole_mismatch, should_fail=.true.), &
      new_unittest("load-error", test_load_error, should_fail=.true.) &
      ]

end subroutine collect_wavefunction_restart

subroutine test_write_read(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(wavefunction_type) :: wfn1, wfn2
   integer, parameter :: nat = 5, nsh = 12, nao = 20, nspin = 2

   call new_wavefunction(wfn1, nat, nsh, nao, nspin, 0.0_wp)
   call random_number(wfn1%emo)
   call random_number(wfn1%density)

   call save_wavefunction(wfn1, ".test-wavefunction-restart.npz", error)
   if (.not. allocated(error)) then
      call load_wavefunction(wfn2, ".test-wavefunction-restart.npz", error)
      call delete_file(".test-wavefunction-restart.npz")
   end if
   if (allocated(error)) return

   call check(error, all(abs(wfn1%emo - wfn2%emo) < thr), &
      & "Orbital energies are not equal")
   if (allocated(error)) return

   call check(error, all(abs(wfn1%density - wfn2%density) < thr), &
      & "Density matrices are not equal")
end subroutine test_write_read

subroutine test_write_twice(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(wavefunction_type) :: wfn1, wfn2
   integer, parameter :: nat = 5, nsh = 12, nao = 20, nspin = 2

   call new_wavefunction(wfn1, nat, nsh, nao, nspin, 0.0_wp)
   call random_number(wfn1%emo)
   call random_number(wfn1%density)

   call save_wavefunction(wfn1, ".test-wavefunction-restart2.npz", error)
   if (.not. allocated(error)) then
      call save_wavefunction(wfn1, ".test-wavefunction-restart2.npz", error)
   end if
   if (.not. allocated(error)) then
      call load_wavefunction(wfn2, ".test-wavefunction-restart2.npz", error)
      call delete_file(".test-wavefunction-restart2.npz")
   end if
   if (allocated(error)) return

   call check(error, all(abs(wfn1%emo - wfn2%emo) < thr), &
      & "Orbital energies are not equal")
   if (allocated(error)) return

   call check(error, all(abs(wfn1%density - wfn2%density) < thr), &
      & "Density matrices are not equal")
end subroutine test_write_twice

subroutine test_file_not_found(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(wavefunction_type) :: wfn
   character(len=*), parameter :: filename = ".test-file-not-found.npz"

   call load_wavefunction(wfn, filename, error)
end subroutine test_file_not_found

subroutine test_nat_mismatch(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(wavefunction_type) :: wfn
   integer, parameter :: nat = 5, nsh = 12, nao = 20, nspin = 2
   character(len=*), parameter :: prefix = "tblite0_", filename = ".test-nat-mismatch.npz"
   integer :: stat

   call new_wavefunction(wfn, nat, nsh, nao, nspin, 0.0_wp)

   call save_wavefunction(wfn, filename, error)
   call load_wavefunction(wfn, filename, error, nat=10)
   call delete_file(filename)
end subroutine test_nat_mismatch

subroutine test_nsh_mismatch(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(wavefunction_type) :: wfn
   integer, parameter :: nat = 5, nsh = 12, nao = 20, nspin = 2
   character(len=*), parameter :: prefix = "tblite0_", filename = ".test-nsh-mismatch.npz"
   integer :: stat

   call new_wavefunction(wfn, nat, nsh, nao, nspin, 0.0_wp)

   call save_wavefunction(wfn, filename, error)
   call load_wavefunction(wfn, filename, error, nsh=10)
   call delete_file(filename)
end subroutine test_nsh_mismatch

subroutine test_nao_mismatch(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(wavefunction_type) :: wfn
   integer, parameter :: nat = 5, nsh = 12, nao = 20, nspin = 2
   character(len=*), parameter :: prefix = "tblite0_", filename = ".test-nao-mismatch.npz"
   integer :: stat

   call new_wavefunction(wfn, nat, nsh, nao, nspin, 0.0_wp)

   call save_wavefunction(wfn, filename, error)
   call load_wavefunction(wfn, filename, error, nao=10)
   call delete_file(filename)
end subroutine test_nao_mismatch

subroutine test_nspin_mismatch(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(wavefunction_type) :: wfn
   integer, parameter :: nat = 5, nsh = 12, nao = 20, nspin = 2
   character(len=*), parameter :: prefix = "tblite0_", filename = ".test-nspin-mismatch.npz"
   integer :: stat

   call new_wavefunction(wfn, nat, nsh, nao, nspin, 0.0_wp)

   call save_wavefunction(wfn, filename, error)
   call load_wavefunction(wfn, filename, error, nspin=1)
   call delete_file(filename)
end subroutine test_nspin_mismatch

subroutine test_kt_mismatch(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(wavefunction_type) :: wfn
   integer, parameter :: nat = 5, nsh = 12, nao = 20, nspin = 2
   character(len=*), parameter :: prefix = "tblite0_", filename = ".test-kt-mismatch.npz"
   integer :: stat

   call new_wavefunction(wfn, nat, nsh, nao, nspin, 0.0_wp)

   call save_npz(filename, prefix//"kt", [wfn%kt, wfn%kt])

   call save_npz(filename, prefix//"nel", wfn%nel)
   call save_npz(filename, prefix//"n0at", wfn%n0at)
   call save_npz(filename, prefix//"n0sh", wfn%n0sh)

   call save_npz(filename, prefix//"density", wfn%density)
   call save_npz(filename, prefix//"coeff", wfn%coeff)
   call save_npz(filename, prefix//"emo", wfn%emo)
   call save_npz(filename, prefix//"focc", wfn%focc)

   call save_npz(filename, prefix//"qat", wfn%qat)
   call save_npz(filename, prefix//"qsh", wfn%qsh)

   call save_npz(filename, prefix//"dpat", wfn%dpat)
   call save_npz(filename, prefix//"qpat", wfn%qpat)

   call load_wavefunction(wfn, filename, error)
   call delete_file(filename)
end subroutine test_kt_mismatch

subroutine test_dipole_mismatch(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(wavefunction_type) :: wfn
   integer, parameter :: nat = 5, nsh = 12, nao = 20, nspin = 2
   character(len=*), parameter :: prefix = "tblite0_", filename = ".test-dipole-mismatch.npz"
   integer :: stat

   call new_wavefunction(wfn, nat, nsh, nao, nspin, 0.0_wp)

   call save_npz(filename, prefix//"kt", [wfn%kt])

   call save_npz(filename, prefix//"nel", wfn%nel)
   call save_npz(filename, prefix//"n0at", wfn%n0at)
   call save_npz(filename, prefix//"n0sh", wfn%n0sh)

   call save_npz(filename, prefix//"density", wfn%density)
   call save_npz(filename, prefix//"coeff", wfn%coeff)
   call save_npz(filename, prefix//"emo", wfn%emo)
   call save_npz(filename, prefix//"focc", wfn%focc)

   call save_npz(filename, prefix//"qat", wfn%qat)
   call save_npz(filename, prefix//"qsh", wfn%qsh)

   call save_npz(filename, prefix//"dpat", wfn%qpat)
   call save_npz(filename, prefix//"qpat", wfn%qpat)

   call load_wavefunction(wfn, filename, error)
   call delete_file(filename)
end subroutine test_dipole_mismatch

subroutine test_quadrupole_mismatch(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(wavefunction_type) :: wfn
   integer, parameter :: nat = 5, nsh = 12, nao = 20, nspin = 2
   character(len=*), parameter :: prefix = "tblite0_", filename = ".test-quadrpole-mismatch.npz"
   integer :: stat

   call new_wavefunction(wfn, nat, nsh, nao, nspin, 0.0_wp)

   call save_npz(filename, prefix//"kt", [wfn%kt])

   call save_npz(filename, prefix//"nel", wfn%nel)
   call save_npz(filename, prefix//"n0at", wfn%n0at)
   call save_npz(filename, prefix//"n0sh", wfn%n0sh)

   call save_npz(filename, prefix//"density", wfn%density)
   call save_npz(filename, prefix//"coeff", wfn%coeff)
   call save_npz(filename, prefix//"emo", wfn%emo)
   call save_npz(filename, prefix//"focc", wfn%focc)

   call save_npz(filename, prefix//"qat", wfn%qat)
   call save_npz(filename, prefix//"qsh", wfn%qsh)

   call save_npz(filename, prefix//"dpat", wfn%dpat)
   call save_npz(filename, prefix//"qpat", wfn%dpat)

   call load_wavefunction(wfn, filename, error)
   call delete_file(filename)
end subroutine test_quadrupole_mismatch

subroutine test_load_error(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(wavefunction_type) :: wfn
   integer, parameter :: nat = 5, nsh = 12, nao = 20, nspin = 2
   character(len=*), parameter :: prefix = "tblite0_", filename = ".test-load-error.npz"
   integer :: stat

   call new_wavefunction(wfn, nat, nsh, nao, nspin, 0.0_wp)

   call save_npz(filename, prefix//"kt", [1, 1, 1])

   call save_npz(filename, prefix//"nel", wfn%nel)
   call save_npz(filename, prefix//"n0at", wfn%n0at)
   call save_npz(filename, prefix//"n0sh", wfn%n0sh)

   call save_npz(filename, prefix//"density", wfn%density)
   call save_npz(filename, prefix//"coeff", wfn%coeff)
   call save_npz(filename, prefix//"emo", wfn%emo)
   call save_npz(filename, prefix//"focc", wfn%focc)

   call save_npz(filename, prefix//"qat", wfn%qat)
   call save_npz(filename, prefix//"qsh", wfn%qsh)

   call save_npz(filename, prefix//"dpat", wfn%dpat)
   call save_npz(filename, prefix//"qpat", wfn%dpat)

   call load_wavefunction(wfn, filename, error)
   call delete_file(filename)
end subroutine test_load_error

subroutine delete_file(filename)
   character(len=*), intent(in) :: filename

   integer :: io

   open(newunit=io, file=filename)
   close(io, status="delete")
end subroutine delete_file

end module test_wavefunction_restart