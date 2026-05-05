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

module test_qvszp
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use multicharge, only : get_eeqbc_charges
   use tblite_basis_cache, only : basis_cache, cgto_cache
   use tblite_basis_qvszp, only : qvszp_basis_type, qvszp_cgto_type, &
      & new_qvszp_cgto, new_qvszp_basis, maxg
   use tblite_basis_type, only : get_cutoff, cgto_container
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_overlap, only : overlap_cgto, msao, get_overlap
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction

   implicit none
   private

   public :: collect_qvszp

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains

!> Collect all exported unit tests
subroutine collect_qvszp(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("norm-qvszp-mb01", test_qvszp_norm_mb01), &
      new_unittest("norm-qvszp-mb02", test_qvszp_norm_mb02), &
      new_unittest("norm-qvszp-accl3", test_qvszp_norm_accl3), &
      new_unittest("norm-qvszp-ce2", test_qvszp_norm_ce2), &
      new_unittest("norm-grad-qvszp-mb01", test_qvszp_norm_grad_mb01), &
      new_unittest("norm-grad-qvszp-mb02", test_qvszp_norm_grad_mb02), &
      new_unittest("norm-grad-qvszp-accl3", test_qvszp_norm_grad_accl3), &
      new_unittest("norm-grad-qvszp-ce2", test_qvszp_norm_grad_ce2), &
      new_unittest("coeffs-qvszp-mb01", test_qvszp_coeffs_mb01), &
      new_unittest("coeffs-qvszp-mb02", test_qvszp_coeffs_mb02), &
      new_unittest("coeffs-qvszp-accl3", test_qvszp_coeffs_accl3), &
      new_unittest("coeffs-qvszp-ce2", test_qvszp_coeffs_ce2), &
      new_unittest("coeffs-grad-qvszp-mb01", test_qvszp_coeffs_grad_mb01), &
      new_unittest("coeffs-grad-qvszp-mb02", test_qvszp_coeffs_grad_mb02), &
      new_unittest("coeffs-grad-qvszp-accl3", test_qvszp_coeffs_grad_accl3), &
      new_unittest("coeffs-grad-qvszp-ce2", test_qvszp_coeffs_grad_ce2), &
      new_unittest("overlap-qvszp-lih", test_qvszp_overlap_lih), &
      new_unittest("overlap-qvszp-sih4", test_qvszp_overlap_sih4), &
      new_unittest("overlap-qvszp-accl3", test_qvszp_overlap_accl3) &
      & ]

end subroutine collect_qvszp


subroutine make_qvszp_basis(bas, mol, error)
   !> Basis set information
   type(qvszp_basis_type), intent(out) :: bas
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Parameter: Number of shells selected from the q-vSZP basis set
   integer, parameter :: pa_nshell(103) = [ &
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, & !1-20
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, & !21-40
      & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 4, 4, 4, & !41-60
      & 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 2, & !61-80
      & 3, 3, 3, 3, 3, 3, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, & !81-100
      & 4, 4, 4]

   integer :: isp, izp, ish
   integer, allocatable :: nsh_id(:)
   type(cgto_container), allocatable :: cgto(:, :)
   type(qvszp_cgto_type), allocatable :: cgto_qvszp

   ! Setup q-vSZP CGTOs and basis 
   nsh_id = pa_nshell(mol%num)
   allocate(cgto(maxval(nsh_id), mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, nsh_id(isp)
         ! Get the q-vSZP basis set for the current shell
         allocate(cgto_qvszp)
         call new_qvszp_cgto(cgto_qvszp, izp, ish, .true., error)
         if (allocated(error)) return
         call move_alloc(cgto_qvszp, cgto(ish, isp)%raw)
      end do
   end do
   call new_qvszp_basis(bas, mol, nsh_id, cgto, error, accuracy=1.0_wp)

end subroutine make_qvszp_basis


subroutine test_qvszp_norm(error, mol, qdep, thr_in)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Flag for charge adaptation
   logical, intent(in) :: qdep
   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iao
   type(qvszp_basis_type) :: bas
   type(basis_cache) :: bcache
   type(wavefunction_type), allocatable :: wfn_aux
   real(wp) :: cutoff, thr_
   real(wp), allocatable :: lattr(:, :)
   real(wp), allocatable :: overlap(:,:)

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Setup q-vSZP CGTOs and basis
   call make_qvszp_basis(bas, mol, error)
   if (allocated(error)) return

   if (qdep) then
      ! Obtain EEQBC charges for charge adaptation
      allocate(wfn_aux)
      call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
      call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
      if (allocated(error)) return
   end if

   ! Update the basis cache
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   allocate(overlap(bas%nao, bas%nao))
   call get_overlap(mol, lattr, cutoff, bas, bcache, overlap)

   do iao = 1, bas%nao
      call check(error, overlap(iao, iao), 1.0_wp, thr=thr_)
      if (allocated(error)) then
         write(*,*) 'Diagonal error in the overlap:'
         print*, (overlap(iao, iao) - 1.0_wp)
         return
      end if
   end do

end subroutine test_qvszp_norm


subroutine test_qvszp_norm_grad(error, mol, qdep, thr_in)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Flag for charge adaptation
   logical, intent(in) :: qdep
   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, isp, ish, max_nsh
   type(qvszp_basis_type) :: bas
   type(basis_cache) :: bcache
   type(wavefunction_type), allocatable :: wfn_aux
   real(wp) :: thr_, normr, norml
   real(wp), allocatable :: numgrad(:, :)
   
   real(wp), parameter :: step = 1.0e-5_wp

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Setup q-vSZP CGTOs and basis
   call make_qvszp_basis(bas, mol, error)
   if (allocated(error)) return

   if (qdep) then
      ! Obtain EEQBC charges for charge adaptation
      allocate(wfn_aux)
      call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
      call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
      if (allocated(error)) return
   end if

   ! Update the basis cache
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   ! Numeric Gradient Calculation
   max_nsh = maxval(bas%nsh_at)
   allocate(numgrad(max_nsh, mol%nat), source=0.0_wp)
   do iat = 1, mol%nat
      isp = mol%id(iat)
      do ish = 1, bas%nsh_id(isp)
         normr = 0.0_wp
         norml = 0.0_wp
         ! Right hand side
         bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff + step
         ! Get the normalization constants
         call bas%cgto(ish, isp)%raw%get_normalization(bcache%cgto(ish, iat), .false.)
         normr = bcache%cgto(ish, iat)%norm

         ! Left hand side
         bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff - 2*step
         ! Get the normalization constants
         call bas%cgto(ish, isp)%raw%get_normalization(bcache%cgto(ish, iat), .false.)
         norml = bcache%cgto(ish, iat)%norm

         bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff + step
         numgrad(ish, iat) = 0.5_wp*(normr - norml)/step
      end do 
   end do 

   if (qdep) then
      ! Update charges with derivatives
      call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .true.)
      call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1), wfn_aux%dqatdr(:, :, :, 1), &
        & wfn_aux%dqatdL(:, :, :, 1))
      if (allocated(error)) return
   end if

   ! Update the basis cache
   call bas%update(mol, bcache, .true., wfn_aux=wfn_aux)

   if (any(abs(bcache%cgto(:, :)%dnorm - numgrad) > thr_)) then
      call test_failed(error, "Gradient of normalization constant does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', bcache%cgto(:, :)%dnorm
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', numgrad
      write(*,*) 'Difference:'
      print'(3es21.14)', bcache%cgto(:, :)%dnorm-numgrad
   end if

end subroutine test_qvszp_norm_grad


subroutine test_qvszp_coeffs(error, mol, qdep, ref, thr_in)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Flag for charge adaptation
   logical, intent(in) :: qdep
   !> Reference coefficients
   real(wp), intent(in) :: ref(:, :, :)
   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, isp, ish, max_nsh
   type(qvszp_basis_type) :: bas
   type(basis_cache) :: bcache
   type(wavefunction_type), allocatable :: wfn_aux
   real(wp) :: thr_
   real(wp), allocatable :: coeffs(:, :, :)

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Setup q-vSZP CGTOs and basis 
   call make_qvszp_basis(bas, mol, error)
   if (allocated(error)) return

   if (qdep) then
      ! Obtain EEQBC charges for charge adaptation
      allocate(wfn_aux)
      call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
      call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
      if (allocated(error)) return
   end if

   ! Initialize basis cache
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   ! Collect coefficient
   max_nsh = maxval(bas%nsh_at)
   allocate(coeffs(maxg, max_nsh, mol%nat), source=0.0_wp)
   do iat = 1, mol%nat
      isp = mol%id(iat)
      do ish = 1, bas%nsh_id(isp)
         associate(cgto => bas%cgto(ish, isp)%raw)
            call cgto%get_coeffs(bcache%cgto(ish, iat), &
               & coeffs(1:cgto%nprim, ish, iat))
         end associate
      end do
   end do

   if (any(abs(coeffs - ref) > thr_)) then
      call test_failed(error, "Coefficients do not match.")
      write(*,*) 'Reference:'
      print'(3es21.14)', ref
      write(*,*) 'Coefficients:'
      print'(3es21.14)', coeffs
      write(*,*) 'Difference:'
      print'(3es21.14)', coeffs-ref
   end if

end subroutine test_qvszp_coeffs


subroutine test_qvszp_coeffs_grad(error, mol, qdep, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Flag for charge adaptation
   logical, intent(in) :: qdep
   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, isp, ish, max_nsh, nprim
   type(qvszp_basis_type) :: bas
   type(basis_cache) :: bcache
   type(wavefunction_type), allocatable :: wfn_aux
   real(wp) :: thr_
   real(wp), allocatable :: numgrad(:, :, :), dCdqeff(:, :, :), cr(:), cl(:)
   
   real(wp), parameter :: step = 1.0e-5_wp

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Setup q-vSZP CGTOs and basis 
   call make_qvszp_basis(bas, mol, error)
   if (allocated(error)) return

   if (qdep) then
      ! Update charges with derivatives
      allocate(wfn_aux)
      call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .true.)
      call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1), wfn_aux%dqatdr(:, :, :, 1), &
        & wfn_aux%dqatdL(:, :, :, 1))
      if (allocated(error)) return
   end if

   ! Update the basis cache
   call bas%update(mol, bcache, .true., wfn_aux=wfn_aux)

   ! Gradient calculation (numerical and analytical)
   max_nsh = maxval(bas%nsh_at)
   allocate(numgrad(maxg, max_nsh, mol%nat), dCdqeff(maxg, max_nsh, mol%nat), &
      & cr(maxg), cl(maxg), source=0.0_wp)
   do iat = 1, mol%nat
      isp = mol%id(iat)
      do ish = 1, bas%nsh_id(isp)
         nprim = bas%cgto(ish, isp)%raw%nprim
         cr = 0.0_wp 
         cl = 0.0_wp

         ! 1. Right hand side
         bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff + step
         associate(cgto => bas%cgto(ish, isp)%raw)
            ! Update normalization for perturbed charge
            call cgto%get_normalization(bcache%cgto(ish, iat), .false.)
            ! Get perturbed coefficients
            call bas%cgto(ish, isp)%raw%get_coeffs(bcache%cgto(ish, iat), &
               & cr(1:cgto%nprim))
         end associate

         ! Left hand side
         bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff - 2*step
         associate(cgto => bas%cgto(ish, isp)%raw)
            ! Update normalization for perturbed charge
            call cgto%get_normalization(bcache%cgto(ish, iat), .false.)
            ! Get perturbed coefficients
            call bas%cgto(ish, isp)%raw%get_coeffs(bcache%cgto(ish, iat), &
               & cl(1:cgto%nprim))
         end associate

         bcache%cgto(ish, iat)%qeff = bcache%cgto(ish, iat)%qeff + step
         numgrad(1:nprim, ish, iat) = 0.5_wp * (cr(1:nprim) - cl(1:nprim)) / step

         ! Update the basis set cache
         call bas%update(mol, bcache, .true., wfn_aux=wfn_aux)
         associate(cgto => bas%cgto(ish, isp)%raw)
            ! Get analytic coefficient derivatives
            call bas%cgto(ish, isp)%raw%get_coeff_derivs(bcache%cgto(ish, iat), &
               & dCdqeff(1:cgto%nprim, ish, iat))
         end associate
      end do 
   end do 

   if (any(abs(dCdqeff - numgrad) > thr_)) then
      call test_failed(error, "Gradient of CGTO coefficients does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', dCdqeff
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', numgrad
      write(*,*) 'Difference:'
      print'(3es21.14)', dCdqeff-numgrad
   end if

end subroutine test_qvszp_coeffs_grad


subroutine test_qvszp_overlap(error, mol, qdep, ref, thr_in)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Flag for charge adaptation
   logical, intent(in) :: qdep
   !> Reference overlap matrix
   real(wp), intent(in) :: ref(:, :)
   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, isp, ish, max_nsh
   type(qvszp_basis_type) :: bas
   type(basis_cache) :: bcache
   type(wavefunction_type), allocatable :: wfn_aux
   real(wp) :: thr_
   real(wp), allocatable :: overlap(:,:), lattr(:, :)
   real(wp) :: cutoff

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Setup q-vSZP CGTOs and basis 
   call make_qvszp_basis(bas, mol, error)
   if (allocated(error)) return
   call check(error, bas%nao, size(ref, 1))
   if (allocated(error)) return

   if (qdep) then
      ! Obtain EEQBC charges for charge adaptation
      allocate(wfn_aux)
      call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
      call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
      if (allocated(error)) return
   end if

   ! Initialize basis cache
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   
   ! Get overlap matrix
   allocate(overlap(bas%nao, bas%nao))
   call get_overlap(mol, lattr, cutoff, bas, bcache, overlap)

   if (any(abs(overlap - ref) > thr_)) then
      call test_failed(error, "Coefficients do not match.")
      write(*,*) 'Reference:'
      print'(3es21.14)', ref
      write(*,*) 'Overlap:'
      print'(3es21.14)', overlap
      write(*,*) 'Difference:'
      print'(3es21.14)', overlap-ref
   end if

end subroutine test_qvszp_overlap


subroutine test_qvszp_norm_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_qvszp_norm(error, mol, .true., thr_in=thr)

end subroutine test_qvszp_norm_mb01

subroutine test_qvszp_norm_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_qvszp_norm(error, mol, .false., thr_in=thr)

end subroutine test_qvszp_norm_mb02

subroutine test_qvszp_norm_accl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "AcCl3")
   call test_qvszp_norm(error, mol, .true., thr_in=thr)

end subroutine test_qvszp_norm_accl3

subroutine test_qvszp_norm_ce2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "Ce2")
   call test_qvszp_norm(error, mol, .false., thr_in=thr)

end subroutine test_qvszp_norm_ce2


subroutine test_qvszp_norm_grad_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_qvszp_norm_grad(error, mol, .true., thr_in=thr1*100)

end subroutine test_qvszp_norm_grad_mb01

subroutine test_qvszp_norm_grad_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_qvszp_norm_grad(error, mol, .true., thr_in=thr1*10)

end subroutine test_qvszp_norm_grad_mb02

subroutine test_qvszp_norm_grad_accl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "AcCl3")
   call test_qvszp_norm_grad(error, mol, .true., thr_in=thr1*10)

end subroutine test_qvszp_norm_grad_accl3

subroutine test_qvszp_norm_grad_ce2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "Ce2")
   call test_qvszp_norm_grad(error, mol, .true., thr_in=thr1)

end subroutine test_qvszp_norm_grad_ce2


subroutine test_qvszp_coeffs_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   integer, parameter :: max_nsh = 3
   integer, parameter :: nat = 16
   real(wp), parameter :: coeffs(maxg, max_nsh, nat) = reshape([&
      & 5.78231778504263E-02_wp,-2.83825852749142E-01_wp, 9.90885477161646E-02_wp, &
      & 9.24361149893659E-03_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.46872829539007E-01_wp, 6.42693389305345E-02_wp, 5.20235024840289E-02_wp, &
      & 3.70736072326586E-03_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.22039577921690E-02_wp, 5.86140971684357E-02_wp, 1.01729322095785E-01_wp, &
      & 1.58596821810171E-01_wp, 1.91324695508849E-01_wp, 1.20872033083347E-01_wp, &
      & 8.47266497965465E-02_wp, 1.19353497297525E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.30878624567548E+00_wp, 1.11087228189312E+00_wp, 2.93160316956843E-01_wp, &
      & 8.11013269863575E-02_wp, 3.93116810672270E-02_wp, 1.30672685290969E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.45959284270498E+00_wp, 1.77781135416278E+00_wp, 1.28688356938413E+00_wp, &
      & 6.27177066211598E-01_wp, 1.82157413918656E-01_wp, 2.04454633055948E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.11224095611268E-02_wp, 5.76089148965212E-02_wp, 9.87763809435179E-02_wp, &
      & 1.55074905631894E-01_wp, 1.85516128281786E-01_wp, 1.19643965491189E-01_wp, &
      & 8.38670268534455E-02_wp, 1.39416689507856E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.04366474800031E+00_wp, 5.36275019085109E-01_wp, 3.49590820032480E-01_wp, &
      & 1.30955754526052E-01_wp, 1.30499431711682E-01_wp, 2.38994700818764E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 2.28242811794748E+00_wp, 2.67789308926715E+00_wp, 1.94578409878920E+00_wp, &
      & 1.02822189132134E+00_wp, 3.12632009568882E-01_wp, 4.03578657004155E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.65466341794698E-02_wp, 6.25440904960307E-02_wp, 1.13535055981685E-01_wp, &
      & 1.72539506800578E-01_wp, 2.14579569906112E-01_wp, 1.25491311875952E-01_wp, &
      & 8.79596146782660E-02_wp, 3.61448591586499E-03_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.03160644370461E-02_wp, 5.68530644156860E-02_wp, 9.65717253815038E-02_wp, &
      & 1.52437088753209E-01_wp, 1.81181443752145E-01_wp, 1.18709456096912E-01_wp, &
      & 8.32128621908047E-02_wp, 1.54213827823462E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.31125079779487E+00_wp, 1.11319825009748E+00_wp, 2.93967885020420E-01_wp, &
      & 8.12694083174010E-02_wp, 3.93968266221745E-02_wp, 1.26581517547700E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.46839971166148E+00_wp, 1.78734142452929E+00_wp, 1.29662411182965E+00_wp, &
      & 6.30757728525424E-01_wp, 1.82702181216329E-01_wp, 1.96791759810448E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.12232084953397E+00_wp, 9.26314123847090E-01_wp, 2.33143456592824E-01_wp, &
      & 8.90781961099526E-02_wp, 3.09827426052464E-02_wp, 1.11394239080102E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 9.16099390651602E-01_wp, 1.18266572197101E+00_wp, 8.87208496204766E-01_wp, &
      & 4.68853378024867E-01_wp, 1.22392434445639E-01_wp, 1.41112566630377E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.77210818719151E-02_wp, 6.35755813814786E-02_wp, 1.16712751697994E-01_wp, &
      & 1.76251490702988E-01_wp, 2.20848617647926E-01_wp, 1.26648425456669E-01_wp, &
      & 8.87693229559139E-02_wp, 1.28582865426124E-03_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.65966067275163E-02_wp, 6.25882661977842E-02_wp, 1.13670403988952E-01_wp, &
      & 1.72697986685960E-01_wp, 2.14846499937608E-01_wp, 1.25541387272827E-01_wp, &
      & 8.79946570811966E-02_wp, 3.51611645398002E-03_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 9.65615231920409E-02_wp,-4.27584012954246E-01_wp, 2.30026527824260E-01_wp, &
      & 1.84570650040638E-01_wp, 2.40675920781590E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-4.88486997320055E-01_wp, 3.17986635539090E-01_wp, 3.44727960957614E-01_wp, &
      & 1.27152412242587E-01_wp, 1.62227567465251E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.47087487453936E+00_wp, 3.90885636065272E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-7.11339708756932E-01_wp, 4.61755223237068E-01_wp, 2.04479318166316E-01_wp, &
      & 9.63588423911522E-02_wp,-6.73804126212010E-02_wp, 2.09249668502405E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.77466256015808E-01_wp, 4.45527106003421E-01_wp, 3.37676092043305E-01_wp, &
      & 1.12404294173381E-01_wp, 1.53008354105313E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-9.76128979597185E-01_wp, 6.36162649052069E-01_wp, 2.82667980924178E-01_wp, &
      & 1.34533778351486E-01_wp,-1.54550598807571E-01_wp, 3.03018270556527E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 4.15102986122589E-01_wp, 5.35453741021593E-01_wp, 3.84110024115255E-01_wp, &
      & 1.27843038999618E-01_wp, 3.34816394264130E-03_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.14339255013510E+00_wp, 9.43962946697837E-01_wp, 2.37600615517664E-01_wp, &
      & 9.07767373007212E-02_wp, 3.15876233492005E-02_wp, 8.62743129406114E-03_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 9.54841944180201E-01_wp, 1.21397248451509E+00_wp, 9.22755306541706E-01_wp, &
      & 4.82631432780964E-01_wp, 1.25159221696361E-01_wp, 1.07399243886539E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.02449119722697E-01_wp,-3.64903943231804E-01_wp, 5.07367994529731E-02_wp, &
      & 1.75325139616380E-01_wp, 1.24107869042705E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-2.16180247496466E-01_wp, 1.77698908164938E-01_wp, 1.01831199842622E-01_wp, &
      & 3.32262226908301E-04_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.22884197436823E-01_wp, 4.98040576830414E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp],&
      & shape(coeffs))

   call get_structure(mol, "MB16-43", "01")
   call test_qvszp_coeffs(error, mol, .true., coeffs, thr_in=thr1)

end subroutine test_qvszp_coeffs_mb01

subroutine test_qvszp_coeffs_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   integer, parameter :: max_nsh = 3
   integer, parameter :: nat = 16
   real(wp), parameter :: coeffs(maxg, max_nsh, nat) = reshape([&
      & 3.09198448290606E-02_wp, 5.74195463765993E-02_wp, 9.82227878717171E-02_wp, &
      & 1.54413208275303E-01_wp, 1.84427524708545E-01_wp, 1.19410706876252E-01_wp, &
      & 8.37037458761371E-02_wp, 1.43146705278220E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.11912519716861E-01_wp,-3.99363669309258E-01_wp, 1.89376581697049E-01_wp, &
      & 1.48432771751156E-01_wp, 1.82423059323969E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-4.65392657049056E-01_wp, 4.45994379700401E-01_wp, 2.13970350152015E-01_wp, &
      & 7.15187051054873E-02_wp, 1.25628040440038E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 8.59467249531832E-01_wp, 2.40333548081508E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-7.32506234542098E-01_wp, 4.75691864406586E-01_wp, 2.10725389909986E-01_wp, &
      & 9.94058643292783E-02_wp,-7.42253276609760E-02_wp, 2.16713759865377E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.81376180582447E-01_wp, 4.54368877798388E-01_wp, 3.42359568091845E-01_wp, &
      & 1.13961639252036E-01_wp, 1.42073223436629E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.31657604072179E+00_wp, 1.11822715158310E+00_wp, 2.95716210627004E-01_wp, &
      & 8.16327909711588E-02_wp, 3.95809508712867E-02_wp, 1.17684126671039E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.48734350826409E+00_wp, 1.80780487829380E+00_wp, 1.31763563793538E+00_wp, &
      & 6.38453705086415E-01_wp, 1.83857330977302E-01_wp, 1.80041374793366E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 9.18136319069972E-02_wp,-2.91732926834885E-01_wp, 1.35378429021833E-01_wp, &
      & 2.06819772600119E-02_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.36588148242061E-01_wp, 1.29155408899155E-01_wp, 1.18965176587476E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 2.33979261753881E-01_wp, 2.52308973149882E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.62925738750279E-02_wp, 6.23191151386609E-02_wp, 1.12846761217219E-01_wp, &
      & 1.71733072062377E-01_wp, 2.13222250983907E-01_wp, 1.25235593592078E-01_wp, &
      & 8.77806630507687E-02_wp, 4.11363096327537E-03_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.21793718777976E-02_wp, 5.85913590770360E-02_wp, 1.01662249229691E-01_wp, &
      & 1.58516970620712E-01_wp, 1.91192725937990E-01_wp, 1.20844445239769E-01_wp, &
      & 8.47073393161427E-02_wp, 1.19812371907194E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.18229683464308E-02_wp, 5.82611536392856E-02_wp, 1.00689661125612E-01_wp, &
      & 1.57358325745577E-01_wp, 1.89279285374572E-01_wp, 1.20442795601827E-01_wp, &
      & 8.44261968131542E-02_wp, 1.26449629371037E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 8.72383220549090E-02_wp,-3.29329440769257E-01_wp, 1.28832059576612E-01_wp, &
      & 8.35434613399693E-02_wp, 4.08665309993546E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.82859311749713E-01_wp, 1.99370606658476E-01_wp, 1.04577642340278E-01_wp, &
      & 1.03530726029163E-02_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 5.19596905379966E-01_wp, 1.54773896897002E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.57861444113624E-02_wp, 6.18687549054559E-02_wp, 1.11473835597803E-01_wp, &
      & 1.70121998415949E-01_wp, 2.10515426762305E-01_wp, 1.24720254150136E-01_wp, &
      & 8.74200196676238E-02_wp, 5.10383318168844E-03_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-7.08801833930433E-01_wp, 4.60084269437038E-01_wp, 2.03730456872988E-01_wp, &
      & 9.59935535161101E-02_wp,-6.65610577247126E-02_wp, 2.08355067017201E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.76984380177493E-01_wp, 4.44443621891247E-01_wp, 3.37100623466775E-01_wp, &
      & 1.12212938064388E-01_wp, 1.54337641301781E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 2.73196544572315E-02_wp,-2.59792866808796E-01_wp, 5.50273776792902E-02_wp, &
      & 6.73446670180127E-02_wp, 1.18596147649306E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 4.97107620136939E-02_wp, 1.31356008556578E-01_wp, 1.33508441120512E-01_wp, &
      & 6.97998914585093E-02_wp, 1.38575027970997E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.03410736586767E+00_wp, 5.30092767074584E-01_wp, 3.46108299603462E-01_wp, &
      & 1.29737324275681E-01_wp, 1.29319969042254E-01_wp, 2.60365565989879E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 2.24767440457712E+00_wp, 2.63631048956501E+00_wp, 1.91516937994751E+00_wp, &
      & 1.01233716356374E+00_wp, 3.07916675958105E-01_wp, 4.40407463206866E-02_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.22460596820579E-02_wp, 5.86530225568324E-02_wp, 1.01844174439499E-01_wp, &
      & 1.58733539254611E-01_wp, 1.91550677585503E-01_wp, 1.20919239689585E-01_wp, &
      & 8.47596926413675E-02_wp, 1.18567397045022E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 3.09469026959438E-02_wp, 5.74448615066049E-02_wp, 9.82967444468721E-02_wp, &
      & 1.54501632821707E-01_wp, 1.84572949237111E-01_wp, 1.19441923402381E-01_wp, &
      & 8.37255975187334E-02_wp, 1.42648962110969E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.11210770553683E-01_wp,-3.97451725784224E-01_wp, 1.87281140825621E-01_wp, &
      & 1.47705064672842E-01_wp, 1.91016455273785E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-4.58209817486145E-01_wp, 4.38636381784568E-01_wp, 2.10830522563045E-01_wp, &
      & 7.04562808601214E-02_wp, 1.37094225750490E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 8.48929042715517E-01_wp, 2.40825232525925E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp],&
      & shape(coeffs))

   call get_structure(mol, "MB16-43", "02")
   call test_qvszp_coeffs(error, mol, .true., coeffs, thr_in=thr1)

end subroutine test_qvszp_coeffs_mb02

subroutine test_qvszp_coeffs_accl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   integer, parameter :: max_nsh = 4
   integer, parameter :: nat = 4
   real(wp), parameter :: coeffs(maxg, max_nsh, nat) = reshape([&
      &-4.45942062245529E-01_wp, 3.34402817034430E-01_wp, 2.95404320366263E-01_wp, &
      &-2.37392419129959E-01_wp,-8.79261388754130E-02_wp, 7.63064577379641E-02_wp, &
      & 7.02758450239785E-03_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 7.77610624569431E-01_wp,-1.56248846005694E-01_wp,-7.40249168328618E-02_wp, &
      & 1.60008825974348E-02_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 2.70306382498416E+00_wp,-4.50214579677167E-01_wp,-1.74255563675224E+00_wp, &
      & 1.31740517888646E-01_wp, 3.58995581412661E-02_wp, 1.41481850538834E-03_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 2.04382306093911E+00_wp, 1.39308820542945E+00_wp, 1.67718826376502E+00_wp, &
      & 3.60100024374824E-01_wp, 3.10851619211348E-02_wp, 5.28785303570289E-03_wp, &
      & 1.20374788587376E-04_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 9.53328417110122E-02_wp,-4.23605499463639E-01_wp, 2.25971137203251E-01_wp, &
      & 1.82795320819728E-01_wp, 2.57971318975002E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-4.77335901731840E-01_wp, 3.10066877215319E-01_wp, 3.37178241233296E-01_wp, &
      & 1.24370249222245E-01_wp, 1.83544170934418E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.40627122017885E+00_wp, 3.93965361200963E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 9.53304499825054E-02_wp,-4.23597747295585E-01_wp, 2.25963248972874E-01_wp, &
      & 1.82791861991029E-01_wp, 2.58004879120182E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-4.77314332742814E-01_wp, 3.10051570065073E-01_wp, 3.37163632573930E-01_wp, &
      & 1.24364865708002E-01_wp, 1.83584962212038E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.40614406065904E+00_wp, 3.93971375684371E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 9.53333622114728E-02_wp,-4.23607186527050E-01_wp, 2.25972853884212E-01_wp, &
      & 1.82796073546572E-01_wp, 2.57964015386103E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-4.77340595761050E-01_wp, 3.10070208495071E-01_wp, 3.37181420492775E-01_wp, &
      & 1.24371420827929E-01_wp, 1.83535293363906E-02_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.40629889256866E+00_wp, 3.93964052307564E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp],&
      & shape(coeffs))

   call get_structure(mol, "f-block", "AcCl3")
   call test_qvszp_coeffs(error, mol, .true., coeffs, thr_in=thr1)

end subroutine test_qvszp_coeffs_accl3

subroutine test_qvszp_coeffs_ce2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   integer, parameter :: max_nsh = 4
   integer, parameter :: nat = 2
   real(wp), parameter :: coeffs(maxg, max_nsh, nat) = reshape([&
      &-4.86142688949204E-01_wp, 4.26564930231890E-01_wp, 2.34699884700147E-01_wp, &
      &-2.46168615922551E-01_wp,-5.52038782183826E-02_wp, 7.18065414269633E-02_wp, &
      & 7.79936444554429E-03_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 6.17175949108133E-01_wp,-1.12793579650491E-01_wp,-7.25536446170685E-02_wp, &
      & 2.00226362954444E-02_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 2.87446689976590E+00_wp,-2.34812746288182E+00_wp,-3.18867323164170E+00_wp, &
      &-3.30133433431023E-01_wp, 2.07406672362190E-01_wp, 3.92002043751221E-02_wp, &
      & 1.17365956425021E-03_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.57876360180012E+02_wp, 6.32047964168888E+01_wp, 1.40816452918020E+01_wp, &
      & 1.97209010727878E+00_wp, 1.85569847312590E-01_wp, 9.89268664157448E-03_wp, &
      & 5.50788903963886E-04_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-4.86142688949204E-01_wp, 4.26564930231890E-01_wp, 2.34699884700147E-01_wp, &
      &-2.46168615922551E-01_wp,-5.52038782183826E-02_wp, 7.18065414269633E-02_wp, &
      & 7.79936444554429E-03_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 6.17175949108133E-01_wp,-1.12793579650491E-01_wp,-7.25536446170685E-02_wp, &
      & 2.00226362954444E-02_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 2.87446689976590E+00_wp,-2.34812746288182E+00_wp,-3.18867323164170E+00_wp, &
      &-3.30133433431023E-01_wp, 2.07406672362190E-01_wp, 3.92002043751221E-02_wp, &
      & 1.17365956425021E-03_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.57876360180012E+02_wp, 6.32047964168888E+01_wp, 1.40816452918020E+01_wp, &
      & 1.97209010727878E+00_wp, 1.85569847312590E-01_wp, 9.89268664157448E-03_wp, &
      & 5.50788903963886E-04_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp],&
      & shape(coeffs))

   call get_structure(mol, "f-block", "Ce2")
   call test_qvszp_coeffs(error, mol, .true., coeffs, thr_in=thr1)

end subroutine test_qvszp_coeffs_ce2


subroutine test_qvszp_coeffs_grad_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_qvszp_coeffs_grad(error, mol, .true., thr_in=thr1*10)

end subroutine test_qvszp_coeffs_grad_mb01

subroutine test_qvszp_coeffs_grad_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_qvszp_coeffs_grad(error, mol, .true., thr_in=thr1*10)

end subroutine test_qvszp_coeffs_grad_mb02

subroutine test_qvszp_coeffs_grad_accl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "AcCl3")
   call test_qvszp_coeffs_grad(error, mol, .true., thr_in=thr1*10)

end subroutine test_qvszp_coeffs_grad_accl3

subroutine test_qvszp_coeffs_grad_ce2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "Ce2")
   call test_qvszp_coeffs_grad(error, mol, .true., thr_in=thr1*100)

end subroutine test_qvszp_coeffs_grad_ce2


subroutine test_qvszp_overlap_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 5
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp,-5.35147025308507E-17_wp, &
      & 0.00000000000000E+00_wp, 4.59550974647630E-01_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp,-5.35147025308507E-17_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 5.36198714710137E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 4.59550974647630E-01_wp, &
      & 0.00000000000000E+00_wp, 5.36198714710137E-01_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp], shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_qvszp_overlap(error, mol, .true., overlap, thr_in=thr1)

end subroutine test_qvszp_overlap_lih

subroutine test_qvszp_overlap_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 13
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 4.24354635991397E-01_wp, 4.24354635991397E-01_wp, 4.24354635991397E-01_wp, &
      & 4.24354635991397E-01_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 2.99012109170759E-01_wp,-2.99012109170759E-01_wp, &
      &-2.99012109170759E-01_wp, 2.99012109170759E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,-2.99012109170759E-01_wp, &
      &-2.99012109170759E-01_wp, 2.99012109170759E-01_wp, 2.99012109170759E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 2.99012109170759E-01_wp,-2.99012109170759E-01_wp, 2.99012109170759E-01_wp, &
      &-2.99012109170759E-01_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 2.07220673209862E-01_wp, 2.07220673209862E-01_wp, &
      &-2.07220673209862E-01_wp,-2.07220673209862E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,-2.07220673209862E-01_wp, &
      & 2.07220673209862E-01_wp,-2.07220673209862E-01_wp, 2.07220673209862E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp,-2.07220673209862E-01_wp, 2.07220673209862E-01_wp, &
      & 2.07220673209862E-01_wp,-2.07220673209862E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 4.24354635991397E-01_wp, 2.99012109170759E-01_wp,-2.99012109170759E-01_wp, &
      & 2.99012109170759E-01_wp, 2.07220673209862E-01_wp,-2.07220673209862E-01_wp, &
      & 0.00000000000000E+00_wp,-2.07220673209862E-01_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 8.47025274324054E-02_wp, 8.47025274324054E-02_wp, &
      & 8.47025274324054E-02_wp, 4.24354635991397E-01_wp,-2.99012109170759E-01_wp, &
      &-2.99012109170759E-01_wp,-2.99012109170759E-01_wp, 2.07220673209862E-01_wp, &
      & 2.07220673209862E-01_wp, 0.00000000000000E+00_wp, 2.07220673209862E-01_wp, &
      & 0.00000000000000E+00_wp, 8.47025274324054E-02_wp, 1.00000000000000E+00_wp, &
      & 8.47025274324054E-02_wp, 8.47025274324054E-02_wp, 4.24354635991397E-01_wp, &
      &-2.99012109170759E-01_wp, 2.99012109170759E-01_wp, 2.99012109170759E-01_wp, &
      &-2.07220673209862E-01_wp,-2.07220673209862E-01_wp, 0.00000000000000E+00_wp, &
      & 2.07220673209862E-01_wp, 0.00000000000000E+00_wp, 8.47025274324054E-02_wp, &
      & 8.47025274324054E-02_wp, 1.00000000000000E+00_wp, 8.47025274324054E-02_wp, &
      & 4.24354635991397E-01_wp, 2.99012109170759E-01_wp, 2.99012109170759E-01_wp, &
      &-2.99012109170759E-01_wp,-2.07220673209862E-01_wp, 2.07220673209862E-01_wp, &
      & 0.00000000000000E+00_wp,-2.07220673209862E-01_wp, 0.00000000000000E+00_wp, &
      & 8.47025274324054E-02_wp, 8.47025274324054E-02_wp, 8.47025274324054E-02_wp, &
      & 1.00000000000000E+00_wp], shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_qvszp_overlap(error, mol, .true., overlap, thr_in=thr1)

end subroutine test_qvszp_overlap_sih4


subroutine test_qvszp_overlap_accl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 43
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000000E+00_wp, 3.80787903596223E-19_wp,-1.39103434520588E-18_wp, &
      & 1.06470847525002E-20_wp,-3.01783389277602E-38_wp,-1.07627484423508E-36_wp, &
      & 0.00000000000000E+00_wp,-6.67477505145243E-38_wp, 0.00000000000000E+00_wp, &
      &-1.96904359169296E-35_wp,-1.50608272372749E-55_wp,-1.96904359169296E-35_wp, &
      & 3.15046974670873E-34_wp, 4.92260897923239E-36_wp, 0.00000000000000E+00_wp, &
      &-9.84521795846479E-36_wp, 2.54097308108428E-01_wp, 8.84015748872831E-02_wp, &
      & 1.18568805663372E-01_wp, 6.13839836625273E-02_wp, 6.03579005564660E-03_wp, &
      & 1.16586831846457E-02_wp, 5.30901581911619E-03_wp, 8.09551661319805E-03_wp, &
      &-2.25063789807228E-03_wp, 2.54015002903256E-01_wp,-1.46891904967046E-01_wp, &
      &-1.49666892064391E-02_wp, 6.20359068065403E-02_wp,-1.01525739787462E-02_wp, &
      & 2.44939466846409E-03_wp,-8.03334969409318E-03_wp,-1.03443698561445E-03_wp, &
      &-9.87606572009307E-03_wp, 2.54100770949792E-01_wp, 8.63967922051214E-03_wp, &
      &-1.52597050519540E-02_wp,-1.59163304337980E-01_wp,-1.52943032490713E-03_wp, &
      &-1.46633394881259E-04_wp,-8.00806831828508E-03_wp, 2.70133358657425E-03_wp, &
      & 1.40463499393994E-02_wp, 3.80787903596223E-19_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.45308601905858E-38_wp, 1.51962460181940E-20_wp, &
      & 8.78415155804615E-19_wp, 3.44998159597460E-20_wp,-1.55806156526104E-55_wp, &
      & 5.97554340940556E-20_wp,-2.25679311363995E-17_wp, 7.10603277090496E-38_wp, &
      &-2.25679311363995E-17_wp, 8.60945840649161E-37_wp, 1.18387963466327E-38_wp, &
      & 1.11147630094025E-36_wp, 4.58514610896495E-38_wp,-1.96944531188248E-01_wp, &
      & 1.11099714046488E-01_wp,-5.99769020303631E-02_wp,-3.10505039986088E-02_wp, &
      & 9.21528007593678E-03_wp, 1.78001603555801E-02_wp,-4.74680107941960E-03_wp, &
      & 2.28529789455802E-03_wp,-1.14528472297217E-02_wp, 3.27207402638203E-01_wp, &
      & 3.22941795761639E-02_wp,-1.25835373046284E-02_wp, 5.21579045812418E-02_wp, &
      & 1.23189251392594E-02_wp,-2.97204527841683E-03_wp, 1.41196984223314E-02_wp, &
      & 4.81518800022261E-04_wp, 2.25763429422497E-02_wp,-1.92480863136031E-02_wp, &
      & 1.55388722976938E-01_wp, 7.54386858235250E-04_wp, 7.86848138263948E-03_wp, &
      &-1.95185887127872E-02_wp,-1.87133528061882E-03_wp,-8.31407752597991E-04_wp, &
      & 7.45576160756131E-05_wp,-6.69530328480171E-04_wp,-1.39103434520588E-18_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.55806156526104E-55_wp,-5.97554340940557E-20_wp, 1.01430645332808E-18_wp, &
      & 1.51962460181940E-20_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-3.74375612102395E-38_wp,-1.40591867050797E-36_wp, 2.25679311363995E-17_wp, &
      & 8.98849947354287E-38_wp, 0.00000000000000E+00_wp, 8.48911217938689E-54_wp, &
      &-2.64152283199680E-01_wp,-5.99769020303630E-02_wp, 7.53726584380687E-02_wp, &
      &-4.16465563996464E-02_wp, 2.28529789455802E-03_wp, 1.52317735269676E-02_wp, &
      & 1.87636847823807E-02_wp, 1.05765868823361E-02_wp,-8.52146612532555E-04_wp, &
      & 3.33388793782123E-02_wp,-1.25835373046284E-02_wp, 1.54514302091994E-01_wp, &
      & 5.31432380635052E-03_wp, 4.81518800022266E-04_wp,-1.80953157569141E-02_wp, &
      &-1.73426461800191E-03_wp, 7.64207750034068E-03_wp, 4.68404497660954E-04_wp, &
      & 3.39966464568259E-02_wp, 7.54386858235250E-04_wp, 1.54483414153759E-01_wp, &
      &-1.38975883295297E-02_wp, 7.45576160756131E-05_wp, 1.06436196936872E-03_wp, &
      &-1.76577636714531E-03_wp,-1.96080622593259E-02_wp,-6.84740160431334E-04_wp, &
      & 1.06470847525002E-20_wp, 1.45308601905858E-38_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp,-5.97554340940557E-20_wp,-1.55806156526104E-55_wp, &
      &-8.77355672927606E-21_wp, 8.78415155804615E-19_wp, 1.51962460181940E-20_wp, &
      &-4.58514610896495E-38_wp,-1.11147630094025E-36_wp, 1.18387963466327E-38_wp, &
      &-5.50430931586381E-38_wp,-2.25679311363995E-17_wp, 7.10603277090496E-38_wp, &
      & 2.25679311363995E-17_wp,-1.36753670964550E-01_wp,-3.10505039986088E-02_wp, &
      &-4.16465563996465E-02_wp, 1.34256066540272E-01_wp, 1.20006269402234E-02_wp, &
      & 2.28529789455801E-03_wp,-3.29606752232506E-03_wp, 1.60958671305149E-02_wp, &
      & 7.07026322024399E-03_wp,-1.38187383035352E-01_wp, 5.21579045812417E-02_wp, &
      & 5.31432380635052E-03_wp, 1.33768917813984E-01_wp,-1.99750078696926E-02_wp, &
      & 4.81518800022263E-04_wp,-5.96308078148157E-03_wp,-2.03523628309494E-03_wp, &
      & 5.65151124424488E-03_wp, 3.54595227630930E-01_wp, 7.86848138263948E-03_wp, &
      &-1.38975883295297E-02_wp, 1.08598206891990E-02_wp, 1.83487211019493E-03_wp, &
      & 7.45576160756129E-05_wp, 1.53164951820822E-02_wp,-3.24081560148148E-03_wp, &
      &-2.66184213957890E-02_wp,-3.01783389277602E-38_wp, 1.51962460181940E-20_wp, &
      &-1.55806156526104E-55_wp,-5.97554340940557E-20_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 6.48243550947880E-39_wp, &
      & 0.00000000000000E+00_wp, 5.75303560479036E-21_wp, 1.15619026221186E-18_wp, &
      &-1.48542740583277E-21_wp, 1.82210360890911E-55_wp, 2.81203870812377E-20_wp, &
      & 0.00000000000000E+00_wp, 1.08909790854540E-19_wp, 9.39269140638547E-02_wp, &
      & 1.57879134513454E-02_wp, 1.13128049853127E-01_wp,-4.01647340285818E-02_wp, &
      & 2.03088780126363E-03_wp, 3.21222676800050E-02_wp, 6.23539859418376E-02_wp, &
      &-4.73536065838154E-03_wp,-1.32004971106857E-02_wp,-1.57596313406539E-01_wp, &
      & 1.65918515991503E-01_wp, 2.39592007111032E-02_wp, 6.46200339972795E-02_wp, &
      &-2.50063922595361E-03_wp,-1.94343664150374E-02_wp, 2.65362127613629E-02_wp, &
      &-7.80159109383848E-04_wp, 9.69870137039866E-02_wp,-2.38036029993225E-02_wp, &
      & 1.75681033204998E-01_wp, 3.68979292051442E-03_wp, 2.88359568781052E-02_wp, &
      &-9.99777769247366E-02_wp,-1.18826943267016E-02_wp, 3.99357968717228E-03_wp, &
      &-3.35887729018429E-03_wp,-2.08808975488039E-02_wp,-1.07627484423508E-36_wp, &
      & 8.78415155804615E-19_wp,-5.97554340940557E-20_wp,-1.55806156526104E-55_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 3.74263588640198E-39_wp, &
      & 0.00000000000000E+00_wp,-6.48243550947880E-39_wp, 0.00000000000000E+00_wp, &
      & 4.69733390126680E-21_wp, 1.46247785483871E-18_wp, 6.88805997185843E-20_wp, &
      &-1.11580602507777E-55_wp, 8.89244718622854E-20_wp, 1.44049938424137E-55_wp, &
      & 1.81428466445328E-01_wp, 3.04958057485193E-02_wp, 1.19785228301668E-01_wp, &
      & 1.13128049853127E-01_wp, 3.21222676800050E-02_wp, 4.74480545609624E-02_wp, &
      & 2.99992901960328E-02_wp, 6.46825514802640E-02_wp, 2.67194117879095E-02_wp, &
      & 3.80214486134918E-02_wp,-4.00292506432198E-02_wp, 1.58148927579536E-01_wp, &
      & 2.39592007111032E-02_wp,-1.94343664150374E-02_wp,-7.83660577007973E-02_wp, &
      &-2.53485649076460E-02_wp, 4.28895439268176E-02_wp,-1.24601861450841E-02_wp, &
      &-2.28215895902830E-03_wp, 1.68433343419291E-02_wp,-9.29593376498922E-03_wp, &
      & 3.68979292051442E-03_wp,-1.18826943267016E-02_wp, 2.28230442519527E-02_wp, &
      & 1.52057505474754E-03_wp, 6.46610665419610E-03_wp,-2.65879307064732E-03_wp, &
      & 0.00000000000000E+00_wp, 3.44998159597460E-20_wp, 1.01430645332808E-18_wp, &
      &-8.77355672927606E-21_wp, 0.00000000000000E+00_wp, 3.74263588640198E-39_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 6.72093437846902E-35_wp, 0.00000000000000E+00_wp,-9.74118783064143E-20_wp, &
      & 1.55119201273741E-18_wp, 5.14567147571516E-21_wp, 0.00000000000000E+00_wp, &
      & 1.05014599663578E-35_wp, 8.26171003312602E-02_wp, 1.31191984817152E-01_wp, &
      &-5.34044364485053E-02_wp, 9.10966423724774E-02_wp, 6.23539859418376E-02_wp, &
      & 2.99992901960328E-02_wp,-6.08947107781740E-02_wp, 2.08308046732050E-02_wp, &
      &-2.32506834337747E-02_wp,-1.24700031612179E-01_wp, 9.14204973299121E-02_wp, &
      & 3.82445305554984E-02_wp,-3.86090265072014E-02_wp, 2.65362127613629E-02_wp, &
      &-2.53485649076459E-02_wp, 4.30726416673874E-02_wp, 1.07052952349079E-02_wp, &
      & 2.58134914103780E-02_wp,-1.24635216090334E-01_wp,-5.36709154215204E-03_wp, &
      & 3.89999669170718E-02_wp, 9.88745071119326E-02_wp, 3.99357968717228E-03_wp, &
      & 1.52057505474754E-03_wp, 4.29241979763334E-02_wp,-2.80125851933167E-02_wp, &
      &-3.66771973089427E-02_wp,-6.67477505145243E-38_wp,-1.55806156526104E-55_wp, &
      & 1.51962460181940E-20_wp, 8.78415155804615E-19_wp, 6.48243550947880E-39_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp,-1.44049938424137E-55_wp,-8.89244718622854E-20_wp, &
      &-1.11580602507777E-55_wp,-3.63853919423637E-21_wp, 1.46247785483871E-18_wp, &
      & 4.69733390126681E-21_wp, 0.00000000000000E+00_wp, 1.25979678918588E-01_wp, &
      & 1.13128049853127E-01_wp, 8.31760577394356E-02_wp,-5.38710373613605E-02_wp, &
      &-4.73536065838155E-03_wp, 6.46825514802641E-02_wp, 2.08308046732049E-02_wp, &
      &-7.89872988986472E-04_wp,-5.39637251474037E-02_wp,-1.60573521281883E-02_wp, &
      & 2.39592007111032E-02_wp,-6.67900122547896E-02_wp, 6.58407939881903E-03_wp, &
      &-7.80159109383848E-04_wp, 4.28895439268176E-02_wp, 1.07052952349079E-02_wp, &
      & 5.07681213946361E-03_wp, 1.45016162485314E-02_wp, 4.20427601155703E-02_wp, &
      & 3.68979292051442E-03_wp, 1.71253063589435E-01_wp,-5.09310803815546E-02_wp, &
      &-3.35887729018430E-03_wp, 6.46610665419611E-03_wp,-2.80125851933167E-02_wp, &
      &-9.59469256266461E-02_wp, 2.47798997497724E-02_wp, 0.00000000000000E+00_wp, &
      & 5.97554340940556E-20_wp, 0.00000000000000E+00_wp, 1.51962460181940E-20_wp, &
      & 0.00000000000000E+00_wp,-6.48243550947880E-39_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp,-1.08909790854540E-19_wp, &
      & 0.00000000000000E+00_wp,-2.81203870812376E-20_wp, 0.00000000000000E+00_wp, &
      &-1.48542740583275E-21_wp, 1.15619026221186E-18_wp, 5.75303560479037E-21_wp, &
      &-3.50236622699164E-02_wp, 6.72812209188627E-02_wp,-4.21834215549399E-02_wp, &
      &-9.03958686999726E-02_wp,-1.32004971106857E-02_wp, 2.67194117879095E-02_wp, &
      &-2.32506834337745E-02_wp,-5.39637251474036E-02_wp,-2.84481487341263E-02_wp, &
      &-1.53304132696364E-01_wp, 6.48160526463661E-02_wp, 2.33066650210194E-02_wp, &
      &-1.65835732705046E-01_wp, 9.69870137039865E-02_wp,-1.24601861450841E-02_wp, &
      & 2.58134914103780E-02_wp, 1.45016162485314E-02_wp,-7.85753275275793E-03_wp, &
      & 2.18613252334541E-01_wp, 2.88358154032139E-02_wp,-3.38872073617419E-02_wp, &
      &-1.75683639506761E-01_wp,-2.08808975488039E-02_wp,-2.65879307064732E-03_wp, &
      &-3.66771973089427E-02_wp, 2.47798997497724E-02_wp, 8.95196235017786E-02_wp, &
      &-1.96904359169296E-35_wp,-2.25679311363995E-17_wp, 0.00000000000000E+00_wp, &
      &-4.58514610896495E-38_wp, 5.75303560479036E-21_wp, 0.00000000000000E+00_wp, &
      & 6.72093437846902E-35_wp,-1.44049938424137E-55_wp,-1.08909790854540E-19_wp, &
      & 1.00000000000000E+00_wp, 3.59788212745732E-37_wp,-1.53355933114843E-16_wp, &
      & 1.31853921310412E-52_wp, 7.60535815940642E-38_wp,-1.10951581687897E-36_wp, &
      & 8.24087008190073E-53_wp,-6.90829699178357E-03_wp,-2.35965157129680E-02_wp, &
      &-1.00686468390867E-02_wp, 3.79367235854106E-02_wp, 2.39145846637979E-03_wp, &
      &-2.35718098176066E-02_wp,-9.57453075945590E-03_wp, 3.00999157525351E-02_wp, &
      & 1.26504268561520E-02_wp,-3.29679821667578E-02_wp,-1.08925883548512E-02_wp, &
      & 6.06459295822026E-03_wp,-9.75228411749248E-02_wp, 7.27989483631098E-02_wp, &
      &-2.21369737805701E-03_wp, 7.81770089972111E-03_wp, 1.19554387789143E-02_wp, &
      &-3.65086815583784E-02_wp,-1.46929221355490E-02_wp, 9.88694820965624E-02_wp, &
      & 2.75607121079114E-03_wp, 1.78113863857822E-02_wp,-7.97563207247439E-02_wp, &
      &-1.07945544482274E-02_wp, 3.46636784544117E-03_wp,-3.63018922687454E-03_wp, &
      &-2.06584952521769E-02_wp,-1.50608272372749E-55_wp, 7.10603277090496E-38_wp, &
      &-3.74375612102395E-38_wp,-1.11147630094025E-36_wp, 1.15619026221186E-18_wp, &
      & 4.69733390126680E-21_wp, 0.00000000000000E+00_wp,-8.89244718622854E-20_wp, &
      & 0.00000000000000E+00_wp, 3.59788212745732E-37_wp, 1.00000000000000E+00_wp, &
      & 2.78690751225190E-37_wp, 1.52107163188128E-37_wp,-8.59427256225192E-37_wp, &
      & 0.00000000000000E+00_wp, 1.10951581687897E-36_wp,-7.05965539825512E-02_wp, &
      &-2.94596995776699E-02_wp,-6.76611859786402E-02_wp, 1.47842448386819E-02_wp, &
      &-7.94708085696002E-03_wp,-4.78952506397010E-02_wp,-4.56434544158290E-02_wp, &
      &-1.09822144323712E-02_wp, 1.93495163335749E-02_wp,-1.49480638222916E-02_wp, &
      & 2.09658246501309E-02_wp,-5.63527703772376E-02_wp, 2.86143317062409E-03_wp, &
      & 3.21447323224210E-03_wp, 4.81490264803972E-02_wp, 1.45916756667056E-02_wp, &
      &-8.80174082039844E-04_wp, 1.79387619533889E-02_wp,-2.30262119462451E-03_wp, &
      & 1.55254730408111E-02_wp,-8.49668395531560E-03_wp, 3.64904446403644E-03_wp, &
      &-1.51100766908978E-02_wp, 2.51138516221906E-02_wp, 2.24585339359639E-03_wp, &
      & 8.10425997898489E-03_wp,-3.93928241184114E-03_wp,-1.96904359169296E-35_wp, &
      &-2.25679311363995E-17_wp,-1.40591867050797E-36_wp, 1.18387963466327E-38_wp, &
      &-1.48542740583277E-21_wp, 1.46247785483871E-18_wp,-9.74118783064143E-20_wp, &
      &-1.11580602507777E-55_wp,-2.81203870812376E-20_wp,-1.53355933114843E-16_wp, &
      & 2.78690751225190E-37_wp, 1.00000000000000E+00_wp,-4.43806326751587E-37_wp, &
      &-1.17821701973295E-37_wp, 8.59427256225191E-37_wp,-7.60535815940642E-38_wp, &
      &-6.84799155786955E-02_wp,-4.46212021655553E-02_wp,-1.37272102926845E-02_wp, &
      &-6.28122403337026E-02_wp,-4.09105202500581E-02_wp,-4.26278351333780E-02_wp, &
      & 1.59808958342320E-02_wp,-3.89223133329991E-02_wp,-9.97144257955998E-03_wp, &
      &-6.24421611509400E-02_wp, 4.33251578953506E-02_wp, 2.95228364104495E-02_wp, &
      &-2.89208448206675E-02_wp, 2.95587486501468E-02_wp,-2.69476339326544E-02_wp, &
      & 3.90584941414175E-02_wp, 1.38995536511192E-02_wp, 1.77545188302525E-02_wp, &
      & 3.66899251388403E-03_wp,-2.48915244600366E-02_wp,-1.77101926740922E-03_wp, &
      &-4.35489790425279E-03_wp, 1.98701368923876E-02_wp, 5.96903307305780E-03_wp, &
      &-2.28634033709060E-03_wp, 2.13969759880309E-03_wp, 4.14067960720862E-03_wp, &
      & 3.15046974670873E-34_wp, 8.60945840649161E-37_wp, 2.25679311363995E-17_wp, &
      &-5.50430931586381E-38_wp, 1.82210360890911E-55_wp, 6.88805997185843E-20_wp, &
      & 1.55119201273741E-18_wp,-3.63853919423637E-21_wp, 0.00000000000000E+00_wp, &
      & 1.31853921310412E-52_wp, 1.52107163188128E-37_wp,-4.43806326751587E-37_wp, &
      & 1.00000000000000E+00_wp, 1.43915285098293E-37_wp, 0.00000000000000E+00_wp, &
      & 4.94452204914044E-53_wp, 1.11377180110526E-02_wp,-4.06104166340001E-02_wp, &
      & 5.78090962167377E-02_wp,-2.81989224102442E-02_wp,-3.85019952853237E-02_wp, &
      & 9.01015352853153E-03_wp, 4.96883636365857E-02_wp, 6.25643963580342E-03_wp, &
      & 1.43567037523914E-02_wp,-1.60584621291964E-02_wp, 1.79475165704033E-02_wp, &
      &-5.97886245721637E-02_wp,-7.57965842719622E-03_wp, 9.93010263126267E-03_wp, &
      & 4.33330289303085E-02_wp, 1.89823707176358E-02_wp,-1.83005574403083E-02_wp, &
      & 9.65965344344436E-03_wp,-1.63768200326637E-02_wp,-1.07618052470052E-03_wp, &
      &-5.96497855969987E-02_wp, 1.98257879724107E-02_wp, 1.52776930829613E-03_wp, &
      &-2.53958702119006E-03_wp, 1.93441892736988E-02_wp, 4.67851932496283E-02_wp, &
      &-1.40310950956885E-02_wp, 4.92260897923239E-36_wp, 1.18387963466327E-38_wp, &
      & 8.98849947354287E-38_wp,-2.25679311363995E-17_wp, 2.81203870812377E-20_wp, &
      &-1.11580602507777E-55_wp, 5.14567147571516E-21_wp, 1.46247785483871E-18_wp, &
      &-1.48542740583275E-21_wp, 7.60535815940642E-38_wp,-8.59427256225192E-37_wp, &
      &-1.17821701973295E-37_wp, 1.43915285098293E-37_wp, 1.00000000000000E+00_wp, &
      & 2.78690751225191E-37_wp, 6.35108148699737E-17_wp,-4.75508499079760E-02_wp, &
      &-6.28122403337026E-02_wp,-9.53185340207594E-03_wp, 2.22191439566478E-03_wp, &
      &-1.07800651365268E-02_wp,-3.89223133329991E-02_wp, 1.10967598716634E-02_wp, &
      &-1.36009695001246E-02_wp, 4.03490316805468E-02_wp, 2.63707934812862E-02_wp, &
      &-2.89208448206675E-02_wp,-1.24681882819496E-02_wp,-1.29411967704638E-02_wp, &
      & 2.85563588822670E-03_wp, 1.38995536511192E-02_wp,-1.64953208490581E-02_wp, &
      & 9.43567242540067E-05_wp,-2.32665999172373E-02_wp,-6.75915108878867E-02_wp, &
      &-4.35489790425279E-03_wp, 3.26263593187400E-02_wp, 5.50995923856506E-02_wp, &
      & 5.00961158918477E-03_wp, 2.13969759880309E-03_wp, 4.21197909788809E-02_wp, &
      &-3.33331054393402E-02_wp,-3.58779170386596E-02_wp, 0.00000000000000E+00_wp, &
      & 1.11147630094025E-36_wp, 0.00000000000000E+00_wp, 7.10603277090496E-38_wp, &
      & 0.00000000000000E+00_wp, 8.89244718622854E-20_wp, 0.00000000000000E+00_wp, &
      & 4.69733390126681E-21_wp, 1.15619026221186E-18_wp,-1.10951581687897E-36_wp, &
      & 0.00000000000000E+00_wp, 8.59427256225191E-37_wp, 0.00000000000000E+00_wp, &
      & 2.78690751225191E-37_wp, 1.00000000000000E+00_wp, 3.59788212745732E-37_wp, &
      & 2.63241892778876E-02_wp,-3.94472919947311E-02_wp, 2.52296431764687E-02_wp, &
      & 6.71168089133414E-02_wp, 1.93495163335749E-02_wp,-1.40185210728789E-02_wp, &
      & 1.70196258253043E-02_wp, 5.00035908751888E-02_wp, 3.67296165679995E-02_wp, &
      &-1.45409490249614E-02_wp, 1.19937250746915E-02_wp,-5.48179865441023E-02_wp, &
      &-1.71090429600988E-02_wp, 1.79387619533889E-02_wp, 3.28875100040146E-02_wp, &
      & 1.41942672028148E-02_wp,-3.38881391734682E-02_wp, 2.22365978086306E-03_wp, &
      & 2.11473661472861E-02_wp, 3.10192369717647E-03_wp, 7.80338889697969E-02_wp, &
      &-2.56047302565900E-02_wp,-3.93928241184114E-03_wp, 6.12494819781462E-03_wp, &
      &-2.06260083674999E-02_wp,-6.15774511801703E-02_wp, 2.06395253003514E-02_wp, &
      &-9.84521795846479E-36_wp, 4.58514610896495E-38_wp, 8.48911217938689E-54_wp, &
      & 2.25679311363995E-17_wp, 1.08909790854540E-19_wp, 1.44049938424137E-55_wp, &
      & 1.05014599663578E-35_wp, 0.00000000000000E+00_wp, 5.75303560479037E-21_wp, &
      & 8.24087008190073E-53_wp, 1.10951581687897E-36_wp,-7.60535815940642E-38_wp, &
      & 4.94452204914044E-53_wp, 6.35108148699737E-17_wp, 3.59788212745732E-37_wp, &
      & 1.00000000000000E+00_wp, 2.70517750245542E-02_wp,-1.37535295553788E-02_wp, &
      & 3.94271944904242E-02_wp, 4.32214705417513E-03_wp, 9.12736563336260E-03_wp, &
      & 1.49771210540312E-03_wp, 3.74923157441737E-02_wp, 1.29065608770154E-02_wp, &
      & 5.48727535960099E-05_wp, 8.44993729212913E-02_wp,-8.01724527297483E-02_wp, &
      &-1.55439996115142E-02_wp,-5.98525729418569E-03_wp,-2.81718581267872E-02_wp, &
      & 1.73384368231226E-02_wp,-2.00373447295399E-02_wp,-3.02269317344520E-03_wp, &
      &-6.92781117208129E-02_wp, 8.95164752724059E-02_wp, 2.04421151448564E-02_wp, &
      &-1.67913351825952E-02_wp,-7.47087628118558E-02_wp,-2.11528101890322E-02_wp, &
      &-2.80263201427949E-03_wp,-2.11188100405667E-02_wp, 1.83948733471546E-02_wp, &
      & 7.06498716627562E-02_wp, 2.54097308108428E-01_wp,-1.96944531188248E-01_wp, &
      &-2.64152283199680E-01_wp,-1.36753670964550E-01_wp, 9.39269140638547E-02_wp, &
      & 1.81428466445328E-01_wp, 8.26171003312602E-02_wp, 1.25979678918588E-01_wp, &
      &-3.50236622699164E-02_wp,-6.90829699178357E-03_wp,-7.05965539825512E-02_wp, &
      &-6.84799155786955E-02_wp, 1.11377180110526E-02_wp,-4.75508499079760E-02_wp, &
      & 2.63241892778876E-02_wp, 2.70517750245542E-02_wp, 1.00000000000000E+00_wp, &
      &-2.12970904433604E-16_wp,-2.36587368571639E-16_wp,-8.74653785794152E-17_wp, &
      & 2.82675664048246E-32_wp, 2.68801499323733E-31_wp, 0.00000000000000E+00_wp, &
      & 4.99109716129684E-32_wp, 0.00000000000000E+00_wp, 1.29382455777710E-03_wp, &
      &-4.70408009048437E-03_wp,-2.66961059912368E-03_wp, 1.31089245615731E-05_wp, &
      &-1.79174260621499E-06_wp, 3.64885390100942E-04_wp,-6.60523487340886E-05_wp, &
      &-1.01683112541574E-06_wp,-3.21476991791121E-04_wp, 1.31949180791217E-03_wp, &
      &-1.62520053543003E-03_wp,-2.72684695870236E-03_wp,-4.49381885326583E-03_wp, &
      & 2.09708745158415E-04_wp, 1.27251158228804E-04_wp,-6.60162410569535E-05_wp, &
      & 3.51860365217760E-04_wp, 2.52010492935596E-04_wp, 8.84015748872831E-02_wp, &
      & 1.11099714046488E-01_wp,-5.99769020303630E-02_wp,-3.10505039986088E-02_wp, &
      & 1.57879134513454E-02_wp, 3.04958057485193E-02_wp, 1.31191984817152E-01_wp, &
      & 1.13128049853127E-01_wp, 6.72812209188627E-02_wp,-2.35965157129680E-02_wp, &
      &-2.94596995776699E-02_wp,-4.46212021655553E-02_wp,-4.06104166340001E-02_wp, &
      &-6.28122403337026E-02_wp,-3.94472919947311E-02_wp,-1.37535295553788E-02_wp, &
      &-2.12970904433604E-16_wp, 1.00000000000000E+00_wp, 5.21981842528518E-32_wp, &
      & 3.34732613046338E-32_wp, 1.42104882963912E-16_wp, 4.96931221098021E-16_wp, &
      & 1.26001336902374E-16_wp,-6.79321935979898E-48_wp, 2.18240717336515E-16_wp, &
      & 4.70355454401235E-03_wp,-1.18691154225944E-02_wp,-7.88547357633682E-03_wp, &
      & 3.87210323027677E-05_wp,-9.35337759900350E-06_wp, 1.90479973079573E-03_wp, &
      &-1.56636584659408E-04_wp,-5.93515175398985E-06_wp,-1.47995712366862E-03_wp, &
      & 1.62524006318938E-03_wp, 4.35275317136063E-04_wp,-2.72869099621408E-03_wp, &
      &-4.49685781022301E-03_wp, 3.35340135917542E-05_wp, 2.03484221241593E-05_wp, &
      &-5.09220431996983E-05_wp, 6.95206092026806E-04_wp, 6.35643019044151E-04_wp, &
      & 1.18568805663372E-01_wp,-5.99769020303631E-02_wp, 7.53726584380687E-02_wp, &
      &-4.16465563996465E-02_wp, 1.13128049853127E-01_wp, 1.19785228301668E-01_wp, &
      &-5.34044364485053E-02_wp, 8.31760577394356E-02_wp,-4.21834215549399E-02_wp, &
      &-1.00686468390867E-02_wp,-6.76611859786402E-02_wp,-1.37272102926845E-02_wp, &
      & 5.78090962167377E-02_wp,-9.53185340207594E-03_wp, 2.52296431764687E-02_wp, &
      & 3.94271944904242E-02_wp,-2.36587368571639E-16_wp, 5.21981842528518E-32_wp, &
      & 1.00000000000000E+00_wp, 8.20365510201205E-32_wp,-6.79321935979898E-48_wp, &
      &-2.18240717336515E-16_wp, 5.73806748539344E-16_wp, 1.42104882963912E-16_wp, &
      & 0.00000000000000E+00_wp, 2.66931234645683E-03_wp,-7.88547357633682E-03_wp, &
      &-2.44932472562660E-03_wp, 2.19745574599337E-05_wp,-5.93515175398985E-06_wp, &
      & 8.12208603593496E-04_wp,-4.78610168691184E-04_wp,-2.26339286491003E-06_wp, &
      &-1.06489331954156E-03_wp, 2.72691328045650E-03_wp,-2.72869099621409E-03_wp, &
      &-2.51676640870696E-03_wp,-7.54506460969073E-03_wp, 6.95206092026806E-04_wp, &
      & 2.84130284308594E-04_wp,-4.85672551966368E-04_wp, 7.85644602357878E-04_wp, &
      & 8.35440743356501E-04_wp, 6.13839836625273E-02_wp,-3.10505039986088E-02_wp, &
      &-4.16465563996464E-02_wp, 1.34256066540272E-01_wp,-4.01647340285818E-02_wp, &
      & 1.13128049853127E-01_wp, 9.10966423724774E-02_wp,-5.38710373613605E-02_wp, &
      &-9.03958686999726E-02_wp, 3.79367235854106E-02_wp, 1.47842448386819E-02_wp, &
      &-6.28122403337026E-02_wp,-2.81989224102442E-02_wp, 2.22191439566478E-03_wp, &
      & 6.71168089133414E-02_wp, 4.32214705417513E-03_wp,-8.74653785794152E-17_wp, &
      & 3.34732613046338E-32_wp, 8.20365510201205E-32_wp, 1.00000000000000E+00_wp, &
      &-2.18240717336515E-16_wp,-6.79321935979898E-48_wp,-8.20442924323749E-17_wp, &
      & 4.96931221098021E-16_wp, 1.42104882963912E-16_wp,-1.31074600140053E-05_wp, &
      & 3.87210323027668E-05_wp, 2.19745574599332E-05_wp, 2.02564930595702E-03_wp, &
      &-3.96445997120969E-04_wp,-5.93515175398994E-06_wp, 4.36501320637859E-07_wp, &
      &-2.24986908287379E-04_wp, 6.33394226917790E-06_wp, 4.49392815091020E-03_wp, &
      &-4.49685781022301E-03_wp,-7.54506460969072E-03_wp,-1.03726225848758E-02_wp, &
      & 1.00797286939984E-03_wp, 6.95206092026806E-04_wp,-1.40803816383847E-04_wp, &
      & 1.69122990883716E-03_wp, 9.95989952642546E-04_wp, 6.03579005564660E-03_wp, &
      & 9.21528007593678E-03_wp, 2.28529789455802E-03_wp, 1.20006269402234E-02_wp, &
      & 2.03088780126363E-03_wp, 3.21222676800050E-02_wp, 6.23539859418376E-02_wp, &
      &-4.73536065838155E-03_wp,-1.32004971106857E-02_wp, 2.39145846637979E-03_wp, &
      &-7.94708085696002E-03_wp,-4.09105202500581E-02_wp,-3.85019952853237E-02_wp, &
      &-1.07800651365268E-02_wp, 1.93495163335749E-02_wp, 9.12736563336260E-03_wp, &
      & 2.82675664048246E-32_wp, 1.42104882963912E-16_wp,-6.79321935979898E-48_wp, &
      &-2.18240717336515E-16_wp, 1.00000000000000E+00_wp, 8.79747389274545E-32_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.79194394128006E-06_wp, 9.35524508620751E-06_wp, 5.93633950594922E-06_wp, &
      & 3.96526886714902E-04_wp,-6.12285793081369E-07_wp,-2.52966764276612E-08_wp, &
      & 2.71248570573907E-09_wp,-3.61311370947521E-07_wp, 2.31745400404824E-08_wp, &
      & 2.09703620833611E-04_wp,-3.35328867003048E-05_wp,-6.95175841802704E-04_wp, &
      &-1.00792913004316E-03_wp, 1.14546171462770E-06_wp, 7.26736963500946E-07_wp, &
      &-3.15498276364994E-07_wp, 2.88230422225898E-06_wp, 2.15843791949335E-06_wp, &
      & 1.16586831846457E-02_wp, 1.78001603555801E-02_wp, 1.52317735269676E-02_wp, &
      & 2.28529789455801E-03_wp, 3.21222676800050E-02_wp, 4.74480545609624E-02_wp, &
      & 2.99992901960328E-02_wp, 6.46825514802641E-02_wp, 2.67194117879095E-02_wp, &
      &-2.35718098176066E-02_wp,-4.78952506397010E-02_wp,-4.26278351333780E-02_wp, &
      & 9.01015352853153E-03_wp,-3.89223133329991E-02_wp,-1.40185210728789E-02_wp, &
      & 1.49771210540312E-03_wp, 2.68801499323733E-31_wp, 4.96931221098021E-16_wp, &
      &-2.18240717336515E-16_wp,-6.79321935979898E-48_wp, 8.79747389274545E-32_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 3.64926391650736E-04_wp,-1.90518004144670E-03_wp, &
      &-8.12369591870538E-04_wp, 5.93633950594930E-06_wp,-2.52966764276617E-08_wp, &
      & 4.53921590635269E-06_wp,-1.17832979789977E-06_wp,-1.31530791529738E-08_wp, &
      &-4.35807171297131E-06_wp, 1.27248048791160E-04_wp,-2.03477383270100E-05_wp, &
      &-2.84118041224261E-04_wp,-6.95175841802704E-04_wp, 7.26736963500946E-07_wp, &
      & 3.88789458191834E-07_wp,-4.18923418785128E-07_wp, 1.61223564290462E-06_wp, &
      & 1.44107420785261E-06_wp, 5.30901581911619E-03_wp,-4.74680107941960E-03_wp, &
      & 1.87636847823807E-02_wp,-3.29606752232506E-03_wp, 6.23539859418376E-02_wp, &
      & 2.99992901960328E-02_wp,-6.08947107781740E-02_wp, 2.08308046732049E-02_wp, &
      &-2.32506834337745E-02_wp,-9.57453075945590E-03_wp,-4.56434544158290E-02_wp, &
      & 1.59808958342320E-02_wp, 4.96883636365857E-02_wp, 1.10967598716634E-02_wp, &
      & 1.70196258253043E-02_wp, 3.74923157441737E-02_wp, 0.00000000000000E+00_wp, &
      & 1.26001336902374E-16_wp, 5.73806748539344E-16_wp,-8.20442924323749E-17_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      & 5.07922392016529E-32_wp, 0.00000000000000E+00_wp,-6.60597709240123E-05_wp, &
      & 1.56667034861073E-04_wp, 4.78706966073588E-04_wp,-4.36586176632834E-07_wp, &
      & 2.71248570573917E-09_wp,-1.17832979789977E-06_wp,-2.85802716870683E-07_wp, &
      & 3.28366782286897E-09_wp, 4.86678020566489E-07_wp,-6.60146279212627E-05_wp, &
      & 5.09198968253041E-05_wp, 4.85651186285602E-04_wp, 1.40797881474581E-04_wp, &
      &-3.15498276364994E-07_wp,-4.18923418785128E-07_wp,-3.15321837140630E-07_wp, &
      &-1.15835917868009E-06_wp,-3.79139535106234E-07_wp, 8.09551661319805E-03_wp, &
      & 2.28529789455802E-03_wp, 1.05765868823361E-02_wp, 1.60958671305149E-02_wp, &
      &-4.73536065838154E-03_wp, 6.46825514802640E-02_wp, 2.08308046732050E-02_wp, &
      &-7.89872988986472E-04_wp,-5.39637251474036E-02_wp, 3.00999157525351E-02_wp, &
      &-1.09822144323712E-02_wp,-3.89223133329991E-02_wp, 6.25643963580342E-03_wp, &
      &-1.36009695001246E-02_wp, 5.00035908751888E-02_wp, 1.29065608770154E-02_wp, &
      & 4.99109716129684E-32_wp,-6.79321935979898E-48_wp, 1.42104882963912E-16_wp, &
      & 4.96931221098021E-16_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 5.07922392016529E-32_wp, 1.00000000000000E+00_wp, 8.79747389274545E-32_wp, &
      &-1.01694538499750E-06_wp, 5.93633950594922E-06_wp, 2.26384149315157E-06_wp, &
      & 2.25032813908280E-04_wp,-3.61311370947521E-07_wp,-1.31530791529735E-08_wp, &
      & 3.28366782286887E-09_wp,-1.80672211752091E-07_wp, 1.41588495185786E-08_wp, &
      & 3.51851767356019E-04_wp,-6.95175841802703E-04_wp,-7.85610749179761E-04_wp, &
      &-1.69115652064343E-03_wp, 2.88230422225898E-06_wp, 1.61223564290462E-06_wp, &
      &-1.15835917868009E-06_wp, 4.26369032405522E-06_wp, 3.25838762232018E-06_wp, &
      &-2.25063789807228E-03_wp,-1.14528472297217E-02_wp,-8.52146612532555E-04_wp, &
      & 7.07026322024399E-03_wp,-1.32004971106857E-02_wp, 2.67194117879095E-02_wp, &
      &-2.32506834337747E-02_wp,-5.39637251474037E-02_wp,-2.84481487341263E-02_wp, &
      & 1.26504268561520E-02_wp, 1.93495163335749E-02_wp,-9.97144257955998E-03_wp, &
      & 1.43567037523914E-02_wp, 4.03490316805468E-02_wp, 3.67296165679995E-02_wp, &
      & 5.48727535960099E-05_wp, 0.00000000000000E+00_wp, 2.18240717336515E-16_wp, &
      & 0.00000000000000E+00_wp, 1.42104882963912E-16_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 8.79747389274545E-32_wp, &
      & 1.00000000000000E+00_wp,-3.21513115613133E-04_wp, 1.48025174283598E-03_wp, &
      & 1.06510642767749E-03_wp,-6.33521415324573E-06_wp, 2.31745400404829E-08_wp, &
      &-4.35807171297131E-06_wp, 4.86678020566489E-07_wp, 1.41588495185789E-08_wp, &
      & 3.54559378534955E-06_wp, 2.52004334949106E-04_wp,-6.35615240401303E-04_wp, &
      &-8.35404391158209E-04_wp,-9.95946946791490E-04_wp, 2.15843791949334E-06_wp, &
      & 1.44107420785261E-06_wp,-3.79139535106234E-07_wp, 3.25838762232018E-06_wp, &
      & 1.94316371017185E-06_wp, 2.54015002903256E-01_wp, 3.27207402638203E-01_wp, &
      & 3.33388793782123E-02_wp,-1.38187383035352E-01_wp,-1.57596313406539E-01_wp, &
      & 3.80214486134918E-02_wp,-1.24700031612179E-01_wp,-1.60573521281883E-02_wp, &
      &-1.53304132696364E-01_wp,-3.29679821667578E-02_wp,-1.49480638222916E-02_wp, &
      &-6.24421611509400E-02_wp,-1.60584621291964E-02_wp, 2.63707934812862E-02_wp, &
      &-1.45409490249614E-02_wp, 8.44993729212913E-02_wp, 1.29382455777710E-03_wp, &
      & 4.70355454401235E-03_wp, 2.66931234645683E-03_wp,-1.31074600140053E-05_wp, &
      &-1.79194394128006E-06_wp, 3.64926391650736E-04_wp,-6.60597709240123E-05_wp, &
      &-1.01694538499750E-06_wp,-3.21513115613133E-04_wp, 1.00000000000000E+00_wp, &
      & 7.89097014481192E-16_wp,-2.57700489830573E-17_wp,-1.66137508501279E-16_wp, &
      &-1.21868953836835E-31_wp, 6.84685277291211E-32_wp, 0.00000000000000E+00_wp, &
      &-1.69769255313159E-32_wp, 0.00000000000000E+00_wp, 1.30055599019737E-03_wp, &
      & 3.12486060309987E-03_wp,-5.87037410579942E-06_wp,-4.44405940448674E-03_wp, &
      &-4.02364515623123E-04_wp,-5.31502848774173E-07_wp,-2.46860584227914E-04_wp, &
      & 7.55883392450601E-07_wp, 1.44651631713556E-04_wp,-1.46891904967046E-01_wp, &
      & 3.22941795761639E-02_wp,-1.25835373046284E-02_wp, 5.21579045812417E-02_wp, &
      & 1.65918515991503E-01_wp,-4.00292506432198E-02_wp, 9.14204973299121E-02_wp, &
      & 2.39592007111032E-02_wp, 6.48160526463661E-02_wp,-1.08925883548512E-02_wp, &
      & 2.09658246501309E-02_wp, 4.33251578953506E-02_wp, 1.79475165704033E-02_wp, &
      &-2.89208448206675E-02_wp, 1.19937250746915E-02_wp,-8.01724527297483E-02_wp, &
      &-4.70408009048437E-03_wp,-1.18691154225944E-02_wp,-7.88547357633682E-03_wp, &
      & 3.87210323027668E-05_wp, 9.35524508620751E-06_wp,-1.90518004144670E-03_wp, &
      & 1.56667034861073E-04_wp, 5.93633950594922E-06_wp, 1.48025174283598E-03_wp, &
      & 7.89097014481192E-16_wp, 1.00000000000000E+00_wp, 1.46714904086106E-32_wp, &
      &-2.07148771992988E-32_wp, 1.27754620656300E-16_wp,-1.89119697801711E-17_wp, &
      &-4.83950076282793E-16_wp, 1.02807148607052E-47_wp,-8.38226120448632E-16_wp, &
      &-3.12528576924840E-03_wp,-4.06591770985219E-03_wp, 1.14615352525398E-05_wp, &
      & 8.67674574241973E-03_wp, 1.17681382581479E-03_wp, 1.55451059079795E-06_wp, &
      & 7.99881727160953E-04_wp,-2.91563193175748E-06_wp,-8.21786212367834E-04_wp, &
      &-1.49666892064391E-02_wp,-1.25835373046284E-02_wp, 1.54514302091994E-01_wp, &
      & 5.31432380635052E-03_wp, 2.39592007111032E-02_wp, 1.58148927579536E-01_wp, &
      & 3.82445305554984E-02_wp,-6.67900122547896E-02_wp, 2.33066650210194E-02_wp, &
      & 6.06459295822026E-03_wp,-5.63527703772376E-02_wp, 2.95228364104495E-02_wp, &
      &-5.97886245721637E-02_wp,-1.24681882819496E-02_wp,-5.48179865441023E-02_wp, &
      &-1.55439996115142E-02_wp,-2.66961059912368E-03_wp,-7.88547357633682E-03_wp, &
      &-2.44932472562660E-03_wp, 2.19745574599332E-05_wp, 5.93633950594922E-06_wp, &
      &-8.12369591870538E-04_wp, 4.78706966073588E-04_wp, 2.26384149315157E-06_wp, &
      & 1.06510642767749E-03_wp,-2.57700489830573E-17_wp, 1.46714904086106E-32_wp, &
      & 1.00000000000000E+00_wp,-2.19414423356693E-33_wp, 1.02807148607052E-47_wp, &
      & 8.38226120448632E-16_wp,-2.18376616869757E-17_wp, 1.27754620656300E-16_wp, &
      & 0.00000000000000E+00_wp, 5.87117282441953E-06_wp, 1.14615352525397E-05_wp, &
      & 2.03515403058432E-03_wp,-1.63001650308424E-05_wp,-2.91563193175748E-06_wp, &
      & 2.63824717478037E-04_wp,-2.36111571061122E-06_wp,-3.75201606011243E-04_wp, &
      & 1.04818118901892E-06_wp, 6.20359068065403E-02_wp, 5.21579045812418E-02_wp, &
      & 5.31432380635052E-03_wp, 1.33768917813984E-01_wp, 6.46200339972795E-02_wp, &
      & 2.39592007111032E-02_wp,-3.86090265072014E-02_wp, 6.58407939881903E-03_wp, &
      &-1.65835732705046E-01_wp,-9.75228411749248E-02_wp, 2.86143317062409E-03_wp, &
      &-2.89208448206675E-02_wp,-7.57965842719622E-03_wp,-1.29411967704638E-02_wp, &
      &-1.71090429600988E-02_wp,-5.98525729418569E-03_wp, 1.31089245615731E-05_wp, &
      & 3.87210323027677E-05_wp, 2.19745574599337E-05_wp, 2.02564930595702E-03_wp, &
      & 3.96526886714902E-04_wp, 5.93633950594930E-06_wp,-4.36586176632834E-07_wp, &
      & 2.25032813908280E-04_wp,-6.33521415324573E-06_wp,-1.66137508501279E-16_wp, &
      &-2.07148771992988E-32_wp,-2.19414423356693E-33_wp, 1.00000000000000E+00_wp, &
      & 8.38226120448632E-16_wp, 1.02807148607052E-47_wp,-7.37591646261335E-17_wp, &
      &-1.89119697801711E-17_wp, 1.27754620656300E-16_wp, 4.44466405981728E-03_wp, &
      & 8.67674574241973E-03_wp,-1.63001650308425E-05_wp,-1.03045664021254E-02_wp, &
      &-1.94339728637851E-03_wp,-2.91563193175748E-06_wp,-1.13756175508771E-03_wp, &
      & 3.65087296883579E-06_wp, 4.18299324733108E-04_wp,-1.01525739787462E-02_wp, &
      & 1.23189251392594E-02_wp, 4.81518800022266E-04_wp,-1.99750078696926E-02_wp, &
      &-2.50063922595361E-03_wp,-1.94343664150374E-02_wp, 2.65362127613629E-02_wp, &
      &-7.80159109383848E-04_wp, 9.69870137039865E-02_wp, 7.27989483631098E-02_wp, &
      & 3.21447323224210E-03_wp, 2.95587486501468E-02_wp, 9.93010263126267E-03_wp, &
      & 2.85563588822670E-03_wp, 1.79387619533889E-02_wp,-2.81718581267872E-02_wp, &
      &-1.79174260621499E-06_wp,-9.35337759900350E-06_wp,-5.93515175398985E-06_wp, &
      &-3.96445997120969E-04_wp,-6.12285793081369E-07_wp,-2.52966764276617E-08_wp, &
      & 2.71248570573917E-09_wp,-3.61311370947521E-07_wp, 2.31745400404829E-08_wp, &
      &-1.21868953836835E-31_wp, 1.27754620656300E-16_wp, 1.02807148607052E-47_wp, &
      & 8.38226120448632E-16_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 4.06307374054129E-31_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-4.02309480439717E-04_wp,-1.17652899287768E-03_wp, 2.91492178236114E-06_wp, &
      & 1.94292520555368E-03_wp, 5.70674555842992E-06_wp, 7.87968341926250E-09_wp, &
      & 3.54611710428893E-06_wp,-1.17502876050916E-08_wp,-2.35045189264994E-06_wp, &
      & 2.44939466846409E-03_wp,-2.97204527841683E-03_wp,-1.80953157569141E-02_wp, &
      & 4.81518800022263E-04_wp,-1.94343664150374E-02_wp,-7.83660577007973E-02_wp, &
      &-2.53485649076459E-02_wp, 4.28895439268176E-02_wp,-1.24601861450841E-02_wp, &
      &-2.21369737805701E-03_wp, 4.81490264803972E-02_wp,-2.69476339326544E-02_wp, &
      & 4.33330289303085E-02_wp, 1.38995536511192E-02_wp, 3.28875100040146E-02_wp, &
      & 1.73384368231226E-02_wp, 3.64885390100942E-04_wp, 1.90479973079573E-03_wp, &
      & 8.12208603593496E-04_wp,-5.93515175398994E-06_wp,-2.52966764276612E-08_wp, &
      & 4.53921590635269E-06_wp,-1.17832979789977E-06_wp,-1.31530791529735E-08_wp, &
      &-4.35807171297131E-06_wp, 6.84685277291211E-32_wp,-1.89119697801711E-17_wp, &
      & 8.38226120448632E-16_wp, 1.02807148607052E-47_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp,-3.51872507675822E-31_wp, &
      & 0.00000000000000E+00_wp,-5.31430150125044E-07_wp,-1.55413434112477E-06_wp, &
      &-2.63759193622194E-04_wp, 2.91492178236113E-06_wp, 7.87968341926241E-09_wp, &
      &-2.58414322342849E-07_wp, 5.60586440485866E-09_wp, 4.02801096758346E-07_wp, &
      &-3.63692714880030E-09_wp,-8.03334969409318E-03_wp, 1.41196984223314E-02_wp, &
      &-1.73426461800191E-03_wp,-5.96308078148157E-03_wp, 2.65362127613629E-02_wp, &
      &-2.53485649076460E-02_wp, 4.30726416673874E-02_wp, 1.07052952349079E-02_wp, &
      & 2.58134914103780E-02_wp, 7.81770089972111E-03_wp, 1.45916756667056E-02_wp, &
      & 3.90584941414175E-02_wp, 1.89823707176358E-02_wp,-1.64953208490581E-02_wp, &
      & 1.41942672028148E-02_wp,-2.00373447295399E-02_wp,-6.60523487340886E-05_wp, &
      &-1.56636584659408E-04_wp,-4.78610168691184E-04_wp, 4.36501320637859E-07_wp, &
      & 2.71248570573907E-09_wp,-1.17832979789977E-06_wp,-2.85802716870683E-07_wp, &
      & 3.28366782286887E-09_wp, 4.86678020566489E-07_wp, 0.00000000000000E+00_wp, &
      &-4.83950076282793E-16_wp,-2.18376616869757E-17_wp,-7.37591646261335E-17_wp, &
      & 4.06307374054129E-31_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,-2.46826818781407E-04_wp, &
      &-7.99687633395845E-04_wp, 2.36053787840104E-06_wp, 1.13728572222362E-03_wp, &
      & 3.54611710428893E-06_wp, 5.60586440485871E-09_wp, 2.20044667218825E-06_wp, &
      &-7.97244984431508E-09_wp,-1.27484309740475E-06_wp,-1.03443698561445E-03_wp, &
      & 4.81518800022261E-04_wp, 7.64207750034068E-03_wp,-2.03523628309494E-03_wp, &
      &-7.80159109383848E-04_wp, 4.28895439268176E-02_wp, 1.07052952349079E-02_wp, &
      & 5.07681213946361E-03_wp, 1.45016162485314E-02_wp, 1.19554387789143E-02_wp, &
      &-8.80174082039844E-04_wp, 1.38995536511192E-02_wp,-1.83005574403083E-02_wp, &
      & 9.43567242540067E-05_wp,-3.38881391734682E-02_wp,-3.02269317344520E-03_wp, &
      &-1.01683112541574E-06_wp,-5.93515175398985E-06_wp,-2.26339286491003E-06_wp, &
      &-2.24986908287379E-04_wp,-3.61311370947521E-07_wp,-1.31530791529738E-08_wp, &
      & 3.28366782286897E-09_wp,-1.80672211752091E-07_wp, 1.41588495185789E-08_wp, &
      &-1.69769255313159E-32_wp, 1.02807148607052E-47_wp, 1.27754620656300E-16_wp, &
      &-1.89119697801711E-17_wp, 0.00000000000000E+00_wp,-3.51872507675822E-31_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 7.55780003161801E-07_wp, 2.91492178236114E-06_wp, 3.75108420444024E-04_wp, &
      &-3.64998611613997E-06_wp,-1.17502876050915E-08_wp, 4.02801096758346E-07_wp, &
      &-7.97244984431500E-09_wp,-5.48031491053023E-07_wp, 3.65883314752584E-09_wp, &
      &-9.87606572009307E-03_wp, 2.25763429422497E-02_wp, 4.68404497660954E-04_wp, &
      & 5.65151124424488E-03_wp, 9.69870137039866E-02_wp,-1.24601861450841E-02_wp, &
      & 2.58134914103780E-02_wp, 1.45016162485314E-02_wp,-7.85753275275793E-03_wp, &
      &-3.65086815583784E-02_wp, 1.79387619533889E-02_wp, 1.77545188302525E-02_wp, &
      & 9.65965344344436E-03_wp,-2.32665999172373E-02_wp, 2.22365978086306E-03_wp, &
      &-6.92781117208129E-02_wp,-3.21476991791121E-04_wp,-1.47995712366862E-03_wp, &
      &-1.06489331954156E-03_wp, 6.33394226917790E-06_wp, 2.31745400404824E-08_wp, &
      &-4.35807171297131E-06_wp, 4.86678020566489E-07_wp, 1.41588495185786E-08_wp, &
      & 3.54559378534955E-06_wp, 0.00000000000000E+00_wp,-8.38226120448632E-16_wp, &
      & 0.00000000000000E+00_wp, 1.27754620656300E-16_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 1.44631846347869E-04_wp, 8.21584787949762E-04_wp, &
      &-1.04792588750760E-06_wp,-4.18199240292262E-04_wp,-2.35045189264994E-06_wp, &
      &-3.63692714880034E-09_wp,-1.27484309740475E-06_wp, 3.65883314752589E-09_wp, &
      & 1.37003353935441E-08_wp, 2.54100770949792E-01_wp,-1.92480863136031E-02_wp, &
      & 3.39966464568259E-02_wp, 3.54595227630930E-01_wp,-2.38036029993225E-02_wp, &
      &-2.28215895902830E-03_wp,-1.24635216090334E-01_wp, 4.20427601155703E-02_wp, &
      & 2.18613252334541E-01_wp,-1.46929221355490E-02_wp,-2.30262119462451E-03_wp, &
      & 3.66899251388403E-03_wp,-1.63768200326637E-02_wp,-6.75915108878867E-02_wp, &
      & 2.11473661472861E-02_wp, 8.95164752724059E-02_wp, 1.31949180791217E-03_wp, &
      & 1.62524006318938E-03_wp, 2.72691328045650E-03_wp, 4.49392815091020E-03_wp, &
      & 2.09703620833611E-04_wp, 1.27248048791160E-04_wp,-6.60146279212627E-05_wp, &
      & 3.51851767356019E-04_wp, 2.52004334949106E-04_wp, 1.30055599019737E-03_wp, &
      &-3.12528576924840E-03_wp, 5.87117282441953E-06_wp, 4.44466405981728E-03_wp, &
      &-4.02309480439717E-04_wp,-5.31430150125044E-07_wp,-2.46826818781407E-04_wp, &
      & 7.55780003161801E-07_wp, 1.44631846347869E-04_wp, 1.00000000000000E+00_wp, &
      &-3.39826415555430E-17_wp, 1.71756912647956E-17_wp, 1.32489189884931E-16_wp, &
      &-1.41339118481328E-32_wp, 1.15788362790797E-33_wp, 0.00000000000000E+00_wp, &
      & 2.56643374397691E-32_wp, 0.00000000000000E+00_wp, 8.63967922051214E-03_wp, &
      & 1.55388722976938E-01_wp, 7.54386858235250E-04_wp, 7.86848138263948E-03_wp, &
      & 1.75681033204998E-01_wp, 1.68433343419291E-02_wp,-5.36709154215204E-03_wp, &
      & 3.68979292051442E-03_wp, 2.88358154032139E-02_wp, 9.88694820965624E-02_wp, &
      & 1.55254730408111E-02_wp,-2.48915244600366E-02_wp,-1.07618052470052E-03_wp, &
      &-4.35489790425279E-03_wp, 3.10192369717647E-03_wp, 2.04421151448564E-02_wp, &
      &-1.62520053543003E-03_wp, 4.35275317136063E-04_wp,-2.72869099621409E-03_wp, &
      &-4.49685781022301E-03_wp,-3.35328867003048E-05_wp,-2.03477383270100E-05_wp, &
      & 5.09198968253041E-05_wp,-6.95175841802703E-04_wp,-6.35615240401303E-04_wp, &
      & 3.12486060309987E-03_wp,-4.06591770985219E-03_wp, 1.14615352525397E-05_wp, &
      & 8.67674574241973E-03_wp,-1.17652899287768E-03_wp,-1.55413434112477E-06_wp, &
      &-7.99687633395845E-04_wp, 2.91492178236114E-06_wp, 8.21584787949762E-04_wp, &
      &-3.39826415555430E-17_wp, 1.00000000000000E+00_wp, 1.23625096714642E-33_wp, &
      &-6.58317486018941E-33_wp,-4.32688697434234E-16_wp,-1.05751403193283E-17_wp, &
      &-1.61845529916613E-17_wp, 2.04001046366190E-48_wp,-2.80324680793483E-17_wp, &
      &-1.52597050519540E-02_wp, 7.54386858235250E-04_wp, 1.54483414153759E-01_wp, &
      &-1.38975883295297E-02_wp, 3.68979292051442E-03_wp,-9.29593376498922E-03_wp, &
      & 3.89999669170718E-02_wp, 1.71253063589435E-01_wp,-3.38872073617419E-02_wp, &
      & 2.75607121079114E-03_wp,-8.49668395531560E-03_wp,-1.77101926740922E-03_wp, &
      &-5.96497855969987E-02_wp, 3.26263593187400E-02_wp, 7.80338889697969E-02_wp, &
      &-1.67913351825952E-02_wp,-2.72684695870236E-03_wp,-2.72869099621408E-03_wp, &
      &-2.51676640870696E-03_wp,-7.54506460969072E-03_wp,-6.95175841802704E-04_wp, &
      &-2.84118041224261E-04_wp, 4.85651186285602E-04_wp,-7.85610749179761E-04_wp, &
      &-8.35404391158209E-04_wp,-5.87037410579942E-06_wp, 1.14615352525398E-05_wp, &
      & 2.03515403058432E-03_wp,-1.63001650308425E-05_wp, 2.91492178236114E-06_wp, &
      &-2.63759193622194E-04_wp, 2.36053787840104E-06_wp, 3.75108420444024E-04_wp, &
      &-1.04792588750760E-06_wp, 1.71756912647956E-17_wp, 1.23625096714642E-33_wp, &
      & 1.00000000000000E+00_wp, 3.36605816337719E-33_wp, 2.04001046366190E-48_wp, &
      & 2.80324680793483E-17_wp,-1.22111202201645E-17_wp,-4.32688697434234E-16_wp, &
      & 0.00000000000000E+00_wp,-1.59163304337980E-01_wp, 7.86848138263948E-03_wp, &
      &-1.38975883295297E-02_wp, 1.08598206891990E-02_wp, 2.88359568781052E-02_wp, &
      & 3.68979292051442E-03_wp, 9.88745071119326E-02_wp,-5.09310803815546E-02_wp, &
      &-1.75683639506761E-01_wp, 1.78113863857822E-02_wp, 3.64904446403644E-03_wp, &
      &-4.35489790425279E-03_wp, 1.98257879724107E-02_wp, 5.50995923856506E-02_wp, &
      &-2.56047302565900E-02_wp,-7.47087628118558E-02_wp,-4.49381885326583E-03_wp, &
      &-4.49685781022301E-03_wp,-7.54506460969073E-03_wp,-1.03726225848758E-02_wp, &
      &-1.00792913004316E-03_wp,-6.95175841802704E-04_wp, 1.40797881474581E-04_wp, &
      &-1.69115652064343E-03_wp,-9.95946946791490E-04_wp,-4.44405940448674E-03_wp, &
      & 8.67674574241973E-03_wp,-1.63001650308424E-05_wp,-1.03045664021254E-02_wp, &
      & 1.94292520555368E-03_wp, 2.91492178236113E-06_wp, 1.13728572222362E-03_wp, &
      &-3.64998611613997E-06_wp,-4.18199240292262E-04_wp, 1.32489189884931E-16_wp, &
      &-6.58317486018941E-33_wp, 3.36605816337719E-33_wp, 1.00000000000000E+00_wp, &
      & 2.80324680793483E-17_wp, 2.04001046366190E-48_wp, 2.49812935938963E-16_wp, &
      &-1.05751403193283E-17_wp,-4.32688697434233E-16_wp,-1.52943032490713E-03_wp, &
      &-1.95185887127872E-02_wp, 7.45576160756131E-05_wp, 1.83487211019493E-03_wp, &
      &-9.99777769247366E-02_wp,-1.18826943267016E-02_wp, 3.99357968717228E-03_wp, &
      &-3.35887729018430E-03_wp,-2.08808975488039E-02_wp,-7.97563207247439E-02_wp, &
      &-1.51100766908978E-02_wp, 1.98701368923876E-02_wp, 1.52776930829613E-03_wp, &
      & 5.00961158918477E-03_wp,-3.93928241184114E-03_wp,-2.11528101890322E-02_wp, &
      & 2.09708745158415E-04_wp, 3.35340135917542E-05_wp, 6.95206092026806E-04_wp, &
      & 1.00797286939984E-03_wp, 1.14546171462770E-06_wp, 7.26736963500946E-07_wp, &
      &-3.15498276364994E-07_wp, 2.88230422225898E-06_wp, 2.15843791949334E-06_wp, &
      &-4.02364515623123E-04_wp, 1.17681382581479E-03_wp,-2.91563193175748E-06_wp, &
      &-1.94339728637851E-03_wp, 5.70674555842992E-06_wp, 7.87968341926241E-09_wp, &
      & 3.54611710428893E-06_wp,-1.17502876050915E-08_wp,-2.35045189264994E-06_wp, &
      &-1.41339118481328E-32_wp,-4.32688697434234E-16_wp, 2.04001046366190E-48_wp, &
      & 2.80324680793483E-17_wp, 1.00000000000000E+00_wp,-4.39880888958495E-32_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.46633394881259E-04_wp,-1.87133528061882E-03_wp, 1.06436196936872E-03_wp, &
      & 7.45576160756129E-05_wp,-1.18826943267016E-02_wp, 2.28230442519527E-02_wp, &
      & 1.52057505474754E-03_wp, 6.46610665419611E-03_wp,-2.65879307064732E-03_wp, &
      &-1.07945544482274E-02_wp, 2.51138516221906E-02_wp, 5.96903307305780E-03_wp, &
      &-2.53958702119006E-03_wp, 2.13969759880309E-03_wp, 6.12494819781462E-03_wp, &
      &-2.80263201427949E-03_wp, 1.27251158228804E-04_wp, 2.03484221241593E-05_wp, &
      & 2.84130284308594E-04_wp, 6.95206092026806E-04_wp, 7.26736963500946E-07_wp, &
      & 3.88789458191834E-07_wp,-4.18923418785128E-07_wp, 1.61223564290462E-06_wp, &
      & 1.44107420785261E-06_wp,-5.31502848774173E-07_wp, 1.55451059079795E-06_wp, &
      & 2.63824717478037E-04_wp,-2.91563193175748E-06_wp, 7.87968341926250E-09_wp, &
      &-2.58414322342849E-07_wp, 5.60586440485871E-09_wp, 4.02801096758346E-07_wp, &
      &-3.63692714880034E-09_wp, 1.15788362790797E-33_wp,-1.05751403193283E-17_wp, &
      & 2.80324680793483E-17_wp, 2.04001046366190E-48_wp,-4.39880888958495E-32_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp,-8.00806831828508E-03_wp,-8.31407752597991E-04_wp, &
      &-1.76577636714531E-03_wp, 1.53164951820822E-02_wp, 3.99357968717228E-03_wp, &
      & 1.52057505474754E-03_wp, 4.29241979763334E-02_wp,-2.80125851933167E-02_wp, &
      &-3.66771973089427E-02_wp, 3.46636784544117E-03_wp, 2.24585339359639E-03_wp, &
      &-2.28634033709060E-03_wp, 1.93441892736988E-02_wp, 4.21197909788809E-02_wp, &
      &-2.06260083674999E-02_wp,-2.11188100405667E-02_wp,-6.60162410569535E-05_wp, &
      &-5.09220431996983E-05_wp,-4.85672551966368E-04_wp,-1.40803816383847E-04_wp, &
      &-3.15498276364994E-07_wp,-4.18923418785128E-07_wp,-3.15321837140630E-07_wp, &
      &-1.15835917868009E-06_wp,-3.79139535106234E-07_wp,-2.46860584227914E-04_wp, &
      & 7.99881727160953E-04_wp,-2.36111571061122E-06_wp,-1.13756175508771E-03_wp, &
      & 3.54611710428893E-06_wp, 5.60586440485866E-09_wp, 2.20044667218825E-06_wp, &
      &-7.97244984431500E-09_wp,-1.27484309740475E-06_wp, 0.00000000000000E+00_wp, &
      &-1.61845529916613E-17_wp,-1.22111202201645E-17_wp, 2.49812935938963E-16_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      &-2.53965349651559E-32_wp, 0.00000000000000E+00_wp, 2.70133358657425E-03_wp, &
      & 7.45576160756131E-05_wp,-1.96080622593259E-02_wp,-3.24081560148148E-03_wp, &
      &-3.35887729018429E-03_wp, 6.46610665419610E-03_wp,-2.80125851933167E-02_wp, &
      &-9.59469256266461E-02_wp, 2.47798997497724E-02_wp,-3.63018922687454E-03_wp, &
      & 8.10425997898489E-03_wp, 2.13969759880309E-03_wp, 4.67851932496283E-02_wp, &
      &-3.33331054393402E-02_wp,-6.15774511801703E-02_wp, 1.83948733471546E-02_wp, &
      & 3.51860365217760E-04_wp, 6.95206092026806E-04_wp, 7.85644602357878E-04_wp, &
      & 1.69122990883716E-03_wp, 2.88230422225898E-06_wp, 1.61223564290462E-06_wp, &
      &-1.15835917868009E-06_wp, 4.26369032405522E-06_wp, 3.25838762232018E-06_wp, &
      & 7.55883392450601E-07_wp,-2.91563193175748E-06_wp,-3.75201606011243E-04_wp, &
      & 3.65087296883579E-06_wp,-1.17502876050916E-08_wp, 4.02801096758346E-07_wp, &
      &-7.97244984431508E-09_wp,-5.48031491053023E-07_wp, 3.65883314752589E-09_wp, &
      & 2.56643374397691E-32_wp, 2.04001046366190E-48_wp,-4.32688697434234E-16_wp, &
      &-1.05751403193283E-17_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-2.53965349651559E-32_wp, 1.00000000000000E+00_wp,-4.39880888958495E-32_wp, &
      & 1.40463499393994E-02_wp,-6.69530328480171E-04_wp,-6.84740160431334E-04_wp, &
      &-2.66184213957890E-02_wp,-2.08808975488039E-02_wp,-2.65879307064732E-03_wp, &
      &-3.66771973089427E-02_wp, 2.47798997497724E-02_wp, 8.95196235017786E-02_wp, &
      &-2.06584952521769E-02_wp,-3.93928241184114E-03_wp, 4.14067960720862E-03_wp, &
      &-1.40310950956885E-02_wp,-3.58779170386596E-02_wp, 2.06395253003514E-02_wp, &
      & 7.06498716627562E-02_wp, 2.52010492935596E-04_wp, 6.35643019044151E-04_wp, &
      & 8.35440743356501E-04_wp, 9.95989952642546E-04_wp, 2.15843791949335E-06_wp, &
      & 1.44107420785261E-06_wp,-3.79139535106234E-07_wp, 3.25838762232018E-06_wp, &
      & 1.94316371017185E-06_wp, 1.44651631713556E-04_wp,-8.21786212367834E-04_wp, &
      & 1.04818118901892E-06_wp, 4.18299324733108E-04_wp,-2.35045189264994E-06_wp, &
      &-3.63692714880030E-09_wp,-1.27484309740475E-06_wp, 3.65883314752584E-09_wp, &
      & 1.37003353935441E-08_wp, 0.00000000000000E+00_wp,-2.80324680793483E-17_wp, &
      & 0.00000000000000E+00_wp,-4.32688697434233E-16_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,-4.39880888958495E-32_wp, &
      & 1.00000000000000E+00_wp],shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "AcCl3")
   call test_qvszp_overlap(error, mol, .true., overlap, thr_in=thr1)

end subroutine test_qvszp_overlap_accl3

end module test_qvszp