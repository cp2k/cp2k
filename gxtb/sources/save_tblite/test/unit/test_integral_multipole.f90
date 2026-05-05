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

module test_integral_multipole
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use multicharge, only : get_eeqbc_charges
   use tblite_basis_cache, only : basis_cache, cgto_cache
   use tblite_basis_qvszp, only : qvszp_basis_type, qvszp_cgto_type, &
      & new_qvszp_cgto, new_qvszp_basis
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_basis_type, only : basis_type, new_basis, cgto_container, cgto_type, &
      & get_cutoff 
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_dipole, only : dipole_cgto, dipole_grad_cgto, msao, &
      & get_dipole_integrals
   use tblite_integral_multipole, only : multipole_cgto, multipole_grad_cgto, &
      & get_multipole_integrals
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   implicit none
   private

   public :: collect_integral_multipole

   !real(wp), parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine basis_maker(bas, mol, error)
         import :: basis_type, structure_type, error_type
         class(basis_type), allocatable, intent(out) :: bas
         type(structure_type), intent(in) :: mol
         type(error_type), allocatable, intent(out) :: error
      end subroutine basis_maker
   end interface

contains


!> Collect all exported unit tests
subroutine collect_integral_multipole(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("dipole-trans-ss", test_dipole_ss), &
      new_unittest("dipole-trans-pp", test_dipole_pp), &
      new_unittest("dipole-trans-dd", test_dipole_dd), &
      new_unittest("dipole-grad-ss", test_dipole_grad_ss), &
      new_unittest("dipole-grad-qeff-ss", test_dipole_grad_qeff_ss), &
      new_unittest("dipole-grad-qeff-sp", test_dipole_grad_qeff_sp), &
      new_unittest("dipole-grad-qeff-pp", test_dipole_grad_qeff_pp), &
      new_unittest("dipole-grad-qeff-sd", test_dipole_grad_qeff_sd), &
      new_unittest("dipole-grad-qeff-pd", test_dipole_grad_qeff_pd), &
      new_unittest("dipole-grad-qeff-dd", test_dipole_grad_qeff_dd), &
      new_unittest("dipole-grad-qeff-sf", test_dipole_grad_qeff_sf), &
      new_unittest("dipole-grad-qeff-pf", test_dipole_grad_qeff_pf), &
      new_unittest("dipole-grad-qeff-df", test_dipole_grad_qeff_df), &
      new_unittest("dipole-grad-qeff-ff", test_dipole_grad_qeff_ff), &
      new_unittest("multipole-grad-qeff-ss", test_multipole_grad_qeff_ss), &
      new_unittest("multipole-grad-qeff-sp", test_multipole_grad_qeff_sp), &
      new_unittest("multipole-grad-qeff-pp", test_multipole_grad_qeff_pp), &
      new_unittest("multipole-grad-qeff-sd", test_multipole_grad_qeff_sd), &
      new_unittest("multipole-grad-qeff-pd", test_multipole_grad_qeff_pd), &
      new_unittest("multipole-grad-qeff-dd", test_multipole_grad_qeff_dd), &
      new_unittest("multipole-grad-qeff-sf", test_multipole_grad_qeff_sf), &
      new_unittest("multipole-grad-qeff-pf", test_multipole_grad_qeff_pf), &
      new_unittest("multipole-grad-qeff-df", test_multipole_grad_qeff_df), &
      new_unittest("multipole-grad-qeff-ff", test_multipole_grad_qeff_ff), &
      new_unittest("overlap-dipole-diat-alh3", test_overlap_dipole_diat_alh3), &
      new_unittest("overlap-multipole-diat-alh3", test_overlap_multipole_diat_alh3), &
      new_unittest("overlap-dipole-qvszp-alh3", test_overlap_dipole_qvszp_alh3), &
      new_unittest("overlap-multipole-qvszp-alh3", test_overlap_multipole_qvszp_alh3) &
      ]

end subroutine collect_integral_multipole


subroutine make_gen_basis(bas, mol, error)
   class(basis_type), allocatable, intent(out) :: bas
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   integer, parameter :: nsh(20) = [&
      & 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]
   integer, parameter :: lsh(3, 20) = reshape([&
      & 0, 0, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, &
      & 0, 1, 0,  0, 1, 0,  0, 1, 2,  0, 1, 0,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
      & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 0,  0, 1, 2], &
      & shape(lsh))
   integer, parameter :: pqn(3, 20) = reshape([&
      & 1, 0, 0,  1, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0, &
      & 2, 2, 0,  2, 2, 0,  2, 2, 3,  3, 3, 0,  3, 3, 3,  3, 3, 3,  3, 3, 3, &
      & 3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  4, 4, 0,  4, 4, 3], &
      & shape(pqn))
   real(wp), parameter :: zeta(3, 20) = reshape([&
      & 1.230000_wp, 0.000000_wp, 0.000000_wp, 1.669667_wp, 1.500000_wp, 0.000000_wp, &
      & 0.750060_wp, 0.557848_wp, 0.000000_wp, 1.034720_wp, 0.949332_wp, 0.000000_wp, &
      & 1.479444_wp, 1.479805_wp, 0.000000_wp, 2.096432_wp, 1.800000_wp, 0.000000_wp, &
      & 2.339881_wp, 2.014332_wp, 0.000000_wp, 2.439742_wp, 2.137023_wp, 0.000000_wp, &
      & 2.416361_wp, 2.308399_wp, 0.000000_wp, 3.084104_wp, 2.312051_wp, 2.815609_wp, &
      & 0.763787_wp, 0.573553_wp, 0.000000_wp, 1.184203_wp, 0.717769_wp, 1.300000_wp, &
      & 1.352531_wp, 1.391201_wp, 1.000000_wp, 1.773917_wp, 1.718996_wp, 1.250000_wp, &
      & 1.816945_wp, 1.903247_wp, 1.167533_wp, 1.981333_wp, 2.025643_wp, 1.702555_wp, &
      & 2.485265_wp, 2.199650_wp, 2.476089_wp, 2.329679_wp, 2.149419_wp, 1.950531_wp, &
      & 0.875961_wp, 0.631694_wp, 0.000000_wp, 1.267130_wp, 0.786247_wp, 1.380000_wp],&
      & shape(zeta))

   integer :: isp, izp, ish, stat
   integer, allocatable :: nshell(:)
   type(cgto_container), allocatable :: cgto(:, :)

   nshell = nsh(mol%num)
   allocate(cgto(maxval(nshell), mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, nshell(isp)
         allocate(cgto_type :: cgto(ish, isp)%raw)
         call slater_to_gauss(ng, pqn(ish, izp), lsh(ish, izp), zeta(ish, izp), &
            & cgto(ish, isp)%raw, .true., stat)
      end do
   end do

   allocate(bas)
   call new_basis(bas, mol, nshell, cgto, accuracy=1.0_wp)

end subroutine make_gen_basis

subroutine make_qvszp_basis(bas, mol, error)
   class(basis_type), allocatable, intent(out) :: bas
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error

   !> Parameter: Number of shells selected from the q-vSZP basis set
   integer, parameter :: pa_nshell(60) = [ &
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, & !1-20
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, & !21-40
      & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 4, 4, 4] !41-60

   integer :: isp, izp, ish
   integer, allocatable :: nshell(:)
   type(cgto_container), allocatable :: cgto(:, :)
   type(cgto_container), allocatable :: cgto_h0(:, :)
   type(qvszp_cgto_type), allocatable :: cgto_qvszp
   type(qvszp_basis_type), allocatable :: basis_qvzsp
   real(wp) :: alpha(12), k0, k2, k3

   nshell = pa_nshell(mol%num)
   allocate(cgto(maxval(nshell), mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, nshell(isp)
         allocate(cgto_qvszp)
         call new_qvszp_cgto(cgto_qvszp, izp, ish, .true., error)
         if (allocated(error)) return
         call move_alloc(cgto_qvszp, cgto(ish, isp)%raw)
      end do
   end do

   allocate(basis_qvzsp)
   call new_qvszp_basis(basis_qvzsp, mol, nshell, cgto, error, &
      & accuracy=0.1_wp, cgto_h0=cgto_h0)
   if (allocated(error)) return
   call move_alloc(basis_qvzsp, bas)

end subroutine make_qvszp_basis

subroutine test_dipole_ss(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat, i
   type(cgto_type) :: cgtoi, cgtoj
   type(cgto_cache) :: icache, jcache
   real(wp) :: vec(3), r2
   real(wp) :: overlap(1, 1), dipolei(3, 1, 1), dipolej(3, 1, 1)

   call slater_to_gauss(ng, 2, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 1, 0, 1.0_wp, cgtoj, .true., stat)

   call cgtoi%update(icache, .false.)
   call cgtoj%update(jcache, .false.)

   call random_number(vec)
   vec = vec - 0.5_wp
   r2 = sum(vec**2)

   call dipole_cgto(cgtoi, cgtoj, icache, jcache, r2, vec, 100.0_wp, overlap, dipolej)

   vec(:) = -vec

   call dipole_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, overlap, dipolei)

   do i = 1, 3
      call check(error, dipolei(i, 1, 1) + vec(i) * overlap(1, 1), dipolej(i, 1, 1), thr=thr)
      if (allocated(error)) return
   end do

end subroutine test_dipole_ss


subroutine test_dipole_pp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6, n = 2

   integer :: stat, i, j, k
   type(cgto_type) :: cgtoi, cgtoj
   type(cgto_cache) :: icache, jcache
   real(wp) :: vec(3), r2
   real(wp) :: overlap(3, 3), dipolei(3, 3, 3), dipolej(3, 3, 3)

   call slater_to_gauss(ng, 2, 1, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 1, 1.0_wp, cgtoj, .true., stat)

   call cgtoi%update(icache, .false.)
   call cgtoj%update(jcache, .false.)

   call random_number(vec)
   vec = vec - 0.5_wp
   r2 = sum(vec**2)

   call dipole_cgto(cgtoi, cgtoj, icache, jcache, r2, vec, 100.0_wp, overlap, dipolej)

   vec(:) = -vec

   call dipole_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, overlap, dipolei)

   do i = 1, 3
      do j = 1, 3
         do k = 1, 3
            call check(error, dipolei(k, j, i) + vec(k) * overlap(j, i), dipolej(k, i, j), thr=thr)
            if (allocated(error)) return
         end do
      end do
   end do

end subroutine test_dipole_pp


subroutine test_dipole_dd(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6, n = 2

   integer :: stat, i, j, k
   type(cgto_type) :: cgtoi, cgtoj
   type(cgto_cache) :: icache, jcache
   real(wp) :: vec(3), r2
   real(wp) :: overlap(5, 5), dipolei(3, 5, 5), dipolej(3, 5, 5)

   call slater_to_gauss(ng, 3, 2, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 4, 2, 1.0_wp, cgtoj, .true., stat)

   call random_number(vec)
   vec = vec - 0.5_wp
   r2 = sum(vec**2)

   call dipole_cgto(cgtoi, cgtoj, icache, jcache, r2, vec, 100.0_wp, overlap, dipolej)

   vec(:) = -vec

   call dipole_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, overlap, dipolei)

   do i = 1, 5
      do j = 1, 5
         do k = 1, 3
            call check(error, dipolei(k, j, i) + vec(k) * overlap(j, i), dipolej(k, i, j), thr=thr)
            if (allocated(error)) return
         end do
      end do
   end do

end subroutine test_dipole_dd


subroutine test_dipole_grad_ss(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat, i, j
   type(cgto_type) :: cgtoi, cgtoj
   type(cgto_cache) :: icache, jcache
   real(wp) :: vec(3), r2, zero(3)
   real(wp) :: overlap(1, 1), doverlapi(3, 1, 1)
   real(wp) :: dipole(3, 1, 1), ddipolei(3, 3, 1, 1), ddipolej(3, 3, 1, 1)
   real(wp) :: quadrupole(6, 1, 1), dquadrupolei(3, 6, 1, 1), dquadrupolej(3, 6, 1, 1)
   real(wp) :: sr(1, 1), sl(1, 1), dr(3, 1, 1), dl(3, 1, 1), qr(6, 1, 1), ql(6, 1, 1)
   real(wp), parameter :: step = 1.0e-5_wp

   call slater_to_gauss(ng, 2, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 1, 0, 1.0_wp, cgtoj, .true., stat)

   call cgtoi%update(icache, .true.)
   call cgtoj%update(jcache, .true.)

   call random_number(vec)
   vec = vec - 0.5_wp
   zero = 0
   r2 = sum(vec**2)

   call multipole_grad_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & overlap, dipole, quadrupole, doverlapi, ddipolej, dquadrupolej, &
      & ddipolei, dquadrupolei)

   do i = 1, 3
      vec(i) = vec(i) + step
      r2 = sum(vec**2)
      call multipole_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, sr, dr, qr)
      vec(i) = vec(i) - 2*step
      r2 = sum(vec**2)
      call multipole_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, sl, dl, ql)
      vec(i) = vec(i) + step
      ddipolej(i, :, :, :) = 0.5_wp * (dr - dl) / step
      dquadrupolej(i, :, :, :) = 0.5_wp * (qr - ql) / step
   end do

   num: do i = 1, 3
      do j = 1, 3
         call check(error, ddipolei(i, j, 1, 1), ddipolej(i, j, 1, 1), thr=thr1*10)
         if (allocated(error)) exit num
      end do
      do j = 1, 6
         call check(error, dquadrupolei(i, j, 1, 1), dquadrupolej(i, j, 1, 1), thr=thr1*10)
         if (allocated(error)) exit num
      end do
   end do num
   if (allocated(error)) return

end subroutine test_dipole_grad_ss


!> General routine to test charge gradients of dipole integrals for q-vSZP CGTO pairs
subroutine test_dipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in)
   !> Vector between the two centers
   real(wp), intent(inout) :: vec(3)
   !> CGTO for center i
   type(qvszp_cgto_type), intent(in) :: cgtoi
   !> CGTO for center j
   type(qvszp_cgto_type), intent(in) :: cgtoj
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: naoi, naoj
   real(wp) :: r2, cn, q, qeffi, qeffj, thr_
   type(cgto_cache) :: icache, jcache
   real(wp), allocatable :: overlap(:, :), dpint(:, :, :)
   real(wp), allocatable :: doverlap(:, :, :), ddpint(:, :, :, :)
   real(wp), allocatable :: doverlapdqeffi(:, :), doverlapdqeffj(:, :)
   real(wp), allocatable :: ddpintdqeffi(:, :, :), ddpintdqeffj(:, :, :)
   real(wp), allocatable :: sr(:, :), sl(:, :)
   real(wp), allocatable :: dr(:, :, :), dl(:, :, :)
   real(wp), allocatable :: num_dSdqi(:, :), num_dSdqj(:, :)
   real(wp), allocatable :: num_dDdqi(:, :, :), num_dDdqj(:, :, :)
   real(wp), parameter :: step = 1.0e-5_wp

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Setup CN and charge for a neutral dimer
   cn = 1.0_wp
   q = -0.4_wp
   call cgtoi%update(icache, .true., cn, q)
   
   cn = 1.0_wp
   q = 0.4_wp
   call cgtoj%update(jcache, .true., cn, q)

   ! Analytic gradient calculation
   r2 = sum(vec**2)
   naoi = msao(cgtoi%ang)
   naoj = msao(cgtoj%ang)
   allocate(overlap(naoj, naoi), dpint(3, naoj, naoi), doverlap(3, naoj, naoi), &
      & ddpint(3, 3, naoj, naoi), doverlapdqeffi(naoj, naoi), doverlapdqeffj(naoj, naoi), &
      & ddpintdqeffi(3, naoj, naoi), ddpintdqeffj(3, naoj, naoi))
   call dipole_grad_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & overlap, dpint, doverlap, ddpint, &
      & doverlapdqeffj, doverlapdqeffi, ddpintdqeffj, ddpintdqeffi)

   ! Numerical gradient calculation
   allocate(sr(naoj, naoi), sl(naoj, naoi), dr(3, naoj, naoi), dl(3, naoj, naoi), &
      & num_dDdqi(3, naoj, naoi), num_dDdqj(3, naoj, naoi), &
      & num_dSdqi(naoj, naoi), num_dSdqj(naoj, naoi))
   
   ! Right hand side
   sr = 0.0_wp
   dr = 0.0_wp
   qeffi = icache%qeff
   icache%qeff = qeffi + step
   call cgtoi%get_normalization(icache, .false.)
   call dipole_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & sr, dr)

   ! Left hand side
   sl = 0.0_wp
   dl = 0.0_wp
   icache%qeff = qeffi - step
   call cgtoi%get_normalization(icache, .false.)
   call dipole_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & sl, dl)

   icache%qeff = qeffi
   call cgtoi%get_normalization(icache, .false.)
   num_dSdqi = 0.5_wp * (sr - sl) / step
   num_dDdqi = 0.5_wp * (dr - dl) / step

   ! Right hand side
   sr = 0.0_wp
   dr = 0.0_wp
   qeffj = jcache%qeff
   jcache%qeff = qeffj + step
   call cgtoj%get_normalization(jcache, .false.)
   call dipole_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & sr, dr)

   ! Left hand side
   sl = 0.0_wp
   dl = 0.0_wp
   jcache%qeff = qeffj - step
   call cgtoj%get_normalization(jcache, .false.)
   call dipole_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & sl, dl)

   jcache%qeff = qeffj
   call cgtoj%get_normalization(jcache, .false.)
   num_dSdqj = 0.5_wp * (sr - sl) / step
   num_dDdqj = 0.5_wp * (dr - dl) / step

   if (any(abs(ddpintdqeffi - num_dDdqi) > thr_)) then
      call test_failed(error, "Charge gradient of dipole integrals at center i does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', ddpintdqeffi
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dDdqi
      write(*,*) 'Difference:'
      print'(3es21.14)', ddpintdqeffi - num_dDdqi
   end if

   if (any(abs(ddpintdqeffj - num_dDdqj) > thr_)) then
      call test_failed(error, "Charge gradient of dipole integrals at center j does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', ddpintdqeffj
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dDdqj
      write(*,*) 'Difference:'
      print'(3es21.14)', ddpintdqeffj - num_dDdqj
   end if

   if (any(abs(doverlapdqeffi - num_dSdqi) > thr_)) then
      call test_failed(error, "Charge gradient of overlap at center i does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', doverlapdqeffi
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dSdqi
      write(*,*) 'Difference:'
      print'(3es21.14)', doverlapdqeffi - num_dSdqi
   end if

   if (any(abs(doverlapdqeffj - num_dSdqj) > thr_)) then
      call test_failed(error, "Charge gradient of overlap at center j does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', doverlapdqeffj
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dSdqj
      write(*,*) 'Difference:'
      print'(3es21.14)', doverlapdqeffj - num_dSdqj
   end if

end subroutine test_dipole_grad_qeff_gen


subroutine test_dipole_grad_qeff_ss(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 1, 1, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 2, 1, .true., error)
   if (allocated(error)) return

   call test_dipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_dipole_grad_qeff_ss

subroutine test_dipole_grad_qeff_sp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 3, 1, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 4, 2, .true., error)
   if (allocated(error)) return

   call test_dipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1*10)
end subroutine test_dipole_grad_qeff_sp

subroutine test_dipole_grad_qeff_pp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 5, 2, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 6, 2, .true., error)
   if (allocated(error)) return

   call test_dipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_dipole_grad_qeff_pp

subroutine test_dipole_grad_qeff_sd(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 7, 1, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 8, 3, .true., error)
   if (allocated(error)) return

   call test_dipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_dipole_grad_qeff_sd

subroutine test_dipole_grad_qeff_pd(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 9, 2, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 10, 3, .true., error)
   if (allocated(error)) return

   call test_dipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_dipole_grad_qeff_pd

subroutine test_dipole_grad_qeff_dd(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 11, 3, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 12, 3, .true., error)
   if (allocated(error)) return

   call test_dipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1*10)
end subroutine test_dipole_grad_qeff_dd

subroutine test_dipole_grad_qeff_sf(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 13, 1, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 59, 4, .true., error)
   if (allocated(error)) return

   call test_dipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_dipole_grad_qeff_sf

subroutine test_dipole_grad_qeff_pf(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 14, 2, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 60, 4, .true., error)
   if (allocated(error)) return

   call test_dipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_dipole_grad_qeff_pf

subroutine test_dipole_grad_qeff_df(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 15, 3, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 61, 4, .true., error)
   if (allocated(error)) return

   call test_dipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_dipole_grad_qeff_df

subroutine test_dipole_grad_qeff_ff(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 62, 4, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 63, 4, .true., error)
   if (allocated(error)) return

   call test_dipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_dipole_grad_qeff_ff


!> General routine to test charge gradients of multipole integrals for q-vSZP CGTO pairs
subroutine test_multipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in)
   !> Vector between the two centers
   real(wp), intent(inout) :: vec(3)
   !> CGTO for center i
   type(qvszp_cgto_type), intent(in) :: cgtoi
   !> CGTO for center j
   type(qvszp_cgto_type), intent(in) :: cgtoj
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: naoi, naoj, iao, jao
   real(wp) :: r2, cn, q, qeffi, qeffj, thr_
   type(cgto_cache) :: icache, jcache
   real(wp), allocatable :: overlap(:, :), dpint(:, :, :), qpint(:, :, :)
   real(wp), allocatable :: doverlap(:, :, :), ddpintj(:, :, :, :), ddpinti(:, :, :, :)
   real(wp), allocatable :: dqpintj(:, :, :, :), dqpinti(:, :, :, :)
   real(wp), allocatable :: doverlapdqeffi(:, :), doverlapdqeffj(:, :)
   real(wp), allocatable :: ddpintdqeffi(:, :, :), ddpintdqeffj(:, :, :)
   real(wp), allocatable :: dqpintdqeffi(:, :, :), dqpintdqeffj(:, :, :)
   real(wp), allocatable :: sr(:, :), sl(:, :)
   real(wp), allocatable :: dr(:, :, :), dl(:, :, :), dtmp(:, :, :)
   real(wp), allocatable :: qr(:, :, :), ql(:, :, :), qtmp(:, :, :)
   real(wp), allocatable :: num_dSdqi(:, :), num_dSdqj(:, :)
   real(wp), allocatable :: num_dDdqi(:, :, :), num_dDdqj(:, :, :)
   real(wp), allocatable :: num_dQdqi(:, :, :), num_dQdqj(:, :, :)
   real(wp), parameter :: step = 1.0e-5_wp

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Setup CN and charge for a neutral dimer
   cn = 1.0_wp
   q = -0.4_wp
   call cgtoi%update(icache, .true., cn, q)
   
   cn = 1.0_wp
   q = 0.4_wp
   call cgtoj%update(jcache, .true., cn, q)

   ! Analytic gradient calculation
   r2 = sum(vec**2)
   naoi = msao(cgtoi%ang)
   naoj = msao(cgtoj%ang)
   allocate(overlap(naoj, naoi), dpint(3, naoj, naoi), qpint(6, naoj, naoi), &
      & doverlap(3, naoj, naoi), ddpintj(3, 3, naoj, naoi), ddpinti(3, 3, naoj, naoi), &
      & dqpintj(3, 6, naoj, naoi), dqpinti(3, 6, naoj, naoi), &
      & doverlapdqeffi(naoj, naoi), doverlapdqeffj(naoj, naoi), &
      & ddpintdqeffi(3, naoj, naoi), ddpintdqeffj(3, naoj, naoi), &
      & dqpintdqeffi(6, naoj, naoi), dqpintdqeffj(6, naoj, naoi))
   call multipole_grad_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & overlap, dpint, qpint, doverlap, ddpintj, dqpintj, ddpinti, dqpinti, &
      & doverlapdqeffj, doverlapdqeffi, ddpintdqeffj, ddpintdqeffi, &
      & dqpintdqeffj, dqpintdqeffi)

   ! Numerical gradient calculation
   allocate(sr(naoj, naoi), sl(naoj, naoi), dr(3, naoj, naoi), dl(3, naoj, naoi), &
      & dtmp(3, naoj, naoi), qr(6, naoj, naoi), ql(6, naoj, naoi), qtmp(6, naoj, naoi), &
      & num_dSdqi(naoj, naoi), num_dSdqj(naoj, naoi), num_dDdqi(3, naoj, naoi), &
      & num_dDdqj(3, naoj, naoi), num_dQdqi(6, naoj, naoi), num_dQdqj(6, naoj, naoi))
   
   ! Right hand side
   sr = 0.0_wp
   dr = 0.0_wp
   dtmp = 0.0_wp
   qr = 0.0_wp
   qtmp = 0.0_wp
   qeffi = icache%qeff
   icache%qeff = qeffi + step
   call cgtoi%get_normalization(icache, .false.)
   call multipole_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & sr, dtmp, qtmp)
   ! Shift the atom i centered multipole integrals to atom j to align with gradient
   do jao = 1, naoj
      do iao = 1, naoi
         call shift_operator(vec, sr(jao, iao), dtmp(:, jao, iao), qtmp(:, jao, iao), &
            & dr(:, jao, iao), qr(:, jao, iao))
      end do
   end do

   ! Left hand side
   sl = 0.0_wp
   dl = 0.0_wp
   dtmp = 0.0_wp
   ql = 0.0_wp
   qtmp = 0.0_wp
   icache%qeff = qeffi - step
   call cgtoi%get_normalization(icache, .false.)
   call multipole_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & sl, dtmp, qtmp)
   ! Shift the atom i centered multipole integrals to atom j to align with gradient
   do jao = 1, naoj
      do iao = 1, naoi
         call shift_operator(vec, sl(jao, iao), dtmp(:, jao, iao), qtmp(:, jao, iao), &
            & dl(:, jao, iao), ql(:, jao, iao))
      end do
   end do

   icache%qeff = qeffi
   call cgtoi%get_normalization(icache, .false.)
   num_dSdqi = 0.5_wp * (sr - sl) / step
   num_dDdqi = 0.5_wp * (dr - dl) / step
   num_dQdqi = 0.5_wp * (qr - ql) / step

   ! Right hand side
   sr = 0.0_wp
   dr = 0.0_wp
   dtmp = 0.0_wp
   qr = 0.0_wp
   qtmp = 0.0_wp
   qeffj = jcache%qeff
   jcache%qeff = qeffj + step
   call cgtoj%get_normalization(jcache, .false.)
   call multipole_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & sr, dtmp, qtmp)
   ! Shift the atom i centered multipole integrals to atom j to align with gradient
   do jao = 1, naoj
      do iao = 1, naoi
         call shift_operator(vec, sr(jao, iao), dtmp(:, jao, iao), qtmp(:, jao, iao), &
            & dr(:, jao, iao), qr(:, jao, iao))
      end do
   end do   

   ! Left hand side
   sl = 0.0_wp
   dl = 0.0_wp
   dtmp = 0.0_wp
   ql = 0.0_wp
   qtmp = 0.0_wp
   jcache%qeff = qeffj - step
   call cgtoj%get_normalization(jcache, .false.)
   call multipole_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & sl, dtmp, qtmp)
   ! Shift the atom i centered multipole integrals to atom j to align with gradient
   do jao = 1, naoj
      do iao = 1, naoi
         call shift_operator(vec, sl(jao, iao), dtmp(:, jao, iao), qtmp(:, jao, iao), &
            & dl(:, jao, iao), ql(:, jao, iao))
      end do
   end do

   jcache%qeff = qeffj
   call cgtoj%get_normalization(jcache, .false.)
   num_dSdqj = 0.5_wp * (sr - sl) / step
   num_dDdqj = 0.5_wp * (dr - dl) / step
   num_dQdqj = 0.5_wp * (qr - ql) / step

   if (any(abs(doverlapdqeffi - num_dSdqi) > thr_)) then
      write(*,*) "here 1"
      call test_failed(error, "Charge gradient of overlap at center i does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', doverlapdqeffi
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dSdqi
      write(*,*) 'Difference:'
      print'(3es21.14)', doverlapdqeffi - num_dSdqi
   end if

   if (any(abs(doverlapdqeffj - num_dSdqj) > thr_)) then
      write(*,*) "here 2"
      call test_failed(error, "Charge gradient of overlap at center j does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', doverlapdqeffj
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dSdqj
      write(*,*) 'Difference:'
      print'(3es21.14)', doverlapdqeffj - num_dSdqj
   end if

   if (any(abs(ddpintdqeffi - num_dDdqi) > thr_)) then
      write(*,*) "here 3"
      call test_failed(error, "Charge gradient of dipole integrals at center i does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', ddpintdqeffi
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dDdqi
      write(*,*) 'Difference:'
      print'(3es21.14)', ddpintdqeffi - num_dDdqi
   end if

   if (any(abs(ddpintdqeffj - num_dDdqj) > thr_)) then
      write(*,*) "here 4"
      call test_failed(error, "Charge gradient of dipole integrals at center j does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', ddpintdqeffj
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dDdqj
      write(*,*) 'Difference:'
      print'(3es21.14)', ddpintdqeffj - num_dDdqj
   end if

   if (any(abs(dqpintdqeffi - num_dQdqi) > thr_)) then
      write(*,*) "here 5"
      call test_failed(error, "Charge gradient of quadrupole integrals at center i does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', dqpintdqeffi
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dQdqi
      write(*,*) 'Difference:'
      print'(3es21.14)', dqpintdqeffi - num_dQdqi
   end if

   if (any(abs(dqpintdqeffj - num_dQdqj) > thr_)) then
      write(*,*) "here 6"
      call test_failed(error, "Charge gradient of quadrupole integrals at center j does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', dqpintdqeffj
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dQdqj
      write(*,*) 'Difference:'
      print'(3es21.14)', dqpintdqeffj - num_dQdqj
   end if

end subroutine test_multipole_grad_qeff_gen


!> Shift multipole operator from Ket function (center i) to Bra function (center j),
!> the multipole operator on the Bra function can be assembled from the lower moments
!> on the Ket function and the displacement vector using horizontal shift rules.
pure subroutine shift_operator(vec, s, di, qi, dj, qj)
   !> Displacement vector of center i and j
   real(wp),intent(in) :: vec(:)
   !> Overlap integral between basis functions
   real(wp),intent(in) :: s
   !> Dipole integral with operator on Ket function (center i)
   real(wp),intent(in) :: di(:)
   !> Quadrupole integral with operator on Ket function (center i)
   real(wp),intent(in) :: qi(:)
   !> Dipole integral with operator on Bra function (center j)
   real(wp),intent(out) :: dj(:)
   !> Quadrupole integral with operator on Bra function (center j)
   real(wp),intent(out) :: qj(:)

   real(wp) :: tr

   ! Create dipole operator on Bra function from Ket function and shift contribution
   ! due to monopol displacement
   dj(1) = di(1) + vec(1)*s
   dj(2) = di(2) + vec(2)*s
   dj(3) = di(3) + vec(3)*s

   ! For the quadrupole operator on the Bra function we first construct the shift
   ! contribution from the dipole and monopol displacement, since we have to remove
   ! the trace contribution from the shift and the moment integral on the Ket function
   ! is already traceless
   qj(1) = 2*vec(1)*di(1) + vec(1)**2*s
   qj(3) = 2*vec(2)*di(2) + vec(2)**2*s
   qj(6) = 2*vec(3)*di(3) + vec(3)**2*s
   qj(2) = vec(1)*di(2) + vec(2)*di(1) + vec(1)*vec(2)*s
   qj(4) = vec(1)*di(3) + vec(3)*di(1) + vec(1)*vec(3)*s
   qj(5) = vec(2)*di(3) + vec(3)*di(2) + vec(2)*vec(3)*s
   ! Now collect the trace of the shift contribution
   tr = 0.5_wp * (qj(1) + qj(3) + qj(6))

   ! Finally, assemble the quadrupole operator on the Bra function from the operator
   ! on the Ket function and the traceless shift contribution
   qj(1) = qi(1) + 1.5_wp * qj(1) - tr
   qj(2) = qi(2) + 1.5_wp * qj(2)
   qj(3) = qi(3) + 1.5_wp * qj(3) - tr
   qj(4) = qi(4) + 1.5_wp * qj(4)
   qj(5) = qi(5) + 1.5_wp * qj(5)
   qj(6) = qi(6) + 1.5_wp * qj(6) - tr
end subroutine shift_operator


subroutine test_multipole_grad_qeff_ss(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 1, 1, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 2, 1, .true., error)
   if (allocated(error)) return

   call test_multipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1*10)
end subroutine test_multipole_grad_qeff_ss

subroutine test_multipole_grad_qeff_sp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 3, 1, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 4, 2, .true., error)
   if (allocated(error)) return

   call test_multipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1*10)
end subroutine test_multipole_grad_qeff_sp

subroutine test_multipole_grad_qeff_pp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 5, 2, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 6, 2, .true., error)
   if (allocated(error)) return

   call test_multipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1*100)
end subroutine test_multipole_grad_qeff_pp

subroutine test_multipole_grad_qeff_sd(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 7, 1, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 8, 3, .true., error)
   if (allocated(error)) return

   call test_multipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1*10)
end subroutine test_multipole_grad_qeff_sd

subroutine test_multipole_grad_qeff_pd(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 9, 2, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 10, 3, .true., error)
   if (allocated(error)) return

   call test_multipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1*10)
end subroutine test_multipole_grad_qeff_pd

subroutine test_multipole_grad_qeff_dd(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 11, 3, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 12, 3, .true., error)
   if (allocated(error)) return

   call test_multipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1*100)
end subroutine test_multipole_grad_qeff_dd

subroutine test_multipole_grad_qeff_sf(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 13, 1, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 59, 4, .true., error)
   if (allocated(error)) return

   call test_multipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1*10)
end subroutine test_multipole_grad_qeff_sf

subroutine test_multipole_grad_qeff_pf(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 14, 2, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 60, 4, .true., error)
   if (allocated(error)) return

   call test_multipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1*10)
end subroutine test_multipole_grad_qeff_pf

subroutine test_multipole_grad_qeff_df(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 15, 3, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 61, 4, .true., error)
   if (allocated(error)) return

   call test_multipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1*10)
end subroutine test_multipole_grad_qeff_df

subroutine test_multipole_grad_qeff_ff(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   
   type(qvszp_cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)

   ! Randomly oriented vector of length 1.5 Bohr
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 1.5_wp * vec / sqrt(sum(vec**2))

   call new_qvszp_cgto(cgtoi, 62, 4, .true., error)
   if (allocated(error)) return
   call new_qvszp_cgto(cgtoj, 63, 4, .true., error)
   if (allocated(error)) return

   call test_multipole_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1*10)
end subroutine test_multipole_grad_qeff_ff


subroutine test_overlap_dipole_mol(error, mol, make_basis, &
   & qdep, diat_scale, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Factory to create new basis objects
   procedure(basis_maker) :: make_basis
   !> Flag whether the basis is dependent on the charge
   logical, intent(in) :: qdep
   !> Flag whether the diatomic scaling is applied
   logical, intent(in) :: diat_scale
   !> Reference value to check against
   real(wp), intent(in) :: ref(:, :)

   class(basis_type), allocatable :: bas
   type(basis_cache) :: bcache
   type(wavefunction_type), allocatable :: wfn_aux
   real(wp), allocatable :: lattr(:, :), overlap(:, :), overlap_diat(:, :)
   real(wp), allocatable :: dipole(:, :, :)
   real(wp) :: cutoff
   integer :: ii, jj
   real(wp), allocatable :: ksig(:,:), kpi(:,:), kdel(:,:) 

   allocate(ksig(mol%nid, mol%nid), kpi(mol%nid, mol%nid), kdel(mol%nid, mol%nid))
   ksig = 1.2_wp
   kpi = 1.2_wp
   kdel = 1.2_wp

   ! Setup basis set
   call make_basis(bas, mol, error)
   call check(error, bas%nao, size(ref, 1))
   if (allocated(error)) return

   if(qdep) then
      ! Obtain EEQBC charges for charge adaptation
      allocate(wfn_aux)
      call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
      call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
      if (allocated(error)) return
   end if

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   ! Update basis cache
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   allocate(overlap(bas%nao, bas%nao), overlap_diat(bas%nao, bas%nao))
   allocate(dipole(3, bas%nao, bas%nao))
   if(diat_scale) then 
      call get_dipole_integrals(mol, lattr, cutoff, bas, bcache, &
         & ksig, kpi, kdel, overlap, overlap_diat, dipole)

      do ii = 1, size(overlap_diat, 2)
         do jj = 1, size(overlap_diat, 1)
            call check(error, overlap_diat(jj, ii), ref(jj, ii), thr=thr)
            if (allocated(error)) return
         end do
      end do
   else
      call get_dipole_integrals(mol, lattr, cutoff, bas, bcache, &
         & overlap, dipole)

      do ii = 1, size(overlap, 2)
         do jj = 1, size(overlap, 1)
            call check(error, overlap(jj, ii), ref(jj, ii), thr=thr)
            if (allocated(error)) return
         end do
      end do
   end if

end subroutine test_overlap_dipole_mol

subroutine test_overlap_dipole_diat_alh3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 12
   real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      & 9.99999999869333E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.83197454971925E-1_wp, 4.83197454971925E-1_wp, 4.83197454971926E-1_wp,&
      & 0.00000000000000E+0_wp, 9.99999999998060E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 5.19330455705926E-1_wp,-5.19330455705926E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999998060E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999998060E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.99835578400188E-1_wp,-2.99835578400188E-1_wp, 5.99671156800377E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.08220347867244E-1_wp, 4.08220347867244E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999830206E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp,-2.72146898578163E-1_wp,-2.72146898578163E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp,&
      &-2.35686127729836E-1_wp,-2.35686127729836E-1_wp, 4.71372255459672E-1_wp,&
      & 4.83197454971925E-1_wp, 5.19330455705926E-1_wp, 0.00000000000000E+0_wp,&
      &-2.99835578400188E-1_wp,-4.08220347867244E-1_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp, 0.00000000000000E+0_wp,-2.35686127729836E-1_wp,&
      & 9.99999999881495E-1_wp, 4.25520423996964E-2_wp, 4.25520423996966E-2_wp,&
      & 4.83197454971925E-1_wp,-5.19330455705926E-1_wp, 0.00000000000000E+0_wp,&
      &-2.99835578400188E-1_wp, 4.08220347867244E-1_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp, 0.00000000000000E+0_wp,-2.35686127729836E-1_wp,&
      & 4.25520423996964E-2_wp, 9.99999999881495E-1_wp, 4.25520423996966E-2_wp,&
      & 4.83197454971927E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 5.99671156800377E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp, 0.00000000000000E+0_wp, 4.71372255459672E-1_wp,&
      & 4.25520423996966E-2_wp, 4.25520423996966E-2_wp, 9.99999999881495E-1_wp],&
      & shape(overlap_diat))
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "AlH3")
   call test_overlap_dipole_mol(error, mol, make_gen_basis, &
      & .false., .true., overlap_diat)

end subroutine test_overlap_dipole_diat_alh3

subroutine test_overlap_dipole_qvszp_alh3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 12
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 4.22357985345720E-01_wp, 4.22357985345720E-01_wp, 4.22357985345721E-01_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 4.56235837162591E-01_wp,-4.56235837162591E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-2.63407883399776E-01_wp,-2.63407883399776E-01_wp, 5.26815766799553E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-3.08781350958959E-01_wp, 3.08781350958959E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-2.05854233972639E-01_wp,-2.05854233972639E-01_wp,-2.05854233972640E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      &-1.78274996096891E-01_wp,-1.78274996096891E-01_wp, 3.56549992193783E-01_wp, &
      & 4.22357985345720E-01_wp, 4.56235837162591E-01_wp, 0.00000000000000E+00_wp, &
      &-2.63407883399776E-01_wp,-3.08781350958959E-01_wp, 0.00000000000000E+00_wp, &
      &-2.05854233972639E-01_wp, 0.00000000000000E+00_wp,-1.78274996096891E-01_wp, &
      & 1.00000000000000E+00_wp, 6.29386757686522E-02_wp, 6.29386757686522E-02_wp, &
      & 4.22357985345720E-01_wp,-4.56235837162591E-01_wp, 0.00000000000000E+00_wp, &
      &-2.63407883399776E-01_wp, 3.08781350958959E-01_wp, 0.00000000000000E+00_wp, &
      &-2.05854233972639E-01_wp, 0.00000000000000E+00_wp,-1.78274996096891E-01_wp, &
      & 6.29386757686522E-02_wp, 1.00000000000000E+00_wp, 6.29386757686522E-02_wp, &
      & 4.22357985345721E-01_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 5.26815766799553E-01_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-2.05854233972640E-01_wp, 0.00000000000000E+00_wp, 3.56549992193783E-01_wp, &
      & 6.29386757686522E-02_wp, 6.29386757686522E-02_wp, 1.00000000000000E+00_wp],&
      & shape(overlap))
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "AlH3")
   call test_overlap_dipole_mol(error, mol, make_qvszp_basis, &
      & .true., .false., overlap)

end subroutine test_overlap_dipole_qvszp_alh3

subroutine test_overlap_multipole_mol(error, mol, make_basis, &
   & qdep, diat_scale, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Factory to create new basis objects
   procedure(basis_maker) :: make_basis
   !> Flag whether the basis is dependent on the charge
   logical, intent(in) :: qdep
   !> Flag whether the diatomic scaling is applied
   logical, intent(in) :: diat_scale
   !> Reference value to check against
   real(wp), intent(in) :: ref(:, :)

   class(basis_type), allocatable :: bas
   type(basis_cache) :: bcache
   type(wavefunction_type), allocatable :: wfn_aux
   real(wp), allocatable :: lattr(:, :), overlap(:, :), overlap_diat(:, :)
   real(wp), allocatable :: dipole(:, :, :), quadrupole(:, :, :)
   real(wp) :: cutoff
   integer :: ii, jj
   real(wp), allocatable :: ksig(:,:), kpi(:,:), kdel(:,:) 

   allocate(ksig(mol%nid, mol%nid), kpi(mol%nid, mol%nid), kdel(mol%nid, mol%nid))
   ksig = 1.2_wp
   kpi = 1.2_wp
   kdel = 1.2_wp

   ! Setup basis set
   call make_basis(bas, mol, error)
   call check(error, bas%nao, size(ref, 1))
   if (allocated(error)) return

   if(qdep) then
      ! Obtain EEQBC charges for charge adaptation
      allocate(wfn_aux)
      call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
      call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
      if (allocated(error)) return
   end if

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   ! Update basis cache
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   allocate(overlap(bas%nao, bas%nao), overlap_diat(bas%nao, bas%nao))
   allocate(dipole(3, bas%nao, bas%nao), quadrupole(6, bas%nao, bas%nao))
   if(diat_scale) then 
      call get_multipole_integrals(mol, lattr, cutoff, bas, bcache, &
         & ksig, kpi, kdel, overlap, overlap_diat, dipole, quadrupole)

      do ii = 1, size(overlap_diat, 2)
         do jj = 1, size(overlap_diat, 1)
            call check(error, overlap_diat(jj, ii), ref(jj, ii), thr=thr)
            if (allocated(error)) return
         end do
      end do
   else
      call get_multipole_integrals(mol, lattr, cutoff, bas, bcache, &
         & overlap, dipole, quadrupole)

      do ii = 1, size(overlap, 2)
         do jj = 1, size(overlap, 1)
            call check(error, overlap(jj, ii), ref(jj, ii), thr=thr)
            if (allocated(error)) return
         end do
      end do
   end if

end subroutine test_overlap_multipole_mol

subroutine test_overlap_multipole_diat_alh3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 12
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 9.99999999869333E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.83197454971925E-1_wp, 4.83197454971925E-1_wp, 4.83197454971926E-1_wp,&
      & 0.00000000000000E+0_wp, 9.99999999998060E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 5.19330455705926E-1_wp,-5.19330455705926E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999998060E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999998060E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.99835578400188E-1_wp,-2.99835578400188E-1_wp, 5.99671156800377E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.08220347867244E-1_wp, 4.08220347867244E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999830206E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp,-2.72146898578163E-1_wp,-2.72146898578163E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp,&
      &-2.35686127729836E-1_wp,-2.35686127729836E-1_wp, 4.71372255459672E-1_wp,&
      & 4.83197454971925E-1_wp, 5.19330455705926E-1_wp, 0.00000000000000E+0_wp,&
      &-2.99835578400188E-1_wp,-4.08220347867244E-1_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp, 0.00000000000000E+0_wp,-2.35686127729836E-1_wp,&
      & 9.99999999881495E-1_wp, 4.25520423996964E-2_wp, 4.25520423996966E-2_wp,&
      & 4.83197454971925E-1_wp,-5.19330455705926E-1_wp, 0.00000000000000E+0_wp,&
      &-2.99835578400188E-1_wp, 4.08220347867244E-1_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp, 0.00000000000000E+0_wp,-2.35686127729836E-1_wp,&
      & 4.25520423996964E-2_wp, 9.99999999881495E-1_wp, 4.25520423996966E-2_wp,&
      & 4.83197454971927E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 5.99671156800377E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.72146898578163E-1_wp, 0.00000000000000E+0_wp, 4.71372255459672E-1_wp,&
      & 4.25520423996966E-2_wp, 4.25520423996966E-2_wp, 9.99999999881495E-1_wp],&
      & shape(overlap))
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "AlH3")
   call test_overlap_multipole_mol(error, mol, make_gen_basis, &
      & .false., .true., overlap)

end subroutine test_overlap_multipole_diat_alh3


subroutine test_overlap_multipole_qvszp_alh3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 12
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 4.22357985345720E-01_wp, 4.22357985345720E-01_wp, 4.22357985345721E-01_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 4.56235837162591E-01_wp,-4.56235837162591E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-2.63407883399776E-01_wp,-2.63407883399776E-01_wp, 5.26815766799553E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-3.08781350958959E-01_wp, 3.08781350958959E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-2.05854233972639E-01_wp,-2.05854233972639E-01_wp,-2.05854233972640E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      &-1.78274996096891E-01_wp,-1.78274996096891E-01_wp, 3.56549992193783E-01_wp, &
      & 4.22357985345720E-01_wp, 4.56235837162591E-01_wp, 0.00000000000000E+00_wp, &
      &-2.63407883399776E-01_wp,-3.08781350958959E-01_wp, 0.00000000000000E+00_wp, &
      &-2.05854233972639E-01_wp, 0.00000000000000E+00_wp,-1.78274996096891E-01_wp, &
      & 1.00000000000000E+00_wp, 6.29386757686522E-02_wp, 6.29386757686522E-02_wp, &
      & 4.22357985345720E-01_wp,-4.56235837162591E-01_wp, 0.00000000000000E+00_wp, &
      &-2.63407883399776E-01_wp, 3.08781350958959E-01_wp, 0.00000000000000E+00_wp, &
      &-2.05854233972639E-01_wp, 0.00000000000000E+00_wp,-1.78274996096891E-01_wp, &
      & 6.29386757686522E-02_wp, 1.00000000000000E+00_wp, 6.29386757686522E-02_wp, &
      & 4.22357985345721E-01_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 5.26815766799553E-01_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-2.05854233972640E-01_wp, 0.00000000000000E+00_wp, 3.56549992193783E-01_wp, &
      & 6.29386757686522E-02_wp, 6.29386757686522E-02_wp, 1.00000000000000E+00_wp],&
      & shape(overlap))
   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "AlH3")
   call test_overlap_multipole_mol(error, mol, make_qvszp_basis, &
      & .true., .false., overlap)

end subroutine test_overlap_multipole_qvszp_alh3

end module test_integral_multipole
