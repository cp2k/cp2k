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

module test_integral_overlap
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use multicharge, only : get_eeqbc_charges
   use mstore, only : get_structure
   use tblite_basis_cache, only : basis_cache, cgto_cache
   use tblite_basis_qvszp, only : qvszp_basis_type, qvszp_cgto_type, &
      & new_qvszp_cgto, new_qvszp_basis
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_basis_type, only : basis_type, new_basis, cgto_type, cgto_container, &
      & get_cutoff
   use tblite_cutoff, only : get_lattice_points
   use tblite_integral_overlap
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction

   implicit none
   private

   public :: collect_integral_overlap

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
subroutine collect_integral_overlap(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("overlap-alh3", test_overlap_alh3), &
      new_unittest("overlap-diat-alh3", test_overlap_diat_alh3), &
      new_unittest("overlap-qvszp-alh3", test_overlap_qvszp_alh3), &
      new_unittest("overlap-qvszp-diat-alh3", test_overlap_qvszp_diat_alh3), &
      new_unittest("overlap-bh3", test_overlap_bh3), &
      new_unittest("overlap-beh2", test_overlap_beh2), &
      new_unittest("overlap-ch4", test_overlap_ch4), &
      new_unittest("overlap-cl2", test_overlap_cl2), &
      new_unittest("overlap-diat-cl2", test_overlap_diat_cl2), &
      new_unittest("overlap-f2", test_overlap_f2), &
      new_unittest("overlap-h2", test_overlap_h2), &
      new_unittest("overlap-lih", test_overlap_lih), &
      new_unittest("overlap-qvszp-cecl3", test_overlap_qvszp_cecl3), &
      new_unittest("overlap-diat-cecl3", test_overlap_diat_cecl3), &
      new_unittest("overlap-grad-ss", test_overlap_grad_ss), &
      new_unittest("overlap-grad-pp", test_overlap_grad_pp), &
      new_unittest("overlap-grad-dd", test_overlap_grad_dd), &
      new_unittest("overlap-grad-ff", test_overlap_grad_ff), &
      new_unittest("overlap-grad-qeff-ss", test_overlap_grad_qeff_ss), &
      new_unittest("overlap-grad-qeff-sp", test_overlap_grad_qeff_sp), &
      new_unittest("overlap-grad-qeff-pp", test_overlap_grad_qeff_pp), &
      new_unittest("overlap-grad-qeff-sd", test_overlap_grad_qeff_sd), &
      new_unittest("overlap-grad-qeff-pd", test_overlap_grad_qeff_pd), &
      new_unittest("overlap-grad-qeff-dd", test_overlap_grad_qeff_dd), &
      new_unittest("overlap-grad-qeff-sf", test_overlap_grad_qeff_sf), &
      new_unittest("overlap-grad-qeff-pf", test_overlap_grad_qeff_pf), &
      new_unittest("overlap-grad-qeff-df", test_overlap_grad_qeff_df), &
      new_unittest("overlap-grad-qeff-ff", test_overlap_grad_qeff_ff), &
      new_unittest("overlap-diat-grad-ss", test_overlap_diat_grad_ss), &
      new_unittest("overlap-diat-grad-ss_z", test_overlap_diat_grad_ss_z), &
      new_unittest("overlap-diat-grad-sp", test_overlap_diat_grad_sp), &
      new_unittest("overlap-diat-grad-sp_z", test_overlap_diat_grad_sp_z), &
      new_unittest("overlap-diat-grad-pp", test_overlap_diat_grad_pp), &
      new_unittest("overlap-diat-grad-pp_z", test_overlap_diat_grad_pp_z), &
      new_unittest("overlap-diat-grad-sd", test_overlap_diat_grad_sd), &
      new_unittest("overlap-diat-grad-sd_z", test_overlap_diat_grad_sd_z), &
      new_unittest("overlap-diat-grad-pd", test_overlap_diat_grad_pd), &
      new_unittest("overlap-diat-grad-pd_z", test_overlap_diat_grad_pd_z), &
      new_unittest("overlap-diat-grad-dd", test_overlap_diat_grad_dd), &
      new_unittest("overlap-diat-grad-dd_z", test_overlap_diat_grad_dd_z), &
      new_unittest("overlap-diat-grad-sf", test_overlap_diat_grad_sf), &
      new_unittest("overlap-diat-grad-sf_z", test_overlap_diat_grad_sf_z), &
      new_unittest("overlap-diat-grad-pf", test_overlap_diat_grad_pf), &
      new_unittest("overlap-diat-grad-pf_z", test_overlap_diat_grad_pf_z), &
      new_unittest("overlap-diat-grad-df", test_overlap_diat_grad_df), &
      new_unittest("overlap-diat-grad-df_z", test_overlap_diat_grad_df_z), &
      new_unittest("overlap-diat-grad-ff", test_overlap_diat_grad_ff), &
      new_unittest("overlap-diat-grad-ff_z", test_overlap_diat_grad_ff_z) &
      ]

end subroutine collect_integral_overlap


subroutine make_gen_basis(bas, mol, error)
   class(basis_type), allocatable, intent(out) :: bas
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   integer, parameter :: nsh(60) = [&
      & 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, & ! 1-20
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, & ! 21-40
      & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 4, 3, 3] ! 41-60

   integer, parameter :: lsh(4, 60) = reshape([&
      & 0, 0, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0, & ! 1-6
      & 0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 2, 0,  0, 1, 0, 0,  0, 1, 2, 0, & ! 7-12
      & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 13-18
      & 0, 1, 0, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 19-24
      & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 0, 0, & ! 25-30
      & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 31-36
      & 0, 1, 0, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 37-42
      & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 0, 0, & ! 43-48
      & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 49-54
      & 0, 1, 0, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 3,  0, 1, 2, 0,  0, 1, 2, 0],& ! 55-60
      & shape(lsh))

   integer, parameter :: pqn(4, 60) = reshape([&
      & 1, 0, 0, 0,  1, 2, 0, 0,  2, 2, 0, 0,  2, 2, 0, 0,  2, 2, 0, 0,  2, 2, 0, 0, & ! 1-6
      & 2, 2, 0, 0,  2, 2, 0, 0,  2, 2, 0, 0,  2, 2, 3, 0,  3, 3, 0, 0,  3, 3, 3, 0, & ! 7-12
      & 3, 3, 3, 0,  3, 3, 3, 0,  3, 3, 3, 0,  3, 3, 3, 0,  3, 3, 3, 0,  3, 3, 3, 0, & ! 13-18
      & 4, 4, 0, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0, & ! 19-24
      & 4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 0, 0, & ! 25-30
      & 4, 4, 4, 0,  4, 4, 4, 0,  4, 4, 4, 0,  4, 4, 4, 0,  4, 4, 4, 0,  4, 4, 4, 0, & ! 31-36
      & 5, 5, 0, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0, & ! 37-42
      & 5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 0, 0, & ! 43-48
      & 5, 5, 5, 0,  5, 5, 5, 0,  5, 5, 5, 0,  5, 5, 5, 0,  5, 5, 5, 0,  5, 5, 5, 0, & ! 49-54
      & 6, 6, 0, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 4,  6, 6, 5, 0,  6, 6, 5, 0],& ! 55-60
      & shape(pqn))

   real(wp), parameter :: zeta(4, 60) = reshape([&
      & 1.230000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, & !1
      & 1.669667_wp, 1.500000_wp, 0.000000_wp, 0.000000_wp, & !2
      & 0.750060_wp, 0.557848_wp, 0.000000_wp, 0.000000_wp, & !3
      & 1.034720_wp, 0.949332_wp, 0.000000_wp, 0.000000_wp, & !4
      & 1.479444_wp, 1.479805_wp, 0.000000_wp, 0.000000_wp, & !5
      & 2.096432_wp, 1.800000_wp, 0.000000_wp, 0.000000_wp, & !6
      & 2.339881_wp, 2.014332_wp, 0.000000_wp, 0.000000_wp, & !7
      & 2.439742_wp, 2.137023_wp, 0.000000_wp, 0.000000_wp, & !8
      & 2.416361_wp, 2.308399_wp, 0.000000_wp, 0.000000_wp, & !9
      & 3.084104_wp, 2.312051_wp, 2.815609_wp, 0.000000_wp, & !10
      & 0.763787_wp, 0.573553_wp, 0.000000_wp, 0.000000_wp, & !11
      & 1.184203_wp, 0.717769_wp, 1.300000_wp, 0.000000_wp, & !12
      & 1.352531_wp, 1.391201_wp, 1.000000_wp, 0.000000_wp, & !13
      & 1.773917_wp, 1.718996_wp, 1.250000_wp, 0.000000_wp, & !14
      & 1.816945_wp, 1.903247_wp, 1.167533_wp, 0.000000_wp, & !15
      & 1.981333_wp, 2.025643_wp, 1.702555_wp, 0.000000_wp, & !16
      & 2.485265_wp, 2.199650_wp, 2.476089_wp, 0.000000_wp, & !17
      & 2.329679_wp, 2.149419_wp, 1.950531_wp, 0.000000_wp, & !18
      & 0.875961_wp, 0.631694_wp, 0.000000_wp, 0.000000_wp, & !19
      & 1.267130_wp, 0.786247_wp, 1.380000_wp, 0.000000_wp, & !20
      & 2.224492_wp, 1.554183_wp, 2.009535_wp, 0.000000_wp, & !21
      & 2.588796_wp, 0.994410_wp, 1.885617_wp, 0.000000_wp, & !22
      & 3.043706_wp, 4.030076_wp, 1.663291_wp, 0.000000_wp, & !23
      & 2.250127_wp, 2.706815_wp, 1.675019_wp, 0.000000_wp, & !24
      & 2.206053_wp, 2.820197_wp, 1.861022_wp, 0.000000_wp, & !25
      & 1.572970_wp, 1.986214_wp, 2.837906_wp, 0.000000_wp, & !26
      & 1.808266_wp, 1.736758_wp, 2.797674_wp, 0.000000_wp, & !27
      & 2.007589_wp, 2.250756_wp, 2.982916_wp, 0.000000_wp, & !28
      & 2.181599_wp, 2.384590_wp, 3.095025_wp, 0.000000_wp, & !29
      & 2.263767_wp, 2.203629_wp, 0.000000_wp, 0.000000_wp, & !30
      & 2.638221_wp, 2.067523_wp, 2.113616_wp, 0.000000_wp, & !31
      & 2.528919_wp, 2.194417_wp, 1.776619_wp, 0.000000_wp, & !32
      & 3.556676_wp, 2.420754_wp, 1.465797_wp, 0.000000_wp, & !33
      & 2.896526_wp, 2.454218_wp, 2.278836_wp, 0.000000_wp, & !34
      & 3.289210_wp, 2.565269_wp, 1.645016_wp, 0.000000_wp, & !35
      & 5.209881_wp, 2.843367_wp, 2.758388_wp, 0.000000_wp, & !36
      & 1.269729_wp, 1.887305_wp, 0.000000_wp, 0.000000_wp, & !37
      & 1.868807_wp, 1.785463_wp, 2.160122_wp, 0.000000_wp, & !38
      & 0.920018_wp, 1.457324_wp, 2.229013_wp, 0.000000_wp, & !39
      & 6.506473_wp, 1.432023_wp, 2.119714_wp, 0.000000_wp, & !40
      & 2.109733_wp, 2.799447_wp, 2.018973_wp, 0.000000_wp, & !41
      & 2.584133_wp, 3.027953_wp, 2.087336_wp, 0.000000_wp, & !42
      & 2.621415_wp, 3.134876_wp, 2.132598_wp, 0.000000_wp, & !43
      & 2.739844_wp, 2.181678_wp, 2.546096_wp, 0.000000_wp, & !44
      & 1.840571_wp, 2.974826_wp, 3.106937_wp, 0.000000_wp, & !45
      & 1.756228_wp, 3.394247_wp, 3.202653_wp, 0.000000_wp, & !46
      & 3.050188_wp, 2.349519_wp, 3.353329_wp, 0.000000_wp, & !47
      & 2.419991_wp, 2.288929_wp, 0.000000_wp, 0.000000_wp, & !48
      & 2.878139_wp, 2.446597_wp, 2.757735_wp, 0.000000_wp, & !49
      & 3.038232_wp, 2.320821_wp, 1.775133_wp, 0.000000_wp, & !50
      & 2.687507_wp, 2.385653_wp, 2.125961_wp, 0.000000_wp, & !51
      & 2.810717_wp, 2.452747_wp, 2.018718_wp, 0.000000_wp, & !52
      & 2.906869_wp, 2.493771_wp, 1.900737_wp, 0.000000_wp, & !53
      & 4.175313_wp, 2.869379_wp, 2.968948_wp, 0.000000_wp, & !54
      & 1.242993_wp, 1.991420_wp, 0.000000_wp, 0.000000_wp, & !55
      & 1.314003_wp, 1.164384_wp, 2.127596_wp, 0.000000_wp, & !56
      & 2.817373_wp, 1.698633_wp, 2.273697_wp, 0.000000_wp, & !57
      & 2.845039_wp, 1.460181_wp, 2.534989_wp, 2.534989_wp, & !58
      & 2.816971_wp, 1.475453_wp, 2.543502_wp, 0.000000_wp, & !59
      & 2.788903_wp, 1.490724_wp, 2.552016_wp, 0.000000_wp],& !60
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


subroutine test_overlap_mol(error, mol, make_basis, qdep, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Factory to create new basis objects
   procedure(basis_maker) :: make_basis
   !> Flag whether the basis is dependent on the charge
   logical, intent(in) :: qdep
   !> Reference value to check against
   real(wp), intent(in) :: ref(:, :)

   class(basis_type), allocatable :: bas
   type(basis_cache) :: bcache
   type(wavefunction_type), allocatable :: wfn_aux
   real(wp), allocatable :: lattr(:, :), overlap(:, :)
   real(wp) :: cutoff
   integer :: ii, jj

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

   allocate(overlap(bas%nao, bas%nao))
   call get_overlap(mol, lattr, cutoff, bas, bcache, overlap)

   ! where(abs(overlap) < thr) overlap = 0.0_wp
   ! print '(*(6x,"&", 3(es20.14e1, "_wp":, ","), "&", /))', overlap

   do ii = 1, size(overlap, 2)
      do jj = 1, size(overlap, 1)
         call check(error, overlap(jj, ii), ref(jj, ii), thr=thr1)
         if (allocated(error)) return
      end do
   end do


end subroutine test_overlap_mol

subroutine test_overlap_diat_mol(error, mol, make_basis, qdep, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Factory to create new basis objects
   procedure(basis_maker) :: make_basis
   !> Flag whether the basis is dependent on the charge
   logical, intent(in) :: qdep
   !> Reference value to check against
   real(wp), intent(in) :: ref(:, :)

   class(basis_type), allocatable :: bas
   type(basis_cache) :: bcache
   type(wavefunction_type), allocatable :: wfn_aux
   real(wp), allocatable :: lattr(:, :), overlap(:, :), overlap_diat(:, :)
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
   call get_overlap(mol, lattr, cutoff, bas, bcache, ksig, kpi, kdel, &
      & overlap, overlap_diat)

   do ii = 1, size(overlap_diat, 2)
      do jj = 1, size(overlap_diat, 1)
         call check(error, overlap_diat(jj, ii), ref(jj, ii), thr=thr1*10)
         if (allocated(error)) return
      end do
   end do

end subroutine test_overlap_diat_mol

subroutine test_overlap_alh3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 12
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 9.99999999869333E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.02664545809937E-1_wp, 4.02664545809937E-1_wp, 4.02664545809939E-1_wp,&
      & 0.00000000000000E+0_wp, 9.99999999998060E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.32775379754938E-1_wp,-4.32775379754938E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999998060E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999998060E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.49862982000156E-1_wp,-2.49862982000156E-1_wp, 4.99725964000314E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-3.40183623222704E-1_wp, 3.40183623222704E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999830206E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.26789082148469E-1_wp,-2.26789082148469E-1_wp,-2.26789082148469E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830206E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp,&
      &-1.96405106441530E-1_wp,-1.96405106441530E-1_wp, 3.92810212883060E-1_wp,&
      & 4.02664545809937E-1_wp, 4.32775379754938E-1_wp, 0.00000000000000E+0_wp,&
      &-2.49862982000156E-1_wp,-3.40183623222704E-1_wp, 0.00000000000000E+0_wp,&
      &-2.26789082148469E-1_wp, 0.00000000000000E+0_wp,-1.96405106441530E-1_wp,&
      & 9.99999999881495E-1_wp, 3.54600353330803E-2_wp, 3.54600353330805E-2_wp,&
      & 4.02664545809937E-1_wp,-4.32775379754938E-1_wp, 0.00000000000000E+0_wp,&
      &-2.49862982000156E-1_wp, 3.40183623222704E-1_wp, 0.00000000000000E+0_wp,&
      &-2.26789082148469E-1_wp, 0.00000000000000E+0_wp,-1.96405106441530E-1_wp,&
      & 3.54600353330803E-2_wp, 9.99999999881495E-1_wp, 3.54600353330805E-2_wp,&
      & 4.02664545809939E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.99725964000314E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.26789082148469E-1_wp, 0.00000000000000E+0_wp, 3.92810212883060E-1_wp,&
      & 3.54600353330805E-2_wp, 3.54600353330805E-2_wp, 9.99999999881495E-1_wp],&
      & shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "AlH3")
   call test_overlap_mol(error, mol, make_gen_basis, .false., overlap)

end subroutine test_overlap_alh3

subroutine test_overlap_diat_alh3(error)

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
   call test_overlap_diat_mol(error, mol, make_gen_basis, .false., overlap_diat)

end subroutine test_overlap_diat_alh3

subroutine test_overlap_qvszp_alh3(error)

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
   call test_overlap_mol(error, mol, make_qvszp_basis, .true., overlap)

end subroutine test_overlap_qvszp_alh3

subroutine test_overlap_qvszp_diat_alh3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 12
   real(wp), parameter :: overlap_diat(nao, nao) = reshape([&
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 5.06829582414865E-01_wp, 5.06829582414865E-01_wp, 5.06829582414866E-01_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 5.47483004595109E-01_wp,-5.47483004595109E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-3.16089460079732E-01_wp,-3.16089460079732E-01_wp, 6.32178920159464E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-3.70537621150751E-01_wp, 3.70537621150751E-01_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-2.47025080767167E-01_wp,-2.47025080767167E-01_wp,-2.47025080767168E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      &-2.13929995316269E-01_wp,-2.13929995316269E-01_wp, 4.27859990632540E-01_wp, &
      & 5.06829582414865E-01_wp, 5.47483004595109E-01_wp, 0.00000000000000E+00_wp, &
      &-3.16089460079732E-01_wp,-3.70537621150751E-01_wp, 0.00000000000000E+00_wp, &
      &-2.47025080767167E-01_wp, 0.00000000000000E+00_wp,-2.13929995316269E-01_wp, &
      & 1.00000000000000E+00_wp, 7.55264109223826E-02_wp, 7.55264109223826E-02_wp, &
      & 5.06829582414865E-01_wp,-5.47483004595109E-01_wp, 0.00000000000000E+00_wp, &
      &-3.16089460079732E-01_wp, 3.70537621150751E-01_wp, 0.00000000000000E+00_wp, &
      &-2.47025080767167E-01_wp, 0.00000000000000E+00_wp,-2.13929995316269E-01_wp, &
      & 7.55264109223826E-02_wp, 1.00000000000000E+00_wp, 7.55264109223826E-02_wp, &
      & 5.06829582414866E-01_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 6.32178920159464E-01_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-2.47025080767168E-01_wp, 0.00000000000000E+00_wp, 4.27859990632540E-01_wp, &
      & 7.55264109223826E-02_wp, 7.55264109223826E-02_wp, 1.00000000000000E+00_wp],&
      & shape(overlap_diat))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "AlH3")
   call test_overlap_diat_mol(error, mol, make_qvszp_basis, .true., overlap_diat)

end subroutine test_overlap_qvszp_diat_alh3

subroutine test_overlap_bh3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 7
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000600E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 4.80899928496809E-1_wp, 4.80899928496809E-1_wp,&
      & 4.80899928496810E-1_wp, 0.00000000000000E+0_wp, 9.99999999925689E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 4.26420603078617E-1_wp,&
      &-4.26420603078617E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999925689E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp,-2.46194049975443E-1_wp,-2.46194049975443E-1_wp,&
      & 4.92388099950885E-1_wp, 4.80899928496809E-1_wp, 4.26420603078617E-1_wp,&
      & 0.00000000000000E+0_wp,-2.46194049975443E-1_wp, 9.99999999881495E-1_wp,&
      & 1.10583333710332E-1_wp, 1.10583333710331E-1_wp, 4.80899928496809E-1_wp,&
      &-4.26420603078617E-1_wp, 0.00000000000000E+0_wp,-2.46194049975443E-1_wp,&
      & 1.10583333710332E-1_wp, 9.99999999881495E-1_wp, 1.10583333710331E-1_wp,&
      & 4.80899928496810E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.92388099950885E-1_wp, 1.10583333710331E-1_wp, 1.10583333710331E-1_wp,&
      & 9.99999999881495E-1_wp],shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "BH3")
   call test_overlap_mol(error, mol, make_gen_basis, .false., overlap)

end subroutine test_overlap_bh3

subroutine test_overlap_beh2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 6
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000600E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 4.75332491800271E-1_wp, 4.75332491800271E-1_wp,&
      & 0.00000000000000E+0_wp, 9.99999999925689E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999925689E-1_wp,&
      & 0.00000000000000E+0_wp,-5.55623013601629E-1_wp, 5.55623013601629E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 4.75332491800271E-1_wp, 0.00000000000000E+0_wp,-5.55623013601629E-1_wp,&
      & 0.00000000000000E+0_wp, 9.99999999881495E-1_wp, 4.00645037406879E-2_wp,&
      & 4.75332491800271E-1_wp, 0.00000000000000E+0_wp, 5.55623013601629E-1_wp,&
      & 0.00000000000000E+0_wp, 4.00645037406879E-2_wp, 9.99999999881495E-1_wp],&
      & shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "BeH2")
   call test_overlap_mol(error, mol, make_gen_basis, .false., overlap)

end subroutine test_overlap_beh2

subroutine test_overlap_ch4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 8
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000600E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 4.22144465267144E-1_wp, 4.22144465267144E-1_wp,&
      & 4.22144465267144E-1_wp, 4.22144465267144E-1_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 2.57682345348632E-1_wp,-2.57682345348632E-1_wp,-2.57682345348632E-1_wp,&
      & 2.57682345348632E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 2.57682345348632E-1_wp,&
      & 2.57682345348632E-1_wp,-2.57682345348632E-1_wp,-2.57682345348632E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp,-2.57682345348632E-1_wp, 2.57682345348632E-1_wp,&
      &-2.57682345348632E-1_wp, 2.57682345348632E-1_wp, 4.22144465267144E-1_wp,&
      & 2.57682345348632E-1_wp, 2.57682345348632E-1_wp,-2.57682345348632E-1_wp,&
      & 9.99999999881495E-1_wp, 1.71953274424808E-1_wp, 1.71953274424808E-1_wp,&
      & 1.71953274424808E-1_wp, 4.22144465267144E-1_wp,-2.57682345348632E-1_wp,&
      & 2.57682345348632E-1_wp, 2.57682345348632E-1_wp, 1.71953274424808E-1_wp,&
      & 9.99999999881495E-1_wp, 1.71953274424808E-1_wp, 1.71953274424808E-1_wp,&
      & 4.22144465267144E-1_wp,-2.57682345348632E-1_wp,-2.57682345348632E-1_wp,&
      &-2.57682345348632E-1_wp, 1.71953274424808E-1_wp, 1.71953274424808E-1_wp,&
      & 9.99999999881495E-1_wp, 1.71953274424808E-1_wp, 4.22144465267144E-1_wp,&
      & 2.57682345348632E-1_wp,-2.57682345348632E-1_wp, 2.57682345348632E-1_wp,&
      & 1.71953274424808E-1_wp, 1.71953274424808E-1_wp, 1.71953274424808E-1_wp,&
      & 9.99999999881495E-1_wp],shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "CH4")
   call test_overlap_mol(error, mol, make_gen_basis, .false., overlap)

end subroutine test_overlap_ch4

subroutine test_overlap_cl2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 18
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 9.99999999869332E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.17210205959034E-2_wp, 0.00000000000000E+0_wp,-1.66880110378166E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.08719248478375E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999998060E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 6.66413853324784E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-7.28795473739356E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999998060E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.66880110378166E-1_wp, 0.00000000000000E+0_wp,-2.43857348208300E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.54878516243462E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999998060E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 6.66413853324784E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-7.28795473739356E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830205E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.65407586749691E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830205E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 7.28795473739356E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-8.83217466681701E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999830206E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.08719248478375E-1_wp, 0.00000000000000E+0_wp,-1.54878516243462E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.20587356016693E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830205E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 7.28795473739356E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-8.83217466681701E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.65407586749691E-2_wp,&
      & 9.17210205959034E-2_wp, 0.00000000000000E+0_wp, 1.66880110378166E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.08719248478375E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999869332E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 6.66413853324784E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 7.28795473739356E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999998060E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.66880110378166E-1_wp, 0.00000000000000E+0_wp,-2.43857348208300E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.54878516243462E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999998060E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 6.66413853324784E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 7.28795473739356E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999998060E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.65407586749691E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830205E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-7.28795473739356E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-8.83217466681701E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830205E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.08719248478375E-1_wp, 0.00000000000000E+0_wp, 1.54878516243462E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.20587356016693E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999830206E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-7.28795473739356E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-8.83217466681701E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830205E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.65407586749691E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp],&
      & shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "Cl2")
   call test_overlap_mol(error, mol, make_gen_basis, .false., overlap)

end subroutine test_overlap_cl2

subroutine test_overlap_diat_cl2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 18
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 9.99999999869332E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.10065224715084E-1_wp, 0.00000000000000E+0_wp,-2.00256132453799E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.30463098174050E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999998060E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 7.99696623989741E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-8.74554568487227E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999998060E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 2.00256132453799E-1_wp, 0.00000000000000E+0_wp,-2.92628817849960E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.85854219492155E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999998060E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 7.99696623989741E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-8.74554568487227E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830205E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.98489104099629E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830205E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 8.74554568487227E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.05986096001804E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999830206E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.30463098174050E-1_wp, 0.00000000000000E+0_wp,-1.85854219492155E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.44704827220032E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830205E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 8.74554568487227E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.05986096001804E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.98489104099629E-2_wp,&
      & 1.10065224715084E-1_wp, 0.00000000000000E+0_wp, 2.00256132453799E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.30463098174050E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999869332E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 7.99696623989741E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 8.74554568487227E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999998060E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.00256132453799E-1_wp, 0.00000000000000E+0_wp,-2.92628817849960E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.85854219492155E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999998060E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 7.99696623989741E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 8.74554568487227E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999998060E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.98489104099629E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830205E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-8.74554568487227E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.05986096001804E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830205E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.30463098174050E-1_wp, 0.00000000000000E+0_wp, 1.85854219492155E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.44704827220032E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999830206E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-8.74554568487227E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.05986096001804E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 9.99999999830205E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.98489104099629E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 9.99999999830206E-1_wp],&
      & shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "Cl2")
   call test_overlap_diat_mol(error, mol, make_gen_basis, .false., overlap)

end subroutine test_overlap_diat_cl2

subroutine test_overlap_f2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 8
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000600E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.47958427103274E-1_wp, 0.00000000000000E+0_wp,&
      &-2.00620208590218E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 7.94087596423658E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 2.00620208590218E-1_wp,&
      & 0.00000000000000E+0_wp,-2.36369347922990E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 7.94087596423658E-2_wp, 1.47958427103274E-1_wp,&
      & 0.00000000000000E+0_wp, 2.00620208590218E-1_wp, 0.00000000000000E+0_wp,&
      & 1.00000000000600E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 7.94087596423658E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.00620208590218E-1_wp, 0.00000000000000E+0_wp,-2.36369347922990E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 7.94087596423658E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp],shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "F2")
   call test_overlap_mol(error, mol, make_gen_basis, .false., overlap)

end subroutine test_overlap_f2

subroutine test_overlap_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 2
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 9.99999999881495E-1_wp, 6.61346655776026E-1_wp, 6.61346655776026E-1_wp,&
      & 9.99999999881495E-1_wp],shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_overlap_mol(error, mol, make_gen_basis, .false., overlap)

end subroutine test_overlap_h2

subroutine test_overlap_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 5
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000600E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 3.99089038384911E-1_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 4.65790780903622E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 9.99999999925689E-1_wp, 0.00000000000000E+0_wp, 3.99089038384911E-1_wp,&
      & 0.00000000000000E+0_wp, 4.65790780903622E-1_wp, 0.00000000000000E+0_wp,&
      & 9.99999999881495E-1_wp],shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_overlap_mol(error, mol, make_gen_basis, .false., overlap)

end subroutine test_overlap_lih

subroutine test_overlap_qvszp_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 43
   real(wp), parameter :: overlap(nao, nao) = reshape([&
      & 1.00000000000000E+00_wp,-1.69197570393807E-17_wp,-5.99581626875016E-17_wp, &
      & 4.80516738344421E-17_wp, 5.62681572696870E-33_wp,-9.58319782095767E-33_wp, &
      & 0.00000000000000E+00_wp,-3.27240132584904E-33_wp, 0.00000000000000E+00_wp, &
      & 1.85014596230672E-33_wp, 3.49139538903901E-49_wp, 1.85014596230672E-33_wp, &
      &-1.97348902646050E-32_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 2.71640292351996E-01_wp, 1.12463267268395E-01_wp, &
      & 7.03448491839345E-02_wp, 7.65751470885278E-02_wp, 5.48274140070622E-03_wp, &
      & 5.03665525449538E-03_wp,-1.58329216687767E-03_wp, 3.42940967584335E-03_wp, &
      &-2.15958353166773E-03_wp, 2.71520165109221E-01_wp,-1.16844510294259E-01_wp, &
      &-5.87523732166724E-02_wp, 8.00151496696516E-02_wp,-6.00032867731247E-03_wp, &
      & 4.40583503658527E-03_wp,-2.43654952707911E-03_wp,-3.01711692731095E-03_wp, &
      &-2.32656622881648E-03_wp, 2.71960093496385E-01_wp, 3.69310700214585E-02_wp, &
      &-5.82987274843623E-02_wp,-1.36620107880504E-01_wp,-3.15805596256058E-03_wp, &
      &-1.34761014903253E-03_wp,-2.39073219428602E-03_wp, 4.98525073426546E-03_wp, &
      & 5.41449904305805E-03_wp,-1.69197570393807E-17_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 5.14277548501503E-34_wp,-1.12796983051776E-17_wp, &
      &-4.29449030052634E-18_wp, 3.73851650792758E-18_wp, 0.00000000000000E+00_wp, &
      & 6.47530053666554E-18_wp, 0.00000000000000E+00_wp, 1.10622338413461E-32_wp, &
      & 0.00000000000000E+00_wp,-1.03536604457902E-32_wp, 7.18742570233806E-34_wp, &
      &-1.33665181596106E-32_wp, 2.78367800472584E-33_wp,-2.86186129356754E-01_wp, &
      & 7.65422377302736E-02_wp,-6.28658493087070E-02_wp,-6.84337476518306E-02_wp, &
      & 1.46212313795836E-02_wp, 1.34316205111712E-02_wp,-1.16485318540639E-02_wp, &
      & 8.38206255027617E-04_wp,-2.00334127727306E-02_wp, 2.97076450678070E-01_wp, &
      & 6.86437027853556E-02_wp,-5.45295352799983E-02_wp, 7.42640456540307E-02_wp, &
      & 1.53233105239048E-02_wp,-1.12513800515564E-02_wp, 1.22854141999639E-02_wp, &
      & 7.28954899753139E-04_wp, 2.08214407972141E-02_wp,-9.40524924437796E-02_wp, &
      & 1.66395373175104E-01_wp, 1.71049702610842E-02_wp, 4.00846293427067E-02_wp, &
      &-2.39637014091444E-02_wp,-1.02258248777719E-02_wp,-3.90605105295778E-03_wp, &
      & 4.36203320089825E-04_wp,-5.92939585074226E-03_wp,-5.99581626875016E-17_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp,-6.47530053666554E-18_wp,-4.95885026208224E-18_wp, &
      &-1.12796983051776E-17_wp, 0.00000000000000E+00_wp,-1.07971242649150E-48_wp, &
      &-2.27286357326236E-33_wp, 1.69074567081486E-32_wp, 0.00000000000000E+00_wp, &
      & 1.39927419792190E-32_wp, 0.00000000000000E+00_wp, 1.07971242649150E-48_wp, &
      &-1.79007071349708E-01_wp,-6.28658493087069E-02_wp, 1.37726441565078E-01_wp, &
      &-4.28047466038034E-02_wp, 8.38206255027629E-04_wp, 2.02755819182050E-02_wp, &
      & 1.38459566041498E-02_wp, 1.38054469286112E-02_wp,-3.30159001164187E-04_wp, &
      & 1.49377548505845E-01_wp,-5.45295352799983E-02_wp, 1.49671163139505E-01_wp, &
      & 3.73418393029971E-02_wp, 7.28954899753139E-04_wp,-2.07945739169201E-02_wp, &
      &-1.14668148754686E-02_wp, 1.42401293829613E-02_wp, 2.82644825525715E-04_wp, &
      & 1.48469584634808E-01_wp, 1.71049702610842E-02_wp, 1.50229426898005E-01_wp, &
      &-6.32768799009697E-02_wp, 4.36203320089823E-04_wp, 6.58929517447636E-03_wp, &
      &-1.13413923700238E-02_wp,-2.43759581585471E-02_wp,-7.47872262938034E-04_wp, &
      & 4.80516738344421E-17_wp, 5.14277548501503E-34_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp,-6.47530053666554E-18_wp, 0.00000000000000E+00_wp, &
      & 6.51233685287207E-18_wp,-4.29449030052634E-18_wp,-1.12796983051776E-17_wp, &
      &-2.78367800472584E-33_wp, 1.33665181596106E-32_wp, 7.18742570233806E-34_wp, &
      &-8.56876948787715E-33_wp, 0.00000000000000E+00_wp, 1.10622338413462E-32_wp, &
      & 0.00000000000000E+00_wp,-1.94861357690152E-01_wp,-6.84337476518306E-02_wp, &
      &-4.28047466038034E-02_wp, 1.30452638983255E-01_wp, 2.04180181502975E-02_wp, &
      & 8.38206255027626E-04_wp,-7.93137227608380E-03_wp, 1.27713025088428E-02_wp, &
      & 1.29217553049478E-02_wp,-2.03438027207878E-01_wp, 7.42640456540307E-02_wp, &
      & 3.73418393029971E-02_wp, 1.26233954895127E-01_wp,-2.12520947060716E-02_wp, &
      & 7.28954899753138E-04_wp,-8.41305469540806E-03_wp,-1.06860904005050E-02_wp, &
      & 1.34886570383591E-02_wp, 3.47930933401967E-01_wp, 4.00846293427067E-02_wp, &
      &-6.32768799009697E-02_wp, 2.89448721134034E-02_wp, 7.42537820441489E-03_wp, &
      & 4.36203320089827E-04_wp, 1.44497604843774E-02_wp,-1.17215694036479E-02_wp, &
      &-2.54399760655455E-02_wp, 5.62681572696870E-33_wp,-1.12796983051776E-17_wp, &
      & 0.00000000000000E+00_wp,-6.47530053666554E-18_wp, 1.00000000000000E+00_wp, &
      & 2.68345343283949E-34_wp, 1.31995475669199E-33_wp, 1.22511181451474E-33_wp, &
      & 0.00000000000000E+00_wp,-7.06028367843927E-18_wp, 1.95844915502933E-17_wp, &
      & 1.82295740707302E-18_wp,-2.69012192599270E-48_wp, 2.34329103691514E-17_wp, &
      & 0.00000000000000E+00_wp, 9.07552716128947E-17_wp, 1.70704186868856E-01_wp, &
      & 9.98283434556062E-02_wp, 1.23166359697416E-01_wp,-8.50733399053298E-03_wp, &
      & 2.89571170460196E-02_wp, 7.00627156995697E-02_wp, 2.75252469034757E-02_wp, &
      & 2.22211300508841E-02_wp,-4.39114994972613E-02_wp,-1.84778662615733E-01_wp, &
      & 1.20053023357136E-01_wp, 1.11282240180291E-01_wp,-3.68657002091701E-03_wp, &
      & 3.93059495632419E-02_wp,-6.78524889755777E-02_wp,-1.15115579619026E-02_wp, &
      & 2.46537787898557E-02_wp, 5.06831870845617E-02_wp,-1.00336855750964E-01_wp, &
      & 1.35742399194017E-01_wp, 6.00159024456601E-02_wp, 9.36731736887267E-02_wp, &
      &-5.36255231615542E-02_wp,-3.16478193052413E-02_wp,-6.03726229855204E-03_wp, &
      &-4.77276088174359E-02_wp,-6.59477421558427E-02_wp,-9.58319782095767E-33_wp, &
      &-4.29449030052634E-18_wp,-6.47530053666554E-18_wp, 0.00000000000000E+00_wp, &
      & 2.68345343283949E-34_wp, 1.00000000000000E+00_wp, 7.07318635897479E-34_wp, &
      &-1.14311435114137E-33_wp,-1.22511181451474E-33_wp, 3.12720012523314E-33_wp, &
      &-5.76469748382549E-18_wp, 2.47726400460995E-17_wp, 5.73986735927941E-17_wp, &
      & 1.64735651613881E-48_wp, 7.41013689730956E-17_wp,-2.12672811742396E-48_wp, &
      & 1.56815373354403E-01_wp, 9.17061217858072E-02_wp,-2.94369704500862E-02_wp, &
      & 1.23166359697416E-01_wp, 7.00627156995698E-02_wp, 1.70513728755861E-02_wp, &
      &-5.70042503664237E-02_wp, 1.23397864595546E-02_wp, 7.17137135298872E-03_wp, &
      & 1.35676618656563E-01_wp,-8.81508072307406E-02_wp, 6.61585443518252E-02_wp, &
      & 1.11282240180291E-01_wp,-6.78524889755777E-02_wp,-3.28094341533597E-03_wp, &
      &-6.26875333115744E-02_wp, 7.67644703164060E-03_wp, 3.85781955243907E-03_wp, &
      &-4.28158863348252E-02_wp, 5.79241903805787E-02_wp,-2.13609749530813E-02_wp, &
      & 6.00159024456601E-02_wp,-3.16478193052413E-02_wp, 7.03473190075172E-03_wp, &
      & 1.99265001040424E-02_wp, 4.53565142340132E-03_wp,-4.11332714987785E-02_wp, &
      & 0.00000000000000E+00_wp, 3.73851650792758E-18_wp,-4.95885026208224E-18_wp, &
      & 6.51233685287207E-18_wp, 1.31995475669199E-33_wp, 7.07318635897479E-34_wp, &
      & 1.00000000000000E+00_wp, 1.54929256180770E-34_wp, 0.00000000000000E+00_wp, &
      & 2.50176010018651E-32_wp, 0.00000000000000E+00_wp,-8.11739826571558E-17_wp, &
      & 2.62753526467356E-17_wp,-6.31490969816898E-18_wp, 0.00000000000000E+00_wp, &
      & 3.90900015654143E-34_wp,-4.92955224712704E-02_wp, 2.54563919960000E-02_wp, &
      &-1.38548535508843E-01_wp, 1.73330102244392E-02_wp, 2.75252469034757E-02_wp, &
      &-5.70042503664237E-02_wp,-3.87375227679548E-02_wp,-3.88136407780398E-02_wp, &
      &-1.08418518353207E-02_wp,-7.50329502303753E-02_wp, 4.49644282170747E-03_wp, &
      & 1.31043068444519E-01_wp,-3.07916516106645E-03_wp,-1.15115579619026E-02_wp, &
      &-6.26875333115744E-02_wp,-1.73980857687373E-02_wp, 4.29284383811856E-02_wp, &
      &-4.46348915793432E-03_wp,-7.59576632463371E-02_wp,-1.66254425274986E-03_wp, &
      & 1.31051968499253E-01_wp, 6.15029499645749E-03_wp,-6.03726229855206E-03_wp, &
      & 1.99265001040424E-02_wp,-1.68184321679367E-02_wp,-7.37146416909488E-02_wp, &
      & 1.03509093333787E-02_wp,-3.27240132584904E-33_wp, 0.00000000000000E+00_wp, &
      &-1.12796983051776E-17_wp,-4.29449030052634E-18_wp, 1.22511181451474E-33_wp, &
      &-1.14311435114137E-33_wp, 1.54929256180770E-34_wp, 1.00000000000000E+00_wp, &
      & 2.68345343283949E-34_wp, 2.12672811742396E-48_wp,-7.41013689730957E-17_wp, &
      & 1.64735651613881E-48_wp, 4.46531547015598E-18_wp, 2.47726400460995E-17_wp, &
      &-5.76469748382549E-18_wp,-3.12720012523314E-33_wp, 1.06774065630676E-01_wp, &
      & 1.23166359697416E-01_wp,-2.00433474574100E-02_wp,-5.32126747743511E-03_wp, &
      & 2.22211300508841E-02_wp, 1.23397864595546E-02_wp,-3.88136407780398E-02_wp, &
      & 7.33039657022230E-03_wp,-5.98154810786167E-02_wp,-9.29113821534032E-02_wp, &
      & 1.11282240180291E-01_wp,-4.53053876036286E-02_wp,-1.85370059074957E-03_wp, &
      & 2.46537787898556E-02_wp, 7.67644703164059E-03_wp, 4.29284383811855E-02_wp, &
      & 2.67198462630366E-03_wp, 5.36114109247963E-02_wp, 1.58389968302148E-01_wp, &
      & 6.00159024456601E-02_wp, 7.90212333633180E-02_wp,-1.47870798823357E-01_wp, &
      &-4.77276088174359E-02_wp, 4.53565142340133E-03_wp,-7.37146416909488E-02_wp, &
      &-8.51805157287561E-03_wp, 5.60424300043609E-02_wp, 0.00000000000000E+00_wp, &
      & 6.47530053666554E-18_wp, 0.00000000000000E+00_wp,-1.12796983051776E-17_wp, &
      & 0.00000000000000E+00_wp,-1.22511181451474E-33_wp, 0.00000000000000E+00_wp, &
      & 2.68345343283949E-34_wp, 1.00000000000000E+00_wp,-9.07552716128947E-17_wp, &
      & 0.00000000000000E+00_wp,-2.34329103691514E-17_wp, 0.00000000000000E+00_wp, &
      & 1.82295740707302E-18_wp, 1.95844915502933E-17_wp,-7.06028367843926E-18_wp, &
      &-6.72382525101829E-02_wp, 6.50214007740149E-02_wp,-4.85136946316211E-02_wp, &
      &-1.49893345269856E-01_wp,-4.39114994972613E-02_wp, 7.17137135298871E-03_wp, &
      &-1.08418518353207E-02_wp,-5.98154810786167E-02_wp,-6.52290168122297E-02_wp, &
      &-7.16460413032241E-02_wp,-6.20570469522165E-02_wp, 4.31485533199955E-02_wp, &
      &-1.60025140456901E-01_wp, 5.06831870845617E-02_wp, 3.85781955243902E-03_wp, &
      &-4.46348915793430E-03_wp, 5.36114109247964E-02_wp,-7.17565927775726E-02_wp, &
      & 1.72027923471804E-01_wp, 1.12154537526214E-01_wp,-1.02897494601968E-01_wp, &
      &-6.73737892461001E-02_wp,-6.59477421558427E-02_wp,-4.11332714987785E-02_wp, &
      & 1.03509093333786E-02_wp, 5.60424300043609E-02_wp, 2.09775110980944E-02_wp, &
      & 1.85014596230672E-33_wp, 0.00000000000000E+00_wp,-1.07971242649150E-48_wp, &
      &-2.78367800472584E-33_wp,-7.06028367843927E-18_wp, 3.12720012523314E-33_wp, &
      & 2.50176010018651E-32_wp, 2.12672811742396E-48_wp,-9.07552716128947E-17_wp, &
      & 1.00000000000000E+00_wp,-1.15800615204402E-32_wp,-1.34520575938223E-16_wp, &
      &-2.21079215219747E-48_wp,-1.60919326636305E-32_wp,-9.90795173428150E-33_wp, &
      & 3.97942587395545E-47_wp,-8.87271221601863E-03_wp,-3.09654780788219E-02_wp, &
      &-7.43465509581947E-03_wp, 4.03455271857084E-02_wp, 6.33223591146903E-03_wp, &
      &-1.90767803600434E-02_wp,-3.52809157214824E-03_wp, 1.80847699207693E-02_wp, &
      & 2.84103623255625E-02_wp, 1.03166611962926E-02_wp,-3.46732517765560E-02_wp, &
      &-7.21339766751295E-03_wp,-4.26017054342602E-02_wp, 6.79736167259399E-03_wp, &
      & 1.80558583064038E-02_wp, 2.29483867913033E-03_wp, 1.56108390830926E-02_wp, &
      &-3.37969956237433E-02_wp,-3.23280067060678E-02_wp, 3.46135196140600E-02_wp, &
      & 2.24627183919275E-02_wp, 2.41519984031383E-02_wp,-1.20517061069994E-02_wp, &
      &-1.11470873876836E-02_wp,-7.06986288266467E-03_wp,-2.27065694199979E-02_wp, &
      &-2.84653330969374E-02_wp, 3.49139538903901E-49_wp, 1.10622338413461E-32_wp, &
      &-2.27286357326236E-33_wp, 1.33665181596106E-32_wp, 1.95844915502933E-17_wp, &
      &-5.76469748382549E-18_wp, 0.00000000000000E+00_wp,-7.41013689730957E-17_wp, &
      & 0.00000000000000E+00_wp,-1.15800615204402E-32_wp, 1.00000000000000E+00_wp, &
      &-8.96987708334444E-33_wp,-3.21838653272609E-32_wp,-7.67466641237984E-33_wp, &
      & 0.00000000000000E+00_wp, 9.90795173428151E-33_wp,-4.73660932709948E-02_wp, &
      &-3.87145441954640E-02_wp,-1.39172471953937E-04_wp,-6.87222204880001E-03_wp, &
      &-2.08260457665859E-02_wp,-1.56776548042670E-02_wp, 1.58174413977799E-02_wp, &
      &-8.58333743425472E-03_wp, 2.15747204908677E-02_wp,-4.27747085919899E-02_wp, &
      & 3.79563333930383E-02_wp,-1.28973512282152E-02_wp,-9.30137356639771E-03_wp, &
      & 2.23784051718079E-02_wp,-4.62094080745720E-03_wp, 2.17475896955325E-02_wp, &
      & 4.08897524068078E-03_wp, 2.07659262335779E-02_wp,-2.30896339987546E-02_wp, &
      & 2.65554267780228E-02_wp,-7.21700332353314E-03_wp, 2.76714656194347E-02_wp, &
      &-1.55601171838614E-02_wp, 2.75700778549133E-03_wp, 1.18768916878722E-02_wp, &
      &-2.47390282508854E-03_wp,-2.68986291447388E-02_wp, 1.85014596230672E-33_wp, &
      & 0.00000000000000E+00_wp, 1.69074567081486E-32_wp, 7.18742570233806E-34_wp, &
      & 1.82295740707302E-18_wp, 2.47726400460995E-17_wp,-8.11739826571558E-17_wp, &
      & 1.64735651613881E-48_wp,-2.34329103691514E-17_wp,-1.34520575938223E-16_wp, &
      &-8.96987708334444E-33_wp, 1.00000000000000E+00_wp,-3.96318069371260E-33_wp, &
      & 2.49295148858128E-32_wp, 7.67466641237983E-33_wp, 1.60919326636305E-32_wp, &
      &-1.78219480839546E-03_wp,-1.98249860284091E-02_wp, 4.44635747164258E-02_wp, &
      &-1.41324121957292E-02_wp,-1.86596872921103E-02_wp, 1.49962751770109E-02_wp, &
      & 2.36951984957053E-02_wp, 1.76904901292571E-02_wp, 3.44816543864968E-03_wp, &
      &-8.98833973673464E-03_wp,-1.17908674490106E-02_wp, 4.60414462317829E-02_wp, &
      & 4.97717098610996E-03_wp,-1.07989162672097E-02_wp,-2.18241822493265E-02_wp, &
      &-1.46101566731171E-02_wp, 2.27818699288189E-02_wp,-3.49031703133335E-03_wp, &
      & 2.94912936553633E-03_wp,-5.38016899490450E-03_wp,-1.46043900971410E-02_wp, &
      & 2.55351514429919E-03_wp, 9.38158733271537E-04_wp, 8.15193580412787E-03_wp, &
      & 4.54596058457116E-03_wp, 1.24342367450678E-02_wp,-3.54028932533078E-03_wp, &
      &-1.97348902646050E-32_wp,-1.03536604457902E-32_wp, 0.00000000000000E+00_wp, &
      &-8.56876948787715E-33_wp,-2.69012192599270E-48_wp, 5.73986735927941E-17_wp, &
      & 2.62753526467356E-17_wp, 4.46531547015598E-18_wp, 0.00000000000000E+00_wp, &
      &-2.21079215219747E-48_wp,-3.21838653272609E-32_wp,-3.96318069371260E-33_wp, &
      & 1.00000000000000E+00_wp,-4.63202460817608E-33_wp, 0.00000000000000E+00_wp, &
      & 1.65809411414810E-48_wp, 3.24058962055940E-02_wp, 1.52689570210902E-02_wp, &
      & 1.18305883730963E-02_wp, 1.03964846316267E-02_wp, 6.48042478460700E-03_wp, &
      & 1.90127212748087E-02_wp,-1.26905584426564E-02_wp, 1.29455773741406E-02_wp, &
      &-2.55255858706855E-03_wp,-3.14589601176654E-02_wp, 1.93989802055032E-02_wp, &
      &-1.32426681021874E-03_wp,-1.32844264627656E-02_wp, 1.14264510979431E-02_wp, &
      &-1.07140661502377E-02_wp, 2.01491235766777E-02_wp, 7.33699516068696E-03_wp, &
      & 4.43048983970076E-03_wp,-3.14829673682181E-02_wp,-6.16924951477627E-03_wp, &
      &-1.74920521618782E-03_wp, 2.28220718695871E-02_wp, 6.25251359789685E-03_wp, &
      & 3.31254163982971E-03_wp, 2.04280163304608E-02_wp,-1.22541750328176E-02_wp, &
      &-1.07199585105105E-02_wp, 0.00000000000000E+00_wp, 7.18742570233806E-34_wp, &
      & 1.39927419792190E-32_wp, 0.00000000000000E+00_wp, 2.34329103691514E-17_wp, &
      & 1.64735651613881E-48_wp,-6.31490969816898E-18_wp, 2.47726400460995E-17_wp, &
      & 1.82295740707302E-18_wp,-1.60919326636305E-32_wp,-7.67466641237984E-33_wp, &
      & 2.49295148858128E-32_wp,-4.63202460817608E-33_wp, 1.00000000000000E+00_wp, &
      &-8.96987708334445E-33_wp, 4.48401919794077E-17_wp,-1.21347914664084E-03_wp, &
      &-1.41324121957292E-02_wp, 3.02748164506590E-02_wp,-8.69182383450863E-03_wp, &
      &-9.84542294601399E-03_wp, 1.76904901292571E-02_wp, 1.61338306646406E-02_wp, &
      & 1.06014784284417E-03_wp, 9.60821199164087E-03_wp, 6.15521728410922E-03_wp, &
      & 4.97717098610996E-03_wp,-3.15291938146309E-02_wp,-7.93117290423926E-03_wp, &
      & 6.89125459326437E-03_wp, 2.27818699288189E-02_wp, 1.00050389184092E-02_wp, &
      &-4.15730696428098E-03_wp, 3.68961344026113E-03_wp,-1.09097941608253E-02_wp, &
      & 2.55351514429918E-03_wp, 5.40264159538580E-02_wp,-1.41361919844603E-02_wp, &
      &-6.80798170999700E-03_wp, 1.24342367450678E-02_wp,-1.68169951513391E-02_wp, &
      &-3.44851572401250E-02_wp, 1.11501042348473E-02_wp, 0.00000000000000E+00_wp, &
      &-1.33665181596106E-32_wp, 0.00000000000000E+00_wp, 1.10622338413462E-32_wp, &
      & 0.00000000000000E+00_wp, 7.41013689730956E-17_wp, 0.00000000000000E+00_wp, &
      &-5.76469748382549E-18_wp, 1.95844915502933E-17_wp,-9.90795173428150E-33_wp, &
      & 0.00000000000000E+00_wp, 7.67466641237983E-33_wp, 0.00000000000000E+00_wp, &
      &-8.96987708334445E-33_wp, 1.00000000000000E+00_wp,-1.15800615204402E-32_wp, &
      & 1.86569140347025E-02_wp,-1.13388893495405E-02_wp, 5.48183028392548E-05_wp, &
      & 4.17558546060307E-02_wp, 2.15747204908677E-02_wp, 3.32185778425088E-03_wp, &
      &-6.23029310690556E-03_wp, 7.57152293456893E-03_wp, 2.54497483960324E-02_wp, &
      &-1.65854568657006E-02_wp,-8.36777003079139E-03_wp,-5.00081635897833E-03_wp, &
      &-3.73170014599802E-02_wp, 2.07659262335779E-02_wp,-3.07044231541223E-03_wp, &
      & 8.43240603387193E-03_wp,-2.81830113107483E-04_wp,-2.31261633888675E-02_wp, &
      & 3.95872659233849E-02_wp, 2.73507326766188E-02_wp, 1.23735798390771E-02_wp, &
      &-2.77419228427968E-02_wp,-2.68986291447388E-02_wp,-2.55315572948441E-04_wp, &
      &-2.03629762869517E-02_wp, 5.45027214313180E-03_wp, 1.48688146646377E-02_wp, &
      & 0.00000000000000E+00_wp, 2.78367800472584E-33_wp, 1.07971242649150E-48_wp, &
      & 0.00000000000000E+00_wp, 9.07552716128947E-17_wp,-2.12672811742396E-48_wp, &
      & 3.90900015654143E-34_wp,-3.12720012523314E-33_wp,-7.06028367843926E-18_wp, &
      & 3.97942587395545E-47_wp, 9.90795173428151E-33_wp, 1.60919326636305E-32_wp, &
      & 1.65809411414810E-48_wp, 4.48401919794077E-17_wp,-1.15800615204402E-32_wp, &
      & 1.00000000000000E+00_wp, 3.92063216163170E-02_wp, 4.08306349686239E-03_wp, &
      & 3.28519027436555E-02_wp, 1.66821454033021E-02_wp, 2.09689647083699E-02_wp, &
      & 1.71466011097604E-02_wp, 1.55897643811477E-02_wp, 1.87072473534496E-02_wp, &
      &-1.26546050185148E-03_wp, 4.39506017098200E-02_wp,-8.68941884558286E-03_wp, &
      &-3.07302102712587E-02_wp, 2.15241086275338E-02_wp,-2.62400842432184E-02_wp, &
      & 1.82997764559348E-02_wp, 9.77637424121205E-03_wp,-1.90984663399790E-02_wp, &
      &-1.62238688233457E-03_wp, 3.19021256333458E-02_wp, 4.25304632316621E-02_wp, &
      &-2.21668001594183E-02_wp,-3.10357127665506E-03_wp,-3.08724369198555E-02_wp, &
      &-2.12187080366463E-02_wp, 6.97672627774367E-03_wp, 1.36980314278988E-02_wp, &
      & 3.14704172073689E-03_wp, 2.71640292351996E-01_wp,-2.86186129356754E-01_wp, &
      &-1.79007071349708E-01_wp,-1.94861357690152E-01_wp, 1.70704186868856E-01_wp, &
      & 1.56815373354403E-01_wp,-4.92955224712704E-02_wp, 1.06774065630676E-01_wp, &
      &-6.72382525101829E-02_wp,-8.87271221601863E-03_wp,-4.73660932709948E-02_wp, &
      &-1.78219480839546E-03_wp, 3.24058962055940E-02_wp,-1.21347914664084E-03_wp, &
      & 1.86569140347025E-02_wp, 3.92063216163170E-02_wp, 1.00000000000000E+00_wp, &
      &-3.24870256971979E-16_wp,-5.05897538405822E-17_wp,-2.12011316571604E-16_wp, &
      &-1.61061419113769E-32_wp,-1.29140636663106E-31_wp, 0.00000000000000E+00_wp, &
      & 1.20922074604819E-32_wp, 0.00000000000000E+00_wp, 1.71162831672690E-03_wp, &
      &-6.05377992113713E-03_wp,-3.40829260949412E-03_wp, 8.97569608799371E-05_wp, &
      &-1.30502848732678E-05_wp, 4.95551420739942E-04_wp,-9.30672388945664E-05_wp, &
      &-7.34734167161412E-06_wp,-4.40000332916343E-04_wp, 1.90418230344655E-03_wp, &
      &-2.21142038639408E-03_wp,-3.76575523760944E-03_wp,-6.24060272675448E-03_wp, &
      & 3.10876048741532E-04_wp, 1.87591352959682E-04_wp,-1.00621858844341E-04_wp, &
      & 5.29380626134436E-04_wp, 3.83563276743010E-04_wp, 1.12463267268395E-01_wp, &
      & 7.65422377302736E-02_wp,-6.28658493087069E-02_wp,-6.84337476518306E-02_wp, &
      & 9.98283434556062E-02_wp, 9.17061217858072E-02_wp, 2.54563919960000E-02_wp, &
      & 1.23166359697416E-01_wp, 6.50214007740149E-02_wp,-3.09654780788219E-02_wp, &
      &-3.87145441954640E-02_wp,-1.98249860284091E-02_wp, 1.52689570210902E-02_wp, &
      &-1.41324121957292E-02_wp,-1.13388893495405E-02_wp, 4.08306349686239E-03_wp, &
      &-3.24870256971979E-16_wp, 1.00000000000000E+00_wp, 6.43479108906961E-32_wp, &
      &-1.51427681500289E-32_wp, 4.89204577231437E-17_wp,-4.20964955515460E-17_wp, &
      &-1.53550765102059E-16_wp,-2.86239545746042E-47_wp,-2.65957726697841E-16_wp, &
      & 6.04732062392753E-03_wp,-1.47041207437914E-02_wp,-9.75933422332976E-03_wp, &
      & 2.57010849848264E-04_wp,-6.51497161910206E-05_wp, 2.47389499407760E-03_wp, &
      &-2.07149974319160E-04_wp,-4.12364049545137E-05_wp,-1.92355606367149E-03_wp, &
      & 2.21107213999637E-03_wp, 8.11755796279226E-04_wp,-3.54265554044613E-03_wp, &
      &-5.87088231461763E-03_wp,-8.50486442997838E-06_wp,-5.13207444451533E-06_wp, &
      &-6.62141698033521E-05_wp, 9.71741044910374E-04_wp, 9.09305302729587E-04_wp, &
      & 7.03448491839345E-02_wp,-6.28658493087070E-02_wp, 1.37726441565078E-01_wp, &
      &-4.28047466038034E-02_wp, 1.23166359697416E-01_wp,-2.94369704500862E-02_wp, &
      &-1.38548535508843E-01_wp,-2.00433474574100E-02_wp,-4.85136946316211E-02_wp, &
      &-7.43465509581947E-03_wp,-1.39172471953937E-04_wp, 4.44635747164258E-02_wp, &
      & 1.18305883730963E-02_wp, 3.02748164506590E-02_wp, 5.48183028392548E-05_wp, &
      & 3.28519027436555E-02_wp,-5.05897538405822E-17_wp, 6.43479108906961E-32_wp, &
      & 1.00000000000000E+00_wp, 1.74525095580671E-32_wp,-2.86239545746042E-47_wp, &
      & 2.65957726697841E-16_wp,-4.86088460772499E-17_wp, 4.89204577231437E-17_wp, &
      & 0.00000000000000E+00_wp, 3.40465600967899E-03_wp,-9.75933422332976E-03_wp, &
      &-2.86420408011397E-03_wp, 1.44697724646240E-04_wp,-4.12364049545137E-05_wp, &
      & 1.01993488185775E-03_wp,-6.48972176334236E-04_wp,-1.51221644135675E-05_wp, &
      &-1.39031692292213E-03_wp, 3.76516222024187E-03_wp,-3.54265554044613E-03_wp, &
      &-3.14050988116331E-03_wp,-9.99733291855449E-03_wp, 9.71741044910373E-04_wp, &
      & 3.81146552296546E-04_wp,-7.18068099254989E-04_wp, 1.07559115769634E-03_wp, &
      & 1.19894787919603E-03_wp, 7.65751470885278E-02_wp,-6.84337476518306E-02_wp, &
      &-4.28047466038034E-02_wp, 1.30452638983255E-01_wp,-8.50733399053298E-03_wp, &
      & 1.23166359697416E-01_wp, 1.73330102244392E-02_wp,-5.32126747743511E-03_wp, &
      &-1.49893345269856E-01_wp, 4.03455271857084E-02_wp,-6.87222204880001E-03_wp, &
      &-1.41324121957292E-02_wp, 1.03964846316267E-02_wp,-8.69182383450863E-03_wp, &
      & 4.17558546060307E-02_wp, 1.66821454033021E-02_wp,-2.12011316571604E-16_wp, &
      &-1.51427681500289E-32_wp, 1.74525095580671E-32_wp, 1.00000000000000E+00_wp, &
      & 2.65957726697841E-16_wp,-2.86239545746042E-47_wp,-2.82442394353367E-17_wp, &
      &-4.20964955515460E-17_wp, 4.89204577231436E-17_wp,-8.96611914772649E-05_wp, &
      & 2.57010849848263E-04_wp, 1.44697724646240E-04_wp, 2.62651388901871E-03_wp, &
      &-5.44826901506361E-04_wp,-4.12364049545137E-05_wp, 3.07132938155310E-06_wp, &
      &-3.06738851766667E-04_wp, 4.47078532418239E-05_wp, 6.23961997945212E-03_wp, &
      &-5.87088231461763E-03_wp,-9.99733291855449E-03_wp,-1.36753998857564E-02_wp, &
      & 1.40513816130653E-03_wp, 9.71741044910373E-04_wp,-1.86855620562661E-04_wp, &
      & 2.39276368394753E-03_wp, 1.40773945062066E-03_wp, 5.48274140070622E-03_wp, &
      & 1.46212313795836E-02_wp, 8.38206255027629E-04_wp, 2.04180181502975E-02_wp, &
      & 2.89571170460196E-02_wp, 7.00627156995698E-02_wp, 2.75252469034757E-02_wp, &
      & 2.22211300508841E-02_wp,-4.39114994972613E-02_wp, 6.33223591146903E-03_wp, &
      &-2.08260457665859E-02_wp,-1.86596872921103E-02_wp, 6.48042478460700E-03_wp, &
      &-9.84542294601399E-03_wp, 2.15747204908677E-02_wp, 2.09689647083699E-02_wp, &
      &-1.61061419113769E-32_wp, 4.89204577231437E-17_wp,-2.86239545746042E-47_wp, &
      & 2.65957726697841E-16_wp, 1.00000000000000E+00_wp, 8.75870038752334E-32_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.30642898219517E-05_wp, 6.52724778871420E-05_wp, 4.13143676691053E-05_wp, &
      & 5.45884768170620E-04_wp,-1.44197580318789E-06_wp,-3.00514505778023E-07_wp, &
      & 3.29424435896269E-08_wp,-8.45903480310013E-07_wp, 2.78024017274903E-07_wp, &
      & 3.10924905667185E-04_wp, 8.51213619928546E-06_wp,-9.72009902682531E-04_wp, &
      &-1.40552518683217E-03_wp, 3.03164987712131E-06_wp, 1.91875419651281E-06_wp, &
      &-8.68877459996167E-07_wp, 8.20652768913813E-06_wp, 6.23653211157604E-06_wp, &
      & 5.03665525449538E-03_wp, 1.34316205111712E-02_wp, 2.02755819182050E-02_wp, &
      & 8.38206255027626E-04_wp, 7.00627156995697E-02_wp, 1.70513728755861E-02_wp, &
      &-5.70042503664237E-02_wp, 1.23397864595546E-02_wp, 7.17137135298871E-03_wp, &
      &-1.90767803600434E-02_wp,-1.56776548042670E-02_wp, 1.49962751770109E-02_wp, &
      & 1.90127212748087E-02_wp, 1.76904901292571E-02_wp, 3.32185778425088E-03_wp, &
      & 1.71466011097604E-02_wp,-1.29140636663106E-31_wp,-4.20964955515460E-17_wp, &
      & 2.65957726697841E-16_wp,-2.86239545746042E-47_wp, 8.75870038752334E-32_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 4.96083223094062E-04_wp,-2.47855655767690E-03_wp, &
      &-1.02183539811886E-03_wp, 4.13143676691057E-05_wp,-3.00514505778022E-07_wp, &
      & 9.96138524408819E-06_wp,-2.72409461592086E-06_wp,-1.53891219610808E-07_wp, &
      &-9.70670982165104E-06_wp, 1.87620834602998E-04_wp, 5.13646243467491E-06_wp, &
      &-3.81250263554302E-04_wp,-9.72009902682532E-04_wp, 1.91875419651281E-06_wp, &
      & 1.00972417245683E-06_wp,-1.21871184009155E-06_wp, 4.52957647647255E-06_wp, &
      & 4.16421483137780E-06_wp,-1.58329216687767E-03_wp,-1.16485318540639E-02_wp, &
      & 1.38459566041498E-02_wp,-7.93137227608380E-03_wp, 2.75252469034757E-02_wp, &
      &-5.70042503664237E-02_wp,-3.87375227679548E-02_wp,-3.88136407780398E-02_wp, &
      &-1.08418518353207E-02_wp,-3.52809157214824E-03_wp, 1.58174413977799E-02_wp, &
      & 2.36951984957053E-02_wp,-1.26905584426564E-02_wp, 1.61338306646406E-02_wp, &
      &-6.23029310690556E-03_wp, 1.55897643811477E-02_wp, 0.00000000000000E+00_wp, &
      &-1.53550765102059E-16_wp,-4.86088460772499E-17_wp,-2.82442394353367E-17_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      & 5.05683802648788E-32_wp, 0.00000000000000E+00_wp,-9.31671142549513E-05_wp, &
      & 2.07525566721164E-04_wp, 6.50217215235351E-04_wp,-3.07689813908525E-06_wp, &
      & 3.29424435896268E-08_wp,-2.72409461592087E-06_wp,-6.61949337710925E-07_wp, &
      & 4.03890556081740E-08_wp, 1.11067967383642E-06_wp,-1.00637672461044E-04_wp, &
      & 6.62314833876290E-05_wp, 7.18270198948683E-04_wp, 1.86904479297035E-04_wp, &
      &-8.68877459996168E-07_wp,-1.21871184009155E-06_wp,-9.96852139011951E-07_wp, &
      &-3.43919070258945E-06_wp,-1.07203332966111E-06_wp, 3.42940967584335E-03_wp, &
      & 8.38206255027617E-04_wp, 1.38054469286112E-02_wp, 1.27713025088428E-02_wp, &
      & 2.22211300508841E-02_wp, 1.23397864595546E-02_wp,-3.88136407780398E-02_wp, &
      & 7.33039657022230E-03_wp,-5.98154810786167E-02_wp, 1.80847699207693E-02_wp, &
      &-8.58333743425472E-03_wp, 1.76904901292571E-02_wp, 1.29455773741406E-02_wp, &
      & 1.06014784284417E-03_wp, 7.57152293456893E-03_wp, 1.87072473534496E-02_wp, &
      & 1.20922074604819E-32_wp,-2.86239545746042E-47_wp, 4.89204577231437E-17_wp, &
      &-4.20964955515460E-17_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 5.05683802648788E-32_wp, 1.00000000000000E+00_wp, 8.75870038752334E-32_wp, &
      &-7.35522649129989E-06_wp, 4.13143676691053E-05_wp, 1.51503426040404E-05_wp, &
      & 3.07334433234873E-04_wp,-8.45903480310013E-07_wp,-1.53891219610809E-07_wp, &
      & 4.03890556081742E-08_wp,-4.15735028643728E-07_wp, 1.69138892771752E-07_wp, &
      & 5.29463822990539E-04_wp,-9.72009902682532E-04_wp,-1.07588382966497E-03_wp, &
      &-2.39342273702004E-03_wp, 8.20652768913813E-06_wp, 4.52957647647255E-06_wp, &
      &-3.43919070258945E-06_wp, 1.21870355324475E-05_wp, 9.48860588650961E-06_wp, &
      &-2.15958353166773E-03_wp,-2.00334127727306E-02_wp,-3.30159001164187E-04_wp, &
      & 1.29217553049478E-02_wp,-4.39114994972613E-02_wp, 7.17137135298872E-03_wp, &
      &-1.08418518353207E-02_wp,-5.98154810786167E-02_wp,-6.52290168122297E-02_wp, &
      & 2.84103623255625E-02_wp, 2.15747204908677E-02_wp, 3.44816543864968E-03_wp, &
      &-2.55255858706855E-03_wp, 9.60821199164087E-03_wp, 2.54497483960324E-02_wp, &
      &-1.26546050185148E-03_wp, 0.00000000000000E+00_wp,-2.65957726697841E-16_wp, &
      & 0.00000000000000E+00_wp, 4.89204577231436E-17_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 8.75870038752334E-32_wp, &
      & 1.00000000000000E+00_wp,-4.40472520469573E-04_wp, 1.92716499172006E-03_wp, &
      & 1.39294549545589E-03_wp,-4.47927913764923E-05_wp, 2.78024017274902E-07_wp, &
      &-9.70670982165104E-06_wp, 1.11067967383642E-06_wp, 1.69138892771751E-07_wp, &
      & 7.92357040657421E-06_wp, 3.83623557110607E-04_wp,-9.09558628989837E-04_wp, &
      &-1.19927959972736E-03_wp,-1.40812402016892E-03_wp, 6.23653211157604E-06_wp, &
      & 4.16421483137780E-06_wp,-1.07203332966111E-06_wp, 9.48860588650961E-06_wp, &
      & 5.67169465412100E-06_wp, 2.71520165109221E-01_wp, 2.97076450678070E-01_wp, &
      & 1.49377548505845E-01_wp,-2.03438027207878E-01_wp,-1.84778662615733E-01_wp, &
      & 1.35676618656563E-01_wp,-7.50329502303753E-02_wp,-9.29113821534032E-02_wp, &
      &-7.16460413032241E-02_wp, 1.03166611962926E-02_wp,-4.27747085919899E-02_wp, &
      &-8.98833973673464E-03_wp,-3.14589601176654E-02_wp, 6.15521728410922E-03_wp, &
      &-1.65854568657006E-02_wp, 4.39506017098200E-02_wp, 1.71162831672690E-03_wp, &
      & 6.04732062392753E-03_wp, 3.40465600967899E-03_wp,-8.96611914772649E-05_wp, &
      &-1.30642898219517E-05_wp, 4.96083223094062E-04_wp,-9.31671142549513E-05_wp, &
      &-7.35522649129989E-06_wp,-4.40472520469573E-04_wp, 1.00000000000000E+00_wp, &
      & 5.68342833001024E-16_wp, 9.65402139259676E-17_wp,-9.64820142606403E-17_wp, &
      & 0.00000000000000E+00_wp,-1.60745268934820E-32_wp, 0.00000000000000E+00_wp, &
      &-3.13738942933159E-34_wp, 0.00000000000000E+00_wp, 1.56288685030650E-03_wp, &
      & 3.71154089755697E-03_wp, 1.04851555431008E-05_wp,-5.22939971122128E-03_wp, &
      &-4.94668272806928E-04_wp, 9.91829668611550E-07_wp,-3.02545648489902E-04_wp, &
      &-1.39744486879612E-06_wp, 1.72938846732419E-04_wp,-1.16844510294259E-01_wp, &
      & 6.86437027853556E-02_wp,-5.45295352799983E-02_wp, 7.42640456540307E-02_wp, &
      & 1.20053023357136E-01_wp,-8.81508072307406E-02_wp, 4.49644282170747E-03_wp, &
      & 1.11282240180291E-01_wp,-6.20570469522165E-02_wp,-3.46732517765560E-02_wp, &
      & 3.79563333930383E-02_wp,-1.17908674490106E-02_wp, 1.93989802055032E-02_wp, &
      & 4.97717098610996E-03_wp,-8.36777003079139E-03_wp,-8.68941884558286E-03_wp, &
      &-6.05377992113713E-03_wp,-1.47041207437914E-02_wp,-9.75933422332976E-03_wp, &
      & 2.57010849848263E-04_wp, 6.52724778871420E-05_wp,-2.47855655767690E-03_wp, &
      & 2.07525566721164E-04_wp, 4.13143676691053E-05_wp, 1.92716499172006E-03_wp, &
      & 5.68342833001024E-16_wp, 1.00000000000000E+00_wp, 2.60081076457839E-32_wp, &
      &-7.62088734896664E-32_wp, 2.51922571253414E-18_wp, 7.06006807073242E-17_wp, &
      &-1.51399551488900E-17_wp,-2.85936298470185E-47_wp,-2.62231715421916E-17_wp, &
      &-3.71492005494146E-03_wp,-4.72820131138755E-03_wp,-2.02072092878370E-05_wp, &
      & 1.00782076126601E-02_wp, 1.41301030138934E-03_wp,-2.83314216013732E-06_wp, &
      & 9.57447578053522E-04_wp, 5.29666419068912E-06_wp,-9.83316616938061E-04_wp, &
      &-5.87523732166724E-02_wp,-5.45295352799983E-02_wp, 1.49671163139505E-01_wp, &
      & 3.73418393029971E-02_wp, 1.11282240180291E-01_wp, 6.61585443518252E-02_wp, &
      & 1.31043068444519E-01_wp,-4.53053876036286E-02_wp, 4.31485533199955E-02_wp, &
      &-7.21339766751295E-03_wp,-1.28973512282152E-02_wp, 4.60414462317829E-02_wp, &
      &-1.32426681021874E-03_wp,-3.15291938146309E-02_wp,-5.00081635897833E-03_wp, &
      &-3.07302102712587E-02_wp,-3.40829260949412E-03_wp,-9.75933422332976E-03_wp, &
      &-2.86420408011397E-03_wp, 1.44697724646240E-04_wp, 4.13143676691053E-05_wp, &
      &-1.02183539811886E-03_wp, 6.50217215235351E-04_wp, 1.51503426040404E-05_wp, &
      & 1.39294549545589E-03_wp, 9.65402139259676E-17_wp, 2.60081076457839E-32_wp, &
      & 1.00000000000000E+00_wp,-1.02551007975096E-32_wp,-2.85936298470185E-47_wp, &
      & 2.62231715421916E-17_wp, 8.15226440226889E-17_wp, 2.51922571253414E-18_wp, &
      & 0.00000000000000E+00_wp,-1.04947017105172E-05_wp,-2.02072092878372E-05_wp, &
      & 2.42470022040784E-03_wp, 2.84710790830729E-05_wp, 5.29666419068914E-06_wp, &
      & 3.27824684470463E-04_wp, 4.30892340221626E-06_wp,-4.61890723453825E-04_wp, &
      &-1.85174398080753E-06_wp, 8.00151496696516E-02_wp, 7.42640456540307E-02_wp, &
      & 3.73418393029971E-02_wp, 1.26233954895127E-01_wp,-3.68657002091701E-03_wp, &
      & 1.11282240180291E-01_wp,-3.07916516106645E-03_wp,-1.85370059074957E-03_wp, &
      &-1.60025140456901E-01_wp,-4.26017054342602E-02_wp,-9.30137356639771E-03_wp, &
      & 4.97717098610996E-03_wp,-1.32844264627656E-02_wp,-7.93117290423926E-03_wp, &
      &-3.73170014599802E-02_wp, 2.15241086275338E-02_wp, 8.97569608799371E-05_wp, &
      & 2.57010849848264E-04_wp, 1.44697724646240E-04_wp, 2.62651388901871E-03_wp, &
      & 5.45884768170620E-04_wp, 4.13143676691057E-05_wp,-3.07689813908525E-06_wp, &
      & 3.07334433234873E-04_wp,-4.47927913764923E-05_wp,-9.64820142606403E-17_wp, &
      &-7.62088734896664E-32_wp,-1.02551007975096E-32_wp, 1.00000000000000E+00_wp, &
      & 2.62231715421916E-17_wp,-2.85936298470185E-47_wp,-1.45447564328101E-18_wp, &
      & 7.06006807073242E-17_wp, 2.51922571253414E-18_wp, 5.23416079701783E-03_wp, &
      & 1.00782076126601E-02_wp, 2.84710790830726E-05_wp,-1.17749989229400E-02_wp, &
      &-2.31383978324007E-03_wp, 5.29666419068912E-06_wp,-1.34900199846330E-03_wp, &
      &-6.53663012175235E-06_wp, 4.61638976007214E-04_wp,-6.00032867731247E-03_wp, &
      & 1.53233105239048E-02_wp, 7.28954899753139E-04_wp,-2.12520947060716E-02_wp, &
      & 3.93059495632419E-02_wp,-6.78524889755777E-02_wp,-1.15115579619026E-02_wp, &
      & 2.46537787898556E-02_wp, 5.06831870845617E-02_wp, 6.79736167259399E-03_wp, &
      & 2.23784051718079E-02_wp,-1.07989162672097E-02_wp, 1.14264510979431E-02_wp, &
      & 6.89125459326437E-03_wp, 2.07659262335779E-02_wp,-2.62400842432184E-02_wp, &
      &-1.30502848732678E-05_wp,-6.51497161910206E-05_wp,-4.12364049545137E-05_wp, &
      &-5.44826901506361E-04_wp,-1.44197580318789E-06_wp,-3.00514505778022E-07_wp, &
      & 3.29424435896268E-08_wp,-8.45903480310013E-07_wp, 2.78024017274902E-07_wp, &
      & 0.00000000000000E+00_wp, 2.51922571253414E-18_wp,-2.85936298470185E-47_wp, &
      & 2.62231715421916E-17_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 8.75221073988281E-32_wp, 0.00000000000000E+00_wp, &
      &-4.94213978435561E-04_wp,-1.41074225126450E-03_wp,-5.28809507065689E-06_wp, &
      & 2.31010912667401E-03_wp, 9.61420010084685E-06_wp,-2.01878790936736E-08_wp, &
      & 5.95206283931439E-06_wp, 2.98328163074629E-08_wp,-3.86639863179754E-06_wp, &
      & 4.40583503658527E-03_wp,-1.12513800515564E-02_wp,-2.07945739169201E-02_wp, &
      & 7.28954899753138E-04_wp,-6.78524889755777E-02_wp,-3.28094341533597E-03_wp, &
      &-6.26875333115744E-02_wp, 7.67644703164059E-03_wp, 3.85781955243902E-03_wp, &
      & 1.80558583064038E-02_wp,-4.62094080745720E-03_wp,-2.18241822493265E-02_wp, &
      &-1.07140661502377E-02_wp, 2.27818699288189E-02_wp,-3.07044231541223E-03_wp, &
      & 1.82997764559348E-02_wp, 4.95551420739942E-04_wp, 2.47389499407760E-03_wp, &
      & 1.01993488185775E-03_wp,-4.12364049545137E-05_wp,-3.00514505778023E-07_wp, &
      & 9.96138524408819E-06_wp,-2.72409461592087E-06_wp,-1.53891219610809E-07_wp, &
      &-9.70670982165104E-06_wp,-1.60745268934820E-32_wp, 7.06006807073242E-17_wp, &
      & 2.62231715421916E-17_wp,-2.85936298470185E-47_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 5.05309122667567E-32_wp, 0.00000000000000E+00_wp, &
      &-8.75221073988281E-32_wp, 9.90918790229741E-07_wp, 2.82859462893849E-06_wp, &
      &-3.27281568130850E-04_wp,-5.28809507065713E-06_wp,-2.01878790936732E-08_wp, &
      &-4.54326231998641E-07_wp,-1.43761384587137E-08_wp, 7.03114203587377E-07_wp, &
      & 9.16218156987502E-09_wp,-2.43654952707911E-03_wp, 1.22854141999639E-02_wp, &
      &-1.14668148754686E-02_wp,-8.41305469540806E-03_wp,-1.15115579619026E-02_wp, &
      &-6.26875333115744E-02_wp,-1.73980857687373E-02_wp, 4.29284383811855E-02_wp, &
      &-4.46348915793430E-03_wp, 2.29483867913033E-03_wp, 2.17475896955325E-02_wp, &
      &-1.46101566731171E-02_wp, 2.01491235766777E-02_wp, 1.00050389184092E-02_wp, &
      & 8.43240603387193E-03_wp, 9.77637424121205E-03_wp,-9.30672388945664E-05_wp, &
      &-2.07149974319160E-04_wp,-6.48972176334236E-04_wp, 3.07132938155310E-06_wp, &
      & 3.29424435896269E-08_wp,-2.72409461592086E-06_wp,-6.61949337710925E-07_wp, &
      & 4.03890556081742E-08_wp, 1.11067967383642E-06_wp, 0.00000000000000E+00_wp, &
      &-1.51399551488900E-17_wp, 8.15226440226889E-17_wp,-1.45447564328101E-18_wp, &
      & 0.00000000000000E+00_wp, 5.05309122667567E-32_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp,-3.02267796052731E-04_wp, &
      &-9.55905948998607E-04_wp,-4.30191068905611E-06_wp, 1.34682991016976E-03_wp, &
      & 5.95206283931439E-06_wp,-1.43761384587139E-08_wp, 3.68506189455136E-06_wp, &
      & 2.02553538757878E-08_wp,-2.08087508274799E-06_wp,-3.01711692731095E-03_wp, &
      & 7.28954899753139E-04_wp, 1.42401293829613E-02_wp,-1.06860904005050E-02_wp, &
      & 2.46537787898557E-02_wp, 7.67644703164060E-03_wp, 4.29284383811856E-02_wp, &
      & 2.67198462630366E-03_wp, 5.36114109247964E-02_wp, 1.56108390830926E-02_wp, &
      & 4.08897524068078E-03_wp, 2.27818699288189E-02_wp, 7.33699516068696E-03_wp, &
      &-4.15730696428098E-03_wp,-2.81830113107483E-04_wp,-1.90984663399790E-02_wp, &
      &-7.34734167161412E-06_wp,-4.12364049545137E-05_wp,-1.51221644135675E-05_wp, &
      &-3.06738851766667E-04_wp,-8.45903480310013E-07_wp,-1.53891219610808E-07_wp, &
      & 4.03890556081740E-08_wp,-4.15735028643728E-07_wp, 1.69138892771751E-07_wp, &
      &-3.13738942933159E-34_wp,-2.85936298470185E-47_wp, 2.51922571253414E-18_wp, &
      & 7.06006807073242E-17_wp, 8.75221073988281E-32_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.39616148076989E-06_wp,-5.28809507065713E-06_wp, 4.61125496151224E-04_wp, &
      & 6.52609096417557E-06_wp, 2.98328163074624E-08_wp, 7.03114203587377E-07_wp, &
      & 2.02553538757876E-08_wp,-9.45951695663841E-07_wp,-8.93614629390930E-09_wp, &
      &-2.32656622881648E-03_wp, 2.08214407972141E-02_wp, 2.82644825525715E-04_wp, &
      & 1.34886570383591E-02_wp, 5.06831870845617E-02_wp, 3.85781955243907E-03_wp, &
      &-4.46348915793432E-03_wp, 5.36114109247963E-02_wp,-7.17565927775726E-02_wp, &
      &-3.37969956237433E-02_wp, 2.07659262335779E-02_wp,-3.49031703133335E-03_wp, &
      & 4.43048983970076E-03_wp, 3.68961344026113E-03_wp,-2.31261633888675E-02_wp, &
      &-1.62238688233457E-03_wp,-4.40000332916343E-04_wp,-1.92355606367149E-03_wp, &
      &-1.39031692292213E-03_wp, 4.47078532418239E-05_wp, 2.78024017274903E-07_wp, &
      &-9.70670982165104E-06_wp, 1.11067967383642E-06_wp, 1.69138892771752E-07_wp, &
      & 7.92357040657421E-06_wp, 0.00000000000000E+00_wp,-2.62231715421916E-17_wp, &
      & 0.00000000000000E+00_wp, 2.51922571253414E-18_wp, 0.00000000000000E+00_wp, &
      &-8.75221073988281E-32_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 1.72780022831676E-04_wp, 9.81713023881929E-04_wp, &
      & 1.84874816761846E-06_wp,-4.60910086151720E-04_wp,-3.86639863179754E-06_wp, &
      & 9.16218156987517E-09_wp,-2.08087508274799E-06_wp,-8.93614629390945E-09_wp, &
      &-9.33973133689439E-08_wp, 2.71960093496385E-01_wp,-9.40524924437796E-02_wp, &
      & 1.48469584634808E-01_wp, 3.47930933401967E-01_wp,-1.00336855750964E-01_wp, &
      &-4.28158863348252E-02_wp,-7.59576632463371E-02_wp, 1.58389968302148E-01_wp, &
      & 1.72027923471804E-01_wp,-3.23280067060678E-02_wp,-2.30896339987546E-02_wp, &
      & 2.94912936553633E-03_wp,-3.14829673682181E-02_wp,-1.09097941608253E-02_wp, &
      & 3.95872659233849E-02_wp, 3.19021256333458E-02_wp, 1.90418230344655E-03_wp, &
      & 2.21107213999637E-03_wp, 3.76516222024187E-03_wp, 6.23961997945212E-03_wp, &
      & 3.10924905667185E-04_wp, 1.87620834602998E-04_wp,-1.00637672461044E-04_wp, &
      & 5.29463822990539E-04_wp, 3.83623557110607E-04_wp, 1.56288685030650E-03_wp, &
      &-3.71492005494146E-03_wp,-1.04947017105172E-05_wp, 5.23416079701783E-03_wp, &
      &-4.94213978435561E-04_wp, 9.90918790229741E-07_wp,-3.02267796052731E-04_wp, &
      &-1.39616148076989E-06_wp, 1.72780022831676E-04_wp, 1.00000000000000E+00_wp, &
      &-5.00893794320038E-17_wp, 2.86116926205749E-17_wp, 2.97256731543367E-16_wp, &
      &-2.10264728873578E-33_wp,-8.44195070578347E-33_wp, 0.00000000000000E+00_wp, &
      & 2.32374763135264E-32_wp, 0.00000000000000E+00_wp, 3.69310700214585E-02_wp, &
      & 1.66395373175104E-01_wp, 1.71049702610842E-02_wp, 4.00846293427067E-02_wp, &
      & 1.35742399194017E-01_wp, 5.79241903805787E-02_wp,-1.66254425274986E-03_wp, &
      & 6.00159024456601E-02_wp, 1.12154537526214E-01_wp, 3.46135196140600E-02_wp, &
      & 2.65554267780228E-02_wp,-5.38016899490450E-03_wp,-6.16924951477627E-03_wp, &
      & 2.55351514429918E-03_wp, 2.73507326766188E-02_wp, 4.25304632316621E-02_wp, &
      &-2.21142038639408E-03_wp, 8.11755796279226E-04_wp,-3.54265554044613E-03_wp, &
      &-5.87088231461763E-03_wp, 8.51213619928546E-06_wp, 5.13646243467491E-06_wp, &
      & 6.62314833876290E-05_wp,-9.72009902682532E-04_wp,-9.09558628989837E-04_wp, &
      & 3.71154089755697E-03_wp,-4.72820131138755E-03_wp,-2.02072092878372E-05_wp, &
      & 1.00782076126601E-02_wp,-1.41074225126450E-03_wp, 2.82859462893849E-06_wp, &
      &-9.55905948998607E-04_wp,-5.28809507065713E-06_wp, 9.81713023881929E-04_wp, &
      &-5.00893794320038E-17_wp, 1.00000000000000E+00_wp,-1.95703371397740E-33_wp, &
      &-1.56562697118192E-32_wp, 1.90615473333556E-16_wp,-2.06110993504776E-17_wp, &
      & 6.07532093347858E-18_wp,-8.55771842455788E-49_wp, 1.05227645290717E-17_wp, &
      &-5.82987274843623E-02_wp, 1.71049702610842E-02_wp, 1.50229426898005E-01_wp, &
      &-6.32768799009697E-02_wp, 6.00159024456601E-02_wp,-2.13609749530813E-02_wp, &
      & 1.31051968499253E-01_wp, 7.90212333633180E-02_wp,-1.02897494601968E-01_wp, &
      & 2.24627183919275E-02_wp,-7.21700332353314E-03_wp,-1.46043900971410E-02_wp, &
      &-1.74920521618782E-03_wp, 5.40264159538580E-02_wp, 1.23735798390771E-02_wp, &
      &-2.21668001594183E-02_wp,-3.76575523760944E-03_wp,-3.54265554044613E-03_wp, &
      &-3.14050988116331E-03_wp,-9.99733291855449E-03_wp,-9.72009902682531E-04_wp, &
      &-3.81250263554302E-04_wp, 7.18270198948683E-04_wp,-1.07588382966497E-03_wp, &
      &-1.19927959972736E-03_wp, 1.04851555431008E-05_wp,-2.02072092878370E-05_wp, &
      & 2.42470022040784E-03_wp, 2.84710790830726E-05_wp,-5.28809507065689E-06_wp, &
      &-3.27281568130850E-04_wp,-4.30191068905611E-06_wp, 4.61125496151224E-04_wp, &
      & 1.84874816761846E-06_wp, 2.86116926205749E-17_wp,-1.95703371397740E-33_wp, &
      & 1.00000000000000E+00_wp, 2.35968671848661E-32_wp,-8.55771842455788E-49_wp, &
      &-1.05227645290717E-17_wp,-2.37996475165848E-17_wp, 1.90615473333556E-16_wp, &
      & 0.00000000000000E+00_wp,-1.36620107880504E-01_wp, 4.00846293427067E-02_wp, &
      &-6.32768799009697E-02_wp, 2.89448721134034E-02_wp, 9.36731736887267E-02_wp, &
      & 6.00159024456601E-02_wp, 6.15029499645749E-03_wp,-1.47870798823357E-01_wp, &
      &-6.73737892461001E-02_wp, 2.41519984031383E-02_wp, 2.76714656194347E-02_wp, &
      & 2.55351514429919E-03_wp, 2.28220718695871E-02_wp,-1.41361919844603E-02_wp, &
      &-2.77419228427968E-02_wp,-3.10357127665506E-03_wp,-6.24060272675448E-03_wp, &
      &-5.87088231461763E-03_wp,-9.99733291855449E-03_wp,-1.36753998857564E-02_wp, &
      &-1.40552518683217E-03_wp,-9.72009902682532E-04_wp, 1.86904479297035E-04_wp, &
      &-2.39342273702004E-03_wp,-1.40812402016892E-03_wp,-5.22939971122128E-03_wp, &
      & 1.00782076126601E-02_wp, 2.84710790830729E-05_wp,-1.17749989229400E-02_wp, &
      & 2.31010912667401E-03_wp,-5.28809507065713E-06_wp, 1.34682991016976E-03_wp, &
      & 6.52609096417557E-06_wp,-4.60910086151720E-04_wp, 2.97256731543367E-16_wp, &
      &-1.56562697118192E-32_wp, 2.35968671848661E-32_wp, 1.00000000000000E+00_wp, &
      &-1.05227645290717E-17_wp,-8.55771842455788E-49_wp,-1.10051894840836E-16_wp, &
      &-2.06110993504776E-17_wp, 1.90615473333556E-16_wp,-3.15805596256058E-03_wp, &
      &-2.39637014091444E-02_wp, 4.36203320089823E-04_wp, 7.42537820441489E-03_wp, &
      &-5.36255231615542E-02_wp,-3.16478193052413E-02_wp,-6.03726229855206E-03_wp, &
      &-4.77276088174359E-02_wp,-6.59477421558427E-02_wp,-1.20517061069994E-02_wp, &
      &-1.55601171838614E-02_wp, 9.38158733271537E-04_wp, 6.25251359789685E-03_wp, &
      &-6.80798170999700E-03_wp,-2.68986291447388E-02_wp,-3.08724369198555E-02_wp, &
      & 3.10876048741532E-04_wp,-8.50486442997838E-06_wp, 9.71741044910373E-04_wp, &
      & 1.40513816130653E-03_wp, 3.03164987712131E-06_wp, 1.91875419651281E-06_wp, &
      &-8.68877459996168E-07_wp, 8.20652768913813E-06_wp, 6.23653211157604E-06_wp, &
      &-4.94668272806928E-04_wp, 1.41301030138934E-03_wp, 5.29666419068914E-06_wp, &
      &-2.31383978324007E-03_wp, 9.61420010084685E-06_wp,-2.01878790936732E-08_wp, &
      & 5.95206283931439E-06_wp, 2.98328163074624E-08_wp,-3.86639863179754E-06_wp, &
      &-2.10264728873578E-33_wp, 1.90615473333556E-16_wp,-8.55771842455788E-49_wp, &
      &-1.05227645290717E-17_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      &-1.34761014903253E-03_wp,-1.02258248777719E-02_wp, 6.58929517447636E-03_wp, &
      & 4.36203320089827E-04_wp,-3.16478193052413E-02_wp, 7.03473190075172E-03_wp, &
      & 1.99265001040424E-02_wp, 4.53565142340133E-03_wp,-4.11332714987785E-02_wp, &
      &-1.11470873876836E-02_wp, 2.75700778549133E-03_wp, 8.15193580412787E-03_wp, &
      & 3.31254163982971E-03_wp, 1.24342367450678E-02_wp,-2.55315572948441E-04_wp, &
      &-2.12187080366463E-02_wp, 1.87591352959682E-04_wp,-5.13207444451533E-06_wp, &
      & 3.81146552296546E-04_wp, 9.71741044910373E-04_wp, 1.91875419651281E-06_wp, &
      & 1.00972417245683E-06_wp,-1.21871184009155E-06_wp, 4.52957647647255E-06_wp, &
      & 4.16421483137780E-06_wp, 9.91829668611550E-07_wp,-2.83314216013732E-06_wp, &
      & 3.27824684470463E-04_wp, 5.29666419068912E-06_wp,-2.01878790936736E-08_wp, &
      &-4.54326231998641E-07_wp,-1.43761384587139E-08_wp, 7.03114203587377E-07_wp, &
      & 9.16218156987517E-09_wp,-8.44195070578347E-33_wp,-2.06110993504776E-17_wp, &
      &-1.05227645290717E-17_wp,-8.55771842455788E-49_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp,-2.39073219428602E-03_wp,-3.90605105295778E-03_wp, &
      &-1.13413923700238E-02_wp, 1.44497604843774E-02_wp,-6.03726229855204E-03_wp, &
      & 1.99265001040424E-02_wp,-1.68184321679367E-02_wp,-7.37146416909488E-02_wp, &
      & 1.03509093333786E-02_wp,-7.06986288266467E-03_wp, 1.18768916878722E-02_wp, &
      & 4.54596058457116E-03_wp, 2.04280163304608E-02_wp,-1.68169951513391E-02_wp, &
      &-2.03629762869517E-02_wp, 6.97672627774367E-03_wp,-1.00621858844341E-04_wp, &
      &-6.62141698033521E-05_wp,-7.18068099254989E-04_wp,-1.86855620562661E-04_wp, &
      &-8.68877459996167E-07_wp,-1.21871184009155E-06_wp,-9.96852139011951E-07_wp, &
      &-3.43919070258945E-06_wp,-1.07203332966111E-06_wp,-3.02545648489902E-04_wp, &
      & 9.57447578053522E-04_wp, 4.30892340221626E-06_wp,-1.34900199846330E-03_wp, &
      & 5.95206283931439E-06_wp,-1.43761384587137E-08_wp, 3.68506189455136E-06_wp, &
      & 2.02553538757876E-08_wp,-2.08087508274799E-06_wp, 0.00000000000000E+00_wp, &
      & 6.07532093347858E-18_wp,-2.37996475165848E-17_wp,-1.10051894840836E-16_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 4.98525073426546E-03_wp, &
      & 4.36203320089825E-04_wp,-2.43759581585471E-02_wp,-1.17215694036479E-02_wp, &
      &-4.77276088174359E-02_wp, 4.53565142340132E-03_wp,-7.37146416909488E-02_wp, &
      &-8.51805157287561E-03_wp, 5.60424300043609E-02_wp,-2.27065694199979E-02_wp, &
      &-2.47390282508854E-03_wp, 1.24342367450678E-02_wp,-1.22541750328176E-02_wp, &
      &-3.44851572401250E-02_wp, 5.45027214313180E-03_wp, 1.36980314278988E-02_wp, &
      & 5.29380626134436E-04_wp, 9.71741044910374E-04_wp, 1.07559115769634E-03_wp, &
      & 2.39276368394753E-03_wp, 8.20652768913813E-06_wp, 4.52957647647255E-06_wp, &
      &-3.43919070258945E-06_wp, 1.21870355324475E-05_wp, 9.48860588650961E-06_wp, &
      &-1.39744486879612E-06_wp, 5.29666419068912E-06_wp,-4.61890723453825E-04_wp, &
      &-6.53663012175235E-06_wp, 2.98328163074629E-08_wp, 7.03114203587377E-07_wp, &
      & 2.02553538757878E-08_wp,-9.45951695663841E-07_wp,-8.93614629390945E-09_wp, &
      & 2.32374763135264E-32_wp,-8.55771842455788E-49_wp, 1.90615473333556E-16_wp, &
      &-2.06110993504776E-17_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 1.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 5.41449904305805E-03_wp,-5.92939585074226E-03_wp,-7.47872262938034E-04_wp, &
      &-2.54399760655455E-02_wp,-6.59477421558427E-02_wp,-4.11332714987785E-02_wp, &
      & 1.03509093333787E-02_wp, 5.60424300043609E-02_wp, 2.09775110980944E-02_wp, &
      &-2.84653330969374E-02_wp,-2.68986291447388E-02_wp,-3.54028932533078E-03_wp, &
      &-1.07199585105105E-02_wp, 1.11501042348473E-02_wp, 1.48688146646377E-02_wp, &
      & 3.14704172073689E-03_wp, 3.83563276743010E-04_wp, 9.09305302729587E-04_wp, &
      & 1.19894787919603E-03_wp, 1.40773945062066E-03_wp, 6.23653211157604E-06_wp, &
      & 4.16421483137780E-06_wp,-1.07203332966111E-06_wp, 9.48860588650961E-06_wp, &
      & 5.67169465412100E-06_wp, 1.72938846732419E-04_wp,-9.83316616938061E-04_wp, &
      &-1.85174398080753E-06_wp, 4.61638976007214E-04_wp,-3.86639863179754E-06_wp, &
      & 9.16218156987502E-09_wp,-2.08087508274799E-06_wp,-8.93614629390930E-09_wp, &
      &-9.33973133689439E-08_wp, 0.00000000000000E+00_wp, 1.05227645290717E-17_wp, &
      & 0.00000000000000E+00_wp, 1.90615473333556E-16_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 1.00000000000000E+00_wp],&
      & shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "CeCl3")
   call test_overlap_mol(error, mol, make_qvszp_basis, .true., overlap)

end subroutine test_overlap_qvszp_cecl3

subroutine test_overlap_diat_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 43
   real(wp), parameter :: overlap(nao, nao) = reshape([&
   &  1.00000000017687E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,&
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  8.69412276329442E-02_wp,  1.09663227357883E-01_wp, &
   &  6.85934472373425E-02_wp,  7.46686270912962E-02_wp,  5.39873838622696E-02_wp, &
   &  4.95948688317373E-02_wp, -1.55903398924374E-02_wp,  3.37686648082453E-02_wp, &
   & -2.12649579080578E-02_wp,  8.67139986554086E-02_wp, -1.13584328085048E-01_wp, &
   & -5.71130711953149E-02_wp,  7.77825760830674E-02_wp, -5.83831463741125E-02_wp, &
   &  4.28687369766487E-02_wp, -2.37076059225061E-02_wp, -2.93565217287235E-02_wp, &
   & -2.26374693772444E-02_wp,  8.73543363380361E-02_wp,  3.61731541134505E-02_wp, &
   & -5.71022949696438E-02_wp, -1.33816329028274E-01_wp, -3.17825201816373E-02_wp, &
   & -1.35622823871297E-02_wp, -2.40602114448212E-02_wp,  5.01713186690479E-02_wp, &
   &  5.44912525773983E-02_wp,  0.00000000000000E+00_wp,  1.00000000002709E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  2.77555756156289E-17_wp,  0.00000000000000E+00_wp, &
   &  5.55111512312578E-17_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -3.40201809683547E-01_wp, &
   & -4.14280857819152E-02_wp, -1.19418786053541E-01_wp, -1.29995461121435E-01_wp, &
   & -4.72851062125490E-03_wp, -4.34379010897963E-03_wp, -1.71048248105265E-02_wp, &
   & -2.36191595818926E-02_wp, -3.36400249951174E-02_wp,  3.52869972519298E-01_wp, &
   & -5.63511977455039E-02_wp, -1.03417073813248E-01_wp,  1.40844227841612E-01_wp, &
   & -8.05377110703101E-03_wp,  5.91360720857850E-03_wp,  1.17908981123966E-02_wp, &
   & -2.13786714107248E-02_wp,  3.38405303098483E-02_wp, -1.11927689146485E-01_wp, &
   &  1.29184921393787E-01_wp,  3.25455034102347E-02_wp,  7.62687348214290E-02_wp, &
   & -5.17690901910778E-02_wp, -2.20909800759548E-02_wp, -3.71443619535530E-03_wp, &
   & -1.14661853993300E-02_wp, -2.84111485232627E-02_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  1.00000000002709E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.55111512312578E-17_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   & -2.12793435363904E-01_wp, -1.19418786053541E-01_wp,  7.47962781286959E-02_wp, &
   & -8.13110923174573E-02_wp, -2.36191595818926E-02_wp,  2.68161297325115E-02_wp, &
   &  4.18598977604196E-02_wp,  1.82588424513004E-02_wp,  9.30329270286544E-03_wp, &
   &  1.77431941562336E-01_wp, -1.03417073813248E-01_wp,  9.73200705922398E-02_wp, &
   &  7.08200378325435E-02_wp, -2.13786714107248E-02_wp, -3.46284778301005E-02_wp, &
   & -3.79011924083843E-02_wp,  2.37135902185714E-02_wp, -8.28936173267027E-03_wp, &
   &  1.76687050868385E-01_wp,  3.25455034102347E-02_wp,  9.84260921377083E-02_wp, &
   & -1.20396462500221E-01_wp, -1.14661853993300E-02_wp,  1.10648156497081E-02_wp, &
   & -3.77677054156251E-02_wp, -4.09323723050175E-02_wp,  1.96588187822549E-02_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  1.00000000002709E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  5.55111512312578E-17_wp,  0.00000000000000E+00_wp, &
   & -2.77555756156289E-17_wp, -2.31640110135961E-01_wp, -1.29995461121435E-01_wp, &
   & -8.13110923174573E-02_wp,  6.09790994241015E-02_wp,  2.28025303637907E-02_wp, &
   & -2.36191595818926E-02_wp, -1.16465091900958E-02_wp,  1.42627952967513E-02_wp, &
   &  4.31596987206192E-02_wp, -2.41645444821936E-01_wp,  1.40844227841612E-01_wp, &
   &  7.08200378325435E-02_wp,  5.28706271385305E-02_wp, -2.12103821177911E-02_wp, &
   & -2.13786714107248E-02_wp, -8.07440995582121E-03_wp, -1.06651162567620E-02_wp, &
   &  4.57526496102248E-02_wp,  4.14057132845683E-01_wp,  7.62687348214290E-02_wp, &
   & -1.20396462500221E-01_wp, -1.32341154522262E-01_wp, -1.09127406617264E-02_wp, &
   & -1.14661853993300E-02_wp,  1.37409144503486E-02_wp,  1.72266574885550E-02_wp, &
   & -1.29632440844008E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  9.99999999778610E-01_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  8.20155010090769E-02_wp, &
   &  9.28979457969490E-02_wp,  8.00569085802895E-02_wp,  3.56083920092022E-02_wp, &
   &  3.58620481214061E-02_wp,  6.49720464580809E-02_wp,  9.89663399650733E-03_wp, &
   &  2.86445384905872E-02_wp, -3.63460078236818E-02_wp, -8.86514552275003E-02_wp, &
   &  1.07106300572797E-01_wp,  7.22242203094390E-02_wp, -4.50177541095921E-02_wp, &
   &  4.63175554736653E-02_wp, -6.21718876923976E-02_wp,  4.41008058157470E-03_wp, &
   &  2.92431802609242E-02_wp,  4.19093396272880E-02_wp, -4.83242438420612E-02_wp, &
   &  3.82018892333236E-02_wp,  3.90978820690418E-02_wp,  7.46019922122533E-02_wp, &
   & -3.05116894672468E-02_wp, -1.58559058503109E-02_wp,  2.58807011549745E-03_wp, &
   & -4.23724796435138E-02_wp, -5.46710758908302E-02_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  9.99999999778610E-01_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  7.53425656833336E-02_wp,  8.53395943076514E-02_wp,  2.20043365124130E-02_wp, &
   &  8.00569085802895E-02_wp,  6.49720464580809E-02_wp,  2.48213604227393E-02_wp, &
   & -4.12640911091508E-02_wp,  2.13736344720395E-02_wp, -4.31605738168390E-03_wp, &
   &  6.50937154430233E-02_wp, -7.86444738413934E-02_wp,  3.13009539865031E-04_wp, &
   &  7.22242203094390E-02_wp, -6.21718876923976E-02_wp,  7.29603146759239E-03_wp, &
   & -4.67218548895194E-02_wp, -5.71517669530143E-03_wp, -5.66727024069538E-03_wp, &
   & -2.06209903237694E-02_wp,  1.63015647136606E-02_wp, -3.38005258954584E-04_wp, &
   &  3.90978820690418E-02_wp, -1.58559058503109E-02_wp, -1.20233934496212E-04_wp, &
   &  1.48991414842628E-02_wp, -2.81564588294365E-03_wp, -3.12937259411050E-02_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  9.99999999778609E-01_wp,  0.00000000000000E+00_wp,  4.80740671595891E-17_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp, -2.36842285309118E-02_wp, -7.20467680720084E-03_wp, &
   & -6.03429612836281E-02_wp, -4.90559450776114E-03_wp,  9.89663399650732E-03_wp, &
   & -4.12640911091508E-02_wp, -2.76120606165884E-02_wp, -2.80963191173944E-02_wp, &
   & -3.89816083520694E-03_wp, -3.59986382289629E-02_wp,  2.75278884074025E-02_wp, &
   &  6.03006889788791E-02_wp, -1.88511048183665E-02_wp,  4.41008058157471E-03_wp, &
   & -4.67218548895194E-02_wp, -1.18825152389580E-02_wp,  3.19951378324389E-02_wp, &
   &  1.70996375352673E-03_wp, -3.65827353559851E-02_wp, -8.92221837554665E-03_wp, &
   &  6.06254946164389E-02_wp,  3.30062041606793E-02_wp,  2.58807011549746E-03_wp, &
   &  1.48991414842628E-02_wp, -1.14945896141096E-02_wp, -5.51167977457504E-02_wp, &
   & -4.43725612524172E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  9.99999999778610E-01_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  5.13000216813875E-02_wp, &
   &  8.00569085802895E-02_wp,  1.49825391521149E-02_wp,  2.22727564867189E-02_wp, &
   &  2.86445384905872E-02_wp,  2.13736344720395E-02_wp, -2.80963191173944E-02_wp, &
   &  7.98374374734485E-03_wp, -4.25294966518297E-02_wp, -4.45761925024143E-02_wp, &
   &  7.22242203094390E-02_wp, -2.14349010640342E-04_wp, -2.26360646654454E-02_wp, &
   &  2.92431802609242E-02_wp, -5.71517669530144E-03_wp,  3.19951378324389E-02_wp, &
   &  2.86402855372428E-03_wp,  3.82652025183376E-02_wp,  7.62837881761689E-02_wp, &
   &  3.90978820690418E-02_wp,  1.25039201181341E-03_wp, -1.17765372388225E-01_wp, &
   & -4.23724796435138E-02_wp, -2.81564588294364E-03_wp, -5.51167977457504E-02_wp, &
   &  9.53463736577520E-03_wp,  5.68398395712986E-02_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  4.80740671595891E-17_wp, &
   &  0.00000000000000E+00_wp,  9.99999999778609E-01_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   & -3.23048840672788E-02_wp,  1.12518921460422E-03_wp, -3.15334189104662E-02_wp, &
   & -6.94186677160610E-02_wp, -3.63460078236819E-02_wp, -4.31605738168388E-03_wp, &
   & -3.89816083520693E-03_wp, -4.25294966518297E-02_wp, -4.20968008654698E-02_wp, &
   & -3.43736973355472E-02_wp,  2.34893769725024E-03_wp,  2.80042045879740E-02_wp, &
   & -7.46695782592206E-02_wp,  4.19093396272880E-02_wp, -5.66727024069537E-03_wp, &
   &  1.70996375352671E-03_wp,  3.82652025183376E-02_wp, -4.55187727956170E-02_wp, &
   &  8.28521011474394E-02_wp,  5.94862619285365E-02_wp, -6.70334685509428E-02_wp, &
   & -9.41199197927978E-02_wp, -5.46710758908302E-02_wp, -3.12937259411051E-02_wp, &
   & -4.43725612524171E-03_wp,  5.68398395712986E-02_wp,  3.13346757662612E-02_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  9.99999999954551E-01_wp,  0.00000000000000E+00_wp,  1.11022302462516E-16_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp, -6.37703812927541E-03_wp, -2.24936604198012E-02_wp, &
   & -7.22571944109632E-03_wp,  1.99127981578240E-02_wp,  5.81512727763616E-03_wp, &
   & -1.72393727662213E-02_wp, -2.14384122602978E-03_wp,  1.29355636537831E-02_wp, &
   &  2.79559437246025E-02_wp,  7.40406770098628E-03_wp, -2.55601123198589E-02_wp, &
   & -7.00152224711212E-03_wp, -2.04737585704892E-02_wp,  5.71402985740158E-03_wp, &
   &  1.64261267820631E-02_wp,  6.55758393614611E-04_wp,  1.09396897005001E-02_wp, &
   & -3.25759068310625E-02_wp, -2.32970808181820E-02_wp,  1.42279516174685E-02_wp, &
   &  2.18832031916388E-02_wp,  3.48981234262437E-02_wp, -1.02121852619804E-02_wp, &
   & -6.78822012440303E-03_wp, -1.93575430061701E-03_wp, -2.57527981119962E-02_wp, &
   & -2.95098567978740E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  9.99999999954551E-01_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -3.40431849326325E-02_wp, &
   & -4.74827302104553E-02_wp, -1.58927481291258E-02_wp, -2.11545204817318E-02_wp, &
   & -2.63941704687171E-02_wp, -2.11135912416699E-02_wp,  1.60699196580794E-02_wp, &
   & -1.07444245139649E-02_wp,  2.19190965349443E-02_wp, -3.06985789568106E-02_wp, &
   &  4.54124683874554E-02_wp,  4.52717213649874E-03_wp, -2.15442376617703E-02_wp, &
   &  2.77251308119991E-02_wp, -9.78968090562189E-03_wp,  2.20763678154631E-02_wp, &
   &  5.41731475723263E-03_wp,  2.10798636529092E-02_wp, -1.66394753076522E-02_wp, &
   &  1.12164550839581E-02_wp,  2.25211882161450E-03_wp,  3.09187955546773E-02_wp, &
   & -1.08132873910535E-02_wp, -8.10584518942799E-04_wp,  1.20817835115180E-02_wp, &
   & -6.24737443153229E-03_wp, -2.73626654703935E-02_wp,  0.00000000000000E+00_wp, &
   &  2.77555756156289E-17_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.11022302462516E-16_wp, &
   &  0.00000000000000E+00_wp,  9.99999999954550E-01_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   & -1.28090757033867E-03_wp, -1.23204005913030E-02_wp,  2.49038884936724E-02_wp, &
   & -8.75229503253804E-03_wp, -1.60649283629892E-02_wp,  1.43945681826134E-02_wp, &
   &  2.29720108902036E-02_wp,  1.37607056483274E-02_wp,  4.08008511942591E-03_wp, &
   & -6.45075714555441E-03_wp, -1.77211746261234E-03_wp,  2.88573809092836E-02_wp, &
   & -5.59343295542887E-04_wp, -6.89049708634549E-03_wp, -2.15630999197224E-02_wp, &
   & -1.47401248948007E-02_wp,  1.93252201033651E-02_wp, -3.00247894938941E-03_wp, &
   &  2.12528120885591E-03_wp, -2.57616878551625E-03_wp, -9.21701378652493E-03_wp, &
   & -4.47883417102252E-04_wp,  1.72123004852782E-03_wp,  3.87219849948905E-03_wp, &
   &  4.60884129645277E-03_wp,  1.05925682350414E-02_wp, -1.67051092576148E-03_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  9.99999999954551E-01_wp,  0.00000000000000E+00_wp,  5.55111512312578E-17_wp, &
   &  0.00000000000000E+00_wp,  2.32909205984756E-02_wp,  2.60524728274164E-02_wp, &
   &  1.76031258469667E-02_wp,  1.77388758768520E-02_wp,  1.34622945224013E-02_wp, &
   &  1.98904123574332E-02_wp, -9.39476671196776E-03_wp,  1.35431887132267E-02_wp, &
   & -5.30263009400632E-03_wp, -2.25774856886379E-02_wp,  2.85240889524626E-02_wp, &
   &  8.00111168452364E-03_wp, -1.95333031990455E-02_wp,  1.78063068515174E-02_wp, &
   & -1.19714716143595E-02_wp,  1.69654313272323E-02_wp,  8.19806673481372E-03_wp, &
   &  6.90421382037361E-03_wp, -2.26881057605045E-02_wp, -9.07844125129397E-03_wp, &
   &  7.72420256658115E-03_wp,  3.35841236773789E-02_wp,  9.71809772688050E-03_wp, &
   &  3.71972593593363E-03_wp,  1.72717431621451E-02_wp, -1.37604829309815E-02_wp, &
   & -1.66617157727234E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  2.77555756156289E-17_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  9.99999999954550E-01_wp, &
   &  0.00000000000000E+00_wp, -1.66533453693773E-16_wp, -8.72157531858014E-04_wp, &
   & -8.75229503253805E-03_wp,  1.69568159524313E-02_wp, -5.42556064348706E-03_wp, &
   & -9.29097453288089E-03_wp,  1.37607056483274E-02_wp,  1.56414192434809E-02_wp, &
   &  3.55422714503404E-03_wp,  6.96070362833565E-03_wp,  4.41748009542106E-03_wp, &
   & -5.59343295542888E-04_wp, -1.97615416138547E-02_wp, -2.20587677537982E-03_wp, &
   &  4.95776309493031E-03_wp,  1.93252201033651E-02_wp,  1.00940411888978E-02_wp, &
   & -6.57678604430529E-03_wp,  1.43932135734599E-03_wp, -7.86211035482032E-03_wp, &
   & -4.47883417102248E-04_wp,  3.40967488112168E-02_wp, -1.04037311751583E-03_wp, &
   & -4.15510671194854E-03_wp,  1.05925682350414E-02_wp, -1.70496114723901E-02_wp, &
   & -3.24497991168268E-02_wp,  7.47010245689619E-03_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  9.99999999954550E-01_wp,  0.00000000000000E+00_wp, &
   &  1.34091864220632E-02_wp,  3.45519482653732E-03_wp,  6.25995549018286E-03_wp, &
   &  3.07262139246520E-02_wp,  2.19190965349443E-02_wp,  3.36169334446591E-03_wp, &
   & -6.32974114816765E-03_wp,  1.15088725091117E-02_wp,  2.06202786830572E-02_wp, &
   & -1.19030608012565E-02_wp,  4.39408073360094E-03_wp,  1.75536480937154E-03_wp, &
   & -2.76499015081222E-02_wp,  2.10798636529093E-02_wp, -2.01631062617650E-03_wp, &
   &  8.55988639565553E-03_wp,  4.69913239680737E-03_wp, -1.84673717546551E-02_wp, &
   &  2.85284441435993E-02_wp,  2.26838862024821E-02_wp, -3.86126634519725E-03_wp, &
   & -4.16800783537248E-02_wp, -2.73626654703935E-02_wp, -3.96213612481711E-03_wp, &
   & -2.07142641033214E-02_wp,  9.26442849949040E-03_wp,  2.01405830947493E-02_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp, -1.11022302462516E-16_wp,  0.00000000000000E+00_wp, &
   &  9.99999999954551E-01_wp,  2.81785548509627E-02_wp,  2.32672185743682E-02_wp, &
   &  3.19286677421413E-02_wp,  2.38149248312371E-02_wp,  2.22115769413314E-02_wp, &
   &  2.28581711018306E-02_wp,  9.47310434010206E-03_wp,  2.11477642204028E-02_wp, &
   & -1.90384369106144E-03_wp,  3.15424946469630E-02_wp, -2.93108178909011E-02_wp, &
   & -2.98275876071079E-02_wp,  2.89865919149946E-02_wp, -2.83688082507035E-02_wp, &
   &  2.34835834975460E-02_wp,  2.79363404761069E-03_wp, -2.12898936281053E-02_wp, &
   & -2.83300770467368E-03_wp,  2.29901709038050E-02_wp,  3.00639987805416E-02_wp, &
   & -2.15949193473986E-02_wp, -2.25160768370592E-02_wp, -3.08515669011558E-02_wp, &
   & -1.89302474912947E-02_wp,  1.91025315773589E-03_wp,  1.84855092937346E-02_wp, &
   &  5.24876072680659E-03_wp,  8.69412276329442E-02_wp, -3.40201809683547E-01_wp, &
   & -2.12793435363904E-01_wp, -2.31640110135961E-01_wp,  8.20155010090769E-02_wp, &
   &  7.53425656833336E-02_wp, -2.36842285309118E-02_wp,  5.13000216813875E-02_wp, &
   & -3.23048840672788E-02_wp, -6.37703812927541E-03_wp, -3.40431849326325E-02_wp, &
   & -1.28090757033867E-03_wp,  2.32909205984756E-02_wp, -8.72157531858011E-04_wp, &
   &  1.34091864220632E-02_wp,  2.81785548509627E-02_wp,  9.99999999869332E-01_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  9.52137814488968E-05_wp, &
   & -3.36928188670201E-04_wp, -1.89691378004237E-04_wp,  4.99549878650670E-06_wp, &
   & -3.01084193947362E-06_wp,  1.14329075205538E-04_wp, -2.14716191083883E-05_wp, &
   & -1.69511123039559E-06_wp, -1.01512838117467E-04_wp,  1.12758480338272E-04_wp, &
   & -1.31040324396565E-04_wp, -2.23144270067549E-04_wp, -3.69794278272637E-04_wp, &
   &  7.55446878381371E-05_wp,  4.55857897636091E-05_wp, -2.44516968961120E-05_wp, &
   &  1.28642570924251E-04_wp,  9.32081069771156E-05_wp,  1.09663227357883E-01_wp, &
   & -4.14280857819152E-02_wp, -1.19418786053541E-01_wp, -1.29995461121435E-01_wp, &
   &  9.28979457969490E-02_wp,  8.53395943076514E-02_wp, -7.20467680720085E-03_wp, &
   &  8.00569085802895E-02_wp,  1.12518921460423E-03_wp, -2.24936604198012E-02_wp, &
   & -4.74827302104553E-02_wp, -1.23204005913030E-02_wp,  2.60524728274164E-02_wp, &
   & -8.75229503253804E-03_wp,  3.45519482653731E-03_wp,  2.32672185743682E-02_wp, &
   &  0.00000000000000E+00_wp,  9.99999999998060E-01_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  3.36928188670201E-04_wp, -9.41758669062715E-04_wp, -5.84225011720663E-04_wp, &
   &  1.53854928347468E-05_wp, -1.07377169272165E-05_wp,  4.07737526840130E-04_wp, &
   & -4.66396749143034E-05_wp, -6.57520964570179E-06_wp, -3.30286014739083E-04_wp, &
   &  1.31040324396565E-04_wp, -2.04150976546664E-05_wp, -2.26089711420990E-04_wp, &
   & -3.74675458323373E-04_wp,  2.52051845112027E-05_wp,  1.52095173725851E-05_wp, &
   & -1.65836799186663E-05_wp,  1.63404821497647E-04_wp,  1.43467330644811E-04_wp, &
   &  6.85934472373425E-02_wp, -1.19418786053541E-01_wp,  7.47962781286959E-02_wp, &
   & -8.13110923174573E-02_wp,  8.00569085802895E-02_wp,  2.20043365124130E-02_wp, &
   & -6.03429612836281E-02_wp,  1.49825391521149E-02_wp, -3.15334189104662E-02_wp, &
   & -7.22571944109633E-03_wp, -1.58927481291258E-02_wp,  2.49038884936724E-02_wp, &
   &  1.76031258469667E-02_wp,  1.69568159524313E-02_wp,  6.25995549018286E-03_wp, &
   &  3.19286677421413E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  9.99999999998060E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  1.89691378004237E-04_wp, -5.84225011720663E-04_wp, &
   & -2.32983327085689E-04_wp,  8.66206935257109E-06_wp, -6.57520964570179E-06_wp, &
   &  1.86202195665234E-04_wp, -8.81555136696621E-05_wp, -2.76074508981215E-06_wp, &
   & -2.21688220693919E-04_wp,  2.23144270067549E-04_wp, -2.26089711420990E-04_wp, &
   & -2.72645831427766E-04_wp, -6.38022547981310E-04_wp,  1.63404821497647E-04_wp, &
   &  7.35308666900574E-05_wp, -1.02189065578875E-04_wp,  2.07503254464807E-04_wp, &
   &  2.01611185625156E-04_wp,  7.46686270912962E-02_wp, -1.29995461121435E-01_wp, &
   & -8.13110923174573E-02_wp,  6.09790994241015E-02_wp,  3.56083920092022E-02_wp, &
   &  8.00569085802895E-02_wp, -4.90559450776115E-03_wp,  2.22727564867189E-02_wp, &
   & -6.94186677160610E-02_wp,  1.99127981578240E-02_wp, -2.11545204817318E-02_wp, &
   & -8.75229503253804E-03_wp,  1.77388758768520E-02_wp, -5.42556064348706E-03_wp, &
   &  3.07262139246520E-02_wp,  2.38149248312371E-02_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  9.99999999998060E-01_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -4.99549878650671E-06_wp, &
   &  1.53854928347468E-05_wp,  8.66206935257109E-06_wp,  9.57086406720181E-05_wp, &
   & -6.33015324737798E-05_wp, -6.57520964570179E-06_wp,  6.91507707790898E-07_wp, &
   & -3.56389145477078E-05_wp,  6.77924535368498E-06_wp,  3.69794278272637E-04_wp, &
   & -3.74675458323373E-04_wp, -6.38022547981310E-04_wp, -9.44974883405076E-04_wp, &
   &  2.45721973530458E-04_wp,  1.63404821497647E-04_wp, -4.67989527259473E-05_wp, &
   &  4.18431888622898E-04_wp,  2.63356233815959E-04_wp,  5.39873838622696E-02_wp, &
   & -4.72851062125489E-03_wp, -2.36191595818926E-02_wp,  2.28025303637907E-02_wp, &
   &  3.58620481214061E-02_wp,  6.49720464580809E-02_wp,  9.89663399650732E-03_wp, &
   &  2.86445384905872E-02_wp, -3.63460078236818E-02_wp,  5.81512727763616E-03_wp, &
   & -2.63941704687171E-02_wp, -1.60649283629892E-02_wp,  1.34622945224013E-02_wp, &
   & -9.29097453288089E-03_wp,  2.19190965349443E-02_wp,  2.22115769413314E-02_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  9.99999999830205E-01_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   & -3.01084193947362E-06_wp,  1.07377169272165E-05_wp,  6.57520964570179E-06_wp, &
   &  6.33015324737798E-05_wp, -3.73996719315490E-05_wp, -4.64599908927517E-06_wp, &
   &  2.45641225072744E-07_wp, -2.26199417185387E-05_wp,  4.42393841717257E-06_wp, &
   &  7.55446878381371E-05_wp, -2.52051845112026E-05_wp, -1.63404821497647E-04_wp, &
   & -2.45721973530458E-04_wp,  2.47089644023801E-05_wp,  1.61116888311556E-05_wp, &
   & -4.94256032779233E-06_wp,  1.09891385286233E-04_wp,  8.63251486754047E-05_wp, &
   &  4.95948688317373E-02_wp, -4.34379010897964E-03_wp,  2.68161297325115E-02_wp, &
   & -2.36191595818926E-02_wp,  6.49720464580809E-02_wp,  2.48213604227393E-02_wp, &
   & -4.12640911091509E-02_wp,  2.13736344720395E-02_wp, -4.31605738168389E-03_wp, &
   & -1.72393727662213E-02_wp, -2.11135912416699E-02_wp,  1.43945681826134E-02_wp, &
   &  1.98904123574332E-02_wp,  1.37607056483274E-02_wp,  3.36169334446591E-03_wp, &
   &  2.28581711018306E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  9.99999999830205E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  1.14329075205538E-04_wp, -4.07737526840130E-04_wp, &
   & -1.86202195665234E-04_wp,  6.57520964570179E-06_wp, -4.64599908927517E-06_wp, &
   &  1.38897990982117E-04_wp, -4.86344418783074E-05_wp, -2.20750442352880E-06_wp, &
   & -1.45294009909431E-04_wp,  4.55857897636091E-05_wp, -1.52095173725851E-05_wp, &
   & -7.35308666900573E-05_wp, -1.63404821497647E-04_wp,  1.61116888311556E-05_wp, &
   &  7.73094926501920E-06_wp, -1.90066643814781E-05_wp,  5.65625287713878E-05_wp, &
   &  6.13425902955612E-05_wp, -1.55903398924374E-02_wp, -1.71048248105265E-02_wp, &
   &  4.18598977604196E-02_wp, -1.16465091900958E-02_wp,  9.89663399650731E-03_wp, &
   & -4.12640911091508E-02_wp, -2.76120606165884E-02_wp, -2.80963191173944E-02_wp, &
   & -3.89816083520695E-03_wp, -2.14384122602977E-03_wp,  1.60699196580794E-02_wp, &
   &  2.29720108902036E-02_wp, -9.39476671196776E-03_wp,  1.56414192434809E-02_wp, &
   & -6.32974114816765E-03_wp,  9.47310434010206E-03_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  9.99999999830206E-01_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.14716191083883E-05_wp, &
   &  4.66396749143034E-05_wp,  8.81555136696621E-05_wp, -6.91507707790898E-07_wp, &
   &  2.45641225072744E-07_wp, -4.86344418783074E-05_wp, -2.10150359335817E-05_wp, &
   &  7.21083315540951E-07_wp,  8.28198172373845E-06_wp, -2.44516968961120E-05_wp, &
   &  1.65836799186663E-05_wp,  1.02189065578875E-04_wp,  4.67989527259473E-05_wp, &
   & -4.94256032779232E-06_wp, -1.90066643814781E-05_wp, -2.64669742429983E-05_wp, &
   & -5.36365868268805E-05_wp, -6.09820101131120E-06_wp,  3.37686648082453E-02_wp, &
   & -2.36191595818926E-02_wp,  1.82588424513004E-02_wp,  1.42627952967512E-02_wp, &
   &  2.86445384905872E-02_wp,  2.13736344720395E-02_wp, -2.80963191173944E-02_wp, &
   &  7.98374374734485E-03_wp, -4.25294966518297E-02_wp,  1.29355636537831E-02_wp, &
   & -1.07444245139649E-02_wp,  1.37607056483274E-02_wp,  1.35431887132267E-02_wp, &
   &  3.55422714503404E-03_wp,  1.15088725091117E-02_wp,  2.11477642204028E-02_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  9.99999999830205E-01_wp,  0.00000000000000E+00_wp, &
   & -1.69511123039559E-06_wp,  6.57520964570179E-06_wp,  2.76074508981216E-06_wp, &
   &  3.56389145477078E-05_wp, -2.26199417185387E-05_wp, -2.20750442352880E-06_wp, &
   &  7.21083315540952E-07_wp, -9.95740806845498E-06_wp,  2.82716000311223E-06_wp, &
   &  1.28642570924251E-04_wp, -1.63404821497647E-04_wp, -2.07503254464807E-04_wp, &
   & -4.18431888622897E-04_wp,  1.09891385286233E-04_wp,  5.65625287713878E-05_wp, &
   & -5.36365868268805E-05_wp,  1.47306249222992E-04_wp,  1.20892448538428E-04_wp, &
   & -2.12649579080578E-02_wp, -3.36400249951174E-02_wp,  9.30329270286542E-03_wp, &
   &  4.31596987206192E-02_wp, -3.63460078236819E-02_wp, -4.31605738168392E-03_wp, &
   & -3.89816083520695E-03_wp, -4.25294966518297E-02_wp, -4.20968008654698E-02_wp, &
   &  2.79559437246025E-02_wp,  2.19190965349444E-02_wp,  4.08008511942590E-03_wp, &
   & -5.30263009400631E-03_wp,  6.96070362833566E-03_wp,  2.06202786830572E-02_wp, &
   & -1.90384369106144E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  9.99999999830206E-01_wp, -1.01512838117467E-04_wp,  3.30286014739083E-04_wp, &
   &  2.21688220693919E-04_wp, -6.77924535368497E-06_wp,  4.42393841717257E-06_wp, &
   & -1.45294009909431E-04_wp,  8.28198172373844E-06_wp,  2.82716000311223E-06_wp, &
   &  1.11625581646469E-04_wp,  9.32081069771156E-05_wp, -1.43467330644811E-04_wp, &
   & -2.01611185625156E-04_wp, -2.63356233815959E-04_wp,  8.63251486754046E-05_wp, &
   &  6.13425902955612E-05_wp, -6.09820101131119E-06_wp,  1.20892448538428E-04_wp, &
   &  6.12520710210726E-05_wp,  8.67139986554086E-02_wp,  3.52869972519298E-01_wp, &
   &  1.77431941562336E-01_wp, -2.41645444821936E-01_wp, -8.86514552275003E-02_wp, &
   &  6.50937154430233E-02_wp, -3.59986382289629E-02_wp, -4.45761925024143E-02_wp, &
   & -3.43736973355472E-02_wp,  7.40406770098628E-03_wp, -3.06985789568106E-02_wp, &
   & -6.45075714555441E-03_wp, -2.25774856886379E-02_wp,  4.41748009542107E-03_wp, &
   & -1.19030608012565E-02_wp,  3.15424946469630E-02_wp,  9.52137814488968E-05_wp, &
   &  3.36928188670201E-04_wp,  1.89691378004237E-04_wp, -4.99549878650670E-06_wp, &
   & -3.01084193947362E-06_wp,  1.14329075205538E-04_wp, -2.14716191083883E-05_wp, &
   & -1.69511123039559E-06_wp, -1.01512838117467E-04_wp,  9.99999999869332E-01_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  8.24984479851205E-05_wp, &
   &  1.96296194806693E-04_wp,  5.54539527348801E-07_wp, -2.76572801639242E-04_wp, &
   & -1.09131970490034E-04_wp,  2.18813964986787E-07_wp, -6.67465544040966E-05_wp, &
   & -3.08299360534131E-07_wp,  3.81531587038920E-05_wp, -1.13584328085048E-01_wp, &
   & -5.63511977455039E-02_wp, -1.03417073813248E-01_wp,  1.40844227841612E-01_wp, &
   &  1.07106300572798E-01_wp, -7.86444738413934E-02_wp,  2.75278884074025E-02_wp, &
   &  7.22242203094389E-02_wp,  2.34893769725024E-03_wp, -2.55601123198589E-02_wp, &
   &  4.54124683874554E-02_wp, -1.77211746261233E-03_wp,  2.85240889524626E-02_wp, &
   & -5.59343295542886E-04_wp,  4.39408073360095E-03_wp, -2.93108178909011E-02_wp, &
   & -3.36928188670201E-04_wp, -9.41758669062715E-04_wp, -5.84225011720663E-04_wp, &
   &  1.53854928347468E-05_wp,  1.07377169272165E-05_wp, -4.07737526840130E-04_wp, &
   &  4.66396749143034E-05_wp,  6.57520964570179E-06_wp,  3.30286014739083E-04_wp, &
   &  0.00000000000000E+00_wp,  9.99999999998060E-01_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   & -1.96296194806693E-04_wp, -3.22041577165287E-04_wp, -1.14683436360897E-06_wp, &
   &  5.71975805721029E-04_wp,  2.33926970269965E-04_wp, -4.69032930059591E-07_wp, &
   &  1.53453012266193E-04_wp,  8.06130068702386E-07_wp, -1.36261835578487E-04_wp, &
   & -5.71130711953149E-02_wp, -1.03417073813248E-01_wp,  9.73200705922398E-02_wp, &
   &  7.08200378325435E-02_wp,  7.22242203094389E-02_wp,  3.13009539865045E-04_wp, &
   &  6.03006889788791E-02_wp, -2.14349010640352E-04_wp,  2.80042045879740E-02_wp, &
   & -7.00152224711211E-03_wp,  4.52717213649874E-03_wp,  2.88573809092836E-02_wp, &
   &  8.00111168452365E-03_wp, -1.97615416138547E-02_wp,  1.75536480937154E-03_wp, &
   & -2.98275876071079E-02_wp, -1.89691378004237E-04_wp, -5.84225011720663E-04_wp, &
   & -2.32983327085689E-04_wp,  8.66206935257109E-06_wp,  6.57520964570179E-06_wp, &
   & -1.86202195665234E-04_wp,  8.81555136696621E-05_wp,  2.76074508981215E-06_wp, &
   &  2.21688220693919E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  9.99999999998060E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp, -5.54539527348801E-07_wp, -1.14683436360897E-06_wp, &
   &  8.39122168161397E-05_wp,  1.61583974295499E-06_wp,  8.06130068702386E-07_wp, &
   &  3.64987134453440E-05_wp,  6.12105527417005E-07_wp, -5.14250999299680E-05_wp, &
   & -2.81827665248565E-07_wp,  7.77825760830674E-02_wp,  1.40844227841612E-01_wp, &
   &  7.08200378325435E-02_wp,  5.28706271385305E-02_wp, -4.50177541095921E-02_wp, &
   &  7.22242203094389E-02_wp, -1.88511048183665E-02_wp, -2.26360646654454E-02_wp, &
   & -7.46695782592206E-02_wp, -2.04737585704892E-02_wp, -2.15442376617703E-02_wp, &
   & -5.59343295542893E-04_wp, -1.95333031990455E-02_wp, -2.20587677537982E-03_wp, &
   & -2.76499015081222E-02_wp,  2.89865919149946E-02_wp,  4.99549878650671E-06_wp, &
   &  1.53854928347468E-05_wp,  8.66206935257109E-06_wp,  9.57086406720181E-05_wp, &
   &  6.33015324737798E-05_wp,  6.57520964570179E-06_wp, -6.91507707790898E-07_wp, &
   &  3.56389145477078E-05_wp, -6.77924535368498E-06_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  9.99999999998060E-01_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.76572801639242E-04_wp, &
   &  5.71975805721029E-04_wp,  1.61583974295499E-06_wp, -7.21973578621338E-04_wp, &
   & -3.65551535952679E-04_wp,  8.06130068702386E-07_wp, -2.16208620672635E-04_wp, &
   & -1.03268826055672E-06_wp,  8.91322458394565E-05_wp, -5.83831463741125E-02_wp, &
   & -8.05377110703100E-03_wp, -2.13786714107248E-02_wp, -2.12103821177911E-02_wp, &
   &  4.63175554736653E-02_wp, -6.21718876923976E-02_wp,  4.41008058157470E-03_wp, &
   &  2.92431802609242E-02_wp,  4.19093396272880E-02_wp,  5.71402985740159E-03_wp, &
   &  2.77251308119991E-02_wp, -6.89049708634549E-03_wp,  1.78063068515174E-02_wp, &
   &  4.95776309493032E-03_wp,  2.10798636529093E-02_wp, -2.83688082507035E-02_wp, &
   & -3.01084193947362E-06_wp, -1.07377169272165E-05_wp, -6.57520964570179E-06_wp, &
   & -6.33015324737798E-05_wp, -3.73996719315490E-05_wp, -4.64599908927517E-06_wp, &
   &  2.45641225072744E-07_wp, -2.26199417185387E-05_wp,  4.42393841717257E-06_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  9.99999999830205E-01_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   & -1.09131970490034E-04_wp, -2.33926970269965E-04_wp, -8.06130068702387E-07_wp, &
   &  3.65551535952679E-04_wp,  1.55708719431824E-04_wp, -3.37779234265429E-07_wp, &
   &  9.68084081921723E-05_wp,  5.17901618054598E-07_wp, -6.93663114292701E-05_wp, &
   &  4.28687369766487E-02_wp,  5.91360720857849E-03_wp, -3.46284778301005E-02_wp, &
   & -2.13786714107248E-02_wp, -6.21718876923976E-02_wp,  7.29603146759238E-03_wp, &
   & -4.67218548895194E-02_wp, -5.71517669530143E-03_wp, -5.66727024069538E-03_wp, &
   &  1.64261267820631E-02_wp, -9.78968090562189E-03_wp, -2.15630999197224E-02_wp, &
   & -1.19714716143595E-02_wp,  1.93252201033651E-02_wp, -2.01631062617651E-03_wp, &
   &  2.34835834975460E-02_wp,  1.14329075205538E-04_wp,  4.07737526840130E-04_wp, &
   &  1.86202195665234E-04_wp, -6.57520964570179E-06_wp, -4.64599908927517E-06_wp, &
   &  1.38897990982117E-04_wp, -4.86344418783074E-05_wp, -2.20750442352880E-06_wp, &
   & -1.45294009909431E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  9.99999999830205E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  2.18813964986787E-07_wp,  4.69032930059591E-07_wp, &
   & -3.64987134453440E-05_wp, -8.06130068702386E-07_wp, -3.37779234265429E-07_wp, &
   & -1.27556892913121E-05_wp, -2.67921067336703E-07_wp,  2.12542538776661E-05_wp, &
   &  1.81700110721510E-07_wp, -2.37076059225061E-02_wp,  1.17908981123966E-02_wp, &
   & -3.79011924083843E-02_wp, -8.07440995582129E-03_wp,  4.41008058157471E-03_wp, &
   & -4.67218548895194E-02_wp, -1.18825152389580E-02_wp,  3.19951378324389E-02_wp, &
   &  1.70996375352672E-03_wp,  6.55758393614615E-04_wp,  2.20763678154631E-02_wp, &
   & -1.47401248948007E-02_wp,  1.69654313272323E-02_wp,  1.00940411888978E-02_wp, &
   &  8.55988639565554E-03_wp,  2.79363404761070E-03_wp, -2.14716191083883E-05_wp, &
   & -4.66396749143034E-05_wp, -8.81555136696621E-05_wp,  6.91507707790898E-07_wp, &
   &  2.45641225072744E-07_wp, -4.86344418783074E-05_wp, -2.10150359335817E-05_wp, &
   &  7.21083315540951E-07_wp,  8.28198172373845E-06_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  9.99999999830206E-01_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -6.67465544040966E-05_wp, &
   & -1.53453012266193E-04_wp, -6.12105527417006E-07_wp,  2.16208620672635E-04_wp, &
   &  9.68084081921723E-05_wp, -2.67921067336703E-07_wp,  6.15385941191151E-05_wp, &
   &  3.77489132096825E-07_wp, -3.38447711064139E-05_wp, -2.93565217287235E-02_wp, &
   & -2.13786714107248E-02_wp,  2.37135902185714E-02_wp, -1.06651162567620E-02_wp, &
   &  2.92431802609242E-02_wp, -5.71517669530143E-03_wp,  3.19951378324389E-02_wp, &
   &  2.86402855372428E-03_wp,  3.82652025183376E-02_wp,  1.09396897005001E-02_wp, &
   &  5.41731475723262E-03_wp,  1.93252201033651E-02_wp,  8.19806673481372E-03_wp, &
   & -6.57678604430529E-03_wp,  4.69913239680737E-03_wp, -2.12898936281053E-02_wp, &
   & -1.69511123039559E-06_wp, -6.57520964570179E-06_wp, -2.76074508981216E-06_wp, &
   & -3.56389145477078E-05_wp, -2.26199417185387E-05_wp, -2.20750442352880E-06_wp, &
   &  7.21083315540952E-07_wp, -9.95740806845498E-06_wp,  2.82716000311223E-06_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  9.99999999830205E-01_wp,  0.00000000000000E+00_wp, &
   & -3.08299360534131E-07_wp, -8.06130068702386E-07_wp,  5.14250999299680E-05_wp, &
   &  1.03268826055672E-06_wp,  5.17901618054598E-07_wp,  2.12542538776661E-05_wp, &
   &  3.77489132096825E-07_wp, -2.76169058108372E-05_wp, -1.35914068297143E-07_wp, &
   & -2.26374693772444E-02_wp,  3.38405303098482E-02_wp, -8.28936173267023E-03_wp, &
   &  4.57526496102248E-02_wp,  4.19093396272880E-02_wp, -5.66727024069539E-03_wp, &
   &  1.70996375352672E-03_wp,  3.82652025183376E-02_wp, -4.55187727956170E-02_wp, &
   & -3.25759068310625E-02_wp,  2.10798636529092E-02_wp, -3.00247894938940E-03_wp, &
   &  6.90421382037361E-03_wp,  1.43932135734599E-03_wp, -1.84673717546551E-02_wp, &
   & -2.83300770467367E-03_wp, -1.01512838117467E-04_wp, -3.30286014739083E-04_wp, &
   & -2.21688220693919E-04_wp,  6.77924535368497E-06_wp,  4.42393841717257E-06_wp, &
   & -1.45294009909431E-04_wp,  8.28198172373844E-06_wp,  2.82716000311223E-06_wp, &
   &  1.11625581646469E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  9.99999999830206E-01_wp,  3.81531587038920E-05_wp,  1.36261835578487E-04_wp, &
   &  2.81827665248565E-07_wp, -8.91322458394565E-05_wp, -6.93663114292701E-05_wp, &
   &  1.81700110721510E-07_wp, -3.38447711064139E-05_wp, -1.35914068297143E-07_wp, &
   & -1.84534116858889E-05_wp,  8.73543363380361E-02_wp, -1.11927689146485E-01_wp, &
   &  1.76687050868385E-01_wp,  4.14057132845683E-01_wp, -4.83242438420612E-02_wp, &
   & -2.06209903237694E-02_wp, -3.65827353559851E-02_wp,  7.62837881761689E-02_wp, &
   &  8.28521011474394E-02_wp, -2.32970808181820E-02_wp, -1.66394753076522E-02_wp, &
   &  2.12528120885591E-03_wp, -2.26881057605045E-02_wp, -7.86211035482033E-03_wp, &
   &  2.85284441435993E-02_wp,  2.29901709038050E-02_wp,  1.12758480338272E-04_wp, &
   &  1.31040324396565E-04_wp,  2.23144270067549E-04_wp,  3.69794278272637E-04_wp, &
   &  7.55446878381371E-05_wp,  4.55857897636091E-05_wp, -2.44516968961120E-05_wp, &
   &  1.28642570924251E-04_wp,  9.32081069771156E-05_wp,  8.24984479851205E-05_wp, &
   & -1.96296194806693E-04_wp, -5.54539527348801E-07_wp,  2.76572801639242E-04_wp, &
   & -1.09131970490034E-04_wp,  2.18813964986787E-07_wp, -6.67465544040966E-05_wp, &
   & -3.08299360534131E-07_wp,  3.81531587038920E-05_wp,  9.99999999869332E-01_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  3.61731541134506E-02_wp, &
   &  1.29184921393787E-01_wp,  3.25455034102347E-02_wp,  7.62687348214290E-02_wp, &
   &  3.82018892333236E-02_wp,  1.63015647136606E-02_wp, -8.92221837554666E-03_wp, &
   &  3.90978820690418E-02_wp,  5.94862619285365E-02_wp,  1.42279516174685E-02_wp, &
   &  1.12164550839581E-02_wp, -2.57616878551625E-03_wp, -9.07844125129397E-03_wp, &
   & -4.47883417102248E-04_wp,  2.26838862024821E-02_wp,  3.00639987805416E-02_wp, &
   & -1.31040324396565E-04_wp, -2.04150976546664E-05_wp, -2.26089711420990E-04_wp, &
   & -3.74675458323373E-04_wp, -2.52051845112027E-05_wp, -1.52095173725851E-05_wp, &
   &  1.65836799186663E-05_wp, -1.63404821497647E-04_wp, -1.43467330644811E-04_wp, &
   &  1.96296194806693E-04_wp, -3.22041577165287E-04_wp, -1.14683436360897E-06_wp, &
   &  5.71975805721029E-04_wp, -2.33926970269965E-04_wp,  4.69032930059591E-07_wp, &
   & -1.53453012266193E-04_wp, -8.06130068702386E-07_wp,  1.36261835578487E-04_wp, &
   &  0.00000000000000E+00_wp,  9.99999999998060E-01_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   & -5.71022949696438E-02_wp,  3.25455034102347E-02_wp,  9.84260921377083E-02_wp, &
   & -1.20396462500221E-01_wp,  3.90978820690418E-02_wp, -3.38005258954589E-04_wp, &
   &  6.06254946164389E-02_wp,  1.25039201181343E-03_wp, -6.70334685509428E-02_wp, &
   &  2.18832031916388E-02_wp,  2.25211882161450E-03_wp, -9.21701378652493E-03_wp, &
   &  7.72420256658115E-03_wp,  3.40967488112168E-02_wp, -3.86126634519724E-03_wp, &
   & -2.15949193473986E-02_wp, -2.23144270067549E-04_wp, -2.26089711420990E-04_wp, &
   & -2.72645831427766E-04_wp, -6.38022547981310E-04_wp, -1.63404821497647E-04_wp, &
   & -7.35308666900574E-05_wp,  1.02189065578875E-04_wp, -2.07503254464807E-04_wp, &
   & -2.01611185625156E-04_wp,  5.54539527348801E-07_wp, -1.14683436360897E-06_wp, &
   &  8.39122168161397E-05_wp,  1.61583974295499E-06_wp, -8.06130068702386E-07_wp, &
   & -3.64987134453440E-05_wp, -6.12105527417005E-07_wp,  5.14250999299680E-05_wp, &
   &  2.81827665248565E-07_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  9.99999999998060E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp, -1.33816329028274E-01_wp,  7.62687348214290E-02_wp, &
   & -1.20396462500221E-01_wp, -1.32341154522262E-01_wp,  7.46019922122534E-02_wp, &
   &  3.90978820690418E-02_wp,  3.30062041606792E-02_wp, -1.17765372388225E-01_wp, &
   & -9.41199197927978E-02_wp,  3.48981234262437E-02_wp,  3.09187955546773E-02_wp, &
   & -4.47883417102249E-04_wp,  3.35841236773789E-02_wp, -1.04037311751581E-03_wp, &
   & -4.16800783537248E-02_wp, -2.25160768370592E-02_wp, -3.69794278272637E-04_wp, &
   & -3.74675458323373E-04_wp, -6.38022547981310E-04_wp, -9.44974883405076E-04_wp, &
   & -2.45721973530458E-04_wp, -1.63404821497647E-04_wp,  4.67989527259473E-05_wp, &
   & -4.18431888622898E-04_wp, -2.63356233815959E-04_wp, -2.76572801639242E-04_wp, &
   &  5.71975805721029E-04_wp,  1.61583974295499E-06_wp, -7.21973578621338E-04_wp, &
   &  3.65551535952679E-04_wp, -8.06130068702386E-07_wp,  2.16208620672635E-04_wp, &
   &  1.03268826055672E-06_wp, -8.91322458394565E-05_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  9.99999999998060E-01_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -3.17825201816373E-02_wp, &
   & -5.17690901910778E-02_wp, -1.14661853993300E-02_wp, -1.09127406617264E-02_wp, &
   & -3.05116894672468E-02_wp, -1.58559058503109E-02_wp,  2.58807011549745E-03_wp, &
   & -4.23724796435138E-02_wp, -5.46710758908302E-02_wp, -1.02121852619804E-02_wp, &
   & -1.08132873910535E-02_wp,  1.72123004852782E-03_wp,  9.71809772688051E-03_wp, &
   & -4.15510671194854E-03_wp, -2.73626654703935E-02_wp, -3.08515669011558E-02_wp, &
   &  7.55446878381371E-05_wp,  2.52051845112026E-05_wp,  1.63404821497647E-04_wp, &
   &  2.45721973530458E-04_wp,  2.47089644023801E-05_wp,  1.61116888311556E-05_wp, &
   & -4.94256032779233E-06_wp,  1.09891385286233E-04_wp,  8.63251486754047E-05_wp, &
   & -1.09131970490034E-04_wp,  2.33926970269965E-04_wp,  8.06130068702387E-07_wp, &
   & -3.65551535952679E-04_wp,  1.55708719431824E-04_wp, -3.37779234265429E-07_wp, &
   &  9.68084081921723E-05_wp,  5.17901618054598E-07_wp, -6.93663114292701E-05_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  9.99999999830205E-01_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   & -1.35622823871297E-02_wp, -2.20909800759548E-02_wp,  1.10648156497081E-02_wp, &
   & -1.14661853993300E-02_wp, -1.58559058503109E-02_wp, -1.20233934496220E-04_wp, &
   &  1.48991414842628E-02_wp, -2.81564588294364E-03_wp, -3.12937259411051E-02_wp, &
   & -6.78822012440303E-03_wp, -8.10584518942798E-04_wp,  3.87219849948905E-03_wp, &
   &  3.71972593593363E-03_wp,  1.05925682350414E-02_wp, -3.96213612481712E-03_wp, &
   & -1.89302474912948E-02_wp,  4.55857897636091E-05_wp,  1.52095173725851E-05_wp, &
   &  7.35308666900573E-05_wp,  1.63404821497647E-04_wp,  1.61116888311556E-05_wp, &
   &  7.73094926501920E-06_wp, -1.90066643814781E-05_wp,  5.65625287713878E-05_wp, &
   &  6.13425902955612E-05_wp,  2.18813964986787E-07_wp, -4.69032930059591E-07_wp, &
   &  3.64987134453440E-05_wp,  8.06130068702386E-07_wp, -3.37779234265429E-07_wp, &
   & -1.27556892913121E-05_wp, -2.67921067336703E-07_wp,  2.12542538776661E-05_wp, &
   &  1.81700110721510E-07_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  9.99999999830205E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp, -2.40602114448212E-02_wp, -3.71443619535532E-03_wp, &
   & -3.77677054156251E-02_wp,  1.37409144503486E-02_wp,  2.58807011549745E-03_wp, &
   &  1.48991414842628E-02_wp, -1.14945896141096E-02_wp, -5.51167977457504E-02_wp, &
   & -4.43725612524170E-03_wp, -1.93575430061701E-03_wp,  1.20817835115180E-02_wp, &
   &  4.60884129645277E-03_wp,  1.72717431621451E-02_wp, -1.70496114723901E-02_wp, &
   & -2.07142641033214E-02_wp,  1.91025315773589E-03_wp, -2.44516968961120E-05_wp, &
   & -1.65836799186663E-05_wp, -1.02189065578875E-04_wp, -4.67989527259473E-05_wp, &
   & -4.94256032779232E-06_wp, -1.90066643814781E-05_wp, -2.64669742429983E-05_wp, &
   & -5.36365868268805E-05_wp, -6.09820101131120E-06_wp, -6.67465544040966E-05_wp, &
   &  1.53453012266193E-04_wp,  6.12105527417006E-07_wp, -2.16208620672635E-04_wp, &
   &  9.68084081921723E-05_wp, -2.67921067336703E-07_wp,  6.15385941191151E-05_wp, &
   &  3.77489132096825E-07_wp, -3.38447711064139E-05_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  9.99999999830206E-01_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  5.01713186690479E-02_wp, &
   & -1.14661853993300E-02_wp, -4.09323723050174E-02_wp,  1.72266574885550E-02_wp, &
   & -4.23724796435138E-02_wp, -2.81564588294364E-03_wp, -5.51167977457504E-02_wp, &
   &  9.53463736577521E-03_wp,  5.68398395712986E-02_wp, -2.57527981119962E-02_wp, &
   & -6.24737443153230E-03_wp,  1.05925682350414E-02_wp, -1.37604829309815E-02_wp, &
   & -3.24497991168268E-02_wp,  9.26442849949040E-03_wp,  1.84855092937346E-02_wp, &
   &  1.28642570924251E-04_wp,  1.63404821497647E-04_wp,  2.07503254464807E-04_wp, &
   &  4.18431888622897E-04_wp,  1.09891385286233E-04_wp,  5.65625287713878E-05_wp, &
   & -5.36365868268805E-05_wp,  1.47306249222992E-04_wp,  1.20892448538428E-04_wp, &
   & -3.08299360534131E-07_wp,  8.06130068702386E-07_wp, -5.14250999299680E-05_wp, &
   & -1.03268826055672E-06_wp,  5.17901618054598E-07_wp,  2.12542538776661E-05_wp, &
   &  3.77489132096825E-07_wp, -2.76169058108372E-05_wp, -1.35914068297143E-07_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  9.99999999830205E-01_wp,  0.00000000000000E+00_wp, &
   &  5.44912525773983E-02_wp, -2.84111485232627E-02_wp,  1.96588187822549E-02_wp, &
   & -1.29632440844007E-02_wp, -5.46710758908302E-02_wp, -3.12937259411051E-02_wp, &
   & -4.43725612524171E-03_wp,  5.68398395712986E-02_wp,  3.13346757662612E-02_wp, &
   & -2.95098567978740E-02_wp, -2.73626654703935E-02_wp, -1.67051092576148E-03_wp, &
   & -1.66617157727234E-02_wp,  7.47010245689617E-03_wp,  2.01405830947493E-02_wp, &
   &  5.24876072680661E-03_wp,  9.32081069771156E-05_wp,  1.43467330644811E-04_wp, &
   &  2.01611185625156E-04_wp,  2.63356233815959E-04_wp,  8.63251486754046E-05_wp, &
   &  6.13425902955612E-05_wp, -6.09820101131119E-06_wp,  1.20892448538428E-04_wp, &
   &  6.12520710210726E-05_wp,  3.81531587038920E-05_wp, -1.36261835578487E-04_wp, &
   & -2.81827665248565E-07_wp,  8.91322458394565E-05_wp, -6.93663114292701E-05_wp, &
   &  1.81700110721510E-07_wp, -3.38447711064139E-05_wp, -1.35914068297143E-07_wp, &
   & -1.84534116858889E-05_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
   &  9.99999999830206E-01_wp],&
   & shape(overlap))

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "CeCl3")
   call test_overlap_diat_mol(error, mol, make_gen_basis, .false., overlap)

end subroutine test_overlap_diat_cecl3

subroutine test_overlap_grad_ss(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat, i
   type(cgto_type) :: cgtoi, cgtoj
   type(cgto_cache) :: icache, jcache
   real(wp) :: vec(3), r2
   real(wp) :: overlap(1, 1), doverlapi(3, 1, 1), doverlapj(3, 1, 1), sr(1, 1), sl(1, 1)
   real(wp), parameter :: step = 1.0e-5_wp

   call slater_to_gauss(ng, 2, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 1, 0, 1.0_wp, cgtoj, .true., stat)

   call cgtoi%update(icache, .true.)
   call cgtoj%update(jcache, .true.)

   ! Randomly oriented vector of length 0.5
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 0.5_wp * vec / sqrt(sum(vec**2))
   r2 = sum(vec**2)

   call overlap_grad_cgto(cgtoi, cgtoj, icache, jcache, r2, vec, 100.0_wp, &
      & overlap, doverlapj)

   vec(:) = -vec

   call overlap_grad_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & overlap, doverlapi)

   do i = 1, 3
      call check(error, doverlapi(i, 1, 1), -doverlapj(i, 1, 1), thr=thr1*10)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return

   do i = 1, 3
      vec(i) = vec(i) + step
      r2 = sum(vec**2)
      call overlap_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, sr)
      vec(i) = vec(i) - 2*step
      r2 = sum(vec**2)
      call overlap_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, sl)
      vec(i) = vec(i) + step
      doverlapj(i, :, :) = 0.5_wp * (sr - sl) / step
   end do

   do i = 1, 3
      call check(error, doverlapi(i, 1, 1), doverlapj(i, 1, 1), thr=thr1*10)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return

end subroutine test_overlap_grad_ss


subroutine test_overlap_grad_pp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat, i, j, k
   type(cgto_type) :: cgtoi, cgtoj
   type(cgto_cache) :: icache, jcache
   real(wp) :: vec(3), r2
   real(wp) :: overlap(3, 3), doverlapi(3, 3, 3), doverlapj(3, 3, 3), sr(3, 3), sl(3, 3)
   real(wp), parameter :: step = 1.0e-5_wp

   call slater_to_gauss(ng, 3, 1, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 2, 1, 1.0_wp, cgtoj, .true., stat)

   call cgtoi%update(icache, .true.)
   call cgtoj%update(jcache, .true.)

   ! Randomly oriented vector of length 0.5
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 0.5_wp * vec / sqrt(sum(vec**2))
   r2 = sum(vec**2)

   call overlap_grad_cgto(cgtoi, cgtoj, icache, jcache, r2, vec, 100.0_wp, &
      & overlap, doverlapj)

   vec(:) = -vec

   call overlap_grad_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & overlap, doverlapi)

   lp: do i = 1, 3
      do j = 1, 3
         do k = 1, 3
            call check(error, doverlapi(i, j, k), -doverlapj(i, k, j), thr=thr1*10)
            if (allocated(error)) exit
         end do
      end do
   end do lp
   if (allocated(error)) return

   do i = 1, 3
      vec(i) = vec(i) + step
      r2 = sum(vec**2)
      call overlap_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, sr)
      vec(i) = vec(i) - 2*step
      r2 = sum(vec**2)
      call overlap_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, sl)
      vec(i) = vec(i) + step
      doverlapj(i, :, :) = 0.5_wp * (sr - sl) / step
   end do

   num: do i = 1, 3
      do j = 1, 3
         do k = 1, 3
            call check(error, doverlapi(i, j, k), doverlapj(i, j, k), thr=thr1*10)
            if (allocated(error)) exit num
         end do
      end do
   end do num
   if (allocated(error)) return

end subroutine test_overlap_grad_pp


subroutine test_overlap_grad_dd(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat, i, j, k
   type(cgto_type) :: cgtoi, cgtoj
   type(cgto_cache) :: icache, jcache
   real(wp) :: vec(3), r2
   real(wp) :: overlap(5, 5), doverlapi(3, 5, 5), doverlapj(3, 5, 5), sl(5, 5), sr(5, 5)
   real(wp), parameter :: step = 1.0e-5_wp

   call slater_to_gauss(ng, 4, 2, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 2, 1.0_wp, cgtoj, .true., stat)

   call cgtoi%update(icache, .true.)
   call cgtoj%update(jcache, .true.)

   ! Randomly oriented vector of length 0.5
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 0.5_wp * vec / sqrt(sum(vec**2))
   r2 = sum(vec**2)

   call overlap_grad_cgto(cgtoi, cgtoj, icache, jcache, r2, vec, 100.0_wp, &
      & overlap, doverlapj)

   vec(:) = -vec

   call overlap_grad_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & overlap, doverlapi)

   lp: do i = 1, 3
      do j = 1, 5
         do k = 1, 5
            call check(error, doverlapi(i, j, k), -doverlapj(i, k, j), thr=thr1*10)
            if (allocated(error)) exit lp
         end do
      end do
   end do lp
   if (allocated(error)) return

   do i = 1, 3
      vec(i) = vec(i) + step
      r2 = sum(vec**2)
      call overlap_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, sr)
      vec(i) = vec(i) - 2*step
      r2 = sum(vec**2)
      call overlap_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, sl)
      vec(i) = vec(i) + step
      doverlapj(i, :, :) = 0.5_wp * (sr - sl) / step
   end do

   num: do i = 1, 3
      do j = 1, 5
         do k = 1, 5
            call check(error, doverlapi(i, j, k), doverlapj(i, j, k), thr=thr1*10)
            if (allocated(error)) exit num
         end do
      end do
   end do num
   if (allocated(error)) return

end subroutine test_overlap_grad_dd


subroutine test_overlap_grad_ff(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6
   integer :: stat, i, j, k
   type(cgto_type) :: cgtoi, cgtoj
   type(cgto_cache) :: icache, jcache
   real(wp) :: vec(3), r2
   real(wp) :: overlap(7, 7), doverlapi(3, 7, 7), doverlapj(3, 7, 7), sl(7, 7), sr(7, 7)
   real(wp), parameter :: step = 1.0e-5_wp

   call slater_to_gauss(ng, 5, 3, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 4, 3, 1.0_wp, cgtoj, .true., stat)

   call cgtoi%update(icache, .true.)
   call cgtoj%update(jcache, .true.)

   ! Randomly oriented vector of length 0.5
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 0.5_wp * vec / sqrt(sum(vec**2))
   r2 = sum(vec**2)

   call overlap_grad_cgto(cgtoi, cgtoj, icache, jcache, r2, vec, 100.0_wp, &
      & overlap, doverlapj)

   vec(:) = -vec

   call overlap_grad_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & overlap, doverlapi)

   lp: do i = 1, 3
      do j = 1, 7
         do k = 1, 7
            call check(error, doverlapi(i, j, k), -doverlapj(i, k, j), thr=thr1*10)
            if (allocated(error)) exit lp
         end do
      end do
   end do lp
   if (allocated(error)) return

   do i = 1, 3
      vec(i) = vec(i) + step
      r2 = sum(vec**2)
      call overlap_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, sr)
      vec(i) = vec(i) - 2*step
      r2 = sum(vec**2)
      call overlap_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, sl)
      vec(i) = vec(i) + step
      doverlapj(i, :, :) = 0.5_wp * (sr - sl) / step
   end do

   num: do i = 1, 3
      do j = 1, 7
         do k = 1, 7
            call check(error, doverlapi(i, j, k), doverlapj(i, j, k), thr=thr1*10)
            if (allocated(error)) exit num
         end do
      end do
   end do num
   if (allocated(error)) return

end subroutine test_overlap_grad_ff


!> General routine to test charge gradients of overlap integrals for q-vSZP CGTO pairs
subroutine test_overlap_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in)
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
   real(wp), allocatable :: overlap(:, :), doverlap(:, :, :)
   real(wp), allocatable :: doverlapdqeffi(:, :), doverlapdqeffj(:, :)
   real(wp), allocatable :: sr(:, :), sl(:, :)
   real(wp), allocatable :: num_dSdqi(:, :), num_dSdqj(:, :)
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

   ! Analytical gradient calculation
   r2 = sum(vec**2)
   naoi = msao(cgtoi%ang)
   naoj = msao(cgtoj%ang)
   allocate(overlap(naoj, naoi), doverlap(3, naoj, naoi), doverlapdqeffi(naoj, naoi), &
      & doverlapdqeffj(naoj, naoi))
   call overlap_grad_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & overlap, doverlap, doverlapdqeffj, doverlapdqeffi)

   ! Numerical gradient calculation
   allocate(sr(naoj, naoi), sl(naoj, naoi), num_dSdqi(naoj, naoi), num_dSdqj(naoj, naoi))

   ! Right hand side
   sr = 0.0_wp
   qeffi = icache%qeff
   icache%qeff = qeffi + step
   call cgtoi%get_normalization(icache, .false.)
   call overlap_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, sr)
   
   ! Left hand side
   sl = 0.0_wp
   icache%qeff = qeffi - step
   call cgtoi%get_normalization(icache, .false.)
   call overlap_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, sl)
   
   icache%qeff = qeffi
   call cgtoi%get_normalization(icache, .false.)
   num_dSdqi = 0.5_wp * (sr - sl) / step

   ! Right hand side
   sr = 0.0_wp
   qeffj = jcache%qeff
   jcache%qeff = qeffj + step
   call cgtoj%get_normalization(jcache, .false.)
   call overlap_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, sr)
   
   ! Left hand side
   sl = 0.0_wp
   jcache%qeff = qeffj - step
   call cgtoj%get_normalization(jcache, .false.)
   call overlap_cgto(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, sl)
   
   jcache%qeff = qeffj
   call cgtoj%get_normalization(jcache, .false.)
   num_dSdqj = 0.5_wp * (sr - sl) / step

   if (any(abs(doverlapdqeffi - num_dSdqi) > thr_)) then
      call test_failed(error, "Charge gradient of CGTO at center i does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', doverlapdqeffi
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dSdqi
      write(*,*) 'Difference:'
      print'(3es21.14)', doverlapdqeffi-num_dSdqi
   end if

   if (any(abs(doverlapdqeffj - num_dSdqj) > thr_)) then
      call test_failed(error, "Charge gradient of CGTO at center j does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', doverlapdqeffj
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dSdqj
      write(*,*) 'Difference:'
      print'(3es21.14)', doverlapdqeffj-num_dSdqj
   end if

end subroutine test_overlap_grad_qeff_gen


subroutine test_overlap_grad_qeff_ss(error)
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

   call test_overlap_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_overlap_grad_qeff_ss

subroutine test_overlap_grad_qeff_sp(error)
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

   call test_overlap_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_overlap_grad_qeff_sp

subroutine test_overlap_grad_qeff_pp(error)
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

   call test_overlap_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_overlap_grad_qeff_pp

subroutine test_overlap_grad_qeff_sd(error)
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

   call test_overlap_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_overlap_grad_qeff_sd

subroutine test_overlap_grad_qeff_pd(error)
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

   call test_overlap_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_overlap_grad_qeff_pd

subroutine test_overlap_grad_qeff_dd(error)
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

   call test_overlap_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1*100)
end subroutine test_overlap_grad_qeff_dd

subroutine test_overlap_grad_qeff_sf(error)
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

   call test_overlap_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_overlap_grad_qeff_sf

subroutine test_overlap_grad_qeff_pf(error)
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

   call test_overlap_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_overlap_grad_qeff_pf

subroutine test_overlap_grad_qeff_df(error)
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

   call test_overlap_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_overlap_grad_qeff_df

subroutine test_overlap_grad_qeff_ff(error)
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

   call test_overlap_grad_qeff_gen(vec, cgtoi, cgtoj, error, thr_in=thr1)
end subroutine test_overlap_grad_qeff_ff


subroutine test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)
   !> Vector between the two centers
   real(wp), intent(inout) :: vec(3)
   !> Scaling factors for the diatomic frame overlap
   real(wp), intent(in) :: ksig, kpi, kdel
   !> CGTOs for the tested diatomic frame overlap calculation
   type(cgto_type), intent(in) :: cgtoi, cgtoj
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(cgto_cache) :: icache
   type(cgto_cache) :: jcache
   integer :: i, j, k
   real(wp) :: r2
   real(wp) :: overlap(msao(cgtoi%ang), msao(cgtoj%ang)), overlap_diat(msao(cgtoi%ang), msao(cgtoj%ang)), &
      & sl(msao(cgtoj%ang), msao(cgtoi%ang)), sr(msao(cgtoj%ang), msao(cgtoi%ang)), &
      & sl_diat(msao(cgtoj%ang), msao(cgtoi%ang)), sr_diat(msao(cgtoj%ang), msao(cgtoi%ang))
   real(wp) :: doverlapi(3, msao(cgtoi%ang), msao(cgtoj%ang)), doverlapi_diat(3, msao(cgtoi%ang), msao(cgtoj%ang)), &
      & doverlapj(3, msao(cgtoj%ang), msao(cgtoi%ang)), doverlapj_diat(3, msao(cgtoj%ang), msao(cgtoi%ang)), &
      & doverlaptmp(3, msao(cgtoj%ang), msao(cgtoi%ang)), doverlaptmp_diat(3, msao(cgtoj%ang), msao(cgtoi%ang))

   real(wp), parameter :: step = 1.0e-5_wp

   r2 = sum(vec**2)

   call cgtoi%update(icache, .true.)
   call cgtoj%update(jcache, .true.)

   ! Test antisymmetry w.r.t. the exchange of the two centers
   call overlap_grad_cgto_diat(cgtoi, cgtoj, icache, jcache, r2, vec, 100.0_wp, &
      & ksig, kpi, kdel, overlap, doverlapi, overlap_diat, doverlapi_diat)

   vec(:) = -vec

   call overlap_grad_cgto_diat(cgtoj, cgtoi, jcache, icache, r2, vec, 100.0_wp, &
      & ksig, kpi, kdel, overlap, doverlapj, overlap_diat, doverlapj_diat)

   lp: do i = 1, 3
      do j = 1, msao(cgtoi%ang)
         do k = 1, msao(cgtoj%ang)
            call check(error, doverlapi(i, j, k), -doverlapj(i, k, j), thr=thr2)
            if (allocated(error)) exit lp
            call check(error, doverlapi_diat(i, j, k), -doverlapj_diat(i, k, j), thr=thr2)
            if (allocated(error)) exit lp
         end do
      end do
   end do lp
   if (allocated(error)) return

   ! Test the analytical against the numerical gradient
   do i = 1, 3
      vec(i) = vec(i) + step
      r2 = sum(vec**2)
      call overlap_cgto_diat(cgtoj, cgtoi, jcache, icache, r2, vec, &
         & 100.0_wp, ksig, kpi, kdel, sr, sr_diat)

      vec(i) = vec(i) - 2*step
      r2 = sum(vec**2)
      call overlap_cgto_diat(cgtoj, cgtoi, jcache, icache, r2, vec, &
         & 100.0_wp, ksig, kpi, kdel, sl, sl_diat)

      vec(i) = vec(i) + step
      doverlaptmp(i, :, :) = 0.5_wp * (sr - sl) / step
      doverlaptmp_diat(i, :, :) = 0.5_wp * (sr_diat - sl_diat) / step
   end do

   num: do i = 1, 3
      do j = 1, msao(cgtoi%ang)
         do k = 1, msao(cgtoj%ang)
            call check(error, doverlapj(i, k, j), doverlaptmp(i, k, j), thr=thr2)
            if (allocated(error)) then
               write(*,*) "Error", doverlapj(i, k, j), doverlaptmp(i, k, j), &
                  & doverlapj(i, k, j) - doverlaptmp(i, k, j) 
               exit num
            end if
            call check(error, doverlapj_diat(i, k, j), doverlaptmp_diat(i, k, j), thr=thr2)
            if (allocated(error)) then
               write(*,*) "Error", doverlapj_diat(i, k, j), doverlaptmp_diat(i, k, j), &
                  & doverlapj_diat(i, k, j) - doverlaptmp_diat(i, k, j)
               exit num
            end if
         end do
      end do
   end do num
   if (allocated(error)) return

end subroutine test_overlap_diat_grad_gen

subroutine test_overlap_diat_grad_ss(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 2, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 1, 0, 2.0_wp, cgtoj, .true., stat)

   ! Randomly oriented vector of length 0.5
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 0.5_wp * vec / sqrt(sum(vec**2))

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_ss

subroutine test_overlap_diat_grad_ss_z(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 2, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 1, 0, 2.0_wp, cgtoj, .true., stat)

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   ! Vector along the x-axis
   vec(1) = 0.5_wp
   vec(2) = 0.0_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the y-axis
   vec(1) = 0.0_wp
   vec(2) = 0.5_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = -0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector nearly along the z-axis to test the ill-defined gradient
   vec(1) = 1e-7_wp
   vec(2) = 1e-7_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_ss_z

subroutine test_overlap_diat_grad_sp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 4, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 1, 2.0_wp, cgtoj, .true., stat)

   ! Randomly oriented vector of length 0.5
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 0.5_wp * vec / sqrt(sum(vec**2))

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_sp

subroutine test_overlap_diat_grad_sp_z(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 4, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 1, 2.0_wp, cgtoj, .true., stat)

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   ! Vector along the x-axis
   vec(1) = 0.5_wp
   vec(2) = 0.0_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the y-axis
   vec(1) = 0.0_wp
   vec(2) = 0.5_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = -0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector nearly along the z-axis to test the ill-defined gradient
   vec(1) = 1e-7_wp
   vec(2) = 1e-7_wp
   vec(3) = -0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_sp_z

subroutine test_overlap_diat_grad_pp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 3, 1, 1.5_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 2, 1, 1.0_wp, cgtoj, .true., stat)

   ! Randomly oriented vector of length 0.5
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 0.5_wp * vec / sqrt(sum(vec**2))

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_pp

subroutine test_overlap_diat_grad_pp_z(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 3, 1, 1.5_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 2, 1, 1.0_wp, cgtoj, .true., stat)

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   ! Vector along the x-axis
   vec(1) = 0.5_wp
   vec(2) = 0.0_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the y-axis
   vec(1) = 0.0_wp
   vec(2) = 0.5_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = -0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector nearly along the z-axis to test the ill-defined gradient
   vec(1) = 1e-7_wp
   vec(2) = 1e-7_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_pp_z

subroutine test_overlap_diat_grad_sd(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat


   call slater_to_gauss(ng, 4, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 2, 1.8_wp, cgtoj, .true., stat)

   ! Randomly oriented vector of length 0.5
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 0.5_wp * vec / sqrt(sum(vec**2))

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_sd

subroutine test_overlap_diat_grad_sd_z(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 4, 0, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 2, 1.8_wp, cgtoj, .true., stat)

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   ! Vector along the x-axis
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the y-axis
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = -0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector nearly along the z-axis to test the ill-defined gradient
   vec(1) = 1e-7_wp
   vec(2) = 1e-7_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_sd_z

subroutine test_overlap_diat_grad_pd(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 4, 1, 1.5_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 2, 1.5_wp, cgtoj, .true., stat)

   ! Randomly oriented vector of length 0.5
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 0.5_wp * vec / sqrt(sum(vec**2))

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_pd

subroutine test_overlap_diat_grad_pd_z(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 4, 1, 1.5_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 2, 1.5_wp, cgtoj, .true., stat)

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   ! Vector along the x-axis
   vec(1) = 0.5_wp
   vec(2) = 0.0_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the y-axis
   vec(1) = 0.0_wp
   vec(2) = 0.5_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = -0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector nearly along the z-axis to test the ill-defined gradient
   vec(1) = 1e-7_wp
   vec(2) = 1e-7_wp
   vec(3) = -0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_pd_z

subroutine test_overlap_diat_grad_dd(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 4, 2, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 2, 1.3_wp, cgtoj, .true., stat)

   ! Randomly oriented vector of length 0.7
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 0.7_wp * vec / sqrt(sum(vec**2))

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_dd

subroutine test_overlap_diat_grad_dd_z(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 4, 2, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 3, 2, 1.3_wp, cgtoj, .true., stat)

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   ! Vector along the x-axis
   vec(1) = 0.5_wp
   vec(2) = 0.0_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the y-axis
   vec(1) = 0.0_wp
   vec(2) = 0.5_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = -0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector nearly along the z-axis to test the ill-defined gradient
   vec(1) = -1e-7_wp
   vec(2) = -1e-7_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_dd_z

subroutine test_overlap_diat_grad_sf(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 5, 0, 2.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 4, 3, 1.0_wp, cgtoj, .true., stat)

   ! Randomly oriented vector of length 0.5
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 0.5_wp * vec / sqrt(sum(vec**2))

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_sf

subroutine test_overlap_diat_grad_sf_z(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 5, 0, 2.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 4, 3, 1.0_wp, cgtoj, .true., stat)

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   ! Vector along the x-axis
   vec(1) = 0.5_wp
   vec(2) = 0.0_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the y-axis
   vec(1) = 0.0_wp
   vec(2) = 0.5_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = -0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector nearly along the z-axis to test the ill-defined gradient
   vec(1) = -1e-7_wp
   vec(2) = 1e-7_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_sf_z

subroutine test_overlap_diat_grad_pf(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 5, 1, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 4, 3, 2.0_wp, cgtoj, .true., stat)


   ! Randomly oriented vector of length 0.5
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 0.5_wp * vec / sqrt(sum(vec**2))

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_pf

subroutine test_overlap_diat_grad_pf_z(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 5, 1, 1.0_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 4, 3, 2.0_wp, cgtoj, .true., stat)

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   ! Vector along the x-axis
   vec(1) = 0.5_wp
   vec(2) = 0.0_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the y-axis
   vec(1) = 0.0_wp
   vec(2) = 0.5_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = -0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector nearly along the z-axis to test the ill-defined gradient
   vec(1) = -1e-7_wp
   vec(2) = 1e-7_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_pf_z

subroutine test_overlap_diat_grad_df(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 5, 2, 1.3_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 4, 3, 1.4_wp, cgtoj, .true., stat)

   ! Randomly oriented vector of length 0.5
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 0.5_wp * vec / sqrt(sum(vec**2))

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_df

subroutine test_overlap_diat_grad_df_z(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 5, 2, 1.3_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 4, 3, 1.4_wp, cgtoj, .true., stat)

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   ! Vector along the x-axis
   vec(1) = 0.5_wp
   vec(2) = 0.0_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the y-axis
   vec(1) = 0.0_wp
   vec(2) = 0.5_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = -0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector nearly along the z-axis to test the ill-defined gradient
   vec(1) = 1e-7_wp
   vec(2) = -1e-7_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_df_z

subroutine test_overlap_diat_grad_ff(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 5, 3, 1.6_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 4, 3, 1.0_wp, cgtoj, .true., stat)

   ! Randomly oriented vector of length 0.5
   call random_number(vec)
   vec = vec - 0.5_wp
   vec = 0.5_wp * vec / sqrt(sum(vec**2))

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_ff

subroutine test_overlap_diat_grad_ff_z(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: ng = 6

   type(cgto_type) :: cgtoi, cgtoj
   real(wp) :: vec(3)
   real(wp) :: ksig, kpi, kdel
   integer :: stat

   call slater_to_gauss(ng, 5, 3, 1.6_wp, cgtoi, .true., stat)
   call slater_to_gauss(ng, 4, 3, 1.0_wp, cgtoj, .true., stat)

   ksig = 0.1_wp
   kpi = 0.2_wp
   kdel = 0.5_wp

   ! Vector along the x-axis
   vec(1) = 0.5_wp
   vec(2) = 0.0_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the y-axis
   vec(1) = 0.0_wp
   vec(2) = 0.5_wp
   vec(3) = 0.0_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector along the z-axis to test the ill-defined gradient
   vec(1) = 0.0_wp
   vec(2) = 0.0_wp
   vec(3) = -0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

   ! Vector nearly along the z-axis to test the ill-defined gradient
   vec(1) = -1e-7_wp
   vec(2) = -1e-7_wp
   vec(3) = 0.5_wp

   call test_overlap_diat_grad_gen(vec, ksig, kpi, kdel, cgtoi, cgtoj, error)

end subroutine test_overlap_diat_grad_ff_z

end module test_integral_overlap
