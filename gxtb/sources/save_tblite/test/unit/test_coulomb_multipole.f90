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

module test_coulomb_multipole
   use mctc_data_vdwrad, only : get_vdw_rad
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa 
   use mstore, only : get_structure
   use mctc_ncoord, only : cn_count
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_multipole, only : gfn2_multipole, &
      & new_gfn2_multipole, gxtb_multipole, new_gxtb_multipole
   use tblite_scf_potential, only : potential_type
   use tblite_utils_average, only : average_type, new_average, average_id
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_coulomb, only : tb_coulomb, new_coulomb
   implicit none
   private

   public :: collect_coulomb_multipole

   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine multipole_maker(coulomb, mol, error)
         import :: tb_coulomb, structure_type, error_type
         type(tb_coulomb), intent(out) :: coulomb
         type(structure_type), intent(in) :: mol
         type(error_type), allocatable, intent(out) :: error
      end subroutine multipole_maker
   end interface

   interface
      module subroutine test_e_aes_gfn2_co2(error)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine test_e_aes_gfn2_co2
      module subroutine test_e_aes_gfn2_co2_dp(error)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine test_e_aes_gfn2_co2_dp
      module subroutine test_e_aes_gfn2_co2_qp(error)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine test_e_aes_gfn2_co2_qp
      module subroutine test_e_aes_gfn2_co2_sc(error)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine test_e_aes_gfn2_co2_sc
      module subroutine test_e_aes_gfn2_co2_sc_dp(error)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine test_e_aes_gfn2_co2_sc_dp
      module subroutine test_e_aes_gfn2_co2_sc_qp(error)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine test_e_aes_gfn2_co2_sc_qp
   end interface

contains


!> Collect all exported unit tests
subroutine collect_coulomb_multipole(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-gfn2-m01", test_e_aes_gfn2_m01), &
      new_unittest("energy-gfn2-m02", test_e_aes_gfn2_m02), &
      new_unittest("energy-gfn2-co2-pbc", test_e_aes_gfn2_co2), &
      new_unittest("energy-gfn2-co2-pbc-qp", test_e_aes_gfn2_co2_qp), &
      new_unittest("energy-gfn2-co2-pbc-dp", test_e_aes_gfn2_co2_dp), &
      !new_unittest("energy-gfn2-co2-sc", test_e_aes_gfn2_co2_sc), &
      !new_unittest("energy-gfn2-co2-sc-qp", test_e_aes_gfn2_co2_sc_qp), &
      !new_unittest("energy-gfn2-co2-sc-dp", test_e_aes_gfn2_co2_sc_dp), &
      new_unittest("energy-gxtb-m01", test_e_aes_gxtb_m01), &
      new_unittest("energy-gxtb-m02", test_e_aes_gxtb_m02), &
      new_unittest("energy-gxtb-cecl3", test_e_aes_gxtb_cecl3), &
      new_unittest("potential-gfn2-m07", test_p_aes_gfn2_m07), &
      new_unittest("potential-gfn2-m08", test_p_aes_gfn2_m08), &
      new_unittest("potential-gxtb-m07", test_p_aes_gxtb_m07), &
      new_unittest("potential-gxtb-m08", test_p_aes_gxtb_m08), &
      new_unittest("potential-gxtb-cecl3", test_p_aes_gxtb_cecl3), &
      new_unittest("gradient-gfn2-m03", test_g_aes_gfn2_m03), &
      new_unittest("gradient-gfn2-m04", test_g_aes_gfn2_m04), &
      new_unittest("gradient-gfn2-urea-pbc", test_g_aes_gfn2_urea), &
      new_unittest("gradient-gxtb-m03", test_g_aes_gxtb_m03), &
      new_unittest("gradient-gxtb-m04", test_g_aes_gxtb_m04), &
      new_unittest("gradient-gxtb-cecl3", test_g_aes_gxtb_cecl3), &
      new_unittest("sigma-gfn2-m05", test_s_aes_gfn2_m05), &
      new_unittest("sigma-gfn2-m06", test_s_aes_gfn2_m06), &
      !new_unittest("sigma-gfn2-oxacb-pbc", test_s_aes_gfn2_oxacb) &
      new_unittest("sigma-gxtb-m05", test_s_aes_gxtb_m05), &
      new_unittest("sigma-gxtb-m06", test_s_aes_gxtb_m06), &
      new_unittest("sigma-gxtb-cecl3", test_s_aes_gxtb_cecl3) &
      ]

end subroutine collect_coulomb_multipole


!> Factory to create electrostatic objects based on GFN2-xTB values
subroutine make_multipole_gfn2(coulomb, mol, error)

   !> New electrostatic object
   type(tb_coulomb), intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: kexp3 = 3.0_wp, kexp5 = 3.0_wp
   real(wp), parameter :: shift = 1.2_wp, kradexp = 4.0_wp, rmax = 5.0_wp
   !> Dipole exchange-correlation kernel
   real(wp), parameter :: p_dkernel(20) = 0.01_wp * [&
      & 5.563889_wp,-1.000000_wp,-0.500000_wp,-0.613341_wp,-0.481186_wp, &
      &-0.411674_wp, 3.521273_wp,-4.935670_wp,-8.339183_wp,10.000000_wp, &
      & 0.000000_wp,-0.082005_wp, 2.633341_wp,-0.025750_wp, 2.110225_wp, &
      &-0.151117_wp,-2.536958_wp,-2.077329_wp,-0.103383_wp,-0.236675_wp]
   !> Quadrupole exchange-correlation kernel
   real(wp), parameter :: p_qkernel(20) = 0.01_wp * [&
      & 0.027431_wp,-0.337528_wp, 0.020000_wp,-0.058586_wp,-0.058228_wp, &
      & 0.213583_wp, 2.026786_wp,-0.310828_wp,-0.245955_wp,-0.500000_wp, &
      & 0.020000_wp,-0.005516_wp,-0.021887_wp,-0.080000_wp, 0.028679_wp, &
      & 0.442859_wp, 0.122783_wp,-1.083404_wp, 0.025000_wp, 0.010000_wp]
   real(wp), parameter :: p_rad(20) = [&
      & 1.4_wp, 3.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.0_wp, 1.9_wp, 1.8_wp, 2.4_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 3.9_wp, 2.1_wp, 3.1_wp, 2.5_wp, 5.0_wp, 5.0_wp, 5.0_wp]
   real(wp), parameter :: p_vcn(20) = [&
      & 1.0_wp, 1.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 2.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 2.0_wp]
   real(wp), allocatable :: dkernel(:), qkernel(:), rad(:), vcn(:)
   type(gfn2_multipole), allocatable :: aes2
   type(average_type), allocatable :: average

   dkernel = p_dkernel(mol%num)
   qkernel = p_qkernel(mol%num)
   rad = p_rad(mol%num)
   vcn = p_vcn(mol%num)

   ! Setup coulomb interaction collection 
   call new_coulomb(coulomb, mol, error, cn_count_type=cn_count%dexp)
   if (allocated(error)) return

   ! Setup average
   allocate(average)
   call new_average(average, average_id%arithmetic)

   ! Setup damped multipole electrostatics
   allocate(aes2)
   call new_gfn2_multipole(aes2, mol, dkernel, qkernel, shift, kradexp, rmax, &
      & rad, average, vcn, kdmp3=1.0_wp, kdmp5=1.0_wp, kexp3=kexp3, kexp5=kexp5)
   call move_alloc(aes2, coulomb%aes2)

end subroutine make_multipole_gfn2


!> Factory to create electrostatic objects based on g-xTB values
subroutine make_multipole_gxtb(coulomb, mol, error)

   !> New electrostatic object
   type(tb_coulomb), intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Radii in the coordination number in g-xTB
   real(wp), parameter :: p_cn_rcov(60) = 0.5_wp * [&
      &  0.9646730831_wp,  1.0887413026_wp,  2.7319904126_wp,  2.0642646118_wp, & !1-4
      &  1.8403302383_wp,  1.8469537700_wp,  1.7342318861_wp,  1.6797510318_wp, & !5-8
      &  1.2145258133_wp,  1.6166175254_wp,  3.1780229637_wp,  2.8181820446_wp, & !9-12
      &  2.6406018705_wp,  2.5882455353_wp,  2.6392688502_wp,  2.5174870894_wp, & !13-16
      &  2.1965435815_wp,  2.3352396673_wp,  3.8631034402_wp,  3.2237828672_wp, & !17-20
      &  3.0306140794_wp,  2.8070147869_wp,  2.7092132145_wp,  2.6725665527_wp, & !21-24
      &  2.4709027569_wp,  2.6864662999_wp,  2.5048118705_wp,  2.5534778256_wp, & !25-28
      &  2.5338150899_wp,  2.4255689950_wp,  2.6654702197_wp,  2.4751265075_wp, & !29-32
      &  2.4955198981_wp,  2.3976809191_wp,  2.4115345626_wp,  2.8803089804_wp, & !33-36
      &  4.1126610230_wp,  3.3773012295_wp,  2.9145452327_wp,  2.9593967529_wp, & !37-40
      &  2.6794258023_wp,  2.9716991394_wp,  2.6760216235_wp,  2.6320075085_wp, & !41-44
      &  2.6203520127_wp,  2.6918392535_wp,  2.5246085446_wp,  3.0987734452_wp, & !45-48
      &  2.7763347045_wp,  2.7985474546_wp,  2.8399171705_wp,  2.7057365260_wp, & !49-52
      &  2.7285345912_wp,  2.9622121765_wp,  4.4063616833_wp,  3.9116265616_wp, & !53-56
      &  3.2727781933_wp,  3.1827185065_wp,  3.1075853258_wp,  3.0554430377_wp] !57-60

   !> CN-dependent scaling factor for the atomic dipole moments
   real(wp), parameter :: p_aes_dip_scale(60) = [&
      &  5.7861581101_wp, 54.1660578473_wp,  1.3720398902_wp,  2.3283264813_wp, & !1-4
      &  1.0860367872_wp,  0.5103590165_wp,  1.1791825944_wp, -0.6697227790_wp, & !5-8
      & -5.2914673350_wp, 22.7319498845_wp, -0.2485040084_wp,  2.3320119727_wp, & !9-12
      &  1.8525079181_wp,  1.4084234621_wp,  0.4691763855_wp,  0.1724210314_wp, & !13-16
      & -0.9964379564_wp,  1.7872872247_wp, -1.4488064206_wp,  0.2135932485_wp, & !17-20
      &  0.5294779860_wp,  0.2107226865_wp,  0.2168718334_wp,  0.6206648872_wp, & !21-24
      &  0.1702421269_wp, -0.0225723599_wp, -0.3090043243_wp,  0.1222218815_wp, & !25-28
      &  1.2203073535_wp,  5.1993700933_wp,  1.3718278272_wp,  0.3261878349_wp, & !29-32
      &  0.0487258034_wp,  0.3191700250_wp,  0.0773617537_wp,  2.6938090758_wp, & !33-36
      &  0.1409781543_wp,  0.3807891469_wp,  1.8482722568_wp,  0.3825742006_wp, & !37-40
      &  0.3029974134_wp,  0.4992854964_wp,  0.2928175531_wp,  0.0000000000_wp, & !41-44
      &  0.3238105654_wp,  0.2520229533_wp,  0.0400361622_wp,  3.2198856481_wp, & !45-48
      &  2.4180335921_wp,  1.3421356654_wp,  0.9667177416_wp,  0.4000000000_wp, & !49-52
      &  0.3000000000_wp,  0.3000000000_wp,  0.7817916047_wp,  0.2971597943_wp, & !53-56
      &  0.9000000000_wp,  0.2512392363_wp,  0.0000000000_wp,  0.0000000000_wp] !57-60

   !> van-der-Waals radii scaling
   real(wp), parameter :: p_rvdw_scale(60) = autoaa * [&
      &  1.0000000000_wp,  1.0000000000_wp,  1.0600000000_wp,  1.0000000000_wp, & !1-4
      &  1.0200000000_wp,  1.0000000000_wp,  1.0500000000_wp,  1.0800000000_wp, & !5-8
      &  1.1500000000_wp,  0.8000000000_wp,  1.2000000000_wp,  1.1500000000_wp, & !9-12
      &  1.0700000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !13-16
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !17-20
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !21-24
      &  1.0000000000_wp,  1.1000000000_wp,  1.1000000000_wp,  1.1000000000_wp, & !25-28
      &  1.0000000000_wp,  1.1500000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !29-32
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !33-36
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !37-40
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !41-44
      &  1.1000000000_wp,  1.1000000000_wp,  1.1000000000_wp,  1.1000000000_wp, & !45-48
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !49-52
      &  1.0000000000_wp,  1.0000000000_wp,  0.6500000000_wp,  1.0000000000_wp, & !53-56
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp] !57-60

   !> Damping prefactor of multipole electrostatics
   real(wp), parameter :: p_aes_dmp_r2 = 0.3405910191_wp
   real(wp), parameter :: p_aes_dmp_r3 = 0.1691310614_wp
   real(wp), parameter :: p_aes_dmp_r4 = 0.0740343390_wp
   real(wp), parameter :: p_aes_dmp_r5 = -0.0200000000_wp
   !> Exponent in multipole damping
   real(wp), parameter :: p_aes_exp_r2 = 0.5000000000_wp
   real(wp), parameter :: p_aes_exp_r3 = 1.0000000000_wp
   real(wp), parameter :: p_aes_exp_r4 = 1.0000000000_wp
   real(wp), parameter :: p_aes_exp_r5 = 1.0000000000_wp
   !> Coordination number exponent
   real(wp), parameter :: p_cn_exp = 2.0680000000_wp

   type(gxtb_multipole), allocatable :: aes2
   type(average_type), allocatable :: average
   real(wp), allocatable :: tb_cn_rcov(:), rvdw(:, :), aes_dip_scale(:)
   integer :: isp, jsp, izp, jzp

   ! Setup coulomb interaction collection 
   tb_cn_rcov = p_cn_rcov(mol%num)
   call new_coulomb(coulomb, mol, error, cn_count_type=cn_count%erf, &
      & cn_rcov=tb_cn_rcov, cn_exp=p_cn_exp)
   if (allocated(error)) return

   ! Setup average
   allocate(average)
   call new_average(average, average_id%arithmetic)

   ! Setup damped multipole electrostatics
   allocate(aes2)
   call new_average(average, average_id%arithmetic)
   aes_dip_scale = p_aes_dip_scale(mol%num)
   allocate(rvdw(mol%nid, mol%nid)) 
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, isp
         jzp = mol%num(jsp)
         rvdw(jsp, isp) = get_vdw_rad(jzp, izp) &
            & * average%value(p_rvdw_scale(izp), p_rvdw_scale(jzp))
         rvdw(isp, jsp) = rvdw(jsp, isp)
      end do
   end do

   call new_gxtb_multipole(aes2, mol, rvdw, aes_dip_scale, p_aes_dmp_r2, &
      & p_aes_dmp_r3, p_aes_dmp_r4, p_aes_dmp_r5, p_aes_exp_r2, &
      & p_aes_exp_r3, p_aes_exp_r4, p_aes_exp_r5)
   call move_alloc(aes2, coulomb%aes2)

end subroutine make_multipole_gxtb


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(coulomb_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(coulomb_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

!> Return reference to container cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(coulomb_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(coulomb_cache)
      ptr => target
   end select
end subroutine view


subroutine test_generic(error, mol, qat, dpat, qpat, make_multipole, ref, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), contiguous, intent(in) :: qat(:)

   !> Atomic dipole moments for this structure
   real(wp), contiguous, intent(in) :: dpat(:, :)

   !> Atomic quadrupole moments for this structure
   real(wp), contiguous, intent(in) :: qpat(:, :)

   !> Factory to create new electrostatic objects
   procedure(multipole_maker) :: make_multipole

   !> Reference value to check against
   real(wp), intent(in) :: ref

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   type(tb_coulomb) :: coulomb
   type(container_cache) :: cache
   type(coulomb_cache), pointer :: ccache
   type(wavefunction_type) :: wfn
   real(wp) :: energy(mol%nat)
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   energy = 0.0_wp
   wfn%qat = reshape(qat, [size(qat), 1])
   wfn%dpat = reshape(dpat, [shape(dpat), 1])
   wfn%qpat = reshape(qpat, [shape(qpat), 1])
   call taint(cache, ccache)
   call ccache%update(mol)
   call make_multipole(coulomb, mol, error)
   if(allocated(error)) return
   call coulomb%update(mol, cache)
   call coulomb%aes2%get_energy(mol, cache, wfn, energy)

   call check(error, sum(energy), ref, thr=thr_)
   if (allocated(error)) then
      print *, ref, sum(energy), ref - sum(energy)
   end if

end subroutine test_generic


subroutine test_numpot(error, mol, qat, dpat, qpat, make_multipole, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), contiguous, intent(in) :: qat(:)

   !> Atomic dipole moments for this structure
   real(wp), contiguous, intent(in) :: dpat(:, :)

   !> Atomic quadrupole moments for this structure
   real(wp), contiguous, intent(in) :: qpat(:, :)

   !> Factory to create new electrostatic objects
   procedure(multipole_maker) :: make_multipole

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic
   type(tb_coulomb) :: coulomb
   type(container_cache) :: cache
   type(coulomb_cache), pointer :: ccache
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   real(wp) :: er(mol%nat), el(mol%nat)
   real(wp), allocatable :: vat(:), vdp(:, :), vqp(:, :)
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate(pot%vat(mol%nat, 1), pot%vsh(mol%nat, 1), pot%vao(mol%nat, 1), &
      & pot%vdp(3, mol%nat, 1), pot%vqp(6, mol%nat, 1), pot%kao(mol%nat, mol%nat, 1))
   allocate(vat(mol%nat), vdp(3, mol%nat), vqp(6, mol%nat))
   wfn%qat = reshape(qat, [size(qat), 1])
   wfn%dpat = reshape(dpat, [shape(dpat), 1])
   wfn%qpat = reshape(qpat, [shape(qpat), 1])
   call taint(cache, ccache)
   call ccache%update(mol)
   call make_multipole(coulomb, mol, error)
   if(allocated(error)) return
   call pot%reset

   call coulomb%update(mol, cache)
   call coulomb%aes2%get_potential(mol, cache, wfn, pot)

   do iat = 1, mol%nat
      er = 0.0_wp
      el = 0.0_wp
      wfn%qat(iat, 1) = wfn%qat(iat, 1) + step
      call coulomb%aes2%get_energy(mol, cache, wfn, er)
      wfn%qat(iat, 1) = wfn%qat(iat, 1) - 2*step
      call coulomb%aes2%get_energy(mol, cache, wfn, el)
      wfn%qat(iat, 1) = wfn%qat(iat, 1) + step
      vat(iat) = 0.5_wp*(sum(er) - sum(el))/step
   end do

   if (any(abs(pot%vat(:, 1) - vat) > thr_)) then
      call test_failed(error, "Charge-dependent potential does not match")
      write(*,*) "numerical potential:"
      print'(3es21.14)', vat
      write(*,*) "analytical potential:"
      print'(3es21.14)', pot%vat(: ,1)
      write(*,*) "difference:"
      print'(3es21.14)', pot%vat(: ,1) - vat
      return
   end if

   do iat = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         wfn%dpat(ic, iat, 1) = wfn%dpat(ic, iat, 1) + step
         call coulomb%aes2%get_energy(mol, cache, wfn, er)
         wfn%dpat(ic, iat, 1) = wfn%dpat(ic, iat, 1) - 2*step
         call coulomb%aes2%get_energy(mol, cache, wfn, el)
         wfn%dpat(ic, iat, 1) = wfn%dpat(ic, iat, 1) + step
         vdp(ic, iat) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   if (any(abs(pot%vdp(:, :, 1) - vdp) > thr_)) then
      call test_failed(error, "Dipole-dependent potential does not match")
      write(*,*) "numerical potential:"
      print'(3es21.14)', vdp
      write(*,*) "analytical potential:"
      print'(3es21.14)', pot%vdp
      write(*,*) "difference:"
      print'(3es21.14)', pot%vdp(:, :, 1) - vdp
      return
   end if

   do iat = 1, mol%nat
      do ic = 1, 6
         er = 0.0_wp
         el = 0.0_wp
         wfn%qpat(ic, iat, 1) = wfn%qpat(ic, iat, 1) + step
         call coulomb%aes2%get_energy(mol, cache, wfn, er)
         wfn%qpat(ic, iat, 1) = wfn%qpat(ic, iat, 1) - 2*step
         call coulomb%aes2%get_energy(mol, cache, wfn, el)
         wfn%qpat(ic, iat, 1) = wfn%qpat(ic, iat, 1) + step
         vqp(ic, iat) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   if (any(abs(pot%vqp(:, :, 1) - vqp) > thr_)) then
      call test_failed(error, "Quadrupole-dependent potential does not match")
      write(*,*) "numerical potential:"
      print'(3es21.14)', vqp
      write(*,*) "analytical potential:"
      print'(3es21.14)', pot%vqp
      write(*,*) "difference:"
      print'(3es21.14)', pot%vqp(:, :, 1) - vqp
      return
   end if

end subroutine test_numpot


subroutine test_numgrad(error, mol, qat, dpat, qpat, make_multipole, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), contiguous, intent(in) :: qat(:)

   !> Atomic dipole moments for this structure
   real(wp), contiguous, intent(in) :: dpat(:, :)

   !> Atomic quadrupole moments for this structure
   real(wp), contiguous, intent(in) :: qpat(:, :)

   !> Factory to create new electrostatic objects
   procedure(multipole_maker) :: make_multipole

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic
   type(tb_coulomb) :: coulomb
   type(container_cache) :: cache
   type(coulomb_cache), pointer :: ccache
   type(wavefunction_type) :: wfn
   real(wp) :: er(mol%nat), el(mol%nat), sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   wfn%qat = reshape(qat, [size(qat), 1])
   wfn%dpat = reshape(dpat, [shape(dpat), 1])
   wfn%qpat = reshape(qpat, [shape(qpat), 1])
   call taint(cache, ccache)
   call ccache%update(mol)
   call make_multipole(coulomb, mol, error)
   if(allocated(error)) return
   if (any(mol%periodic)) deallocate(coulomb%ncoord)

   do iat = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call ccache%update(mol)
         call coulomb%update(mol, cache)
         call coulomb%aes2%get_energy(mol, cache, wfn, er)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call ccache%update(mol)
         call coulomb%update(mol, cache)
         call coulomb%aes2%get_energy(mol, cache, wfn, el)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call ccache%update(mol)
   call coulomb%update(mol, cache)
   call coulomb%aes2%get_gradient(mol, cache, wfn, gradient, sigma)

   if (any(abs(gradient - numgrad) > thr_)) then
      call test_failed(error, "Gradient of energy does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', gradient
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', numgrad
      write(*,*) 'Difference:'
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_numgrad


subroutine test_numsigma(error, mol, qat, dpat, qpat, make_multipole, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), contiguous, intent(in) :: qat(:)

   !> Atomic dipole moments for this structure
   real(wp), contiguous, intent(in) :: dpat(:, :)

   !> Atomic quadrupole moments for this structure
   real(wp), contiguous, intent(in) :: qpat(:, :)

   !> Factory to create new electrostatic objects
   procedure(multipole_maker) :: make_multipole

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: ic, jc
   type(tb_coulomb) :: coulomb
   type(container_cache) :: cache
   type(coulomb_cache), pointer :: ccache
   type(wavefunction_type) :: wfn
   real(wp) :: er(mol%nat), el(mol%nat), sigma(3, 3), eps(3, 3), numsigma(3, 3)
   real(wp), allocatable :: gradient(:, :), xyz(:, :), lattice(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate(gradient(3, mol%nat), xyz(3, mol%nat))
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   wfn%qat = reshape(qat, [size(qat), 1])
   wfn%dpat = reshape(dpat, [shape(dpat), 1])
   wfn%qpat = reshape(qpat, [shape(qpat), 1])
   call taint(cache, ccache)
   call ccache%update(mol)
   call make_multipole(coulomb, mol, error)
   if(allocated(error)) return
   deallocate(coulomb%ncoord)

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   if (any(mol%periodic)) lattice = mol%lattice
   do ic = 1, 3
      do jc = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call ccache%update(mol)
         call coulomb%update(mol, cache)
         call coulomb%aes2%get_energy(mol, cache, wfn, er)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call ccache%update(mol)
         call coulomb%update(mol, cache)
         call coulomb%aes2%get_energy(mol, cache, wfn, el)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (allocated(lattice)) mol%lattice = lattice
         numsigma(jc, ic) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do
   numsigma = (numsigma + transpose(numsigma)) * 0.5_wp

   call ccache%update(mol)
   call coulomb%update(mol, cache)
   call coulomb%aes2%get_gradient(mol, cache, wfn, gradient, sigma)

   if (any(abs(sigma - numsigma) > thr_)) then
      call test_failed(error, "Strain derivatives do not match")
      write(*,*) 'Analytic sigma:'
      print'(3es21.14)', sigma
      write(*,*) 'Numeric sigma:'
      print'(3es21.14)', numsigma
      write(*,*) 'Difference:'
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_numsigma


subroutine test_e_aes_gfn2_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & 7.58819842725900E-1_wp,-8.68635970899368E-2_wp,-6.70762542651770E-1_wp, &
      &-1.27992265085229E-1_wp,-3.81264979031709E-1_wp, 1.80684122000111E-1_wp, &
      &-7.59652786607576E-2_wp,-4.78522650032789E-1_wp,-4.06862261737401E-1_wp, &
      & 3.18752726096977E-1_wp, 1.96659232729544E-1_wp,-1.54224745650353E-1_wp, &
      & 2.02418264533699E-1_wp, 4.96634503881606E-1_wp,-1.57256622208408E-1_wp, &
      & 3.85746250180502E-1_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      & 1.33576480531463E-1_wp,-2.34704592296522E-1_wp, 7.58166211639165E-2_wp, &
      & 1.31369672464812E-1_wp,-1.44919264929297E-2_wp,-6.34286998225388E-2_wp, &
      & 8.80079688209580E-3_wp,-9.48667597271333E-2_wp, 2.34462844725161E-2_wp, &
      &-9.06986024987756E-2_wp, 6.07865008837948E-2_wp,-5.99215264431395E-2_wp, &
      & 6.52955943309645E-2_wp, 4.03303831812593E-2_wp,-1.24827722880402E-1_wp, &
      &-9.22646680263864E-2_wp,-5.11755678142903E-2_wp,-5.33636860363439E-2_wp, &
      & 2.05026115065220E-2_wp, 4.12542711922320E-2_wp, 7.39929615141817E-2_wp, &
      & 7.90233076843428E-2_wp, 1.32084349195163E-2_wp,-4.42059556105020E-3_wp, &
      &-5.58279293106194E-2_wp,-1.63930555394780E-1_wp, 1.46549500224945E-1_wp, &
      & 4.76562297991379E-2_wp, 2.26979113448579E-2_wp,-2.56352017663062E-2_wp, &
      &-9.54894985148776E-3_wp, 2.28532337248730E-2_wp, 1.22246873138085E-1_wp, &
      & 6.74305363896158E-2_wp, 1.62691110966761E-1_wp, 6.91192523936835E-2_wp, &
      &-2.02780444228464E-1_wp,-2.76407251793201E-1_wp, 1.93484336903394E-1_wp, &
      &-4.01636414648176E-2_wp,-1.04814892924233E-1_wp, 1.73265414643709E-1_wp, &
      &-1.41626754632893E-1_wp, 1.45363019994928E-1_wp, 4.45314656877912E-3_wp, &
      &-2.22200250002484E-2_wp,-5.65142804921671E-2_wp,-6.34943004642575E-2_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      & 1.60685373895528E-1_wp,-3.76875594001669E-1_wp,-7.01783999199528E-1_wp, &
      & 4.91634045940829E-1_wp,-1.29989919109300E+0_wp, 5.41098625304000E-1_wp, &
      & 3.39388835592412E-2_wp,-1.22984974125562E-2_wp,-2.46376891719136E-2_wp, &
      &-5.09780038654557E-2_wp, 7.35078669697673E-3_wp,-9.30119438732752E-3_wp, &
      & 4.96274425923467E-2_wp, 7.47263609781933E-3_wp,-6.07305485759896E-2_wp, &
      & 4.21312628965061E-2_wp, 5.07302009717703E-2_wp, 1.11031059836425E-2_wp, &
      & 5.66468203038845E-3_wp,-1.83819423960093E-2_wp,-1.79160002262014E-2_wp, &
      & 4.26836801043380E-2_wp,-3.47136600519245E-2_wp, 1.22513181958127E-2_wp, &
      & 1.66237098815483E-3_wp,-1.23221078617715E-2_wp, 1.37287022404641E-2_wp, &
      & 8.77505169945025E-3_wp, 1.49173643908511E-2_wp,-1.53910732286182E-2_wp, &
      & 3.44584198267615E-2_wp, 4.77028193434869E-2_wp,-1.02095775834434E-2_wp, &
      & 3.99364926503320E-2_wp,-8.49395407815870E-3_wp,-2.42488422433182E-2_wp, &
      &-2.21353075815794E-2_wp, 1.88130230606430E-2_wp,-1.17206235259210E-4_wp, &
      & 1.38823690862043E-2_wp, 2.11531701659741E-2_wp, 2.22525138168387E-2_wp, &
      & 2.18908648962617E-2_wp, 2.25552935249801E-2_wp,-8.94647405917353E-3_wp, &
      &-3.04722997206538E-2_wp,-4.88052013444195E-3_wp,-1.29443908370892E-2_wp, &
      &-2.17668024013291E-1_wp, 7.02010342047839E-2_wp, 2.94951515530247E-1_wp, &
      &-5.85390403253958E-2_wp,-2.33818107592977E-1_wp,-7.72834915169562E-2_wp, &
      & 1.87552620722739E-1_wp, 5.90464884091131E-2_wp,-8.61994103787640E-2_wp, &
      &-4.64101615910787E-2_wp, 1.44452379133017E-3_wp,-1.01353210343975E-1_wp, &
      &-4.91677347567815E-2_wp, 1.29386733937149E-2_wp,-4.90324828178447E-3_wp, &
      &-2.38451923784347E-2_wp,-2.50679873439096E-4_wp, 5.40709830385660E-2_wp, &
      &-5.32562226964524E-2_wp, 1.09053183041917E-1_wp, 1.26844938114474E-1_wp, &
      &-4.96892366383249E-2_wp,-8.83023874945798E-2_wp,-7.35887154180222E-2_wp, &
      &-5.58921557991512E-1_wp, 6.51789405625317E-1_wp, 8.74949162077760E-1_wp, &
      &-6.21411138856650E-2_wp,-1.91057044554999E-1_wp,-3.16027604086248E-1_wp, &
      & 3.97918433455119E-1_wp, 1.67065427236830E-1_wp,-4.13032184704058E-1_wp, &
      & 2.80705220050717E-1_wp, 1.00281824404333E-1_wp, 1.51137512489375E-2_wp, &
      &-1.47514581652971E-1_wp, 4.24542092366428E-2_wp, 2.13980461856916E-1_wp, &
      & 1.60707682857530E-3_wp,-1.36988324556310E-1_wp,-6.64658802039458E-2_wp, &
      &-6.43374490513316E-1_wp, 7.40795233323125E-1_wp,-1.31732686037617E-1_wp, &
      &-1.15620551853270E+0_wp, 5.00370790867134E-1_wp, 7.75107176550934E-1_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, qat, dpat, qpat, make_multipole_gfn2, &
      & 1.1581304035711445E-2_wp)

end subroutine test_e_aes_gfn2_m01

subroutine test_e_aes_gfn2_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      &-1.68711131737618E-1_wp,-4.48538415935917E-1_wp, 1.26254102056646E-1_wp, &
      &-5.74045143966882E-1_wp, 3.79969821017488E-1_wp, 3.18528551778335E-1_wp, &
      &-1.51724744520024E-1_wp,-1.00445819399251E-1_wp, 4.89362856273315E-1_wp, &
      & 3.39809037472862E-1_wp, 2.35421594085126E-1_wp, 8.34540187873449E-1_wp, &
      &-4.96783950410221E-1_wp,-1.16050449947446E-1_wp,-1.11871005911778E-1_wp, &
      &-5.55715488728085E-1_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      &-2.09632500665695E-2_wp,-5.74601210929132E-2_wp,-1.14051499358462E-1_wp, &
      & 3.39982793143241E-1_wp,-1.60111036617454E-1_wp, 9.33173174274599E-2_wp, &
      &-1.67100063842922E-1_wp,-8.53438505465638E-2_wp,-1.83955088868127E-1_wp, &
      & 8.83056827990798E-2_wp,-6.32454175088047E-2_wp, 7.82604974598679E-2_wp, &
      & 3.21779694137560E-1_wp,-2.39336583970671E-1_wp, 2.55568763301003E-1_wp, &
      &-1.01276712912062E-2_wp,-5.92904371726635E-2_wp, 2.86383530855915E-2_wp, &
      & 3.17457819415172E-3_wp,-8.44412954452076E-2_wp, 1.16119097515870E-1_wp, &
      &-2.46513826370515E-2_wp, 7.07070993626731E-2_wp, 1.58819480451315E-1_wp, &
      &-2.91087994039460E-1_wp,-2.75530063526090E-1_wp,-6.01919158397077E-2_wp, &
      & 6.83469134649230E-2_wp, 2.18191369950846E-2_wp, 1.18228617563197E-2_wp, &
      & 5.37184499308458E-2_wp, 3.33505341589876E-1_wp, 1.38573916212116E-2_wp, &
      & 4.28947077958733E-1_wp,-1.00443418403241E-1_wp, 5.58642437169894E-1_wp, &
      & 1.79542559111173E-2_wp,-1.08744287416268E-2_wp, 1.35012076972140E-2_wp, &
      & 9.68016091220287E-2_wp,-8.92708119621574E-2_wp, 5.47845429661079E-2_wp, &
      & 3.73176478089200E-2_wp,-3.99071267437073E-2_wp,-4.05241248540806E-2_wp, &
      &-4.55400778065022E-2_wp,-8.52515845714788E-2_wp, 3.68564345900222E-1_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      &-1.51381335069229E-2_wp, 1.24599227476325E-2_wp,-8.43782018956266E-3_wp, &
      & 5.34699177321545E-2_wp, 2.50789791091612E-2_wp, 2.35759536964857E-2_wp, &
      &-4.33237636624836E-3_wp, 1.42282082219790E-2_wp, 1.89195986600509E-1_wp, &
      & 6.40630752093287E-2_wp, 7.17324341305498E-2_wp,-1.84863610234258E-1_wp, &
      &-2.84827574506977E-2_wp, 4.06630127109619E-1_wp, 3.23867272709838E-1_wp, &
      & 2.75716558566391E-1_wp,-4.16717188952446E-1_wp,-2.95384515259141E-1_wp, &
      &-5.01622557220743E-3_wp, 1.67965525359685E-2_wp,-2.48354185280116E-3_wp, &
      & 2.47562557716199E-2_wp,-3.15261176149967E-2_wp, 7.49976742500893E-3_wp, &
      & 3.24735522868949E-2_wp,-6.44123804759036E-1_wp, 9.89005820115999E-1_wp, &
      &-6.23704627178151E-2_wp,-5.53400810141925E-1_wp,-1.02147937240289E+0_wp, &
      &-5.14413522187525E-2_wp, 1.13318569048068E-3_wp, 1.07886831001909E-1_wp, &
      & 4.53573037606549E-3_wp,-6.33053742367396E-2_wp,-5.64454787831562E-2_wp, &
      &-2.92819023396152E-2_wp, 1.20585635707230E-2_wp, 1.11481468543835E-2_wp, &
      &-2.22493344477193E-2_wp,-3.33281460914114E-2_wp, 1.81337554852315E-2_wp, &
      & 1.02263502549609E-2_wp, 3.08679170954981E-3_wp,-2.98109958754787E-2_wp, &
      & 8.19907208675088E-3_wp, 2.48916478960544E-2_wp, 1.95846456205177E-2_wp, &
      & 8.73330263130388E-1_wp, 1.52460248675013E+0_wp,-1.21817595319730E-1_wp, &
      & 7.15313702126765E-1_wp, 2.20734329264357E-1_wp,-7.51512667810658E-1_wp, &
      & 1.30251939126722E-1_wp, 4.19166720807981E-2_wp,-5.15759761673018E-2_wp, &
      & 2.39946760617410E-2_wp,-2.39327564586919E-3_wp,-7.86759629594206E-2_wp, &
      & 4.19994503761355E-1_wp,-5.16164556208536E-2_wp, 3.69285564425158E-1_wp, &
      &-5.62275163431651E-2_wp, 2.92457165299565E-1_wp,-7.89280068186512E-1_wp, &
      & 6.14848279261306E-1_wp, 6.86966927211852E-2_wp,-1.04355228043612E+0_wp, &
      & 1.52035617405654E+0_wp,-2.70976038806502E-1_wp, 4.28704001174811E-1_wp, &
      &-2.41577641945364E-3_wp,-2.05245667804739E-3_wp, 3.57638830556812E-3_wp, &
      &-2.59500163165625E-3_wp,-2.02000418974092E-3_wp,-1.16061188611552E-3_wp, &
      & 2.55283427718014E-2_wp,-2.81941215788829E-3_wp,-2.31497442345624E-2_wp, &
      & 1.70710389135005E-2_wp,-5.66625357053174E-3_wp,-2.37859853723870E-3_wp, &
      & 1.29605658553924E-3_wp,-2.50340394706075E-3_wp, 1.45227175625805E-2_wp, &
      &-1.41399081671130E-2_wp, 3.56580399853426E-3_wp,-1.58187741481197E-2_wp, &
      &-5.02865562562276E-2_wp,-5.79945446793003E-2_wp, 2.72941505907548E-2_wp, &
      &-1.67690649973993E-1_wp,-1.87434204231407E-1_wp, 2.29924056654760E-2_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, qat, dpat, qpat, make_multipole_gfn2, &
      & -1.2918812895272689E-2_wp)

end subroutine test_e_aes_gfn2_m02


subroutine test_e_aes_gxtb_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      &  9.11522203271788E-1_wp, -7.82097836418159E-2_wp, -8.08504330144482E-1_wp, &
      & -1.37000889495809E-1_wp, -4.52267097738380E-1_wp,  2.36679420005181E-1_wp, &
      & -1.05648002487249E-1_wp, -7.54297465802927E-1_wp, -6.08441793086180E-1_wp, &
      &  3.25260352094806E-1_wp,  2.59020905785225E-1_wp, -3.24571359648248E-1_wp, &
      &  9.58852841103432E-2_wp,  7.10489048954234E-1_wp, -3.46833126525611E-1_wp, &
      &  1.07691663384632E+0_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      &  2.19676687117578E-1_wp, -2.38628805351101E-1_wp, -6.67926136598954E-2_wp, &
      &  9.08789134832573E-2_wp, -2.47404666933561E-2_wp, -3.80886784161556E-2_wp, &
      & -9.62541097237183E-5_wp, -2.05104062652084E-1_wp,  8.88293051440884E-2_wp, &
      & -7.04778680677458E-2_wp,  2.59319380611968E-2_wp, -4.17796806172496E-2_wp, &
      &  6.34285740786809E-2_wp,  3.27689719952741E-3_wp, -1.09374593302160E-1_wp, &
      & -6.81897640368189E-2_wp, -2.57584452593298E-2_wp, -3.99075013202982E-2_wp, &
      &  1.45900861513800E-2_wp,  4.27077899831544E-2_wp,  2.07619474694645E-2_wp, &
      &  8.25924891571033E-2_wp,  1.26574204510314E-2_wp, -2.16902323044048E-2_wp, &
      & -1.52965155178979E-1_wp, -3.41035990387539E-1_wp,  3.55742579151667E-1_wp, &
      &  4.05295361420255E-2_wp,  1.72033606861238E-2_wp, -1.94953016778706E-2_wp, &
      & -9.17912741528184E-3_wp,  2.47610296849307E-2_wp,  8.33890955235088E-2_wp, &
      &  2.52740581396424E-3_wp, -8.41056876326470E-2_wp, -2.01864567745985E-1_wp, &
      & -2.05833047812537E-1_wp, -4.12857817751401E-1_wp,  1.55593264465156E-1_wp, &
      & -8.77674152971623E-2_wp, -7.06126689541731E-2_wp,  1.96540235611247E-1_wp, &
      & -1.03520226673464E-1_wp,  4.18571029740440E-1_wp, -1.19912437872244E-1_wp, &
      &  2.25562683765269E-1_wp, -2.56052665302714E-1_wp,  2.81354902974319E-2_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      &  2.65296725923411E-1_wp, -5.27073655193522E-1_wp, -5.42486557496542E-1_wp, &
      &  8.28167421750275E-3_wp, -5.46317247544994E-1_wp,  2.77189831573131E-1_wp, &
      &  6.41217499868503E-2_wp, -7.46930900037705E-3_wp, -3.37792994905333E-2_wp, &
      & -1.01035770791397E-1_wp,  5.98601491898783E-3_wp, -3.03424504963170E-2_wp, &
      & -2.02521175901261E-1_wp, -8.31951636773058E-2_wp,  3.48280054488898E-1_wp, &
      & -2.37006624703918E-2_wp, -5.74808636528217E-2_wp, -1.45758878587639E-1_wp, &
      &  2.88277584629268E-2_wp, -5.67732142087247E-2_wp,  2.30062844067880E-3_wp, &
      &  7.98514811249413E-2_wp, -7.28375241620895E-2_wp, -3.11283869036054E-2_wp, &
      & -8.75989514290190E-2_wp, -3.67144378951517E-2_wp, -1.57995426441371E-2_wp, &
      & -1.74983251069399E-1_wp, -1.97780953556410E-2_wp,  1.03398494073156E-1_wp, &
      &  1.48842699421561E-2_wp,  5.03373213745870E-2_wp, -3.33839814996795E-3_wp, &
      &  5.00663019123519E-2_wp, -5.36581134615616E-3_wp, -1.15458717921884E-2_wp, &
      & -7.22359393334886E-2_wp,  8.21021420249764E-3_wp, -1.08978846742727E-2_wp, &
      &  4.30336105634189E-2_wp,  9.33634006465275E-2_wp,  8.31338240077615E-2_wp, &
      &  7.05017170184208E-2_wp,  1.30944393094117E-2_wp, -1.03071446122418E-1_wp, &
      &  5.42154377145829E-2_wp, -4.51145513179275E-2_wp,  3.25697291039990E-2_wp, &
      & -9.69015554152437E-1_wp, -2.42482895873653E-2_wp,  1.46649393275501E+0_wp, &
      & -4.08484244902551E-1_wp, -7.64311630529712E-1_wp, -4.97478378602577E-1_wp, &
      &  1.35077145658060E-1_wp,  5.19846064604624E-2_wp, -6.39768839713368E-2_wp, &
      & -4.55230060343506E-2_wp, -8.01566327759179E-3_wp, -7.11002616867236E-2_wp, &
      & -5.81705358803799E-2_wp,  1.91784706493377E-2_wp, -1.02478599946261E-2_wp, &
      & -1.08662100500335E-2_wp,  1.49974769259695E-3_wp,  6.84183958750060E-2_wp, &
      & -2.92302567753971E-1_wp,  4.13934393013674E-1_wp,  7.04691272247647E-1_wp, &
      & -1.00401329645679E-1_wp, -9.99802936598431E-2_wp, -4.12388704493676E-1_wp, &
      & -8.09789442486292E-1_wp,  1.03041814968889E+0_wp,  1.47726400715949E+0_wp, &
      & -1.34158594470813E-1_wp, -2.84513670676282E-1_wp, -6.67474564673192E-1_wp, &
      &  8.54128547055306E-2_wp,  2.02950688916942E-1_wp, -2.66812602258740E-1_wp, &
      &  3.69136740703624E-3_wp, -1.58863050432295E-1_wp,  1.81399747553209E-1_wp, &
      & -5.81724242258914E-1_wp,  2.03168912567440E-1_wp,  8.37135341886930E-1_wp, &
      & -1.22888546582942E-1_wp, -6.14167400315125E-1_wp, -2.55411099628017E-1_wp, &
      & -1.54540914767487E-1_wp,  5.04103315328362E-1_wp, -1.26606308897257E-1_wp, &
      &  7.07498138688553E-2_wp,  2.20126196204851E-1_wp,  2.81147223664743E-1_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, qat, dpat, qpat, make_multipole_gxtb, &
      & -2.362232239046687E-3_wp)

end subroutine test_e_aes_gxtb_m01

subroutine test_e_aes_gxtb_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & -1.09728877539416E-1_wp, -4.73502390915595E-1_wp, -1.70738958364902E-1_wp, &
      & -6.43349175946093E-1_wp,  8.96525477884417E-1_wp,  3.35126728257308E-1_wp, &
      & -6.65721051857640E-2_wp, -5.04759927966996E-2_wp,  5.14669970222235E-1_wp, &
      &  3.45527006432455E-1_wp,  8.36530865692038E-2_wp,  7.84009389499457E-1_wp, &
      & -6.23680788779179E-1_wp, -7.42958151399180E-2_wp, -6.46587171605577E-2_wp, &
      & -6.82508837469360E-1_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      & -5.72100924896129E-3_wp, -1.14889030429441E-2_wp, -7.91477885010296E-2_wp, &
      &  1.92738210001794E-1_wp, -1.64437984690317E-1_wp,  1.81151035308866E-1_wp, &
      & -2.88872625413932E-1_wp, -3.28203280711641E-1_wp, -2.97054491047612E-1_wp, &
      &  1.06777060279608E-1_wp, -5.93394942685946E-2_wp,  1.36123689767769E-1_wp, &
      & -1.12724584379820E-2_wp,  4.92898216997059E-1_wp,  1.25644046007025E-1_wp, &
      &  2.45475335437662E-3_wp, -2.33869757195913E-2_wp,  1.92252859488222E-2_wp, &
      &  6.64295098033643E-3_wp, -3.75119141201802E-2_wp,  8.79678277641348E-2_wp, &
      & -1.88683675778822E-2_wp,  4.89190607473817E-2_wp,  1.13099172228354E-1_wp, &
      & -8.05941344735003E-2_wp, -3.93689702054228E-2_wp, -3.73850296192678E-1_wp, &
      &  3.35464498222940E-2_wp,  8.97596631899762E-3_wp,  1.31975445678180E-2_wp, &
      &  6.04710609000185E-2_wp,  5.98675202420446E-1_wp,  3.01151053376063E-2_wp, &
      &  4.33810425395489E-2_wp, -5.79215468958236E-2_wp, -3.08013265886630E-2_wp, &
      & -4.10516704096329E-2_wp,  1.82317647520109E-2_wp, -3.21786100523874E-2_wp, &
      &  9.16050573410953E-2_wp, -8.61201612785467E-2_wp,  4.38021159178940E-2_wp, &
      & -7.61167516656987E-3_wp, -3.01626871343206E-2_wp, -4.46796759617521E-2_wp, &
      &  1.04974151235393E-1_wp, -1.32547295220027E-1_wp, -3.30458521872800E-2_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      & -5.09550484125263E-2_wp,  1.72809463697124E-2_wp, -5.35404589798327E-3_wp, &
      &  1.02541956545932E-1_wp,  6.89148732087383E-2_wp,  5.63090943105098E-2_wp, &
      &  4.13194060848729E-1_wp, -2.59439531758353E-1_wp,  3.77546338376341E-1_wp, &
      &  4.92340900547148E-1_wp,  6.34867362405169E-1_wp, -7.90740399225070E-1_wp, &
      &  4.36802951174498E-1_wp,  4.43809834278748E-1_wp, -6.30115658068373E-2_wp, &
      &  5.87182124690058E-1_wp,  1.35442897312958E-2_wp, -3.73791385367662E-1_wp, &
      & -6.06395843433418E-2_wp,  9.16714786287606E-2_wp, -3.49318265474903E-2_wp, &
      &  5.86787439911592E-2_wp, -4.30506012090299E-2_wp,  9.55714108908321E-2_wp, &
      & -4.93660118456911E-1_wp, -1.40583679764684E-2_wp,  1.15120137566704E+0_wp, &
      &  2.67763073761009E-1_wp, -2.82707135461787E-1_wp, -6.57541257210127E-1_wp, &
      & -5.79993303937843E-2_wp,  2.29588142067121E-2_wp,  1.13065045373809E-1_wp, &
      & -5.04590900735779E-3_wp, -7.07456012620653E-2_wp, -5.50657149800245E-2_wp, &
      & -5.94880195051712E-2_wp,  1.01939554244829E-2_wp,  3.59830721818075E-2_wp, &
      & -5.02609796953404E-2_wp, -7.24201080681576E-2_wp,  2.35049473233637E-2_wp, &
      &  1.43703024044761E-2_wp, -3.51914644941350E-3_wp, -5.93985927655152E-2_wp, &
      & -6.50156586718339E-3_wp,  7.83656826263669E-2_wp,  4.50282903610388E-2_wp, &
      &  1.87601660696053E-1_wp,  1.47565624486511E+0_wp, -3.66017021297139E-1_wp, &
      & -7.63361258071281E-2_wp, -3.23650800957690E-1_wp,  1.78415360601088E-1_wp, &
      &  1.33822417316290E-1_wp,  6.42003843696540E-2_wp, -4.89711205637757E-2_wp, &
      &  2.46314123034987E-2_wp,  3.56570152085797E-4_wp, -8.48512967525145E-2_wp, &
      &  5.12809808185255E-1_wp, -1.95377370907338E-1_wp,  7.67894608745372E-1_wp, &
      & -1.01071870467466E-1_wp,  4.84464513956034E-1_wp, -1.28070441693063E+0_wp, &
      &  2.52326782996651E-1_wp, -2.17916916655596E-1_wp,  3.69466829839597E-1_wp, &
      &  1.33708445504522E-1_wp,  2.01094755757356E-1_wp, -6.21793612836248E-1_wp, &
      & -5.54406365864200E-2_wp, -1.38471109322497E-2_wp,  7.55292411653041E-2_wp, &
      & -3.67483712055698E-2_wp, -5.12834488056551E-2_wp, -2.00886045788824E-2_wp, &
      &  4.50503755654972E-2_wp, -2.03301739459057E-2_wp, -4.00778925928376E-2_wp, &
      &  2.37680767403076E-2_wp, -3.19225993322992E-2_wp, -4.97248297265962E-3_wp, &
      &  9.50230209063980E-2_wp,  5.72664475143504E-2_wp, -1.65466923698361E-2_wp, &
      & -1.49706000292990E-2_wp, -7.04684149934698E-4_wp, -7.84763285365622E-2_wp, &
      & -2.31037117470676E-1_wp,  5.30708208278669E-1_wp, -3.11924373343842E-1_wp, &
      & -4.22549359625799E-1_wp, -6.93996919248957E-1_wp,  5.42961490814520E-1_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, qat, dpat, qpat, make_multipole_gxtb, &
      & -1.235061618042014E-3_wp)

end subroutine test_e_aes_gxtb_m02

subroutine test_e_aes_gxtb_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(4) = [&
      &  1.76167728179214E+0_wp, -5.86846576518797E-1_wp, -5.88294165118524E-1_wp, &
      & -5.86536540306923E-1_wp]
   real(wp), parameter :: dpat(3, 4) = reshape([&
      &  1.03102227915076E-1_wp,  1.70815542489309E-1_wp, -2.59639799151207E-1_wp, &
      & -1.87602572329602E-1_wp, -2.74737606941817E-1_wp, -1.49657158899311E-1_wp, &
      & -2.03472803217469E-1_wp,  2.74583536877613E-1_wp,  1.60701200271541E-1_wp, &
      &  3.20057753764996E-1_wp, -9.88144341108652E-2_wp,  1.55919293968284E-1_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 4) = reshape([&
      & -7.51267637005315E-1_wp,  3.31653256315932E-1_wp, -2.26838871462153E-1_wp, &
      & -9.32402326710803E-1_wp, -1.22767701815665E+0_wp,  9.78106508467468E-1_wp, &
      &  2.61814638396274E-3_wp, -2.16516052249151E-2_wp, -9.92760600427545E-3_wp, &
      & -2.24232173276705E-2_wp, -3.20385899497410E-2_wp,  7.30945962031626E-3_wp, &
      & -1.81977682871093E-3_wp,  3.29176108412411E-2_wp, -2.14330340960966E-2_wp, &
      &  5.73574812947845E-3_wp, -3.15197707756747E-2_wp,  2.32528109248094E-2_wp, &
      & -4.08481125421503E-2_wp,  1.54832741157508E-2_wp,  1.99270553528201E-2_wp, &
      & -2.74175062221744E-2_wp, -3.37992023179355E-3_wp,  2.09210571893266E-2_wp],&
      & shape(qpat))

   call get_structure(mol, "f-block", "CeCl3")
   call test_generic(error, mol, qat, dpat, qpat, make_multipole_gxtb, &
      & -1.142557089747410E-2_wp)

end subroutine test_e_aes_gxtb_cecl3


subroutine test_p_aes_gfn2_m07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      &-1.57324874193807E-1_wp, 1.65226165141353E-1_wp, 3.22324668351963E-1_wp, &
      & 3.63568850622050E-2_wp, 4.85675858799958E-2_wp,-3.59193103231093E-1_wp, &
      &-1.93844716126808E-1_wp,-3.86495668155725E-1_wp, 3.10105112155659E-1_wp, &
      & 8.34817252722515E-2_wp,-3.62673839838247E-1_wp, 3.64142535953715E-1_wp, &
      & 3.34643871790136E-1_wp,-4.69889612438212E-1_wp,-1.89224402240103E-1_wp, &
      & 4.53797666616760E-1_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      &-5.34076664252634E-2_wp,-7.07490405401130E-2_wp, 4.71123933536844E-3_wp, &
      &-3.62784359684118E-2_wp,-5.18758692875438E-2_wp, 7.93631053543726E-2_wp, &
      & 1.11938356444161E-1_wp, 1.56044569684314E-1_wp,-2.12215998018174E-1_wp, &
      &-3.92085646913333E-2_wp,-5.69893491841943E-2_wp, 1.28953955392941E-1_wp, &
      &-5.46548686385054E-2_wp,-6.46313930704017E-2_wp,-1.24929838203386E-1_wp, &
      &-2.33507749261263E-1_wp,-2.26409002121572E-1_wp, 2.06262816770210E-1_wp, &
      & 9.19783215014347E-2_wp, 2.54616916700422E-2_wp, 7.72507018523714E-2_wp, &
      &-1.39230892466345E-2_wp,-2.45837616220100E-2_wp, 8.53686058667152E-2_wp, &
      & 7.50097581489031E-3_wp, 2.62577733138995E-2_wp,-2.62783970963764E-1_wp, &
      &-7.99726651126221E-2_wp, 6.21043258370843E-2_wp, 8.31692921200217E-2_wp, &
      &-2.15463989266538E-1_wp, 1.34268047537320E-1_wp,-2.06897976969809E-1_wp, &
      &-3.96442967082123E-2_wp,-3.49442615996672E-2_wp,-2.63618691754457E-2_wp, &
      & 2.29820002870976E-2_wp, 6.60537345477610E-2_wp,-2.64213473449450E-2_wp, &
      &-1.86780038403075E-2_wp, 4.93641633812145E-2_wp,-9.89237975918140E-2_wp, &
      &-1.09709889162893E-1_wp,-1.65530934097439E-1_wp, 1.30980466379114E-1_wp, &
      & 1.37523724251365E-1_wp, 2.02435409108701E-1_wp, 1.60308459308605E-1_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      & 3.12010817257954E-1_wp,-5.72272801501363E-1_wp, 1.47090811311428E-2_wp, &
      &-3.52134738089932E-2_wp, 7.59119224271274E-3_wp,-3.26719898389098E-1_wp, &
      &-6.33288760707047E-2_wp, 3.31876114636840E-2_wp,-4.16806665423449E-2_wp, &
      &-3.00726289406259E-2_wp,-5.23087892909432E-2_wp, 1.05009542613049E-1_wp, &
      & 7.90626551753495E-2_wp,-5.59939019741557E-1_wp,-1.96613381144146E-1_wp, &
      &-1.90698681108742E-1_wp,-2.26666036370009E-1_wp, 1.17550725968797E-1_wp, &
      &-3.48658897797145E-2_wp, 3.97002159088856E-2_wp,-1.12930377647894E-2_wp, &
      &-2.32804989743644E-2_wp,-3.57057607813369E-2_wp, 4.61589275445038E-2_wp, &
      &-2.89670733084137E-2_wp, 4.01458351819485E-2_wp,-1.00338857618633E-2_wp, &
      & 2.30442398366667E-2_wp, 2.48446985131125E-2_wp, 3.90009590702768E-2_wp, &
      &-6.35919535338223E-2_wp, 4.28616775464788E-2_wp,-5.24795523285702E-2_wp, &
      &-4.74846575211650E-2_wp,-7.84863345937851E-2_wp, 1.16071505862394E-1_wp, &
      & 1.38813055815913E-2_wp, 3.67663646279253E-2_wp,-1.92557024583515E-2_wp, &
      & 6.55823370693129E-2_wp, 2.67104420466093E-2_wp, 5.37439687675907E-3_wp, &
      & 5.72363245854551E-2_wp,-9.17823355037531E-2_wp, 7.09074696257253E-3_wp, &
      &-3.36815808228537E-2_wp,-2.99818754969835E-2_wp,-6.43270715480298E-2_wp, &
      & 1.48546633546123E-2_wp, 3.78576983704521E-1_wp,-2.68648501666339E-2_wp, &
      & 2.54076030338518E-1_wp,-5.45865702181845E-2_wp, 1.20101868120234E-2_wp, &
      &-2.34554939533551E-2_wp,-4.71685480620286E-3_wp, 1.08088845822484E-2_wp, &
      &-2.95235812978304E-3_wp,-8.79348738658490E-3_wp, 1.26466093711069E-2_wp, &
      &-8.46372331176201E-2_wp,-4.53540864160115E-1_wp, 3.70587487482896E-1_wp, &
      & 3.56920411656286E-1_wp,-1.73479973145399E-1_wp,-2.85950254365274E-1_wp, &
      &-3.10036341823726E-2_wp, 4.58609916654156E-2_wp,-2.25258036224922E-2_wp, &
      & 7.15172736735734E-2_wp, 6.33957443767007E-2_wp, 5.35294378048649E-2_wp, &
      &-6.17710777868406E-2_wp, 6.50338226589586E-2_wp, 6.95799630181924E-2_wp, &
      &-3.81897973953993E-2_wp,-1.05412473958445E-1_wp,-7.80888523135179E-3_wp, &
      &-1.72750707129033E-2_wp, 2.58426297775894E-2_wp,-2.11922953328261E-3_wp, &
      & 1.08555732061051E-2_wp,-1.09741618008241E-2_wp, 1.93943002461841E-2_wp, &
      & 9.28747570418450E-4_wp, 3.77457812126494E-2_wp, 6.30415751637603E-3_wp, &
      &-1.49438136169359E-2_wp,-3.58731193456367E-2_wp,-7.23290508679160E-3_wp, &
      & 5.30696905643067E-2_wp,-8.24293499514649E-2_wp,-2.49519149355775E-1_wp, &
      &-7.45110987815727E-2_wp, 3.99113295183259E-1_wp, 1.96449458791466E-1_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "07")
   call test_numpot(error, mol, qat, dpat, qpat, make_multipole_gfn2, thr_in=thr1)

end subroutine test_p_aes_gfn2_m07

subroutine test_p_aes_gfn2_m08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      &-2.05667849633573E-1_wp,-3.99553782585639E-1_wp, 3.29243002578226E-1_wp, &
      &-3.11738017870329E-1_wp, 3.58862021565309E-2_wp, 3.21876529970816E-1_wp, &
      & 4.14749012704168E-2_wp, 2.95740260643397E-2_wp,-5.06348061187294E-1_wp, &
      & 3.43067917816204E-1_wp, 6.88373844005491E-1_wp, 7.03353178466825E-2_wp, &
      &-9.62419349868522E-2_wp,-1.32208182278391E-1_wp, 9.79071404471024E-2_wp, &
      &-3.05981053613728E-1_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      & 1.71354209686771E-1_wp,-5.59983867115148E-2_wp, 3.52751805328171E-2_wp, &
      &-8.07286634442154E-2_wp,-1.83810893296477E-2_wp, 1.14804020908522E-1_wp, &
      & 5.43571075658673E-1_wp,-6.30207096252022E-2_wp,-2.90104224098243E-1_wp, &
      &-6.02729537693486E-2_wp,-9.68718121551570E-2_wp, 5.67655999554924E-2_wp, &
      &-5.27791237942470E-2_wp,-6.13064550549191E-3_wp, 1.21704867523378E-1_wp, &
      & 2.69073848219728E-2_wp,-1.28803273925071E-1_wp, 2.35411279721039E-1_wp, &
      & 5.61846402437907E-2_wp, 9.89145696519802E-2_wp, 5.82766954213542E-2_wp, &
      & 6.46336803630342E-2_wp,-1.14523166471961E-1_wp, 1.34901332862348E-2_wp, &
      & 6.93956165710249E-3_wp, 2.80441545396097E-2_wp, 1.18436233762165E-2_wp, &
      & 2.51716438022498E-1_wp, 5.24823997035039E-2_wp,-2.29894641437262E-1_wp, &
      &-3.13902114644814E-2_wp, 3.44669877654266E-1_wp, 1.17062210876413E-1_wp, &
      &-2.14618673822334E-2_wp, 2.78345191993562E-2_wp,-2.37648333478748E-1_wp, &
      &-1.28803705368037E-1_wp,-8.41493196595313E-2_wp, 6.28432899666357E-2_wp, &
      &-2.42778543915217E-2_wp,-4.80798033618682E-2_wp,-1.86217780469857E-2_wp, &
      & 5.09565618678096E-2_wp,-2.80193475804829E-1_wp, 3.45120669963784E-3_wp, &
      &-6.66307827903392E-3_wp,-1.64671876783801E-2_wp, 3.60197289031932E-2_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      & 2.06461741736693E-1_wp, 7.28541667428188E-2_wp,-1.33835840983553E-1_wp, &
      &-6.95983307006687E-2_wp, 4.82075712059372E-2_wp,-7.26259007531408E-2_wp, &
      &-3.74564698563484E-2_wp,-2.03724977582733E-2_wp,-9.57439295106430E-2_wp, &
      &-2.83219382726871E-2_wp,-4.63569953734966E-2_wp, 1.33200399366992E-1_wp, &
      &-2.86835135385895E-1_wp,-1.60537799699520E-1_wp, 4.56277204936500E-1_wp, &
      & 3.14838988713061E-1_wp,-5.46513710350918E-1_wp,-1.69442069550605E-1_wp, &
      & 1.45266340994970E-2_wp, 3.69793732462622E-3_wp, 6.30127369127499E-4_wp, &
      & 7.60426093033395E-3_wp,-7.08695614481962E-3_wp,-1.51567614686242E-2_wp, &
      &-2.17297304116788E-2_wp, 1.02011244860144E-2_wp,-5.55372569258381E-2_wp, &
      &-8.74911034050772E-4_wp, 6.41460422550064E-3_wp, 7.72669873375170E-2_wp, &
      &-1.72479380870008E+0_wp,-9.83892436920738E-1_wp, 4.59174454545800E-1_wp, &
      & 2.08601212023478E-1_wp, 5.21400954414695E-2_wp, 1.26561935415428E+0_wp, &
      &-1.52836481375808E-2_wp, 6.16726765036988E-3_wp, 1.88428715947021E-2_wp, &
      & 1.41119893526534E-2_wp, 2.53574127410936E-2_wp,-3.55922345712144E-3_wp, &
      &-2.16443510745452E-2_wp,-4.74176382224577E-3_wp, 2.11465971221881E-2_wp, &
      & 1.37203264366346E-2_wp,-8.14837599748600E-3_wp, 4.97753952357263E-4_wp, &
      & 7.43773200889999E-3_wp,-3.57910400477298E-2_wp, 9.09088526507622E-3_wp, &
      & 7.36392117765244E-2_wp,-4.16236233710783E-2_wp,-1.65286172739764E-2_wp, &
      &-3.42134325902312E-2_wp,-7.61060502488899E-1_wp,-2.45488811046729E-1_wp, &
      & 5.47673026793362E-1_wp,-3.36894101675850E-1_wp, 2.79702243636961E-1_wp, &
      &-5.21019443867959E-1_wp,-5.78752381345166E-1_wp, 4.79453134960792E-2_wp, &
      &-4.45615097566108E-1_wp, 5.28740776457770E-1_wp, 4.73074130371880E-1_wp, &
      &-1.21614048807072E-1_wp, 7.25837353275695E-3_wp,-3.41565915732748E-1_wp, &
      &-8.16092467802712E-2_wp,-1.28488782348307E-1_wp, 4.63179964539819E-1_wp, &
      &-1.00288109979275E-2_wp, 2.80676999575447E-2_wp, 8.84900791116329E-3_wp, &
      &-3.01420755561142E-2_wp, 3.45178585866044E-2_wp, 1.17980308676424E-3_wp, &
      & 2.85902815666733E-2_wp, 1.72947238727397E-2_wp, 7.70556210678494E-3_wp, &
      & 1.44970209246944E-2_wp,-1.22912847635515E-2_wp,-3.62958436734579E-2_wp, &
      & 1.73148533024057E-1_wp, 2.51424980765172E-1_wp, 5.56773723377576E-2_wp, &
      & 6.10240887950364E-1_wp, 9.08078399932844E-1_wp,-2.28825905361813E-1_wp, &
      & 2.10725818094498E-2_wp, 5.40916928372065E-3_wp,-1.48423865948397E-2_wp, &
      &-4.61098763513862E-3_wp, 2.64280891570932E-2_wp,-6.23019521460990E-3_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "08")
   call test_numpot(error, mol, qat, dpat, qpat, make_multipole_gfn2, thr_in=thr1)

end subroutine test_p_aes_gfn2_m08


subroutine test_p_aes_gxtb_m07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & -3.11013868908128E-1_wp,  2.17279822694310E-1_wp,  3.78045010818655E-1_wp, &
      &  1.06152524717814E-1_wp,  1.17752833761274E-1_wp, -4.36996222849876E-1_wp, &
      & -3.73978892346850E-1_wp, -6.56877699379819E-1_wp,  4.60878146336219E-1_wp, &
      &  1.09768991207884E-1_wp, -4.78208859838842E-1_wp,  3.67220607553746E-1_wp, &
      &  3.47529247713129E-1_wp, -5.85051124309072E-1_wp, -2.81700734745851E-1_wp, &
      &  1.01920021708489E+0_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      &  9.02924681578633E-2_wp,  1.19107308187557E-1_wp, -2.13187913176873E-2_wp, &
      & -2.66556543750565E-2_wp, -3.83271922464750E-2_wp,  5.72378946320154E-2_wp, &
      &  7.19622953150420E-2_wp,  1.02976685729981E-1_wp, -1.90036111981539E-1_wp, &
      & -2.07800921567854E-2_wp, -3.07077352516037E-2_wp,  7.86810804710750E-2_wp, &
      & -3.15657484387209E-2_wp, -3.69623526997453E-2_wp, -7.94304658163339E-2_wp, &
      & -3.69962922805186E-2_wp, -3.48239181138134E-2_wp, -2.90204962219600E-2_wp, &
      &  1.33668685256546E-1_wp,  8.37625598380338E-2_wp,  1.32919855369916E-1_wp, &
      & -1.79830676347187E-1_wp, -2.30941283858552E-1_wp, -1.55684401800542E-1_wp, &
      &  1.17159889798300E-1_wp,  2.43252885644797E-2_wp, -5.10545026285145E-1_wp, &
      & -5.93238772916325E-2_wp,  4.04204801352332E-2_wp,  7.67094365525818E-2_wp, &
      &  1.57915637958136E-3_wp,  2.24145576108900E-1_wp,  7.29297765912292E-2_wp, &
      & -2.64355235821374E-2_wp, -2.16987092922138E-2_wp, -1.91293967914836E-2_wp, &
      &  1.10550208582573E-2_wp,  4.03761057137495E-2_wp, -1.64163937423298E-2_wp, &
      & -4.13194952881613E-2_wp, -8.57538939468032E-3_wp, -1.33957906033093E-1_wp, &
      & -1.05896684019073E-1_wp, -1.81952103131200E-1_wp,  1.23317919177764E-1_wp, &
      & -9.69051024639739E-2_wp,  2.04717513917997E-1_wp, -2.13674160628888E-1_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      &  3.55452843145141E-1_wp, -3.11348810689210E-1_wp,  2.03906065029034E-1_wp, &
      & -7.47554439327699E-2_wp, -4.04465241816248E-2_wp, -5.59358908174175E-1_wp, &
      & -5.34861282290284E-2_wp,  2.97113944947855E-2_wp, -3.42628176338874E-2_wp, &
      & -2.89261638129338E-2_wp, -4.68594509119050E-2_wp,  8.77489458629158E-2_wp, &
      & -7.72041815480637E-2_wp, -1.47501151370978E-1_wp, -1.51816000184381E-1_wp, &
      & -8.34881975266728E-2_wp, -1.16480527313187E-1_wp,  2.29020181732445E-1_wp, &
      & -5.55526886007858E-2_wp,  5.15010580235321E-2_wp, -2.41083122186270E-2_wp, &
      & -3.99148687671438E-2_wp, -6.01407430080340E-2_wp,  7.96610008194130E-2_wp, &
      & -4.50809327407917E-2_wp,  5.57198210449159E-2_wp, -1.99097626354447E-2_wp, &
      &  4.16854551787123E-2_wp,  4.63098415504627E-2_wp,  6.49906953762365E-2_wp, &
      & -3.47205773309090E-1_wp,  4.16445376529872E-1_wp, -2.53077009629518E-1_wp, &
      & -1.09656069729926E-1_wp, -2.47852783400634E-1_wp,  6.00282782938606E-1_wp, &
      &  9.43481981625274E-2_wp,  8.66098942527645E-2_wp, -1.45350626735982E-1_wp, &
      &  3.01824231278657E-1_wp,  7.69349293175310E-2_wp,  5.10024285734536E-2_wp, &
      & -1.82928405597298E-1_wp,  2.22006324198212E-1_wp, -9.59444588221334E-2_wp, &
      &  1.62795887781162E-1_wp,  1.95564375635441E-1_wp,  2.78872864419435E-1_wp, &
      & -2.53169025688482E-2_wp,  6.32972021850510E-1_wp,  1.53900878858293E-1_wp, &
      &  5.79152558950814E-1_wp, -4.81119531335919E-1_wp, -1.28583976289446E-1_wp, &
      & -1.45288097202633E-2_wp, -1.12468707154057E-2_wp,  3.01680425274717E-2_wp, &
      &  1.16486559866170E-2_wp, -2.85011121204970E-2_wp, -1.56392328072086E-2_wp, &
      & -1.22218777061942E-2_wp, -1.35187830527870E+0_wp,  8.37956772433868E-1_wp, &
      &  4.38831795772269E-1_wp, -7.30950854426522E-1_wp, -8.25734894727677E-1_wp, &
      & -2.74918352113803E-2_wp,  5.63466484550578E-2_wp, -6.66582221375100E-3_wp, &
      &  6.93871904667963E-2_wp,  6.03877329557861E-2_wp,  3.41576574251313E-2_wp, &
      & -5.80957532242772E-2_wp,  6.96599591632022E-2_wp,  7.52673607975716E-2_wp, &
      & -3.47127231886956E-2_wp, -9.70374909045714E-2_wp, -1.71716075732944E-2_wp, &
      &  1.65858283563898E-2_wp,  1.80181972439994E-1_wp,  7.48158974709501E-3_wp, &
      &  6.54371951693999E-2_wp,  5.99274899622614E-2_wp, -2.40674181034830E-2_wp, &
      & -6.73155655140176E-2_wp,  1.87544738938634E-1_wp,  1.20484336791164E-1_wp, &
      & -1.48449064379718E-1_wp, -2.58368229827748E-1_wp, -5.31687712771465E-2_wp, &
      & -7.95135455721368E-1_wp, -7.46807592681798E-1_wp,  3.64945971144273E-1_wp, &
      & -5.36594622000848E-1_wp,  1.12612056564190E-1_wp,  4.30189484577094E-1_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "07")
   call test_numpot(error, mol, qat, dpat, qpat, make_multipole_gxtb, thr_in=thr1)

end subroutine test_p_aes_gxtb_m07

subroutine test_p_aes_gxtb_m08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & -3.79780215699870E-1_wp, -5.43777954602357E-1_wp,  4.47153291016247E-1_wp, &
      & -4.00575300169204E-1_wp,  1.22242636026864E-1_wp,  5.41648259537217E-1_wp, &
      &  1.10996311243541E-1_wp,  1.05126121304397E-1_wp, -6.61141370442608E-1_wp, &
      &  4.71849017393574E-1_wp,  8.51393966111946E-1_wp,  4.97683684943671E-2_wp, &
      & -9.84080662874769E-2_wp, -1.41392484648204E-1_wp, -2.25459831072563E-2_wp, &
      & -4.52556596719091E-1_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      &  3.49706701357985E-1_wp,  1.90764397811718E-3_wp, -2.47956693357178E-1_wp, &
      & -2.59945173786500E-1_wp,  1.53320621658160E-2_wp,  2.34715057483132E-1_wp, &
      &  6.06638941497731E-1_wp,  6.91813614479975E-2_wp, -1.33358654566062E-1_wp, &
      & -6.18427132837973E-2_wp, -1.13155740813063E-1_wp,  5.28959349236889E-2_wp, &
      & -4.00959198640804E-2_wp,  1.17665398630840E-3_wp,  7.30444712184876E-2_wp, &
      &  8.66198264870558E-2_wp, -2.41553541727722E-1_wp,  7.78580339839647E-1_wp, &
      &  2.63215526218328E-2_wp,  7.41933775154638E-2_wp,  3.68503370810625E-2_wp, &
      &  3.31990665029995E-2_wp, -7.95296106886445E-2_wp,  8.34178011027989E-3_wp, &
      & -5.15539660644317E-2_wp,  1.53014136031747E-1_wp,  1.44445403783895E-1_wp, &
      &  4.17131872407917E-1_wp,  1.54779488478117E-1_wp, -2.96869373889402E-1_wp, &
      &  7.40272149182254E-2_wp,  4.09121770754506E-2_wp, -5.09583488723214E-2_wp, &
      & -2.94530478271492E-1_wp,  8.05737709462253E-2_wp, -5.31356895194351E-1_wp, &
      & -1.01945378414306E-1_wp, -5.72106414972446E-2_wp,  4.98967028060630E-2_wp, &
      & -4.51701050358509E-4_wp, -3.31782149239246E-2_wp,  7.35252504912148E-3_wp, &
      &  2.69369681546772E-1_wp, -5.07799177236883E-1_wp,  1.39238712300427E-1_wp, &
      & -1.52645750268585E-2_wp,  3.89672088985521E-2_wp, -1.03888804047019E-2_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      &  5.12371694908371E-1_wp,  1.64247461756723E-1_wp, -5.40462850295631E-2_wp, &
      & -4.72218426310372E-2_wp, -1.58069944289599E-1_wp, -4.58325409878807E-1_wp, &
      & -1.36627378519990E-1_wp, -9.62919849763605E-2_wp, -4.31297534739556E-1_wp, &
      &  1.51774477810966E-2_wp, -2.49495992244568E-1_wp,  5.67924913259545E-1_wp, &
      & -3.24181001796519E-1_wp, -7.22222841647461E-2_wp,  4.51959615994358E-1_wp, &
      & -9.35111888674837E-2_wp, -1.12569393474457E+0_wp, -1.27778614197837E-1_wp, &
      & -5.84028952602633E-2_wp,  1.60711189591921E-1_wp,  1.50343297212535E-1_wp, &
      & -3.59209750116391E-2_wp, -1.15053227816128E-1_wp, -9.19404019522734E-2_wp, &
      & -2.66333320547464E-2_wp,  7.93450033307377E-3_wp, -6.67713261905338E-2_wp, &
      & -6.33215054906878E-3_wp, -4.90837596474728E-3_wp,  9.34046582452802E-2_wp, &
      & -3.90903678225012E-1_wp, -3.90541327399119E-1_wp,  5.33334965411967E-1_wp, &
      &  5.12431760634557E-1_wp,  1.40345014644446E-1_wp, -1.42431287186957E-1_wp, &
      & -1.71136563650398E-3_wp,  2.08816810118458E-2_wp,  1.80800010684163E-2_wp, &
      &  1.61898637588993E-2_wp,  2.86069497160145E-2_wp, -1.63686354319124E-2_wp, &
      & -5.22356545297120E-3_wp, -3.01616523562512E-2_wp,  3.51865276778731E-2_wp, &
      &  1.01514077269351E-2_wp, -1.27480823076761E-2_wp, -2.99629622249019E-2_wp, &
      & -1.45996560396537E-1_wp,  4.48319485628847E-2_wp,  4.99177071528605E-2_wp, &
      & -2.05966541191756E-1_wp,  1.33987609296270E-1_wp,  9.60788532436760E-2_wp, &
      & -3.54195857359102E-1_wp, -8.18273948831220E-1_wp, -5.07276983527944E-2_wp, &
      &  4.40150359584884E-1_wp,  1.42680074473957E-1_wp,  4.04923555711894E-1_wp, &
      & -1.44553630126509E-1_wp, -1.15568619322583E-1_wp, -6.40721632835537E-2_wp, &
      & -1.31847113734340E-1_wp,  8.65381118944525E-2_wp,  2.08625793410063E-1_wp, &
      &  1.89396196098601E-1_wp,  4.49960606536607E-2_wp, -1.12104665347949E+0_wp, &
      & -3.29938727939534E-1_wp, -3.78955847141035E-1_wp,  9.31650457380885E-1_wp, &
      & -3.47073915484226E-2_wp,  5.65334320589241E-2_wp,  4.05973879029644E-2_wp, &
      & -4.67455811743250E-2_wp,  5.60056241332738E-2_wp, -5.88999635454179E-3_wp, &
      & -6.29969540748327E-2_wp,  2.04136927308526E-2_wp, -3.89849024245117E-2_wp, &
      & -6.86321408342724E-2_wp,  2.29079655912482E-2_wp,  1.01981856499344E-1_wp, &
      & -9.66570310024295E-2_wp,  7.82331211868712E-2_wp,  2.10186483117440E-1_wp, &
      &  6.84802142903622E-1_wp,  7.39081486282856E-1_wp, -1.13529452115011E-1_wp, &
      & -4.43366464957071E-2_wp, -1.05564780591803E-1_wp,  4.08830515436467E-2_wp, &
      & -2.16029639977796E-2_wp,  8.90913392741227E-3_wp,  3.45359495206132E-3_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "08")
   call test_numpot(error, mol, qat, dpat, qpat, make_multipole_gxtb, thr_in=thr1)

end subroutine test_p_aes_gxtb_m08

subroutine test_p_aes_gxtb_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(4) = [&
      &  1.76167728179214E+0_wp, -5.86846576518797E-1_wp, -5.88294165118524E-1_wp, &
      & -5.86536540306923E-1_wp]
   real(wp), parameter :: dpat(3, 4) = reshape([&
      &  1.03102227915076E-1_wp,  1.70815542489309E-1_wp, -2.59639799151207E-1_wp, &
      & -1.87602572329602E-1_wp, -2.74737606941817E-1_wp, -1.49657158899311E-1_wp, &
      & -2.03472803217469E-1_wp,  2.74583536877613E-1_wp,  1.60701200271541E-1_wp, &
      &  3.20057753764996E-1_wp, -9.88144341108652E-2_wp,  1.55919293968284E-1_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 4) = reshape([&
      & -7.51267637005315E-1_wp,  3.31653256315932E-1_wp, -2.26838871462153E-1_wp, &
      & -9.32402326710803E-1_wp, -1.22767701815665E+0_wp,  9.78106508467468E-1_wp, &
      &  2.61814638396274E-3_wp, -2.16516052249151E-2_wp, -9.92760600427545E-3_wp, &
      & -2.24232173276705E-2_wp, -3.20385899497410E-2_wp,  7.30945962031626E-3_wp, &
      & -1.81977682871093E-3_wp,  3.29176108412411E-2_wp, -2.14330340960966E-2_wp, &
      &  5.73574812947845E-3_wp, -3.15197707756747E-2_wp,  2.32528109248094E-2_wp, &
      & -4.08481125421503E-2_wp,  1.54832741157508E-2_wp,  1.99270553528201E-2_wp, &
      & -2.74175062221744E-2_wp, -3.37992023179355E-3_wp,  2.09210571893266E-2_wp],&
      & shape(qpat))

   call get_structure(mol, "f-block", "CeCl3")
   call test_numpot(error, mol, qat, dpat, qpat, make_multipole_gxtb, thr_in=thr1)

end subroutine test_p_aes_gxtb_cecl3


subroutine test_g_aes_gfn2_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & 3.91783754152349E-3_wp,-8.39000740130402E-1_wp, 1.69682418230267E-2_wp, &
      & 5.91743613403359E-1_wp, 4.09885570880249E-1_wp, 1.60227715582643E-1_wp, &
      &-2.65240763854753E-1_wp,-1.71713502599602E-3_wp, 3.76177400795062E-2_wp, &
      &-4.24669774348641E-4_wp,-5.26879406504080E-1_wp,-3.99087528731630E-1_wp, &
      &-2.54713398335190E-1_wp, 2.33569743280495E-1_wp, 7.82966864852826E-1_wp, &
      & 5.01663149127674E-2_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      &-1.15230315070558E-2_wp, 1.15809130093675E-1_wp,-8.38709243085290E-2_wp, &
      &-2.40842660369881E-2_wp,-6.58555907012431E-2_wp, 1.03542140714310E-1_wp, &
      &-1.90775939412024E-2_wp,-4.79686231986468E-2_wp, 1.02101958055925E-1_wp, &
      & 3.11704676333466E-1_wp, 3.67757428498038E-1_wp,-2.68626030351259E-1_wp, &
      &-2.45303173973439E-1_wp, 2.15168641748650E-1_wp,-2.88169192295142E-1_wp, &
      & 7.46076101279135E-2_wp, 1.01260850934685E-1_wp,-6.59712753367644E-2_wp, &
      & 5.44645627050000E-2_wp,-1.33120135078562E-1_wp,-8.43224527361052E-2_wp, &
      & 3.37208023534908E-2_wp,-1.26192546432335E-1_wp,-1.64536604069415E-2_wp, &
      &-7.47507639838809E-2_wp, 8.15901291652326E-2_wp, 4.84083420254223E-3_wp, &
      & 4.04966160003665E-2_wp,-6.42622155152806E-2_wp,-1.37876178539178E-1_wp, &
      &-3.49009950008589E-2_wp,-1.64999060694037E-3_wp,-1.87784938817728E-2_wp, &
      & 3.38546931980014E-2_wp,-1.99801737922777E-1_wp, 1.21276453010155E-1_wp, &
      & 4.02492772972945E-2_wp, 9.30412809437698E-2_wp, 1.15656280124059E-1_wp, &
      &-1.60471833617121E-2_wp,-3.61079604858184E-2_wp, 5.12703027843782E-2_wp, &
      &-4.20241072314464E-1_wp, 2.40785717844855E-1_wp,-4.50917286857342E-1_wp, &
      & 8.43990446029860E-2_wp, 6.00782470724405E-2_wp, 3.73648402390034E-2_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      & 3.50526349479599E-1_wp, 7.66559042784789E-2_wp, 6.29436173506202E-2_wp, &
      & 1.11124071277413E-1_wp, 7.69742389175388E-2_wp,-4.13469966830219E-1_wp, &
      &-2.75381290899550E-2_wp, 2.72140206636799E-2_wp,-1.14813764067649E-3_wp, &
      &-3.11816795304511E-2_wp,-5.65659880409654E-2_wp, 2.86862667306302E-2_wp, &
      &-4.12709051889617E-2_wp, 1.53683666184660E-2_wp,-1.06041720944760E-2_wp, &
      &-1.83787750456468E-2_wp,-8.33865522179942E-2_wp, 5.18750772834380E-2_wp, &
      &-1.48976685337176E+0_wp, 1.11162195553861E+0_wp, 1.24124790055020E+0_wp, &
      &-7.07035872970699E-1_wp,-8.94732166949645E-1_wp, 2.48518952821565E-1_wp, &
      &-4.62949255267837E-1_wp, 1.99898431180569E-1_wp, 1.35742321098480E-1_wp, &
      &-2.22933568410016E-1_wp, 3.88208326804613E-1_wp, 3.27206934169355E-1_wp, &
      &-4.45229595152185E-1_wp, 3.50705916422203E-1_wp,-1.03664744506441E+0_wp, &
      &-2.15004289804945E+0_wp, 6.34956704177776E-2_wp, 1.48187704021659E+0_wp, &
      & 2.39478118929564E-1_wp, 1.26626807037684E-1_wp,-8.01809526788819E-2_wp, &
      & 1.36276853411382E-1_wp,-2.87966714802874E-2_wp,-1.59297166250682E-1_wp, &
      &-4.60959432062878E-2_wp,-3.97855259662190E-2_wp, 5.81094544183368E-2_wp, &
      &-7.39210019727102E-3_wp,-1.22487758651412E-2_wp,-1.20135112120492E-2_wp, &
      & 3.27707816948706E-2_wp,-8.29825872707352E-2_wp, 3.63750158737722E-2_wp, &
      &-2.13915374389003E-3_wp,-6.48618307640439E-3_wp,-6.91457975686428E-2_wp, &
      &-6.05850234707494E-2_wp,-2.94714753333237E-2_wp, 5.26564201125322E-2_wp, &
      &-4.18830617530320E-3_wp,-1.11315020843992E-2_wp, 7.92860335821687E-3_wp, &
      &-8.10266863560071E-3_wp,-7.97388103984468E-4_wp, 6.51328083982761E-3_wp, &
      &-6.81436512005934E-3_wp,-3.09292341837090E-3_wp, 1.58938779577390E-3_wp, &
      &-1.73587214521756E-1_wp,-2.36522554009785E-1_wp, 3.39223911104960E-1_wp, &
      &-1.51417740117755E-3_wp,-3.56156837783491E-2_wp,-1.65636696583206E-1_wp, &
      & 3.17626683737956E-2_wp,-4.20518147311984E-2_wp,-9.87840768707964E-3_wp, &
      &-5.80847745566733E-2_wp,-7.86844927166163E-2_wp,-2.18842606867161E-2_wp, &
      &-9.33977563713688E-2_wp, 2.43147293048641E-2_wp,-1.03488400945705E-2_wp, &
      &-4.21041056943522E-2_wp,-1.30729562759062E-1_wp, 1.03746596465939E-1_wp, &
      & 5.37212160477299E-1_wp, 5.98187241115427E-2_wp,-2.28495276491597E-1_wp, &
      & 7.13146767579691E-1_wp,-6.97470281887052E-1_wp,-3.08716883985702E-1_wp, &
      & 7.18967989352684E-2_wp, 7.07364711615708E-2_wp,-1.69845236886560E-2_wp, &
      & 2.85053043764643E-2_wp, 2.10540135304005E-2_wp,-5.49122752466128E-2_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "03")
   call test_numgrad(error, mol, qat, dpat, qpat, make_multipole_gfn2, thr_in=thr1)

end subroutine test_g_aes_gfn2_m03

subroutine test_g_aes_gfn2_m04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      &-8.99891426876454E-2_wp, 9.42167151049675E-2_wp,-1.49390480695042E-1_wp, &
      &-2.99114077043291E-1_wp, 4.85528544756154E-1_wp,-6.83158233623334E-2_wp, &
      & 1.49989659279676E-2_wp, 2.79361802448324E-1_wp,-1.24070454865924E-1_wp, &
      &-9.36742767297382E-2_wp,-2.19062009156706E-1_wp, 2.14544824574956E-1_wp, &
      & 3.06160853875672E-1_wp,-3.86107183037218E-1_wp,-1.51404848137219E-3_wp, &
      & 3.64257893712230E-2_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      &-3.64731729820407E-2_wp,-2.13686275288388E-2_wp, 6.68806426144540E-2_wp, &
      & 4.04240419990876E-2_wp,-1.49843414114636E-1_wp,-5.58836609439025E-3_wp, &
      &-6.05242232810174E-2_wp, 9.15234565804751E-2_wp,-7.40096212350601E-2_wp, &
      &-2.00610828069480E-2_wp, 8.79914756075611E-2_wp, 1.17474526864424E-1_wp, &
      &-1.23227254106242E-1_wp, 5.17071359855638E-2_wp,-3.45193834725206E-1_wp, &
      & 8.38659891802097E-2_wp,-1.85766281220651E-3_wp, 2.16613946001174E-3_wp, &
      & 7.28212028877821E-2_wp,-1.11656708078832E-1_wp, 4.44987301144173E-2_wp, &
      & 5.21199968878805E-1_wp, 2.42009354648746E-1_wp,-2.59373648270648E-1_wp, &
      & 2.42771833481858E-3_wp, 1.33481324585323E-1_wp, 4.70393507351052E-2_wp, &
      & 6.59447581952794E-2_wp,-4.77190059636476E-2_wp, 1.41459747150361E-1_wp, &
      &-5.74635892974854E-3_wp,-9.43452985224502E-2_wp,-4.14652054192348E-2_wp, &
      & 1.18721461568371E-1_wp,-8.07370651006255E-2_wp, 1.47715903211519E-1_wp, &
      & 1.80051594470901E-1_wp, 2.31074137154501E-1_wp, 4.61275280422502E-1_wp, &
      & 7.00280160506826E-2_wp,-6.56712228336652E-2_wp, 5.83187156889666E-2_wp, &
      &-7.73514384205769E-2_wp,-6.11970211951996E-3_wp,-1.18441889940799E-1_wp, &
      &-2.63615278011550E-2_wp,-8.98662961570448E-2_wp, 9.70285404014370E-2_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      & 2.07660448410257E-2_wp, 4.77562367062922E-2_wp,-6.38136462827382E-2_wp, &
      &-6.08076485231580E-2_wp,-2.89106396726373E-3_wp, 4.30476014417125E-2_wp, &
      & 5.32310785129656E-2_wp, 5.56201199842301E-2_wp,-1.94185137737050E-1_wp, &
      &-1.69435262877121E-2_wp, 4.06941497321145E-2_wp, 1.40954059224083E-1_wp, &
      &-1.77841513693408E-2_wp,-2.75394740661722E-2_wp, 8.23579661010067E-3_wp, &
      & 3.34177535391165E-2_wp,-5.71071859403073E-3_wp, 9.54835475924012E-3_wp, &
      & 2.42936845781245E-2_wp, 5.48198008679115E-3_wp, 9.79238907789204E-4_wp, &
      & 5.02629146832239E-3_wp,-3.26133215380267E-2_wp,-2.52729234859127E-2_wp, &
      & 6.13790503829809E-1_wp, 5.32657270579180E-1_wp, 1.54187752434240E-1_wp, &
      &-5.16490059386104E-1_wp, 5.67205701727938E-2_wp,-7.67978256264049E-1_wp, &
      & 7.83646726355242E-2_wp, 7.93496900696582E-3_wp,-6.64025879132266E-2_wp, &
      &-2.39938346396998E-2_wp, 2.27466901433374E-2_wp,-1.19620847222975E-2_wp, &
      &-2.05612814844119E-2_wp,-4.95270994769136E-2_wp, 3.87946677272715E-2_wp, &
      & 2.09792192994460E-2_wp, 1.57023828153490E-2_wp,-1.82333862428595E-2_wp, &
      & 2.06990008614840E-1_wp,-9.40931225186795E-1_wp,-9.81260233409964E-1_wp, &
      & 1.81683320051775E-1_wp,-2.29434219507657E-1_wp, 7.74270224795126E-1_wp, &
      &-4.49851593643217E-2_wp,-1.43120792142963E-2_wp, 8.76827127392763E-2_wp, &
      & 7.95606845308971E-5_wp, 1.49470897135840E-2_wp,-4.26975533749547E-2_wp, &
      &-3.98095736465717E-2_wp,-4.34657048996160E-2_wp, 1.60460563069101E-2_wp, &
      & 1.10714431698586E-2_wp, 2.00213367980184E-2_wp, 2.37635173396616E-2_wp, &
      & 2.31027463622808E-2_wp,-3.26295375143408E-2_wp,-1.33550726444072E-1_wp, &
      &-1.53667922610297E-1_wp, 3.82022829991620E-2_wp, 1.10447980081791E-1_wp, &
      &-4.50996897167389E-1_wp,-1.02330031032349E+0_wp, 8.40160246681131E-2_wp, &
      &-3.56277305140736E-1_wp,-6.25766606417323E-1_wp, 3.66980872499278E-1_wp, &
      & 2.08713890794846E-1_wp, 8.09936390516330E-1_wp,-1.01293604136016E+0_wp, &
      & 3.35301352962707E-1_wp, 3.38610231741285E-1_wp, 8.04222150565318E-1_wp, &
      & 2.87096625886200E-2_wp, 3.36901296815917E-2_wp,-2.68793522207332E-3_wp, &
      &-2.05780563238218E-2_wp,-6.21642900499840E-3_wp,-2.60217273665459E-2_wp, &
      &-1.50087159351410E-2_wp,-3.69461704433247E-3_wp,-1.82860627486028E-2_wp, &
      & 4.40210284335156E-2_wp, 2.60057075458535E-2_wp, 3.32947786837439E-2_wp, &
      & 4.27111264045819E-1_wp, 4.85129382065812E-2_wp,-6.56980004019461E-1_wp, &
      &-1.30593579650094E-1_wp,-2.07993009790954E-1_wp, 2.29868739973643E-1_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol, qat, dpat, qpat, make_multipole_gfn2, thr_in=thr1)

end subroutine test_g_aes_gfn2_m04

subroutine test_g_aes_gfn2_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & 2.14480738661017E-1_wp, 2.14527495272549E-1_wp, 2.14471598997661E-1_wp, &
      & 2.14516468108033E-1_wp, 2.45821990520972E-1_wp, 2.45784832080333E-1_wp, &
      & 2.45806360746891E-1_wp, 2.45774007559891E-1_wp, 3.19474245000812E-1_wp, &
      & 3.19579496746895E-1_wp,-2.94626505330352E-1_wp,-2.94543754718096E-1_wp, &
      &-2.94659848864979E-1_wp,-2.94567600483171E-1_wp,-6.50780911474268E-1_wp, &
      &-6.51058612824167E-1_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      & 8.94827368009637E-2_wp, 8.94995909146529E-2_wp, 7.30634961756909E-2_wp, &
      & 8.95014203078334E-2_wp,-8.94859532646144E-2_wp,-7.30378298851327E-2_wp, &
      &-8.94873352015873E-2_wp,-8.95044727986339E-2_wp, 7.30601341120621E-2_wp, &
      &-8.95090716572443E-2_wp, 8.94905477532190E-2_wp,-7.30362173059338E-2_wp, &
      &-9.01917250610402E-3_wp,-9.00201247838739E-3_wp,-1.47850060857521E-1_wp, &
      &-9.00995133125046E-3_wp, 9.02740339798291E-3_wp, 1.47848725725542E-1_wp, &
      & 9.01739468689571E-3_wp, 9.00092652137544E-3_wp,-1.47851490297087E-1_wp, &
      & 9.00924018813641E-3_wp,-9.02542526993692E-3_wp, 1.47850396429899E-1_wp, &
      & 3.58533790400338E-6_wp, 1.76804223012015E-5_wp, 7.90811791429655E-2_wp, &
      & 1.90955806626401E-5_wp, 1.99177597125356E-6_wp,-7.93788543388031E-2_wp, &
      &-4.56180855347671E-2_wp,-4.55680071125238E-2_wp, 5.92409368467603E-2_wp, &
      &-4.56101826239023E-2_wp, 4.56639049054392E-2_wp,-5.91934296227654E-2_wp, &
      & 4.55978078230863E-2_wp, 4.55693703731601E-2_wp, 5.92330972112299E-2_wp, &
      & 4.56164302678149E-2_wp,-4.56412184067195E-2_wp,-5.91891814963370E-2_wp, &
      & 6.11781066870311E-6_wp, 9.15928265940359E-7_wp,-2.04395992132628E-1_wp, &
      & 1.78950595055862E-6_wp,-4.23972317385840E-6_wp, 2.04360605943730E-1_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      &-1.23698956475140E-2_wp, 1.79751870775559E-2_wp,-1.23691904599371E-2_wp, &
      & 4.67729267242994E-3_wp, 4.68004166312949E-3_wp, 2.47390861074511E-2_wp, &
      &-1.23654488813711E-2_wp,-1.79817963662729E-2_wp,-1.23706010436822E-2_wp, &
      &-4.69179866811513E-3_wp, 4.68829261745078E-3_wp, 2.47360499250533E-2_wp, &
      &-1.23743195999797E-2_wp, 1.79683161428587E-2_wp,-1.23718141043582E-2_wp, &
      &-4.68144974224393E-3_wp,-4.68203929764922E-3_wp, 2.47461337043379E-2_wp, &
      &-1.23715411335551E-2_wp,-1.79777626390457E-2_wp,-1.23696001364202E-2_wp, &
      & 4.69398264928409E-3_wp,-4.69414370569245E-3_wp, 2.47411412699752E-2_wp, &
      &-7.32193778340168E-3_wp, 3.45599884203983E-2_wp,-7.32434931893968E-3_wp, &
      &-6.43942937553289E-3_wp,-6.44505606234306E-3_wp, 1.46462871023413E-2_wp, &
      &-7.32508525229408E-3_wp,-3.45618871061652E-2_wp,-7.32335313431533E-3_wp, &
      & 6.45162772813276E-3_wp,-6.44644470074471E-3_wp, 1.46484383866092E-2_wp, &
      &-7.32905564328019E-3_wp, 3.45586433766952E-2_wp,-7.32543310737930E-3_wp, &
      & 6.44113915814194E-3_wp, 6.44697435437880E-3_wp, 1.46544887506594E-2_wp, &
      &-7.32605090318348E-3_wp,-3.45610747417528E-2_wp,-7.32899993430305E-3_wp, &
      &-6.45372642877205E-3_wp, 6.44746289062807E-3_wp, 1.46550508374864E-2_wp, &
      & 1.11951790334487E-1_wp,-4.02377811299945E-1_wp, 1.12018181924718E-1_wp, &
      &-1.04989650284058E-6_wp,-1.61793643693583E-5_wp,-2.23969972259205E-1_wp, &
      & 1.11996170927313E-1_wp, 4.02383573647660E-1_wp, 1.11929776104217E-1_wp, &
      & 1.69716696800307E-5_wp, 3.87072152745758E-6_wp,-2.23925947031530E-1_wp, &
      & 3.55892547851532E-2_wp,-1.40674101590403E-3_wp, 3.55796422243525E-2_wp, &
      &-4.51437869371043E-2_wp,-4.51403945941909E-2_wp,-7.11688970095042E-2_wp, &
      & 3.56007063955610E-2_wp, 1.39242785502690E-3_wp, 3.56047482812895E-2_wp, &
      & 4.51288756967405E-2_wp,-4.51327168181670E-2_wp,-7.12054546768490E-2_wp, &
      & 3.55791799999949E-2_wp,-1.38863474463774E-3_wp, 3.55709732649653E-2_wp, &
      & 4.51424181512024E-2_wp, 4.51358771366434E-2_wp,-7.11501532649592E-2_wp, &
      & 3.55872643357090E-2_wp, 1.37913694676448E-3_wp, 3.56010570423592E-2_wp, &
      &-4.51277501739264E-2_wp, 4.51338445551804E-2_wp,-7.11883213780670E-2_wp, &
      &-2.34031908350780E-2_wp,-5.97048165007334E-2_wp,-2.34051674127262E-2_wp, &
      &-1.53933429016999E-6_wp, 9.26688498685427E-7_wp, 4.68083582478040E-2_wp, &
      &-2.34383117547964E-2_wp, 5.96710054013645E-2_wp,-2.34363320491785E-2_wp, &
      &-2.13796886588954E-6_wp,-2.36716531036422E-6_wp, 4.68746438039748E-2_wp],&
      & shape(qpat))

   call get_structure(mol, "X23", "urea")
   call test_numgrad(error, mol, qat, dpat, qpat, make_multipole_gfn2, thr_in=thr1)

end subroutine test_g_aes_gfn2_urea


subroutine test_g_aes_gxtb_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & -6.72232589509782E-2_wp, -9.34334647769571E-1_wp,  7.42318469479949E-2_wp, &
      &  6.00334278157932E-1_wp,  8.15555915271801E-1_wp,  1.11908749754635E-1_wp, &
      & -4.00078932773686E-1_wp,  8.01675686522246E-2_wp,  8.00529956280133E-2_wp, &
      &  1.00638981680234E-1_wp, -6.70400663640075E-1_wp, -4.32671353819512E-1_wp, &
      & -5.26347736891704E-1_wp,  2.69243405687243E-1_wp,  8.19036648298260E-1_wp, &
      &  7.98862032739356E-2_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      & -1.94227933220811E-1_wp,  2.43336860340003E-1_wp,  1.86135458816528E-1_wp, &
      & -3.41540997137134E-2_wp, -4.53783242583224E-2_wp,  6.49297195316456E-2_wp, &
      & -1.09904296007690E-2_wp, -3.07211207018219E-2_wp,  3.91362229926026E-2_wp, &
      &  1.58619265813120E-1_wp, -1.60898423006623E-1_wp,  1.30694231257800E-1_wp, &
      &  4.78407364655061E-1_wp,  3.88191026914963E-1_wp, -3.41277133661884E-2_wp, &
      &  6.83026550146984E-1_wp,  6.20707987982817E-1_wp, -1.75174999945350E-1_wp, &
      & -5.36638502783345E-2_wp,  1.60568022993560E-1_wp, -3.79772809906278E-1_wp, &
      &  2.30730194779660E-2_wp, -6.88119315791599E-2_wp, -1.87837649754506E-2_wp, &
      & -4.61158388280005E-2_wp,  3.64753029218826E-2_wp, -1.01050094183562E-2_wp, &
      &  2.04951416768851E-2_wp, -3.94296737450275E-2_wp, -9.49332984399296E-2_wp, &
      &  5.71090370126919E-2_wp,  1.08772152349226E-2_wp,  2.71143855136837E-2_wp, &
      & -4.93966139471334E-2_wp, -2.18076225995155E-1_wp, -3.81968718965382E-2_wp, &
      & -3.98701501982369E-2_wp, -2.20206495445272E-2_wp,  4.14281401312083E-5_wp, &
      & -5.29490161853578E-3_wp, -1.26101724854378E-2_wp,  2.14620831911543E-2_wp, &
      &  4.62772421352236E-2_wp,  1.33545247622996E-1_wp, -1.02323590844664E-1_wp, &
      &  2.70093985269488E-2_wp,  3.20959513717059E-2_wp,  3.55554757232086E-2_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      &  3.74428400131661E-2_wp, -2.91513236600262E-1_wp,  4.68310997912392E-1_wp, &
      & -1.00092568690606E-1_wp, -3.64976081689563E-1_wp, -5.05753837925559E-1_wp, &
      & -9.03389076133383E-2_wp,  2.35415513131662E-2_wp, -1.44281760469589E-2_wp, &
      &  7.96985085096959E-3_wp, -1.96822284690929E-1_wp,  1.04767083660301E-1_wp, &
      & -7.61475413149074E-2_wp,  4.23132648712838E-2_wp, -2.40582726063123E-2_wp, &
      & -2.72793820825515E-2_wp, -1.30368111336366E-1_wp,  1.00205813921220E-1_wp, &
      & -8.30155278808316E-2_wp,  4.11533409128402E-1_wp, -1.85703044875811E-1_wp, &
      &  9.22475114839877E-2_wp,  4.53573199701716E-1_wp,  2.68718572756643E-1_wp, &
      & -8.46066783469599E-2_wp,  6.03746141427654E-2_wp, -3.42015389599566E-1_wp, &
      &  4.83931550466230E-1_wp,  3.11769207439250E-1_wp,  4.26622067946527E-1_wp, &
      & -8.08140048136396E-1_wp, -3.97381960460994E-1_wp, -3.98355918697151E-1_wp, &
      & -9.86110632820629E-1_wp,  7.68956819956569E-1_wp,  1.20649596683355E+0_wp, &
      &  4.81228272379237E-2_wp, -1.20517823853204E-1_wp, -9.21456520305997E-2_wp, &
      & -3.28868050418220E-2_wp, -4.86071679226612E-1_wp,  4.40228247926751E-2_wp, &
      & -7.29251442652190E-2_wp, -4.95354547808566E-2_wp,  1.03567857031379E-1_wp, &
      & -2.02533940885859E-2_wp, -1.21110943064115E-2_wp, -3.06427127661604E-2_wp, &
      &  5.54865158499682E-2_wp, -1.20949169949538E-1_wp,  4.67508090487394E-2_wp, &
      &  1.57598547749525E-2_wp,  4.66785183570980E-3_wp, -1.02237324898708E-1_wp, &
      & -7.38954141220622E-2_wp, -4.35753104219451E-2_wp,  5.47968373533771E-2_wp, &
      & -1.77414489117915E-2_wp, -2.27470376197845E-2_wp,  1.90985767686849E-2_wp, &
      &  4.76602906802981E-3_wp,  1.46977567614165E-2_wp,  1.48782740158886E-2_wp, &
      & -3.58333181171945E-2_wp, -1.16229158881186E-4_wp, -1.96443030839202E-2_wp, &
      & -8.54659753146029E-1_wp, -8.39610255129518E-1_wp,  1.48688365589127E+0_wp, &
      & -3.01371250571339E-1_wp, -4.30737799462750E-1_wp, -6.32223902745247E-1_wp, &
      & -7.28741798745070E-2_wp,  1.54198956710758E-1_wp, -2.27812636470315E-2_wp, &
      &  2.04054057701659E-1_wp,  2.26435982953830E-1_wp,  9.56554435215367E-2_wp, &
      & -8.59721877698228E-2_wp,  3.11934701285617E-2_wp, -9.10074517458215E-3_wp, &
      & -4.59104795552940E-2_wp, -1.29435873370815E-1_wp,  9.50729329444049E-2_wp, &
      &  3.34239240829591E-1_wp, -1.02941144084035E-1_wp, -2.67919313396491E-1_wp, &
      &  4.92653475236327E-1_wp, -3.12245402295760E-1_wp, -6.63199274330999E-2_wp, &
      &  1.17901614434976E-1_wp,  1.22710689853952E-1_wp, -1.34237152075345E-2_wp, &
      &  6.20197054432496E-2_wp,  1.69941132047634E-2_wp, -1.04477899227441E-1_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "03")
   call test_numgrad(error, mol, qat, dpat, qpat, make_multipole_gxtb, thr_in=thr1)

end subroutine test_g_aes_gxtb_m03

subroutine test_g_aes_gxtb_m04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & -6.41747650934630E-2_wp, -2.25666193168372E-1_wp, -3.45862671449368E-2_wp, &
      & -3.86060092404173E-1_wp,  7.23899147269319E-1_wp, -5.77708449165084E-2_wp, &
      &  1.19612041167111E-1_wp, -5.46111760044127E-3_wp, -1.79451438301903E-2_wp, &
      & -1.21151427608654E-1_wp, -3.49373668911567E-1_wp,  9.24468471315697E-1_wp, &
      &  1.39403005585106E-1_wp, -6.97541217849458E-1_wp,  1.14935657028527E-1_wp, &
      & -6.25875842906434E-2_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      & -3.57421953327594E-2_wp,  1.25509246271253E-2_wp,  2.37526154738796E-2_wp, &
      & -5.17775860693660E-2_wp, -1.47299624317204E-1_wp, -1.03044617845191E-1_wp, &
      & -6.55007227511516E-2_wp,  8.62513013027920E-2_wp, -7.40695647029688E-2_wp, &
      & -2.38460729006625E-2_wp,  7.62860392723061E-2_wp,  8.13999540001495E-2_wp, &
      & -1.79362908209604E-1_wp,  8.53481540574998E-2_wp, -5.02178529783988E-1_wp, &
      &  3.88411230580228E-2_wp,  1.37013298380537E-2_wp, -1.26252576885668E-2_wp, &
      &  4.53012010688574E-2_wp, -6.48069907146587E-2_wp,  4.62760975241764E-2_wp, &
      &  3.78680714343827E-1_wp,  1.19317413830196E-1_wp, -5.08257518966121E-1_wp, &
      & -1.72635262584643E-3_wp,  1.09885891220385E-1_wp,  3.33183309773249E-2_wp, &
      &  6.17707820193476E-2_wp, -4.16566596458103E-2_wp,  1.44480373398104E-1_wp, &
      & -7.36027621360482E-2_wp,  1.95274279149139E-1_wp,  2.20023954919709E-1_wp, &
      &  4.27483506970933E-1_wp,  1.35628129313216E-1_wp,  4.90386843827805E-1_wp, &
      &  3.07694279524434E-2_wp,  6.42050804506663E-2_wp,  2.09246738801646E-1_wp, &
      &  8.42589862787275E-2_wp, -1.32376737276813E-1_wp,  1.70143500230254E-2_wp, &
      & -6.18055726324971E-2_wp,  9.28086124840346E-3_wp, -6.29412993260402E-2_wp, &
      & -1.31625836962117E-1_wp, -2.80633086809525E-1_wp,  5.57290139902586E-2_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      &  1.74919970405689E-2_wp,  1.00989591675384E-1_wp, -1.60551660896662E-1_wp, &
      & -9.58729857096069E-2_wp, -8.55550676972757E-2_wp,  1.43059663856092E-1_wp, &
      & -3.45079241330293E-2_wp,  1.83110249317935E-1_wp, -3.59176992291322E-1_wp, &
      &  2.75323883571215E-1_wp,  1.87650363262677E-1_wp,  3.93684916424350E-1_wp, &
      & -4.40436732999174E-2_wp, -3.55020617697798E-2_wp,  4.40520032100995E-2_wp, &
      &  5.94812748848310E-2_wp, -2.82322880169285E-2_wp, -8.32991018206997E-6_wp, &
      & -1.32391504068836E-1_wp, -9.27006770868815E-2_wp,  3.12286633811363E-2_wp, &
      & -4.37042705709019E-3_wp,  1.39755144003930E-1_wp,  1.01162840687699E-1_wp, &
      &  2.99318531326219E-1_wp,  2.74801519487720E-1_wp,  3.84594332611328E-1_wp, &
      & -6.36759964899869E-1_wp,  4.96751933862765E-1_wp, -6.83912863937548E-1_wp, &
      &  2.25415288048392E-1_wp, -3.38927290417661E-2_wp, -1.90353528021403E-1_wp, &
      & -3.48714233136383E-3_wp,  9.93905542122104E-4_wp, -3.50617600269882E-2_wp, &
      & -3.87327985944474E-2_wp, -6.91174139519020E-2_wp,  5.49419584154653E-2_wp, &
      &  2.17665409988861E-2_wp,  1.34164608145337E-2_wp, -1.62091598210179E-2_wp, &
      & -1.50326484753311E-1_wp, -5.32948770320869E-2_wp, -5.88727278060565E-1_wp, &
      & -1.07806833740847E+0_wp, -2.50683784340449E-1_wp,  7.39053762813878E-1_wp, &
      & -8.32528646057871E-2_wp, -9.71601178047237E-3_wp,  1.37155697721370E-1_wp, &
      &  8.92639600684630E-3_wp,  5.02842549842410E-2_wp, -5.39028331155829E-2_wp, &
      & -3.74851140591499E-2_wp, -6.77339647622242E-2_wp,  5.74803705621472E-2_wp, &
      & -7.94727661704886E-3_wp,  4.13419432361995E-2_wp, -1.99952565029973E-2_wp, &
      & -3.39878814477004E-1_wp, -8.98908743783573E-2_wp,  6.07231719489638E-2_wp, &
      & -5.38311434058213E-1_wp,  6.78493406110892E-1_wp,  2.79155642528040E-1_wp, &
      &  2.97910585250318E-1_wp, -4.81269125949961E-1_wp, -5.43604732937434E-1_wp, &
      & -4.35058245092393E-1_wp, -6.33509997861268E-1_wp,  2.45694147687118E-1_wp, &
      & -2.01089321484726E-1_wp,  9.13836260032944E-1_wp,  4.68453410810543E-2_wp, &
      &  7.55393040583326E-3_wp, -4.24914432086360E-2_wp,  1.54243980403667E-1_wp, &
      & -8.28769100177347E-2_wp, -3.09911485849456E-1_wp,  1.90499909705355E-1_wp, &
      &  1.71898946913901E-1_wp, -1.17032143519923E-1_wp, -1.07622999687620E-1_wp, &
      & -1.92022483522591E-2_wp, -4.83056876496173E-3_wp, -3.55014219429624E-2_wp, &
      &  4.79270951486242E-2_wp,  4.02023252289722E-2_wp,  5.47036702952215E-2_wp, &
      &  5.30751552125004E-1_wp,  3.12181884074276E-1_wp, -2.81075233594611E-1_wp, &
      & -1.37230587154048E-1_wp,  3.22352214632118E-2_wp, -2.49676318530392E-1_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol, qat, dpat, qpat, make_multipole_gxtb, thr_in=thr1)

end subroutine test_g_aes_gxtb_m04

subroutine test_g_aes_gxtb_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(4) = [&
      &  1.76167728179214E+0_wp, -5.86846576518797E-1_wp, -5.88294165118524E-1_wp, &
      & -5.86536540306923E-1_wp]
   real(wp), parameter :: dpat(3, 4) = reshape([&
      &  1.03102227915076E-1_wp,  1.70815542489309E-1_wp, -2.59639799151207E-1_wp, &
      & -1.87602572329602E-1_wp, -2.74737606941817E-1_wp, -1.49657158899311E-1_wp, &
      & -2.03472803217469E-1_wp,  2.74583536877613E-1_wp,  1.60701200271541E-1_wp, &
      &  3.20057753764996E-1_wp, -9.88144341108652E-2_wp,  1.55919293968284E-1_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 4) = reshape([&
      & -7.51267637005315E-1_wp,  3.31653256315932E-1_wp, -2.26838871462153E-1_wp, &
      & -9.32402326710803E-1_wp, -1.22767701815665E+0_wp,  9.78106508467468E-1_wp, &
      &  2.61814638396274E-3_wp, -2.16516052249151E-2_wp, -9.92760600427545E-3_wp, &
      & -2.24232173276705E-2_wp, -3.20385899497410E-2_wp,  7.30945962031626E-3_wp, &
      & -1.81977682871093E-3_wp,  3.29176108412411E-2_wp, -2.14330340960966E-2_wp, &
      &  5.73574812947845E-3_wp, -3.15197707756747E-2_wp,  2.32528109248094E-2_wp, &
      & -4.08481125421503E-2_wp,  1.54832741157508E-2_wp,  1.99270553528201E-2_wp, &
      & -2.74175062221744E-2_wp, -3.37992023179355E-3_wp,  2.09210571893266E-2_wp],&
      & shape(qpat))

   call get_structure(mol, "f-block", "CeCl3")
   call test_numgrad(error, mol, qat, dpat, qpat, make_multipole_gxtb, thr_in=thr1)

end subroutine test_g_aes_gxtb_cecl3


subroutine test_s_aes_gfn2_m05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & 2.29216160742194E-1_wp, 4.71025340158282E-2_wp,-3.52853843840331E-2_wp, &
      & 3.09109494571158E-2_wp, 3.28899921234171E-1_wp,-2.01752645265381E-1_wp, &
      & 5.46543133415123E-2_wp,-1.09690925298550E-1_wp,-3.47346448335524E-1_wp, &
      & 1.84561311715549E-1_wp,-2.07552642508510E-1_wp, 4.67143385624953E-1_wp, &
      &-1.84214559781700E-2_wp,-1.05015851128102E-1_wp, 6.52406633071011E-2_wp, &
      &-3.82663886540142E-1_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      &-1.36866725567795E-1_wp, 8.98162734808879E-2_wp, 1.93992067964322E-1_wp, &
      &-7.94862211844545E-2_wp, 1.72156341648282E-1_wp, 1.09913940611171E-1_wp, &
      & 4.98989463722096E-3_wp, 5.29170302861464E-2_wp,-2.76614889872771E-2_wp, &
      & 9.88532719472691E-2_wp,-1.40034980230081E-2_wp, 5.34276261339990E-3_wp, &
      &-1.43509540880904E-1_wp,-1.53177661284083E-1_wp,-1.44819803333425E-1_wp, &
      &-2.04156814828549E-2_wp,-1.36595645633772E-1_wp, 9.21268600746920E-2_wp, &
      & 4.82467726581279E-2_wp,-6.83638402837403E-2_wp, 2.25756519899868E-2_wp, &
      &-1.07827737715754E-1_wp,-1.57029392143043E-2_wp,-4.52817748323062E-3_wp, &
      &-5.75565163095502E-2_wp, 8.04433988598797E-2_wp,-6.18156660422769E-2_wp, &
      &-5.40875924342744E-2_wp, 8.84669124477697E-2_wp, 5.79993467185712E-2_wp, &
      & 1.17647171474739E-1_wp, 1.47697439492523E-2_wp, 1.01604279159870E-1_wp, &
      &-5.04764399493363E-1_wp,-9.83938684357523E-3_wp,-3.12571551890408E-1_wp, &
      & 3.55860816794651E-2_wp,-1.31041235123739E-1_wp,-8.75524850448995E-2_wp, &
      & 7.29800329344290E-2_wp,-9.87358811335129E-2_wp,-9.28538506341673E-2_wp, &
      &-8.73682094598224E-2_wp,-2.05679212716268E-1_wp,-1.33550747139953E-1_wp, &
      & 8.30672094549445E-2_wp, 6.79087546790233E-2_wp, 7.61092730005571E-3_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      &-3.21923529804590E-2_wp, 1.01247911778500E-1_wp,-1.60153050884742E-1_wp, &
      & 2.09688971314230E-1_wp,-5.53870278646362E-1_wp, 1.92345403865201E-1_wp, &
      & 1.42745048693680E-1_wp, 3.96650238849656E-1_wp,-1.85556055670220E-2_wp, &
      & 2.09624026315852E-1_wp,-1.00919824206596E-1_wp,-1.24189443126656E-1_wp, &
      &-4.76658863743303E-4_wp, 2.90471575684891E-2_wp, 1.07961071436757E-1_wp, &
      &-9.53915574675821E-2_wp, 3.36258209533697E-2_wp,-1.07484412573014E-1_wp, &
      & 1.27092159973006E-1_wp, 1.21655068472347E-2_wp,-5.40941977072504E-2_wp, &
      &-1.59997403617732E-2_wp,-3.07740637834542E-6_wp,-7.29979622657559E-2_wp, &
      & 5.53303994778954E-1_wp,-2.99239292851686E-1_wp,-2.51102054122008E-1_wp, &
      &-6.26827633889481E-1_wp,-8.42332212885064E-1_wp,-3.02201940656946E-1_wp, &
      & 2.19857432174074E-1_wp, 2.27516201446880E-1_wp,-2.10413952758986E-1_wp, &
      &-4.07388945048188E-1_wp,-1.26641550036282E-1_wp,-9.44347941508889E-3_wp, &
      & 2.93948297678871E-2_wp,-3.43251555808793E-2_wp,-2.11422152270078E-2_wp, &
      & 5.41372742690424E-2_wp,-3.39231131165861E-2_wp,-8.25261454087947E-3_wp, &
      & 2.37070130746082E-1_wp, 3.85988938198882E-3_wp,-1.72140997530014E-1_wp, &
      &-1.84923389243440E-1_wp,-2.83535528712180E-3_wp,-6.49291332160717E-2_wp, &
      &-5.91277580612433E-3_wp, 2.78962136409584E-2_wp,-4.44146350092125E-2_wp, &
      & 4.21395120616379E-2_wp, 1.77901916517195E-2_wp, 5.03274108153378E-2_wp, &
      &-1.69793065817916E-2_wp,-4.35671755332885E-2_wp,-8.26169504201363E-3_wp, &
      &-3.08726223882655E-2_wp, 7.19218791920170E-2_wp, 2.52410016238050E-2_wp, &
      &-1.16489039702432E-1_wp, 1.21709469870565E-1_wp, 1.61694380019930E-1_wp, &
      &-1.78938261587284E-1_wp,-2.40115227627245E-1_wp,-4.52053403174976E-2_wp, &
      & 9.96286547963425E-1_wp,-1.19749514806617E-1_wp,-1.35617319657923E+0_wp, &
      & 8.55455083538778E-1_wp,-9.90077412493467E-1_wp, 3.59886648615801E-1_wp, &
      & 5.22916600071566E-2_wp,-3.91715496559764E-2_wp, 6.83298122566691E-3_wp, &
      &-2.53231481131226E-2_wp,-5.89822395445022E-2_wp,-5.91246412328236E-2_wp, &
      & 4.75850655614168E-3_wp,-1.05070993090875E-2_wp,-9.80745234500818E-3_wp, &
      &-1.41037515552590E-2_wp, 3.45290550320865E-2_wp, 5.04894578886597E-3_wp, &
      & 3.84056205679418E-1_wp, 4.90928417189913E-1_wp,-4.20009072583480E-1_wp, &
      & 1.27378622840086E-1_wp, 2.70322396743841E-1_wp, 3.59528669040627E-2_wp, &
      & 1.02687867732048E-1_wp,-9.41266714173517E-2_wp,-1.57143312769043E-2_wp, &
      &-1.16797766346470E-1_wp,-2.25346145494165E-1_wp,-8.69735364551442E-2_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "05")
   call test_numsigma(error, mol, qat, dpat, qpat, make_multipole_gfn2, thr_in=thr1)

end subroutine test_s_aes_gfn2_m05

subroutine test_s_aes_gfn2_m06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & 2.08158722831082E-1_wp,-3.78015439235964E-1_wp, 3.36498612504260E-2_wp, &
      &-4.11549915371290E-1_wp, 8.14985792195240E-2_wp,-2.00895089502942E-1_wp, &
      & 2.44763776148029E-1_wp, 2.54570821381252E-2_wp, 2.59846290341396E-1_wp, &
      & 4.21681978758186E-1_wp, 1.37097577114182E-1_wp, 4.06943998730574E-2_wp, &
      &-1.10962379048735E-1_wp,-6.44071741016567E-2_wp,-1.91525573827717E-1_wp, &
      &-9.54926965857112E-2_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      &-1.38757863526336E-1_wp, 9.44281294694277E-2_wp,-2.20905494810822E-2_wp, &
      &-6.76271661399369E-2_wp, 1.28863877510786E-1_wp, 6.64008692108340E-3_wp, &
      & 1.05845053789867E-1_wp, 3.20756596369100E-2_wp,-2.61616134638101E-2_wp, &
      &-1.54692589602203E-2_wp,-1.71325991463100E-2_wp, 4.16583321796550E-3_wp, &
      &-5.10286551685211E-2_wp,-1.24258737415072E-1_wp,-1.68925560110076E-1_wp, &
      & 3.89780935454439E-2_wp, 2.74052865410239E-3_wp,-2.21400648822050E-2_wp, &
      & 2.49874310362634E-2_wp, 3.27419756444858E-2_wp, 6.18542730650856E-2_wp, &
      &-7.81265682720928E-2_wp, 8.35615333909530E-2_wp,-4.04047772976296E-2_wp, &
      &-7.60911516191494E-2_wp,-1.44657184795252E-1_wp, 2.74164556915873E-2_wp, &
      &-3.35721006884757E-1_wp,-3.04442343036050E-1_wp, 4.42474331429686E-1_wp, &
      &-8.42480506238376E-2_wp, 6.30134856815068E-2_wp,-7.97174973506382E-3_wp, &
      & 5.28799020569297E-3_wp, 9.67408586419737E-3_wp, 1.28093823252530E-1_wp, &
      & 2.47890384259527E-2_wp, 1.32250383459731E-2_wp, 1.15433780917482E-1_wp, &
      &-3.80693744454969E-4_wp,-7.01021372703028E-2_wp, 3.11636557277436E-2_wp, &
      & 3.22305990386226E-2_wp, 4.60403860383308E-2_wp, 5.62539559627448E-3_wp, &
      & 7.02569888935242E-2_wp,-2.90643759577685E-2_wp,-3.67626155472377E-2_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      &-3.13024148290062E-1_wp, 2.61584538243456E-1_wp,-4.61651985952017E-1_wp, &
      &-2.91003023546858E-1_wp,-1.22014295662189E-1_wp, 7.74676134242079E-1_wp, &
      &-1.34443825936356E-3_wp,-3.99956235299210E-2_wp,-7.47799317738422E-2_wp, &
      &-4.06754103770418E-2_wp,-5.00955621101496E-2_wp, 7.61243700332053E-2_wp, &
      & 7.99137125232521E-2_wp, 3.60405700863902E-2_wp,-2.33868262236604E-2_wp, &
      &-5.05609303399179E-2_wp,-2.41867928499022E-2_wp,-5.65268862995917E-2_wp, &
      &-2.28047139606244E-2_wp,-2.53620428449374E-1_wp,-1.40935360788743E-1_wp, &
      & 5.42468357821061E-2_wp, 6.77550434821008E-2_wp, 1.63740074749368E-1_wp, &
      & 9.84138860796349E-2_wp, 2.05401786481362E-1_wp, 5.08930057773724E-2_wp, &
      & 6.18355817777164E-2_wp,-8.79408298722754E-2_wp,-1.49306891857007E-1_wp, &
      &-3.67037066569778E-2_wp, 3.04612985342592E-2_wp,-3.62087034398719E-4_wp, &
      & 1.52537104446144E-2_wp, 3.53779868504119E-2_wp, 3.70657936913765E-2_wp, &
      & 1.40253565352209E-2_wp,-1.24078640470215E+0_wp,-7.56120808182631E-1_wp, &
      & 1.50109582306222E+0_wp, 3.57103379871705E-1_wp, 7.42095451647409E-1_wp, &
      & 1.87510950260991E-2_wp,-5.86003042127686E-2_wp, 2.54817747989697E-2_wp, &
      & 4.17308181084498E-2_wp,-3.90138382941339E-2_wp,-4.42328698250688E-2_wp, &
      &-7.51608441513982E-3_wp,-9.66374364137592E-1_wp,-6.22437474815178E-1_wp, &
      & 1.81800799781093E-1_wp, 2.73636244652070E-1_wp, 6.29953559230319E-1_wp, &
      & 3.87552943603469E-1_wp, 2.27870646813054E-2_wp,-4.34009118694632E-1_wp, &
      &-6.15147736496459E-1_wp,-3.50708689472064E-1_wp, 4.64561750911630E-2_wp, &
      & 4.19561927014193E-2_wp,-7.87542156734094E-2_wp, 3.00053452114131E-2_wp, &
      & 1.91766065522965E-2_wp, 3.42797457855182E-3_wp,-7.19615379128325E-2_wp, &
      &-5.22181543711780E-2_wp,-1.10008395468616E-2_wp,-5.37722495332387E-2_wp, &
      &-2.73878079383160E-3_wp, 1.20687681297706E-2_wp, 1.05990403904417E-1_wp, &
      &-4.64327979725137E-2_wp, 1.72801416897558E-2_wp,-3.14880932493609E-2_wp, &
      & 6.19928751082781E-2_wp, 3.85665458330814E-2_wp, 7.79208912218747E-2_wp, &
      &-6.90767772589938E-3_wp, 5.79645601278631E-2_wp, 3.15679788385037E-2_wp, &
      & 1.02634959072458E-2_wp, 5.40387161648890E-3_wp,-2.46603011126042E-2_wp, &
      &-4.37662521291874E-3_wp, 3.83793053174837E-3_wp, 3.07085662485559E-2_wp, &
      &-6.05511439994602E-3_wp, 9.50209584313457E-3_wp,-2.63319410356369E-2_wp, &
      & 6.56851486385975E-2_wp,-2.51206026224842E-2_wp,-5.57383484000740E-2_wp, &
      &-4.45029464021450E-2_wp, 3.10148202231310E-3_wp,-9.94680023852338E-3_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "06")
   call test_numsigma(error, mol, qat, dpat, qpat, make_multipole_gfn2, thr_in=thr1)

end subroutine test_s_aes_gfn2_m06

subroutine test_s_aes_gfn2_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      & 3.83148728601709E-1_wp, 3.83118715934305E-1_wp, 3.83080191213008E-1_wp, &
      & 3.82184417033339E-1_wp, 3.47311885772574E-1_wp, 3.45281176165792E-1_wp, &
      & 3.47791244554054E-1_wp, 3.47952947802927E-1_wp,-3.92870149599810E-1_wp, &
      &-3.92136115450962E-1_wp,-3.93531741241129E-1_wp,-3.92928999213039E-1_wp, &
      &-3.37416498051037E-1_wp,-3.36550419758369E-1_wp,-3.37512222809576E-1_wp, &
      &-3.36923160953788E-1_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      & 4.09688027943853E-2_wp, 8.33483477302867E-3_wp,-1.09648576234461E-2_wp, &
      &-4.09474996408295E-2_wp, 8.30233124676478E-3_wp, 1.09874908009843E-2_wp, &
      &-4.09278272213102E-2_wp,-8.32793028126641E-3_wp, 1.09615760238990E-2_wp, &
      & 4.11086171124463E-2_wp,-8.39441030068160E-3_wp,-1.05971938523367E-2_wp, &
      & 5.96300641349022E-2_wp,-4.34579903309867E-2_wp, 6.26025852647896E-2_wp, &
      &-5.98049356163273E-2_wp,-4.34306642465664E-2_wp,-6.22704459500676E-2_wp, &
      &-6.00922329406000E-2_wp, 4.36215944530855E-2_wp,-6.21183058238615E-2_wp, &
      & 5.93589247389314E-2_wp, 4.37917586457069E-2_wp, 6.24175529675174E-2_wp, &
      &-9.47993007577910E-2_wp,-9.21933679240633E-2_wp, 1.41162420287618E-1_wp, &
      & 9.47847920056643E-2_wp,-9.22795166591021E-2_wp,-1.41330245385699E-1_wp, &
      & 9.47957155116519E-2_wp, 9.20217761595336E-2_wp,-1.41075167719715E-1_wp, &
      &-9.47612703050034E-2_wp, 9.20970946041004E-2_wp, 1.41174078870679E-1_wp, &
      & 3.07186908919029E-3_wp, 9.78091090173267E-2_wp,-1.59405890565398E-1_wp, &
      &-2.89501267712337E-3_wp, 9.81266227352185E-2_wp, 1.59671635941258E-1_wp, &
      &-3.07984597314511E-3_wp,-9.78272060157212E-2_wp, 1.59285689458381E-1_wp, &
      & 3.45622260398958E-3_wp,-9.78261483786836E-2_wp,-1.59383710436135E-1_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      & 6.56379864949194E-2_wp,-5.73520679256482E-4_wp,-3.52084956587255E-2_wp, &
      & 7.71130683289638E-3_wp,-1.49139358369174E-2_wp,-3.04294908361939E-2_wp, &
      & 6.54814038446620E-2_wp, 5.17836513451199E-4_wp,-3.51478340536698E-2_wp, &
      & 7.70830546882409E-3_wp, 1.49356522560081E-2_wp,-3.03335697909922E-2_wp, &
      & 6.55759132522838E-2_wp,-5.57374828410120E-4_wp,-3.51772230507483E-2_wp, &
      & 7.71886247591620E-3_wp,-1.49101894592223E-2_wp,-3.03986902015354E-2_wp, &
      & 6.44858926966245E-2_wp, 6.21105139273002E-4_wp,-3.60657223374486E-2_wp, &
      & 8.53829013403492E-3_wp, 1.44562000883513E-2_wp,-2.84201703591760E-2_wp, &
      &-2.16992043497083E-1_wp,-6.14525993351003E-2_wp, 2.61869112357575E-1_wp, &
      & 7.43403756029076E-2_wp, 3.29057424398825E-1_wp,-4.48770688604918E-2_wp, &
      &-2.17253885199018E-1_wp, 6.16371296152827E-2_wp, 2.61441601565784E-1_wp, &
      & 7.36098779877098E-2_wp,-3.29050085539625E-1_wp,-4.41877163667672E-2_wp, &
      &-2.17095729603948E-1_wp,-6.13633493424920E-2_wp, 2.61949523181417E-1_wp, &
      & 7.41593953178476E-2_wp, 3.29100144825087E-1_wp,-4.48537935774700E-2_wp, &
      &-2.16773170522421E-1_wp, 6.14170456026348E-2_wp, 2.61763772584928E-1_wp, &
      & 7.44868464063092E-2_wp,-3.29037954066144E-1_wp,-4.49906020625066E-2_wp, &
      &-7.84671583301996E-2_wp,-2.93177074689629E-3_wp, 9.45239916587787E-2_wp, &
      &-1.35552573154701E-3_wp, 1.17958224427535E-1_wp,-1.60568333285798E-2_wp, &
      &-7.86660963425362E-2_wp, 3.27173353036481E-3_wp, 9.46888561547714E-2_wp, &
      &-1.70110485836640E-3_wp,-1.18041188489037E-1_wp,-1.60227598122358E-2_wp, &
      &-7.83224933313823E-2_wp,-2.95614679715344E-3_wp, 9.43219272902749E-2_wp, &
      &-1.39062937509669E-3_wp, 1.17927803395099E-1_wp,-1.59994339588925E-2_wp, &
      &-7.84941989805875E-2_wp, 2.98202100370017E-3_wp, 9.44901433808538E-2_wp, &
      &-1.32436849481419E-3_wp,-1.18022315802278E-1_wp,-1.59959444002657E-2_wp, &
      &-4.89935283276599E-3_wp,-1.35938726283863E-2_wp, 8.75281455685861E-3_wp, &
      & 1.70009780794440E-2_wp, 1.14577178485397E-2_wp,-3.85346172409284E-3_wp, &
      &-4.97349323995833E-3_wp, 1.37313962587995E-2_wp, 8.85061916230923E-3_wp, &
      & 1.71435461297655E-2_wp,-1.12771254648400E-2_wp,-3.87712592235134E-3_wp, &
      &-4.80233859626600E-3_wp,-1.36263939310827E-2_wp, 8.80915005635918E-3_wp, &
      & 1.70027332793215E-2_wp, 1.14775835735027E-2_wp,-4.00681146009285E-3_wp, &
      &-4.69609623172096E-3_wp, 1.35220039889045E-2_wp, 8.69541313382693E-3_wp, &
      & 1.69993038311265E-2_wp,-1.13773019767008E-2_wp,-3.99931690210575E-3_wp],&
      & shape(qpat))

   call get_structure(mol, "X23", "oxacb")
   call test_numsigma(error, mol, qat, dpat, qpat, make_multipole_gfn2, thr_in=thr1)

end subroutine test_s_aes_gfn2_oxacb


subroutine test_s_aes_gxtb_m05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      &  1.00607193048140E-1_wp,  2.93204616982138E-1_wp, -6.31923900996152E-2_wp, &
      &  9.94718617265660E-2_wp,  3.99673325355180E-1_wp, -2.36928769836361E-1_wp, &
      &  4.01761574338139E-2_wp, -2.15912463333496E-1_wp, -4.88276189648645E-1_wp, &
      &  2.69439849576041E-1_wp, -2.47704971849146E-1_wp,  5.11327796912406E-1_wp, &
      & -4.55515768969466E-2_wp, -4.95754223978477E-2_wp,  1.03276773849253E-1_wp, &
      & -4.70035791316972E-1_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      & -2.26079408431879E-1_wp,  1.97005405139997E-1_wp,  3.20065169895204E-1_wp, &
      & -2.51548120328085E-1_wp,  4.62555635243222E-1_wp,  3.43864128150969E-1_wp, &
      & -1.57240378662275E-2_wp,  9.96473595487471E-2_wp, -6.23873255588818E-2_wp, &
      &  1.09366851044351E-1_wp, -1.86697868192343E-2_wp,  1.20675722706818E-2_wp, &
      & -2.85583506192746E-2_wp, -6.23879183695405E-2_wp, -2.65882546839332E-2_wp, &
      & -1.23582189300156E-1_wp, -2.04105836975350E-1_wp,  3.77015851790974E-1_wp, &
      &  5.00232764033768E-2_wp, -8.01616637433991E-2_wp,  1.49873530375689E-2_wp, &
      &  3.80273142974064E-2_wp,  1.54181698060758E-2_wp, -2.05437765038788E-2_wp, &
      & -8.05101442021778E-2_wp,  6.12986400881258E-2_wp, -3.86026859221855E-1_wp, &
      & -3.57889660705829E-2_wp,  7.16292123993904E-2_wp,  3.59105262311269E-2_wp, &
      &  1.87036556910200E-1_wp, -2.36470988156689E-1_wp,  1.13980405816514E-1_wp, &
      & -4.42639769424700E-1_wp,  9.81524608173087E-2_wp, -5.31374900298224E-2_wp, &
      &  3.27349414346133E-2_wp, -1.01759211218765E-1_wp, -6.18311934544862E-2_wp, &
      &  4.01826072161882E-2_wp, -6.84183946605261E-2_wp, -5.52618929308968E-2_wp, &
      & -2.53065567321631E-1_wp, -3.94281607735358E-1_wp, -2.37619483516441E-1_wp, &
      &  4.21483077653447E-1_wp,  1.22005705964453E-1_wp, -1.66074886547056E-1_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      &  1.08336907259083E-1_wp,  4.48319808918437E-1_wp, -5.94348421061186E-2_wp, &
      &  1.94287075957641E-1_wp, -7.75871626677914E-1_wp, -4.89020651529648E-2_wp, &
      &  9.26684307241988E-1_wp,  6.43736381549019E-1_wp, -9.74226568941271E-1_wp, &
      &  1.96960180694922E-1_wp, -3.01079680950474E-1_wp,  4.75422616992827E-2_wp, &
      & -4.51926138090557E-2_wp, -1.55806335788140E-2_wp,  1.02545952895093E-1_wp, &
      & -1.56672483797956E-3_wp, -6.64632119870171E-2_wp, -5.73533390860377E-2_wp, &
      &  1.12030477698873E-1_wp,  1.47957137237649E-2_wp, -4.61056931662970E-2_wp, &
      & -1.69236708071314E-2_wp,  6.16192209975602E-3_wp, -6.59247845325759E-2_wp, &
      &  3.61261056367041E-1_wp, -6.97514645500418E-2_wp, -1.46285367901093E-1_wp, &
      & -5.24458024781472E-1_wp, -3.14532462348105E-1_wp, -2.14975688465946E-1_wp, &
      & -9.69155527556580E-1_wp,  5.51303323868311E-1_wp,  2.26448080549945E-1_wp, &
      & -6.96686590951129E-1_wp,  3.77722833052181E-1_wp,  7.42707447006641E-1_wp, &
      &  2.28419895479698E-2_wp, -3.24221879376675E-2_wp, -2.05626329854465E-2_wp, &
      &  4.25314191702436E-2_wp, -2.93069502336335E-2_wp, -2.27935656252332E-3_wp, &
      &  9.73296769409362E-1_wp,  1.03572902649684E-1_wp, -6.70757165721952E-1_wp, &
      & -6.94649031290375E-1_wp, -7.61364568073428E-2_wp, -3.02539603687412E-1_wp, &
      &  4.17475567765974E-1_wp,  3.95706296656136E-1_wp, -5.49685272886270E-1_wp, &
      &  2.94246648369059E-1_wp,  9.76007015958937E-2_wp,  1.32209705120296E-1_wp, &
      &  1.13151809883176E-2_wp, -3.61994324949722E-2_wp, -1.66474016709828E-2_wp, &
      & -2.67906507458427E-2_wp,  5.48507187214782E-2_wp,  5.33222068266492E-3_wp, &
      &  3.58692422482463E-2_wp,  2.69422126380810E-1_wp,  6.40797051653688E-1_wp, &
      & -3.09040347062161E-1_wp, -8.41402694538926E-1_wp, -6.76666293901935E-1_wp, &
      &  1.05411127417475E+0_wp,  4.94926730507545E-1_wp, -1.34721062336771E+0_wp, &
      &  4.72109954772651E-1_wp, -8.71336612544720E-1_wp,  2.93099349192962E-1_wp, &
      & -1.51089594353926E-3_wp, -5.55591879945522E-2_wp,  6.75028153915662E-2_wp, &
      & -4.92251918013597E-2_wp,  2.77228604581291E-2_wp, -6.59919194480270E-2_wp, &
      &  8.54614367310513E-3_wp, -2.92266352353522E-2_wp, -2.12414707359616E-2_wp, &
      & -4.36968531744766E-2_wp,  7.90924575841896E-2_wp,  1.26953270628563E-2_wp, &
      &  5.78229883864172E-1_wp,  9.13446422751601E-1_wp, -9.75920448407603E-1_wp, &
      & -4.52604908192712E-1_wp,  2.18431458177776E-1_wp,  3.97690564543428E-1_wp, &
      &  4.30681855042516E-1_wp,  4.12587322166438E-1_wp, -4.76264912486020E-2_wp, &
      & -3.13367294588800E-1_wp,  2.25640400422268E-3_wp, -3.83055363793915E-1_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "05")
   call test_numsigma(error, mol, qat, dpat, qpat, make_multipole_gxtb, thr_in=thr1)

end subroutine test_s_aes_gxtb_m05

subroutine test_s_aes_gxtb_m06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(16) = [&
      &  1.37740593302187E-1_wp, -5.81394515563046E-1_wp,  9.31696088480394E-2_wp, &
      & -5.09165955394015E-1_wp, -2.47190537711075E-1_wp, -2.36870944027405E-1_wp, &
      &  7.17339716830641E-1_wp,  9.05519812337205E-2_wp,  2.25155968874018E-1_wp, &
      &  5.02835647961319E-1_wp,  2.38321215929232E-1_wp,  1.14711929684633E-1_wp, &
      & -3.53674184951065E-2_wp, -6.21300098385260E-2_wp, -3.87247005729922E-1_wp, &
      & -6.04602763036139E-2_wp]
   real(wp), parameter :: dpat(3, 16) = reshape([&
      & -1.95522038867900E-1_wp,  3.57473018393676E-2_wp,  4.84009219878046E-2_wp, &
      &  2.26106133034266E-1_wp,  9.40926853359786E-2_wp,  6.14967504010130E-2_wp, &
      &  6.60530642223028E-2_wp,  8.48677921085219E-3_wp, -1.94168533072856E-2_wp, &
      &  2.90051453117600E-2_wp,  3.73104785632608E-2_wp, -8.68785149044501E-3_wp, &
      & -2.10554673919791E-1_wp,  4.90892187629505E-2_wp, -2.01084649960431E-1_wp, &
      &  5.26552978606727E-2_wp, -2.28431763624763E-2_wp, -6.92294693504492E-2_wp, &
      &  5.63870743450493E-3_wp,  2.55785177860538E-1_wp,  1.91283520547011E-1_wp, &
      & -4.91435603305543E-2_wp,  4.12239229763236E-2_wp, -2.70398871634670E-2_wp, &
      & -3.27440520624705E-1_wp, -4.54406545368392E-1_wp,  9.92897855653986E-2_wp, &
      & -1.29257434792820E-1_wp, -1.37580396923045E-1_wp,  5.14517565381086E-1_wp, &
      & -6.38106082350805E-2_wp,  5.36139960379346E-2_wp, -4.88845303443103E-3_wp, &
      &  3.18157317451078E-3_wp, -1.57620411640091E-3_wp,  7.11369499480335E-2_wp, &
      &  1.59087629356785E-2_wp,  1.74397905413834E-2_wp,  7.63798084603338E-2_wp, &
      &  4.87296812675384E-2_wp, -3.73606759574232E-2_wp,  2.86691375487502E-2_wp, &
      & -3.87052261640813E-2_wp, -2.04684951059333E-1_wp, -6.58263602569087E-2_wp, &
      &  3.34427532504582E-2_wp, -2.38478986911650E-2_wp, -1.21602146680676E-2_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 16) = reshape([&
      & -2.44352317377263E-1_wp,  3.54971910078462E-1_wp, -1.53902927754928E-2_wp, &
      &  3.95612341083278E-2_wp,  2.03259033327176E-2_wp,  2.59742610152756E-1_wp, &
      &  3.47816442568286E-1_wp, -4.88800661090738E-3_wp,  5.34538940299081E-2_wp, &
      &  1.37311778795911E-1_wp,  1.39012590712571E-1_wp, -4.01270336598193E-1_wp, &
      &  1.07074309057369E-1_wp,  5.21658697340422E-2_wp, -3.60261607423933E-2_wp, &
      & -7.54797059094252E-2_wp, -1.72898702392894E-2_wp, -7.10481483149761E-2_wp, &
      &  4.00188953492187E-3_wp,  9.14511229882559E-2_wp,  5.46443852193823E-2_wp, &
      & -2.08941515487900E-2_wp, -2.78814828613285E-2_wp, -5.86462747543024E-2_wp, &
      &  2.37857185187833E-1_wp, -4.50300621052163E-2_wp,  2.25131210542704E-1_wp, &
      &  2.35956396879372E-1_wp, -1.02126556733498E-1_wp, -4.62988395730536E-1_wp, &
      & -5.88840912126312E-2_wp,  6.21321317767934E-2_wp,  5.06025800742305E-3_wp, &
      &  2.72035956312882E-2_wp,  4.06725943342423E-2_wp,  5.38238332052079E-2_wp, &
      &  1.32644417792449E-1_wp, -5.92855102419151E-1_wp, -2.40597255998820E-1_wp, &
      &  9.62238671208510E-1_wp,  1.87284920676935E-1_wp,  1.07952838206371E-1_wp, &
      &  1.69665478919738E-2_wp, -9.18320638850835E-2_wp,  3.92297934395233E-2_wp, &
      &  5.88402935507671E-2_wp, -5.24667682572655E-2_wp, -5.61963413314968E-2_wp, &
      &  6.50867568749236E-2_wp, -4.06711260299059E-1_wp, -3.30436331151086E-1_wp, &
      &  7.00082636152224E-2_wp,  1.30647410120592E-1_wp,  2.65349574276161E-1_wp, &
      & -2.63794166402986E-1_wp, -4.74800609208770E-1_wp, -9.09975412069866E-2_wp, &
      & -4.68700474589242E-1_wp, -8.69696406721567E-2_wp,  3.54791707609976E-1_wp, &
      &  4.30949822160176E-2_wp, -5.01918910339194E-2_wp,  1.51853168593317E-2_wp, &
      &  2.14868807692945E-2_wp,  2.32863617853478E-3_wp, -5.82802990753490E-2_wp, &
      & -8.68683448616701E-2_wp, -5.93050256463196E-3_wp, -7.14638312589508E-2_wp, &
      & -3.34700199420226E-3_wp,  2.30797189549555E-2_wp,  1.58332176120621E-1_wp, &
      & -7.77941889351208E-2_wp,  1.90926078158555E-2_wp, -7.22230016477059E-2_wp, &
      &  8.10420605514724E-2_wp,  4.20299203967978E-2_wp,  1.50017190582826E-1_wp, &
      & -5.91145888962679E-2_wp,  1.33207367054315E-1_wp,  1.22628368079946E-1_wp, &
      & -2.86811588062588E-2_wp, -2.71104008133068E-2_wp, -6.35137791836782E-2_wp, &
      & -2.19489556681615E-1_wp,  1.34459386282128E-1_wp,  5.45881480054307E-1_wp, &
      &  1.00045235551687E-1_wp,  2.66486532407144E-1_wp, -3.26391923372692E-1_wp, &
      &  1.31033489381302E-1_wp, -4.91549468099577E-2_wp, -9.00399823418538E-2_wp, &
      & -8.80774028114536E-2_wp,  4.02287690004152E-3_wp, -4.09935070394480E-2_wp],&
      & shape(qpat))

   call get_structure(mol, "MB16-43", "06")
   call test_numsigma(error, mol, qat, dpat, qpat, make_multipole_gxtb, thr_in=thr1)

end subroutine test_s_aes_gxtb_m06

subroutine test_s_aes_gxtb_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(4) = [&
      &  1.76167728179214E+0_wp, -5.86846576518797E-1_wp, -5.88294165118524E-1_wp, &
      & -5.86536540306923E-1_wp]
   real(wp), parameter :: dpat(3, 4) = reshape([&
      &  1.03102227915076E-1_wp,  1.70815542489309E-1_wp, -2.59639799151207E-1_wp, &
      & -1.87602572329602E-1_wp, -2.74737606941817E-1_wp, -1.49657158899311E-1_wp, &
      & -2.03472803217469E-1_wp,  2.74583536877613E-1_wp,  1.60701200271541E-1_wp, &
      &  3.20057753764996E-1_wp, -9.88144341108652E-2_wp,  1.55919293968284E-1_wp],&
      & shape(dpat))
   real(wp), parameter :: qpat(6, 4) = reshape([&
      & -7.51267637005315E-1_wp,  3.31653256315932E-1_wp, -2.26838871462153E-1_wp, &
      & -9.32402326710803E-1_wp, -1.22767701815665E+0_wp,  9.78106508467468E-1_wp, &
      &  2.61814638396274E-3_wp, -2.16516052249151E-2_wp, -9.92760600427545E-3_wp, &
      & -2.24232173276705E-2_wp, -3.20385899497410E-2_wp,  7.30945962031626E-3_wp, &
      & -1.81977682871093E-3_wp,  3.29176108412411E-2_wp, -2.14330340960966E-2_wp, &
      &  5.73574812947845E-3_wp, -3.15197707756747E-2_wp,  2.32528109248094E-2_wp, &
      & -4.08481125421503E-2_wp,  1.54832741157508E-2_wp,  1.99270553528201E-2_wp, &
      & -2.74175062221744E-2_wp, -3.37992023179355E-3_wp,  2.09210571893266E-2_wp],&
      & shape(qpat))

   call get_structure(mol, "f-block", "CeCl3")
   call test_numgrad(error, mol, qat, dpat, qpat, make_multipole_gxtb, thr_in=thr1)

end subroutine test_s_aes_gxtb_cecl3

end module test_coulomb_multipole


submodule (test_coulomb_multipole) test_supercell_scaling
   implicit none


   integer, parameter :: nat = 12
   real(wp), parameter :: qat1(nat) = [&
      &  4.95105332967126E-01_wp,  4.95110445149787E-01_wp,  4.95109208803526E-01_wp, &
      &  4.95110553060372E-01_wp, -2.47570208775520E-01_wp, -2.47372219387123E-01_wp, &
      & -2.47367867558844E-01_wp, -2.47541478966424E-01_wp, -2.47535153228645E-01_wp, &
      & -2.47738212772185E-01_wp, -2.47741337427010E-01_wp, -2.47569061865040E-01_wp]
   real(wp), parameter :: dpat1(3, nat) = reshape([&
      & -3.26654706213550E-04_wp, -6.13403589170031E-05_wp,  2.83929462708564E-04_wp, &
      & -2.68534201975523E-04_wp, -1.57447758402670E-05_wp,  2.68727495879782E-04_wp, &
      & -2.74196985286096E-04_wp, -4.56547709680684E-05_wp,  2.33174149569738E-04_wp, &
      & -3.32314146822759E-04_wp, -6.08326819539977E-08_wp,  3.19473368568051E-04_wp, &
      & -5.42513865726882E-02_wp, -5.43318111661449E-02_wp, -5.44179905961818E-02_wp, &
      & -5.43660117536740E-02_wp, -5.44385268257028E-02_wp,  5.43650907998299E-02_wp, &
      & -5.43539636657199E-02_wp,  5.44469674651156E-02_wp,  5.43590516972356E-02_wp, &
      &  5.44387502680685E-02_wp, -5.43324554152907E-02_wp,  5.42581516232861E-02_wp, &
      &  5.44165886089117E-02_wp,  5.43323030762035E-02_wp,  5.42418925720472E-02_wp, &
      &  5.43211854223397E-02_wp,  5.42458343314803E-02_wp, -5.43133771609431E-02_wp, &
      &  5.43327026728458E-02_wp, -5.42359656240572E-02_wp, -5.43211640102117E-02_wp, &
      & -5.42499802018258E-02_wp,  5.43502803322248E-02_wp, -5.44204207635553E-02_wp],&
      & shape(dpat1))
   real(wp), parameter :: qpat1(6, nat) = reshape([&
      &  1.18066656167315E-05_wp, -4.34798463050927E-01_wp,  1.44462097284581E-07_wp, &
      & -4.34798514781465E-01_wp, -4.34797971420059E-01_wp, -1.19511277154594E-05_wp, &
      &  2.25596069447498E-05_wp, -4.34783690153835E-01_wp, -2.11985471754161E-05_wp, &
      &  4.34768022081199E-01_wp,  4.34785018224248E-01_wp, -1.36105976999978E-06_wp, &
      &  2.25103495342660E-05_wp,  4.34767817555497E-01_wp,  1.07631391967900E-05_wp, &
      &  4.34783902623050E-01_wp, -4.34783120715119E-01_wp, -3.32734887309449E-05_wp, &
      & -9.48473787443227E-06_wp,  4.34783953044164E-01_wp,  1.08281396449250E-05_wp, &
      & -4.34783789399495E-01_wp,  4.34768955922080E-01_wp, -1.34340177049275E-06_wp, &
      & -9.66531291224371E-05_wp, -1.18702391374322E-01_wp, -1.37547479104327E-05_wp, &
      & -1.18645267151075E-01_wp, -1.18600873572414E-01_wp,  1.10407877032648E-04_wp, &
      & -2.26931493848559E-05_wp, -1.18637258774011E-01_wp,  5.06129802909649E-05_wp, &
      &  1.18680813186751E-01_wp,  1.18647169151941E-01_wp, -2.79198309067752E-05_wp, &
      & -3.70191853311663E-05_wp,  1.18641364188384E-01_wp,  6.65002198458886E-05_wp, &
      &  1.18692897340129E-01_wp, -1.18644737187210E-01_wp, -2.94810345136121E-05_wp, &
      &  1.06241518233796E-04_wp,  1.18591745614562E-01_wp,  9.27256724547743E-06_wp, &
      & -1.18653931140339E-01_wp,  1.18709157591422E-01_wp, -1.15514085479163E-04_wp, &
      &  9.75875617358346E-05_wp, -1.18606840689317E-01_wp,  1.39412337278877E-05_wp, &
      & -1.18668677913009E-01_wp, -1.18717811830031E-01_wp, -1.11528795463611E-04_wp, &
      &  3.00943076451121E-05_wp, -1.18658461646011E-01_wp, -6.19964990261623E-05_wp, &
      &  1.18611269561002E-01_wp,  1.18657681228184E-01_wp,  3.19021913810502E-05_wp, &
      &  4.34882210389453E-05_wp,  1.18645813671023E-01_wp, -6.06140176001579E-05_wp, &
      &  1.18607408054076E-01_wp, -1.18660243635834E-01_wp,  1.71257965608795E-05_wp, &
      & -1.15717831159601E-04_wp,  1.18704043515894E-01_wp, -3.90833585128814E-06_wp, &
      & -1.18646463864190E-01_wp,  1.18587233297690E-01_wp,  1.19626167010667E-04_wp],&
      & shape(qpat1))

   real(wp), parameter :: e02 = 1.5016169607148633E-2_wp, e11 = -3.5486703953320990E-3_wp, &
      & e01 = 1.5706712259423678E-2_wp - e02 - e11

contains

module subroutine test_e_aes_gfn2_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = qat1
   real(wp), parameter :: dpat(3, nat) = dpat1
   real(wp), parameter :: qpat(6, nat) = qpat1

   call get_structure(mol, "X23", "CO2")
   call test_generic(error, mol, qat, dpat, qpat, make_multipole_gfn2, e01+e02+e11)
end subroutine test_e_aes_gfn2_co2

module subroutine test_e_aes_gfn2_co2_dp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat0(nat) = 0.0_wp
   real(wp), parameter :: dpat(3, nat) = dpat1
   real(wp), parameter :: qpat0(6, nat) = 0.0_wp

   call get_structure(mol, "X23", "CO2")
   call test_generic(error, mol, qat0, dpat, qpat0, make_multipole_gfn2, e11)
end subroutine test_e_aes_gfn2_co2_dp

module subroutine test_e_aes_gfn2_co2_qp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(nat) = qat1
   real(wp), parameter :: dpat0(3, nat) = 0.0_wp
   real(wp), parameter :: qpat(6, nat) = qpat1

   call get_structure(mol, "X23", "CO2")
   call test_generic(error, mol, qat, dpat0, qpat, make_multipole_gfn2, e02)
end subroutine test_e_aes_gfn2_co2_qp


module subroutine test_e_aes_gfn2_co2_sc(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   integer, parameter :: supercell(3) = [2, 2, 2]
   integer, parameter :: nsc = product(supercell)
   real(wp), parameter :: qat(nsc*nat) = [spread(qat1, 2, nsc)]
   real(wp), parameter :: dpat(3, nsc*nat) = &
      & reshape([dpat1, dpat1, dpat1, dpat1, dpat1, dpat1, dpat1, dpat1], shape(dpat))
   real(wp), parameter :: qpat(6, nsc*nat) = &
      & reshape([qpat1, qpat1, qpat1, qpat1, qpat1, qpat1, qpat1, qpat1], shape(qpat))

   call get_structure(mol, "X23", "CO2")
   call make_supercell(mol, supercell)
   call test_generic(error, mol, qat, dpat, qpat, make_multipole_gfn2, product(supercell)*(e01+e02+e11))
end subroutine test_e_aes_gfn2_co2_sc


module subroutine test_e_aes_gfn2_co2_sc_dp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   integer, parameter :: supercell(3) = [2, 2, 2]
   integer, parameter :: nsc = product(supercell)
   real(wp), parameter :: qat0(nsc*nat) = 0.0_wp
   real(wp), parameter :: dpat(3, nsc*nat) = &
      & reshape([dpat1, dpat1, dpat1, dpat1, dpat1, dpat1, dpat1, dpat1], shape(dpat))
   real(wp), parameter :: qpat0(6, nsc*nat) = 0.0_wp

   call get_structure(mol, "X23", "CO2")
   call make_supercell(mol, supercell)
   call test_generic(error, mol, qat0, dpat, qpat0, make_multipole_gfn2, product(supercell)*e11)
end subroutine test_e_aes_gfn2_co2_sc_dp


module subroutine test_e_aes_gfn2_co2_sc_qp(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   integer, parameter :: supercell(3) = [2, 2, 2]
   integer, parameter :: nsc = product(supercell)
   real(wp), parameter :: qat(nsc*nat) = [spread(qat1, 2, nsc)]
   real(wp), parameter :: dpat0(3, nsc*nat) = 0.0_wp
   real(wp), parameter :: qpat(6, nsc*nat) = &
      & reshape([qpat1, qpat1, qpat1, qpat1, qpat1, qpat1, qpat1, qpat1], shape(qpat))

   call get_structure(mol, "X23", "CO2")
   call make_supercell(mol, supercell)
   call test_generic(error, mol, qat, dpat0, qpat, make_multipole_gfn2, product(supercell)*e02)
end subroutine test_e_aes_gfn2_co2_sc_qp


subroutine make_supercell(mol, rep)
   use mctc_io, only : new
   type(structure_type), intent(inout) :: mol
   integer, intent(in) :: rep(3)

   real(wp), allocatable :: xyz(:, :), lattice(:, :)
   integer, allocatable :: num(:)
   integer :: i, j, k, c

   num = reshape(spread([mol%num(mol%id)], 2, product(rep)), [product(rep)*mol%nat])
   lattice = reshape(&
      [rep(1)*mol%lattice(:, 1), rep(2)*mol%lattice(:, 2), rep(3)*mol%lattice(:, 3)], &
      shape(mol%lattice))
   allocate(xyz(3, product(rep)*mol%nat))
   c = 0
   do i = 0, rep(1)-1
      do j = 0, rep(2)-1
         do k = 0, rep(3)-1
            xyz(:, c+1:c+mol%nat) = mol%xyz &
               & + spread(matmul(mol%lattice, [real(wp):: i, j, k]), 2, mol%nat)
            c = c + mol%nat
         end do
      end do
   end do

   call new(mol, num, xyz, lattice=lattice)
end subroutine make_supercell


end submodule
