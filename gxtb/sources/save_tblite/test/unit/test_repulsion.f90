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

module test_repulsion
   use mctc_data_vdwrad, only : get_vdw_rad
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mctc_io_convert, only : autoaa
   use mstore, only : get_structure
   use multicharge, only : get_eeqbc_charges
   use tblite_basis_type, only : basis_type
   use tblite_blas, only : gemv
   use tblite_container_cache, only : container_cache
   use tblite_repulsion, only : repulsion_type, gfn_repulsion, gxtb_repulsion, &
      & new_repulsion_gfn, new_repulsion_gxtb
   use tblite_scf_potential, only : potential_type, new_potential
   use tblite_utils_average, only : average_type, new_average, average_id
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   implicit none
   private

   public :: collect_repulsion

   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine repulsion_maker(rep, mol, error)
         import :: repulsion_type, structure_type, error_type
         class(repulsion_type), allocatable, intent(out) :: rep
         type(structure_type), intent(in) :: mol
         type(error_type), allocatable, intent(out) :: error
      end subroutine repulsion_maker
   end interface

contains


!> Collect all exported unit tests
subroutine collect_repulsion(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-gfn1-m01", test_e_gfn1_m01), &
      new_unittest("energy-gfn2-m02", test_e_gfn2_m02), &
      new_unittest("energy-gxtb-m03", test_e_gxtb_m03), &
      new_unittest("energy-gxtb-m04", test_e_gxtb_m04), &
      new_unittest("energy-pbc-gfn2-uracil", test_e_gfn2_uracil), &
      new_unittest("potential-gxtb-m01", test_p_gxtb_m01), &
      new_unittest("potential-gxtb-m02", test_p_gxtb_m02), &
      new_unittest("gradient-gfn1-m03", test_g_gfn1_m03), &
      new_unittest("gradient-gfn2-m04", test_g_gfn2_m04), &
      new_unittest("gradient-gxtb-m05", test_g_gxtb_m05), &
      new_unittest("gradient-pbc-gfn2-urea", test_g_gfn2_urea), &
      new_unittest("sigma-gfn1-m05", test_s_gfn1_m05), &
      new_unittest("sigma-gfn2-m06", test_s_gfn2_m06), &
      new_unittest("sigma-gxtb-m07", test_s_gxtb_m07), &
      new_unittest("sigma-pbc-gfn2-succinic", test_s_gfn2_succinic) &
      ]

end subroutine collect_repulsion


!> Factory to create repulsion objects based on GFN1-xTB
subroutine make_repulsion_gfn1(rep, mol, error)
   !> Repulsion object
   class(repulsion_type), allocatable, intent(out) :: rep
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: alpha_gfn1(1:20) = [&
      & 2.209700_wp, 1.382907_wp, 0.671797_wp, 0.865377_wp, 1.093544_wp, &
      & 1.281954_wp, 1.727773_wp, 2.004253_wp, 2.507078_wp, 3.038727_wp, &
      & 0.704472_wp, 0.862629_wp, 0.929219_wp, 0.948165_wp, 1.067197_wp, &
      & 1.200803_wp, 1.404155_wp, 1.323756_wp, 0.581529_wp, 0.665588_wp]
   real(wp), parameter :: zeff_gfn1(1:20) = [&
      &  1.116244_wp,  0.440231_wp,  2.747587_wp,  4.076830_wp,  4.458376_wp, &
      &  4.428763_wp,  5.498808_wp,  5.171786_wp,  6.931741_wp,  9.102523_wp, &
      & 10.591259_wp, 15.238107_wp, 16.283595_wp, 16.898359_wp, 15.249559_wp, &
      & 15.100323_wp, 17.000000_wp, 17.153132_wp, 20.831436_wp, 19.840212_wp]
   real(wp), allocatable :: alpha(:), zeff(:)
   real(wp) :: kexp, rexp

   type(gfn_repulsion), allocatable :: tmp

   alpha = alpha_gfn1(mol%num)
   zeff = zeff_gfn1(mol%num)
   kexp = 1.5_wp
   rexp = 1.0_wp
   allocate(tmp)
   call new_repulsion_gfn(tmp, mol, alpha=alpha, zeff=zeff, &
      & kexp=kexp, kexp_light=kexp, rexp=rexp)
   call move_alloc(tmp, rep)
end subroutine make_repulsion_gfn1

!> Factory to create repulsion objects based on GFN2-xTB
subroutine make_repulsion_gfn2(rep, mol, error)
   !> Repulsion object
   class(repulsion_type), allocatable, intent(out) :: rep
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: alpha_gfn2(1:20) = [&
      & 2.213717_wp, 3.604670_wp, 0.475307_wp, 0.939696_wp, 1.373856_wp, &
      & 1.247655_wp, 1.682689_wp, 2.165712_wp, 2.421394_wp, 3.318479_wp, &
      & 0.572728_wp, 0.917975_wp, 0.876623_wp, 1.187323_wp, 1.143343_wp, &
      & 1.214553_wp, 1.577144_wp, 0.896198_wp, 0.482206_wp, 0.683051_wp]
   real(wp), parameter :: zeff_gfn2(1:20) = [&
      &  1.105388_wp,  1.094283_wp,  1.289367_wp,  4.221216_wp,  7.192431_wp, &
      &  4.231078_wp,  5.242592_wp,  5.784415_wp,  7.021486_wp, 11.041068_wp, &
      &  5.244917_wp, 18.083164_wp, 17.867328_wp, 40.001111_wp, 19.683502_wp, &
      & 14.995090_wp, 17.353134_wp,  7.266606_wp, 10.439482_wp, 14.786701_wp]
   real(wp), allocatable :: alpha(:), zeff(:)
   real(wp) :: kexp, kexp_light, rexp

   type(gfn_repulsion), allocatable :: tmp

   alpha = alpha_gfn2(mol%num)
   zeff = zeff_gfn2(mol%num)
   kexp = 1.5_wp
   kexp_light = 1.0_wp
   rexp = 1.0_wp
   allocate(tmp)
   call new_repulsion_gfn(tmp, mol, alpha=alpha, zeff=zeff, &
      & kexp=kexp, kexp_light=kexp_light, rexp=rexp)
   call move_alloc(tmp, rep)
end subroutine make_repulsion_gfn2

!> Factory to create repulsion objects based on g-xTB
subroutine make_repulsion_gxtb(rep, mol, error)
   !> Repulsion object
   class(repulsion_type), allocatable, intent(out) :: rep
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Damping exponent in g-xTB repulsion
   real(wp), parameter :: alpha_gxtb(1:20) = [&
      &  2.3064151254_wp,  1.7009718420_wp,  1.0458011914_wp,  1.4067559016_wp, & !1-4
      &  1.3617823162_wp,  1.3504447002_wp,  1.5681331528_wp,  1.6076347346_wp, & !5-8
      &  1.9750426600_wp,  1.7245863480_wp,  1.0327533208_wp,  1.0099829332_wp, & !9-12
      &  1.0348717084_wp,  0.9990576652_wp,  0.9818397342_wp,  1.0908412780_wp, & !13-16
      &  1.2017500102_wp,  1.3421422540_wp,  0.9233240796_wp,  0.9379985182_wp] !17-20

   !> Environment-independent effective nuclear charge in g-xTB repulsion
   real(wp), parameter :: zeff_gxtb(1:20) = [&
      &  1.1651363510_wp,  1.6147404156_wp,  2.4088644573_wp,  2.9910757114_wp, & !1-4
      &  2.9653240876_wp,  3.0013735364_wp,  3.1702771654_wp,  2.7937777512_wp, & !5-8
      &  2.4618005088_wp,  2.8424842044_wp,  5.3104161736_wp,  3.7135952381_wp, & !9-12
      &  4.0415845090_wp,  3.6542325309_wp,  3.2819133459_wp,  3.4919910081_wp, & !13-16
      &  3.7183123254_wp,  4.7127657589_wp,  8.2615266398_wp,  5.6725795921_wp] !17-20

   !> CN-dependence of the damping exponent in g-xTB repulsion
   real(wp), parameter :: kcn_gxtb(1:20) = [&
      & -0.0531679812_wp,  0.1522345594_wp,  0.0670220679_wp,  0.0070419276_wp, & !1-4
      &  0.1269377930_wp,  0.0005415406_wp, -0.0157761568_wp, -0.0404982270_wp, & !5-8
      & -0.1192200271_wp, -0.3759724092_wp,  0.1098475549_wp,  0.1136516455_wp, & !9-12
      &  0.0056233749_wp,  0.0501531120_wp,  0.0327507128_wp, -0.0791156425_wp, & !13-16
      & -0.0543143398_wp,  0.1107306518_wp,  0.1542241501_wp,  0.1113109564_wp] !17-20

   !> Charge-dependence of the effective nuclear charge in g-xTB repulsion
   real(wp), parameter :: kq_gxtb(1:20) = [&
      & -0.0081209590_wp,  0.2790638838_wp, -0.0829596185_wp, -0.0183031721_wp, & !1-4
      & -0.0146182371_wp,  0.0030745124_wp,  0.0011233853_wp, -0.0193806701_wp, & !5-8
      &  0.1838885885_wp,  0.0973313826_wp, -0.1688817066_wp, -0.0104042318_wp, & !9-12
      &  0.0041519279_wp,  0.0147646445_wp, -0.0020241599_wp,  0.0429527086_wp, & !13-16
      & -0.0670716286_wp,  0.1211869919_wp, -0.1283368421_wp, -0.0100450923_wp] !17-20

   !> Offset atomic radii in g-xTB repulsion
   real(wp), parameter :: roffset_gxtb(1:20) = [&
      &  0.2605095070_wp,  0.9302084918_wp,  0.2311088371_wp,  0.3075571987_wp, & !1-4
      &  0.3593779508_wp,  0.4976308757_wp,  0.4809703425_wp,  0.4672401499_wp, & !5-8
      &  0.2998780561_wp,  0.3330185987_wp,  0.3428410987_wp,  0.2041348539_wp, & !9-12
      &  0.4459445436_wp,  0.3481122618_wp,  0.3642906897_wp,  0.3969689179_wp, & !13-16
      &  0.3606684818_wp,  0.1331283171_wp,  0.2093700751_wp,  0.2705625792_wp] !17-20

   !> First-order expansion coeffient in g-xTB repulsion
   real(wp), parameter :: k1_gxtb(1:20) = [&
      & -0.3590022732_wp,  3.2448304666_wp, -0.1000972069_wp,  0.0747194459_wp, & !1-4
      &  0.0827507883_wp,  0.1159291480_wp,  0.0096119015_wp,  0.0557337069_wp, & !5-8
      &  0.1311191951_wp,  0.2715035680_wp, -1.7977844176_wp,  0.9161847112_wp, & !9-12
      &  0.6137328343_wp,  0.1242287340_wp,  0.1697366618_wp,  0.0939307333_wp, & !13-16
      &  0.3167667611_wp,  1.4014300137_wp, -1.4845782485_wp,  0.4494744288_wp] !17-20

   !> Radii in the coordination number in g-xTB
   real(wp), parameter :: cn_rcov_gxtb(1:20) = 0.5_wp * [&
      &  0.9646730831_wp,  1.0887413026_wp,  2.7319904126_wp,  2.0642646118_wp, & !1-4
      &  1.8403302383_wp,  1.8469537700_wp,  1.7342318861_wp,  1.6797510318_wp, & !5-8
      &  1.2145258133_wp,  1.6166175254_wp,  3.1780229637_wp,  2.8181820446_wp, & !9-12
      &  2.6406018705_wp,  2.5882455353_wp,  2.6392688502_wp,  2.5174870894_wp, & !13-16
      &  2.1965435815_wp,  2.3352396673_wp,  3.8631034402_wp,  3.2237828672_wp] !17-20

   !> van-der-Waals radii scaling
   real(wp), parameter :: rvdw_scale(1:20) = autoaa * [&
      &  1.0000000000_wp,  1.0000000000_wp,  1.0600000000_wp,  1.0000000000_wp, & !1-4
      &  1.0200000000_wp,  1.0000000000_wp,  1.0500000000_wp,  1.0800000000_wp, & !5-8
      &  1.1500000000_wp,  0.8000000000_wp,  1.2000000000_wp,  1.1500000000_wp, & !9-12
      &  1.0700000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !13-16
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp] !17-20

   !> Coordination number exponent
   real(wp), parameter :: p_cn_exp = 2.0680000000_wp

   !> Distance damping exponent repulsion
   real(wp), parameter :: p_kexp = 1.5000000000_wp

   !> Expansion coefficients repulsion
   real(wp), parameter :: p_k2_light = 0.1704681501_wp / 3.75372519958925_wp**2
   real(wp), parameter :: p_k2 = 0.1203923899_wp / 3.75372519958925_wp**2
   real(wp), parameter :: p_k3 = 0.5868618150_wp / 3.75372519958925_wp**3
   real(wp), parameter :: p_k4 = 2.3046216529_wp / 3.75372519958925_wp**4

   !> Short-range-repulsion correction prefactor, 
   !> damping exponent scaling and distance damping exponent 
   real(wp), parameter :: p_short = 0.0046511298_wp
   real(wp), parameter :: p_short_alpha = 0.7300000000_wp
   real(wp), parameter :: p_short_exp = 2.0000000000_wp

   integer :: isp, jsp, izp, jzp
   real(wp), allocatable :: alpha(:), zeff(:), roffset(:), kq(:), k1(:)
   real(wp), allocatable :: kcn(:), cn_rcov(:), rvdw(:, :)
   type(average_type) :: average

   type(gxtb_repulsion), allocatable :: tmp

   alpha = alpha_gxtb(mol%num)
   zeff = zeff_gxtb(mol%num)
   kcn = kcn_gxtb(mol%num)
   kq = kq_gxtb(mol%num)
   roffset = roffset_gxtb(mol%num)
   k1 = k1_gxtb(mol%num)
   cn_rcov = cn_rcov_gxtb(mol%num)

   call new_average(average, average_id%arithmetic)
   allocate(rvdw(mol%nid, mol%nid)) 
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do jsp = 1, isp
         jzp = mol%num(jsp)
         rvdw(jsp, isp) = get_vdw_rad(jzp, izp) &
            & * average%value(rvdw_scale(izp), rvdw_scale(jzp))
         rvdw(isp, jsp) = rvdw(jsp, isp)
      end do
   end do

   allocate(tmp)
   call new_repulsion_gxtb(tmp, mol, alpha=alpha, zeff=zeff, kcn=kcn, kq=kq, &
      & roffset=roffset, kexp=p_kexp, k1=k1, k2=p_k2, k2_light=p_k2_light, &
      & k3=p_k3, k4=p_k4, kshort=p_short, kshort_alpha=p_short_alpha, &
      & kshort_exp=p_short_exp, rvdw=rvdw, cn_rcov=cn_rcov, cn_exp=p_cn_exp, &
      & error=error)
   if (allocated(error)) return
   call move_alloc(tmp, rep)

end subroutine make_repulsion_gxtb

subroutine test_generic(error, mol, make_repulsion, ref, qat, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Factory to produce repulsion objects
   procedure(repulsion_maker) :: make_repulsion

   !> Reference value to compare against
   real(wp), intent(in) :: ref

   !> Atomic partial charges for this structure
   real(wp), intent(in), optional :: qat(:)

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   class(repulsion_type), allocatable :: rep
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   real(wp) :: energy(mol%nat)
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   energy = 0.0_wp
   call make_repulsion(rep, mol, error)
   if(allocated(error)) return
   call rep%update(mol, cache)
   call rep%get_engrad(mol, cache, energy)
   if (present(qat)) then
      wfn%qat = reshape(qat, [size(qat), 1])
      call rep%get_energy(mol, cache, wfn, energy)
   end if

   call check(error, sum(energy), ref, thr=thr_)
   if (allocated(error)) then
      print*,sum(energy), ref, sum(energy) - ref
   end if

end subroutine test_generic


subroutine test_numpot(error, mol, make_repulsion, qat, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Factory to create a repulsion object
   procedure(repulsion_maker) :: make_repulsion

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ish
   type(basis_type) :: bas
   class(repulsion_type), allocatable :: rep
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   real(wp), allocatable :: numvat(:), numvsh(:)
   real(wp) :: er(mol%nat), el(mol%nat)
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   wfn%qat = reshape(qat, [size(qat), 1])
   bas%nsh = size(qat)
   bas%nao = size(qat)

   ! Setup potential with dummy basis set
   call new_potential(pot, mol, bas, 1, .false.)

   ! Setup repulsion object
   call make_repulsion(rep, mol, error)
   if(allocated(error)) return

   ! Update repulsion cache
   call rep%update(mol, cache)

   ! Numerical atomic potential
   allocate(numvat(mol%nat), source=0.0_wp)
   do iat = 1, mol%nat
      er = 0.0_wp
      el = 0.0_wp
      ! Right hand side
      wfn%qat(iat, 1) = wfn%qat(iat, 1) + step
      call rep%get_energy(mol, cache, wfn, er)
      ! Left hand side
      wfn%qat(iat, 1) = wfn%qat(iat, 1) - 2*step
      call rep%get_energy(mol, cache, wfn, el)

      wfn%qat(iat, 1) = wfn%qat(iat, 1) + step
      numvat(iat) = 0.5_wp*(sum(er) - sum(el))/step
   end do

   ! Analytic potentials
   call pot%reset()
   call rep%get_potential(mol, cache, wfn, pot)

   if (any(abs(pot%vat(:, 1) - numvat) > thr_)) then
      call test_failed(error, "Atom-resolved potential does not match")
      write(*,*) "numerical potential:"
      print'(3es21.14)', numvat
      write(*,*) "analytical potential:"
      print'(3es21.14)', pot%vat(: ,1)
      write(*,*) "difference:"
      print'(3es21.14)', pot%vat(: ,1) - numvat
   end if

end subroutine test_numpot


subroutine test_numgrad(error, mol, make_repulsion, qat, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Factory to create a repulsion object
   procedure(repulsion_maker) :: make_repulsion

   !> Atomic partial charges for this structure
   real(wp), intent(in), optional :: qat(:)

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic
   class(repulsion_type), allocatable :: rep
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat), sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call make_repulsion(rep, mol, error)
   if(allocated(error)) return

   if (present(qat)) then
      wfn%qat = reshape(qat, [size(qat), 1])
   end if

   do iat = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp

         ! Right hand
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         ! Update and calculate repulsion
         call rep%update(mol, cache)
         call rep%get_engrad(mol, cache, er)
         if (present(qat)) then
            call rep%get_energy(mol, cache, wfn, er)
         end if
         
         ! Left hand
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         ! Update and calculate repulsion
         call rep%update(mol, cache)
         call rep%get_engrad(mol, cache, el)
         if (present(qat)) then
            call rep%get_energy(mol, cache, wfn, el)
         end if
      
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call rep%update(mol, cache)
   call rep%get_engrad(mol, cache, energy, gradient, sigma)
   if (present(qat)) then
      call rep%get_gradient(mol, cache, wfn, gradient, sigma)
   end if

   if (any(abs(gradient - numgrad) > thr_)) then
      call test_failed(error, "Gradient of repulsion energy does not match")
      write(*,*) "numerical gradient:"
      print'(3es21.14)', numgrad
      write(*,*) "analytical gradient:"
      print'(3es21.14)', gradient
      write(*,*) "difference:"
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_numgrad


subroutine test_numsigma(error, mol, make_repulsion, qat, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Factory to create a repulsion object
   procedure(repulsion_maker) :: make_repulsion

   !> Atomic partial charges for this structure
   real(wp), intent(in), optional :: qat(:)

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: ic, jc
   class(repulsion_type), allocatable :: rep
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)
   real(wp) :: sigma(3, 3), eps(3, 3), numsigma(3, 3)
   real(wp), allocatable :: gradient(:, :), xyz(:, :), lat(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate(gradient(3, mol%nat), xyz(3, mol%nat), lat(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call make_repulsion(rep, mol, error)
   if(allocated(error)) return
   
   if (present(qat)) then
      wfn%qat = reshape(qat, [size(qat), 1])
   end if

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   if (any(mol%periodic)) &
   lat(:, :) = mol%lattice
   do ic = 1, 3
      do jc = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         
         ! Right hand
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (any(mol%periodic)) mol%lattice(:, :) = matmul(eps, lat)
         ! Update and calculate repulsion
         call rep%update(mol, cache)
         call rep%get_engrad(mol, cache, er)
         if (present(qat)) then
            call rep%get_energy(mol, cache, wfn, er)
         end if

         ! Left hand
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (any(mol%periodic)) mol%lattice(:, :) = matmul(eps, lat)
         ! Update and calculate repulsion
         call rep%update(mol, cache)
         call rep%get_engrad(mol, cache, el)
         if (present(qat)) then
            call rep%get_energy(mol, cache, wfn, el)
         end if

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (any(mol%periodic)) mol%lattice(:, :) = lat
         numsigma(jc, ic) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call rep%update(mol, cache)
   call rep%get_engrad(mol, cache, energy, gradient, sigma)
   if (present(qat)) then
      call rep%get_gradient(mol, cache, wfn, gradient, sigma)
   end if

   if (any(abs(sigma - numsigma) > thr_)) then
      call test_failed(error, "Strain derivatives do not match")
      write(*,*) "numerical sigma:"
      print'(3es21.14)', numsigma
      write(*,*) "analytical sigma:"
      print'(3es21.14)', sigma
      write(*,*) "difference:"
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_numsigma


subroutine test_e_gfn1_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, make_repulsion_gfn1, 0.16777923624986593_wp)

end subroutine test_e_gfn1_m01

subroutine test_e_gfn2_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, make_repulsion_gfn2, 0.10745931926703985_wp)

end subroutine test_e_gfn2_m02

subroutine test_e_gxtb_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qat(16) = [ &
      & -6.722325895098202E-2_wp, -9.343346477695722E-1_wp,  7.423184694799456E-2_wp, &
      &  6.003342781579326E-1_wp,  8.155559152718005E-1_wp,  1.119087497546364E-1_wp, &
      & -4.000789327736840E-1_wp,  8.016756865222574E-2_wp,  8.005299562801238E-2_wp, &
      &  1.006389816802352E-1_wp, -6.704006636400752E-1_wp, -4.326713538195071E-1_wp, &
      & -5.263477368917056E-1_wp,  2.692434056872415E-1_wp,  8.190366482982609E-1_wp, &
      &  7.988620327393425E-2_wp]

   call get_structure(mol, "MB16-43", "03")
   call test_generic(error, mol, make_repulsion_gxtb, 1.83264694725369_wp, qat, thr_in=thr1)

end subroutine test_e_gxtb_m03

subroutine test_e_gxtb_m04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qat(16) = [ &
      & -6.417476509346409E-2_wp, -2.256661931683622E-1_wp, -3.458626714493751E-2_wp, &
      & -3.860600924041810E-1_wp,  7.238991472693167E-1_wp, -5.777084491650908E-2_wp, &
      &  1.196120411671119E-1_wp, -5.461117600438992E-3_wp, -1.794514383018986E-2_wp, &
      & -1.211514276086560E-1_wp, -3.493736689115619E-1_wp,  9.244684713156893E-1_wp, &
      &  1.394030055851094E-1_wp, -6.975412178494735E-1_wp,  1.149356570285273E-1_wp, &
      & -6.258758429064248E-2_wp]

   call get_structure(mol, "MB16-43", "04")
   call test_generic(error, mol, make_repulsion_gxtb, 1.91678411629162_wp, qat, thr_in=thr1)

end subroutine test_e_gxtb_m04

subroutine test_e_gfn2_uracil(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "uracil")
   call test_generic(error, mol, make_repulsion_gfn2, 1.0401472262740301_wp)

end subroutine test_e_gfn2_uracil


subroutine test_p_gxtb_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qat(16) = [ &
      &  9.115222032717870E-1_wp, -7.820978364181719E-2_wp, -8.085043301444808E-1_wp, &
      & -1.370008894958108E-1_wp, -4.522670977383800E-1_wp,  2.366794200051805E-1_wp, &
      & -1.056480024872486E-1_wp, -7.542974658029240E-1_wp, -6.084417930861774E-1_wp, &
      &  3.252603520948056E-1_wp,  2.590209057852254E-1_wp, -3.245713596482531E-1_wp, &
      &  9.588528411034369E-2_wp,  7.104890489542336E-1_wp, -3.468331265256146E-1_wp, &
      &  1.076916633846320E+0_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_numpot(error, mol, make_repulsion_gxtb, qat, thr_in=thr1*10)

end subroutine test_p_gxtb_m01

subroutine test_p_gxtb_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qat(16) = [ &
      & -1.097288775394150E-1_wp, -4.735023909155970E-1_wp, -1.707389583649015E-1_wp, &
      & -6.433491759460952E-1_wp,  8.965254778844154E-1_wp,  3.351267282573084E-1_wp, &
      & -6.657210518576462E-2_wp, -5.047599279669956E-2_wp,  5.146699702222322E-1_wp, &
      &  3.455270064324552E-1_wp,  8.365308656920356E-2_wp,  7.840093894994580E-1_wp, &
      & -6.236807887791764E-1_wp, -7.429581513991601E-2_wp, -6.465871716055793E-2_wp, &
      & -6.825088374693600E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_numpot(error, mol, make_repulsion_gxtb, qat, thr_in=thr1*10)

end subroutine test_p_gxtb_m02


subroutine test_g_gfn1_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   call test_numgrad(error, mol, make_repulsion_gfn1, thr_in=thr1*10)

end subroutine test_g_gfn1_m03

subroutine test_g_gfn2_m04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol, make_repulsion_gfn2, thr_in=thr1*10)

end subroutine test_g_gfn2_m04

subroutine test_g_gxtb_m05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: qat(16) = [ &
      &  1.006071930481401E-1_wp,  2.932046169821381E-1_wp, -6.319239009961519E-2_wp, &
      &  9.947186172656597E-2_wp,  3.996733253551799E-1_wp, -2.369287698363611E-1_wp, &
      &  4.017615743381386E-2_wp, -2.159124633334960E-1_wp, -4.882761896486452E-1_wp, &
      &  2.694398495760408E-1_wp, -2.477049718491459E-1_wp,  5.113277969124058E-1_wp, &
      & -4.555157689694656E-2_wp, -4.957542239784773E-2_wp,  1.032767738492532E-1_wp, &
      & -4.700357913169724E-1_wp]

   call get_structure(mol, "MB16-43", "05")
   call test_numgrad(error, mol, make_repulsion_gxtb, qat, thr_in=thr1*10)

end subroutine test_g_gxtb_m05

subroutine test_g_gfn2_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "urea")
   call test_numgrad(error, mol, make_repulsion_gfn2, thr_in=thr1*10)

end subroutine test_g_gfn2_urea


subroutine test_s_gfn1_m05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numsigma(error, mol, make_repulsion_gfn1, thr_in=thr1*100)

end subroutine test_s_gfn1_m05

subroutine test_s_gfn2_m06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "06")
   call test_numsigma(error, mol, make_repulsion_gfn2, thr_in=thr1*100)

end subroutine test_s_gfn2_m06

subroutine test_s_gxtb_m07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   
   real(wp), parameter :: qat(16) = [ &
      & -3.110138689081248E-1_wp,  2.172798226943077E-1_wp,  3.780450108186539E-1_wp, &
      &  1.061525247178133E-1_wp,  1.177528337612729E-1_wp, -4.369962228498777E-1_wp, &
      & -3.739788923468512E-1_wp, -6.568776993798167E-1_wp,  4.608781463362221E-1_wp, &
      &  1.097689912078816E-1_wp, -4.782088598388384E-1_wp,  3.672206075537463E-1_wp, &
      &  3.475292477131291E-1_wp, -5.850511243090697E-1_wp, -2.817007347458540E-1_wp, &
      &  1.019200217084888E+0_wp]

   call get_structure(mol, "MB16-43", "07")
   call test_numsigma(error, mol, make_repulsion_gxtb, qat, thr_in=thr1*100)

end subroutine test_s_gxtb_m07

subroutine test_s_gfn2_succinic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "succinic")
   call test_numsigma(error, mol, make_repulsion_gfn2, thr_in=thr1*100)

end subroutine test_s_gfn2_succinic

end module test_repulsion