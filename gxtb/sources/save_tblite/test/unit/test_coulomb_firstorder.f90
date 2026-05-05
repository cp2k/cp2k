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

module test_coulomb_firstorder
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use mstore, only : get_structure
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_basis_type
   use tblite_ceh_ceh, only : get_effective_qat
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_firstorder, only : onsite_firstorder, new_onsite_firstorder
   use tblite_cutoff, only : get_lattice_points
   use tblite_scf, only: new_potential, potential_type
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_coulomb, only : tb_coulomb, new_coulomb
   implicit none
   private

   public :: collect_coulomb_firstorder

   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine coulomb_maker(coulomb, mol, shell, error)
         import :: tb_coulomb, structure_type, error_type
         type(tb_coulomb), allocatable, intent(out) :: coulomb
         type(structure_type), intent(in) :: mol
         logical, intent(in) :: shell
         type(error_type), allocatable, intent(out) :: error
      end subroutine coulomb_maker
   end interface

   abstract interface
      subroutine charge_maker(wfn, mol, nshell, error)
         import :: wavefunction_type, structure_type, error_type
         type(wavefunction_type), intent(inout) :: wfn
         type(structure_type), intent(in) :: mol
         integer, optional, intent(in) :: nshell(:)
         type(error_type), allocatable, intent(out) :: error
      end subroutine charge_maker
   end interface

contains


!> Collect all exported unit tests
subroutine collect_coulomb_firstorder(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-atom-gxtb-m01", test_e_gxtb_m01_atom), &
      new_unittest("energy-shell-gxtb-m01", test_e_gxtb_m01_shell), &
      new_unittest("energy-atom-gxtb-m02", test_e_gxtb_m02_atom), &
      new_unittest("energy-shell-gxtb-m02", test_e_gxtb_m02_shell), &
      new_unittest("energy-shell-gxtb-oxacb-pbc", test_e_gxtb_oxacb), &
      new_unittest("energy-atom-gxtb-oxacb-sc", test_e_gxtb_oxacb_sc), &
      new_unittest("potential-atom-gxtb-m01", test_p_gxtb_m01_atom), &
      new_unittest("potential-shell-gxtb-m01", test_p_gxtb_m01_shell), &
      new_unittest("potential-atom-gxtb-m02", test_p_gxtb_m02_atom), &
      new_unittest("potential-shell-gxtb-m02", test_p_gxtb_m02_shell), &
      new_unittest("gradient-atom-gxtb-mb03", test_g_gxtb_m03_atom), &
      new_unittest("gradient-shell-gxtb-mb03", test_g_gxtb_m03_shell), &
      new_unittest("gradient-atom-gxtb-mb04", test_g_gxtb_m04_atom), &
      new_unittest("gradient-shell-gxtb-mb04", test_g_gxtb_m04_shell), &
      new_unittest("gradient-atom-gxtb-urea-pbc", test_g_gxtb_urea_atom), &
      new_unittest("sigma-atom-gxtb-mb05", test_s_gxtb_m05_atom), &
      new_unittest("sigma-shell-gxtb-mb05", test_s_gxtb_m05_shell), &
      new_unittest("sigma-atom-gxtb-mb06", test_s_gxtb_m06_atom), &
      new_unittest("sigma-shell-gxtb-mb06", test_s_gxtb_m06_shell), &
      new_unittest("sigma-atom-gxtb-pyrazine-pbc", test_s_gxtb_pyrazine_atom), &
      new_unittest("potential-gradient-atom-gxtb-effceh-lih", test_pg_gxtb_ceh_lih_atom), &
      new_unittest("potential-gradient-shell-gxtb-effceh-lih", test_pg_gxtb_ceh_lih_shell), &
      new_unittest("potential-gradient-atom-gxtb-effceh-m15", test_pg_gxtb_ceh_m15_atom), &
      new_unittest("potential-gradient-shell-gxtb-effceh-m15", test_pg_gxtb_ceh_m15_shell), &
      new_unittest("potential-gradient-atom-gxtb-effceh-m16", test_pg_gxtb_ceh_m16_atom), &
      new_unittest("potential-gradient-shell-gxtb-effceh-m16", test_pg_gxtb_ceh_m16_shell), &
      new_unittest("potential-sigma-atom-gxtb-effceh-lih", test_ps_gxtb_ceh_lih_atom), &
      new_unittest("potential-sigma-shell-gxtb-effceh-lih", test_ps_gxtb_ceh_lih_shell), &
      new_unittest("potential-sigma-atom-gxtb-effceh-co2", test_ps_gxtb_ceh_co2_atom), &
      new_unittest("potential-sigma-shell-gxtb-effceh-co2", test_ps_gxtb_ceh_co2_shell), &
      new_unittest("potential-sigma-atom-gxtb-effceh-m05", test_ps_gxtb_ceh_mb05_atom), &
      new_unittest("potential-sigma-shell-gxtb-effceh-m05", test_ps_gxtb_ceh_mb05_shell), &
      new_unittest("potential-sigma-atom-gxtb-effceh-m17", test_ps_gxtb_ceh_mb17_atom), &
      new_unittest("potential-sigma-shell-gxtb-effceh-m17", test_ps_gxtb_ceh_mb17_shell) &
      ]

end subroutine collect_coulomb_firstorder


!> Factory to setup the CEH basis set for testing of the potential gradient
subroutine make_basis(bas, mol, ng)
   type(basis_type), intent(out) :: bas
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: ng

   integer, parameter :: nsh(20) = [&
   & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]
   integer, parameter :: lsh(3, 20) = reshape([&
   & 0, 0, 0,  0, 0, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0, &
   & 0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 0,  0, 1, 2, &
   & 0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2, &
   & 0, 1, 0,  0, 1, 2], shape(lsh))

   integer, parameter :: pqn(3, 20) = reshape([&
   & 1, 0, 0,  1, 0, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0, &
   & 2, 2, 0,  2, 2, 0,  2, 2, 0,  2, 2, 0,  3, 3, 0,  3, 3, 3, &
   & 3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3,  3, 3, 3, &
   & 4, 4, 0,  4, 4, 3], shape(pqn))

   real(wp), parameter :: zeta(3, 20) = reshape([&
   & 1.23363166_wp, 0.00000000_wp, 0.00000000_wp, 2.27004605_wp, 0.00000000_wp, 0.00000000_wp, &
   & 0.86185456_wp, 1.42017184_wp, 0.00000000_wp, 1.76817995_wp, 1.44095844_wp, 0.00000000_wp, &
   & 2.06339837_wp, 1.52051807_wp, 0.00000000_wp, 2.56058582_wp, 1.86484737_wp, 0.00000000_wp, &
   & 2.71233631_wp, 2.19848968_wp, 0.00000000_wp, 3.21585650_wp, 2.41309737_wp, 0.00000000_wp, &
   & 3.82146807_wp, 2.63063636_wp, 0.00000000_wp, 4.62721228_wp, 2.53599954_wp, 0.00000000_wp, &
   & 0.93221172_wp, 1.55333839_wp, 0.00000000_wp, 1.77220557_wp, 1.59942632_wp, 2.98596647_wp, &
   & 2.26040231_wp, 1.78718151_wp, 2.00990188_wp, 1.85259089_wp, 1.81733349_wp, 1.65269988_wp, &
   & 2.65701241_wp, 2.03189759_wp, 2.03883661_wp, 2.60609998_wp, 2.16530440_wp, 2.41888232_wp, &
   & 2.78818934_wp, 2.24732894_wp, 1.99081182_wp, 2.55424399_wp, 2.20946190_wp, 1.93619550_wp, &
   & 1.73713827_wp, 1.33788617_wp, 0.00000000_wp, 2.47982574_wp, 1.07250770_wp, 2.11920764_wp],&
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

   call new_basis(bas, mol, nshell, cgto, accuracy=1.0_wp)

end subroutine make_basis


!> Factory to create onsite firstorder objects based on g-xTB values
subroutine make_coulomb_gxtb(coulomb, mol, shell, error)

   !> New coulomb object
   type(tb_coulomb), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Parameter: Number of shells selected from the q-vSZP basis set
   integer, parameter :: max_elem = 60
   integer, parameter :: p_nshell(max_elem) = [ &
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, & !1-20
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, & !21-40
      & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 4, 4, 4] !41-60

   !> Exponent, slope and offset of IP and EA split
   real(wp), parameter :: p_split_exp = 1.0000000000_wp
   real(wp), parameter :: p_split_slope = 0.0250000000_wp
   real(wp), parameter :: p_split_offset = 0.6666666660_wp

   !> Coordination number exponent
   real(wp), parameter :: p_cn_exp = 2.0680000000_wp

   !> IP/EA for first-order tight-binding
   real(wp), parameter :: p_ipea(4, 60) = reshape([&
      & -0.2236400997_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !1
      & -1.1075334263_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !2
      & -0.0492550777_wp, -0.2305159060_wp,  0.0000000000_wp,  0.0000000000_wp, & !3
      & -0.1009144987_wp, -0.2329069247_wp,  0.0000000000_wp,  0.0000000000_wp, & !4
      & -0.0818410518_wp, -0.3271309601_wp,  0.0000000000_wp,  0.0000000000_wp, & !5
      & -0.1778126924_wp, -0.4137050056_wp,  0.0000000000_wp,  0.0000000000_wp, & !6
      & -0.3463826358_wp, -0.4337759736_wp,  0.0000000000_wp,  0.0000000000_wp, & !7
      & -0.3459474621_wp, -0.4222766259_wp,  0.0000000000_wp,  0.0000000000_wp, & !8
      & -0.6834425857_wp, -0.4828894953_wp,  0.0000000000_wp,  0.0000000000_wp, & !9
      & -0.7163662160_wp, -0.4479512887_wp,  0.0000000000_wp,  0.0000000000_wp, & !10
      & -0.0508350279_wp, -0.2355591746_wp,  0.0000000000_wp,  0.0000000000_wp, & !11
      & -0.1773136481_wp, -0.2545024159_wp, -0.1501328541_wp,  0.0000000000_wp, & !12
      & -0.2283159759_wp, -0.3064793968_wp, -0.2560000000_wp,  0.0000000000_wp, & !13
      & -0.2515580615_wp, -0.2914315049_wp, -0.2755234416_wp,  0.0000000000_wp, & !14
      & -0.3005595368_wp, -0.2465951614_wp, -0.2189727363_wp,  0.0000000000_wp, & !15
      & -0.2457685577_wp, -0.2702870926_wp, -0.2357223088_wp,  0.0000000000_wp, & !16
      & -0.2817668461_wp, -0.2717163334_wp, -0.2828866286_wp,  0.0000000000_wp, & !17
      & -0.2325241630_wp, -0.3840452930_wp, -0.2928248171_wp,  0.0000000000_wp, & !18
      & -0.0060226065_wp, -0.1871655070_wp,  0.0000000000_wp,  0.0000000000_wp, & !19
      & -0.0770137434_wp, -0.1515379257_wp, -0.0478074694_wp,  0.0000000000_wp, & !20
      & -0.2469355319_wp, -0.1215351154_wp, -0.1770773067_wp,  0.0000000000_wp, & !21
      & -0.2028768917_wp, -0.0392612780_wp, -0.2215183861_wp,  0.0000000000_wp, & !22
      & -0.1629580112_wp, -0.1325469534_wp, -0.2159219247_wp,  0.0000000000_wp, & !23
      & -0.1394072312_wp, -0.2120951802_wp, -0.1896027217_wp,  0.0000000000_wp, & !24
      & -0.1754639178_wp, -0.1751567823_wp, -0.1938121594_wp,  0.0000000000_wp, & !25
      & -0.1547856717_wp, -0.1553160844_wp, -0.1618817164_wp,  0.0000000000_wp, & !26
      & -0.1990300049_wp, -0.1925579773_wp, -0.1850989193_wp,  0.0000000000_wp, & !27
      & -0.1211017411_wp, -0.3257904674_wp, -0.2096679474_wp,  0.0000000000_wp, & !28
      & -0.1402553079_wp, -0.2281424215_wp, -0.2501727542_wp,  0.0000000000_wp, & !29
      & -0.1205889582_wp, -0.0506286556_wp,  0.0000000000_wp,  0.0000000000_wp, & !30
      & -0.2147187913_wp, -0.2517410182_wp, -0.3044988154_wp,  0.0000000000_wp, & !31
      & -0.2910649696_wp, -0.2836820301_wp, -0.3300697248_wp,  0.0000000000_wp, & !32
      & -0.3219796985_wp, -0.2233293047_wp, -0.3358771314_wp,  0.0000000000_wp, & !33
      & -0.2861354178_wp, -0.2783561512_wp, -0.3169690889_wp,  0.0000000000_wp, & !34
      & -0.1205238420_wp, -0.2171649042_wp, -0.1856106723_wp,  0.0000000000_wp, & !35
      & -0.3248436816_wp, -0.3137685530_wp, -0.2038739084_wp,  0.0000000000_wp, & !36
      & -0.0179263928_wp, -0.1528218030_wp,  0.0000000000_wp,  0.0000000000_wp, & !37
      & -0.1156173327_wp, -0.1294757846_wp, -0.0860441641_wp,  0.0000000000_wp, & !38
      & -0.1292287046_wp, -0.0920945481_wp, -0.1784847378_wp,  0.0000000000_wp, & !39
      & -0.1752814706_wp, -0.0749896803_wp, -0.1684808948_wp,  0.0000000000_wp, & !40
      & -0.1691597656_wp, -0.0840127675_wp, -0.1791904532_wp,  0.0000000000_wp, & !41
      & -0.1722251014_wp, -0.0933448262_wp, -0.1499786907_wp,  0.0000000000_wp, & !42
      & -0.1041107268_wp, -0.1887955442_wp, -0.1797036857_wp,  0.0000000000_wp, & !43
      & -0.1437341003_wp, -0.1883735496_wp, -0.1841951068_wp,  0.0000000000_wp, & !44
      & -0.1057676951_wp, -0.1347699894_wp, -0.1606046242_wp,  0.0000000000_wp, & !45
      & -0.1715676919_wp, -0.2017276899_wp, -0.1525956319_wp,  0.0000000000_wp, & !46
      & -0.1638338009_wp, -0.2401222311_wp, -0.1032313774_wp,  0.0000000000_wp, & !47
      & -0.0988363496_wp, -0.0369765190_wp,  0.0000000000_wp,  0.0000000000_wp, & !48
      & -0.2180514628_wp, -0.1654159521_wp, -0.3457837668_wp,  0.0000000000_wp, & !49
      & -0.2412412699_wp, -0.2067782636_wp, -0.3084253381_wp,  0.0000000000_wp, & !50
      & -0.3207592010_wp, -0.2019467169_wp, -0.3253009951_wp,  0.0000000000_wp, & !51
      & -0.2744605924_wp, -0.1680824350_wp, -0.3783446430_wp,  0.0000000000_wp, & !52
      & -0.2377393210_wp, -0.2121912141_wp, -0.3253560916_wp,  0.0000000000_wp, & !53
      & -0.2987386026_wp, -0.2030374838_wp, -0.4253309567_wp,  0.0000000000_wp, & !54
      & -0.0677642843_wp, -0.2007465216_wp,  0.0000000000_wp,  0.0000000000_wp, & !55
      & -0.0509377496_wp, -0.0732018209_wp, -0.0540130497_wp,  0.0000000000_wp, & !56
      & -0.1218611287_wp, -0.1065790422_wp, -0.1081467941_wp,  0.0000000000_wp, & !57
      & -0.1710753191_wp, -0.1270853009_wp, -0.0624907619_wp, -0.3760017932_wp, & !58
      & -0.1644692595_wp, -0.1289050868_wp, -0.0427084369_wp, -0.4452902549_wp, & !59
      & -0.1340309951_wp, -0.1114222354_wp, -0.0593760960_wp, -0.4155633191_wp],& !60
      & shape(p_ipea))

   !> CN-dependence of the IP/EA for first-order tight-binding
   real(wp), parameter :: p_ipea_cn(60) = [&
      &  0.8067865070_wp,  0.2636069004_wp,  0.8646552118_wp,  0.5564736628_wp, & !1-4
      &  0.1467607257_wp,  0.1146711024_wp,  0.0251686527_wp,  0.1396723495_wp, & !5-8
      & -0.3378956948_wp,  0.5094660384_wp,  0.6453866485_wp,  0.5131213442_wp, & !9-12
      &  0.2905415760_wp,  0.0926453042_wp,  0.0974992421_wp,  0.0754281035_wp, & !13-16
      &  0.4499953879_wp, -0.2542185480_wp,  0.9622719304_wp,  0.3359066499_wp, & !17-20
      &  0.1215411848_wp,  0.1573069763_wp,  0.0596716713_wp,  0.4134271556_wp, & !21-24
      &  0.1455695987_wp,  0.0069530383_wp,  0.0948240914_wp,  0.0172331407_wp, & !25-28
      &  0.1384955329_wp,  1.0067843809_wp,  0.2732955043_wp,  0.1383288804_wp, & !29-32
      &  0.0218915887_wp,  0.0843631633_wp,  0.3191912129_wp,  0.0731610406_wp, & !33-36
      &  0.2958977960_wp,  0.3754678992_wp,  0.2452548742_wp,  0.2498401539_wp, & !37-40
      &  0.2348466385_wp,  0.3170989699_wp,  0.0567316849_wp,  0.0778681332_wp, & !41-44
      &  0.1000445263_wp,  0.1120457255_wp,  0.0241971449_wp,  1.2281256408_wp, & !45-48
      &  0.4013897951_wp,  0.2599035978_wp,  0.1156283118_wp,  0.3003489452_wp, & !49-52
      &  0.2690292677_wp,  0.0485166560_wp,  0.1988069756_wp,  0.1413607929_wp, & !53-56
      &  0.7858825497_wp,  0.0218956999_wp,  0.0225704851_wp,  0.0189202167_wp] !57-60
   
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

   real(wp), allocatable :: ipea(:, :), ipea_cn(:), tb_cn_rcov(:)
   integer, allocatable :: nshell(:)
   integer :: isp, izp, ish
   type(onsite_firstorder), allocatable :: tmp

   ! Setup coulomb interaction collection 
   tb_cn_rcov = p_cn_rcov(mol%num)
   allocate(coulomb)
   call new_coulomb(coulomb, mol, error, cn_count_type=cn_count%erf, &
      & cn_rcov=tb_cn_rcov, cn_exp=p_cn_exp)
   if (allocated(error)) return

   ! Setup onsite firstorder object
   allocate(tmp)
   ipea_cn = p_ipea_cn(mol%num)
   if (shell) then
      allocate(ipea(4, mol%nid))
      ipea(:, :) = 0.0_wp
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, 4
            ipea(ish, isp) = p_ipea(ish, izp) 
         end do
      end do
      nshell = p_nshell(mol%num)
   else
      allocate(ipea(1, mol%nid))
      ipea(:, :) = 0.0_wp
      do isp = 1, mol%nid
         izp = mol%num(isp)
         ipea(1, isp) = p_ipea(1, izp) 
      end do
   end if 
   call new_onsite_firstorder(tmp, mol, ipea, ipea_cn, p_split_exp, &
      & p_split_slope, p_split_offset, nshell)
   if (allocated(error)) return
   call move_alloc(tmp, coulomb%es1)

end subroutine make_coulomb_gxtb


!> Procedure to create CN based effective charges and gradients from CEH
subroutine get_charges_effceh(wfn, mol, nshell, error)

   !> New wavefunction object
   type(wavefunction_type), intent(inout) :: wfn

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return shell-resolved charges
   integer, optional, intent(in) :: nshell(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: ceh_cov_radii(20) = 0.5 * [&
   &  2.4040551903_wp,  1.8947380542_wp,  3.4227634078_wp,  3.5225408137_wp, &
   &  3.6150631704_wp,  2.8649682108_wp,  2.4695867541_wp,  2.3533691180_wp, &
   &  2.4992147462_wp,  3.3442607441_wp,  4.4665909451_wp,  4.3877250907_wp, &
   &  4.6647077385_wp,  4.2086223530_wp,  4.4750280107_wp,  4.2847281423_wp, &
   &  3.8560304959_wp,  3.9017061017_wp,  5.2392192639_wp,  5.1872031383_wp]

   real(wp), parameter :: pauling_en_ceh(20) = (1.0_wp/3.98_wp) * [ &
   &  1.9435211923_wp,  3.6116085622_wp,  2.4630915335_wp,  2.0658837656_wp, &
   &  2.3619778807_wp,  2.9484294262_wp,  3.8753937411_wp,  4.6235054741_wp, &
   &  3.9800000000_wp,  3.5124865506_wp,  2.3578254072_wp,  2.4225832022_wp, &
   &  2.1120078826_wp,  2.4607564741_wp,  2.7410779326_wp,  3.3517034720_wp, &
   &  4.1093492601_wp,  3.7979559518_wp,  2.4147937668_wp,  2.1974781961_wp]
   
   real(wp), allocatable :: lattr(:, :)
   real(wp), allocatable :: cn_en(:), dcn_endr(:, :, :), dcn_endL(:, :, :)
   class(ncoord_type), allocatable :: ncoord_en
   integer :: iat, ii, ish

   allocate(cn_en(mol%nat), dcn_endr(3, mol%nat, mol%nat), dcn_endL(3, 3, mol%nat))

   ! Get electronegativity-weighted coordination number
   call new_ncoord(ncoord_en, mol, cn_count%erf_en, error, &
      & rcov=ceh_cov_radii(mol%num), en=pauling_en_ceh(mol%num))
   if (allocated(error)) return
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call ncoord_en%get_coordination_number(mol, lattr, cn_en, dcn_endr, dcn_endL)

   ! Get effective charges and their gradients
   call get_effective_qat(mol, cn_en, wfn%qat, &
      & dcn_endr, dcn_endL, wfn%dqatdr, wfn%dqatdL)

   ! If shell-resolved charges are requested, partition atomic charges on shells
   if(present(nshell)) then
      ii = 0
      do iat = 1, mol%nat
         do ish = 1, nshell(iat)
            wfn%qsh(ii+ish, :) = wfn%qat(iat, :) / real(nshell(iat), wp)
            wfn%dqshdr(:, :, ii+ish, :) = wfn%dqatdr(:, :, iat, :) &
               & / real(nshell(iat), wp)
            wfn%dqshdL(:, :, ii+ish, :) = wfn%dqatdL(:, :, iat, :) &
               & / real(nshell(iat), wp)
         end do 
         ii = ii + nshell(iat)
      end do 
   end if

end subroutine get_charges_effceh


subroutine test_generic(error, mol, qat, qsh, make_coulomb, ref, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges for this structure
   real(wp), intent(in), optional :: qsh(:)

   !> Factory to create new onsite firstorder objects
   procedure(coulomb_maker) :: make_coulomb

   !> Reference value to check against
   real(wp), intent(in) :: ref

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   type(tb_coulomb), allocatable :: coulomb
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   real(wp) :: energy(mol%nat)
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   energy = 0.0_wp
   if (present(qsh)) then
      wfn%qsh = reshape(qsh, [size(qsh), 1])
   else
      wfn%qsh = reshape(qat, [size(qat), 1])
   end if
   wfn%qat = reshape(qat, [size(qat), 1])
   call make_coulomb(coulomb, mol, present(qsh), error)
   if (allocated(error)) return
   call coulomb%update(mol, cache)
   call coulomb%es1%get_energy(mol, cache, wfn, energy)

   call check(error, sum(energy), ref, thr=thr_)
   if (allocated(error)) then
      print*,ref, sum(energy), abs(sum(energy) - ref)
   end if

end subroutine test_generic


subroutine test_numpot(error, mol, qat, qsh, make_coulomb, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in), optional :: qsh(:)

   !> Factory to create new coulomb objects
   procedure(coulomb_maker) :: make_coulomb

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ish
   type(basis_type) :: bas
   type(tb_coulomb), allocatable :: coulomb
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot
   real(wp), allocatable :: numvat(:), numvsh(:)
   real(wp) :: er(mol%nat), el(mol%nat)
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   if (present(qsh)) then
      wfn%qsh = reshape(qsh, [size(qsh), 1])
      bas%nsh = size(qsh)
      bas%nao = size(qsh)
   else
      wfn%qsh = reshape(qat, [size(qat), 1])
      bas%nsh = size(qat)
      bas%nao = size(qat)
   end if
   wfn%qat = reshape(qat, [size(qat), 1])

   ! Setup potential with dummy basis set
   call new_potential(pot, mol, bas, 1, .true.)

   ! Setup coulomb object
   call make_coulomb(coulomb, mol, present(qsh), error)
   if (allocated(error)) return
   call coulomb%update(mol, cache)

   ! Numerical atomic potential
   allocate(numvat(mol%nat), source=0.0_wp)
   do iat = 1, mol%nat
      er = 0.0_wp
      el = 0.0_wp
      ! Right hand side
      wfn%qat(iat, 1) = wfn%qat(iat, 1) + step
      call coulomb%es1%get_energy(mol, cache, wfn, er)
      ! Left hand side
      wfn%qat(iat, 1) = wfn%qat(iat, 1) - 2*step
      call coulomb%es1%get_energy(mol, cache, wfn, el)

      wfn%qat(iat, 1) = wfn%qat(iat, 1) + step
      numvat(iat) = 0.5_wp*(sum(er) - sum(el))/step
   end do

   if (present(qsh)) then
      ! Numerical shell potential
      allocate(numvsh(bas%nsh), source=0.0_wp)   
      do ish = 1, bas%nsh
         er = 0.0_wp
         el = 0.0_wp
         ! Right hand side
         wfn%qsh(ish, 1) = wfn%qsh(ish, 1) + step
         call coulomb%es1%get_energy(mol, cache, wfn, er)
         ! Left hand side
         wfn%qsh(ish, 1) = wfn%qsh(ish, 1) - 2*step
         call coulomb%es1%get_energy(mol, cache, wfn, el)

         wfn%qsh(ish, 1) = wfn%qsh(ish, 1) + step
         numvsh(ish) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end if

   ! Analytic potentials
   call pot%reset()
   call coulomb%es1%get_potential(mol, cache, wfn, pot)

   if (any(abs(pot%vat(:, 1) - numvat) > thr_)) then
      call test_failed(error, "Atom-resolved potential does not match")
      write(*,*) "numerical potential:"
      print'(3es21.14)', numvat
      write(*,*) "analytical potential:"
      print'(3es21.14)', pot%vat(: ,1)
      write(*,*) "difference:"
      print'(3es21.14)', pot%vat(: ,1) - numvat
   end if

   if (present(qsh)) then
      if (any(abs(pot%vsh(:,1) - numvsh) > thr_)) then
         call test_failed(error, "Shell-resolved potential does not match")
         write(*,*) "numerical potential:"
         print'(3es21.14)', numvsh
         write(*,*) "analytical potential:"
         print'(3es21.14)', pot%vsh(: ,1)
         write(*,*) "difference:"
         print'(3es21.14)', pot%vsh(: ,1) - numvsh
      end if
   end if

end subroutine test_numpot


subroutine test_numgrad(error, mol, qat, qsh, make_coulomb, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges for this structure
   real(wp), intent(in), optional :: qsh(:)

   !> Factory to create new onsite firstorder objects
   procedure(coulomb_maker) :: make_coulomb

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic
   type(tb_coulomb), allocatable :: coulomb
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
   if (present(qsh)) then
      wfn%qsh = reshape(qsh, [size(qsh), 1])
   else
      wfn%qsh = reshape(qat, [size(qat), 1])
   end if
   wfn%qat = reshape(qat, [size(qat), 1])
   call make_coulomb(coulomb, mol, present(qsh), error)
   if (allocated(error)) return

   do iat = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call coulomb%update(mol, cache)
         call coulomb%es1%get_energy(mol, cache, wfn, er)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call coulomb%update(mol, cache)
         call coulomb%es1%get_energy(mol, cache, wfn, el)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call coulomb%update(mol, cache)
   call coulomb%es1%get_gradient(mol, cache, wfn, gradient, sigma)

   if (any(abs(gradient - numgrad) > thr_)) then
      call test_failed(error, "Gradient of energy does not match")
      write(*,*) "numerical gradient:"
      print'(3es21.14)', numgrad
      write(*,*) "analytical gradient:"
      print'(3es21.14)', gradient
      write(*,*) "difference:"
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_numgrad


subroutine test_numpotgrad(error, mol, get_charges, make_coulomb, shell, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Procedure to provide atomic charges (and their gradients)
   procedure(charge_maker) :: get_charges

   !> Factory to create new onsite firstorder objects
   procedure(coulomb_maker) :: make_coulomb

   !> Test shell-resolved
   logical, intent(in) :: shell

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic
   type(basis_type) :: bas
   type(potential_type) :: potl, potr
   type(tb_coulomb), allocatable :: coulomb
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   real(wp), allocatable :: numvatgrad(:, :, :), numvshgrad(:, :, :)
   integer, allocatable :: nshell(:)
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Setup potentials and wavefunction with dummy basis set 
   call make_basis(bas, mol, 6)
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, 1, 0.0_wp, .true.)
   call new_potential(potr, mol, bas, 1, .true.)
   call new_potential(potl, mol, bas, 1, .true.)
   call make_coulomb(coulomb, mol, shell, error)
   if (allocated(error)) return

   if(shell) then
      allocate(numvshgrad(3, mol%nat, bas%nsh), source=0.0_wp)
      allocate(nshell(mol%nat))
      nshell = bas%nsh_at
   end if
   allocate(numvatgrad(3, mol%nat, mol%nat), source=0.0_wp)

   do iat = 1, mol%nat
      do ic = 1, 3
         call potr%reset
         call potl%reset
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_charges(wfn, mol, nshell, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%es1%get_potential(mol, cache, wfn, potr)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_charges(wfn, mol, nshell, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%es1%get_potential(mol, cache, wfn, potl)
         
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         if(shell) then
            numvshgrad(ic, iat, :) = 0.5_wp*(potr%vsh(:,1) - potl%vsh(:,1))/step
         end if 
         numvatgrad(ic, iat, :) = 0.5_wp*(potr%vat(:,1) - potl%vat(:,1))/step
      end do
   end do

   call get_charges(wfn, mol, nshell, error)
   if (allocated(error)) return
   call coulomb%update(mol, cache)
   call coulomb%es1%get_potential(mol, cache, wfn, potl)
   call coulomb%es1%get_potential_gradient(mol, cache, wfn, potl)

   if(shell) then
      if (any(abs(potl%dvshdr(:,:,:,1) - numvshgrad) > thr_)) then
         call test_failed(error, "Gradient of shell-resolved potential does not match")
         write(*,*) "numerical potential gradient:"
         print'(3es21.14)', numvshgrad
         write(*,*) "analytical potential gradient:"
         print'(3es21.14)', potl%dvshdr(:,:,:,1)
         write(*,*) "difference:"
         print'(3es21.14)', potl%dvshdr(:,:,:,1) - numvshgrad
      end if
   end if
   if (any(abs(potl%dvatdr(:,:,:,1) - numvatgrad) > thr_)) then
      call test_failed(error, "Gradient of atom-resolved potential does not match")
      write(*,*) "numerical potential gradient:"
      print'(3es21.14)', numvatgrad
      write(*,*) "analytical potential gradient:"
      print'(3es21.14)', potl%dvatdr(:,:,:,1)
      write(*,*) "difference:"
      print'(3es21.14)', potl%dvatdr(:,:,:,1) - numvatgrad
   end if

end subroutine test_numpotgrad


subroutine test_numsigma(error, mol, qat, qsh, make_coulomb, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges for this structure
   real(wp), intent(in), optional :: qsh(:)

   !> Factory to create new electrostatic objects
   procedure(coulomb_maker) :: make_coulomb

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: ic, jc
   type(tb_coulomb), allocatable :: coulomb
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   real(wp) :: energy(mol%nat), er(mol%nat), el(mol%nat)
   real(wp) :: sigma(3, 3), eps(3, 3), numsigma(3, 3)
   real(wp), allocatable :: gradient(:, :), xyz(:, :), lattice(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate(gradient(3, mol%nat), xyz(3, mol%nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   if (present(qsh)) then
      wfn%qsh = reshape(qsh, [size(qsh), 1])
   else
      wfn%qsh = reshape(qat, [size(qat), 1])
   end if
   wfn%qat = reshape(qat, [size(qat), 1])
   call make_coulomb(coulomb, mol, present(qsh), error)
   if (allocated(error)) return

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
         call coulomb%update(mol, cache)
         call coulomb%es1%get_energy(mol, cache, wfn, er)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call coulomb%update(mol, cache)
         call coulomb%es1%get_energy(mol, cache, wfn, el)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (allocated(lattice)) mol%lattice = lattice
         numsigma(jc, ic) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call coulomb%update(mol, cache)
   call coulomb%es1%get_gradient(mol, cache, wfn, gradient, sigma)

   if (any(abs(sigma - numsigma) > thr_)) then
      call test_failed(error, "Strain derivatives do not match")
      write(*,*) "numerical gradient:"
      print'(3es21.14)', numsigma
      write(*,*) "analytical gradient:"
      print'(3es21.14)', sigma
      write(*,*) "difference:"
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_numsigma


subroutine test_numpotsigma(error, mol, get_charges, make_coulomb, shell, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Procedure to provide atomic charges (and their gradients)
   procedure(charge_maker) :: get_charges

   !> Factory to create new electrostatic objects
   procedure(coulomb_maker) :: make_coulomb

   !> Test shell-resolved
   logical, intent(in) :: shell

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: ic, jc
   type(basis_type) :: bas
   type(potential_type) :: potl, potr
   type(tb_coulomb), allocatable :: coulomb
   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   real(wp) :: eps(3, 3)
   real(wp), allocatable :: numvatsigma(:, :, :), numvshsigma(:, :, :)
   real(wp), allocatable :: xyz(:, :), lattice(:, :)
   integer, allocatable :: nshell(:)
   real(wp), parameter :: unity(3, 3) = reshape(&
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   allocate(xyz(3, mol%nat), source=0.0_wp)

   ! Setup potentials and wavefunction with dummy basis set 
   call make_basis(bas, mol, 6)
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, 1, 0.0_wp, .true.)
   call new_potential(potr, mol, bas, 1, .true.)
   call new_potential(potl, mol, bas, 1, .true.)
   call make_coulomb(coulomb, mol, shell, error)
   if (allocated(error)) return

   if(shell) then
      allocate(numvshsigma(3, 3, bas%nsh), source=0.0_wp)
      allocate(nshell(mol%nat))
      nshell = bas%nsh_at
   end if
   allocate(numvatsigma(3, 3, mol%nat), source=0.0_wp)

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   if (any(mol%periodic)) lattice = mol%lattice
   do ic = 1, 3
      do jc = 1, 3
         call potr%reset
         call potl%reset
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call get_charges(wfn, mol, nshell, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%es1%get_potential(mol, cache, wfn, potr)

         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call get_charges(wfn, mol, nshell, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%es1%get_potential(mol, cache, wfn, potl)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (allocated(lattice)) mol%lattice = lattice
         if(shell) then
            numvshsigma(jc, ic, :) = 0.5_wp*(potr%vsh(:,1) - potl%vsh(:,1))/step
         end if
         numvatsigma(jc, ic, :) = 0.5_wp*(potr%vat(:,1) - potl%vat(:,1))/step
      end do
   end do

   call get_charges(wfn, mol, nshell, error)
   if (allocated(error)) return
   call coulomb%update(mol, cache)
   call coulomb%es1%get_potential(mol, cache, wfn, potl)
   call coulomb%es1%get_potential_gradient(mol, cache, wfn, potl)

   if(shell) then
      if (any(abs(potl%dvshdL(:,:,:,1) - numvshsigma) > thr_)) then
         call test_failed(error, "Sigma of shell-resolved potential does not match")
         write(*,*) "numerical potential sigma:"
         print'(3es21.14)', numvshsigma
         write(*,*) "analytical potential sigma:"
         print'(3es21.14)', potl%dvshdL(:,:,:,1)
         write(*,*) "difference:"
         print'(3es21.14)', potl%dvshdL(:,:,:,1) - numvshsigma
      end if
   end if
   if (any(abs(potl%dvatdL(:,:,:,1) - numvatsigma) > thr_)) then
      call test_failed(error, "Sigma of atom-resolved potential does not match")
      write(*,*) "numerical potential sigma:"
      print'(3es21.14)', numvatsigma
      write(*,*) "analytical potential sigma:"
      print'(3es21.14)', potl%dvatdL(:,:,:,1)
      write(*,*) "difference:"
      print'(3es21.14)', potl%dvatdL(:,:,:,1) - numvatsigma
   end if

end subroutine test_numpotsigma


subroutine test_e_gxtb_m01_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  9.115222032717878E-1_wp, -7.820978364181586E-2_wp, -8.085043301444819E-1_wp, &
      & -1.370008894958090E-1_wp, -4.522670977383796E-1_wp,  2.366794200051806E-1_wp, &
      & -1.056480024872493E-1_wp, -7.542974658029267E-1_wp, -6.084417930861796E-1_wp, &
      &  3.252603520948062E-1_wp,  2.590209057852253E-1_wp, -3.245713596482476E-1_wp, &
      &  9.588528411034325E-2_wp,  7.104890489542338E-1_wp, -3.468331265256113E-1_wp, &
      &  1.076916633846319E+0_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, qat, qsh, make_coulomb_gxtb, 0.72838286138592401_wp)

end subroutine test_e_gxtb_m01_atom

subroutine test_e_gxtb_m01_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  9.115222032717878E-1_wp, -7.820978364181586E-2_wp, -8.085043301444819E-1_wp, &
      & -1.370008894958090E-1_wp, -4.522670977383796E-1_wp,  2.366794200051806E-1_wp, &
      & -1.056480024872493E-1_wp, -7.542974658029267E-1_wp, -6.084417930861796E-1_wp, &
      &  3.252603520948062E-1_wp,  2.590209057852253E-1_wp, -3.245713596482476E-1_wp, &
      &  9.588528411034325E-2_wp,  7.104890489542338E-1_wp, -3.468331265256113E-1_wp, &
      &  1.076916633846319E+0_wp]
   real(wp), parameter :: qsh(*) = [ &
      &  4.452001114394941E-1_wp,  4.663220918322937E-1_wp, -7.820978364181586E-2_wp, &
      & -1.925657752689702E-1_wp, -6.159385548755116E-1_wp, -1.370008894958090E-1_wp, &
      & -9.716756357792056E-2_wp, -3.550995341604590E-1_wp,  2.366794200051806E-1_wp, &
      & -1.056480024872493E-1_wp, -1.692417149620185E-1_wp, -5.850557508409082E-1_wp, &
      & -1.817505622195756E-1_wp, -4.266912308666040E-1_wp,  3.252603520948062E-1_wp, &
      &  2.590209057852253E-1_wp, -7.093896295382240E-2_wp, -2.932293779939634E-1_wp, &
      &  3.959698129953824E-2_wp,  1.473378889672090E-1_wp, -5.145260485686576E-2_wp, &
      &  3.544877620709267E-1_wp,  3.560012868833071E-1_wp, -6.786424039218564E-2_wp, &
      & -2.789688861334256E-1_wp,  9.044625085911038E-1_wp,  5.497811147662567E-1_wp, &
      & -3.773269895110409E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, qat, qsh, make_coulomb_gxtb, 0.599763979699109_wp)

end subroutine test_e_gxtb_m01_shell

subroutine test_e_gxtb_m02_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      & -1.097288775394156E-1_wp, -4.735023909155955E-1_wp, -1.707389583649019E-1_wp, &
      & -6.433491759460934E-1_wp,  8.965254778844174E-1_wp,  3.351267282573080E-1_wp, &
      & -6.657210518576395E-2_wp, -5.047599279669956E-2_wp,  5.146699702222350E-1_wp, &
      &  3.455270064324552E-1_wp,  8.365308656920378E-2_wp,  7.840093894994569E-1_wp, &
      & -6.236807887791791E-1_wp, -7.429581513991801E-2_wp, -6.465871716055771E-2_wp, &
      & -6.825088374693604E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, qat, qsh, make_coulomb_gxtb, 0.46027485527623380_wp)

end subroutine test_e_gxtb_m02_atom

subroutine test_e_gxtb_m02_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      & -1.097288775394156E-1_wp, -4.735023909155955E-1_wp, -1.707389583649019E-1_wp, &
      & -6.433491759460934E-1_wp,  8.965254778844174E-1_wp,  3.351267282573080E-1_wp, &
      & -6.657210518576395E-2_wp, -5.047599279669956E-2_wp,  5.146699702222350E-1_wp, &
      &  3.455270064324552E-1_wp,  8.365308656920378E-2_wp,  7.840093894994569E-1_wp, &
      & -6.236807887791791E-1_wp, -7.429581513991801E-2_wp, -6.465871716055771E-2_wp, &
      & -6.825088374693604E-1_wp]
   real(wp), parameter :: qsh(*) = [ &
      & -1.097288775394156E-1_wp, -1.649240032771995E-1_wp, -4.263046913300617E-1_wp, &
      &  1.177263036916656E-1_wp,  8.265314281056635E-2_wp, -2.533921011754683E-1_wp, &
      & -1.316164635096304E-1_wp, -5.117327124364630E-1_wp,  6.670371585258130E-1_wp, &
      &  2.689765878292192E-1_wp, -3.948826847061487E-2_wp,  3.351267282573080E-1_wp, &
      & -6.657210518576395E-2_wp, -5.047599279669956E-2_wp,  1.895960251239481E-1_wp, &
      &  2.395294429005155E-1_wp,  8.554450219777143E-2_wp,  3.455270064324552E-1_wp, &
      &  1.534569746334520E-1_wp, -6.980388806424820E-2_wp,  3.986364504419448E-1_wp, &
      &  3.853729390575120E-1_wp, -1.606712300918751E-1_wp, -4.630095586873040E-1_wp, &
      & -7.429581513991801E-2_wp, -6.465871716055771E-2_wp, -2.013953379594855E-1_wp, &
      & -6.189894955055433E-1_wp,  1.378759959956683E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, qat, qsh, make_coulomb_gxtb, 0.404135326140822_wp)

end subroutine test_e_gxtb_m02_shell

subroutine test_e_gxtb_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! Currently a regression test based on GFN2 charges
   real(wp), parameter :: qat(*) = [&
      & 0.43142759318744994_wp, 0.43150967285357533_wp, 0.43146165224958366_wp, &
      & 0.43147407117646708_wp, 0.90726023139894296_wp, 0.90728336242768870_wp, &
      & 0.90737069516079072_wp, 0.90712544694131525_wp,-0.77945631396614101_wp, &
      &-0.78097308159910250_wp,-0.78002089184600454_wp,-0.77924355344945684_wp, &
      &-0.55876993965122157_wp,-0.55823164158343652_wp,-0.55927229556741787_wp, &
      &-0.55894500773303069_wp]   
   real(wp), parameter :: qsh(*)=[&
      & 0.43142759318744994_wp, 0.43150967285357533_wp, 0.43146165224958366_wp, &
      & 0.43147407117646708_wp, 0.12631625958529813_wp, 0.78094397181364483_wp, &
      & 0.12627921822540056_wp, 0.78100414420228814_wp, 0.12624063356686954_wp, &
      & 0.78113006159392118_wp, 0.12626556398329924_wp, 0.78085988295801601_wp, &
      & 0.19373236130795646_wp,-0.97318867527409747_wp, 0.19324264343503161_wp, &
      &-0.97421572503413412_wp, 0.19355130642179486_wp,-0.97357219826779939_wp, &
      & 0.19382295575165664_wp,-0.97306650920111348_wp, 0.23748355725433878_wp, &
      &-0.79625349690556035_wp, 0.23769535856264201_wp,-0.79592700014607853_wp, &
      & 0.23731240854425084_wp,-0.79658470411166871_wp, 0.23739675443019159_wp, &
      &-0.79634176216322228_wp]

   call get_structure(mol, "X23", "oxacb")
   call test_generic(error, mol, qat, qsh, make_coulomb_gxtb, 0.52870369969517261_wp)

end subroutine test_e_gxtb_oxacb

subroutine test_e_gxtb_oxacb_sc(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! Currently a regression test based on GFN2 charges
   real(wp), parameter :: qat1(*) = [&
      & 3.41731844312030E-1_wp, 3.41716020106239E-1_wp, 3.41730526585671E-1_wp,&
      & 3.41714427217954E-1_wp, 3.80996046757999E-1_wp, 3.80989821246195E-1_wp,&
      & 3.81000747720282E-1_wp, 3.80990494183703E-1_wp,-3.70406587264474E-1_wp,&
      &-3.70407565207006E-1_wp,-3.70417590212352E-1_wp,-3.70399716470705E-1_wp,&
      &-3.52322260586075E-1_wp,-3.52304269439196E-1_wp,-3.52313440903261E-1_wp,&
      &-3.52298498047004E-1_wp]
   integer, parameter :: supercell(*) = [2, 2, 2]
   real(wp), parameter :: qat(*) = [spread(qat1, 2, product(supercell))]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "oxacb")
   call make_supercell(mol, supercell)
   call test_generic(error, mol, qat, qsh, make_coulomb_gxtb, &
      & 0.3826334951637162_wp*product(supercell), 1.0e-7_wp)

end subroutine test_e_gxtb_oxacb_sc

subroutine make_supercell(mol, rep)
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


subroutine test_p_gxtb_m01_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  9.115222032717878E-1_wp, -7.820978364181586E-2_wp, -8.085043301444819E-1_wp, &
      & -1.370008894958090E-1_wp, -4.522670977383796E-1_wp,  2.366794200051806E-1_wp, &
      & -1.056480024872493E-1_wp, -7.542974658029267E-1_wp, -6.084417930861796E-1_wp, &
      &  3.252603520948062E-1_wp,  2.590209057852253E-1_wp, -3.245713596482476E-1_wp, &
      &  9.588528411034325E-2_wp,  7.104890489542338E-1_wp, -3.468331265256113E-1_wp, &
      &  1.076916633846319E+0_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "01")
   call test_numpot(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_p_gxtb_m01_atom

subroutine test_p_gxtb_m01_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  9.115222032717878E-1_wp, -7.820978364181586E-2_wp, -8.085043301444819E-1_wp, &
      & -1.370008894958090E-1_wp, -4.522670977383796E-1_wp,  2.366794200051806E-1_wp, &
      & -1.056480024872493E-1_wp, -7.542974658029267E-1_wp, -6.084417930861796E-1_wp, &
      &  3.252603520948062E-1_wp,  2.590209057852253E-1_wp, -3.245713596482476E-1_wp, &
      &  9.588528411034325E-2_wp,  7.104890489542338E-1_wp, -3.468331265256113E-1_wp, &
      &  1.076916633846319E+0_wp]
   real(wp), parameter :: qsh(*) = [ &
      &  4.452001114394941E-1_wp,  4.663220918322937E-1_wp, -7.820978364181586E-2_wp, &
      & -1.925657752689702E-1_wp, -6.159385548755116E-1_wp, -1.370008894958090E-1_wp, &
      & -9.716756357792056E-2_wp, -3.550995341604590E-1_wp,  2.366794200051806E-1_wp, &
      & -1.056480024872493E-1_wp, -1.692417149620185E-1_wp, -5.850557508409082E-1_wp, &
      & -1.817505622195756E-1_wp, -4.266912308666040E-1_wp,  3.252603520948062E-1_wp, &
      &  2.590209057852253E-1_wp, -7.093896295382240E-2_wp, -2.932293779939634E-1_wp, &
      &  3.959698129953824E-2_wp,  1.473378889672090E-1_wp, -5.145260485686576E-2_wp, &
      &  3.544877620709267E-1_wp,  3.560012868833071E-1_wp, -6.786424039218564E-2_wp, &
      & -2.789688861334256E-1_wp,  9.044625085911038E-1_wp,  5.497811147662567E-1_wp, &
      & -3.773269895110409E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_numpot(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_p_gxtb_m01_shell

subroutine test_p_gxtb_m02_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      & -1.097288775394156E-1_wp, -4.735023909155955E-1_wp, -1.707389583649019E-1_wp, &
      & -6.433491759460934E-1_wp,  8.965254778844174E-1_wp,  3.351267282573080E-1_wp, &
      & -6.657210518576395E-2_wp, -5.047599279669956E-2_wp,  5.146699702222350E-1_wp, &
      &  3.455270064324552E-1_wp,  8.365308656920378E-2_wp,  7.840093894994569E-1_wp, &
      & -6.236807887791791E-1_wp, -7.429581513991801E-2_wp, -6.465871716055771E-2_wp, &
      & -6.825088374693604E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "02")
   call test_numpot(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_p_gxtb_m02_atom

subroutine test_p_gxtb_m02_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      & -1.097288775394156E-1_wp, -4.735023909155955E-1_wp, -1.707389583649019E-1_wp, &
      & -6.433491759460934E-1_wp,  8.965254778844174E-1_wp,  3.351267282573080E-1_wp, &
      & -6.657210518576395E-2_wp, -5.047599279669956E-2_wp,  5.146699702222350E-1_wp, &
      &  3.455270064324552E-1_wp,  8.365308656920378E-2_wp,  7.840093894994569E-1_wp, &
      & -6.236807887791791E-1_wp, -7.429581513991801E-2_wp, -6.465871716055771E-2_wp, &
      & -6.825088374693604E-1_wp]
   real(wp), parameter :: qsh(*) = [ &
      & -1.097288775394156E-1_wp, -1.649240032771995E-1_wp, -4.263046913300617E-1_wp, &
      &  1.177263036916656E-1_wp,  8.265314281056635E-2_wp, -2.533921011754683E-1_wp, &
      & -1.316164635096304E-1_wp, -5.117327124364630E-1_wp,  6.670371585258130E-1_wp, &
      &  2.689765878292192E-1_wp, -3.948826847061487E-2_wp,  3.351267282573080E-1_wp, &
      & -6.657210518576395E-2_wp, -5.047599279669956E-2_wp,  1.895960251239481E-1_wp, &
      &  2.395294429005155E-1_wp,  8.554450219777143E-2_wp,  3.455270064324552E-1_wp, &
      &  1.534569746334520E-1_wp, -6.980388806424820E-2_wp,  3.986364504419448E-1_wp, &
      &  3.853729390575120E-1_wp, -1.606712300918751E-1_wp, -4.630095586873040E-1_wp, &
      & -7.429581513991801E-2_wp, -6.465871716055771E-2_wp, -2.013953379594855E-1_wp, &
      & -6.189894955055433E-1_wp,  1.378759959956683E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_numpot(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_p_gxtb_m02_shell


subroutine test_g_gxtb_m03_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & -6.722325895097825E-2_wp, -9.343346477695711E-1_wp,  7.423184694799490E-2_wp, &
      &  6.003342781579319E-1_wp,  8.155559152718013E-1_wp,  1.119087497546355E-1_wp, &
      & -4.000789327736856E-1_wp,  8.016756865222463E-2_wp,  8.005299562801327E-2_wp, &
      &  1.006389816802343E-1_wp, -6.704006636400752E-1_wp, -4.326713538195121E-1_wp, &
      & -5.263477368917042E-1_wp,  2.692434056872427E-1_wp,  8.190366482982604E-1_wp, &
      &  7.988620327393559E-2_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "03")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_g_gxtb_m03_atom

subroutine test_g_gxtb_m03_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & -6.722325895097825E-2_wp, -9.343346477695711E-1_wp,  7.423184694799490E-2_wp, &
      &  6.003342781579319E-1_wp,  8.155559152718013E-1_wp,  1.119087497546355E-1_wp, &
      & -4.000789327736856E-1_wp,  8.016756865222463E-2_wp,  8.005299562801327E-2_wp, &
      &  1.006389816802343E-1_wp, -6.704006636400752E-1_wp, -4.326713538195121E-1_wp, &
      & -5.263477368917042E-1_wp,  2.692434056872427E-1_wp,  8.190366482982604E-1_wp, &
      &  7.988620327393559E-2_wp]
   real(wp), parameter :: qsh(*) = [&
      & -1.800091594906625E-1_wp,  1.127859005396843E-1_wp, -1.895026460120119E-1_wp, &
      & -7.448320017575591E-1_wp,  7.423184694799490E-2_wp,  3.531328646188461E-1_wp, &
      &  2.472014135390858E-1_wp,  7.256572899281329E-1_wp,  1.929197853611723E-1_wp, &
      & -1.030211600175039E-1_wp, -4.389454578059868E-1_wp,  4.294231434737540E-1_wp, &
      &  1.214310640868683E-1_wp, -2.459935249292498E-1_wp, -1.540854078444358E-1_wp, &
      &  8.016756865222463E-2_wp,  8.005299562801327E-2_wp,  1.006389816802343E-1_wp, &
      & -1.682560705402540E-1_wp, -5.021445930998212E-1_wp, -1.913448527993926E-1_wp, &
      & -3.445558639142448E-1_wp,  1.032293628941253E-1_wp, -2.292258440600861E-1_wp, &
      & -2.971218928316182E-1_wp,  2.692434056872427E-1_wp,  3.989086310954189E-1_wp, &
      &  4.201280172028415E-1_wp,  7.988620327393559E-2_wp]

   call get_structure(mol, "MB16-43", "03")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_g_gxtb_m03_shell

subroutine test_g_gxtb_m04_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & -6.417476509346298E-2_wp, -2.256661931683722E-1_wp, -3.458626714493684E-2_wp, &
      & -3.860600924041735E-1_wp,  7.238991472693194E-1_wp, -5.777084491650841E-2_wp, &
      &  1.196120411671105E-1_wp, -5.461117600441268E-3_wp, -1.794514383019030E-2_wp, &
      & -1.211514276086538E-1_wp, -3.493736689115672E-1_wp,  9.244684713156974E-1_wp, &
      &  1.394030055851056E-1_wp, -6.975412178494578E-1_wp,  1.149356570285270E-1_wp, &
      & -6.258758429064337E-2_wp]
   real(wp), allocatable :: qsh(:) 

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_g_gxtb_m04_atom

subroutine test_g_gxtb_m04_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & -6.417476509346298E-2_wp, -2.256661931683722E-1_wp, -3.458626714493684E-2_wp, &
      & -3.860600924041735E-1_wp,  7.238991472693194E-1_wp, -5.777084491650841E-2_wp, &
      &  1.196120411671105E-1_wp, -5.461117600441268E-3_wp, -1.794514383019030E-2_wp, &
      & -1.211514276086538E-1_wp, -3.493736689115672E-1_wp,  9.244684713156974E-1_wp, &
      &  1.394030055851056E-1_wp, -6.975412178494578E-1_wp,  1.149356570285270E-1_wp, &
      & -6.258758429064337E-2_wp]
   real(wp), parameter :: qsh(*) = [&
      & -6.417476509346298E-2_wp,  1.351102135612946E-1_wp, -3.607764067296668E-1_wp, &
      & -3.458626714493684E-2_wp, -9.911483231505502E-2_wp, -2.869452600891185E-1_wp, &
      &  2.880944438849421E-1_wp,  4.358047033843773E-1_wp, -5.777084491650841E-2_wp, &
      &  1.196120411671105E-1_wp, -3.939325366594362E-1_wp,  1.711280521507259E-1_wp, &
      &  2.173433669082690E-1_wp, -1.794514383019030E-2_wp, -1.211514276086538E-1_wp, &
      & -1.811065357961501E-1_wp, -1.682671331154171E-1_wp,  5.809438805678405E-1_wp, &
      &  5.754056566852046E-1_wp, -2.318810659373478E-1_wp,  1.140488675203044E-1_wp, &
      & -1.540659926577748E-1_wp,  1.794201307225760E-1_wp, -1.738439396866633E-1_wp, &
      & -5.236972781627944E-1_wp,  1.149356570285270E-1_wp, -1.328352578199055E-2_wp, &
      & -4.930405850865283E-2_wp]

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_g_gxtb_m04_shell

subroutine test_g_gxtb_urea_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   ! Currently a regression test based on GFN1 charges
   real(wp), parameter :: qat(*) = [&
      & 5.55723890858218E-1_wp, 5.55765354442035E-1_wp, 2.50200231242017E-1_wp,&
      & 2.50282053284422E-1_wp, 2.39786980460652E-1_wp, 2.39895142481200E-1_wp,&
      & 2.50103678240412E-1_wp, 2.50425041601730E-1_wp, 2.39464477136495E-1_wp,&
      & 2.40360053062669E-1_wp,-4.38369096728919E-1_wp,-4.38451412936599E-1_wp,&
      &-4.38310020776279E-1_wp,-4.38617373848238E-1_wp,-6.59141030224988E-1_wp,&
      &-6.59117968294813E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "urea")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_g_gxtb_urea_atom


subroutine test_s_gxtb_m05_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &  1.006071930481400E-1_wp,  2.932046169821376E-1_wp, -6.319239009961786E-2_wp, &
      &  9.947186172656641E-2_wp,  3.996733253551813E-1_wp, -2.369287698363657E-1_wp, &
      &  4.017615743381375E-2_wp, -2.159124633334978E-1_wp, -4.882761896486447E-1_wp, &
      &  2.694398495760397E-1_wp, -2.477049718491462E-1_wp,  5.113277969124146E-1_wp, &
      & -4.555157689694500E-2_wp, -4.957542239784729E-2_wp,  1.032767738492523E-1_wp, &
      & -4.700357913169708E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "05")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_s_gxtb_m05_atom

subroutine test_s_gxtb_m05_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &  1.006071930481400E-1_wp,  2.932046169821376E-1_wp, -6.319239009961786E-2_wp, &
      &  9.947186172656641E-2_wp,  3.996733253551813E-1_wp, -2.369287698363657E-1_wp, &
      &  4.017615743381375E-2_wp, -2.159124633334978E-1_wp, -4.882761896486447E-1_wp, &
      &  2.694398495760397E-1_wp, -2.477049718491462E-1_wp,  5.113277969124146E-1_wp, &
      & -4.555157689694500E-2_wp, -4.957542239784729E-2_wp,  1.032767738492523E-1_wp, &
      & -4.700357913169708E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      &  1.239103318667341E-1_wp, -2.330313881859403E-2_wp, -2.277289244926683E-1_wp, &
      &  3.286240428012324E-1_wp,  1.923094986735736E-1_wp, -6.319239009961786E-2_wp, &
      &  9.947186172656641E-2_wp,  1.327520370431289E-1_wp,  2.669212883120524E-1_wp, &
      & -1.498303762898072E-1_wp, -4.572492914243425E-1_wp,  3.701508978777839E-1_wp, &
      &  4.017615743381375E-2_wp, -7.538154807238073E-2_wp, -1.825334695562404E-1_wp, &
      &  4.200255429512334E-2_wp, -7.970923512887529E-2_wp, -4.085669545197694E-1_wp, &
      &  2.694398495760397E-1_wp, -1.935153183172251E-1_wp, -4.267132252332066E-1_wp, &
      &  3.725235717012856E-1_wp,  6.936572403651375E-2_wp,  3.506170360257168E-1_wp, &
      &  9.134503685018397E-2_wp, -4.555157689694500E-2_wp, -4.957542239784729E-2_wp, &
      & -1.804516816371773E-1_wp,  5.496093141504543E-2_wp,  2.287675240713841E-1_wp, &
      & -1.327734865281642E-1_wp, -3.372623047888066E-1_wp]

   call get_structure(mol, "MB16-43", "05")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_s_gxtb_m05_shell

subroutine test_s_gxtb_m06_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &  1.377405933021885E-1_wp, -5.813945155630489E-1_wp,  9.316960884803893E-2_wp, &
      & -5.091659553940089E-1_wp, -2.471905377110744E-1_wp, -2.368709440274053E-1_wp, &
      &  7.173397168306400E-1_wp,  9.055198123372066E-2_wp,  2.251559688740128E-1_wp, &
      &  5.028356479613176E-1_wp,  2.383212159292323E-1_wp,  1.147119296846322E-1_wp, &
      & -3.536741849510672E-2_wp, -6.213000983852623E-2_wp, -3.872470057299220E-1_wp, &
      & -6.046027630361439E-2_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "06")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_s_gxtb_m06_atom

subroutine test_s_gxtb_m06_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &  1.377405933021885E-1_wp, -5.813945155630489E-1_wp,  9.316960884803893E-2_wp, &
      & -5.091659553940089E-1_wp, -2.471905377110744E-1_wp, -2.368709440274053E-1_wp, &
      &  7.173397168306400E-1_wp,  9.055198123372066E-2_wp,  2.251559688740128E-1_wp, &
      &  5.028356479613176E-1_wp,  2.383212159292323E-1_wp,  1.147119296846322E-1_wp, &
      & -3.536741849510672E-2_wp, -6.213000983852623E-2_wp, -3.872470057299220E-1_wp, &
      & -6.046027630361439E-2_wp]
   real(wp), parameter :: qsh(*) = [&
      &  1.338902997227369E-1_wp,  3.850293579451547E-3_wp, -1.002115992895833E-1_wp, &
      & -4.811829162734655E-1_wp,  9.316960884803893E-2_wp, -2.344402290465235E-1_wp, &
      & -2.747257263474854E-1_wp,  1.227409477745669E-1_wp, -3.699314854856413E-1_wp, &
      & -2.368709440274053E-1_wp,  5.393875951602345E-1_wp,  2.782261577536143E-1_wp, &
      & -1.002740360832087E-1_wp,  9.055198123372066E-2_wp, -5.156077569100348E-2_wp, &
      &  2.767167445650163E-1_wp, -2.049048981039336E-2_wp,  2.441204310828009E-1_wp, &
      &  2.792057066889100E-1_wp,  2.383212159292323E-1_wp,  1.147119296846322E-1_wp, &
      & -3.536741849510672E-2_wp, -6.213000983852623E-2_wp, -1.679654378577344E-1_wp, &
      & -2.192815678721876E-1_wp, -6.046027630361439E-2_wp]

   call get_structure(mol, "MB16-43", "06")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_s_gxtb_m06_shell

subroutine test_s_gxtb_pyrazine_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   ! Currently a regression test based on GFN2 charges
   real(wp), parameter :: qat(*) = [&
      & 1.23705026512437E-1_wp, 1.22537765989959E-1_wp, 1.23932626831231E-1_wp,&
      & 1.22418959326822E-1_wp, 1.23788033684569E-1_wp, 1.23643389058068E-1_wp,&
      & 1.22389811551880E-1_wp, 1.22402837718155E-1_wp, 3.17174594622166E-2_wp,&
      & 3.15585817789390E-2_wp, 3.14290005061543E-2_wp, 3.10358506314526E-2_wp,&
      & 3.11084225356749E-2_wp, 3.11187325528499E-2_wp, 3.10707965296438E-2_wp,&
      & 3.13508824430484E-2_wp,-3.09107070033859E-1_wp,-3.09144123845085E-1_wp,&
      &-3.08473365847409E-1_wp,-3.08483617386745E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "pyrazine")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_gxtb, thr_in=thr1)

end subroutine test_s_gxtb_pyrazine_atom


subroutine test_pg_gxtb_ceh_lih_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .false., thr_in=thr1)

end subroutine test_pg_gxtb_ceh_lih_atom

subroutine test_pg_gxtb_ceh_lih_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .true., thr_in=thr1)

end subroutine test_pg_gxtb_ceh_lih_shell

subroutine test_pg_gxtb_ceh_m15_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "15")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .false., thr_in=thr1)

end subroutine test_pg_gxtb_ceh_m15_atom

subroutine test_pg_gxtb_ceh_m15_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "15")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .true., thr_in=thr1)

end subroutine test_pg_gxtb_ceh_m15_shell

subroutine test_pg_gxtb_ceh_m16_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "16")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .false., thr_in=thr1)

end subroutine test_pg_gxtb_ceh_m16_atom

subroutine test_pg_gxtb_ceh_m16_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "16")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .true., thr_in=thr1)

end subroutine test_pg_gxtb_ceh_m16_shell


subroutine test_ps_gxtb_ceh_lih_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .false., thr_in=thr1)

end subroutine test_ps_gxtb_ceh_lih_atom

subroutine test_ps_gxtb_ceh_lih_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .true., thr_in=thr1)

end subroutine test_ps_gxtb_ceh_lih_shell

subroutine test_ps_gxtb_ceh_co2_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "CO2")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .false., thr_in=thr1)

end subroutine test_ps_gxtb_ceh_co2_atom

subroutine test_ps_gxtb_ceh_co2_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "CO2")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .true., thr_in=thr1)

end subroutine test_ps_gxtb_ceh_co2_shell

subroutine test_ps_gxtb_ceh_mb05_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .false., thr_in=thr1*10)

end subroutine test_ps_gxtb_ceh_mb05_atom

subroutine test_ps_gxtb_ceh_mb05_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .true., thr_in=thr1*10)

end subroutine test_ps_gxtb_ceh_mb05_shell

subroutine test_ps_gxtb_ceh_mb17_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "17")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .false., thr_in=thr1)

end subroutine test_ps_gxtb_ceh_mb17_atom

subroutine test_ps_gxtb_ceh_mb17_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "17")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_gxtb, &
      & .true., thr_in=thr1)

end subroutine test_ps_gxtb_ceh_mb17_shell

end module test_coulomb_firstorder
