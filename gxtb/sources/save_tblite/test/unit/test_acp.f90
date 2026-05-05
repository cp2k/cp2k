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

module test_acp
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type
   use mstore, only : get_structure
   use multicharge, only : get_eeqbc_charges
   use tblite_acp_cache, only : acp_cache
   use tblite_acp_type, only : acp_type, new_acp, get_acp, get_acp_gradient
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_cache, only : basis_cache
   use tblite_basis_qvszp, only : qvszp_basis_type, qvszp_cgto_type, &
      & new_qvszp_cgto, new_qvszp_basis
   use tblite_basis_type, only : basis_type, cgto_container, cgto_type, &
      & new_cgto, get_cutoff
   use tblite_blas, only : gemv
   use tblite_cutoff, only : get_lattice_points
   use tblite_scf_iterator, only : get_electronic_energy
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   implicit none
   private

   public :: collect_acp

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine basis_maker(bas, mol, error, scale_h0_basis, ncoord_bas)
         import :: basis_type, structure_type, error_type, ncoord_type
         class(basis_type), allocatable, intent(out) :: bas
         type(structure_type), intent(in) :: mol
         type(error_type), allocatable, intent(out) :: error
         logical, intent(in), optional :: scale_h0_basis
         class(ncoord_type), intent(out), allocatable, optional :: ncoord_bas
      end subroutine basis_maker
   end interface

contains


!> Collect all exported unit tests
subroutine collect_acp(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("acp-gxtb-h2", test_acp_gxtb_h2), &
      new_unittest("acp-gxtb-lih", test_acp_gxtb_lih), &
      new_unittest("acp-gxtb-s2", test_acp_gxtb_s2), &
      new_unittest("acp-gxtb-sih4", test_acp_gxtb_sih4), &
      new_unittest("acp-gxtb-cecl3", test_acp_gxtb_cecl3), &
      new_unittest("acp-gxtb-ce2", test_acp_gxtb_ce2), &
      new_unittest("qeff-gradient-gxtb-h2", test_qeff_grad_acp_gxtb_h2), &
      new_unittest("qeff-gradient-gxtb-lih", test_qeff_grad_acp_gxtb_lih), &
      new_unittest("qeff-gradient-gxtb-no", test_qeff_grad_acp_gxtb_no), &
      new_unittest("qeff-gradient-gxtb-s2", test_qeff_grad_acp_gxtb_s2), &
      new_unittest("qeff-gradient-gxtb-sih4", test_qeff_grad_acp_gxtb_sih4), &
      new_unittest("qeff-gradient-gxtb-cecl3", test_qeff_grad_acp_gxtb_cecl3), &
      new_unittest("gradient-gxtb-h2", test_g_acp_gxtb_h2), &
      new_unittest("gradient-gxtb-lih", test_g_acp_gxtb_lih), &
      new_unittest("gradient-gxtb-no", test_g_acp_gxtb_no), &
      new_unittest("gradient-gxtb-s2", test_g_acp_gxtb_s2), &
      new_unittest("gradient-gxtb-sih4", test_g_acp_gxtb_sih4), &
      new_unittest("gradient-gxtb-cecl3", test_g_acp_gxtb_cecl3) &
      ]

end subroutine collect_acp


subroutine make_qvszp_basis(bas, mol, error, scale_h0_basis, ncoord_bas)
   !> Basis set information
   class(basis_type), allocatable, intent(out) :: bas
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Scale basis set for H0
   logical, intent(in), optional :: scale_h0_basis
   !> Optional coordination number used in the basis set
   class(ncoord_type), intent(out), allocatable, optional :: ncoord_bas

   !> Parameter: Number of shells selected from the q-vSZP basis set
   integer, parameter :: pa_nshell(60) = [ &
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, & !1-20
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, & !21-40
      & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 4, 4, 4] !41-60

   !> Parameter: Scaling factor for the exponents in the q-vSZP basis set for H0
   real(wp), parameter :: ps_h0_qvszp_exp_scal(4, 60) = reshape([&
      &  1.0620846674_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !1
      &  0.8163974124_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !2
      &  1.1014998730_wp,  1.0499004352_wp,  0.0000000000_wp,  0.0000000000_wp, & !3
      &  1.1532244901_wp,  1.2602742203_wp,  0.0000000000_wp,  0.0000000000_wp, & !4
      &  1.0988296196_wp,  1.3187829969_wp,  0.0000000000_wp,  0.0000000000_wp, & !5
      &  1.1130996849_wp,  1.3138642534_wp,  0.0000000000_wp,  0.0000000000_wp, & !6
      &  1.0859858922_wp,  1.3257240965_wp,  0.0000000000_wp,  0.0000000000_wp, & !7
      &  1.0279753385_wp,  1.2723840699_wp,  0.0000000000_wp,  0.0000000000_wp, & !8
      &  0.8013994888_wp,  1.1664223021_wp,  0.0000000000_wp,  0.0000000000_wp, & !9
      &  0.7911486642_wp,  1.1226241564_wp,  0.0000000000_wp,  0.0000000000_wp, & !10
      &  1.1684677989_wp,  1.0384234320_wp,  0.0000000000_wp,  0.0000000000_wp, & !11
      &  1.0835068321_wp,  1.1280262032_wp,  0.9996533790_wp,  0.0000000000_wp, & !12
      &  0.9291714103_wp,  1.1513410877_wp,  0.8950841734_wp,  0.0000000000_wp, & !13
      &  1.0910392740_wp,  1.2010270023_wp,  1.1678636759_wp,  0.0000000000_wp, & !14
      &  1.1181829599_wp,  1.2019678153_wp,  1.2481922351_wp,  0.0000000000_wp, & !15
      &  0.9725756488_wp,  1.1742057774_wp,  1.1879062205_wp,  0.0000000000_wp, & !16
      &  1.1120954193_wp,  1.1399372972_wp,  1.2782888382_wp,  0.0000000000_wp, & !17
      &  0.8070769651_wp,  0.9964584149_wp,  1.1825178739_wp,  0.0000000000_wp, & !18
      &  1.2600071523_wp,  1.1010384925_wp,  0.0000000000_wp,  0.0000000000_wp, & !19
      &  1.0527809748_wp,  1.0342129559_wp,  0.9926415358_wp,  0.0000000000_wp, & !20
      &  1.1135641399_wp,  1.0112089602_wp,  1.1005123829_wp,  0.0000000000_wp, & !21
      &  1.0824155476_wp,  0.9549459968_wp,  1.6557926027_wp,  0.0000000000_wp, & !22
      &  1.0497530411_wp,  1.0045935806_wp,  1.8219073034_wp,  0.0000000000_wp, & !23
      &  0.9891295760_wp,  1.0717234428_wp,  1.6544442197_wp,  0.0000000000_wp, & !24
      &  1.0408276438_wp,  0.9586434423_wp,  1.3659523184_wp,  0.0000000000_wp, & !25
      &  1.1042352123_wp,  0.9290804856_wp,  1.2545416485_wp,  0.0000000000_wp, & !26
      &  1.0808396387_wp,  1.0215146898_wp,  1.2048044970_wp,  0.0000000000_wp, & !27
      &  1.0679425789_wp,  0.9667082917_wp,  1.4565888504_wp,  0.0000000000_wp, & !28
      &  1.0930771910_wp,  0.8769888633_wp,  1.3906906336_wp,  0.0000000000_wp, & !29
      &  1.1322328919_wp,  1.1592893464_wp,  0.0000000000_wp,  0.0000000000_wp, & !30
      &  0.9508858180_wp,  1.1435821682_wp,  0.9536339733_wp,  0.0000000000_wp, & !31
      &  0.9933770024_wp,  1.1717176209_wp,  1.0343302681_wp,  0.0000000000_wp, & !32
      &  0.9801330011_wp,  1.1542395263_wp,  1.2585023868_wp,  0.0000000000_wp, & !33
      &  1.0906903682_wp,  1.1729545470_wp,  1.2734763745_wp,  0.0000000000_wp, & !34
      &  1.0608907285_wp,  1.1775870650_wp,  1.1543606866_wp,  0.0000000000_wp, & !35
      &  0.8980779646_wp,  1.1311445432_wp,  1.1605830960_wp,  0.0000000000_wp, & !36
      &  1.2039370131_wp,  1.1762952233_wp,  0.0000000000_wp,  0.0000000000_wp, & !37
      &  0.9134154055_wp,  1.0357797174_wp,  0.8745452968_wp,  0.0000000000_wp, & !38
      &  1.0497974799_wp,  0.9898502205_wp,  0.9998510478_wp,  0.0000000000_wp, & !39
      &  1.1364241186_wp,  0.9622662134_wp,  1.3614020634_wp,  0.0000000000_wp, & !40
      &  1.0734316008_wp,  1.0138236031_wp,  1.4874165102_wp,  0.0000000000_wp, & !41
      &  0.9704792044_wp,  1.0398619534_wp,  1.9177590688_wp,  0.0000000000_wp, & !42
      &  1.1434064252_wp,  1.0655883593_wp,  1.7283632903_wp,  0.0000000000_wp, & !43
      &  1.0487575935_wp,  1.0126758902_wp,  1.6265342450_wp,  0.0000000000_wp, & !44
      &  1.0692292092_wp,  1.0026796095_wp,  1.6333197781_wp,  0.0000000000_wp, & !45
      &  1.0729544626_wp,  0.9847336304_wp,  1.4519986839_wp,  0.0000000000_wp, & !46
      &  1.2103370858_wp,  1.0373177212_wp,  1.4601053538_wp,  0.0000000000_wp, & !47
      &  1.1670159301_wp,  1.1402442510_wp,  0.0000000000_wp,  0.0000000000_wp, & !48
      &  0.9792785527_wp,  1.1549317947_wp,  0.9856053019_wp,  0.0000000000_wp, & !49
      &  0.9265049593_wp,  1.1628552897_wp,  1.0327513899_wp,  0.0000000000_wp, & !50
      &  0.9091528979_wp,  1.1722450396_wp,  1.0816566766_wp,  0.0000000000_wp, & !51
      &  0.9811550136_wp,  1.1078703898_wp,  1.1597119979_wp,  0.0000000000_wp, & !52
      &  0.9594210334_wp,  1.1371770500_wp,  1.0581229671_wp,  0.0000000000_wp, & !53
      &  0.8868060549_wp,  1.1198864731_wp,  1.0824504600_wp,  0.0000000000_wp, & !54
      &  1.3746277416_wp,  1.0457727034_wp,  0.0000000000_wp,  0.0000000000_wp, & !55
      &  0.9345919097_wp,  1.0529418454_wp,  0.9195245358_wp,  0.0000000000_wp, & !56
      &  1.0195116315_wp,  1.0069184760_wp,  0.9391047717_wp,  0.0000000000_wp, & !57
      &  1.0753714450_wp,  1.1077438646_wp,  1.0260945802_wp,  1.6635363040_wp, & !58
      &  1.2799870696_wp,  1.1021925863_wp,  1.0228191807_wp,  1.5175499791_wp, & !59
      &  1.1544985695_wp,  0.9966709187_wp,  1.2326992368_wp,  1.4961183536_wp],& !60
      & shape(ps_h0_qvszp_exp_scal))

   !> Parameter: Scaling of the k0 parameter in q-vSZP for H0
   real(wp), parameter :: pa_h0_qvszp_k0_scal(60) = [&
      &  0.7401773802_wp,  0.8297893344_wp,  0.8688182361_wp,  0.4826574818_wp, & !1-4
      &  0.8269570159_wp,  0.8045707711_wp,  0.6172097287_wp,  0.7768616480_wp, & !5-8
      &  1.0680591133_wp,  0.8831106195_wp,  0.9648691211_wp,  0.7584091416_wp, & !9-12
      &  1.0064287409_wp,  0.8409060168_wp,  0.5996029735_wp,  0.6427980929_wp, & !13-16
      &  0.7602455915_wp,  0.0266180416_wp,  0.9962245388_wp,  1.1609352422_wp, & !17-20
      &  1.1904230396_wp,  1.3561400801_wp,  1.3698508205_wp,  1.6149252821_wp, & !21-24
      &  1.1217511879_wp,  1.1611074890_wp,  1.0876577477_wp,  0.1203999287_wp, & !25-28
      &  1.3297855006_wp,  1.0467788477_wp,  1.0744457428_wp,  0.7686842287_wp, & !29-32
      &  0.5874122119_wp,  0.6158231523_wp,  0.8063827609_wp,  0.6083377968_wp, & !33-36
      &  1.0722884393_wp,  1.0710584434_wp,  0.7013584281_wp,  0.9903729238_wp, & !37-40
      &  1.4373221426_wp,  0.9779781520_wp,  1.1018404188_wp,  1.1521476276_wp, & !41-44
      &  1.3960383727_wp,  0.6819168284_wp,  0.9365263998_wp,  1.1132166924_wp, & !45-48
      &  0.9786191810_wp,  0.5149746806_wp,  0.3222169183_wp,  0.6536595967_wp, & !49-52
      &  0.6223719136_wp,  0.9540774226_wp,  0.8639208774_wp,  0.9600595068_wp, & !53-56
      &  0.9105314363_wp,  0.8195362185_wp,  0.6860825124_wp,  0.1614716024_wp] !57-60

   !> Parameter: Scaling of the k2 parameter in q-vSZP for H0
   real(wp), parameter :: pa_h0_qvszp_k2_scal(60) = [&
      &  0.6129511778_wp,  0.7897783231_wp,  0.9829848653_wp,  0.8759370242_wp, & !1-4
      &  0.7840913450_wp,  0.7905355159_wp,  0.7540170345_wp,  0.5528879874_wp, & !5-8
      & -0.2295655984_wp,  0.2129544706_wp,  1.0340178348_wp,  1.0235343906_wp, & !9-12
      &  1.0486873645_wp,  0.8253693377_wp,  0.6227306404_wp,  1.1647945829_wp, & !13-16
      &  1.1664494026_wp,  0.3907378512_wp,  1.0086685177_wp,  0.9916989871_wp, & !17-20
      &  1.4957174827_wp,  1.1908451336_wp,  0.9668065885_wp,  0.9826826550_wp, & !21-24
      &  0.9643044000_wp,  0.9707348786_wp,  0.9986296428_wp,  0.1020914169_wp, & !25-28
      &  1.5416582135_wp,  0.9237324130_wp,  1.0058201440_wp,  1.0209414874_wp, & !29-32
      &  0.6154754704_wp,  0.5174049071_wp,  1.5308335537_wp,  1.9523046259_wp, & !33-36
      &  0.9967596540_wp,  0.9925354309_wp,  1.0425025333_wp,  1.0298781463_wp, & !37-40
      &  1.0408850481_wp,  1.0228199733_wp,  0.9744252388_wp,  0.9991072925_wp, & !41-44
      &  1.0026006169_wp,  0.0571473677_wp,  0.8459767948_wp,  0.9749847762_wp, & !45-48
      &  1.0105949467_wp,  1.0492305424_wp,  0.5382106310_wp,  0.9108840337_wp, & !49-52
      &  1.1700253101_wp,  0.9793239579_wp,  1.0828797286_wp,  1.0440865904_wp, & !53-56
      &  0.9797774075_wp,  1.6013493194_wp,  2.4175437074_wp,  1.4605346600_wp] !57-60

   !> Parameter: Scaling of the k3 parameter in q-vSZP for H0
   real(wp), parameter :: pa_h0_qvszp_k3_scal(60) = [&
      & -0.0019370870_wp,-13.6514723719_wp,  0.0219518919_wp,  0.9833983588_wp, & !1-4
      &  1.7839253028_wp, -0.0597381241_wp,  0.5279772860_wp, -0.7374267949_wp, & !5-8
      & -0.6999826070_wp, -0.0965206689_wp,  1.1750061305_wp, -0.4537139062_wp, & !9-12
      &  2.5018663890_wp,  3.4206602512_wp,  1.3476190261_wp,  1.6825849259_wp, & !13-16
      &  0.2011624543_wp, 20.6614159582_wp,  1.8183979662_wp,  0.2769432546_wp, & !17-20
      &  0.7862470210_wp,  0.5554056029_wp,  1.4317558829_wp,  1.6202867535_wp, & !21-24
      &  0.1499030579_wp,  1.2392829398_wp,  1.1882525718_wp,  0.0624229608_wp, & !25-28
      & -2.0083915226_wp,  1.1987516124_wp, -1.0023885159_wp, -7.1604532284_wp, & !29-32
      &  0.8974741441_wp, -6.4644482733_wp,  0.1165796763_wp,  1.7159909027_wp, & !33-36
      &  0.6674865463_wp,  1.0188309720_wp,  1.0411334603_wp,  0.2097893861_wp, & !37-40
      & -0.4914723425_wp, -0.2543783835_wp,  0.7191733407_wp,  0.4499075208_wp, & !41-44
      &  0.4209110509_wp,  0.5986627624_wp,  0.9490974423_wp,  0.0939965449_wp, & !45-48
      & 11.2969034989_wp, -0.0533916904_wp,  5.6600822608_wp, -0.1705453433_wp, & !49-52
      &  2.6442605504_wp,  1.4250005488_wp,  0.5015235440_wp,  0.5985555238_wp, & !53-56
      &  1.0229471751_wp,  0.0446861462_wp,  0.0172994827_wp, -0.0204035724_wp] !57-60

   integer :: isp, izp, ish, stat, nprim
   integer, allocatable :: nshell(:)
   type(cgto_container), allocatable :: cgto(:, :)
   type(cgto_container), allocatable :: cgto_h0(:, :)
   type(qvszp_cgto_type), allocatable :: cgto_qvszp
   type(qvszp_basis_type), allocatable :: qvszp_bas
   real(wp) :: alpha(12), k0, k2, k3

   nshell = pa_nshell(mol%num)
   allocate(cgto(maxval(nshell), mol%nid))
   if (scale_h0_basis) then
      allocate(cgto_h0(maxval(nshell), mol%nid))
   end if
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, nshell(isp)
         allocate(cgto_qvszp)
         call new_qvszp_cgto(cgto_qvszp, izp, ish, .true., error)
         if (allocated(error)) return
         call move_alloc(cgto_qvszp, cgto(ish, isp)%raw)

         if (scale_h0_basis) then
            select type(p_cgto => cgto(ish, isp)%raw)
            type is (qvszp_cgto_type)
               nprim = p_cgto%nprim
               alpha(1:nprim) = p_cgto%alpha(1:nprim) &
                  & * ps_h0_qvszp_exp_scal(ish, izp)
               k0 = p_cgto%k0 * pa_h0_qvszp_k0_scal(izp)
               k2 = p_cgto%k2 * pa_h0_qvszp_k2_scal(izp)
               k3 = p_cgto%k3 * pa_h0_qvszp_k3_scal(izp)
            end select

            allocate(cgto_qvszp)
            call new_qvszp_cgto(cgto_qvszp, izp, ish, .true., error, &
               & expos=alpha, k0=k0, k2=k2, k3=k3)
            if (allocated(error)) return
            call move_alloc(cgto_qvszp, cgto_h0(ish, isp)%raw)
         end if

      end do
   end do

   allocate(qvszp_bas)
   call new_qvszp_basis(qvszp_bas, mol, nshell, cgto, error, &
      & accuracy=0.1_wp, cgto_h0=cgto_h0)
   if (allocated(error)) return
   ! Save the coordination number object
   if (present(ncoord_bas)) then
      ncoord_bas = qvszp_bas%ncoord
   end if
   call move_alloc(qvszp_bas, bas)

end subroutine make_qvszp_basis


subroutine make_acp(acp, mol)
   !> Atomic correction potential
   type(acp_type), intent(out) :: acp
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Number of atomic correction potential projectors
   integer, parameter :: p_nacp(60) = [&
      & 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 3, 4, & !1-20
      & 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, & !21-40
      & 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4] !41-60

   !> Angular momentum of atomic correction potential projectors
   integer, parameter :: p_l_acp(4, 60) = reshape([&
      & 0, 1, 2, 1, 0, 1, 2, 0, 0, 1, 2, 1, 0, 1, 2, 1, 0, 1, 2, 1, & !1-5
      & 0, 1, 2, 1, 0, 1, 2, 1, 0, 1, 2, 1, 0, 1, 2, 1, 0, 1, 2, 1, & !6-10
      & 0, 1, 2, 0, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & !11-15
      & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 0, 0, 1, 2, 3, & !16-20
      & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & !21-25
      & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & !26-30
      & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & !31-35
      & 0, 1, 2, 3, 0, 1, 2, 0, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & !36-40
      & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & !41-45
      & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & !46-50
      & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 0, & !51-55
      & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],& !56-60
      & shape(p_l_acp))

   !> Atomic correction potential level for the projectors
   real(wp), parameter :: p_acp_level(4, 60) = reshape([&
      & -0.1454617634_wp, -0.0374942112_wp, -0.0022824328_wp, -0.0033848998_wp, & !1
      & -0.4915474215_wp,  0.0281599235_wp, -0.0413950073_wp,  0.0000000000_wp, & !2
      &  0.0066299538_wp, -0.0670874205_wp, -0.0621048472_wp, -0.0008812649_wp, & !3
      & -0.0380822059_wp, -0.1093247889_wp, -0.0364004361_wp, -0.0019587898_wp, & !4
      & -0.1748514687_wp, -0.1460924817_wp, -0.0685038109_wp, -0.0037953126_wp, & !5
      & -0.3035649222_wp, -0.2085407560_wp, -0.0669114590_wp, -0.0063286708_wp, & !6
      & -0.3686208196_wp, -0.1869055488_wp, -0.0666227821_wp, -0.0075923156_wp, & !7
      & -0.3461467666_wp, -0.2566836350_wp, -0.7143271768_wp, -0.0165945334_wp, & !8
      & -0.7195861086_wp, -0.3141288061_wp,  0.0510057262_wp, -0.0016637314_wp, & !9
      & -0.9152350003_wp, -0.2845040970_wp, -0.0854540987_wp,  0.0030921426_wp, & !10
      &  0.1381279292_wp, -0.0602380954_wp, -0.0427757019_wp,  0.0000000000_wp, & !11
      & -0.0584805400_wp, -0.0631712887_wp, -0.0199624707_wp, -0.0503233082_wp, & !12
      & -0.1095303670_wp, -0.0651295179_wp, -0.0747754369_wp, -0.1039499848_wp, & !13
      & -0.1711469006_wp, -0.0954958739_wp, -0.0688230927_wp, -0.1354100980_wp, & !14
      & -0.2515899493_wp, -0.1098397595_wp, -0.0644176084_wp, -0.0697593023_wp, & !15
      & -0.4731747274_wp, -0.1024829796_wp, -0.0556818155_wp, -0.2105296460_wp, & !16
      & -0.5662777833_wp, -0.1882739645_wp, -0.0243104252_wp, -0.0032020242_wp, & !17
      & -0.6089084679_wp, -0.3342216270_wp, -0.0794422538_wp,  0.0254845544_wp, & !18
      &  0.0525063856_wp, -0.0918776201_wp, -0.0408481736_wp,  0.0000000000_wp, & !19
      &  0.0182408034_wp, -0.0020410867_wp, -0.0318880520_wp, -0.0668297020_wp, & !20
      &  0.0006499494_wp, -0.0017582019_wp, -0.1000811685_wp, -0.0211237246_wp, & !21
      & -0.0010829003_wp, -0.0054720091_wp, -0.0679843405_wp, -0.0032623720_wp, & !22
      & -0.0095663666_wp, -0.0190493059_wp, -0.0526698222_wp, -0.0022652288_wp, & !23
      & -0.0283600781_wp, -0.0035873155_wp, -0.0654733191_wp, -0.0021735487_wp, & !24
      & -0.0163516474_wp, -0.0337519148_wp, -0.0710823355_wp, -0.0020097747_wp, & !25
      & -0.0435929692_wp, -0.0209327669_wp, -0.0592823897_wp, -0.0001133359_wp, & !26
      & -0.0246714569_wp, -0.0266370104_wp, -0.0498142191_wp, -0.0018233728_wp, & !27
      & -0.0614209792_wp, -0.0834964564_wp, -0.0119503484_wp, -0.0034797366_wp, & !28
      & -0.0278135238_wp, -0.1315325538_wp, -0.0192658941_wp, -0.0018300423_wp, & !29
      & -0.0712696364_wp, -0.0846554541_wp, -0.0054216047_wp, -0.0061194330_wp, & !30
      & -0.1621764885_wp, -0.0866206892_wp, -0.0387456814_wp, -0.0287028250_wp, & !31
      & -0.1981507481_wp, -0.0712423873_wp, -0.0337600302_wp, -0.1482828378_wp, & !32
      & -0.3825023416_wp, -0.0830428510_wp, -0.0257027454_wp, -0.0667680542_wp, & !33
      & -0.3871471949_wp, -0.1251530288_wp, -0.0638564076_wp, -0.0269950057_wp, & !34
      & -0.4996752870_wp, -0.1164755419_wp, -0.0146311832_wp, -0.1666611885_wp, & !35
      & -0.6100649261_wp, -0.2470998486_wp, -0.0276722307_wp, -0.0016789015_wp, & !36
      &  0.0512749784_wp, -0.0683804795_wp, -0.0100092749_wp,  0.0000000000_wp, & !37
      &  0.0088946771_wp,  0.0681974612_wp, -0.0801685872_wp, -0.0476643511_wp, & !38
      & -0.0039427093_wp,  0.0056510647_wp, -0.1448105549_wp, -0.0298700941_wp, & !39
      &  0.0029956163_wp,  0.0003767109_wp, -0.0989309033_wp, -0.0100348649_wp, & !40
      &  0.0048007639_wp, -0.0011953295_wp, -0.0447431041_wp, -0.0018389225_wp, & !41
      & -0.0089982746_wp, -0.0080274862_wp, -0.0738949595_wp, -0.0010869979_wp, & !42
      & -0.0146683536_wp, -0.0186094420_wp, -0.0584097779_wp, -0.0020297977_wp, & !43
      & -0.0634130687_wp, -0.0078267430_wp, -0.0859998379_wp, -0.0030063828_wp, & !44
      & -0.0652105641_wp, -0.0241538699_wp, -0.0610913453_wp, -0.0020226291_wp, & !45
      & -0.0482050491_wp, -0.0529220928_wp, -0.0008601133_wp, -0.0035558280_wp, & !46
      & -0.0587181875_wp, -0.0945264892_wp, -0.0139422609_wp, -0.0024292348_wp, & !47
      & -0.0582120988_wp, -0.0872329670_wp, -0.0072648601_wp, -0.0029963682_wp, & !48
      & -0.1509675332_wp, -0.0570533977_wp, -0.0259138331_wp, -0.0450679521_wp, & !49
      & -0.2779660146_wp, -0.0303245240_wp, -0.0263124111_wp, -0.1964445235_wp, & !50
      & -0.3723689410_wp, -0.0790543776_wp, -0.0224141256_wp, -0.1396938451_wp, & !51
      & -0.4916536536_wp, -0.0804833789_wp, -0.0134706492_wp, -0.1800422276_wp, & !52
      & -0.5404689875_wp, -0.1069349036_wp, -0.0181481100_wp, -0.2188501589_wp, & !53
      & -0.4947601781_wp, -0.1492769138_wp, -0.0592818264_wp, -0.1200662433_wp, & !54
      &  0.0784940227_wp, -0.0004567508_wp, -0.0259614907_wp,  0.0000000000_wp, & !55
      & -0.0150965820_wp,  0.0694314766_wp, -0.1175067262_wp, -0.0232222233_wp, & !56
      &  0.0072897139_wp,  0.0647156816_wp, -0.0878552331_wp, -0.0349598510_wp, & !57
      & -0.0063418143_wp, -0.0028679307_wp, -0.0220549840_wp, -0.0020365145_wp, & !58
      & -0.0127442947_wp, -0.0011426164_wp, -0.0057410045_wp, -0.0015739788_wp, & !59
      & -0.0027602069_wp, -0.0041233996_wp, -0.0399775790_wp, -0.0014031092_wp],& !60
      & shape(p_acp_level))

   !> Atomic correction potential exponents for the projectors
   real(wp), parameter :: p_acp_exp(4, 60) = reshape([&
      &  0.5327194644_wp,  1.0861418701_wp,  0.3820887161_wp,  0.2074179438_wp, & !1
      &  0.8370716133_wp,  0.3352413470_wp,  0.3964778826_wp,  0.0000000000_wp, & !2
      &  0.3193903557_wp,  0.3381668370_wp,  0.2979817420_wp,  0.2241574660_wp, & !3
      &  0.1271891396_wp,  0.4512929798_wp,  0.2793237280_wp,  0.0488136169_wp, & !4
      &  0.1896471980_wp,  0.6132174105_wp,  0.4127478446_wp,  0.1146556963_wp, & !5
      &  0.2401597446_wp,  0.7746090267_wp,  0.4210714202_wp,  0.1146636923_wp, & !6
      &  0.2315842367_wp,  0.8165844095_wp,  0.6726706507_wp,  0.1678891922_wp, & !7
      &  0.2817367898_wp,  0.8179152162_wp,  2.6590490707_wp,  0.1799675540_wp, & !8
      &  0.6546188993_wp,  0.8383919277_wp,  0.8166686414_wp,  0.1833270302_wp, & !9
      &  1.1631103182_wp,  1.1063225763_wp,  0.5025779870_wp,  0.2900000000_wp, & !10
      &  0.3646665344_wp,  0.2410474383_wp,  0.2496912493_wp,  0.0000000000_wp, & !11
      &  0.3741373622_wp,  0.2312569367_wp,  0.2357767694_wp,  0.2875002912_wp, & !12
      &  0.1922877331_wp,  0.2569065389_wp,  0.3701661639_wp,  0.5947548592_wp, & !13
      &  0.1799036201_wp,  0.3557830451_wp,  0.2150016432_wp,  0.6634019151_wp, & !14
      &  0.1528213277_wp,  0.3954994646_wp,  0.2089405784_wp,  0.5430061479_wp, & !15
      &  0.2621221715_wp,  0.3870942403_wp,  0.2342010464_wp,  0.7514873896_wp, & !16
      &  0.2581341402_wp,  0.3828550836_wp,  0.2218865476_wp,  0.8235068168_wp, & !17
      &  0.2944791115_wp,  0.4212668863_wp,  0.2970622780_wp,  0.2625745582_wp, & !18
      &  0.0684078811_wp,  0.2170723131_wp,  0.1839318536_wp,  0.0000000000_wp, & !19
      &  0.0810208465_wp,  0.2037374614_wp,  0.2056459567_wp,  0.2991591813_wp, & !20
      &  0.2518188761_wp,  0.2221327310_wp,  0.3369411449_wp,  0.1656994710_wp, & !21
      &  0.3117232621_wp,  0.1741021079_wp,  0.2313537763_wp,  0.2216918178_wp, & !22
      &  0.1291195934_wp,  0.1691630909_wp,  0.2229798289_wp,  0.2984028231_wp, & !23
      &  0.0814712602_wp,  0.1474973925_wp,  0.2288188694_wp,  0.2613209109_wp, & !24
      &  0.1467426136_wp,  0.1305039572_wp,  0.2640359310_wp,  0.2899512855_wp, & !25
      &  0.0971319913_wp,  0.1357994416_wp,  0.2699425875_wp,  0.6150588582_wp, & !26
      &  0.1208144429_wp,  0.1609477570_wp,  0.2907096531_wp,  0.4031105848_wp, & !27
      &  0.1064873863_wp,  0.2630301643_wp,  0.2494200189_wp,  0.3648719310_wp, & !28
      &  0.1903277970_wp,  0.2709475201_wp,  0.3618956301_wp,  0.4994027890_wp, & !29
      &  0.1980792209_wp,  0.2027933018_wp,  0.1154125767_wp,  0.3858792191_wp, & !30
      &  0.1925919823_wp,  0.3760134963_wp,  0.1722113021_wp,  0.4355463373_wp, & !31
      &  0.1815767281_wp,  0.3334121728_wp,  0.1817309994_wp,  0.5441676930_wp, & !32
      &  0.2329578726_wp,  0.3243967255_wp,  0.1558464345_wp,  0.3937752729_wp, & !33
      &  0.1619411034_wp,  0.3867544092_wp,  0.7423849005_wp,  0.2200491487_wp, & !34
      &  0.2471164021_wp,  0.3688857049_wp,  0.1382485819_wp,  0.5647161739_wp, & !35
      &  0.2992714487_wp,  0.3569005238_wp,  0.1948768999_wp,  0.1721807465_wp, & !36
      &  0.0888093055_wp,  0.1864975635_wp,  0.3584109076_wp,  0.0000000000_wp, & !37
      &  0.1830175467_wp,  0.1442595409_wp,  0.2535251813_wp,  0.1812859805_wp, & !38
      &  0.8373677850_wp,  0.0778760089_wp,  0.3767218896_wp,  0.1396933443_wp, & !39
      &  0.2648983273_wp,  0.1850312721_wp,  0.2594323729_wp,  0.1068143161_wp, & !40
      &  0.2483409752_wp,  0.3374175071_wp,  0.1745934149_wp,  0.2946680952_wp, & !41
      &  0.0741692624_wp,  0.1225427312_wp,  0.1853030583_wp,  0.2355684447_wp, & !42
      &  0.0699017304_wp,  0.0936212509_wp,  0.2109482185_wp,  0.2032300530_wp, & !43
      &  0.0696721794_wp,  0.0965552790_wp,  0.2573783232_wp,  0.2690368743_wp, & !44
      &  0.0897440303_wp,  0.1365913082_wp,  0.2793670924_wp,  0.3509398327_wp, & !45
      &  0.1201311409_wp,  0.1901069709_wp,  0.2022851634_wp,  0.2856415010_wp, & !46
      &  0.1538005691_wp,  0.2289732305_wp,  0.2152245222_wp,  0.3243444924_wp, & !47
      &  0.1528842982_wp,  0.1718749423_wp,  0.1152391798_wp,  0.2689796113_wp, & !48
      &  0.1936765244_wp,  0.3509404825_wp,  0.1386412468_wp,  0.3062571913_wp, & !49
      &  0.1854949576_wp,  0.2688700662_wp,  0.1265019386_wp,  0.4645385480_wp, & !50
      &  0.1931323954_wp,  0.2499227224_wp,  0.1001634557_wp,  0.3762632494_wp, & !51
      &  0.2229723716_wp,  0.2173427473_wp,  0.0921964149_wp,  0.3923077384_wp, & !52
      &  0.2426704625_wp,  0.2841747839_wp,  0.1121929658_wp,  0.4459375000_wp, & !53
      &  0.4091709168_wp,  0.2326517711_wp,  0.2382389263_wp,  0.3541404616_wp, & !54
      &  0.0828524658_wp,  0.0727942966_wp,  0.1817561678_wp,  0.0000000000_wp, & !55
      &  0.5947453098_wp,  0.1392333913_wp,  0.3262591022_wp,  0.1415989655_wp, & !56
      &  0.1378169685_wp,  0.1687266498_wp,  0.3045791774_wp,  0.1694070041_wp, & !57
      &  0.5690779020_wp,  0.5996451145_wp,  0.1901160496_wp,  0.5000000000_wp, & !58
      &  1.1366496791_wp,  0.4636739029_wp,  0.2730623477_wp,  0.4500000000_wp, & !59
      &  0.7790100475_wp,  0.6141190465_wp,  0.2715680143_wp,  0.3523819477_wp],& !60
      & shape(p_acp_exp))

   integer :: isp, izp, iproj, lproj, mproj
   type(cgto_container), allocatable :: cgtp(:, :)
   real(wp) :: alpha(1), coeff(1)
   real(wp), allocatable :: levels(:, :)
   integer, allocatable :: nproj(:)

   nproj = p_nacp(mol%num)
   mproj = maxval(nproj)
   coeff(1) = 1.0_wp
   allocate(cgtp(mproj, mol%nid))
   allocate(levels(mproj, mol%nid), source=0.0_wp)
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do iproj = 1, nproj(isp)
         lproj = p_l_acp(iproj, izp)
         alpha(1) = p_acp_exp(iproj, izp)
         levels(iproj, isp) = p_acp_level(iproj, izp)
         allocate(cgto_type :: cgtp(iproj, isp)%raw)
         call new_cgto(cgtp(iproj, isp)%raw, 1, lproj, alpha, coeff, .true.)
      end do
   end do

   call new_acp(acp, mol, nproj, cgtp, levels, accuracy=0.1_wp)

end subroutine make_acp


subroutine test_acp_gen(error, mol, make_basis, ref, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Factory to create new basis set objects
   procedure(basis_maker) :: make_basis

   !> Reference Hamiltonian matrix
   real(wp), intent(in) :: ref(:, :)

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   class(basis_type), allocatable :: bas
   type(basis_cache) :: bcache
   type(acp_type) :: acp
   type(acp_cache) :: acache
   type(wavefunction_type) :: wfn_aux
   type(adjacency_list) :: list
   real(wp), allocatable :: lattr(:, :)
   real(wp), allocatable :: hamiltonian(:, :)
   real(wp) :: cutoff, thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Setup basis set
   call make_basis(bas, mol, error=error, scale_h0_basis=.false.)
   if (allocated(error)) return
   call check(error, bas%nao, size(ref, 1))
   if (allocated(error)) return

   ! Obtain EEQBC charges for charge adaptation
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)

   ! Update basis set with charge/CN scaling
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   ! Create ACP object
   call make_acp(acp, mol)

   ! Update ACP cache
   call acp%update(mol, acache)

   ! Build ACP Hamiltonian contribution
   allocate(hamiltonian(bas%nao, bas%nao), source=0.0_wp)
   call get_acp(mol, lattr, list, bas, bcache, acp, acache, hamiltonian)

   ! where(abs(hamiltonian) < thr) hamiltonian = 0.0_wp
   ! print '(*(6x,"&", 3(es20.14e1, "_wp":, ","), "&", /))', hamiltonian

   if (any(abs(hamiltonian - ref) > thr_)) then
      call test_failed(error, "ACP does not match.")
      write(*,*) 'Reference:'
      print'(3es21.14)', ref
      write(*,*) 'ACP:'
      print'(3es21.14)', hamiltonian
      write(*,*) 'Difference:'
      print'(3es21.14)', hamiltonian-ref
   end if

end subroutine test_acp_gen


subroutine test_acp_qeff_numgrad(error, mol, density, make_basis, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Density matrix for this structure
   real(wp), intent(in) :: density(:, :, :)

   !> Factory to create new basis set objects
   procedure(basis_maker) :: make_basis

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, isp, ish, nspin
   class(basis_type), allocatable :: bas
   type(basis_cache) :: bcache
   type(acp_type) :: acp
   type(acp_cache) :: acache
   class(ncoord_type), allocatable :: ncoord_bas
   type(wavefunction_type) :: wfn_aux, wfn
   type(adjacency_list) :: list
   real(wp), parameter :: cn_cutoff = 30.0_wp
   real(wp), allocatable :: cnbas(:)
   real(wp), allocatable :: lattr(:, :)
   real(wp) :: thr_, cutoff
   real(wp), allocatable :: hamiltonian(:, :), dEdqbas(:), dEdcnbas(:)
   real(wp), allocatable :: er(:), el(:), num_dEdqbas(:), num_dEdcnbas(:)
   real(wp), allocatable :: gradient(:, :)
   real(wp) :: sigma(3, 3)
   
   real(wp), parameter :: step = 1.0e-5_wp

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in
   nspin = size(density, dim=3)

   ! Setup vasis set
   call make_basis(bas, mol, error=error, scale_h0_basis=.false., &
      & ncoord_bas=ncoord_bas)
   if (allocated(error)) return

   ! Obtain EEQBC charges for charge adaptation
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return

   ! Calculate coordination number
   if (allocated(ncoord_bas)) then
      allocate(cnbas(mol%nat), source=0.0_wp)
      call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
      call ncoord_bas%get_coordination_number(mol, lattr, cnbas)
   end if

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)

   ! Update basis set with charge/CN scaling
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   ! Setup wavefunction
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, nspin, 0.0_wp, .false.)
   wfn%density = density

   ! Create ACP object
   call make_acp(acp, mol)

   ! Update ACP cache
   call acp%update(mol, acache)

   ! Calculate numerical gradient of the ACP energy for the coordination number
   allocate(num_dEdcnbas(mol%nat), num_dEdqbas(mol%nat), er(bas%nao), el(bas%nao), &
      & hamiltonian(bas%nao, bas%nao), source=0.0_wp)
   do iat = 1, mol%nat
      isp = mol%id(iat)

      ! Right hand side
      cnbas(iat) = cnbas(iat) + step
      do ish = 1, bas%nsh_id(isp)
         associate(cgto => bas%cgto(ish, isp)%raw)
            ! Perturb effective charge by coordination number
            call cgto%get_qeff(wfn_aux%qat(iat, 1), cnbas(iat), bcache%cgto(ish, iat)%qeff)
            ! Update normalization
            call cgto%get_normalization(bcache%cgto(ish, iat), .false.)
         end associate
      end do
      hamiltonian = 0.0_wp
      er = 0.0_wp
      ! Build ACP Hamiltonian contribution
      call get_acp(mol, lattr, list, bas, bcache, acp, acache, hamiltonian)
      ! Calculate electronic energy
      call get_electronic_energy(hamiltonian, wfn%density, er)

      ! Left hand side
      cnbas(iat) = cnbas(iat) - 2*step
      do ish = 1, bas%nsh_id(isp)
         associate(cgto => bas%cgto(ish, isp)%raw)
            ! Perturb effective charge by coordination number
            call cgto%get_qeff(wfn_aux%qat(iat, 1), cnbas(iat), bcache%cgto(ish, iat)%qeff)
            ! Update normalization
            call cgto%get_normalization(bcache%cgto(ish, iat), .false.)
         end associate
      end do
      hamiltonian = 0.0_wp
      el = 0.0_wp
      ! Build ACP Hamiltonian contribution
      call get_acp(mol, lattr, list, bas, bcache, acp, acache, hamiltonian)
      ! Calculate electronic energy
      call get_electronic_energy(hamiltonian, wfn%density, el)

      cnbas(iat) = cnbas(iat) + step
      do ish = 1, bas%nsh_id(isp)
         associate(cgto => bas%cgto(ish, isp)%raw)
            ! Perturb effective charge by coordination number
            call cgto%get_qeff(wfn_aux%qat(iat, 1), cnbas(iat), bcache%cgto(ish, iat)%qeff)
            ! Update normalization
            call cgto%get_normalization(bcache%cgto(ish, iat), .false.)
         end associate
      end do
      num_dEdcnbas(iat) = 0.5_wp * (sum(er) - sum(el)) / step
   end do

   ! Calculate numerical gradient of the ACP energy for the atomic charges
   do iat = 1, mol%nat
      isp = mol%id(iat)

      ! Right hand side
      wfn_aux%qat(iat, 1) = wfn_aux%qat(iat, 1) + step
      call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)
      hamiltonian = 0.0_wp
      er = 0.0_wp
      ! Build ACP Hamiltonian contribution
      call get_acp(mol, lattr, list, bas, bcache, acp, acache, hamiltonian)
      ! Calculate electronic energy
      call get_electronic_energy(hamiltonian, wfn%density, er)

      ! Left hand side
      wfn_aux%qat(iat, 1) = wfn_aux%qat(iat, 1) - 2*step
      call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)
      hamiltonian = 0.0_wp
      el = 0.0_wp
      ! Build ACP Hamiltonian contribution
      call get_acp(mol, lattr, list, bas, bcache, acp, acache, hamiltonian)
      ! Calculate electronic energy
      call get_electronic_energy(hamiltonian, wfn%density, el)

      wfn_aux%qat(iat, 1) = wfn_aux%qat(iat, 1) + step
      call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)
      num_dEdqbas(iat) = 0.5_wp * (sum(er) - sum(el)) / step
   end do

   ! Update charges with derivatives
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .true.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1), wfn_aux%dqatdr(:, :, :, 1), &
     & wfn_aux%dqatdL(:, :, :, 1))
   if (allocated(error)) return
   call bas%update(mol, bcache, .true., wfn_aux=wfn_aux)

   ! Prepare ACP with updated pv_overlap
   call acp%update(mol, acache)
   call get_acp(mol, lattr, list, bas, bcache, acp, acache, & 
      & hamiltonian)

   ! Analytic effective charge gradient of the ACP energy
   allocate(dEdcnbas(mol%nat), dEdqbas(mol%nat), gradient(3, mol%nat), source=0.0_wp)
   sigma = 0.0_wp
   call get_acp_gradient(mol, lattr, list, bas, bcache, acp, acache, & 
      & wfn, dEdcnbas, dEdqbas, gradient, sigma)

   if (any(abs(dEdcnbas - num_dEdcnbas) > thr_)) then
      call test_failed(error, "Gradient of ACP energy w.r.t. coordination number does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', dEdcnbas
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dEdcnbas
      write(*,*) 'Difference:'
      print'(3es21.14)', dEdcnbas-num_dEdcnbas
   end if

   if (any(abs(dEdqbas - num_dEdqbas) > thr_)) then
      call test_failed(error, "Gradient of ACP energy w.r.t. atomic charge does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', dEdqbas
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dEdqbas
      write(*,*) 'Difference:'
      print'(3es21.14)', dEdqbas-num_dEdqbas
   end if

end subroutine test_acp_qeff_numgrad


subroutine test_acp_numgrad(error, mol, density, make_basis, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Density matrix for this structure
   real(wp), intent(in) :: density(:, :, :)

   !> Factory to create new basis set objects
   procedure(basis_maker) :: make_basis

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic, ish, nspin
   class(basis_type), allocatable :: bas
   type(basis_cache) :: bcache
   type(acp_type) :: acp
   type(acp_cache) :: acache
   type(wavefunction_type) :: wfn_aux, wfn
   type(adjacency_list) :: list
   real(wp), allocatable :: lattr(:, :)
   real(wp) :: thr_, cutoff, sigma(3, 3)
   real(wp), allocatable :: hamiltonian(:, :), er(:), el(:), numgrad(:, :)
   real(wp), allocatable :: dEdq(:), dEdcnbas(:), dEdqbas(:), gradient(:, :)
   real(wp), parameter :: step = 1.0e-5_wp

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in
   nspin = size(density, dim=3)

   ! Setup basis set
   call make_basis(bas, mol, error=error, scale_h0_basis=.false.)
   if (allocated(error)) return

   ! Obtain EEQBC charges for charge adaptation
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)

   ! Update basis set with charge/CN scaling
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   ! Setup wavefunction
   call new_wavefunction(wfn, mol%nat, bas%nsh, bas%nao, nspin, &
      & 0.0_wp, .false.)
   wfn%density(:, :, :) = density

   ! Create ACP object
   call make_acp(acp, mol)

   ! Update ACP cache
   call acp%update(mol, acache)

   ! Numerical gradient of the ACP energy
   allocate(numgrad(3, mol%nat), hamiltonian(bas%nao, bas%nao), er(bas%nao), &
      & el(bas%nao), source=0.0_wp)
   do iat = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp

         ! Right hand side
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         hamiltonian = 0.0_wp
         ! Update the basis set
         call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
         if (allocated(error)) return
         call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)
         ! Update ACP
         call acp%update(mol, acache)
         ! Build ACP Hamiltonian contribution
         call get_acp(mol, lattr, list, bas, bcache, acp, acache, & 
            & hamiltonian)
         ! Obtain the electronic energy contribution
         call get_electronic_energy(hamiltonian, wfn%density, er)

         ! Left hand side
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         hamiltonian = 0.0_wp
         ! Update the basis set
         call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
         if (allocated(error)) return
         call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)
         ! Update ACP
         call acp%update(mol, acache)
         ! Build ACP Hamiltonian contribution
         call get_acp(mol, lattr, list, bas, bcache, acp, acache, & 
            & hamiltonian)
         ! Obtain the electronic energy contribution
         call get_electronic_energy(hamiltonian, wfn%density, el)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   ! Update charges with derivatives
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .true.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1), wfn_aux%dqatdr(:, :, :, 1), &
     & wfn_aux%dqatdL(:, :, :, 1))
   if (allocated(error)) return
   call bas%update(mol, bcache, .true., wfn_aux=wfn_aux)

   ! Prepare ACP with updated pv_overlap
   call acp%update(mol, acache)
   call get_acp(mol, lattr, list, bas, bcache, acp, acache, & 
      & hamiltonian)

   ! Analytic gradient of the ACP energy 
   allocate(dEdq(mol%nat), dEdcnbas(mol%nat), dEdqbas(mol%nat), gradient(3, mol%nat), &
      & source=0.0_wp)
   sigma(:, :) = 0.0_wp
   call get_acp_gradient(mol, lattr, list, bas, bcache, acp, acache, & 
      & wfn, dEdcnbas, dEdqbas, gradient, sigma)

   ! Add basis set gradient due to the charge and CN dependence
   if (bas%charge_dependent) then 
      call bas%get_basis_gradient(mol, dEdcnbas, dEdqbas, &
         & dEdq, gradient, sigma)
   end if

   ! Add gradient contribution due to external charge model derivatives
   if (allocated(wfn_aux%dqatdr)) then
      call gemv(wfn_aux%dqatdr(:, :, :, 1), dEdq, gradient, beta=1.0_wp)
   end if
   if (allocated(wfn_aux%dqatdL)) then
      call gemv(wfn_aux%dqatdL(:, :, :, 1), dEdq, sigma, beta=1.0_wp)
   end if

   if (any(abs(gradient - numgrad) > thr_)) then
      call test_failed(error, "Gradient of ACP does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', gradient
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', numgrad
      write(*,*) 'Difference:'
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_acp_numgrad


subroutine test_acp_gxtb_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 2
   real(wp), parameter :: acp(nao, nao) = reshape([&
      & -1.93539607281311E-01_wp, -1.70281555812325E-01_wp, -1.70281555812325E-01_wp, &
      & -1.93539607281311E-01_wp], shape(acp))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_acp_gen(error, mol, make_qvszp_basis, acp, thr_in=thr1)

end subroutine test_acp_gxtb_h2

subroutine test_acp_gxtb_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 5
   real(wp), parameter :: acp(nao, nao) = reshape([&
      & -1.04894205003368E-02_wp,  0.00000000000000E+00_wp, -1.79213752727101E-02_wp, &
      &  0.00000000000000E+00_wp, -3.75783994396044E-02_wp,  0.00000000000000E+00_wp, &
      & -3.83556246681492E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.79213752727101E-02_wp,  0.00000000000000E+00_wp, &
      & -6.52987506789424E-02_wp,  0.00000000000000E+00_wp, -7.68689393529899E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.83556246681492E-02_wp,  0.00000000000000E+00_wp, -3.75783994396044E-02_wp, &
      &  0.00000000000000E+00_wp, -7.68689393529899E-02_wp,  0.00000000000000E+00_wp, &
      & -1.44506610375521E-01_wp],shape(acp))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_acp_gen(error, mol, make_qvszp_basis, acp, thr_in=thr1)

end subroutine test_acp_gxtb_lih

subroutine test_acp_gxtb_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 18
   real(wp), parameter :: acp(nao, nao) = reshape([&
      & -4.89052651366402E-01_wp,  0.00000000000000E+00_wp, -4.96132069235079E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.46073969377099E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.95851971748565E-01_wp,  0.00000000000000E+00_wp,  1.88079311363139E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.27778281695495E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.04991066714579E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.22463255035933E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.59105643213628E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  3.43959063015462E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -4.96132069235079E-02_wp,  0.00000000000000E+00_wp, -1.70848987419594E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.34582686157549E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.88079311363139E-01_wp,  0.00000000000000E+00_wp,  6.55510460029589E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -4.06630616840534E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.04991066714579E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.22463255035933E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.59105643213628E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  3.43959063015462E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.37079185446647E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.27821829490348E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.22463255035933E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.85928429915537E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -3.43959063015462E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  3.49872422969472E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.46073969377099E-02_wp,  0.00000000000000E+00_wp, -5.34582686157549E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -8.10814261663483E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.27778281695495E-01_wp,  0.00000000000000E+00_wp,  4.06630616840534E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -8.94557944320327E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.22463255035933E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.85928429915537E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.43959063015462E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  3.49872422969472E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -4.37079185446646E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.27821829490348E-02_wp, &
      & -1.95851971748565E-01_wp,  0.00000000000000E+00_wp, -1.88079311363139E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.27778281695495E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -4.89052651366402E-01_wp,  0.00000000000000E+00_wp,  4.96132069235079E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.46073969377099E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.59105643213628E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -3.43959063015462E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.04991066714579E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.22463255035933E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.88079311363139E-01_wp,  0.00000000000000E+00_wp,  6.55510460029589E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  4.06630616840534E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  4.96132069235079E-02_wp,  0.00000000000000E+00_wp, -1.70848987419594E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  5.34582686157549E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.59105643213628E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -3.43959063015462E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.04991066714579E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.22463255035933E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.27821829490348E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.37079185446647E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  3.43959063015462E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  3.49872422969472E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.22463255035933E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.85928429915537E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.27778281695495E-01_wp,  0.00000000000000E+00_wp, -4.06630616840534E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -8.94557944320327E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.46073969377099E-02_wp,  0.00000000000000E+00_wp,  5.34582686157549E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -8.10814261663483E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  3.43959063015462E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  3.49872422969472E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.22463255035933E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.85928429915537E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.27821829490348E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -4.37079185446646E-02_wp],&
      & shape(acp))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "S2")
   call test_acp_gen(error, mol, make_qvszp_basis, acp, thr_in=thr1)

end subroutine test_acp_gxtb_s2

subroutine test_acp_gxtb_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 13
   real(wp), parameter :: acp(nao, nao) = reshape([&
      & -2.19558716718933E-01_wp, -6.12533616268649E-20_wp, -1.81627593093385E-18_wp, &
      &  3.53070031358048E-18_wp,  1.82815105314704E-18_wp,  3.35473652523001E-18_wp, &
      &  3.00369861789520E-20_wp,  7.40528452554017E-20_wp,  5.59128352808655E-21_wp, &
      & -1.12822566489047E-01_wp, -1.12822566489047E-01_wp, -1.12822566489047E-01_wp, &
      & -1.12822566489047E-01_wp, -6.12533616268649E-20_wp, -1.26435936037915E-01_wp, &
      & -2.72108237756479E-18_wp, -2.50836265373244E-18_wp,  3.36729226374660E-18_wp, &
      & -1.64278976529652E-18_wp, -5.13171156572083E-20_wp,  3.86788348521085E-02_wp, &
      &  1.60723274315891E-19_wp, -5.56863254485502E-02_wp,  5.56863254485502E-02_wp, &
      &  5.56863254485502E-02_wp, -5.56863254485502E-02_wp, -1.81627593093385E-18_wp, &
      & -2.72108237756479E-18_wp, -1.26435936037915E-01_wp,  9.59253847275846E-19_wp, &
      &  3.86788348521084E-02_wp,  9.29722795561503E-20_wp, -2.14053789319158E-19_wp, &
      &  3.32226456377319E-18_wp,  2.50703088147650E-20_wp,  5.56863254485502E-02_wp, &
      &  5.56863254485502E-02_wp, -5.56863254485502E-02_wp, -5.56863254485502E-02_wp, &
      &  3.53070031358048E-18_wp, -2.50836265373244E-18_wp,  9.59253847275846E-19_wp, &
      & -1.26435936037915E-01_wp,  1.02154688207015E-19_wp,  3.86788348521085E-02_wp, &
      &  5.13171156572083E-20_wp, -1.68344734676473E-18_wp,  2.73085288618749E-19_wp, &
      & -5.56863254485502E-02_wp,  5.56863254485502E-02_wp, -5.56863254485502E-02_wp, &
      &  5.56863254485502E-02_wp,  1.82815105314704E-18_wp,  3.36729226374660E-18_wp, &
      &  3.86788348521084E-02_wp,  1.02154688207015E-19_wp, -9.64335932602174E-02_wp, &
      &  2.58089550037090E-19_wp, -5.20084261335142E-19_wp,  6.70790940818228E-18_wp, &
      &  1.06590766475164E-20_wp, -4.90158691040716E-02_wp, -4.90158691040716E-02_wp, &
      &  4.90158691040716E-02_wp,  4.90158691040716E-02_wp,  3.35473652523001E-18_wp, &
      & -1.64278976529652E-18_wp,  9.29722795561503E-20_wp,  3.86788348521085E-02_wp, &
      &  2.58089550037090E-19_wp, -9.64335932602174E-02_wp, -8.93115661727134E-20_wp, &
      & -3.21307757287911E-18_wp, -1.32631821052136E-19_wp,  4.90158691040716E-02_wp, &
      & -4.90158691040716E-02_wp,  4.90158691040716E-02_wp, -4.90158691040716E-02_wp, &
      &  3.00369861789520E-20_wp, -5.13171156572083E-20_wp, -2.14053789319158E-19_wp, &
      &  5.13171156572083E-20_wp, -5.20084261335142E-19_wp, -8.93115661727134E-20_wp, &
      & -6.32375714785341E-02_wp, -8.68712868561810E-20_wp,  8.16699840135520E-21_wp, &
      &  4.85318612398898E-20_wp,  4.85318612398898E-20_wp,  7.36781018422772E-20_wp, &
      &  0.00000000000000E+00_wp,  7.40528452554017E-20_wp,  3.86788348521085E-02_wp, &
      &  3.32226456377319E-18_wp, -1.68344734676473E-18_wp,  6.70790940818228E-18_wp, &
      & -3.21307757287911E-18_wp, -8.68712868561810E-20_wp, -9.64335932602174E-02_wp, &
      & -6.37231455652560E-20_wp,  4.90158691040716E-02_wp, -4.90158691040716E-02_wp, &
      & -4.90158691040716E-02_wp,  4.90158691040716E-02_wp,  5.59128352808655E-21_wp, &
      &  1.60723274315891E-19_wp,  2.50703088147650E-20_wp,  2.73085288618749E-19_wp, &
      &  1.06590766475164E-20_wp, -1.32631821052136E-19_wp,  8.16699840135520E-21_wp, &
      & -6.37231455652560E-20_wp, -6.32375714785342E-02_wp, -1.94945919560809E-21_wp, &
      & -7.14731447325590E-22_wp, -7.14731447325590E-22_wp,  0.00000000000000E+00_wp, &
      & -1.12822566489047E-01_wp, -5.56863254485502E-02_wp,  5.56863254485502E-02_wp, &
      & -5.56863254485502E-02_wp, -4.90158691040716E-02_wp,  4.90158691040716E-02_wp, &
      &  4.85318612398898E-20_wp,  4.90158691040716E-02_wp, -1.94945919560809E-21_wp, &
      & -1.91835946538684E-01_wp, -2.87572278366319E-02_wp, -2.87572278366319E-02_wp, &
      & -2.87572278366319E-02_wp, -1.12822566489047E-01_wp,  5.56863254485502E-02_wp, &
      &  5.56863254485502E-02_wp,  5.56863254485502E-02_wp, -4.90158691040716E-02_wp, &
      & -4.90158691040716E-02_wp,  4.85318612398898E-20_wp, -4.90158691040716E-02_wp, &
      & -7.14731447325590E-22_wp, -2.87572278366319E-02_wp, -1.91835946538684E-01_wp, &
      & -2.87572278366319E-02_wp, -2.87572278366319E-02_wp, -1.12822566489047E-01_wp, &
      &  5.56863254485502E-02_wp, -5.56863254485502E-02_wp, -5.56863254485502E-02_wp, &
      &  4.90158691040716E-02_wp,  4.90158691040716E-02_wp,  7.36781018422772E-20_wp, &
      & -4.90158691040716E-02_wp, -7.14731447325590E-22_wp, -2.87572278366319E-02_wp, &
      & -2.87572278366319E-02_wp, -1.91835946538684E-01_wp, -2.87572278366319E-02_wp, &
      & -1.12822566489047E-01_wp, -5.56863254485502E-02_wp, -5.56863254485502E-02_wp, &
      &  5.56863254485502E-02_wp,  4.90158691040716E-02_wp, -4.90158691040716E-02_wp, &
      &  0.00000000000000E+00_wp,  4.90158691040716E-02_wp,  0.00000000000000E+00_wp, &
      & -2.87572278366319E-02_wp, -2.87572278366319E-02_wp, -2.87572278366319E-02_wp, &
      & -1.91835946538684E-01_wp], shape(acp))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_acp_gen(error, mol, make_qvszp_basis, acp, thr_in=thr1)

end subroutine test_acp_gxtb_sih4

subroutine test_acp_gxtb_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 43
   real(wp), parameter :: acp(nao, nao) = reshape([&
      & -1.25828652907136E-01_wp,  1.24385765686107E-02_wp, -1.77136928976386E-02_wp, &
      &  7.42971053606532E-03_wp,  1.97217497712655E-02_wp, -4.30229044215817E-02_wp, &
      &  3.45343904532766E-02_wp, -2.97331744855494E-02_wp, -5.75053457822863E-03_wp, &
      &  5.54580571228155E-03_wp,  2.03103843447941E-02_wp,  1.40169576966735E-03_wp, &
      &  5.47968043558296E-03_wp,  1.07278097705703E-03_wp, -7.48017735275227E-03_wp, &
      & -2.06397041668984E-02_wp, -1.48559589744239E-01_wp, -1.52628625504330E-02_wp, &
      & -9.78572901136424E-03_wp, -1.04613450568844E-02_wp,  1.22188631365986E-04_wp, &
      &  1.05164283044589E-04_wp, -3.33193594682712E-05_wp,  6.97101537810912E-05_wp, &
      & -5.03190637760557E-05_wp, -1.48432467164886E-01_wp,  1.59789345408535E-02_wp, &
      &  7.86090622740261E-03_wp, -1.07790943066180E-02_wp, -1.27665197584264E-04_wp, &
      &  9.09045560501680E-05_wp, -4.86414400685599E-05_wp, -6.81623384664728E-05_wp, &
      & -4.76424207175327E-05_wp, -1.48717408713262E-01_wp, -4.86158090263790E-03_wp, &
      &  7.85398826305840E-03_wp,  1.87095261785214E-02_wp, -7.63846793567684E-05_wp, &
      & -3.71182533473672E-05_wp, -5.22708318669255E-05_wp,  1.12778748110421E-04_wp, &
      &  1.23534957245937E-04_wp,  1.24385765686107E-02_wp, -1.03125078064864E-01_wp, &
      & -4.40065641085201E-02_wp,  2.01472907605830E-02_wp,  5.22883100814775E-02_wp, &
      & -6.86528702046933E-04_wp,  1.68378827675317E-04_wp,  4.01877002667501E-02_wp, &
      &  9.32181686759718E-03_wp, -4.90293810557277E-03_wp, -2.42831588212534E-03_wp, &
      &  2.27643597709559E-03_wp,  9.42227688265525E-03_wp, -1.91107916858107E-03_wp, &
      &  8.76379992310948E-03_wp, -1.27635303758120E-06_wp,  1.57871737258077E-01_wp, &
      & -7.80735266568450E-03_wp,  7.22727666559655E-03_wp,  8.37930611992706E-03_wp, &
      & -8.56685031261012E-04_wp, -7.94253918522852E-04_wp,  6.60215897087348E-04_wp, &
      & -7.12321790163897E-05_wp,  1.13998032583339E-03_wp, -1.63766016547088E-01_wp, &
      & -7.04167861250398E-03_wp,  6.07118885182019E-03_wp, -8.60065585989530E-03_wp, &
      & -9.06402463743951E-04_wp,  6.71013904079389E-04_wp, -7.05282356721445E-04_wp, &
      & -6.91061363050009E-05_wp, -1.18133903391043E-03_wp,  5.19384684456390E-02_wp, &
      & -1.85357592391555E-02_wp, -2.46398559310288E-03_wp, -4.87029200019459E-03_wp, &
      &  1.35108882652281E-03_wp,  5.74648679544558E-04_wp,  2.25525625276539E-04_wp, &
      & -2.81128495367743E-05_wp,  3.17929622266108E-04_wp, -1.77136928976386E-02_wp, &
      & -4.40065641085201E-02_wp, -4.78223159557428E-02_wp, -3.02717737462205E-02_wp, &
      &  4.01877002667501E-02_wp,  8.48908477617474E-03_wp,  6.37132453183005E-03_wp, &
      &  5.62167756631225E-03_wp, -1.47809386836961E-02_wp,  9.10005187947799E-04_wp, &
      &  1.24827439890850E-03_wp, -1.10936319942318E-03_wp,  9.49755780159458E-03_wp, &
      & -7.26067369823746E-04_wp, -1.82827927482483E-04_wp, -2.33673007575579E-03_wp, &
      &  9.86516695809701E-02_wp,  6.97161013380999E-03_wp, -1.48556373913106E-02_wp, &
      &  4.71828156485312E-03_wp, -7.91839346167452E-05_wp, -1.15634621595344E-03_wp, &
      & -7.57820556200692E-04_wp, -7.88063333265930E-04_wp,  3.05027484424016E-05_wp, &
      & -8.24141497514839E-02_wp,  6.26413883650035E-03_wp, -1.60883864897687E-02_wp, &
      & -4.72850772381024E-03_wp, -6.41600226156657E-05_wp,  1.17784079475200E-03_wp, &
      &  6.26371532742956E-04_wp, -8.03926132406255E-04_wp, -3.18324443953198E-05_wp, &
      & -8.19135677852889E-02_wp, -2.40084660452909E-03_wp, -1.61978251961245E-02_wp, &
      &  7.37698144704994E-03_wp, -3.10768840010638E-05_wp, -3.69738454192088E-04_wp, &
      &  6.19853361974670E-04_wp,  1.37916365829163E-03_wp,  6.79885084012242E-05_wp, &
      &  7.42971053606532E-03_wp,  2.01472907605830E-02_wp, -3.02717737462205E-02_wp, &
      & -1.14668520935123E-01_wp,  1.75192609593194E-02_wp,  4.01877002667501E-02_wp, &
      &  6.19148128617498E-04_wp, -3.02484060694392E-02_wp, -4.77390285313419E-02_wp, &
      &  6.97153239996482E-03_wp, -6.53475154935186E-03_wp, -1.91107916858107E-03_wp, &
      &  6.53322121931173E-03_wp,  3.66690494229969E-03_wp, -8.15211223434213E-03_wp, &
      &  2.80929363497976E-03_wp,  1.07464329218713E-01_wp,  8.30563136751999E-03_wp, &
      &  4.84627947279294E-03_wp, -1.44542530038716E-02_wp, -1.17176664909735E-03_wp, &
      & -7.34193101819050E-05_wp,  4.50111932224473E-04_wp, -7.39453736940471E-04_wp, &
      & -7.17014908185187E-04_wp,  1.12212689149786E-01_wp, -8.78227953469313E-03_wp, &
      & -4.68770050091010E-03_wp, -1.37758758940700E-02_wp,  1.22445808457141E-03_wp, &
      & -6.17865918406354E-05_wp,  4.85201490650374E-04_wp,  6.11932520430462E-04_wp, &
      & -7.41690819306007E-04_wp, -1.91841763102795E-01_wp, -4.61547832061718E-03_wp, &
      &  7.20830889488899E-03_wp, -2.41994025059087E-03_wp, -4.46183006989497E-04_wp, &
      & -3.76499538606373E-05_wp, -8.28982338224282E-04_wp,  7.06812891483817E-04_wp, &
      &  1.47192601575797E-03_wp,  1.97217497712655E-02_wp,  5.22883100814775E-02_wp, &
      &  4.01877002667501E-02_wp,  1.75192609593194E-02_wp, -8.06217206026291E-02_wp, &
      & -5.80825674682523E-03_wp, -9.89538422295905E-03_wp, -1.53616113922316E-02_wp, &
      &  9.95961560349825E-03_wp,  2.72056667116756E-04_wp, -2.30276967664078E-03_wp, &
      & -1.64375338953638E-03_wp, -1.10570461441990E-02_wp, -1.05408863219153E-04_wp, &
      & -1.42736845948091E-03_wp,  2.59074722693840E-03_wp, -1.02847833133567E-01_wp, &
      & -2.08994563679602E-02_wp, -2.18263312628109E-02_wp, -3.11905791900308E-03_wp, &
      & -1.29180707287719E-03_wp, -4.19761610355791E-03_wp, -2.31799411400137E-03_wp, &
      & -9.84243160227329E-04_wp,  2.82533063802975E-03_wp,  1.11497753814434E-01_wp, &
      & -2.50960905014986E-02_wp, -1.99721330976610E-02_wp,  6.18070892689290E-03_wp, &
      & -1.88669284060823E-03_wp,  4.10301463748736E-03_wp,  1.32888425998985E-03_wp, &
      & -1.21608001627479E-03_wp, -3.25380930264990E-03_wp,  6.05862979021574E-02_wp, &
      & -1.80352012943086E-02_wp, -1.06756185124730E-02_wp, -1.87257954244266E-02_wp, &
      &  4.08949521173900E-03_wp,  2.46855923047715E-03_wp,  7.06493130153477E-04_wp, &
      &  2.94472814877927E-03_wp,  4.23605560773337E-03_wp, -4.30229044215817E-02_wp, &
      & -6.86528702046933E-04_wp,  8.48908477617474E-03_wp,  4.01877002667501E-02_wp, &
      & -5.80825674682523E-03_wp, -5.79808253301321E-02_wp,  7.79737619592565E-03_wp, &
      &  2.99193424794528E-03_wp,  2.05803813098816E-02_wp, -1.44941102813819E-03_wp, &
      &  9.49222073278328E-03_wp,  1.08880066238596E-03_wp, -1.54558072974832E-03_wp, &
      &  5.16572890897420E-04_wp,  6.76951820004811E-04_wp, -8.29191001564372E-03_wp, &
      & -9.47472663004035E-02_wp, -2.02806970805528E-02_wp,  1.16416920603337E-04_wp, &
      & -2.23179396782681E-02_wp, -4.21441236835231E-03_wp, -5.86516156535148E-04_wp, &
      &  3.91693362204840E-03_wp, -3.20574242868777E-04_wp, -8.85178950927925E-04_wp, &
      & -8.20109191272605E-02_wp,  1.91447457326663E-02_wp, -6.45123301936370E-03_wp, &
      & -2.02431039231869E-02_wp,  4.11569522960024E-03_wp,  6.88999474022867E-04_wp, &
      &  4.25138510531403E-03_wp, -9.99216266148978E-04_wp, -6.16239968034470E-04_wp, &
      &  2.56311359234808E-02_wp, -8.04688333086930E-03_wp,  2.55469425915402E-03_wp, &
      & -9.94639913468045E-03_wp,  2.47736866628175E-03_wp, -6.42098560522297E-04_wp, &
      & -1.34070192851108E-03_wp, -5.79976700283422E-04_wp,  2.74971980733658E-03_wp, &
      &  3.45343904532766E-02_wp,  1.68378827675317E-04_wp,  6.37132453183005E-03_wp, &
      &  6.19148128617498E-04_wp, -9.89538422295905E-03_wp,  7.79737619592565E-03_wp, &
      & -3.71226945570069E-02_wp,  5.50359194596325E-03_wp,  3.25234079395653E-03_wp, &
      & -1.93136502609584E-03_wp, -4.56760986183595E-03_wp,  8.97074285440113E-05_wp, &
      & -2.07396138882269E-03_wp, -6.27154089895719E-06_wp,  1.67733823165555E-03_wp, &
      &  7.03971476153940E-03_wp,  2.99046801880864E-02_wp, -1.26740740840466E-03_wp, &
      &  2.18786769454792E-02_wp, -8.32323794306251E-04_wp, -2.30709328096433E-03_wp, &
      &  3.92555859673508E-03_wp,  2.73915880650123E-03_wp,  2.67415049357420E-03_wp, &
      &  9.10403575227964E-04_wp,  4.53850314769518E-02_wp, -4.21782272230626E-03_wp, &
      & -2.09148271449537E-02_wp,  2.83073985067648E-03_wp,  1.32573287532450E-03_wp, &
      &  4.25590467462933E-03_wp,  1.28748004724623E-03_wp, -2.91212401178862E-03_wp, &
      &  5.09071225369334E-04_wp,  4.59573086010083E-02_wp,  1.31473819328767E-03_wp, &
      & -2.09683094066007E-02_wp, -5.11296154237426E-03_wp,  7.11359213551861E-04_wp, &
      & -1.34803332410466E-03_wp,  1.24630969923380E-03_wp,  4.99691238869995E-03_wp, &
      & -1.20802612207777E-03_wp, -2.97331744855494E-02_wp,  4.01877002667501E-02_wp, &
      &  5.62167756631225E-03_wp, -3.02484060694392E-02_wp, -1.53616113922316E-02_wp, &
      &  2.99193424794528E-03_wp,  5.50359194596325E-03_wp, -5.89258856925753E-02_wp, &
      & -1.91304061719755E-02_wp,  5.96102708186324E-03_wp,  3.78309926828190E-03_wp, &
      &  5.16572890897420E-04_wp, -9.77292271465542E-04_wp,  6.95689619568503E-04_wp, &
      & -7.48665765736813E-03_wp, -3.80232209766427E-03_wp, -6.45277323661452E-02_wp, &
      & -2.23574438679134E-02_wp,  3.60494438675301E-05_wp, -2.72035265480656E-03_wp, &
      & -1.00415395221762E-03_wp, -3.22290077597349E-04_wp,  2.66683529399740E-03_wp, &
      & -3.34509372056147E-04_wp,  4.13411728764280E-03_wp,  5.59086836779910E-02_wp, &
      & -1.95337410242953E-02_wp,  4.83712284118301E-03_wp,  2.49604078425354E-03_wp, &
      & -1.20237553574870E-03_wp, -1.00902322331104E-03_wp, -2.90231460734150E-03_wp, &
      & -8.07477642097572E-05_wp, -3.70359812207770E-03_wp, -9.56773554417690E-02_wp, &
      & -1.07070837962044E-02_wp, -7.96711611455442E-03_wp,  2.97198107743107E-02_wp, &
      &  2.94623206270063E-03_wp, -5.67818610146154E-04_wp,  4.98846655687548E-03_wp, &
      &  1.30621192663252E-03_wp, -3.16675114990010E-03_wp, -5.75053457822863E-03_wp, &
      &  9.32181686759718E-03_wp, -1.47809386836961E-02_wp, -4.77390285313419E-02_wp, &
      &  9.95961560349825E-03_wp,  2.05803813098816E-02_wp,  3.25234079395653E-03_wp, &
      & -1.91304061719755E-02_wp, -5.86912462448705E-02_wp,  3.87465639185594E-03_wp, &
      & -1.42736845948091E-03_wp, -1.07539723744054E-03_wp,  4.07161612962649E-03_wp, &
      &  2.43182227708096E-03_wp, -5.18448544188334E-03_wp,  3.48086167753146E-04_wp, &
      &  4.05137620115022E-02_wp, -6.82044720016571E-03_wp,  8.60555102193064E-03_wp, &
      &  2.33602382114053E-02_wp,  2.82231112983009E-03_wp, -8.99081510915722E-04_wp, &
      &  9.13004186632757E-04_wp,  4.13799215953495E-03_wp,  4.76064067087658E-03_wp, &
      &  4.31718884395347E-02_wp,  6.22576158241734E-03_wp, -7.89226083399055E-03_wp, &
      &  2.50602052746388E-02_wp, -3.25567830088934E-03_wp, -6.25358787170301E-04_wp, &
      &  5.19740015161591E-04_wp, -3.70621769882035E-03_wp,  5.26406528302970E-03_wp, &
      & -1.03758188121652E-01_wp, -1.87303883668798E-02_wp,  1.85641133217241E-02_wp, &
      &  1.79712750515064E-02_wp,  4.24061930791193E-03_wp,  2.76738016020276E-03_wp, &
      & -1.21964200598010E-03_wp, -3.15299445695016E-03_wp, -6.87492098453706E-04_wp, &
      &  5.54580571228155E-03_wp, -4.90293810557277E-03_wp,  9.10005187947799E-04_wp, &
      &  6.97153239996482E-03_wp,  2.72056667116756E-04_wp, -1.44941102813819E-03_wp, &
      & -1.93136502609584E-03_wp,  5.96102708186324E-03_wp,  3.87465639185594E-03_wp, &
      & -3.36052670656921E-03_wp, -8.04561899253917E-04_wp,  2.68608195897778E-04_wp, &
      & -3.02141929664058E-04_wp, -5.17911884825695E-04_wp,  3.92150145560411E-04_wp, &
      &  5.38818778251740E-04_wp,  5.37264316346838E-03_wp,  5.17871156628720E-03_wp, &
      &  1.44616420289608E-03_wp, -5.83910768860982E-03_wp, -2.45740780952653E-04_wp, &
      &  7.12463330965630E-04_wp,  1.80261448015424E-04_wp, -8.30802228313411E-04_wp, &
      & -9.74854593975121E-04_wp, -6.19562905821017E-03_wp,  5.72938304298449E-03_wp, &
      &  1.31370824880511E-03_wp,  6.35220285606474E-03_wp, -2.91280985594157E-04_wp, &
      & -6.67697221533828E-04_wp, -1.52525015161597E-04_wp, -7.31684150583113E-04_wp, &
      &  1.19212714032177E-03_wp,  1.95069941633896E-02_wp, -4.89467044125057E-03_wp, &
      & -4.22866821911025E-03_wp, -5.48604211085857E-03_wp,  5.08592883764901E-04_wp, &
      &  5.68971260800863E-04_wp,  4.71547751468779E-04_wp,  6.04130139874219E-04_wp, &
      &  8.80222799798718E-04_wp,  2.03103843447941E-02_wp, -2.42831588212534E-03_wp, &
      &  1.24827439890850E-03_wp, -6.53475154935186E-03_wp, -2.30276967664078E-03_wp, &
      &  9.49222073278328E-03_wp, -4.56760986183595E-03_wp,  3.78309926828190E-03_wp, &
      & -1.42736845948091E-03_wp, -8.04561899253917E-04_wp, -5.21753749244060E-03_wp, &
      & -2.87967073629453E-04_wp, -5.30717012523834E-04_wp, -1.65147186915758E-04_wp, &
      &  7.66563481990073E-04_wp,  3.15657471315477E-03_wp,  2.86461993455965E-02_wp, &
      &  8.30523164560073E-03_wp,  1.49861704048234E-03_wp,  2.60300972608257E-03_wp, &
      &  4.17220898351571E-04_wp,  2.54978265763439E-04_wp, -5.09962171845208E-04_wp, &
      &  1.66342198195187E-04_wp, -6.95552059523375E-04_wp,  2.58684080549418E-02_wp, &
      & -8.08932227998177E-03_wp,  9.09159946186731E-04_wp,  2.85339220985403E-03_wp, &
      & -4.79850820334331E-04_wp, -9.80566540976684E-05_wp, -7.02262147674337E-04_wp, &
      & -5.31636132362328E-05_wp, -6.73240722701929E-04_wp,  1.39989725878205E-02_wp, &
      & -3.75196695283167E-03_wp,  4.97369373169149E-04_wp, -5.84672123530152E-03_wp, &
      &  7.42244987502995E-04_wp, -2.32593503386459E-04_wp, -3.83932213400042E-04_wp, &
      & -9.81408204540011E-05_wp,  8.68076184347637E-04_wp,  1.40169576966735E-03_wp, &
      &  2.27643597709559E-03_wp, -1.10936319942318E-03_wp, -1.91107916858107E-03_wp, &
      & -1.64375338953638E-03_wp,  1.08880066238596E-03_wp,  8.97074285440113E-05_wp, &
      &  5.16572890897420E-04_wp, -1.07539723744054E-03_wp,  2.68608195897778E-04_wp, &
      & -2.87967073629453E-04_wp, -1.96975145540478E-03_wp, -9.17378014923665E-05_wp, &
      &  1.88383842869318E-04_wp,  4.72235029015719E-05_wp,  2.79454274793133E-04_wp, &
      &  1.07963883617155E-03_wp,  3.20276562779747E-03_wp, -6.85619208876269E-03_wp, &
      &  2.22498238582897E-03_wp,  7.39335574925860E-04_wp, -5.36915132944082E-04_wp, &
      & -8.36160404249023E-04_wp, -7.66955316111852E-04_wp, -9.27727532771920E-05_wp, &
      &  5.42433907494828E-03_wp,  1.42611170790771E-03_wp, -7.36499265690783E-03_wp, &
      & -4.88435399046622E-04_wp,  5.41658068020377E-04_wp,  7.44041314818760E-04_wp, &
      &  4.97896553885498E-04_wp, -9.11735143042169E-04_wp,  1.44759300901617E-04_wp, &
      & -1.78330982432131E-03_wp,  8.04639859241270E-04_wp,  2.32385736417165E-03_wp, &
      & -2.27943884761570E-04_wp, -1.67320135499335E-06_wp, -4.55649357916293E-04_wp, &
      & -1.52361598653761E-04_wp, -4.95378318767215E-04_wp,  2.00393724987738E-04_wp, &
      &  5.47968043558296E-03_wp,  9.42227688265525E-03_wp,  9.49755780159458E-03_wp, &
      &  6.53322121931173E-03_wp, -1.10570461441990E-02_wp, -1.54558072974832E-03_wp, &
      & -2.07396138882269E-03_wp, -9.77292271465542E-04_wp,  4.07161612962649E-03_wp, &
      & -3.02141929664058E-04_wp, -5.30717012523834E-04_wp, -9.17378014923665E-05_wp, &
      & -3.52572969536796E-03_wp, -7.31701695974809E-05_wp,  1.31469147246490E-04_wp, &
      &  8.18956457842527E-04_wp, -1.95120133236498E-02_wp, -3.66220264036705E-03_wp, &
      & -2.60212134365512E-03_wp, -2.48410546752564E-03_wp,  1.12352474585631E-04_wp, &
      & -5.67923214754790E-04_wp,  5.54957628376847E-04_wp, -3.86311955660163E-04_wp, &
      & -4.38823513467838E-05_wp,  1.89870055319542E-02_wp, -4.41843917631192E-03_wp, &
      & -4.58886089322819E-04_wp,  3.15279532772002E-03_wp, -8.21055955018999E-05_wp, &
      &  2.76899969210025E-04_wp, -7.99561565752778E-04_wp, -1.91261386348345E-04_wp, &
      & -2.83030026582609E-05_wp,  1.89963978669526E-02_wp,  1.53068603695454E-03_wp, &
      & -3.80821127377831E-04_wp, -5.22929960303891E-03_wp, -4.84198784867568E-05_wp, &
      & -8.63122520046658E-05_wp, -8.07421860714000E-04_wp,  3.12758283417442E-04_wp, &
      &  7.55965934758805E-05_wp,  1.07278097705703E-03_wp, -1.91107916858107E-03_wp, &
      & -7.26067369823746E-04_wp,  3.66690494229969E-03_wp, -1.05408863219153E-04_wp, &
      &  5.16572890897420E-04_wp, -6.27154089895719E-06_wp,  6.95689619568503E-04_wp, &
      &  2.43182227708096E-03_wp, -5.17911884825695E-04_wp, -1.65147186915758E-04_wp, &
      &  1.88383842869318E-04_wp, -7.31701695974809E-05_wp, -2.07967450789482E-03_wp, &
      &  1.57695902195876E-04_wp,  1.51260887793367E-05_wp,  7.37880931550436E-04_wp, &
      &  2.23188910975574E-03_wp, -4.66316030671915E-03_wp,  1.45432900462871E-03_wp, &
      &  3.58004835715111E-04_wp, -7.66595814168157E-04_wp, -5.69555866987976E-04_wp, &
      &  6.74785502036212E-05_wp, -4.32516812518343E-04_wp, -3.71779562676939E-03_wp, &
      & -4.81654918573338E-04_wp,  5.02485435363909E-03_wp,  1.04462637504234E-03_wp, &
      & -3.22102405552695E-04_wp, -9.11517062014260E-04_wp, -3.39366838079932E-04_wp, &
      &  3.67496148663196E-05_wp, -2.21518236705814E-04_wp,  6.58433744077556E-03_wp, &
      & -2.36773492516553E-04_wp, -8.65683524767902E-03_wp,  1.61616869778505E-03_wp, &
      &  3.48893625053997E-04_wp, -4.95519435157707E-04_wp,  5.69311471106019E-04_wp, &
      &  1.24372303668858E-03_wp, -5.46213445312269E-04_wp, -7.48017735275227E-03_wp, &
      &  8.76379992310948E-03_wp, -1.82827927482483E-04_wp, -8.15211223434213E-03_wp, &
      & -1.42736845948091E-03_wp,  6.76951820004811E-04_wp,  1.67733823165555E-03_wp, &
      & -7.48665765736813E-03_wp, -5.18448544188334E-03_wp,  3.92150145560411E-04_wp, &
      &  7.66563481990073E-04_wp,  4.72235029015719E-05_wp,  1.31469147246490E-04_wp, &
      &  1.57695902195876E-04_wp, -3.51702346703292E-03_wp, -1.11890867288079E-03_wp, &
      & -1.12826358557136E-02_wp,  9.14556952050571E-04_wp, -5.87976161018276E-04_wp, &
      & -7.16660870806842E-03_wp, -6.95994194829474E-04_wp, -8.97420001800235E-05_wp, &
      &  2.00874780020662E-04_wp, -8.09507511064301E-05_wp, -1.07583796836397E-03_wp, &
      &  9.97448258350704E-03_wp,  6.68246736287945E-04_wp,  4.06258828399510E-04_wp, &
      &  6.16928059515968E-03_wp, -6.65808452996053E-04_wp,  1.31190741707232E-04_wp, &
      & -2.69670773652115E-04_wp,  2.27204629200503E-04_wp,  9.86859904428150E-04_wp, &
      & -2.39000220401368E-02_wp, -4.73706456943321E-03_wp, -9.45802605153256E-04_wp, &
      &  6.47962570529719E-03_wp,  8.60035138174569E-04_wp, -1.58982541838826E-04_wp, &
      &  6.53396386314785E-04_wp,  1.78224438556927E-05_wp, -2.36311205706913E-04_wp, &
      & -2.06397041668984E-02_wp, -1.27635303758120E-06_wp, -2.33673007575579E-03_wp, &
      &  2.80929363497976E-03_wp,  2.59074722693840E-03_wp, -8.29191001564372E-03_wp, &
      &  7.03971476153940E-03_wp, -3.80232209766427E-03_wp,  3.48086167753146E-04_wp, &
      &  5.38818778251740E-04_wp,  3.15657471315477E-03_wp,  2.79454274793133E-04_wp, &
      &  8.18956457842527E-04_wp,  1.51260887793367E-05_wp, -1.11890867288079E-03_wp, &
      & -5.25881886891614E-03_wp, -2.37104327383569E-02_wp, -2.56339518965953E-03_wp, &
      & -6.31723296520610E-03_wp, -3.86424802340436E-03_wp, -6.10647191805846E-04_wp, &
      & -3.01735414716944E-04_wp, -7.97436904239392E-04_wp, -5.02968082214876E-04_wp, &
      & -1.76080352598168E-06_wp, -2.65792397263913E-02_wp,  3.57898856564963E-03_wp, &
      &  5.86425786116436E-03_wp, -4.84868654522685E-03_wp,  7.44274306189312E-04_wp, &
      & -3.63847703996112E-04_wp, -6.44782905349123E-04_wp,  5.28187001336908E-04_wp, &
      & -1.90259945547730E-05_wp, -1.93171535151488E-02_wp, -7.12057804834572E-03_wp, &
      &  4.27732710879221E-03_wp,  2.49094343283048E-03_wp,  1.02652897622931E-03_wp, &
      &  8.05058288924409E-04_wp, -4.63456011701027E-04_wp, -2.32061344714121E-04_wp, &
      &  1.95318300369122E-05_wp, -1.48559589744239E-01_wp,  1.57871737258077E-01_wp, &
      &  9.86516695809701E-02_wp,  1.07464329218713E-01_wp, -1.02847833133567E-01_wp, &
      & -9.47472663004035E-02_wp,  2.99046801880864E-02_wp, -6.45277323661452E-02_wp, &
      &  4.05137620115022E-02_wp,  5.37264316346838E-03_wp,  2.86461993455965E-02_wp, &
      &  1.07963883617155E-03_wp, -1.95120133236498E-02_wp,  7.37880931550436E-04_wp, &
      & -1.12826358557136E-02_wp, -2.37104327383569E-02_wp, -5.53793951859975E-01_wp, &
      & -1.11918254394686E-03_wp, -7.02064371813275E-04_wp, -7.63036194486640E-04_wp, &
      & -5.64339014514890E-04_wp, -5.18640589867192E-04_wp,  1.63029518139009E-04_wp, &
      & -3.53212928115925E-04_wp,  2.22196893088208E-04_wp, -7.40114343879413E-04_wp, &
      &  2.06873533904093E-03_wp,  1.12162746065931E-03_wp, -5.03899555244931E-04_wp, &
      &  1.16539266895151E-04_wp,  1.36585648709915E-04_wp,  1.52799935208122E-04_wp, &
      &  1.22055995569097E-05_wp, -4.09952475150115E-04_wp, -8.07483219780851E-04_wp, &
      &  3.89132241418222E-04_wp,  1.27237497190952E-03_wp,  2.31089284908727E-03_wp, &
      &  3.02499921490600E-04_wp,  2.99903267180984E-05_wp,  1.43465241471736E-04_wp, &
      &  1.24532319723597E-04_wp,  2.54765110784266E-04_wp, -1.52628625504330E-02_wp, &
      & -7.80735266568450E-03_wp,  6.97161013380999E-03_wp,  8.30563136751999E-03_wp, &
      & -2.08994563679602E-02_wp, -2.02806970805528E-02_wp, -1.26740740840466E-03_wp, &
      & -2.23574438679134E-02_wp, -6.82044720016571E-03_wp,  5.17871156628720E-03_wp, &
      &  8.30523164560073E-03_wp,  3.20276562779747E-03_wp, -3.66220264036705E-03_wp, &
      &  2.23188910975574E-03_wp,  9.14556952050571E-04_wp, -2.56339518965953E-03_wp, &
      & -1.11918254394686E-03_wp, -1.80150630002225E-01_wp, -3.63224106167918E-04_wp, &
      & -3.84447657447443E-04_wp, -4.43215492230896E-04_wp, -4.08043370600211E-04_wp, &
      &  2.83393132902750E-04_wp, -1.03839369198027E-04_wp,  4.73384203402255E-04_wp, &
      & -2.03800981000527E-03_wp,  9.69761865202517E-04_wp,  4.67254661951171E-04_wp, &
      & -6.52924948725766E-05_wp, -6.37373091802286E-05_wp,  2.41929172083578E-05_wp, &
      &  1.03241547261647E-04_wp, -1.80028030863896E-04_wp, -1.18914642450772E-04_wp, &
      & -1.14056920873264E-03_wp, -8.00699133260832E-04_wp,  1.15240489414810E-04_wp, &
      &  7.31771144838515E-04_wp,  5.80135768131265E-04_wp,  1.42166530058440E-04_wp, &
      &  2.34002249906474E-04_wp,  4.87456854387733E-05_wp,  2.61570234221497E-05_wp, &
      & -9.78572901136424E-03_wp,  7.22727666559655E-03_wp, -1.48556373913106E-02_wp, &
      &  4.84627947279294E-03_wp, -2.18263312628109E-02_wp,  1.16416920603337E-04_wp, &
      &  2.18786769454792E-02_wp,  3.60494438675301E-05_wp,  8.60555102193064E-03_wp, &
      &  1.44616420289608E-03_wp,  1.49861704048234E-03_wp, -6.85619208876269E-03_wp, &
      & -2.60212134365512E-03_wp, -4.66316030671915E-03_wp, -5.87976161018276E-04_wp, &
      & -6.31723296520610E-03_wp, -7.02064371813275E-04_wp, -3.63224106167918E-04_wp, &
      & -1.79804732269566E-01_wp, -2.49537614260684E-04_wp, -1.03839369198027E-04_wp, &
      & -5.03489979961417E-04_wp, -2.64368248313250E-04_wp, -3.42969119181269E-04_wp, &
      &  4.07395236248837E-05_wp, -1.19046440199972E-03_wp,  4.86727329054344E-04_wp, &
      &  3.88422248374962E-04_wp, -1.47068427016465E-04_wp, -9.94249783669650E-05_wp, &
      & -5.55654692782165E-05_wp,  1.75119752134599E-05_wp,  1.31390113518260E-04_wp, &
      & -1.71580088279704E-04_wp, -1.33297514956569E-03_wp,  1.09793216130541E-04_wp, &
      &  4.08519928269525E-04_wp,  5.86107550140136E-04_wp, -8.69977306092166E-06_wp, &
      &  8.51516735792755E-05_wp,  2.70034982579483E-05_wp, -9.39921062338574E-05_wp, &
      &  1.88854498309724E-04_wp, -1.04613450568844E-02_wp,  8.37930611992706E-03_wp, &
      &  4.71828156485312E-03_wp, -1.44542530038716E-02_wp, -3.11905791900308E-03_wp, &
      & -2.23179396782681E-02_wp, -8.32323794306251E-04_wp, -2.72035265480656E-03_wp, &
      &  2.33602382114053E-02_wp, -5.83910768860982E-03_wp,  2.60300972608257E-03_wp, &
      &  2.22498238582897E-03_wp, -2.48410546752564E-03_wp,  1.45432900462871E-03_wp, &
      & -7.16660870806842E-03_wp, -3.86424802340436E-03_wp, -7.63036194486640E-04_wp, &
      & -3.84447657447443E-04_wp, -2.49537614260684E-04_wp, -1.79852829573764E-01_wp, &
      & -5.20957081262844E-04_wp, -1.03839369198027E-04_wp,  1.93002097875518E-04_wp, &
      & -3.26564323350443E-04_wp, -2.34043066438151E-04_wp, -4.42875710728222E-04_wp, &
      &  1.08504035030278E-04_wp, -4.89318505993571E-05_wp, -1.06907196538622E-03_wp, &
      &  5.20830775091095E-04_wp,  7.34966413016542E-05_wp,  2.21845557823726E-04_wp, &
      &  2.14620690114397E-04_wp, -3.52467126491866E-04_wp, -2.00684704012956E-03_wp, &
      &  7.18665143275904E-04_wp,  5.79161912118600E-04_wp,  8.32033228600458E-04_wp, &
      & -2.52954185608582E-04_wp, -2.59569178391652E-04_wp, -6.31901251564753E-08_wp, &
      &  4.14683904709787E-05_wp,  1.25565553689502E-04_wp,  1.22188631365986E-04_wp, &
      & -8.56685031261012E-04_wp, -7.91839346167452E-05_wp, -1.17176664909735E-03_wp, &
      & -1.29180707287719E-03_wp, -4.21441236835231E-03_wp, -2.30709328096433E-03_wp, &
      & -1.00415395221762E-03_wp,  2.82231112983009E-03_wp, -2.45740780952653E-04_wp, &
      &  4.17220898351571E-04_wp,  7.39335574925860E-04_wp,  1.12352474585631E-04_wp, &
      &  3.58004835715111E-04_wp, -6.95994194829474E-04_wp, -6.10647191805846E-04_wp, &
      & -5.64339014514890E-04_wp, -4.43215492230896E-04_wp, -1.03839369198027E-04_wp, &
      & -5.20957081262844E-04_wp, -1.30046092302377E-02_wp, -9.76856133957780E-05_wp, &
      &  1.84480124100187E-04_wp, -1.45644993259553E-04_wp, -1.21395161495332E-06_wp, &
      & -9.10284138485143E-05_wp, -4.77048328964455E-05_wp, -1.93938384851233E-04_wp, &
      & -4.67568520710679E-04_wp,  2.19929064886028E-04_wp,  1.33016814104858E-04_wp, &
      &  1.44462502759741E-04_wp,  1.89048164681098E-05_wp, -2.16609483564533E-04_wp, &
      &  2.26111221181509E-04_wp, -1.06376043470459E-04_wp, -1.91109962842684E-04_wp, &
      & -1.14278551480164E-04_wp,  7.01664168355671E-05_wp, -5.94774590154050E-05_wp, &
      &  5.03175909601935E-05_wp,  1.74197214046754E-04_wp,  1.22717877592350E-04_wp, &
      &  1.05164283044589E-04_wp, -7.94253918522852E-04_wp, -1.15634621595344E-03_wp, &
      & -7.34193101819050E-05_wp, -4.19761610355791E-03_wp, -5.86516156535148E-04_wp, &
      &  3.92555859673508E-03_wp, -3.22290077597349E-04_wp, -8.99081510915722E-04_wp, &
      &  7.12463330965630E-04_wp,  2.54978265763439E-04_wp, -5.36915132944082E-04_wp, &
      & -5.67923214754790E-04_wp, -7.66595814168157E-04_wp, -8.97420001800235E-05_wp, &
      & -3.01735414716944E-04_wp, -5.18640589867192E-04_wp, -4.08043370600211E-04_wp, &
      & -5.03489979961417E-04_wp, -1.03839369198027E-04_wp, -9.76856133957780E-05_wp, &
      & -1.29880687997255E-02_wp, -8.58892089079018E-05_wp, -1.58827681120098E-04_wp, &
      &  1.46410945302004E-04_wp,  1.24149336402112E-04_wp, -2.59023898595789E-05_wp, &
      &  3.27571645562904E-05_wp,  7.02519389599410E-06_wp, -9.54597779408670E-05_wp, &
      &  1.96338531478262E-05_wp,  5.33540875159566E-06_wp,  4.44923166720558E-06_wp, &
      & -1.18705354827859E-04_wp, -2.59615426211534E-05_wp, -3.21000382075668E-04_wp, &
      &  9.59724181677427E-05_wp, -2.31127485983110E-05_wp,  1.86432654235490E-04_wp, &
      &  1.15848837570902E-04_wp,  5.70569150448062E-06_wp, -4.04058331363471E-05_wp, &
      &  1.52180609814273E-04_wp, -3.33193594682712E-05_wp,  6.60215897087348E-04_wp, &
      & -7.57820556200692E-04_wp,  4.50111932224473E-04_wp, -2.31799411400137E-03_wp, &
      &  3.91693362204840E-03_wp,  2.73915880650123E-03_wp,  2.66683529399740E-03_wp, &
      &  9.13004186632757E-04_wp,  1.80261448015424E-04_wp, -5.09962171845208E-04_wp, &
      & -8.36160404249023E-04_wp,  5.54957628376847E-04_wp, -5.69555866987976E-04_wp, &
      &  2.00874780020662E-04_wp, -7.97436904239392E-04_wp,  1.63029518139009E-04_wp, &
      &  2.83393132902750E-04_wp, -2.64368248313250E-04_wp,  1.93002097875518E-04_wp, &
      &  1.84480124100187E-04_wp, -8.58892089079018E-05_wp, -1.29002622237164E-02_wp, &
      & -5.84771163103343E-05_wp, -7.26582142215733E-05_wp,  1.37908468629568E-04_wp, &
      & -1.51732984494235E-04_wp,  4.15790697134960E-05_wp,  1.40872631925668E-04_wp, &
      & -1.31220450920107E-04_wp, -2.01402171296678E-05_wp, -8.64747156169791E-05_wp, &
      &  6.17504374997657E-05_wp, -2.25104432070926E-06_wp,  1.35653516694249E-04_wp, &
      &  8.82953817794205E-05_wp,  2.79961019637383E-05_wp, -1.77686280092715E-04_wp, &
      & -1.04252999265368E-04_wp,  5.04740298657699E-05_wp, -7.77665017724439E-05_wp, &
      & -3.50180487991005E-05_wp,  7.66217690306814E-05_wp,  6.97101537810912E-05_wp, &
      & -7.12321790163897E-05_wp, -7.88063333265930E-04_wp, -7.39453736940471E-04_wp, &
      & -9.84243160227329E-04_wp, -3.20574242868777E-04_wp,  2.67415049357420E-03_wp, &
      & -3.34509372056147E-04_wp,  4.13799215953495E-03_wp, -8.30802228313411E-04_wp, &
      &  1.66342198195187E-04_wp, -7.66955316111852E-04_wp, -3.86311955660163E-04_wp, &
      &  6.74785502036212E-05_wp, -8.09507511064301E-05_wp, -5.02968082214876E-04_wp, &
      & -3.53212928115925E-04_wp, -1.03839369198027E-04_wp, -3.42969119181269E-04_wp, &
      & -3.26564323350443E-04_wp, -1.45644993259553E-04_wp, -1.58827681120098E-04_wp, &
      & -5.84771163103343E-05_wp, -1.28629880096263E-02_wp, -1.01197352683435E-04_wp, &
      & -1.04431524664554E-04_wp, -1.12359803560323E-04_wp,  1.08134412507555E-04_wp, &
      & -3.41213365342892E-04_wp,  1.06675956726340E-04_wp, -1.54260416860735E-05_wp, &
      &  3.28872408592862E-06_wp,  1.41968094941072E-04_wp, -2.08505370512871E-04_wp, &
      &  1.39139734809984E-04_wp,  1.65719226559504E-04_wp,  4.69201968911274E-06_wp, &
      & -2.69841363270746E-06_wp, -1.00772340047709E-04_wp, -4.89572880202854E-05_wp, &
      &  4.44639318123372E-06_wp,  2.67181858282058E-05_wp,  9.43211898934151E-05_wp, &
      & -5.03190637760557E-05_wp,  1.13998032583339E-03_wp,  3.05027484424016E-05_wp, &
      & -7.17014908185187E-04_wp,  2.82533063802975E-03_wp, -8.85178950927925E-04_wp, &
      &  9.10403575227964E-04_wp,  4.13411728764280E-03_wp,  4.76064067087658E-03_wp, &
      & -9.74854593975121E-04_wp, -6.95552059523375E-04_wp, -9.27727532771920E-05_wp, &
      & -4.38823513467838E-05_wp, -4.32516812518343E-04_wp, -1.07583796836397E-03_wp, &
      & -1.76080352598168E-06_wp,  2.22196893088208E-04_wp,  4.73384203402255E-04_wp, &
      &  4.07395236248837E-05_wp, -2.34043066438151E-04_wp, -1.21395161495332E-06_wp, &
      &  1.46410945302004E-04_wp, -7.26582142215733E-05_wp, -1.01197352683435E-04_wp, &
      & -1.30071965141969E-02_wp, -4.12678042752354E-04_wp,  1.55111723535769E-04_wp, &
      &  8.24802573331507E-05_wp, -3.76896608118955E-04_wp,  2.48191750190466E-04_wp, &
      & -8.06540525782525E-05_wp,  6.79814886977871E-05_wp,  1.77666580716642E-04_wp, &
      & -1.12345546941984E-05_wp,  3.21712323184730E-04_wp,  5.56755372337909E-04_wp, &
      &  1.10082610129905E-04_wp, -2.28178751641412E-04_wp, -3.68757719580196E-04_wp, &
      & -1.49306817466030E-04_wp, -1.49762625073490E-04_wp, -8.58543413206954E-06_wp, &
      &  6.34735644322809E-05_wp, -1.48432467164886E-01_wp, -1.63766016547088E-01_wp, &
      & -8.24141497514839E-02_wp,  1.12212689149786E-01_wp,  1.11497753814434E-01_wp, &
      & -8.20109191272605E-02_wp,  4.53850314769518E-02_wp,  5.59086836779910E-02_wp, &
      &  4.31718884395347E-02_wp, -6.19562905821017E-03_wp,  2.58684080549418E-02_wp, &
      &  5.42433907494828E-03_wp,  1.89870055319542E-02_wp, -3.71779562676939E-03_wp, &
      &  9.97448258350704E-03_wp, -2.65792397263913E-02_wp, -7.40114343879413E-04_wp, &
      & -2.03800981000527E-03_wp, -1.19046440199972E-03_wp, -4.42875710728222E-04_wp, &
      & -9.10284138485143E-05_wp,  1.24149336402112E-04_wp,  1.37908468629568E-04_wp, &
      & -1.04431524664554E-04_wp, -4.12678042752354E-04_wp, -5.53792622930797E-01_wp, &
      &  1.15987427301119E-03_wp,  5.82061907975335E-04_wp, -7.92853779360320E-04_wp, &
      &  6.10539348685428E-04_wp, -4.48395785182379E-04_wp,  2.48018803026518E-04_wp, &
      &  3.06876876513905E-04_wp,  2.36798499647136E-04_wp, -7.01929220519856E-04_wp, &
      & -9.60758728436176E-04_wp,  4.11958642562427E-04_wp,  1.95546097947158E-03_wp, &
      & -3.02791495369278E-04_wp, -1.70543423858042E-04_wp, -2.25321777298567E-04_wp, &
      & -2.91324079951097E-04_wp, -1.67830034552342E-05_wp,  1.59789345408535E-02_wp, &
      & -7.04167861250398E-03_wp,  6.26413883650035E-03_wp, -8.78227953469313E-03_wp, &
      & -2.50960905014986E-02_wp,  1.91447457326663E-02_wp, -4.21782272230626E-03_wp, &
      & -1.95337410242953E-02_wp,  6.22576158241734E-03_wp,  5.72938304298449E-03_wp, &
      & -8.08932227998177E-03_wp,  1.42611170790771E-03_wp, -4.41843917631192E-03_wp, &
      & -4.81654918573338E-04_wp,  6.68246736287945E-04_wp,  3.57898856564963E-03_wp, &
      &  2.06873533904093E-03_wp,  9.69761865202517E-04_wp,  4.86727329054344E-04_wp, &
      &  1.08504035030278E-04_wp, -4.77048328964455E-05_wp, -2.59023898595789E-05_wp, &
      & -1.51732984494235E-04_wp, -1.12359803560323E-04_wp,  1.55111723535769E-04_wp, &
      &  1.15987427301119E-03_wp, -1.80170892223689E-01_wp, -3.11921138310343E-04_wp, &
      &  4.19784650397748E-04_wp, -4.76154937747586E-04_wp,  3.50145168522627E-04_wp, &
      & -3.19871621020701E-04_wp, -9.39135476357906E-05_wp, -4.94864011712153E-04_wp, &
      &  1.53047453828104E-03_wp,  1.35567286307258E-04_wp, -5.80934424806554E-04_wp, &
      & -8.01154801552224E-04_wp,  1.04187954686724E-04_wp,  6.44561587683816E-05_wp, &
      &  2.52579789461132E-04_wp,  3.65096087868755E-04_wp,  8.34663571047527E-05_wp, &
      &  7.86090622740261E-03_wp,  6.07118885182019E-03_wp, -1.60883864897687E-02_wp, &
      & -4.68770050091010E-03_wp, -1.99721330976610E-02_wp, -6.45123301936370E-03_wp, &
      & -2.09148271449537E-02_wp,  4.83712284118301E-03_wp, -7.89226083399055E-03_wp, &
      &  1.31370824880511E-03_wp,  9.09159946186731E-04_wp, -7.36499265690783E-03_wp, &
      & -4.58886089322819E-04_wp,  5.02485435363909E-03_wp,  4.06258828399510E-04_wp, &
      &  5.86425786116436E-03_wp,  1.12162746065931E-03_wp,  4.67254661951171E-04_wp, &
      &  3.88422248374962E-04_wp, -4.89318505993571E-05_wp, -1.93938384851233E-04_wp, &
      &  3.27571645562904E-05_wp,  4.15790697134960E-05_wp,  1.08134412507555E-04_wp, &
      &  8.24802573331507E-05_wp,  5.82061907975335E-04_wp, -3.11921138310343E-04_wp, &
      & -1.79702568561196E-01_wp,  2.05280008700651E-04_wp, -9.39135476357907E-05_wp, &
      &  4.91288757338998E-04_wp,  2.06675128813526E-04_wp, -3.36071342566094E-04_wp, &
      & -3.68741443676301E-05_wp,  4.18898288020449E-04_wp, -4.13645402313573E-04_wp, &
      & -5.59037279757227E-04_wp, -4.73036088619893E-04_wp,  2.81779188724190E-04_wp, &
      &  4.48366532573882E-05_wp, -6.39039615159346E-06_wp,  4.00831729080692E-04_wp, &
      &  2.81883813580586E-04_wp, -1.07790943066180E-02_wp, -8.60065585989530E-03_wp, &
      & -4.72850772381024E-03_wp, -1.37758758940700E-02_wp,  6.18070892689290E-03_wp, &
      & -2.02431039231869E-02_wp,  2.83073985067648E-03_wp,  2.49604078425354E-03_wp, &
      &  2.50602052746388E-02_wp,  6.35220285606474E-03_wp,  2.85339220985403E-03_wp, &
      & -4.88435399046622E-04_wp,  3.15279532772002E-03_wp,  1.04462637504234E-03_wp, &
      &  6.16928059515968E-03_wp, -4.84868654522685E-03_wp, -5.03899555244931E-04_wp, &
      & -6.52924948725766E-05_wp, -1.47068427016465E-04_wp, -1.06907196538622E-03_wp, &
      & -4.67568520710679E-04_wp,  7.02519389599410E-06_wp,  1.40872631925668E-04_wp, &
      & -3.41213365342892E-04_wp, -3.76896608118955E-04_wp, -7.92853779360320E-04_wp, &
      &  4.19784650397748E-04_wp,  2.05280008700651E-04_wp, -1.79840664583753E-01_wp, &
      &  5.50458645134115E-04_wp, -9.39135476357907E-05_wp,  2.19163405064857E-04_wp, &
      &  2.76396879787367E-04_wp, -2.39518557550638E-04_wp, -1.55471638648371E-03_wp, &
      & -5.53292418186806E-04_wp, -2.36632916520753E-04_wp,  6.09124162910329E-04_wp, &
      & -3.75362383395788E-05_wp, -2.03337453504803E-04_wp,  4.59221399345166E-05_wp, &
      &  1.67180345932289E-04_wp,  1.59750178698940E-04_wp, -1.27665197584264E-04_wp, &
      & -9.06402463743951E-04_wp, -6.41600226156657E-05_wp,  1.22445808457141E-03_wp, &
      & -1.88669284060823E-03_wp,  4.11569522960024E-03_wp,  1.32573287532450E-03_wp, &
      & -1.20237553574870E-03_wp, -3.25567830088934E-03_wp, -2.91280985594157E-04_wp, &
      & -4.79850820334331E-04_wp,  5.41658068020377E-04_wp, -8.21055955018999E-05_wp, &
      & -3.22102405552695E-04_wp, -6.65808452996053E-04_wp,  7.44274306189312E-04_wp, &
      &  1.16539266895151E-04_wp, -6.37373091802286E-05_wp, -9.94249783669650E-05_wp, &
      &  5.20830775091095E-04_wp,  2.19929064886028E-04_wp, -9.54597779408670E-05_wp, &
      & -1.31220450920107E-04_wp,  1.06675956726340E-04_wp,  2.48191750190466E-04_wp, &
      &  6.10539348685428E-04_wp, -4.76154937747586E-04_wp, -9.39135476357907E-05_wp, &
      &  5.50458645134115E-04_wp, -1.30340655785626E-02_wp,  8.52715529897453E-05_wp, &
      & -1.99091824833555E-04_wp, -1.25956457946074E-04_wp,  1.00194068832734E-06_wp, &
      & -2.26847551232843E-04_wp, -1.11987136516778E-04_wp,  4.37406360677530E-05_wp, &
      &  7.07420244329337E-05_wp,  1.30210820842156E-04_wp,  1.20760079179363E-04_wp, &
      &  8.27485398535997E-05_wp, -3.92252629812274E-05_wp, -8.15014124206673E-05_wp, &
      &  9.09045560501680E-05_wp,  6.71013904079389E-04_wp,  1.17784079475200E-03_wp, &
      & -6.17865918406354E-05_wp,  4.10301463748736E-03_wp,  6.88999474022867E-04_wp, &
      &  4.25590467462933E-03_wp, -1.00902322331104E-03_wp, -6.25358787170301E-04_wp, &
      & -6.67697221533828E-04_wp, -9.80566540976684E-05_wp,  7.44041314818760E-04_wp, &
      &  2.76899969210025E-04_wp, -9.11517062014260E-04_wp,  1.31190741707232E-04_wp, &
      & -3.63847703996112E-04_wp,  1.36585648709915E-04_wp,  2.41929172083578E-05_wp, &
      & -5.55654692782165E-05_wp,  7.34966413016542E-05_wp,  1.33016814104858E-04_wp, &
      &  1.96338531478262E-05_wp, -2.01402171296678E-05_wp, -1.54260416860735E-05_wp, &
      & -8.06540525782525E-05_wp, -4.48395785182379E-04_wp,  3.50145168522627E-04_wp, &
      &  4.91288757338998E-04_wp, -9.39135476357907E-05_wp,  8.52715529897453E-05_wp, &
      & -1.29805546557839E-02_wp, -7.41913132212086E-05_wp,  1.72337282053953E-04_wp, &
      &  1.26577758057487E-04_wp, -3.29810591225137E-04_wp,  2.07096536554246E-04_wp, &
      &  3.93156931611294E-04_wp,  3.14016845264834E-04_wp, -1.36118247439483E-04_wp, &
      & -2.27003362745431E-05_wp, -3.93667800219882E-05_wp, -2.93406929073157E-04_wp, &
      & -1.73309086839298E-04_wp, -4.86414400685599E-05_wp, -7.05282356721445E-04_wp, &
      &  6.26371532742956E-04_wp,  4.85201490650374E-04_wp,  1.32888425998985E-03_wp, &
      &  4.25138510531403E-03_wp,  1.28748004724623E-03_wp, -2.90231460734150E-03_wp, &
      &  5.19740015161591E-04_wp, -1.52525015161597E-04_wp, -7.02262147674337E-04_wp, &
      &  4.97896553885498E-04_wp, -7.99561565752778E-04_wp, -3.39366838079932E-04_wp, &
      & -2.69670773652115E-04_wp, -6.44782905349123E-04_wp,  1.52799935208122E-04_wp, &
      &  1.03241547261647E-04_wp,  1.75119752134599E-05_wp,  2.21845557823726E-04_wp, &
      &  1.44462502759741E-04_wp,  5.33540875159566E-06_wp, -8.64747156169791E-05_wp, &
      &  3.28872408592862E-06_wp,  6.79814886977871E-05_wp,  2.48018803026518E-04_wp, &
      & -3.19871621020701E-04_wp,  2.06675128813526E-04_wp,  2.19163405064857E-04_wp, &
      & -1.99091824833555E-04_wp, -7.41913132212086E-05_wp, -1.28736435150089E-02_wp, &
      &  5.08155300300249E-05_wp, -7.71773940436147E-05_wp, -2.25802427606423E-04_wp, &
      &  1.27860568480973E-04_wp, -5.76243885410779E-06_wp,  2.23071756314345E-04_wp, &
      & -3.82249450093695E-05_wp,  2.90340489996125E-06_wp,  1.19565038471065E-04_wp, &
      & -4.22643021920042E-05_wp, -1.79297232198718E-04_wp, -6.81623384664728E-05_wp, &
      & -6.91061363050009E-05_wp, -8.03926132406255E-04_wp,  6.11932520430462E-04_wp, &
      & -1.21608001627479E-03_wp, -9.99216266148978E-04_wp, -2.91212401178862E-03_wp, &
      & -8.07477642097572E-05_wp, -3.70621769882035E-03_wp, -7.31684150583113E-04_wp, &
      & -5.31636132362328E-05_wp, -9.11735143042169E-04_wp, -1.91261386348345E-04_wp, &
      &  3.67496148663196E-05_wp,  2.27204629200503E-04_wp,  5.28187001336908E-04_wp, &
      &  1.22055995569097E-05_wp, -1.80028030863896E-04_wp,  1.31390113518260E-04_wp, &
      &  2.14620690114397E-04_wp,  1.89048164681098E-05_wp,  4.44923166720558E-06_wp, &
      &  6.17504374997657E-05_wp,  1.41968094941072E-04_wp,  1.77666580716642E-04_wp, &
      &  3.06876876513905E-04_wp, -9.39135476357906E-05_wp, -3.36071342566094E-04_wp, &
      &  2.76396879787367E-04_wp, -1.25956457946074E-04_wp,  1.72337282053953E-04_wp, &
      &  5.08155300300249E-05_wp, -1.28468305442905E-02_wp,  8.76141009797303E-05_wp, &
      & -6.60975683924875E-05_wp, -2.52398951221758E-04_wp, -9.01132734772037E-05_wp, &
      &  2.56892502776247E-05_wp,  1.84787897688217E-04_wp,  6.81990600707138E-05_wp, &
      &  1.69443037068557E-05_wp,  5.61998503299879E-05_wp,  2.71098448232089E-05_wp, &
      & -4.76424207175327E-05_wp, -1.18133903391043E-03_wp, -3.18324443953198E-05_wp, &
      & -7.41690819306007E-04_wp, -3.25380930264990E-03_wp, -6.16239968034470E-04_wp, &
      &  5.09071225369334E-04_wp, -3.70359812207770E-03_wp,  5.26406528302970E-03_wp, &
      &  1.19212714032177E-03_wp, -6.73240722701929E-04_wp,  1.44759300901617E-04_wp, &
      & -2.83030026582609E-05_wp, -2.21518236705814E-04_wp,  9.86859904428150E-04_wp, &
      & -1.90259945547730E-05_wp, -4.09952475150115E-04_wp, -1.18914642450772E-04_wp, &
      & -1.71580088279704E-04_wp, -3.52467126491866E-04_wp, -2.16609483564533E-04_wp, &
      & -1.18705354827859E-04_wp, -2.25104432070926E-06_wp, -2.08505370512871E-04_wp, &
      & -1.12345546941984E-05_wp,  2.36798499647136E-04_wp, -4.94864011712153E-04_wp, &
      & -3.68741443676301E-05_wp, -2.39518557550638E-04_wp,  1.00194068832734E-06_wp, &
      &  1.26577758057487E-04_wp, -7.71773940436147E-05_wp,  8.76141009797303E-05_wp, &
      & -1.30362106832147E-02_wp,  2.01165203854045E-04_wp, -1.38730019743377E-04_wp, &
      & -3.95569062796219E-04_wp, -8.90364084467836E-05_wp,  3.76747633109434E-05_wp, &
      & -5.62726044513473E-05_wp,  1.63390683549491E-04_wp,  2.53053704317884E-04_wp, &
      &  7.49922749173758E-05_wp, -1.48717408713262E-01_wp,  5.19384684456390E-02_wp, &
      & -8.19135677852889E-02_wp, -1.91841763102795E-01_wp,  6.05862979021574E-02_wp, &
      &  2.56311359234808E-02_wp,  4.59573086010083E-02_wp, -9.56773554417690E-02_wp, &
      & -1.03758188121652E-01_wp,  1.95069941633896E-02_wp,  1.39989725878205E-02_wp, &
      & -1.78330982432131E-03_wp,  1.89963978669526E-02_wp,  6.58433744077556E-03_wp, &
      & -2.39000220401368E-02_wp, -1.93171535151488E-02_wp, -8.07483219780851E-04_wp, &
      & -1.14056920873264E-03_wp, -1.33297514956569E-03_wp, -2.00684704012956E-03_wp, &
      &  2.26111221181509E-04_wp, -2.59615426211534E-05_wp,  1.35653516694249E-04_wp, &
      &  1.39139734809984E-04_wp,  3.21712323184730E-04_wp, -7.01929220519856E-04_wp, &
      &  1.53047453828104E-03_wp,  4.18898288020449E-04_wp, -1.55471638648371E-03_wp, &
      & -2.26847551232843E-04_wp, -3.29810591225137E-04_wp, -2.25802427606423E-04_wp, &
      & -6.60975683924875E-05_wp,  2.01165203854045E-04_wp, -5.53801685267588E-01_wp, &
      & -3.67370071688180E-04_wp,  5.82458611794306E-04_wp,  1.36656275308419E-03_wp, &
      &  3.32029384500713E-04_wp,  1.41546225512443E-04_wp,  2.51544504085765E-04_wp, &
      & -5.24441732338870E-04_wp, -5.69553870423969E-04_wp, -4.86158090263790E-03_wp, &
      & -1.85357592391555E-02_wp, -2.40084660452909E-03_wp, -4.61547832061718E-03_wp, &
      & -1.80352012943086E-02_wp, -8.04688333086930E-03_wp,  1.31473819328767E-03_wp, &
      & -1.07070837962044E-02_wp, -1.87303883668798E-02_wp, -4.89467044125057E-03_wp, &
      & -3.75196695283167E-03_wp,  8.04639859241270E-04_wp,  1.53068603695454E-03_wp, &
      & -2.36773492516553E-04_wp, -4.73706456943321E-03_wp, -7.12057804834572E-03_wp, &
      &  3.89132241418222E-04_wp, -8.00699133260832E-04_wp,  1.09793216130541E-04_wp, &
      &  7.18665143275904E-04_wp, -1.06376043470459E-04_wp, -3.21000382075668E-04_wp, &
      &  8.82953817794205E-05_wp,  1.65719226559504E-04_wp,  5.56755372337909E-04_wp, &
      & -9.60758728436176E-04_wp,  1.35567286307258E-04_wp, -4.13645402313573E-04_wp, &
      & -5.53292418186806E-04_wp, -1.11987136516778E-04_wp,  2.07096536554246E-04_wp, &
      &  1.27860568480973E-04_wp, -2.52398951221758E-04_wp, -1.38730019743377E-04_wp, &
      & -3.67370071688180E-04_wp, -1.79644313328851E-01_wp,  8.89876760142140E-05_wp, &
      &  2.23671497413935E-04_wp,  5.31397248426354E-04_wp,  2.26564801834993E-04_wp, &
      &  1.02054573278130E-04_wp, -4.95035735966572E-05_wp,  8.07806855518498E-05_wp, &
      &  7.85398826305840E-03_wp, -2.46398559310288E-03_wp, -1.61978251961245E-02_wp, &
      &  7.20830889488899E-03_wp, -1.06756185124730E-02_wp,  2.55469425915402E-03_wp, &
      & -2.09683094066007E-02_wp, -7.96711611455442E-03_wp,  1.85641133217241E-02_wp, &
      & -4.22866821911025E-03_wp,  4.97369373169149E-04_wp,  2.32385736417165E-03_wp, &
      & -3.80821127377831E-04_wp, -8.65683524767902E-03_wp, -9.45802605153256E-04_wp, &
      &  4.27732710879221E-03_wp,  1.27237497190952E-03_wp,  1.15240489414810E-04_wp, &
      &  4.08519928269525E-04_wp,  5.79161912118600E-04_wp, -1.91109962842684E-04_wp, &
      &  9.59724181677427E-05_wp,  2.79961019637383E-05_wp,  4.69201968911274E-06_wp, &
      &  1.10082610129905E-04_wp,  4.11958642562427E-04_wp, -5.80934424806554E-04_wp, &
      & -5.59037279757227E-04_wp, -2.36632916520753E-04_wp,  4.37406360677530E-05_wp, &
      &  3.93156931611294E-04_wp, -5.76243885410779E-06_wp, -9.01132734772037E-05_wp, &
      & -3.95569062796219E-04_wp,  5.82458611794306E-04_wp,  8.89876760142140E-05_wp, &
      & -1.79730654584363E-01_wp, -3.62315061415856E-04_wp, -4.95035735966572E-05_wp, &
      & -1.55933615620841E-04_wp,  2.07986907463626E-04_wp,  5.78446720207333E-04_wp, &
      &  8.60225386342648E-05_wp,  1.87095261785214E-02_wp, -4.87029200019459E-03_wp, &
      &  7.37698144704994E-03_wp, -2.41994025059087E-03_wp, -1.87257954244266E-02_wp, &
      & -9.94639913468045E-03_wp, -5.11296154237426E-03_wp,  2.97198107743107E-02_wp, &
      &  1.79712750515064E-02_wp, -5.48604211085857E-03_wp, -5.84672123530152E-03_wp, &
      & -2.27943884761570E-04_wp, -5.22929960303891E-03_wp,  1.61616869778505E-03_wp, &
      &  6.47962570529719E-03_wp,  2.49094343283048E-03_wp,  2.31089284908727E-03_wp, &
      &  7.31771144838515E-04_wp,  5.86107550140136E-04_wp,  8.32033228600458E-04_wp, &
      & -1.14278551480164E-04_wp, -2.31127485983110E-05_wp, -1.77686280092715E-04_wp, &
      & -2.69841363270746E-06_wp, -2.28178751641412E-04_wp,  1.95546097947158E-03_wp, &
      & -8.01154801552224E-04_wp, -4.73036088619893E-04_wp,  6.09124162910329E-04_wp, &
      &  7.07420244329337E-05_wp,  3.14016845264834E-04_wp,  2.23071756314345E-04_wp, &
      &  2.56892502776247E-05_wp, -8.90364084467836E-05_wp,  1.36656275308419E-03_wp, &
      &  2.23671497413935E-04_wp, -3.62315061415856E-04_wp, -1.80429171066847E-01_wp, &
      & -2.51916636131474E-04_wp, -4.95035735966572E-05_wp, -3.77298703194386E-04_wp, &
      &  3.98609879103522E-04_wp,  7.00549995343506E-04_wp, -7.63846793567684E-05_wp, &
      &  1.35108882652281E-03_wp, -3.10768840010638E-05_wp, -4.46183006989497E-04_wp, &
      &  4.08949521173900E-03_wp,  2.47736866628175E-03_wp,  7.11359213551861E-04_wp, &
      &  2.94623206270063E-03_wp,  4.24061930791193E-03_wp,  5.08592883764901E-04_wp, &
      &  7.42244987502995E-04_wp, -1.67320135499335E-06_wp, -4.84198784867568E-05_wp, &
      &  3.48893625053997E-04_wp,  8.60035138174569E-04_wp,  1.02652897622931E-03_wp, &
      &  3.02499921490600E-04_wp,  5.80135768131265E-04_wp, -8.69977306092166E-06_wp, &
      & -2.52954185608582E-04_wp,  7.01664168355671E-05_wp,  1.86432654235490E-04_wp, &
      & -1.04252999265368E-04_wp, -1.00772340047709E-04_wp, -3.68757719580196E-04_wp, &
      & -3.02791495369278E-04_wp,  1.04187954686724E-04_wp,  2.81779188724190E-04_wp, &
      & -3.75362383395788E-05_wp,  1.30210820842156E-04_wp, -1.36118247439483E-04_wp, &
      & -3.82249450093695E-05_wp,  1.84787897688217E-04_wp,  3.76747633109434E-05_wp, &
      &  3.32029384500713E-04_wp,  5.31397248426354E-04_wp, -4.95035735966572E-05_wp, &
      & -2.51916636131474E-04_wp, -1.30380548473772E-02_wp, -1.49051426720838E-04_wp, &
      & -1.08358385538101E-04_wp,  3.79727615859641E-05_wp, -2.65603191426876E-06_wp, &
      & -3.71182533473672E-05_wp,  5.74648679544558E-04_wp, -3.69738454192088E-04_wp, &
      & -3.76499538606373E-05_wp,  2.46855923047715E-03_wp, -6.42098560522297E-04_wp, &
      & -1.34803332410466E-03_wp, -5.67818610146154E-04_wp,  2.76738016020276E-03_wp, &
      &  5.68971260800863E-04_wp, -2.32593503386459E-04_wp, -4.55649357916293E-04_wp, &
      & -8.63122520046658E-05_wp, -4.95519435157707E-04_wp, -1.58982541838826E-04_wp, &
      &  8.05058288924409E-04_wp,  2.99903267180984E-05_wp,  1.42166530058440E-04_wp, &
      &  8.51516735792755E-05_wp, -2.59569178391652E-04_wp, -5.94774590154050E-05_wp, &
      &  1.15848837570902E-04_wp,  5.04740298657699E-05_wp, -4.89572880202854E-05_wp, &
      & -1.49306817466030E-04_wp, -1.70543423858042E-04_wp,  6.44561587683816E-05_wp, &
      &  4.48366532573882E-05_wp, -2.03337453504803E-04_wp,  1.20760079179363E-04_wp, &
      & -2.27003362745431E-05_wp,  2.90340489996125E-06_wp,  6.81990600707138E-05_wp, &
      & -5.62726044513473E-05_wp,  1.41546225512443E-04_wp,  2.26564801834993E-04_wp, &
      & -1.55933615620841E-04_wp, -4.95035735966572E-05_wp, -1.49051426720838E-04_wp, &
      & -1.27523108476812E-02_wp,  2.38644103832816E-05_wp,  9.37490197627848E-05_wp, &
      & -4.16134092031767E-05_wp, -5.22708318669255E-05_wp,  2.25525625276539E-04_wp, &
      &  6.19853361974670E-04_wp, -8.28982338224282E-04_wp,  7.06493130153477E-04_wp, &
      & -1.34070192851108E-03_wp,  1.24630969923380E-03_wp,  4.98846655687548E-03_wp, &
      & -1.21964200598010E-03_wp,  4.71547751468779E-04_wp, -3.83932213400042E-04_wp, &
      & -1.52361598653761E-04_wp, -8.07421860714000E-04_wp,  5.69311471106019E-04_wp, &
      &  6.53396386314785E-04_wp, -4.63456011701027E-04_wp,  1.43465241471736E-04_wp, &
      &  2.34002249906474E-04_wp,  2.70034982579483E-05_wp, -6.31901251564753E-08_wp, &
      &  5.03175909601935E-05_wp,  5.70569150448062E-06_wp, -7.77665017724439E-05_wp, &
      &  4.44639318123372E-06_wp, -1.49762625073490E-04_wp, -2.25321777298567E-04_wp, &
      &  2.52579789461132E-04_wp, -6.39039615159346E-06_wp,  4.59221399345166E-05_wp, &
      &  8.27485398535997E-05_wp, -3.93667800219882E-05_wp,  1.19565038471065E-04_wp, &
      &  1.69443037068557E-05_wp,  1.63390683549491E-04_wp,  2.51544504085765E-04_wp, &
      &  1.02054573278130E-04_wp,  2.07986907463626E-04_wp, -3.77298703194386E-04_wp, &
      & -1.08358385538101E-04_wp,  2.38644103832816E-05_wp, -1.28728570956844E-02_wp, &
      & -8.82402366526675E-05_wp,  1.85758046632624E-04_wp,  1.12778748110421E-04_wp, &
      & -2.81128495367743E-05_wp,  1.37916365829163E-03_wp,  7.06812891483817E-04_wp, &
      &  2.94472814877927E-03_wp, -5.79976700283422E-04_wp,  4.99691238869995E-03_wp, &
      &  1.30621192663252E-03_wp, -3.15299445695016E-03_wp,  6.04130139874219E-04_wp, &
      & -9.81408204540011E-05_wp, -4.95378318767215E-04_wp,  3.12758283417442E-04_wp, &
      &  1.24372303668858E-03_wp,  1.78224438556927E-05_wp, -2.32061344714121E-04_wp, &
      &  1.24532319723597E-04_wp,  4.87456854387733E-05_wp, -9.39921062338574E-05_wp, &
      &  4.14683904709787E-05_wp,  1.74197214046754E-04_wp, -4.04058331363471E-05_wp, &
      & -3.50180487991005E-05_wp,  2.67181858282058E-05_wp, -8.58543413206954E-06_wp, &
      & -2.91324079951097E-04_wp,  3.65096087868755E-04_wp,  4.00831729080692E-04_wp, &
      &  1.67180345932289E-04_wp, -3.92252629812274E-05_wp, -2.93406929073157E-04_wp, &
      & -4.22643021920042E-05_wp,  5.61998503299879E-05_wp,  2.53053704317884E-04_wp, &
      & -5.24441732338870E-04_wp, -4.95035735966572E-05_wp,  5.78446720207333E-04_wp, &
      &  3.98609879103522E-04_wp,  3.79727615859641E-05_wp,  9.37490197627848E-05_wp, &
      & -8.82402366526675E-05_wp, -1.30739972630176E-02_wp, -1.45598047446246E-04_wp, &
      &  1.23534957245937E-04_wp,  3.17929622266108E-04_wp,  6.79885084012242E-05_wp, &
      &  1.47192601575797E-03_wp,  4.23605560773337E-03_wp,  2.74971980733658E-03_wp, &
      & -1.20802612207777E-03_wp, -3.16675114990010E-03_wp, -6.87492098453706E-04_wp, &
      &  8.80222799798718E-04_wp,  8.68076184347637E-04_wp,  2.00393724987738E-04_wp, &
      &  7.55965934758805E-05_wp, -5.46213445312269E-04_wp, -2.36311205706913E-04_wp, &
      &  1.95318300369122E-05_wp,  2.54765110784266E-04_wp,  2.61570234221497E-05_wp, &
      &  1.88854498309724E-04_wp,  1.25565553689502E-04_wp,  1.22717877592350E-04_wp, &
      &  1.52180609814273E-04_wp,  7.66217690306814E-05_wp,  9.43211898934151E-05_wp, &
      &  6.34735644322809E-05_wp, -1.67830034552342E-05_wp,  8.34663571047527E-05_wp, &
      &  2.81883813580586E-04_wp,  1.59750178698940E-04_wp, -8.15014124206673E-05_wp, &
      & -1.73309086839298E-04_wp, -1.79297232198718E-04_wp,  2.71098448232089E-05_wp, &
      &  7.49922749173758E-05_wp, -5.69553870423969E-04_wp,  8.07806855518498E-05_wp, &
      &  8.60225386342648E-05_wp,  7.00549995343506E-04_wp, -2.65603191426876E-06_wp, &
      & -4.16134092031767E-05_wp,  1.85758046632624E-04_wp, -1.45598047446246E-04_wp, &
      & -1.30350273737420E-02_wp], shape(acp))

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "CeCl3")
   call test_acp_gen(error, mol, make_qvszp_basis, acp, thr_in=thr1)

end subroutine test_acp_gxtb_cecl3

subroutine test_acp_gxtb_ce2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 32
   real(wp), parameter :: acp(nao, nao) = reshape([&
      & -1.94261827957380E-04_wp,  0.00000000000000E+00_wp, -2.53251893234292E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.92288807951860E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -6.83088912817165E-05_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.27487967800631E-04_wp,  0.00000000000000E+00_wp, &
      & -1.16821235485607E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.42593250313637E-04_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.98584769957678E-05_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.74637222880388E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -7.37750306978122E-04_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  3.30777206003856E-20_wp, &
      &  0.00000000000000E+00_wp, -2.50419191345714E-04_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  3.01320354010766E-05_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.16600078702005E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.12808025549363E-20_wp,  0.00000000000000E+00_wp,  9.63351866794148E-06_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.53251893234292E-04_wp,  0.00000000000000E+00_wp, &
      & -7.37358828314758E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.31701253795925E-04_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.81071038494121E-05_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.16821235485607E-04_wp, &
      &  0.00000000000000E+00_wp, -9.64136254433220E-06_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.61595620378633E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.05198281984097E-05_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.74637222880388E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -7.37750306978122E-04_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.50419191345714E-04_wp,  0.00000000000000E+00_wp, &
      & -3.30777206003856E-20_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  3.01320354010766E-05_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.16600078702005E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  9.63351866794148E-06_wp, &
      &  0.00000000000000E+00_wp, -2.12808025549363E-20_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.99616477960280E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.94080796222895E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.84235364597170E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.60995406446780E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -7.37750306978122E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.18538891440886E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  4.86777811925436E-22_wp,  0.00000000000000E+00_wp, &
      & -7.44132128902110E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.16600078702005E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.30200657249153E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.11766827385847E-19_wp, &
      &  0.00000000000000E+00_wp, -2.32097767078994E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.92288807951860E-04_wp,  0.00000000000000E+00_wp,  1.31701253795925E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.04416794603269E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.59099483202302E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.42593250313637E-04_wp,  0.00000000000000E+00_wp, &
      & -2.61595620378633E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -6.30013419318605E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.08792272561842E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -7.37750306978122E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.18538891440886E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -7.44132128902110E-04_wp,  0.00000000000000E+00_wp, -4.86777811925436E-22_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.16600078702005E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.30200657249153E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.32097767078994E-03_wp,  0.00000000000000E+00_wp, &
      & -2.11766827385847E-19_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.99616477960280E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.15809925790646E-20_wp,  0.00000000000000E+00_wp, &
      & -1.94080796222895E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.84235364597170E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -9.40391683609779E-20_wp, &
      &  0.00000000000000E+00_wp,  1.60995406446779E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  3.30777206003856E-20_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  4.86777811925436E-22_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -8.96315397540680E-04_wp,  0.00000000000000E+00_wp,  6.01970189327411E-20_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.12808025549363E-20_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  5.88340510892458E-21_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.74893445341902E-05_wp,  0.00000000000000E+00_wp, &
      &  1.87633036946880E-20_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.94080796222895E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.02125320281756E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.60995406446780E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  8.07774603188906E-05_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.50419191345714E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -7.44132128902110E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  6.01970189327411E-20_wp,  0.00000000000000E+00_wp, &
      & -1.14993954806604E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  9.63351866794171E-06_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.32097767078994E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.90975742235966E-20_wp, &
      &  0.00000000000000E+00_wp, -1.40713283317480E-04_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -6.83088912817165E-05_wp,  0.00000000000000E+00_wp,  2.81071038494121E-05_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.59099483202302E-04_wp,  0.00000000000000E+00_wp,  1.15809925790646E-20_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -9.84158896808183E-04_wp,  0.00000000000000E+00_wp,  7.44431390585514E-21_wp, &
      &  0.00000000000000E+00_wp,  4.98584769957678E-05_wp,  0.00000000000000E+00_wp, &
      & -5.05198281984101E-05_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.08792272561842E-03_wp,  0.00000000000000E+00_wp, &
      &  9.40391683609779E-20_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.48602978957833E-04_wp,  0.00000000000000E+00_wp, &
      &  2.80772590782603E-36_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.50419191345714E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -7.44132128902110E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.14993954806604E-03_wp,  0.00000000000000E+00_wp, -6.01970189327411E-20_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  9.63351866794171E-06_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.32097767078994E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.40713283317480E-04_wp,  0.00000000000000E+00_wp, &
      & -2.90975742235966E-20_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.94080796222895E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  7.44431390585514E-21_wp,  0.00000000000000E+00_wp, &
      & -1.02125320281756E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.60995406446779E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.44732942036055E-36_wp, &
      &  0.00000000000000E+00_wp,  8.07774603188905E-05_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.30777206003856E-20_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.86777811925436E-22_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -6.01970189327411E-20_wp,  0.00000000000000E+00_wp, &
      & -8.96315397540680E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.12808025549363E-20_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.88340510892458E-21_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.87633036946880E-20_wp, &
      &  0.00000000000000E+00_wp, -1.74893445341902E-05_wp,  1.27487967800631E-04_wp, &
      &  0.00000000000000E+00_wp,  1.16821235485607E-04_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.42593250313637E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  4.98584769957678E-05_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.94261827957380E-04_wp,  0.00000000000000E+00_wp,  2.53251893234292E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.92288807951860E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  6.83088912817165E-05_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  3.01320354010766E-05_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.16600078702005E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.12808025549363E-20_wp,  0.00000000000000E+00_wp, &
      &  9.63351866794171E-06_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.74637222880388E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  7.37750306978118E-04_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  5.39834532408105E-20_wp, &
      &  0.00000000000000E+00_wp, -2.50419191345714E-04_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.16821235485607E-04_wp,  0.00000000000000E+00_wp, -9.64136254433220E-06_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.61595620378633E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.05198281984101E-05_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.53251893234292E-04_wp,  0.00000000000000E+00_wp, &
      & -7.37358828314758E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.31701253795929E-04_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.81071038494123E-05_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  3.01320354010766E-05_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.16600078702005E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  9.63351866794171E-06_wp,  0.00000000000000E+00_wp, -2.12808025549363E-20_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.74637222880388E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  7.37750306978118E-04_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.50419191345714E-04_wp,  0.00000000000000E+00_wp, &
      & -5.39834532408105E-20_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -4.84235364597170E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.60995406446780E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.99616477960280E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.94080796222894E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.16600078702005E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.30200657249153E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  5.88340510892458E-21_wp,  0.00000000000000E+00_wp,  2.32097767078994E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  7.37750306978118E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.18538891440886E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -6.21977190069225E-20_wp,  0.00000000000000E+00_wp, &
      &  7.44132128902108E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.42593250313637E-04_wp, &
      &  0.00000000000000E+00_wp,  2.61595620378633E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -6.30013419318605E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.08792272561842E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.92288807951860E-04_wp,  0.00000000000000E+00_wp, -1.31701253795929E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.04416794603269E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.59099483202300E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.16600078702005E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.30200657249153E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.32097767078994E-03_wp, &
      &  0.00000000000000E+00_wp, -5.88340510892458E-21_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  7.37750306978118E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.18538891440886E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  7.44132128902108E-04_wp,  0.00000000000000E+00_wp,  6.21977190069225E-20_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -4.84235364597170E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  9.40391683609779E-20_wp,  0.00000000000000E+00_wp, -1.60995406446779E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.99616477960280E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.15809925790646E-20_wp,  0.00000000000000E+00_wp, &
      &  1.94080796222894E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.12808025549363E-20_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.11766827385847E-19_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.74893445341902E-05_wp, &
      &  0.00000000000000E+00_wp,  2.90975742235966E-20_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  5.39834532408105E-20_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -6.21977190069225E-20_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -8.96315397540680E-04_wp,  0.00000000000000E+00_wp,  8.08579155603246E-20_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.60995406446780E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  8.07774603188906E-05_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.94080796222894E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.02125320281756E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  9.63351866794148E-06_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.32097767078994E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.87633036946880E-20_wp,  0.00000000000000E+00_wp, -1.40713283317480E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.50419191345714E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  7.44132128902108E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  8.08579155603246E-20_wp,  0.00000000000000E+00_wp, &
      & -1.14993954806604E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -4.98584769957678E-05_wp, &
      &  0.00000000000000E+00_wp, -5.05198281984097E-05_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.08792272561842E-03_wp, &
      &  0.00000000000000E+00_wp, -9.40391683609779E-20_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.48602978957833E-04_wp, &
      &  0.00000000000000E+00_wp, -2.44732942036055E-36_wp,  0.00000000000000E+00_wp, &
      &  6.83088912817165E-05_wp,  0.00000000000000E+00_wp,  2.81071038494123E-05_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.59099483202300E-04_wp,  0.00000000000000E+00_wp, -1.15809925790646E-20_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -9.84158896808183E-04_wp,  0.00000000000000E+00_wp,  7.44431390585514E-21_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  9.63351866794148E-06_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.32097767078994E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.40713283317480E-04_wp, &
      &  0.00000000000000E+00_wp, -1.87633036946880E-20_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.50419191345714E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  7.44132128902108E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.14993954806604E-03_wp,  0.00000000000000E+00_wp, -8.08579155603246E-20_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.60995406446779E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.80772590782603E-36_wp,  0.00000000000000E+00_wp,  8.07774603188905E-05_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.94080796222894E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  7.44431390585514E-21_wp,  0.00000000000000E+00_wp, &
      & -1.02125320281756E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.12808025549363E-20_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.11766827385847E-19_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.90975742235966E-20_wp,  0.00000000000000E+00_wp, -1.74893445341902E-05_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.39834532408105E-20_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  6.21977190069225E-20_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -8.08579155603246E-20_wp,  0.00000000000000E+00_wp, &
      & -8.96315397540680E-04_wp], shape(acp))

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "Ce2")
   call test_acp_gen(error, mol, make_qvszp_basis, acp, thr_in=thr1)

end subroutine test_acp_gxtb_ce2


subroutine test_qeff_grad_acp_gxtb_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(2, 2, 1) = reshape([&
      & 5.93683766916992E-01_wp, 5.93683766916992E-01_wp, 5.93683766916992E-01_wp, &
      & 5.93683766916992E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "H2")
   call test_acp_qeff_numgrad(error, mol, density, make_qvszp_basis, thr_in=thr1*10)

end subroutine test_qeff_grad_acp_gxtb_h2

subroutine test_qeff_grad_acp_gxtb_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(5, 5, 1) = reshape([&
      & 7.43138968868807E-02_wp, 6.30732585440217E-45_wp, 1.15038033099932E-01_wp, &
      & 0.00000000000000E+00_wp, 2.77067464359637E-01_wp, 6.30732585440217E-45_wp, &
      & 5.35328667990135E-88_wp, 9.76375066853567E-45_wp, 0.00000000000000E+00_wp, &
      & 2.35158544306899E-44_wp, 1.15038033099932E-01_wp, 9.76375066853567E-45_wp, &
      & 1.78079062111965E-01_wp, 0.00000000000000E+00_wp, 4.28900884910329E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 2.77067464359637E-01_wp, &
      & 2.35158544306899E-44_wp, 4.28900884910329E-01_wp, 0.00000000000000E+00_wp, &
      & 1.03300167293785E+00_wp], shape(density))

   call get_structure(mol, "MB16-43", "LiH")
   call test_acp_qeff_numgrad(error, mol, density, make_qvszp_basis, thr_in=thr1)

end subroutine test_qeff_grad_acp_gxtb_lih

subroutine test_qeff_grad_acp_gxtb_no(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(8, 8, 2) = reshape([&
      & 9.42009046415725E-01_wp, 8.68128306151435E-17_wp,-3.13530611390051E-01_wp, &
      & 2.83821071259847E-17_wp,-1.91210322106378E-01_wp,-2.95777592933361E-16_wp, &
      &-2.58330200255829E-02_wp,-1.88965537110059E-16_wp, 8.68128306151435E-17_wp, &
      & 7.55421320827554E-01_wp, 4.13035452754259E-17_wp,-4.05630015646006E-01_wp, &
      &-7.56771273088709E-17_wp,-3.58191420308672E-02_wp, 4.29283702464554E-17_wp, &
      & 3.33220482034872E-01_wp,-3.13530611390051E-01_wp, 4.13035452754259E-17_wp, &
      & 4.49934652749067E-01_wp,-1.37354302670472E-16_wp,-1.71924644995845E-03_wp, &
      &-8.68821483997484E-17_wp,-4.07076133805002E-01_wp, 1.81509948244029E-16_wp, &
      & 2.83821071259847E-17_wp,-4.05630015646006E-01_wp,-1.37354302670472E-16_wp, &
      & 5.92359279622604E-01_wp, 1.27369180291386E-16_wp, 3.33220482034872E-01_wp, &
      & 1.10259833116572E-17_wp, 9.81344853538432E-02_wp,-1.91210322106378E-01_wp, &
      &-7.56771273088709E-17_wp,-1.71924644995845E-03_wp, 1.27369180291386E-16_wp, &
      & 1.01808686826282E+00_wp, 1.30483506117222E-16_wp, 2.62723878284128E-01_wp, &
      &-2.46091681044612E-16_wp,-2.95777592933361E-16_wp,-3.58191420308672E-02_wp, &
      &-8.68821483997484E-17_wp, 3.33220482034872E-01_wp, 1.30483506117222E-16_wp, &
      & 8.63442034138714E-01_wp,-5.57245658663308E-17_wp,-2.73736867008512E-01_wp, &
      &-2.58330200255829E-02_wp, 4.29283702464554E-17_wp,-4.07076133805002E-01_wp, &
      & 1.10259833116572E-17_wp, 2.62723878284128E-01_wp,-5.57245658663308E-17_wp, &
      & 5.33778495173348E-01_wp,-1.58550889927256E-17_wp,-1.88965537110059E-16_wp, &
      & 3.33220482034872E-01_wp, 1.81509948244029E-16_wp, 9.81344853538432E-02_wp, &
      &-2.46091681044612E-16_wp,-2.73736867008512E-01_wp,-1.58550889927256E-17_wp, &
      & 7.53400640342027E-01_wp, 9.54305653170233E-01_wp, 1.09464057770057E-17_wp, &
      &-3.11307650277096E-01_wp,-7.36206352015003E-17_wp,-2.03773626333043E-01_wp, &
      & 1.65290626496932E-16_wp,-1.70067070429939E-02_wp,-8.93509590142745E-17_wp, &
      & 1.09464057770057E-17_wp, 2.40282587827900E-01_wp, 4.33933558995002E-17_wp, &
      & 1.85040392733642E-05_wp, 6.74251196637872E-17_wp, 3.64459920248371E-01_wp, &
      & 9.51214979595760E-17_wp, 6.68285106103067E-06_wp,-3.11307650277096E-01_wp, &
      & 4.33933558995002E-17_wp, 4.38731609559950E-01_wp, 2.61613977605231E-17_wp, &
      &-4.10640194707355E-03_wp, 1.10151806812206E-16_wp,-4.09548671902720E-01_wp, &
      & 2.83066560924208E-19_wp,-7.36206352015003E-17_wp, 1.85040392733642E-05_wp, &
      & 2.61613977605231E-17_wp, 2.40290026395803E-01_wp, 3.09134882803496E-17_wp, &
      & 6.68285106153466E-06_wp,-3.66627068485342E-18_wp, 3.64462606734314E-01_wp, &
      &-2.03773626333043E-01_wp, 6.74251196637872E-17_wp,-4.10640194707355E-03_wp, &
      & 3.09134882803496E-17_wp, 1.03092149514987E+00_wp,-3.17330346614473E-16_wp, &
      & 2.53665476333746E-01_wp,-1.60068098192912E-17_wp, 1.65290626496932E-16_wp, &
      & 3.64459920248371E-01_wp, 1.10151806812206E-16_wp, 6.68285106153466E-06_wp, &
      &-3.17330346614473E-16_wp, 5.52811733573677E-01_wp,-3.44555050846798E-17_wp, &
      &-2.22979264938861E-05_wp,-1.70067070429939E-02_wp, 9.51214979595760E-17_wp, &
      &-4.09548671902720E-01_wp,-3.66627068485342E-18_wp, 2.53665476333746E-01_wp, &
      &-3.44555050846798E-17_wp, 5.38687781287961E-01_wp, 2.26914599114829E-17_wp, &
      &-8.93509590142745E-17_wp, 6.68285106103067E-06_wp, 2.83066560924208E-19_wp, &
      & 3.64462606734314E-01_wp,-1.60068098192912E-17_wp,-2.22979264938861E-05_wp, &
      & 2.26914599114829E-17_wp, 5.52802769874570E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "NO")
   call test_acp_qeff_numgrad(error, mol, density, make_qvszp_basis, thr_in=thr1*10)

end subroutine test_qeff_grad_acp_gxtb_no

subroutine test_qeff_grad_acp_gxtb_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(18, 18, 1) = reshape([&
      & 2.01512274410725E+00_wp, 3.05183510489231E-16_wp,-3.99929219814921E-01_wp, &
      & 7.69080296474742E-16_wp, 1.45947463264597E-04_wp,-4.42069429901860E-17_wp, &
      &-3.51352512626506E-02_wp,-1.08853924469132E-16_wp, 8.19811912191895E-05_wp, &
      &-2.61461503009236E-01_wp,-2.53847513922895E-16_wp, 8.46120257485079E-02_wp, &
      &-1.05417595234967E-15_wp,-2.59372907294474E-05_wp,-4.35889228793221E-17_wp, &
      & 1.18677041085995E-02_wp,-9.23447662047026E-17_wp,-1.45694206907493E-05_wp, &
      & 3.05183510489231E-16_wp, 1.11646346539116E+00_wp, 1.51923309365924E-15_wp, &
      & 5.22665236979214E-01_wp,-2.23623283632863E-17_wp, 2.23688961934518E-02_wp, &
      &-3.01660087631464E-16_wp,-3.85852024970845E-02_wp,-2.43149330627879E-18_wp, &
      &-6.95644355001582E-16_wp, 5.02302414177719E-01_wp, 5.74361021043659E-16_wp, &
      &-5.26740092615372E-01_wp, 7.76809871729795E-17_wp,-7.25879700203994E-02_wp, &
      & 3.46931693965287E-16_wp,-4.72231756427874E-02_wp, 9.75451031609286E-17_wp, &
      &-3.99929219814921E-01_wp, 1.51923309365924E-15_wp, 8.16076844347648E-01_wp, &
      & 1.34007033367206E-15_wp,-2.39395566884248E-05_wp,-7.79787012198157E-17_wp, &
      & 5.39351481157056E-02_wp, 1.48581899014128E-18_wp,-1.34472592440345E-05_wp, &
      &-8.46120257485160E-02_wp,-1.44825907786549E-15_wp,-7.72403989719101E-01_wp, &
      &-1.76921592471880E-15_wp,-1.32742422451715E-07_wp,-1.15813374568686E-16_wp, &
      & 4.74250266390859E-02_wp,-1.82310535840248E-17_wp,-7.45636935369285E-08_wp, &
      & 7.69080296474742E-16_wp, 5.22665236979214E-01_wp, 1.34007033367206E-15_wp, &
      & 1.70364350645588E+00_wp,-8.02983753144685E-18_wp,-3.85852024970846E-02_wp, &
      &-4.70255064846748E-16_wp,-2.09790427620324E-02_wp,-1.96205026639568E-17_wp, &
      &-6.57821241716489E-16_wp,-5.26740092615371E-01_wp, 1.54875272040670E-15_wp, &
      &-8.94554593703680E-02_wp,-1.47010689589666E-16_wp,-4.72231756427873E-02_wp, &
      & 5.33348516140237E-16_wp,-1.25640104065888E-01_wp, 1.48730635567615E-16_wp, &
      & 1.45947463264597E-04_wp,-2.23623283632863E-17_wp,-2.39395566884248E-05_wp, &
      &-8.02983753144685E-18_wp, 1.06235239203714E-08_wp,-6.27312114532195E-19_wp, &
      &-2.27401943747449E-06_wp, 6.00156399076560E-19_wp, 5.96741544153476E-09_wp, &
      &-2.59372907296047E-05_wp,-1.26125559712104E-17_wp, 1.32742422317945E-07_wp, &
      & 8.58050673700442E-18_wp,-2.35397458118598E-09_wp, 1.42171549599656E-18_wp, &
      & 1.27475826808189E-06_wp, 7.55995255055379E-19_wp,-1.32226786236991E-09_wp, &
      &-4.42069429901860E-17_wp, 2.23688961934518E-02_wp,-7.79787012198157E-17_wp, &
      &-3.85852024970846E-02_wp,-6.27312114532195E-19_wp, 4.86266676876034E-03_wp, &
      & 3.34959389043418E-17_wp, 4.04941245032744E-03_wp, 6.35609235074649E-19_wp, &
      & 6.31298797212952E-17_wp, 7.25879700204016E-02_wp,-7.29477821938508E-17_wp, &
      & 4.72231756427944E-02_wp, 1.03467976744611E-17_wp,-7.56324676480364E-04_wp, &
      &-1.77577932859663E-17_wp, 2.96701633209598E-03_wp,-6.46245268795547E-18_wp, &
      &-3.51352512626506E-02_wp,-3.01660087631464E-16_wp, 5.39351481157056E-02_wp, &
      &-4.70255064846748E-16_wp,-2.27401943747449E-06_wp, 3.34959389043418E-17_wp, &
      & 3.73731424127451E-03_wp, 7.17979015797382E-17_wp,-1.27735568790833E-06_wp, &
      & 1.18677041086054E-02_wp, 3.72567464311869E-16_wp,-4.74250266390819E-02_wp, &
      & 8.06318001452668E-16_wp, 1.27475826807370E-06_wp, 2.50819963630132E-17_wp, &
      & 2.76687889875194E-03_wp, 4.75201425698955E-17_wp, 7.16053564725540E-07_wp, &
      &-1.08853924469132E-16_wp,-3.85852024970845E-02_wp, 1.48581899014128E-18_wp, &
      &-2.09790427620324E-02_wp, 6.00156399076560E-19_wp, 4.04941245032744E-03_wp, &
      & 7.17979015797382E-17_wp, 9.41191550512845E-03_wp, 2.27476146379299E-19_wp, &
      &-3.86583716090084E-17_wp, 4.72231756427946E-02_wp,-3.15405023331482E-16_wp, &
      & 1.25640104065898E-01_wp, 2.15478167900987E-18_wp, 2.96701633209597E-03_wp, &
      &-1.89192522936108E-17_wp, 2.57692316259317E-03_wp,-1.20436763189445E-17_wp, &
      & 8.19811912191895E-05_wp,-2.43149330627879E-18_wp,-1.34472592440345E-05_wp, &
      &-1.96205026639568E-17_wp, 5.96741544153476E-09_wp, 6.35609235074649E-19_wp, &
      &-1.27735568790833E-06_wp, 2.27476146379299E-19_wp, 3.35199951718306E-09_wp, &
      &-1.45694206906697E-05_wp, 9.40823506771232E-18_wp, 7.45636930694748E-08_wp, &
      & 5.82344702349143E-19_wp,-1.32226786235450E-09_wp, 3.30923532167611E-19_wp, &
      & 7.16053564753441E-07_wp, 1.42332895146119E-18_wp,-7.42740518019966E-10_wp, &
      &-2.61461503009236E-01_wp,-6.95644355001582E-16_wp,-8.46120257485160E-02_wp, &
      &-6.57821241716489E-16_wp,-2.59372907296047E-05_wp, 6.31298797212952E-17_wp, &
      & 1.18677041086054E-02_wp,-3.86583716090084E-17_wp,-1.45694206906697E-05_wp, &
      & 2.01512274410725E+00_wp, 6.37025174552147E-16_wp, 3.99929219814931E-01_wp, &
      & 8.33671871420983E-16_wp, 1.45947463264675E-04_wp, 5.98740820853847E-17_wp, &
      &-3.51352512626559E-02_wp, 6.97851591164381E-17_wp, 8.19811912194447E-05_wp, &
      &-2.53847513922895E-16_wp, 5.02302414177719E-01_wp,-1.44825907786549E-15_wp, &
      &-5.26740092615371E-01_wp,-1.26125559712104E-17_wp, 7.25879700204016E-02_wp, &
      & 3.72567464311869E-16_wp, 4.72231756427946E-02_wp, 9.40823506771232E-18_wp, &
      & 6.37025174552147E-16_wp, 1.11646346539112E+00_wp,-3.83819455107282E-16_wp, &
      & 5.22665236979232E-01_wp, 1.65728863182037E-16_wp,-2.23688961934538E-02_wp, &
      &-2.35707388589142E-16_wp, 3.85852024970803E-02_wp,-7.60309218403261E-17_wp, &
      & 8.46120257485079E-02_wp, 5.74361021043659E-16_wp,-7.72403989719101E-01_wp, &
      & 1.54875272040670E-15_wp, 1.32742422317945E-07_wp,-7.29477821938508E-17_wp, &
      &-4.74250266390819E-02_wp,-3.15405023331482E-16_wp, 7.45636930694748E-08_wp, &
      & 3.99929219814931E-01_wp,-3.83819455107282E-16_wp, 8.16076844347611E-01_wp, &
      &-2.03826552101056E-15_wp, 2.39395566885481E-05_wp,-4.64841795882909E-17_wp, &
      &-5.39351481157069E-02_wp,-2.20759476011753E-16_wp, 1.34472592445486E-05_wp, &
      &-1.05417595234967E-15_wp,-5.26740092615372E-01_wp,-1.76921592471880E-15_wp, &
      &-8.94554593703680E-02_wp, 8.58050673700442E-18_wp, 4.72231756427944E-02_wp, &
      & 8.06318001452668E-16_wp, 1.25640104065898E-01_wp, 5.82344702349143E-19_wp, &
      & 8.33671871420983E-16_wp, 5.22665236979232E-01_wp,-2.03826552101056E-15_wp, &
      & 1.70364350645586E+00_wp, 3.39236205107499E-18_wp, 3.85852024970803E-02_wp, &
      &-3.39364583169101E-16_wp, 2.09790427620258E-02_wp,-1.47879398954886E-16_wp, &
      &-2.59372907294474E-05_wp, 7.76809871729795E-17_wp,-1.32742422451715E-07_wp, &
      &-1.47010689589666E-16_wp,-2.35397458118598E-09_wp, 1.03467976744611E-17_wp, &
      & 1.27475826807370E-06_wp, 2.15478167900987E-18_wp,-1.32226786235450E-09_wp, &
      & 1.45947463264675E-04_wp, 1.65728863182037E-16_wp, 2.39395566885481E-05_wp, &
      & 3.39236205107499E-18_wp, 1.06235239203804E-08_wp,-3.14583296736217E-18_wp, &
      &-2.27401943748182E-06_wp, 1.01384858020498E-17_wp, 5.96741544155305E-09_wp, &
      &-4.35889228793221E-17_wp,-7.25879700203994E-02_wp,-1.15813374568686E-16_wp, &
      &-4.72231756427873E-02_wp, 1.42171549599656E-18_wp,-7.56324676480364E-04_wp, &
      & 2.50819963630132E-17_wp, 2.96701633209597E-03_wp, 3.30923532167611E-19_wp, &
      & 5.98740820853847E-17_wp,-2.23688961934538E-02_wp,-4.64841795882909E-17_wp, &
      & 3.85852024970803E-02_wp,-3.14583296736217E-18_wp, 4.86266676875980E-03_wp, &
      &-2.75950952491621E-17_wp, 4.04941245032639E-03_wp,-7.72981644699408E-18_wp, &
      & 1.18677041085995E-02_wp, 3.46931693965287E-16_wp, 4.74250266390859E-02_wp, &
      & 5.33348516140237E-16_wp, 1.27475826808189E-06_wp,-1.77577932859663E-17_wp, &
      & 2.76687889875194E-03_wp,-1.89192522936108E-17_wp, 7.16053564753441E-07_wp, &
      &-3.51352512626559E-02_wp,-2.35707388589142E-16_wp,-5.39351481157069E-02_wp, &
      &-3.39364583169101E-16_wp,-2.27401943748182E-06_wp,-2.75950952491621E-17_wp, &
      & 3.73731424127476E-03_wp,-3.67074117602398E-17_wp,-1.27735568794358E-06_wp, &
      &-9.23447662047026E-17_wp,-4.72231756427874E-02_wp,-1.82310535840248E-17_wp, &
      &-1.25640104065888E-01_wp, 7.55995255055379E-19_wp, 2.96701633209598E-03_wp, &
      & 4.75201425698955E-17_wp, 2.57692316259317E-03_wp, 1.42332895146119E-18_wp, &
      & 6.97851591164381E-17_wp, 3.85852024970803E-02_wp,-2.20759476011753E-16_wp, &
      & 2.09790427620258E-02_wp, 1.01384858020498E-17_wp, 4.04941245032639E-03_wp, &
      &-3.67074117602398E-17_wp, 9.41191550512677E-03_wp,-1.22124983920611E-17_wp, &
      &-1.45694206907493E-05_wp, 9.75451031609286E-17_wp,-7.45636935369285E-08_wp, &
      & 1.48730635567615E-16_wp,-1.32226786236991E-09_wp,-6.46245268795547E-18_wp, &
      & 7.16053564725540E-07_wp,-1.20436763189445E-17_wp,-7.42740518019966E-10_wp, &
      & 8.19811912194447E-05_wp,-7.60309218403261E-17_wp, 1.34472592445486E-05_wp, &
      &-1.47879398954886E-16_wp, 5.96741544155305E-09_wp,-7.72981644699408E-18_wp, &
      &-1.27735568794358E-06_wp,-1.22124983920611E-17_wp, 3.35199951720075E-09_wp],&
      & shape(density))

   call get_structure(mol, "MB16-43", "S2")
   call test_acp_qeff_numgrad(error, mol, density, make_qvszp_basis, thr_in=thr1*10)

end subroutine test_qeff_grad_acp_gxtb_s2

subroutine test_qeff_grad_acp_gxtb_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(13, 13, 1) = reshape([&
      & 6.67076596984081E-01_wp, 1.17812091373782E-15_wp, 1.16003518499611E-15_wp, &
      &-9.41177582170035E-16_wp,-5.18269183303055E-17_wp, 2.81268188491212E-17_wp, &
      & 1.09658020136346E-17_wp,-2.39914571721069E-16_wp, 2.14126324730246E-19_wp, &
      & 2.51982366144073E-01_wp, 2.51982366144073E-01_wp, 2.51982366144073E-01_wp, &
      & 2.51982366144077E-01_wp, 1.17812091373782E-15_wp, 4.65632832545944E-01_wp, &
      &-3.08027819846221E-16_wp,-5.01335488096171E-17_wp, 5.51922671972401E-16_wp, &
      &-4.45360639961878E-16_wp,-7.86723569094682E-18_wp,-5.58630716839440E-02_wp, &
      &-7.00845070605205E-17_wp, 3.05808254914735E-01_wp,-3.05808254914736E-01_wp, &
      &-3.05808254914736E-01_wp, 3.05808254914735E-01_wp, 1.16003518499611E-15_wp, &
      &-3.08027819846221E-16_wp, 4.65632832545945E-01_wp,-3.98576615792397E-17_wp, &
      &-5.58630716839442E-02_wp,-3.56963383436534E-16_wp, 9.89192849362209E-17_wp, &
      & 3.11962003002779E-16_wp, 1.36197924766342E-16_wp,-3.05808254914736E-01_wp, &
      &-3.05808254914736E-01_wp, 3.05808254914735E-01_wp, 3.05808254914735E-01_wp, &
      &-9.41177582170035E-16_wp,-5.01335488096171E-17_wp,-3.98576615792397E-17_wp, &
      & 4.65632832545944E-01_wp,-4.07263775340787E-16_wp,-5.58630716839439E-02_wp, &
      & 5.44186509945490E-17_wp,-5.04840605943758E-16_wp, 1.59146606601200E-17_wp, &
      & 3.05808254914736E-01_wp,-3.05808254914735E-01_wp, 3.05808254914736E-01_wp, &
      &-3.05808254914735E-01_wp,-5.18269183303055E-17_wp, 5.51922671972401E-16_wp, &
      &-5.58630716839442E-02_wp,-4.07263775340787E-16_wp, 6.70202477111097E-03_wp, &
      & 9.22921619598948E-17_wp,-1.18675804605582E-17_wp,-9.93121625154311E-17_wp, &
      &-1.63399870082742E-17_wp, 3.66885392776897E-02_wp, 3.66885392776896E-02_wp, &
      &-3.66885392776901E-02_wp,-3.66885392776889E-02_wp, 2.81268188491212E-17_wp, &
      &-4.45360639961878E-16_wp,-3.56963383436534E-16_wp,-5.58630716839439E-02_wp, &
      & 9.22921619598948E-17_wp, 6.70202477111090E-03_wp,-6.52873420637083E-18_wp, &
      & 1.14654961676060E-16_wp,-1.90931946190504E-18_wp,-3.66885392776896E-02_wp, &
      & 3.66885392776898E-02_wp,-3.66885392776895E-02_wp, 3.66885392776888E-02_wp, &
      & 1.09658020136346E-17_wp,-7.86723569094682E-18_wp, 9.89192849362209E-17_wp, &
      & 5.44186509945490E-17_wp,-1.18675804605582E-17_wp,-6.52873420637083E-18_wp, &
      & 2.76875750299656E-32_wp, 9.43850864112936E-19_wp, 3.19815635497338E-32_wp, &
      &-3.02508000498458E-17_wp,-9.13968622807606E-17_wp, 1.10015080012159E-16_wp, &
      & 2.82015265779562E-17_wp,-2.39914571721069E-16_wp,-5.58630716839440E-02_wp, &
      & 3.11962003002779E-16_wp,-5.04840605943758E-16_wp,-9.93121625154311E-17_wp, &
      & 1.14654961676060E-16_wp, 9.43850864112936E-19_wp, 6.70202477111093E-03_wp, &
      & 8.40820399293783E-18_wp,-3.66885392776899E-02_wp, 3.66885392776897E-02_wp, &
      & 3.66885392776894E-02_wp,-3.66885392776889E-02_wp, 2.14126324730246E-19_wp, &
      &-7.00845070605205E-17_wp, 1.36197924766342E-16_wp, 1.59146606601200E-17_wp, &
      &-1.63399870082742E-17_wp,-1.90931946190504E-18_wp, 3.19815635497338E-32_wp, &
      & 8.40820399293783E-18_wp, 5.09307325669261E-32_wp,-1.24944740611828E-16_wp, &
      &-5.37917417731081E-17_wp, 1.46010682147439E-16_wp, 3.30493376435646E-17_wp, &
      & 2.51982366144073E-01_wp, 3.05808254914735E-01_wp,-3.05808254914736E-01_wp, &
      & 3.05808254914736E-01_wp, 3.66885392776897E-02_wp,-3.66885392776896E-02_wp, &
      &-3.02508000498458E-17_wp,-3.66885392776899E-02_wp,-1.24944740611828E-16_wp, &
      & 6.97710523806227E-01_wp,-1.05657986637187E-01_wp,-1.05657986637186E-01_wp, &
      &-1.05657986637189E-01_wp, 2.51982366144073E-01_wp,-3.05808254914736E-01_wp, &
      &-3.05808254914736E-01_wp,-3.05808254914735E-01_wp, 3.66885392776896E-02_wp, &
      & 3.66885392776898E-02_wp,-9.13968622807606E-17_wp, 3.66885392776897E-02_wp, &
      &-5.37917417731081E-17_wp,-1.05657986637187E-01_wp, 6.97710523806228E-01_wp, &
      &-1.05657986637187E-01_wp,-1.05657986637188E-01_wp, 2.51982366144073E-01_wp, &
      &-3.05808254914736E-01_wp, 3.05808254914735E-01_wp, 3.05808254914736E-01_wp, &
      &-3.66885392776901E-02_wp,-3.66885392776895E-02_wp, 1.10015080012159E-16_wp, &
      & 3.66885392776894E-02_wp, 1.46010682147439E-16_wp,-1.05657986637186E-01_wp, &
      &-1.05657986637187E-01_wp, 6.97710523806227E-01_wp,-1.05657986637189E-01_wp, &
      & 2.51982366144077E-01_wp, 3.05808254914735E-01_wp, 3.05808254914735E-01_wp, &
      &-3.05808254914735E-01_wp,-3.66885392776889E-02_wp, 3.66885392776888E-02_wp, &
      & 2.82015265779562E-17_wp,-3.66885392776889E-02_wp, 3.30493376435646E-17_wp, &
      &-1.05657986637189E-01_wp,-1.05657986637188E-01_wp,-1.05657986637189E-01_wp, &
      & 6.97710523806226E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "SiH4")
   call test_acp_qeff_numgrad(error, mol, density, make_qvszp_basis, thr_in=thr1*10)

end subroutine test_qeff_grad_acp_gxtb_sih4

subroutine test_qeff_grad_acp_gxtb_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(43, 43, 2) = reshape([&
      & 6.91475853693640E-03_wp,-4.66439474613685E-04_wp, 7.43124629353357E-04_wp, &
      &-3.09904488027013E-04_wp,-2.72087829478411E-03_wp, 6.74488694981665E-03_wp, &
      &-5.28266720158637E-03_wp, 4.80689218346076E-03_wp, 9.22601160310286E-04_wp, &
      &-2.71966921742355E-03_wp,-1.01894320855533E-02_wp,-2.27063129413855E-03_wp, &
      &-2.82181425185299E-03_wp,-2.14613395036928E-03_wp, 3.61217559584011E-03_wp, &
      & 1.12028617141975E-02_wp,-8.64033744991843E-04_wp, 3.10371374558179E-02_wp, &
      & 2.30333249960052E-02_wp, 2.04795146556047E-02_wp, 1.17566904494905E-03_wp, &
      & 1.12383332749237E-03_wp,-3.37161261865593E-04_wp, 7.55695280200463E-04_wp, &
      &-4.91653601592628E-04_wp,-1.02695166967481E-03_wp,-3.45289149212513E-02_wp, &
      &-1.34121866164477E-02_wp, 2.17842958083843E-02_wp,-1.32239449246803E-03_wp, &
      & 9.23763188397653E-04_wp,-5.87060101040320E-04_wp,-6.06266471532816E-04_wp, &
      &-5.25957119780333E-04_wp,-9.44887800903170E-04_wp, 9.71328579477840E-03_wp, &
      &-1.32264514234976E-02_wp,-4.00506954182355E-02_wp,-7.30975427202840E-04_wp, &
      &-2.86085329176629E-04_wp,-5.83829373771101E-04_wp, 1.06473924707788E-03_wp, &
      & 1.20864032141500E-03_wp,-4.66439474613685E-04_wp, 2.37196189260911E-03_wp, &
      & 3.76291798633030E-04_wp,-1.89524702707700E-04_wp,-2.70383567603791E-03_wp, &
      & 7.51687695827497E-04_wp, 1.57852557823706E-05_wp,-1.70421269429076E-03_wp, &
      &-1.24711040742963E-04_wp, 7.89473512438832E-04_wp,-3.57098061627412E-04_wp, &
      & 4.01348859464832E-03_wp,-3.38878489221253E-03_wp, 6.34879336473785E-03_wp, &
      &-2.30645681562275E-03_wp,-6.70882068405559E-04_wp, 6.96944723684135E-03_wp, &
      &-1.08436093353191E-02_wp,-1.77208852888238E-02_wp,-1.93565751770573E-02_wp, &
      &-8.75401591254688E-04_wp,-5.65296230452197E-04_wp, 1.94865564589505E-04_wp, &
      &-6.09542770971772E-04_wp,-6.86637147631842E-05_wp,-7.77156635669613E-03_wp, &
      &-1.12023415897402E-02_wp,-1.78971328952172E-02_wp, 2.23980508202744E-02_wp, &
      &-8.58573131874772E-04_wp, 6.50886246444007E-04_wp,-1.23589021153629E-04_wp, &
      &-7.14602097204805E-04_wp, 1.61848868701550E-04_wp, 2.47681579669712E-03_wp, &
      & 1.51505155319784E-02_wp, 4.13303227392449E-03_wp, 1.19937863232238E-02_wp, &
      &-4.00974483671558E-04_wp,-1.44661987626867E-04_wp,-2.20905355403101E-05_wp, &
      &-4.45102163842215E-04_wp,-5.21718641127864E-04_wp, 7.43124629353357E-04_wp, &
      & 3.76291798633030E-04_wp, 1.92419514812669E-03_wp, 2.43411943306366E-04_wp, &
      &-1.38558534277301E-03_wp,-1.29700970256260E-03_wp, 3.20274034780636E-04_wp, &
      &-6.00771636380427E-04_wp, 4.53126654298869E-04_wp, 6.89363513179920E-04_wp, &
      & 2.85458472275171E-03_wp,-8.42769400830715E-03_wp,-3.49285659188215E-03_wp, &
      &-8.73740458812932E-03_wp,-2.03180597856797E-03_wp, 3.57585309496399E-03_wp, &
      & 4.50241444891910E-03_wp,-1.94202483147024E-02_wp, 1.08132873675751E-02_wp, &
      &-1.31156713245789E-02_wp,-6.35109734804627E-04_wp,-1.41358099309022E-04_wp, &
      & 4.86112152718710E-04_wp,-6.63799363359947E-05_wp, 2.63148272706973E-04_wp, &
      &-3.69338959899031E-03_wp,-1.62131073853426E-02_wp, 1.12706765571304E-02_wp, &
      & 1.05297118737512E-02_wp,-6.83652635175080E-04_wp,-3.67543219007418E-06_wp, &
      &-5.64274993106786E-04_wp,-4.99490489386992E-05_wp,-1.70869611258782E-04_wp, &
      &-3.50216085250586E-03_wp, 4.54245265617874E-03_wp, 1.19055309005284E-02_wp, &
      &-1.89552860619114E-02_wp,-4.21614601949817E-04_wp,-8.16969303097783E-05_wp, &
      &-5.05453807145344E-04_wp, 3.80106625787291E-05_wp, 4.70419178269535E-04_wp, &
      &-3.09904488027013E-04_wp,-1.89524702707700E-04_wp, 2.43411943306366E-04_wp, &
      & 2.47928402744502E-03_wp,-1.13559737079954E-03_wp,-1.51157026310338E-03_wp, &
      &-1.74075334951878E-05_wp, 1.78516116698151E-03_wp, 2.08500104756173E-03_wp, &
      &-2.62166403396250E-03_wp, 1.57889816852864E-03_wp, 5.04088894829648E-03_wp, &
      &-2.33600703356005E-03_wp, 1.82408462816258E-03_wp, 3.82256307533336E-03_wp, &
      &-1.43702087113090E-03_wp, 4.61273866776342E-03_wp,-1.97237097978908E-02_wp, &
      &-1.20991357049096E-02_wp, 5.08336438917607E-03_wp,-2.83097626903721E-04_wp, &
      &-5.91992885184355E-04_wp, 1.38111314514985E-04_wp,-9.22681999257671E-05_wp, &
      & 7.37417254726881E-04_wp, 5.31578919345793E-03_wp, 2.10073897134136E-02_wp, &
      & 1.09429414011269E-02_wp, 3.20167044114457E-03_wp, 3.34046150737528E-04_wp, &
      &-7.19136030754418E-04_wp, 4.42987781466214E-05_wp, 2.30642331218448E-04_wp, &
      & 6.80365098404692E-04_wp,-8.84877256219468E-03_wp, 1.35788864585994E-02_wp, &
      &-2.03601998000453E-02_wp,-2.27788852326930E-02_wp,-6.78210598851248E-04_wp, &
      &-4.08503618156882E-04_wp,-1.59367877771823E-04_wp, 1.00614012061151E-03_wp, &
      & 6.42501172500844E-04_wp,-2.72087829478411E-03_wp,-2.70383567603791E-03_wp, &
      &-1.38558534277301E-03_wp,-1.13559737079954E-03_wp, 4.40223717129586E-02_wp, &
      & 5.57984043828882E-03_wp, 2.63086726948487E-03_wp, 7.70835864593437E-03_wp, &
      &-3.97661022352564E-04_wp, 2.26637457646800E-04_wp, 7.74500745277695E-03_wp, &
      &-1.79565826907083E-03_wp, 1.67084160854124E-02_wp,-6.61485609356795E-03_wp, &
      &-7.67995123149824E-04_wp, 2.18809239728138E-03_wp, 2.24397055345458E-02_wp, &
      & 6.64549541657225E-02_wp, 8.85871421287740E-02_wp, 6.96854518071820E-04_wp, &
      & 9.90846322150787E-04_wp, 2.46975374061967E-03_wp, 3.12047632637008E-04_wp, &
      & 1.29596421478787E-03_wp,-1.29311438090641E-03_wp,-2.40673203813859E-02_wp, &
      & 7.93390435947526E-02_wp, 7.02903965683136E-02_wp,-1.26672314141543E-02_wp, &
      & 1.72311948156705E-03_wp,-2.33142012668874E-03_wp, 1.41048483944506E-04_wp, &
      & 1.23754087718086E-03_wp, 1.32235237491202E-03_wp,-1.35903661755520E-02_wp, &
      & 7.65937840551549E-02_wp, 3.28686383454786E-02_wp, 6.44067474820537E-02_wp, &
      &-1.00233070466989E-03_wp,-6.88819688386836E-04_wp, 3.53857312558554E-04_wp, &
      &-1.40694498188145E-03_wp,-2.36540545013052E-03_wp, 6.74488694981665E-03_wp, &
      & 7.51687695827497E-04_wp,-1.29700970256260E-03_wp,-1.51157026310338E-03_wp, &
      & 5.57984043828882E-03_wp, 3.14464127057532E-02_wp, 5.42414731300381E-03_wp, &
      &-4.96716029907383E-03_wp,-7.32631149808548E-03_wp, 3.52486701113476E-03_wp, &
      &-1.72853100046916E-02_wp, 6.66597965239516E-03_wp, 1.41280628812287E-03_wp, &
      & 1.69880913473563E-03_wp, 6.38096027411938E-04_wp, 9.57106430383583E-03_wp, &
      & 2.07478739672486E-02_wp, 6.75325136801770E-02_wp,-2.15309610649779E-02_wp, &
      & 7.89233538672891E-02_wp, 2.81701658165517E-03_wp, 5.88757207551287E-04_wp, &
      &-1.82030938336181E-03_wp, 5.97326378130184E-04_wp,-5.00617661066083E-04_wp, &
      & 1.84804474365555E-02_wp,-6.44185602640536E-02_wp, 3.79320344196822E-02_wp, &
      & 6.87221148902045E-02_wp,-2.66981711043062E-03_wp, 3.40407076461401E-04_wp, &
      &-2.12977682744976E-03_wp,-3.25448102001268E-04_wp,-6.09023766965693E-04_wp, &
      &-5.22282280974187E-03_wp, 3.75844989308166E-02_wp,-1.69469206157941E-02_wp, &
      & 4.36814617080786E-02_wp,-7.35276572443060E-04_wp,-6.06690581189437E-04_wp, &
      & 6.06034207884786E-04_wp, 1.10660155876927E-04_wp,-1.13219077295543E-03_wp, &
      &-5.28266720158637E-03_wp, 1.57852557823706E-05_wp, 3.20274034780636E-04_wp, &
      &-1.74075334951878E-05_wp, 2.63086726948487E-03_wp, 5.42414731300381E-03_wp, &
      & 2.93589728175265E-02_wp, 3.67709738831447E-03_wp,-7.47063038603620E-04_wp, &
      & 3.90836704630479E-03_wp, 7.09621952656464E-05_wp, 2.07760552467874E-03_wp, &
      &-1.29405483046609E-04_wp, 2.89718143673136E-03_wp, 3.16800401176324E-04_wp, &
      &-1.69606053549689E-02_wp,-7.54471811196526E-03_wp, 1.55318190585747E-02_wp, &
      &-9.64076933204909E-02_wp, 1.10515854057054E-02_wp, 1.45947727313188E-04_wp, &
      &-2.16103631752938E-03_wp,-1.41298732458204E-03_wp,-1.48077115538771E-03_wp, &
      &-3.68872484730471E-05_wp,-9.73129495503449E-03_wp, 3.62507400525617E-04_wp, &
      & 8.64415228683305E-02_wp,-5.52582904665705E-03_wp, 2.36297575249206E-04_wp, &
      &-1.94458556002641E-03_wp,-1.28072235876254E-03_wp, 1.67826970775080E-03_wp, &
      &-2.64680965053990E-04_wp,-9.91183630949829E-03_wp,-5.57175718616085E-03_wp, &
      & 8.67153633545441E-02_wp, 2.90410287476276E-03_wp, 3.78053496299328E-04_wp, &
      & 9.38622543054151E-04_wp,-1.22901953593365E-03_wp,-2.36862461809501E-03_wp, &
      & 3.42392873712068E-05_wp, 4.80689218346076E-03_wp,-1.70421269429076E-03_wp, &
      &-6.00771636380427E-04_wp, 1.78516116698151E-03_wp, 7.70835864593437E-03_wp, &
      &-4.96716029907383E-03_wp, 3.67709738831447E-03_wp, 3.44193006387433E-02_wp, &
      & 4.76228788479840E-03_wp,-1.12594099043911E-02_wp,-4.22637413690153E-03_wp, &
      &-1.10988128220400E-02_wp, 9.39053014236173E-04_wp,-1.16060082172845E-03_wp, &
      & 1.07624737460322E-02_wp, 5.60381954566253E-03_wp, 1.41316668863599E-02_wp, &
      & 7.88193671551557E-02_wp,-1.33489350096091E-02_wp, 5.13853337547530E-03_wp, &
      & 1.55817983246274E-03_wp, 5.12251036343349E-04_wp,-1.38894336044170E-03_wp, &
      & 4.98434304022621E-06_wp,-1.55155659357652E-03_wp,-1.19435034444835E-02_wp, &
      & 7.87707633106604E-02_wp,-3.00272463633091E-02_wp,-2.09804044214836E-03_wp, &
      & 1.14598547458529E-03_wp,-1.69273880680804E-04_wp, 1.52074117567880E-03_wp, &
      &-6.69143343856796E-04_wp, 1.75734486260590E-03_wp, 2.14421083889722E-02_wp, &
      & 3.28723251997694E-02_wp, 4.78972976105153E-02_wp,-1.04712017981954E-01_wp, &
      &-1.59454356614859E-03_wp,-7.29860323367932E-05_wp,-2.39582887843716E-03_wp, &
      & 4.78040749240586E-04_wp, 2.59547000700204E-03_wp, 9.22601160310286E-04_wp, &
      &-1.24711040742963E-04_wp, 4.53126654298869E-04_wp, 2.08500104756173E-03_wp, &
      &-3.97661022352564E-04_wp,-7.32631149808548E-03_wp,-7.47063038603620E-04_wp, &
      & 4.76228788479840E-03_wp, 4.33364941711815E-02_wp,-3.15071146362812E-03_wp, &
      &-6.67108487319752E-04_wp, 6.99176018534552E-03_wp,-6.22633014911112E-03_wp, &
      &-3.38055775952816E-03_wp, 7.69017399473797E-03_wp,-1.76961128934806E-03_wp, &
      &-8.85374052511765E-03_wp, 3.47268154939371E-02_wp,-3.45265330160401E-02_wp, &
      &-8.92468896301965E-02_wp,-1.30909318966849E-03_wp,-3.70938201881427E-04_wp, &
      &-1.10750429754854E-04_wp,-1.31927484389760E-03_wp,-1.78842876359272E-03_wp, &
      &-8.63307150232750E-03_wp,-3.79182506231862E-02_wp, 3.52704972597678E-02_wp, &
      &-9.70603235678339E-02_wp, 1.09039808896316E-03_wp,-3.59947499872617E-04_wp, &
      &-3.60824594442789E-04_wp, 1.85106237168900E-03_wp,-2.29662810010202E-03_wp, &
      & 2.20191489255154E-02_wp, 7.33216018913181E-02_wp,-7.19861577004959E-02_wp, &
      &-4.69373915752494E-02_wp,-2.10917504939881E-03_wp,-1.45840439450965E-03_wp, &
      & 1.08123286917808E-04_wp, 2.31640389941183E-03_wp, 5.16614467009417E-04_wp, &
      &-2.71966921742355E-03_wp, 7.89473512438832E-04_wp, 6.89363513179920E-04_wp, &
      &-2.62166403396250E-03_wp, 2.26637457646800E-04_wp, 3.52486701113476E-03_wp, &
      & 3.90836704630479E-03_wp,-1.12594099043911E-02_wp,-3.15071146362812E-03_wp, &
      & 2.34151430731864E-02_wp, 1.85415808477836E-02_wp,-4.24734078204873E-02_wp, &
      & 2.35575691065080E-03_wp,-3.33058390368733E-02_wp,-5.57795006981639E-03_wp, &
      & 1.01108773593490E-02_wp, 4.01102854511713E-04_wp,-4.56200941512203E-02_wp, &
      &-6.60382781523261E-03_wp, 4.54256036670506E-02_wp, 8.13365306994022E-04_wp, &
      &-1.45138153891137E-03_wp,-4.61944048276185E-04_wp, 2.51862372218516E-04_wp, &
      & 2.07065310902566E-03_wp, 3.41550940627179E-04_wp,-5.80400261978655E-02_wp, &
      &-1.74917522456011E-02_wp,-5.30554321998536E-02_wp, 5.62332461341726E-04_wp, &
      & 1.78951892334922E-03_wp, 5.49306638141359E-04_wp, 1.18354420618126E-04_wp, &
      &-2.60744135287044E-03_wp,-1.84935397407614E-04_wp, 4.00823788730297E-02_wp, &
      & 5.07256361601547E-02_wp, 6.01728319805672E-02_wp,-1.08074030232713E-03_wp, &
      &-5.41519705134412E-04_wp, 1.74128670059472E-04_wp,-2.34898043486087E-03_wp, &
      &-2.72927683099535E-03_wp,-1.01894320855533E-02_wp,-3.57098061627412E-04_wp, &
      & 2.85458472275171E-03_wp, 1.57889816852864E-03_wp, 7.74500745277695E-03_wp, &
      &-1.72853100046916E-02_wp, 7.09621952656464E-05_wp,-4.22637413690153E-03_wp, &
      &-6.67108487319752E-04_wp, 1.85415808477836E-02_wp, 7.81937062031174E-02_wp, &
      &-1.56952868345472E-01_wp, 4.63674964411300E-03_wp,-1.52201455602624E-01_wp, &
      &-3.58517678616803E-02_wp, 3.07585591512265E-02_wp, 4.52883325072763E-04_wp, &
      &-9.32147819443877E-02_wp,-5.48349175380994E-03_wp,-3.42959429035476E-02_wp, &
      &-1.17137609070403E-03_wp,-4.30791961546821E-03_wp,-8.16503856575622E-04_wp, &
      &-1.99648246109626E-03_wp, 1.66375144799501E-03_wp, 7.55268864056282E-05_wp, &
      & 7.36998183172769E-02_wp, 8.59929205900434E-03_wp,-3.29569302511059E-02_wp, &
      & 2.35295984835294E-03_wp,-7.98432581416592E-04_wp, 2.29759184132251E-03_wp, &
      &-6.70858028794763E-04_wp, 2.19545673464568E-03_wp, 9.16685372743052E-04_wp, &
      & 2.57126867846738E-02_wp, 1.14966740670106E-02_wp, 3.88665260783668E-02_wp, &
      &-9.20983705270224E-04_wp,-1.88651521181885E-03_wp, 2.10285595596077E-03_wp, &
      & 5.97064769894395E-04_wp,-3.32760114542134E-03_wp,-2.27063129413855E-03_wp, &
      & 4.01348859464832E-03_wp,-8.42769400830715E-03_wp, 5.04088894829648E-03_wp, &
      &-1.79565826907083E-03_wp, 6.66597965239516E-03_wp, 2.07760552467874E-03_wp, &
      &-1.10988128220400E-02_wp, 6.99176018534552E-03_wp,-4.24734078204873E-02_wp, &
      &-1.56952868345472E-01_wp, 4.62813275841273E-01_wp,-9.10464660774819E-03_wp, &
      & 4.38047333432606E-01_wp, 9.59547188118657E-02_wp,-1.25382771228594E-01_wp, &
      &-1.38374238093702E-03_wp,-7.41273093138561E-03_wp,-1.10468068244529E-02_wp, &
      &-1.56187513126319E-03_wp,-4.41907081296685E-03_wp, 6.66432321276951E-03_wp, &
      & 6.93367015797582E-03_wp, 3.50242290466801E-03_wp, 1.53391697941402E-03_wp, &
      & 5.43354709716936E-04_wp, 1.40237003812439E-02_wp, 3.28804253492689E-02_wp, &
      & 1.21189122316512E-02_wp,-6.66975469219293E-05_wp,-3.13907873268688E-03_wp, &
      &-3.25713971009019E-03_wp, 4.31017066296889E-03_wp,-1.47887531130849E-03_wp, &
      &-1.84087838689988E-04_wp, 1.41827982066621E-02_wp,-6.62740370087469E-02_wp, &
      & 2.39879765902584E-02_wp, 1.05566543856893E-03_wp, 3.60304872454716E-03_wp, &
      &-3.54139161540637E-03_wp,-2.95929713093146E-03_wp, 4.28984898113645E-03_wp, &
      &-2.82181425185299E-03_wp,-3.38878489221253E-03_wp,-3.49285659188215E-03_wp, &
      &-2.33600703356005E-03_wp, 1.67084160854124E-02_wp, 1.41280628812287E-03_wp, &
      &-1.29405483046609E-04_wp, 9.39053014236173E-04_wp,-6.22633014911112E-03_wp, &
      & 2.35575691065080E-03_wp, 4.63674964411300E-03_wp,-9.10464660774819E-03_wp, &
      & 1.55693316213506E-02_wp,-8.69456073030934E-03_wp,-2.07001602169140E-03_wp, &
      &-1.64641127449032E-03_wp, 7.40412926616573E-04_wp, 5.19987722271890E-02_wp, &
      & 3.50963776224895E-02_wp, 3.53821139878432E-02_wp, 2.03766503438216E-03_wp, &
      & 1.38182723045576E-03_wp,-8.56882698100729E-04_wp, 9.67443679055029E-04_wp, &
      &-7.89383163152661E-04_wp,-3.12016989034458E-04_wp, 4.89673716424246E-02_wp, &
      & 6.85566681173991E-03_wp,-4.33726604714417E-02_wp, 2.32152383555904E-03_wp, &
      &-7.19828442743159E-04_wp, 1.28190858547835E-03_wp, 7.83301665643090E-04_wp, &
      & 1.95937371861310E-04_wp,-2.19901739923401E-04_wp,-2.48692793280382E-02_wp, &
      & 6.56157086460467E-03_wp, 6.03966283303314E-02_wp, 1.66757073509023E-03_wp, &
      & 4.66111281355355E-04_wp, 1.29891259805378E-03_wp,-8.75032383974094E-04_wp, &
      &-1.61471980448914E-03_wp,-2.14613395036928E-03_wp, 6.34879336473785E-03_wp, &
      &-8.73740458812932E-03_wp, 1.82408462816258E-03_wp,-6.61485609356795E-03_wp, &
      & 1.69880913473563E-03_wp, 2.89718143673136E-03_wp,-1.16060082172845E-03_wp, &
      &-3.38055775952816E-03_wp,-3.33058390368733E-02_wp,-1.52201455602624E-01_wp, &
      & 4.38047333432606E-01_wp,-8.69456073030934E-03_wp, 4.36466806645170E-01_wp, &
      & 9.34741416386527E-02_wp,-1.19354104143835E-01_wp,-1.37121127231969E-03_wp, &
      & 5.30583940781882E-04_wp,-2.67897186114664E-02_wp, 4.59737348138373E-03_wp, &
      &-4.02177217080790E-03_wp, 6.15256712973667E-03_wp, 6.28787735416212E-03_wp, &
      & 3.19802362534811E-03_wp, 1.32955044990478E-03_wp, 1.31200798097231E-03_wp, &
      & 1.24744952830266E-02_wp,-7.32813987752802E-02_wp, 3.72566964773942E-03_wp, &
      & 1.02640693286939E-04_wp,-1.47429100000553E-05_wp,-1.17086414461730E-03_wp, &
      & 2.24601344723718E-03_wp,-2.00559116494408E-03_wp,-1.12586025628635E-03_wp, &
      & 2.39815493205939E-02_wp, 2.95058137735713E-02_wp, 2.61728970657019E-02_wp, &
      & 5.32051158886537E-04_wp, 4.14351970336048E-03_wp,-5.00522584023764E-03_wp, &
      &-6.01050825620884E-03_wp, 3.53071480912169E-03_wp, 3.61217559584011E-03_wp, &
      &-2.30645681562275E-03_wp,-2.03180597856797E-03_wp, 3.82256307533336E-03_wp, &
      &-7.67995123149824E-04_wp, 6.38096027411938E-04_wp, 3.16800401176324E-04_wp, &
      & 1.07624737460322E-02_wp, 7.69017399473797E-03_wp,-5.57795006981639E-03_wp, &
      &-3.58517678616803E-02_wp, 9.59547188118657E-02_wp,-2.07001602169140E-03_wp, &
      & 9.34741416386527E-02_wp, 3.90520218771133E-02_wp,-2.06735101101431E-02_wp, &
      &-2.07679312335869E-04_wp,-3.36427527890217E-03_wp,-2.68963265017502E-03_wp, &
      & 7.67642844476133E-02_wp, 1.35665085458970E-03_wp, 1.35542875743456E-03_wp, &
      & 7.30947572371426E-04_wp, 2.17212419192107E-03_wp, 2.27422974116808E-03_wp, &
      &-5.83309848297215E-04_wp, 1.33899406090693E-02_wp,-2.03214366844796E-03_wp, &
      &-5.90050199661660E-02_wp, 1.83884922532074E-03_wp,-5.27133217598405E-04_wp, &
      & 3.10277119084791E-04_wp, 1.78806056495109E-03_wp,-1.42996134361120E-03_wp, &
      & 9.38729389402143E-05_wp, 4.82654826901009E-02_wp,-1.14992339504637E-02_wp, &
      &-6.93910045600238E-02_wp,-1.97145781885208E-03_wp, 3.56893172299223E-05_wp, &
      &-2.34151302886779E-03_wp, 4.72722357628950E-04_wp, 2.70282186006250E-03_wp, &
      & 1.12028617141975E-02_wp,-6.70882068405559E-04_wp, 3.57585309496399E-03_wp, &
      &-1.43702087113090E-03_wp, 2.18809239728138E-03_wp, 9.57106430383583E-03_wp, &
      &-1.69606053549689E-02_wp, 5.60381954566253E-03_wp,-1.76961128934806E-03_wp, &
      & 1.01108773593490E-02_wp, 3.07585591512265E-02_wp,-1.25382771228594E-01_wp, &
      &-1.64641127449032E-03_wp,-1.19354104143835E-01_wp,-2.06735101101431E-02_wp, &
      & 5.64165987847063E-02_wp,-3.60508414656166E-04_wp, 2.36891004263969E-02_wp, &
      & 7.45384160365143E-02_wp, 3.27358666275927E-02_wp, 2.50498630191137E-03_wp, &
      & 3.12765956419526E-04_wp,-1.49064981294504E-03_wp, 7.64346376836081E-04_wp, &
      &-3.81834224001756E-04_wp,-3.89914338100669E-04_wp,-5.16374969918261E-02_wp, &
      &-4.56666343914478E-02_wp, 4.41302193900175E-02_wp,-2.26312525914848E-03_wp, &
      & 2.64354771088751E-03_wp, 1.76499596098071E-04_wp,-2.64863626414187E-03_wp, &
      &-2.83279606761995E-05_wp,-2.41011357456463E-04_wp, 5.49357483690476E-02_wp, &
      &-2.23923554739493E-02_wp,-3.95054115593926E-02_wp,-2.59058745650320E-03_wp, &
      &-2.26084037290193E-03_wp, 4.94952180347453E-04_wp, 2.52712586784912E-03_wp, &
      &-5.02041178188871E-04_wp,-8.64033744991843E-04_wp, 6.96944723684135E-03_wp, &
      & 4.50241444891910E-03_wp, 4.61273866776342E-03_wp, 2.24397055345458E-02_wp, &
      & 2.07478739672486E-02_wp,-7.54471811196526E-03_wp, 1.41316668863599E-02_wp, &
      &-8.85374052511765E-03_wp, 4.01102854511713E-04_wp, 4.52883325072763E-04_wp, &
      &-1.38374238093702E-03_wp, 7.40412926616573E-04_wp,-1.37121127231969E-03_wp, &
      &-2.07679312335869E-04_wp,-3.60508414656166E-04_wp, 9.79385233275982E-01_wp, &
      &-6.73932117684745E-02_wp,-4.24221242518477E-02_wp,-4.54829773946338E-02_wp, &
      &-3.44767227988328E-03_wp,-3.20186908954779E-03_wp, 1.08518088939978E-03_wp, &
      &-2.17906881067328E-03_wp, 1.36803923011760E-03_wp,-1.74429112109329E-03_wp, &
      &-4.60039554544568E-03_wp,-4.62867164240676E-03_wp,-2.29461193380960E-02_wp, &
      & 1.08781582530057E-03_wp, 8.21284940311572E-04_wp, 8.63411726541678E-04_wp, &
      & 3.71234417787815E-04_wp,-2.26885600405926E-03_wp,-4.67052712125803E-04_wp, &
      &-2.18177810703489E-02_wp,-4.85160230141416E-03_wp, 2.44761467320601E-03_wp, &
      & 2.10650456051256E-03_wp, 5.13537734104285E-04_wp, 8.26572953278356E-04_wp, &
      & 7.14219643152562E-04_wp, 1.14929051465422E-03_wp, 3.10371374558179E-02_wp, &
      &-1.08436093353191E-02_wp,-1.94202483147024E-02_wp,-1.97237097978908E-02_wp, &
      & 6.64549541657225E-02_wp, 6.75325136801770E-02_wp, 1.55318190585747E-02_wp, &
      & 7.88193671551557E-02_wp, 3.47268154939371E-02_wp,-4.56200941512203E-02_wp, &
      &-9.32147819443877E-02_wp,-7.41273093138561E-03_wp, 5.19987722271890E-02_wp, &
      & 5.30583940781882E-04_wp,-3.36427527890217E-03_wp, 2.36891004263969E-02_wp, &
      &-6.73932117684745E-02_wp, 8.89614738721827E-01_wp,-1.28950331967082E-02_wp, &
      &-1.39586534564777E-02_wp, 1.53963342226492E-02_wp, 1.29909255552899E-02_wp, &
      &-1.31338339197786E-02_wp,-1.58767803862221E-05_wp,-2.41453153734097E-02_wp, &
      & 5.83630129305121E-03_wp, 4.74494717503912E-03_wp,-1.32063220712849E-02_wp, &
      &-1.19859393613802E-02_wp,-6.98434276122301E-04_wp, 2.29788467999063E-03_wp, &
      & 3.50868101513821E-04_wp,-9.19684784734447E-04_wp,-1.01534697330281E-03_wp, &
      &-1.79249305334720E-02_wp,-4.44215180079783E-02_wp,-4.45142330247158E-03_wp, &
      & 3.22808677054423E-02_wp, 5.51511024939449E-03_wp, 2.23194762525990E-03_wp, &
      & 2.14969208354629E-03_wp, 2.26334563113351E-04_wp,-2.69348432555279E-04_wp, &
      & 2.30333249960052E-02_wp,-1.77208852888238E-02_wp, 1.08132873675751E-02_wp, &
      &-1.20991357049096E-02_wp, 8.85871421287740E-02_wp,-2.15309610649779E-02_wp, &
      &-9.64076933204909E-02_wp,-1.33489350096091E-02_wp,-3.45265330160401E-02_wp, &
      &-6.60382781523261E-03_wp,-5.48349175380994E-03_wp,-1.10468068244529E-02_wp, &
      & 3.50963776224895E-02_wp,-2.67897186114664E-02_wp,-2.68963265017502E-03_wp, &
      & 7.45384160365143E-02_wp,-4.24221242518477E-02_wp,-1.28950331967082E-02_wp, &
      & 9.00965756211415E-01_wp,-9.14415453328002E-03_wp, 1.33841889817848E-03_wp, &
      & 2.15951288119083E-02_wp, 1.24101348943288E-02_wp, 1.49920354704816E-02_wp, &
      &-4.40469099717627E-04_wp, 9.85577803669907E-04_wp,-1.22201335362130E-02_wp, &
      & 2.11341370625129E-02_wp,-1.39791945191686E-02_wp,-7.01721881174691E-04_wp, &
      & 1.08687409943675E-04_wp,-1.44491706265718E-03_wp, 4.54075284812082E-04_wp, &
      &-1.98967421502637E-03_wp, 1.69947916255744E-03_wp,-1.66229817006404E-02_wp, &
      & 2.23035839446687E-02_wp,-9.27699182882334E-03_wp, 4.79918109489262E-04_wp, &
      & 1.40000483985680E-04_wp,-1.14187690595404E-03_wp, 2.76159075758330E-04_wp, &
      & 1.54213587965840E-03_wp, 2.04795146556047E-02_wp,-1.93565751770573E-02_wp, &
      &-1.31156713245789E-02_wp, 5.08336438917607E-03_wp, 6.96854518071820E-04_wp, &
      & 7.89233538672891E-02_wp, 1.10515854057054E-02_wp, 5.13853337547530E-03_wp, &
      &-8.92468896301965E-02_wp, 4.54256036670506E-02_wp,-3.42959429035476E-02_wp, &
      &-1.56187513126319E-03_wp, 3.53821139878432E-02_wp, 4.59737348138373E-03_wp, &
      & 7.67642844476133E-02_wp, 3.27358666275927E-02_wp,-4.54829773946338E-02_wp, &
      &-1.39586534564777E-02_wp,-9.14415453328002E-03_wp, 9.01197413335927E-01_wp, &
      & 2.37424843742069E-02_wp, 2.21729555499775E-04_wp,-8.75069972743066E-03_wp, &
      & 1.29168909226540E-02_wp, 1.70306110002357E-02_wp,-2.30415156877423E-02_wp, &
      & 1.66719242999173E-02_wp, 2.25021243603302E-03_wp,-5.05076972762294E-02_wp, &
      & 5.43204781193947E-03_wp,-9.57154239876852E-05_wp, 2.35591709642209E-03_wp, &
      & 2.63084615326078E-03_wp,-2.80062748714972E-03_wp, 1.34259756367200E-02_wp, &
      & 9.67852731635338E-03_wp,-1.30717775428052E-02_wp,-3.47522554510675E-03_wp, &
      &-2.41917211488738E-03_wp,-1.02717741613644E-03_wp,-7.37519467376673E-04_wp, &
      & 2.37177085330210E-03_wp, 1.59044196104408E-03_wp, 1.17566904494905E-03_wp, &
      &-8.75401591254688E-04_wp,-6.35109734804627E-04_wp,-2.83097626903721E-04_wp, &
      & 9.90846322150787E-04_wp, 2.81701658165517E-03_wp, 1.45947727313188E-04_wp, &
      & 1.55817983246274E-03_wp,-1.30909318966849E-03_wp, 8.13365306994022E-04_wp, &
      &-1.17137609070403E-03_wp,-4.41907081296685E-03_wp, 2.03766503438216E-03_wp, &
      &-4.02177217080790E-03_wp, 1.35665085458970E-03_wp, 2.50498630191137E-03_wp, &
      &-3.44767227988328E-03_wp, 1.53963342226492E-02_wp, 1.33841889817848E-03_wp, &
      & 2.37424843742069E-02_wp, 9.71080523127880E-04_wp, 2.24537932149475E-04_wp, &
      &-5.12935597022903E-04_wp, 3.54826393034930E-04_wp, 1.78175387132417E-05_wp, &
      &-9.92864289497205E-04_wp,-4.93805000442397E-04_wp,-1.05118342011321E-03_wp, &
      &-5.15532684544786E-03_wp, 2.01568567478042E-04_wp, 8.47489533143351E-05_wp, &
      & 1.28915036185723E-04_wp, 5.96716941142845E-05_wp,-1.61436989305193E-04_wp, &
      & 9.05992698760134E-04_wp,-1.86454186027253E-03_wp,-1.56221941951145E-03_wp, &
      &-8.01779151303163E-04_wp, 5.29754230196242E-05_wp,-1.67720095140052E-05_wp, &
      & 6.88550780570842E-05_wp, 1.49459436389064E-04_wp, 3.48745958380269E-05_wp, &
      & 1.12383332749237E-03_wp,-5.65296230452197E-04_wp,-1.41358099309022E-04_wp, &
      &-5.91992885184355E-04_wp, 2.46975374061967E-03_wp, 5.88757207551287E-04_wp, &
      &-2.16103631752938E-03_wp, 5.12251036343349E-04_wp,-3.70938201881427E-04_wp, &
      &-1.45138153891137E-03_wp,-4.30791961546821E-03_wp, 6.66432321276951E-03_wp, &
      & 1.38182723045576E-03_wp, 6.15256712973667E-03_wp, 1.35542875743456E-03_wp, &
      & 3.12765956419526E-04_wp,-3.20186908954779E-03_wp, 1.29909255552899E-02_wp, &
      & 2.15951288119083E-02_wp, 2.21729555499775E-04_wp, 2.24537932149475E-04_wp, &
      & 8.41821697249052E-04_wp, 2.01357871246874E-04_wp, 4.33545604698566E-04_wp, &
      &-3.54326607641045E-04_wp, 6.92419916187062E-04_wp,-2.26880794834199E-03_wp, &
      &-2.60472690055268E-04_wp,-6.40001479013474E-04_wp,-5.97084418745616E-05_wp, &
      & 4.30725140314870E-05_wp,-9.13448849034343E-05_wp, 5.35251386105152E-05_wp, &
      &-1.48760158878772E-04_wp,-3.65628894307319E-04_wp,-3.34642237400698E-03_wp, &
      &-3.20723614998436E-04_wp,-3.02716907708290E-04_wp, 1.66642160459250E-04_wp, &
      & 1.29518345210917E-04_wp,-6.06645445284161E-05_wp,-3.58157141410844E-05_wp, &
      & 1.47884993200396E-04_wp,-3.37161261865593E-04_wp, 1.94865564589505E-04_wp, &
      & 4.86112152718710E-04_wp, 1.38111314514985E-04_wp, 3.12047632637008E-04_wp, &
      &-1.82030938336181E-03_wp,-1.41298732458204E-03_wp,-1.38894336044170E-03_wp, &
      &-1.10750429754854E-04_wp,-4.61944048276185E-04_wp,-8.16503856575622E-04_wp, &
      & 6.93367015797582E-03_wp,-8.56882698100729E-04_wp, 6.28787735416212E-03_wp, &
      & 7.30947572371426E-04_wp,-1.49064981294504E-03_wp, 1.08518088939978E-03_wp, &
      &-1.31338339197786E-02_wp, 1.24101348943288E-02_wp,-8.75069972743066E-03_wp, &
      &-5.12935597022903E-04_wp, 2.01357871246874E-04_wp, 5.57704342368946E-04_wp, &
      & 1.30088685792905E-04_wp, 2.05807133618905E-04_wp, 7.04869867155257E-04_wp, &
      &-4.71965410417172E-04_wp, 1.54322398343485E-03_wp, 1.57912740212245E-03_wp, &
      &-8.60224070526115E-05_wp,-8.26444834398727E-05_wp,-1.22481720956570E-04_wp, &
      & 5.33348705135931E-05_wp,-3.16133131415677E-06_wp, 6.83482385411022E-04_wp, &
      & 1.43421135623700E-03_wp, 1.29954197474567E-03_wp,-6.84034941755676E-04_wp, &
      &-6.77514813142387E-05_wp, 3.57114046682632E-05_wp,-1.37991308366808E-04_wp, &
      &-1.10652382973498E-04_wp, 8.18569440156331E-05_wp, 7.55695280200463E-04_wp, &
      &-6.09542770971772E-04_wp,-6.63799363359947E-05_wp,-9.22681999257671E-05_wp, &
      & 1.29596421478787E-03_wp, 5.97326378130184E-04_wp,-1.48077115538771E-03_wp, &
      & 4.98434304022621E-06_wp,-1.31927484389760E-03_wp, 2.51862372218516E-04_wp, &
      &-1.99648246109626E-03_wp, 3.50242290466801E-03_wp, 9.67443679055029E-04_wp, &
      & 3.19802362534811E-03_wp, 2.17212419192107E-03_wp, 7.64346376836081E-04_wp, &
      &-2.17906881067328E-03_wp,-1.58767803862221E-05_wp, 1.49920354704816E-02_wp, &
      & 1.29168909226540E-02_wp, 3.54826393034930E-04_wp, 4.33545604698566E-04_wp, &
      & 1.30088685792905E-04_wp, 4.82820109678503E-04_wp, 2.51108863692765E-04_wp, &
      &-7.79172865409633E-04_wp,-7.66705658334690E-04_wp, 2.26017817007804E-04_wp, &
      &-3.30674201547255E-03_wp, 1.03484025349429E-04_wp,-1.06917681623253E-06_wp, &
      &-1.97825617795629E-07_wp, 1.08987968326367E-04_wp,-1.54513083585891E-04_wp, &
      & 9.18499763895862E-04_wp, 2.49064579621402E-04_wp,-3.01515131463446E-04_wp, &
      &-1.87144089306515E-03_wp,-4.15809966344942E-05_wp, 1.50166180257202E-05_wp, &
      &-9.17945277868815E-05_wp, 3.21446802816235E-05_wp, 1.29877997003160E-04_wp, &
      &-4.91653601592628E-04_wp,-6.86637147631842E-05_wp, 2.63148272706973E-04_wp, &
      & 7.37417254726881E-04_wp,-1.29311438090641E-03_wp,-5.00617661066083E-04_wp, &
      &-3.68872484730471E-05_wp,-1.55155659357652E-03_wp,-1.78842876359272E-03_wp, &
      & 2.07065310902566E-03_wp, 1.66375144799501E-03_wp, 1.53391697941402E-03_wp, &
      &-7.89383163152661E-04_wp, 1.32955044990478E-03_wp, 2.27422974116808E-03_wp, &
      &-3.81834224001756E-04_wp, 1.36803923011760E-03_wp,-2.41453153734097E-02_wp, &
      &-4.40469099717627E-04_wp, 1.70306110002357E-02_wp, 1.78175387132417E-05_wp, &
      &-3.54326607641045E-04_wp, 2.05807133618905E-04_wp, 2.51108863692765E-04_wp, &
      & 1.00296300267918E-03_wp,-2.26339054665830E-03_wp, 1.22673110868470E-03_wp, &
      & 1.42194057769701E-03_wp,-2.96462717477253E-03_wp, 1.86342724132078E-04_wp, &
      &-1.01495495974086E-04_wp, 4.60297941197738E-05_wp, 1.31730529441688E-04_wp, &
      &-5.02008128433043E-05_wp, 2.19806549831163E-03_wp, 5.36800952916396E-03_wp, &
      &-1.11865135939340E-04_wp,-2.11155401443376E-03_wp,-3.09423947555397E-04_wp, &
      &-1.21929197699361E-04_wp,-1.22852682598191E-04_wp, 4.05393279050851E-05_wp, &
      & 4.47136744942552E-05_wp,-1.02695166967481E-03_wp,-7.77156635669613E-03_wp, &
      &-3.69338959899031E-03_wp, 5.31578919345793E-03_wp,-2.40673203813859E-02_wp, &
      & 1.84804474365555E-02_wp,-9.73129495503449E-03_wp,-1.19435034444835E-02_wp, &
      &-8.63307150232750E-03_wp, 3.41550940627179E-04_wp, 7.55268864056282E-05_wp, &
      & 5.43354709716936E-04_wp,-3.12016989034458E-04_wp, 1.31200798097231E-03_wp, &
      &-5.83309848297215E-04_wp,-3.89914338100669E-04_wp,-1.74429112109329E-03_wp, &
      & 5.83630129305121E-03_wp, 9.85577803669907E-04_wp,-2.30415156877423E-02_wp, &
      &-9.92864289497205E-04_wp, 6.92419916187062E-04_wp, 7.04869867155257E-04_wp, &
      &-7.79172865409633E-04_wp,-2.26339054665830E-03_wp, 9.80319653876711E-01_wp, &
      & 6.96148081109755E-02_wp, 3.39853186201215E-02_wp,-4.74581076068379E-02_wp, &
      & 3.77653932652044E-03_wp,-2.82491840284670E-03_wp, 1.53601953767155E-03_wp, &
      & 1.88938374714825E-03_wp, 1.39873724097239E-03_wp,-3.03607735046782E-03_wp, &
      & 1.59310863120059E-02_wp, 1.87885921811414E-02_wp, 5.08467535689386E-03_wp, &
      &-1.69023816854932E-03_wp,-6.52194452941343E-04_wp,-1.24575925065902E-03_wp, &
      &-1.84891013912983E-03_wp,-3.95093353420261E-04_wp,-3.45289149212513E-02_wp, &
      &-1.12023415897402E-02_wp,-1.62131073853426E-02_wp, 2.10073897134136E-02_wp, &
      & 7.93390435947526E-02_wp,-6.44185602640536E-02_wp, 3.62507400525617E-04_wp, &
      & 7.87707633106604E-02_wp,-3.79182506231862E-02_wp,-5.80400261978655E-02_wp, &
      & 7.36998183172769E-02_wp, 1.40237003812439E-02_wp, 4.89673716424246E-02_wp, &
      & 1.24744952830266E-02_wp, 1.33899406090693E-02_wp,-5.16374969918261E-02_wp, &
      &-4.60039554544568E-03_wp, 4.74494717503912E-03_wp,-1.22201335362130E-02_wp, &
      & 1.66719242999173E-02_wp,-4.93805000442397E-04_wp,-2.26880794834199E-03_wp, &
      &-4.71965410417172E-04_wp,-7.66705658334690E-04_wp, 1.22673110868470E-03_wp, &
      & 6.96148081109755E-02_wp, 8.84483522477175E-01_wp,-6.27087331795843E-03_wp, &
      & 1.39890422776412E-02_wp, 1.70641049794269E-02_wp,-1.41030186631462E-02_wp, &
      & 1.27499204873158E-02_wp, 3.86236907846822E-04_wp, 2.29766788919789E-02_wp, &
      & 9.93178929003970E-03_wp,-2.68014119418027E-03_wp,-3.97830522408287E-02_wp, &
      &-1.26975405204845E-02_wp, 1.43642412751507E-03_wp,-5.15623780963257E-04_wp, &
      & 1.97029419194627E-03_wp, 3.39979714316641E-03_wp, 1.57372281682921E-04_wp, &
      &-1.34121866164477E-02_wp,-1.78971328952172E-02_wp, 1.12706765571304E-02_wp, &
      & 1.09429414011269E-02_wp, 7.02903965683136E-02_wp, 3.79320344196822E-02_wp, &
      & 8.64415228683305E-02_wp,-3.00272463633091E-02_wp, 3.52704972597678E-02_wp, &
      &-1.74917522456011E-02_wp, 8.59929205900434E-03_wp, 3.28804253492689E-02_wp, &
      & 6.85566681173991E-03_wp,-7.32813987752802E-02_wp,-2.03214366844796E-03_wp, &
      &-4.56666343914478E-02_wp,-4.62867164240676E-03_wp,-1.32063220712849E-02_wp, &
      & 2.11341370625129E-02_wp, 2.25021243603302E-03_wp,-1.05118342011321E-03_wp, &
      &-2.60472690055268E-04_wp, 1.54322398343485E-03_wp, 2.26017817007804E-04_wp, &
      & 1.42194057769701E-03_wp, 3.39853186201215E-02_wp,-6.27087331795843E-03_wp, &
      & 9.01725438066412E-01_wp, 8.33818371867964E-03_wp,-3.41304494398283E-04_wp, &
      &-2.29885786612005E-02_wp,-1.59420543751472E-02_wp, 1.67323852681583E-02_wp, &
      & 5.90543543298967E-04_wp, 1.89070265565719E-02_wp,-2.18204966861747E-02_wp, &
      &-1.64116449632240E-02_wp,-3.50049447183280E-02_wp, 1.92048599024866E-03_wp, &
      & 2.73294147835100E-04_wp, 2.40805559729263E-04_wp, 3.98704788703111E-03_wp, &
      & 3.83007723350443E-03_wp, 2.17842958083843E-02_wp, 2.23980508202744E-02_wp, &
      & 1.05297118737512E-02_wp, 3.20167044114457E-03_wp,-1.26672314141543E-02_wp, &
      & 6.87221148902045E-02_wp,-5.52582904665705E-03_wp,-2.09804044214836E-03_wp, &
      &-9.70603235678339E-02_wp,-5.30554321998536E-02_wp,-3.29569302511059E-02_wp, &
      & 1.21189122316512E-02_wp,-4.33726604714417E-02_wp, 3.72566964773942E-03_wp, &
      &-5.90050199661660E-02_wp, 4.41302193900175E-02_wp,-2.29461193380960E-02_wp, &
      &-1.19859393613802E-02_wp,-1.39791945191686E-02_wp,-5.05076972762294E-02_wp, &
      &-5.15532684544786E-03_wp,-6.40001479013474E-04_wp, 1.57912740212245E-03_wp, &
      &-3.30674201547255E-03_wp,-2.96462717477253E-03_wp,-4.74581076068379E-02_wp, &
      & 1.39890422776412E-02_wp, 8.33818371867964E-03_wp, 9.00709592745119E-01_wp, &
      &-2.32277295829346E-02_wp,-5.91310561328952E-04_wp,-1.03765797970833E-02_wp, &
      &-1.37570262673529E-02_wp, 1.79170712797340E-02_wp, 1.33699745251540E-02_wp, &
      & 1.12806043159415E-02_wp,-9.44493581999975E-03_wp,-1.90246054652711E-03_wp, &
      &-1.78419681380527E-03_wp,-9.11286541924675E-04_wp,-5.17440346512460E-04_wp, &
      & 2.69307620534820E-03_wp, 1.58235914557377E-03_wp,-1.32239449246803E-03_wp, &
      &-8.58573131874772E-04_wp,-6.83652635175080E-04_wp, 3.34046150737528E-04_wp, &
      & 1.72311948156705E-03_wp,-2.66981711043062E-03_wp, 2.36297575249206E-04_wp, &
      & 1.14598547458529E-03_wp, 1.09039808896316E-03_wp, 5.62332461341726E-04_wp, &
      & 2.35295984835294E-03_wp,-6.66975469219293E-05_wp, 2.32152383555904E-03_wp, &
      & 1.02640693286939E-04_wp, 1.83884922532074E-03_wp,-2.26312525914848E-03_wp, &
      & 1.08781582530057E-03_wp,-6.98434276122301E-04_wp,-7.01721881174691E-04_wp, &
      & 5.43204781193947E-03_wp, 2.01568567478042E-04_wp,-5.97084418745616E-05_wp, &
      &-8.60224070526115E-05_wp, 1.03484025349429E-04_wp, 1.86342724132078E-04_wp, &
      & 3.77653932652044E-03_wp, 1.70641049794269E-02_wp,-3.41304494398283E-04_wp, &
      &-2.32277295829346E-02_wp, 9.71954797249471E-04_wp,-2.66463312336126E-04_wp, &
      & 5.36546546854042E-04_wp, 3.78261896705893E-04_wp,-2.20377842351566E-05_wp, &
      &-1.08352826422062E-03_wp,-1.82074210886211E-03_wp,-1.01687330465836E-03_wp, &
      & 1.84168810429607E-03_wp, 1.24242288971349E-04_wp, 2.60909585417293E-05_wp, &
      & 9.82423307905643E-05_wp,-1.33251711431608E-05_wp,-7.80121783928019E-05_wp, &
      & 9.23763188397653E-04_wp, 6.50886246444007E-04_wp,-3.67543219007418E-06_wp, &
      &-7.19136030754418E-04_wp,-2.33142012668874E-03_wp, 3.40407076461401E-04_wp, &
      &-1.94458556002641E-03_wp,-1.69273880680804E-04_wp,-3.59947499872617E-04_wp, &
      & 1.78951892334922E-03_wp,-7.98432581416592E-04_wp,-3.13907873268688E-03_wp, &
      &-7.19828442743159E-04_wp,-1.47429100000553E-05_wp,-5.27133217598405E-04_wp, &
      & 2.64354771088751E-03_wp, 8.21284940311572E-04_wp, 2.29788467999063E-03_wp, &
      & 1.08687409943675E-04_wp,-9.57154239876852E-05_wp, 8.47489533143351E-05_wp, &
      & 4.30725140314870E-05_wp,-8.26444834398727E-05_wp,-1.06917681623253E-06_wp, &
      &-1.01495495974086E-04_wp,-2.82491840284670E-03_wp,-1.41030186631462E-02_wp, &
      &-2.29885786612005E-02_wp,-5.91310561328952E-04_wp,-2.66463312336126E-04_wp, &
      & 8.50416540478395E-04_wp, 2.20268366240725E-04_wp,-4.57895934640666E-04_wp, &
      &-3.97935121573021E-04_wp,-1.96526810803974E-03_wp, 3.13903002693450E-03_wp, &
      & 3.85187428686917E-03_wp, 2.11734626341255E-03_wp,-1.34264453635969E-04_wp, &
      &-2.51838077167750E-05_wp,-5.85408576229357E-05_wp,-2.31098244614495E-04_wp, &
      &-1.74472472563477E-04_wp,-5.87060101040320E-04_wp,-1.23589021153629E-04_wp, &
      &-5.64274993106786E-04_wp, 4.42987781466214E-05_wp, 1.41048483944506E-04_wp, &
      &-2.12977682744976E-03_wp,-1.28072235876254E-03_wp, 1.52074117567880E-03_wp, &
      &-3.60824594442789E-04_wp, 5.49306638141359E-04_wp, 2.29759184132251E-03_wp, &
      &-3.25713971009019E-03_wp, 1.28190858547835E-03_wp,-1.17086414461730E-03_wp, &
      & 3.10277119084791E-04_wp, 1.76499596098071E-04_wp, 8.63411726541678E-04_wp, &
      & 3.50868101513821E-04_wp,-1.44491706265718E-03_wp, 2.35591709642209E-03_wp, &
      & 1.28915036185723E-04_wp,-9.13448849034343E-05_wp,-1.22481720956570E-04_wp, &
      &-1.97825617795629E-07_wp, 4.60297941197738E-05_wp, 1.53601953767155E-03_wp, &
      & 1.27499204873158E-02_wp,-1.59420543751472E-02_wp,-1.03765797970833E-02_wp, &
      & 5.36546546854042E-04_wp, 2.20268366240725E-04_wp, 6.09265712051661E-04_wp, &
      &-1.47519943592145E-04_wp, 1.22238728619821E-04_wp,-1.25095952083766E-03_wp, &
      & 2.27846314193697E-04_wp, 1.88274865293780E-04_wp, 2.05809775537036E-03_wp, &
      & 1.14606708087605E-05_wp,-2.53947944812789E-05_wp, 7.94171880809830E-05_wp, &
      &-5.67477088128098E-05_wp,-1.57956030127115E-04_wp,-6.06266471532816E-04_wp, &
      &-7.14602097204805E-04_wp,-4.99490489386992E-05_wp, 2.30642331218448E-04_wp, &
      & 1.23754087718086E-03_wp,-3.25448102001268E-04_wp, 1.67826970775080E-03_wp, &
      &-6.69143343856796E-04_wp, 1.85106237168900E-03_wp, 1.18354420618126E-04_wp, &
      &-6.70858028794763E-04_wp, 4.31017066296889E-03_wp, 7.83301665643090E-04_wp, &
      & 2.24601344723718E-03_wp, 1.78806056495109E-03_wp,-2.64863626414187E-03_wp, &
      & 3.71234417787815E-04_wp,-9.19684784734447E-04_wp, 4.54075284812082E-04_wp, &
      & 2.63084615326078E-03_wp, 5.96716941142845E-05_wp, 5.35251386105152E-05_wp, &
      & 5.33348705135931E-05_wp, 1.08987968326367E-04_wp, 1.31730529441688E-04_wp, &
      & 1.88938374714825E-03_wp, 3.86236907846822E-04_wp, 1.67323852681583E-02_wp, &
      &-1.37570262673529E-02_wp, 3.78261896705893E-04_wp,-4.57895934640666E-04_wp, &
      &-1.47519943592145E-04_wp, 5.67070951883082E-04_wp,-2.66365841731747E-04_wp, &
      &-2.17685295273928E-05_wp,-2.21004109845667E-03_wp,-1.02849359837745E-03_wp, &
      &-9.02442359022939E-04_wp, 1.12455065949128E-04_wp, 6.84526285745388E-05_wp, &
      &-1.17897739761924E-05_wp, 2.91415636451971E-05_wp, 1.16342256034976E-04_wp, &
      &-5.25957119780333E-04_wp, 1.61848868701550E-04_wp,-1.70869611258782E-04_wp, &
      & 6.80365098404692E-04_wp, 1.32235237491202E-03_wp,-6.09023766965693E-04_wp, &
      &-2.64680965053990E-04_wp, 1.75734486260590E-03_wp,-2.29662810010202E-03_wp, &
      &-2.60744135287044E-03_wp, 2.19545673464568E-03_wp,-1.47887531130849E-03_wp, &
      & 1.95937371861310E-04_wp,-2.00559116494408E-03_wp,-1.42996134361120E-03_wp, &
      &-2.83279606761995E-05_wp,-2.26885600405926E-03_wp,-1.01534697330281E-03_wp, &
      &-1.98967421502637E-03_wp,-2.80062748714972E-03_wp,-1.61436989305193E-04_wp, &
      &-1.48760158878772E-04_wp,-3.16133131415677E-06_wp,-1.54513083585891E-04_wp, &
      &-5.02008128433043E-05_wp, 1.39873724097239E-03_wp, 2.29766788919789E-02_wp, &
      & 5.90543543298967E-04_wp, 1.79170712797340E-02_wp,-2.20377842351566E-05_wp, &
      &-3.97935121573021E-04_wp, 1.22238728619821E-04_wp,-2.66365841731747E-04_wp, &
      & 9.76918208350585E-04_wp, 1.37310447814700E-03_wp,-3.85050874153481E-04_wp, &
      &-4.29155366732576E-03_wp,-8.98951630115442E-04_wp, 7.11466472205285E-06_wp, &
      &-7.22055400547691E-05_wp, 1.11585039122059E-04_wp, 2.53643008088297E-04_wp, &
      & 2.63123070309016E-05_wp,-9.44887800903170E-04_wp, 2.47681579669712E-03_wp, &
      &-3.50216085250586E-03_wp,-8.84877256219468E-03_wp,-1.35903661755520E-02_wp, &
      &-5.22282280974187E-03_wp,-9.91183630949829E-03_wp, 2.14421083889722E-02_wp, &
      & 2.20191489255154E-02_wp,-1.84935397407614E-04_wp, 9.16685372743052E-04_wp, &
      &-1.84087838689988E-04_wp,-2.19901739923401E-04_wp,-1.12586025628635E-03_wp, &
      & 9.38729389402143E-05_wp,-2.41011357456463E-04_wp,-4.67052712125803E-04_wp, &
      &-1.79249305334720E-02_wp, 1.69947916255744E-03_wp, 1.34259756367200E-02_wp, &
      & 9.05992698760134E-04_wp,-3.65628894307319E-04_wp, 6.83482385411022E-04_wp, &
      & 9.18499763895862E-04_wp, 2.19806549831163E-03_wp,-3.03607735046782E-03_wp, &
      & 9.93178929003970E-03_wp, 1.89070265565719E-02_wp, 1.33699745251540E-02_wp, &
      &-1.08352826422062E-03_wp,-1.96526810803974E-03_wp,-1.25095952083766E-03_wp, &
      &-2.17685295273928E-05_wp, 1.37310447814700E-03_wp, 9.79754791069694E-01_wp, &
      &-2.23989504876204E-02_wp, 3.40636152248749E-02_wp, 8.18979621687276E-02_wp, &
      & 2.09012487006927E-03_wp, 8.40503347739348E-04_wp, 1.56038857676965E-03_wp, &
      &-3.26760910350976E-03_wp,-3.45211447044762E-03_wp, 9.71328579477840E-03_wp, &
      & 1.51505155319784E-02_wp, 4.54245265617874E-03_wp, 1.35788864585994E-02_wp, &
      & 7.65937840551549E-02_wp, 3.75844989308166E-02_wp,-5.57175718616085E-03_wp, &
      & 3.28723251997694E-02_wp, 7.33216018913181E-02_wp, 4.00823788730297E-02_wp, &
      & 2.57126867846738E-02_wp, 1.41827982066621E-02_wp,-2.48692793280382E-02_wp, &
      & 2.39815493205939E-02_wp, 4.82654826901009E-02_wp, 5.49357483690476E-02_wp, &
      &-2.18177810703489E-02_wp,-4.44215180079783E-02_wp,-1.66229817006404E-02_wp, &
      & 9.67852731635338E-03_wp,-1.86454186027253E-03_wp,-3.34642237400698E-03_wp, &
      & 1.43421135623700E-03_wp, 2.49064579621402E-04_wp, 5.36800952916396E-03_wp, &
      & 1.59310863120059E-02_wp,-2.68014119418027E-03_wp,-2.18204966861747E-02_wp, &
      & 1.12806043159415E-02_wp,-1.82074210886211E-03_wp, 3.13903002693450E-03_wp, &
      & 2.27846314193697E-04_wp,-2.21004109845667E-03_wp,-3.85050874153481E-04_wp, &
      &-2.23989504876204E-02_wp, 9.07216790441336E-01_wp, 4.71099439338665E-03_wp, &
      & 6.98669050592126E-03_wp,-2.76452181520636E-02_wp,-1.32617256121538E-02_wp, &
      &-5.84858535890091E-03_wp,-9.24621410530014E-04_wp,-8.64305835005504E-03_wp, &
      &-1.32264514234976E-02_wp, 4.13303227392449E-03_wp, 1.19055309005284E-02_wp, &
      &-2.03601998000453E-02_wp, 3.28686383454786E-02_wp,-1.69469206157941E-02_wp, &
      & 8.67153633545441E-02_wp, 4.78972976105153E-02_wp,-7.19861577004959E-02_wp, &
      & 5.07256361601547E-02_wp, 1.14966740670106E-02_wp,-6.62740370087469E-02_wp, &
      & 6.56157086460467E-03_wp, 2.95058137735713E-02_wp,-1.14992339504637E-02_wp, &
      &-2.23923554739493E-02_wp,-4.85160230141416E-03_wp,-4.45142330247158E-03_wp, &
      & 2.23035839446687E-02_wp,-1.30717775428052E-02_wp,-1.56221941951145E-03_wp, &
      &-3.20723614998436E-04_wp, 1.29954197474567E-03_wp,-3.01515131463446E-04_wp, &
      &-1.11865135939340E-04_wp, 1.87885921811414E-02_wp,-3.97830522408287E-02_wp, &
      &-1.64116449632240E-02_wp,-9.44493581999975E-03_wp,-1.01687330465836E-03_wp, &
      & 3.85187428686917E-03_wp, 1.88274865293780E-04_wp,-1.02849359837745E-03_wp, &
      &-4.29155366732576E-03_wp, 3.40636152248749E-02_wp, 4.71099439338665E-03_wp, &
      & 9.03653999747069E-01_wp,-1.07477493791768E-02_wp,-6.79035487611753E-04_wp, &
      & 7.65027028291042E-03_wp,-1.49468871725635E-02_wp,-2.63446016222559E-02_wp, &
      &-1.04417300492389E-03_wp,-4.00506954182355E-02_wp, 1.19937863232238E-02_wp, &
      &-1.89552860619114E-02_wp,-2.27788852326930E-02_wp, 6.44067474820537E-02_wp, &
      & 4.36814617080786E-02_wp, 2.90410287476276E-03_wp,-1.04712017981954E-01_wp, &
      &-4.69373915752494E-02_wp, 6.01728319805672E-02_wp, 3.88665260783668E-02_wp, &
      & 2.39879765902584E-02_wp, 6.03966283303314E-02_wp, 2.61728970657019E-02_wp, &
      &-6.93910045600238E-02_wp,-3.95054115593926E-02_wp, 2.44761467320601E-03_wp, &
      & 3.22808677054423E-02_wp,-9.27699182882334E-03_wp,-3.47522554510675E-03_wp, &
      &-8.01779151303163E-04_wp,-3.02716907708290E-04_wp,-6.84034941755676E-04_wp, &
      &-1.87144089306515E-03_wp,-2.11155401443376E-03_wp, 5.08467535689386E-03_wp, &
      &-1.26975405204845E-02_wp,-3.50049447183280E-02_wp,-1.90246054652711E-03_wp, &
      & 1.84168810429607E-03_wp, 2.11734626341255E-03_wp, 2.05809775537036E-03_wp, &
      &-9.02442359022939E-04_wp,-8.98951630115442E-04_wp, 8.18979621687276E-02_wp, &
      & 6.98669050592126E-03_wp,-1.07477493791768E-02_wp, 8.78220908101747E-01_wp, &
      & 8.78198962890551E-03_wp, 7.29973779520088E-04_wp, 1.48135170964552E-02_wp, &
      &-1.41265088250659E-02_wp,-2.63165868327368E-02_wp,-7.30975427202840E-04_wp, &
      &-4.00974483671558E-04_wp,-4.21614601949817E-04_wp,-6.78210598851248E-04_wp, &
      &-1.00233070466989E-03_wp,-7.35276572443060E-04_wp, 3.78053496299328E-04_wp, &
      &-1.59454356614859E-03_wp,-2.10917504939881E-03_wp,-1.08074030232713E-03_wp, &
      &-9.20983705270224E-04_wp, 1.05566543856893E-03_wp, 1.66757073509023E-03_wp, &
      & 5.32051158886537E-04_wp,-1.97145781885208E-03_wp,-2.59058745650320E-03_wp, &
      & 2.10650456051256E-03_wp, 5.51511024939449E-03_wp, 4.79918109489262E-04_wp, &
      &-2.41917211488738E-03_wp, 5.29754230196242E-05_wp, 1.66642160459250E-04_wp, &
      &-6.77514813142387E-05_wp,-4.15809966344942E-05_wp,-3.09423947555397E-04_wp, &
      &-1.69023816854932E-03_wp, 1.43642412751507E-03_wp, 1.92048599024866E-03_wp, &
      &-1.78419681380527E-03_wp, 1.24242288971349E-04_wp,-1.34264453635969E-04_wp, &
      & 1.14606708087605E-05_wp, 1.12455065949128E-04_wp, 7.11466472205285E-06_wp, &
      & 2.09012487006927E-03_wp,-2.76452181520636E-02_wp,-6.79035487611753E-04_wp, &
      & 8.78198962890551E-03_wp, 9.72829001384354E-04_wp, 4.31613119157028E-04_wp, &
      & 3.39329532587659E-04_wp,-1.10362577227386E-04_wp, 8.20860388579246E-06_wp, &
      &-2.86085329176629E-04_wp,-1.44661987626867E-04_wp,-8.16969303097783E-05_wp, &
      &-4.08503618156882E-04_wp,-6.88819688386836E-04_wp,-6.06690581189437E-04_wp, &
      & 9.38622543054151E-04_wp,-7.29860323367932E-05_wp,-1.45840439450965E-03_wp, &
      &-5.41519705134412E-04_wp,-1.88651521181885E-03_wp, 3.60304872454716E-03_wp, &
      & 4.66111281355355E-04_wp, 4.14351970336048E-03_wp, 3.56893172299223E-05_wp, &
      &-2.26084037290193E-03_wp, 5.13537734104285E-04_wp, 2.23194762525990E-03_wp, &
      & 1.40000483985680E-04_wp,-1.02717741613644E-03_wp,-1.67720095140052E-05_wp, &
      & 1.29518345210917E-04_wp, 3.57114046682632E-05_wp, 1.50166180257202E-05_wp, &
      &-1.21929197699361E-04_wp,-6.52194452941343E-04_wp,-5.15623780963257E-04_wp, &
      & 2.73294147835100E-04_wp,-9.11286541924675E-04_wp, 2.60909585417293E-05_wp, &
      &-2.51838077167750E-05_wp,-2.53947944812789E-05_wp, 6.84526285745388E-05_wp, &
      &-7.22055400547691E-05_wp, 8.40503347739348E-04_wp,-1.32617256121538E-02_wp, &
      & 7.65027028291042E-03_wp, 7.29973779520088E-04_wp, 4.31613119157028E-04_wp, &
      & 3.07330023331222E-04_wp,-6.96308160685538E-05_wp,-2.72899832827896E-04_wp, &
      & 1.35850391213250E-04_wp,-5.83829373771101E-04_wp,-2.20905355403101E-05_wp, &
      &-5.05453807145344E-04_wp,-1.59367877771823E-04_wp, 3.53857312558554E-04_wp, &
      & 6.06034207884786E-04_wp,-1.22901953593365E-03_wp,-2.39582887843716E-03_wp, &
      & 1.08123286917808E-04_wp, 1.74128670059472E-04_wp, 2.10285595596077E-03_wp, &
      &-3.54139161540637E-03_wp, 1.29891259805378E-03_wp,-5.00522584023764E-03_wp, &
      &-2.34151302886779E-03_wp, 4.94952180347453E-04_wp, 8.26572953278356E-04_wp, &
      & 2.14969208354629E-03_wp,-1.14187690595404E-03_wp,-7.37519467376673E-04_wp, &
      & 6.88550780570842E-05_wp,-6.06645445284161E-05_wp,-1.37991308366808E-04_wp, &
      &-9.17945277868815E-05_wp,-1.22852682598191E-04_wp,-1.24575925065902E-03_wp, &
      & 1.97029419194627E-03_wp, 2.40805559729263E-04_wp,-5.17440346512460E-04_wp, &
      & 9.82423307905643E-05_wp,-5.85408576229357E-05_wp, 7.94171880809830E-05_wp, &
      &-1.17897739761924E-05_wp, 1.11585039122059E-04_wp, 1.56038857676965E-03_wp, &
      &-5.84858535890091E-03_wp,-1.49468871725635E-02_wp, 1.48135170964552E-02_wp, &
      & 3.39329532587659E-04_wp,-6.96308160685538E-05_wp, 5.91760718214581E-04_wp, &
      & 2.52912255874505E-04_wp,-4.23264327989375E-04_wp, 1.06473924707788E-03_wp, &
      &-4.45102163842215E-04_wp, 3.80106625787291E-05_wp, 1.00614012061151E-03_wp, &
      &-1.40694498188145E-03_wp, 1.10660155876927E-04_wp,-2.36862461809501E-03_wp, &
      & 4.78040749240586E-04_wp, 2.31640389941183E-03_wp,-2.34898043486087E-03_wp, &
      & 5.97064769894395E-04_wp,-2.95929713093146E-03_wp,-8.75032383974094E-04_wp, &
      &-6.01050825620884E-03_wp, 4.72722357628950E-04_wp, 2.52712586784912E-03_wp, &
      & 7.14219643152562E-04_wp, 2.26334563113351E-04_wp, 2.76159075758330E-04_wp, &
      & 2.37177085330210E-03_wp, 1.49459436389064E-04_wp,-3.58157141410844E-05_wp, &
      &-1.10652382973498E-04_wp, 3.21446802816235E-05_wp, 4.05393279050851E-05_wp, &
      &-1.84891013912983E-03_wp, 3.39979714316641E-03_wp, 3.98704788703111E-03_wp, &
      & 2.69307620534820E-03_wp,-1.33251711431608E-05_wp,-2.31098244614495E-04_wp, &
      &-5.67477088128098E-05_wp, 2.91415636451971E-05_wp, 2.53643008088297E-04_wp, &
      &-3.26760910350976E-03_wp,-9.24621410530014E-04_wp,-2.63446016222559E-02_wp, &
      &-1.41265088250659E-02_wp,-1.10362577227386E-04_wp,-2.72899832827896E-04_wp, &
      & 2.52912255874505E-04_wp, 1.08310028976737E-03_wp, 4.38482863597593E-04_wp, &
      & 1.20864032141500E-03_wp,-5.21718641127864E-04_wp, 4.70419178269535E-04_wp, &
      & 6.42501172500844E-04_wp,-2.36540545013052E-03_wp,-1.13219077295543E-03_wp, &
      & 3.42392873712068E-05_wp, 2.59547000700204E-03_wp, 5.16614467009417E-04_wp, &
      &-2.72927683099535E-03_wp,-3.32760114542134E-03_wp, 4.28984898113645E-03_wp, &
      &-1.61471980448914E-03_wp, 3.53071480912169E-03_wp, 2.70282186006250E-03_wp, &
      &-5.02041178188871E-04_wp, 1.14929051465422E-03_wp,-2.69348432555279E-04_wp, &
      & 1.54213587965840E-03_wp, 1.59044196104408E-03_wp, 3.48745958380269E-05_wp, &
      & 1.47884993200396E-04_wp, 8.18569440156331E-05_wp, 1.29877997003160E-04_wp, &
      & 4.47136744942552E-05_wp,-3.95093353420261E-04_wp, 1.57372281682921E-04_wp, &
      & 3.83007723350443E-03_wp, 1.58235914557377E-03_wp,-7.80121783928019E-05_wp, &
      &-1.74472472563477E-04_wp,-1.57956030127115E-04_wp, 1.16342256034976E-04_wp, &
      & 2.63123070309016E-05_wp,-3.45211447044762E-03_wp,-8.64305835005504E-03_wp, &
      &-1.04417300492389E-03_wp,-2.63165868327368E-02_wp, 8.20860388579246E-06_wp, &
      & 1.35850391213250E-04_wp,-4.23264327989375E-04_wp, 4.38482863597593E-04_wp, &
      & 9.40813877254780E-04_wp, 7.11952239175513E-03_wp,-4.48759433563561E-04_wp, &
      & 7.04758808977056E-04_wp,-2.94624945675376E-04_wp,-2.74166581356306E-03_wp, &
      & 6.84153465016098E-03_wp,-5.31279761822322E-03_wp, 4.85658416213561E-03_wp, &
      & 9.50019121810167E-04_wp,-2.43299151738308E-03_wp,-9.06679526463874E-03_wp, &
      &-3.80607531989858E-04_wp,-2.37874651841337E-03_wp,-3.06390494210790E-04_wp, &
      & 3.36294433121553E-03_wp, 8.99503672589712E-03_wp,-1.11688855965218E-03_wp, &
      & 3.20046292353588E-02_wp, 2.32606695048053E-02_wp, 2.10966188523540E-02_wp, &
      & 1.14006433275913E-03_wp, 1.14147229303145E-03_wp,-3.02272998088046E-04_wp, &
      & 7.61032562708332E-04_wp,-4.79144450852178E-04_wp,-1.30082860861281E-03_wp, &
      &-3.51696354973964E-02_wp,-1.39407590739093E-02_wp, 2.24022500427106E-02_wp, &
      &-1.30001097214208E-03_wp, 9.03676498870593E-04_wp,-5.83171540616389E-04_wp, &
      &-5.86376985234641E-04_wp,-5.20172683811391E-04_wp,-1.20825923162660E-03_wp, &
      & 1.01661826271254E-02_wp,-1.38165545028928E-02_wp,-4.08742499801280E-02_wp, &
      &-7.18067691148896E-04_wp,-2.63641797653755E-04_wp,-5.93966433083412E-04_wp, &
      & 1.02979125404226E-03_wp, 1.20885516832012E-03_wp,-4.48759433563561E-04_wp, &
      & 2.28687129134323E-03_wp, 5.41396449946994E-04_wp,-2.59430166947037E-04_wp, &
      &-2.64919936937235E-03_wp, 5.77530559995348E-04_wp,-5.26976013457427E-05_wp, &
      &-1.70360657293966E-03_wp,-2.12815610709220E-04_wp, 1.06766405478892E-03_wp, &
      & 1.23297231618413E-03_wp,-1.24027234557119E-03_wp,-2.62975652749730E-03_wp, &
      & 6.81819221360930E-04_wp,-2.77984913399897E-03_wp, 6.15223438265065E-04_wp, &
      & 7.17489050802549E-03_wp,-1.15284447884868E-02_wp,-1.66995151518806E-02_wp, &
      &-1.96285390992293E-02_wp,-8.09390428732359E-04_wp,-6.27932511020678E-04_wp, &
      & 1.21450766756529E-04_wp,-6.28894926244606E-04_wp,-7.62484937973928E-05_wp, &
      &-8.00507020680877E-03_wp,-1.18379012605914E-02_wp,-1.74855485990935E-02_wp, &
      & 2.22060923985087E-02_wp,-8.41884382085824E-04_wp, 6.54018323036404E-04_wp, &
      &-9.79735922630597E-05_wp,-7.38014920537391E-04_wp, 1.67408952638943E-04_wp, &
      & 2.55413860729764E-03_wp, 1.44077180987680E-02_wp, 4.81975145073898E-03_wp, &
      & 1.15402494798059E-02_wp,-3.91472860045389E-04_wp,-1.85447400076914E-04_wp, &
      & 2.69372961745180E-05_wp,-3.83960703078053E-04_wp,-5.52056921452567E-04_wp, &
      & 7.04758808977056E-04_wp, 5.41396449946994E-04_wp, 1.62716519681097E-03_wp, &
      & 3.61628176171920E-04_wp,-1.63409127414125E-03_wp,-1.07250498654219E-03_wp, &
      & 3.31847658400096E-04_wp,-6.83598102195532E-04_wp, 5.75709000782720E-04_wp, &
      &-1.94547758751811E-04_wp,-6.02151183304568E-04_wp, 1.49179954780989E-03_wp, &
      &-2.93152544359381E-03_wp, 1.04205333697116E-03_wp, 1.48481072357888E-04_wp, &
      & 6.08528200623539E-04_wp, 4.57690583399911E-03_wp,-1.88737263402893E-02_wp, &
      & 8.75664181962324E-03_wp,-1.26486079338153E-02_wp,-7.10043290946803E-04_wp, &
      & 2.89620397193752E-06_wp, 6.06131016453163E-04_wp, 2.82293389952361E-06_wp, &
      & 2.86985578991511E-04_wp,-3.80262208043012E-03_wp,-1.56986199671825E-02_wp, &
      & 9.95289435445029E-03_wp, 1.08109379474349E-02_wp,-6.58426823850518E-04_wp, &
      &-2.72580344710504E-05_wp,-5.84203711930999E-04_wp, 1.67615675302901E-05_wp, &
      &-1.95672888720868E-04_wp,-3.63428057490533E-03_wp, 5.13544987508667E-03_wp, &
      & 1.02100793239369E-02_wp,-1.81555166175934E-02_wp,-3.91264029795358E-04_wp, &
      & 9.88365504780387E-06_wp,-5.75792938733272E-04_wp,-4.71112289725868E-05_wp, &
      & 5.46257757174947E-04_wp,-2.94624945675376E-04_wp,-2.59430166947037E-04_wp, &
      & 3.61628176171920E-04_wp, 2.43622668168914E-03_wp,-1.06123973390260E-03_wp, &
      &-1.65324946777566E-03_wp,-7.25296306806351E-05_wp, 1.80304070856658E-03_wp, &
      & 2.13517632377420E-03_wp,-1.81460080347714E-03_wp, 2.37543992326715E-03_wp, &
      & 6.76725101117992E-04_wp,-1.81030095889485E-03_wp,-1.71651898481651E-03_wp, &
      & 2.34703941229088E-03_wp,-2.97757430265402E-04_wp, 4.75605525497978E-03_wp, &
      &-1.99507322125130E-02_wp,-1.13575027794173E-02_wp, 4.49366063601718E-03_wp, &
      &-2.46172687351389E-04_wp,-6.32133131699997E-04_wp, 8.27214390433451E-05_wp, &
      &-1.20209486448711E-04_wp, 7.01975281929006E-04_wp, 5.46887348151783E-03_wp, &
      & 2.07465098916450E-02_wp, 1.12870324819613E-02_wp, 2.62360032993394E-03_wp, &
      & 3.30170019630099E-04_wp,-6.85792113938533E-04_wp, 6.13175091245618E-05_wp, &
      & 2.02457181360978E-04_wp, 6.65644753860844E-04_wp,-9.08181098068503E-03_wp, &
      & 1.33513181820170E-02_wp,-1.99268786228122E-02_wp,-2.34217202320808E-02_wp, &
      &-6.69813653944928E-04_wp,-4.38712775661282E-04_wp,-1.20006019568062E-04_wp, &
      & 1.02053886193372E-03_wp, 5.94007597988317E-04_wp,-2.74166581356306E-03_wp, &
      &-2.64919936937235E-03_wp,-1.63409127414125E-03_wp,-1.06123973390260E-03_wp, &
      & 4.26660636479651E-02_wp, 5.48867639267497E-03_wp, 2.54176254686534E-03_wp, &
      & 7.36832296027678E-03_wp,-4.03975654473191E-04_wp,-3.26942203750127E-04_wp, &
      & 3.91589546660145E-03_wp, 4.64342996080958E-03_wp, 1.33992841815537E-02_wp, &
      & 7.49785230431801E-04_wp, 8.17240065739435E-04_wp,-3.19238355446981E-04_wp, &
      & 2.23646244137205E-02_wp, 6.59548065857667E-02_wp, 8.68754738499316E-02_wp, &
      & 1.25040456174467E-03_wp, 8.43040832431269E-04_wp, 2.57139294590359E-03_wp, &
      & 4.69721320177459E-04_wp, 1.34881283254879E-03_wp,-1.20147515573048E-03_wp, &
      &-2.39915293849405E-02_wp, 7.91814684776250E-02_wp, 6.88501216429399E-02_wp, &
      &-1.27801136787652E-02_wp, 1.67318715110993E-03_wp,-2.26525736225408E-03_wp, &
      & 1.03312589060382E-04_wp, 1.25511859423495E-03_wp, 1.23392237280349E-03_wp, &
      &-1.35620731019085E-02_wp, 7.56618991833905E-02_wp, 3.15084854540211E-02_wp, &
      & 6.45053437127237E-02_wp,-9.38359904069108E-04_wp,-6.34043613155897E-04_wp, &
      & 3.14862311605270E-04_wp,-1.37787803489033E-03_wp,-2.24225399100554E-03_wp, &
      & 6.84153465016098E-03_wp, 5.77530559995348E-04_wp,-1.07250498654219E-03_wp, &
      &-1.65324946777566E-03_wp, 5.48867639267497E-03_wp, 3.04783663710111E-02_wp, &
      & 5.04606442096020E-03_wp,-4.66371208138313E-03_wp,-7.19052804954172E-03_wp, &
      & 3.26384647666875E-03_wp,-1.16648070070560E-02_wp,-1.74839240093656E-03_wp, &
      & 1.30572623508845E-03_wp,-5.03897420550892E-03_wp,-8.96763548845411E-04_wp, &
      & 9.84759792548497E-03_wp, 2.07224740193060E-02_wp, 6.70688097195503E-02_wp, &
      &-1.96459370841130E-02_wp, 7.81809041418509E-02_wp, 2.84420597616065E-03_wp, &
      & 3.81213248123500E-04_wp,-1.93389141426929E-03_wp, 4.82420477979741E-04_wp, &
      &-5.32810176688762E-04_wp, 1.83968711176355E-02_wp,-6.46585082512025E-02_wp, &
      & 3.76360397992066E-02_wp, 6.76726646968307E-02_wp,-2.56333530435510E-03_wp, &
      & 3.57609160918523E-04_wp,-2.02312098196774E-03_wp,-3.47311396599134E-04_wp, &
      &-5.79280405069415E-04_wp,-5.21056540703168E-03_wp, 3.66511396790259E-02_wp, &
      &-1.61959501492248E-02_wp, 4.22610487435226E-02_wp,-7.23150166636224E-04_wp, &
      &-6.96300015777415E-04_wp, 6.92000783727072E-04_wp, 2.14568368478460E-04_wp, &
      &-1.19565402301603E-03_wp,-5.31279761822322E-03_wp,-5.26976013457427E-05_wp, &
      & 3.31847658400096E-04_wp,-7.25296306806351E-05_wp, 2.54176254686534E-03_wp, &
      & 5.04606442096020E-03_wp, 2.84679587556318E-02_wp, 3.54894557191521E-03_wp, &
      &-7.72046463324211E-04_wp, 3.39812061949517E-03_wp, 1.80295732067388E-03_wp, &
      &-1.72699146503798E-03_wp, 2.08068962045090E-04_wp,-9.57025303579991E-04_wp, &
      &-6.72245961471031E-04_wp,-1.25000483023965E-02_wp,-7.48190813025657E-03_wp, &
      & 1.50194621101850E-02_wp,-9.51746257217108E-02_wp, 1.07105823451796E-02_wp, &
      & 2.74276798196848E-04_wp,-2.28477526111477E-03_wp,-1.54670546994761E-03_wp, &
      &-1.52509696471355E-03_wp,-8.34604782021266E-05_wp,-9.70631584445687E-03_wp, &
      & 5.89106335176432E-04_wp, 8.54634077370308E-02_wp,-6.23156571113132E-03_wp, &
      & 2.39237960359503E-04_wp,-1.86157192226641E-03_wp,-1.22914115729558E-03_wp, &
      & 1.62380339070375E-03_wp,-2.57899109827892E-04_wp,-9.88014923090613E-03_wp, &
      &-6.25144814937098E-03_wp, 8.60612099656805E-02_wp, 3.14099372304278E-03_wp, &
      & 3.78532516914335E-04_wp, 9.38872653962537E-04_wp,-1.19990637124712E-03_wp, &
      &-2.29785457537160E-03_wp, 4.74724238071550E-05_wp, 4.85658416213561E-03_wp, &
      &-1.70360657293966E-03_wp,-6.83598102195532E-04_wp, 1.80304070856658E-03_wp, &
      & 7.36832296027678E-03_wp,-4.66371208138313E-03_wp, 3.54894557191521E-03_wp, &
      & 3.34393781938852E-02_wp, 4.73561694057961E-03_wp,-9.20066557038704E-03_wp, &
      &-4.78340485526786E-03_wp,-5.05261443662806E-03_wp, 6.84791172570740E-04_wp, &
      & 2.03099825524059E-03_wp, 9.45917682510810E-03_wp, 3.78879388099217E-03_wp, &
      & 1.40822121394522E-02_wp, 7.85903417261118E-02_wp,-1.34573967461195E-02_wp, &
      & 5.53985700082034E-03_wp, 1.51665855913071E-03_wp, 4.92005695111929E-04_wp, &
      &-1.33814279600240E-03_wp, 9.60493988837844E-06_wp,-1.49543254126618E-03_wp, &
      &-1.19162210963848E-02_wp, 7.77031024200890E-02_wp,-3.00536616447891E-02_wp, &
      &-2.00732274111886E-03_wp, 1.10329669950783E-03_wp,-1.68670894926453E-04_wp, &
      & 1.46237940335675E-03_wp,-6.46986175423066E-04_wp, 1.68257298559314E-03_wp, &
      & 2.13600011676214E-02_wp, 3.25281161052191E-02_wp, 4.67381503377974E-02_wp, &
      &-1.03858523700383E-01_wp,-1.50740974288112E-03_wp, 3.82651994863863E-05_wp, &
      &-2.40313487407486E-03_wp, 3.71631530997041E-04_wp, 2.60157115829100E-03_wp, &
      & 9.50019121810167E-04_wp,-2.12815610709220E-04_wp, 5.75709000782720E-04_wp, &
      & 2.13517632377420E-03_wp,-4.03975654473191E-04_wp,-7.19052804954172E-03_wp, &
      &-7.72046463324211E-04_wp, 4.73561694057961E-03_wp, 4.18770813681427E-02_wp, &
      &-2.50290440454551E-03_wp, 7.65727418567875E-04_wp, 2.12624345393360E-03_wp, &
      &-4.97027877917978E-03_wp,-5.83472669844681E-03_wp, 5.51732011957964E-03_wp, &
      &-4.27826077946260E-04_wp,-8.82403162605877E-03_wp, 3.37743473393135E-02_wp, &
      &-3.36885096134562E-02_wp,-8.78793958800207E-02_wp,-1.20767003906759E-03_wp, &
      &-4.48556039118914E-04_wp,-1.81086224545010E-04_wp,-1.30225149564302E-03_wp, &
      &-1.72564923333566E-03_wp,-8.62878718054269E-03_wp,-3.71725125852225E-02_wp, &
      & 3.53696185434499E-02_wp,-9.55348684789787E-02_wp, 1.04090454647082E-03_wp, &
      &-3.59221400042112E-04_wp,-3.57354903215494E-04_wp, 1.79429741173164E-03_wp, &
      &-2.19262341202585E-03_wp, 2.19832369216930E-02_wp, 7.20780368313856E-02_wp, &
      &-7.10320337351804E-02_wp,-4.72446672147909E-02_wp,-2.04798561894034E-03_wp, &
      &-1.49699282601156E-03_wp, 1.86593826628326E-04_wp, 2.32423304833032E-03_wp, &
      & 4.22160666293566E-04_wp,-2.43299151738308E-03_wp, 1.06766405478892E-03_wp, &
      &-1.94547758751811E-04_wp,-1.81460080347714E-03_wp,-3.26942203750127E-04_wp, &
      & 3.26384647666875E-03_wp, 3.39812061949517E-03_wp,-9.20066557038704E-03_wp, &
      &-2.50290440454551E-03_wp, 1.18831170923388E-02_wp, 3.11417561356817E-03_wp, &
      &-1.79633667274163E-03_wp, 1.08428680718754E-03_wp, 2.95504975538613E-03_wp, &
      & 9.40583767877227E-04_wp,-8.18360385983116E-04_wp, 2.70527565881747E-05_wp, &
      &-3.42955794705075E-02_wp,-1.00765000114561E-02_wp, 3.53152003768071E-02_wp, &
      & 3.49311376473727E-04_wp,-7.42699380840104E-04_wp, 4.39431782527856E-06_wp, &
      & 3.27030474840866E-04_wp, 1.58749721680075E-03_wp, 5.59939284193674E-04_wp, &
      &-4.33600687024676E-02_wp,-1.57001763493095E-02_wp,-3.92643295028619E-02_wp, &
      & 4.03441445098763E-04_wp, 1.24040844923251E-03_wp, 2.73735765525048E-04_wp, &
      & 2.82584056419282E-04_wp,-2.02260087963652E-03_wp,-8.83664108597079E-04_wp, &
      & 3.13114544459643E-02_wp, 3.67817089890721E-02_wp, 5.10642880072443E-02_wp, &
      &-6.93371177426940E-04_wp,-1.52128659470246E-04_wp,-6.60971351200672E-05_wp, &
      &-2.04580162271832E-03_wp,-1.86857557022780E-03_wp,-9.06679526463874E-03_wp, &
      & 1.23297231618413E-03_wp,-6.02151183304568E-04_wp, 2.37543992326715E-03_wp, &
      & 3.91589546660145E-03_wp,-1.16648070070560E-02_wp, 1.80295732067388E-03_wp, &
      &-4.78340485526786E-03_wp, 7.65727418567875E-04_wp, 3.11417561356817E-03_wp, &
      & 1.53506396343270E-02_wp, 4.93458019575497E-04_wp, 1.16823409155740E-03_wp, &
      & 3.13052003010725E-04_wp,-1.86532900379925E-03_wp,-8.96110846963260E-03_wp, &
      &-9.72183490348135E-04_wp,-6.96677285675340E-02_wp,-2.42640812741529E-02_wp, &
      &-2.43105967018250E-02_wp,-1.88909058022181E-03_wp,-1.84620872265780E-03_wp, &
      & 8.40605401675121E-04_wp,-8.31609991216517E-04_wp, 1.57654392910010E-03_wp, &
      &-6.32238839105142E-04_wp, 6.52028457208331E-02_wp, 1.09811119074071E-04_wp, &
      &-2.52064987261408E-02_wp, 1.91383514382578E-03_wp,-1.05856014409679E-03_wp, &
      & 1.27412333131771E-03_wp, 3.48350771197466E-04_wp, 1.31355340280619E-03_wp, &
      & 5.25854327644141E-05_wp, 2.45048495855217E-02_wp,-2.56655009577627E-03_wp, &
      & 4.06754362288473E-02_wp,-4.24473085273045E-04_wp,-4.73872134127807E-04_wp, &
      & 6.75497357338132E-04_wp,-5.61485634191517E-04_wp,-1.60968883271809E-03_wp, &
      &-3.80607531989858E-04_wp,-1.24027234557119E-03_wp, 1.49179954780989E-03_wp, &
      & 6.76725101117992E-04_wp, 4.64342996080958E-03_wp,-1.74839240093656E-03_wp, &
      &-1.72699146503798E-03_wp,-5.05261443662806E-03_wp, 2.12624345393360E-03_wp, &
      &-1.79633667274163E-03_wp, 4.93458019575497E-04_wp, 5.57899823593500E-03_wp, &
      &-1.29897186515772E-04_wp,-7.80032096534039E-04_wp,-1.84892158652720E-04_wp, &
      &-9.62092557672180E-04_wp, 3.41667470650106E-05_wp,-2.14598915207327E-02_wp, &
      & 4.22313704080375E-02_wp,-1.51574111619925E-02_wp,-7.82803317697533E-04_wp, &
      & 8.27111347825122E-04_wp, 1.17930469715194E-03_wp, 5.54949096391098E-04_wp, &
      & 3.31242618197427E-04_wp,-2.19959469457473E-04_wp,-3.87427809378408E-03_wp, &
      & 4.92359767360958E-02_wp, 3.76026077555199E-04_wp,-1.50246094084959E-04_wp, &
      &-1.27777504135598E-03_wp,-1.02227177755387E-03_wp, 9.55873747395274E-04_wp, &
      &-4.76745213910508E-05_wp, 5.49263377585457E-04_wp,-6.37450174956073E-03_wp, &
      &-1.46067658017479E-02_wp,-4.95983018633469E-03_wp, 1.86494551430902E-04_wp, &
      &-2.30162533491548E-05_wp, 1.12528083885465E-04_wp, 5.83275218807750E-04_wp, &
      & 4.19280315359991E-04_wp,-2.37874651841337E-03_wp,-2.62975652749730E-03_wp, &
      &-2.93152544359381E-03_wp,-1.81030095889485E-03_wp, 1.33992841815537E-02_wp, &
      & 1.30572623508845E-03_wp, 2.08068962045090E-04_wp, 6.84791172570740E-04_wp, &
      &-4.97027877917978E-03_wp, 1.08428680718754E-03_wp, 1.16823409155740E-03_wp, &
      &-1.29897186515772E-04_wp, 9.94428221843759E-03_wp,-7.19846587011024E-05_wp, &
      &-1.93960408497296E-04_wp,-2.74856880867185E-03_wp, 1.27031368525096E-03_wp, &
      & 4.18745285288000E-02_wp, 2.72966013888841E-02_wp, 2.85915868985861E-02_wp, &
      & 1.52421575467950E-03_wp, 1.17133429597759E-03_wp,-5.70259377186292E-04_wp, &
      & 8.00501008899383E-04_wp,-5.93648368390666E-04_wp,-8.89167773592477E-04_wp, &
      & 4.01388520677539E-02_wp, 6.11538319289190E-03_wp,-3.43827040289076E-02_wp, &
      & 1.81283369635923E-03_wp,-6.03244529432152E-04_wp, 9.60330421290320E-04_wp, &
      & 6.66961474578595E-04_wp, 1.38916670200248E-04_wp,-8.40001714013650E-04_wp, &
      &-1.91461291146221E-02_wp, 5.62002262438172E-03_wp, 4.94693462764965E-02_wp, &
      & 1.30712345342775E-03_wp, 4.10116747608010E-04_wp, 9.58776573157464E-04_wp, &
      &-7.52025031397998E-04_wp,-1.22301380113751E-03_wp,-3.06390494210790E-04_wp, &
      & 6.81819221360930E-04_wp, 1.04205333697116E-03_wp,-1.71651898481651E-03_wp, &
      & 7.49785230431801E-04_wp,-5.03897420550892E-03_wp,-9.57025303579991E-04_wp, &
      & 2.03099825524059E-03_wp,-5.83472669844681E-03_wp, 2.95504975538613E-03_wp, &
      & 3.13052003010725E-04_wp,-7.80032096534039E-04_wp,-7.19846587011024E-05_wp, &
      & 6.02218778761440E-03_wp,-1.96369388859724E-04_wp, 4.00523412011083E-04_wp, &
      & 2.43415200594010E-05_wp,-1.51845049947398E-02_wp, 2.87327242622507E-02_wp, &
      &-9.49249765705680E-03_wp,-5.48713980681419E-04_wp, 5.54022846554013E-04_wp, &
      & 8.01687607586436E-04_wp, 3.77777202539557E-04_wp, 1.95523826272956E-04_wp, &
      & 6.68905535754955E-04_wp,-5.70282419272623E-03_wp,-3.25394854693152E-02_wp, &
      &-4.86038808135561E-03_wp,-6.16005794351641E-05_wp, 1.02574808923772E-03_wp, &
      & 4.82756606241754E-04_wp,-5.53868614429023E-04_wp,-4.61704646590128E-04_wp, &
      &-4.42109441645594E-04_wp, 1.09535528038150E-03_wp, 5.78076943435337E-02_wp, &
      &-1.68594875769994E-03_wp,-1.50326172865466E-04_wp, 5.22664015462132E-04_wp, &
      &-1.09131629766287E-03_wp,-1.82188558275468E-03_wp,-3.25003079162081E-05_wp, &
      & 3.36294433121553E-03_wp,-2.77984913399897E-03_wp, 1.48481072357888E-04_wp, &
      & 2.34703941229088E-03_wp, 8.17240065739435E-04_wp,-8.96763548845411E-04_wp, &
      &-6.72245961471031E-04_wp, 9.45917682510810E-03_wp, 5.51732011957964E-03_wp, &
      & 9.40583767877227E-04_wp,-1.86532900379925E-03_wp,-1.84892158652720E-04_wp, &
      &-1.93960408497296E-04_wp,-1.96369388859724E-04_wp, 1.14272521505678E-02_wp, &
      & 3.81891246340942E-03_wp, 4.46942267787944E-04_wp,-4.11175285560280E-03_wp, &
      & 9.67588448002750E-03_wp, 5.64940155335306E-02_wp, 1.60503540214145E-03_wp, &
      & 1.52655440518358E-04_wp,-3.37131360832177E-04_wp, 1.19451943984374E-03_wp, &
      & 1.48935194146030E-03_wp,-9.31069490872293E-04_wp, 7.87169218444636E-03_wp, &
      & 3.89347577592281E-03_wp,-4.72243223724727E-02_wp, 1.35155308784772E-03_wp, &
      &-2.00017008301476E-04_wp, 5.23748289028596E-04_wp, 8.65609593211022E-04_wp, &
      &-8.12160921509295E-04_wp, 1.03952145445761E-03_wp, 3.42404939671829E-02_wp, &
      &-2.44882129383441E-03_wp,-6.15913706401642E-02_wp,-1.62689399094805E-03_wp, &
      &-5.53268037298782E-04_wp,-1.20341681703358E-03_wp, 1.03007344334728E-03_wp, &
      & 1.53200056222919E-03_wp, 8.99503672589712E-03_wp, 6.15223438265065E-04_wp, &
      & 6.08528200623539E-04_wp,-2.97757430265402E-04_wp,-3.19238355446981E-04_wp, &
      & 9.84759792548497E-03_wp,-1.25000483023965E-02_wp, 3.78879388099217E-03_wp, &
      &-4.27826077946260E-04_wp,-8.18360385983116E-04_wp,-8.96110846963260E-03_wp, &
      &-9.62092557672180E-04_wp,-2.74856880867185E-03_wp, 4.00523412011083E-04_wp, &
      & 3.81891246340942E-03_wp, 1.48383546641407E-02_wp, 2.14569885077566E-04_wp, &
      & 2.60228009071659E-02_wp, 4.55480551778864E-02_wp, 3.11414586934289E-02_wp, &
      & 1.24763488181776E-03_wp, 1.53320635726693E-03_wp, 7.48314603836746E-08_wp, &
      & 1.22592371170191E-03_wp,-1.03037430243239E-04_wp, 6.85864626229303E-04_wp, &
      &-3.93198699151587E-02_wp,-4.33959140312128E-02_wp, 3.84544626969406E-02_wp, &
      &-1.77253159605672E-03_wp, 1.76953567734716E-03_wp,-2.92690161508619E-04_wp, &
      &-1.42187034413862E-03_wp,-3.56619680821680E-04_wp, 2.40006331974931E-04_wp, &
      & 4.74589429392484E-02_wp,-2.86025601155193E-02_wp,-2.74869037414869E-02_wp, &
      &-1.81297319011913E-03_wp,-9.90759348859760E-04_wp,-3.97806973349798E-04_wp, &
      & 1.22820555419493E-03_wp, 4.56941288207337E-04_wp,-1.11688855965218E-03_wp, &
      & 7.17489050802549E-03_wp, 4.57690583399911E-03_wp, 4.75605525497978E-03_wp, &
      & 2.23646244137205E-02_wp, 2.07224740193060E-02_wp,-7.48190813025657E-03_wp, &
      & 1.40822121394522E-02_wp,-8.82403162605877E-03_wp, 2.70527565881747E-05_wp, &
      &-9.72183490348135E-04_wp, 3.41667470650106E-05_wp, 1.27031368525096E-03_wp, &
      & 2.43415200594010E-05_wp, 4.46942267787944E-04_wp, 2.14569885077566E-04_wp, &
      & 9.79953185022663E-01_wp,-6.56913535165527E-02_wp,-4.15746208633432E-02_wp, &
      &-4.42906788242035E-02_wp,-3.26497777864028E-03_wp,-3.02957050163200E-03_wp, &
      & 1.03160173383718E-03_wp,-2.06224126752306E-03_wp, 1.29590186401169E-03_wp, &
      &-1.73578590097797E-03_wp,-4.78755981358030E-03_wp,-4.80964039210149E-03_wp, &
      &-2.26550292609638E-02_wp, 1.03244385712820E-03_wp, 8.28645456652198E-04_wp, &
      & 8.36183698869690E-04_wp, 3.56016852022734E-04_wp,-2.26721397739033E-03_wp, &
      &-4.67936049598523E-04_wp,-2.15316812457995E-02_wp,-5.04588128389799E-03_wp, &
      & 2.25098483451279E-03_wp, 2.05812072051300E-03_wp, 5.05766248283278E-04_wp, &
      & 7.91454294590990E-04_wp, 7.16816433450695E-04_wp, 1.18902455312510E-03_wp, &
      & 3.20046292353588E-02_wp,-1.15284447884868E-02_wp,-1.88737263402893E-02_wp, &
      &-1.99507322125130E-02_wp, 6.59548065857667E-02_wp, 6.70688097195503E-02_wp, &
      & 1.50194621101850E-02_wp, 7.85903417261118E-02_wp, 3.37743473393135E-02_wp, &
      &-3.42955794705075E-02_wp,-6.96677285675340E-02_wp,-2.14598915207327E-02_wp, &
      & 4.18745285288000E-02_wp,-1.51845049947398E-02_wp,-4.11175285560280E-03_wp, &
      & 2.60228009071659E-02_wp,-6.56913535165527E-02_wp, 8.97832343104439E-01_wp, &
      &-9.01190569007224E-03_wp,-1.45144645636852E-02_wp, 1.55005079703965E-02_wp, &
      & 1.20367234816796E-02_wp,-1.34556376466509E-02_wp,-3.44948217516065E-04_wp, &
      &-2.38477695614042E-02_wp, 5.97412493903286E-03_wp, 3.58652467089850E-03_wp, &
      &-1.13563354350285E-02_wp,-9.59168996535537E-03_wp,-6.71551375240948E-04_wp, &
      & 2.30409327479318E-03_wp, 4.33495391693732E-04_wp,-1.03517926424578E-03_wp, &
      &-9.32392910000505E-04_wp,-1.75636665246337E-02_wp,-4.80394095409307E-02_wp, &
      &-4.04081526462186E-03_wp, 2.91067062791004E-02_wp, 5.41082241160409E-03_wp, &
      & 2.05070858833908E-03_wp, 2.29339909979139E-03_wp, 4.07757821937570E-04_wp, &
      &-4.27488088446933E-04_wp, 2.32606695048053E-02_wp,-1.66995151518806E-02_wp, &
      & 8.75664181962324E-03_wp,-1.13575027794173E-02_wp, 8.68754738499316E-02_wp, &
      &-1.96459370841130E-02_wp,-9.51746257217108E-02_wp,-1.34573967461195E-02_wp, &
      &-3.36885096134562E-02_wp,-1.00765000114561E-02_wp,-2.42640812741529E-02_wp, &
      & 4.22313704080375E-02_wp, 2.72966013888841E-02_wp, 2.87327242622507E-02_wp, &
      & 9.67588448002750E-03_wp, 4.55480551778864E-02_wp,-4.15746208633432E-02_wp, &
      &-9.01190569007224E-03_wp, 8.98949522434250E-01_wp,-5.85249818311979E-03_wp, &
      &-4.20396123694783E-05_wp, 2.32696964109612E-02_wp, 1.40699101827117E-02_wp, &
      & 1.57206108208981E-02_wp, 4.82043997286318E-05_wp, 1.29832452985130E-03_wp, &
      &-1.01189621165397E-02_wp, 1.64501682042587E-02_wp,-1.13646125614395E-02_wp, &
      &-6.58322439410475E-04_wp,-9.95021813207843E-05_wp,-1.69154717257235E-03_wp, &
      & 8.98861322776682E-04_wp,-2.17156562704802E-03_wp, 1.74725131414714E-03_wp, &
      &-1.23060391666245E-02_wp, 1.53478173145747E-02_wp,-5.75999068005425E-03_wp, &
      & 6.00827861785124E-04_wp, 6.76325363903056E-04_wp,-1.66014282991114E-03_wp, &
      &-2.79625792411043E-04_wp, 2.06014446932550E-03_wp, 2.10966188523540E-02_wp, &
      &-1.96285390992293E-02_wp,-1.26486079338153E-02_wp, 4.49366063601718E-03_wp, &
      & 1.25040456174467E-03_wp, 7.81809041418509E-02_wp, 1.07105823451796E-02_wp, &
      & 5.53985700082034E-03_wp,-8.78793958800207E-02_wp, 3.53152003768071E-02_wp, &
      &-2.43105967018250E-02_wp,-1.51574111619925E-02_wp, 2.85915868985861E-02_wp, &
      &-9.49249765705680E-03_wp, 5.64940155335306E-02_wp, 3.11414586934289E-02_wp, &
      &-4.42906788242035E-02_wp,-1.45144645636852E-02_wp,-5.85249818311979E-03_wp, &
      & 9.09568901510487E-01_wp, 2.36486322778074E-02_wp,-3.50001068677793E-04_wp, &
      &-9.07751827589855E-03_wp, 1.23464513582589E-02_wp, 1.66020669497448E-02_wp, &
      &-2.28312252496096E-02_wp, 1.33634035581138E-02_wp, 1.41764609986243E-03_wp, &
      &-5.45979916231767E-02_wp, 5.34400854366441E-03_wp,-4.96044569603787E-06_wp, &
      & 2.40201841492500E-03_wp, 2.48198884031881E-03_wp,-2.72121814524801E-03_wp, &
      & 1.34763066806292E-02_wp, 1.29062067504150E-02_wp,-9.80190748894087E-03_wp, &
      &-5.04907733624699E-03_wp,-2.40085902556875E-03_wp,-1.16119924081091E-03_wp, &
      &-5.81093006975118E-04_wp, 2.45027490200433E-03_wp, 1.39966952488796E-03_wp, &
      & 1.14006433275913E-03_wp,-8.09390428732359E-04_wp,-7.10043290946803E-04_wp, &
      &-2.46172687351389E-04_wp, 8.43040832431269E-04_wp, 2.84420597616065E-03_wp, &
      & 2.74276798196848E-04_wp, 1.51665855913071E-03_wp,-1.20767003906759E-03_wp, &
      & 3.49311376473727E-04_wp,-1.88909058022181E-03_wp,-7.82803317697533E-04_wp, &
      & 1.52421575467950E-03_wp,-5.48713980681419E-04_wp, 1.60503540214145E-03_wp, &
      & 1.24763488181776E-03_wp,-3.26497777864028E-03_wp, 1.55005079703965E-02_wp, &
      &-4.20396123694783E-05_wp, 2.36486322778074E-02_wp, 9.18511119490073E-04_wp, &
      & 2.19078907283795E-04_wp,-4.80302880422056E-04_wp, 3.36883951800646E-04_wp, &
      & 1.96443072182944E-05_wp,-9.42305899811192E-04_wp,-2.89537305206576E-04_wp, &
      &-1.28022878279291E-03_wp,-5.00721188868637E-03_wp, 1.98157365462310E-04_wp, &
      & 6.86048439806489E-05_wp, 1.12897189556710E-04_wp, 8.04397669198610E-05_wp, &
      &-1.60993997929097E-04_wp, 9.24348836836812E-04_wp,-1.64629006596351E-03_wp, &
      &-1.88908493067724E-03_wp,-5.38795740743591E-04_wp, 5.99810932074870E-05_wp, &
      & 1.05424830110441E-05_wp, 4.34016120608557E-05_wp, 1.16049685662805E-04_wp, &
      & 5.66604383383437E-05_wp, 1.14147229303145E-03_wp,-6.27932511020678E-04_wp, &
      & 2.89620397193752E-06_wp,-6.32133131699997E-04_wp, 2.57139294590359E-03_wp, &
      & 3.81213248123500E-04_wp,-2.28477526111477E-03_wp, 4.92005695111929E-04_wp, &
      &-4.48556039118914E-04_wp,-7.42699380840104E-04_wp,-1.84620872265780E-03_wp, &
      & 8.27111347825122E-04_wp, 1.17133429597759E-03_wp, 5.54022846554013E-04_wp, &
      & 1.52655440518358E-04_wp, 1.53320635726693E-03_wp,-3.02957050163200E-03_wp, &
      & 1.20367234816796E-02_wp, 2.32696964109612E-02_wp,-3.50001068677793E-04_wp, &
      & 2.19078907283795E-04_wp, 7.86474390598036E-04_wp, 1.81027823094628E-04_wp, &
      & 4.09324704116321E-04_wp,-3.38483160385702E-04_wp, 7.13341494118634E-04_wp, &
      &-2.49338759770268E-03_wp, 1.68735866975485E-04_wp,-7.41309576022951E-04_wp, &
      &-5.88299190152288E-05_wp, 5.59832212371159E-05_wp,-6.75477375787827E-05_wp, &
      & 1.76396292110411E-05_wp,-1.27813548615559E-04_wp,-3.12515170166521E-04_wp, &
      &-3.56717280611371E-03_wp, 3.00321001380880E-04_wp,-6.77069550327594E-04_wp, &
      & 1.50498405466485E-04_wp, 8.35670166748950E-05_wp,-1.78507090937897E-05_wp, &
      & 5.03854201101717E-06_wp, 9.88365521261470E-05_wp,-3.02272998088046E-04_wp, &
      & 1.21450766756529E-04_wp, 6.06131016453163E-04_wp, 8.27214390433451E-05_wp, &
      & 4.69721320177459E-04_wp,-1.93389141426929E-03_wp,-1.54670546994761E-03_wp, &
      &-1.33814279600240E-03_wp,-1.81086224545010E-04_wp, 4.39431782527856E-06_wp, &
      & 8.40605401675121E-04_wp, 1.17930469715194E-03_wp,-5.70259377186292E-04_wp, &
      & 8.01687607586436E-04_wp,-3.37131360832177E-04_wp, 7.48314603836746E-08_wp, &
      & 1.03160173383718E-03_wp,-1.34556376466509E-02_wp, 1.40699101827117E-02_wp, &
      &-9.07751827589855E-03_wp,-4.80302880422056E-04_wp, 1.81027823094628E-04_wp, &
      & 5.18247582373794E-04_wp, 1.22109199673484E-04_wp, 1.93563736986084E-04_wp, &
      & 6.94512792359044E-04_wp,-7.43041114473141E-04_wp, 1.90087784391794E-03_wp, &
      & 1.41762342890863E-03_wp,-8.72768819207063E-05_wp,-6.27592023472381E-05_wp, &
      &-9.94397095236854E-05_wp, 1.69999823758605E-05_wp, 1.02867249071339E-05_wp, &
      & 6.86263434466501E-04_wp, 1.12323534619301E-03_wp, 1.84122040819062E-03_wp, &
      &-1.05685564533458E-03_wp,-7.41528292047566E-05_wp,-5.64959811538514E-06_wp, &
      &-9.38849407025514E-05_wp,-6.44076461307671E-05_wp, 4.02013059080886E-05_wp, &
      & 7.61032562708332E-04_wp,-6.28894926244606E-04_wp, 2.82293389952361E-06_wp, &
      &-1.20209486448711E-04_wp, 1.34881283254879E-03_wp, 4.82420477979741E-04_wp, &
      &-1.52509696471355E-03_wp, 9.60493988837844E-06_wp,-1.30225149564302E-03_wp, &
      & 3.27030474840866E-04_wp,-8.31609991216517E-04_wp, 5.54949096391098E-04_wp, &
      & 8.00501008899383E-04_wp, 3.77777202539557E-04_wp, 1.19451943984374E-03_wp, &
      & 1.22592371170191E-03_wp,-2.06224126752306E-03_wp,-3.44948217516065E-04_wp, &
      & 1.57206108208981E-02_wp, 1.23464513582589E-02_wp, 3.36883951800646E-04_wp, &
      & 4.09324704116321E-04_wp, 1.22109199673484E-04_wp, 4.58251142523757E-04_wp, &
      & 2.38692760788488E-04_wp,-7.45009963602784E-04_wp,-8.81995629091355E-04_wp, &
      & 4.45232108061955E-04_wp,-3.31802835554324E-03_wp, 9.72796042636698E-05_wp, &
      & 4.59931850098103E-06_wp, 8.59895681225748E-06_wp, 8.83186095660335E-05_wp, &
      &-1.38666929038812E-04_wp, 9.28623086750384E-04_wp, 1.19596404451922E-04_wp, &
      & 5.87522199748577E-06_wp,-2.03335470252935E-03_wp,-4.17552021346788E-05_wp, &
      &-5.09405249893958E-06_wp,-6.69138499163433E-05_wp, 5.12364190514789E-05_wp, &
      & 1.03787374695025E-04_wp,-4.79144450852178E-04_wp,-7.62484937973928E-05_wp, &
      & 2.86985578991511E-04_wp, 7.01975281929006E-04_wp,-1.20147515573048E-03_wp, &
      &-5.32810176688762E-04_wp,-8.34604782021266E-05_wp,-1.49543254126618E-03_wp, &
      &-1.72564923333566E-03_wp, 1.58749721680075E-03_wp, 1.57654392910010E-03_wp, &
      & 3.31242618197427E-04_wp,-5.93648368390666E-04_wp, 1.95523826272956E-04_wp, &
      & 1.48935194146030E-03_wp,-1.03037430243239E-04_wp, 1.29590186401169E-03_wp, &
      &-2.38477695614042E-02_wp, 4.82043997286318E-05_wp, 1.66020669497448E-02_wp, &
      & 1.96443072182944E-05_wp,-3.38483160385702E-04_wp, 1.93563736986084E-04_wp, &
      & 2.38692760788488E-04_wp, 9.55996160124409E-04_wp,-2.25331678846492E-03_wp, &
      & 1.15282575115815E-03_wp, 1.46723806272380E-03_wp,-2.97238116407506E-03_wp, &
      & 1.74422156274755E-04_wp,-9.53909049759322E-05_wp, 4.56054718637480E-05_wp, &
      & 1.19023687658630E-04_wp,-4.22374644344222E-05_wp, 2.15549509066777E-03_wp, &
      & 5.25101958893515E-03_wp, 1.29756344050563E-05_wp,-2.17978866795573E-03_wp, &
      &-2.95874145259961E-04_wp,-1.25680696991265E-04_wp,-1.08710889283059E-04_wp, &
      & 5.08755648891243E-05_wp, 3.70646570213923E-05_wp,-1.30082860861281E-03_wp, &
      &-8.00507020680877E-03_wp,-3.80262208043012E-03_wp, 5.46887348151783E-03_wp, &
      &-2.39915293849405E-02_wp, 1.83968711176355E-02_wp,-9.70631584445687E-03_wp, &
      &-1.19162210963848E-02_wp,-8.62878718054269E-03_wp, 5.59939284193674E-04_wp, &
      &-6.32238839105142E-04_wp,-2.19959469457473E-04_wp,-8.89167773592477E-04_wp, &
      & 6.68905535754955E-04_wp,-9.31069490872293E-04_wp, 6.85864626229303E-04_wp, &
      &-1.73578590097797E-03_wp, 5.97412493903286E-03_wp, 1.29832452985130E-03_wp, &
      &-2.28312252496096E-02_wp,-9.42305899811192E-04_wp, 7.13341494118634E-04_wp, &
      & 6.94512792359044E-04_wp,-7.45009963602784E-04_wp,-2.25331678846492E-03_wp, &
      & 9.80974511856617E-01_wp, 6.77159511266575E-02_wp, 3.31547266836938E-02_wp, &
      &-4.62531679479647E-02_wp, 3.56741815719079E-03_wp,-2.66711332001603E-03_wp, &
      & 1.45224189985259E-03_wp, 1.78354903062082E-03_wp, 1.32108871221265E-03_wp, &
      &-3.02557120515403E-03_wp, 1.59240225965459E-02_wp, 1.87604065635367E-02_wp, &
      & 4.69837569301324E-03_wp,-1.69113411385192E-03_wp,-6.61974620886243E-04_wp, &
      &-1.25193134896855E-03_wp,-1.80059796257016E-03_wp,-3.51490801006542E-04_wp, &
      &-3.51696354973964E-02_wp,-1.18379012605914E-02_wp,-1.56986199671825E-02_wp, &
      & 2.07465098916450E-02_wp, 7.91814684776250E-02_wp,-6.46585082512025E-02_wp, &
      & 5.89106335176432E-04_wp, 7.77031024200890E-02_wp,-3.71725125852225E-02_wp, &
      &-4.33600687024676E-02_wp, 6.52028457208331E-02_wp,-3.87427809378408E-03_wp, &
      & 4.01388520677539E-02_wp,-5.70282419272623E-03_wp, 7.87169218444636E-03_wp, &
      &-3.93198699151587E-02_wp,-4.78755981358030E-03_wp, 3.58652467089850E-03_wp, &
      &-1.01189621165397E-02_wp, 1.33634035581138E-02_wp,-2.89537305206576E-04_wp, &
      &-2.49338759770268E-03_wp,-7.43041114473141E-04_wp,-8.81995629091355E-04_wp, &
      & 1.15282575115815E-03_wp, 6.77159511266575E-02_wp, 8.93041610503373E-01_wp, &
      &-4.10359825688933E-03_wp, 1.29601523129638E-02_wp, 1.67095113530953E-02_wp, &
      &-1.36647926268428E-02_wp, 1.27251154417897E-02_wp,-3.25733892456355E-05_wp, &
      & 2.26849602998123E-02_wp, 9.66633023541558E-03_wp,-5.19297868799919E-03_wp, &
      &-3.84633749324155E-02_wp,-1.24613040121676E-02_wp, 1.38676485413577E-03_wp, &
      &-6.88153886369985E-04_wp, 2.13804102298870E-03_wp, 3.54546470887379E-03_wp, &
      &-6.75694770587501E-06_wp,-1.39407590739093E-02_wp,-1.74855485990935E-02_wp, &
      & 9.95289435445029E-03_wp, 1.12870324819613E-02_wp, 6.88501216429399E-02_wp, &
      & 3.76360397992066E-02_wp, 8.54634077370308E-02_wp,-3.00536616447891E-02_wp, &
      & 3.53696185434499E-02_wp,-1.57001763493095E-02_wp, 1.09811119074071E-04_wp, &
      & 4.92359767360958E-02_wp, 6.11538319289190E-03_wp,-3.25394854693152E-02_wp, &
      & 3.89347577592281E-03_wp,-4.33959140312128E-02_wp,-4.80964039210149E-03_wp, &
      &-1.13563354350285E-02_wp, 1.64501682042587E-02_wp, 1.41764609986243E-03_wp, &
      &-1.28022878279291E-03_wp, 1.68735866975485E-04_wp, 1.90087784391794E-03_wp, &
      & 4.45232108061955E-04_wp, 1.46723806272380E-03_wp, 3.31547266836938E-02_wp, &
      &-4.10359825688933E-03_wp, 9.08117769603495E-01_wp, 9.12873283600103E-03_wp, &
      &-3.41667105887394E-04_wp,-2.27240427085580E-02_wp,-1.60656013699057E-02_wp, &
      & 1.70173883505189E-02_wp, 3.01970480128187E-04_wp, 1.87221098337115E-02_wp, &
      &-2.32165734653952E-02_wp,-2.18459102876546E-02_wp,-3.26536302874007E-02_wp, &
      & 1.96920444748069E-03_wp, 5.26310641005283E-04_wp,-1.60834820233167E-05_wp, &
      & 3.65950247565885E-03_wp, 4.03502323498763E-03_wp, 2.24022500427106E-02_wp, &
      & 2.22060923985087E-02_wp, 1.08109379474349E-02_wp, 2.62360032993394E-03_wp, &
      &-1.27801136787652E-02_wp, 6.76726646968307E-02_wp,-6.23156571113132E-03_wp, &
      &-2.00732274111886E-03_wp,-9.55348684789787E-02_wp,-3.92643295028619E-02_wp, &
      &-2.52064987261408E-02_wp, 3.76026077555199E-04_wp,-3.43827040289076E-02_wp, &
      &-4.86038808135561E-03_wp,-4.72243223724727E-02_wp, 3.84544626969406E-02_wp, &
      &-2.26550292609638E-02_wp,-9.59168996535537E-03_wp,-1.13646125614395E-02_wp, &
      &-5.45979916231767E-02_wp,-5.00721188868637E-03_wp,-7.41309576022951E-04_wp, &
      & 1.41762342890863E-03_wp,-3.31802835554324E-03_wp,-2.97238116407506E-03_wp, &
      &-4.62531679479647E-02_wp, 1.29601523129638E-02_wp, 9.12873283600103E-03_wp, &
      & 9.09645248969580E-01_wp,-2.27265262407864E-02_wp,-5.32299199307859E-04_wp, &
      &-1.00495621815063E-02_wp,-1.37190880213787E-02_wp, 1.76759686357001E-02_wp, &
      & 1.35404942972661E-02_wp, 8.14988231623697E-03_wp,-1.12416390562587E-02_wp, &
      &-3.84377541807435E-03_wp,-1.75012611402543E-03_wp,-9.81189201589842E-04_wp, &
      &-3.98266659773968E-04_wp, 2.71148530315096E-03_wp, 1.45446242844250E-03_wp, &
      &-1.30001097214208E-03_wp,-8.41884382085824E-04_wp,-6.58426823850518E-04_wp, &
      & 3.30170019630099E-04_wp, 1.67318715110993E-03_wp,-2.56333530435510E-03_wp, &
      & 2.39237960359503E-04_wp, 1.10329669950783E-03_wp, 1.04090454647082E-03_wp, &
      & 4.03441445098763E-04_wp, 1.91383514382578E-03_wp,-1.50246094084959E-04_wp, &
      & 1.81283369635923E-03_wp,-6.16005794351641E-05_wp, 1.35155308784772E-03_wp, &
      &-1.77253159605672E-03_wp, 1.03244385712820E-03_wp,-6.71551375240948E-04_wp, &
      &-6.58322439410475E-04_wp, 5.34400854366441E-03_wp, 1.98157365462310E-04_wp, &
      &-5.88299190152288E-05_wp,-8.72768819207063E-05_wp, 9.72796042636698E-05_wp, &
      & 1.74422156274755E-04_wp, 3.56741815719079E-03_wp, 1.67095113530953E-02_wp, &
      &-3.41667105887394E-04_wp,-2.27265262407864E-02_wp, 9.20468008197760E-04_wp, &
      &-2.49874089186393E-04_wp, 5.11060765003657E-04_wp, 3.55781691260968E-04_wp, &
      &-2.00686966103915E-05_wp,-1.11069365429560E-03_wp,-1.77595705940923E-03_wp, &
      &-9.69766591163357E-04_wp, 1.79487344769638E-03_wp, 1.18227059103687E-04_wp, &
      & 2.23628169815579E-05_wp, 9.59540695118697E-05_wp,-6.27125251361419E-06_wp, &
      &-7.33810769916618E-05_wp, 9.03676498870593E-04_wp, 6.54018323036404E-04_wp, &
      &-2.72580344710504E-05_wp,-6.85792113938533E-04_wp,-2.26525736225408E-03_wp, &
      & 3.57609160918523E-04_wp,-1.86157192226641E-03_wp,-1.68670894926453E-04_wp, &
      &-3.59221400042112E-04_wp, 1.24040844923251E-03_wp,-1.05856014409679E-03_wp, &
      &-1.27777504135598E-03_wp,-6.03244529432152E-04_wp, 1.02574808923772E-03_wp, &
      &-2.00017008301476E-04_wp, 1.76953567734716E-03_wp, 8.28645456652198E-04_wp, &
      & 2.30409327479318E-03_wp,-9.95021813207843E-05_wp,-4.96044569603787E-06_wp, &
      & 6.86048439806489E-05_wp, 5.59832212371159E-05_wp,-6.27592023472381E-05_wp, &
      & 4.59931850098103E-06_wp,-9.53909049759322E-05_wp,-2.66711332001603E-03_wp, &
      &-1.36647926268428E-02_wp,-2.27240427085580E-02_wp,-5.32299199307859E-04_wp, &
      &-2.49874089186393E-04_wp, 8.05224008513454E-04_wp, 2.12701081422343E-04_wp, &
      &-4.33732331900280E-04_wp,-3.76373327444553E-04_wp,-1.92622491079819E-03_wp, &
      & 3.14448950646073E-03_wp, 3.64600111382412E-03_wp, 2.20063557970678E-03_wp, &
      &-1.22466148179130E-04_wp,-1.14817581188053E-05_wp,-6.53524280632076E-05_wp, &
      &-2.32450086452956E-04_wp,-1.56880573246973E-04_wp,-5.83171540616389E-04_wp, &
      &-9.79735922630597E-05_wp,-5.84203711930999E-04_wp, 6.13175091245618E-05_wp, &
      & 1.03312589060382E-04_wp,-2.02312098196774E-03_wp,-1.22914115729558E-03_wp, &
      & 1.46237940335675E-03_wp,-3.57354903215494E-04_wp, 2.73735765525048E-04_wp, &
      & 1.27412333131771E-03_wp,-1.02227177755387E-03_wp, 9.60330421290320E-04_wp, &
      & 4.82756606241754E-04_wp, 5.23748289028596E-04_wp,-2.92690161508619E-04_wp, &
      & 8.36183698869690E-04_wp, 4.33495391693732E-04_wp,-1.69154717257235E-03_wp, &
      & 2.40201841492500E-03_wp, 1.12897189556710E-04_wp,-6.75477375787827E-05_wp, &
      &-9.94397095236854E-05_wp, 8.59895681225748E-06_wp, 4.56054718637480E-05_wp, &
      & 1.45224189985259E-03_wp, 1.27251154417897E-02_wp,-1.60656013699057E-02_wp, &
      &-1.00495621815063E-02_wp, 5.11060765003657E-04_wp, 2.12701081422343E-04_wp, &
      & 5.85429402822957E-04_wp,-1.44929109080051E-04_wp, 1.19686671170248E-04_wp, &
      &-1.26135270204616E-03_wp, 3.43298532897802E-04_wp, 3.83817970980841E-06_wp, &
      & 2.16308617823632E-03_wp, 1.39871480504945E-05_wp,-1.14181893240403E-05_wp, &
      & 6.32963494893182E-05_wp,-6.58156831460961E-05_wp,-1.36918692981879E-04_wp, &
      &-5.86376985234641E-04_wp,-7.38014920537391E-04_wp, 1.67615675302901E-05_wp, &
      & 2.02457181360978E-04_wp, 1.25511859423495E-03_wp,-3.47311396599134E-04_wp, &
      & 1.62380339070375E-03_wp,-6.46986175423066E-04_wp, 1.79429741173164E-03_wp, &
      & 2.82584056419282E-04_wp, 3.48350771197466E-04_wp, 9.55873747395274E-04_wp, &
      & 6.66961474578595E-04_wp,-5.53868614429023E-04_wp, 8.65609593211022E-04_wp, &
      &-1.42187034413862E-03_wp, 3.56016852022734E-04_wp,-1.03517926424578E-03_wp, &
      & 8.98861322776682E-04_wp, 2.48198884031881E-03_wp, 8.04397669198610E-05_wp, &
      & 1.76396292110411E-05_wp, 1.69999823758605E-05_wp, 8.83186095660335E-05_wp, &
      & 1.19023687658630E-04_wp, 1.78354903062082E-03_wp,-3.25733892456355E-05_wp, &
      & 1.70173883505189E-02_wp,-1.37190880213787E-02_wp, 3.55781691260968E-04_wp, &
      &-4.33732331900280E-04_wp,-1.44929109080051E-04_wp, 5.39623595018820E-04_wp, &
      &-2.59213800285309E-04_wp,-4.02492245769748E-05_wp,-2.31816357724689E-03_wp, &
      &-6.97563672919024E-04_wp,-1.10536488762762E-03_wp, 1.02416470270356E-04_wp, &
      & 4.30064827780183E-05_wp, 1.06480846519756E-05_wp, 5.28437067393878E-05_wp, &
      & 9.07005439053863E-05_wp,-5.20172683811391E-04_wp, 1.67408952638943E-04_wp, &
      &-1.95672888720868E-04_wp, 6.65644753860844E-04_wp, 1.23392237280349E-03_wp, &
      &-5.79280405069415E-04_wp,-2.57899109827892E-04_wp, 1.68257298559314E-03_wp, &
      &-2.19262341202585E-03_wp,-2.02260087963652E-03_wp, 1.31355340280619E-03_wp, &
      &-4.76745213910508E-05_wp, 1.38916670200248E-04_wp,-4.61704646590128E-04_wp, &
      &-8.12160921509295E-04_wp,-3.56619680821680E-04_wp,-2.26721397739033E-03_wp, &
      &-9.32392910000505E-04_wp,-2.17156562704802E-03_wp,-2.72121814524801E-03_wp, &
      &-1.60993997929097E-04_wp,-1.27813548615559E-04_wp, 1.02867249071339E-05_wp, &
      &-1.38666929038812E-04_wp,-4.22374644344222E-05_wp, 1.32108871221265E-03_wp, &
      & 2.26849602998123E-02_wp, 3.01970480128187E-04_wp, 1.76759686357001E-02_wp, &
      &-2.00686966103915E-05_wp,-3.76373327444553E-04_wp, 1.19686671170248E-04_wp, &
      &-2.59213800285309E-04_wp, 9.32898777817762E-04_wp, 1.33666846700373E-03_wp, &
      &-3.22194483116891E-04_wp,-4.38986218965655E-03_wp,-8.07406407670555E-04_wp, &
      & 5.55356868900987E-06_wp,-6.21201449865174E-05_wp, 9.70425331312997E-05_wp, &
      & 2.32661759562063E-04_wp, 3.40426800813876E-05_wp,-1.20825923162660E-03_wp, &
      & 2.55413860729764E-03_wp,-3.63428057490533E-03_wp,-9.08181098068503E-03_wp, &
      &-1.35620731019085E-02_wp,-5.21056540703168E-03_wp,-9.88014923090613E-03_wp, &
      & 2.13600011676214E-02_wp, 2.19832369216930E-02_wp,-8.83664108597079E-04_wp, &
      & 5.25854327644141E-05_wp, 5.49263377585457E-04_wp,-8.40001714013650E-04_wp, &
      &-4.42109441645594E-04_wp, 1.03952145445761E-03_wp, 2.40006331974931E-04_wp, &
      &-4.67936049598523E-04_wp,-1.75636665246337E-02_wp, 1.74725131414714E-03_wp, &
      & 1.34763066806292E-02_wp, 9.24348836836812E-04_wp,-3.12515170166521E-04_wp, &
      & 6.86263434466501E-04_wp, 9.28623086750384E-04_wp, 2.15549509066777E-03_wp, &
      &-3.02557120515403E-03_wp, 9.66633023541558E-03_wp, 1.87221098337115E-02_wp, &
      & 1.35404942972661E-02_wp,-1.11069365429560E-03_wp,-1.92622491079819E-03_wp, &
      &-1.26135270204616E-03_wp,-4.02492245769748E-05_wp, 1.33666846700373E-03_wp, &
      & 9.80364656316398E-01_wp,-2.18477518470186E-02_wp, 3.31124366005085E-02_wp, &
      & 7.98537612496825E-02_wp, 1.97816816411445E-03_wp, 7.95287884239188E-04_wp, &
      & 1.47614735664917E-03_wp,-3.09112108744030E-03_wp,-3.26315136987568E-03_wp, &
      & 1.01661826271254E-02_wp, 1.44077180987680E-02_wp, 5.13544987508667E-03_wp, &
      & 1.33513181820170E-02_wp, 7.56618991833905E-02_wp, 3.66511396790259E-02_wp, &
      &-6.25144814937098E-03_wp, 3.25281161052191E-02_wp, 7.20780368313856E-02_wp, &
      & 3.13114544459643E-02_wp, 2.45048495855217E-02_wp,-6.37450174956073E-03_wp, &
      &-1.91461291146221E-02_wp, 1.09535528038150E-03_wp, 3.42404939671829E-02_wp, &
      & 4.74589429392484E-02_wp,-2.15316812457995E-02_wp,-4.80394095409307E-02_wp, &
      &-1.23060391666245E-02_wp, 1.29062067504150E-02_wp,-1.64629006596351E-03_wp, &
      &-3.56717280611371E-03_wp, 1.12323534619301E-03_wp, 1.19596404451922E-04_wp, &
      & 5.25101958893515E-03_wp, 1.59240225965459E-02_wp,-5.19297868799919E-03_wp, &
      &-2.32165734653952E-02_wp, 8.14988231623697E-03_wp,-1.77595705940923E-03_wp, &
      & 3.14448950646073E-03_wp, 3.43298532897802E-04_wp,-2.31816357724689E-03_wp, &
      &-3.22194483116891E-04_wp,-2.18477518470186E-02_wp, 9.14810407256850E-01_wp, &
      & 7.48033986208809E-03_wp, 5.09953559406022E-03_wp,-2.71827543681182E-02_wp, &
      &-1.35396965598019E-02_wp,-5.28749471869153E-03_wp,-4.55353025077479E-04_wp, &
      &-8.97158346346022E-03_wp,-1.38165545028928E-02_wp, 4.81975145073898E-03_wp, &
      & 1.02100793239369E-02_wp,-1.99268786228122E-02_wp, 3.15084854540211E-02_wp, &
      &-1.61959501492248E-02_wp, 8.60612099656805E-02_wp, 4.67381503377974E-02_wp, &
      &-7.10320337351804E-02_wp, 3.67817089890721E-02_wp,-2.56655009577627E-03_wp, &
      &-1.46067658017479E-02_wp, 5.62002262438172E-03_wp, 5.78076943435337E-02_wp, &
      &-2.44882129383441E-03_wp,-2.86025601155193E-02_wp,-5.04588128389799E-03_wp, &
      &-4.04081526462186E-03_wp, 1.53478173145747E-02_wp,-9.80190748894087E-03_wp, &
      &-1.88908493067724E-03_wp, 3.00321001380880E-04_wp, 1.84122040819062E-03_wp, &
      & 5.87522199748577E-06_wp, 1.29756344050563E-05_wp, 1.87604065635367E-02_wp, &
      &-3.84633749324155E-02_wp,-2.18459102876546E-02_wp,-1.12416390562587E-02_wp, &
      &-9.69766591163357E-04_wp, 3.64600111382412E-03_wp, 3.83817970980841E-06_wp, &
      &-6.97563672919024E-04_wp,-4.38986218965655E-03_wp, 3.31124366005085E-02_wp, &
      & 7.48033986208809E-03_wp, 9.07045088096433E-01_wp,-6.65538482288324E-03_wp, &
      &-4.69539950191276E-04_wp, 8.56203889791546E-03_wp,-1.56700021141413E-02_wp, &
      &-2.67566253419294E-02_wp,-5.14238323166957E-05_wp,-4.08742499801280E-02_wp, &
      & 1.15402494798059E-02_wp,-1.81555166175934E-02_wp,-2.34217202320808E-02_wp, &
      & 6.45053437127237E-02_wp, 4.22610487435226E-02_wp, 3.14099372304278E-03_wp, &
      &-1.03858523700383E-01_wp,-4.72446672147909E-02_wp, 5.10642880072443E-02_wp, &
      & 4.06754362288473E-02_wp,-4.95983018633469E-03_wp, 4.94693462764965E-02_wp, &
      &-1.68594875769994E-03_wp,-6.15913706401642E-02_wp,-2.74869037414869E-02_wp, &
      & 2.25098483451279E-03_wp, 2.91067062791004E-02_wp,-5.75999068005425E-03_wp, &
      &-5.04907733624699E-03_wp,-5.38795740743591E-04_wp,-6.77069550327594E-04_wp, &
      &-1.05685564533458E-03_wp,-2.03335470252935E-03_wp,-2.17978866795573E-03_wp, &
      & 4.69837569301324E-03_wp,-1.24613040121676E-02_wp,-3.26536302874007E-02_wp, &
      &-3.84377541807435E-03_wp, 1.79487344769638E-03_wp, 2.20063557970678E-03_wp, &
      & 2.16308617823632E-03_wp,-1.10536488762762E-03_wp,-8.07406407670555E-04_wp, &
      & 7.98537612496825E-02_wp, 5.09953559406022E-03_wp,-6.65538482288324E-03_wp, &
      & 8.85720560609369E-01_wp, 8.47405040928226E-03_wp, 4.32609508439594E-05_wp, &
      & 1.51332117219235E-02_wp,-1.32002681370673E-02_wp,-2.63946101633413E-02_wp, &
      &-7.18067691148896E-04_wp,-3.91472860045389E-04_wp,-3.91264029795358E-04_wp, &
      &-6.69813653944928E-04_wp,-9.38359904069108E-04_wp,-7.23150166636224E-04_wp, &
      & 3.78532516914335E-04_wp,-1.50740974288112E-03_wp,-2.04798561894034E-03_wp, &
      &-6.93371177426940E-04_wp,-4.24473085273045E-04_wp, 1.86494551430902E-04_wp, &
      & 1.30712345342775E-03_wp,-1.50326172865466E-04_wp,-1.62689399094805E-03_wp, &
      &-1.81297319011913E-03_wp, 2.05812072051300E-03_wp, 5.41082241160409E-03_wp, &
      & 6.00827861785124E-04_wp,-2.40085902556875E-03_wp, 5.99810932074870E-05_wp, &
      & 1.50498405466485E-04_wp,-7.41528292047566E-05_wp,-4.17552021346788E-05_wp, &
      &-2.95874145259961E-04_wp,-1.69113411385192E-03_wp, 1.38676485413577E-03_wp, &
      & 1.96920444748069E-03_wp,-1.75012611402543E-03_wp, 1.18227059103687E-04_wp, &
      &-1.22466148179130E-04_wp, 1.39871480504945E-05_wp, 1.02416470270356E-04_wp, &
      & 5.55356868900987E-06_wp, 1.97816816411445E-03_wp,-2.71827543681182E-02_wp, &
      &-4.69539950191276E-04_wp, 8.47405040928226E-03_wp, 9.24530237513602E-04_wp, &
      & 4.12508019867173E-04_wp, 3.21310649550766E-04_wp,-1.04518495185120E-04_wp, &
      & 1.07469239920213E-05_wp,-2.63641797653755E-04_wp,-1.85447400076914E-04_wp, &
      & 9.88365504780387E-06_wp,-4.38712775661282E-04_wp,-6.34043613155897E-04_wp, &
      &-6.96300015777415E-04_wp, 9.38872653962537E-04_wp, 3.82651994863863E-05_wp, &
      &-1.49699282601156E-03_wp,-1.52128659470246E-04_wp,-4.73872134127807E-04_wp, &
      &-2.30162533491548E-05_wp, 4.10116747608010E-04_wp, 5.22664015462132E-04_wp, &
      &-5.53268037298782E-04_wp,-9.90759348859760E-04_wp, 5.05766248283278E-04_wp, &
      & 2.05070858833908E-03_wp, 6.76325363903056E-04_wp,-1.16119924081091E-03_wp, &
      & 1.05424830110441E-05_wp, 8.35670166748950E-05_wp,-5.64959811538514E-06_wp, &
      &-5.09405249893958E-06_wp,-1.25680696991265E-04_wp,-6.61974620886243E-04_wp, &
      &-6.88153886369985E-04_wp, 5.26310641005283E-04_wp,-9.81189201589842E-04_wp, &
      & 2.23628169815579E-05_wp,-1.14817581188053E-05_wp,-1.14181893240403E-05_wp, &
      & 4.30064827780183E-05_wp,-6.21201449865174E-05_wp, 7.95287884239188E-04_wp, &
      &-1.35396965598019E-02_wp, 8.56203889791546E-03_wp, 4.32609508439594E-05_wp, &
      & 4.12508019867173E-04_wp, 2.88003393911525E-04_wp,-6.50182407583915E-05_wp, &
      &-2.53634831707603E-04_wp, 1.27233371966009E-04_wp,-5.93966433083412E-04_wp, &
      & 2.69372961745180E-05_wp,-5.75792938733272E-04_wp,-1.20006019568062E-04_wp, &
      & 3.14862311605270E-04_wp, 6.92000783727072E-04_wp,-1.19990637124712E-03_wp, &
      &-2.40313487407486E-03_wp, 1.86593826628326E-04_wp,-6.60971351200672E-05_wp, &
      & 6.75497357338132E-04_wp, 1.12528083885465E-04_wp, 9.58776573157464E-04_wp, &
      &-1.09131629766287E-03_wp,-1.20341681703358E-03_wp,-3.97806973349798E-04_wp, &
      & 7.91454294590990E-04_wp, 2.29339909979139E-03_wp,-1.66014282991114E-03_wp, &
      &-5.81093006975118E-04_wp, 4.34016120608557E-05_wp,-1.78507090937897E-05_wp, &
      &-9.38849407025514E-05_wp,-6.69138499163433E-05_wp,-1.08710889283059E-04_wp, &
      &-1.25193134896855E-03_wp, 2.13804102298870E-03_wp,-1.60834820233167E-05_wp, &
      &-3.98266659773968E-04_wp, 9.59540695118697E-05_wp,-6.53524280632076E-05_wp, &
      & 6.32963494893182E-05_wp, 1.06480846519756E-05_wp, 9.70425331312997E-05_wp, &
      & 1.47614735664917E-03_wp,-5.28749471869153E-03_wp,-1.56700021141413E-02_wp, &
      & 1.51332117219235E-02_wp, 3.21310649550766E-04_wp,-6.50182407583915E-05_wp, &
      & 5.65537994689618E-04_wp, 2.39816673197981E-04_wp,-4.01145908943552E-04_wp, &
      & 1.02979125404226E-03_wp,-3.83960703078053E-04_wp,-4.71112289725868E-05_wp, &
      & 1.02053886193372E-03_wp,-1.37787803489033E-03_wp, 2.14568368478460E-04_wp, &
      &-2.29785457537160E-03_wp, 3.71631530997041E-04_wp, 2.32423304833032E-03_wp, &
      &-2.04580162271832E-03_wp,-5.61485634191517E-04_wp, 5.83275218807750E-04_wp, &
      &-7.52025031397998E-04_wp,-1.82188558275468E-03_wp, 1.03007344334728E-03_wp, &
      & 1.22820555419493E-03_wp, 7.16816433450695E-04_wp, 4.07757821937570E-04_wp, &
      &-2.79625792411043E-04_wp, 2.45027490200433E-03_wp, 1.16049685662805E-04_wp, &
      & 5.03854201101717E-06_wp,-6.44076461307671E-05_wp, 5.12364190514789E-05_wp, &
      & 5.08755648891243E-05_wp,-1.80059796257016E-03_wp, 3.54546470887379E-03_wp, &
      & 3.65950247565885E-03_wp, 2.71148530315096E-03_wp,-6.27125251361419E-06_wp, &
      &-2.32450086452956E-04_wp,-6.58156831460961E-05_wp, 5.28437067393878E-05_wp, &
      & 2.32661759562063E-04_wp,-3.09112108744030E-03_wp,-4.55353025077479E-04_wp, &
      &-2.67566253419294E-02_wp,-1.32002681370673E-02_wp,-1.04518495185120E-04_wp, &
      &-2.53634831707603E-04_wp, 2.39816673197981E-04_wp, 1.01972520181179E-03_wp, &
      & 4.18800952996620E-04_wp, 1.20885516832012E-03_wp,-5.52056921452567E-04_wp, &
      & 5.46257757174947E-04_wp, 5.94007597988317E-04_wp,-2.24225399100554E-03_wp, &
      &-1.19565402301603E-03_wp, 4.74724238071550E-05_wp, 2.60157115829100E-03_wp, &
      & 4.22160666293566E-04_wp,-1.86857557022780E-03_wp,-1.60968883271809E-03_wp, &
      & 4.19280315359991E-04_wp,-1.22301380113751E-03_wp,-3.25003079162081E-05_wp, &
      & 1.53200056222919E-03_wp, 4.56941288207337E-04_wp, 1.18902455312510E-03_wp, &
      &-4.27488088446933E-04_wp, 2.06014446932550E-03_wp, 1.39966952488796E-03_wp, &
      & 5.66604383383437E-05_wp, 9.88365521261470E-05_wp, 4.02013059080886E-05_wp, &
      & 1.03787374695025E-04_wp, 3.70646570213923E-05_wp,-3.51490801006542E-04_wp, &
      &-6.75694770587501E-06_wp, 4.03502323498763E-03_wp, 1.45446242844250E-03_wp, &
      &-7.33810769916618E-05_wp,-1.56880573246973E-04_wp,-1.36918692981879E-04_wp, &
      & 9.07005439053863E-05_wp, 3.40426800813876E-05_wp,-3.26315136987568E-03_wp, &
      &-8.97158346346022E-03_wp,-5.14238323166957E-05_wp,-2.63946101633413E-02_wp, &
      & 1.07469239920213E-05_wp, 1.27233371966009E-04_wp,-4.01145908943552E-04_wp, &
      & 4.18800952996620E-04_wp, 8.91699965008989E-04_wp], shape(density))

   call get_structure(mol, "f-block", "CeCl3")
   call test_acp_qeff_numgrad(error, mol, density, make_qvszp_basis, thr_in=thr1*10)

end subroutine test_qeff_grad_acp_gxtb_cecl3


subroutine test_g_acp_gxtb_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(2, 2, 1) = reshape([&
      & 5.93683766916992E-01_wp, 5.93683766916992E-01_wp, 5.93683766916992E-01_wp, &
      & 5.93683766916992E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "H2")
   call test_acp_numgrad(error, mol, density, make_qvszp_basis, thr_in=thr1)

end subroutine test_g_acp_gxtb_h2

subroutine test_g_acp_gxtb_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(5, 5, 1) = reshape([&
      & 7.43138968868807E-02_wp, 6.30732585440217E-45_wp, 1.15038033099932E-01_wp, &
      & 0.00000000000000E+00_wp, 2.77067464359637E-01_wp, 6.30732585440217E-45_wp, &
      & 5.35328667990135E-88_wp, 9.76375066853567E-45_wp, 0.00000000000000E+00_wp, &
      & 2.35158544306899E-44_wp, 1.15038033099932E-01_wp, 9.76375066853567E-45_wp, &
      & 1.78079062111965E-01_wp, 0.00000000000000E+00_wp, 4.28900884910329E-01_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, &
      & 0.00000000000000E+00_wp, 0.00000000000000E+00_wp, 2.77067464359637E-01_wp, &
      & 2.35158544306899E-44_wp, 4.28900884910329E-01_wp, 0.00000000000000E+00_wp, &
      & 1.03300167293785E+00_wp], shape(density))

   call get_structure(mol, "MB16-43", "LiH")
   call test_acp_numgrad(error, mol, density, make_qvszp_basis, thr_in=thr1)

end subroutine test_g_acp_gxtb_lih

subroutine test_g_acp_gxtb_no(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(8, 8, 2) = reshape([&
      & 9.42009046415725E-01_wp, 8.68128306151435E-17_wp,-3.13530611390051E-01_wp, &
      & 2.83821071259847E-17_wp,-1.91210322106378E-01_wp,-2.95777592933361E-16_wp, &
      &-2.58330200255829E-02_wp,-1.88965537110059E-16_wp, 8.68128306151435E-17_wp, &
      & 7.55421320827554E-01_wp, 4.13035452754259E-17_wp,-4.05630015646006E-01_wp, &
      &-7.56771273088709E-17_wp,-3.58191420308672E-02_wp, 4.29283702464554E-17_wp, &
      & 3.33220482034872E-01_wp,-3.13530611390051E-01_wp, 4.13035452754259E-17_wp, &
      & 4.49934652749067E-01_wp,-1.37354302670472E-16_wp,-1.71924644995845E-03_wp, &
      &-8.68821483997484E-17_wp,-4.07076133805002E-01_wp, 1.81509948244029E-16_wp, &
      & 2.83821071259847E-17_wp,-4.05630015646006E-01_wp,-1.37354302670472E-16_wp, &
      & 5.92359279622604E-01_wp, 1.27369180291386E-16_wp, 3.33220482034872E-01_wp, &
      & 1.10259833116572E-17_wp, 9.81344853538432E-02_wp,-1.91210322106378E-01_wp, &
      &-7.56771273088709E-17_wp,-1.71924644995845E-03_wp, 1.27369180291386E-16_wp, &
      & 1.01808686826282E+00_wp, 1.30483506117222E-16_wp, 2.62723878284128E-01_wp, &
      &-2.46091681044612E-16_wp,-2.95777592933361E-16_wp,-3.58191420308672E-02_wp, &
      &-8.68821483997484E-17_wp, 3.33220482034872E-01_wp, 1.30483506117222E-16_wp, &
      & 8.63442034138714E-01_wp,-5.57245658663308E-17_wp,-2.73736867008512E-01_wp, &
      &-2.58330200255829E-02_wp, 4.29283702464554E-17_wp,-4.07076133805002E-01_wp, &
      & 1.10259833116572E-17_wp, 2.62723878284128E-01_wp,-5.57245658663308E-17_wp, &
      & 5.33778495173348E-01_wp,-1.58550889927256E-17_wp,-1.88965537110059E-16_wp, &
      & 3.33220482034872E-01_wp, 1.81509948244029E-16_wp, 9.81344853538432E-02_wp, &
      &-2.46091681044612E-16_wp,-2.73736867008512E-01_wp,-1.58550889927256E-17_wp, &
      & 7.53400640342027E-01_wp, 9.54305653170233E-01_wp, 1.09464057770057E-17_wp, &
      &-3.11307650277096E-01_wp,-7.36206352015003E-17_wp,-2.03773626333043E-01_wp, &
      & 1.65290626496932E-16_wp,-1.70067070429939E-02_wp,-8.93509590142745E-17_wp, &
      & 1.09464057770057E-17_wp, 2.40282587827900E-01_wp, 4.33933558995002E-17_wp, &
      & 1.85040392733642E-05_wp, 6.74251196637872E-17_wp, 3.64459920248371E-01_wp, &
      & 9.51214979595760E-17_wp, 6.68285106103067E-06_wp,-3.11307650277096E-01_wp, &
      & 4.33933558995002E-17_wp, 4.38731609559950E-01_wp, 2.61613977605231E-17_wp, &
      &-4.10640194707355E-03_wp, 1.10151806812206E-16_wp,-4.09548671902720E-01_wp, &
      & 2.83066560924208E-19_wp,-7.36206352015003E-17_wp, 1.85040392733642E-05_wp, &
      & 2.61613977605231E-17_wp, 2.40290026395803E-01_wp, 3.09134882803496E-17_wp, &
      & 6.68285106153466E-06_wp,-3.66627068485342E-18_wp, 3.64462606734314E-01_wp, &
      &-2.03773626333043E-01_wp, 6.74251196637872E-17_wp,-4.10640194707355E-03_wp, &
      & 3.09134882803496E-17_wp, 1.03092149514987E+00_wp,-3.17330346614473E-16_wp, &
      & 2.53665476333746E-01_wp,-1.60068098192912E-17_wp, 1.65290626496932E-16_wp, &
      & 3.64459920248371E-01_wp, 1.10151806812206E-16_wp, 6.68285106153466E-06_wp, &
      &-3.17330346614473E-16_wp, 5.52811733573677E-01_wp,-3.44555050846798E-17_wp, &
      &-2.22979264938861E-05_wp,-1.70067070429939E-02_wp, 9.51214979595760E-17_wp, &
      &-4.09548671902720E-01_wp,-3.66627068485342E-18_wp, 2.53665476333746E-01_wp, &
      &-3.44555050846798E-17_wp, 5.38687781287961E-01_wp, 2.26914599114829E-17_wp, &
      &-8.93509590142745E-17_wp, 6.68285106103067E-06_wp, 2.83066560924208E-19_wp, &
      & 3.64462606734314E-01_wp,-1.60068098192912E-17_wp,-2.22979264938861E-05_wp, &
      & 2.26914599114829E-17_wp, 5.52802769874570E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "NO")
   call test_acp_numgrad(error, mol, density, make_qvszp_basis, thr_in=thr1*10)

end subroutine test_g_acp_gxtb_no

subroutine test_g_acp_gxtb_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(18, 18, 1) = reshape([&
      & 2.01512274410725E+00_wp, 3.05183510489231E-16_wp,-3.99929219814921E-01_wp, &
      & 7.69080296474742E-16_wp, 1.45947463264597E-04_wp,-4.42069429901860E-17_wp, &
      &-3.51352512626506E-02_wp,-1.08853924469132E-16_wp, 8.19811912191895E-05_wp, &
      &-2.61461503009236E-01_wp,-2.53847513922895E-16_wp, 8.46120257485079E-02_wp, &
      &-1.05417595234967E-15_wp,-2.59372907294474E-05_wp,-4.35889228793221E-17_wp, &
      & 1.18677041085995E-02_wp,-9.23447662047026E-17_wp,-1.45694206907493E-05_wp, &
      & 3.05183510489231E-16_wp, 1.11646346539116E+00_wp, 1.51923309365924E-15_wp, &
      & 5.22665236979214E-01_wp,-2.23623283632863E-17_wp, 2.23688961934518E-02_wp, &
      &-3.01660087631464E-16_wp,-3.85852024970845E-02_wp,-2.43149330627879E-18_wp, &
      &-6.95644355001582E-16_wp, 5.02302414177719E-01_wp, 5.74361021043659E-16_wp, &
      &-5.26740092615372E-01_wp, 7.76809871729795E-17_wp,-7.25879700203994E-02_wp, &
      & 3.46931693965287E-16_wp,-4.72231756427874E-02_wp, 9.75451031609286E-17_wp, &
      &-3.99929219814921E-01_wp, 1.51923309365924E-15_wp, 8.16076844347648E-01_wp, &
      & 1.34007033367206E-15_wp,-2.39395566884248E-05_wp,-7.79787012198157E-17_wp, &
      & 5.39351481157056E-02_wp, 1.48581899014128E-18_wp,-1.34472592440345E-05_wp, &
      &-8.46120257485160E-02_wp,-1.44825907786549E-15_wp,-7.72403989719101E-01_wp, &
      &-1.76921592471880E-15_wp,-1.32742422451715E-07_wp,-1.15813374568686E-16_wp, &
      & 4.74250266390859E-02_wp,-1.82310535840248E-17_wp,-7.45636935369285E-08_wp, &
      & 7.69080296474742E-16_wp, 5.22665236979214E-01_wp, 1.34007033367206E-15_wp, &
      & 1.70364350645588E+00_wp,-8.02983753144685E-18_wp,-3.85852024970846E-02_wp, &
      &-4.70255064846748E-16_wp,-2.09790427620324E-02_wp,-1.96205026639568E-17_wp, &
      &-6.57821241716489E-16_wp,-5.26740092615371E-01_wp, 1.54875272040670E-15_wp, &
      &-8.94554593703680E-02_wp,-1.47010689589666E-16_wp,-4.72231756427873E-02_wp, &
      & 5.33348516140237E-16_wp,-1.25640104065888E-01_wp, 1.48730635567615E-16_wp, &
      & 1.45947463264597E-04_wp,-2.23623283632863E-17_wp,-2.39395566884248E-05_wp, &
      &-8.02983753144685E-18_wp, 1.06235239203714E-08_wp,-6.27312114532195E-19_wp, &
      &-2.27401943747449E-06_wp, 6.00156399076560E-19_wp, 5.96741544153476E-09_wp, &
      &-2.59372907296047E-05_wp,-1.26125559712104E-17_wp, 1.32742422317945E-07_wp, &
      & 8.58050673700442E-18_wp,-2.35397458118598E-09_wp, 1.42171549599656E-18_wp, &
      & 1.27475826808189E-06_wp, 7.55995255055379E-19_wp,-1.32226786236991E-09_wp, &
      &-4.42069429901860E-17_wp, 2.23688961934518E-02_wp,-7.79787012198157E-17_wp, &
      &-3.85852024970846E-02_wp,-6.27312114532195E-19_wp, 4.86266676876034E-03_wp, &
      & 3.34959389043418E-17_wp, 4.04941245032744E-03_wp, 6.35609235074649E-19_wp, &
      & 6.31298797212952E-17_wp, 7.25879700204016E-02_wp,-7.29477821938508E-17_wp, &
      & 4.72231756427944E-02_wp, 1.03467976744611E-17_wp,-7.56324676480364E-04_wp, &
      &-1.77577932859663E-17_wp, 2.96701633209598E-03_wp,-6.46245268795547E-18_wp, &
      &-3.51352512626506E-02_wp,-3.01660087631464E-16_wp, 5.39351481157056E-02_wp, &
      &-4.70255064846748E-16_wp,-2.27401943747449E-06_wp, 3.34959389043418E-17_wp, &
      & 3.73731424127451E-03_wp, 7.17979015797382E-17_wp,-1.27735568790833E-06_wp, &
      & 1.18677041086054E-02_wp, 3.72567464311869E-16_wp,-4.74250266390819E-02_wp, &
      & 8.06318001452668E-16_wp, 1.27475826807370E-06_wp, 2.50819963630132E-17_wp, &
      & 2.76687889875194E-03_wp, 4.75201425698955E-17_wp, 7.16053564725540E-07_wp, &
      &-1.08853924469132E-16_wp,-3.85852024970845E-02_wp, 1.48581899014128E-18_wp, &
      &-2.09790427620324E-02_wp, 6.00156399076560E-19_wp, 4.04941245032744E-03_wp, &
      & 7.17979015797382E-17_wp, 9.41191550512845E-03_wp, 2.27476146379299E-19_wp, &
      &-3.86583716090084E-17_wp, 4.72231756427946E-02_wp,-3.15405023331482E-16_wp, &
      & 1.25640104065898E-01_wp, 2.15478167900987E-18_wp, 2.96701633209597E-03_wp, &
      &-1.89192522936108E-17_wp, 2.57692316259317E-03_wp,-1.20436763189445E-17_wp, &
      & 8.19811912191895E-05_wp,-2.43149330627879E-18_wp,-1.34472592440345E-05_wp, &
      &-1.96205026639568E-17_wp, 5.96741544153476E-09_wp, 6.35609235074649E-19_wp, &
      &-1.27735568790833E-06_wp, 2.27476146379299E-19_wp, 3.35199951718306E-09_wp, &
      &-1.45694206906697E-05_wp, 9.40823506771232E-18_wp, 7.45636930694748E-08_wp, &
      & 5.82344702349143E-19_wp,-1.32226786235450E-09_wp, 3.30923532167611E-19_wp, &
      & 7.16053564753441E-07_wp, 1.42332895146119E-18_wp,-7.42740518019966E-10_wp, &
      &-2.61461503009236E-01_wp,-6.95644355001582E-16_wp,-8.46120257485160E-02_wp, &
      &-6.57821241716489E-16_wp,-2.59372907296047E-05_wp, 6.31298797212952E-17_wp, &
      & 1.18677041086054E-02_wp,-3.86583716090084E-17_wp,-1.45694206906697E-05_wp, &
      & 2.01512274410725E+00_wp, 6.37025174552147E-16_wp, 3.99929219814931E-01_wp, &
      & 8.33671871420983E-16_wp, 1.45947463264675E-04_wp, 5.98740820853847E-17_wp, &
      &-3.51352512626559E-02_wp, 6.97851591164381E-17_wp, 8.19811912194447E-05_wp, &
      &-2.53847513922895E-16_wp, 5.02302414177719E-01_wp,-1.44825907786549E-15_wp, &
      &-5.26740092615371E-01_wp,-1.26125559712104E-17_wp, 7.25879700204016E-02_wp, &
      & 3.72567464311869E-16_wp, 4.72231756427946E-02_wp, 9.40823506771232E-18_wp, &
      & 6.37025174552147E-16_wp, 1.11646346539112E+00_wp,-3.83819455107282E-16_wp, &
      & 5.22665236979232E-01_wp, 1.65728863182037E-16_wp,-2.23688961934538E-02_wp, &
      &-2.35707388589142E-16_wp, 3.85852024970803E-02_wp,-7.60309218403261E-17_wp, &
      & 8.46120257485079E-02_wp, 5.74361021043659E-16_wp,-7.72403989719101E-01_wp, &
      & 1.54875272040670E-15_wp, 1.32742422317945E-07_wp,-7.29477821938508E-17_wp, &
      &-4.74250266390819E-02_wp,-3.15405023331482E-16_wp, 7.45636930694748E-08_wp, &
      & 3.99929219814931E-01_wp,-3.83819455107282E-16_wp, 8.16076844347611E-01_wp, &
      &-2.03826552101056E-15_wp, 2.39395566885481E-05_wp,-4.64841795882909E-17_wp, &
      &-5.39351481157069E-02_wp,-2.20759476011753E-16_wp, 1.34472592445486E-05_wp, &
      &-1.05417595234967E-15_wp,-5.26740092615372E-01_wp,-1.76921592471880E-15_wp, &
      &-8.94554593703680E-02_wp, 8.58050673700442E-18_wp, 4.72231756427944E-02_wp, &
      & 8.06318001452668E-16_wp, 1.25640104065898E-01_wp, 5.82344702349143E-19_wp, &
      & 8.33671871420983E-16_wp, 5.22665236979232E-01_wp,-2.03826552101056E-15_wp, &
      & 1.70364350645586E+00_wp, 3.39236205107499E-18_wp, 3.85852024970803E-02_wp, &
      &-3.39364583169101E-16_wp, 2.09790427620258E-02_wp,-1.47879398954886E-16_wp, &
      &-2.59372907294474E-05_wp, 7.76809871729795E-17_wp,-1.32742422451715E-07_wp, &
      &-1.47010689589666E-16_wp,-2.35397458118598E-09_wp, 1.03467976744611E-17_wp, &
      & 1.27475826807370E-06_wp, 2.15478167900987E-18_wp,-1.32226786235450E-09_wp, &
      & 1.45947463264675E-04_wp, 1.65728863182037E-16_wp, 2.39395566885481E-05_wp, &
      & 3.39236205107499E-18_wp, 1.06235239203804E-08_wp,-3.14583296736217E-18_wp, &
      &-2.27401943748182E-06_wp, 1.01384858020498E-17_wp, 5.96741544155305E-09_wp, &
      &-4.35889228793221E-17_wp,-7.25879700203994E-02_wp,-1.15813374568686E-16_wp, &
      &-4.72231756427873E-02_wp, 1.42171549599656E-18_wp,-7.56324676480364E-04_wp, &
      & 2.50819963630132E-17_wp, 2.96701633209597E-03_wp, 3.30923532167611E-19_wp, &
      & 5.98740820853847E-17_wp,-2.23688961934538E-02_wp,-4.64841795882909E-17_wp, &
      & 3.85852024970803E-02_wp,-3.14583296736217E-18_wp, 4.86266676875980E-03_wp, &
      &-2.75950952491621E-17_wp, 4.04941245032639E-03_wp,-7.72981644699408E-18_wp, &
      & 1.18677041085995E-02_wp, 3.46931693965287E-16_wp, 4.74250266390859E-02_wp, &
      & 5.33348516140237E-16_wp, 1.27475826808189E-06_wp,-1.77577932859663E-17_wp, &
      & 2.76687889875194E-03_wp,-1.89192522936108E-17_wp, 7.16053564753441E-07_wp, &
      &-3.51352512626559E-02_wp,-2.35707388589142E-16_wp,-5.39351481157069E-02_wp, &
      &-3.39364583169101E-16_wp,-2.27401943748182E-06_wp,-2.75950952491621E-17_wp, &
      & 3.73731424127476E-03_wp,-3.67074117602398E-17_wp,-1.27735568794358E-06_wp, &
      &-9.23447662047026E-17_wp,-4.72231756427874E-02_wp,-1.82310535840248E-17_wp, &
      &-1.25640104065888E-01_wp, 7.55995255055379E-19_wp, 2.96701633209598E-03_wp, &
      & 4.75201425698955E-17_wp, 2.57692316259317E-03_wp, 1.42332895146119E-18_wp, &
      & 6.97851591164381E-17_wp, 3.85852024970803E-02_wp,-2.20759476011753E-16_wp, &
      & 2.09790427620258E-02_wp, 1.01384858020498E-17_wp, 4.04941245032639E-03_wp, &
      &-3.67074117602398E-17_wp, 9.41191550512677E-03_wp,-1.22124983920611E-17_wp, &
      &-1.45694206907493E-05_wp, 9.75451031609286E-17_wp,-7.45636935369285E-08_wp, &
      & 1.48730635567615E-16_wp,-1.32226786236991E-09_wp,-6.46245268795547E-18_wp, &
      & 7.16053564725540E-07_wp,-1.20436763189445E-17_wp,-7.42740518019966E-10_wp, &
      & 8.19811912194447E-05_wp,-7.60309218403261E-17_wp, 1.34472592445486E-05_wp, &
      &-1.47879398954886E-16_wp, 5.96741544155305E-09_wp,-7.72981644699408E-18_wp, &
      &-1.27735568794358E-06_wp,-1.22124983920611E-17_wp, 3.35199951720075E-09_wp],&
      & shape(density))

   call get_structure(mol, "MB16-43", "S2")
   call test_acp_numgrad(error, mol, density, make_qvszp_basis, thr_in=thr1*10)

end subroutine test_g_acp_gxtb_s2

subroutine test_g_acp_gxtb_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(13, 13, 1) = reshape([&
      & 6.67076596984081E-01_wp, 1.17812091373782E-15_wp, 1.16003518499611E-15_wp, &
      &-9.41177582170035E-16_wp,-5.18269183303055E-17_wp, 2.81268188491212E-17_wp, &
      & 1.09658020136346E-17_wp,-2.39914571721069E-16_wp, 2.14126324730246E-19_wp, &
      & 2.51982366144073E-01_wp, 2.51982366144073E-01_wp, 2.51982366144073E-01_wp, &
      & 2.51982366144077E-01_wp, 1.17812091373782E-15_wp, 4.65632832545944E-01_wp, &
      &-3.08027819846221E-16_wp,-5.01335488096171E-17_wp, 5.51922671972401E-16_wp, &
      &-4.45360639961878E-16_wp,-7.86723569094682E-18_wp,-5.58630716839440E-02_wp, &
      &-7.00845070605205E-17_wp, 3.05808254914735E-01_wp,-3.05808254914736E-01_wp, &
      &-3.05808254914736E-01_wp, 3.05808254914735E-01_wp, 1.16003518499611E-15_wp, &
      &-3.08027819846221E-16_wp, 4.65632832545945E-01_wp,-3.98576615792397E-17_wp, &
      &-5.58630716839442E-02_wp,-3.56963383436534E-16_wp, 9.89192849362209E-17_wp, &
      & 3.11962003002779E-16_wp, 1.36197924766342E-16_wp,-3.05808254914736E-01_wp, &
      &-3.05808254914736E-01_wp, 3.05808254914735E-01_wp, 3.05808254914735E-01_wp, &
      &-9.41177582170035E-16_wp,-5.01335488096171E-17_wp,-3.98576615792397E-17_wp, &
      & 4.65632832545944E-01_wp,-4.07263775340787E-16_wp,-5.58630716839439E-02_wp, &
      & 5.44186509945490E-17_wp,-5.04840605943758E-16_wp, 1.59146606601200E-17_wp, &
      & 3.05808254914736E-01_wp,-3.05808254914735E-01_wp, 3.05808254914736E-01_wp, &
      &-3.05808254914735E-01_wp,-5.18269183303055E-17_wp, 5.51922671972401E-16_wp, &
      &-5.58630716839442E-02_wp,-4.07263775340787E-16_wp, 6.70202477111097E-03_wp, &
      & 9.22921619598948E-17_wp,-1.18675804605582E-17_wp,-9.93121625154311E-17_wp, &
      &-1.63399870082742E-17_wp, 3.66885392776897E-02_wp, 3.66885392776896E-02_wp, &
      &-3.66885392776901E-02_wp,-3.66885392776889E-02_wp, 2.81268188491212E-17_wp, &
      &-4.45360639961878E-16_wp,-3.56963383436534E-16_wp,-5.58630716839439E-02_wp, &
      & 9.22921619598948E-17_wp, 6.70202477111090E-03_wp,-6.52873420637083E-18_wp, &
      & 1.14654961676060E-16_wp,-1.90931946190504E-18_wp,-3.66885392776896E-02_wp, &
      & 3.66885392776898E-02_wp,-3.66885392776895E-02_wp, 3.66885392776888E-02_wp, &
      & 1.09658020136346E-17_wp,-7.86723569094682E-18_wp, 9.89192849362209E-17_wp, &
      & 5.44186509945490E-17_wp,-1.18675804605582E-17_wp,-6.52873420637083E-18_wp, &
      & 2.76875750299656E-32_wp, 9.43850864112936E-19_wp, 3.19815635497338E-32_wp, &
      &-3.02508000498458E-17_wp,-9.13968622807606E-17_wp, 1.10015080012159E-16_wp, &
      & 2.82015265779562E-17_wp,-2.39914571721069E-16_wp,-5.58630716839440E-02_wp, &
      & 3.11962003002779E-16_wp,-5.04840605943758E-16_wp,-9.93121625154311E-17_wp, &
      & 1.14654961676060E-16_wp, 9.43850864112936E-19_wp, 6.70202477111093E-03_wp, &
      & 8.40820399293783E-18_wp,-3.66885392776899E-02_wp, 3.66885392776897E-02_wp, &
      & 3.66885392776894E-02_wp,-3.66885392776889E-02_wp, 2.14126324730246E-19_wp, &
      &-7.00845070605205E-17_wp, 1.36197924766342E-16_wp, 1.59146606601200E-17_wp, &
      &-1.63399870082742E-17_wp,-1.90931946190504E-18_wp, 3.19815635497338E-32_wp, &
      & 8.40820399293783E-18_wp, 5.09307325669261E-32_wp,-1.24944740611828E-16_wp, &
      &-5.37917417731081E-17_wp, 1.46010682147439E-16_wp, 3.30493376435646E-17_wp, &
      & 2.51982366144073E-01_wp, 3.05808254914735E-01_wp,-3.05808254914736E-01_wp, &
      & 3.05808254914736E-01_wp, 3.66885392776897E-02_wp,-3.66885392776896E-02_wp, &
      &-3.02508000498458E-17_wp,-3.66885392776899E-02_wp,-1.24944740611828E-16_wp, &
      & 6.97710523806227E-01_wp,-1.05657986637187E-01_wp,-1.05657986637186E-01_wp, &
      &-1.05657986637189E-01_wp, 2.51982366144073E-01_wp,-3.05808254914736E-01_wp, &
      &-3.05808254914736E-01_wp,-3.05808254914735E-01_wp, 3.66885392776896E-02_wp, &
      & 3.66885392776898E-02_wp,-9.13968622807606E-17_wp, 3.66885392776897E-02_wp, &
      &-5.37917417731081E-17_wp,-1.05657986637187E-01_wp, 6.97710523806228E-01_wp, &
      &-1.05657986637187E-01_wp,-1.05657986637188E-01_wp, 2.51982366144073E-01_wp, &
      &-3.05808254914736E-01_wp, 3.05808254914735E-01_wp, 3.05808254914736E-01_wp, &
      &-3.66885392776901E-02_wp,-3.66885392776895E-02_wp, 1.10015080012159E-16_wp, &
      & 3.66885392776894E-02_wp, 1.46010682147439E-16_wp,-1.05657986637186E-01_wp, &
      &-1.05657986637187E-01_wp, 6.97710523806227E-01_wp,-1.05657986637189E-01_wp, &
      & 2.51982366144077E-01_wp, 3.05808254914735E-01_wp, 3.05808254914735E-01_wp, &
      &-3.05808254914735E-01_wp,-3.66885392776889E-02_wp, 3.66885392776888E-02_wp, &
      & 2.82015265779562E-17_wp,-3.66885392776889E-02_wp, 3.30493376435646E-17_wp, &
      &-1.05657986637189E-01_wp,-1.05657986637188E-01_wp,-1.05657986637189E-01_wp, &
      & 6.97710523806226E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "SiH4")
   call test_acp_numgrad(error, mol, density, make_qvszp_basis, thr_in=thr1*10)

end subroutine test_g_acp_gxtb_sih4

subroutine test_g_acp_gxtb_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: density(43, 43, 2) = reshape([&
      & 6.91475853693640E-03_wp,-4.66439474613685E-04_wp, 7.43124629353357E-04_wp, &
      &-3.09904488027013E-04_wp,-2.72087829478411E-03_wp, 6.74488694981665E-03_wp, &
      &-5.28266720158637E-03_wp, 4.80689218346076E-03_wp, 9.22601160310286E-04_wp, &
      &-2.71966921742355E-03_wp,-1.01894320855533E-02_wp,-2.27063129413855E-03_wp, &
      &-2.82181425185299E-03_wp,-2.14613395036928E-03_wp, 3.61217559584011E-03_wp, &
      & 1.12028617141975E-02_wp,-8.64033744991843E-04_wp, 3.10371374558179E-02_wp, &
      & 2.30333249960052E-02_wp, 2.04795146556047E-02_wp, 1.17566904494905E-03_wp, &
      & 1.12383332749237E-03_wp,-3.37161261865593E-04_wp, 7.55695280200463E-04_wp, &
      &-4.91653601592628E-04_wp,-1.02695166967481E-03_wp,-3.45289149212513E-02_wp, &
      &-1.34121866164477E-02_wp, 2.17842958083843E-02_wp,-1.32239449246803E-03_wp, &
      & 9.23763188397653E-04_wp,-5.87060101040320E-04_wp,-6.06266471532816E-04_wp, &
      &-5.25957119780333E-04_wp,-9.44887800903170E-04_wp, 9.71328579477840E-03_wp, &
      &-1.32264514234976E-02_wp,-4.00506954182355E-02_wp,-7.30975427202840E-04_wp, &
      &-2.86085329176629E-04_wp,-5.83829373771101E-04_wp, 1.06473924707788E-03_wp, &
      & 1.20864032141500E-03_wp,-4.66439474613685E-04_wp, 2.37196189260911E-03_wp, &
      & 3.76291798633030E-04_wp,-1.89524702707700E-04_wp,-2.70383567603791E-03_wp, &
      & 7.51687695827497E-04_wp, 1.57852557823706E-05_wp,-1.70421269429076E-03_wp, &
      &-1.24711040742963E-04_wp, 7.89473512438832E-04_wp,-3.57098061627412E-04_wp, &
      & 4.01348859464832E-03_wp,-3.38878489221253E-03_wp, 6.34879336473785E-03_wp, &
      &-2.30645681562275E-03_wp,-6.70882068405559E-04_wp, 6.96944723684135E-03_wp, &
      &-1.08436093353191E-02_wp,-1.77208852888238E-02_wp,-1.93565751770573E-02_wp, &
      &-8.75401591254688E-04_wp,-5.65296230452197E-04_wp, 1.94865564589505E-04_wp, &
      &-6.09542770971772E-04_wp,-6.86637147631842E-05_wp,-7.77156635669613E-03_wp, &
      &-1.12023415897402E-02_wp,-1.78971328952172E-02_wp, 2.23980508202744E-02_wp, &
      &-8.58573131874772E-04_wp, 6.50886246444007E-04_wp,-1.23589021153629E-04_wp, &
      &-7.14602097204805E-04_wp, 1.61848868701550E-04_wp, 2.47681579669712E-03_wp, &
      & 1.51505155319784E-02_wp, 4.13303227392449E-03_wp, 1.19937863232238E-02_wp, &
      &-4.00974483671558E-04_wp,-1.44661987626867E-04_wp,-2.20905355403101E-05_wp, &
      &-4.45102163842215E-04_wp,-5.21718641127864E-04_wp, 7.43124629353357E-04_wp, &
      & 3.76291798633030E-04_wp, 1.92419514812669E-03_wp, 2.43411943306366E-04_wp, &
      &-1.38558534277301E-03_wp,-1.29700970256260E-03_wp, 3.20274034780636E-04_wp, &
      &-6.00771636380427E-04_wp, 4.53126654298869E-04_wp, 6.89363513179920E-04_wp, &
      & 2.85458472275171E-03_wp,-8.42769400830715E-03_wp,-3.49285659188215E-03_wp, &
      &-8.73740458812932E-03_wp,-2.03180597856797E-03_wp, 3.57585309496399E-03_wp, &
      & 4.50241444891910E-03_wp,-1.94202483147024E-02_wp, 1.08132873675751E-02_wp, &
      &-1.31156713245789E-02_wp,-6.35109734804627E-04_wp,-1.41358099309022E-04_wp, &
      & 4.86112152718710E-04_wp,-6.63799363359947E-05_wp, 2.63148272706973E-04_wp, &
      &-3.69338959899031E-03_wp,-1.62131073853426E-02_wp, 1.12706765571304E-02_wp, &
      & 1.05297118737512E-02_wp,-6.83652635175080E-04_wp,-3.67543219007418E-06_wp, &
      &-5.64274993106786E-04_wp,-4.99490489386992E-05_wp,-1.70869611258782E-04_wp, &
      &-3.50216085250586E-03_wp, 4.54245265617874E-03_wp, 1.19055309005284E-02_wp, &
      &-1.89552860619114E-02_wp,-4.21614601949817E-04_wp,-8.16969303097783E-05_wp, &
      &-5.05453807145344E-04_wp, 3.80106625787291E-05_wp, 4.70419178269535E-04_wp, &
      &-3.09904488027013E-04_wp,-1.89524702707700E-04_wp, 2.43411943306366E-04_wp, &
      & 2.47928402744502E-03_wp,-1.13559737079954E-03_wp,-1.51157026310338E-03_wp, &
      &-1.74075334951878E-05_wp, 1.78516116698151E-03_wp, 2.08500104756173E-03_wp, &
      &-2.62166403396250E-03_wp, 1.57889816852864E-03_wp, 5.04088894829648E-03_wp, &
      &-2.33600703356005E-03_wp, 1.82408462816258E-03_wp, 3.82256307533336E-03_wp, &
      &-1.43702087113090E-03_wp, 4.61273866776342E-03_wp,-1.97237097978908E-02_wp, &
      &-1.20991357049096E-02_wp, 5.08336438917607E-03_wp,-2.83097626903721E-04_wp, &
      &-5.91992885184355E-04_wp, 1.38111314514985E-04_wp,-9.22681999257671E-05_wp, &
      & 7.37417254726881E-04_wp, 5.31578919345793E-03_wp, 2.10073897134136E-02_wp, &
      & 1.09429414011269E-02_wp, 3.20167044114457E-03_wp, 3.34046150737528E-04_wp, &
      &-7.19136030754418E-04_wp, 4.42987781466214E-05_wp, 2.30642331218448E-04_wp, &
      & 6.80365098404692E-04_wp,-8.84877256219468E-03_wp, 1.35788864585994E-02_wp, &
      &-2.03601998000453E-02_wp,-2.27788852326930E-02_wp,-6.78210598851248E-04_wp, &
      &-4.08503618156882E-04_wp,-1.59367877771823E-04_wp, 1.00614012061151E-03_wp, &
      & 6.42501172500844E-04_wp,-2.72087829478411E-03_wp,-2.70383567603791E-03_wp, &
      &-1.38558534277301E-03_wp,-1.13559737079954E-03_wp, 4.40223717129586E-02_wp, &
      & 5.57984043828882E-03_wp, 2.63086726948487E-03_wp, 7.70835864593437E-03_wp, &
      &-3.97661022352564E-04_wp, 2.26637457646800E-04_wp, 7.74500745277695E-03_wp, &
      &-1.79565826907083E-03_wp, 1.67084160854124E-02_wp,-6.61485609356795E-03_wp, &
      &-7.67995123149824E-04_wp, 2.18809239728138E-03_wp, 2.24397055345458E-02_wp, &
      & 6.64549541657225E-02_wp, 8.85871421287740E-02_wp, 6.96854518071820E-04_wp, &
      & 9.90846322150787E-04_wp, 2.46975374061967E-03_wp, 3.12047632637008E-04_wp, &
      & 1.29596421478787E-03_wp,-1.29311438090641E-03_wp,-2.40673203813859E-02_wp, &
      & 7.93390435947526E-02_wp, 7.02903965683136E-02_wp,-1.26672314141543E-02_wp, &
      & 1.72311948156705E-03_wp,-2.33142012668874E-03_wp, 1.41048483944506E-04_wp, &
      & 1.23754087718086E-03_wp, 1.32235237491202E-03_wp,-1.35903661755520E-02_wp, &
      & 7.65937840551549E-02_wp, 3.28686383454786E-02_wp, 6.44067474820537E-02_wp, &
      &-1.00233070466989E-03_wp,-6.88819688386836E-04_wp, 3.53857312558554E-04_wp, &
      &-1.40694498188145E-03_wp,-2.36540545013052E-03_wp, 6.74488694981665E-03_wp, &
      & 7.51687695827497E-04_wp,-1.29700970256260E-03_wp,-1.51157026310338E-03_wp, &
      & 5.57984043828882E-03_wp, 3.14464127057532E-02_wp, 5.42414731300381E-03_wp, &
      &-4.96716029907383E-03_wp,-7.32631149808548E-03_wp, 3.52486701113476E-03_wp, &
      &-1.72853100046916E-02_wp, 6.66597965239516E-03_wp, 1.41280628812287E-03_wp, &
      & 1.69880913473563E-03_wp, 6.38096027411938E-04_wp, 9.57106430383583E-03_wp, &
      & 2.07478739672486E-02_wp, 6.75325136801770E-02_wp,-2.15309610649779E-02_wp, &
      & 7.89233538672891E-02_wp, 2.81701658165517E-03_wp, 5.88757207551287E-04_wp, &
      &-1.82030938336181E-03_wp, 5.97326378130184E-04_wp,-5.00617661066083E-04_wp, &
      & 1.84804474365555E-02_wp,-6.44185602640536E-02_wp, 3.79320344196822E-02_wp, &
      & 6.87221148902045E-02_wp,-2.66981711043062E-03_wp, 3.40407076461401E-04_wp, &
      &-2.12977682744976E-03_wp,-3.25448102001268E-04_wp,-6.09023766965693E-04_wp, &
      &-5.22282280974187E-03_wp, 3.75844989308166E-02_wp,-1.69469206157941E-02_wp, &
      & 4.36814617080786E-02_wp,-7.35276572443060E-04_wp,-6.06690581189437E-04_wp, &
      & 6.06034207884786E-04_wp, 1.10660155876927E-04_wp,-1.13219077295543E-03_wp, &
      &-5.28266720158637E-03_wp, 1.57852557823706E-05_wp, 3.20274034780636E-04_wp, &
      &-1.74075334951878E-05_wp, 2.63086726948487E-03_wp, 5.42414731300381E-03_wp, &
      & 2.93589728175265E-02_wp, 3.67709738831447E-03_wp,-7.47063038603620E-04_wp, &
      & 3.90836704630479E-03_wp, 7.09621952656464E-05_wp, 2.07760552467874E-03_wp, &
      &-1.29405483046609E-04_wp, 2.89718143673136E-03_wp, 3.16800401176324E-04_wp, &
      &-1.69606053549689E-02_wp,-7.54471811196526E-03_wp, 1.55318190585747E-02_wp, &
      &-9.64076933204909E-02_wp, 1.10515854057054E-02_wp, 1.45947727313188E-04_wp, &
      &-2.16103631752938E-03_wp,-1.41298732458204E-03_wp,-1.48077115538771E-03_wp, &
      &-3.68872484730471E-05_wp,-9.73129495503449E-03_wp, 3.62507400525617E-04_wp, &
      & 8.64415228683305E-02_wp,-5.52582904665705E-03_wp, 2.36297575249206E-04_wp, &
      &-1.94458556002641E-03_wp,-1.28072235876254E-03_wp, 1.67826970775080E-03_wp, &
      &-2.64680965053990E-04_wp,-9.91183630949829E-03_wp,-5.57175718616085E-03_wp, &
      & 8.67153633545441E-02_wp, 2.90410287476276E-03_wp, 3.78053496299328E-04_wp, &
      & 9.38622543054151E-04_wp,-1.22901953593365E-03_wp,-2.36862461809501E-03_wp, &
      & 3.42392873712068E-05_wp, 4.80689218346076E-03_wp,-1.70421269429076E-03_wp, &
      &-6.00771636380427E-04_wp, 1.78516116698151E-03_wp, 7.70835864593437E-03_wp, &
      &-4.96716029907383E-03_wp, 3.67709738831447E-03_wp, 3.44193006387433E-02_wp, &
      & 4.76228788479840E-03_wp,-1.12594099043911E-02_wp,-4.22637413690153E-03_wp, &
      &-1.10988128220400E-02_wp, 9.39053014236173E-04_wp,-1.16060082172845E-03_wp, &
      & 1.07624737460322E-02_wp, 5.60381954566253E-03_wp, 1.41316668863599E-02_wp, &
      & 7.88193671551557E-02_wp,-1.33489350096091E-02_wp, 5.13853337547530E-03_wp, &
      & 1.55817983246274E-03_wp, 5.12251036343349E-04_wp,-1.38894336044170E-03_wp, &
      & 4.98434304022621E-06_wp,-1.55155659357652E-03_wp,-1.19435034444835E-02_wp, &
      & 7.87707633106604E-02_wp,-3.00272463633091E-02_wp,-2.09804044214836E-03_wp, &
      & 1.14598547458529E-03_wp,-1.69273880680804E-04_wp, 1.52074117567880E-03_wp, &
      &-6.69143343856796E-04_wp, 1.75734486260590E-03_wp, 2.14421083889722E-02_wp, &
      & 3.28723251997694E-02_wp, 4.78972976105153E-02_wp,-1.04712017981954E-01_wp, &
      &-1.59454356614859E-03_wp,-7.29860323367932E-05_wp,-2.39582887843716E-03_wp, &
      & 4.78040749240586E-04_wp, 2.59547000700204E-03_wp, 9.22601160310286E-04_wp, &
      &-1.24711040742963E-04_wp, 4.53126654298869E-04_wp, 2.08500104756173E-03_wp, &
      &-3.97661022352564E-04_wp,-7.32631149808548E-03_wp,-7.47063038603620E-04_wp, &
      & 4.76228788479840E-03_wp, 4.33364941711815E-02_wp,-3.15071146362812E-03_wp, &
      &-6.67108487319752E-04_wp, 6.99176018534552E-03_wp,-6.22633014911112E-03_wp, &
      &-3.38055775952816E-03_wp, 7.69017399473797E-03_wp,-1.76961128934806E-03_wp, &
      &-8.85374052511765E-03_wp, 3.47268154939371E-02_wp,-3.45265330160401E-02_wp, &
      &-8.92468896301965E-02_wp,-1.30909318966849E-03_wp,-3.70938201881427E-04_wp, &
      &-1.10750429754854E-04_wp,-1.31927484389760E-03_wp,-1.78842876359272E-03_wp, &
      &-8.63307150232750E-03_wp,-3.79182506231862E-02_wp, 3.52704972597678E-02_wp, &
      &-9.70603235678339E-02_wp, 1.09039808896316E-03_wp,-3.59947499872617E-04_wp, &
      &-3.60824594442789E-04_wp, 1.85106237168900E-03_wp,-2.29662810010202E-03_wp, &
      & 2.20191489255154E-02_wp, 7.33216018913181E-02_wp,-7.19861577004959E-02_wp, &
      &-4.69373915752494E-02_wp,-2.10917504939881E-03_wp,-1.45840439450965E-03_wp, &
      & 1.08123286917808E-04_wp, 2.31640389941183E-03_wp, 5.16614467009417E-04_wp, &
      &-2.71966921742355E-03_wp, 7.89473512438832E-04_wp, 6.89363513179920E-04_wp, &
      &-2.62166403396250E-03_wp, 2.26637457646800E-04_wp, 3.52486701113476E-03_wp, &
      & 3.90836704630479E-03_wp,-1.12594099043911E-02_wp,-3.15071146362812E-03_wp, &
      & 2.34151430731864E-02_wp, 1.85415808477836E-02_wp,-4.24734078204873E-02_wp, &
      & 2.35575691065080E-03_wp,-3.33058390368733E-02_wp,-5.57795006981639E-03_wp, &
      & 1.01108773593490E-02_wp, 4.01102854511713E-04_wp,-4.56200941512203E-02_wp, &
      &-6.60382781523261E-03_wp, 4.54256036670506E-02_wp, 8.13365306994022E-04_wp, &
      &-1.45138153891137E-03_wp,-4.61944048276185E-04_wp, 2.51862372218516E-04_wp, &
      & 2.07065310902566E-03_wp, 3.41550940627179E-04_wp,-5.80400261978655E-02_wp, &
      &-1.74917522456011E-02_wp,-5.30554321998536E-02_wp, 5.62332461341726E-04_wp, &
      & 1.78951892334922E-03_wp, 5.49306638141359E-04_wp, 1.18354420618126E-04_wp, &
      &-2.60744135287044E-03_wp,-1.84935397407614E-04_wp, 4.00823788730297E-02_wp, &
      & 5.07256361601547E-02_wp, 6.01728319805672E-02_wp,-1.08074030232713E-03_wp, &
      &-5.41519705134412E-04_wp, 1.74128670059472E-04_wp,-2.34898043486087E-03_wp, &
      &-2.72927683099535E-03_wp,-1.01894320855533E-02_wp,-3.57098061627412E-04_wp, &
      & 2.85458472275171E-03_wp, 1.57889816852864E-03_wp, 7.74500745277695E-03_wp, &
      &-1.72853100046916E-02_wp, 7.09621952656464E-05_wp,-4.22637413690153E-03_wp, &
      &-6.67108487319752E-04_wp, 1.85415808477836E-02_wp, 7.81937062031174E-02_wp, &
      &-1.56952868345472E-01_wp, 4.63674964411300E-03_wp,-1.52201455602624E-01_wp, &
      &-3.58517678616803E-02_wp, 3.07585591512265E-02_wp, 4.52883325072763E-04_wp, &
      &-9.32147819443877E-02_wp,-5.48349175380994E-03_wp,-3.42959429035476E-02_wp, &
      &-1.17137609070403E-03_wp,-4.30791961546821E-03_wp,-8.16503856575622E-04_wp, &
      &-1.99648246109626E-03_wp, 1.66375144799501E-03_wp, 7.55268864056282E-05_wp, &
      & 7.36998183172769E-02_wp, 8.59929205900434E-03_wp,-3.29569302511059E-02_wp, &
      & 2.35295984835294E-03_wp,-7.98432581416592E-04_wp, 2.29759184132251E-03_wp, &
      &-6.70858028794763E-04_wp, 2.19545673464568E-03_wp, 9.16685372743052E-04_wp, &
      & 2.57126867846738E-02_wp, 1.14966740670106E-02_wp, 3.88665260783668E-02_wp, &
      &-9.20983705270224E-04_wp,-1.88651521181885E-03_wp, 2.10285595596077E-03_wp, &
      & 5.97064769894395E-04_wp,-3.32760114542134E-03_wp,-2.27063129413855E-03_wp, &
      & 4.01348859464832E-03_wp,-8.42769400830715E-03_wp, 5.04088894829648E-03_wp, &
      &-1.79565826907083E-03_wp, 6.66597965239516E-03_wp, 2.07760552467874E-03_wp, &
      &-1.10988128220400E-02_wp, 6.99176018534552E-03_wp,-4.24734078204873E-02_wp, &
      &-1.56952868345472E-01_wp, 4.62813275841273E-01_wp,-9.10464660774819E-03_wp, &
      & 4.38047333432606E-01_wp, 9.59547188118657E-02_wp,-1.25382771228594E-01_wp, &
      &-1.38374238093702E-03_wp,-7.41273093138561E-03_wp,-1.10468068244529E-02_wp, &
      &-1.56187513126319E-03_wp,-4.41907081296685E-03_wp, 6.66432321276951E-03_wp, &
      & 6.93367015797582E-03_wp, 3.50242290466801E-03_wp, 1.53391697941402E-03_wp, &
      & 5.43354709716936E-04_wp, 1.40237003812439E-02_wp, 3.28804253492689E-02_wp, &
      & 1.21189122316512E-02_wp,-6.66975469219293E-05_wp,-3.13907873268688E-03_wp, &
      &-3.25713971009019E-03_wp, 4.31017066296889E-03_wp,-1.47887531130849E-03_wp, &
      &-1.84087838689988E-04_wp, 1.41827982066621E-02_wp,-6.62740370087469E-02_wp, &
      & 2.39879765902584E-02_wp, 1.05566543856893E-03_wp, 3.60304872454716E-03_wp, &
      &-3.54139161540637E-03_wp,-2.95929713093146E-03_wp, 4.28984898113645E-03_wp, &
      &-2.82181425185299E-03_wp,-3.38878489221253E-03_wp,-3.49285659188215E-03_wp, &
      &-2.33600703356005E-03_wp, 1.67084160854124E-02_wp, 1.41280628812287E-03_wp, &
      &-1.29405483046609E-04_wp, 9.39053014236173E-04_wp,-6.22633014911112E-03_wp, &
      & 2.35575691065080E-03_wp, 4.63674964411300E-03_wp,-9.10464660774819E-03_wp, &
      & 1.55693316213506E-02_wp,-8.69456073030934E-03_wp,-2.07001602169140E-03_wp, &
      &-1.64641127449032E-03_wp, 7.40412926616573E-04_wp, 5.19987722271890E-02_wp, &
      & 3.50963776224895E-02_wp, 3.53821139878432E-02_wp, 2.03766503438216E-03_wp, &
      & 1.38182723045576E-03_wp,-8.56882698100729E-04_wp, 9.67443679055029E-04_wp, &
      &-7.89383163152661E-04_wp,-3.12016989034458E-04_wp, 4.89673716424246E-02_wp, &
      & 6.85566681173991E-03_wp,-4.33726604714417E-02_wp, 2.32152383555904E-03_wp, &
      &-7.19828442743159E-04_wp, 1.28190858547835E-03_wp, 7.83301665643090E-04_wp, &
      & 1.95937371861310E-04_wp,-2.19901739923401E-04_wp,-2.48692793280382E-02_wp, &
      & 6.56157086460467E-03_wp, 6.03966283303314E-02_wp, 1.66757073509023E-03_wp, &
      & 4.66111281355355E-04_wp, 1.29891259805378E-03_wp,-8.75032383974094E-04_wp, &
      &-1.61471980448914E-03_wp,-2.14613395036928E-03_wp, 6.34879336473785E-03_wp, &
      &-8.73740458812932E-03_wp, 1.82408462816258E-03_wp,-6.61485609356795E-03_wp, &
      & 1.69880913473563E-03_wp, 2.89718143673136E-03_wp,-1.16060082172845E-03_wp, &
      &-3.38055775952816E-03_wp,-3.33058390368733E-02_wp,-1.52201455602624E-01_wp, &
      & 4.38047333432606E-01_wp,-8.69456073030934E-03_wp, 4.36466806645170E-01_wp, &
      & 9.34741416386527E-02_wp,-1.19354104143835E-01_wp,-1.37121127231969E-03_wp, &
      & 5.30583940781882E-04_wp,-2.67897186114664E-02_wp, 4.59737348138373E-03_wp, &
      &-4.02177217080790E-03_wp, 6.15256712973667E-03_wp, 6.28787735416212E-03_wp, &
      & 3.19802362534811E-03_wp, 1.32955044990478E-03_wp, 1.31200798097231E-03_wp, &
      & 1.24744952830266E-02_wp,-7.32813987752802E-02_wp, 3.72566964773942E-03_wp, &
      & 1.02640693286939E-04_wp,-1.47429100000553E-05_wp,-1.17086414461730E-03_wp, &
      & 2.24601344723718E-03_wp,-2.00559116494408E-03_wp,-1.12586025628635E-03_wp, &
      & 2.39815493205939E-02_wp, 2.95058137735713E-02_wp, 2.61728970657019E-02_wp, &
      & 5.32051158886537E-04_wp, 4.14351970336048E-03_wp,-5.00522584023764E-03_wp, &
      &-6.01050825620884E-03_wp, 3.53071480912169E-03_wp, 3.61217559584011E-03_wp, &
      &-2.30645681562275E-03_wp,-2.03180597856797E-03_wp, 3.82256307533336E-03_wp, &
      &-7.67995123149824E-04_wp, 6.38096027411938E-04_wp, 3.16800401176324E-04_wp, &
      & 1.07624737460322E-02_wp, 7.69017399473797E-03_wp,-5.57795006981639E-03_wp, &
      &-3.58517678616803E-02_wp, 9.59547188118657E-02_wp,-2.07001602169140E-03_wp, &
      & 9.34741416386527E-02_wp, 3.90520218771133E-02_wp,-2.06735101101431E-02_wp, &
      &-2.07679312335869E-04_wp,-3.36427527890217E-03_wp,-2.68963265017502E-03_wp, &
      & 7.67642844476133E-02_wp, 1.35665085458970E-03_wp, 1.35542875743456E-03_wp, &
      & 7.30947572371426E-04_wp, 2.17212419192107E-03_wp, 2.27422974116808E-03_wp, &
      &-5.83309848297215E-04_wp, 1.33899406090693E-02_wp,-2.03214366844796E-03_wp, &
      &-5.90050199661660E-02_wp, 1.83884922532074E-03_wp,-5.27133217598405E-04_wp, &
      & 3.10277119084791E-04_wp, 1.78806056495109E-03_wp,-1.42996134361120E-03_wp, &
      & 9.38729389402143E-05_wp, 4.82654826901009E-02_wp,-1.14992339504637E-02_wp, &
      &-6.93910045600238E-02_wp,-1.97145781885208E-03_wp, 3.56893172299223E-05_wp, &
      &-2.34151302886779E-03_wp, 4.72722357628950E-04_wp, 2.70282186006250E-03_wp, &
      & 1.12028617141975E-02_wp,-6.70882068405559E-04_wp, 3.57585309496399E-03_wp, &
      &-1.43702087113090E-03_wp, 2.18809239728138E-03_wp, 9.57106430383583E-03_wp, &
      &-1.69606053549689E-02_wp, 5.60381954566253E-03_wp,-1.76961128934806E-03_wp, &
      & 1.01108773593490E-02_wp, 3.07585591512265E-02_wp,-1.25382771228594E-01_wp, &
      &-1.64641127449032E-03_wp,-1.19354104143835E-01_wp,-2.06735101101431E-02_wp, &
      & 5.64165987847063E-02_wp,-3.60508414656166E-04_wp, 2.36891004263969E-02_wp, &
      & 7.45384160365143E-02_wp, 3.27358666275927E-02_wp, 2.50498630191137E-03_wp, &
      & 3.12765956419526E-04_wp,-1.49064981294504E-03_wp, 7.64346376836081E-04_wp, &
      &-3.81834224001756E-04_wp,-3.89914338100669E-04_wp,-5.16374969918261E-02_wp, &
      &-4.56666343914478E-02_wp, 4.41302193900175E-02_wp,-2.26312525914848E-03_wp, &
      & 2.64354771088751E-03_wp, 1.76499596098071E-04_wp,-2.64863626414187E-03_wp, &
      &-2.83279606761995E-05_wp,-2.41011357456463E-04_wp, 5.49357483690476E-02_wp, &
      &-2.23923554739493E-02_wp,-3.95054115593926E-02_wp,-2.59058745650320E-03_wp, &
      &-2.26084037290193E-03_wp, 4.94952180347453E-04_wp, 2.52712586784912E-03_wp, &
      &-5.02041178188871E-04_wp,-8.64033744991843E-04_wp, 6.96944723684135E-03_wp, &
      & 4.50241444891910E-03_wp, 4.61273866776342E-03_wp, 2.24397055345458E-02_wp, &
      & 2.07478739672486E-02_wp,-7.54471811196526E-03_wp, 1.41316668863599E-02_wp, &
      &-8.85374052511765E-03_wp, 4.01102854511713E-04_wp, 4.52883325072763E-04_wp, &
      &-1.38374238093702E-03_wp, 7.40412926616573E-04_wp,-1.37121127231969E-03_wp, &
      &-2.07679312335869E-04_wp,-3.60508414656166E-04_wp, 9.79385233275982E-01_wp, &
      &-6.73932117684745E-02_wp,-4.24221242518477E-02_wp,-4.54829773946338E-02_wp, &
      &-3.44767227988328E-03_wp,-3.20186908954779E-03_wp, 1.08518088939978E-03_wp, &
      &-2.17906881067328E-03_wp, 1.36803923011760E-03_wp,-1.74429112109329E-03_wp, &
      &-4.60039554544568E-03_wp,-4.62867164240676E-03_wp,-2.29461193380960E-02_wp, &
      & 1.08781582530057E-03_wp, 8.21284940311572E-04_wp, 8.63411726541678E-04_wp, &
      & 3.71234417787815E-04_wp,-2.26885600405926E-03_wp,-4.67052712125803E-04_wp, &
      &-2.18177810703489E-02_wp,-4.85160230141416E-03_wp, 2.44761467320601E-03_wp, &
      & 2.10650456051256E-03_wp, 5.13537734104285E-04_wp, 8.26572953278356E-04_wp, &
      & 7.14219643152562E-04_wp, 1.14929051465422E-03_wp, 3.10371374558179E-02_wp, &
      &-1.08436093353191E-02_wp,-1.94202483147024E-02_wp,-1.97237097978908E-02_wp, &
      & 6.64549541657225E-02_wp, 6.75325136801770E-02_wp, 1.55318190585747E-02_wp, &
      & 7.88193671551557E-02_wp, 3.47268154939371E-02_wp,-4.56200941512203E-02_wp, &
      &-9.32147819443877E-02_wp,-7.41273093138561E-03_wp, 5.19987722271890E-02_wp, &
      & 5.30583940781882E-04_wp,-3.36427527890217E-03_wp, 2.36891004263969E-02_wp, &
      &-6.73932117684745E-02_wp, 8.89614738721827E-01_wp,-1.28950331967082E-02_wp, &
      &-1.39586534564777E-02_wp, 1.53963342226492E-02_wp, 1.29909255552899E-02_wp, &
      &-1.31338339197786E-02_wp,-1.58767803862221E-05_wp,-2.41453153734097E-02_wp, &
      & 5.83630129305121E-03_wp, 4.74494717503912E-03_wp,-1.32063220712849E-02_wp, &
      &-1.19859393613802E-02_wp,-6.98434276122301E-04_wp, 2.29788467999063E-03_wp, &
      & 3.50868101513821E-04_wp,-9.19684784734447E-04_wp,-1.01534697330281E-03_wp, &
      &-1.79249305334720E-02_wp,-4.44215180079783E-02_wp,-4.45142330247158E-03_wp, &
      & 3.22808677054423E-02_wp, 5.51511024939449E-03_wp, 2.23194762525990E-03_wp, &
      & 2.14969208354629E-03_wp, 2.26334563113351E-04_wp,-2.69348432555279E-04_wp, &
      & 2.30333249960052E-02_wp,-1.77208852888238E-02_wp, 1.08132873675751E-02_wp, &
      &-1.20991357049096E-02_wp, 8.85871421287740E-02_wp,-2.15309610649779E-02_wp, &
      &-9.64076933204909E-02_wp,-1.33489350096091E-02_wp,-3.45265330160401E-02_wp, &
      &-6.60382781523261E-03_wp,-5.48349175380994E-03_wp,-1.10468068244529E-02_wp, &
      & 3.50963776224895E-02_wp,-2.67897186114664E-02_wp,-2.68963265017502E-03_wp, &
      & 7.45384160365143E-02_wp,-4.24221242518477E-02_wp,-1.28950331967082E-02_wp, &
      & 9.00965756211415E-01_wp,-9.14415453328002E-03_wp, 1.33841889817848E-03_wp, &
      & 2.15951288119083E-02_wp, 1.24101348943288E-02_wp, 1.49920354704816E-02_wp, &
      &-4.40469099717627E-04_wp, 9.85577803669907E-04_wp,-1.22201335362130E-02_wp, &
      & 2.11341370625129E-02_wp,-1.39791945191686E-02_wp,-7.01721881174691E-04_wp, &
      & 1.08687409943675E-04_wp,-1.44491706265718E-03_wp, 4.54075284812082E-04_wp, &
      &-1.98967421502637E-03_wp, 1.69947916255744E-03_wp,-1.66229817006404E-02_wp, &
      & 2.23035839446687E-02_wp,-9.27699182882334E-03_wp, 4.79918109489262E-04_wp, &
      & 1.40000483985680E-04_wp,-1.14187690595404E-03_wp, 2.76159075758330E-04_wp, &
      & 1.54213587965840E-03_wp, 2.04795146556047E-02_wp,-1.93565751770573E-02_wp, &
      &-1.31156713245789E-02_wp, 5.08336438917607E-03_wp, 6.96854518071820E-04_wp, &
      & 7.89233538672891E-02_wp, 1.10515854057054E-02_wp, 5.13853337547530E-03_wp, &
      &-8.92468896301965E-02_wp, 4.54256036670506E-02_wp,-3.42959429035476E-02_wp, &
      &-1.56187513126319E-03_wp, 3.53821139878432E-02_wp, 4.59737348138373E-03_wp, &
      & 7.67642844476133E-02_wp, 3.27358666275927E-02_wp,-4.54829773946338E-02_wp, &
      &-1.39586534564777E-02_wp,-9.14415453328002E-03_wp, 9.01197413335927E-01_wp, &
      & 2.37424843742069E-02_wp, 2.21729555499775E-04_wp,-8.75069972743066E-03_wp, &
      & 1.29168909226540E-02_wp, 1.70306110002357E-02_wp,-2.30415156877423E-02_wp, &
      & 1.66719242999173E-02_wp, 2.25021243603302E-03_wp,-5.05076972762294E-02_wp, &
      & 5.43204781193947E-03_wp,-9.57154239876852E-05_wp, 2.35591709642209E-03_wp, &
      & 2.63084615326078E-03_wp,-2.80062748714972E-03_wp, 1.34259756367200E-02_wp, &
      & 9.67852731635338E-03_wp,-1.30717775428052E-02_wp,-3.47522554510675E-03_wp, &
      &-2.41917211488738E-03_wp,-1.02717741613644E-03_wp,-7.37519467376673E-04_wp, &
      & 2.37177085330210E-03_wp, 1.59044196104408E-03_wp, 1.17566904494905E-03_wp, &
      &-8.75401591254688E-04_wp,-6.35109734804627E-04_wp,-2.83097626903721E-04_wp, &
      & 9.90846322150787E-04_wp, 2.81701658165517E-03_wp, 1.45947727313188E-04_wp, &
      & 1.55817983246274E-03_wp,-1.30909318966849E-03_wp, 8.13365306994022E-04_wp, &
      &-1.17137609070403E-03_wp,-4.41907081296685E-03_wp, 2.03766503438216E-03_wp, &
      &-4.02177217080790E-03_wp, 1.35665085458970E-03_wp, 2.50498630191137E-03_wp, &
      &-3.44767227988328E-03_wp, 1.53963342226492E-02_wp, 1.33841889817848E-03_wp, &
      & 2.37424843742069E-02_wp, 9.71080523127880E-04_wp, 2.24537932149475E-04_wp, &
      &-5.12935597022903E-04_wp, 3.54826393034930E-04_wp, 1.78175387132417E-05_wp, &
      &-9.92864289497205E-04_wp,-4.93805000442397E-04_wp,-1.05118342011321E-03_wp, &
      &-5.15532684544786E-03_wp, 2.01568567478042E-04_wp, 8.47489533143351E-05_wp, &
      & 1.28915036185723E-04_wp, 5.96716941142845E-05_wp,-1.61436989305193E-04_wp, &
      & 9.05992698760134E-04_wp,-1.86454186027253E-03_wp,-1.56221941951145E-03_wp, &
      &-8.01779151303163E-04_wp, 5.29754230196242E-05_wp,-1.67720095140052E-05_wp, &
      & 6.88550780570842E-05_wp, 1.49459436389064E-04_wp, 3.48745958380269E-05_wp, &
      & 1.12383332749237E-03_wp,-5.65296230452197E-04_wp,-1.41358099309022E-04_wp, &
      &-5.91992885184355E-04_wp, 2.46975374061967E-03_wp, 5.88757207551287E-04_wp, &
      &-2.16103631752938E-03_wp, 5.12251036343349E-04_wp,-3.70938201881427E-04_wp, &
      &-1.45138153891137E-03_wp,-4.30791961546821E-03_wp, 6.66432321276951E-03_wp, &
      & 1.38182723045576E-03_wp, 6.15256712973667E-03_wp, 1.35542875743456E-03_wp, &
      & 3.12765956419526E-04_wp,-3.20186908954779E-03_wp, 1.29909255552899E-02_wp, &
      & 2.15951288119083E-02_wp, 2.21729555499775E-04_wp, 2.24537932149475E-04_wp, &
      & 8.41821697249052E-04_wp, 2.01357871246874E-04_wp, 4.33545604698566E-04_wp, &
      &-3.54326607641045E-04_wp, 6.92419916187062E-04_wp,-2.26880794834199E-03_wp, &
      &-2.60472690055268E-04_wp,-6.40001479013474E-04_wp,-5.97084418745616E-05_wp, &
      & 4.30725140314870E-05_wp,-9.13448849034343E-05_wp, 5.35251386105152E-05_wp, &
      &-1.48760158878772E-04_wp,-3.65628894307319E-04_wp,-3.34642237400698E-03_wp, &
      &-3.20723614998436E-04_wp,-3.02716907708290E-04_wp, 1.66642160459250E-04_wp, &
      & 1.29518345210917E-04_wp,-6.06645445284161E-05_wp,-3.58157141410844E-05_wp, &
      & 1.47884993200396E-04_wp,-3.37161261865593E-04_wp, 1.94865564589505E-04_wp, &
      & 4.86112152718710E-04_wp, 1.38111314514985E-04_wp, 3.12047632637008E-04_wp, &
      &-1.82030938336181E-03_wp,-1.41298732458204E-03_wp,-1.38894336044170E-03_wp, &
      &-1.10750429754854E-04_wp,-4.61944048276185E-04_wp,-8.16503856575622E-04_wp, &
      & 6.93367015797582E-03_wp,-8.56882698100729E-04_wp, 6.28787735416212E-03_wp, &
      & 7.30947572371426E-04_wp,-1.49064981294504E-03_wp, 1.08518088939978E-03_wp, &
      &-1.31338339197786E-02_wp, 1.24101348943288E-02_wp,-8.75069972743066E-03_wp, &
      &-5.12935597022903E-04_wp, 2.01357871246874E-04_wp, 5.57704342368946E-04_wp, &
      & 1.30088685792905E-04_wp, 2.05807133618905E-04_wp, 7.04869867155257E-04_wp, &
      &-4.71965410417172E-04_wp, 1.54322398343485E-03_wp, 1.57912740212245E-03_wp, &
      &-8.60224070526115E-05_wp,-8.26444834398727E-05_wp,-1.22481720956570E-04_wp, &
      & 5.33348705135931E-05_wp,-3.16133131415677E-06_wp, 6.83482385411022E-04_wp, &
      & 1.43421135623700E-03_wp, 1.29954197474567E-03_wp,-6.84034941755676E-04_wp, &
      &-6.77514813142387E-05_wp, 3.57114046682632E-05_wp,-1.37991308366808E-04_wp, &
      &-1.10652382973498E-04_wp, 8.18569440156331E-05_wp, 7.55695280200463E-04_wp, &
      &-6.09542770971772E-04_wp,-6.63799363359947E-05_wp,-9.22681999257671E-05_wp, &
      & 1.29596421478787E-03_wp, 5.97326378130184E-04_wp,-1.48077115538771E-03_wp, &
      & 4.98434304022621E-06_wp,-1.31927484389760E-03_wp, 2.51862372218516E-04_wp, &
      &-1.99648246109626E-03_wp, 3.50242290466801E-03_wp, 9.67443679055029E-04_wp, &
      & 3.19802362534811E-03_wp, 2.17212419192107E-03_wp, 7.64346376836081E-04_wp, &
      &-2.17906881067328E-03_wp,-1.58767803862221E-05_wp, 1.49920354704816E-02_wp, &
      & 1.29168909226540E-02_wp, 3.54826393034930E-04_wp, 4.33545604698566E-04_wp, &
      & 1.30088685792905E-04_wp, 4.82820109678503E-04_wp, 2.51108863692765E-04_wp, &
      &-7.79172865409633E-04_wp,-7.66705658334690E-04_wp, 2.26017817007804E-04_wp, &
      &-3.30674201547255E-03_wp, 1.03484025349429E-04_wp,-1.06917681623253E-06_wp, &
      &-1.97825617795629E-07_wp, 1.08987968326367E-04_wp,-1.54513083585891E-04_wp, &
      & 9.18499763895862E-04_wp, 2.49064579621402E-04_wp,-3.01515131463446E-04_wp, &
      &-1.87144089306515E-03_wp,-4.15809966344942E-05_wp, 1.50166180257202E-05_wp, &
      &-9.17945277868815E-05_wp, 3.21446802816235E-05_wp, 1.29877997003160E-04_wp, &
      &-4.91653601592628E-04_wp,-6.86637147631842E-05_wp, 2.63148272706973E-04_wp, &
      & 7.37417254726881E-04_wp,-1.29311438090641E-03_wp,-5.00617661066083E-04_wp, &
      &-3.68872484730471E-05_wp,-1.55155659357652E-03_wp,-1.78842876359272E-03_wp, &
      & 2.07065310902566E-03_wp, 1.66375144799501E-03_wp, 1.53391697941402E-03_wp, &
      &-7.89383163152661E-04_wp, 1.32955044990478E-03_wp, 2.27422974116808E-03_wp, &
      &-3.81834224001756E-04_wp, 1.36803923011760E-03_wp,-2.41453153734097E-02_wp, &
      &-4.40469099717627E-04_wp, 1.70306110002357E-02_wp, 1.78175387132417E-05_wp, &
      &-3.54326607641045E-04_wp, 2.05807133618905E-04_wp, 2.51108863692765E-04_wp, &
      & 1.00296300267918E-03_wp,-2.26339054665830E-03_wp, 1.22673110868470E-03_wp, &
      & 1.42194057769701E-03_wp,-2.96462717477253E-03_wp, 1.86342724132078E-04_wp, &
      &-1.01495495974086E-04_wp, 4.60297941197738E-05_wp, 1.31730529441688E-04_wp, &
      &-5.02008128433043E-05_wp, 2.19806549831163E-03_wp, 5.36800952916396E-03_wp, &
      &-1.11865135939340E-04_wp,-2.11155401443376E-03_wp,-3.09423947555397E-04_wp, &
      &-1.21929197699361E-04_wp,-1.22852682598191E-04_wp, 4.05393279050851E-05_wp, &
      & 4.47136744942552E-05_wp,-1.02695166967481E-03_wp,-7.77156635669613E-03_wp, &
      &-3.69338959899031E-03_wp, 5.31578919345793E-03_wp,-2.40673203813859E-02_wp, &
      & 1.84804474365555E-02_wp,-9.73129495503449E-03_wp,-1.19435034444835E-02_wp, &
      &-8.63307150232750E-03_wp, 3.41550940627179E-04_wp, 7.55268864056282E-05_wp, &
      & 5.43354709716936E-04_wp,-3.12016989034458E-04_wp, 1.31200798097231E-03_wp, &
      &-5.83309848297215E-04_wp,-3.89914338100669E-04_wp,-1.74429112109329E-03_wp, &
      & 5.83630129305121E-03_wp, 9.85577803669907E-04_wp,-2.30415156877423E-02_wp, &
      &-9.92864289497205E-04_wp, 6.92419916187062E-04_wp, 7.04869867155257E-04_wp, &
      &-7.79172865409633E-04_wp,-2.26339054665830E-03_wp, 9.80319653876711E-01_wp, &
      & 6.96148081109755E-02_wp, 3.39853186201215E-02_wp,-4.74581076068379E-02_wp, &
      & 3.77653932652044E-03_wp,-2.82491840284670E-03_wp, 1.53601953767155E-03_wp, &
      & 1.88938374714825E-03_wp, 1.39873724097239E-03_wp,-3.03607735046782E-03_wp, &
      & 1.59310863120059E-02_wp, 1.87885921811414E-02_wp, 5.08467535689386E-03_wp, &
      &-1.69023816854932E-03_wp,-6.52194452941343E-04_wp,-1.24575925065902E-03_wp, &
      &-1.84891013912983E-03_wp,-3.95093353420261E-04_wp,-3.45289149212513E-02_wp, &
      &-1.12023415897402E-02_wp,-1.62131073853426E-02_wp, 2.10073897134136E-02_wp, &
      & 7.93390435947526E-02_wp,-6.44185602640536E-02_wp, 3.62507400525617E-04_wp, &
      & 7.87707633106604E-02_wp,-3.79182506231862E-02_wp,-5.80400261978655E-02_wp, &
      & 7.36998183172769E-02_wp, 1.40237003812439E-02_wp, 4.89673716424246E-02_wp, &
      & 1.24744952830266E-02_wp, 1.33899406090693E-02_wp,-5.16374969918261E-02_wp, &
      &-4.60039554544568E-03_wp, 4.74494717503912E-03_wp,-1.22201335362130E-02_wp, &
      & 1.66719242999173E-02_wp,-4.93805000442397E-04_wp,-2.26880794834199E-03_wp, &
      &-4.71965410417172E-04_wp,-7.66705658334690E-04_wp, 1.22673110868470E-03_wp, &
      & 6.96148081109755E-02_wp, 8.84483522477175E-01_wp,-6.27087331795843E-03_wp, &
      & 1.39890422776412E-02_wp, 1.70641049794269E-02_wp,-1.41030186631462E-02_wp, &
      & 1.27499204873158E-02_wp, 3.86236907846822E-04_wp, 2.29766788919789E-02_wp, &
      & 9.93178929003970E-03_wp,-2.68014119418027E-03_wp,-3.97830522408287E-02_wp, &
      &-1.26975405204845E-02_wp, 1.43642412751507E-03_wp,-5.15623780963257E-04_wp, &
      & 1.97029419194627E-03_wp, 3.39979714316641E-03_wp, 1.57372281682921E-04_wp, &
      &-1.34121866164477E-02_wp,-1.78971328952172E-02_wp, 1.12706765571304E-02_wp, &
      & 1.09429414011269E-02_wp, 7.02903965683136E-02_wp, 3.79320344196822E-02_wp, &
      & 8.64415228683305E-02_wp,-3.00272463633091E-02_wp, 3.52704972597678E-02_wp, &
      &-1.74917522456011E-02_wp, 8.59929205900434E-03_wp, 3.28804253492689E-02_wp, &
      & 6.85566681173991E-03_wp,-7.32813987752802E-02_wp,-2.03214366844796E-03_wp, &
      &-4.56666343914478E-02_wp,-4.62867164240676E-03_wp,-1.32063220712849E-02_wp, &
      & 2.11341370625129E-02_wp, 2.25021243603302E-03_wp,-1.05118342011321E-03_wp, &
      &-2.60472690055268E-04_wp, 1.54322398343485E-03_wp, 2.26017817007804E-04_wp, &
      & 1.42194057769701E-03_wp, 3.39853186201215E-02_wp,-6.27087331795843E-03_wp, &
      & 9.01725438066412E-01_wp, 8.33818371867964E-03_wp,-3.41304494398283E-04_wp, &
      &-2.29885786612005E-02_wp,-1.59420543751472E-02_wp, 1.67323852681583E-02_wp, &
      & 5.90543543298967E-04_wp, 1.89070265565719E-02_wp,-2.18204966861747E-02_wp, &
      &-1.64116449632240E-02_wp,-3.50049447183280E-02_wp, 1.92048599024866E-03_wp, &
      & 2.73294147835100E-04_wp, 2.40805559729263E-04_wp, 3.98704788703111E-03_wp, &
      & 3.83007723350443E-03_wp, 2.17842958083843E-02_wp, 2.23980508202744E-02_wp, &
      & 1.05297118737512E-02_wp, 3.20167044114457E-03_wp,-1.26672314141543E-02_wp, &
      & 6.87221148902045E-02_wp,-5.52582904665705E-03_wp,-2.09804044214836E-03_wp, &
      &-9.70603235678339E-02_wp,-5.30554321998536E-02_wp,-3.29569302511059E-02_wp, &
      & 1.21189122316512E-02_wp,-4.33726604714417E-02_wp, 3.72566964773942E-03_wp, &
      &-5.90050199661660E-02_wp, 4.41302193900175E-02_wp,-2.29461193380960E-02_wp, &
      &-1.19859393613802E-02_wp,-1.39791945191686E-02_wp,-5.05076972762294E-02_wp, &
      &-5.15532684544786E-03_wp,-6.40001479013474E-04_wp, 1.57912740212245E-03_wp, &
      &-3.30674201547255E-03_wp,-2.96462717477253E-03_wp,-4.74581076068379E-02_wp, &
      & 1.39890422776412E-02_wp, 8.33818371867964E-03_wp, 9.00709592745119E-01_wp, &
      &-2.32277295829346E-02_wp,-5.91310561328952E-04_wp,-1.03765797970833E-02_wp, &
      &-1.37570262673529E-02_wp, 1.79170712797340E-02_wp, 1.33699745251540E-02_wp, &
      & 1.12806043159415E-02_wp,-9.44493581999975E-03_wp,-1.90246054652711E-03_wp, &
      &-1.78419681380527E-03_wp,-9.11286541924675E-04_wp,-5.17440346512460E-04_wp, &
      & 2.69307620534820E-03_wp, 1.58235914557377E-03_wp,-1.32239449246803E-03_wp, &
      &-8.58573131874772E-04_wp,-6.83652635175080E-04_wp, 3.34046150737528E-04_wp, &
      & 1.72311948156705E-03_wp,-2.66981711043062E-03_wp, 2.36297575249206E-04_wp, &
      & 1.14598547458529E-03_wp, 1.09039808896316E-03_wp, 5.62332461341726E-04_wp, &
      & 2.35295984835294E-03_wp,-6.66975469219293E-05_wp, 2.32152383555904E-03_wp, &
      & 1.02640693286939E-04_wp, 1.83884922532074E-03_wp,-2.26312525914848E-03_wp, &
      & 1.08781582530057E-03_wp,-6.98434276122301E-04_wp,-7.01721881174691E-04_wp, &
      & 5.43204781193947E-03_wp, 2.01568567478042E-04_wp,-5.97084418745616E-05_wp, &
      &-8.60224070526115E-05_wp, 1.03484025349429E-04_wp, 1.86342724132078E-04_wp, &
      & 3.77653932652044E-03_wp, 1.70641049794269E-02_wp,-3.41304494398283E-04_wp, &
      &-2.32277295829346E-02_wp, 9.71954797249471E-04_wp,-2.66463312336126E-04_wp, &
      & 5.36546546854042E-04_wp, 3.78261896705893E-04_wp,-2.20377842351566E-05_wp, &
      &-1.08352826422062E-03_wp,-1.82074210886211E-03_wp,-1.01687330465836E-03_wp, &
      & 1.84168810429607E-03_wp, 1.24242288971349E-04_wp, 2.60909585417293E-05_wp, &
      & 9.82423307905643E-05_wp,-1.33251711431608E-05_wp,-7.80121783928019E-05_wp, &
      & 9.23763188397653E-04_wp, 6.50886246444007E-04_wp,-3.67543219007418E-06_wp, &
      &-7.19136030754418E-04_wp,-2.33142012668874E-03_wp, 3.40407076461401E-04_wp, &
      &-1.94458556002641E-03_wp,-1.69273880680804E-04_wp,-3.59947499872617E-04_wp, &
      & 1.78951892334922E-03_wp,-7.98432581416592E-04_wp,-3.13907873268688E-03_wp, &
      &-7.19828442743159E-04_wp,-1.47429100000553E-05_wp,-5.27133217598405E-04_wp, &
      & 2.64354771088751E-03_wp, 8.21284940311572E-04_wp, 2.29788467999063E-03_wp, &
      & 1.08687409943675E-04_wp,-9.57154239876852E-05_wp, 8.47489533143351E-05_wp, &
      & 4.30725140314870E-05_wp,-8.26444834398727E-05_wp,-1.06917681623253E-06_wp, &
      &-1.01495495974086E-04_wp,-2.82491840284670E-03_wp,-1.41030186631462E-02_wp, &
      &-2.29885786612005E-02_wp,-5.91310561328952E-04_wp,-2.66463312336126E-04_wp, &
      & 8.50416540478395E-04_wp, 2.20268366240725E-04_wp,-4.57895934640666E-04_wp, &
      &-3.97935121573021E-04_wp,-1.96526810803974E-03_wp, 3.13903002693450E-03_wp, &
      & 3.85187428686917E-03_wp, 2.11734626341255E-03_wp,-1.34264453635969E-04_wp, &
      &-2.51838077167750E-05_wp,-5.85408576229357E-05_wp,-2.31098244614495E-04_wp, &
      &-1.74472472563477E-04_wp,-5.87060101040320E-04_wp,-1.23589021153629E-04_wp, &
      &-5.64274993106786E-04_wp, 4.42987781466214E-05_wp, 1.41048483944506E-04_wp, &
      &-2.12977682744976E-03_wp,-1.28072235876254E-03_wp, 1.52074117567880E-03_wp, &
      &-3.60824594442789E-04_wp, 5.49306638141359E-04_wp, 2.29759184132251E-03_wp, &
      &-3.25713971009019E-03_wp, 1.28190858547835E-03_wp,-1.17086414461730E-03_wp, &
      & 3.10277119084791E-04_wp, 1.76499596098071E-04_wp, 8.63411726541678E-04_wp, &
      & 3.50868101513821E-04_wp,-1.44491706265718E-03_wp, 2.35591709642209E-03_wp, &
      & 1.28915036185723E-04_wp,-9.13448849034343E-05_wp,-1.22481720956570E-04_wp, &
      &-1.97825617795629E-07_wp, 4.60297941197738E-05_wp, 1.53601953767155E-03_wp, &
      & 1.27499204873158E-02_wp,-1.59420543751472E-02_wp,-1.03765797970833E-02_wp, &
      & 5.36546546854042E-04_wp, 2.20268366240725E-04_wp, 6.09265712051661E-04_wp, &
      &-1.47519943592145E-04_wp, 1.22238728619821E-04_wp,-1.25095952083766E-03_wp, &
      & 2.27846314193697E-04_wp, 1.88274865293780E-04_wp, 2.05809775537036E-03_wp, &
      & 1.14606708087605E-05_wp,-2.53947944812789E-05_wp, 7.94171880809830E-05_wp, &
      &-5.67477088128098E-05_wp,-1.57956030127115E-04_wp,-6.06266471532816E-04_wp, &
      &-7.14602097204805E-04_wp,-4.99490489386992E-05_wp, 2.30642331218448E-04_wp, &
      & 1.23754087718086E-03_wp,-3.25448102001268E-04_wp, 1.67826970775080E-03_wp, &
      &-6.69143343856796E-04_wp, 1.85106237168900E-03_wp, 1.18354420618126E-04_wp, &
      &-6.70858028794763E-04_wp, 4.31017066296889E-03_wp, 7.83301665643090E-04_wp, &
      & 2.24601344723718E-03_wp, 1.78806056495109E-03_wp,-2.64863626414187E-03_wp, &
      & 3.71234417787815E-04_wp,-9.19684784734447E-04_wp, 4.54075284812082E-04_wp, &
      & 2.63084615326078E-03_wp, 5.96716941142845E-05_wp, 5.35251386105152E-05_wp, &
      & 5.33348705135931E-05_wp, 1.08987968326367E-04_wp, 1.31730529441688E-04_wp, &
      & 1.88938374714825E-03_wp, 3.86236907846822E-04_wp, 1.67323852681583E-02_wp, &
      &-1.37570262673529E-02_wp, 3.78261896705893E-04_wp,-4.57895934640666E-04_wp, &
      &-1.47519943592145E-04_wp, 5.67070951883082E-04_wp,-2.66365841731747E-04_wp, &
      &-2.17685295273928E-05_wp,-2.21004109845667E-03_wp,-1.02849359837745E-03_wp, &
      &-9.02442359022939E-04_wp, 1.12455065949128E-04_wp, 6.84526285745388E-05_wp, &
      &-1.17897739761924E-05_wp, 2.91415636451971E-05_wp, 1.16342256034976E-04_wp, &
      &-5.25957119780333E-04_wp, 1.61848868701550E-04_wp,-1.70869611258782E-04_wp, &
      & 6.80365098404692E-04_wp, 1.32235237491202E-03_wp,-6.09023766965693E-04_wp, &
      &-2.64680965053990E-04_wp, 1.75734486260590E-03_wp,-2.29662810010202E-03_wp, &
      &-2.60744135287044E-03_wp, 2.19545673464568E-03_wp,-1.47887531130849E-03_wp, &
      & 1.95937371861310E-04_wp,-2.00559116494408E-03_wp,-1.42996134361120E-03_wp, &
      &-2.83279606761995E-05_wp,-2.26885600405926E-03_wp,-1.01534697330281E-03_wp, &
      &-1.98967421502637E-03_wp,-2.80062748714972E-03_wp,-1.61436989305193E-04_wp, &
      &-1.48760158878772E-04_wp,-3.16133131415677E-06_wp,-1.54513083585891E-04_wp, &
      &-5.02008128433043E-05_wp, 1.39873724097239E-03_wp, 2.29766788919789E-02_wp, &
      & 5.90543543298967E-04_wp, 1.79170712797340E-02_wp,-2.20377842351566E-05_wp, &
      &-3.97935121573021E-04_wp, 1.22238728619821E-04_wp,-2.66365841731747E-04_wp, &
      & 9.76918208350585E-04_wp, 1.37310447814700E-03_wp,-3.85050874153481E-04_wp, &
      &-4.29155366732576E-03_wp,-8.98951630115442E-04_wp, 7.11466472205285E-06_wp, &
      &-7.22055400547691E-05_wp, 1.11585039122059E-04_wp, 2.53643008088297E-04_wp, &
      & 2.63123070309016E-05_wp,-9.44887800903170E-04_wp, 2.47681579669712E-03_wp, &
      &-3.50216085250586E-03_wp,-8.84877256219468E-03_wp,-1.35903661755520E-02_wp, &
      &-5.22282280974187E-03_wp,-9.91183630949829E-03_wp, 2.14421083889722E-02_wp, &
      & 2.20191489255154E-02_wp,-1.84935397407614E-04_wp, 9.16685372743052E-04_wp, &
      &-1.84087838689988E-04_wp,-2.19901739923401E-04_wp,-1.12586025628635E-03_wp, &
      & 9.38729389402143E-05_wp,-2.41011357456463E-04_wp,-4.67052712125803E-04_wp, &
      &-1.79249305334720E-02_wp, 1.69947916255744E-03_wp, 1.34259756367200E-02_wp, &
      & 9.05992698760134E-04_wp,-3.65628894307319E-04_wp, 6.83482385411022E-04_wp, &
      & 9.18499763895862E-04_wp, 2.19806549831163E-03_wp,-3.03607735046782E-03_wp, &
      & 9.93178929003970E-03_wp, 1.89070265565719E-02_wp, 1.33699745251540E-02_wp, &
      &-1.08352826422062E-03_wp,-1.96526810803974E-03_wp,-1.25095952083766E-03_wp, &
      &-2.17685295273928E-05_wp, 1.37310447814700E-03_wp, 9.79754791069694E-01_wp, &
      &-2.23989504876204E-02_wp, 3.40636152248749E-02_wp, 8.18979621687276E-02_wp, &
      & 2.09012487006927E-03_wp, 8.40503347739348E-04_wp, 1.56038857676965E-03_wp, &
      &-3.26760910350976E-03_wp,-3.45211447044762E-03_wp, 9.71328579477840E-03_wp, &
      & 1.51505155319784E-02_wp, 4.54245265617874E-03_wp, 1.35788864585994E-02_wp, &
      & 7.65937840551549E-02_wp, 3.75844989308166E-02_wp,-5.57175718616085E-03_wp, &
      & 3.28723251997694E-02_wp, 7.33216018913181E-02_wp, 4.00823788730297E-02_wp, &
      & 2.57126867846738E-02_wp, 1.41827982066621E-02_wp,-2.48692793280382E-02_wp, &
      & 2.39815493205939E-02_wp, 4.82654826901009E-02_wp, 5.49357483690476E-02_wp, &
      &-2.18177810703489E-02_wp,-4.44215180079783E-02_wp,-1.66229817006404E-02_wp, &
      & 9.67852731635338E-03_wp,-1.86454186027253E-03_wp,-3.34642237400698E-03_wp, &
      & 1.43421135623700E-03_wp, 2.49064579621402E-04_wp, 5.36800952916396E-03_wp, &
      & 1.59310863120059E-02_wp,-2.68014119418027E-03_wp,-2.18204966861747E-02_wp, &
      & 1.12806043159415E-02_wp,-1.82074210886211E-03_wp, 3.13903002693450E-03_wp, &
      & 2.27846314193697E-04_wp,-2.21004109845667E-03_wp,-3.85050874153481E-04_wp, &
      &-2.23989504876204E-02_wp, 9.07216790441336E-01_wp, 4.71099439338665E-03_wp, &
      & 6.98669050592126E-03_wp,-2.76452181520636E-02_wp,-1.32617256121538E-02_wp, &
      &-5.84858535890091E-03_wp,-9.24621410530014E-04_wp,-8.64305835005504E-03_wp, &
      &-1.32264514234976E-02_wp, 4.13303227392449E-03_wp, 1.19055309005284E-02_wp, &
      &-2.03601998000453E-02_wp, 3.28686383454786E-02_wp,-1.69469206157941E-02_wp, &
      & 8.67153633545441E-02_wp, 4.78972976105153E-02_wp,-7.19861577004959E-02_wp, &
      & 5.07256361601547E-02_wp, 1.14966740670106E-02_wp,-6.62740370087469E-02_wp, &
      & 6.56157086460467E-03_wp, 2.95058137735713E-02_wp,-1.14992339504637E-02_wp, &
      &-2.23923554739493E-02_wp,-4.85160230141416E-03_wp,-4.45142330247158E-03_wp, &
      & 2.23035839446687E-02_wp,-1.30717775428052E-02_wp,-1.56221941951145E-03_wp, &
      &-3.20723614998436E-04_wp, 1.29954197474567E-03_wp,-3.01515131463446E-04_wp, &
      &-1.11865135939340E-04_wp, 1.87885921811414E-02_wp,-3.97830522408287E-02_wp, &
      &-1.64116449632240E-02_wp,-9.44493581999975E-03_wp,-1.01687330465836E-03_wp, &
      & 3.85187428686917E-03_wp, 1.88274865293780E-04_wp,-1.02849359837745E-03_wp, &
      &-4.29155366732576E-03_wp, 3.40636152248749E-02_wp, 4.71099439338665E-03_wp, &
      & 9.03653999747069E-01_wp,-1.07477493791768E-02_wp,-6.79035487611753E-04_wp, &
      & 7.65027028291042E-03_wp,-1.49468871725635E-02_wp,-2.63446016222559E-02_wp, &
      &-1.04417300492389E-03_wp,-4.00506954182355E-02_wp, 1.19937863232238E-02_wp, &
      &-1.89552860619114E-02_wp,-2.27788852326930E-02_wp, 6.44067474820537E-02_wp, &
      & 4.36814617080786E-02_wp, 2.90410287476276E-03_wp,-1.04712017981954E-01_wp, &
      &-4.69373915752494E-02_wp, 6.01728319805672E-02_wp, 3.88665260783668E-02_wp, &
      & 2.39879765902584E-02_wp, 6.03966283303314E-02_wp, 2.61728970657019E-02_wp, &
      &-6.93910045600238E-02_wp,-3.95054115593926E-02_wp, 2.44761467320601E-03_wp, &
      & 3.22808677054423E-02_wp,-9.27699182882334E-03_wp,-3.47522554510675E-03_wp, &
      &-8.01779151303163E-04_wp,-3.02716907708290E-04_wp,-6.84034941755676E-04_wp, &
      &-1.87144089306515E-03_wp,-2.11155401443376E-03_wp, 5.08467535689386E-03_wp, &
      &-1.26975405204845E-02_wp,-3.50049447183280E-02_wp,-1.90246054652711E-03_wp, &
      & 1.84168810429607E-03_wp, 2.11734626341255E-03_wp, 2.05809775537036E-03_wp, &
      &-9.02442359022939E-04_wp,-8.98951630115442E-04_wp, 8.18979621687276E-02_wp, &
      & 6.98669050592126E-03_wp,-1.07477493791768E-02_wp, 8.78220908101747E-01_wp, &
      & 8.78198962890551E-03_wp, 7.29973779520088E-04_wp, 1.48135170964552E-02_wp, &
      &-1.41265088250659E-02_wp,-2.63165868327368E-02_wp,-7.30975427202840E-04_wp, &
      &-4.00974483671558E-04_wp,-4.21614601949817E-04_wp,-6.78210598851248E-04_wp, &
      &-1.00233070466989E-03_wp,-7.35276572443060E-04_wp, 3.78053496299328E-04_wp, &
      &-1.59454356614859E-03_wp,-2.10917504939881E-03_wp,-1.08074030232713E-03_wp, &
      &-9.20983705270224E-04_wp, 1.05566543856893E-03_wp, 1.66757073509023E-03_wp, &
      & 5.32051158886537E-04_wp,-1.97145781885208E-03_wp,-2.59058745650320E-03_wp, &
      & 2.10650456051256E-03_wp, 5.51511024939449E-03_wp, 4.79918109489262E-04_wp, &
      &-2.41917211488738E-03_wp, 5.29754230196242E-05_wp, 1.66642160459250E-04_wp, &
      &-6.77514813142387E-05_wp,-4.15809966344942E-05_wp,-3.09423947555397E-04_wp, &
      &-1.69023816854932E-03_wp, 1.43642412751507E-03_wp, 1.92048599024866E-03_wp, &
      &-1.78419681380527E-03_wp, 1.24242288971349E-04_wp,-1.34264453635969E-04_wp, &
      & 1.14606708087605E-05_wp, 1.12455065949128E-04_wp, 7.11466472205285E-06_wp, &
      & 2.09012487006927E-03_wp,-2.76452181520636E-02_wp,-6.79035487611753E-04_wp, &
      & 8.78198962890551E-03_wp, 9.72829001384354E-04_wp, 4.31613119157028E-04_wp, &
      & 3.39329532587659E-04_wp,-1.10362577227386E-04_wp, 8.20860388579246E-06_wp, &
      &-2.86085329176629E-04_wp,-1.44661987626867E-04_wp,-8.16969303097783E-05_wp, &
      &-4.08503618156882E-04_wp,-6.88819688386836E-04_wp,-6.06690581189437E-04_wp, &
      & 9.38622543054151E-04_wp,-7.29860323367932E-05_wp,-1.45840439450965E-03_wp, &
      &-5.41519705134412E-04_wp,-1.88651521181885E-03_wp, 3.60304872454716E-03_wp, &
      & 4.66111281355355E-04_wp, 4.14351970336048E-03_wp, 3.56893172299223E-05_wp, &
      &-2.26084037290193E-03_wp, 5.13537734104285E-04_wp, 2.23194762525990E-03_wp, &
      & 1.40000483985680E-04_wp,-1.02717741613644E-03_wp,-1.67720095140052E-05_wp, &
      & 1.29518345210917E-04_wp, 3.57114046682632E-05_wp, 1.50166180257202E-05_wp, &
      &-1.21929197699361E-04_wp,-6.52194452941343E-04_wp,-5.15623780963257E-04_wp, &
      & 2.73294147835100E-04_wp,-9.11286541924675E-04_wp, 2.60909585417293E-05_wp, &
      &-2.51838077167750E-05_wp,-2.53947944812789E-05_wp, 6.84526285745388E-05_wp, &
      &-7.22055400547691E-05_wp, 8.40503347739348E-04_wp,-1.32617256121538E-02_wp, &
      & 7.65027028291042E-03_wp, 7.29973779520088E-04_wp, 4.31613119157028E-04_wp, &
      & 3.07330023331222E-04_wp,-6.96308160685538E-05_wp,-2.72899832827896E-04_wp, &
      & 1.35850391213250E-04_wp,-5.83829373771101E-04_wp,-2.20905355403101E-05_wp, &
      &-5.05453807145344E-04_wp,-1.59367877771823E-04_wp, 3.53857312558554E-04_wp, &
      & 6.06034207884786E-04_wp,-1.22901953593365E-03_wp,-2.39582887843716E-03_wp, &
      & 1.08123286917808E-04_wp, 1.74128670059472E-04_wp, 2.10285595596077E-03_wp, &
      &-3.54139161540637E-03_wp, 1.29891259805378E-03_wp,-5.00522584023764E-03_wp, &
      &-2.34151302886779E-03_wp, 4.94952180347453E-04_wp, 8.26572953278356E-04_wp, &
      & 2.14969208354629E-03_wp,-1.14187690595404E-03_wp,-7.37519467376673E-04_wp, &
      & 6.88550780570842E-05_wp,-6.06645445284161E-05_wp,-1.37991308366808E-04_wp, &
      &-9.17945277868815E-05_wp,-1.22852682598191E-04_wp,-1.24575925065902E-03_wp, &
      & 1.97029419194627E-03_wp, 2.40805559729263E-04_wp,-5.17440346512460E-04_wp, &
      & 9.82423307905643E-05_wp,-5.85408576229357E-05_wp, 7.94171880809830E-05_wp, &
      &-1.17897739761924E-05_wp, 1.11585039122059E-04_wp, 1.56038857676965E-03_wp, &
      &-5.84858535890091E-03_wp,-1.49468871725635E-02_wp, 1.48135170964552E-02_wp, &
      & 3.39329532587659E-04_wp,-6.96308160685538E-05_wp, 5.91760718214581E-04_wp, &
      & 2.52912255874505E-04_wp,-4.23264327989375E-04_wp, 1.06473924707788E-03_wp, &
      &-4.45102163842215E-04_wp, 3.80106625787291E-05_wp, 1.00614012061151E-03_wp, &
      &-1.40694498188145E-03_wp, 1.10660155876927E-04_wp,-2.36862461809501E-03_wp, &
      & 4.78040749240586E-04_wp, 2.31640389941183E-03_wp,-2.34898043486087E-03_wp, &
      & 5.97064769894395E-04_wp,-2.95929713093146E-03_wp,-8.75032383974094E-04_wp, &
      &-6.01050825620884E-03_wp, 4.72722357628950E-04_wp, 2.52712586784912E-03_wp, &
      & 7.14219643152562E-04_wp, 2.26334563113351E-04_wp, 2.76159075758330E-04_wp, &
      & 2.37177085330210E-03_wp, 1.49459436389064E-04_wp,-3.58157141410844E-05_wp, &
      &-1.10652382973498E-04_wp, 3.21446802816235E-05_wp, 4.05393279050851E-05_wp, &
      &-1.84891013912983E-03_wp, 3.39979714316641E-03_wp, 3.98704788703111E-03_wp, &
      & 2.69307620534820E-03_wp,-1.33251711431608E-05_wp,-2.31098244614495E-04_wp, &
      &-5.67477088128098E-05_wp, 2.91415636451971E-05_wp, 2.53643008088297E-04_wp, &
      &-3.26760910350976E-03_wp,-9.24621410530014E-04_wp,-2.63446016222559E-02_wp, &
      &-1.41265088250659E-02_wp,-1.10362577227386E-04_wp,-2.72899832827896E-04_wp, &
      & 2.52912255874505E-04_wp, 1.08310028976737E-03_wp, 4.38482863597593E-04_wp, &
      & 1.20864032141500E-03_wp,-5.21718641127864E-04_wp, 4.70419178269535E-04_wp, &
      & 6.42501172500844E-04_wp,-2.36540545013052E-03_wp,-1.13219077295543E-03_wp, &
      & 3.42392873712068E-05_wp, 2.59547000700204E-03_wp, 5.16614467009417E-04_wp, &
      &-2.72927683099535E-03_wp,-3.32760114542134E-03_wp, 4.28984898113645E-03_wp, &
      &-1.61471980448914E-03_wp, 3.53071480912169E-03_wp, 2.70282186006250E-03_wp, &
      &-5.02041178188871E-04_wp, 1.14929051465422E-03_wp,-2.69348432555279E-04_wp, &
      & 1.54213587965840E-03_wp, 1.59044196104408E-03_wp, 3.48745958380269E-05_wp, &
      & 1.47884993200396E-04_wp, 8.18569440156331E-05_wp, 1.29877997003160E-04_wp, &
      & 4.47136744942552E-05_wp,-3.95093353420261E-04_wp, 1.57372281682921E-04_wp, &
      & 3.83007723350443E-03_wp, 1.58235914557377E-03_wp,-7.80121783928019E-05_wp, &
      &-1.74472472563477E-04_wp,-1.57956030127115E-04_wp, 1.16342256034976E-04_wp, &
      & 2.63123070309016E-05_wp,-3.45211447044762E-03_wp,-8.64305835005504E-03_wp, &
      &-1.04417300492389E-03_wp,-2.63165868327368E-02_wp, 8.20860388579246E-06_wp, &
      & 1.35850391213250E-04_wp,-4.23264327989375E-04_wp, 4.38482863597593E-04_wp, &
      & 9.40813877254780E-04_wp, 7.11952239175513E-03_wp,-4.48759433563561E-04_wp, &
      & 7.04758808977056E-04_wp,-2.94624945675376E-04_wp,-2.74166581356306E-03_wp, &
      & 6.84153465016098E-03_wp,-5.31279761822322E-03_wp, 4.85658416213561E-03_wp, &
      & 9.50019121810167E-04_wp,-2.43299151738308E-03_wp,-9.06679526463874E-03_wp, &
      &-3.80607531989858E-04_wp,-2.37874651841337E-03_wp,-3.06390494210790E-04_wp, &
      & 3.36294433121553E-03_wp, 8.99503672589712E-03_wp,-1.11688855965218E-03_wp, &
      & 3.20046292353588E-02_wp, 2.32606695048053E-02_wp, 2.10966188523540E-02_wp, &
      & 1.14006433275913E-03_wp, 1.14147229303145E-03_wp,-3.02272998088046E-04_wp, &
      & 7.61032562708332E-04_wp,-4.79144450852178E-04_wp,-1.30082860861281E-03_wp, &
      &-3.51696354973964E-02_wp,-1.39407590739093E-02_wp, 2.24022500427106E-02_wp, &
      &-1.30001097214208E-03_wp, 9.03676498870593E-04_wp,-5.83171540616389E-04_wp, &
      &-5.86376985234641E-04_wp,-5.20172683811391E-04_wp,-1.20825923162660E-03_wp, &
      & 1.01661826271254E-02_wp,-1.38165545028928E-02_wp,-4.08742499801280E-02_wp, &
      &-7.18067691148896E-04_wp,-2.63641797653755E-04_wp,-5.93966433083412E-04_wp, &
      & 1.02979125404226E-03_wp, 1.20885516832012E-03_wp,-4.48759433563561E-04_wp, &
      & 2.28687129134323E-03_wp, 5.41396449946994E-04_wp,-2.59430166947037E-04_wp, &
      &-2.64919936937235E-03_wp, 5.77530559995348E-04_wp,-5.26976013457427E-05_wp, &
      &-1.70360657293966E-03_wp,-2.12815610709220E-04_wp, 1.06766405478892E-03_wp, &
      & 1.23297231618413E-03_wp,-1.24027234557119E-03_wp,-2.62975652749730E-03_wp, &
      & 6.81819221360930E-04_wp,-2.77984913399897E-03_wp, 6.15223438265065E-04_wp, &
      & 7.17489050802549E-03_wp,-1.15284447884868E-02_wp,-1.66995151518806E-02_wp, &
      &-1.96285390992293E-02_wp,-8.09390428732359E-04_wp,-6.27932511020678E-04_wp, &
      & 1.21450766756529E-04_wp,-6.28894926244606E-04_wp,-7.62484937973928E-05_wp, &
      &-8.00507020680877E-03_wp,-1.18379012605914E-02_wp,-1.74855485990935E-02_wp, &
      & 2.22060923985087E-02_wp,-8.41884382085824E-04_wp, 6.54018323036404E-04_wp, &
      &-9.79735922630597E-05_wp,-7.38014920537391E-04_wp, 1.67408952638943E-04_wp, &
      & 2.55413860729764E-03_wp, 1.44077180987680E-02_wp, 4.81975145073898E-03_wp, &
      & 1.15402494798059E-02_wp,-3.91472860045389E-04_wp,-1.85447400076914E-04_wp, &
      & 2.69372961745180E-05_wp,-3.83960703078053E-04_wp,-5.52056921452567E-04_wp, &
      & 7.04758808977056E-04_wp, 5.41396449946994E-04_wp, 1.62716519681097E-03_wp, &
      & 3.61628176171920E-04_wp,-1.63409127414125E-03_wp,-1.07250498654219E-03_wp, &
      & 3.31847658400096E-04_wp,-6.83598102195532E-04_wp, 5.75709000782720E-04_wp, &
      &-1.94547758751811E-04_wp,-6.02151183304568E-04_wp, 1.49179954780989E-03_wp, &
      &-2.93152544359381E-03_wp, 1.04205333697116E-03_wp, 1.48481072357888E-04_wp, &
      & 6.08528200623539E-04_wp, 4.57690583399911E-03_wp,-1.88737263402893E-02_wp, &
      & 8.75664181962324E-03_wp,-1.26486079338153E-02_wp,-7.10043290946803E-04_wp, &
      & 2.89620397193752E-06_wp, 6.06131016453163E-04_wp, 2.82293389952361E-06_wp, &
      & 2.86985578991511E-04_wp,-3.80262208043012E-03_wp,-1.56986199671825E-02_wp, &
      & 9.95289435445029E-03_wp, 1.08109379474349E-02_wp,-6.58426823850518E-04_wp, &
      &-2.72580344710504E-05_wp,-5.84203711930999E-04_wp, 1.67615675302901E-05_wp, &
      &-1.95672888720868E-04_wp,-3.63428057490533E-03_wp, 5.13544987508667E-03_wp, &
      & 1.02100793239369E-02_wp,-1.81555166175934E-02_wp,-3.91264029795358E-04_wp, &
      & 9.88365504780387E-06_wp,-5.75792938733272E-04_wp,-4.71112289725868E-05_wp, &
      & 5.46257757174947E-04_wp,-2.94624945675376E-04_wp,-2.59430166947037E-04_wp, &
      & 3.61628176171920E-04_wp, 2.43622668168914E-03_wp,-1.06123973390260E-03_wp, &
      &-1.65324946777566E-03_wp,-7.25296306806351E-05_wp, 1.80304070856658E-03_wp, &
      & 2.13517632377420E-03_wp,-1.81460080347714E-03_wp, 2.37543992326715E-03_wp, &
      & 6.76725101117992E-04_wp,-1.81030095889485E-03_wp,-1.71651898481651E-03_wp, &
      & 2.34703941229088E-03_wp,-2.97757430265402E-04_wp, 4.75605525497978E-03_wp, &
      &-1.99507322125130E-02_wp,-1.13575027794173E-02_wp, 4.49366063601718E-03_wp, &
      &-2.46172687351389E-04_wp,-6.32133131699997E-04_wp, 8.27214390433451E-05_wp, &
      &-1.20209486448711E-04_wp, 7.01975281929006E-04_wp, 5.46887348151783E-03_wp, &
      & 2.07465098916450E-02_wp, 1.12870324819613E-02_wp, 2.62360032993394E-03_wp, &
      & 3.30170019630099E-04_wp,-6.85792113938533E-04_wp, 6.13175091245618E-05_wp, &
      & 2.02457181360978E-04_wp, 6.65644753860844E-04_wp,-9.08181098068503E-03_wp, &
      & 1.33513181820170E-02_wp,-1.99268786228122E-02_wp,-2.34217202320808E-02_wp, &
      &-6.69813653944928E-04_wp,-4.38712775661282E-04_wp,-1.20006019568062E-04_wp, &
      & 1.02053886193372E-03_wp, 5.94007597988317E-04_wp,-2.74166581356306E-03_wp, &
      &-2.64919936937235E-03_wp,-1.63409127414125E-03_wp,-1.06123973390260E-03_wp, &
      & 4.26660636479651E-02_wp, 5.48867639267497E-03_wp, 2.54176254686534E-03_wp, &
      & 7.36832296027678E-03_wp,-4.03975654473191E-04_wp,-3.26942203750127E-04_wp, &
      & 3.91589546660145E-03_wp, 4.64342996080958E-03_wp, 1.33992841815537E-02_wp, &
      & 7.49785230431801E-04_wp, 8.17240065739435E-04_wp,-3.19238355446981E-04_wp, &
      & 2.23646244137205E-02_wp, 6.59548065857667E-02_wp, 8.68754738499316E-02_wp, &
      & 1.25040456174467E-03_wp, 8.43040832431269E-04_wp, 2.57139294590359E-03_wp, &
      & 4.69721320177459E-04_wp, 1.34881283254879E-03_wp,-1.20147515573048E-03_wp, &
      &-2.39915293849405E-02_wp, 7.91814684776250E-02_wp, 6.88501216429399E-02_wp, &
      &-1.27801136787652E-02_wp, 1.67318715110993E-03_wp,-2.26525736225408E-03_wp, &
      & 1.03312589060382E-04_wp, 1.25511859423495E-03_wp, 1.23392237280349E-03_wp, &
      &-1.35620731019085E-02_wp, 7.56618991833905E-02_wp, 3.15084854540211E-02_wp, &
      & 6.45053437127237E-02_wp,-9.38359904069108E-04_wp,-6.34043613155897E-04_wp, &
      & 3.14862311605270E-04_wp,-1.37787803489033E-03_wp,-2.24225399100554E-03_wp, &
      & 6.84153465016098E-03_wp, 5.77530559995348E-04_wp,-1.07250498654219E-03_wp, &
      &-1.65324946777566E-03_wp, 5.48867639267497E-03_wp, 3.04783663710111E-02_wp, &
      & 5.04606442096020E-03_wp,-4.66371208138313E-03_wp,-7.19052804954172E-03_wp, &
      & 3.26384647666875E-03_wp,-1.16648070070560E-02_wp,-1.74839240093656E-03_wp, &
      & 1.30572623508845E-03_wp,-5.03897420550892E-03_wp,-8.96763548845411E-04_wp, &
      & 9.84759792548497E-03_wp, 2.07224740193060E-02_wp, 6.70688097195503E-02_wp, &
      &-1.96459370841130E-02_wp, 7.81809041418509E-02_wp, 2.84420597616065E-03_wp, &
      & 3.81213248123500E-04_wp,-1.93389141426929E-03_wp, 4.82420477979741E-04_wp, &
      &-5.32810176688762E-04_wp, 1.83968711176355E-02_wp,-6.46585082512025E-02_wp, &
      & 3.76360397992066E-02_wp, 6.76726646968307E-02_wp,-2.56333530435510E-03_wp, &
      & 3.57609160918523E-04_wp,-2.02312098196774E-03_wp,-3.47311396599134E-04_wp, &
      &-5.79280405069415E-04_wp,-5.21056540703168E-03_wp, 3.66511396790259E-02_wp, &
      &-1.61959501492248E-02_wp, 4.22610487435226E-02_wp,-7.23150166636224E-04_wp, &
      &-6.96300015777415E-04_wp, 6.92000783727072E-04_wp, 2.14568368478460E-04_wp, &
      &-1.19565402301603E-03_wp,-5.31279761822322E-03_wp,-5.26976013457427E-05_wp, &
      & 3.31847658400096E-04_wp,-7.25296306806351E-05_wp, 2.54176254686534E-03_wp, &
      & 5.04606442096020E-03_wp, 2.84679587556318E-02_wp, 3.54894557191521E-03_wp, &
      &-7.72046463324211E-04_wp, 3.39812061949517E-03_wp, 1.80295732067388E-03_wp, &
      &-1.72699146503798E-03_wp, 2.08068962045090E-04_wp,-9.57025303579991E-04_wp, &
      &-6.72245961471031E-04_wp,-1.25000483023965E-02_wp,-7.48190813025657E-03_wp, &
      & 1.50194621101850E-02_wp,-9.51746257217108E-02_wp, 1.07105823451796E-02_wp, &
      & 2.74276798196848E-04_wp,-2.28477526111477E-03_wp,-1.54670546994761E-03_wp, &
      &-1.52509696471355E-03_wp,-8.34604782021266E-05_wp,-9.70631584445687E-03_wp, &
      & 5.89106335176432E-04_wp, 8.54634077370308E-02_wp,-6.23156571113132E-03_wp, &
      & 2.39237960359503E-04_wp,-1.86157192226641E-03_wp,-1.22914115729558E-03_wp, &
      & 1.62380339070375E-03_wp,-2.57899109827892E-04_wp,-9.88014923090613E-03_wp, &
      &-6.25144814937098E-03_wp, 8.60612099656805E-02_wp, 3.14099372304278E-03_wp, &
      & 3.78532516914335E-04_wp, 9.38872653962537E-04_wp,-1.19990637124712E-03_wp, &
      &-2.29785457537160E-03_wp, 4.74724238071550E-05_wp, 4.85658416213561E-03_wp, &
      &-1.70360657293966E-03_wp,-6.83598102195532E-04_wp, 1.80304070856658E-03_wp, &
      & 7.36832296027678E-03_wp,-4.66371208138313E-03_wp, 3.54894557191521E-03_wp, &
      & 3.34393781938852E-02_wp, 4.73561694057961E-03_wp,-9.20066557038704E-03_wp, &
      &-4.78340485526786E-03_wp,-5.05261443662806E-03_wp, 6.84791172570740E-04_wp, &
      & 2.03099825524059E-03_wp, 9.45917682510810E-03_wp, 3.78879388099217E-03_wp, &
      & 1.40822121394522E-02_wp, 7.85903417261118E-02_wp,-1.34573967461195E-02_wp, &
      & 5.53985700082034E-03_wp, 1.51665855913071E-03_wp, 4.92005695111929E-04_wp, &
      &-1.33814279600240E-03_wp, 9.60493988837844E-06_wp,-1.49543254126618E-03_wp, &
      &-1.19162210963848E-02_wp, 7.77031024200890E-02_wp,-3.00536616447891E-02_wp, &
      &-2.00732274111886E-03_wp, 1.10329669950783E-03_wp,-1.68670894926453E-04_wp, &
      & 1.46237940335675E-03_wp,-6.46986175423066E-04_wp, 1.68257298559314E-03_wp, &
      & 2.13600011676214E-02_wp, 3.25281161052191E-02_wp, 4.67381503377974E-02_wp, &
      &-1.03858523700383E-01_wp,-1.50740974288112E-03_wp, 3.82651994863863E-05_wp, &
      &-2.40313487407486E-03_wp, 3.71631530997041E-04_wp, 2.60157115829100E-03_wp, &
      & 9.50019121810167E-04_wp,-2.12815610709220E-04_wp, 5.75709000782720E-04_wp, &
      & 2.13517632377420E-03_wp,-4.03975654473191E-04_wp,-7.19052804954172E-03_wp, &
      &-7.72046463324211E-04_wp, 4.73561694057961E-03_wp, 4.18770813681427E-02_wp, &
      &-2.50290440454551E-03_wp, 7.65727418567875E-04_wp, 2.12624345393360E-03_wp, &
      &-4.97027877917978E-03_wp,-5.83472669844681E-03_wp, 5.51732011957964E-03_wp, &
      &-4.27826077946260E-04_wp,-8.82403162605877E-03_wp, 3.37743473393135E-02_wp, &
      &-3.36885096134562E-02_wp,-8.78793958800207E-02_wp,-1.20767003906759E-03_wp, &
      &-4.48556039118914E-04_wp,-1.81086224545010E-04_wp,-1.30225149564302E-03_wp, &
      &-1.72564923333566E-03_wp,-8.62878718054269E-03_wp,-3.71725125852225E-02_wp, &
      & 3.53696185434499E-02_wp,-9.55348684789787E-02_wp, 1.04090454647082E-03_wp, &
      &-3.59221400042112E-04_wp,-3.57354903215494E-04_wp, 1.79429741173164E-03_wp, &
      &-2.19262341202585E-03_wp, 2.19832369216930E-02_wp, 7.20780368313856E-02_wp, &
      &-7.10320337351804E-02_wp,-4.72446672147909E-02_wp,-2.04798561894034E-03_wp, &
      &-1.49699282601156E-03_wp, 1.86593826628326E-04_wp, 2.32423304833032E-03_wp, &
      & 4.22160666293566E-04_wp,-2.43299151738308E-03_wp, 1.06766405478892E-03_wp, &
      &-1.94547758751811E-04_wp,-1.81460080347714E-03_wp,-3.26942203750127E-04_wp, &
      & 3.26384647666875E-03_wp, 3.39812061949517E-03_wp,-9.20066557038704E-03_wp, &
      &-2.50290440454551E-03_wp, 1.18831170923388E-02_wp, 3.11417561356817E-03_wp, &
      &-1.79633667274163E-03_wp, 1.08428680718754E-03_wp, 2.95504975538613E-03_wp, &
      & 9.40583767877227E-04_wp,-8.18360385983116E-04_wp, 2.70527565881747E-05_wp, &
      &-3.42955794705075E-02_wp,-1.00765000114561E-02_wp, 3.53152003768071E-02_wp, &
      & 3.49311376473727E-04_wp,-7.42699380840104E-04_wp, 4.39431782527856E-06_wp, &
      & 3.27030474840866E-04_wp, 1.58749721680075E-03_wp, 5.59939284193674E-04_wp, &
      &-4.33600687024676E-02_wp,-1.57001763493095E-02_wp,-3.92643295028619E-02_wp, &
      & 4.03441445098763E-04_wp, 1.24040844923251E-03_wp, 2.73735765525048E-04_wp, &
      & 2.82584056419282E-04_wp,-2.02260087963652E-03_wp,-8.83664108597079E-04_wp, &
      & 3.13114544459643E-02_wp, 3.67817089890721E-02_wp, 5.10642880072443E-02_wp, &
      &-6.93371177426940E-04_wp,-1.52128659470246E-04_wp,-6.60971351200672E-05_wp, &
      &-2.04580162271832E-03_wp,-1.86857557022780E-03_wp,-9.06679526463874E-03_wp, &
      & 1.23297231618413E-03_wp,-6.02151183304568E-04_wp, 2.37543992326715E-03_wp, &
      & 3.91589546660145E-03_wp,-1.16648070070560E-02_wp, 1.80295732067388E-03_wp, &
      &-4.78340485526786E-03_wp, 7.65727418567875E-04_wp, 3.11417561356817E-03_wp, &
      & 1.53506396343270E-02_wp, 4.93458019575497E-04_wp, 1.16823409155740E-03_wp, &
      & 3.13052003010725E-04_wp,-1.86532900379925E-03_wp,-8.96110846963260E-03_wp, &
      &-9.72183490348135E-04_wp,-6.96677285675340E-02_wp,-2.42640812741529E-02_wp, &
      &-2.43105967018250E-02_wp,-1.88909058022181E-03_wp,-1.84620872265780E-03_wp, &
      & 8.40605401675121E-04_wp,-8.31609991216517E-04_wp, 1.57654392910010E-03_wp, &
      &-6.32238839105142E-04_wp, 6.52028457208331E-02_wp, 1.09811119074071E-04_wp, &
      &-2.52064987261408E-02_wp, 1.91383514382578E-03_wp,-1.05856014409679E-03_wp, &
      & 1.27412333131771E-03_wp, 3.48350771197466E-04_wp, 1.31355340280619E-03_wp, &
      & 5.25854327644141E-05_wp, 2.45048495855217E-02_wp,-2.56655009577627E-03_wp, &
      & 4.06754362288473E-02_wp,-4.24473085273045E-04_wp,-4.73872134127807E-04_wp, &
      & 6.75497357338132E-04_wp,-5.61485634191517E-04_wp,-1.60968883271809E-03_wp, &
      &-3.80607531989858E-04_wp,-1.24027234557119E-03_wp, 1.49179954780989E-03_wp, &
      & 6.76725101117992E-04_wp, 4.64342996080958E-03_wp,-1.74839240093656E-03_wp, &
      &-1.72699146503798E-03_wp,-5.05261443662806E-03_wp, 2.12624345393360E-03_wp, &
      &-1.79633667274163E-03_wp, 4.93458019575497E-04_wp, 5.57899823593500E-03_wp, &
      &-1.29897186515772E-04_wp,-7.80032096534039E-04_wp,-1.84892158652720E-04_wp, &
      &-9.62092557672180E-04_wp, 3.41667470650106E-05_wp,-2.14598915207327E-02_wp, &
      & 4.22313704080375E-02_wp,-1.51574111619925E-02_wp,-7.82803317697533E-04_wp, &
      & 8.27111347825122E-04_wp, 1.17930469715194E-03_wp, 5.54949096391098E-04_wp, &
      & 3.31242618197427E-04_wp,-2.19959469457473E-04_wp,-3.87427809378408E-03_wp, &
      & 4.92359767360958E-02_wp, 3.76026077555199E-04_wp,-1.50246094084959E-04_wp, &
      &-1.27777504135598E-03_wp,-1.02227177755387E-03_wp, 9.55873747395274E-04_wp, &
      &-4.76745213910508E-05_wp, 5.49263377585457E-04_wp,-6.37450174956073E-03_wp, &
      &-1.46067658017479E-02_wp,-4.95983018633469E-03_wp, 1.86494551430902E-04_wp, &
      &-2.30162533491548E-05_wp, 1.12528083885465E-04_wp, 5.83275218807750E-04_wp, &
      & 4.19280315359991E-04_wp,-2.37874651841337E-03_wp,-2.62975652749730E-03_wp, &
      &-2.93152544359381E-03_wp,-1.81030095889485E-03_wp, 1.33992841815537E-02_wp, &
      & 1.30572623508845E-03_wp, 2.08068962045090E-04_wp, 6.84791172570740E-04_wp, &
      &-4.97027877917978E-03_wp, 1.08428680718754E-03_wp, 1.16823409155740E-03_wp, &
      &-1.29897186515772E-04_wp, 9.94428221843759E-03_wp,-7.19846587011024E-05_wp, &
      &-1.93960408497296E-04_wp,-2.74856880867185E-03_wp, 1.27031368525096E-03_wp, &
      & 4.18745285288000E-02_wp, 2.72966013888841E-02_wp, 2.85915868985861E-02_wp, &
      & 1.52421575467950E-03_wp, 1.17133429597759E-03_wp,-5.70259377186292E-04_wp, &
      & 8.00501008899383E-04_wp,-5.93648368390666E-04_wp,-8.89167773592477E-04_wp, &
      & 4.01388520677539E-02_wp, 6.11538319289190E-03_wp,-3.43827040289076E-02_wp, &
      & 1.81283369635923E-03_wp,-6.03244529432152E-04_wp, 9.60330421290320E-04_wp, &
      & 6.66961474578595E-04_wp, 1.38916670200248E-04_wp,-8.40001714013650E-04_wp, &
      &-1.91461291146221E-02_wp, 5.62002262438172E-03_wp, 4.94693462764965E-02_wp, &
      & 1.30712345342775E-03_wp, 4.10116747608010E-04_wp, 9.58776573157464E-04_wp, &
      &-7.52025031397998E-04_wp,-1.22301380113751E-03_wp,-3.06390494210790E-04_wp, &
      & 6.81819221360930E-04_wp, 1.04205333697116E-03_wp,-1.71651898481651E-03_wp, &
      & 7.49785230431801E-04_wp,-5.03897420550892E-03_wp,-9.57025303579991E-04_wp, &
      & 2.03099825524059E-03_wp,-5.83472669844681E-03_wp, 2.95504975538613E-03_wp, &
      & 3.13052003010725E-04_wp,-7.80032096534039E-04_wp,-7.19846587011024E-05_wp, &
      & 6.02218778761440E-03_wp,-1.96369388859724E-04_wp, 4.00523412011083E-04_wp, &
      & 2.43415200594010E-05_wp,-1.51845049947398E-02_wp, 2.87327242622507E-02_wp, &
      &-9.49249765705680E-03_wp,-5.48713980681419E-04_wp, 5.54022846554013E-04_wp, &
      & 8.01687607586436E-04_wp, 3.77777202539557E-04_wp, 1.95523826272956E-04_wp, &
      & 6.68905535754955E-04_wp,-5.70282419272623E-03_wp,-3.25394854693152E-02_wp, &
      &-4.86038808135561E-03_wp,-6.16005794351641E-05_wp, 1.02574808923772E-03_wp, &
      & 4.82756606241754E-04_wp,-5.53868614429023E-04_wp,-4.61704646590128E-04_wp, &
      &-4.42109441645594E-04_wp, 1.09535528038150E-03_wp, 5.78076943435337E-02_wp, &
      &-1.68594875769994E-03_wp,-1.50326172865466E-04_wp, 5.22664015462132E-04_wp, &
      &-1.09131629766287E-03_wp,-1.82188558275468E-03_wp,-3.25003079162081E-05_wp, &
      & 3.36294433121553E-03_wp,-2.77984913399897E-03_wp, 1.48481072357888E-04_wp, &
      & 2.34703941229088E-03_wp, 8.17240065739435E-04_wp,-8.96763548845411E-04_wp, &
      &-6.72245961471031E-04_wp, 9.45917682510810E-03_wp, 5.51732011957964E-03_wp, &
      & 9.40583767877227E-04_wp,-1.86532900379925E-03_wp,-1.84892158652720E-04_wp, &
      &-1.93960408497296E-04_wp,-1.96369388859724E-04_wp, 1.14272521505678E-02_wp, &
      & 3.81891246340942E-03_wp, 4.46942267787944E-04_wp,-4.11175285560280E-03_wp, &
      & 9.67588448002750E-03_wp, 5.64940155335306E-02_wp, 1.60503540214145E-03_wp, &
      & 1.52655440518358E-04_wp,-3.37131360832177E-04_wp, 1.19451943984374E-03_wp, &
      & 1.48935194146030E-03_wp,-9.31069490872293E-04_wp, 7.87169218444636E-03_wp, &
      & 3.89347577592281E-03_wp,-4.72243223724727E-02_wp, 1.35155308784772E-03_wp, &
      &-2.00017008301476E-04_wp, 5.23748289028596E-04_wp, 8.65609593211022E-04_wp, &
      &-8.12160921509295E-04_wp, 1.03952145445761E-03_wp, 3.42404939671829E-02_wp, &
      &-2.44882129383441E-03_wp,-6.15913706401642E-02_wp,-1.62689399094805E-03_wp, &
      &-5.53268037298782E-04_wp,-1.20341681703358E-03_wp, 1.03007344334728E-03_wp, &
      & 1.53200056222919E-03_wp, 8.99503672589712E-03_wp, 6.15223438265065E-04_wp, &
      & 6.08528200623539E-04_wp,-2.97757430265402E-04_wp,-3.19238355446981E-04_wp, &
      & 9.84759792548497E-03_wp,-1.25000483023965E-02_wp, 3.78879388099217E-03_wp, &
      &-4.27826077946260E-04_wp,-8.18360385983116E-04_wp,-8.96110846963260E-03_wp, &
      &-9.62092557672180E-04_wp,-2.74856880867185E-03_wp, 4.00523412011083E-04_wp, &
      & 3.81891246340942E-03_wp, 1.48383546641407E-02_wp, 2.14569885077566E-04_wp, &
      & 2.60228009071659E-02_wp, 4.55480551778864E-02_wp, 3.11414586934289E-02_wp, &
      & 1.24763488181776E-03_wp, 1.53320635726693E-03_wp, 7.48314603836746E-08_wp, &
      & 1.22592371170191E-03_wp,-1.03037430243239E-04_wp, 6.85864626229303E-04_wp, &
      &-3.93198699151587E-02_wp,-4.33959140312128E-02_wp, 3.84544626969406E-02_wp, &
      &-1.77253159605672E-03_wp, 1.76953567734716E-03_wp,-2.92690161508619E-04_wp, &
      &-1.42187034413862E-03_wp,-3.56619680821680E-04_wp, 2.40006331974931E-04_wp, &
      & 4.74589429392484E-02_wp,-2.86025601155193E-02_wp,-2.74869037414869E-02_wp, &
      &-1.81297319011913E-03_wp,-9.90759348859760E-04_wp,-3.97806973349798E-04_wp, &
      & 1.22820555419493E-03_wp, 4.56941288207337E-04_wp,-1.11688855965218E-03_wp, &
      & 7.17489050802549E-03_wp, 4.57690583399911E-03_wp, 4.75605525497978E-03_wp, &
      & 2.23646244137205E-02_wp, 2.07224740193060E-02_wp,-7.48190813025657E-03_wp, &
      & 1.40822121394522E-02_wp,-8.82403162605877E-03_wp, 2.70527565881747E-05_wp, &
      &-9.72183490348135E-04_wp, 3.41667470650106E-05_wp, 1.27031368525096E-03_wp, &
      & 2.43415200594010E-05_wp, 4.46942267787944E-04_wp, 2.14569885077566E-04_wp, &
      & 9.79953185022663E-01_wp,-6.56913535165527E-02_wp,-4.15746208633432E-02_wp, &
      &-4.42906788242035E-02_wp,-3.26497777864028E-03_wp,-3.02957050163200E-03_wp, &
      & 1.03160173383718E-03_wp,-2.06224126752306E-03_wp, 1.29590186401169E-03_wp, &
      &-1.73578590097797E-03_wp,-4.78755981358030E-03_wp,-4.80964039210149E-03_wp, &
      &-2.26550292609638E-02_wp, 1.03244385712820E-03_wp, 8.28645456652198E-04_wp, &
      & 8.36183698869690E-04_wp, 3.56016852022734E-04_wp,-2.26721397739033E-03_wp, &
      &-4.67936049598523E-04_wp,-2.15316812457995E-02_wp,-5.04588128389799E-03_wp, &
      & 2.25098483451279E-03_wp, 2.05812072051300E-03_wp, 5.05766248283278E-04_wp, &
      & 7.91454294590990E-04_wp, 7.16816433450695E-04_wp, 1.18902455312510E-03_wp, &
      & 3.20046292353588E-02_wp,-1.15284447884868E-02_wp,-1.88737263402893E-02_wp, &
      &-1.99507322125130E-02_wp, 6.59548065857667E-02_wp, 6.70688097195503E-02_wp, &
      & 1.50194621101850E-02_wp, 7.85903417261118E-02_wp, 3.37743473393135E-02_wp, &
      &-3.42955794705075E-02_wp,-6.96677285675340E-02_wp,-2.14598915207327E-02_wp, &
      & 4.18745285288000E-02_wp,-1.51845049947398E-02_wp,-4.11175285560280E-03_wp, &
      & 2.60228009071659E-02_wp,-6.56913535165527E-02_wp, 8.97832343104439E-01_wp, &
      &-9.01190569007224E-03_wp,-1.45144645636852E-02_wp, 1.55005079703965E-02_wp, &
      & 1.20367234816796E-02_wp,-1.34556376466509E-02_wp,-3.44948217516065E-04_wp, &
      &-2.38477695614042E-02_wp, 5.97412493903286E-03_wp, 3.58652467089850E-03_wp, &
      &-1.13563354350285E-02_wp,-9.59168996535537E-03_wp,-6.71551375240948E-04_wp, &
      & 2.30409327479318E-03_wp, 4.33495391693732E-04_wp,-1.03517926424578E-03_wp, &
      &-9.32392910000505E-04_wp,-1.75636665246337E-02_wp,-4.80394095409307E-02_wp, &
      &-4.04081526462186E-03_wp, 2.91067062791004E-02_wp, 5.41082241160409E-03_wp, &
      & 2.05070858833908E-03_wp, 2.29339909979139E-03_wp, 4.07757821937570E-04_wp, &
      &-4.27488088446933E-04_wp, 2.32606695048053E-02_wp,-1.66995151518806E-02_wp, &
      & 8.75664181962324E-03_wp,-1.13575027794173E-02_wp, 8.68754738499316E-02_wp, &
      &-1.96459370841130E-02_wp,-9.51746257217108E-02_wp,-1.34573967461195E-02_wp, &
      &-3.36885096134562E-02_wp,-1.00765000114561E-02_wp,-2.42640812741529E-02_wp, &
      & 4.22313704080375E-02_wp, 2.72966013888841E-02_wp, 2.87327242622507E-02_wp, &
      & 9.67588448002750E-03_wp, 4.55480551778864E-02_wp,-4.15746208633432E-02_wp, &
      &-9.01190569007224E-03_wp, 8.98949522434250E-01_wp,-5.85249818311979E-03_wp, &
      &-4.20396123694783E-05_wp, 2.32696964109612E-02_wp, 1.40699101827117E-02_wp, &
      & 1.57206108208981E-02_wp, 4.82043997286318E-05_wp, 1.29832452985130E-03_wp, &
      &-1.01189621165397E-02_wp, 1.64501682042587E-02_wp,-1.13646125614395E-02_wp, &
      &-6.58322439410475E-04_wp,-9.95021813207843E-05_wp,-1.69154717257235E-03_wp, &
      & 8.98861322776682E-04_wp,-2.17156562704802E-03_wp, 1.74725131414714E-03_wp, &
      &-1.23060391666245E-02_wp, 1.53478173145747E-02_wp,-5.75999068005425E-03_wp, &
      & 6.00827861785124E-04_wp, 6.76325363903056E-04_wp,-1.66014282991114E-03_wp, &
      &-2.79625792411043E-04_wp, 2.06014446932550E-03_wp, 2.10966188523540E-02_wp, &
      &-1.96285390992293E-02_wp,-1.26486079338153E-02_wp, 4.49366063601718E-03_wp, &
      & 1.25040456174467E-03_wp, 7.81809041418509E-02_wp, 1.07105823451796E-02_wp, &
      & 5.53985700082034E-03_wp,-8.78793958800207E-02_wp, 3.53152003768071E-02_wp, &
      &-2.43105967018250E-02_wp,-1.51574111619925E-02_wp, 2.85915868985861E-02_wp, &
      &-9.49249765705680E-03_wp, 5.64940155335306E-02_wp, 3.11414586934289E-02_wp, &
      &-4.42906788242035E-02_wp,-1.45144645636852E-02_wp,-5.85249818311979E-03_wp, &
      & 9.09568901510487E-01_wp, 2.36486322778074E-02_wp,-3.50001068677793E-04_wp, &
      &-9.07751827589855E-03_wp, 1.23464513582589E-02_wp, 1.66020669497448E-02_wp, &
      &-2.28312252496096E-02_wp, 1.33634035581138E-02_wp, 1.41764609986243E-03_wp, &
      &-5.45979916231767E-02_wp, 5.34400854366441E-03_wp,-4.96044569603787E-06_wp, &
      & 2.40201841492500E-03_wp, 2.48198884031881E-03_wp,-2.72121814524801E-03_wp, &
      & 1.34763066806292E-02_wp, 1.29062067504150E-02_wp,-9.80190748894087E-03_wp, &
      &-5.04907733624699E-03_wp,-2.40085902556875E-03_wp,-1.16119924081091E-03_wp, &
      &-5.81093006975118E-04_wp, 2.45027490200433E-03_wp, 1.39966952488796E-03_wp, &
      & 1.14006433275913E-03_wp,-8.09390428732359E-04_wp,-7.10043290946803E-04_wp, &
      &-2.46172687351389E-04_wp, 8.43040832431269E-04_wp, 2.84420597616065E-03_wp, &
      & 2.74276798196848E-04_wp, 1.51665855913071E-03_wp,-1.20767003906759E-03_wp, &
      & 3.49311376473727E-04_wp,-1.88909058022181E-03_wp,-7.82803317697533E-04_wp, &
      & 1.52421575467950E-03_wp,-5.48713980681419E-04_wp, 1.60503540214145E-03_wp, &
      & 1.24763488181776E-03_wp,-3.26497777864028E-03_wp, 1.55005079703965E-02_wp, &
      &-4.20396123694783E-05_wp, 2.36486322778074E-02_wp, 9.18511119490073E-04_wp, &
      & 2.19078907283795E-04_wp,-4.80302880422056E-04_wp, 3.36883951800646E-04_wp, &
      & 1.96443072182944E-05_wp,-9.42305899811192E-04_wp,-2.89537305206576E-04_wp, &
      &-1.28022878279291E-03_wp,-5.00721188868637E-03_wp, 1.98157365462310E-04_wp, &
      & 6.86048439806489E-05_wp, 1.12897189556710E-04_wp, 8.04397669198610E-05_wp, &
      &-1.60993997929097E-04_wp, 9.24348836836812E-04_wp,-1.64629006596351E-03_wp, &
      &-1.88908493067724E-03_wp,-5.38795740743591E-04_wp, 5.99810932074870E-05_wp, &
      & 1.05424830110441E-05_wp, 4.34016120608557E-05_wp, 1.16049685662805E-04_wp, &
      & 5.66604383383437E-05_wp, 1.14147229303145E-03_wp,-6.27932511020678E-04_wp, &
      & 2.89620397193752E-06_wp,-6.32133131699997E-04_wp, 2.57139294590359E-03_wp, &
      & 3.81213248123500E-04_wp,-2.28477526111477E-03_wp, 4.92005695111929E-04_wp, &
      &-4.48556039118914E-04_wp,-7.42699380840104E-04_wp,-1.84620872265780E-03_wp, &
      & 8.27111347825122E-04_wp, 1.17133429597759E-03_wp, 5.54022846554013E-04_wp, &
      & 1.52655440518358E-04_wp, 1.53320635726693E-03_wp,-3.02957050163200E-03_wp, &
      & 1.20367234816796E-02_wp, 2.32696964109612E-02_wp,-3.50001068677793E-04_wp, &
      & 2.19078907283795E-04_wp, 7.86474390598036E-04_wp, 1.81027823094628E-04_wp, &
      & 4.09324704116321E-04_wp,-3.38483160385702E-04_wp, 7.13341494118634E-04_wp, &
      &-2.49338759770268E-03_wp, 1.68735866975485E-04_wp,-7.41309576022951E-04_wp, &
      &-5.88299190152288E-05_wp, 5.59832212371159E-05_wp,-6.75477375787827E-05_wp, &
      & 1.76396292110411E-05_wp,-1.27813548615559E-04_wp,-3.12515170166521E-04_wp, &
      &-3.56717280611371E-03_wp, 3.00321001380880E-04_wp,-6.77069550327594E-04_wp, &
      & 1.50498405466485E-04_wp, 8.35670166748950E-05_wp,-1.78507090937897E-05_wp, &
      & 5.03854201101717E-06_wp, 9.88365521261470E-05_wp,-3.02272998088046E-04_wp, &
      & 1.21450766756529E-04_wp, 6.06131016453163E-04_wp, 8.27214390433451E-05_wp, &
      & 4.69721320177459E-04_wp,-1.93389141426929E-03_wp,-1.54670546994761E-03_wp, &
      &-1.33814279600240E-03_wp,-1.81086224545010E-04_wp, 4.39431782527856E-06_wp, &
      & 8.40605401675121E-04_wp, 1.17930469715194E-03_wp,-5.70259377186292E-04_wp, &
      & 8.01687607586436E-04_wp,-3.37131360832177E-04_wp, 7.48314603836746E-08_wp, &
      & 1.03160173383718E-03_wp,-1.34556376466509E-02_wp, 1.40699101827117E-02_wp, &
      &-9.07751827589855E-03_wp,-4.80302880422056E-04_wp, 1.81027823094628E-04_wp, &
      & 5.18247582373794E-04_wp, 1.22109199673484E-04_wp, 1.93563736986084E-04_wp, &
      & 6.94512792359044E-04_wp,-7.43041114473141E-04_wp, 1.90087784391794E-03_wp, &
      & 1.41762342890863E-03_wp,-8.72768819207063E-05_wp,-6.27592023472381E-05_wp, &
      &-9.94397095236854E-05_wp, 1.69999823758605E-05_wp, 1.02867249071339E-05_wp, &
      & 6.86263434466501E-04_wp, 1.12323534619301E-03_wp, 1.84122040819062E-03_wp, &
      &-1.05685564533458E-03_wp,-7.41528292047566E-05_wp,-5.64959811538514E-06_wp, &
      &-9.38849407025514E-05_wp,-6.44076461307671E-05_wp, 4.02013059080886E-05_wp, &
      & 7.61032562708332E-04_wp,-6.28894926244606E-04_wp, 2.82293389952361E-06_wp, &
      &-1.20209486448711E-04_wp, 1.34881283254879E-03_wp, 4.82420477979741E-04_wp, &
      &-1.52509696471355E-03_wp, 9.60493988837844E-06_wp,-1.30225149564302E-03_wp, &
      & 3.27030474840866E-04_wp,-8.31609991216517E-04_wp, 5.54949096391098E-04_wp, &
      & 8.00501008899383E-04_wp, 3.77777202539557E-04_wp, 1.19451943984374E-03_wp, &
      & 1.22592371170191E-03_wp,-2.06224126752306E-03_wp,-3.44948217516065E-04_wp, &
      & 1.57206108208981E-02_wp, 1.23464513582589E-02_wp, 3.36883951800646E-04_wp, &
      & 4.09324704116321E-04_wp, 1.22109199673484E-04_wp, 4.58251142523757E-04_wp, &
      & 2.38692760788488E-04_wp,-7.45009963602784E-04_wp,-8.81995629091355E-04_wp, &
      & 4.45232108061955E-04_wp,-3.31802835554324E-03_wp, 9.72796042636698E-05_wp, &
      & 4.59931850098103E-06_wp, 8.59895681225748E-06_wp, 8.83186095660335E-05_wp, &
      &-1.38666929038812E-04_wp, 9.28623086750384E-04_wp, 1.19596404451922E-04_wp, &
      & 5.87522199748577E-06_wp,-2.03335470252935E-03_wp,-4.17552021346788E-05_wp, &
      &-5.09405249893958E-06_wp,-6.69138499163433E-05_wp, 5.12364190514789E-05_wp, &
      & 1.03787374695025E-04_wp,-4.79144450852178E-04_wp,-7.62484937973928E-05_wp, &
      & 2.86985578991511E-04_wp, 7.01975281929006E-04_wp,-1.20147515573048E-03_wp, &
      &-5.32810176688762E-04_wp,-8.34604782021266E-05_wp,-1.49543254126618E-03_wp, &
      &-1.72564923333566E-03_wp, 1.58749721680075E-03_wp, 1.57654392910010E-03_wp, &
      & 3.31242618197427E-04_wp,-5.93648368390666E-04_wp, 1.95523826272956E-04_wp, &
      & 1.48935194146030E-03_wp,-1.03037430243239E-04_wp, 1.29590186401169E-03_wp, &
      &-2.38477695614042E-02_wp, 4.82043997286318E-05_wp, 1.66020669497448E-02_wp, &
      & 1.96443072182944E-05_wp,-3.38483160385702E-04_wp, 1.93563736986084E-04_wp, &
      & 2.38692760788488E-04_wp, 9.55996160124409E-04_wp,-2.25331678846492E-03_wp, &
      & 1.15282575115815E-03_wp, 1.46723806272380E-03_wp,-2.97238116407506E-03_wp, &
      & 1.74422156274755E-04_wp,-9.53909049759322E-05_wp, 4.56054718637480E-05_wp, &
      & 1.19023687658630E-04_wp,-4.22374644344222E-05_wp, 2.15549509066777E-03_wp, &
      & 5.25101958893515E-03_wp, 1.29756344050563E-05_wp,-2.17978866795573E-03_wp, &
      &-2.95874145259961E-04_wp,-1.25680696991265E-04_wp,-1.08710889283059E-04_wp, &
      & 5.08755648891243E-05_wp, 3.70646570213923E-05_wp,-1.30082860861281E-03_wp, &
      &-8.00507020680877E-03_wp,-3.80262208043012E-03_wp, 5.46887348151783E-03_wp, &
      &-2.39915293849405E-02_wp, 1.83968711176355E-02_wp,-9.70631584445687E-03_wp, &
      &-1.19162210963848E-02_wp,-8.62878718054269E-03_wp, 5.59939284193674E-04_wp, &
      &-6.32238839105142E-04_wp,-2.19959469457473E-04_wp,-8.89167773592477E-04_wp, &
      & 6.68905535754955E-04_wp,-9.31069490872293E-04_wp, 6.85864626229303E-04_wp, &
      &-1.73578590097797E-03_wp, 5.97412493903286E-03_wp, 1.29832452985130E-03_wp, &
      &-2.28312252496096E-02_wp,-9.42305899811192E-04_wp, 7.13341494118634E-04_wp, &
      & 6.94512792359044E-04_wp,-7.45009963602784E-04_wp,-2.25331678846492E-03_wp, &
      & 9.80974511856617E-01_wp, 6.77159511266575E-02_wp, 3.31547266836938E-02_wp, &
      &-4.62531679479647E-02_wp, 3.56741815719079E-03_wp,-2.66711332001603E-03_wp, &
      & 1.45224189985259E-03_wp, 1.78354903062082E-03_wp, 1.32108871221265E-03_wp, &
      &-3.02557120515403E-03_wp, 1.59240225965459E-02_wp, 1.87604065635367E-02_wp, &
      & 4.69837569301324E-03_wp,-1.69113411385192E-03_wp,-6.61974620886243E-04_wp, &
      &-1.25193134896855E-03_wp,-1.80059796257016E-03_wp,-3.51490801006542E-04_wp, &
      &-3.51696354973964E-02_wp,-1.18379012605914E-02_wp,-1.56986199671825E-02_wp, &
      & 2.07465098916450E-02_wp, 7.91814684776250E-02_wp,-6.46585082512025E-02_wp, &
      & 5.89106335176432E-04_wp, 7.77031024200890E-02_wp,-3.71725125852225E-02_wp, &
      &-4.33600687024676E-02_wp, 6.52028457208331E-02_wp,-3.87427809378408E-03_wp, &
      & 4.01388520677539E-02_wp,-5.70282419272623E-03_wp, 7.87169218444636E-03_wp, &
      &-3.93198699151587E-02_wp,-4.78755981358030E-03_wp, 3.58652467089850E-03_wp, &
      &-1.01189621165397E-02_wp, 1.33634035581138E-02_wp,-2.89537305206576E-04_wp, &
      &-2.49338759770268E-03_wp,-7.43041114473141E-04_wp,-8.81995629091355E-04_wp, &
      & 1.15282575115815E-03_wp, 6.77159511266575E-02_wp, 8.93041610503373E-01_wp, &
      &-4.10359825688933E-03_wp, 1.29601523129638E-02_wp, 1.67095113530953E-02_wp, &
      &-1.36647926268428E-02_wp, 1.27251154417897E-02_wp,-3.25733892456355E-05_wp, &
      & 2.26849602998123E-02_wp, 9.66633023541558E-03_wp,-5.19297868799919E-03_wp, &
      &-3.84633749324155E-02_wp,-1.24613040121676E-02_wp, 1.38676485413577E-03_wp, &
      &-6.88153886369985E-04_wp, 2.13804102298870E-03_wp, 3.54546470887379E-03_wp, &
      &-6.75694770587501E-06_wp,-1.39407590739093E-02_wp,-1.74855485990935E-02_wp, &
      & 9.95289435445029E-03_wp, 1.12870324819613E-02_wp, 6.88501216429399E-02_wp, &
      & 3.76360397992066E-02_wp, 8.54634077370308E-02_wp,-3.00536616447891E-02_wp, &
      & 3.53696185434499E-02_wp,-1.57001763493095E-02_wp, 1.09811119074071E-04_wp, &
      & 4.92359767360958E-02_wp, 6.11538319289190E-03_wp,-3.25394854693152E-02_wp, &
      & 3.89347577592281E-03_wp,-4.33959140312128E-02_wp,-4.80964039210149E-03_wp, &
      &-1.13563354350285E-02_wp, 1.64501682042587E-02_wp, 1.41764609986243E-03_wp, &
      &-1.28022878279291E-03_wp, 1.68735866975485E-04_wp, 1.90087784391794E-03_wp, &
      & 4.45232108061955E-04_wp, 1.46723806272380E-03_wp, 3.31547266836938E-02_wp, &
      &-4.10359825688933E-03_wp, 9.08117769603495E-01_wp, 9.12873283600103E-03_wp, &
      &-3.41667105887394E-04_wp,-2.27240427085580E-02_wp,-1.60656013699057E-02_wp, &
      & 1.70173883505189E-02_wp, 3.01970480128187E-04_wp, 1.87221098337115E-02_wp, &
      &-2.32165734653952E-02_wp,-2.18459102876546E-02_wp,-3.26536302874007E-02_wp, &
      & 1.96920444748069E-03_wp, 5.26310641005283E-04_wp,-1.60834820233167E-05_wp, &
      & 3.65950247565885E-03_wp, 4.03502323498763E-03_wp, 2.24022500427106E-02_wp, &
      & 2.22060923985087E-02_wp, 1.08109379474349E-02_wp, 2.62360032993394E-03_wp, &
      &-1.27801136787652E-02_wp, 6.76726646968307E-02_wp,-6.23156571113132E-03_wp, &
      &-2.00732274111886E-03_wp,-9.55348684789787E-02_wp,-3.92643295028619E-02_wp, &
      &-2.52064987261408E-02_wp, 3.76026077555199E-04_wp,-3.43827040289076E-02_wp, &
      &-4.86038808135561E-03_wp,-4.72243223724727E-02_wp, 3.84544626969406E-02_wp, &
      &-2.26550292609638E-02_wp,-9.59168996535537E-03_wp,-1.13646125614395E-02_wp, &
      &-5.45979916231767E-02_wp,-5.00721188868637E-03_wp,-7.41309576022951E-04_wp, &
      & 1.41762342890863E-03_wp,-3.31802835554324E-03_wp,-2.97238116407506E-03_wp, &
      &-4.62531679479647E-02_wp, 1.29601523129638E-02_wp, 9.12873283600103E-03_wp, &
      & 9.09645248969580E-01_wp,-2.27265262407864E-02_wp,-5.32299199307859E-04_wp, &
      &-1.00495621815063E-02_wp,-1.37190880213787E-02_wp, 1.76759686357001E-02_wp, &
      & 1.35404942972661E-02_wp, 8.14988231623697E-03_wp,-1.12416390562587E-02_wp, &
      &-3.84377541807435E-03_wp,-1.75012611402543E-03_wp,-9.81189201589842E-04_wp, &
      &-3.98266659773968E-04_wp, 2.71148530315096E-03_wp, 1.45446242844250E-03_wp, &
      &-1.30001097214208E-03_wp,-8.41884382085824E-04_wp,-6.58426823850518E-04_wp, &
      & 3.30170019630099E-04_wp, 1.67318715110993E-03_wp,-2.56333530435510E-03_wp, &
      & 2.39237960359503E-04_wp, 1.10329669950783E-03_wp, 1.04090454647082E-03_wp, &
      & 4.03441445098763E-04_wp, 1.91383514382578E-03_wp,-1.50246094084959E-04_wp, &
      & 1.81283369635923E-03_wp,-6.16005794351641E-05_wp, 1.35155308784772E-03_wp, &
      &-1.77253159605672E-03_wp, 1.03244385712820E-03_wp,-6.71551375240948E-04_wp, &
      &-6.58322439410475E-04_wp, 5.34400854366441E-03_wp, 1.98157365462310E-04_wp, &
      &-5.88299190152288E-05_wp,-8.72768819207063E-05_wp, 9.72796042636698E-05_wp, &
      & 1.74422156274755E-04_wp, 3.56741815719079E-03_wp, 1.67095113530953E-02_wp, &
      &-3.41667105887394E-04_wp,-2.27265262407864E-02_wp, 9.20468008197760E-04_wp, &
      &-2.49874089186393E-04_wp, 5.11060765003657E-04_wp, 3.55781691260968E-04_wp, &
      &-2.00686966103915E-05_wp,-1.11069365429560E-03_wp,-1.77595705940923E-03_wp, &
      &-9.69766591163357E-04_wp, 1.79487344769638E-03_wp, 1.18227059103687E-04_wp, &
      & 2.23628169815579E-05_wp, 9.59540695118697E-05_wp,-6.27125251361419E-06_wp, &
      &-7.33810769916618E-05_wp, 9.03676498870593E-04_wp, 6.54018323036404E-04_wp, &
      &-2.72580344710504E-05_wp,-6.85792113938533E-04_wp,-2.26525736225408E-03_wp, &
      & 3.57609160918523E-04_wp,-1.86157192226641E-03_wp,-1.68670894926453E-04_wp, &
      &-3.59221400042112E-04_wp, 1.24040844923251E-03_wp,-1.05856014409679E-03_wp, &
      &-1.27777504135598E-03_wp,-6.03244529432152E-04_wp, 1.02574808923772E-03_wp, &
      &-2.00017008301476E-04_wp, 1.76953567734716E-03_wp, 8.28645456652198E-04_wp, &
      & 2.30409327479318E-03_wp,-9.95021813207843E-05_wp,-4.96044569603787E-06_wp, &
      & 6.86048439806489E-05_wp, 5.59832212371159E-05_wp,-6.27592023472381E-05_wp, &
      & 4.59931850098103E-06_wp,-9.53909049759322E-05_wp,-2.66711332001603E-03_wp, &
      &-1.36647926268428E-02_wp,-2.27240427085580E-02_wp,-5.32299199307859E-04_wp, &
      &-2.49874089186393E-04_wp, 8.05224008513454E-04_wp, 2.12701081422343E-04_wp, &
      &-4.33732331900280E-04_wp,-3.76373327444553E-04_wp,-1.92622491079819E-03_wp, &
      & 3.14448950646073E-03_wp, 3.64600111382412E-03_wp, 2.20063557970678E-03_wp, &
      &-1.22466148179130E-04_wp,-1.14817581188053E-05_wp,-6.53524280632076E-05_wp, &
      &-2.32450086452956E-04_wp,-1.56880573246973E-04_wp,-5.83171540616389E-04_wp, &
      &-9.79735922630597E-05_wp,-5.84203711930999E-04_wp, 6.13175091245618E-05_wp, &
      & 1.03312589060382E-04_wp,-2.02312098196774E-03_wp,-1.22914115729558E-03_wp, &
      & 1.46237940335675E-03_wp,-3.57354903215494E-04_wp, 2.73735765525048E-04_wp, &
      & 1.27412333131771E-03_wp,-1.02227177755387E-03_wp, 9.60330421290320E-04_wp, &
      & 4.82756606241754E-04_wp, 5.23748289028596E-04_wp,-2.92690161508619E-04_wp, &
      & 8.36183698869690E-04_wp, 4.33495391693732E-04_wp,-1.69154717257235E-03_wp, &
      & 2.40201841492500E-03_wp, 1.12897189556710E-04_wp,-6.75477375787827E-05_wp, &
      &-9.94397095236854E-05_wp, 8.59895681225748E-06_wp, 4.56054718637480E-05_wp, &
      & 1.45224189985259E-03_wp, 1.27251154417897E-02_wp,-1.60656013699057E-02_wp, &
      &-1.00495621815063E-02_wp, 5.11060765003657E-04_wp, 2.12701081422343E-04_wp, &
      & 5.85429402822957E-04_wp,-1.44929109080051E-04_wp, 1.19686671170248E-04_wp, &
      &-1.26135270204616E-03_wp, 3.43298532897802E-04_wp, 3.83817970980841E-06_wp, &
      & 2.16308617823632E-03_wp, 1.39871480504945E-05_wp,-1.14181893240403E-05_wp, &
      & 6.32963494893182E-05_wp,-6.58156831460961E-05_wp,-1.36918692981879E-04_wp, &
      &-5.86376985234641E-04_wp,-7.38014920537391E-04_wp, 1.67615675302901E-05_wp, &
      & 2.02457181360978E-04_wp, 1.25511859423495E-03_wp,-3.47311396599134E-04_wp, &
      & 1.62380339070375E-03_wp,-6.46986175423066E-04_wp, 1.79429741173164E-03_wp, &
      & 2.82584056419282E-04_wp, 3.48350771197466E-04_wp, 9.55873747395274E-04_wp, &
      & 6.66961474578595E-04_wp,-5.53868614429023E-04_wp, 8.65609593211022E-04_wp, &
      &-1.42187034413862E-03_wp, 3.56016852022734E-04_wp,-1.03517926424578E-03_wp, &
      & 8.98861322776682E-04_wp, 2.48198884031881E-03_wp, 8.04397669198610E-05_wp, &
      & 1.76396292110411E-05_wp, 1.69999823758605E-05_wp, 8.83186095660335E-05_wp, &
      & 1.19023687658630E-04_wp, 1.78354903062082E-03_wp,-3.25733892456355E-05_wp, &
      & 1.70173883505189E-02_wp,-1.37190880213787E-02_wp, 3.55781691260968E-04_wp, &
      &-4.33732331900280E-04_wp,-1.44929109080051E-04_wp, 5.39623595018820E-04_wp, &
      &-2.59213800285309E-04_wp,-4.02492245769748E-05_wp,-2.31816357724689E-03_wp, &
      &-6.97563672919024E-04_wp,-1.10536488762762E-03_wp, 1.02416470270356E-04_wp, &
      & 4.30064827780183E-05_wp, 1.06480846519756E-05_wp, 5.28437067393878E-05_wp, &
      & 9.07005439053863E-05_wp,-5.20172683811391E-04_wp, 1.67408952638943E-04_wp, &
      &-1.95672888720868E-04_wp, 6.65644753860844E-04_wp, 1.23392237280349E-03_wp, &
      &-5.79280405069415E-04_wp,-2.57899109827892E-04_wp, 1.68257298559314E-03_wp, &
      &-2.19262341202585E-03_wp,-2.02260087963652E-03_wp, 1.31355340280619E-03_wp, &
      &-4.76745213910508E-05_wp, 1.38916670200248E-04_wp,-4.61704646590128E-04_wp, &
      &-8.12160921509295E-04_wp,-3.56619680821680E-04_wp,-2.26721397739033E-03_wp, &
      &-9.32392910000505E-04_wp,-2.17156562704802E-03_wp,-2.72121814524801E-03_wp, &
      &-1.60993997929097E-04_wp,-1.27813548615559E-04_wp, 1.02867249071339E-05_wp, &
      &-1.38666929038812E-04_wp,-4.22374644344222E-05_wp, 1.32108871221265E-03_wp, &
      & 2.26849602998123E-02_wp, 3.01970480128187E-04_wp, 1.76759686357001E-02_wp, &
      &-2.00686966103915E-05_wp,-3.76373327444553E-04_wp, 1.19686671170248E-04_wp, &
      &-2.59213800285309E-04_wp, 9.32898777817762E-04_wp, 1.33666846700373E-03_wp, &
      &-3.22194483116891E-04_wp,-4.38986218965655E-03_wp,-8.07406407670555E-04_wp, &
      & 5.55356868900987E-06_wp,-6.21201449865174E-05_wp, 9.70425331312997E-05_wp, &
      & 2.32661759562063E-04_wp, 3.40426800813876E-05_wp,-1.20825923162660E-03_wp, &
      & 2.55413860729764E-03_wp,-3.63428057490533E-03_wp,-9.08181098068503E-03_wp, &
      &-1.35620731019085E-02_wp,-5.21056540703168E-03_wp,-9.88014923090613E-03_wp, &
      & 2.13600011676214E-02_wp, 2.19832369216930E-02_wp,-8.83664108597079E-04_wp, &
      & 5.25854327644141E-05_wp, 5.49263377585457E-04_wp,-8.40001714013650E-04_wp, &
      &-4.42109441645594E-04_wp, 1.03952145445761E-03_wp, 2.40006331974931E-04_wp, &
      &-4.67936049598523E-04_wp,-1.75636665246337E-02_wp, 1.74725131414714E-03_wp, &
      & 1.34763066806292E-02_wp, 9.24348836836812E-04_wp,-3.12515170166521E-04_wp, &
      & 6.86263434466501E-04_wp, 9.28623086750384E-04_wp, 2.15549509066777E-03_wp, &
      &-3.02557120515403E-03_wp, 9.66633023541558E-03_wp, 1.87221098337115E-02_wp, &
      & 1.35404942972661E-02_wp,-1.11069365429560E-03_wp,-1.92622491079819E-03_wp, &
      &-1.26135270204616E-03_wp,-4.02492245769748E-05_wp, 1.33666846700373E-03_wp, &
      & 9.80364656316398E-01_wp,-2.18477518470186E-02_wp, 3.31124366005085E-02_wp, &
      & 7.98537612496825E-02_wp, 1.97816816411445E-03_wp, 7.95287884239188E-04_wp, &
      & 1.47614735664917E-03_wp,-3.09112108744030E-03_wp,-3.26315136987568E-03_wp, &
      & 1.01661826271254E-02_wp, 1.44077180987680E-02_wp, 5.13544987508667E-03_wp, &
      & 1.33513181820170E-02_wp, 7.56618991833905E-02_wp, 3.66511396790259E-02_wp, &
      &-6.25144814937098E-03_wp, 3.25281161052191E-02_wp, 7.20780368313856E-02_wp, &
      & 3.13114544459643E-02_wp, 2.45048495855217E-02_wp,-6.37450174956073E-03_wp, &
      &-1.91461291146221E-02_wp, 1.09535528038150E-03_wp, 3.42404939671829E-02_wp, &
      & 4.74589429392484E-02_wp,-2.15316812457995E-02_wp,-4.80394095409307E-02_wp, &
      &-1.23060391666245E-02_wp, 1.29062067504150E-02_wp,-1.64629006596351E-03_wp, &
      &-3.56717280611371E-03_wp, 1.12323534619301E-03_wp, 1.19596404451922E-04_wp, &
      & 5.25101958893515E-03_wp, 1.59240225965459E-02_wp,-5.19297868799919E-03_wp, &
      &-2.32165734653952E-02_wp, 8.14988231623697E-03_wp,-1.77595705940923E-03_wp, &
      & 3.14448950646073E-03_wp, 3.43298532897802E-04_wp,-2.31816357724689E-03_wp, &
      &-3.22194483116891E-04_wp,-2.18477518470186E-02_wp, 9.14810407256850E-01_wp, &
      & 7.48033986208809E-03_wp, 5.09953559406022E-03_wp,-2.71827543681182E-02_wp, &
      &-1.35396965598019E-02_wp,-5.28749471869153E-03_wp,-4.55353025077479E-04_wp, &
      &-8.97158346346022E-03_wp,-1.38165545028928E-02_wp, 4.81975145073898E-03_wp, &
      & 1.02100793239369E-02_wp,-1.99268786228122E-02_wp, 3.15084854540211E-02_wp, &
      &-1.61959501492248E-02_wp, 8.60612099656805E-02_wp, 4.67381503377974E-02_wp, &
      &-7.10320337351804E-02_wp, 3.67817089890721E-02_wp,-2.56655009577627E-03_wp, &
      &-1.46067658017479E-02_wp, 5.62002262438172E-03_wp, 5.78076943435337E-02_wp, &
      &-2.44882129383441E-03_wp,-2.86025601155193E-02_wp,-5.04588128389799E-03_wp, &
      &-4.04081526462186E-03_wp, 1.53478173145747E-02_wp,-9.80190748894087E-03_wp, &
      &-1.88908493067724E-03_wp, 3.00321001380880E-04_wp, 1.84122040819062E-03_wp, &
      & 5.87522199748577E-06_wp, 1.29756344050563E-05_wp, 1.87604065635367E-02_wp, &
      &-3.84633749324155E-02_wp,-2.18459102876546E-02_wp,-1.12416390562587E-02_wp, &
      &-9.69766591163357E-04_wp, 3.64600111382412E-03_wp, 3.83817970980841E-06_wp, &
      &-6.97563672919024E-04_wp,-4.38986218965655E-03_wp, 3.31124366005085E-02_wp, &
      & 7.48033986208809E-03_wp, 9.07045088096433E-01_wp,-6.65538482288324E-03_wp, &
      &-4.69539950191276E-04_wp, 8.56203889791546E-03_wp,-1.56700021141413E-02_wp, &
      &-2.67566253419294E-02_wp,-5.14238323166957E-05_wp,-4.08742499801280E-02_wp, &
      & 1.15402494798059E-02_wp,-1.81555166175934E-02_wp,-2.34217202320808E-02_wp, &
      & 6.45053437127237E-02_wp, 4.22610487435226E-02_wp, 3.14099372304278E-03_wp, &
      &-1.03858523700383E-01_wp,-4.72446672147909E-02_wp, 5.10642880072443E-02_wp, &
      & 4.06754362288473E-02_wp,-4.95983018633469E-03_wp, 4.94693462764965E-02_wp, &
      &-1.68594875769994E-03_wp,-6.15913706401642E-02_wp,-2.74869037414869E-02_wp, &
      & 2.25098483451279E-03_wp, 2.91067062791004E-02_wp,-5.75999068005425E-03_wp, &
      &-5.04907733624699E-03_wp,-5.38795740743591E-04_wp,-6.77069550327594E-04_wp, &
      &-1.05685564533458E-03_wp,-2.03335470252935E-03_wp,-2.17978866795573E-03_wp, &
      & 4.69837569301324E-03_wp,-1.24613040121676E-02_wp,-3.26536302874007E-02_wp, &
      &-3.84377541807435E-03_wp, 1.79487344769638E-03_wp, 2.20063557970678E-03_wp, &
      & 2.16308617823632E-03_wp,-1.10536488762762E-03_wp,-8.07406407670555E-04_wp, &
      & 7.98537612496825E-02_wp, 5.09953559406022E-03_wp,-6.65538482288324E-03_wp, &
      & 8.85720560609369E-01_wp, 8.47405040928226E-03_wp, 4.32609508439594E-05_wp, &
      & 1.51332117219235E-02_wp,-1.32002681370673E-02_wp,-2.63946101633413E-02_wp, &
      &-7.18067691148896E-04_wp,-3.91472860045389E-04_wp,-3.91264029795358E-04_wp, &
      &-6.69813653944928E-04_wp,-9.38359904069108E-04_wp,-7.23150166636224E-04_wp, &
      & 3.78532516914335E-04_wp,-1.50740974288112E-03_wp,-2.04798561894034E-03_wp, &
      &-6.93371177426940E-04_wp,-4.24473085273045E-04_wp, 1.86494551430902E-04_wp, &
      & 1.30712345342775E-03_wp,-1.50326172865466E-04_wp,-1.62689399094805E-03_wp, &
      &-1.81297319011913E-03_wp, 2.05812072051300E-03_wp, 5.41082241160409E-03_wp, &
      & 6.00827861785124E-04_wp,-2.40085902556875E-03_wp, 5.99810932074870E-05_wp, &
      & 1.50498405466485E-04_wp,-7.41528292047566E-05_wp,-4.17552021346788E-05_wp, &
      &-2.95874145259961E-04_wp,-1.69113411385192E-03_wp, 1.38676485413577E-03_wp, &
      & 1.96920444748069E-03_wp,-1.75012611402543E-03_wp, 1.18227059103687E-04_wp, &
      &-1.22466148179130E-04_wp, 1.39871480504945E-05_wp, 1.02416470270356E-04_wp, &
      & 5.55356868900987E-06_wp, 1.97816816411445E-03_wp,-2.71827543681182E-02_wp, &
      &-4.69539950191276E-04_wp, 8.47405040928226E-03_wp, 9.24530237513602E-04_wp, &
      & 4.12508019867173E-04_wp, 3.21310649550766E-04_wp,-1.04518495185120E-04_wp, &
      & 1.07469239920213E-05_wp,-2.63641797653755E-04_wp,-1.85447400076914E-04_wp, &
      & 9.88365504780387E-06_wp,-4.38712775661282E-04_wp,-6.34043613155897E-04_wp, &
      &-6.96300015777415E-04_wp, 9.38872653962537E-04_wp, 3.82651994863863E-05_wp, &
      &-1.49699282601156E-03_wp,-1.52128659470246E-04_wp,-4.73872134127807E-04_wp, &
      &-2.30162533491548E-05_wp, 4.10116747608010E-04_wp, 5.22664015462132E-04_wp, &
      &-5.53268037298782E-04_wp,-9.90759348859760E-04_wp, 5.05766248283278E-04_wp, &
      & 2.05070858833908E-03_wp, 6.76325363903056E-04_wp,-1.16119924081091E-03_wp, &
      & 1.05424830110441E-05_wp, 8.35670166748950E-05_wp,-5.64959811538514E-06_wp, &
      &-5.09405249893958E-06_wp,-1.25680696991265E-04_wp,-6.61974620886243E-04_wp, &
      &-6.88153886369985E-04_wp, 5.26310641005283E-04_wp,-9.81189201589842E-04_wp, &
      & 2.23628169815579E-05_wp,-1.14817581188053E-05_wp,-1.14181893240403E-05_wp, &
      & 4.30064827780183E-05_wp,-6.21201449865174E-05_wp, 7.95287884239188E-04_wp, &
      &-1.35396965598019E-02_wp, 8.56203889791546E-03_wp, 4.32609508439594E-05_wp, &
      & 4.12508019867173E-04_wp, 2.88003393911525E-04_wp,-6.50182407583915E-05_wp, &
      &-2.53634831707603E-04_wp, 1.27233371966009E-04_wp,-5.93966433083412E-04_wp, &
      & 2.69372961745180E-05_wp,-5.75792938733272E-04_wp,-1.20006019568062E-04_wp, &
      & 3.14862311605270E-04_wp, 6.92000783727072E-04_wp,-1.19990637124712E-03_wp, &
      &-2.40313487407486E-03_wp, 1.86593826628326E-04_wp,-6.60971351200672E-05_wp, &
      & 6.75497357338132E-04_wp, 1.12528083885465E-04_wp, 9.58776573157464E-04_wp, &
      &-1.09131629766287E-03_wp,-1.20341681703358E-03_wp,-3.97806973349798E-04_wp, &
      & 7.91454294590990E-04_wp, 2.29339909979139E-03_wp,-1.66014282991114E-03_wp, &
      &-5.81093006975118E-04_wp, 4.34016120608557E-05_wp,-1.78507090937897E-05_wp, &
      &-9.38849407025514E-05_wp,-6.69138499163433E-05_wp,-1.08710889283059E-04_wp, &
      &-1.25193134896855E-03_wp, 2.13804102298870E-03_wp,-1.60834820233167E-05_wp, &
      &-3.98266659773968E-04_wp, 9.59540695118697E-05_wp,-6.53524280632076E-05_wp, &
      & 6.32963494893182E-05_wp, 1.06480846519756E-05_wp, 9.70425331312997E-05_wp, &
      & 1.47614735664917E-03_wp,-5.28749471869153E-03_wp,-1.56700021141413E-02_wp, &
      & 1.51332117219235E-02_wp, 3.21310649550766E-04_wp,-6.50182407583915E-05_wp, &
      & 5.65537994689618E-04_wp, 2.39816673197981E-04_wp,-4.01145908943552E-04_wp, &
      & 1.02979125404226E-03_wp,-3.83960703078053E-04_wp,-4.71112289725868E-05_wp, &
      & 1.02053886193372E-03_wp,-1.37787803489033E-03_wp, 2.14568368478460E-04_wp, &
      &-2.29785457537160E-03_wp, 3.71631530997041E-04_wp, 2.32423304833032E-03_wp, &
      &-2.04580162271832E-03_wp,-5.61485634191517E-04_wp, 5.83275218807750E-04_wp, &
      &-7.52025031397998E-04_wp,-1.82188558275468E-03_wp, 1.03007344334728E-03_wp, &
      & 1.22820555419493E-03_wp, 7.16816433450695E-04_wp, 4.07757821937570E-04_wp, &
      &-2.79625792411043E-04_wp, 2.45027490200433E-03_wp, 1.16049685662805E-04_wp, &
      & 5.03854201101717E-06_wp,-6.44076461307671E-05_wp, 5.12364190514789E-05_wp, &
      & 5.08755648891243E-05_wp,-1.80059796257016E-03_wp, 3.54546470887379E-03_wp, &
      & 3.65950247565885E-03_wp, 2.71148530315096E-03_wp,-6.27125251361419E-06_wp, &
      &-2.32450086452956E-04_wp,-6.58156831460961E-05_wp, 5.28437067393878E-05_wp, &
      & 2.32661759562063E-04_wp,-3.09112108744030E-03_wp,-4.55353025077479E-04_wp, &
      &-2.67566253419294E-02_wp,-1.32002681370673E-02_wp,-1.04518495185120E-04_wp, &
      &-2.53634831707603E-04_wp, 2.39816673197981E-04_wp, 1.01972520181179E-03_wp, &
      & 4.18800952996620E-04_wp, 1.20885516832012E-03_wp,-5.52056921452567E-04_wp, &
      & 5.46257757174947E-04_wp, 5.94007597988317E-04_wp,-2.24225399100554E-03_wp, &
      &-1.19565402301603E-03_wp, 4.74724238071550E-05_wp, 2.60157115829100E-03_wp, &
      & 4.22160666293566E-04_wp,-1.86857557022780E-03_wp,-1.60968883271809E-03_wp, &
      & 4.19280315359991E-04_wp,-1.22301380113751E-03_wp,-3.25003079162081E-05_wp, &
      & 1.53200056222919E-03_wp, 4.56941288207337E-04_wp, 1.18902455312510E-03_wp, &
      &-4.27488088446933E-04_wp, 2.06014446932550E-03_wp, 1.39966952488796E-03_wp, &
      & 5.66604383383437E-05_wp, 9.88365521261470E-05_wp, 4.02013059080886E-05_wp, &
      & 1.03787374695025E-04_wp, 3.70646570213923E-05_wp,-3.51490801006542E-04_wp, &
      &-6.75694770587501E-06_wp, 4.03502323498763E-03_wp, 1.45446242844250E-03_wp, &
      &-7.33810769916618E-05_wp,-1.56880573246973E-04_wp,-1.36918692981879E-04_wp, &
      & 9.07005439053863E-05_wp, 3.40426800813876E-05_wp,-3.26315136987568E-03_wp, &
      &-8.97158346346022E-03_wp,-5.14238323166957E-05_wp,-2.63946101633413E-02_wp, &
      & 1.07469239920213E-05_wp, 1.27233371966009E-04_wp,-4.01145908943552E-04_wp, &
      & 4.18800952996620E-04_wp, 8.91699965008989E-04_wp], shape(density))

   call get_structure(mol, "f-block", "CeCl3")
   call test_acp_numgrad(error, mol, density, make_qvszp_basis, thr_in=thr1*10)

end subroutine test_g_acp_gxtb_cecl3

end module test_acp