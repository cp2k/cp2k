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

module test_hamiltonian
   use mctc_data_covrad, only : get_covalent_rad
   use mctc_data_paulingen, only : get_pauling_en
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use mstore, only : get_structure
   use multicharge, only : get_eeqbc_charges
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_cache, only : basis_cache
   use tblite_basis_qvszp, only : qvszp_basis_type, qvszp_cgto_type, &
      & new_qvszp_cgto, new_qvszp_basis
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_basis_type, only : basis_type, new_basis, cgto_container, &
      & cgto_type, get_cutoff
   use tblite_blas, only : gemv
   use tblite_ceh_ceh, only : ceh_h0spec
   use tblite_cutoff, only : get_lattice_points
   use tblite_scf_iterator, only : get_electronic_energy
   use tblite_scf_potential, only : potential_type, new_potential
   use tblite_wavefunction_spin, only : updown_to_magnet
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_gfn2, only : gfn2_h0spec
   use tblite_xtb_gxtb, only : gxtb_h0spec
   use tblite_xtb_h0, only : new_hamiltonian, tb_hamiltonian, &
      & get_hamiltonian, get_hamiltonian_gradient, get_selfenergy, &
      & get_anisotropy, get_anisotropy_gradient
   use tblite_xtb_spec, only : tb_h0spec
   implicit none
   private

   public :: collect_hamiltonian

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1e5*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine basis_maker(bas, mol, error, ng, scale_h0_basis, ncoord_bas)
         import :: basis_type, structure_type, error_type, ncoord_type
         class(basis_type), allocatable, intent(out) :: bas
         type(structure_type), intent(in) :: mol
         type(error_type), allocatable, intent(out) :: error
         integer, intent(in), optional :: ng
         logical, intent(in), optional :: scale_h0_basis
         class(ncoord_type), intent(out), allocatable, optional :: ncoord_bas
      end subroutine basis_maker

      subroutine ncoord_maker(ncoord, ncoord_en, mol, cn_cutoff, error)
         import :: ncoord_type, structure_type, error_type, wp
         class(ncoord_type), allocatable, intent(out) :: ncoord
         class(ncoord_type), allocatable, intent(out) :: ncoord_en
         type(structure_type), intent(in) :: mol
         real(wp), intent(in) :: cn_cutoff
         type(error_type), allocatable, intent(out) :: error
      end subroutine ncoord_maker
   end interface

contains


!> Collect all exported unit tests
subroutine collect_hamiltonian(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("selfenergy-gfn2-h2", test_selfenergy_gfn2_h2), &
      new_unittest("selfenergy-gfn2-lih", test_selfenergy_gfn2_lih), &
      new_unittest("selfenergy-gfn2-s2", test_selfenergy_gfn2_s2), &
      new_unittest("selfenergy-gfn2-sih4", test_selfenergy_gfn2_sih4), &
      new_unittest("selfenergy-gxtb-h2", test_selfenergy_gxtb_h2), &
      new_unittest("selfenergy-gxtb-lih", test_selfenergy_gxtb_lih), &
      new_unittest("selfenergy-gxtb-s2", test_selfenergy_gxtb_s2), &
      new_unittest("selfenergy-gxtb-sih4", test_selfenergy_gxtb_sih4), &
      new_unittest("selfenergy-gxtb-cecl3", test_selfenergy_gxtb_cecl3), &
      new_unittest("selfenergy-gxtb-ce2", test_selfenergy_gxtb_ce2), &
      new_unittest("selfenergy-ceh-h2", test_selfenergy_ceh_h2), &
      new_unittest("selfenergy-ceh-lih", test_selfenergy_ceh_lih), &
      new_unittest("selfenergy-ceh-s2", test_selfenergy_ceh_s2), &
      new_unittest("selfenergy-ceh-sih4", test_selfenergy_ceh_sih4), &
      new_unittest("selfenergy-ceh-accl6", test_selfenergy_ceh_accl6), &
      new_unittest("anisotropy-gxtb-h2", test_anisotropy_gxtb_h2), &
      new_unittest("anisotropy-gxtb-lih", test_anisotropy_gxtb_lih), &
      new_unittest("anisotropy-gxtb-sih4", test_anisotropy_gxtb_sih4), &
      new_unittest("anisotropy-gxtb-cecl3", test_anisotropy_gxtb_cecl3), &
      new_unittest("anisotropy-gxtb-m01", test_anisotropy_gxtb_m01), &
      new_unittest("anisotropy-gxtb-m02", test_anisotropy_gxtb_m02), &
      new_unittest("hamiltonian-gfn2-h2", test_hamiltonian_gfn2_h2), &
      new_unittest("hamiltonian-gfn2-lih", test_hamiltonian_gfn2_lih), &
      new_unittest("hamiltonian-gfn2-s2", test_hamiltonian_gfn2_s2), &
      new_unittest("hamiltonian-gfn2-sih4", test_hamiltonian_gfn2_sih4), &
      new_unittest("hamiltonian-gxtb-h2", test_hamiltonian_gxtb_h2), &
      new_unittest("hamiltonian-gxtb-lih", test_hamiltonian_gxtb_lih), &
      new_unittest("hamiltonian-gxtb-s2", test_hamiltonian_gxtb_s2), &
      new_unittest("hamiltonian-gxtb-sih4", test_hamiltonian_gxtb_sih4), &
      new_unittest("hamiltonian-gxtb-cecl3", test_hamiltonian_gxtb_cecl3), &
      new_unittest("hamiltonian-gxtb-ce2", test_hamiltonian_gxtb_ce2), &
      new_unittest("hamiltonian-ceh-h2", test_hamiltonian_ceh_h2), &
      new_unittest("hamiltonian-ceh-lih", test_hamiltonian_ceh_lih), &
      new_unittest("hamiltonian-ceh-s2", test_hamiltonian_ceh_s2), &
      new_unittest("hamiltonian-ceh-sih4", test_hamiltonian_ceh_sih4), &
      new_unittest("hamiltonian-ceh-panp", test_hamiltonian_ceh_panp), &
      new_unittest("anisotropy-gradient-gxtb-h2", test_anisotropy_gradient_gxtb_h2), &
      new_unittest("anisotropy-gradient-gxtb-lih", test_anisotropy_gradient_gxtb_lih), &
      new_unittest("anisotropy-gradient-gxtb-sih4", test_anisotropy_gradient_gxtb_sih4), &
      new_unittest("anisotropy-gradient-gxtb-cecl3", test_anisotropy_gradient_gxtb_cecl3), &
      new_unittest("anisotropy-gradient-gxtb-m01", test_anisotropy_gradient_gxtb_m01), &
      new_unittest("anisotropy-gradient-gxtb-m02", test_anisotropy_gradient_gxtb_m02), &
      new_unittest("energy-aniso-grad-gxtb-h2", test_g_energy_aniso_gxtb_h2), &
      new_unittest("energy-aniso-grad-gxtb-lih", test_g_energy_aniso_gxtb_lih), &
      new_unittest("energy-aniso-grad-gxtb-no", test_g_energy_aniso_gxtb_no), &
      new_unittest("energy-aniso-grad-gxtb-s2", test_g_energy_aniso_gxtb_s2), &
      new_unittest("energy-aniso-grad-gxtb-sih4", test_g_energy_aniso_gxtb_sih4), &
      new_unittest("energy-aniso-grad-gxtb-cecl3", test_g_energy_aniso_gxtb_cecl3), &
      new_unittest("energy-qeff-grad-gxtb-h2", test_qg_hamiltonian_gxtb_h2), &
      new_unittest("energy-qeff-grad-gxtb-lih", test_qg_hamiltonian_gxtb_lih), &
      new_unittest("energy-qeff-grad-gxtb-no", test_qg_hamiltonian_gxtb_no), &
      new_unittest("energy-qeff-grad-gxtb-s2", test_qg_hamiltonian_gxtb_s2), &
      new_unittest("energy-qeff-grad-gxtb-sih4", test_qg_hamiltonian_gxtb_sih4), &
      new_unittest("energy-qeff-grad-gxtb-cecl3", test_qg_hamiltonian_gxtb_cecl3), &
      new_unittest("hamiltonian-gradient-gfn2-h2", test_g_hamiltonian_gfn2_h2), &
      new_unittest("hamiltonian-gradient-gfn2-lih", test_g_hamiltonian_gfn2_lih), &
      new_unittest("hamiltonian-gradient-gfn2-s2", test_g_hamiltonian_gfn2_s2), &
      new_unittest("hamiltonian-gradient-gfn2-no", test_g_hamiltonian_gfn2_no), &
      new_unittest("hamiltonian-gradient-gfn2-sih4", test_g_hamiltonian_gfn2_sih4), &
      new_unittest("hamiltonian-gradient-gxtb-h2", test_g_hamiltonian_gxtb_h2), &
      new_unittest("hamiltonian-gradient-gxtb-lih", test_g_hamiltonian_gxtb_lih), &
      new_unittest("hamiltonian-gradient-gxtb-no", test_g_hamiltonian_gxtb_no), &
      new_unittest("hamiltonian-gradient-gxtb-s2", test_g_hamiltonian_gxtb_s2), &
      new_unittest("hamiltonian-gradient-gxtb-sih4", test_g_hamiltonian_gxtb_sih4), &
      new_unittest("hamiltonian-gradient-gxtb-cecl3", test_g_hamiltonian_gxtb_cecl3) &
      ]

end subroutine collect_hamiltonian


subroutine make_gfn2_basis(bas, mol, error, ng, scale_h0_basis, ncoord_bas)
   !> Basis set information
   class(basis_type), allocatable, intent(out) :: bas
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Number of Gaussians per Slater function
   integer, intent(in), optional :: ng
   !> Scale basis set for H0 (unused)
   logical, intent(in), optional :: scale_h0_basis
   !> Optional coordination number used in the basis set (unused)
   class(ncoord_type), intent(out), allocatable, optional :: ncoord_bas
   
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

   integer :: isp, izp, ish, stat, local_ng
   integer, allocatable :: nshell(:)
   type(cgto_container), allocatable :: cgto(:, :)

   if (present(ng)) then
      local_ng = ng
   else
      local_ng = 6
   end if

   nshell = nsh(mol%num)
   allocate(cgto(maxval(nshell), mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, nshell(isp)
         allocate(cgto_type :: cgto(ish, isp)%raw)
         call slater_to_gauss(local_ng, pqn(ish, izp), lsh(ish, izp), &
            & zeta(ish, izp), cgto(ish, isp)%raw, .true., stat)
      end do
   end do

   allocate(bas)
   call new_basis(bas, mol, nshell, cgto, accuracy=1.0_wp)

end subroutine make_gfn2_basis


subroutine make_qvszp_basis(bas, mol, error, ng, scale_h0_basis, ncoord_bas)
   !> Basis set information
   class(basis_type), allocatable, intent(out) :: bas
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Number of Gaussians per Slater function (unused)
   integer, intent(in), optional :: ng
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


subroutine make_ceh_basis(bas, mol, error, ng, scale_h0_basis, ncoord_bas)
   !> Basis set information
   class(basis_type), allocatable, intent(out) :: bas
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Number of Gaussians per Slater function
   integer, intent(in), optional :: ng
   !> Scale basis set for H0
   logical, intent(in), optional :: scale_h0_basis
   !> Optional coordination number used in the basis set (unused)
   class(ncoord_type), intent(out), allocatable, optional :: ncoord_bas

   integer, parameter :: nsh(103) = [&
   & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, & ! 1-20
   & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, & ! 21-40
   & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, & ! 41-60
   & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, & ! 61-80
   & 3, 3, 3, 3, 3, 3, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, & ! 81-100
   & 4, 4, 4]

   integer, parameter :: lsh(4, 103) = reshape([&
   & 0, 0, 0, 0,  0, 0, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0, & ! 1-6
   & 0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 2, 0, & ! 7-12
   & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 13-18
   & 0, 1, 0, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 19-24
   & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 0, 0, & ! 25-30
   & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 31-36
   & 0, 1, 0, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 37-42
   & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 0, 0, & ! 43-48
   & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 49-54
   & 0, 1, 0, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 55-60
   & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 61-66
   & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 67-72
   & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 73-78
   & 0, 1, 2, 0,  0, 1, 0, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 2, 0, & ! 79-84
   & 0, 1, 2, 0,  0, 1, 2, 0,  0, 1, 0, 0,  0, 1, 2, 0,  0, 1, 2, 3,  0, 1, 2, 3, & ! 85-90
   & 0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3, & ! 91-96
   & 0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3, & ! 97-102
   & 0, 1, 2, 3], shape(lsh))

   integer, parameter :: pqn(4, 103) = reshape([&
   & 1, 0, 0, 0,  1, 0, 0, 0,  2, 2, 0, 0,  2, 2, 0, 0,  2, 2, 0, 0,  2, 2, 0, 0, & ! 1-6
   & 2, 2, 0, 0,  2, 2, 0, 0,  2, 2, 0, 0,  2, 2, 0, 0,  3, 3, 0, 0,  3, 3, 3, 0, & ! 7-12
   & 3, 3, 3, 0,  3, 3, 3, 0,  3, 3, 3, 0,  3, 3, 3, 0,  3, 3, 3, 0,  3, 3, 3, 0, & ! 13-18
   & 4, 4, 0, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0, & ! 19-24
   & 4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 3, 0,  4, 4, 0, 0, & ! 25-30
   & 4, 4, 4, 0,  4, 4, 4, 0,  4, 4, 4, 0,  4, 4, 4, 0,  4, 4, 4, 0,  4, 4, 4, 0, & ! 31-36
   & 5, 5, 0, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0, & ! 37-42
   & 5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 4, 0,  5, 5, 0, 0, & ! 43-48
   & 5, 5, 5, 0,  5, 5, 5, 0,  5, 5, 5, 0,  5, 5, 5, 0,  5, 5, 5, 0,  5, 5, 5, 0, & ! 49-54
   & 6, 6, 0, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0, & ! 55-60
   & 6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0, & ! 61-66
   & 6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0, & ! 67-72
   & 6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0, & ! 73-78
   & 6, 6, 5, 0,  6, 6, 0, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 5, 0, & ! 79-84
   & 6, 6, 5, 0,  6, 6, 5, 0,  6, 6, 0, 0,  6, 6, 5, 0,  6, 6, 5, 5,  6, 6, 5, 5, & ! 85-90
   & 6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5, & ! 91-96
   & 6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5,  6, 6, 5, 5, & ! 97-102
   & 6, 6, 5, 5], shape(pqn))

   real(wp), parameter :: zeta(4, 103) = reshape([&
   & 1.23363166_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, & ! 1
   & 2.27004605_wp, 0.00000000_wp, 0.00000000_wp, 0.00000000_wp, & ! 2
   & 0.86185456_wp, 1.42017184_wp, 0.00000000_wp, 0.00000000_wp, & ! 3
   & 1.76817995_wp, 1.44095844_wp, 0.00000000_wp, 0.00000000_wp, & ! 4
   & 2.06339837_wp, 1.52051807_wp, 0.00000000_wp, 0.00000000_wp, & ! 5
   & 2.56058582_wp, 1.86484737_wp, 0.00000000_wp, 0.00000000_wp, & ! 6
   & 2.71233631_wp, 2.19848968_wp, 0.00000000_wp, 0.00000000_wp, & ! 7
   & 3.21585650_wp, 2.41309737_wp, 0.00000000_wp, 0.00000000_wp, & ! 8
   & 3.82146807_wp, 2.63063636_wp, 0.00000000_wp, 0.00000000_wp, & ! 9
   & 4.62721228_wp, 2.53599954_wp, 0.00000000_wp, 0.00000000_wp, & ! 10
   & 0.93221172_wp, 1.55333839_wp, 0.00000000_wp, 0.00000000_wp, & ! 11
   & 1.77220557_wp, 1.59942632_wp, 2.98596647_wp, 0.00000000_wp, & ! 12
   & 2.26040231_wp, 1.78718151_wp, 2.00990188_wp, 0.00000000_wp, & ! 13
   & 1.85259089_wp, 1.81733349_wp, 1.65269988_wp, 0.00000000_wp, & ! 14
   & 2.65701241_wp, 2.03189759_wp, 2.03883661_wp, 0.00000000_wp, & ! 15
   & 2.60609998_wp, 2.16530440_wp, 2.41888232_wp, 0.00000000_wp, & ! 16
   & 2.78818934_wp, 2.24732894_wp, 1.99081182_wp, 0.00000000_wp, & ! 17
   & 2.55424399_wp, 2.20946190_wp, 1.93619550_wp, 0.00000000_wp, & ! 18
   & 1.73713827_wp, 1.33788617_wp, 0.00000000_wp, 0.00000000_wp, & ! 19
   & 2.47982574_wp, 1.07250770_wp, 2.11920764_wp, 0.00000000_wp, & ! 20
   & 2.22449249_wp, 1.55418319_wp, 2.00953578_wp, 0.00000000_wp, & ! 21
   & 2.58879616_wp, 0.99441077_wp, 1.88561781_wp, 0.00000000_wp, & ! 22
   & 3.04370654_wp, 4.03007600_wp, 1.66329169_wp, 0.00000000_wp, & ! 23
   & 2.25012727_wp, 2.70681556_wp, 1.67501904_wp, 0.00000000_wp, & ! 24
   & 2.20605319_wp, 2.82019792_wp, 1.86102254_wp, 0.00000000_wp, & ! 25
   & 1.57297015_wp, 1.98621494_wp, 2.83790684_wp, 0.00000000_wp, & ! 26
   & 1.80826602_wp, 1.73675835_wp, 2.79767448_wp, 0.00000000_wp, & ! 27
   & 2.00758945_wp, 2.25075692_wp, 2.98291663_wp, 0.00000000_wp, & ! 28
   & 2.18159986_wp, 2.38459096_wp, 3.09502522_wp, 0.00000000_wp, & ! 29
   & 2.26376756_wp, 2.20362977_wp, 0.00000000_wp, 0.00000000_wp, & ! 30
   & 2.63822153_wp, 2.06752328_wp, 2.11361643_wp, 0.00000000_wp, & ! 31
   & 2.52891955_wp, 2.19441794_wp, 1.77661998_wp, 0.00000000_wp, & ! 32
   & 3.55667605_wp, 2.42075463_wp, 1.46579772_wp, 0.00000000_wp, & ! 33
   & 2.89652631_wp, 2.45421858_wp, 2.27883625_wp, 0.00000000_wp, & ! 34
   & 3.28921099_wp, 2.56526915_wp, 1.64501640_wp, 0.00000000_wp, & ! 35
   & 5.20988189_wp, 2.84336725_wp, 2.75838814_wp, 0.00000000_wp, & ! 36
   & 1.26972917_wp, 1.88730596_wp, 0.00000000_wp, 0.00000000_wp, & ! 37
   & 1.86880714_wp, 1.78546342_wp, 2.16012236_wp, 0.00000000_wp, & ! 38
   & 0.92001877_wp, 1.45732462_wp, 2.22901395_wp, 0.00000000_wp, & ! 39
   & 6.50647305_wp, 1.43202338_wp, 2.11971490_wp, 0.00000000_wp, & ! 40
   & 2.10973371_wp, 2.79944781_wp, 2.01897369_wp, 0.00000000_wp, & ! 41
   & 2.58413333_wp, 3.02795359_wp, 2.08733665_wp, 0.00000000_wp, & ! 42
   & 2.62141555_wp, 3.13487625_wp, 2.13259872_wp, 0.00000000_wp, & ! 43
   & 2.73984499_wp, 2.18167834_wp, 2.54609647_wp, 0.00000000_wp, & ! 44
   & 1.84057176_wp, 2.97482636_wp, 3.10693700_wp, 0.00000000_wp, & ! 45
   & 1.75622839_wp, 3.39424756_wp, 3.20265306_wp, 0.00000000_wp, & ! 46
   & 3.05018811_wp, 2.34951987_wp, 3.35332952_wp, 0.00000000_wp, & ! 47
   & 2.41999128_wp, 2.28892954_wp, 0.00000000_wp, 0.00000000_wp, & ! 48
   & 2.87813961_wp, 2.44659724_wp, 2.75773502_wp, 0.00000000_wp, & ! 49
   & 3.03823214_wp, 2.32082155_wp, 1.77513328_wp, 0.00000000_wp, & ! 50
   & 2.68750711_wp, 2.38565373_wp, 2.12596190_wp, 0.00000000_wp, & ! 51
   & 2.81071790_wp, 2.45274786_wp, 2.01871821_wp, 0.00000000_wp, & ! 52
   & 2.90686956_wp, 2.49377102_wp, 1.90073732_wp, 0.00000000_wp, & ! 53
   & 4.17531340_wp, 2.86937955_wp, 2.96894812_wp, 0.00000000_wp, & ! 54
   & 1.24299361_wp, 1.99142040_wp, 0.00000000_wp, 0.00000000_wp, & ! 55
   & 1.31400366_wp, 1.16438481_wp, 2.12759606_wp, 0.00000000_wp, & ! 56
   & 2.81737350_wp, 1.69863323_wp, 2.27369715_wp, 0.00000000_wp, & ! 57
   & 2.84503901_wp, 1.46018192_wp, 2.53498936_wp, 0.00000000_wp, & ! 58
   & 2.81697107_wp, 1.47545307_wp, 2.54350275_wp, 0.00000000_wp, & ! 59
   & 2.78890313_wp, 1.49072422_wp, 2.55201615_wp, 0.00000000_wp, & ! 60
   & 2.76083520_wp, 1.50599537_wp, 2.56052955_wp, 0.00000000_wp, & ! 61
   & 2.73276726_wp, 1.52126652_wp, 2.56904294_wp, 0.00000000_wp, & ! 62
   & 2.70469932_wp, 1.53653767_wp, 2.57755634_wp, 0.00000000_wp, & ! 63
   & 2.67663138_wp, 1.55180881_wp, 2.58606974_wp, 0.00000000_wp, & ! 64
   & 2.64856345_wp, 1.56707996_wp, 2.59458313_wp, 0.00000000_wp, & ! 65
   & 2.62049551_wp, 1.58235111_wp, 2.60309653_wp, 0.00000000_wp, & ! 66
   & 2.59242757_wp, 1.59762226_wp, 2.61160992_wp, 0.00000000_wp, & ! 67
   & 2.56435964_wp, 1.61289341_wp, 2.62012332_wp, 0.00000000_wp, & ! 68
   & 2.53629170_wp, 1.62816456_wp, 2.62863672_wp, 0.00000000_wp, & ! 69
   & 2.50822376_wp, 1.64343571_wp, 2.63715011_wp, 0.00000000_wp, & ! 70
   & 2.48015583_wp, 1.65870685_wp, 2.64566351_wp, 0.00000000_wp, & ! 71
   & 3.19537752_wp, 2.24853837_wp, 2.41492177_wp, 0.00000000_wp, & ! 72
   & 3.14122020_wp, 2.48723489_wp, 2.21933576_wp, 0.00000000_wp, & ! 73
   & 3.17661283_wp, 3.39538568_wp, 2.37502789_wp, 0.00000000_wp, & ! 74
   & 3.14538352_wp, 2.58361113_wp, 2.47139347_wp, 0.00000000_wp, & ! 75
   & 1.81565647_wp, 2.48106221_wp, 3.18585355_wp, 0.00000000_wp, & ! 76
   & 2.11798490_wp, 2.85857032_wp, 3.47048400_wp, 0.00000000_wp, & ! 77
   & 2.71241232_wp, 3.37886078_wp, 3.64124964_wp, 0.00000000_wp, & ! 78
   & 2.80572458_wp, 2.82570220_wp, 3.72064445_wp, 0.00000000_wp, & ! 79
   & 2.61951362_wp, 2.69607886_wp, 0.00000000_wp, 0.00000000_wp, & ! 80
   & 3.05383193_wp, 2.61683803_wp, 3.32179612_wp, 0.00000000_wp, & ! 81
   & 3.02135073_wp, 2.59250246_wp, 4.24674489_wp, 0.00000000_wp, & ! 82
   & 3.16405210_wp, 2.63238785_wp, 3.04625573_wp, 0.00000000_wp, & ! 83
   & 2.96133467_wp, 2.71388453_wp, 2.31022562_wp, 0.00000000_wp, & ! 84
   & 2.98240599_wp, 2.95960758_wp, 2.43778345_wp, 0.00000000_wp, & ! 85
   & 3.07936232_wp, 2.68589775_wp, 2.10311395_wp, 0.00000000_wp, & ! 86
   & 1.81913220_wp, 3.23064408_wp, 0.00000000_wp, 0.00000000_wp, & ! 87
   & 2.43263729_wp, 2.47485608_wp, 2.09113715_wp, 0.00000000_wp, & ! 88
   & 3.65108887_wp, 3.45440279_wp, 1.97314608_wp, 1.98901892_wp, & ! 89
   & 3.35816295_wp, 2.81245896_wp, 2.05947820_wp, 2.04247660_wp, & ! 90
   & 3.08262439_wp, 2.24936413_wp, 2.13696560_wp, 2.09705269_wp, & ! 91
   & 2.82447317_wp, 1.76511830_wp, 2.20560830_wp, 2.15274719_wp, & ! 92
   & 2.58370931_wp, 1.35972146_wp, 2.26540630_wp, 2.20956009_wp, & ! 93
   & 2.36033280_wp, 1.03317362_wp, 2.31635959_wp, 2.26749140_wp, & ! 94
   & 2.15434364_wp, 0.78547477_wp, 2.35846817_wp, 2.32654112_wp, & ! 95
   & 1.96574183_wp, 0.61662492_wp, 2.39173204_wp, 2.38670925_wp, & ! 96
   & 1.79452738_wp, 0.52662407_wp, 2.41615121_wp, 2.44799578_wp, & ! 97
   & 1.64070027_wp, 0.51547221_wp, 2.43172568_wp, 2.51040072_wp, & ! 98
   & 1.50426052_wp, 0.58316935_wp, 2.43845543_wp, 2.57392407_wp, & ! 99
   & 1.38520812_wp, 0.72971549_wp, 2.43634048_wp, 2.63856583_wp, & ! 100
   & 1.28354307_wp, 0.95511062_wp, 2.42538083_wp, 2.70432599_wp, & ! 101
   & 1.19926537_wp, 1.25935475_wp, 2.40557646_wp, 2.77120456_wp, & ! 102
   & 1.13237503_wp, 1.64244787_wp, 2.37692740_wp, 2.83920154_wp],& ! 103
   & shape(zeta))

   integer :: isp, izp, ish, stat, local_ng
   integer, allocatable :: nshell(:)
   type(cgto_container), allocatable :: cgto(:, :)

   if (present(ng)) then
      local_ng = ng
   else
      local_ng = 6
   end if

   nshell = nsh(mol%num)
   allocate(cgto(maxval(nshell), mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      do ish = 1, nshell(isp)
         allocate(cgto_type :: cgto(ish, isp)%raw)
         call slater_to_gauss(local_ng, pqn(ish, izp), lsh(ish, izp), &
            & zeta(ish, izp), cgto(ish, isp)%raw, .true., stat)
      end do
   end do

   allocate(bas)
   call new_basis(bas, mol, nshell, cgto, accuracy=1.0_wp)

end subroutine make_ceh_basis


subroutine make_gfn2_ncoord(ncoord, ncoord_en, mol, cn_cutoff, error)
   !> Coordination number
   class(ncoord_type), allocatable, intent(out) :: ncoord
   !> Electronegativity scaled coordination number
   class(ncoord_type), allocatable, intent(out) :: ncoord_en
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Real space cutoff
   real(wp), intent(in) :: cn_cutoff
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), allocatable :: rcov(:)

   allocate(rcov(mol%nid))
   rcov(:) = get_covalent_rad(mol%num)
   call new_ncoord(ncoord, mol, cn_count%dexp, error, cutoff=cn_cutoff, rcov=rcov)

end subroutine make_gfn2_ncoord


subroutine make_gxtb_ncoord(ncoord, ncoord_en, mol, cn_cutoff, error)
   !> Coordination number
   class(ncoord_type), allocatable, intent(out) :: ncoord
   !> Electronegativity scaled coordination number
   class(ncoord_type), allocatable, intent(out) :: ncoord_en
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Real space cutoff
   real(wp), intent(in) :: cn_cutoff
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

   !> Coordination number exponent
   real(wp), parameter :: p_kcn = 2.0680000000_wp

   call new_ncoord(ncoord, mol, cn_count_type=cn_count%erf, error=error, &
      & cutoff=cn_cutoff, kcn=p_kcn, rcov=p_cn_rcov(mol%num))

end subroutine make_gxtb_ncoord


subroutine make_ceh_ncoord(ncoord, ncoord_en, mol, cn_cutoff, error)
   !> Coordination number
   class(ncoord_type), allocatable, intent(out) :: ncoord
   !> Electronegativity scaled coordination number
   class(ncoord_type), allocatable, intent(out) :: ncoord_en
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Real space cutoff
   real(wp), intent(in) :: cn_cutoff
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), allocatable :: rcov(:), en(:)

   allocate(rcov(mol%nid), en(mol%nid))
   rcov(:) = get_covalent_rad(mol%num)
   en(:) = get_pauling_en(mol%num)
   call new_ncoord(ncoord, mol, cn_count%erf, error, &
      & cutoff=cn_cutoff, rcov=rcov)
   if (allocated(error)) return
   call new_ncoord(ncoord_en, mol, cn_count%erf_en, error, &
      & cutoff=cn_cutoff, rcov=rcov)
   if (allocated(error)) return

end subroutine make_ceh_ncoord


subroutine test_selfenergy_gen(error, mol, make_basis, spec, make_ncoord, ref, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Factory to create new basis set objects
   procedure(basis_maker) :: make_basis

   !> Hamiltonian specification
   class(tb_h0spec), intent(in) :: spec

   !> Factory to create new coordination number objects
   procedure(ncoord_maker) :: make_ncoord

   !> Reference self-energies
   real(wp), intent(in) :: ref(:)

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   class(basis_type), allocatable :: bas
   type(tb_hamiltonian) :: h0
   class(ncoord_type), allocatable :: ncoord, ncoord_en
   real(wp), parameter :: cn_cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :), cn(:), cn_en(:)
   real(wp), allocatable :: selfenergy(:)
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Setup basis set
   call make_basis(bas, mol, error=error, ng=6, scale_h0_basis=.true.)
   if (allocated(error)) return

   ! Setup Hamiltonian
   call new_hamiltonian(h0, mol, bas, spec)

   ! Obtain coordination number
   call make_ncoord(ncoord, ncoord_en, mol, cn_cutoff, error)
   if (allocated(error)) return
   call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
   if (allocated(ncoord)) then
      allocate(cn(mol%nat))
      call ncoord%get_coordination_number(mol, lattr, cn)
   end if
   if (allocated(ncoord_en)) then
      allocate(cn_en(mol%nat))
      call ncoord_en%get_coordination_number(mol, lattr, cn_en)
   end if

   ! Calculate self-energies
   allocate(selfenergy(bas%nsh))
   call get_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
      & selfenergy=selfenergy)

   ! where(abs(selfenergy) < thr) selfenergy = 0.0_wp
   ! print '(*(6x,"&", 3(es20.14e1, "_wp":, ","), "&", /))', selfenergy

   if (any(abs(selfenergy - ref) > thr_)) then
      call test_failed(error, "Self energy does not match.")
      write(*,*) 'Reference:'
      print'(3es21.14)', ref
      write(*,*) 'Selfenergy:'
      print'(3es21.14)', selfenergy
      write(*,*) 'Difference:'
      print'(3es21.14)', selfenergy-ref
   end if

end subroutine test_selfenergy_gen


subroutine test_anisotropy_gen(error, mol, make_basis, spec, ref, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Factory to create new basis set objects
   procedure(basis_maker) :: make_basis

   !> Hamiltonian specification
   class(tb_h0spec), intent(in) :: spec

   !> Reference Hamiltonian matrix
   real(wp), intent(in) :: ref(:, :)

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   class(basis_type), allocatable :: bas
   type(tb_hamiltonian) :: h0
   real(wp), allocatable :: aniso_dip(:, :)
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Setup basis set
   call make_basis(bas, mol, error=error, ng=6, scale_h0_basis=.true.)
   if (allocated(error)) return

   ! Setup Hamiltonian
   call new_hamiltonian(h0, mol, bas, spec)

   ! Calculate the anisotropic dipole correction
   allocate(aniso_dip(3, mol%nat))
   call get_anisotropy(h0, mol, aniso_dip)

   ! where(abs(aniso_dip) < thr) aniso_dip = 0.0_wp
   ! print '(*(6x,"&", 3(es20.14e1, "_wp":, ","), "&", /))', aniso_dip

   if (any(abs(aniso_dip - ref) > thr_)) then
      call test_failed(error, "Dipol anisotropy does not match.")
      write(*,*) 'reference:'
      print'(3es21.14)', ref
      write(*,*) 'Dipole anisotropy:'
      print'(3es21.14)', aniso_dip
      write(*,*) 'Difference:'
      print'(3es21.14)', aniso_dip-ref
   end if

end subroutine test_anisotropy_gen


subroutine test_anisotropy_grad(error, mol, make_basis, spec, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Factory to create new basis set objects
   procedure(basis_maker) :: make_basis

   !> Hamiltonian specification
   class(tb_h0spec), intent(in) :: spec

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   class(basis_type), allocatable :: bas
   type(tb_hamiltonian) :: h0
   integer :: iat, ic
   real(wp) :: er, el, sigma(3, 3)
   real(wp), allocatable :: dEdad(:, :), aniso_dip(:, :)
   real(wp), allocatable :: num_grad(:, :), gradient(:, :)
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Setup basis set
   call make_basis(bas, mol, error=error, ng=6, scale_h0_basis=.true.)
   if (allocated(error)) return

   ! Setup Hamiltonian
   call new_hamiltonian(h0, mol, bas, spec)

   ! Generate random energy derivatives
   allocate(dEdad(3, mol%nat), source=0.0_wp)
   call random_number(dEdad)
   dEdad = dEdad - 0.5_wp
   dEdad = 0.1_wp

   ! Calculate the numerical gradient of the anisotropic dipole term
   allocate(aniso_dip(3, mol%nat), num_grad(3, mol%nat))
   num_grad = 0.0_wp
   do iat = 1, mol%nat
      do ic = 1, 3
         ! Right hand side
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_anisotropy(h0, mol, aniso_dip)
         er = sum(aniso_dip * dEdad)
         ! Left hand side
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_anisotropy(h0, mol, aniso_dip)
         el = sum(aniso_dip * dEdad)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         num_grad(ic, iat) = 0.5_wp * (er - el) / step
      end do
   end do

   allocate(gradient(3, mol%nat), source=0.0_wp)
   call get_anisotropy_gradient(h0, mol, dEdad, gradient, sigma)

   if (any(abs(gradient - num_grad) > thr_)) then
      call test_failed(error, "Gradient of hamiltonian anisotropy does not match")
      write(*,*) "numerical gradient:"
      print'(3es21.14)', num_grad
      write(*,*) "analytical gradient:"
      print'(3es21.14)', gradient
      write(*,*) "difference:"
      print'(3es21.14)', gradient - num_grad
   end if

end subroutine test_anisotropy_grad

subroutine test_energy_anisotropy_gradient(error, mol, density, make_basis, spec, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Density matrix for this structure
   real(wp), allocatable, intent(inout) :: density(:, :, :)

   !> Factory to create new basis set objects
   procedure(basis_maker) :: make_basis

   !> Hamiltonian specification
   class(tb_h0spec), intent(in) :: spec

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic, nspin
   class(basis_type), allocatable :: bas
   type(basis_cache) :: bcache
   type(tb_hamiltonian) :: h0
   type(wavefunction_type) :: wfn_aux
   type(potential_type) :: pot
   type(adjacency_list) :: list
   real(wp) :: thr_, cutoff, sigma(3, 3)
   real(wp), allocatable :: lattr(:, :), selfenergy(:), dsedcn(:)
   real(wp), allocatable :: aniso_dip(:, :), overlap(:, :), hamiltonian(:, :)
   real(wp), allocatable :: dpint(:, :, :), qpint(:, :, :)
   real(wp), allocatable :: dEdcn(:), dEdad(:, :)
   real(wp), allocatable :: dEdqbas(:), dEdcnbas(:)
   real(wp), allocatable :: wdensity(:, :, :), gradient(:, :)
   real(wp), allocatable :: num_dEdad(:, :)
   real(wp), allocatable :: er(:), el(:)
   real(wp), parameter :: step = 1.0e-5_wp

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in
   nspin = size(density, dim=3)

   ! Setup basis set
   call make_basis(bas, mol, error=error, scale_h0_basis=.true.)
   if (allocated(error)) return

   ! Setup Hamiltonian
   call new_hamiltonian(h0, mol, bas, spec)

   ! Obtain EEQBC charges for charge adaptation
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)

   ! Update basis set with charge/CN scaling
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)
   
   ! Calculate selfenergies
   allocate(selfenergy(bas%nsh), dsedcn(bas%nsh))
   call get_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, &
      & selfenergy=selfenergy, dsedcn=dsedcn)

   ! Calculate anisotropic dipole corrections
   allocate(aniso_dip(3, mol%nat))
   call get_anisotropy(h0, mol, aniso_dip)

   ! Calculate numerical gradient w.r.t. anisotropic dipole term
   allocate(num_dEdad(3, mol%nat), er(bas%nao), el(bas%nao))
   allocate(overlap(bas%nao, bas%nao), dpint(3, bas%nao, bas%nao), &
      & qpint(6, bas%nao, bas%nao), hamiltonian(bas%nao, bas%nao))
   num_dEdad = 0.0_wp
   do iat = 1, mol%nat
      ! Dipole anisotropy derivatives
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         ! Right hand side + step
         aniso_dip(ic, iat) = aniso_dip(ic, iat) + step
         call get_hamiltonian(mol, lattr, list, bas, bcache, h0, selfenergy, &
            & aniso_dip, overlap, dpint, qpint, hamiltonian)
         call get_electronic_energy(hamiltonian, density, er)
         
         ! Left hand side - step
         aniso_dip(ic, iat) = aniso_dip(ic, iat) - 2.0_wp*step
         call get_hamiltonian(mol, lattr, list, bas, bcache, h0, selfenergy, &
            & aniso_dip, overlap, dpint, qpint, hamiltonian)
         call get_electronic_energy(hamiltonian, density, el)
      
         aniso_dip(ic, iat) = aniso_dip(ic, iat) + step   
         num_dEdad(ic, iat) = 0.5_wp * (sum(er) - sum(el)) / step
      end do
   end do

   ! Analytical gradient w.r.t. anisotropic dipole term
   allocate(wdensity(size(density, 1), size(density, 2), size(density, 3)), &
      & dEdcn(mol%nat), dEdad(3, mol%nat), & 
      & dEdcnbas(mol%nat), dEdqbas(mol%nat), gradient(3, mol%nat), source=0.0_wp)
   
   ! Setup dummy potential
   call new_potential(pot, mol, bas, nspin, .false.)
   call pot%reset()

   ! Convert the density and wdensity to total density and magnetization density
   call updown_to_magnet(density)
   call updown_to_magnet(wdensity)

   ! Build energy derivatives
   call get_hamiltonian_gradient(mol, lattr, list, bas, bcache, h0, selfenergy, &
      & dsedcn, aniso_dip, pot, density, wdensity, &
      & dEdcn, dEdad, dEdcnbas, dEdqbas, gradient, sigma)

   if (any(abs(dEdad - num_dEdad) > thr_)) then
      call test_failed(error, "Gradient of energy w.r.t the dipole aniso vector does not match")
      write(*,*) "Numerical dEdad:"
      print'(3es21.14)', num_dEdad
      write(*,*) "Analytical dEdad:"
      print'(3es21.14)', dEdad
      write(*,*) "Difference:"
      print'(3es21.14)', dEdad - num_dEdad
   end if

end subroutine test_energy_anisotropy_gradient


subroutine test_hamiltonian_gen(error, mol, make_basis, spec, make_ncoord, ref, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Factory to create new basis set objects
   procedure(basis_maker) :: make_basis

   !> Hamiltonian specification
   class(tb_h0spec), intent(in) :: spec

   !> Factory to create new coordination number objects
   procedure(ncoord_maker) :: make_ncoord

   !> Reference Hamiltonian matrix
   real(wp), intent(in) :: ref(:, :)

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   class(basis_type), allocatable :: bas
   type(basis_cache) :: bcache
   type(tb_hamiltonian) :: h0
   type(wavefunction_type) :: wfn_aux
   class(ncoord_type), allocatable :: ncoord, ncoord_en
   type(adjacency_list) :: list
   real(wp), parameter :: cn_cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :), cn(:), cn_en(:)
   real(wp), allocatable :: selfenergy(:), aniso_dip(:, :)
   real(wp), allocatable :: overlap(:, :), hamiltonian(:, :)
   real(wp), allocatable :: dpint(:, :, :), qpint(:, :, :)
   real(wp) :: cutoff, thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Setup basis set
   call make_basis(bas, mol, error=error, ng=6, scale_h0_basis=.true.)
   if (allocated(error)) return
   call check(error, bas%nao, size(ref, 1))
   if (allocated(error)) return

   ! Setup Hamiltonian
   call new_hamiltonian(h0, mol, bas, spec)

   ! Obtain coordination number
   call make_ncoord(ncoord, ncoord_en, mol, cn_cutoff, error)
   if (allocated(error)) return
   call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
   if (allocated(ncoord)) then
      allocate(cn(mol%nat))
      call ncoord%get_coordination_number(mol, lattr, cn)
   end if
   if (allocated(ncoord_en)) then
      allocate(cn_en(mol%nat))
      call ncoord_en%get_coordination_number(mol, lattr, cn_en)
   end if

   ! Obtain EEQBC charges for charge adaptation
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)

   ! Calculate self-energies
   allocate(selfenergy(bas%nsh))
   call get_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
      & selfenergy=selfenergy)

   ! Calculate the anisotropic dipole correction
   allocate(aniso_dip(3, mol%nat))
   call get_anisotropy(h0, mol, aniso_dip)

   ! Update basis set with charge/CN scaling
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   ! Build Hamiltonian matrix
   allocate(overlap(bas%nao, bas%nao), dpint(3, bas%nao, bas%nao), &
      & qpint(6, bas%nao, bas%nao), hamiltonian(bas%nao, bas%nao))
   call get_hamiltonian(mol, lattr, list, bas, bcache, h0, selfenergy, &
      & aniso_dip, overlap, dpint, qpint, hamiltonian)

   ! where(abs(hamiltonian) < thr) hamiltonian = 0.0_wp
   ! print '(*(6x,"&", 3(es20.14e1, "_wp":, ","), "&", /))', hamiltonian

   if (any(abs(hamiltonian - ref) > thr_)) then
      call test_failed(error, "Hamiltonian does not match.")
      write(*,*) 'Reference:'
      print'(3es21.14)', ref
      write(*,*) 'Hamiltonian:'
      print'(3es21.14)', hamiltonian
      write(*,*) 'Difference:'
      print'(3es21.14)', hamiltonian-ref
   end if

   !allocate(eigval(bas%nao))
   !call sygvd%solve(hamiltonian, overlap, eigval, error)
   !if (allocated(error)) return

   !print '(*("&", 3(es20.14e1, "_wp":, ","), "&", /))', eigval

end subroutine test_hamiltonian_gen


subroutine test_hamiltonian_qeff_numgrad(error, mol, density, make_basis, spec, make_ncoord, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Density matrix for this structure
   real(wp), allocatable, intent(inout) :: density(:, :, :)

   !> Factory to create new basis set objects
   procedure(basis_maker) :: make_basis

   !> Hamiltonian specification
   class(tb_h0spec), intent(in) :: spec

   !> Factory to create new coordination number objects
   procedure(ncoord_maker) :: make_ncoord

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ish, isp, nspin
   class(basis_type), allocatable :: bas
   type(basis_cache) :: bcache
   type(tb_hamiltonian) :: h0
   type(wavefunction_type) :: wfn_aux
   class(ncoord_type), allocatable :: ncoord, ncoord_en, ncoord_bas
   type(potential_type) :: pot
   type(adjacency_list) :: list
   real(wp), parameter :: cn_cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :), cn(:), cn_en(:), cnbas(:)
   real(wp), allocatable :: selfenergy(:), dsedcn(:)
   real(wp), allocatable :: aniso_dip(:, :), dEdad(:, :)
   real(wp) :: thr_, cutoff, sigma(3, 3)
   real(wp), allocatable :: overlap(:, :), hamiltonian(:, :), wdensity(:, :, :)
   real(wp), allocatable :: dpint(:, :, :), qpint(:, :, :)
   real(wp), allocatable :: er(:), el(:)
   real(wp), allocatable :: dEdcn(:), dEdqbas(:), dEdcnbas(:), gradient(:, :)
   real(wp), allocatable :: num_dEdqbas(:), num_dEdcnbas(:)
   real(wp), parameter :: step = 1.0e-5_wp

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in
   nspin = size(density, dim=3)

   ! Setup basis set and extract the coordination number
   call make_basis(bas, mol, error=error, scale_h0_basis=.true., &
      & ncoord_bas=ncoord_bas)
   if (allocated(error)) return

   ! Setup Hamiltonian
   call new_hamiltonian(h0, mol, bas, spec)

   ! Setup coordination number
   call make_ncoord(ncoord, ncoord_en, mol, cn_cutoff, error)
   if (allocated(error)) return
   call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
   if (allocated(ncoord)) allocate(cn(mol%nat))
   if (allocated(ncoord_en)) allocate(cn_en(mol%nat))
   if (allocated(ncoord)) then
      call ncoord%get_coordination_number(mol, lattr, cn)
   end if
   if (allocated(ncoord_en)) then
      call ncoord_en%get_coordination_number(mol, lattr, cn_en)
   end if
   if (allocated(ncoord_bas)) then
      allocate(cnbas(mol%nat))
      call ncoord_bas%get_coordination_number(mol, lattr, cnbas)
   end if

   ! Obtain EEQBC charges for charge adaptation
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)

   ! Update basis set with charge/CN scaling
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   ! Calculate selfenergies
   allocate(selfenergy(bas%nsh), dsedcn(bas%nsh))
   call get_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
      & selfenergy=selfenergy, dsedcn=dsedcn)

   ! Calculate the anisotropic dipole correction
   allocate(aniso_dip(3, mol%nat))
   call get_anisotropy(h0, mol, aniso_dip)

   ! Calculate numerical gradient
   allocate(num_dEdcnbas(mol%nat), num_dEdqbas(mol%nat), er(bas%nao), el(bas%nao))
   allocate(overlap(bas%nao, bas%nao), dpint(3, bas%nao, bas%nao), &
      & qpint(6, bas%nao, bas%nao), hamiltonian(bas%nao, bas%nao))
   num_dEdcnbas(:) = 0.0_wp
   num_dEdqbas(:) = 0.0_wp
   do iat = 1, mol%nat
      isp = mol%id(iat)
      ! Calculate numerical gradient for the coordination number
      ! Right hand side
      cnbas(iat) = cnbas(iat) + step
      do ish = 1, bas%nsh_id(isp)
         associate(cgto => bas%cgto(ish, isp)%raw)
            ! Perturb effective charge by coordination number
            call cgto%get_qeff(wfn_aux%qat(iat, 1), cnbas(iat), bcache%cgto(ish, iat)%qeff)
            ! Update normalization
            call cgto%get_normalization(bcache%cgto(ish, iat), .false.)
         end associate
         associate(cgto_h0 => bas%cgto_h0(ish, isp)%raw)
            ! Perturb effective charge by coordination number for scaled H0
            call cgto_h0%get_qeff(wfn_aux%qat(iat, 1), cnbas(iat), bcache%cgto_h0(ish, iat)%qeff)
            ! Update normalization of scaled H0
            call cgto_h0%get_normalization(bcache%cgto_h0(ish, iat), .false.)
         end associate
      end do
      er = 0.0_wp
      ! Build Hamiltonian matrix
      call get_hamiltonian(mol, lattr, list, bas, bcache, h0, selfenergy, &
         & aniso_dip, overlap, dpint, qpint, hamiltonian)
      ! Obtain the electronic energy contribution
      call get_electronic_energy(hamiltonian, density, er)

      ! Left hand side
      cnbas(iat) = cnbas(iat) - 2*step
      do ish = 1, bas%nsh_id(isp)
         associate(cgto => bas%cgto(ish, isp)%raw)
            ! Perturb effective charge by coordination number
            call cgto%get_qeff(wfn_aux%qat(iat, 1), cnbas(iat), bcache%cgto(ish, iat)%qeff)
            ! Update normalization
            call cgto%get_normalization(bcache%cgto(ish, iat), .false.)
         end associate
         associate(cgto_h0 => bas%cgto_h0(ish, isp)%raw)
            ! Perturb effective charge by coordination number for scaled H0
            call cgto_h0%get_qeff(wfn_aux%qat(iat, 1), cnbas(iat), bcache%cgto_h0(ish, iat)%qeff)
            ! Update normalization of scaled H0
            call cgto_h0%get_normalization(bcache%cgto_h0(ish, iat), .false.)
         end associate
      end do
      el = 0.0_wp
      ! Build Hamiltonian matrix
      call get_hamiltonian(mol, lattr, list, bas, bcache, h0, selfenergy, &
         & aniso_dip, overlap, dpint, qpint, hamiltonian)
      ! Obtain the electronic energy contribution
      call get_electronic_energy(hamiltonian, density, el)

      cnbas(iat) = cnbas(iat) + step
      do ish = 1, bas%nsh_id(isp)
         associate(cgto => bas%cgto(ish, isp)%raw)
            ! Perturb effective charge by coordination number
            call cgto%get_qeff(wfn_aux%qat(iat, 1), cnbas(iat), bcache%cgto(ish, iat)%qeff)
            ! Update normalization
            call cgto%get_normalization(bcache%cgto(ish, iat), .false.)
         end associate
         associate(cgto_h0 => bas%cgto_h0(ish, isp)%raw)
            ! Perturb effective charge by coordination number for scaled H0
            call cgto_h0%get_qeff(wfn_aux%qat(iat, 1), cnbas(iat), bcache%cgto_h0(ish, iat)%qeff)
            ! Update normalization of scaled H0
            call cgto_h0%get_normalization(bcache%cgto_h0(ish, iat), .false.)
         end associate
      end do
      num_dEdcnbas(iat) = 0.5_wp*(sum(er) - sum(el))/step


      ! Calculate numerical gradient for the atomic charges
      ! Right hand side
      wfn_aux%qat(iat, 1) = wfn_aux%qat(iat, 1) + step
      call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)
      er = 0.0_wp
      ! Build Hamiltonian matrix
      call get_hamiltonian(mol, lattr, list, bas, bcache, h0, selfenergy, &
         & aniso_dip, overlap, dpint, qpint, hamiltonian)
      ! Obtain the electronic energy contribution
      call get_electronic_energy(hamiltonian, density, er)

      ! Left hand side
      wfn_aux%qat(iat, 1) = wfn_aux%qat(iat, 1) - 2*step
      call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)
      el = 0.0_wp
      ! Build Hamiltonian matrix
      call get_hamiltonian(mol, lattr, list, bas, bcache, h0, selfenergy, &
         & aniso_dip, overlap, dpint, qpint, hamiltonian)
      ! Obtain the electronic energy contribution
      call get_electronic_energy(hamiltonian, density, el)

      wfn_aux%qat(iat, 1) = wfn_aux%qat(iat, 1) + step
      call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)
      num_dEdqbas(iat) = 0.5_wp*(sum(er) - sum(el))/step

   end do

   ! Update charges with derivatives
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .true.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1), wfn_aux%dqatdr(:, :, :, 1), &
     & wfn_aux%dqatdL(:, :, :, 1))
   if (allocated(error)) return
   call bas%update(mol, bcache, .true., wfn_aux=wfn_aux)
   
   ! Setup dummy potential
   call new_potential(pot, mol, bas, nspin, .false.)
   call pot%reset()

   ! Prepare analytic electronic energy gradient
   allocate(wdensity(size(density, 1), size(density, 2), size(density, 3)), &
      & dEdcn(mol%nat), dEdad(3, mol%nat), &
      & dEdcnbas(mol%nat), dEdqbas(mol%nat), gradient(3, mol%nat), source=0.0_wp)
   ! Convert the density and wdensity to total density and magnetization density
   call updown_to_magnet(density)
   call updown_to_magnet(wdensity)

   ! Analytic gradient of the electronic energy
   call get_hamiltonian_gradient(mol, lattr, list, bas, bcache, h0, selfenergy, &
      & dsedcn, aniso_dip, pot, density, wdensity, &
      & dEdcn, dEdad, dEdcnbas, dEdqbas, gradient, sigma)

   if (any(abs(dEdcnbas - num_dEdcnbas) > thr_)) then
      call test_failed(error, "Gradient of electronic energy w.r.t. coordination number does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', dEdcnbas
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dEdcnbas
      write(*,*) 'Difference:'
      print'(3es21.14)', dEdcnbas-num_dEdcnbas
   end if

   if (any(abs(dEdqbas - num_dEdqbas) > thr_)) then
      call test_failed(error, "Gradient of electronic energy w.r.t. atomic charge does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', dEdqbas
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', num_dEdqbas
      write(*,*) 'Difference:'
      print'(3es21.14)', dEdqbas-num_dEdqbas
   end if

end subroutine test_hamiltonian_qeff_numgrad


subroutine test_hamiltonian_numgrad(error, mol, density, make_basis, spec, make_ncoord, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Density matrix for this structure
   real(wp), allocatable, intent(inout) :: density(:, :, :)

   !> Factory to create new basis set objects
   procedure(basis_maker) :: make_basis

   !> Hamiltonian specification
   class(tb_h0spec), intent(in) :: spec

   !> Factory to create new coordination number objects
   procedure(ncoord_maker) :: make_ncoord

   !> Test threshold
   real(wp), intent(in), optional :: thr_in

   integer :: iat, ic, nspin
   class(basis_type), allocatable :: bas
   type(basis_cache) :: bcache
   type(tb_hamiltonian) :: h0
   type(wavefunction_type) :: wfn_aux
   class(ncoord_type), allocatable :: ncoord, ncoord_en
   type(potential_type) :: pot
   type(adjacency_list) :: list
   real(wp), parameter :: cn_cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :), cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: cn_en(:), dcn_endr(:, :, :), dcn_endL(:, :, :)
   real(wp), allocatable :: selfenergy(:), dsedcn(:)
   real(wp), allocatable :: aniso_dip(:, :), dEdad(:, :)
   real(wp) :: thr_, cutoff, sigma(3, 3)
   real(wp), allocatable :: overlap(:, :), hamiltonian(:, :), wdensity(:, :, :)
   real(wp), allocatable :: dpint(:, :, :), qpint(:, :, :)
   real(wp), allocatable :: er(:), el(:), numgrad(:, :)
   real(wp), allocatable :: dEdcn(:), dEdq(:), dEdqbas(:), dEdcnbas(:), gradient(:, :)
   real(wp), parameter :: step = 1.0e-5_wp

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in
   nspin = size(density, dim=3)

   ! Setup basis set
   call make_basis(bas, mol, error=error, scale_h0_basis=.true.)
   if (allocated(error)) return

   ! Setup Hamiltonian
   call new_hamiltonian(h0, mol, bas, spec)

   ! Setup coordination number
   call make_ncoord(ncoord, ncoord_en, mol, cn_cutoff, error)
   if (allocated(error)) return
   call get_lattice_points(mol%periodic, mol%lattice, cn_cutoff, lattr)
   if (allocated(ncoord)) allocate(cn(mol%nat))
   if (allocated(ncoord_en)) allocate(cn_en(mol%nat))

   ! Obtain EEQBC charges for charge adaptation
   call new_wavefunction(wfn_aux, mol%nat, 0, 0, 1, 0.0_wp, .false.)
   call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
   if (allocated(error)) return

   cutoff = get_cutoff(bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)

   ! Update basis set with charge/CN scaling
   call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)

   ! Calculate numerical gradient
   allocate(numgrad(3, mol%nat), er(bas%nao), el(bas%nao))
   allocate(selfenergy(bas%nsh), aniso_dip(3, mol%nat))
   allocate(overlap(bas%nao, bas%nao), dpint(3, bas%nao, bas%nao), &
      & qpint(6, bas%nao, bas%nao), hamiltonian(bas%nao, bas%nao))
   numgrad(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   do iat = 1, mol%nat
      do ic = 1, 3
         ! Right hand side
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         er = 0.0_wp
         ! Update coordination numbers
         if (allocated(ncoord)) then
            call ncoord%get_coordination_number(mol, lattr, cn)
         end if
         if (allocated(ncoord_en)) then
            call ncoord_en%get_coordination_number(mol, lattr, cn_en)
         end if
         ! Calculate selfenergies
         call get_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
            & selfenergy=selfenergy)
         ! Calculate the anisotropic dipole correction
         call get_anisotropy(h0, mol, aniso_dip)
         ! Update the basis set
         call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
         if (allocated(error)) return
         call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)
         ! Build Hamiltonian matrix
         call get_hamiltonian(mol, lattr, list, bas, bcache, h0, selfenergy, &
            & aniso_dip, overlap, dpint, qpint, hamiltonian)
         ! Obtain the electronic energy contribution
         call get_electronic_energy(hamiltonian, density, er)

         ! Left hand side
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         el = 0.0_wp
         ! Update coordination numbers
         if (allocated(ncoord)) then
            call ncoord%get_coordination_number(mol, lattr, cn)
         end if
         if (allocated(ncoord_en)) then
            call ncoord_en%get_coordination_number(mol, lattr, cn_en)
         end if
         ! Calculate selfenergies
         call get_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
            & selfenergy=selfenergy)
         ! Calculate the anisotropic dipole correction
         call get_anisotropy(h0, mol, aniso_dip)
         ! Update the basis set
         call get_eeqbc_charges(mol, error, wfn_aux%qat(:, 1))
         if (allocated(error)) return
         call bas%update(mol, bcache, .false., wfn_aux=wfn_aux)
         ! Build Hamiltonian matrix
         call get_hamiltonian(mol, lattr, list, bas, bcache, h0, selfenergy, &
            & aniso_dip, overlap, dpint, qpint, hamiltonian)
         ! Obtain the electronic energy contribution
         call get_electronic_energy(hamiltonian, density, el)

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

   ! Update coordination numbers and derivatives
   if (allocated(ncoord)) then
      allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
      call ncoord%get_coordination_number(mol, lattr, cn, dcndr, dcndL)
   end if
   if (allocated(ncoord_en)) then
      allocate(dcn_endr(3, mol%nat, mol%nat), dcn_endL(3, 3, mol%nat))
      call ncoord_en%get_coordination_number(mol, lattr, cn_en, dcn_endr, dcn_endL)
   end if
   
   ! Calculate selfenergies and derivatives
   allocate(dsedcn(bas%nsh))
   call get_selfenergy(h0, mol%id, bas%ish_at, bas%nsh_id, cn=cn, cn_en=cn_en, &
      & selfenergy=selfenergy, dsedcn=dsedcn)
   
   ! Calculate the anisotropic dipole correction
   call get_anisotropy(h0, mol, aniso_dip)
   
   ! Setup dummy potential
   call new_potential(pot, mol, bas, nspin, .false.)
   call pot%reset()

   ! Prepare analytic electronic energy gradient
   allocate(wdensity(size(density, 1), size(density, 2), size(density, 3)), &
      & dEdcn(mol%nat), dEdQ(mol%nat), dEdad(3, mol%nat), &
      & dEdcnbas(mol%nat), dEdqbas(mol%nat), gradient(3, mol%nat), source=0.0_wp)
   ! Convert the density and wdensity to total density and magnetization density
   call updown_to_magnet(density)
   call updown_to_magnet(wdensity)

   ! Analytic gradient of the electronic energy
   call get_hamiltonian_gradient(mol, lattr, list, bas, bcache, h0, selfenergy, &
      & dsedcn, aniso_dip, pot, density, wdensity, &
      & dEdcn, dEdad, dEdcnbas, dEdqbas, gradient, sigma)

   ! Add contribution from CN derivatives of the selfenergy
   if (allocated(dcndr)) then
      call gemv(dcndr, dEdcn, gradient, beta=1.0_wp)
   end if
   if (allocated(dcndL)) then
      call gemv(dcndL, dEdcn, sigma, beta=1.0_wp)
   end if

   ! Add contribution from the anisotropy derivatives
   call get_anisotropy_gradient(h0, mol, dEdad, gradient, sigma)

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
      call test_failed(error, "Gradient of electronic energy does not match")
      write(*,*) 'Analytic gradient:'
      print'(3es21.14)', gradient
      write(*,*) 'Numeric gradient:'
      print'(3es21.14)', numgrad
      write(*,*) 'Difference:'
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_hamiltonian_numgrad


subroutine test_selfenergy_gfn2_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 2
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & -3.91986886330804E-1_wp, -3.91986886330804E-1_wp],&
      & shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_selfenergy_gen(error, mol, make_gfn2_basis, gfn2_h0spec(mol), &
      & make_gfn2_ncoord, selfenergy)

end subroutine test_selfenergy_gfn2_h2

subroutine test_selfenergy_gfn2_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 3
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & -1.85652593278629E-1_wp, -7.93540992090888E-2_wp, -3.91761150304657E-1_wp],&
      & shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_selfenergy_gen(error, mol, make_gfn2_basis, gfn2_h0spec(mol), &
      & make_gfn2_ncoord, selfenergy)

end subroutine test_selfenergy_gfn2_lih

subroutine test_selfenergy_gfn2_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 6
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & -7.35145168851922E-1_wp, -4.17765769244959E-1_wp, -2.27200179659920E-2_wp, &
      & -7.35145168851922E-1_wp, -4.17765769244959E-1_wp, -2.27200179659920E-2_wp],&
      & shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "S2")
   call test_selfenergy_gen(error, mol, make_gfn2_basis, gfn2_h0spec(mol), &
      & make_gfn2_ncoord, selfenergy)

end subroutine test_selfenergy_gfn2_s2

subroutine test_selfenergy_gfn2_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 7
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & -5.52421014778594E-1_wp, -2.35769691509448E-1_wp, -4.13801903148179E-2_wp, &
      & -3.91823589915583E-1_wp, -3.91823589915583E-1_wp, -3.91823589915583E-1_wp, &
      & -3.91823589915583E-1_wp], shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_selfenergy_gen(error, mol, make_gfn2_basis, gfn2_h0spec(mol), &
      & make_gfn2_ncoord, selfenergy)

end subroutine test_selfenergy_gfn2_sih4


subroutine test_selfenergy_gxtb_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 2
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & -2.13219267896436E-01_wp, -2.13219267896436E-01_wp],&
      & shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_selfenergy_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & make_gxtb_ncoord, selfenergy)

end subroutine test_selfenergy_gxtb_h2

subroutine test_selfenergy_gxtb_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 3
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & -1.44548305390441E-01_wp, -1.57093875988090E-01_wp, -2.10484850069579E-01_wp],&
      & shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_selfenergy_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & make_gxtb_ncoord, selfenergy)

end subroutine test_selfenergy_gxtb_lih

subroutine test_selfenergy_gxtb_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 6
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & -2.97250714884559E-01_wp, -3.56811624241647E-01_wp, -6.57362067963437E-02_wp, &
      & -2.97250714884559E-01_wp, -3.56811624241647E-01_wp, -6.57362067963437E-02_wp],&
      & shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "S2")
   call test_selfenergy_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & make_gxtb_ncoord, selfenergy)

end subroutine test_selfenergy_gxtb_s2

subroutine test_selfenergy_gxtb_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 7
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & -3.15985884830748E-01_wp, -2.76904874686258E-01_wp, -5.03259873131728E-02_wp, &
      & -2.11107874809141E-01_wp, -2.11107874809141E-01_wp, -2.11107874809141E-01_wp, &
      & -2.11107874809141E-01_wp], shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_selfenergy_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & make_gxtb_ncoord, selfenergy)

end subroutine test_selfenergy_gxtb_sih4

subroutine test_selfenergy_gxtb_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 13
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & -2.82049302516615E-01_wp, -1.71264401650916E-01_wp, -3.30898966840816E-02_wp, &
      & -3.84756195587905E-01_wp, -3.43061374750722E-01_wp, -3.09507307342230E-01_wp, &
      & -1.31964145071846E-01_wp, -3.43064588788769E-01_wp, -3.09502658142368E-01_wp, &
      & -1.31960989638649E-01_wp, -3.43055500606573E-01_wp, -3.09515804463225E-01_wp, &
      & -1.31969912106704E-01_wp], shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "CeCl3")
   call test_selfenergy_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & make_gxtb_ncoord, selfenergy)

end subroutine test_selfenergy_gxtb_cecl3

subroutine test_selfenergy_gxtb_ce2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 8
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & -2.87625642334199E-01_wp, -1.71872094404136E-01_wp, -3.55267110074857E-02_wp, &
      & -3.86583909524310E-01_wp, -2.87625642334199E-01_wp, -1.71872094404136E-01_wp, &
      & -3.55267110074857E-02_wp, -3.86583909524310E-01_wp], shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "Ce2")
   call test_selfenergy_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & make_gxtb_ncoord, selfenergy)

end subroutine test_selfenergy_gxtb_ce2


subroutine test_selfenergy_ceh_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 2
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & -5.2057326046758E-01_wp, -5.2057326046758E-01_wp],&
      & shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_selfenergy_gen(error, mol, make_ceh_basis, ceh_h0spec(mol), &
      & make_ceh_ncoord, selfenergy)

end subroutine test_selfenergy_ceh_h2

subroutine test_selfenergy_ceh_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 3
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & -5.7614182696741E-02_wp, -1.3057703854461E-01_wp, -3.6985761349230E-01_wp],&
      & shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_selfenergy_gen(error, mol, make_ceh_basis, ceh_h0spec(mol), &
      & make_ceh_ncoord, selfenergy)

end subroutine test_selfenergy_ceh_lih

subroutine test_selfenergy_ceh_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 6
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & -6.9008304496671E-01_wp, -5.6274208401578E-01_wp, -5.7343694597688E-02_wp, & 
      & -6.9008304496671E-01_wp, -5.6274208401578E-01_wp, -5.7343694597688E-02_wp],&
      & shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "S2")
   call test_selfenergy_gen(error, mol, make_ceh_basis, ceh_h0spec(mol), &
      & make_ceh_ncoord, selfenergy)

end subroutine test_selfenergy_ceh_s2

subroutine test_selfenergy_ceh_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 7
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & -7.0849504464403E-01_wp, -4.7605638972741E-01_wp, -1.8541704653682E-01_wp, & 
      & -4.8652644697198E-01_wp, -4.8652644697198E-01_wp, -4.8652644697198E-01_wp, &
      & -4.8652644697198E-01_wp], shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_selfenergy_gen(error, mol, make_ceh_basis, ceh_h0spec(mol), &
      & make_ceh_ncoord, selfenergy)

end subroutine test_selfenergy_ceh_sih4

subroutine test_selfenergy_ceh_accl6(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nsh = 22
   real(wp), parameter :: selfenergy(nsh) = reshape([&
      & 3.47146858519381E-00_wp, -4.35870279866614E-01_wp, -5.8999616375307806E-01_wp, &
      & 5.07577833416763E-02_wp, -3.63354240235948E-01_wp, -1.8456541144073824E-01_wp, &
      & 9.62240909632790E-02_wp, -3.63356385871318E-01_wp, -1.8459106004014270E-01_wp, &
      & 9.62177622670393E-02_wp, -3.63388268461152E-01_wp, -1.8459960442278944E-01_wp, &
      & 9.62157809840490E-02_wp, -3.63390391126379E-01_wp, -1.8461247387922086E-01_wp, &
      & 9.62126097565016E-02_wp, -3.63357326972503E-01_wp, -1.8458205115288362E-01_wp, &
      & 9.62199920825849E-02_wp, -3.63349003068365E-01_wp, -1.8457090877938098E-01_wp, &
      & 9.62227112995023E-02_wp], shape(selfenergy))

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "AcCl6")
   call test_selfenergy_gen(error, mol, make_ceh_basis, ceh_h0spec(mol), &
      & make_ceh_ncoord, selfenergy)

end subroutine test_selfenergy_ceh_accl6


subroutine test_anisotropy_gxtb_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 2
   real(wp), parameter :: aniso_dip(3, nat) = reshape([&
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.78247703054970E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.78247703054970E-03_wp],&
      & shape(aniso_dip))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_anisotropy_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & aniso_dip)

end subroutine test_anisotropy_gxtb_h2

subroutine test_anisotropy_gxtb_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 2
   real(wp), parameter :: aniso_dip(3, nat) = reshape([&
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  5.73035312260568E-04_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.73035312260568E-04_wp],&
      & shape(aniso_dip))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_anisotropy_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & aniso_dip)

end subroutine test_anisotropy_gxtb_lih

subroutine test_anisotropy_gxtb_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 5
   real(wp), parameter :: aniso_dip(3, nat) = reshape([&
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.42403063051974E-03_wp,  1.42403063051974E-03_wp, -1.42403063051974E-03_wp, &
      & -1.42403063051974E-03_wp, -1.42403063051974E-03_wp, -1.42403063051974E-03_wp, &
      &  1.42403063051974E-03_wp, -1.42403063051974E-03_wp,  1.42403063051974E-03_wp, &
      & -1.42403063051974E-03_wp,  1.42403063051974E-03_wp,  1.42403063051974E-03_wp],&
      & shape(aniso_dip))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_anisotropy_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & aniso_dip)

end subroutine test_anisotropy_gxtb_sih4

subroutine test_anisotropy_gxtb_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 4
   real(wp), parameter :: aniso_dip(3, nat) = reshape([&
      & -2.40272071326107E-05_wp, -5.36541530689848E-05_wp,  7.25187122682972E-05_wp, &
      &  1.17417711234263E-04_wp,  1.72447327137560E-04_wp,  1.07864385539449E-04_wp, &
      &  1.21044734519321E-04_wp, -1.76759186066648E-04_wp, -8.88789866388204E-05_wp, &
      & -2.14435238620973E-04_wp,  5.79660119980730E-05_wp, -9.15041111689260E-05_wp],&
      & shape(aniso_dip))

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "CeCl3")
   call test_anisotropy_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & aniso_dip)

end subroutine test_anisotropy_gxtb_cecl3

subroutine test_anisotropy_gxtb_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 16
   real(wp), parameter :: aniso_dip(3, nat) = reshape([&
      & -8.84626352915912E-05_wp,  1.64076823423698E-05_wp,  4.34533831880603E-04_wp, &
      &  2.34336187941528E-03_wp, -2.87955822245846E-04_wp, -1.22622047501200E-03_wp, &
      &  6.72117399017162E-05_wp, -1.58348973119830E-03_wp,  3.96564528721346E-04_wp, &
      & -1.89595063777219E-03_wp,  1.34283071936362E-03_wp, -1.24771243019455E-03_wp, &
      &  3.09662564980785E-03_wp,  1.81907173480367E-03_wp, -5.85584199120327E-03_wp, &
      & -1.68483558036286E-03_wp, -9.90523488368302E-04_wp, -8.63045709234935E-04_wp, &
      &  5.67162306639455E-04_wp,  1.27213205976199E-03_wp,  1.95425458175822E-03_wp, &
      &  1.33984020638103E-03_wp,  5.53137610836614E-04_wp, -6.56768042556105E-04_wp, &
      & -1.28675351499619E-03_wp,  8.39048421708812E-04_wp, -1.08168843191651E-03_wp, &
      & -1.62405229927545E-03_wp, -4.75329516575678E-04_wp,  4.63915038469426E-04_wp, &
      & -2.07905986216747E-04_wp,  1.92950736822851E-04_wp,  2.11747652769870E-03_wp, &
      &  1.31356785830739E-04_wp,  4.19302237070717E-04_wp,  1.77134046313117E-04_wp, &
      & -8.87138322320912E-04_wp, -2.14605042991286E-03_wp,  1.35077459095300E-04_wp, &
      & -4.09267819074569E-03_wp,  1.52655182602118E-03_wp,  6.04268116369991E-03_wp, &
      &  4.60066591691792E-04_wp, -7.66763578046723E-04_wp, -2.32857268172792E-05_wp, &
      &  3.76215200731378E-03_wp, -1.73132046238412E-03_wp, -7.67074370701964E-04_wp],&
      & shape(aniso_dip))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_anisotropy_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & aniso_dip)

end subroutine test_anisotropy_gxtb_m01

subroutine test_anisotropy_gxtb_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 16
   real(wp), parameter :: aniso_dip(3, nat) = reshape([&
      & -6.73779304838932E-04_wp, -1.04760461066428E-03_wp, -2.26861344079635E-03_wp, &
      &  3.99809449346525E-04_wp, -1.86661973264431E-04_wp,  1.43857584426252E-04_wp, &
      & -1.09537428454888E-03_wp, -5.71372948408353E-04_wp, -1.73149848205811E-03_wp, &
      &  1.24969345250539E-03_wp, -8.00535922113103E-04_wp,  1.34126043978537E-03_wp, &
      & -1.10472198113902E-03_wp,  8.04570496004220E-04_wp, -3.73703686556675E-04_wp, &
      & -1.86694694753554E-05_wp,  1.39071567131937E-03_wp, -7.87258315701775E-04_wp, &
      & -1.47976523417233E-04_wp, -1.61090330964705E-03_wp,  2.09371540515845E-03_wp, &
      & -3.27408245200211E-04_wp,  1.08979923852725E-03_wp,  2.39836219288716E-03_wp, &
      &  2.81363247255721E-04_wp,  1.96991766907280E-03_wp, -1.23212576480048E-03_wp, &
      & -1.49322262769944E-03_wp, -2.79547313714176E-04_wp, -3.71281659844228E-04_wp, &
      & -1.38295352251774E-03_wp,  2.69610079208768E-03_wp,  5.84358117168985E-04_wp, &
      &  1.12403640784570E-03_wp,  1.50892786731795E-03_wp,  1.53263697885809E-04_wp, &
      & -3.55555761047778E-04_wp, -2.03303688376572E-03_wp,  1.03123118620759E-04_wp, &
      &  1.76781575702663E-03_wp, -1.38804064696822E-03_wp,  9.64139810302573E-04_wp, &
      &  1.23269203122986E-03_wp, -9.50348470534244E-04_wp, -1.01170176146506E-03_wp, &
      &  5.44251374674761E-04_wp, -5.91979655249699E-04_wp, -5.89725501267750E-06_wp],&
      & shape(aniso_dip))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_anisotropy_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & aniso_dip)

end subroutine test_anisotropy_gxtb_m02


subroutine test_hamiltonian_gfn2_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 2
   ! Updated reference values after fixing Hamiltonian bug in the original xtb
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      &-3.91986886284352E-1_wp,-4.69784177830833E-1_wp,-4.69784177830833E-1_wp,&
      &-3.91986886284352E-1_wp],shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_hamiltonian_gen(error, mol, make_gfn2_basis, gfn2_h0spec(mol), &
      & make_gfn2_ncoord, hamiltonian, thr_in=thr1)

end subroutine test_hamiltonian_gfn2_h2

subroutine test_hamiltonian_gfn2_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 5
   ! Updated reference values after fixing Hamiltonian bug in the original xtb
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      &-1.85652593279743E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.04060203683820E-1_wp, 0.00000000000000E+0_wp,&
      &-7.93540992031919E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-7.93540992031919E-2_wp, 0.00000000000000E+0_wp,-2.64332064837894E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-7.93540992031919E-2_wp, 0.00000000000000E+0_wp,-2.04060203683820E-1_wp,&
      & 0.00000000000000E+0_wp,-2.64332064837894E-1_wp, 0.00000000000000E+0_wp,&
      &-3.91761150258231E-1_wp],shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_hamiltonian_gen(error, mol, make_gfn2_basis, gfn2_h0spec(mol), &
      & make_gfn2_ncoord, hamiltonian, thr_in=thr1)

end subroutine test_hamiltonian_gfn2_lih

subroutine test_hamiltonian_gfn2_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 18
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      &-7.35145168755863E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.92782989368414E-1_wp, 0.00000000000000E+0_wp, 2.36427133955232E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.05951879709695E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-4.17765769244149E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-9.33756600277693E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.20176383632229E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-4.17765769244149E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.36427133955232E-1_wp, 0.00000000000000E+0_wp, 2.58607490728104E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.23733826451155E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.17765769244149E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-9.33756600277693E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.20176383632229E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.27200179621343E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.00344447885664E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.27200179621343E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.20176383632229E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.45142819682488E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.27200179621343E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.05951879709695E-1_wp, 0.00000000000000E+0_wp, 1.23733826451155E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-9.40352474081307E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.27200179621343E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.20176383632229E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 2.45142819682488E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.27200179621343E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.00344447885664E-2_wp,&
      &-1.92782989368414E-1_wp, 0.00000000000000E+0_wp,-2.36427133955232E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.05951879709695E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-7.35145168755863E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-9.33756600277693E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.20176383632229E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-4.17765769244149E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 2.36427133955232E-1_wp, 0.00000000000000E+0_wp, 2.58607490728104E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.23733826451155E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-4.17765769244149E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-9.33756600277693E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.20176383632229E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.17765769244149E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.00344447885664E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.27200179621343E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 1.20176383632229E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 2.45142819682488E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.27200179621343E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.05951879709695E-1_wp, 0.00000000000000E+0_wp,-1.23733826451155E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-9.40352474081307E-3_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.27200179621343E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 1.20176383632229E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 2.45142819682488E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.27200179621343E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-1.00344447885664E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-2.27200179621343E-2_wp],&
      & shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "S2")
   call test_hamiltonian_gen(error, mol, make_gfn2_basis, gfn2_h0spec(mol), &
      & make_gfn2_ncoord, hamiltonian, thr_in=thr1)

end subroutine test_hamiltonian_gfn2_s2

subroutine test_hamiltonian_gfn2_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 13
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      &-5.52421014706411E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-3.36004320989251E-1_wp,-3.36004320989251E-1_wp,-3.36004320989251E-1_wp,&
      &-3.36004320989251E-1_wp, 0.00000000000000E+0_wp,-2.35769691508991E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.53874698572085E-1_wp, 1.53874698572085E-1_wp,&
      & 1.53874698572085E-1_wp,-1.53874698572085E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-2.35769691508991E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.53874698572085E-1_wp,&
      & 1.53874698572085E-1_wp,-1.53874698572085E-1_wp,-1.53874698572085E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-2.35769691508991E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-1.53874698572085E-1_wp, 1.53874698572085E-1_wp,-1.53874698572085E-1_wp,&
      & 1.53874698572085E-1_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-4.13801903077918E-2_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-1.23912378597848E-1_wp,-1.23912378597848E-1_wp,&
      & 1.23912378597848E-1_wp, 1.23912378597848E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-4.13801903077918E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 1.23912378597848E-1_wp,&
      &-1.23912378597848E-1_wp, 1.23912378597848E-1_wp,-1.23912378597848E-1_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-4.13801903077918E-2_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,-4.13801903077918E-2_wp,&
      & 0.00000000000000E+0_wp, 1.23912378597848E-1_wp,-1.23912378597848E-1_wp,&
      &-1.23912378597848E-1_wp, 1.23912378597848E-1_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp,-4.13801903077918E-2_wp, 0.00000000000000E+0_wp,&
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp,&
      &-3.36004320989251E-1_wp,-1.53874698572085E-1_wp, 1.53874698572085E-1_wp,&
      &-1.53874698572085E-1_wp,-1.23912378597848E-1_wp, 1.23912378597848E-1_wp,&
      & 0.00000000000000E+0_wp, 1.23912378597848E-1_wp, 0.00000000000000E+0_wp,&
      &-3.91823589869151E-1_wp,-4.31486730191658E-2_wp,-4.31486730191658E-2_wp,&
      &-4.31486730191658E-2_wp,-3.36004320989251E-1_wp, 1.53874698572085E-1_wp,&
      & 1.53874698572085E-1_wp, 1.53874698572085E-1_wp,-1.23912378597848E-1_wp,&
      &-1.23912378597848E-1_wp, 0.00000000000000E+0_wp,-1.23912378597848E-1_wp,&
      & 0.00000000000000E+0_wp,-4.31486730191658E-2_wp,-3.91823589869151E-1_wp,&
      &-4.31486730191658E-2_wp,-4.31486730191658E-2_wp,-3.36004320989251E-1_wp,&
      & 1.53874698572085E-1_wp,-1.53874698572085E-1_wp,-1.53874698572085E-1_wp,&
      & 1.23912378597848E-1_wp, 1.23912378597848E-1_wp, 0.00000000000000E+0_wp,&
      &-1.23912378597848E-1_wp, 0.00000000000000E+0_wp,-4.31486730191658E-2_wp,&
      &-4.31486730191658E-2_wp,-3.91823589869151E-1_wp,-4.31486730191658E-2_wp,&
      &-3.36004320989251E-1_wp,-1.53874698572085E-1_wp,-1.53874698572085E-1_wp,&
      & 1.53874698572085E-1_wp, 1.23912378597848E-1_wp,-1.23912378597848E-1_wp,&
      & 0.00000000000000E+0_wp, 1.23912378597848E-1_wp, 0.00000000000000E+0_wp,&
      &-4.31486730191658E-2_wp,-4.31486730191658E-2_wp,-4.31486730191658E-2_wp,&
      &-3.91823589869151E-1_wp],shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_hamiltonian_gen(error, mol, make_gfn2_basis, gfn2_h0spec(mol), &
      & make_gfn2_ncoord, hamiltonian, thr_in=thr1)

end subroutine test_hamiltonian_gfn2_sih4


subroutine test_hamiltonian_gxtb_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 2
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -2.13219267896436E-01_wp, -3.44938325597689E-01_wp, -3.44938325597689E-01_wp, &
      & -2.13219267896436E-01_wp], shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_hamiltonian_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & make_gxtb_ncoord, hamiltonian, thr_in=thr1)

end subroutine test_hamiltonian_gxtb_h2

subroutine test_hamiltonian_gxtb_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 5
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -1.44548305390441E-01_wp,  0.00000000000000E+00_wp, -1.14908124967573E-03_wp, &
      &  0.00000000000000E+00_wp, -1.65220670435799E-01_wp,  0.00000000000000E+00_wp, &
      & -1.57093875988090E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.14908124967573E-03_wp,  0.00000000000000E+00_wp, &
      & -1.57093875988090E-01_wp,  0.00000000000000E+00_wp, -2.19328788817602E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.57093875988090E-01_wp,  0.00000000000000E+00_wp, -1.65220670435799E-01_wp, &
      &  0.00000000000000E+00_wp, -2.19328788817602E-01_wp,  0.00000000000000E+00_wp, &
      & -2.10484850069579E-01_wp],shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_hamiltonian_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & make_gxtb_ncoord, hamiltonian, thr_in=thr1)

end subroutine test_hamiltonian_gxtb_lih

subroutine test_hamiltonian_gxtb_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 18
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -2.97250714884559E-01_wp,  0.00000000000000E+00_wp, -3.37367040176380E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.22486203195724E-01_wp,  0.00000000000000E+00_wp,  1.91625974434395E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -8.75394921306863E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -3.56811624241647E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.73935956395111E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.24754361444952E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.00075329638765E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.37367040176380E-03_wp,  0.00000000000000E+00_wp, -3.56811624241647E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.16313996330869E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.91625974434395E-01_wp,  0.00000000000000E+00_wp,  2.24671700962941E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.30315630241442E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.56811624241647E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.73935956395111E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.24754361444952E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.00075329638765E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -6.57362067963437E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.13544543871412E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.73935956395111E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -6.57362067963437E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.00075329638765E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  6.27596935976111E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -3.16313996330869E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -6.57362067963437E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -8.75394921306862E-02_wp,  0.00000000000000E+00_wp,  1.30315630241442E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.27399676452334E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.73935956395111E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -6.57362067963437E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.00075329638765E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  6.27596935976111E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -6.57362067963437E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.13544543871412E-02_wp, &
      & -1.22486203195724E-01_wp,  0.00000000000000E+00_wp, -1.91625974434395E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -8.75394921306862E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.97250714884559E-01_wp,  0.00000000000000E+00_wp,  3.37367040176380E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.24754361444952E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.00075329638765E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -3.56811624241647E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.73935956395111E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.91625974434395E-01_wp,  0.00000000000000E+00_wp,  2.24671700962941E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.30315630241442E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  3.37367040176380E-03_wp,  0.00000000000000E+00_wp, -3.56811624241647E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  3.16313996330869E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.24754361444952E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.00075329638765E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.56811624241647E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.73935956395111E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.13544543871412E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -6.57362067963437E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.00075329638765E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  6.27596935976111E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.73935956395111E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -6.57362067963437E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -8.75394921306863E-02_wp,  0.00000000000000E+00_wp, -1.30315630241442E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.27399676452334E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  3.16313996330869E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -6.57362067963437E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.00075329638765E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  6.27596935976111E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.73935956395111E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -6.57362067963437E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.13544543871412E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -6.57362067963437E-02_wp],&
      & shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "S2")
   call test_hamiltonian_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & make_gxtb_ncoord, hamiltonian, thr_in=thr1)

end subroutine test_hamiltonian_gxtb_s2

subroutine test_hamiltonian_gxtb_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 13
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -3.15985884830748E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.15000538121138E-01_wp, -2.15000538121138E-01_wp, -2.15000538121138E-01_wp, &
      & -2.15000538121138E-01_wp,  0.00000000000000E+00_wp, -2.76904874686258E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.48050962773337E-01_wp,  1.48050962773337E-01_wp, &
      &  1.48050962773337E-01_wp, -1.48050962773337E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.76904874686258E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.48050962773337E-01_wp, &
      &  1.48050962773337E-01_wp, -1.48050962773337E-01_wp, -1.48050962773337E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.76904874686258E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.48050962773337E-01_wp,  1.48050962773337E-01_wp, -1.48050962773337E-01_wp, &
      &  1.48050962773337E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.03259873131728E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -7.52781979076199E-02_wp, -7.52781979076199E-02_wp, &
      &  7.52781979076199E-02_wp,  7.52781979076199E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.03259873131728E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  7.52781979076199E-02_wp, &
      & -7.52781979076199E-02_wp,  7.52781979076199E-02_wp, -7.52781979076199E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.03259873131728E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -7.02192610787061E-18_wp, -7.02192610787061E-18_wp, -7.02192610787061E-18_wp, &
      & -7.02192610787061E-18_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.03259873131728E-02_wp, &
      &  0.00000000000000E+00_wp,  7.52781979076199E-02_wp, -7.52781979076199E-02_wp, &
      & -7.52781979076199E-02_wp,  7.52781979076199E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.03259873131728E-02_wp,  1.26983892219261E-18_wp, &
      &  1.26983892219261E-18_wp,  1.26983892219261E-18_wp,  1.26983892219261E-18_wp, &
      & -2.15000538121138E-01_wp, -1.48050962773337E-01_wp,  1.48050962773337E-01_wp, &
      & -1.48050962773337E-01_wp, -7.52781979076199E-02_wp,  7.52781979076199E-02_wp, &
      & -7.02192610787061E-18_wp,  7.52781979076199E-02_wp,  1.26983892219261E-18_wp, &
      & -2.11107874809141E-01_wp, -5.46460522039508E-02_wp, -5.46460522039508E-02_wp, &
      & -5.46460522039508E-02_wp, -2.15000538121138E-01_wp,  1.48050962773337E-01_wp, &
      &  1.48050962773337E-01_wp,  1.48050962773337E-01_wp, -7.52781979076199E-02_wp, &
      & -7.52781979076199E-02_wp, -7.02192610787061E-18_wp, -7.52781979076199E-02_wp, &
      &  1.26983892219261E-18_wp, -5.46460522039508E-02_wp, -2.11107874809141E-01_wp, &
      & -5.46460522039508E-02_wp, -5.46460522039508E-02_wp, -2.15000538121138E-01_wp, &
      &  1.48050962773337E-01_wp, -1.48050962773337E-01_wp, -1.48050962773337E-01_wp, &
      &  7.52781979076199E-02_wp,  7.52781979076199E-02_wp, -7.02192610787061E-18_wp, &
      & -7.52781979076199E-02_wp,  1.26983892219261E-18_wp, -5.46460522039508E-02_wp, &
      & -5.46460522039508E-02_wp, -2.11107874809141E-01_wp, -5.46460522039508E-02_wp, &
      & -2.15000538121138E-01_wp, -1.48050962773337E-01_wp, -1.48050962773337E-01_wp, &
      &  1.48050962773337E-01_wp,  7.52781979076199E-02_wp, -7.52781979076199E-02_wp, &
      & -7.02192610787061E-18_wp,  7.52781979076199E-02_wp,  1.26983892219261E-18_wp, &
      & -5.46460522039508E-02_wp, -5.46460522039508E-02_wp, -5.46460522039508E-02_wp, &
      & -2.11107874809141E-01_wp], shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_hamiltonian_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & make_gxtb_ncoord, hamiltonian, thr_in=thr1)

end subroutine test_hamiltonian_gxtb_sih4


subroutine test_hamiltonian_gxtb_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 43
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -2.82049302516615E-01_wp,  1.44249297095608E-04_wp, -1.94966701972367E-04_wp, &
      &  6.45971940996523E-05_wp,  7.30868969135059E-23_wp,  3.52696966316430E-21_wp, &
      &  2.32481781012070E-21_wp, -1.18967618614470E-21_wp,  4.33699871562573E-21_wp, &
      & -2.38433810047813E-21_wp, -9.34911860900292E-37_wp, -2.97984598215748E-37_wp, &
      &  1.85878247589460E-36_wp, -2.13549118334625E-21_wp,  5.29429411528367E-37_wp, &
      &  2.13549118334625E-21_wp, -1.22327615420258E-01_wp, -5.50799865621865E-02_wp, &
      & -3.44795609739750E-02_wp, -3.75066680349345E-02_wp, -4.35783008427281E-03_wp, &
      & -4.00713627845440E-03_wp,  1.25536395265239E-03_wp, -2.72874000813000E-03_wp, &
      &  1.71594545101538E-03_wp, -1.22264710871624E-01_wp,  5.72416997864790E-02_wp, &
      &  2.87567955170003E-02_wp, -3.91833859030749E-02_wp,  4.74237605027460E-03_wp, &
      & -3.47983485181397E-03_wp,  1.92821272181122E-03_wp,  2.38169502570973E-03_wp, &
      &  1.84157392408560E-03_wp, -1.22460964229986E-01_wp, -1.80789789843406E-02_wp, &
      &  2.85423127657356E-02_wp,  6.69377814947157E-02_wp,  2.53842716456936E-03_wp, &
      &  1.08208897724546E-03_wp,  1.92548114105724E-03_wp, -4.00772120793916E-03_wp, &
      & -4.35807794043936E-03_wp,  1.44249297095608E-04_wp, -1.71264401650916E-01_wp, &
      &  2.91849642575157E-22_wp,  3.26196065363849E-21_wp,  2.32143355452967E-05_wp, &
      & -7.00653101551784E-05_wp, -2.99292395117194E-05_wp, -1.13328138613986E-37_wp, &
      & -5.18389634661959E-05_wp, -1.67528584467143E-21_wp, -8.72808114945579E-22_wp, &
      & -7.19034672801566E-22_wp,  1.24541508040964E-21_wp, -1.60268886077255E-22_wp, &
      &  1.60782395518071E-21_wp, -6.20718726692424E-22_wp,  1.27854242498941E-01_wp, &
      & -3.00111123699892E-02_wp,  3.45024665091326E-02_wp,  3.75303108914129E-02_wp, &
      & -3.17082096701095E-03_wp, -2.91430593493387E-03_wp,  3.96576726638335E-03_wp, &
      &  1.42880749622715E-03_wp,  7.10720144856875E-03_wp, -1.32703366657036E-01_wp, &
      & -2.56636215505385E-02_wp,  2.98766519495830E-02_wp, -4.07102422249159E-02_wp, &
      & -3.10850629609557E-03_wp,  2.28176804805362E-03_wp, -3.75035959129963E-03_wp, &
      &  1.29667791400357E-03_wp, -7.31261979409684E-03_wp,  4.20157391316526E-02_wp, &
      & -7.92556852535859E-02_wp, -9.37929446104940E-03_wp, -2.19972176845851E-02_wp, &
      &  9.29502413478776E-03_wp,  3.96252195216349E-03_wp,  1.18692149774725E-03_wp, &
      &  6.85553886364135E-04_wp,  3.37301693901983E-03_wp, -1.94966701972367E-04_wp, &
      &  2.91849642575157E-22_wp, -1.71264401650916E-01_wp, -4.53954640887330E-21_wp, &
      & -1.13328138613986E-37_wp,  5.18389634661959E-05_wp, -8.09044513578937E-05_wp, &
      &  2.32143355452967E-05_wp, -8.25660808937685E-38_wp,  8.63705225287944E-38_wp, &
      &  5.06814718062176E-22_wp, -2.03375430998060E-21_wp,  4.37412401264199E-21_wp, &
      & -1.10402464140243E-21_wp, -7.49692080369944E-22_wp,  5.84172583786121E-38_wp, &
      &  7.99617118035541E-02_wp,  3.44692060710472E-02_wp, -6.35473649848168E-02_wp, &
      &  2.34718724175513E-02_wp,  1.42631249047077E-03_wp, -6.69257556262513E-03_wp, &
      & -6.19800318861863E-03_wp, -4.55736928349791E-03_wp, -5.61513696179550E-04_wp, &
      & -6.67360487355299E-02_wp,  2.99090550798106E-02_wp, -7.01120822706820E-02_wp, &
      & -2.04731091749614E-02_wp,  1.29921644154981E-03_wp,  7.36400805627289E-03_wp, &
      &  5.35181103965646E-03_wp, -5.04039164139594E-03_wp,  5.04989555182662E-04_wp, &
      & -6.63241385241726E-02_wp, -9.37803943936767E-03_wp, -7.03908687281913E-02_wp, &
      &  3.47237944532958E-02_wp,  6.85386149539902E-04_wp, -2.33498037015163E-03_wp, &
      &  5.30964575093812E-03_wp,  8.64699576834706E-03_wp, -1.17771780132484E-03_wp, &
      &  6.45971940996523E-05_wp,  3.26196065363849E-21_wp, -4.53954640887330E-21_wp, &
      & -1.71264401650916E-01_wp,  5.18389634661959E-05_wp, -1.13328138613986E-37_wp, &
      & -1.34028028761354E-05_wp, -7.00653101551784E-05_wp,  2.32143355452967E-05_wp, &
      &  6.20718726692424E-22_wp, -1.60782395518071E-21_wp, -1.60268886077255E-22_wp, &
      &  6.76074258723784E-22_wp, -5.56833293700226E-22_wp, -8.72808114945580E-22_wp, &
      &  3.07445655464504E-22_wp,  8.70535260943938E-02_wp,  3.75263899004499E-02_wp, &
      &  2.34920665793948E-02_wp, -5.95719246562874E-02_wp, -6.45228663309892E-03_wp, &
      &  1.42853729713138E-03_wp,  2.70061544651703E-03_wp, -4.03948280660208E-03_wp, &
      & -6.06232790463037E-03_wp,  9.08810651002315E-02_wp, -4.07300329452995E-02_wp, &
      & -2.04608707777862E-02_wp, -5.72568009065204E-02_wp,  6.54751512891106E-03_wp, &
      &  1.29781635285664E-03_wp,  2.56648577780100E-03_wp,  3.28899897889399E-03_wp, &
      & -6.38013719371571E-03_wp, -1.55408826901311E-01_wp, -2.19742878982292E-02_wp, &
      &  3.46922406378281E-02_wp, -3.83302369012886E-03_wp, -1.02076481691057E-03_wp, &
      &  6.84246145643275E-04_wp, -4.39726862536289E-03_wp,  1.61134893933712E-03_wp, &
      &  6.96923546958978E-03_wp,  7.30868969135059E-23_wp,  2.32143355452967E-05_wp, &
      & -1.13328138613986E-37_wp,  5.18389634661959E-05_wp, -3.30898966840816E-02_wp, &
      &  3.71719534051773E-22_wp,  4.46297333397708E-22_wp, -5.56737885815481E-22_wp, &
      &  0.00000000000000E+00_wp,  1.36064515088429E-05_wp, -3.35309640479922E-05_wp, &
      & -3.51317067298182E-06_wp,  8.72052034038496E-37_wp, -7.84511474868008E-06_wp, &
      &  0.00000000000000E+00_wp, -3.03839987707241E-05_wp, -7.97812857625614E-02_wp, &
      & -4.84524422734976E-02_wp, -6.03467200158205E-02_wp,  4.82889756857431E-03_wp, &
      & -8.71650801556639E-03_wp, -1.65027194137899E-02_wp, -9.10037767448247E-03_wp, &
      & -3.88958831715503E-03_wp,  1.10719455673932E-02_wp,  8.63568720380417E-02_wp, &
      & -5.83479881335090E-02_wp, -5.45064290098753E-02_wp,  1.12280142171848E-03_wp, &
      & -1.09947461474657E-02_wp,  1.60903803016326E-02_wp,  5.22815186685464E-03_wp, &
      & -4.73227358642259E-03_wp, -1.27802615861412E-02_wp,  4.68870629578708E-02_wp, &
      & -6.72860880385563E-02_wp, -2.93962056981400E-02_wp, -4.56745992964253E-02_wp, &
      &  1.24385119738042E-02_wp,  9.71299399288886E-03_wp,  2.78692490239979E-03_wp, &
      &  1.15672992586680E-02_wp,  1.66316199240656E-02_wp,  3.52696966316430E-21_wp, &
      & -7.00653101551784E-05_wp,  5.18389634661959E-05_wp, -1.13328138613986E-37_wp, &
      &  3.71719534051773E-22_wp, -3.30898966840816E-02_wp, -3.21432768243631E-22_wp, &
      & -3.86504828363668E-22_wp,  5.56737885815481E-22_wp,  8.17062954224834E-21_wp, &
      &  1.11096211355291E-05_wp, -4.24136874131500E-05_wp, -1.92165281078489E-05_wp, &
      & -5.34020628137626E-37_wp, -2.48084311112085E-05_wp,  6.89417666436082E-37_wp, &
      & -7.32828514503475E-02_wp, -4.45063128648662E-02_wp,  1.50786081940942E-02_wp, &
      & -6.03304570737948E-02_wp, -1.64973231293152E-02_wp, -5.91499794173133E-03_wp, &
      &  1.53633897818081E-02_wp, -1.24811473959793E-03_wp, -3.52293118440440E-03_wp, &
      & -6.34132609908191E-02_wp,  4.28456294498154E-02_wp, -3.30931552610038E-02_wp, &
      & -5.45166480109714E-02_wp,  1.60940501411438E-02_wp, -8.92957506230816E-04_wp, &
      &  1.66594585688245E-02_wp, -3.95496295995601E-03_wp, -2.45202972143322E-03_wp, &
      &  2.00098126116035E-02_wp, -2.87162310555367E-02_wp,  1.06805031152813E-02_wp, &
      & -2.94045683181393E-02_wp,  9.71665080564901E-03_wp, -6.18518502445096E-03_wp, &
      & -5.29715029859336E-03_wp, -2.24309708993980E-03_wp,  1.08435978400652E-02_wp, &
      &  2.32481781012070E-21_wp, -2.99292395117194E-05_wp, -8.09044513578937E-05_wp, &
      & -1.34028028761354E-05_wp,  4.46297333397708E-22_wp, -3.21432768243631E-22_wp, &
      & -3.30898966840816E-02_wp,  2.14612373047833E-22_wp,  8.50269723022520E-22_wp, &
      & -1.51129341109871E-21_wp, -5.03553923197153E-53_wp,  2.71762746718437E-05_wp, &
      & -4.49865089774473E-05_wp,  1.21699802025309E-05_wp,  1.07055984358930E-37_wp, &
      &  6.28440688161235E-37_wp,  2.30448328514989E-02_wp, -1.28455561910549E-02_wp, &
      &  6.83504280096736E-02_wp, -8.74562697366483E-03_wp, -9.09347989988447E-03_wp, &
      &  1.53629073822107E-02_wp,  7.11152645967737E-03_wp,  1.04610161044146E-02_wp, &
      &  3.58184950097580E-03_wp,  3.50622164993731E-02_wp, -1.80614896024256E-03_wp, &
      & -6.45901676781609E-02_wp,  1.23206159814852E-03_wp,  5.23416140050281E-03_wp, &
      &  1.66618117420543E-02_wp,  1.40607882463047E-03_wp, -1.14074143898209E-02_wp, &
      &  2.02859278523153E-03_wp,  3.54875371122529E-02_wp,  6.84419358448112E-04_wp, &
      & -6.45840327694736E-02_wp, -2.54957147564150E-03_wp,  2.79196298847885E-03_wp, &
      & -5.29353832250156E-03_wp,  1.23195416282767E-03_wp,  1.95920615572102E-02_wp, &
      & -4.78485486891422E-03_wp, -1.18967618614470E-21_wp, -1.13328138613986E-37_wp, &
      &  2.32143355452967E-05_wp, -7.00653101551784E-05_wp, -5.56737885815481E-22_wp, &
      & -3.86504828363668E-22_wp,  2.14612373047833E-22_wp, -3.30898966840816E-02_wp, &
      &  3.71719534051772E-22_wp, -6.89417666436082E-37_wp,  2.48084311112085E-05_wp, &
      & -5.34020628137626E-37_wp, -8.60547552811564E-06_wp, -4.24136874131500E-05_wp, &
      &  1.11096211355291E-05_wp, -8.17062954224834E-21_wp, -4.98969803349299E-02_wp, &
      & -6.03284819663359E-02_wp,  1.02666515671617E-02_wp,  3.01926990960727E-03_wp, &
      & -3.88835431023663E-03_wp, -1.24809100651277E-03_wp,  1.04604503850031E-02_wp, &
      & -4.93172921888944E-03_wp,  1.62477905574884E-02_wp,  4.34278552584213E-02_wp, &
      & -5.45249572759087E-02_wp,  2.26640099821330E-02_wp,  5.63887926599605E-04_wp, &
      & -4.73405294240286E-03_wp, -3.95584904622056E-03_wp, -1.14098611341426E-02_wp, &
      & -3.96103665762880E-03_wp, -1.45346758443479E-02_wp, -7.40139380537715E-02_wp, &
      & -2.93949752468011E-02_wp, -3.95039474886829E-02_wp,  7.21000252707905E-02_wp, &
      &  1.15667054229119E-02_wp, -2.24202832374431E-03_wp,  1.95905390960769E-02_wp, &
      &  1.50621445903036E-03_wp, -1.24039529309376E-02_wp,  4.33699871562573E-21_wp, &
      & -5.18389634661959E-05_wp, -8.25660808937685E-38_wp,  2.32143355452967E-05_wp, &
      &  0.00000000000000E+00_wp,  5.56737885815481E-22_wp,  8.50269723022520E-22_wp, &
      &  3.71719534051772E-22_wp, -3.30898966840816E-02_wp,  3.03839987707241E-05_wp, &
      &  0.00000000000000E+00_wp,  7.84511474868008E-06_wp, -5.62820145027376E-37_wp, &
      & -3.51317067298182E-06_wp, -3.35309640479922E-05_wp,  1.36064515088429E-05_wp, &
      &  3.14258887863276E-02_wp, -3.25134707287603E-02_wp,  2.37706475434611E-02_wp, &
      &  7.38802941269571E-02_wp,  1.10719077913541E-02_wp, -3.52376921367520E-03_wp, &
      &  3.58421095837541E-03_wp,  1.62544335785286E-02_wp,  1.50326271779877E-02_wp, &
      &  3.34787949774353E-02_wp,  3.10853266781575E-02_wp, -2.11308415221338E-02_wp, &
      &  7.88563413039282E-02_wp, -1.27799450512849E-02_wp, -2.45391191981019E-03_wp, &
      &  2.02898523708873E-03_wp, -1.45266116562791E-02_wp,  1.70056225505001E-02_wp, &
      & -8.03768537963252E-02_wp, -5.51476974935269E-02_wp,  5.03924938223464E-02_wp, &
      &  3.22052519482362E-02_wp,  1.66312459498249E-02_wp,  1.08379531091101E-02_wp, &
      & -4.78213571004275E-03_wp, -1.23981006173077E-02_wp, -6.38003057473169E-03_wp, &
      & -2.38433810047813E-21_wp, -1.67528584467143E-21_wp,  8.63705225287944E-38_wp, &
      &  6.20718726692424E-22_wp,  1.36064515088429E-05_wp,  8.17062954224834E-21_wp, &
      & -1.51129341109871E-21_wp, -6.89417666436082E-37_wp,  3.03839987707241E-05_wp, &
      & -3.84756195587905E-01_wp,  1.83599815535358E-23_wp, -4.13489202945117E-21_wp, &
      &  9.14607585545216E-40_wp,  1.54503068998479E-22_wp, -3.66176109146505E-22_wp, &
      & -5.09816450787819E-37_wp,  4.75886999668046E-03_wp,  1.69426506412195E-02_wp, &
      &  4.44433537258549E-03_wp, -2.01386131552921E-02_wp, -3.91810195680666E-03_wp, &
      &  1.14593011472230E-02_wp,  1.84631179330279E-03_wp, -9.97250550671195E-03_wp, &
      & -1.77943546323289E-02_wp, -5.53339867232688E-03_wp,  1.90537965326382E-02_wp, &
      &  4.31200999154014E-03_wp,  2.11524379525669E-02_wp, -4.06800608320505E-03_wp, &
      & -1.08747856125109E-02_wp, -1.00734392764357E-03_wp, -8.54903961597109E-03_wp, &
      &  2.09746823885651E-02_wp,  1.73454412855337E-02_wp, -1.66820784761393E-02_wp, &
      & -1.34280708663993E-02_wp, -1.67778877362497E-02_wp,  7.21583090636907E-03_wp, &
      &  5.83108517600569E-03_wp,  3.07583906118868E-03_wp,  1.50076815099539E-02_wp, &
      &  1.80436610877076E-02_wp, -9.34911860900292E-37_wp, -8.72808114945579E-22_wp, &
      &  5.06814718062176E-22_wp, -1.60782395518071E-21_wp, -3.35309640479922E-05_wp, &
      &  1.11096211355291E-05_wp, -5.03553923197153E-53_wp,  2.48084311112085E-05_wp, &
      &  0.00000000000000E+00_wp,  1.83599815535358E-23_wp, -3.84756195587905E-01_wp, &
      &  1.42215805587042E-23_wp,  3.09006137996958E-22_wp, -2.83638794500689E-22_wp, &
      &  3.23278400199631E-38_wp,  3.66176109146505E-22_wp,  2.54039953642054E-02_wp, &
      &  2.51717246396447E-02_wp,  3.33208185045141E-03_wp,  7.09188019591595E-03_wp, &
      &  1.44219746765395E-02_wp,  1.11604549536128E-02_wp, -9.96507167802176E-03_wp, &
      &  5.82865133966052E-03_wp, -1.35917510010593E-02_wp,  2.29386761989074E-02_wp, &
      & -2.44579306701262E-02_wp,  4.19056227844808E-03_wp,  8.14297903923234E-03_wp, &
      & -1.53458073169640E-02_wp,  4.16324998368276E-03_wp, -1.36975722903759E-02_wp, &
      & -2.79152778454753E-03_wp, -1.30802518555029E-02_wp,  1.23897273967202E-02_wp, &
      & -1.28622273407926E-02_wp,  2.40561204825662E-03_wp, -1.73596479127738E-02_wp, &
      &  8.54681563809014E-03_wp, -6.61655303633574E-04_wp, -7.48641242775246E-03_wp, &
      &  2.50307227825363E-03_wp,  1.69543290975158E-02_wp, -2.97984598215748E-37_wp, &
      & -7.19034672801566E-22_wp, -2.03375430998060E-21_wp, -1.60268886077255E-22_wp, &
      & -3.51317067298182E-06_wp, -4.24136874131500E-05_wp,  2.71762746718437E-05_wp, &
      & -5.34020628137626E-37_wp,  7.84511474868008E-06_wp, -4.13489202945117E-21_wp, &
      &  1.42215805587042E-23_wp, -3.84756195587905E-01_wp, -1.46470443658602E-22_wp, &
      & -2.39355125267618E-22_wp,  2.83638794500688E-22_wp, -1.54503068998479E-22_wp, &
      &  9.53412288597674E-04_wp,  1.04164988814151E-02_wp, -2.28035081904679E-02_wp, &
      &  7.41810044871485E-03_wp,  1.10214491110984E-02_wp, -9.31338835945515E-03_wp, &
      & -1.47600115667260E-02_wp, -1.00433992903600E-02_wp, -2.35752253379770E-03_wp, &
      &  4.81777612959216E-03_wp,  5.05173599101536E-03_wp, -2.42505451150253E-02_wp, &
      & -1.86431864299246E-03_wp,  5.76357979927393E-03_wp,  1.36648707944988E-02_wp, &
      &  9.30321768047876E-03_wp, -1.33552863799000E-02_wp,  2.09347539252052E-03_wp, &
      & -1.58185587295662E-03_wp,  2.66756045848352E-03_wp,  7.70095280590281E-03_wp, &
      & -9.23682340861661E-04_wp, -8.32867220971488E-04_wp, -3.90926292479287E-03_wp, &
      & -2.89890593204321E-03_wp, -7.29801731095905E-03_wp,  1.71447115623448E-03_wp, &
      &  1.85878247589460E-36_wp,  1.24541508040964E-21_wp,  4.37412401264199E-21_wp, &
      &  6.76074258723784E-22_wp,  8.72052034038496E-37_wp, -1.92165281078489E-05_wp, &
      & -4.49865089774473E-05_wp, -8.60547552811564E-06_wp, -5.62820145027376E-37_wp, &
      &  9.14607585545216E-40_wp,  3.09006137996958E-22_wp, -1.46470443658602E-22_wp, &
      & -3.84756195587905E-01_wp,  7.34399262141419E-24_wp,  0.00000000000000E+00_wp, &
      & -2.10683974717548E-37_wp, -1.73810551248513E-02_wp, -1.14387005235194E-02_wp, &
      & -8.32928373190366E-03_wp, -7.78883725745195E-03_wp, -5.88236391647153E-03_wp, &
      & -1.20443103945359E-02_wp,  6.98891079238934E-03_wp, -8.20092460021227E-03_wp, &
      &  2.31684940473338E-03_wp,  1.68697649561765E-02_wp, -1.35968392032966E-02_wp, &
      & -1.12626976409064E-03_wp,  9.30960954316449E-03_wp, -8.81814508092532E-03_wp, &
      &  6.94549201767614E-03_wp, -1.17019824580443E-02_wp, -4.75624462216638E-03_wp, &
      & -3.41988583440891E-03_wp,  1.68918016472164E-02_wp,  4.31789471624345E-03_wp, &
      & -8.94443018558470E-04_wp, -1.59788259256774E-02_wp, -4.82108443058118E-03_wp, &
      & -2.15283464842901E-03_wp, -1.18814193255205E-02_wp,  7.96411287819765E-03_wp, &
      &  8.26736085420374E-03_wp, -2.13549118334625E-21_wp, -1.60268886077255E-22_wp, &
      & -1.10402464140243E-21_wp, -5.56833293700226E-22_wp, -7.84511474868008E-06_wp, &
      & -5.34020628137626E-37_wp,  1.21699802025309E-05_wp, -4.24136874131500E-05_wp, &
      & -3.51317067298182E-06_wp,  1.54503068998479E-22_wp, -2.83638794500689E-22_wp, &
      & -2.39355125267618E-22_wp,  7.34399262141419E-24_wp, -3.84756195587905E-01_wp, &
      &  1.42215805587040E-23_wp,  2.99383189497237E-20_wp,  6.49162929004437E-04_wp, &
      &  7.41807983486012E-03_wp, -1.55264932827591E-02_wp,  4.57254405270234E-03_wp, &
      &  6.05042597742747E-03_wp, -1.00431560089288E-02_wp, -1.00499026180057E-02_wp, &
      & -1.40170440652097E-03_wp, -5.29642017191528E-03_wp, -3.29936174029389E-03_wp, &
      & -1.86381319619824E-03_wp,  1.66076967037159E-02_wp,  3.60713294844584E-03_wp, &
      & -3.84491928366204E-03_wp, -1.33565652750205E-02_wp, -6.37089993897873E-03_wp, &
      &  3.30685791597706E-03_wp, -1.69727547016348E-03_wp,  5.85124250919746E-03_wp, &
      & -9.24289111039500E-04_wp, -2.84849806569867E-02_wp,  5.83388156658344E-03_wp, &
      &  3.57554930201469E-03_wp, -7.29652776304610E-03_wp,  1.07238201385890E-02_wp, &
      &  2.11147832638772E-02_wp, -6.05296949082871E-03_wp,  5.29429411528367E-37_wp, &
      &  1.60782395518071E-21_wp, -7.49692080369944E-22_wp, -8.72808114945580E-22_wp, &
      &  0.00000000000000E+00_wp, -2.48084311112085E-05_wp,  1.07055984358930E-37_wp, &
      &  1.11096211355291E-05_wp, -3.35309640479922E-05_wp, -3.66176109146505E-22_wp, &
      &  3.23278400199631E-38_wp,  2.83638794500688E-22_wp,  0.00000000000000E+00_wp, &
      &  1.42215805587040E-23_wp, -3.84756195587905E-01_wp,  1.83599815535362E-23_wp, &
      & -1.00065792413757E-02_wp,  3.79288917767306E-03_wp, -1.31242521959848E-03_wp, &
      & -2.29259981882272E-02_wp, -1.35918260500578E-02_wp, -1.98059989513990E-03_wp, &
      &  3.92528733235931E-03_wp, -5.84329628223042E-03_wp, -1.47314018915257E-02_wp, &
      &  8.89313866184798E-03_wp,  2.41916650667538E-03_wp,  1.62429225136105E-03_wp, &
      &  2.05358603966318E-02_wp, -1.30797735530431E-02_wp,  1.53175138815014E-03_wp, &
      & -5.31028305259376E-03_wp, -1.20274973483545E-03_wp,  1.33154500801721E-02_wp, &
      & -2.12398836114840E-02_wp, -1.55361327930224E-02_wp, -4.12321264555874E-03_wp, &
      &  1.95992221609897E-02_wp,  1.69537717181884E-02_wp,  1.16661862628272E-03_wp, &
      &  1.28337379049229E-02_wp, -4.28268532092944E-03_wp, -1.06333460602266E-02_wp, &
      &  2.13549118334625E-21_wp, -6.20718726692424E-22_wp,  5.84172583786121E-38_wp, &
      &  3.07445655464504E-22_wp, -3.03839987707241E-05_wp,  6.89417666436082E-37_wp, &
      &  6.28440688161235E-37_wp, -8.17062954224834E-21_wp,  1.36064515088429E-05_wp, &
      & -5.09816450787819E-37_wp,  3.66176109146505E-22_wp, -1.54503068998479E-22_wp, &
      & -2.10683974717548E-37_wp,  2.99383189497237E-20_wp,  1.83599815535362E-23_wp, &
      & -3.84756195587905E-01_wp, -2.10295778274710E-02_wp, -6.41899382251196E-03_wp, &
      & -1.96395717561546E-02_wp, -1.15387649328002E-02_wp, -1.33152873429029E-02_wp, &
      & -1.22254756852146E-02_wp, -8.15804427249692E-03_wp, -1.23469239768497E-02_wp, &
      &  8.68482020562939E-04_wp, -2.35665313676380E-02_wp,  9.50278449584669E-03_wp, &
      &  1.83649724797938E-02_wp, -1.45359391681317E-02_wp,  1.68391261707037E-02_wp, &
      & -1.28074391540540E-02_wp, -4.29381852446258E-03_wp,  1.25249548674218E-02_wp, &
      &  1.23547696834148E-03_wp, -1.71144624585088E-02_wp, -2.30865132440670E-02_wp, &
      &  1.32493653760844E-02_wp,  5.86413536429997E-03_wp,  1.93556808322858E-02_wp, &
      &  1.26758264811348E-02_wp, -3.03621334385227E-03_wp, -9.82695047005563E-03_wp, &
      & -2.36406963669832E-03_wp, -1.22327615420258E-01_wp,  1.27854242498941E-01_wp, &
      &  7.99617118035541E-02_wp,  8.70535260943938E-02_wp, -7.97812857625614E-02_wp, &
      & -7.32828514503475E-02_wp,  2.30448328514989E-02_wp, -4.98969803349299E-02_wp, &
      &  3.14258887863276E-02_wp,  4.75886999668046E-03_wp,  2.54039953642054E-02_wp, &
      &  9.53412288597674E-04_wp, -1.73810551248513E-02_wp,  6.49162929004437E-04_wp, &
      & -1.00065792413757E-02_wp, -2.10295778274710E-02_wp, -3.43061374750722E-01_wp, &
      & -1.65418040536170E-04_wp, -1.03467624553790E-04_wp, -1.12631538215260E-04_wp, &
      & -6.33698647955640E-21_wp,  4.09038452692581E-20_wp,  3.72394203862578E-20_wp, &
      &  2.74507948952301E-20_wp,  0.00000000000000E+00_wp, -6.87735962909471E-04_wp, &
      &  2.60655279888896E-03_wp,  1.46746609513964E-03_wp, -3.89734065804676E-05_wp, &
      &  4.36679701675687E-06_wp, -1.63384313160929E-04_wp,  3.06898178248518E-05_wp, &
      &  2.45843203883275E-06_wp,  1.45073157618217E-04_wp, -7.68777982505079E-04_wp, &
      &  9.56612085875917E-04_wp,  1.62955162476691E-03_wp,  2.70067759387231E-03_wp, &
      & -1.03434910211958E-04_wp, -6.24079275925733E-05_wp,  3.35029740378149E-05_wp, &
      & -1.76243674900006E-04_wp, -1.27739980850912E-04_wp, -5.50799865621865E-02_wp, &
      & -3.00111123699892E-02_wp,  3.44692060710472E-02_wp,  3.75263899004499E-02_wp, &
      & -4.84524422734976E-02_wp, -4.45063128648662E-02_wp, -1.28455561910549E-02_wp, &
      & -6.03284819663359E-02_wp, -3.25134707287603E-02_wp,  1.69426506412195E-02_wp, &
      &  2.51717246396447E-02_wp,  1.04164988814151E-02_wp, -1.14387005235194E-02_wp, &
      &  7.41807983486012E-03_wp,  3.79288917767306E-03_wp, -6.41899382251196E-03_wp, &
      & -1.65418040536170E-04_wp, -3.09507307342230E-01_wp, -3.75164769628872E-20_wp, &
      & -1.08147160715032E-20_wp, -8.75382452393096E-05_wp, -8.04159690619524E-05_wp, &
      &  7.42266851442242E-05_wp, -1.55399181923735E-35_wp,  1.28564389947214E-04_wp, &
      & -2.60412522885625E-03_wp,  7.17085576259614E-03_wp,  4.75433714518235E-03_wp, &
      & -1.26339087306690E-04_wp,  2.40089256053084E-05_wp, -8.98657562999095E-04_wp, &
      &  7.92632213479003E-05_wp,  1.51256860607161E-05_wp,  7.02988952396418E-04_wp, &
      & -9.57205693037985E-04_wp, -3.88768686841018E-04_wp,  1.73451176272957E-03_wp, &
      &  2.87464324655446E-03_wp, -4.61727732809598E-06_wp, -2.78663932274328E-06_wp, &
      &  2.58447201736011E-05_wp, -3.55155811788161E-04_wp, -3.29583774267699E-04_wp, &
      & -3.44795609739750E-02_wp,  3.45024665091326E-02_wp, -6.35473649848168E-02_wp, &
      &  2.34920665793948E-02_wp, -6.03467200158205E-02_wp,  1.50786081940942E-02_wp, &
      &  6.83504280096736E-02_wp,  1.02666515671617E-02_wp,  2.37706475434611E-02_wp, &
      &  4.44433537258549E-03_wp,  3.33208185045141E-03_wp, -2.28035081904679E-02_wp, &
      & -8.32928373190366E-03_wp, -1.55264932827591E-02_wp, -1.31242521959848E-03_wp, &
      & -1.96395717561546E-02_wp, -1.03467624553790E-04_wp, -3.75164769628872E-20_wp, &
      & -3.09507307342230E-01_wp, -2.77176759945251E-20_wp, -1.55399181923735E-35_wp, &
      & -1.28564389947214E-04_wp, -9.28563627701256E-05_wp, -8.75382452393096E-05_wp, &
      &  1.70739219691179E-36_wp, -1.46615810879207E-03_wp,  4.75454046870674E-03_wp, &
      &  1.40282309266397E-03_wp, -7.11307167513486E-05_wp,  1.51263117275510E-05_wp, &
      & -3.76268420816934E-04_wp,  2.29763869112618E-04_wp,  5.65776236896460E-06_wp, &
      &  5.02680137199205E-04_wp, -1.62942594588283E-03_wp,  1.73321821956634E-03_wp, &
      &  1.54556977885277E-03_wp,  4.89331043593652E-03_wp, -3.54896780062363E-04_wp, &
      & -1.41957523900783E-04_wp,  2.56951299406999E-04_wp, -4.00882021496380E-04_wp, &
      & -4.38280925153620E-04_wp, -3.75066680349345E-02_wp,  3.75303108914129E-02_wp, &
      &  2.34718724175513E-02_wp, -5.95719246562874E-02_wp,  4.82889756857431E-03_wp, &
      & -6.03304570737948E-02_wp, -8.74562697366483E-03_wp,  3.01926990960727E-03_wp, &
      &  7.38802941269571E-02_wp, -2.01386131552921E-02_wp,  7.09188019591595E-03_wp, &
      &  7.41810044871485E-03_wp, -7.78883725745195E-03_wp,  4.57254405270234E-03_wp, &
      & -2.29259981882272E-02_wp, -1.15387649328002E-02_wp, -1.12631538215260E-04_wp, &
      & -1.08147160715032E-20_wp, -2.77176759945251E-20_wp, -3.09507307342230E-01_wp, &
      & -1.28564389947214E-04_wp, -1.55399181923735E-35_wp,  5.05402294533029E-05_wp, &
      & -8.04159690619524E-05_wp, -8.75382452393096E-05_wp,  3.82834247139429E-05_wp, &
      & -1.24075412973248E-04_wp, -6.98532480783309E-05_wp, -1.27208097659004E-03_wp, &
      &  1.89463328917381E-04_wp,  1.48611784419206E-05_wp, -1.14136446631704E-06_wp, &
      &  1.06664443475968E-04_wp, -1.60537165162568E-05_wp, -2.70009412242538E-03_wp, &
      &  2.87207255225638E-03_wp,  4.89258324876961E-03_wp,  6.70161058209863E-03_wp, &
      & -5.15947192695021E-04_wp, -3.54845178480825E-04_wp,  7.28145525409465E-05_wp, &
      & -8.79111080065194E-04_wp, -5.22477795870369E-04_wp, -4.35783008427281E-03_wp, &
      & -3.17082096701095E-03_wp,  1.42631249047077E-03_wp, -6.45228663309892E-03_wp, &
      & -8.71650801556639E-03_wp, -1.64973231293152E-02_wp, -9.09347989988447E-03_wp, &
      & -3.88835431023663E-03_wp,  1.10719077913541E-02_wp, -3.91810195680666E-03_wp, &
      &  1.44219746765395E-02_wp,  1.10214491110984E-02_wp, -5.88236391647153E-03_wp, &
      &  6.05042597742747E-03_wp, -1.35918260500578E-02_wp, -1.33152873429029E-02_wp, &
      & -6.33698647955640E-21_wp, -8.75382452393096E-05_wp, -1.55399181923735E-35_wp, &
      & -1.28564389947214E-04_wp, -1.31964145071846E-01_wp,  1.32581538481864E-20_wp, &
      & -4.31104573284133E-21_wp,  1.06440543981631E-20_wp,  0.00000000000000E+00_wp, &
      &  4.24281913647378E-06_wp, -2.33644567493396E-05_wp, -1.47173926684956E-05_wp, &
      & -1.89799716481203E-04_wp,  1.77639161911680E-08_wp,  3.88087767133332E-09_wp, &
      & -4.34972667736503E-10_wp,  1.04997034154379E-08_wp, -3.62481338184538E-09_wp, &
      & -1.03556388698046E-04_wp,  4.82147566947588E-06_wp,  3.55334622119110E-04_wp, &
      &  5.16595300284712E-04_wp, -5.01311685363092E-08_wp, -3.13118404150676E-08_wp, &
      &  1.47267310289273E-08_wp, -1.26636043050946E-07_wp, -9.57912697815438E-08_wp, &
      & -4.00713627845440E-03_wp, -2.91430593493387E-03_wp, -6.69257556262513E-03_wp, &
      &  1.42853729713138E-03_wp, -1.65027194137899E-02_wp, -5.91499794173133E-03_wp, &
      &  1.53629073822107E-02_wp, -1.24809100651277E-03_wp, -3.52376921367520E-03_wp, &
      &  1.14593011472230E-02_wp,  1.11604549536128E-02_wp, -9.31338835945515E-03_wp, &
      & -1.20443103945359E-02_wp, -1.00431560089288E-02_wp, -1.98059989513990E-03_wp, &
      & -1.22254756852146E-02_wp,  4.09038452692581E-20_wp, -8.04159690619524E-05_wp, &
      & -1.28564389947214E-04_wp, -1.55399181923735E-35_wp,  1.32581538481864E-20_wp, &
      & -1.31964145071846E-01_wp,  6.14534767204849E-21_wp,  3.73347512151710E-21_wp, &
      & -1.06440543981631E-20_wp, -1.63546678198079E-04_wp,  9.00251736101473E-04_wp, &
      &  3.76937616207027E-04_wp, -1.49830407249654E-05_wp,  4.58504509571370E-09_wp, &
      & -1.42860566246132E-07_wp,  3.74319341451120E-08_wp,  2.37296555422636E-09_wp, &
      &  1.37454215806030E-07_wp, -6.24964462065621E-05_wp,  2.90896957572468E-06_wp, &
      &  1.42164392328653E-04_wp,  3.55386238415830E-04_wp, -3.13490792846406E-08_wp, &
      & -1.71160683090322E-08_wp,  1.83373639019706E-08_wp, -7.07156954555868E-08_wp, &
      & -6.32983605190282E-08_wp,  1.25536395265239E-03_wp,  3.96576726638335E-03_wp, &
      & -6.19800318861863E-03_wp,  2.70061544651703E-03_wp, -9.10037767448247E-03_wp, &
      &  1.53633897818081E-02_wp,  7.11152645967737E-03_wp,  1.04604503850031E-02_wp, &
      &  3.58421095837541E-03_wp,  1.84631179330279E-03_wp, -9.96507167802176E-03_wp, &
      & -1.47600115667260E-02_wp,  6.98891079238934E-03_wp, -1.00499026180057E-02_wp, &
      &  3.92528733235931E-03_wp, -8.15804427249692E-03_wp,  3.72394203862578E-20_wp, &
      &  7.42266851442242E-05_wp, -9.28563627701256E-05_wp,  5.05402294533029E-05_wp, &
      & -4.31104573284133E-21_wp,  6.14534767204849E-21_wp, -1.31964145071846E-01_wp, &
      &  7.65459869320790E-21_wp,  0.00000000000000E+00_wp,  3.07095858900258E-05_wp, &
      & -7.93446795212911E-05_wp, -2.30137055242606E-04_wp,  1.21031707400086E-06_wp, &
      & -5.71152880407719E-10_wp,  3.73787245733402E-08_wp,  7.41722882527973E-09_wp, &
      & -6.18427735470074E-10_wp, -1.69383698561871E-08_wp,  3.34943007108530E-05_wp, &
      & -2.57530509938010E-05_wp, -2.56966254793346E-04_wp, -7.27937944021209E-05_wp, &
      &  1.45400090914039E-08_wp,  1.81378783229037E-08_wp,  1.23615015110711E-08_wp, &
      &  5.14621015767376E-08_wp,  1.80914543507347E-08_wp, -2.72874000813000E-03_wp, &
      &  1.42880749622715E-03_wp, -4.55736928349791E-03_wp, -4.03948280660208E-03_wp, &
      & -3.88958831715503E-03_wp, -1.24811473959793E-03_wp,  1.04610161044146E-02_wp, &
      & -4.93172921888944E-03_wp,  1.62544335785286E-02_wp, -9.97250550671195E-03_wp, &
      &  5.82865133966052E-03_wp, -1.00433992903600E-02_wp, -8.20092460021227E-03_wp, &
      & -1.40170440652097E-03_wp, -5.84329628223042E-03_wp, -1.23469239768497E-02_wp, &
      &  2.74507948952301E-20_wp, -1.55399181923735E-35_wp, -8.75382452393096E-05_wp, &
      & -8.04159690619524E-05_wp,  1.06440543981631E-20_wp,  3.73347512151710E-21_wp, &
      &  7.65459869320790E-21_wp, -1.31964145071846E-01_wp,  1.32581538481864E-20_wp, &
      &  2.38880259091061E-06_wp, -1.47180195500970E-05_wp, -5.50956708759523E-06_wp, &
      & -1.06861567997999E-04_wp,  1.05043808765199E-08_wp,  2.01321765740406E-09_wp, &
      & -4.90759822628855E-10_wp,  5.02295446556431E-09_wp, -2.18318090214474E-09_wp, &
      & -1.76235000855311E-04_wp,  3.55075516525441E-04_wp,  4.00906656710059E-04_wp, &
      &  8.79171512391947E-04_wp, -1.25836528205761E-07_wp, -7.01864955584681E-08_wp, &
      &  5.14706096744002E-08_wp, -1.90963189386393E-07_wp, -1.47396387488401E-07_wp, &
      &  1.71594545101538E-03_wp,  7.10720144856875E-03_wp, -5.61513696179550E-04_wp, &
      & -6.06232790463037E-03_wp,  1.10719455673932E-02_wp, -3.52293118440440E-03_wp, &
      &  3.58184950097580E-03_wp,  1.62477905574884E-02_wp,  1.50326271779877E-02_wp, &
      & -1.77943546323289E-02_wp, -1.35917510010593E-02_wp, -2.35752253379770E-03_wp, &
      &  2.31684940473338E-03_wp, -5.29642017191528E-03_wp, -1.47314018915257E-02_wp, &
      &  8.68482020562939E-04_wp,  0.00000000000000E+00_wp,  1.28564389947214E-04_wp, &
      &  1.70739219691179E-36_wp, -8.75382452393096E-05_wp,  0.00000000000000E+00_wp, &
      & -1.06440543981631E-20_wp,  0.00000000000000E+00_wp,  1.32581538481864E-20_wp, &
      & -1.31964145071846E-01_wp,  1.45209014418172E-04_wp, -7.04197132465522E-04_wp, &
      & -5.03537480938425E-04_wp,  1.60797268274272E-05_wp, -4.17030606962770E-09_wp, &
      &  1.37411394974633E-07_wp, -1.69839222959105E-08_wp, -2.51883908039411E-09_wp, &
      & -1.13529689571677E-07_wp, -1.27648799043105E-04_wp,  3.29466739037622E-04_wp, &
      &  4.38012532774768E-04_wp,  5.22073157288066E-04_wp, -9.50086216996425E-08_wp, &
      & -6.26984303999139E-08_wp,  1.80182645728791E-08_wp, -1.46855880598522E-07_wp, &
      & -9.05158124333577E-08_wp, -1.22264710871624E-01_wp, -1.32703366657036E-01_wp, &
      & -6.67360487355299E-02_wp,  9.08810651002315E-02_wp,  8.63568720380417E-02_wp, &
      & -6.34132609908191E-02_wp,  3.50622164993731E-02_wp,  4.34278552584213E-02_wp, &
      &  3.34787949774353E-02_wp, -5.53339867232688E-03_wp,  2.29386761989074E-02_wp, &
      &  4.81777612959216E-03_wp,  1.68697649561765E-02_wp, -3.29936174029389E-03_wp, &
      &  8.89313866184798E-03_wp, -2.35665313676380E-02_wp, -6.87735962909471E-04_wp, &
      & -2.60412522885625E-03_wp, -1.46615810879207E-03_wp,  3.82834247139429E-05_wp, &
      &  4.24281913647378E-06_wp, -1.63546678198079E-04_wp,  3.07095858900258E-05_wp, &
      &  2.38880259091061E-06_wp,  1.45209014418172E-04_wp, -3.43064588788769E-01_wp, &
      &  1.69580169139235E-04_wp,  8.52691954660391E-05_wp, -1.16128541944409E-04_wp, &
      & -8.91013369739260E-21_wp,  1.65166955809329E-20_wp, -1.22035636042108E-19_wp, &
      &  9.27267678358380E-21_wp, -1.22035636042108E-19_wp, -6.25422081420641E-04_wp, &
      & -1.59131822740711E-03_wp, -4.26826904290270E-06_wp,  2.24242727282653E-03_wp, &
      &  1.61689839914723E-04_wp, -2.95314950354044E-07_wp,  9.89001535888966E-05_wp, &
      &  4.16203740883160E-07_wp, -5.65718702045342E-05_wp,  5.72416997864790E-02_wp, &
      & -2.56636215505385E-02_wp,  2.99090550798106E-02_wp, -4.07300329452995E-02_wp, &
      & -5.83479881335090E-02_wp,  4.28456294498154E-02_wp, -1.80614896024256E-03_wp, &
      & -5.45249572759087E-02_wp,  3.10853266781575E-02_wp,  1.90537965326382E-02_wp, &
      & -2.44579306701262E-02_wp,  5.05173599101536E-03_wp, -1.35968392032966E-02_wp, &
      & -1.86381319619824E-03_wp,  2.41916650667538E-03_wp,  9.50278449584669E-03_wp, &
      &  2.60655279888896E-03_wp,  7.17085576259614E-03_wp,  4.75454046870674E-03_wp, &
      & -1.24075412973248E-04_wp, -2.33644567493396E-05_wp,  9.00251736101473E-04_wp, &
      & -7.93446795212911E-05_wp, -1.47180195500970E-05_wp, -7.04197132465522E-04_wp, &
      &  1.69580169139235E-04_wp, -3.09502658142368E-01_wp, -9.84453710223430E-21_wp, &
      &  7.02454996156433E-21_wp, -9.02476461398641E-05_wp,  6.62657435476439E-05_wp, &
      & -7.60871552207326E-05_wp, -4.49799112719485E-36_wp, -1.31786818645688E-04_wp, &
      &  1.59292180207359E-03_wp,  2.30096460217580E-03_wp,  9.27681631603521E-06_wp, &
      & -4.89108303800000E-03_wp, -5.12483329906770E-04_wp,  9.39308604426324E-07_wp, &
      & -3.45586942456063E-04_wp, -1.73184270840948E-06_wp,  3.47756644904370E-04_wp, &
      &  2.87567955170003E-02_wp,  2.98766519495830E-02_wp, -7.01120822706820E-02_wp, &
      & -2.04608707777862E-02_wp, -5.45064290098753E-02_wp, -3.30931552610038E-02_wp, &
      & -6.45901676781609E-02_wp,  2.26640099821330E-02_wp, -2.11308415221338E-02_wp, &
      &  4.31200999154014E-03_wp,  4.19056227844808E-03_wp, -2.42505451150253E-02_wp, &
      & -1.12626976409064E-03_wp,  1.66076967037159E-02_wp,  1.62429225136105E-03_wp, &
      &  1.83649724797938E-02_wp,  1.46746609513964E-03_wp,  4.75433714518235E-03_wp, &
      &  1.40282309266397E-03_wp, -6.98532480783309E-05_wp, -1.47173926684956E-05_wp, &
      &  3.76937616207027E-04_wp, -2.30137055242606E-04_wp, -5.50956708759523E-06_wp, &
      & -5.03537480938425E-04_wp,  8.52691954660391E-05_wp, -9.84453710223430E-21_wp, &
      & -3.09502658142368E-01_wp,  7.29177902324088E-21_wp, -4.49799112719485E-36_wp, &
      &  1.31786818645688E-04_wp,  7.65170897505659E-05_wp, -9.02476461398641E-05_wp, &
      &  0.00000000000000E+00_wp,  4.72737579798131E-06_wp,  1.03335124097745E-05_wp, &
      & -1.16988676260113E-03_wp, -1.45616318049003E-05_wp, -1.93313870347339E-06_wp, &
      & -1.12817914398629E-04_wp, -1.51673759940898E-06_wp,  1.58997830454910E-04_wp, &
      &  6.76344244920974E-07_wp, -3.91833859030749E-02_wp, -4.07102422249159E-02_wp, &
      & -2.04731091749614E-02_wp, -5.72568009065204E-02_wp,  1.12280142171848E-03_wp, &
      & -5.45166480109714E-02_wp,  1.23206159814852E-03_wp,  5.63887926599605E-04_wp, &
      &  7.88563413039282E-02_wp,  2.11524379525669E-02_wp,  8.14297903923234E-03_wp, &
      & -1.86431864299246E-03_wp,  9.30960954316449E-03_wp,  3.60713294844584E-03_wp, &
      &  2.05358603966318E-02_wp, -1.45359391681317E-02_wp, -3.89734065804676E-05_wp, &
      & -1.26339087306690E-04_wp, -7.11307167513486E-05_wp, -1.27208097659004E-03_wp, &
      & -1.89799716481203E-04_wp, -1.49830407249654E-05_wp,  1.21031707400086E-06_wp, &
      & -1.06861567997999E-04_wp,  1.60797268274272E-05_wp, -1.16128541944409E-04_wp, &
      &  7.02454996156433E-21_wp,  7.29177902324088E-21_wp, -3.09502658142368E-01_wp, &
      &  1.31786818645688E-04_wp, -4.49799112719485E-36_wp,  5.21045027925806E-05_wp, &
      &  6.62657435476439E-05_wp, -9.02476461398641E-05_wp, -2.24402827110634E-03_wp, &
      & -4.88955261261686E-03_wp, -1.30684694907557E-05_wp,  5.72033091687622E-03_wp, &
      &  8.33220930307883E-04_wp, -1.73125816033888E-06_wp,  4.86860503191426E-04_wp, &
      &  2.15042503779012E-06_wp, -1.71992315683826E-04_wp,  4.74237605027460E-03_wp, &
      & -3.10850629609557E-03_wp,  1.29921644154981E-03_wp,  6.54751512891106E-03_wp, &
      & -1.09947461474657E-02_wp,  1.60940501411438E-02_wp,  5.23416140050281E-03_wp, &
      & -4.73405294240286E-03_wp, -1.27799450512849E-02_wp, -4.06800608320505E-03_wp, &
      & -1.53458073169640E-02_wp,  5.76357979927393E-03_wp, -8.81814508092532E-03_wp, &
      & -3.84491928366204E-03_wp, -1.30797735530431E-02_wp,  1.68391261707037E-02_wp, &
      &  4.36679701675687E-06_wp,  2.40089256053084E-05_wp,  1.51263117275510E-05_wp, &
      &  1.89463328917381E-04_wp,  1.77639161911680E-08_wp,  4.58504509571370E-09_wp, &
      & -5.71152880407719E-10_wp,  1.05043808765199E-08_wp, -4.17030606962770E-09_wp, &
      & -8.91013369739260E-21_wp, -9.02476461398641E-05_wp, -4.49799112719485E-36_wp, &
      &  1.31786818645688E-04_wp, -1.31960989638649E-01_wp, -1.19717926066414E-22_wp, &
      &  1.65153134013025E-20_wp,  1.71604262059467E-20_wp, -9.85556005895367E-36_wp, &
      &  1.61576215080994E-04_wp,  5.11776201333169E-04_wp,  1.85805755747294E-06_wp, &
      & -8.32076288834617E-04_wp, -1.27614142239949E-07_wp,  1.26323225508708E-10_wp, &
      & -7.88700816457376E-08_wp, -1.85025274120340E-10_wp,  5.06332576027045E-08_wp, &
      & -3.47983485181397E-03_wp,  2.28176804805362E-03_wp,  7.36400805627289E-03_wp, &
      &  1.29781635285664E-03_wp,  1.60903803016326E-02_wp, -8.92957506230816E-04_wp, &
      &  1.66618117420543E-02_wp, -3.95584904622056E-03_wp, -2.45391191981019E-03_wp, &
      & -1.08747856125109E-02_wp,  4.16324998368276E-03_wp,  1.36648707944988E-02_wp, &
      &  6.94549201767614E-03_wp, -1.33565652750205E-02_wp,  1.53175138815014E-03_wp, &
      & -1.28074391540540E-02_wp, -1.63384313160929E-04_wp, -8.98657562999095E-04_wp, &
      & -3.76268420816934E-04_wp,  1.48611784419206E-05_wp,  3.88087767133332E-09_wp, &
      & -1.42860566246132E-07_wp,  3.73787245733402E-08_wp,  2.01321765740406E-09_wp, &
      &  1.37411394974633E-07_wp,  1.65166955809329E-20_wp,  6.62657435476439E-05_wp, &
      &  1.31786818645688E-04_wp, -4.49799112719485E-36_wp, -1.19717926066414E-22_wp, &
      & -1.31960989638649E-01_wp,  9.90757668941204E-21_wp, -1.43026809569896E-20_wp, &
      & -1.71604262059467E-20_wp, -3.52819224301537E-07_wp, -1.11422599312532E-06_wp, &
      &  1.12681425808760E-04_wp,  2.05960363970976E-06_wp,  4.06121639039401E-10_wp, &
      &  5.16194371281006E-09_wp,  2.63497629891793E-10_wp, -8.14788659713529E-09_wp, &
      & -1.81352835022564E-10_wp,  1.92821272181122E-03_wp, -3.75035959129963E-03_wp, &
      &  5.35181103965646E-03_wp,  2.56648577780100E-03_wp,  5.22815186685464E-03_wp, &
      &  1.66594585688245E-02_wp,  1.40607882463047E-03_wp, -1.14098611341426E-02_wp, &
      &  2.02898523708873E-03_wp, -1.00734392764357E-03_wp, -1.36975722903759E-02_wp, &
      &  9.30321768047876E-03_wp, -1.17019824580443E-02_wp, -6.37089993897873E-03_wp, &
      & -5.31028305259376E-03_wp, -4.29381852446258E-03_wp,  3.06898178248518E-05_wp, &
      &  7.92632213479003E-05_wp,  2.29763869112618E-04_wp, -1.14136446631704E-06_wp, &
      & -4.34972667736503E-10_wp,  3.74319341451120E-08_wp,  7.41722882527973E-09_wp, &
      & -4.90759822628855E-10_wp, -1.69839222959105E-08_wp, -1.22035636042108E-19_wp, &
      & -7.60871552207326E-05_wp,  7.65170897505659E-05_wp,  5.21045027925806E-05_wp, &
      &  1.65153134013025E-20_wp,  9.90757668941204E-21_wp, -1.31960989638649E-01_wp, &
      & -6.91191768412672E-23_wp,  3.24076230367709E-20_wp,  9.88136428776565E-05_wp, &
      &  3.44994863566807E-04_wp,  1.53756540080326E-06_wp, -4.86139135782912E-04_wp, &
      & -7.87794974228051E-08_wp,  1.09096231782629E-10_wp, -4.88212018653482E-08_wp, &
      & -1.54164306323999E-10_wp,  2.76855400590531E-08_wp,  2.38169502570973E-03_wp, &
      &  1.29667791400357E-03_wp, -5.04039164139594E-03_wp,  3.28899897889399E-03_wp, &
      & -4.73227358642259E-03_wp, -3.95496295995601E-03_wp, -1.14074143898209E-02_wp, &
      & -3.96103665762880E-03_wp, -1.45266116562791E-02_wp, -8.54903961597109E-03_wp, &
      & -2.79152778454753E-03_wp, -1.33552863799000E-02_wp, -4.75624462216638E-03_wp, &
      &  3.30685791597706E-03_wp, -1.20274973483545E-03_wp,  1.25249548674218E-02_wp, &
      &  2.45843203883275E-06_wp,  1.51256860607161E-05_wp,  5.65776236896460E-06_wp, &
      &  1.06664443475968E-04_wp,  1.04997034154379E-08_wp,  2.37296555422636E-09_wp, &
      & -6.18427735470074E-10_wp,  5.02295446556431E-09_wp, -2.51883908039411E-09_wp, &
      &  9.27267678358380E-21_wp, -4.49799112719485E-36_wp, -9.02476461398641E-05_wp, &
      &  6.62657435476439E-05_wp,  1.71604262059467E-20_wp, -1.43026809569896E-20_wp, &
      & -6.91191768412672E-23_wp, -1.31960989638649E-01_wp, -1.19717926066415E-22_wp, &
      &  4.96989235046967E-07_wp,  2.05902006007559E-06_wp, -1.58721004179660E-04_wp, &
      & -2.55372433687406E-06_wp, -5.97302624741964E-10_wp, -8.11966989908650E-09_wp, &
      & -3.70804240816506E-10_wp,  1.08491741544299E-08_wp,  1.81846130288217E-10_wp, &
      &  1.84157392408560E-03_wp, -7.31261979409684E-03_wp,  5.04989555182662E-04_wp, &
      & -6.38013719371571E-03_wp, -1.27802615861412E-02_wp, -2.45202972143322E-03_wp, &
      &  2.02859278523153E-03_wp, -1.45346758443479E-02_wp,  1.70056225505001E-02_wp, &
      &  2.09746823885651E-02_wp, -1.30802518555029E-02_wp,  2.09347539252052E-03_wp, &
      & -3.41988583440891E-03_wp, -1.69727547016348E-03_wp,  1.33154500801721E-02_wp, &
      &  1.23547696834148E-03_wp,  1.45073157618217E-04_wp,  7.02988952396418E-04_wp, &
      &  5.02680137199205E-04_wp, -1.60537165162568E-05_wp, -3.62481338184538E-09_wp, &
      &  1.37454215806030E-07_wp, -1.69383698561871E-08_wp, -2.18318090214474E-09_wp, &
      & -1.13529689571677E-07_wp, -1.22035636042108E-19_wp, -1.31786818645688E-04_wp, &
      &  0.00000000000000E+00_wp, -9.02476461398641E-05_wp, -9.85556005895367E-36_wp, &
      & -1.71604262059467E-20_wp,  3.24076230367709E-20_wp, -1.19717926066415E-22_wp, &
      & -1.31960989638649E-01_wp, -5.64438261527678E-05_wp, -3.46917623300954E-04_wp, &
      & -6.49080416742357E-07_wp,  1.71314470754176E-04_wp,  5.02195332420630E-08_wp, &
      & -5.34814212462765E-11_wp,  2.74296510512656E-08_wp,  5.71037372972422E-11_wp, &
      & -1.00558672663130E-09_wp, -1.22460964229986E-01_wp,  4.20157391316526E-02_wp, &
      & -6.63241385241726E-02_wp, -1.55408826901311E-01_wp,  4.68870629578708E-02_wp, &
      &  2.00098126116035E-02_wp,  3.54875371122529E-02_wp, -7.40139380537715E-02_wp, &
      & -8.03768537963252E-02_wp,  1.73454412855337E-02_wp,  1.23897273967202E-02_wp, &
      & -1.58185587295662E-03_wp,  1.68918016472164E-02_wp,  5.85124250919746E-03_wp, &
      & -2.12398836114840E-02_wp, -1.71144624585088E-02_wp, -7.68777982505079E-04_wp, &
      & -9.57205693037985E-04_wp, -1.62942594588283E-03_wp, -2.70009412242538E-03_wp, &
      & -1.03556388698046E-04_wp, -6.24964462065621E-05_wp,  3.34943007108530E-05_wp, &
      & -1.76235000855311E-04_wp, -1.27648799043105E-04_wp, -6.25422081420641E-04_wp, &
      &  1.59292180207359E-03_wp,  4.72737579798131E-06_wp, -2.24402827110634E-03_wp, &
      &  1.61576215080994E-04_wp, -3.52819224301537E-07_wp,  9.88136428776565E-05_wp, &
      &  4.96989235046967E-07_wp, -5.64438261527678E-05_wp, -3.43055500606573E-01_wp, &
      & -5.56044630948233E-05_wp,  8.77762122515226E-05_wp,  2.05699096783231E-04_wp, &
      & -2.29383246087960E-21_wp,  1.92672553066413E-21_wp,  0.00000000000000E+00_wp, &
      &  1.18124164857075E-20_wp,  0.00000000000000E+00_wp, -1.80789789843406E-02_wp, &
      & -7.92556852535859E-02_wp, -9.37803943936767E-03_wp, -2.19742878982292E-02_wp, &
      & -6.72860880385563E-02_wp, -2.87162310555367E-02_wp,  6.84419358448112E-04_wp, &
      & -2.93949752468011E-02_wp, -5.51476974935269E-02_wp, -1.66820784761393E-02_wp, &
      & -1.28622273407926E-02_wp,  2.66756045848352E-03_wp,  4.31789471624345E-03_wp, &
      & -9.24289111039500E-04_wp, -1.55361327930224E-02_wp, -2.30865132440670E-02_wp, &
      &  9.56612085875917E-04_wp, -3.88768686841018E-04_wp,  1.73321821956634E-03_wp, &
      &  2.87207255225638E-03_wp,  4.82147566947588E-06_wp,  2.90896957572468E-06_wp, &
      & -2.57530509938010E-05_wp,  3.55075516525441E-04_wp,  3.29466739037622E-04_wp, &
      & -1.59131822740711E-03_wp,  2.30096460217580E-03_wp,  1.03335124097745E-05_wp, &
      & -4.88955261261686E-03_wp,  5.11776201333169E-04_wp, -1.11422599312532E-06_wp, &
      &  3.44994863566807E-04_wp,  2.05902006007559E-06_wp, -3.46917623300954E-04_wp, &
      & -5.56044630948233E-05_wp, -3.09515804463225E-01_wp,  1.09680426008225E-21_wp, &
      &  1.46438047844144E-20_wp,  1.59868971108584E-04_wp,  6.82195155929510E-05_wp, &
      &  2.49505899856878E-05_wp, -6.32347368154228E-37_wp,  4.32156895340305E-05_wp, &
      &  2.85423127657356E-02_wp, -9.37929446104940E-03_wp, -7.03908687281913E-02_wp, &
      &  3.46922406378281E-02_wp, -2.93962056981400E-02_wp,  1.06805031152813E-02_wp, &
      & -6.45840327694736E-02_wp, -3.95039474886829E-02_wp,  5.03924938223464E-02_wp, &
      & -1.34280708663993E-02_wp,  2.40561204825662E-03_wp,  7.70095280590281E-03_wp, &
      & -8.94443018558470E-04_wp, -2.84849806569867E-02_wp, -4.12321264555874E-03_wp, &
      &  1.32493653760844E-02_wp,  1.62955162476691E-03_wp,  1.73451176272957E-03_wp, &
      &  1.54556977885277E-03_wp,  4.89258324876961E-03_wp,  3.55334622119110E-04_wp, &
      &  1.42164392328653E-04_wp, -2.56966254793346E-04_wp,  4.00906656710059E-04_wp, &
      &  4.38012532774768E-04_wp, -4.26826904290270E-06_wp,  9.27681631603521E-06_wp, &
      & -1.16988676260113E-03_wp, -1.30684694907557E-05_wp,  1.85805755747294E-06_wp, &
      &  1.12681425808760E-04_wp,  1.53756540080326E-06_wp, -1.58721004179660E-04_wp, &
      & -6.49080416742357E-07_wp,  8.77762122515226E-05_wp,  1.09680426008225E-21_wp, &
      & -3.09515804463225E-01_wp, -2.40888451735672E-20_wp, -6.32347368154228E-37_wp, &
      & -4.32156895340305E-05_wp,  7.87731113831522E-05_wp,  1.59868971108584E-04_wp, &
      & -1.60941952528841E-36_wp,  6.69377814947157E-02_wp, -2.19972176845851E-02_wp, &
      &  3.47237944532958E-02_wp, -3.83302369012886E-03_wp, -4.56745992964253E-02_wp, &
      & -2.94045683181393E-02_wp, -2.54957147564150E-03_wp,  7.21000252707905E-02_wp, &
      &  3.22052519482362E-02_wp, -1.67778877362497E-02_wp, -1.73596479127738E-02_wp, &
      & -9.23682340861661E-04_wp, -1.59788259256774E-02_wp,  5.83388156658344E-03_wp, &
      &  1.95992221609897E-02_wp,  5.86413536429997E-03_wp,  2.70067759387231E-03_wp, &
      &  2.87464324655446E-03_wp,  4.89331043593652E-03_wp,  6.70161058209863E-03_wp, &
      &  5.16595300284712E-04_wp,  3.55386238415830E-04_wp, -7.27937944021209E-05_wp, &
      &  8.79171512391947E-04_wp,  5.22073157288066E-04_wp,  2.24242727282653E-03_wp, &
      & -4.89108303800000E-03_wp, -1.45616318049003E-05_wp,  5.72033091687622E-03_wp, &
      & -8.32076288834617E-04_wp,  2.05960363970976E-06_wp, -4.86139135782912E-04_wp, &
      & -2.55372433687406E-06_wp,  1.71314470754176E-04_wp,  2.05699096783231E-04_wp, &
      &  1.46438047844144E-20_wp, -2.40888451735672E-20_wp, -3.09515804463225E-01_wp, &
      & -4.32156895340305E-05_wp, -6.32347368154228E-37_wp, -9.23003935046092E-05_wp, &
      &  6.82195155929510E-05_wp,  1.59868971108584E-04_wp,  2.53842716456936E-03_wp, &
      &  9.29502413478776E-03_wp,  6.85386149539902E-04_wp, -1.02076481691057E-03_wp, &
      &  1.24385119738042E-02_wp,  9.71665080564901E-03_wp,  2.79196298847885E-03_wp, &
      &  1.15667054229119E-02_wp,  1.66312459498249E-02_wp,  7.21583090636907E-03_wp, &
      &  8.54681563809014E-03_wp, -8.32867220971488E-04_wp, -4.82108443058118E-03_wp, &
      &  3.57554930201469E-03_wp,  1.69537717181884E-02_wp,  1.93556808322858E-02_wp, &
      & -1.03434910211958E-04_wp, -4.61727732809598E-06_wp, -3.54896780062363E-04_wp, &
      & -5.15947192695021E-04_wp, -5.01311685363092E-08_wp, -3.13490792846406E-08_wp, &
      &  1.45400090914039E-08_wp, -1.25836528205761E-07_wp, -9.50086216996425E-08_wp, &
      &  1.61689839914723E-04_wp, -5.12483329906770E-04_wp, -1.93313870347339E-06_wp, &
      &  8.33220930307883E-04_wp, -1.27614142239949E-07_wp,  4.06121639039401E-10_wp, &
      & -7.87794974228051E-08_wp, -5.97302624741964E-10_wp,  5.02195332420630E-08_wp, &
      & -2.29383246087960E-21_wp,  1.59868971108584E-04_wp, -6.32347368154228E-37_wp, &
      & -4.32156895340305E-05_wp, -1.31969912106704E-01_wp,  4.14751484954949E-21_wp, &
      &  0.00000000000000E+00_wp, -1.12115386014571E-21_wp,  0.00000000000000E+00_wp, &
      &  1.08208897724546E-03_wp,  3.96252195216349E-03_wp, -2.33498037015163E-03_wp, &
      &  6.84246145643275E-04_wp,  9.71299399288886E-03_wp, -6.18518502445096E-03_wp, &
      & -5.29353832250156E-03_wp, -2.24202832374431E-03_wp,  1.08379531091101E-02_wp, &
      &  5.83108517600569E-03_wp, -6.61655303633574E-04_wp, -3.90926292479287E-03_wp, &
      & -2.15283464842901E-03_wp, -7.29652776304610E-03_wp,  1.16661862628272E-03_wp, &
      &  1.26758264811348E-02_wp, -6.24079275925733E-05_wp, -2.78663932274328E-06_wp, &
      & -1.41957523900783E-04_wp, -3.54845178480825E-04_wp, -3.13118404150676E-08_wp, &
      & -1.71160683090322E-08_wp,  1.81378783229037E-08_wp, -7.01864955584681E-08_wp, &
      & -6.26984303999139E-08_wp, -2.95314950354044E-07_wp,  9.39308604426324E-07_wp, &
      & -1.12817914398629E-04_wp, -1.73125816033888E-06_wp,  1.26323225508708E-10_wp, &
      &  5.16194371281006E-09_wp,  1.09096231782629E-10_wp, -8.11966989908650E-09_wp, &
      & -5.34814212462765E-11_wp,  1.92672553066413E-21_wp,  6.82195155929510E-05_wp, &
      & -4.32156895340305E-05_wp, -6.32347368154228E-37_wp,  4.14751484954949E-21_wp, &
      & -1.31969912106704E-01_wp, -6.47298482958116E-22_wp,  0.00000000000000E+00_wp, &
      &  1.12115386014571E-21_wp,  1.92548114105724E-03_wp,  1.18692149774725E-03_wp, &
      &  5.30964575093812E-03_wp, -4.39726862536289E-03_wp,  2.78692490239979E-03_wp, &
      & -5.29715029859336E-03_wp,  1.23195416282767E-03_wp,  1.95905390960769E-02_wp, &
      & -4.78213571004275E-03_wp,  3.07583906118868E-03_wp, -7.48641242775246E-03_wp, &
      & -2.89890593204321E-03_wp, -1.18814193255205E-02_wp,  1.07238201385890E-02_wp, &
      &  1.28337379049229E-02_wp, -3.03621334385227E-03_wp,  3.35029740378149E-05_wp, &
      &  2.58447201736011E-05_wp,  2.56951299406999E-04_wp,  7.28145525409465E-05_wp, &
      &  1.47267310289273E-08_wp,  1.83373639019706E-08_wp,  1.23615015110711E-08_wp, &
      &  5.14706096744002E-08_wp,  1.80182645728791E-08_wp,  9.89001535888966E-05_wp, &
      & -3.45586942456063E-04_wp, -1.51673759940898E-06_wp,  4.86860503191426E-04_wp, &
      & -7.88700816457376E-08_wp,  2.63497629891793E-10_wp, -4.88212018653482E-08_wp, &
      & -3.70804240816506E-10_wp,  2.74296510512656E-08_wp,  0.00000000000000E+00_wp, &
      &  2.49505899856878E-05_wp,  7.87731113831522E-05_wp, -9.23003935046092E-05_wp, &
      &  0.00000000000000E+00_wp, -6.47298482958116E-22_wp, -1.31969912106704E-01_wp, &
      &  2.39456881485537E-21_wp,  0.00000000000000E+00_wp, -4.00772120793916E-03_wp, &
      &  6.85553886364135E-04_wp,  8.64699576834706E-03_wp,  1.61134893933712E-03_wp, &
      &  1.15672992586680E-02_wp, -2.24309708993980E-03_wp,  1.95920615572102E-02_wp, &
      &  1.50621445903036E-03_wp, -1.23981006173077E-02_wp,  1.50076815099539E-02_wp, &
      &  2.50307227825363E-03_wp, -7.29801731095905E-03_wp,  7.96411287819765E-03_wp, &
      &  2.11147832638772E-02_wp, -4.28268532092944E-03_wp, -9.82695047005563E-03_wp, &
      & -1.76243674900006E-04_wp, -3.55155811788161E-04_wp, -4.00882021496380E-04_wp, &
      & -8.79111080065194E-04_wp, -1.26636043050946E-07_wp, -7.07156954555868E-08_wp, &
      &  5.14621015767376E-08_wp, -1.90963189386393E-07_wp, -1.46855880598522E-07_wp, &
      &  4.16203740883160E-07_wp, -1.73184270840948E-06_wp,  1.58997830454910E-04_wp, &
      &  2.15042503779012E-06_wp, -1.85025274120340E-10_wp, -8.14788659713529E-09_wp, &
      & -1.54164306323999E-10_wp,  1.08491741544299E-08_wp,  5.71037372972422E-11_wp, &
      &  1.18124164857075E-20_wp, -6.32347368154228E-37_wp,  1.59868971108584E-04_wp, &
      &  6.82195155929510E-05_wp, -1.12115386014571E-21_wp,  0.00000000000000E+00_wp, &
      &  2.39456881485537E-21_wp, -1.31969912106704E-01_wp,  4.14751484954949E-21_wp, &
      & -4.35807794043936E-03_wp,  3.37301693901983E-03_wp, -1.17771780132484E-03_wp, &
      &  6.96923546958978E-03_wp,  1.66316199240656E-02_wp,  1.08435978400652E-02_wp, &
      & -4.78485486891422E-03_wp, -1.24039529309376E-02_wp, -6.38003057473169E-03_wp, &
      &  1.80436610877076E-02_wp,  1.69543290975158E-02_wp,  1.71447115623448E-03_wp, &
      &  8.26736085420374E-03_wp, -6.05296949082871E-03_wp, -1.06333460602266E-02_wp, &
      & -2.36406963669832E-03_wp, -1.27739980850912E-04_wp, -3.29583774267699E-04_wp, &
      & -4.38280925153620E-04_wp, -5.22477795870369E-04_wp, -9.57912697815438E-08_wp, &
      & -6.32983605190282E-08_wp,  1.80914543507347E-08_wp, -1.47396387488401E-07_wp, &
      & -9.05158124333577E-08_wp, -5.65718702045342E-05_wp,  3.47756644904370E-04_wp, &
      &  6.76344244920974E-07_wp, -1.71992315683826E-04_wp,  5.06332576027045E-08_wp, &
      & -1.81352835022564E-10_wp,  2.76855400590531E-08_wp,  1.81846130288217E-10_wp, &
      & -1.00558672663130E-09_wp,  0.00000000000000E+00_wp,  4.32156895340305E-05_wp, &
      & -1.60941952528841E-36_wp,  1.59868971108584E-04_wp,  0.00000000000000E+00_wp, &
      &  1.12115386014571E-21_wp,  0.00000000000000E+00_wp,  4.14751484954949E-21_wp, &
      & -1.31969912106704E-01_wp], shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "CeCl3")
   call test_hamiltonian_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & make_gxtb_ncoord, hamiltonian, thr_in=thr1)

end subroutine test_hamiltonian_gxtb_cecl3


subroutine test_hamiltonian_gxtb_ce2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 32
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -2.87625642334199E-01_wp,  0.00000000000000E+00_wp, -1.06392206676407E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.51010910312432E-19_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.23075165890850E-01_wp,  0.00000000000000E+00_wp, &
      &  1.87167066396519E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.05870802283058E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.52102452913589E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.71872094404136E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.11593011135074E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.95153090409821E-19_wp, &
      &  0.00000000000000E+00_wp,  3.65193452759191E-20_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.27572276160301E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  6.48656428418055E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  8.51298096285352E-18_wp,  0.00000000000000E+00_wp, -4.26924026302343E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.06392206676407E-02_wp,  0.00000000000000E+00_wp, &
      & -1.71872094404136E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.75266671550807E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.22114399857983E-19_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.87167066396519E-01_wp, &
      &  0.00000000000000E+00_wp, -3.70623786546151E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  4.30581076031436E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.76867634334542E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.71872094404136E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.11593011135074E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  3.65193452759191E-20_wp,  0.00000000000000E+00_wp, &
      & -1.95153090409821E-19_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.27572276160301E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  6.48656428418055E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -4.26924026302343E-03_wp, &
      &  0.00000000000000E+00_wp, -8.51298096285352E-18_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.55267110074857E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.94545602278596E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.71051711379041E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  4.63785620983120E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -4.11593011135074E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.55267110074857E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.46082884787843E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -6.48656428418055E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  4.10045094066887E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.37365712309899E-19_wp, &
      &  0.00000000000000E+00_wp, -6.81029027818177E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -4.75266671550807E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.55267110074857E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.61010314851148E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.05870802283058E-02_wp,  0.00000000000000E+00_wp, &
      & -4.30581076031437E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.56892474799869E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  3.66210140844275E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -4.11593011135074E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.55267110074857E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.46082884787843E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -6.48656428418055E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  4.10045094066887E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -6.81029027818177E-02_wp,  0.00000000000000E+00_wp, &
      & -1.37365712309899E-19_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.55267110074857E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.94545602278596E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.71051711379041E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -3.43414280774747E-20_wp, &
      &  0.00000000000000E+00_wp,  4.63785620983120E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.95153090409821E-19_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.86583909524310E-01_wp,  0.00000000000000E+00_wp,  8.49892475232076E-19_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.74162681065267E-17_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -7.84443966800847E-03_wp,  0.00000000000000E+00_wp, &
      & -6.61427873242330E-19_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.94545602278596E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.86583909524310E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.63785620983120E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.90880838382641E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  3.65193452759191E-20_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.46082884787843E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  8.49892475232076E-19_wp,  0.00000000000000E+00_wp, &
      & -3.86583909524310E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -4.26924026302340E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  6.81029027818177E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  5.79076948848866E-18_wp, &
      &  0.00000000000000E+00_wp, -3.06105121967253E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.51010910312432E-19_wp,  0.00000000000000E+00_wp,  1.22114399857983E-19_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.61010314851148E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.86583909524310E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.52102452913590E-02_wp,  0.00000000000000E+00_wp, &
      & -1.76867634334542E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -3.66210140844275E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.98236873369862E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  3.65193452759191E-20_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.46082884787843E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.86583909524310E-01_wp,  0.00000000000000E+00_wp, -6.59292117557936E-19_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -4.26924026302340E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  6.81029027818177E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -3.06105121967253E-02_wp,  0.00000000000000E+00_wp, &
      & -5.79076948848866E-18_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.94545602278596E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.86583909524310E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.63785620983120E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.41704402753271E-20_wp, &
      &  0.00000000000000E+00_wp,  2.90880838382641E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.95153090409821E-19_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -6.59292117557936E-19_wp,  0.00000000000000E+00_wp, &
      & -3.86583909524310E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.74162681065267E-17_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  6.61427873242330E-19_wp, &
      &  0.00000000000000E+00_wp, -7.84443966800847E-03_wp, -2.23075165890850E-01_wp, &
      &  0.00000000000000E+00_wp, -1.87167066396519E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.05870802283058E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.52102452913590E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.87625642334199E-01_wp,  0.00000000000000E+00_wp,  1.06392206676407E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  8.50659238294207E-19_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  9.58827975233050E-35_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.27572276160301E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -6.48656428418055E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.74162681065267E-17_wp,  0.00000000000000E+00_wp, &
      & -4.26924026302340E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.71872094404136E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  4.11593011135074E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  3.65193452759193E-20_wp, &
      &  0.00000000000000E+00_wp,  1.95153090409822E-19_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.87167066396519E-01_wp,  0.00000000000000E+00_wp, -3.70623786546151E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -4.30581076031437E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.76867634334542E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.06392206676407E-02_wp,  0.00000000000000E+00_wp, &
      & -1.71872094404136E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  4.75266671550807E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -3.53786835543724E-19_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.27572276160301E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -6.48656428418055E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -4.26924026302340E-03_wp,  0.00000000000000E+00_wp, -1.74162681065267E-17_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.71872094404136E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  4.11593011135074E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.95153090409822E-19_wp,  0.00000000000000E+00_wp, &
      & -3.65193452759193E-20_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -4.71051711379041E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -4.63785620983120E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.55267110074857E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.94545602278596E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  6.48656428418055E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  4.10045094066887E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  6.81029027818177E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  4.11593011135074E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.55267110074857E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.23320443974980E-19_wp,  0.00000000000000E+00_wp, &
      &  2.46082884787843E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.05870802283058E-02_wp, &
      &  0.00000000000000E+00_wp,  4.30581076031436E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.56892474799869E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -3.66210140844275E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  8.50659238294207E-19_wp,  0.00000000000000E+00_wp,  4.75266671550807E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.55267110074857E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.61010314851148E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  6.48656428418055E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  4.10045094066887E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  6.81029027818177E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  4.11593011135074E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.55267110074857E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.46082884787843E-03_wp,  0.00000000000000E+00_wp,  2.23320443974980E-19_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -4.71051711379041E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -4.63785620983120E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.55267110074857E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.94545602278596E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  8.51298096285352E-18_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.37365712309899E-19_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -7.84443966800847E-03_wp, &
      &  0.00000000000000E+00_wp,  5.79076948848866E-18_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  3.65193452759193E-20_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.23320443974980E-19_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.86583909524310E-01_wp,  0.00000000000000E+00_wp, -2.09889486718617E-18_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  4.63785620983120E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.90880838382641E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.94545602278596E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.86583909524310E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.26924026302343E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -6.81029027818177E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -6.61427873242330E-19_wp,  0.00000000000000E+00_wp, -3.06105121967253E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.95153090409822E-19_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.46082884787843E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.09889486718617E-18_wp,  0.00000000000000E+00_wp, &
      & -3.86583909524310E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.52102452913589E-02_wp, &
      &  0.00000000000000E+00_wp, -1.76867634334542E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  3.66210140844275E-02_wp, &
      &  0.00000000000000E+00_wp, -3.43414280774747E-20_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.98236873369862E-02_wp, &
      &  0.00000000000000E+00_wp, -1.41704402753271E-20_wp,  0.00000000000000E+00_wp, &
      &  9.58827975233050E-35_wp,  0.00000000000000E+00_wp, -3.53786835543724E-19_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.61010314851148E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.86583909524310E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.26924026302343E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -6.81029027818177E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -3.06105121967253E-02_wp, &
      &  0.00000000000000E+00_wp,  6.61427873242330E-19_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.95153090409822E-19_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.46082884787843E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.86583909524310E-01_wp,  0.00000000000000E+00_wp,  2.90661545467494E-18_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  4.63785620983120E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.90880838382641E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.94545602278596E-03_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.86583909524310E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -8.51298096285352E-18_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.37365712309899E-19_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.79076948848866E-18_wp,  0.00000000000000E+00_wp, -7.84443966800847E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.65193452759193E-20_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.23320443974980E-19_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.90661545467494E-18_wp,  0.00000000000000E+00_wp, &
      & -3.86583909524310E-01_wp], shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "Ce2")
   call test_hamiltonian_gen(error, mol, make_qvszp_basis, gxtb_h0spec(mol), &
      & make_gxtb_ncoord, hamiltonian)

end subroutine test_hamiltonian_gxtb_ce2


subroutine test_hamiltonian_ceh_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 2
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -5.20573260467584E-01_wp,  -5.12556505508839E-01_wp,  -5.12556505508839E-01_wp, &
      & -5.20573260467584E-01_wp], shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_hamiltonian_gen(error, mol, make_ceh_basis, ceh_h0spec(mol), &
      & make_ceh_ncoord, hamiltonian, thr_in=thr1*10)

end subroutine test_hamiltonian_ceh_h2

subroutine test_hamiltonian_ceh_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 5
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -5.76141826967410E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.64597683251088E-01_wp,  0.00000000000000E+00_wp, &
      & -1.30577038544615E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.30577038544615E-01_wp,  0.00000000000000E+00_wp, -2.17667285277957E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.30577038544615E-01_wp,  0.00000000000000E+00_wp, -1.64597683251088E-01_wp, &
      &  0.00000000000000E+00_wp, -2.17667285277957E-01_wp,  0.00000000000000E+00_wp, &
      & -3.69857613492308E-01_wp], shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_hamiltonian_gen(error, mol, make_ceh_basis, ceh_h0spec(mol), &
      & make_ceh_ncoord, hamiltonian, thr_in=thr1*10)

end subroutine test_hamiltonian_ceh_lih

subroutine test_hamiltonian_ceh_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 18
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -6.90083044966714E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.08003200438262E-01_wp,  0.00000000000000E+00_wp,  2.54418424031753E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.19476852196354E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.62742084015788E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.36072924446956E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  9.36704132556010E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.62742084015788E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.54418424031753E-01_wp,  0.00000000000000E+00_wp,  4.02200104061939E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.56606945484210E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.62742084015788E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.36072924446956E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  9.36704132556010E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.73436945976887E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.31497588937989E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.73436945976887E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -9.36704132556010E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.34814881206664E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.73436945976887E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.19476852196354E-01_wp,  0.00000000000000E+00_wp,  1.56606945484210E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.50747597557354E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.73436945976887E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -9.36704132556010E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.34814881206664E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.73436945976887E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.31497588937989E-03_wp, &
      & -1.08003200438262E-01_wp,  0.00000000000000E+00_wp, -2.54418424031753E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.19476852196354E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -6.90083044966714E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.36072924446956E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -9.36704132556010E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.62742084015788E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.54418424031753E-01_wp,  0.00000000000000E+00_wp,  4.02200104061939E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.56606945484210E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.62742084015788E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.36072924446956E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -9.36704132556010E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.62742084015788E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.31497588937989E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.73436945976887E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  9.36704132556010E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.34814881206664E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.73436945976887E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.19476852196354E-01_wp,  0.00000000000000E+00_wp, -1.56606945484210E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.50747597557354E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -5.73436945976887E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  9.36704132556010E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.34814881206664E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -5.73436945976887E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.31497588937989E-03_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.73436945976887E-02_wp],&
      & shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "S2")
   call test_hamiltonian_gen(error, mol, make_ceh_basis, ceh_h0spec(mol), &
      & make_ceh_ncoord, hamiltonian, thr_in=thr1*10)

end subroutine test_hamiltonian_ceh_s2

subroutine test_hamiltonian_ceh_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nao = 13
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -7.08495044644038E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.38605666294581E-01_wp, -3.38605666294581E-01_wp, -3.38605666294581E-01_wp, &
      & -3.38605666294581E-01_wp,  0.00000000000000E+00_wp, -4.76056389727416E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -2.31521047829169E-01_wp,  2.31521047829169E-01_wp, &
      &  2.31521047829169E-01_wp, -2.31521047829169E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.76056389727416E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  2.31521047829169E-01_wp, &
      &  2.31521047829169E-01_wp, -2.31521047829169E-01_wp, -2.31521047829169E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -4.76056389727416E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.31521047829169E-01_wp,  2.31521047829169E-01_wp, -2.31521047829169E-01_wp, &
      &  2.31521047829169E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.85417046536822E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.59872557670323E-01_wp, -1.59872557670323E-01_wp, &
      &  1.59872557670323E-01_wp,  1.59872557670323E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.85417046536822E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.59872557670323E-01_wp, &
      & -1.59872557670323E-01_wp,  1.59872557670323E-01_wp, -1.59872557670323E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.85417046536822E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.72825675197766E-17_wp, -2.72825675197766E-17_wp, -2.72825675197766E-17_wp, &
      & -2.72825675197766E-17_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -1.85417046536822E-01_wp, &
      &  0.00000000000000E+00_wp,  1.59872557670323E-01_wp, -1.59872557670323E-01_wp, &
      & -1.59872557670323E-01_wp,  1.59872557670323E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -1.85417046536822E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.38605666294581E-01_wp, -2.31521047829169E-01_wp,  2.31521047829169E-01_wp, &
      & -2.31521047829169E-01_wp, -1.59872557670323E-01_wp,  1.59872557670323E-01_wp, &
      & -2.72825675197766E-17_wp,  1.59872557670323E-01_wp,  0.00000000000000E+00_wp, &
      & -4.86526446971988E-01_wp, -4.43220553944223E-02_wp, -4.43220553944223E-02_wp, &
      & -4.43220553944223E-02_wp, -3.38605666294581E-01_wp,  2.31521047829169E-01_wp, &
      &  2.31521047829169E-01_wp,  2.31521047829169E-01_wp, -1.59872557670323E-01_wp, &
      & -1.59872557670323E-01_wp, -2.72825675197766E-17_wp, -1.59872557670323E-01_wp, &
      &  0.00000000000000E+00_wp, -4.43220553944223E-02_wp, -4.86526446971988E-01_wp, &
      & -4.43220553944223E-02_wp, -4.43220553944223E-02_wp, -3.38605666294581E-01_wp, &
      &  2.31521047829169E-01_wp, -2.31521047829169E-01_wp, -2.31521047829169E-01_wp, &
      &  1.59872557670323E-01_wp,  1.59872557670323E-01_wp, -2.72825675197766E-17_wp, &
      & -1.59872557670323E-01_wp,  0.00000000000000E+00_wp, -4.43220553944223E-02_wp, &
      & -4.43220553944223E-02_wp, -4.86526446971988E-01_wp, -4.43220553944223E-02_wp, &
      & -3.38605666294581E-01_wp, -2.31521047829169E-01_wp, -2.31521047829169E-01_wp, &
      &  2.31521047829169E-01_wp,  1.59872557670323E-01_wp, -1.59872557670323E-01_wp, &
      & -2.72825675197766E-17_wp,  1.59872557670323E-01_wp,  0.00000000000000E+00_wp, &
      & -4.43220553944223E-02_wp, -4.43220553944223E-02_wp, -4.43220553944223E-02_wp, &
      & -4.86526446971988E-01_wp], shape(hamiltonian))

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_hamiltonian_gen(error, mol, make_ceh_basis, ceh_h0spec(mol), &
      & make_ceh_ncoord, hamiltonian, thr_in=thr1*10)

end subroutine test_hamiltonian_ceh_sih4

subroutine test_hamiltonian_ceh_panp(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   integer, parameter :: nao = 32
   real(wp), parameter :: hamiltonian(nao, nao) = reshape([&
      & -5.74853679880744E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -1.66451247411204E-1_wp,  0.00000000000000E+0_wp, &
      &  4.47287558776629E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -2.81497449467025E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  1.43046322995670E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -1.98501599971948E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -6.72422997001521E-2_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  7.79400568668741E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.07353930537686E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -1.98501599971948E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -2.92653724182373E-1_wp, &
      &  0.00000000000000E+0_wp,  1.13711621512660E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -2.26026875478745E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  5.71010461356486E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -1.98501599971948E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -6.72422997001521E-2_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  7.79400568668741E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.07353930537686E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.22477733305189E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -7.31025577232244E-2_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  8.28237773528076E-2_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.22477733305189E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -6.09662854995410E-2_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  1.39411396274477E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -1.81621678388459E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.22477733305189E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -2.68834064553585E-1_wp,  0.00000000000000E+0_wp, &
      & -3.66147236116575E-2_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -2.59793482034505E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  1.01762616199382E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.22477733305189E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -6.09662854995410E-2_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  1.39411396274477E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -1.81621678388459E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.22477733305189E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -7.31025577232244E-2_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  8.28237773528076E-2_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.62889659406096E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -4.46979917198434E-2_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.62889659406096E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -8.09796791881055E-2_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  2.04475579369224E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.62889659406096E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -2.96021748823223E-2_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  1.79738538099156E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -2.46023684066937E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.62889659406096E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -9.82268170755756E-2_wp,  0.00000000000000E+0_wp, &
      & -3.66405043298728E-2_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -1.02885773548641E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  1.55511095251092E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.62889659406096E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -2.96021748823223E-2_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  1.79738538099156E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -2.46023684066937E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.62889659406096E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -8.09796791881055E-2_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  2.04475579369224E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.62889659406096E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -4.46979917198434E-2_wp, -1.66451247411204E-1_wp, &
      &  0.00000000000000E+0_wp, -2.92653724182373E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -2.68834064553585E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -9.82268170755756E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.45780799289089E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -6.72422997001521E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -6.09662854995410E-2_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -2.96021748823223E-2_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -1.36147828319041E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  4.47287558776629E-1_wp,  0.00000000000000E+0_wp,  1.13711621512660E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.66147236116575E-2_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.66405043298728E-2_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -1.36147828319041E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -6.72422997001521E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -6.09662854995410E-2_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -2.96021748823223E-2_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -1.36147828319041E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -7.31025577232244E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -8.09796791881055E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.49035738267303E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  7.79400568668741E-2_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  1.39411396274477E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  1.79738538099156E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.49035738267303E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -2.81497449467025E-1_wp, &
      &  0.00000000000000E+0_wp, -2.26026875478745E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -2.59793482034505E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.02885773548641E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.49035738267303E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  7.79400568668741E-2_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  1.39411396274477E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  1.79738538099156E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.49035738267303E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -7.31025577232244E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -8.09796791881055E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -3.49035738267303E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -4.46979917198434E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.37809711810083E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  8.28237773528076E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  2.04475579369224E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.37809711810083E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -1.07353930537686E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.81621678388459E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -2.46023684066937E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.37809711810083E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  1.43046322995670E-1_wp, &
      &  0.00000000000000E+0_wp,  5.71010461356486E-2_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  1.01762616199382E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  1.55511095251092E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.37809711810083E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp, -1.07353930537686E-1_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -1.81621678388459E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -2.46023684066937E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.37809711810083E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  8.28237773528076E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  2.04475579369224E-1_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.37809711810083E-1_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, -4.46979917198434E-2_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      &  0.00000000000000E+0_wp,  0.00000000000000E+0_wp,  0.00000000000000E+0_wp, &
      & -4.37809711810083E-1_wp], shape(hamiltonian))

   call get_structure(mol, "f-block", "PaNp")
   call test_hamiltonian_gen(error, mol, make_ceh_basis, ceh_h0spec(mol), &
      & make_ceh_ncoord, hamiltonian, thr_in=thr1*10)

end subroutine test_hamiltonian_ceh_panp

subroutine test_anisotropy_gradient_gxtb_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "H2")
   call test_anisotropy_grad(error, mol, make_qvszp_basis, gxtb_h0spec(mol), thr_in=thr1)

end subroutine test_anisotropy_gradient_gxtb_h2

subroutine test_anisotropy_gradient_gxtb_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_anisotropy_grad(error, mol, make_qvszp_basis, gxtb_h0spec(mol), thr_in=thr1)

end subroutine test_anisotropy_gradient_gxtb_lih

subroutine test_anisotropy_gradient_gxtb_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "SiH4")
   call test_anisotropy_grad(error, mol, make_qvszp_basis, gxtb_h0spec(mol), thr_in=thr1)

end subroutine test_anisotropy_gradient_gxtb_sih4

subroutine test_anisotropy_gradient_gxtb_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "f-block", "CeCl3")
   call test_anisotropy_grad(error, mol, make_qvszp_basis, gxtb_h0spec(mol), thr_in=thr1)

end subroutine test_anisotropy_gradient_gxtb_cecl3

subroutine test_anisotropy_gradient_gxtb_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_anisotropy_grad(error, mol, make_qvszp_basis, gxtb_h0spec(mol), thr_in=thr1)

end subroutine test_anisotropy_gradient_gxtb_m01

subroutine test_anisotropy_gradient_gxtb_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_anisotropy_grad(error, mol, make_qvszp_basis, gxtb_h0spec(mol), thr_in=thr1)

end subroutine test_anisotropy_gradient_gxtb_m02


subroutine test_g_energy_aniso_gxtb_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :)
   
   allocate(density(2, 2, 1))
   density(:, :, :) = reshape([&
      & 5.93683766916992E-01_wp, 5.93683766916992E-01_wp, 5.93683766916992E-01_wp, &
      & 5.93683766916992E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "H2")
   call test_energy_anisotropy_gradient(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), thr_in=thr1)

end subroutine test_g_energy_aniso_gxtb_h2

subroutine test_g_energy_aniso_gxtb_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :)
   
   allocate(density(5, 5, 1))
   density(:, :, :) = reshape([&
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
   call test_energy_anisotropy_gradient(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), thr_in=thr1)

end subroutine test_g_energy_aniso_gxtb_lih

subroutine test_g_energy_aniso_gxtb_no(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :)
   
   allocate(density(8, 8, 2))
   density(:, :, :) = reshape([&
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
   call test_energy_anisotropy_gradient(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), thr_in=thr1*10)

end subroutine test_g_energy_aniso_gxtb_no

subroutine test_g_energy_aniso_gxtb_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :)
   
   allocate(density(18, 18, 1))
   density(:, :, :) = reshape([&
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
   call test_energy_anisotropy_gradient(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), thr_in=thr1*10)

end subroutine test_g_energy_aniso_gxtb_s2

subroutine test_g_energy_aniso_gxtb_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :)
   
   allocate(density(13, 13, 1))
   density(:, :, :) = reshape([&
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
   call test_energy_anisotropy_gradient(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), thr_in=thr1*10)

end subroutine test_g_energy_aniso_gxtb_sih4

subroutine test_g_energy_aniso_gxtb_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :)
   
   allocate(density(43, 43, 2))
   density(:, :, :) = reshape([&
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
   call test_energy_anisotropy_gradient(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), thr_in=thr1*100)

end subroutine test_g_energy_aniso_gxtb_cecl3


subroutine test_qg_hamiltonian_gxtb_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(2, 2, 1))
   density(:, :, :) = reshape([&
      & 5.93683766916992E-01_wp, 5.93683766916992E-01_wp, 5.93683766916992E-01_wp, &
      & 5.93683766916992E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "H2")
   call test_hamiltonian_qeff_numgrad(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), make_gxtb_ncoord, thr_in=thr1*10)

end subroutine test_qg_hamiltonian_gxtb_h2

subroutine test_qg_hamiltonian_gxtb_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(5, 5, 1))
   density(:, :, :) = reshape([&
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
   call test_hamiltonian_qeff_numgrad(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), make_gxtb_ncoord, thr_in=thr1)

end subroutine test_qg_hamiltonian_gxtb_lih

subroutine test_qg_hamiltonian_gxtb_no(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(8, 8, 2))
   density(:, :, :) = reshape([&
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
   call test_hamiltonian_qeff_numgrad(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), make_gxtb_ncoord, thr_in=thr1*10)

end subroutine test_qg_hamiltonian_gxtb_no

subroutine test_qg_hamiltonian_gxtb_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(18, 18, 1))
   density(:, :, :) = reshape([&
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
   call test_hamiltonian_qeff_numgrad(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), make_gxtb_ncoord, thr_in=thr1*10)

end subroutine test_qg_hamiltonian_gxtb_s2

subroutine test_qg_hamiltonian_gxtb_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(13, 13, 1))
   density(:, :, :) = reshape([&
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
   call test_hamiltonian_qeff_numgrad(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), make_gxtb_ncoord, thr_in=thr1*10)

end subroutine test_qg_hamiltonian_gxtb_sih4

subroutine test_qg_hamiltonian_gxtb_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(43, 43, 2))
   density(:, :, :) = reshape([&
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
   call test_hamiltonian_qeff_numgrad(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), make_gxtb_ncoord, thr_in=thr1*10)

end subroutine test_qg_hamiltonian_gxtb_cecl3


subroutine test_g_hamiltonian_gfn2_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(2, 2, 1))
   density(:, :, :) = reshape([&
      & 6.01869285390526E-1_wp, 6.01869285390526E-1_wp, 6.01869285390526E-1_wp, &
      & 6.01869285390526E-1_wp], shape(density))

   call get_structure(mol, "MB16-43", "H2")
   call test_hamiltonian_numgrad(error, mol, density, make_gfn2_basis, &
      & gfn2_h0spec(mol), make_gfn2_ncoord, thr_in=thr1)

end subroutine test_g_hamiltonian_gfn2_h2

subroutine test_g_hamiltonian_gfn2_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(5, 5, 1))
   density(:, :, :) = reshape([&
      & 2.19539073988330E-1_wp, 0.00000000000000E+0_wp, 1.40923711006316E-1_wp, &
      & 0.00000000000000E+0_wp, 4.75256757831863E-1_wp, 0.00000000000000E+0_wp, &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, &
      & 0.00000000000000E+0_wp, 1.40923711006316E-1_wp, 0.00000000000000E+0_wp, &
      & 9.04599439316547E-2_wp, 0.00000000000000E+0_wp, 3.05070731955700E-1_wp, &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, &
      & 0.00000000000000E+0_wp, 0.00000000000000E+0_wp, 4.75256757831863E-1_wp, &
      & 0.00000000000000E+0_wp, 3.05070731955700E-1_wp, 0.00000000000000E+0_wp, &
      & 1.02883273469971E+0_wp], shape(density))

   call get_structure(mol, "MB16-43", "LiH")
   call test_hamiltonian_numgrad(error, mol, density, make_gfn2_basis, &
      & gfn2_h0spec(mol), make_gfn2_ncoord, thr_in=thr1)

end subroutine test_g_hamiltonian_gfn2_lih

subroutine test_g_hamiltonian_gfn2_no(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(8, 8, 2))
   ! Density matrix using the spin-polarized hamiltonian
   density(:, :, :) = reshape([&
      &  7.74679282212467E-01_wp,  0.00000000000000E+00_wp, -3.61378196533156E-01_wp, &
      &  0.00000000000000E+00_wp,  4.00028571475358E-02_wp,  0.00000000000000E+00_wp, &
      & -7.88582402887757E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  6.67585512041934E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  6.36483804983739E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -3.61378196533156E-01_wp,  0.00000000000000E+00_wp, &
      &  4.83753129547642E-01_wp,  0.00000000000000E+00_wp,  6.51793701014339E-02_wp, &
      &  0.00000000000000E+00_wp, -3.03648635028530E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  6.67585512047539E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  6.36483804938595E-02_wp,  4.00028571475358E-02_wp, &
      &  0.00000000000000E+00_wp,  6.51793701014339E-02_wp,  0.00000000000000E+00_wp, &
      &  8.35070730887832E-01_wp,  0.00000000000000E+00_wp,  2.97112877190300E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  6.36483804983739E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  8.01613055806132E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -7.88582402887757E-02_wp,  0.00000000000000E+00_wp, -3.03648635028530E-01_wp, &
      &  0.00000000000000E+00_wp,  2.97112877190300E-01_wp,  0.00000000000000E+00_wp, &
      &  5.65044622268906E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  6.36483804938595E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  8.01613055809767E-01_wp,  7.74679282212467E-01_wp,  0.00000000000000E+00_wp, &
      & -3.61378196533156E-01_wp,  0.00000000000000E+00_wp,  4.00028571475358E-02_wp, &
      &  0.00000000000000E+00_wp, -7.88582402887757E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.72982732227426E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  3.84309747826930E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -3.61378196533156E-01_wp, &
      &  0.00000000000000E+00_wp,  4.83753129547642E-01_wp,  0.00000000000000E+00_wp, &
      &  6.51793701014339E-02_wp,  0.00000000000000E+00_wp, -3.03648635028530E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  2.72982732227498E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  3.84309747826951E-01_wp, &
      &  4.00028571475358E-02_wp,  0.00000000000000E+00_wp,  6.51793701014339E-02_wp, &
      &  0.00000000000000E+00_wp,  8.35070730887832E-01_wp,  0.00000000000000E+00_wp, &
      &  2.97112877190300E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  3.84309747826930E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  5.41037819753934E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -7.88582402887757E-02_wp,  0.00000000000000E+00_wp, &
      & -3.03648635028530E-01_wp,  0.00000000000000E+00_wp,  2.97112877190300E-01_wp, &
      &  0.00000000000000E+00_wp,  5.65044622268906E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  3.84309747826951E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  5.41037819753851E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "NO")
   call test_hamiltonian_numgrad(error, mol, density, make_gfn2_basis, &
      & gfn2_h0spec(mol), make_gfn2_ncoord, thr_in=thr1*10)

end subroutine test_g_hamiltonian_gfn2_no

subroutine test_g_hamiltonian_gfn2_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(18, 18, 1))
   density(:, :, :) = reshape([&
      &  2.08383407952596E+00_wp,  0.00000000000000E+00_wp, -2.86885457434526E-01_wp, &
      &  0.00000000000000E+00_wp, -1.16562682086441E-12_wp,  0.00000000000000E+00_wp, &
      & -6.69739674887478E-04_wp,  0.00000000000000E+00_wp,  2.78018424605231E-12_wp, &
      & -3.99763276513242E-01_wp,  0.00000000000000E+00_wp,  7.20789576967972E-02_wp, &
      &  0.00000000000000E+00_wp,  1.56518653223086E-13_wp,  0.00000000000000E+00_wp, &
      &  1.98462616863759E-03_wp,  0.00000000000000E+00_wp, -3.72373014274475E-13_wp, &
      &  0.00000000000000E+00_wp,  1.40640424983986E+00_wp,  0.00000000000000E+00_wp, &
      & -7.22875156360156E-11_wp,  0.00000000000000E+00_wp, -4.28929526173378E-02_wp, &
      &  0.00000000000000E+00_wp,  9.41288664490691E-12_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  3.60016476420680E-01_wp,  0.00000000000000E+00_wp, &
      &  7.23104311110505E-11_wp,  0.00000000000000E+00_wp, -9.32017279406466E-02_wp, &
      &  0.00000000000000E+00_wp,  9.45829764333529E-12_wp,  0.00000000000000E+00_wp, &
      & -2.86885457434526E-01_wp,  0.00000000000000E+00_wp,  8.09558374407343E-01_wp, &
      &  0.00000000000000E+00_wp,  1.13679924111833E-13_wp,  0.00000000000000E+00_wp, &
      & -2.27717812278863E-02_wp,  0.00000000000000E+00_wp, -2.71178933012841E-13_wp, &
      & -7.20789576967987E-02_wp,  0.00000000000000E+00_wp, -7.90979746024854E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -2.30013575054105E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -7.22875292213000E-11_wp,  0.00000000000000E+00_wp, &
      &  1.40640425018458E+00_wp,  0.00000000000000E+00_wp,  9.41308207589007E-12_wp, &
      &  0.00000000000000E+00_wp, -4.28929526622254E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  7.23096931461026E-11_wp,  0.00000000000000E+00_wp, &
      &  3.60016476075854E-01_wp,  0.00000000000000E+00_wp,  9.45862783964508E-12_wp, &
      &  0.00000000000000E+00_wp, -9.32017279857499E-02_wp,  0.00000000000000E+00_wp, &
      & -1.16562682086441E-12_wp,  0.00000000000000E+00_wp,  1.13679924111833E-13_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.56162495193329E-13_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -4.28929526173378E-02_wp,  0.00000000000000E+00_wp, &
      &  9.41308296146505E-12_wp,  0.00000000000000E+00_wp,  9.56674657154071E-03_wp, &
      &  0.00000000000000E+00_wp, -1.23270777661548E-12_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  9.32017279406464E-02_wp,  0.00000000000000E+00_wp, &
      & -9.45824745645229E-12_wp,  0.00000000000000E+00_wp,  8.13392095026911E-03_wp, &
      &  0.00000000000000E+00_wp, -1.23012011470127E-12_wp,  0.00000000000000E+00_wp, &
      & -6.69739674887478E-04_wp,  0.00000000000000E+00_wp, -2.27717812278863E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  6.80924065724479E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  1.98462616863702E-03_wp,  0.00000000000000E+00_wp,  2.30013575054114E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  6.78087189604356E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  9.41288946212473E-12_wp,  0.00000000000000E+00_wp, &
      & -4.28929526622254E-02_wp,  0.00000000000000E+00_wp, -1.23270780717507E-12_wp, &
      &  0.00000000000000E+00_wp,  9.56674657741918E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -9.45824812834174E-12_wp,  0.00000000000000E+00_wp, &
      &  9.32017279857503E-02_wp,  0.00000000000000E+00_wp, -1.23012395220030E-12_wp, &
      &  0.00000000000000E+00_wp,  8.13392095613514E-03_wp,  0.00000000000000E+00_wp, &
      &  2.78018424605231E-12_wp,  0.00000000000000E+00_wp, -2.71178933012841E-13_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.72536267101461E-13_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -3.99763276513242E-01_wp,  0.00000000000000E+00_wp, -7.20789576967987E-02_wp, &
      &  0.00000000000000E+00_wp,  1.56162495193329E-13_wp,  0.00000000000000E+00_wp, &
      &  1.98462616863702E-03_wp,  0.00000000000000E+00_wp, -3.72536267101461E-13_wp, &
      &  2.08383407952596E+00_wp,  0.00000000000000E+00_wp,  2.86885457434528E-01_wp, &
      &  0.00000000000000E+00_wp, -1.16561088864467E-12_wp,  0.00000000000000E+00_wp, &
      & -6.69739674887130E-04_wp,  0.00000000000000E+00_wp,  2.77970757907873E-12_wp, &
      &  0.00000000000000E+00_wp,  3.60016476420680E-01_wp,  0.00000000000000E+00_wp, &
      &  7.23096934386063E-11_wp,  0.00000000000000E+00_wp,  9.32017279406464E-02_wp, &
      &  0.00000000000000E+00_wp, -9.45824531112395E-12_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  1.40640424983985E+00_wp,  0.00000000000000E+00_wp, &
      & -7.22883526483221E-11_wp,  0.00000000000000E+00_wp,  4.28929526173370E-02_wp, &
      &  0.00000000000000E+00_wp, -9.41292333286485E-12_wp,  0.00000000000000E+00_wp, &
      &  7.20789576967972E-02_wp,  0.00000000000000E+00_wp, -7.90979746024854E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.30013575054114E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.86885457434528E-01_wp,  0.00000000000000E+00_wp,  8.09558374407336E-01_wp, &
      &  0.00000000000000E+00_wp, -1.13733004181200E-13_wp,  0.00000000000000E+00_wp, &
      &  2.27717812278852E-02_wp,  0.00000000000000E+00_wp,  2.71107064233052E-13_wp, &
      &  0.00000000000000E+00_wp,  7.23104446963354E-11_wp,  0.00000000000000E+00_wp, &
      &  3.60016476075854E-01_wp,  0.00000000000000E+00_wp, -9.45824483615390E-12_wp, &
      &  0.00000000000000E+00_wp,  9.32017279857503E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -7.22883523558188E-11_wp,  0.00000000000000E+00_wp, &
      &  1.40640425018457E+00_wp,  0.00000000000000E+00_wp, -9.41298996900585E-12_wp, &
      &  0.00000000000000E+00_wp,  4.28929526622246E-02_wp,  0.00000000000000E+00_wp, &
      &  1.56518653223086E-13_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -1.16561088864467E-12_wp,  0.00000000000000E+00_wp, -1.13733004181200E-13_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -9.32017279406466E-02_wp,  0.00000000000000E+00_wp, &
      &  9.45862872522002E-12_wp,  0.00000000000000E+00_wp,  8.13392095026911E-03_wp, &
      &  0.00000000000000E+00_wp, -1.23012392164070E-12_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  4.28929526173370E-02_wp,  0.00000000000000E+00_wp, &
      & -9.41299258930421E-12_wp,  0.00000000000000E+00_wp,  9.56674657154062E-03_wp, &
      &  0.00000000000000E+00_wp, -1.23274238171834E-12_wp,  0.00000000000000E+00_wp, &
      &  1.98462616863759E-03_wp,  0.00000000000000E+00_wp, -2.30013575054105E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  6.78087189604356E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -6.69739674887130E-04_wp,  0.00000000000000E+00_wp,  2.27717812278852E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  6.80924065724415E-04_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  9.45829872582955E-12_wp,  0.00000000000000E+00_wp, &
      & -9.32017279857499E-02_wp,  0.00000000000000E+00_wp, -1.23012014526085E-12_wp, &
      &  0.00000000000000E+00_wp,  8.13392095613514E-03_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp, -9.41292441535906E-12_wp,  0.00000000000000E+00_wp, &
      &  4.28929526622246E-02_wp,  0.00000000000000E+00_wp, -1.23274241227792E-12_wp, &
      &  0.00000000000000E+00_wp,  9.56674657741898E-03_wp,  0.00000000000000E+00_wp, &
      & -3.72373014274475E-13_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.77970757907873E-12_wp,  0.00000000000000E+00_wp,  2.71107064233052E-13_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp],&
      & shape(density))

   call get_structure(mol, "MB16-43", "S2")
   call test_hamiltonian_numgrad(error, mol, density, make_gfn2_basis, &
      & gfn2_h0spec(mol), make_gfn2_ncoord, thr_in=thr1*10)

end subroutine test_g_hamiltonian_gfn2_s2

subroutine test_g_hamiltonian_gfn2_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(13, 13, 1))
   density(:, :, :) = reshape([&
      &  1.03168544260937E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.34756550848050E-01_wp,  2.34756550848050E-01_wp,  2.34756550848050E-01_wp, &
      &  2.34756550848050E-01_wp,  0.00000000000000E+00_wp,  3.42715656393563E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -6.96002001260335E-02_wp, &
      &  0.00000000000000E+00_wp,  2.88009541315403E-01_wp, -2.88009541315403E-01_wp, &
      & -2.88009541315403E-01_wp,  2.88009541315403E-01_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  3.42715656393564E-01_wp,  0.00000000000000E+00_wp, &
      & -6.96002001260335E-02_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -2.88009541315403E-01_wp, &
      & -2.88009541315403E-01_wp,  2.88009541315404E-01_wp,  2.88009541315403E-01_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  3.42715656393563E-01_wp,  0.00000000000000E+00_wp, -6.96002001260336E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.88009541315403E-01_wp, -2.88009541315403E-01_wp,  2.88009541315403E-01_wp, &
      & -2.88009541315403E-01_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      & -6.96002001260335E-02_wp,  0.00000000000000E+00_wp,  1.41347142075733E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  5.84902420995309E-02_wp,  5.84902420995308E-02_wp, &
      & -5.84902420995311E-02_wp, -5.84902420995310E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -6.96002001260336E-02_wp, &
      &  0.00000000000000E+00_wp,  1.41347142075733E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -5.84902420995308E-02_wp, &
      &  5.84902420995309E-02_wp, -5.84902420995311E-02_wp,  5.84902420995312E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, -6.96002001260335E-02_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  1.41347142075733E-02_wp, &
      &  0.00000000000000E+00_wp, -5.84902420995308E-02_wp,  5.84902420995309E-02_wp, &
      &  5.84902420995308E-02_wp, -5.84902420995312E-02_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  0.00000000000000E+00_wp,  0.00000000000000E+00_wp,  0.00000000000000E+00_wp, &
      &  2.34756550848050E-01_wp,  2.88009541315403E-01_wp, -2.88009541315403E-01_wp, &
      &  2.88009541315403E-01_wp,  5.84902420995309E-02_wp, -5.84902420995308E-02_wp, &
      &  0.00000000000000E+00_wp, -5.84902420995308E-02_wp,  0.00000000000000E+00_wp, &
      &  7.79525794260635E-01_wp, -1.88617847166152E-01_wp, -1.88617847166152E-01_wp, &
      & -1.88617847166151E-01_wp,  2.34756550848050E-01_wp, -2.88009541315403E-01_wp, &
      & -2.88009541315403E-01_wp, -2.88009541315403E-01_wp,  5.84902420995308E-02_wp, &
      &  5.84902420995309E-02_wp,  0.00000000000000E+00_wp,  5.84902420995309E-02_wp, &
      &  0.00000000000000E+00_wp, -1.88617847166152E-01_wp,  7.79525794260635E-01_wp, &
      & -1.88617847166152E-01_wp, -1.88617847166152E-01_wp,  2.34756550848050E-01_wp, &
      & -2.88009541315403E-01_wp,  2.88009541315404E-01_wp,  2.88009541315403E-01_wp, &
      & -5.84902420995311E-02_wp, -5.84902420995311E-02_wp,  0.00000000000000E+00_wp, &
      &  5.84902420995308E-02_wp,  0.00000000000000E+00_wp, -1.88617847166152E-01_wp, &
      & -1.88617847166152E-01_wp,  7.79525794260635E-01_wp, -1.88617847166151E-01_wp, &
      &  2.34756550848050E-01_wp,  2.88009541315403E-01_wp,  2.88009541315403E-01_wp, &
      & -2.88009541315403E-01_wp, -5.84902420995310E-02_wp,  5.84902420995312E-02_wp, &
      &  0.00000000000000E+00_wp, -5.84902420995312E-02_wp,  0.00000000000000E+00_wp, &
      & -1.88617847166151E-01_wp, -1.88617847166152E-01_wp, -1.88617847166151E-01_wp, &
      &  7.79525794260635E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "SiH4")
   call test_hamiltonian_numgrad(error, mol, density, make_gfn2_basis, &
      & gfn2_h0spec(mol), make_gfn2_ncoord, thr_in=thr1*10)

end subroutine test_g_hamiltonian_gfn2_sih4


subroutine test_g_hamiltonian_gxtb_h2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(2, 2, 1))
   density(:, :, :) = reshape([&
      & 5.93683766916992E-01_wp, 5.93683766916992E-01_wp, 5.93683766916992E-01_wp, &
      & 5.93683766916992E-01_wp], shape(density))

   call get_structure(mol, "MB16-43", "H2")
   call test_hamiltonian_numgrad(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), make_gxtb_ncoord, thr_in=thr1*10)

end subroutine test_g_hamiltonian_gxtb_h2

subroutine test_g_hamiltonian_gxtb_lih(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(5, 5, 1))
   density(:, :, :) = reshape([&
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
   call test_hamiltonian_numgrad(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), make_gxtb_ncoord, thr_in=thr1)

end subroutine test_g_hamiltonian_gxtb_lih

subroutine test_g_hamiltonian_gxtb_no(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(8, 8, 2))
   density(:, :, :) = reshape([&
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
   call test_hamiltonian_numgrad(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), make_gxtb_ncoord, thr_in=thr1*10)

end subroutine test_g_hamiltonian_gxtb_no

subroutine test_g_hamiltonian_gxtb_s2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(18, 18, 1))
   density(:, :, :) = reshape([&
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
   call test_hamiltonian_numgrad(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), make_gxtb_ncoord, thr_in=thr1*10)

end subroutine test_g_hamiltonian_gxtb_s2

subroutine test_g_hamiltonian_gxtb_sih4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(13, 13, 1))
   density(:, :, :) = reshape([&
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
   call test_hamiltonian_numgrad(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), make_gxtb_ncoord, thr_in=thr1*10)

end subroutine test_g_hamiltonian_gxtb_sih4

subroutine test_g_hamiltonian_gxtb_cecl3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: density(:, :, :) 
   
   allocate(density(43, 43, 2))
   density(:, :, :) = reshape([&
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
   call test_hamiltonian_numgrad(error, mol, density, make_qvszp_basis, &
      & gxtb_h0spec(mol), make_gxtb_ncoord, thr_in=thr1*10)

end subroutine test_g_hamiltonian_gxtb_cecl3

end module test_hamiltonian
