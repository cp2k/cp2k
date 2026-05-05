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

module test_coulomb_thirdorder
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use mstore, only : get_structure
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_basis_type
   use tblite_blas, only : gemv
   use tblite_ceh_ceh, only : get_effective_qat
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_thirdorder, only : onsite_thirdorder, new_onsite_thirdorder, &
      & twobody_thirdorder, new_twobody_thirdorder
   use tblite_cutoff, only : get_lattice_points
   use tblite_scf, only: new_potential, potential_type
   use tblite_utils_average, only : average_type, new_average, average_id
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_coulomb, only : tb_coulomb, new_coulomb
   implicit none
   private

   public :: collect_coulomb_thirdorder

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
subroutine collect_coulomb_thirdorder(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-atom-ons-gfn1-m01", test_e_onsite_gfn1_m01_atom), &
      new_unittest("energy-atom-ons-gfn2-m02", test_e_onsite_gfn2_m02_atom), &
      new_unittest("energy-shell-ons-gfn2-m07", test_e_onsite_gfn2_m07_shell), &
      new_unittest("energy-atom-ons-gfn2-oxacb-pbc", test_e_onsite_gfn2_oxacb_shell), &
      new_unittest("energy-atom-ons-gfn2-oxacb-sc", test_e_onsite_gfn2_oxacb_sc_atom), &
      new_unittest("energy-shell-twob-gxtb-lih", test_e_twobody_gxtb_lih_shell), &
      new_unittest("energy-atom-twob-gxtb-m01", test_e_twobody_gxtb_m01_atom), &
      new_unittest("energy-shell-twob-gxtb-m01", test_e_twobody_gxtb_m01_shell), &
      new_unittest("energy-atom-twob-gxtb-m02", test_e_twobody_gxtb_m02_atom), &
      new_unittest("energy-shell-twob-gxtb-m02", test_e_twobody_gxtb_m02_shell), &
      new_unittest("potential-atom-twob-gxtb-lih", test_p_twobody_gxtb_lih_shell), &
      new_unittest("potential-atom-twob-gxtb-m01", test_p_twobody_gxtb_m01_atom), &
      new_unittest("potential-shell-twob-gxtb-m01", test_p_twobody_gxtb_m01_shell), &
      new_unittest("potential-atom-twob-gxtb-m02", test_p_twobody_gxtb_m02_atom), &
      new_unittest("potential-shell-twob-gxtb-m02", test_p_twobody_gxtb_m02_shell), &
      new_unittest("taumat-deriv-shell-twob-gxtb-h2", test_taumat_twobody_gxtb_h2_shell), &
      new_unittest("taumat-deriv-shell-twob-gxtb-lih", test_taumat_twobody_gxtb_lih_shell), &
      new_unittest("taumat-deriv-shell-twob-gxtb-s2", test_taumat_twobody_gxtb_s2_shell), &
      new_unittest("taumat-deriv-shell-twob-gxtb-cecl3", test_taumat_twobody_gxtb_cecl3_shell), &
      new_unittest("gradient-atom-ons-gfn1-m03", test_g_onsite_gfn1_m03_atom), &
      new_unittest("gradient-atom-ons-gfn1-urea-pbc", test_g_onsite_gfn1_urea_atom), &
      new_unittest("gradient-atom-ons-gfn2-m04", test_g_onsite_gfn2_m04_atom), &
      new_unittest("gradient-shell-ons-gfn2-m08", test_g_onsite_gfn2_m08_shell), &
      new_unittest("gradient-shell-twob-gxtb-lih", test_g_twobody_gxtb_lih_shell), &
      new_unittest("gradient-atom-twob-gxtb-m03", test_g_twobody_gxtb_m03_atom), &
      new_unittest("gradient-shell-twob-gxtb-m03", test_g_twobody_gxtb_m03_shell), &
      new_unittest("gradient-atom-twob-gxtb-m04", test_g_twobody_gxtb_m04_atom), &
      new_unittest("gradient-shell-twob-gxtb-m04", test_g_twobody_gxtb_m04_shell), &
      new_unittest("sigma-shell-ons-gfn2-m09", test_s_onsite_gfn2_m09_shell), &
      new_unittest("sigma-atom-ons-gfn2-pyrazine-pbc", test_s_onsite_gfn2_pyrazine_atom), &
      new_unittest("sigma-shell-twob-gxtb-lih", test_s_twobody_gxtb_lih_shell), &
      new_unittest("sigma-atom-twob-gxtb-m05", test_s_twobody_gxtb_m05_atom), &
      new_unittest("sigma-shell-twob-gxtb-m05", test_s_twobody_gxtb_m05_shell), &
      new_unittest("sigma-atom-twob-gxtb-m06", test_s_twobody_gxtb_m06_atom), &
      new_unittest("sigma-shell-twob-gxtb-m06", test_s_twobody_gxtb_m06_shell), &
      new_unittest("potential-gradient-atom-ons-effceh-lih", test_pg_ceh_lih_atom), &
      new_unittest("potential-gradient-atom-ons-effceh-m15", test_pg_ceh_m15_atom), &
      new_unittest("potential-gradient-shell-ons-effceh-m16", test_pg_ceh_m16_shell), &
      new_unittest("potential-sigma-atom-ons-effceh-lih", test_ps_ceh_lih_atom), &
      new_unittest("potential-sigma-atom-ons-effceh-co2", test_ps_ceh_co2_atom), &
      new_unittest("potential-sigma-atom-ons-effceh-m05", test_ps_ceh_m05_atom), &
      new_unittest("potential-sigma-shell-ons-effceh-m17", test_ps_ceh_m17_shell) &
      ]

end subroutine collect_coulomb_thirdorder


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


!> Factory to create onsite thirdorder objects based on GFN1-xTB values
subroutine make_coulomb_o1(coulomb, mol, shell, error)

   !> New coulomb object
   type(tb_coulomb), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: p_hubbard_derivs(20) = 0.1_wp * [&
      & 0.000000_wp, 1.500000_wp, 1.027370_wp, 0.900554_wp, 1.300000_wp, &
      & 1.053856_wp, 0.042507_wp,-0.005102_wp, 1.615037_wp, 1.600000_wp, &
      & 1.200000_wp, 1.100000_wp, 1.200000_wp, 1.500000_wp, 1.500000_wp, &
      & 1.500000_wp, 1.000000_wp, 0.829312_wp, 0.732923_wp, 1.116963_wp] 
   
   real(wp), allocatable :: hubbard_derivs(:, :)
   type(onsite_thirdorder), allocatable :: tmp

   ! Setup coulomb interaction collection 
   allocate(coulomb)
   call new_coulomb(coulomb, mol, error)
   if (allocated(error)) return

   allocate(tmp)
   hubbard_derivs = reshape(p_hubbard_derivs(mol%num), [1, mol%nid])
   call new_onsite_thirdorder(tmp, mol, hubbard_derivs)
   call move_alloc(tmp, coulomb%es3)

end subroutine make_coulomb_o1

!> Factory to create onsite thirdorder objects based on GFN2-xTB values
subroutine make_coulomb_o2(coulomb, mol, shell, error)

   !> New coulomb object
   type(tb_coulomb), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: p_hubbard_derivs(20) = 0.1_wp * [&
      & 0.800000_wp, 2.000000_wp, 1.303821_wp, 0.574239_wp, 0.946104_wp, &
      & 1.500000_wp,-0.639780_wp,-0.517134_wp, 1.426212_wp, 0.500000_wp, &
      & 1.798727_wp, 2.349164_wp, 1.400000_wp, 1.936289_wp, 0.711291_wp, &
      &-0.501722_wp, 1.495483_wp,-0.315455_wp, 2.033085_wp, 2.006898_wp]
   real(wp), parameter :: shell_hubbard_derivs(0:4) = &
      & [1.0_wp, 0.5_wp, spread(0.25_wp, 1, 3)]
   integer, parameter :: shell_count(20) = [&
      & 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]

   real(wp), allocatable :: hubbard_derivs(:, :)
   integer :: isp, izp, ish, il
   type(onsite_thirdorder), allocatable :: tmp

   ! Setup coulomb interaction collection 
   allocate(coulomb)
   call new_coulomb(coulomb, mol, error)
   if (allocated(error)) return

   allocate(tmp)
   if (shell) then
      allocate(hubbard_derivs(3, mol%nid))
      hubbard_derivs(:, :) = 0.0_wp
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            il = ish - 1
            hubbard_derivs(ish, isp) = p_hubbard_derivs(izp) * shell_hubbard_derivs(il)
         end do
      end do
      call new_onsite_thirdorder(tmp, mol, hubbard_derivs, shell_count(mol%num))
   else
      hubbard_derivs = reshape(p_hubbard_derivs(mol%num), [1, mol%nid])
      call new_onsite_thirdorder(tmp, mol, hubbard_derivs)
   end if
   call move_alloc(tmp, coulomb%es3)

end subroutine make_coulomb_o2

!> Factory to create twobody thirdorder objects based on g-xTB values
subroutine make_coulomb_t_gxtb(coulomb, mol, shell, error)

   !> New electrostatic object
   type(tb_coulomb), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Atomic Hubbard parameters or chemical hardness
   real(wp), parameter :: p_hubbard_parameter(60) = [&
      &  0.4725928800_wp,  0.9220339100_wp,  0.1745288800_wp,  0.2570073300_wp, & !1-4
      &  0.3394908600_wp,  0.4219541200_wp,  0.5043819300_wp,  0.5869186300_wp, & !5-8
      &  0.6693135100_wp,  0.7519160700_wp,  0.1796410500_wp,  0.2215727600_wp, & !9-12
      &  0.2634857800_wp,  0.3053964500_wp,  0.3473401400_wp,  0.3892472500_wp, & !13-16
      &  0.4311567000_wp,  0.4730826900_wp,  0.1710546900_wp,  0.2027624400_wp, & !17-20
      &  0.2100732200_wp,  0.2173964700_wp,  0.2247103900_wp,  0.2320150100_wp, & !21-24
      &  0.2393396900_wp,  0.2466563800_wp,  0.2539825500_wp,  0.2612886300_wp, & !25-28
      &  0.2685947600_wp,  0.2759256500_wp,  0.3076299900_wp,  0.3393158000_wp, & !29-32
      &  0.3723598500_wp,  0.4027354900_wp,  0.4344577600_wp,  0.4661170800_wp, & !33-36
      &  0.1558507900_wp,  0.1864932400_wp,  0.1935621000_wp,  0.2006331100_wp, & !37-40
      &  0.2077052200_wp,  0.2147725400_wp,  0.2218461400_wp,  0.2289187200_wp, & !41-44
      &  0.2359862100_wp,  0.2430561200_wp,  0.2501301800_wp,  0.2571993700_wp, & !45-48
      &  0.2878478000_wp,  0.3184867300_wp,  0.3491243100_wp,  0.3797659300_wp, & !49-52
      &  0.4104080800_wp,  0.4410577700_wp,  0.0501933200_wp,  0.0676257000_wp, & !53-56
      &  0.0850444500_wp,  0.1024773600_wp,  0.1199110500_wp,  0.1373277200_wp] !57-60

   !> Shell-scaling for the Hubbard parameters or chemical hardness
   !> for second-order tight-binding
   real(wp), parameter :: p_shell_hubbard(4, 60) = reshape([&
      &  0.9464169217_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !1
      &  1.0469495122_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !2
      &  0.8132756349_wp,  0.8805765138_wp,  0.0000000000_wp,  0.0000000000_wp, & !3
      &  1.0492612127_wp,  1.0475558484_wp,  0.0000000000_wp,  0.0000000000_wp, & !4
      &  0.9555440530_wp,  0.9716757622_wp,  0.0000000000_wp,  0.0000000000_wp, & !5
      &  1.0165874455_wp,  0.8821434392_wp,  0.0000000000_wp,  0.0000000000_wp, & !6
      &  0.9726072270_wp,  0.7921417667_wp,  0.0000000000_wp,  0.0000000000_wp, & !7
      &  1.1369464012_wp,  0.7918709113_wp,  0.0000000000_wp,  0.0000000000_wp, & !8
      &  1.0960193236_wp,  0.8299693123_wp,  0.0000000000_wp,  0.0000000000_wp, & !9
      &  1.6207110415_wp,  0.7932001047_wp,  0.0000000000_wp,  0.0000000000_wp, & !10
      &  0.8662947624_wp,  0.7789473893_wp,  0.0000000000_wp,  0.0000000000_wp, & !11
      &  1.1039691910_wp,  0.9056377427_wp,  1.0816848716_wp,  0.0000000000_wp, & !12
      &  1.2307657286_wp,  0.8989571619_wp,  1.6068779370_wp,  0.0000000000_wp, & !13
      &  1.0396148307_wp,  0.8349303125_wp,  0.9359597426_wp,  0.0000000000_wp, & !14
      &  1.1951393269_wp,  0.8488074903_wp,  0.8472579619_wp,  0.0000000000_wp, & !15
      &  0.8895919837_wp,  0.7986074171_wp,  0.6518407631_wp,  0.0000000000_wp, & !16
      &  1.2766923782_wp,  0.8474086558_wp,  0.9144929669_wp,  0.0000000000_wp, & !17
      &  1.3146553673_wp,  0.7842427139_wp,  0.2000000030_wp,  0.0000000000_wp, & !18
      &  0.8700795322_wp,  0.6906432456_wp,  0.0000000000_wp,  0.0000000000_wp, & !19
      &  1.0256912887_wp,  1.3920647067_wp,  0.7954270392_wp,  0.0000000000_wp, & !20
      &  1.2541646765_wp,  0.6017316709_wp,  1.5481334140_wp,  0.0000000000_wp, & !21
      &  1.3457526976_wp,  0.8171693594_wp,  1.3211861774_wp,  0.0000000000_wp, & !22
      &  1.3807837037_wp,  0.8137856445_wp,  1.3478527812_wp,  0.0000000000_wp, & !23
      &  0.9795461354_wp,  0.5485797738_wp,  1.2682366602_wp,  0.0000000000_wp, & !24
      &  1.1032546994_wp,  0.6582642385_wp,  1.6730571933_wp,  0.0000000000_wp, & !25
      &  1.0439120723_wp,  0.7622240663_wp,  1.7529923287_wp,  0.0000000000_wp, & !26
      &  0.9780845068_wp,  0.7229254198_wp,  1.9313123159_wp,  0.0000000000_wp, & !27
      &  0.9318598687_wp,  0.5098226174_wp,  1.6893501712_wp,  0.0000000000_wp, & !28
      &  1.1018736360_wp,  0.4264834619_wp,  1.8212779865_wp,  0.0000000000_wp, & !29
      &  0.9758721814_wp,  0.7686020673_wp,  0.0000000000_wp,  0.0000000000_wp, & !30
      &  0.9821762645_wp,  0.8219884850_wp,  0.6082808941_wp,  0.0000000000_wp, & !31
      &  0.8553022267_wp,  0.7640261264_wp,  0.5034949661_wp,  0.0000000000_wp, & !32
      &  1.0135146878_wp,  0.7475060062_wp,  0.5608727667_wp,  0.0000000000_wp, & !33
      &  0.8733342408_wp,  0.7173571104_wp,  0.6343295211_wp,  0.0000000000_wp, & !34
      &  0.8762328090_wp,  0.7604833639_wp,  0.7589458711_wp,  0.0000000000_wp, & !35
      &  1.0268812360_wp,  0.6974261439_wp,  0.6627860654_wp,  0.0000000000_wp, & !36
      &  1.0586849132_wp,  0.8876719775_wp,  0.0000000000_wp,  0.0000000000_wp, & !37
      &  1.0125740729_wp,  0.7734175839_wp,  1.1935739696_wp,  0.0000000000_wp, & !38
      &  1.0861908792_wp,  0.8742022235_wp,  1.0183786400_wp,  0.0000000000_wp, & !39
      &  1.2760300829_wp,  1.0377723675_wp,  1.1353505307_wp,  0.0000000000_wp, & !40
      &  1.2990721726_wp,  0.9534014229_wp,  1.1862665906_wp,  0.0000000000_wp, & !41
      &  1.1828449223_wp,  0.9924100339_wp,  1.3028844186_wp,  0.0000000000_wp, & !42
      &  1.2457763061_wp,  0.5758137723_wp,  1.3884996994_wp,  0.0000000000_wp, & !43
      &  1.0864993879_wp,  0.6880088096_wp,  1.5587710441_wp,  0.0000000000_wp, & !44
      &  1.1150512984_wp,  0.6189764277_wp,  1.6268976268_wp,  0.0000000000_wp, & !45
      &  1.1468284449_wp,  0.6084583370_wp,  1.5080051732_wp,  0.0000000000_wp, & !46
      &  1.2974499968_wp,  0.5935388618_wp,  1.3043870339_wp,  0.0000000000_wp, & !47
      &  1.1245569530_wp,  0.6938172206_wp,  0.0000000000_wp,  0.0000000000_wp, & !48
      &  1.0778327303_wp,  0.7596073011_wp,  0.4508068027_wp,  0.0000000000_wp, & !49
      &  0.8971500149_wp,  0.7210165892_wp,  0.4526771713_wp,  0.0000000000_wp, & !50
      &  0.7975020744_wp,  0.6730867257_wp,  0.3907781437_wp,  0.0000000000_wp, & !51
      &  0.9766829398_wp,  0.6536945179_wp,  0.4610940105_wp,  0.0000000000_wp, & !52
      &  0.8814212108_wp,  0.7125909325_wp,  0.8644415285_wp,  0.0000000000_wp, & !53
      &  2.0273685645_wp,  0.6617072240_wp,  0.4249197117_wp,  0.0000000000_wp, & !54
      &  2.0761614803_wp,  4.4214513738_wp,  0.0000000000_wp,  0.0000000000_wp, & !55
      &  2.1787168814_wp,  1.8457906962_wp,  2.6311808555_wp,  0.0000000000_wp, & !56
      &  2.1642824902_wp,  1.3006449626_wp,  2.3535134678_wp,  0.0000000000_wp, & !57
      &  1.8443014730_wp,  0.8313487372_wp,  2.0414076297_wp,  5.3670117853_wp, & !58
      &  1.6171183532_wp,  0.5541990090_wp,  1.9983544977_wp,  4.8618221978_wp, & !59
      &  1.1642363075_wp,  0.6557853312_wp,  1.5579302981_wp,  4.9536261223_wp],& !60
      & shape(p_shell_hubbard))

   !> CN-dependence of the Hubbard parameter for second-order tight-binding
   real(wp), parameter :: pa_tb2_hubbard_cn(60) = [&
      &  1.0227655971_wp,  1.8268631210_wp,  0.8833851822_wp, -0.0697354595_wp, & !1-4
      &  0.0320241836_wp,  0.1592026095_wp,  0.1912729298_wp,  0.0736673324_wp, & !5-8
      & -0.0543136451_wp,  0.2498562442_wp,  0.9039947639_wp, -0.0863905029_wp, & !9-12
      & -0.0159996132_wp,  0.2181798725_wp,  0.0013271661_wp,  0.0808688669_wp, & !13-16
      & -0.0819748983_wp, -0.3436711724_wp,  1.0577464191_wp, -0.3420185483_wp, & !17-20
      & -0.0016736512_wp, -0.0044604372_wp, -0.0108927366_wp,  0.3828669208_wp, & !21-24
      &  0.0654548718_wp,  0.0267336071_wp,  0.0570803443_wp,  0.0863538344_wp, & !25-28
      &  0.2525728488_wp, -0.0670798662_wp,  0.1802824320_wp,  0.3583183378_wp, & !29-32
      & -0.0016231720_wp, -0.0510806043_wp,  0.1862622360_wp, -0.3873812217_wp, & !33-36
      &  0.3069504766_wp,  0.0813168350_wp,  0.0433573692_wp,  0.0199935683_wp, & !37-40
      &  0.1177473201_wp,  0.0260182660_wp,  0.0335642042_wp,  0.0355746241_wp, & !41-44
      &  0.0549381766_wp,  0.0532426334_wp, -0.0037690605_wp, -0.0117487916_wp, & !45-48
      &  0.0489277980_wp,  0.3280168961_wp,  0.2755642779_wp,  0.0188372641_wp, & !49-52
      &  0.2193034068_wp, -0.1700364031_wp,  0.5024281804_wp,  0.1359279318_wp, & !53-56
      &  0.2392737272_wp,  0.0328019344_wp,  0.1333449065_wp,  0.0444103029_wp] !57-60

   !> Parameter: Derivative of Hubbard parameters or chemical hardness 
   !> for third-order tight-binding
   real(wp), parameter :: p_hubbard_derivs(60) = [&
      & -0.4000000000_wp, -0.8500000000_wp, -0.2000000000_wp, -0.4400000000_wp, & !1-4
      & -0.4500000000_wp, -0.4800000000_wp, -0.5000000000_wp, -0.5800000000_wp, & !5-8
      & -0.6600000000_wp, -0.8000000000_wp, -0.2000000000_wp, -0.3000000000_wp, & !9-12
      & -0.3000000000_wp, -0.3200000000_wp, -0.3000000000_wp, -0.5800000000_wp, & !13-16
      & -0.5800000000_wp, -0.6000000000_wp, -0.2000000000_wp, -0.3000000000_wp, & !17-20
      & -0.4000000000_wp, -0.4000000000_wp, -0.4000000000_wp, -0.4000000000_wp, & !21-24
      & -0.3500000000_wp, -0.4500000000_wp, -0.4500000000_wp, -0.4500000000_wp, & !25-28
      & -0.4500000000_wp, -0.4500000000_wp, -0.3000000000_wp, -0.3000000000_wp, & !29-32
      & -0.3500000000_wp, -0.4000000000_wp, -0.4500000000_wp, -0.5500000000_wp, & !33-36
      & -0.2000000000_wp, -0.2000000000_wp, -0.3500000000_wp, -0.4000000000_wp, & !37-40
      & -0.4000000000_wp, -0.3500000000_wp, -0.3500000000_wp, -0.3500000000_wp, & !41-44
      & -0.3500000000_wp, -0.4000000000_wp, -0.4500000000_wp, -0.4500000000_wp, & !45-48
      & -0.3000000000_wp, -0.2500000000_wp, -0.3000000000_wp, -0.4000000000_wp, & !49-52
      & -0.4500000000_wp, -0.5000000000_wp, -0.2000000000_wp, -0.3000000000_wp, & !53-56
      & -0.3000000000_wp, -0.3500000000_wp, -0.3500000000_wp, -0.3500000000_wp] !57-60

   !> Shell-dependence of the third-order tight-binding
   real(wp), parameter :: p_tb3_kshell(0:3) = &
      & [1.0000000000_wp, 1.1500000000_wp, 1.3000000000_wp, 1.4500000000_wp]

   !> Coulomb-kernel exponent for second-order tight-binding
   real(wp), parameter :: p_gexp = 1.0000000000_wp
   
   !> Exponent in radius-dependent hubbard scaling for second-order tight-binding
   real(wp), parameter :: p_hubbard_exp = 0.2946211550_wp

   !> Coordination number exponent
   real(wp), parameter :: p_cn_exp = 2.0680000000_wp

   !> Exponent and expansion coefficient for third-order tight-binding
   real(wp), parameter :: p_texp = 1.3000000000_wp
   real(wp), parameter :: p_onsite = 0.2093327496_wp
   real(wp), parameter :: p_offsite = 2.3000000000_wp

   !> Number of shells selected from the q-vSZP basis set
   integer, parameter :: p_nshell(60) = [ &
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, & !1-20
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, & !21-40
      & 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 4, 4, 4] !41-60
   
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

   integer :: isp, izp, ish
   integer, allocatable :: ish_at(:)
   real(wp), allocatable :: hubbard_derivs(:, :), hubbard(:, :), hardness_cn(:)
   real(wp), allocatable :: tb_cn_rcov(:)
   type(twobody_thirdorder), allocatable :: tmp
   type(average_type), allocatable :: average

   ! Setup coulomb interaction collection 
   tb_cn_rcov = p_cn_rcov(mol%num)
   allocate(coulomb)
   call new_coulomb(coulomb, mol, error, cn_count_type=cn_count%erf, &
      & cn_rcov=tb_cn_rcov, cn_exp=p_cn_exp)
   if (allocated(error)) return

   ! Setup effective coulomb interaction
   allocate(tmp)
   if (shell) then
      allocate(hubbard(4, mol%nid), hubbard_derivs(4, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, p_nshell(izp)
            hubbard(ish, isp) = p_hubbard_parameter(izp) * p_shell_hubbard(ish, izp)
            hubbard_derivs(ish, isp) = p_tb3_kshell(ish-1) * p_hubbard_derivs(izp)
         end do
      end do
      hardness_cn = pa_tb2_hubbard_cn(mol%num)
      allocate(average)
      call new_average(average, average_id%harmonic)
      call new_twobody_thirdorder(tmp, mol, p_texp, p_onsite, p_offsite, &
         & hubbard, average, hubbard_derivs, hardness_cn, p_nshell(mol%num))
   else
      allocate(hubbard(1, mol%nid), hubbard_derivs(1, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         hubbard(1, isp) = p_hubbard_parameter(izp)
         hubbard_derivs(1, isp) = p_hubbard_derivs(izp)
      end do
      hardness_cn = pa_tb2_hubbard_cn(mol%num)
      allocate(average)
      call new_average(average, average_id%harmonic)
      call new_twobody_thirdorder(tmp, mol, p_texp,  p_onsite, p_offsite, &
         & hubbard, average, hubbard_derivs, hardness_cn)
   end if
   call move_alloc(tmp, coulomb%es3)

end subroutine make_coulomb_t_gxtb

!> Factory to create onsite thirdorder objects based on CEH values
subroutine make_coulomb_oceh(coulomb, mol, shell, error)

   !> New coulomb object
   type(tb_coulomb), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: p_hubbard_derivs(20) = [&
      &  0.8936213810_wp, -0.3936567743_wp, -0.7726174171_wp, -0.2849896764_wp, &
      &  0.0126634714_wp, -0.0082561791_wp,  0.0992949802_wp, -0.0267387652_wp, &
      & -0.0632999086_wp, -1.0106414497_wp, -0.3492075197_wp, -0.3191627473_wp, &
      &  0.0467483747_wp, -0.0920002125_wp, -0.0728788864_wp, -0.0213909690_wp, &
      & -0.0206065548_wp, -0.0432378706_wp, -0.0686554093_wp, -0.1671301006_wp] 
   integer, parameter :: shell_count(20) = [&
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]

   real(wp), allocatable :: hubbard_derivs(:, :)
   integer :: isp, izp, ish, il
   type(onsite_thirdorder), allocatable :: tmp

   ! Setup coulomb interaction collection 
   allocate(coulomb)
   call new_coulomb(coulomb, mol, error)
   if (allocated(error)) return

   allocate(tmp)
   if (shell) then
      ! If shell-resolved, we use the atomic parameter for each shell
      allocate(hubbard_derivs(3, mol%nid))
      hubbard_derivs(:, :) = 0.0_wp
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            il = ish - 1
            hubbard_derivs(ish, isp) = p_hubbard_derivs(izp)
         end do
      end do
      call new_onsite_thirdorder(tmp, mol, hubbard_derivs, shell_count(mol%num))
   else
      hubbard_derivs = reshape(p_hubbard_derivs(mol%num), [1, mol%nid])
      call new_onsite_thirdorder(tmp, mol, hubbard_derivs)
   end if
   call move_alloc(tmp, coulomb%es3)

end subroutine make_coulomb_oceh


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

   !> Factory to create new onsite thirdorder objects
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
   call coulomb%es3%get_energy(mol, cache, wfn, energy)

   call check(error, sum(energy), ref, thr=thr_)
   if (allocated(error)) then
      print*, ref, sum(energy), ref - sum(energy)
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
      call coulomb%es3%get_energy(mol, cache, wfn, er)
      ! Left hand side
      wfn%qat(iat, 1) = wfn%qat(iat, 1) - 2*step
      call coulomb%es3%get_energy(mol, cache, wfn, el)

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
         call coulomb%es3%get_energy(mol, cache, wfn, er)
         ! Left hand side
         wfn%qsh(ish, 1) = wfn%qsh(ish, 1) - 2*step
         call coulomb%es3%get_energy(mol, cache, wfn, el)

         wfn%qsh(ish, 1) = wfn%qsh(ish, 1) + step
         numvsh(ish) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end if

   ! Analytic potentials
   call pot%reset()
   call coulomb%es3%get_potential(mol, cache, wfn, pot)

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


subroutine test_taumat_numgrad(error, mol, qat, qsh, make_coulomb, thr_in)

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

   type(tb_coulomb), allocatable, target :: coulomb
   type(container_cache) :: cache
   type(coulomb_cache), pointer :: ptr
   type(twobody_thirdorder), pointer :: es3
   integer :: iat, ic, ndim
   real(wp), allocatable :: dtaudr(:, :, :), dtaudL(:, :, :)
   real(wp), allocatable :: dtaudcn(:, :), qvec(:)
   real(wp), allocatable :: taumatr(:, :), taumatl(:, :), numdtaudr(:, :, :)
   real(wp), allocatable :: dtau(:, :), dtauq(:)
   real(wp), parameter :: step = 1.0e-5_wp
   real(wp) :: thr_

   thr_ = thr
   if (present(thr_in)) thr_ = thr_in

   ! Select charges
   if (present(qsh)) then
      ndim = size(qsh, 1)
      allocate(qvec(ndim))
      qvec(:) = qsh
   else
      ndim = size(qat, 1)
      allocate(qvec(ndim))
      qvec(:) = qat
   end if

   ! Setup coulomb object
   call make_coulomb(coulomb, mol, present(qsh), error)
   if (allocated(error)) return

   ! Check for thirdorder object with tau dependence
   nullify(es3)
   select type(target => coulomb%es3)
   type is(twobody_thirdorder)
      es3 => target
   end select
   if (.not. associated(es3)) then
      call test_failed(error, "Coulomb object is not of twobody thirdorder type")
      return
   end if

   call coulomb%update(mol, cache)
   
   ! Numerical derivative of the Coulomb matrix
   allocate(dtaudr(3, mol%nat, ndim), numdtaudr(3, mol%nat, ndim), &
      & dtaudL(3, 3, ndim), dtaudcn(ndim, ndim), dtau(ndim, ndim), dtauq(ndim))
   dtaudr(:, :, :) = 0.0_wp 
   numdtaudr(:, :, :) = 0.0_wp
   allocate(taumatr(ndim, ndim), taumatl(ndim, ndim), source=0.0_wp)
   do iat = 1, mol%nat
      do ic = 1, 3
         ! Right hand side
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call coulomb%es3%update(mol, cache)
         call view(cache, ptr)
         taumatr(:, :) = ptr%taumat(:, :)

         ! Left hand side
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call coulomb%es3%update(mol, cache)
         call view(cache, ptr)
         taumatl(:, :) = ptr%taumat(:, :)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step

         ! Full matrix derivative since third order is not bilinear
         dtau(:, :) = (taumatr - taumatl) * (0.5_wp/step)
         call gemv(dtau, qvec, dtauq)
         numdtaudr(ic, iat, :) = dtauq(:)
      end do
   end do

   ! Analytic derivative for Coulomb matrix
   call coulomb%es3%update(mol, cache)
   call view(cache, ptr)
   call es3%get_tau_derivs(mol, ptr, qat, qsh, dtaudr, dtaudL, dtaudcn)

   if (any(abs(dtaudr - numdtaudr) > thr_)) then
      call test_failed(error, "Gradient of tau matrix does not match")
      write(*,*) "numerical gradient:"
      print'(3es21.14)', numdtaudr
      write(*,*) "analytical gradient:"
      print'(3es21.14)', dtaudr
      write(*,*) "difference:"
      print'(3es21.14)', dtaudr - numdtaudr
   end if

end subroutine test_taumat_numgrad

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


subroutine test_numgrad(error, mol, qat, qsh, make_coulomb, thr_in)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges for this structure
   real(wp), intent(in), optional :: qsh(:)

   !> Factory to create new onsite thirdorder objects
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

   call coulomb%update(mol, cache)

   do iat = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call coulomb%update(mol, cache)
         call coulomb%es3%get_energy(mol, cache, wfn, er)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call coulomb%update(mol, cache)
         call coulomb%es3%get_energy(mol, cache, wfn, el)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call coulomb%update(mol, cache)
   call coulomb%es3%get_gradient(mol, cache, wfn, gradient, sigma)

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

   !> Factory to create new onsite thirdorder objects
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
   real(wp), allocatable :: numpotgrad(:, :, :)
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
      allocate(numpotgrad(3, mol%nat, bas%nsh), source=0.0_wp)
   else
      allocate(numpotgrad(3, mol%nat, mol%nat), source=0.0_wp)
   end if
   
   do iat = 1, mol%nat
      do ic = 1, 3
         call potr%reset
         call potl%reset
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_charges(wfn, mol, coulomb%es3%nsh_at, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%es3%get_potential(mol, cache, wfn, potr)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_charges(wfn, mol, coulomb%es3%nsh_at, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%es3%get_potential(mol, cache, wfn, potl)
         
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         if(shell) then
            numpotgrad(ic, iat, :) = 0.5_wp*(potr%vsh(:,1) - potl%vsh(:,1))/step
         else
            numpotgrad(ic, iat, :) = 0.5_wp*(potr%vat(:,1) - potl%vat(:,1))/step
         end if
      end do
   end do

   call get_charges(wfn, mol, coulomb%es3%nsh_at, error)
   if (allocated(error)) return
   call coulomb%update(mol, cache)
   call coulomb%es3%get_potential(mol, cache, wfn, potl)
   call coulomb%es3%get_potential_gradient(mol, cache, wfn, potl)

   if(shell) then
      if (any(abs(potl%dvshdr(:,:,:,1) - numpotgrad) > thr_)) then
         call test_failed(error, "Gradient of shell-resolved potential does not match")
         print'(3es21.14)', potl%dvshdr(:,:,:,1) - numpotgrad
      end if
   else
      if (any(abs(potl%dvatdr(:,:,:,1) - numpotgrad) > thr_)) then
         call test_failed(error, "Gradient of atom-resolved potential does not match")
         print'(3es21.14)', potl%dvatdr(:,:,:,1) - numpotgrad
      end if
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

   call coulomb%update(mol, cache)

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
         call coulomb%es3%get_energy(mol, cache, wfn, er)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call coulomb%update(mol, cache)
         call coulomb%es3%get_energy(mol, cache, wfn, el)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (allocated(lattice)) mol%lattice = lattice
         numsigma(jc, ic) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call coulomb%update(mol, cache)
   call coulomb%es3%get_gradient(mol, cache, wfn, gradient, sigma)

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
   real(wp), allocatable :: numpotsigma(:, :, :), xyz(:, :), lattice(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 5.0e-5_wp
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
      allocate(numpotsigma(3, 3, bas%nsh), source=0.0_wp)
   else
      allocate(numpotsigma(3, 3, mol%nat), source=0.0_wp)
   end if

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
         call get_charges(wfn, mol, coulomb%es3%nsh_at, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%es3%get_potential(mol, cache, wfn, potr)

         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call get_charges(wfn, mol, coulomb%es3%nsh_at, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%es3%get_potential(mol, cache, wfn, potl)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (allocated(lattice)) mol%lattice = lattice
         if(shell) then
            numpotsigma(jc, ic, :) = 0.5_wp*(potr%vsh(:,1) - potl%vsh(:,1))/step
         else
            numpotsigma(jc, ic, :) = 0.5_wp*(potr%vat(:,1) - potl%vat(:,1))/step
         end if
      end do
   end do

   call get_charges(wfn, mol, coulomb%es3%nsh_at, error)
   if (allocated(error)) return
   call coulomb%update(mol, cache)
   call coulomb%es3%get_potential(mol, cache, wfn, potl)
   call coulomb%es3%get_potential_gradient(mol, cache, wfn, potl)

   if(shell) then
      if (any(abs(potl%dvshdL(:,:,:,1) - numpotsigma) > thr_)) then
         call test_failed(error, "Sigma of shell-resolved potential does not match")
         print'(3es21.14)', potl%dvshdL(:,:,:,1) - numpotsigma
      end if
   else
      if (any(abs(potl%dvatdL(:,:,:,1) - numpotsigma) > thr_)) then
         call test_failed(error, "Sigma of atom-resolved potential does not match")
         print'(3es21.14)', potl%dvatdL(:,:,:,1) - numpotsigma
      end if
   end if

end subroutine test_numpotsigma


subroutine test_e_onsite_gfn1_m01_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 0.80928533739041852_wp, -0.0985788702919893_wp, -1.1789498512968015_wp, &
      &-0.07630466864269804_wp, -0.5449852206782641_wp,  0.3220779574144627_wp, &
      &-0.02986669797918212_wp, -1.1079458022741191_wp, -0.6616262018005230_wp, &
      & 0.56887324141215356_wp,  0.3231023443889342_wp,  0.0959376082660442_wp, &
      & 0.27886474952257534_wp,  0.8972931001216102_wp, -0.3079177622577322_wp, &
      & 0.71074073670510574_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, qat, qsh, make_coulomb_o1, 5.9181986430716059E-2_wp)

end subroutine test_e_onsite_gfn1_m01_atom

subroutine test_e_onsite_gfn2_m02_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.38394711236234E-2_wp,-1.68354976558608E-1_wp,-3.47642833746823E-1_wp, &
      &-7.05489267186003E-1_wp, 7.73548301641266E-1_wp, 2.30207581365386E-1_wp, &
      & 1.02748501676354E-1_wp, 9.47818107467040E-2_wp, 2.44260351729187E-2_wp, &
      & 2.34984927037408E-1_wp,-3.17839896393030E-1_wp, 6.67112994818879E-1_wp, &
      &-4.78119977010488E-1_wp, 6.57536027459275E-2_wp, 1.08259054549882E-1_wp, &
      &-3.58215329983396E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, qat, qsh, make_coulomb_o2, 4.9292337881703292E-2_wp)

end subroutine test_e_onsite_gfn2_m02_atom

subroutine test_e_onsite_gfn2_m07_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.986820554251878_wp,  0.2133150304286247_wp,  1.2219537941124052_wp, &
      & 0.637526091539623_wp,  0.6288673640786957_wp, -0.4793698954628428_wp, &
      &-0.381215207423406_wp, -1.0863256053671664_wp,  0.8906955087366708_wp, &
      & 0.159224917897852_wp,  0.9937121523714770_wp,  0.4332827430406178_wp, &
      & 0.459237705248708_wp, -0.6369517649821776_wp, -0.3853972301505781_wp, &
      & 0.318264950183372_wp]
   real(wp), parameter :: qsh(*) = [&
      &-0.2372082175033831_wp, -2.7496123367484957_wp,  0.2133150304286247_wp, &
      & 1.2459662589528724_wp, -0.0240124648404671_wp,  0.6375260915396232_wp, &
      & 0.6288673640786957_wp,  0.0427487776207620_wp, -0.4119239254478257_wp, &
      &-0.1101947476357790_wp,  0.1962526334754219_wp, -0.5774678408988283_wp, &
      &-0.0954706513992993_wp, -0.9908549539678670_wp,  0.1217211077680937_wp, &
      & 0.7689744009685770_wp,  0.1592249178978529_wp,  0.2104607013244239_wp, &
      & 0.7889318241389675_wp, -0.0056803730919144_wp,  0.4332827430406178_wp, &
      & 0.4592377052487081_wp,  0.2833371093194456_wp, -0.9202888743016233_wp, &
      & 0.1737487393096157_wp, -0.5591459694601939_wp,  1.0261224370878024_wp, &
      &-0.3777382313598706_wp, -0.3301192555445596_wp]
   real(wp), allocatable :: qsh0(:)

   call get_structure(mol, "MB16-43", "07")
   call test_generic(error, mol, qat, qsh0, make_coulomb_o2, -1.2136915157280399_wp)
   call test_generic(error, mol, qat, qsh, make_coulomb_o2, -0.34488498240315002_wp)

end subroutine test_e_onsite_gfn2_m07_shell

subroutine test_e_onsite_gfn2_oxacb_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
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
   call test_generic(error, mol, qat, qsh, make_coulomb_o2, 0.10439618995869278_wp)

end subroutine test_e_onsite_gfn2_oxacb_shell

subroutine test_e_onsite_gfn2_oxacb_sc_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat1(*) = [&
      & 3.41731844312030E-1_wp, 3.41716020106239E-1_wp, 3.41730526585671E-1_wp, &
      & 3.41714427217954E-1_wp, 3.80996046757999E-1_wp, 3.80989821246195E-1_wp, &
      & 3.81000747720282E-1_wp, 3.80990494183703E-1_wp,-3.70406587264474E-1_wp, &
      &-3.70407565207006E-1_wp,-3.70417590212352E-1_wp,-3.70399716470705E-1_wp, &
      &-3.52322260586075E-1_wp,-3.52304269439196E-1_wp,-3.52313440903261E-1_wp, &
      &-3.52298498047004E-1_wp]
   integer, parameter :: supercell(*) = [2, 2, 2]
   real(wp), parameter :: qat(*) = [spread(qat1, 2, product(supercell))]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "oxacb")
   call make_supercell(mol, supercell)
   call test_generic(error, mol, qat, qsh, make_coulomb_o2, &
      & 0.02183660722324880_wp*product(supercell), 1.0e-7_wp)

end subroutine test_e_onsite_gfn2_oxacb_sc_atom

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


subroutine test_e_twobody_gxtb_lih_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  3.903043569209723E-1_wp, -3.903043569695019E-1_wp]
   real(wp), parameter :: qsh(*) = [ &
      &  1.883240895671248E-1_wp,  2.019802673538475E-1_wp, -3.903043569695019E-1_wp]

   call get_structure(mol, "MB16-43", "LiH")
   call test_generic(error, mol, qat, qsh, make_coulomb_t_gxtb, &
      & -4.518092516045460E-4_wp, thr_in=thr)

end subroutine test_e_twobody_gxtb_lih_shell

subroutine test_e_twobody_gxtb_m01_atom(error)

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
   call test_generic(error, mol, qat, qsh, make_coulomb_t_gxtb, &
      & -1.1878845985344208E-2_wp, thr_in=thr)

end subroutine test_e_twobody_gxtb_m01_atom

subroutine test_e_twobody_gxtb_m01_shell(error)

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
   call test_generic(error, mol, qat, qsh, make_coulomb_t_gxtb, &
      & -4.386569688922437E-3_wp, thr_in=thr)

end subroutine test_e_twobody_gxtb_m01_shell

subroutine test_e_twobody_gxtb_m02_atom(error)

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
   call test_generic(error, mol, qat, qsh, make_coulomb_t_gxtb, &
      & -4.3093926507609508E-3_wp, thr_in=thr)

end subroutine test_e_twobody_gxtb_m02_atom

subroutine test_e_twobody_gxtb_m02_shell(error)

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
   call test_generic(error, mol, qat, qsh, make_coulomb_t_gxtb, &
      & -2.490306308579771E-3_wp, thr_in=thr)

end subroutine test_e_twobody_gxtb_m02_shell


subroutine test_p_twobody_gxtb_lih_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  3.903043569209723E-1_wp, -3.903043569695019E-1_wp]
   real(wp), parameter :: qsh(*) = [ &
      &  1.883240895671248E-1_wp,  2.019802673538475E-1_wp, -3.903043569695019E-1_wp]

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpot(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_p_twobody_gxtb_lih_shell

subroutine test_p_twobody_gxtb_m01_atom(error)

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
   call test_numpot(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_p_twobody_gxtb_m01_atom

subroutine test_p_twobody_gxtb_m01_shell(error)

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
   call test_numpot(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_p_twobody_gxtb_m01_shell

subroutine test_p_twobody_gxtb_m02_atom(error)

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
   call test_numpot(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_p_twobody_gxtb_m02_atom

subroutine test_p_twobody_gxtb_m02_shell(error)

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
   call test_numpot(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_p_twobody_gxtb_m02_shell


subroutine test_taumat_twobody_gxtb_h2_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      & 1.0_wp, 1.0_wp]
   real(wp), parameter :: qsh(*) = [ &
      & 1.0_wp, 1.0_wp]
   
   call get_structure(mol, "MB16-43", "H2")
   call test_taumat_numgrad(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1*10)

end subroutine test_taumat_twobody_gxtb_h2_shell

subroutine test_taumat_twobody_gxtb_lih_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  3.903043569209723E-1_wp, -3.903043569695019E-1_wp]
   real(wp), parameter :: qsh(*) = [ &
      &  1.883240895671248E-1_wp,  2.019802673538475E-1_wp, -3.903043569695019E-1_wp]

   call get_structure(mol, "MB16-43", "LiH")
   call test_taumat_numgrad(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_taumat_twobody_gxtb_lih_shell

subroutine test_taumat_twobody_gxtb_s2_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      & -4.234713968376269E-11_wp, -4.241411388772320E-11_wp]
   real(wp), parameter :: qsh(*) = [ &
      & -1.658294650172361E-1_wp,  8.500381954781311E-2_wp,  8.082564542707581E-2_wp, &
      & -1.658294650172407E-1_wp,  8.500381954775360E-2_wp,  8.082564542707300E-2_wp]

   call get_structure(mol, "MB16-43", "S2")
   call test_taumat_numgrad(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_taumat_twobody_gxtb_s2_shell

subroutine test_taumat_twobody_gxtb_cecl3_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  1.761677281792140E+0_wp, -5.868465765187969E-1_wp, -5.882941651185237E-1_wp, &
      & -5.865365403069231E-1_wp]
   real(wp), parameter :: qsh(*) = [ &
      &  8.561147975329674E-1_wp,  1.080811850518806E-1_wp,  8.401587265651704E-1_wp, &
      & -4.267742735787827E-2_wp, -9.508455719506892E-2_wp, -5.372626500963689E-1_wp, &
      &  4.550063077264088E-2_wp, -9.613010319135307E-2_wp, -5.376894880457614E-1_wp, &
      &  4.552542611859081E-2_wp, -9.547301109487960E-2_wp, -5.365766770437510E-1_wp, &
      &  4.551314783170746E-2_wp]

   call get_structure(mol, "f-block", "CeCl3")
   call test_taumat_numgrad(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_taumat_twobody_gxtb_cecl3_shell


subroutine test_g_onsite_gfn1_m03_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.77788256288236E-1_wp,-8.22943267808161E-1_wp, 4.04578389873281E-2_wp,&
      & 5.79710531992282E-1_wp, 6.99601887637659E-1_wp, 6.84309612639107E-2_wp,&
      &-3.42971414989811E-1_wp, 4.64954031865410E-2_wp, 6.77012204116428E-2_wp,&
      & 8.49931225363225E-2_wp,-5.22285304699699E-1_wp,-2.92515001764712E-1_wp,&
      &-3.98375452377043E-1_wp, 2.09769668297792E-1_wp, 7.23140464830357E-1_wp,&
      & 3.65775987838250E-2_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "03")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_o1)

end subroutine test_g_onsite_gfn1_m03_atom

subroutine test_g_onsite_gfn1_urea_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 5.55723890858218E-1_wp, 5.55765354442035E-1_wp, 2.50200231242017E-1_wp,&
      & 2.50282053284422E-1_wp, 2.39786980460652E-1_wp, 2.39895142481200E-1_wp,&
      & 2.50103678240412E-1_wp, 2.50425041601730E-1_wp, 2.39464477136495E-1_wp,&
      & 2.40360053062669E-1_wp,-4.38369096728919E-1_wp,-4.38451412936599E-1_wp,&
      &-4.38310020776279E-1_wp,-4.38617373848238E-1_wp,-6.59141030224988E-1_wp,&
      &-6.59117968294813E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "urea")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_o1)

end subroutine test_g_onsite_gfn1_urea_atom

subroutine test_g_onsite_gfn2_m04_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 9.33596160193497E-2_wp,-3.41088061922851E-1_wp, 7.32474961830646E-2_wp,&
      &-2.21649975471802E-1_wp, 6.24413528413759E-3_wp, 1.07366683260668E-1_wp,&
      & 1.25982547197317E-1_wp, 9.65935501843890E-2_wp, 1.02704543049803E-1_wp,&
      & 1.45380937882263E-1_wp,-1.55978251071729E-1_wp, 3.42948437914661E-1_wp,&
      & 5.65504846503244E-2_wp,-3.37789986050220E-1_wp, 1.13510089629769E-1_wp,&
      &-2.07382246739143E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_o2)

end subroutine test_g_onsite_gfn2_m04_atom

subroutine test_g_onsite_gfn2_m08_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.11048312695985E-1_wp,-5.02011645803230E-1_wp, 4.15238062649689E-1_wp, &
      &-3.25959600753673E-1_wp, 2.51473641195433E-2_wp, 2.93748490123740E-1_wp, &
      & 2.56736194030896E-2_wp, 2.38762690307426E-2_wp,-6.03118603733083E-1_wp, &
      & 3.91990240426822E-1_wp, 8.97114734113785E-1_wp, 1.93532936362436E-1_wp, &
      &-1.03136268223866E-1_wp,-1.04447608767710E-1_wp,-2.64818891498402E-2_wp, &
      &-3.90117787102468E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      & 9.06259904944829E-1_wp,-1.11730821567902E+0_wp, 2.78017329305492E-1_wp, &
      &-7.80028989546297E-1_wp, 1.11352815063389E+0_wp,-6.98290073981154E-1_wp, &
      & 2.03943236255318E-1_wp,-5.29902840441233E-1_wp, 4.38219939650397E-2_wp, &
      &-1.86746328945826E-2_wp, 4.65996457236599E-2_wp, 4.97590807484258E-1_wp, &
      &-2.50441962186972E-1_wp, 4.83295451755440E-2_wp,-2.26559244782012E-2_wp, &
      & 4.50331992744248E-2_wp,-2.11569328297532E-2_wp, 3.12470620007346E-1_wp, &
      &-9.15589243372491E-1_wp, 1.06394261835743E+0_wp,-6.71952361588756E-1_wp, &
      & 1.82322476598938E+0_wp,-9.26110009158329E-1_wp, 9.78357111140355E-1_wp, &
      &-7.84824170464332E-1_wp,-9.43549308434806E-2_wp,-8.78133979988158E-3_wp, &
      &-7.07783143624696E-2_wp,-3.36692984665466E-2_wp, 6.75375129657761E-1_wp, &
      &-7.01857024422455E-1_wp, 2.11598132242645E-1_wp,-6.01715925641418E-1_wp]

   call get_structure(mol, "MB16-43", "08")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_o2)

end subroutine test_g_onsite_gfn2_m08_shell

subroutine test_g_twobody_gxtb_lih_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  3.903043569209723E-1_wp, -3.903043569695019E-1_wp]
   real(wp), parameter :: qsh(*) = [ &
      &  1.883240895671248E-1_wp,  2.019802673538475E-1_wp, -3.903043569695019E-1_wp]

   call get_structure(mol, "MB16-43", "LiH")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_g_twobody_gxtb_lih_shell

subroutine test_g_twobody_gxtb_m03_atom(error)

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
   call test_numgrad(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_g_twobody_gxtb_m03_atom

subroutine test_g_twobody_gxtb_m03_shell(error)

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
   call test_numgrad(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_g_twobody_gxtb_m03_shell

subroutine test_g_twobody_gxtb_m04_atom(error)

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
   call test_numgrad(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_g_twobody_gxtb_m04_atom

subroutine test_g_twobody_gxtb_m04_shell(error)

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
   call test_numgrad(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_g_twobody_gxtb_m04_shell


subroutine test_s_onsite_gfn2_m09_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.13260038539900E-2_wp, 1.10070523471231E-2_wp,-9.97165474630829E-2_wp, &
      &-8.78527301724521E-2_wp, 2.89049242695863E-1_wp, 3.57284006856323E-2_wp, &
      &-1.73226219187217E-1_wp, 1.61174372420268E-1_wp,-8.89089419183055E-2_wp, &
      & 3.23950178196666E-2_wp, 1.88420688366637E-1_wp, 4.14882523279327E-2_wp, &
      &-2.23498403532295E-1_wp,-3.55334728213004E-1_wp,-7.15753987897201E-2_wp, &
      & 3.52175946466941E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      & 1.40887956776581E-3_wp,-1.27348820058716E-2_wp, 2.63961183739554E-2_wp, &
      &-1.53890176131402E-2_wp,-8.73648546608390E-2_wp,-1.23517435478204E-2_wp, &
      &-8.71021014527735E-2_wp,-7.50559382736492E-4_wp, 7.82044211296174E-1_wp, &
      &-4.92995083533018E-1_wp, 4.84143136555792E-2_wp,-1.26858387490357E-2_wp, &
      & 9.72488073646510E-1_wp,-1.14571448042039E+0_wp, 1.07574874045191E+0_wp, &
      &-9.14574293473561E-1_wp,-7.63358458276189E-2_wp,-1.25730981035572E-2_wp, &
      & 4.44349073468088E-2_wp,-1.20397879426510E-2_wp, 5.20245311277456E-1_wp, &
      & 1.92282483805197E-1_wp,-5.24107355799204E-1_wp, 5.39382871928999E-2_wp, &
      &-1.24499232808976E-2_wp, 7.97368410133983E-2_wp,-3.13209082796440E-1_wp, &
      & 9.97387287057362E-3_wp, 1.94446888375020E-1_wp,-5.49781435696375E-1_wp, &
      &-6.89789344411558E-2_wp,-2.59643153694089E-3_wp, 1.09519797601190E+0_wp, &
      &-7.43022154621128E-1_wp]

   call get_structure(mol, "MB16-43", "09")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_o2)

end subroutine test_s_onsite_gfn2_m09_shell

subroutine test_s_onsite_gfn2_pyrazine_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
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
   call test_numsigma(error, mol, qat, qsh, make_coulomb_o2)

end subroutine test_s_onsite_gfn2_pyrazine_atom

subroutine test_s_twobody_gxtb_lih_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  3.903043569209723E-1_wp, -3.903043569695019E-1_wp]
   real(wp), parameter :: qsh(*) = [ &
      &  1.883240895671248E-1_wp,  2.019802673538475E-1_wp, -3.903043569695019E-1_wp]

   call get_structure(mol, "MB16-43", "LiH")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_s_twobody_gxtb_lih_shell

subroutine test_s_twobody_gxtb_m05_atom(error)

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
   call test_numsigma(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_s_twobody_gxtb_m05_atom

subroutine test_s_twobody_gxtb_m05_shell(error)

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
   call test_numsigma(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_s_twobody_gxtb_m05_shell

subroutine test_s_twobody_gxtb_m06_atom(error)

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
   call test_numsigma(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_s_twobody_gxtb_m06_atom

subroutine test_s_twobody_gxtb_m06_shell(error)

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
   call test_numsigma(error, mol, qat, qsh, make_coulomb_t_gxtb, thr_in=thr1)

end subroutine test_s_twobody_gxtb_m06_shell


subroutine test_pg_ceh_lih_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_oceh, &
      & .false., thr_in=thr1)

end subroutine test_pg_ceh_lih_atom

subroutine test_pg_ceh_m15_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "15")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_oceh, &
      & .false., thr_in=thr1)

end subroutine test_pg_ceh_m15_atom

subroutine test_pg_ceh_m16_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "16")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_oceh, &
      & .true., thr_in=thr1)

end subroutine test_pg_ceh_m16_shell


subroutine test_ps_ceh_lih_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_oceh, &
      & .false., thr_in=thr1)

end subroutine test_ps_ceh_lih_atom

subroutine test_ps_ceh_co2_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "CO2")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_oceh, &
      & .false., thr_in=thr1)

end subroutine test_ps_ceh_co2_atom

subroutine test_ps_ceh_m05_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_oceh, &
      & .false., thr_in=thr1*10)

end subroutine test_ps_ceh_m05_atom

subroutine test_ps_ceh_m17_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "17")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_oceh, &
      & .true., thr_in=thr1*100)

end subroutine test_ps_ceh_m17_shell

end module test_coulomb_thirdorder
