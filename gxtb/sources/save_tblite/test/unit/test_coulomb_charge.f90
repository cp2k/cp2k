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

module test_coulomb_charge
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use mstore, only : get_structure
   use tblite_basis_type, only : basis_type, new_basis, cgto_container, cgto_type
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_ceh_ceh, only : get_effective_qat
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_charge, only : effective_coulomb, new_effective_coulomb, &
      & gamma_coulomb, new_gamma_coulomb
   use tblite_cutoff, only : get_lattice_points
   use tblite_scf_potential, only: new_potential, potential_type
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_utils_average, only : average_type, new_average, average_id
   use tblite_xtb_coulomb, only : tb_coulomb, new_coulomb
   implicit none
   private

   public :: collect_coulomb_charge

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
subroutine collect_coulomb_charge(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-atom-eff-gfn1-m01", test_e_effective_gfn1_m01_atom), &
      new_unittest("energy-shell-eff-gfn1-m07", test_e_effective_gfn1_m07_shell), &
      new_unittest("energy-atom-eff-gfn2-m01", test_e_effective_gfn2_m02_atom), &
      new_unittest("energy-atom-eff-gfn2-oxacb-pbc", test_e_effective_gfn2_oxacb_atom), &
      new_unittest("energy-atom-eff-gfn2-oxacb-sc", test_e_effective_gfn2_oxacb_sc_atom), &
      new_unittest("energy-shell-eff-gxtb-lih", test_e_effective_gxtb_lih_shell), &
      new_unittest("energy-atom-eff-gxtb-m01", test_e_effective_gxtb_m01_atom), &
      new_unittest("energy-shell-eff-gxtb-m01", test_e_effective_gxtb_m01_shell), &
      new_unittest("energy-atom-eff-gxtb-m02", test_e_effective_gxtb_m02_atom), &
      new_unittest("energy-shell-eff-gxtb-m02", test_e_effective_gxtb_m02_shell), &
      new_unittest("energy-atom-gam-gfn1-m10", test_e_gamma_gfn1_m10_atom), &
      new_unittest("energy-shell-gam-gfn2-m13", test_e_gamma_gfn2_m13_shell), &
      new_unittest("potential-atom-eff-gxtb-m01", test_p_effective_gxtb_m01_atom), &
      new_unittest("potential-shell-eff-gxtb-m01", test_p_effective_gxtb_m01_shell), &
      new_unittest("potential-atom-eff-gxtb-m02", test_p_effective_gxtb_m02_atom), &
      new_unittest("potential-shell-eff-gxtb-m02", test_p_effective_gxtb_m02_shell), &
      new_unittest("amat-deriv-shell-eff-gfn1-m07", test_amat_effective_gfn1_m07_shell), &
      new_unittest("amat-deriv-shell-gam-gfn2-m14", test_amat_gamma_gfn2_m14_shell), &
      new_unittest("amat-deriv-shell-eff-gxtb-h2", test_amat_effective_gxtb_h2_shell), &
      new_unittest("amat-deriv-shell-eff-gxtb-lih", test_amat_effective_gxtb_lih_shell), &
      new_unittest("amat-deriv-shell-eff-gxtb-s2", test_amat_effective_gxtb_s2_shell), &
      new_unittest("amat-deriv-shell-eff-gxtb-cecl3", test_amat_effective_gxtb_cecl3_shell), &
      new_unittest("gradient-atom-eff-gfn1-m03", test_g_effective_gfn1_m03_atom), &
      new_unittest("gradient-shell-eff-gfn1-m08", test_g_effective_gfn1_m08_shell), &
      new_unittest("gradient-atom-eff-gfn2-m04", test_g_effective_gfn2_m04_atom), &
      new_unittest("gradient-atom-eff-gfn2-co2-pbc", test_g_effective_gfn2_co2_atom), &
      new_unittest("gradient-shell-eff-gxtb-lih", test_g_effective_gxtb_lih_shell), &
      new_unittest("gradient-atom-eff-gxtb-m03", test_g_effective_gxtb_m03_atom), &
      new_unittest("gradient-shell-eff-gxtb-m03", test_g_effective_gxtb_m03_shell), &
      new_unittest("gradient-atom-eff-gxtb-m04", test_g_effective_gxtb_m04_atom), &
      new_unittest("gradient-shell-eff-gxtb-m04", test_g_effective_gxtb_m04_shell), &
      new_unittest("gradient-atom-gam-gfn1-m11", test_g_gamma_gfn1_m11_atom), &
      new_unittest("gradient-atom-gam-gfn1-urea-pbc", test_g_gamma_gfn1_urea_atom), &
      new_unittest("gradient-shell-gam-gfn2-m14", test_g_gamma_gfn2_m14_shell), &
      new_unittest("sigma-atom-eff-gfn1-m05", test_s_effective_gfn1_m05_atom), &
      new_unittest("sigma-shell-eff-gfn1-m09", test_s_effective_gfn1_m09_shell), &
      new_unittest("sigma-atom-eff-gfn2-m06", test_s_effective_gfn2_m06_atom), &
      new_unittest("sigma-atom-eff-gfn2-ammonia-pbc", test_s_effective_gfn2_ammonia_atom), &
      new_unittest("sigma-atom-eff-gxtb-m05", test_s_effective_gxtb_m05_atom), &
      !new_unittest("sigma-shell-eff-gxtb-m05", test_s_effective_gxtb_m05_shell), & ! not working because of diagonal contributions are not removed from the first-order part
      new_unittest("sigma-atom-eff-gxtb-m06", test_s_effective_gxtb_m06_atom), &
      !new_unittest("sigma-shell-eff-gxtb-m06", test_s_effective_gxtb_m06_shell), & ! not working because of diagonal contributions are not removed from the first-order part
      new_unittest("sigma-atom-gam-gfn1-m12", test_s_gamma_gfn1_m12_atom), &
      new_unittest("sigma-shell-gam-gfn2-m15", test_s_gamma_gfn2_m15_shell), &
      new_unittest("sigma-atom-gam-gfn2-pyrazine-pbc", test_s_gamma_gfn2_pyrazine_atom), &
      new_unittest("potential-gradient-atom-eff-effceh-lih", test_pg_ceh_lih_atom), &
      new_unittest("potential-gradient-atom-eff-effceh-m15", test_pg_ceh_m15_atom), &
      new_unittest("potential-gradient-shell-eff-effceh-m16", test_pg_ceh_m16_shell), &
      new_unittest("potential-sigma-atom-eff-effceh-lih", test_ps_ceh_lih_atom), &
      new_unittest("potential-sigma-atom-eff-effceh-co2", test_ps_ceh_co2_atom), &
      new_unittest("potential-sigma-atom-eff-effceh-m05", test_ps_ceh_m05_atom), &
      new_unittest("potential-sigma-shell-eff-effceh-m17", test_ps_ceh_m17_shell) &
      ]

end subroutine collect_coulomb_charge


!> Factory to setup the CEH basis set for testing of the potential (gradient)
subroutine make_basis(bas, mol, ng)
   !> Basis set information
   type(basis_type), intent(out) :: bas
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of Gaussians per Slater function
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


!> Factory to create electrostatic objects based on GFN1-xTB values
subroutine make_coulomb_e1(coulomb, mol, shell, error)

   !> New electrostatic object
   type(tb_coulomb), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: hubbard_parameter(20) = [&
      & 0.470099_wp, 1.441379_wp, 0.205342_wp, 0.274022_wp, 0.340530_wp, &
      & 0.479988_wp, 0.476106_wp, 0.583349_wp, 0.788194_wp, 0.612878_wp, &
      & 0.165908_wp, 0.354151_wp, 0.221658_wp, 0.438331_wp, 0.798319_wp, &
      & 0.643959_wp, 0.519712_wp, 0.529906_wp, 0.114358_wp, 0.134187_wp]
   integer, parameter :: shell_count(20) = [&
      & 2, 1, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 3, 3, 3, 3, 2, 3]
   real(wp), parameter :: shell_scale(3, 20) = reshape([&
      & 0.0_wp, 0.0000000_wp, 0.0000000_wp,  0.0_wp, 0.0000000_wp, 0.0000000_wp, &
      & 0.0_wp,-0.0772012_wp, 0.0000000_wp,  0.0_wp, 0.1113005_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0165643_wp, 0.0000000_wp,  0.0_wp,-0.0471181_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0315090_wp, 0.0000000_wp,  0.0_wp, 0.0374608_wp, 0.0000000_wp, &
      & 0.0_wp,-0.0827352_wp, 0.0000000_wp,  0.0_wp,-0.3892542_wp, 0.0000000_wp, &
      & 0.0_wp,-0.3004391_wp, 0.0000000_wp,  0.0_wp, 0.0674819_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0503564_wp, 0.0000000_wp,  0.0_wp,-0.5925834_wp, 0.0000000_wp, &
      & 0.0_wp,-0.2530875_wp, 0.0000000_wp,  0.0_wp,-0.1678147_wp, 0.0000000_wp, &
      & 0.0_wp,-0.4481841_wp, 0.0000000_wp,  0.0_wp,-0.1450000_wp, 0.0000000_wp, &
      & 0.0_wp,-0.5332978_wp, 0.0000000_wp,  0.0_wp, 1.1522018_wp, 0.0000000_wp],&
      & shape(shell_scale)) + 1.0_wp
   real(wp), parameter :: gexp = 2.0_wp
   real(wp), allocatable :: hubbard(:, :)
   integer :: isp, izp, ish
   type(effective_coulomb), allocatable :: tmp
   type(average_type), allocatable :: average

   ! Setup coulomb interaction collection 
   allocate(coulomb)
   call new_coulomb(coulomb, mol, error)
   if (allocated(error)) return

   ! Setup effective coulomb interaction
   allocate(tmp)
   if (shell) then
      allocate(hubbard(3, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            hubbard(ish, isp) = hubbard_parameter(izp) * shell_scale(ish, izp)
         end do
      end do
      allocate(average)
      call new_average(average, average_id%harmonic)
      call new_effective_coulomb(tmp, mol, gexp, hubbard, average, &
         & nshell=shell_count(mol%num))
   else
      hubbard = reshape(hubbard_parameter(mol%num), [1, mol%nid])
      allocate(average)
      call new_average(average, average_id%harmonic)
      call new_effective_coulomb(tmp, mol, gexp, hubbard, average)
   end if
   call move_alloc(tmp, coulomb%es2)

end subroutine make_coulomb_e1

!> Factory to create electrostatic objects based on GFN2-xTB values
subroutine make_coulomb_e2(coulomb, mol, shell, error)

   !> New electrostatic object
   type(tb_coulomb), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: hubbard_parameter(20) = [&
      & 0.405771_wp, 0.642029_wp, 0.245006_wp, 0.684789_wp, 0.513556_wp, &
      & 0.538015_wp, 0.461493_wp, 0.451896_wp, 0.531518_wp, 0.850000_wp, &
      & 0.271056_wp, 0.344822_wp, 0.364801_wp, 0.720000_wp, 0.297739_wp, &
      & 0.339971_wp, 0.248514_wp, 0.502376_wp, 0.247602_wp, 0.320378_wp]
   integer, parameter :: shell_count(20) = [&
      & 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]
   real(wp), parameter :: shell_scale(3, 20) = reshape([&
      & 0.0_wp, 0.0000000_wp, 0.0000000_wp, 0.0_wp, 0.0000000_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1972612_wp, 0.0000000_wp, 0.0_wp, 0.9658467_wp, 0.0000000_wp, &
      & 0.0_wp, 0.3994080_wp, 0.0000000_wp, 0.0_wp, 0.1056358_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1164892_wp, 0.0000000_wp, 0.0_wp, 0.1497020_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1677376_wp, 0.0000000_wp, 0.0_wp, 0.1190576_wp,-0.3200000_wp, &
      & 0.0_wp, 0.1018894_wp, 0.0000000_wp, 0.0_wp, 1.4000000_wp,-0.0500000_wp, &
      & 0.0_wp,-0.0603699_wp, 0.2000000_wp, 0.0_wp,-0.5580042_wp,-0.2300000_wp, &
      & 0.0_wp,-0.1558060_wp,-0.3500000_wp, 0.0_wp,-0.1085866_wp,-0.2500000_wp, &
      & 0.0_wp, 0.4989400_wp, 0.5000000_wp, 0.0_wp,-0.0461133_wp,-0.0100000_wp, &
      & 0.0_wp, 0.3483655_wp, 0.0000000_wp, 0.0_wp, 1.5000000_wp,-0.2500000_wp],&
      & shape(shell_scale)) + 1.0_wp
   real(wp), parameter :: gexp = 2.0_wp
   real(wp), allocatable :: hubbard(:, :)
   integer :: isp, izp, ish
   type(effective_coulomb), allocatable :: tmp
   type(average_type), allocatable :: average

   ! Setup coulomb interaction collection 
   allocate(coulomb)
   call new_coulomb(coulomb, mol, error)
   if (allocated(error)) return

   ! Setup effective coulomb interaction
   allocate(tmp)
   if (shell) then
      allocate(hubbard(3, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            hubbard(ish, isp) = hubbard_parameter(izp) * shell_scale(ish, izp)
         end do
      end do
      allocate(average)
      call new_average(average, average_id%arithmetic)
      call new_effective_coulomb(tmp, mol, gexp, hubbard, average, &
         & nshell=shell_count(mol%num))
   else
      hubbard = reshape(hubbard_parameter(mol%num), [1, mol%nid])
      allocate(average)
      call new_average(average, average_id%arithmetic)
      call new_effective_coulomb(tmp, mol, gexp, hubbard, average)
   end if
   call move_alloc(tmp, coulomb%es2)

end subroutine make_coulomb_e2

!> Factory to create electrostatic objects based on GFN2-xTB values
subroutine make_coulomb_e_gxtb(coulomb, mol, shell, error)

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
   real(wp), parameter :: p_hubbard_cn(60) = [&
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

   !> Parameter: Reference shell occupation numbers averaged over the fitset
   real(wp), parameter :: p_refocc(4, 60) = reshape([&
      &  1.000000000000000_wp,  0.000000000000000_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !1
      &  2.000000000000000_wp,  0.000000000000000_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !2
      &  0.389964587324730_wp,  0.610035412626740_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !3
      &  0.866957380679220_wp,  1.133042619270460_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !4
      &  1.034919550460500_wp,  1.965080449490200_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !5
      &  1.035393989459650_wp,  2.964606010490320_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !6
      &  1.374608989337600_wp,  3.625391010612030_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !7
      &  1.673030369494220_wp,  4.326969630455940_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !8
      &  1.855971592974940_wp,  5.144028406975800_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !9
      &  2.024702720311780_wp,  5.975297279650150_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !10
      &  0.449668065502160_wp,  0.550331934445460_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !11
      &  0.905701704906800_wp,  0.747464940854180_wp,  0.346833354189680_wp,  0.000000000000000_wp, & !12
      &  1.189369318979370_wp,  1.401517047075460_wp,  0.409113633893730_wp,  0.000000000000000_wp, & !13
      &  1.387904678714880_wp,  2.324096632620660_wp,  0.287998688612370_wp,  0.000000000000000_wp, & !14
      &  1.569243405127110_wp,  2.991989632466600_wp,  0.438766962357950_wp,  0.000000000000000_wp, & !15
      &  1.757590244716220_wp,  4.089812500791580_wp,  0.152597254449820_wp,  0.000000000000000_wp, & !16
      &  1.875585647425410_wp,  5.067168204623920_wp,  0.057246147899970_wp,  0.000000000000000_wp, & !17
      &  2.032894154961360_wp,  5.927640487540430_wp,  0.039465357452110_wp,  0.000000000000000_wp, & !18
      &  0.549184343841530_wp,  0.450815656109330_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !19
      &  0.734785469716480_wp,  0.463488821669950_wp,  0.801725708560800_wp,  0.000000000000000_wp, & !20
      &  0.398850830889560_wp,  0.522067821265710_wp,  2.079081347792660_wp,  0.000000000000000_wp, & !21
      &  0.466595373323450_wp,  0.411147465207960_wp,  3.122257161416360_wp,  0.000000000000000_wp, & !22
      &  0.547573971664490_wp,  0.417024291972550_wp,  4.035401736314930_wp,  0.000000000000000_wp, & !23
      &  0.527949050985710_wp,  0.416603885183040_wp,  5.055447063783960_wp,  0.000000000000000_wp, & !24
      &  0.468277447309330_wp,  0.534770584581240_wp,  5.996951968057700_wp,  0.000000000000000_wp, & !25
      &  0.697283458577130_wp,  0.575283674370380_wp,  6.727432867004380_wp,  0.000000000000000_wp, & !26
      &  0.647270484720000_wp,  0.569297730542570_wp,  7.783431784687350_wp,  0.000000000000000_wp, & !27
      &  0.642232166509890_wp,  0.422019553437360_wp,  8.935748280005919_wp,  0.000000000000000_wp, & !28
      &  0.631701741931860_wp,  0.291625758227430_wp, 10.076672499791499_wp,  0.000000000000000_wp, & !29
      &  1.198698133534100_wp,  0.801301866414550_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !30
      &  1.343277117190120_wp,  1.389496581419170_wp,  0.267226301340220_wp,  0.000000000000000_wp, & !31
      &  1.527472867049060_wp,  2.240738055951160_wp,  0.231789076949180_wp,  0.000000000000000_wp, & !32
      &  1.785648562165390_wp,  2.967997490018190_wp,  0.246353947766780_wp,  0.000000000000000_wp, & !33
      &  1.847591386890390_wp,  3.950561565539410_wp,  0.201847047521140_wp,  0.000000000000000_wp, & !34
      &  1.898875710042050_wp,  5.021246602344470_wp,  0.079877687562420_wp,  0.000000000000000_wp, & !35
      &  2.053317045258220_wp,  5.854386583477780_wp,  0.092296371211850_wp,  0.000000000000000_wp, & !36
      &  0.490543450661540_wp,  0.509456549286070_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !37
      &  0.707514414863100_wp,  0.248351778014970_wp,  1.044133807071210_wp,  0.000000000000000_wp, & !38
      &  0.320648435337010_wp,  0.276274228431790_wp,  2.403077336186480_wp,  0.000000000000000_wp, & !39
      &  0.464607903069180_wp,  0.385254233164940_wp,  3.150137863714800_wp,  0.000000000000000_wp, & !40
      &  0.470319885212620_wp,  0.348580765978570_wp,  4.181099348759150_wp,  0.000000000000000_wp, & !41
      &  0.526770293776420_wp,  0.377561955749200_wp,  5.095667750423620_wp,  0.000000000000000_wp, & !42
      &  0.495435979579140_wp,  0.307536283587170_wp,  6.197027736782700_wp,  0.000000000000000_wp, & !43
      &  0.579624809297460_wp,  0.480359853710220_wp,  6.940015336940880_wp,  0.000000000000000_wp, & !44
      &  0.475272441433770_wp,  0.419667445643640_wp,  8.105060112872950_wp,  0.000000000000000_wp, & !45
      &  0.433824472468830_wp,  0.251792608076740_wp,  9.314382919404551_wp,  0.000000000000000_wp, & !46
      &  0.581223558603990_wp,  0.239262896010660_wp, 10.179513545335659_wp,  0.000000000000000_wp, & !47
      &  1.269902889516260_wp,  0.730097110434900_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !48
      &  1.500177584702030_wp,  1.257589325737200_wp,  0.242233089511580_wp,  0.000000000000000_wp, & !49
      &  1.586658968457270_wp,  2.147898714645540_wp,  0.265442316846750_wp,  0.000000000000000_wp, & !50
      &  1.913169164328120_wp,  2.780920558155240_wp,  0.305910277465950_wp,  0.000000000000000_wp, & !51
      &  1.884454140098540_wp,  3.919041194364880_wp,  0.196504665485680_wp,  0.000000000000000_wp, & !52
      &  1.926023784512330_wp,  4.980077406759490_wp,  0.093898808679490_wp,  0.000000000000000_wp, & !53
      &  2.059350430129080_wp,  5.781517383702110_wp,  0.159132186124940_wp,  0.000000000000000_wp, & !54
      &  0.433650332022660_wp,  0.566349667931310_wp,  0.000000000000000_wp,  0.000000000000000_wp, & !55
      &  0.649245222582690_wp,  0.075763820649950_wp,  1.274990956716850_wp,  0.000000000000000_wp, & !56
      &  0.320648435337010_wp,  0.276274228431790_wp,  2.403077336186480_wp,  0.000000000000000_wp, & !57
      &  0.908673499222234_wp,  0.140323854143347_wp,  1.710527220137139_wp,  1.240475426497280_wp, & !58
      &  0.879052624944570_wp,  0.148085394305093_wp,  1.610300718980087_wp,  2.362561261770250_wp, & !59
      &  0.909880261220291_wp,  0.165784953202508_wp,  1.621691076486291_wp,  3.302643709090910_wp],& !60
      & shape(p_refocc))

   !> Parameter: Shell-distributed effective valence nuclear charge
   real(wp), parameter :: p_zeffsh(4, 60) = reshape([&
      &  1.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !1
      &  2.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !2
      &  1.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !3
      &  2.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !4
      &  2.0000000000_wp,  1.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !5
      &  2.0000000000_wp,  2.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !6
      &  2.0000000000_wp,  3.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !7
      &  2.0000000000_wp,  4.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !8
      &  2.0000000000_wp,  5.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !9
      &  2.0000000000_wp,  6.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !10
      &  1.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !11
      &  2.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !12
      &  2.0000000000_wp,  1.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !13
      &  2.0000000000_wp,  2.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !14
      &  2.0000000000_wp,  3.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !15
      &  2.0000000000_wp,  4.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !16
      &  2.0000000000_wp,  5.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !17
      &  2.0000000000_wp,  6.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !18
      &  1.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !19
      &  2.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !20
      &  2.0000000000_wp,  0.0000000000_wp,  1.0000000000_wp,  0.0000000000_wp, & !21
      &  2.0000000000_wp,  0.0000000000_wp,  2.0000000000_wp,  0.0000000000_wp, & !22
      &  2.0000000000_wp,  0.0000000000_wp,  3.0000000000_wp,  0.0000000000_wp, & !23
      &  2.0000000000_wp,  0.0000000000_wp,  4.0000000000_wp,  0.0000000000_wp, & !24
      &  2.0000000000_wp,  0.0000000000_wp,  5.0000000000_wp,  0.0000000000_wp, & !25
      &  2.0000000000_wp,  0.0000000000_wp,  6.0000000000_wp,  0.0000000000_wp, & !26
      &  2.0000000000_wp,  0.0000000000_wp,  7.0000000000_wp,  0.0000000000_wp, & !27
      &  2.0000000000_wp,  0.0000000000_wp,  8.0000000000_wp,  0.0000000000_wp, & !28
      &  2.0000000000_wp,  0.0000000000_wp,  9.0000000000_wp,  0.0000000000_wp, & !29
      &  2.0000000000_wp,  0.0000000000_wp, 10.0000000000_wp,  0.0000000000_wp, & !30
      &  2.0000000000_wp,  1.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !31
      &  2.0000000000_wp,  2.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !32
      &  2.0000000000_wp,  3.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !33
      &  2.0000000000_wp,  4.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !34
      &  2.0000000000_wp,  5.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !35
      &  2.0000000000_wp,  6.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !36
      &  1.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !37
      &  2.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !38
      &  2.0000000000_wp,  0.0000000000_wp,  1.0000000000_wp,  0.0000000000_wp, & !39
      &  2.0000000000_wp,  0.0000000000_wp,  2.0000000000_wp,  0.0000000000_wp, & !40
      &  2.0000000000_wp,  0.0000000000_wp,  3.0000000000_wp,  0.0000000000_wp, & !41
      &  2.0000000000_wp,  0.0000000000_wp,  4.0000000000_wp,  0.0000000000_wp, & !42
      &  2.0000000000_wp,  0.0000000000_wp,  5.0000000000_wp,  0.0000000000_wp, & !43
      &  2.0000000000_wp,  0.0000000000_wp,  6.0000000000_wp,  0.0000000000_wp, & !44
      &  2.0000000000_wp,  0.0000000000_wp,  7.0000000000_wp,  0.0000000000_wp, & !45
      &  2.0000000000_wp,  0.0000000000_wp,  8.0000000000_wp,  0.0000000000_wp, & !46
      &  2.0000000000_wp,  0.0000000000_wp,  9.0000000000_wp,  0.0000000000_wp, & !47
      &  2.0000000000_wp,  0.0000000000_wp, 10.0000000000_wp,  0.0000000000_wp, & !48
      &  2.0000000000_wp,  1.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !49
      &  2.0000000000_wp,  2.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !50
      &  2.0000000000_wp,  3.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !51
      &  2.0000000000_wp,  4.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !52
      &  2.0000000000_wp,  5.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !53
      &  2.0000000000_wp,  6.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !54
      &  1.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !55
      &  2.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp,  0.0000000000_wp, & !56
      &  2.0000000000_wp,  0.0000000000_wp,  1.0000000000_wp,  0.0000000000_wp, & !57
      &  2.0000000000_wp,  0.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & !58
      &  2.0000000000_wp,  0.0000000000_wp,  1.0000000000_wp,  2.0000000000_wp, & !59
      &  2.0000000000_wp,  0.0000000000_wp,  1.0000000000_wp,  3.0000000000_wp],& !60
      & shape(p_zeffsh))

   !> Coulomb-kernel exponent for second-order tight-binding
   real(wp), parameter :: p_gexp = 1.0000000000_wp
   
   !> Exponent in radius-dependent hubbard scaling for second-order tight-binding
   real(wp), parameter :: p_hubbard_exp = 0.2946211550_wp

   !> Coordination number exponent
   real(wp), parameter :: p_cn_exp = 2.0680000000_wp

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

   real(wp), allocatable :: hubbard(:, :), hardness_cn(:), refqsh(:), tb_cn_rcov(:)
   integer :: isp, izp, ish, ii, iat, ind
   integer, allocatable :: ish_at(:)
   type(effective_coulomb), allocatable :: tmp
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
      allocate(hubbard(4, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, p_nshell(izp)
            hubbard(ish, isp) = p_hubbard_parameter(izp) * p_shell_hubbard(ish, izp)
         end do
      end do

      allocate(ish_at(mol%nat))
      ind = 0
      do iat = 1, mol%nat
         isp = mol%id(iat)
         izp = mol%num(isp)        
         ish_at(iat) = ind
         ind = ind + p_nshell(izp)
      end do
   
      allocate(refqsh(ind))
      do iat = 1, mol%nat
         isp = mol%id(iat)
         izp = mol%num(isp)
         ii = ish_at(iat)
         do ish = 1, p_nshell(izp)
            refqsh(ii + ish) = p_refocc(ish, izp) - p_zeffsh(ish, izp)
         end do
      end do

      hardness_cn = p_hubbard_cn(mol%num)
      allocate(average)
      call new_average(average, average_id%harmonic)
      call new_effective_coulomb(tmp, mol, p_gexp, hubbard, average, &
         & hardness_cn, p_hubbard_exp, refqsh, nshell=p_nshell(mol%num))
   else
      allocate(hubbard(1, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         hubbard(1, isp) = p_hubbard_parameter(izp)
      end do
      hardness_cn = p_hubbard_cn(mol%num)
      allocate(average)
      call new_average(average, average_id%harmonic)
      call new_effective_coulomb(tmp, mol, p_gexp, hubbard, average, &
         & hardness_cn, p_hubbard_exp, refqsh)
   end if
   call move_alloc(tmp, coulomb%es2)

end subroutine make_coulomb_e_gxtb

!> Factory to create electrostatic objects based on CEH values
subroutine make_coulomb_eceh(coulomb, mol, shell, error)

   !> New electrostatic object
   type(tb_coulomb), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: hubbard_parameter(20) = [&
      & 0.47259288_wp, 0.92203391_wp, 0.17452888_wp, 0.25700733_wp, &
      & 0.33949086_wp, 0.42195412_wp, 0.50438193_wp, 0.58691863_wp, &
      & 0.66931351_wp, 0.75191607_wp, 0.17964105_wp, 0.22157276_wp, &
      & 0.26348578_wp, 0.30539645_wp, 0.34734014_wp, 0.38924725_wp, &
      & 0.43115670_wp, 0.47308269_wp, 0.17105469_wp, 0.20276244_wp]
   integer, parameter :: shell_count(20) = [&
      & 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]

   real(wp), parameter :: gexp = 1.0_wp
   real(wp), allocatable :: hubbard(:, :)
   integer :: isp, izp, ish
   type(effective_coulomb), allocatable :: tmp
   type(average_type), allocatable :: average

   ! Setup coulomb interaction collection 
   allocate(coulomb)
   call new_coulomb(coulomb, mol, error)
   if (allocated(error)) return

   ! Setup effective coulomb interaction
   allocate(tmp)
   if (shell) then
      ! If shell-resolved, we use the atomic parameter for each shell
      allocate(hubbard(3, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            hubbard(ish, isp) = hubbard_parameter(izp)
         end do
      end do
      allocate(average)
      call new_average(average, average_id%arithmetic)
      call new_effective_coulomb(tmp, mol, gexp, hubbard, average, &
         & nshell=shell_count(mol%num))
   else
      hubbard = reshape(hubbard_parameter(mol%num), [1, mol%nid])
      allocate(average)
      call new_average(average, average_id%arithmetic)
      call new_effective_coulomb(tmp, mol, gexp, hubbard, average)
   end if
   call move_alloc(tmp, coulomb%es2)

end subroutine make_coulomb_eceh

!> Factory to create electrostatic objects based on GFN1-xTB values
subroutine make_coulomb_g1(coulomb, mol, shell, error)

   !> New electrostatic object
   type(tb_coulomb), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: hubbard_parameter(20) = [&
      & 0.470099_wp, 1.441379_wp, 0.205342_wp, 0.274022_wp, 0.340530_wp, &
      & 0.479988_wp, 0.476106_wp, 0.583349_wp, 0.788194_wp, 0.612878_wp, &
      & 0.165908_wp, 0.354151_wp, 0.221658_wp, 0.438331_wp, 0.798319_wp, &
      & 0.643959_wp, 0.519712_wp, 0.529906_wp, 0.114358_wp, 0.134187_wp]
   integer, parameter :: shell_count(20) = [&
      & 2, 1, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 3, 3, 3, 3, 3, 3, 2, 3]
   real(wp), parameter :: shell_scale(3, 20) = reshape([&
      & 0.0_wp, 0.0000000_wp, 0.0000000_wp,  0.0_wp, 0.0000000_wp, 0.0000000_wp, &
      & 0.0_wp,-0.0772012_wp, 0.0000000_wp,  0.0_wp, 0.1113005_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0165643_wp, 0.0000000_wp,  0.0_wp,-0.0471181_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0315090_wp, 0.0000000_wp,  0.0_wp, 0.0374608_wp, 0.0000000_wp, &
      & 0.0_wp,-0.0827352_wp, 0.0000000_wp,  0.0_wp,-0.3892542_wp, 0.0000000_wp, &
      & 0.0_wp,-0.3004391_wp, 0.0000000_wp,  0.0_wp, 0.0674819_wp, 0.0000000_wp, &
      & 0.0_wp, 0.0503564_wp, 0.0000000_wp,  0.0_wp,-0.5925834_wp, 0.0000000_wp, &
      & 0.0_wp,-0.2530875_wp, 0.0000000_wp,  0.0_wp,-0.1678147_wp, 0.0000000_wp, &
      & 0.0_wp,-0.4481841_wp, 0.0000000_wp,  0.0_wp,-0.1450000_wp, 0.0000000_wp, &
      & 0.0_wp,-0.5332978_wp, 0.0000000_wp,  0.0_wp, 1.1522018_wp, 0.0000000_wp],&
      & shape(shell_scale)) + 1.0_wp
   real(wp), parameter :: gexp = 2.0_wp
   real(wp), allocatable :: hubbard(:, :)
   integer :: isp, izp, ish
   type(gamma_coulomb), allocatable :: tmp

   ! Setup coulomb interaction collection 
   allocate(coulomb)
   call new_coulomb(coulomb, mol, error)
   if (allocated(error)) return

   ! Setup gamma coulomb interaction
   allocate(tmp)
   if (shell) then
      allocate(hubbard(3, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            hubbard(ish, isp) = hubbard_parameter(izp) * shell_scale(ish, izp)
         end do
      end do
      call new_gamma_coulomb(tmp, mol, hubbard, nshell=shell_count(mol%num))
   else
      hubbard = reshape(hubbard_parameter(mol%num), [1, mol%nid])
      call new_gamma_coulomb(tmp, mol, hubbard)
   end if
   call move_alloc(tmp, coulomb%es2)

end subroutine make_coulomb_g1

!> Factory to create electrostatic objects based on GFN1-xTB values
subroutine make_coulomb_g2(coulomb, mol, shell, error)

   !> New electrostatic object
   type(tb_coulomb), allocatable, intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Return a shell resolved object
   logical, intent(in) :: shell

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: hubbard_parameter(20) = [&
      & 0.405771_wp, 0.642029_wp, 0.245006_wp, 0.684789_wp, 0.513556_wp, &
      & 0.538015_wp, 0.461493_wp, 0.451896_wp, 0.531518_wp, 0.850000_wp, &
      & 0.271056_wp, 0.344822_wp, 0.364801_wp, 0.720000_wp, 0.297739_wp, &
      & 0.339971_wp, 0.248514_wp, 0.502376_wp, 0.247602_wp, 0.320378_wp]
   integer, parameter :: shell_count(20) = [&
      & 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3]
   real(wp), parameter :: shell_scale(3, 20) = reshape([&
      & 0.0_wp, 0.0000000_wp, 0.0000000_wp, 0.0_wp, 0.0000000_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1972612_wp, 0.0000000_wp, 0.0_wp, 0.9658467_wp, 0.0000000_wp, &
      & 0.0_wp, 0.3994080_wp, 0.0000000_wp, 0.0_wp, 0.1056358_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1164892_wp, 0.0000000_wp, 0.0_wp, 0.1497020_wp, 0.0000000_wp, &
      & 0.0_wp, 0.1677376_wp, 0.0000000_wp, 0.0_wp, 0.1190576_wp,-0.3200000_wp, &
      & 0.0_wp, 0.1018894_wp, 0.0000000_wp, 0.0_wp, 1.4000000_wp,-0.0500000_wp, &
      & 0.0_wp,-0.0603699_wp, 0.2000000_wp, 0.0_wp,-0.5580042_wp,-0.2300000_wp, &
      & 0.0_wp,-0.1558060_wp,-0.3500000_wp, 0.0_wp,-0.1085866_wp,-0.2500000_wp, &
      & 0.0_wp, 0.4989400_wp, 0.5000000_wp, 0.0_wp,-0.0461133_wp,-0.0100000_wp, &
      & 0.0_wp, 0.3483655_wp, 0.0000000_wp, 0.0_wp, 1.5000000_wp,-0.2500000_wp],&
      & shape(shell_scale)) + 1.0_wp
   real(wp), parameter :: gexp = 2.0_wp
   real(wp), allocatable :: hubbard(:, :)
   integer :: isp, izp, ish
   type(gamma_coulomb), allocatable :: tmp

   ! Setup coulomb interaction collection 
   allocate(coulomb)
   call new_coulomb(coulomb, mol, error)
   if (allocated(error)) return

   ! Setup gamma coulomb interaction
   allocate(tmp)
   if (shell) then
      allocate(hubbard(3, mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         do ish = 1, shell_count(izp)
            hubbard(ish, isp) = hubbard_parameter(izp) * shell_scale(ish, izp)
         end do
      end do
      call new_gamma_coulomb(tmp, mol, hubbard, nshell=shell_count(mol%num))
   else
      hubbard = reshape(hubbard_parameter(mol%num), [1, mol%nid])
      call new_gamma_coulomb(tmp, mol, hubbard)
   end if
   call move_alloc(tmp, coulomb%es2)

end subroutine make_coulomb_g2


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

   !> Factory to create new electrostatic objects
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
   call coulomb%es2%get_energy(mol, cache, wfn, energy)

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
      call coulomb%es2%get_energy(mol, cache, wfn, er)
      ! Left hand side
      wfn%qat(iat, 1) = wfn%qat(iat, 1) - 2*step
      call coulomb%es2%get_energy(mol, cache, wfn, el)

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
         call coulomb%es2%get_energy(mol, cache, wfn, er)
         ! Left hand side
         wfn%qsh(ish, 1) = wfn%qsh(ish, 1) - 2*step
         call coulomb%es2%get_energy(mol, cache, wfn, el)

         wfn%qsh(ish, 1) = wfn%qsh(ish, 1) + step
         numvsh(ish) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end if

   ! Analytic potentials
   call pot%reset()
   call coulomb%es2%get_potential(mol, cache, wfn, pot)

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


subroutine test_amat_numgrad(error, mol, qat, qsh, make_coulomb, thr_in)

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

   integer :: iat, ic, ndim, is, ish, jsh
   type(tb_coulomb), allocatable :: coulomb
   type(container_cache) :: cache
   type(coulomb_cache), pointer :: ptr
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), datr(:, :)
   real(wp), allocatable :: dadcn(:, :), datrdcn(:), qvec(:)
   real(wp), allocatable :: amatr(:, :), amatl(:, :), numdadr(:, :, :)
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
   call coulomb%update(mol, cache)
   
   ! Numerical derivative of the Coulomb matrix
   allocate(dadr(3, mol%nat, ndim), numdadr(3, mol%nat, ndim), dadL(3, 3, ndim), &
      & datr(3, ndim), dadcn(ndim, ndim), datrdcn(ndim))
   dadr(:, :, :) = 0.0_wp 
   numdadr(:, :, :) = 0.0_wp
   allocate(amatr(ndim, ndim), amatl(ndim, ndim), source=0.0_wp)
   do iat = 1, mol%nat
      do ic = 1, 3
         ! Right hand side
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call coulomb%es2%update(mol, cache)
         call view(cache, ptr)
         amatr(:, :) = ptr%amat(:, :)

         ! Left hand side
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call coulomb%es2%update(mol, cache)
         call view(cache, ptr)
         amatl(:, :) = ptr%amat(:, :)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step

         ! Finite difference only for columns centered on the displaced atom, 
         ! as we calculate the once-contracted potential gradient.
         ! Possible since the second-order energy is bilinear in the charge
         is = coulomb%es2%offset(iat)
         do ish = 1, coulomb%es2%nshell(iat) 
            do jsh = 1, ndim
               numdadr(ic, iat, jsh) = numdadr(ic, iat, jsh) + 0.5_wp * qvec(is+ish) &
                  & * (amatr(jsh, is+ish) - amatl(jsh, is+ish))/step
            end do 
         end do 
      end do
   end do

   ! Analytic derivative for Coulomb matrix
   call coulomb%es2%update(mol, cache)
   call view(cache, ptr)
   call coulomb%es2%get_coulomb_derivs(mol, ptr, qat, qsh, &
      & dadr, dadL, datr, dadcn, datrdcn)

   if (any(abs(dadr - numdadr) > thr_)) then
      call test_failed(error, "Gradient of Coulomb A matrix does not match")
      write(*,*) "numerical gradient:"
      print'(3es21.14)', numdadr
      write(*,*) "analytical gradient:"
      print'(3es21.14)', dadr
      write(*,*) "difference:"
      print'(3es21.14)', dadr - numdadr
   end if

end subroutine test_amat_numgrad

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

   !> Factory to create new electrostatic objects
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
         call coulomb%es2%get_energy(mol, cache, wfn, er)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call coulomb%update(mol, cache)
         call coulomb%es2%get_energy(mol, cache, wfn, el)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call coulomb%update(mol, cache)
   call coulomb%es2%get_gradient(mol, cache, wfn, gradient, sigma)

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

   !> Factory to create new electrostatic objects
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
         call get_charges(wfn, mol, coulomb%es2%nshell, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%es2%get_potential(mol, cache, wfn, potr)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_charges(wfn, mol, coulomb%es2%nshell, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%es2%get_potential(mol, cache, wfn, potl)
         
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         if(shell) then
            numpotgrad(ic, iat, :) = 0.5_wp*(potr%vsh(:,1) - potl%vsh(:,1))/step
         else
            numpotgrad(ic, iat, :) = 0.5_wp*(potr%vat(:,1) - potl%vat(:,1))/step
         end if
      end do
   end do

   call get_charges(wfn, mol, coulomb%es2%nshell, error)
   if (allocated(error)) return
   call coulomb%update(mol, cache)
   call coulomb%es2%get_potential(mol, cache, wfn, potl)
   call coulomb%es2%get_potential_gradient(mol, cache, wfn, potl)
   
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
         call coulomb%es2%get_energy(mol, cache, wfn, er)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call coulomb%update(mol, cache)
         call coulomb%es2%get_energy(mol, cache, wfn, el)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (allocated(lattice)) mol%lattice = lattice
         numsigma(jc, ic) = 0.5_wp*(sum(er) - sum(el))/step
      end do
   end do

   call coulomb%update(mol, cache)
   call coulomb%es2%get_gradient(mol, cache, wfn, gradient, sigma)

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
         call get_charges(wfn, mol, coulomb%es2%nshell, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%es2%get_potential(mol, cache, wfn, potr)

         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call get_charges(wfn, mol, coulomb%es2%nshell, error)
         if (allocated(error)) return
         call coulomb%update(mol, cache)
         call coulomb%es2%get_potential(mol, cache, wfn, potl)
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

   call get_charges(wfn, mol, coulomb%es2%nshell, error)
   if (allocated(error)) return
   call coulomb%update(mol, cache)
   call coulomb%es2%get_potential(mol, cache, wfn, potl)
   call coulomb%es2%get_potential_gradient(mol, cache, wfn, potl)

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


subroutine test_e_effective_gfn1_m01_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.73347900345264E-1_wp, 1.07626888948184E-1_wp,-3.66999593831010E-1_wp,&
      & 4.92833325937897E-2_wp,-1.83332156197733E-1_wp, 2.33302086605469E-1_wp,&
      & 6.61837152062315E-2_wp,-5.43944165050002E-1_wp,-2.70264356583716E-1_wp,&
      & 2.66618968841682E-1_wp, 2.62725033202480E-1_wp,-7.15315510172571E-2_wp,&
      &-3.73300777019193E-1_wp, 3.84585237785621E-2_wp,-5.05851088366940E-1_wp,&
      & 5.17677238544189E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, qat, qsh, make_coulomb_e1, 0.10952019883948200_wp)

end subroutine test_e_effective_gfn1_m01_atom

subroutine test_e_effective_gfn1_m07_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.49712318034775E-1_wp, 2.12665850975202E-1_wp, 3.35977061494489E-1_wp, &
      & 3.16737890491354E-2_wp, 4.12434432866631E-2_wp,-3.21014009885608E-1_wp, &
      &-3.06535419089495E-1_wp,-5.36251066565321E-1_wp, 4.48758364798896E-1_wp, &
      & 6.00309584480896E-2_wp,-2.75470557482709E-1_wp, 3.60263594022551E-1_wp, &
      & 3.77425314022796E-1_wp,-6.30561365518420E-1_wp,-2.50675164323255E-1_wp, &
      & 6.02181524801775E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      & 8.85960229060055E-1_wp,-1.03567241653662E+0_wp, 2.34499192077770E-1_wp, &
      &-2.18333480864186E-2_wp, 1.09026104661485E+0_wp,-7.54283954798938E-1_wp, &
      & 4.12740327203921E-2_wp,-9.60021563849638E-3_wp, 5.17672944681095E-2_wp, &
      &-1.05238375989861E-2_wp, 5.94332546515988E-2_wp,-3.94897989828280E-1_wp, &
      & 1.44506731071946E-2_wp, 1.57870128213110E-1_wp,-4.64405557396352E-1_wp, &
      & 4.78122334280047E-1_wp,-1.01437364107707E+0_wp, 9.10337331767967E-1_wp, &
      &-4.61579000227231E-1_wp, 9.07619848805192E-2_wp,-3.07310018122722E-2_wp, &
      & 1.13955875471381E-1_wp,-3.99913576087036E-1_wp, 1.04872002787662E-2_wp, &
      & 4.12951024314537E-1_wp,-5.26874026571100E-2_wp, 4.04435991881125E-1_wp, &
      &-2.70107073714884E-2_wp, 3.13675308978710E-1_wp,-9.44236655190031E-1_wp, &
      & 1.75329569882602E-1_wp,-4.26004749886597E-1_wp, 1.24860566181157E+0_wp, &
      &-6.46424080267374E-1_wp]
   real(wp), allocatable :: qsh0(:)

   call get_structure(mol, "MB16-43", "07")
   call test_generic(error, mol, qat, qsh0, make_coulomb_e1, 0.13650692645610521_wp)
   call test_generic(error, mol, qat, qsh, make_coulomb_e1, 0.12017418620257683_wp)

end subroutine test_e_effective_gfn1_m07_shell

subroutine test_e_effective_gfn2_m02_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 7.38394711236234E-2_wp,-1.68354976558608E-1_wp,-3.47642833746823E-1_wp,&
      &-7.05489267186003E-1_wp, 7.73548301641266E-1_wp, 2.30207581365386E-1_wp,&
      & 1.02748501676354E-1_wp, 9.47818107467040E-2_wp, 2.44260351729187E-2_wp,&
      & 2.34984927037408E-1_wp,-3.17839896393030E-1_wp, 6.67112994818879E-1_wp,&
      &-4.78119977010488E-1_wp, 6.57536027459275E-2_wp, 1.08259054549882E-1_wp,&
      &-3.58215329983396E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, qat, qsh, make_coulomb_e2, 0.10635843572138280_wp)

end subroutine test_e_effective_gfn2_m02_atom

subroutine test_e_effective_gfn2_oxacb_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 3.41731844312030E-1_wp, 3.41716020106239E-1_wp, 3.41730526585671E-1_wp,&
      & 3.41714427217954E-1_wp, 3.80996046757999E-1_wp, 3.80989821246195E-1_wp,&
      & 3.81000747720282E-1_wp, 3.80990494183703E-1_wp,-3.70406587264474E-1_wp,&
      &-3.70407565207006E-1_wp,-3.70417590212352E-1_wp,-3.70399716470705E-1_wp,&
      &-3.52322260586075E-1_wp,-3.52304269439196E-1_wp,-3.52313440903261E-1_wp,&
      &-3.52298498047004E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "oxacb")
   call test_generic(error, mol, qat, qsh, make_coulomb_e2, 0.10130450083781417_wp)

end subroutine test_e_effective_gfn2_oxacb_atom

subroutine test_e_effective_gfn2_oxacb_sc_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
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
   call test_generic(error, mol, qat, qsh, make_coulomb_e2, &
      & 0.10130450083781417_wp*product(supercell), 1.0e-7_wp)

end subroutine test_e_effective_gfn2_oxacb_sc_atom

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

subroutine test_e_effective_gxtb_lih_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  3.903043569209723E-1_wp, -3.903043569695019E-1_wp]
   real(wp), parameter :: qsh(*) = [ &
      &  1.883240895671248E-1_wp,  2.019802673538475E-1_wp, -3.903043569695019E-1_wp]

   call get_structure(mol, "MB16-43", "LiH")
   call test_generic(error, mol, qat, qsh, make_coulomb_e_gxtb, &
      & 1.923159063644453E-2_wp, thr_in=thr)

end subroutine test_e_effective_gxtb_lih_shell

subroutine test_e_effective_gxtb_m01_atom(error)

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
   call test_generic(error, mol, qat, qsh, make_coulomb_e_gxtb, &
      & 0.31486225110686500_wp, thr_in=thr1)

end subroutine test_e_effective_gxtb_m01_atom

subroutine test_e_effective_gxtb_m01_shell(error)

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
   call test_generic(error, mol, qat, qsh, make_coulomb_e_gxtb, &
      & 0.212397063611163_wp, thr_in=thr1)

end subroutine test_e_effective_gxtb_m01_shell

subroutine test_e_effective_gxtb_m02_atom(error)

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
   call test_generic(error, mol, qat, qsh, make_coulomb_e_gxtb, &
      & 0.19952723917612791_wp, thr_in=thr1)

end subroutine test_e_effective_gxtb_m02_atom

subroutine test_e_effective_gxtb_m02_shell(error)

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
   call test_generic(error, mol, qat, qsh, make_coulomb_e_gxtb, &
      & 0.144986266743238_wp, thr_in=thr1)

end subroutine test_e_effective_gxtb_m02_shell

subroutine test_e_gamma_gfn1_m10_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.68014634286730E-2_wp, 4.42632083957210E-1_wp, 5.27534110370006E-3_wp, &
      &-3.49386920985384E-1_wp,-2.13178401684440E-1_wp, 4.07942205949245E-1_wp, &
      &-4.08514260972236E-1_wp, 1.34978625814380E-1_wp,-2.48254330281858E-1_wp, &
      &-3.65112235756872E-1_wp, 3.19858617682441E-1_wp, 2.96731604233838E-2_wp, &
      & 2.84061022228221E-1_wp,-3.25028474749853E-2_wp,-5.33914616408729E-3_wp, &
      & 1.46685494419848E-2_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "10")
   call test_generic(error, mol, qat, qsh, make_coulomb_g1, &  
      & 7.964758847421499E-2_wp, thr2)

end subroutine test_e_gamma_gfn1_m10_atom

subroutine test_e_gamma_gfn2_m13_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.62608233282119E-1_wp, 3.73633121487967E-1_wp, 1.51424532948944E-1_wp, &
      & 8.11274419840145E-2_wp, 4.55582555217907E-1_wp, 1.89469664895825E-1_wp, &
      &-3.59350817183894E-1_wp,-1.38911850317377E-1_wp,-1.83689392824396E-1_wp, &
      &-1.88906495161279E-1_wp, 5.33440028285669E-2_wp, 1.94112134916556E-1_wp, &
      & 2.02080948377078E-1_wp, 1.74595453525400E-1_wp,-4.46124496927388E-1_wp, &
      &-2.95778570663624E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      & 5.72134559421376E-2_wp,-2.68193499948548E-1_wp,-5.17064903935052E-2_wp, &
      & 3.73632853173886E-1_wp, 1.51424477665324E-1_wp, 8.11205008366953E-2_wp, &
      & 1.05453876337982E+0_wp,-4.64589774617786E-1_wp,-1.34371775944868E-1_wp, &
      & 3.21958772020979E-1_wp,-1.32435004307411E-1_wp, 2.88638899705825E-1_wp, &
      &-6.47972813769995E-1_wp, 6.82118177705109E-2_wp,-1.10631364729565E-1_wp, &
      &-9.64955180671905E-2_wp, 1.27185941911165E-1_wp,-3.10873201558534E-1_wp, &
      & 9.97036415531523E-2_wp,-2.88615729133477E-1_wp,-1.09656595674679E-1_wp, &
      & 1.63000176490660E-1_wp, 1.94112048312228E-1_wp,-4.48012133376332E-2_wp, &
      & 2.46906120671464E-1_wp, 1.74594629792853E-1_wp, 2.06932598673206E-1_wp, &
      &-6.52990922441455E-1_wp,-6.98603488893812E-3_wp,-2.88854759086314E-1_wp]
   real(wp), allocatable :: qsh0(:)

   call get_structure(mol, "MB16-43", "13")
   call test_generic(error, mol, qat, qsh0, make_coulomb_g2, &
      & 4.4263535114062461E-2_wp, thr2)
   if (allocated(error)) return
   call test_generic(error, mol, qat, qsh, make_coulomb_g2, &
      & 6.1195162961497636E-2_wp, thr2)

end subroutine test_e_gamma_gfn2_m13_shell


subroutine test_p_effective_gxtb_m01_atom(error)

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
   call test_numpot(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_p_effective_gxtb_m01_atom

subroutine test_p_effective_gxtb_m01_shell(error)

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
   call test_numpot(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_p_effective_gxtb_m01_shell

subroutine test_p_effective_gxtb_m02_atom(error)

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
   call test_numpot(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_p_effective_gxtb_m02_atom

subroutine test_p_effective_gxtb_m02_shell(error)

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
   call test_numpot(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_p_effective_gxtb_m02_shell


subroutine test_amat_effective_gfn1_m07_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-1.49712318034775E-1_wp, 2.12665850975202E-1_wp, 3.35977061494489E-1_wp, &
      & 3.16737890491354E-2_wp, 4.12434432866631E-2_wp,-3.21014009885608E-1_wp, &
      &-3.06535419089495E-1_wp,-5.36251066565321E-1_wp, 4.48758364798896E-1_wp, &
      & 6.00309584480896E-2_wp,-2.75470557482709E-1_wp, 3.60263594022551E-1_wp, &
      & 3.77425314022796E-1_wp,-6.30561365518420E-1_wp,-2.50675164323255E-1_wp, &
      & 6.02181524801775E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      & 8.85960229060055E-1_wp,-1.03567241653662E+0_wp, 2.34499192077770E-1_wp, &
      &-2.18333480864186E-2_wp, 1.09026104661485E+0_wp,-7.54283954798938E-1_wp, &
      & 4.12740327203921E-2_wp,-9.60021563849638E-3_wp, 5.17672944681095E-2_wp, &
      &-1.05238375989861E-2_wp, 5.94332546515988E-2_wp,-3.94897989828280E-1_wp, &
      & 1.44506731071946E-2_wp, 1.57870128213110E-1_wp,-4.64405557396352E-1_wp, &
      & 4.78122334280047E-1_wp,-1.01437364107707E+0_wp, 9.10337331767967E-1_wp, &
      &-4.61579000227231E-1_wp, 9.07619848805192E-2_wp,-3.07310018122722E-2_wp, &
      & 1.13955875471381E-1_wp,-3.99913576087036E-1_wp, 1.04872002787662E-2_wp, &
      & 4.12951024314537E-1_wp,-5.26874026571100E-2_wp, 4.04435991881125E-1_wp, &
      &-2.70107073714884E-2_wp, 3.13675308978710E-1_wp,-9.44236655190031E-1_wp, &
      & 1.75329569882602E-1_wp,-4.26004749886597E-1_wp, 1.24860566181157E+0_wp, &
      &-6.46424080267374E-1_wp]
   real(wp), allocatable :: qsh0(:)

   call get_structure(mol, "MB16-43", "07")
   call test_amat_numgrad(error, mol, qat, qsh, make_coulomb_e1, thr_in=thr1)

end subroutine test_amat_effective_gfn1_m07_shell

subroutine test_amat_gamma_gfn2_m14_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-5.35371225694038E-1_wp, 1.27905155882876E-1_wp, 2.06910619292535E-1_wp, &
      &-1.93061647443670E-1_wp, 5.46833043573218E-1_wp, 2.98577669101319E-1_wp, &
      &-3.62405585534705E-1_wp, 2.07231134137244E-1_wp, 2.85826164709174E-1_wp, &
      &-1.76518940177473E-1_wp, 9.44972704818130E-2_wp,-1.17451405142691E-1_wp, &
      &-1.41198286268662E-1_wp, 1.05227974737201E-2_wp,-1.31666840327078E-1_wp, &
      &-1.20629924063582E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      & 2.92177048596496E-1_wp,-8.27551283559270E-1_wp, 3.81811514612779E-1_wp, &
      & 4.17784916263666E-1_wp,-6.71683789364610E-1_wp,-5.11576334445887E-2_wp, &
      & 2.68488368058409E-1_wp,-1.04391705441501E-2_wp,-1.93062932974169E-1_wp, &
      & 7.61952673849415E-1_wp,-2.15114295745990E-1_wp, 2.98579359296096E-1_wp, &
      & 2.43594536377056E-2_wp,-3.79828052486015E-1_wp,-6.93861559657517E-3_wp, &
      & 9.40499498626918E-1_wp, 3.55733525506819E-1_wp,-1.08899859046528E+0_wp, &
      & 2.85827068858603E-1_wp,-1.76521898491791E-1_wp, 1.16110895118352E+0_wp, &
      &-1.06660369382576E+0_wp,-1.17453146019547E-1_wp,-1.41199843368407E-1_wp, &
      & 1.11757931464608E+0_wp,-1.10704872438469E+0_wp,-1.31668946428877E-1_wp, &
      &-1.20631076436794E-1_wp]

   call get_structure(mol, "MB16-43", "14")
   call test_amat_numgrad(error, mol, qat, qsh, make_coulomb_g2, thr_in=thr2*10)

end subroutine test_amat_gamma_gfn2_m14_shell

subroutine test_amat_effective_gxtb_h2_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      & 1.0_wp, 1.0_wp]
   real(wp), parameter :: qsh(*) = [ &
      & 1.0_wp, 1.0_wp]
   
   call get_structure(mol, "MB16-43", "H2")
   call test_amat_numgrad(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_amat_effective_gxtb_h2_shell

subroutine test_amat_effective_gxtb_lih_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  3.903043569209723E-1_wp, -3.903043569695019E-1_wp]
   real(wp), parameter :: qsh(*) = [ &
      &  1.883240895671248E-1_wp,  2.019802673538475E-1_wp, -3.903043569695019E-1_wp]

   call get_structure(mol, "MB16-43", "LiH")
   call test_amat_numgrad(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_amat_effective_gxtb_lih_shell

subroutine test_amat_effective_gxtb_s2_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      & -4.234713968376269E-11_wp, -4.241411388772320E-11_wp]
   real(wp), parameter :: qsh(*) = [ &
      & -1.658294650172361E-1_wp,  8.500381954781311E-2_wp,  8.082564542707581E-2_wp, &
      & -1.658294650172407E-1_wp,  8.500381954775360E-2_wp,  8.082564542707300E-2_wp]

   call get_structure(mol, "MB16-43", "S2")
   call test_amat_numgrad(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_amat_effective_gxtb_s2_shell

subroutine test_amat_effective_gxtb_cecl3_shell(error)

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
   call test_amat_numgrad(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_amat_effective_gxtb_cecl3_shell


subroutine test_g_effective_gfn1_m03_atom(error)

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
   call test_numgrad(error, mol, qat, qsh, make_coulomb_e1, thr_in=thr1)

end subroutine test_g_effective_gfn1_m03_atom

subroutine test_g_effective_gfn1_m08_shell(error)

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
   call test_numgrad(error, mol, qat, qsh, make_coulomb_e1, thr_in=thr1)

end subroutine test_g_effective_gfn1_m08_shell

subroutine test_g_effective_gfn2_m04_atom(error)

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
   call test_numgrad(error, mol, qat, qsh, make_coulomb_e2, thr_in=thr1)

end subroutine test_g_effective_gfn2_m04_atom

subroutine test_g_effective_gfn2_co2_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 4.56275672862067E-1_wp, 4.56284770386671E-1_wp, 4.56284770386671E-1_wp,&
      & 4.56284770386671E-1_wp,-2.28127680925611E-1_wp,-2.28138283131909E-1_wp,&
      &-2.28145770512561E-1_wp,-2.28145770512561E-1_wp,-2.28150142163058E-1_wp,&
      &-2.28145770512561E-1_wp,-2.28138283131909E-1_wp,-2.28138283131909E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "CO2")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_e2, thr_in=thr1)

end subroutine test_g_effective_gfn2_co2_atom

subroutine test_g_effective_gxtb_lih_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [ &
      &  3.903043569209723E-1_wp, -3.903043569695019E-1_wp]
   real(wp), parameter :: qsh(*) = [ &
      &  1.883240895671248E-1_wp,  2.019802673538475E-1_wp, -3.903043569695019E-1_wp]

   call get_structure(mol, "MB16-43", "LiH")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_g_effective_gxtb_lih_shell

subroutine test_g_effective_gxtb_m03_atom(error)

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
   call test_numgrad(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_g_effective_gxtb_m03_atom

subroutine test_g_effective_gxtb_m03_shell(error)

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
   call test_numgrad(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_g_effective_gxtb_m03_shell

subroutine test_g_effective_gxtb_m04_atom(error)

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
   call test_numgrad(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_g_effective_gxtb_m04_atom

subroutine test_g_effective_gxtb_m04_shell(error)

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
   call test_numgrad(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_g_effective_gxtb_m04_shell

subroutine test_g_gamma_gfn1_m11_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.90839692052395E-1_wp,-2.66254004019676E-1_wp,-1.14585349859674E-1_wp, &
      &-1.94168606876184E-1_wp, 1.69243400195097E-1_wp, 5.94099012700995E-2_wp, &
      & 3.58048105537982E-1_wp, 3.65662316444953E-2_wp, 3.38204991437465E-1_wp, &
      &-4.07570397699211E-1_wp, 5.27525279437458E-1_wp,-2.17656283937311E-1_wp, &
      &-3.01540791618602E-1_wp, 2.89587569744254E-1_wp,-1.87215698022911E-2_wp, &
      & 3.27512165984881E-2_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "11")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_g1, thr_in=thr1*10)

end subroutine test_g_gamma_gfn1_m11_atom

subroutine test_g_gamma_gfn1_urea_atom(error)

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
   call test_numgrad(error, mol, qat, qsh, make_coulomb_g1, thr_in=thr2*10)

end subroutine test_g_gamma_gfn1_urea_atom

subroutine test_g_gamma_gfn2_m14_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-5.35371225694038E-1_wp, 1.27905155882876E-1_wp, 2.06910619292535E-1_wp, &
      &-1.93061647443670E-1_wp, 5.46833043573218E-1_wp, 2.98577669101319E-1_wp, &
      &-3.62405585534705E-1_wp, 2.07231134137244E-1_wp, 2.85826164709174E-1_wp, &
      &-1.76518940177473E-1_wp, 9.44972704818130E-2_wp,-1.17451405142691E-1_wp, &
      &-1.41198286268662E-1_wp, 1.05227974737201E-2_wp,-1.31666840327078E-1_wp, &
      &-1.20629924063582E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      & 2.92177048596496E-1_wp,-8.27551283559270E-1_wp, 3.81811514612779E-1_wp, &
      & 4.17784916263666E-1_wp,-6.71683789364610E-1_wp,-5.11576334445887E-2_wp, &
      & 2.68488368058409E-1_wp,-1.04391705441501E-2_wp,-1.93062932974169E-1_wp, &
      & 7.61952673849415E-1_wp,-2.15114295745990E-1_wp, 2.98579359296096E-1_wp, &
      & 2.43594536377056E-2_wp,-3.79828052486015E-1_wp,-6.93861559657517E-3_wp, &
      & 9.40499498626918E-1_wp, 3.55733525506819E-1_wp,-1.08899859046528E+0_wp, &
      & 2.85827068858603E-1_wp,-1.76521898491791E-1_wp, 1.16110895118352E+0_wp, &
      &-1.06660369382576E+0_wp,-1.17453146019547E-1_wp,-1.41199843368407E-1_wp, &
      & 1.11757931464608E+0_wp,-1.10704872438469E+0_wp,-1.31668946428877E-1_wp, &
      &-1.20631076436794E-1_wp]

   call get_structure(mol, "MB16-43", "14")
   call test_numgrad(error, mol, qat, qsh, make_coulomb_g2, thr_in=thr2)

end subroutine test_g_gamma_gfn2_m14_shell


subroutine test_s_effective_gfn1_m05_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.01138111283277E-1_wp, 1.30358706339300E-1_wp, 9.38825924720944E-2_wp,&
      & 8.92795900801844E-2_wp, 5.13625440660610E-2_wp,-2.65500121876709E-2_wp,&
      & 9.26496972837658E-2_wp,-9.61095258223972E-2_wp,-4.92845009674246E-1_wp,&
      & 2.66730531684206E-1_wp, 3.37256104303071E-2_wp, 1.63170419985976E-1_wp,&
      & 6.91343155032824E-2_wp, 1.04287482572171E-1_wp, 6.09307909835941E-2_wp,&
      &-3.38869622433350E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "05")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_e1, thr_in=thr1)

end subroutine test_s_effective_gfn1_m05_atom

subroutine test_s_effective_gfn1_m09_shell(error)

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
   call test_numsigma(error, mol, qat, qsh, make_coulomb_e1, thr_in=thr1)

end subroutine test_s_effective_gfn1_m09_shell

subroutine test_s_effective_gfn2_m06_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-2.13983049532933E-1_wp,-5.10521279217923E-1_wp, 7.70190120699491E-2_wp,&
      &-3.68835155548212E-1_wp,-4.08747874260092E-1_wp,-4.09471309598929E-2_wp,&
      & 2.94164204769172E-1_wp, 9.76819709672870E-2_wp,-7.84337476935767E-3_wp,&
      & 7.07702520795024E-1_wp, 2.38774840136381E-1_wp, 1.08934666297455E-1_wp,&
      & 1.10156911889136E-1_wp, 9.25098455002779E-2_wp,-1.96776817442259E-1_wp,&
      & 2.07107093059868E-2_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "06")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_e2, thr_in=thr1)

end subroutine test_s_effective_gfn2_m06_atom

subroutine test_s_effective_gfn2_ammonia_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.95376975876519E-1_wp, 2.95376975876519E-1_wp, 2.95376975876519E-1_wp,&
      & 2.95329109335847E-1_wp, 2.95332441468412E-1_wp, 2.95347202855778E-1_wp,&
      & 2.95347202855779E-1_wp, 2.95329109335848E-1_wp, 2.95332441468411E-1_wp,&
      & 2.95347202855777E-1_wp, 2.95329109335847E-1_wp, 2.95332441468412E-1_wp,&
      &-8.86118742099358E-1_wp,-8.86012815503436E-1_wp,-8.86012815503437E-1_wp,&
      &-8.86012815503434E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "ammonia")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_e2, thr_in=thr2)

end subroutine test_s_effective_gfn2_ammonia_atom

subroutine test_s_effective_gxtb_m05_atom(error)

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
   call test_numsigma(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_s_effective_gxtb_m05_atom

subroutine test_s_effective_gxtb_m05_shell(error)

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
   call test_numsigma(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_s_effective_gxtb_m05_shell

subroutine test_s_effective_gxtb_m06_atom(error)

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
   call test_numsigma(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_s_effective_gxtb_m06_atom

subroutine test_s_effective_gxtb_m06_shell(error)

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
   call test_numsigma(error, mol, qat, qsh, make_coulomb_e_gxtb, thr_in=thr1)

end subroutine test_s_effective_gxtb_m06_shell

subroutine test_s_gamma_gfn1_m12_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 2.57523069224332E-1_wp,-7.13834599743419E-2_wp,-9.43945788149514E-2_wp, &
      & 1.02376579062895E-2_wp, 1.97960723912756E-1_wp, 3.16253846282241E-1_wp, &
      &-3.91548233613895E-1_wp,-1.68829398890385E-1_wp,-3.99798824173873E-1_wp, &
      & 4.22333212859043E-1_wp, 3.16980455307683E-1_wp,-1.09800808615267E-1_wp, &
      &-1.07582049146789E-1_wp,-4.48115660454027E-1_wp, 7.19397904194672E-2_wp, &
      & 1.98224257774052E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "MB16-43", "12")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_g1, thr_in=thr2)

end subroutine test_s_gamma_gfn1_m12_atom

subroutine test_s_gamma_gfn2_m15_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      &-5.41402496268596E-2_wp,-1.33777153976276E-1_wp, 4.14313829600631E-1_wp, &
      &-1.16641170075389E-1_wp, 4.56021377424607E-1_wp,-5.20378766868989E-1_wp, &
      &-1.63965423099635E-1_wp,-7.65345311273482E-2_wp, 2.08304494730413E-1_wp, &
      &-1.71827679329874E-1_wp,-3.30458156481514E-1_wp, 5.58638267294323E-1_wp, &
      &-3.10094747569162E-1_wp, 1.56794592474036E-1_wp, 2.13459748796815E-1_wp, &
      &-1.29714432165776E-1_wp]
   real(wp), parameter :: qsh(*) = [&
      &-5.41503330947710E-2_wp,-1.33779040525986E-1_wp, 6.99211871952128E-2_wp, &
      & 5.19410210372243E-1_wp,-1.75021497979320E-1_wp,-1.16646288193380E-1_wp, &
      & 5.88579806566877E-1_wp,-1.32544517759254E-1_wp, 7.82133136453700E-4_wp, &
      &-5.23133226533954E-1_wp, 1.95314227063571E-3_wp,-1.63971802434084E-1_wp, &
      &-7.65416499768423E-2_wp, 6.26092659280207E-1_wp, 3.00350609071998E-1_wp, &
      &-7.18137148569625E-1_wp,-1.71830433005059E-1_wp,-1.31761941444373E-1_wp, &
      &-1.98713319565700E-1_wp, 8.86264348375974E-2_wp, 7.46929250011706E-1_wp, &
      &-2.76884353276538E-1_wp, 4.27025703238206E-2_wp,-3.52590124894769E-1_wp, &
      &-2.13102162478342E-4_wp, 1.13142747674328E+0_wp,-9.74609683930292E-1_wp, &
      & 1.02802427493832E+0_wp,-8.14556120567196E-1_wp,-1.29715170834719E-1_wp]

   call get_structure(mol, "MB16-43", "15")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_g2, thr_in=thr2)

end subroutine test_s_gamma_gfn2_m15_shell

subroutine test_s_gamma_gfn2_pyrazine_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: qat(*) = [&
      & 1.23705026512437E-1_wp, 1.22537765989959E-1_wp, 1.23932626831231E-1_wp, &
      & 1.22418959326822E-1_wp, 1.23788033684569E-1_wp, 1.23643389058068E-1_wp, &
      & 1.22389811551880E-1_wp, 1.22402837718155E-1_wp, 3.17174594622166E-2_wp, &
      & 3.15585817789390E-2_wp, 3.14290005061543E-2_wp, 3.10358506314526E-2_wp, &
      & 3.11084225356749E-2_wp, 3.11187325528499E-2_wp, 3.10707965296438E-2_wp, &
      & 3.13508824430484E-2_wp,-3.09107070033859E-1_wp,-3.09144123845085E-1_wp, &
      &-3.08473365847409E-1_wp,-3.08483617386745E-1_wp]
   real(wp), allocatable :: qsh(:)

   call get_structure(mol, "X23", "pyrazine")
   call test_numsigma(error, mol, qat, qsh, make_coulomb_g2, thr_in=thr1)

end subroutine test_s_gamma_gfn2_pyrazine_atom


subroutine test_pg_ceh_lih_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_eceh, &
      & .false., thr_in=thr1)

end subroutine test_pg_ceh_lih_atom

subroutine test_pg_ceh_m15_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "15")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_eceh, &
      & .false., thr_in=thr1)

end subroutine test_pg_ceh_m15_atom

subroutine test_pg_ceh_m16_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "16")
   call test_numpotgrad(error, mol, get_charges_effceh, make_coulomb_eceh, &
      & .true., thr_in=thr1)

end subroutine test_pg_ceh_m16_shell


subroutine test_ps_ceh_lih_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "LiH")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_eceh, &
      & .false., thr_in=thr2)

end subroutine test_ps_ceh_lih_atom

subroutine test_ps_ceh_co2_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "CO2")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_eceh, &
      & .false., thr_in=thr2)

end subroutine test_ps_ceh_co2_atom

subroutine test_ps_ceh_m05_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_eceh, &
      & .false., thr_in=thr2)

end subroutine test_ps_ceh_m05_atom

subroutine test_ps_ceh_m17_shell(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "17")
   call test_numpotsigma(error, mol, get_charges_effceh, make_coulomb_eceh, &
      & .true., thr_in=thr2)

end subroutine test_ps_ceh_m17_shell

end module test_coulomb_charge
