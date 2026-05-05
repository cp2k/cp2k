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

module test_xtb_param
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_context_type, only : context_type
   use tblite_param, only : param_record, acp_record, charge_record, dispersion_record, &
      & element_record, exchange_record, firstorder_record, fourthorder_record, &
      & halogen_record, hamiltonian_record, increment_record, multipole_record, &
      & repulsion_record, spin_record, thirdorder_record, param_mask, count
   use tblite_param_molecular_moments, only:  molecular_multipole_record
   use tblite_toml, only : toml_table
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator, new_xtb_calculator
   use tblite_xtb_gfn2, only : export_gfn2_param, new_gfn2_calculator
   use tblite_xtb_gfn1, only : export_gfn1_param, new_gfn1_calculator
   use tblite_xtb_gxtb, only : export_gxtb_param, new_gxtb_calculator
   use tblite_xtb_ipea1, only : export_ipea1_param, new_ipea1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_xtb_param

   real(wp), parameter :: acc = 1.0_wp
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 3.166808578545117e-06_wp

contains


!> Collect all exported unit tests
subroutine collect_xtb_param(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("param-empty", test_param_empty, should_fail=.true.), &
      new_unittest("param-minimal", test_param_minimal), &
      new_unittest("param-invalid", test_param_invalid, should_fail=.true.), &
      new_unittest("element-empty", test_element_empty, should_fail=.true.), &
      new_unittest("acp-empty", test_acp_empty), &
      new_unittest("charge-empty", test_charge_empty, should_fail=.true.), &
      new_unittest("dispersion-empty", test_dispersion_empty, should_fail=.true.), &
      new_unittest("exchange-empty", test_exchange_empty, should_fail=.true.), &
      new_unittest("firstorder-empty", test_firstorder_empty, should_fail=.true.), &
      new_unittest("fourthorder-empty", test_fourthorder_empty, should_fail=.true.), &
      new_unittest("halogen-empty", test_halogen_empty, should_fail=.true.), &
      new_unittest("hamiltonian-empty", test_hamiltonian_empty, should_fail=.true.), &
      new_unittest("increment-empty", test_increment_empty), &
      new_unittest("multipole-empty", test_multipole_empty, should_fail=.true.), &
      new_unittest("repulsion-empty", test_repulsion_empty, should_fail=.true.), &
      new_unittest("spin-empty", test_spin_empty, should_fail=.true.), &
      new_unittest("thirdorder-empty", test_thirdorder_empty, should_fail=.true.), &
      new_unittest("mol-multipole-empty", test_mol_multipole_empty), &
      new_unittest("mask-gfn2", test_mask_gfn2), &
      new_unittest("mask-gfn1", test_mask_gfn1), &
      new_unittest("mask-gxtb", test_mask_gxtb), &
      new_unittest("gfn2-xtb-m02", test_gfn2_mb02), &
      new_unittest("gfn2-xtb-m05", test_gfn2_mb05), &
      new_unittest("gfn2-xtb-file-m07", test_gfn2_round_trip), &
      new_unittest("gfn1-xtb-m01", test_gfn1_mb01), &
      new_unittest("gfn1-xtb-m04", test_gfn1_mb04), &
      new_unittest("gfn1-xtb-file-m08", test_gfn1_round_trip), &
      new_unittest("ipea1-xtb-m03", test_ipea1_mb03), &
      new_unittest("ipea1-xtb-m06", test_ipea1_mb06), &
      new_unittest("ipea1-xtb-file-m09", test_ipea1_round_trip), &
      new_unittest("gxtb-xtb-water", test_gxtb_water), &
      new_unittest("gxtb-xtb-m07", test_gxtb_mb07), &
      new_unittest("gxtb-xtb-m08", test_gxtb_mb08), &
      new_unittest("gxtb-xtb-m10", test_gxtb_mb10), &
      new_unittest("gxtb-xtb-file-m10", test_gxtb_round_trip) &
      ]

end subroutine collect_xtb_param


subroutine test_param_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(param_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_param_empty


subroutine test_param_minimal(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: io
   type(param_record) :: param

   open(newunit=io, status="scratch")
   write(io, '(a)') &
      "[hamiltonian.xtb]", &
      "wexp = 0.5000000000000000", &
      "kpol = 2.0000000000000000", &
      "enscale = 0.0200000000000000", &
      "cn = ""dexp""", &
      "default_etemp = 300.0000000000000000", &
      "default_guess = ""sad""", &
      "charge_model = ""none""", &
      "[hamiltonian.xtb.shell]", &
      "ss = 1.8500000000000001", &
      "pp = 2.2300000000000000", &
      "dd = 2.2300000000000000", &
      "sd = 2.0000000000000000", &
      "pd = 2.0000000000000000", &
      "[element]", &
      "[element.H]", &
      "zeff = 1.1053880000000000", &
      "arep = 2.2137169999999999", &
      "gam = 0.4057710000000000", &
      "gam3 = 0.0800000000000000", &
      "dkernel = 0.0556388900000000", &
      "qkernel = 0.0002743100000000", &
      "mprad = 1.3999999999999999", &
      "mpvcn = 1.0000000000000000", &
      "wll_scale = 1.0000000000000000", &
      "en = 2.2000000000000002", &
      "[element.H.shells]", &
      "[element.H.shells.1s]", &
      "level = -10.7072109999999991", &
      "kcn = -0.0500000000000000", &
      "shpoly = -0.0095361800000000", &
      "refocc = 1.0000000000000000", &
      "lgam = 1.0000000000000000", &
      "[element.H.shells.1s.sto_ng]", &
      "slater = 1.2300000000000000", &
      "ngauss = 3", &
      "[element.C]", &
      "zeff = 4.2310780000000001", &
      "arep = 1.2476550000000000", &
      "gam = 0.5380150000000000", &
      "gam3 = 0.1500000000000000", &
      "dkernel = -0.0041167400000000", &
      "qkernel = 0.0021358300000000", &
      "mprad = 3.0000000000000000", &
      "mpvcn = 3.0000000000000000", &
      "wll_scale = 1.0000000000000000", &
      "en = 2.5499999999999998", &
      "[element.C.shells]", &
      "[element.C.shells.2s]", &
      "level = -13.9709220000000016", &
      "kcn = -0.0102144000000000", &
      "shpoly = -0.0229432100000000", &
      "refocc = 1.0000000000000000", &
      "lgam = 1.0000000000000000", &
      "[element.C.shells.2s.sto_ng]", &
      "slater = 2.0964320000000001", &
      "ngauss = 4", &
      "[element.C.shells.2p]", &
      "level = -10.0632920000000006", &
      "kcn = 0.0161657000000000", &
      "shpoly = -0.0027110200000000", &
      "refocc = 3.0000000000000000", &
      "lgam = 1.1056357999999999", &
      "[element.C.shells.2p.sto_ng]", &
      "slater = 1.8000000000000000", &
      "ngauss = 4", &
      ""
   rewind io

   call param%load(io, error)
   close(io)
end subroutine test_param_minimal


subroutine test_param_invalid(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: io
   type(param_record) :: param

   open(newunit=io, status="scratch")
   write(io, '(a)') &
      "[hamiltonian.xtb]", &
      "wexp = 5.0000000000000000E-01", &
      "kpol = 2.0000000000000000E+00", &
      "enscale = 2.0000000000000000E-02", &
      "cn = ""dexp""", &
      "shell = {ss=1.85, pp=2.23, dd=2.23, sd=2.00, pd=2.00"
   rewind io

   call param%load(io, error)
   close(io)
end subroutine test_param_invalid


subroutine test_element_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(element_record) :: param

   table = toml_table()
   table%key = "Te"
   call param%load(table, error)
end subroutine test_element_empty

subroutine test_mol_multipole_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(molecular_multipole_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_mol_multipole_empty


subroutine test_acp_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(acp_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_acp_empty


subroutine test_charge_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(charge_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_charge_empty


subroutine test_dispersion_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(dispersion_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_dispersion_empty


subroutine test_exchange_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(exchange_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_exchange_empty


subroutine test_firstorder_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(firstorder_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_firstorder_empty


subroutine test_fourthorder_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(fourthorder_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_fourthorder_empty


subroutine test_halogen_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(halogen_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_halogen_empty


subroutine test_hamiltonian_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(hamiltonian_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_hamiltonian_empty


subroutine test_increment_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(increment_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_increment_empty


subroutine test_multipole_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(multipole_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_multipole_empty


subroutine test_repulsion_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(repulsion_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_repulsion_empty


subroutine test_spin_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(spin_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_spin_empty


subroutine test_thirdorder_empty(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table) :: table
   type(thirdorder_record) :: param

   table = toml_table()
   call param%load(table, error)
end subroutine test_thirdorder_empty


subroutine test_mask_gfn2(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(param_record), target :: param, base
   type(param_mask) :: mask1, mask2
   real(wp), allocatable :: array(:)
   type(toml_table) :: table
   integer :: io

   call export_gen_param("gfn2", base)

   mask1%ref => base%record
   open(newunit=io, status="scratch")
   write(io, '(a)') &
      "hamiltonian = {}", &
      "dispersion = {}", &
      "repulsion = {}", &
      "charge = {}", &
      "thirdorder = {}", &
      "multipole = {}", &
      "[element]", &
      "[element.H]", &
      "increment = false", &
      "zeff = false", &
      "k0 = false", &
      "k1 = false", &
      "k2 = false", &
      "k3 = false", &
      "k0_h0_scale = false", &
      "k1_h0_scale = false", &
      "k2_h0_scale = false", &
      "k3_h0_scale = false", &
      "cn_rcov = false", &
      "cn_avg = false", &
      "h0_rad = false", &
      "h0_dip_scale = false", &
      "h0_diat_scale_sig = false", &
      "h0_diat_scale_pi = false", &
      "h0_diat_scale_del = false", &
      "rvdw_scale = false", &
      "rep_cn = false", &
      "rep_q = false", &
      "rep_q2 = false", &
      "rep_roffset = false", &
      "rep_k1 = false", &
      "ipea_cn = false", &
      "gam_cn = false", &
      "gam4 = false", &
      "aes_dip_scale = false", &
      "wll_scale = false", &
      "cscale = false", &
      "crad = false", &
      "[element.H.shells]", &
      "[element.H.shells.1s]", &
      "shpoly2 = false", &
      "shpoly4 = false", &
      "h0_exp_scale = false", &
      "ipea = false", &
      "lgam_fx = false", &
      "avg_exp_fx = false", &
      "[element.H.shells.1s.sto_ng]", &
      "slater = false", &
      "[element.C]", &
      "increment = false", &
      "k0 = false", &
      "k1 = false", &
      "k2 = false", &
      "k3 = false", &
      "k0_h0_scale = false", &
      "k1_h0_scale = false", &
      "k2_h0_scale = false", &
      "k3_h0_scale = false", &
      "cn_rcov = false", &
      "cn_avg = false", &
      "h0_rad = false", &
      "h0_dip_scale = false", &
      "h0_diat_scale_sig = false", &
      "h0_diat_scale_pi = false", &
      "h0_diat_scale_del = false", &
      "rvdw_scale = false", &
      "rep_cn = false", &
      "rep_q = false", &
      "rep_q2 = false", &
      "rep_roffset = false", &
      "rep_k1 = false", &
      "ipea_cn = false", &
      "gam_cn = false", &
      "gam4 = false", &
      "aes_dip_scale = false", &
      "wll_scale = false", &
      "cscale = false", &
      "crad = false", &
      "[element.C.shells]", &
      "[element.C.shells.2s]", &
      "shpoly2 = false", &
      "shpoly4 = false", &
      "h0_exp_scale = false", &
      "ipea = false", &
      "lgam_fx = false", &
      "avg_exp_fx = false", &
      "[element.C.shells.2s.sto_ng]", &
      "[element.C.shells.2p]", &
      "shpoly2 = false", &
      "shpoly4 = false", &
      "h0_exp_scale = false", &
      "ipea = false", &
      "lgam_fx = false", &
      "avg_exp_fx = false", &
      "[element.C.shells.2p.sto_ng]", &
      ""
   rewind io

   call mask1%load(io, error)
   close(io)
   if (allocated(error)) return

   call check(error, count(mask1), 25)
   if (allocated(error)) return

   table = toml_table()
   call mask1%dump(table, error)
   if (allocated(error)) return

   mask2%ref => base%record
   call mask2%load(table, error)
   if (allocated(error)) return

   call check(error, count(mask2), 25)
   if (allocated(error)) return

   allocate(array(count(mask2)))
   call base%dump(array, mask2, error)
   if (allocated(error)) return

   call param%load(array, base, mask2, error)
   if (allocated(error)) return

end subroutine test_mask_gfn2


subroutine test_mask_gfn1(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(param_record), target :: param, base
   type(param_mask) :: mask1, mask2
   real(wp), allocatable :: array(:)
   type(toml_table) :: table
   integer :: io

   call export_gen_param("gfn1", base)

   mask1%ref => base%record
   open(newunit=io, status="scratch")
   write(io, '(a)') &
      "hamiltonian = {}", &
      "dispersion = {}", &
      "repulsion = {}", &
      "charge = {}", &
      "thirdorder = {}", &
      "halogen = {}", &
      "[element]", &
      "[element.H]", &
      "increment = false", &
      "zeff = false", &
      "k0 = false", &
      "k1 = false", &
      "k2 = false", &
      "k3 = false", &
      "k0_h0_scale = false", &
      "k1_h0_scale = false", &
      "k2_h0_scale = false", &
      "k3_h0_scale = false", &
      "cn_rcov = false", &
      "h0_rad = false", &
      "h0_dip_scale = false", &
      "h0_diat_scale_sig = false", &
      "h0_diat_scale_pi = false", &
      "h0_diat_scale_del = false", &
      "rvdw_scale = false", &
      "rep_cn = false", &
      "rep_q = false", &
      "rep_q2 = false", &
      "rep_roffset = false", &
      "rep_k1 = false", &
      "ipea_cn = false", &
      "gam_cn = false", &
      "gam4 = false", &
      "wll_scale = false", &
      "cscale = false", &
      "crad = false", &
      "[element.H.shells]", &
      "[element.H.shells.1s]", &
      "shpoly2 = false", &
      "shpoly4 = false", &
      "h0_exp_scale = false", &
      "ipea = false", &
      "lgam_fx = false", &
      "avg_exp_fx = false", &
      "[element.H.shells.1s.sto_ng]", &
      "slater = false", &
      "[element.H.shells.2s]", &
      "shpoly2 = false", &
      "shpoly4 = false", &
      "h0_exp_scale = false", &
      "ipea = false", &
      "lgam_fx = false", &
      "avg_exp_fx = false", &
      "[element.H.shells.2s.sto_ng]", &
      "[element.C]", &
      "increment = false", &
      "k0 = false", &
      "k1 = false", &
      "k2 = false", &
      "k3 = false", &
      "k0_h0_scale = false", &
      "k1_h0_scale = false", &
      "k2_h0_scale = false", &
      "k3_h0_scale = false", &
      "cn_rcov = false", &
      "h0_rad = false", &
      "h0_dip_scale = false", &
      "h0_diat_scale_sig = false", &
      "h0_diat_scale_pi = false", &
      "h0_diat_scale_del = false", &
      "rvdw_scale = false", &
      "rep_cn = false", &
      "rep_q = false", &
      "rep_q2 = false", &
      "rep_roffset = false", &
      "rep_k1 = false", &
      "ipea_cn = false", &
      "gam_cn = false", &
      "gam4 = false", &
      "wll_scale = false", &
      "cscale = false", &
      "crad = false", &
      "[element.C.shells]", &
      "[element.C.shells.2s]", &
      "shpoly = false", &
      "shpoly2 = false", &
      "shpoly4 = false", &
      "h0_exp_scale = false", &
      "ipea = false", &
      "lgam_fx = false", &
      "avg_exp_fx = false", &
      "[element.C.shells.2s.sto_ng]", &
      "[element.C.shells.2p]", &
      "shpoly2 = false", &
      "shpoly4 = false", &
      "h0_exp_scale = false", &
      "ipea = false", &
      "lgam_fx = false", &
      "avg_exp_fx = false", &
      "[element.C.shells.2p.sto_ng]", &
      ""
   rewind io

   call mask1%load(io, error)
   close(io)
   if (allocated(error)) return

   call check(error, count(mask1), 29)
   if (allocated(error)) return

   table = toml_table()
   call mask1%dump(table, error)
   if (allocated(error)) return

   mask2%ref => base%record
   call mask2%load(table, error)
   if (allocated(error)) return

   call check(error, count(mask2), 29)
   if (allocated(error)) return

   allocate(array(count(mask2)))
   call base%dump(array, mask2, error)
   if (allocated(error)) return

   call param%load(array, base, mask2, error)
   if (allocated(error)) return

end subroutine test_mask_gfn1


subroutine test_mask_gxtb(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(param_record), target :: param, base
   type(param_mask) :: mask1, mask2
   real(wp), allocatable :: array(:)
   type(toml_table) :: table
   integer :: io

   call export_gen_param("gxtb", base)

   mask1%ref => base%record
   open(newunit=io, status="scratch")
   write(io, '(a)') &
      "hamiltonian = {}", &
      "acp = {}", &
      "dispersion = {}", &
      "increment = {}", &
      "repulsion = {}", &
      "firstorder = {}", &
      "charge = {}", &
      "thirdorder = {}", &
      "fourthorder = {}", &
      "multipole = {}", &
      "exchange = {}", &
      "spin = {}", &
      "[element]", &
      "[element.H]", &
      "increment = false", &
      "xbond = false", &
      "gam4 = false", &
      "dkernel = false", &
      "qkernel = false", &
      "mprad = false", &
      "mpvcn = false", &
      "en = false", &
      "[element.H.acp]", &
      "[element.H.shells]", &
      "[element.H.shells.1s]", &
      "shpoly = false", &
      "zeffsh = false", &
      "[element.H.shells.1s.q_vszp]", &
      "expos = [true, false, false, false, false, false, false, false]", &
      "coeffs = [false, true, false, false, false, false, false, false]", &
      "coeffs_env = [false, false, true, false, false, false, false, false]", &
      "[element.C]", &
      "xbond = false", &
      "dkernel = false", &
      "qkernel = false", &
      "mprad = false", &
      "mpvcn = false", &
      "en = false", &
      "[element.C.acp]", &
      "[element.C.shells]", &
      "[element.C.shells.2s]", &
      "shpoly = false", &
      "zeffsh = false", &
      "[element.C.shells.2s.q_vszp]", &
      "expos = [false, false, false, false, false, false]", &
      "coeffs = [false, false, false, false, false, false]", &
      "coeffs_env = [false, false, false, false, false, false]", &
      "[element.C.shells.2p]", &
      "shpoly = false", &
      "[element.C.shells.2p.q_vszp]", &
      "expos = [true, true, true, true, true, true]", &
      "coeffs = [true, true, true, true, true, true]", &
      "coeffs_env = [true, true, true, true, true, true]", &
      ""
   rewind io

   call mask1%load(io, error)
   close(io)
   if (allocated(error)) return

   call check(error, count(mask1), 110)
   if (allocated(error)) return

   table = toml_table()
   call mask1%dump(table, error)
   if (allocated(error)) return

   mask2%ref => base%record
   call mask2%load(table, error)
   if (allocated(error)) return

   call check(error, count(mask2), 110)
   if (allocated(error)) return

   allocate(array(count(mask2)))
   call base%dump(array, mask2, error)
   if (allocated(error)) return

   call param%load(array, base, mask2, error)
   if (allocated(error)) return

end subroutine test_mask_gxtb


subroutine export_gen_param(method, param)
   character(len=*), intent(in) :: method
   type(param_record), intent(out) :: param
   select case(method)
   case("gfn1")
      call export_gfn1_param(param)
   case("gfn2")
      call export_gfn2_param(param)
   case("ipea1")
      call export_ipea1_param(param)
   case("gxtb")
      call export_gxtb_param(param)
   end select
end subroutine export_gen_param


subroutine new_gen_calculator(calc, method, mol, error)
   type(xtb_calculator), intent(out) :: calc
   character(len=*), intent(in) :: method
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error
   select case(method)
   case("gfn1")
      call new_gfn1_calculator(calc, mol, error)
   case("gfn2")
      call new_gfn2_calculator(calc, mol, error)
   case("ipea1")
      call new_ipea1_calculator(calc, mol, error)
   case("gxtb")
      call new_gxtb_calculator(calc, mol, error)
   end select
end subroutine new_gen_calculator


subroutine test_gen(mol, method, error)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Method name
   character(len=*), intent(in) :: method
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(param_record) :: param
   type(wavefunction_type) :: wfn
   integer :: nspin
   real(wp) :: energy1, energy2
   
   nspin = 1
   if (mod(mol%uhf, 2) == 1) nspin = 2

   call export_gen_param(method, param)
   call new_xtb_calculator(calc, mol, param, error)
   if (allocated(error)) return

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & nspin, calc%default_etemp * kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy1, verbosity=0)

   call new_gen_calculator(calc, method, mol, error)
   if (allocated(error)) return

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & nspin, calc%default_etemp * kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy2, verbosity=0)

   call check(error, energy2, energy1, thr=thr)

end subroutine test_gen


subroutine test_round_trip(mol, method, error)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Method name
   character(len=*), intent(in) :: method
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(param_record) :: param1, param2
   type(wavefunction_type) :: wfn
   real(wp) :: energy1, energy2
   integer :: nspin
   character(len=:), allocatable :: io

   nspin = 1
   if (mod(mol%uhf, 2) == 1) nspin = 2

   call export_gen_param(method, param1)

   io = "." // method // "-" // get_name() // ".toml"
   call param1%dump(io, error)
   if (.not.allocated(error)) then
      call param2%load(io, error)
   end if
   call delete_file(io)
   if (allocated(error)) return

   call new_xtb_calculator(calc, mol, param1, error)
   if (allocated(error)) return

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & nspin, calc%default_etemp * kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy1, verbosity=0)

   call new_xtb_calculator(calc, mol, param2, error)
   if (allocated(error)) return

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & nspin, calc%default_etemp * kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy2, verbosity=0)

   call check(error, energy2, energy1, thr=thr)

end subroutine test_round_trip


subroutine test_gfn1_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_gen(mol, "gfn1", error)

end subroutine test_gfn1_mb01


subroutine test_gfn2_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "02")
   call test_gen(mol, "gfn2", error)

end subroutine test_gfn2_mb02


subroutine test_ipea1_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   call test_gen(mol, "ipea1", error)

end subroutine test_ipea1_mb03


subroutine test_gfn1_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_gen(mol, "gfn1", error)

end subroutine test_gfn1_mb04


subroutine test_gfn2_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_gen(mol, "gfn2", error)

end subroutine test_gfn2_mb05


subroutine test_ipea1_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "06")
   call test_gen(mol, "ipea1", error)

end subroutine test_ipea1_mb06


subroutine test_gxtb_water(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "ICE10", "gas")
   call test_gen(mol, "gxtb", error)

end subroutine test_gxtb_water


subroutine test_gxtb_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "07")
   call test_gen(mol, "gxtb", error)

end subroutine test_gxtb_mb07


subroutine test_gxtb_mb08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "08")
   call test_gen(mol, "gxtb", error)

end subroutine test_gxtb_mb08


subroutine test_gxtb_mb10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "10")
   call test_gen(mol, "gxtb", error)

end subroutine test_gxtb_mb10

subroutine test_gfn2_round_trip(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "07")
   call test_round_trip(mol, "gfn2", error)

end subroutine test_gfn2_round_trip


subroutine test_gfn1_round_trip(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "08")
   call test_round_trip(mol, "gfn1", error)

end subroutine test_gfn1_round_trip


subroutine test_ipea1_round_trip(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "09")
   call test_round_trip(mol, "ipea1", error)

end subroutine test_ipea1_round_trip


subroutine test_gxtb_round_trip(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "10")
   call test_round_trip(mol, "gxtb", error)

end subroutine test_gxtb_round_trip


function get_name() result(name)
   character(len=20) :: name
   real :: val

   call random_number(val)
   write(name, '(a, z8.8)') "tblite-test-", int(val*1.0e9)
end function get_name


subroutine delete_file(file)
   character(len=*), intent(in) :: file
   integer :: unit
   logical :: exist
   inquire(file=file, exist=exist)
   if (exist) then
      open(newunit=unit, file=file)
      close(unit, status="delete")
   end if
end subroutine delete_file


end module test_xtb_param
