! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

module dftd4_param_functionals
   implicit none
   public

   private :: new_funcgroup

   enum, bind(C)
      enumerator :: p_invalid, p_default, &
         & p_hf, p_blyp, p_bpbe, p_bp, p_bpw, p_lb94, p_mpwlyp, p_mpwpw, &
         & p_olyp, p_opbe, p_pbe, p_rpbe, p_revpbe, p_pw86pbe, &
         & p_rpw86pbe, p_pw91, p_pwp, p_xlyp, p_b97, p_tpss, p_revtpss, &
         & p_scan, p_rscan, p_r2scan, p_b1lyp, p_b3lyp, p_bhlyp, p_b1p, &
         & p_b3p, p_b1pw, p_b3pw, p_o3lyp, p_revpbe0, p_revpbe38, &
         & p_pbe0, p_pwp1, p_pw1pw, p_mpw1pw, p_mpw1lyp, p_pw6b95, &
         & p_tpssh, p_tpss0, p_x3lyp, p_m06l, p_m06, p_b97d, &
         & p_wb97, p_wb97x_2008, p_b97m, p_wb97m, p_camb3lyp, p_lcblyp, &
         & p_lh07tsvwn, p_lh07ssvwn, p_lh12ctssirpw92, p_lh12ctssifpw92, &
         & p_lh14tcalpbe, p_lh20t, p_b2plyp, p_b2gpplyp, p_mpw2plyp, p_pwpb95, &
         & p_dsdblyp, p_dsdpbe, p_dsdpbeb95, p_dsdpbep86, p_dsdsvwn, &
         & p_dodblyp, p_dodpbe, p_dodpbeb95, p_dodpbep86, p_dodsvwn, &
         & p_pbe0_2, p_pbe0_dh, p_hsesol, p_dftb_3ob, p_dftb_mio, p_dftb_ob2, &
         & p_dftb_matsci, p_dftb_pbc, p_b1b95, p_pbesol, p_hse06, p_mpwb1k, &
         & p_hse03, p_revtpssh, p_mn12sx, p_glyp, p_mpw1b95, &
         & p_revpbe0dh, p_revtpss0, p_revdsdpbep86, p_revdsdpbe, &
         & p_revdsdblyp, p_revdodpbep86, p_am05, p_hse12, p_hse12s, &
         & p_r2scanh, p_r2scan0, p_r2scan50, p_r2scan_3c, p_camqtp01, &
         & p_lcwpbe, p_lcwpbeh, p_wb97x_rev, p_wb97m_rev, &
         & p_wb97x_3c, p_wr2scan, p_r2scan0_dh, p_r2scan_cidh, &
         & p_r2scan_qidh, p_r2scan0_2, p_pr2scan50, p_pr2scan69, &
         & p_kpr2scan50, p_wpr2scan50, p_wb97x, p_last
   end enum
   integer, parameter :: df_enum = kind(p_invalid)


   !> Group different spellings/names of functionals
   type functional_group
      character(len=:), allocatable :: names(:)
   end type functional_group

contains


!> Create a new group of functional names
function new_funcgroup(input_names) result(group)

   !> List of spellings/names of the functional
   character(len=*), intent(in) :: input_names(:)

   !> Functional with possibly different spellings
   type(functional_group) :: group
   
   integer :: n, i, max_len
   n = size(input_names)

   ! Determine the length of the longest name
   max_len = 0
   do i = 1, n
      max_len = max(max_len, len_trim(input_names(i)))
   end do

   ! Allocate based on the longest name's length
   allocate(character(len=max_len) :: group%names(n))
   do i = 1, n
      group%names(i) = trim(input_names(i))
   end do
end function new_funcgroup


!> Collect all supported functionals
subroutine get_functionals(funcs)
   !DEC$ ATTRIBUTES DLLEXPORT :: get_functionals

   !> Collection of functionals with possibly different spellings/names
   type(functional_group), allocatable, intent(out) :: funcs(:)

   allocate(funcs(p_last - 1))

   funcs(p_hf) = new_funcgroup([character(len=20) :: 'hf'])
   funcs(p_am05) = new_funcgroup([character(len=20) :: 'am05'])
   funcs(p_blyp) = new_funcgroup([character(len=20) :: 'b-lyp', 'blyp'])
   funcs(p_bpbe) = new_funcgroup([character(len=20) :: 'bpbe'])
   funcs(p_bp) = new_funcgroup([character(len=20) :: 'b-p', 'bp86', 'bp', 'b-p86'])
   funcs(p_bpw) = new_funcgroup([character(len=20) :: 'bpw', 'b-pw'])
   funcs(p_lb94) = new_funcgroup([character(len=20) :: 'lb94'])
   funcs(p_mpwlyp) = new_funcgroup([character(len=20) :: 'mpwlyp', 'mpw-lyp'])
   funcs(p_mpwpw) = new_funcgroup([character(len=20) :: 'mpwpw', 'mpw-pw', 'mpwpw91'])
   funcs(p_olyp) = new_funcgroup([character(len=20) :: 'o-lyp', 'olyp'])
   funcs(p_opbe) = new_funcgroup([character(len=20) :: 'opbe'])
   funcs(p_pbe) = new_funcgroup([character(len=20) :: 'pbe'])
   funcs(p_rpbe) = new_funcgroup([character(len=20) :: 'rpbe'])
   funcs(p_revpbe) = new_funcgroup([character(len=20) :: 'revpbe'])
   funcs(p_pw86pbe) = new_funcgroup([character(len=20) :: 'pw86pbe'])
   funcs(p_rpw86pbe) = new_funcgroup([character(len=20) :: 'rpw86pbe'])
   funcs(p_pw91) = new_funcgroup([character(len=20) :: 'pw91'])
   funcs(p_pwp) = new_funcgroup([character(len=20) :: 'pwp', 'pw-p', 'pw91p86'])
   funcs(p_xlyp) = new_funcgroup([character(len=20) :: 'x-lyp', 'xlyp'])
   funcs(p_b97) = new_funcgroup([character(len=20) :: 'b97'])
   funcs(p_tpss) = new_funcgroup([character(len=20) :: 'tpss'])
   funcs(p_revtpss) = new_funcgroup([character(len=20) :: 'revtpss'])
   funcs(p_scan) = new_funcgroup([character(len=20) :: 'scan'])
   funcs(p_rscan) = new_funcgroup([character(len=20) :: 'rscan'])
   funcs(p_r2scan) = new_funcgroup([character(len=20) :: 'r2scan', 'r²scan'])
   funcs(p_r2scanh) = new_funcgroup([character(len=20) :: 'r2scanh', 'r²scanh'])
   funcs(p_r2scan0) = new_funcgroup([character(len=20) :: 'r2scan0', 'r²scan0'])
   funcs(p_r2scan50) = new_funcgroup([character(len=20) :: 'r2scan50', 'r²scan50'])
   funcs(p_r2scan_3c) = new_funcgroup([character(len=20) :: 'r2scan-3c', &
      & 'r²scan-3c', 'r2scan_3c', 'r²scan_3c', 'r2scan3c'])
   funcs(p_wr2scan) = new_funcgroup([character(len=20) :: 'wr2scan', 'wr²scan'])
   funcs(p_r2scan0_dh) = new_funcgroup([character(len=20) :: 'r2scan0-dh', &
      & 'r²scan0-dh', 'r2scan0dh', 'r²scan0dh'])
   funcs(p_r2scan_cidh) = new_funcgroup([character(len=20) :: 'r2scan-cidh', &
      & 'r²scan-cidh', 'r2scancidh', 'r²scancidh'])
   funcs(p_r2scan_qidh) = new_funcgroup([character(len=20) :: 'r2scan-qidh', &  
      & 'r²scan-qidh', 'r2scanqidh', 'r²scanqidh'])
    funcs(p_r2scan0_2) = new_funcgroup([character(len=20) :: 'r2scan0-2', &
      & 'r²scan0-2', 'r2scan02', 'r²scan02'])
   funcs(p_pr2scan50) = new_funcgroup([character(len=20) :: 'pr2scan50', &
      & 'pr²scan50', 'pr2scan50', 'pr²scan50'])
   funcs(p_pr2scan69) = new_funcgroup([character(len=20) :: 'pr2scan69', &
      & 'pr²scan69', 'pr2scan69', 'pr²scan69'])
   funcs(p_kpr2scan50) = new_funcgroup([character(len=20) :: 'kpr2scan50', & 
      & 'kpr²scan50', 'kpr2scan50', 'kpr²scan50'])
   funcs(p_wpr2scan50) = new_funcgroup([character(len=20) :: 'wpr2scan50', &
      & 'wpr²scan50', 'wpr2scan50', 'wpr²scan50'])
   funcs(p_b1lyp) = new_funcgroup([character(len=20) :: 'b1lyp', 'b1-lyp'])
   funcs(p_b3lyp) = new_funcgroup([character(len=20) :: 'b3-lyp', 'b3lyp'])
   funcs(p_bhlyp) = new_funcgroup([character(len=20) :: 'bh-lyp', 'bhlyp'])
   funcs(p_b1p) = new_funcgroup([character(len=20) :: 'b1p', 'b1-p', 'b1p86'])
   funcs(p_b3p) = new_funcgroup([character(len=20) :: 'b3p', 'b3-p', 'b3p86'])
   funcs(p_b1pw) = new_funcgroup([character(len=20) :: 'b1pw', 'b1-pw', 'b1pw91'])
   funcs(p_b3pw) = new_funcgroup([character(len=20) :: 'b3pw', 'b3-pw', 'b3pw91'])
   funcs(p_o3lyp) = new_funcgroup([character(len=20) :: 'o3-lyp', 'o3lyp'])
   funcs(p_revpbe0) = new_funcgroup([character(len=20) :: 'revpbe0'])
   funcs(p_revpbe38) = new_funcgroup([character(len=20) :: 'revpbe38'])
   funcs(p_pbe0) = new_funcgroup([character(len=20) :: 'pbe0'])
   funcs(p_pwp1) = new_funcgroup([character(len=20) :: 'pwp1'])
   funcs(p_pw1pw) = new_funcgroup([character(len=20) :: 'pw1pw', 'pw1-pw'])
   funcs(p_mpw1pw) = new_funcgroup([character(len=20) :: 'mpw1pw', 'mpw1-pw', 'mpw1pw91'])
   funcs(p_mpw1lyp) = new_funcgroup([character(len=20) :: 'mpw1lyp', 'mpw1-lyp'])
   funcs(p_pw6b95) = new_funcgroup([character(len=20) :: 'pw6b95'])
   funcs(p_tpssh) = new_funcgroup([character(len=20) :: 'tpssh'])
   funcs(p_tpss0) = new_funcgroup([character(len=20) :: 'tpss0'])
   funcs(p_x3lyp) = new_funcgroup([character(len=20) :: 'x3-lyp', 'x3lyp'])
   funcs(p_m06) = new_funcgroup([character(len=20) :: 'm06'])
   funcs(p_m06l) = new_funcgroup([character(len=20) :: 'm06l'])
   funcs(p_mn12sx) = new_funcgroup([character(len=20) :: 'mn12sx', 'mn12-sx'])
   funcs(p_b97d) = new_funcgroup([character(len=20) :: 'b97d'])
   funcs(p_lh07tsvwn) = new_funcgroup([character(len=20) :: 'lh07tsvwn', 'lh07t-svwn'])
   funcs(p_lh07ssvwn) = new_funcgroup([character(len=20) :: 'lh07ssvwn', 'lh07s-svwn'])
   funcs(p_lh12ctssirpw92) = new_funcgroup([character(len=20) :: 'lh12ctssirpw92', 'lh12ct-ssirpw92'])
   funcs(p_lh12ctssifpw92) = new_funcgroup([character(len=20) :: 'lh12ctssifpw92', 'lh12ct-ssifpw92'])
   funcs(p_lh14tcalpbe) = new_funcgroup([character(len=20) :: 'lh14tcalpbe', 'lh14t-calpbe'])
   funcs(p_lh20t) = new_funcgroup([character(len=20) :: 'lh20t'])
   funcs(p_b2plyp) = new_funcgroup([character(len=20) :: 'b2plyp', 'b2-plyp'])
   funcs(p_b2gpplyp) = new_funcgroup([character(len=20) :: 'b2gpplyp', 'b2gp-plyp'])
   funcs(p_mpw2plyp) = new_funcgroup([character(len=20) :: 'mpw2plyp'])
   funcs(p_pwpb95) = new_funcgroup([character(len=20) :: 'pwpb95'])
   funcs(p_dsdblyp) = new_funcgroup([character(len=20) :: 'dsdblyp', 'dsd-blyp'])
   funcs(p_dsdpbe) = new_funcgroup([character(len=20) :: 'dsdpbe', 'dsd-pbe'])
   funcs(p_dsdpbeb95) = new_funcgroup([character(len=20) :: 'dsdpbeb95', 'dsd-pbeb95'])
   funcs(p_dsdpbep86) = new_funcgroup([character(len=20) :: 'dsdpbep86', 'dsd-pbep86'])
   funcs(p_dsdsvwn) = new_funcgroup([character(len=20) :: 'dsdsvwn', 'dsd-svwn'])
   funcs(p_dodblyp) = new_funcgroup([character(len=20) :: 'dodblyp', 'dod-blyp'])
   funcs(p_dodpbe) = new_funcgroup([character(len=20) :: 'dodpbe', 'dod-pbe'])
   funcs(p_dodpbeb95) = new_funcgroup([character(len=20) :: 'dodpbeb95', 'dod-pbeb95'])
   funcs(p_dodpbep86) = new_funcgroup([character(len=20) :: 'dodpbep86', 'dod-pbep86'])
   funcs(p_dodsvwn) = new_funcgroup([character(len=20) :: 'dodsvwn', 'dod-svwn'])
   funcs(p_pbe0_2) = new_funcgroup([character(len=20) :: 'pbe02', 'pbe0-2'])
   funcs(p_pbe0_dh) = new_funcgroup([character(len=20) :: 'pbe0dh', 'pbe0-dh'])
   funcs(p_dftb_3ob) = new_funcgroup([character(len=20) :: 'dftb3', 'dftb(3ob)'])
   funcs(p_dftb_mio) = new_funcgroup([character(len=20) :: 'dftb(mio)'])
   funcs(p_dftb_pbc) = new_funcgroup([character(len=20) :: 'dftb(pbc)'])
   funcs(p_dftb_matsci) = new_funcgroup([character(len=20) :: 'dftb(matsci)'])
   funcs(p_dftb_ob2) = new_funcgroup([character(len=20) :: 'lc-dftb', 'dftb(ob2)'])
   funcs(p_b1b95) = new_funcgroup([character(len=20) :: 'b1b95'])
   funcs(p_pbesol) = new_funcgroup([character(len=20) :: 'pbesol'])
   funcs(p_mpwb1k) = new_funcgroup([character(len=20) :: 'mpwb1k'])
   funcs(p_mpw1b95) = new_funcgroup([character(len=20) :: 'mpw1b95'])
   funcs(p_hse03) = new_funcgroup([character(len=20) :: 'hse03'])
   funcs(p_hse06) = new_funcgroup([character(len=20) :: 'hse06'])
   funcs(p_hse12) = new_funcgroup([character(len=20) :: 'hse12'])
   funcs(p_hse12s) = new_funcgroup([character(len=20) :: 'hse12s'])
   funcs(p_hsesol) = new_funcgroup([character(len=20) :: 'hsesol'])
   funcs(p_revtpssh) = new_funcgroup([character(len=20) :: 'revtpssh'])
   funcs(p_glyp) = new_funcgroup([character(len=20) :: 'glyp', 'g-lyp'])
   funcs(p_revpbe0dh) = new_funcgroup([character(len=20) :: 'revpbe0dh', 'revpbe0-dh'])
   funcs(p_revtpss0) = new_funcgroup([character(len=20) :: 'revtpss0'])
   funcs(p_revdsdpbep86) = new_funcgroup([character(len=20) :: 'revdsd-pbep86', 'revdsdpbep86'])
   funcs(p_revdsdpbe) = new_funcgroup([character(len=20) :: 'revdsd-pbe', 'revdsd-pbepbe', 'revdsdpbe', 'revdsdpbepbe'])
   funcs(p_revdsdblyp) = new_funcgroup([character(len=20) :: 'revdsd-blyp', 'revdsdblyp'])
   funcs(p_revdodpbep86) = new_funcgroup([character(len=20) :: 'revdod-pbep86', 'revdodpbep86'])
   funcs(p_b97m) = new_funcgroup([character(len=20) :: 'b97m'])
   funcs(p_wb97m) = new_funcgroup([character(len=20) :: 'wb97m', 'ωb97m', 'omegab97m'])
   funcs(p_wb97m_rev) = new_funcgroup([character(len=20) :: 'wb97m-rev', &
      & 'ωb97m-rev', 'omegab97m-rev', 'wb97m_rev', 'ωb97m_rev', 'omegab97m_rev'])
   funcs(p_wb97) = new_funcgroup([character(len=20) :: 'wb97', 'ωb97', 'omegab97'])
   funcs(p_wb97x_2008) = new_funcgroup([character(len=20) :: 'wb97x_2008', &
      & 'ωb97x_2008', 'omegab97x_2008', 'wb97x-2008', 'ωb97x-2008', &
      & 'omegab97x-2008'])
   funcs(p_wb97x) = new_funcgroup([character(len=20) :: 'wb97x', 'ωb97x', &
      & 'omegab97x'])
   funcs(p_wb97x_rev) = new_funcgroup([character(len=20) :: 'wb97x-rev', &
      & 'ωb97x-rev', 'omegab97x-rev', 'wb97x_rev', 'ωb97x_rev', 'omegab97x_rev'])
   funcs(p_wb97x_3c) = new_funcgroup([character(len=20) :: 'wb97x-3c', &
      & 'ωb97x-3c', 'omegab97x-3c', 'wb97x_3c', 'ωb97x_3c', 'omegab97x_3c'])
   funcs(p_camb3lyp) = new_funcgroup([character(len=20) :: 'cam-b3lyp', 'camb3lyp'])
   funcs(p_camqtp01) = new_funcgroup([character(len=20) :: 'cam-qtp01', &
      & 'camqtp01', 'camqtp(01)', 'cam-qtp(01)'])
   funcs(p_lcblyp) = new_funcgroup([character(len=20) :: 'lc-blyp', 'lcblyp'])
   funcs(p_lcwpbe) = new_funcgroup([character(len=20) :: 'lc-wpbe', &
      & 'lcwpbe', 'lc-ωpbe', 'lcωpbe', 'lc-omegapbe', 'lcomegapbe'])
   funcs(p_lcwpbeh) = new_funcgroup([character(len=20) :: 'lc-wpbeh', &
      & 'lcwpbeh', 'lc-ωpbeh', 'lcωpbeh', 'lc-omegapbeh', 'lcomegapbeh'])

end subroutine get_functionals


!> Get the unique identifier for most functionals, returns none if
!> the functional was not known at the time I implemented this mapping
pure function get_functional_id(df) result(num)
   integer(df_enum) :: num
   character(len=*), intent(in) :: df
   select case(df)
   case default
      num = p_invalid
   case('default')
      num = p_default
   case('hf')
      num = p_hf
   case('am05', 'gga_x_am05:gga_c_am05')
      num = p_am05
   case('b-lyp', 'blyp', 'gga_x_b88:gga_c_lyp')
      num = p_blyp
   case('bpbe', 'gga_x_b88:gga_c_pbe')
      num = p_bpbe
   case('b-p', 'bp86', 'bp', 'b-p86', 'gga_x_b88:gga_c_p86')
      num = p_bp
   case('bpw', 'b-pw', 'gga_x_b88:gga_c_pw91')
      num = p_bpw
   case('lb94', 'gga_x_lb') ! no gga_c_lb
      num = p_lb94
   case('mpwlyp', 'mpw-lyp', 'gga_x_mpw91:gga_c_lyp')
      num = p_mpwlyp
   case('mpwpw', 'mpw-pw', 'mpwpw91', 'gga_x_mpw91:gga_c_pw91')
      num = p_mpwpw
   case('o-lyp', 'olyp', 'gga_x_optx:gga_c_lyp')
      num = p_olyp
   case('opbe', 'gga_x_optx:gga_c_pbe')
      num = p_opbe
   case('pbe', 'gga_x_pbe:gga_c_pbe')
      num = p_pbe
   case('rpbe', 'gga_x_rpbe:gga_c_pbe')
      num = p_rpbe
   case('revpbe', 'gga_x_pbe_r:gga_c_pbe')
      num = p_revpbe
   case('pbesol', 'gga_x_pbe_sol:gga_c_pbe_sol')
      num = p_pbesol
   case('pw86pbe', 'gga_x_pw86:gga_c_pbe')
      num = p_pw86pbe
   case('rpw86pbe', 'gga_x_rpw86:gga_c_pbe')
      num = p_rpw86pbe
   case('pw91', 'gga_x_pw91:gga_c_pw91')
      num = p_pw91
   case('pwp', 'pw-p', 'pw91p86', 'gga_x_pw91:gga_c_p86')
      num = p_pwp
   case('x-lyp', 'xlyp', 'gga_xc_xlyp')
      num = p_xlyp
   case('b97', 'hyb_gga_xc_b97')
      num = p_b97
   case('b97d', 'gga_xc_b97_d')
      num = p_b97d
   case('tpss', 'mgga_c_tpss:mgga_x_tpss')
      num = p_tpss
   case('revtpss', 'mgga_c_revtpss:mgga_x_revtpss')
      num = p_revtpss
   case('scan', 'mgga_x_scan:mgga_c_scan')
      num = p_scan
   case('rscan', 'mgga_x_rscan:mgga_c_rscan')
      num = p_rscan
   case('r2scan', 'r²scan', 'mgga_x_r2scan:mgga_c_r2scan')
      num = p_r2scan
   case('r2scanh', 'r²scanh', 'hyb_mgga_xc_r2scanh')
      num = p_r2scanh
   case('r2scan0', 'r²scan0', 'hyb_mgga_xc_r2scan0')
      num = p_r2scan0
   case('r2scan50', 'r²scan50', 'hyb_mgga_xc_r2scan50')
      num = p_r2scan50
   case('r2scan-3c', 'r²scan-3c', 'r2scan_3c', 'r²scan_3c', 'r2scan3c')
      num = p_r2scan_3c
   case('b1lyp', 'b1-lyp', 'hyb_gga_xc_b1lyp')
      num = p_b1lyp
   case('b3-lyp', 'b3lyp', 'hyb_gga_xc_b3lyp', 'hyb_gga_xc_b3lyp3', 'hyb_gga_xc_b3lyp5')
      num = p_b3lyp
   case('bh-lyp', 'bhlyp', 'hyb_gga_xc_bhandh', 'hyb_gga_xc_bhandhlyp')
      num = p_bhlyp
   case('b1p', 'b1-p', 'b1p86') ! 0.75 b88 + 0.25 hf; p86 (nonloc) + pw81 (loc)
      num = p_b1p
   case('b3p', 'b3-p', 'b3p86', 'hyb_gga_xc_b3p86', 'hyb_gga_xc_b3p86_nwchem')
      num = p_b3p
   case('b1pw', 'b1-pw', 'b1pw91', 'hyb_gga_xc_b1pw91')
      num = p_b1pw
   case('b3pw', 'b3-pw', 'b3pw91', 'hyb_gga_xc_b3pw91')
      num = p_b3pw
   case('o3-lyp', 'o3lyp', 'hyb_gga_xc_o3lyp')
      num = p_o3lyp
   case('revpbe0') ! no libxc
      num = p_revpbe0
   case('revpbe38') ! no libxc
      num = p_revpbe38
   case('pbe0', 'hyb_gga_xc_pbeh')
      num = p_pbe0
   case('pwp1') ! no libxc
      num = p_pwp1
   case('pw1pw', 'pw1-pw') ! no libxc
      num = p_pw1pw
   case('mpw1pw', 'mpw1-pw', 'mpw1pw91', 'hyb_gga_xc_mpw1pw')
      num = p_mpw1pw
   case('mpw1lyp', 'mpw1-lyp', 'hyb_gga_xc_mpw1lyp')
      num = p_mpw1lyp
   case('pw6b95', 'hyb_mgga_xc_pw6b95')
      num = p_pw6b95
   case('tpssh', 'hyb_mgga_xc_tpssh')
      num = p_tpssh
   case('tpss0', 'hyb_mgga_xc_tpss0')
      num = p_tpss0
   case('x3-lyp', 'x3lyp', 'hyb_gga_xc_x3lyp')
      num = p_x3lyp
   case('m06', 'mgga_x_m06:mgga_c_m06')
      num = p_m06
   case('m06l', 'mgga_x_m06_l:mgga_c_m06_l')
      num = p_m06l
   case('mn12sx', 'mn12-sx', 'mgga_c_mn12_sx:mgga_c_mn12_sx')
      num = p_mn12sx
   case('cam-b3lyp', 'camb3lyp', 'hyb_gga_xc_cam_b3lyp')
      num = p_camb3lyp
   case('cam-qtp01', 'camqtp01', 'camqtp(01)', 'cam-qtp(01)', &
      & 'hyb_gga_xc_cam_qtp_01')
      num = p_camqtp01
   case('lc-blyp', 'lcblyp', 'hyb_gga_xc_lc_blyp')
      num = p_lcblyp
   case('lc-wpbe', 'lcwpbe', 'lc-ωpbe', 'lcωpbe', 'lc-omegapbe', 'lcomegapbe', &
      & 'hyb_gga_xc_lc_wpbe', 'hyb_gga_xc_lc_wpbe08_whs', &
      & 'hyb_gga_xc_lc_wpbe_whs', 'hyb_gga_xc_lrc_wpbe')
      num = p_lcwpbe
   case('lc-wpbeh', 'lcwpbeh', 'lc-ωpbeh', 'lcωpbeh', 'lc-omegapbeh', &
      & 'lcomegapbeh', 'hyb_gga_xc_lc_wpbeh_whs', 'hyb_gga_xc_lrc_wpbeh')
      num = p_lcwpbeh
   case('lh07tsvwn', 'lh07t-svwn') ! no libxc
      num = p_lh07tsvwn
   case('lh07ssvwn', 'lh07s-svwn') ! no libxc
      num = p_lh07ssvwn
   case('lh12ctssirpw92', 'lh12ct-ssirpw92') ! no libxc
      num = p_lh12ctssirpw92
   case('lh12ctssifpw92', 'lh12ct-ssifpw92') ! no libxc
      num = p_lh12ctssifpw92
   case('lh14tcalpbe', 'lh14t-calpbe') ! no libxc
      num = p_lh14tcalpbe
   case('lh20t') ! no libxc
      num = p_lh20t
   case('b2plyp', 'b2-plyp', 'xc_hyb_gga_xc_b2plyp') ! only in code
      num = p_b2plyp
   case('b2gpplyp', 'b2gp-plyp', 'xc_hyb_gga_xc_b2gpplyp') ! only in code
      num = p_b2gpplyp
   case('mpw2plyp') ! no libxc
      num = p_mpw2plyp
   case('pwpb95') ! no libxc
      num = p_pwpb95
   case('dsdblyp', 'dsd-blyp') ! no libxc
      num = p_dsdblyp
   case('dsdpbe', 'dsd-pbe') ! no libxc
      num = p_dsdpbe
   case('dsdpbeb95', 'dsd-pbeb95') ! no libxc
      num = p_dsdpbeb95
   case('dsdpbep86', 'dsd-pbep86') ! no libxc
      num = p_dsdpbep86
   case('dsdsvwn', 'dsd-svwn') ! no libxc
      num = p_dsdsvwn
   case('dodblyp', 'dod-blyp') ! no libxc
      num = p_dodblyp
   case('dodpbe', 'dod-pbe') ! no libxc
      num = p_dodpbe
   case('dodpbeb95', 'dod-pbeb95') ! no libxc
      num = p_dodpbeb95
   case('dodpbep86', 'dod-pbep86') ! no libxc
      num = p_dodpbep86
   case('dodsvwn', 'dod-svwn') ! no libxc
      num = p_dodsvwn
   case('pbe02', 'pbe0-2') ! no libxc
      num = p_pbe0_2
   case('pbe0dh', 'pbe0-dh') ! no libxc
      num = p_pbe0_dh
   case('dftb3', 'dftb(3ob)') ! no libxc
      num = p_dftb_3ob
   case('dftb(mio)') ! no libxc
      num = p_dftb_mio
   case('dftb(pbc)') ! no libxc
      num = p_dftb_pbc
   case('dftb(matsci)') ! no libxc
      num = p_dftb_matsci
   case('lc-dftb', 'dftb(ob2)') ! no libxc
      num = p_dftb_ob2
   case('b1b95', 'hyb_mgga_xc_b88b95')
      num = p_b1b95
   case('mpwb1k', 'hyb_mgga_xc_mpwb1k')
      num = p_mpwb1k
   case('mpw1b95', 'hyb_mgga_xc_mpw1b95')
      num = p_mpw1b95
   case('hse03', 'hyb_gga_xc_hse03')
      num = p_hse03
   case('hse06', 'hyb_gga_xc_hse06')
      num = p_hse06
   case('hse12', 'hyb_gga_xc_hse12')
      num = p_hse12
   case('hse12s', 'hyb_gga_xc_hse12s')
      num = p_hse12s
   case('hsesol', 'hyb_gga_xc_hse_sol')
      num = p_hsesol
   case('glyp', 'g-lyp', 'gga_x_g96:gga_c_lyp')
      num = p_glyp
   case('revpbe0dh', 'revpbe0-dh') ! no libxc
      num = p_revpbe0dh
   case('revtpssh', 'hyb_mgga_xc_revtpssh')
      num = p_revtpssh
    case('revtpss0') ! no libxc
      num = p_revtpss0
   case('revdsd-pbep86', 'revdsdpbep86') ! no libxc
      num = p_revdsdpbep86
   case('revdsd-pbe', 'revdsd-pbepbe', 'revdsdpbe', 'revdsdpbepbe') ! no libxc
      num = p_revdsdpbe
   case('revdsd-blyp', 'revdsdblyp') ! no libxc
      num = p_revdsdblyp
   case('revdod-pbep86', 'revdodpbep86') ! no libxc
      num = p_revdodpbep86
   case('b97m', 'mgga_xc_b97m_v')
      num = p_b97m
   case('wb97m', 'ωb97m', 'omegab97m', 'hyb_mgga_xc_wb97m_v')
      num = p_wb97m
   case('wb97m-rev', 'ωb97m-rev', 'omegab97m-rev', 'wb97m_rev', 'ωb97m_rev', &
      & 'omegab97m_rev') ! D4 re-parametrization
      num = p_wb97m_rev
   case('wb97', 'ωb97', 'omegab97', 'hyb_gga_xc_wb97')
      num = p_wb97
   case('wb97x-2008', 'ωb97x-2008', 'omegab97x-2008', 'hyb_gga_xc_wb97x', &
      & 'wb97x_2008', 'ωb97x_2008', 'omegab97x_2008')
      num = p_wb97x_2008
   case('wb97x', 'ωb97x', 'omegab97x', 'hyb_gga_xc_wb97x_v')
      num = p_wb97x
   case('wb97x-rev', 'ωb97x-rev', 'omegab97x-rev', 'wb97x_rev', 'ωb97x_rev', &
      & 'omegab97x_rev') ! D4 re-parametrization
      num = p_wb97x_rev
   case('wb97x-3c', 'ωb97x-3c', 'omegab97x-3c', 'wb97x_3c', 'ωb97x_3c', &
      & 'omegab97x_3c') ! no libxc
      num = p_wb97x_3c
   case('wr2scan', 'wr²scan') ! no libxc
      num = p_wr2scan
   case('r2scan0-dh', 'r²scan0-dh', 'r2scan0dh', 'r²scan0dh') ! no libxc
      num = p_r2scan0_dh
   case('r2scan-cidh', 'r²scan-cidh', 'r2scancidh', 'r²scancidh') ! no libxc
      num = p_r2scan_cidh
   case('r2scan-qidh', 'r²scan-qidh', 'r2scanqidh', 'r²scanqidh') ! no libxc
      num = p_r2scan_qidh
   case('r2scan0-2', 'r²scan0-2', 'r2scan02', 'r²scan02') ! no libxc
      num = p_r2scan0_2
   case('pr2scan50', 'pr²scan50') ! no libxc
      num = p_pr2scan50
   case('pr2scan69', 'pr²scan69') ! no libxc
      num = p_pr2scan69
   case('kpr2scan50', 'kpr²scan50') ! no libxc
      num = p_kpr2scan50
   case('wpr2scan50', 'wpr²scan50') ! no libxc
      num = p_wpr2scan50
   end select
end function get_functional_id


end module dftd4_param_functionals
