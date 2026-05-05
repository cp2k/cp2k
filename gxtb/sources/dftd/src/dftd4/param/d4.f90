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

module dftd4_param_d4
   use dftd4_damping_type, only : twobody_damping_function, threebody_damping_function
   use dftd4_param_functionals
   use dftd4_param_type, only : param_type
   use dftd4_utils, only : lowercase
   use mctc_env, only : error_type, fatal_error, wp
   implicit none
   private

   public :: get_damping_d4

contains


!> Retrieve D4 damping parameters from functional ID
subroutine get_damping_d4(error, id, damping_2b, damping_3b, param)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Functional ID for which to retrieve the damping parameters
   integer, intent(in) :: id
   !> Type of two-body damping function to use
   integer, intent(in) :: damping_2b
   !> Type of three-body damping function to use
   integer, intent(in) :: damping_3b
   !> Damping parameters for the functional
   type(param_type), intent(inout) :: param

   logical :: found

   select case(damping_2b)
   case(twobody_damping_function%rational)
      ! Add optional three-body damping parameters
      select case(damping_3b)
      case(threebody_damping_function%rational)
         call fatal_error(error, "No D4 damping parameters available for "// &
            & "three-body rational damping combined with two-body rational damping.")
         return
      case(threebody_damping_function%screened)
         call fatal_error(error, "No D4 damping parameters available for "// &
            & "three-body screened damping combined with two-body rational damping.")
         return
      case(threebody_damping_function%zero)
         call fatal_error(error, "No D4 damping parameters available for "// &
            & "three-body zero damping combined with two-body rational damping.")
         return
      case(threebody_damping_function%zero_avg)
         found = get_d4eeq_bj_atm_0avg_parameter(id, param)
      case(: 0)
         ! No three-body ATM dispersion term
         found = get_d4eeq_bj_parameter(id, param)
         ! If no parameters are available without ATM, we only remove the three-body term
         if (.not. found) then
            found = get_d4eeq_bj_atm_0avg_parameter(id, param)
            deallocate(param%s9)
         end if
      case default
         call fatal_error(error, "Unknown D4 three-body damping function combined "// &
            & "with two-body rational damping.")
         return
      end select
   case(twobody_damping_function%screened)
      call fatal_error(error, "No D4 damping parameters available for two-body "// &
         & "screened damping.")
      return
   case(twobody_damping_function%zero)
      call fatal_error(error, "No D4 damping parameters available for two-body "// &
         & "zero damping.")
      return
   case(twobody_damping_function%mzero)
      call fatal_error(error, "No D4 damping parameters available for two-body "// &
         & "modified zero damping.")
      return
   case(twobody_damping_function%optpower)
      call fatal_error(error, "No D4 damping parameters available for two-body "// &
         & "optimized power damping.")
      return
   case(twobody_damping_function%cso)
      call fatal_error(error, "No D4 damping parameters available for two-body "// &
         & "C-Six-Only damping.")
      return
   case(twobody_damping_function%koide)
      call fatal_error(error, "No D4 damping parameters available for two-body "// &
         & "spherical wave Koide damping.")
      return
   case default
      call fatal_error(error, "Unknown two-body damping function for D4.")
      return
   end select

   if (.not. found) then
      call fatal_error(error, "No D4 damping parameters available for this functional.")
      return
   end if

end subroutine get_damping_d4


function get_d4eeq_bj_parameter(dfnum, param) result(found)
   !> Functional ID for which to retrieve the damping parameters
   integer(df_enum), intent(in) :: dfnum
   !> Damping parameters for the functional
   type(param_type), intent(inout) :: param
   !> Whether parameters were found for the functional
   logical :: found

   found = .true.
   select case(dfnum)
   case(p_default)
      ! No functional specified, use default parameters
      call dftd_param(param)
   case(p_dftb_3ob)
      call dftd_param(param, & ! (SAW191202)
         &  s6=1.0_wp, s8=0.4727337_wp, a1=0.5467502_wp, a2=4.4955068_wp)
   case(p_dftb_matsci)
      call dftd_param(param, & ! (SAW191202)
         &  s6=1.0_wp, s8=2.7711819_wp, a1=0.4681712_wp, a2=5.2918629_wp)
   case(p_dftb_mio)
      call dftd_param(param, & ! (SAW191202)
         &  s6=1.0_wp, s8=1.1948145_wp, a1=0.6074567_wp, a2=4.9336133_wp)
   case(p_dftb_ob2)
      call dftd_param(param, & ! (SAW191202)
         &  s6=1.0_wp, s8=2.7611320_wp, a1=0.6037249_wp, a2=5.3900004_wp)
   case(p_dftb_pbc)
      call dftd_param(param, & ! (SAW191202)
         &  s6=1.0_wp, s8=1.7303734_wp, a1=0.5546548_wp, a2=4.7973454_wp)
   case default
      found = .false.
   end select

contains

   pure subroutine dftd_param(par, s6, s8, a1, a2, alp)
      !> Damping parameters for the functional
      type(param_type), intent(inout) :: par
      !> Required parameters for D4 rational damping
      real(wp), intent(in), optional :: s8, a1, a2
      !> Optional parameters for D4 rational damping with global defaults
      real(wp), intent(in), optional :: s6, alp

      if (.not. allocated(par%s6)) then
         if (present(s6)) then
            par%s6 = s6
         else
            par%s6 = 1.0_wp
         end if
      end if

      ! No default values for s8, only assigned if present
      if (.not. allocated(par%s8)) then
         if (present(s8)) then
            par%s8 = s8
         end if
      end if

      ! No default values for a1, only assigned if present
      if (abs(par%a1) < epsilon(1.0_wp)) then
         if (present(a1)) then
            par%a1 = a1
         end if
      end if

      ! No default values for a2, only assigned if present
      if (abs(par%a2) < epsilon(1.0_wp)) then
         if (present(a2)) then
            par%a2 = a2
         end if
      end if

      if (.not. allocated(par%alp)) then
         if (present(alp)) then
            par%alp = alp
         else
            par%alp = 16.0_wp
         end if
      end if

   end subroutine dftd_param

end function get_d4eeq_bj_parameter

function get_d4eeq_bj_atm_0avg_parameter(dfnum, param) result(found)
   !> Functional ID for which to retrieve the damping parameters
   integer(df_enum), intent(in) :: dfnum
   !> Damping parameters for the functional
   type(param_type), intent(inout) :: param
   !> Whether parameters were found for the functional
   logical :: found

   found = .true.
   select case(dfnum)
   case(p_default)
      ! No functional specified, use default parameters
      call dftd_param(param)
   case(p_b1b95)
      call dftd_param(param, & ! (SAW190107)
         &  s6=1.0000_wp, s8=1.27701162_wp, a1=0.40554715_wp, a2=4.63323074_wp )
      !  Fitset: MD= 0.22852 MAD= 0.35189 RMSD= 0.46982
   case(p_b1lyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.98553711_wp, a1=0.39309040_wp, a2=4.55465145_wp )
      !  Fitset: MD= -0.04797 MAD= 0.25597 RMSD= 0.38778
   case(p_b1p)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.36115015_wp, a1=0.48665293_wp, a2=5.05219572_wp )
      !  Fitset: MD= -0.01406 MAD= 0.27441 RMSD= 0.47328
   case(p_b1pw)
      call dftd_param(param, & ! (SAW190107)
         &  s6=1.0000_wp, s8=3.02227550_wp, a1=0.47396846_wp, a2=4.49845309_wp )
      !  Fitset: MD= 0.10485 MAD= 0.32175 RMSD= 0.48508
   case(p_b2gpplyp)
      call dftd_param(param, & ! (SAW190107)
         &  s6=0.5600_wp, s8=0.94633372_wp, a1=0.42907301_wp, a2=5.18802602_wp )
      !  Fitset: MD= -0.05248 MAD= 0.18110 RMSD= 0.27365
   case(p_b2plyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.6400_wp, s8=1.16888646_wp, a1=0.44154604_wp, a2=4.73114642_wp )
      !  Fitset: MD= -0.03761 MAD= 0.18247 RMSD= 0.27109
   case(p_b3lyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.02929367_wp, a1=0.40868035_wp, a2=4.53807137_wp )
      !  Fitset: MD= -0.05892 MAD= 0.26117 RMSD= 0.40531
   case(p_b3p)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.08822155_wp, a1=0.47324238_wp, a2=4.98682134_wp )
      !  Fitset: MD= -0.02970 MAD= 0.26962 RMSD= 0.46761
   case(p_b3pw)
      call dftd_param(param, & ! (SAW190107)
         &  s6=1.0000_wp, s8=2.88364295_wp, a1=0.46990860_wp, a2=4.51641422_wp )
      !  Fitset: MD= 0.06643 MAD= 0.29151 RMSD= 0.45541
   case(p_b97)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.87854260_wp, a1=0.29319126_wp, a2=4.51647719_wp )
      !  Fitset: MD= -0.13017 MAD= 0.24778 RMSD= 0.36116
   case(p_bhlyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.65281646_wp, a1=0.27263660_wp, a2=5.48634586_wp )
      !  Fitset: MD= -0.15832 MAD= 0.34132 RMSD= 0.57342
   case(p_blyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.34076671_wp, a1=0.44488865_wp, a2=4.09330090_wp )
      !  Fitset: MD= 0.04801 MAD= 0.28161 RMSD= 0.38321
   case(p_bpbe)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.64405246_wp, a1=0.52905620_wp, a2=4.11311891_wp )
      !  Fitset: MD= 0.19316 MAD= 0.41912 RMSD= 0.60452
   case(p_bp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.35497927_wp, a1=0.43645861_wp, a2=4.92406854_wp )
      !  Fitset: MD= 0.08252 MAD= 0.32681 RMSD= 0.47063
   case(p_bpw)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.24571506_wp, a1=0.50050454_wp, a2=4.12346483_wp )
      !  Fitset: MD= 0.20607 MAD= 0.41941 RMSD= 0.59589
   case(p_camb3lyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.66041301_wp, a1=0.40267156_wp, a2=5.17432195_wp )
      !  Fitset: MD= -0.19675 MAD= 0.34901 RMSD= 0.59087
   case(p_camqtp01)
      call dftd_param(param, & ! (10.1021/acs.jctc.3c00717)
         &  s6=1.0000_wp, s8=1.156_wp, a1=0.461_wp, a2=6.375_wp )
   case(p_dodblyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.4700_wp, s8=1.31146043_wp, a1=0.43407294_wp, a2=4.27914360_wp )
      !  Fitset: MD= 0.03323 MAD= 0.13858 RMSD= 0.20861
   case(p_dodpbeb95)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.5600_wp, s8=0.01574635_wp, a1=0.43745720_wp, a2=3.69180763_wp )
      !  Fitset: MD= 0.03704 MAD= 0.13343 RMSD= 0.18278
   case(p_dodpbe)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.4800_wp, s8=0.92051454_wp, a1=0.43037052_wp, a2=4.38067238_wp )
      !  Fitset: MD= 0.01065 MAD= 0.13414 RMSD= 0.21424
   case(p_dodpbep86)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.4600_wp, s8=0.71405681_wp, a1=0.42408665_wp, a2=4.52884439_wp )
      !  Fitset: MD= -0.03740 MAD= 0.12467 RMSD= 0.18127
   case(p_dodsvwn)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.4200_wp, s8=0.94500207_wp, a1=0.47449026_wp, a2=5.05316093_wp )
      !  Fitset: MD= -0.07427 MAD= 0.16970 RMSD= 0.25286
   case(p_dsdblyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.5400_wp, s8=0.63018237_wp, a1=0.47591835_wp, a2=4.73713781_wp )
      !  Fitset: MD= -0.01981 MAD= 0.14823 RMSD= 0.21530
   case(p_dsdpbeb95)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.5400_wp, s8=-0.14668670_wp, a1=0.46394587_wp, a2=3.64913860_wp )
      !  Fitset: MD= 0.02996 MAD= 0.12414 RMSD= 0.16860
   case(p_dsdpbe)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.4500_wp, s8=0.70584116_wp, a1=0.45787085_wp, a2=4.44566742_wp )
      !  Fitset: MD= 0.00866 MAD= 0.13406 RMSD= 0.21380
   case(p_dsdpbep86)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.4700_wp, s8=0.37586675_wp, a1=0.53698768_wp, a2=5.13022435_wp )
      !  Fitset: MD= -0.05273 MAD= 0.14259 RMSD= 0.21271
   case(p_dsdsvwn)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.4100_wp, s8=0.72914436_wp, a1=0.51347412_wp, a2=5.11858541_wp )
      !  Fitset: MD= -0.08974 MAD= 0.32285 RMSD= 0.43146
   case(p_glyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=4.23798924_wp, a1=0.38426465_wp, a2=4.38412863_wp )
      !  Fitset: MD= 0.63466 MAD= 0.89568 RMSD= 1.11309
   case(p_hf)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.61679827_wp, a1=0.44959224_wp, a2=3.35743605_wp )
      !  Fitset: MD= -0.02597 MAD= 0.34732 RMSD= 0.49719
   case(p_lb94)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.59538499_wp, a1=0.42088944_wp, a2=3.28193223_wp )
      !  Fitset: MD= 0.31701 MAD= 0.53196 RMSD= 0.74553
   case(p_lcblyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.60344180_wp, a1=0.45769839_wp, a2=7.86924893_wp )
      !  Fitset: MD= -0.39724 MAD= 0.72327 RMSD= 1.18218
   case(p_lcwpbe)
      call dftd_param(param, & ! (10.1021/acs.jctc.3c00717)
         &  s6=1.0000_wp, s8=1.170_wp, a1=0.378_wp, a2=4.816_wp )
   case(p_lcwpbeh)
      call dftd_param(param, & ! (10.1021/acs.jctc.3c00717)
         &  s6=1.0000_wp, s8=1.318_wp, a1=0.386_wp, a2=5.010_wp )
   case(p_lh07ssvwn)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.16675531_wp, a1=0.35965552_wp, a2=4.31947614_wp )
      !  Fitset: MD= 0.32224 MAD= 0.59006 RMSD= 0.86272
   case(p_lh07tsvwn)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.09333001_wp, a1=0.35025189_wp, a2=4.34166515_wp )
      !  Fitset: MD= 0.24243 MAD= 0.43497 RMSD= 0.61671
   case(p_lh12ctssifpw92)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.68467610_wp, a1=0.34190416_wp, a2=3.91039666_wp )
      !  Fitset: MD= 0.55106 MAD= 0.80783 RMSD= 1.11048
   case(p_lh12ctssirpw92)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.48973402_wp, a1=0.34026075_wp, a2=3.96948081_wp )
      !  Fitset: MD= 0.47785 MAD= 0.71188 RMSD= 0.98422
   case(p_lh14tcalpbe)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.28130770_wp, a1=0.38822021_wp, a2=4.92501211_wp )
      !  Fitset: MD= -0.02105 MAD= 0.22968 RMSD= 0.36045
   case(p_lh20t)
      call dftd_param(param, & ! (10.1021/acs.jctc.0c00498)
         & s6=1.000_wp, s8=0.113_wp, a1=0.479_wp, a2=4.635_wp )
   case(p_m06)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.16366729_wp, a1=0.53456413_wp, a2=6.06192174_wp )
      !  Fitset: MD= 0.01788 MAD= 0.24914 RMSD= 0.38604
   case(p_m06l)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.59493760_wp, a1=0.71422359_wp, a2=6.35314182_wp )
      !  Fitset: MD= 0.08395 MAD= 0.24888 RMSD= 0.34879
   case(p_mn12sx)
      call dftd_param(param, & ! (SAW211021)
         &  s6=1.0000_wp, s8=0.85964873_wp, a1=0.62662681_wp, a2=5.62088906_wp )
      !  Fitset: MD= 0.16131 MAD= 0.34142 RMSD= 0.47113
   case(p_mpw1b95)
      call dftd_param(param, & ! (SAW190107)
         &  s6=1.0000_wp, s8=0.50093024_wp, a1=0.41585097_wp, a2=4.99154869_wp )
      !  Fitset: MD= 0.00585 MAD= 0.15695 RMSD= 0.21297
   case(p_mpw1lyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.15591153_wp, a1=0.25603493_wp, a2=5.32083895_wp )
      !  Fitset: MD= -0.26979 MAD= 0.41542 RMSD= 0.60678
   case(p_mpw1pw)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.80841716_wp, a1=0.42961819_wp, a2=4.68892341_wp )
      !  Fitset: MD= -0.08840 MAD= 0.26815 RMSD= 0.45231
   case(p_mpw2plyp)
      call dftd_param(param, & ! (SAW190107)
         &  s6=0.7500_wp, s8=0.45788846_wp, a1=0.42997704_wp, a2=5.07650682_wp )
      !  Fitset: MD= -0.18921 MAD= 0.30115 RMSD= 0.44049
   case(p_mpwb1k)
      call dftd_param(param, & ! (SAW190107)
         &  s6=1.0000_wp, s8=0.57338313_wp, a1=0.44687975_wp, a2=5.21266777_wp )
      !  Fitset: MD= -0.00870 MAD= 0.17226 RMSD= 0.23614
   case(p_mpwlyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.25842942_wp, a1=0.25773894_wp, a2=5.02319542_wp )
      !  Fitset: MD= -0.24426 MAD= 0.39145 RMSD= 0.54503
   case(p_mpwpw)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.82596836_wp, a1=0.34526745_wp, a2=4.84620734_wp )
      !  Fitset: MD= -0.06278 MAD= 0.27913 RMSD= 0.43988
   case(p_o3lyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.75762508_wp, a1=0.10348980_wp, a2=6.16233282_wp )
      !  Fitset: MD= -0.19268 MAD= 0.38577 RMSD= 0.62168
   case(p_olyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=2.74836820_wp, a1=0.60184498_wp, a2=2.53292167_wp )
      !  Fitset: MD= 0.12352 MAD= 0.37113 RMSD= 0.58291
   case(p_opbe)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=3.06917417_wp, a1=0.68267534_wp, a2=2.22849018_wp )
      !  Fitset: MD= 0.26699 MAD= 0.55308 RMSD= 0.85023
   case(p_pbe0_2)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.5000_wp, s8=0.64299082_wp, a1=0.76542115_wp, a2=5.78578675_wp )
      !  Fitset: MD= -0.04260 MAD= 0.21186 RMSD= 0.34045
   case(p_pbe0)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.20065498_wp, a1=0.40085597_wp, a2=5.02928789_wp )
      !  Fitset: MD= -0.17892 MAD= 0.30557 RMSD= 0.51050
   case(p_pbe0_dh)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.8750_wp, s8=0.96811578_wp, a1=0.47592488_wp, a2=5.08622873_wp )
      !  Fitset: MD= -0.13857 MAD= 0.27919 RMSD= 0.47256
   case(p_pbe)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.95948085_wp, a1=0.38574991_wp, a2=4.80688534_wp )
      !  Fitset: MD= -0.20544 MAD= 0.33635 RMSD= 0.51168
   case(p_pbesol)
      call dftd_param(param, & ! (SAW211021)
         &  s6=1.0000_wp, s8=1.71885698_wp, a1=0.47901421_wp, a2=5.96771589_wp )
      !  Fitset: MD= -0.28899 MAD= 0.52215 RMSD= 0.93584
   case(p_am05)
      call dftd_param(param, & ! (SAW211021)
         &  s6=1.0000_wp, s8=1.71885838_wp, a1=0.47901431_wp, a2=5.96771581_wp )
      !  Fitset: MD= -0.28899 MAD= 0.52215 RMSD= 0.93584
   case(p_pw1pw)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.96850170_wp, a1=0.42427511_wp, a2=5.02060636_wp )
      !  Fitset: MD= -0.27325 MAD= 0.42206 RMSD= 0.64119
   case(p_pw6b95)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=-0.31926054_wp, a1=0.04142919_wp, a2=5.84655608_wp )
      !  Fitset: MD= -0.04767 MAD= 0.14330 RMSD= 0.18958
   case(p_pw86pbe)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.21362856_wp, a1=0.40510366_wp, a2=4.66737724_wp )
      !  Fitset: MD= -0.11505 MAD= 0.24691 RMSD= 0.38101
   case(p_pw91)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.77283111_wp, a1=0.39581542_wp, a2=4.93405761_wp )
      !  Fitset: MD= -0.33019 MAD= 0.48611 RMSD= 0.68110
   case(p_pwp1)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.60492565_wp, a1=0.46855837_wp, a2=5.76921413_wp )
      !  Fitset: MD= -0.35321 MAD= 0.54026 RMSD= 0.86629
   case(p_pwpb95)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.8200_wp, s8=-0.34639127_wp, a1=0.41080636_wp, a2=3.83878274_wp )
      !  Fitset: MD= 0.02143 MAD= 0.13040 RMSD= 0.17599
   case(p_pwp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=0.32801227_wp, a1=0.35874687_wp, a2=6.05861168_wp )
      !  Fitset: MD= -0.42482 MAD= 0.62607 RMSD= 0.91840
   case(p_revpbe0)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.57185414_wp, a1=0.38705966_wp, a2=4.11028876_wp )
      !  Fitset: MD= 0.02724 MAD= 0.21587 RMSD= 0.36040
   case(p_revpbe0dh)
      call dftd_param(param, & ! (SAW190103)
         &  s6=0.8750_wp, s8=1.24456037_wp, a1=0.36730560_wp, a2=4.71126482_wp )
      !  Fitset: MD= -0.01089 MAD= 0.20910 RMSD= 0.33564
   case(p_revpbe38)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.66597472_wp, a1=0.39476833_wp, a2=4.39026628_wp )
      !  Fitset: MD= -0.01326 MAD= 0.22598 RMSD= 0.36210
   case(p_revpbe)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.74676530_wp, a1=0.53634900_wp, a2=3.07261485_wp )
      !  Fitset: MD= 0.05649 MAD= 0.25212 RMSD= 0.40863
   case(p_revtpss0)
      call dftd_param(param, & ! (SAW190107)
         &  s6=1.0000_wp, s8=1.54664499_wp, a1=0.45890964_wp, a2=4.78426405_wp )
      !  Fitset: MD= -0.05298 MAD= 0.19965 RMSD= 0.32081
   case(p_revtpss)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.53089454_wp, a1=0.44880597_wp, a2=4.64042317_wp )
      !  Fitset: MD= -0.01904 MAD= 0.19568 RMSD= 0.29618
   case(p_revtpssh)
      call dftd_param(param, & ! (SAW190107)
         &  s6=1.0000_wp, s8=1.52740307_wp, a1=0.45161957_wp, a2=4.70779483_wp )
      !  Fitset: MD= -0.03731 MAD= 0.19133 RMSD= 0.29091
   case(p_rpbe)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.31183787_wp, a1=0.46169493_wp, a2=3.15711757_wp )
      !  Fitset: MD= -0.07156 MAD= 0.26348 RMSD= 0.38671
   case(p_rpw86pbe)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.12624034_wp, a1=0.38151218_wp, a2=4.75480472_wp )
      !  Fitset: MD= -0.12740 MAD= 0.26294 RMSD= 0.40614
   case(p_scan)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.46126056_wp, a1=0.62930855_wp, a2=6.31284039_wp )
      !  Fitset: MD= -0.13170 MAD= 0.28640 RMSD= 0.51183
   case(p_rscan)
      call dftd_param(param, & ! (10.1063/5.0041008)
         &  s6=1.0000_wp, s8=0.87728975_wp, a1=0.49116966_wp, a2=5.75859346_wp )
   case(p_r2scan)
      call dftd_param(param, & ! (10.1063/5.0041008)
         &  s6=1.0000_wp, s8=0.60187490_wp, a1=0.51559235_wp, a2=5.77342911_wp )
   case(p_r2scanh)
      call dftd_param(param, & ! (10.1063/5.0086040)
         & s6=1.0_wp, s8=0.8324_wp, a1=0.4944_wp, a2=5.9019_wp)
   case(p_r2scan0)
      call dftd_param(param, & ! (10.1063/5.0086040)
         & s6=1.0_wp, s8=0.8992_wp, a1=0.4778_wp, a2=5.8779_wp)
   case(p_r2scan50)
      call dftd_param(param, & ! (10.1063/5.0086040)
         & s6=1.0_wp, s8=1.0471_wp, a1=0.4574_wp, a2=5.8969_wp)
   case(p_r2scan_3c)
      call dftd_param(param, & ! (10.1063/5.0040021)
         & s6=1.0_wp, s8=0.00_wp, a1=0.42_wp, a2=5.65_wp, s9=2.00_wp)
   case(p_tpss0)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.62438102_wp, a1=0.40329022_wp, a2=4.80537871_wp )
      !  Fitset: MD= -0.09569 MAD= 0.26733 RMSD= 0.44767
   case(p_tpss)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.76596355_wp, a1=0.42822303_wp, a2=4.54257102_wp )
      !  Fitset: MD= -0.09296 MAD= 0.27505 RMSD= 0.42537
   case(p_tpssh)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.85897750_wp, a1=0.44286966_wp, a2=4.60230534_wp )
      !  Fitset: MD=  0.02238 MAD= 0.16042 RMSD= 0.33519
   case(p_b97d)
      call dftd_param(param, & ! (SAW201029)
         &  s6=1.0000_wp, s8=1.69460052_wp, a1=0.28904684_wp, a2=4.13407323_wp )
      !  Fitset: MD= -0.09858 MAD= 0.26757 RMSD= 0.42380
   case(p_wb97)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=6.55792598_wp, a1=0.76666802_wp, a2=8.36027334_wp )
      !  Fitset: MD= -0.12779 MAD= 0.36152 RMSD= 0.49991
   case(p_wb97x_2008)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=-0.07519516_wp, a1=0.45094893_wp, a2=6.78425255_wp )
      !  S22x5: MD= 0.05 MAD= 0.16 RMSD= 0.22
      !  S66x8: MD= 0.06 MAD= 0.16 RMSD= 0.21
      !  NCI10: MD= 0.08 MAD= 0.15 RMSD= 0.25
   case(p_wb97x)
      call dftd_param(param, & ! (10.1002/jcc.26411)
         &  s6=1.0000_wp, s8=0.5093_wp, a1=0.0662_wp, a2=5.4487_wp )
   case(p_wb97x_rev)
      call dftd_param(param, & ! (10.1063/5.0133026)
         &  s6=1.0000_wp, s8=0.4485_wp, a1=0.3306_wp, a2=4.279_wp )
   case(p_wb97x_3c)
      call dftd_param(param, & ! (10.1063/5.0133026)
         &  s6=1.0000_wp, s8=0.0_wp, a1=0.2464_wp, a2=4.737_wp )
   case(p_b97m)
      call dftd_param(param, & ! (10.1002/jcc.26411)
         &  s6=1.0000_wp, s8=0.6633_wp, a1=0.4288_wp, a2=3.9935_wp )
      !  S22x5: MD= 0.03 MAD= 0.12 RMSD= 0.18
      !  S66x8: MD= 0.09 MAD= 0.17 RMSD= 0.22
      !  NCI10: MD= 0.09 MAD= 0.15 RMSD= 0.32
   case(p_wb97m)
      call dftd_param(param, & ! (10.1002/jcc.26411)
         &  s6=1.0000_wp, s8=0.7761_wp, a1=0.7514_wp, a2=2.7099_wp )
      !  Fitset: MD= -0.20216 MAD= 0.34696 RMSD= 0.53641
   case(p_wb97m_rev) 
      call dftd_param(param, & ! (10.1021/acs.jctc.3c00717)
         &  s6=1.0000_wp, s8=0.842_wp, a1=0.359_wp, a2=4.668_wp )
   case(p_x3lyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.54701429_wp, a1=0.20318443_wp, a2=5.61852648_wp )
      !  Fitset: MD= -0.15607 MAD= 0.31342 RMSD= 0.49546
   case(p_xlyp)
      call dftd_param(param, & ! (SAW190103)
         &  s6=1.0000_wp, s8=1.62972054_wp, a1=0.11268673_wp, a2=5.40786417_wp )
      !  Fitset: MD= -0.03900 MAD= 0.27562 RMSD= 0.38491
   case(p_revdsdpbep86)
      call dftd_param(param, & ! (WTMAD2)
         &  s6=0.5132_wp, s8=0.00000000_wp, a1=0.44000000_wp, a2=3.60000000_wp )
   case(p_revdsdpbe)
      call dftd_param(param, & ! (WTMAD2)
         &  s6=0.6706_wp, s8=0.00000000_wp, a1=0.40000000_wp, a2=3.60000000_wp )
   case(p_revdsdblyp)
      call dftd_param(param, & !(WTMAD2)
         &  s6=0.6141_wp, s8=0.00000000_wp, a1=0.38000000_wp, a2=3.52000000_wp )
   case(p_revdodpbep86)
      call dftd_param(param, & !(WTMAD2)
         &  s6=0.5552_wp, s8=0.00000000_wp, a1=0.44000000_wp, a2=3.60000000_wp )
   case(p_dftb_3ob)
      call dftd_param(param, & ! (SAW191202)
         &  s6=1.0_wp, s8=0.6635015_wp, a1=0.5523240_wp, a2=4.3537076_wp)
   case(p_dftb_matsci)
      call dftd_param(param, & ! (SAW191202)
         &  s6=1.0_wp, s8=3.3157614_wp, a1=0.4826330_wp, a2=5.3811976_wp)
   case(p_dftb_mio)
      call dftd_param(param, & ! (SAW191202)
         &  s6=1.0_wp, s8=1.2916225_wp, a1=0.5965326_wp, a2=4.8778602_wp)
   case(p_dftb_ob2)
      call dftd_param(param, & ! (SAW191202)
         &  s6=1.0_wp, s8=2.9692689_wp, a1=0.6068916_wp, a2=5.4476789_wp)
   case(p_dftb_pbc)
      call dftd_param(param, & ! (SAW191202)
         &  s6=1.0_wp, s8=2.1667394_wp, a1=0.5646391_wp, a2=4.9576353_wp)
   case(p_hse03)
      call dftd_param(param, & ! (SAW211107)
         &  s6=1.0_wp, s8=1.19812280_wp, a1=0.38662939_wp, a2=5.22925796_wp)
   case(p_hse06)
      call dftd_param(param, & ! (SAW211107)
         &  s6=1.0_wp, s8=1.19528249_wp, a1=0.38663183_wp, a2=5.19133469_wp)
   case(p_hse12)
      call dftd_param(param, & ! (SAW211107)
         &  s6=1.0_wp, s8=1.23500792_wp, a1=0.39226921_wp, a2=5.22036266_wp)
   case(p_hse12s)
      call dftd_param(param, & ! (SAW211107)
         &  s6=1.0_wp, s8=1.23767762_wp, a1=0.39989137_wp, a2=5.34809245_wp)
   case(p_hsesol)
      call dftd_param(param, & ! (SAW211107)
         &  s6=1.0_wp, s8=1.82207807_wp, a1=0.45646268_wp, a2=5.59662251_wp)
   case(p_wr2scan) ! (10.1063/5.0174988)
      call dftd_param(param, &
         & s6=1.0_wp, s8=1.0_wp, a1=0.3834_wp, a2=5.7889_wp)
   case(p_r2scan0_dh) ! (10.1063/5.0174988)
      call dftd_param(param, &
         & s6=0.9424_wp, s8=0.3856_wp, a1=0.4271_wp, a2=5.8565_wp)
   case(p_r2scan_cidh) ! (10.1063/5.0174988)
      call dftd_param(param, &
         & s6=0.8666_wp, s8=0.5336_wp, a1=0.4171_wp, a2=5.9125_wp)
   case(p_r2scan_qidh) ! (10.1063/5.0174988)
      call dftd_param(param, &
         & s6=0.7867_wp, s8=0.2955_wp, a1=0.4001_wp, a2=5.8300_wp)
   case(p_r2scan0_2) ! (10.1063/5.0174988)
      call dftd_param(param, &
         & s6=0.7386_wp, s8=0.0000_wp, a1=0.4030_wp, a2=5.5142_wp)
   case(p_pr2scan50) ! (10.1063/5.0174988)
      call dftd_param(param, &
         & s6=0.7964_wp, s8=0.3421_wp, a1=0.4663_wp, a2=5.7916_wp)
   case(p_pr2scan69) ! (10.1063/5.0174988)
      call dftd_param(param, &
         & s6=0.7167_wp, s8=0.0000_wp, a1=0.4644_wp, a2=5.2563_wp)
   case(p_kpr2scan50) ! (10.1063/5.0174988)
      call dftd_param(param, &
         & s6=0.8402_wp, s8=0.1212_wp, a1=0.4382_wp, a2=5.8232_wp)
   case(p_wpr2scan50) ! (10.1063/5.0174988)
      call dftd_param(param, &
         & s6=0.8143_wp, s8=0.3842_wp, a1=0.4135_wp, a2=5.8773_wp)
   case default
      found = .false.
   end select

contains

   pure subroutine dftd_param(par, s6, s8, a1, a2, alp, s9, rs9)
      !> Damping parameters for the functional
      type(param_type), intent(inout) :: par
      !> Required parameters for D4 rational damping
      real(wp), intent(in), optional :: s8, a1, a2
      !> Optional parameters for D4 rational damping with global defaults
      real(wp), intent(in), optional :: s6, alp, s9, rs9

      if (.not. allocated(par%s6)) then
         if (present(s6)) then
            par%s6 = s6
         else
            par%s6 = 1.0_wp
         end if
      end if

      ! No default values for s8, only assigned if present
      if (.not. allocated(par%s8)) then
         if (present(s8)) then
            par%s8 = s8
         end if
      end if

      ! No default values for a1, only assigned if present
      if (abs(par%a1-1.0_wp) < epsilon(1.0_wp)) then
         if (present(a1)) then
            par%a1 = a1
         end if
      end if

      ! No default values for a2, only assigned if present
      if (abs(par%a2) < epsilon(1.0_wp)) then
         if (present(a2)) then
            par%a2 = a2
         end if
      end if

      if (.not. allocated(par%alp)) then
         if (present(alp)) then
            par%alp = alp
         else
            par%alp = 16.0_wp
         end if
      end if

      if (.not. allocated(par%s9)) then
         if (present(s9)) then
            par%s9 = s9
         else
            par%s9 = 1.0_wp
         end if
      end if

      if (.not. allocated(param%rs9)) then
         if (present(rs9)) then
            par%rs9 = rs9
         else
            par%rs9 = 1.0_wp
         end if
      end if

   end subroutine dftd_param

end function get_d4eeq_bj_atm_0avg_parameter

end module dftd4_param_d4
