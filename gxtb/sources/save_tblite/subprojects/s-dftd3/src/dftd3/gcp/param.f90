! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

module dftd3_gcp_param
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use dftd3_data, only : get_vdw_rad
   implicit none
   private

   public :: gcp_param, get_gcp_param

   !> Parameters for the geometric counter-poise correction
   type :: gcp_param
      !> Basis set superposition error correction
      real(wp), allocatable :: emiss(:)
      !> Number of virtual orbitals
      real(wp), allocatable :: xv(:)
      !> Slater exponents
      real(wp), allocatable :: slater(:)
      !> Van der Waals radii for effective nuclear charges
      real(wp), allocatable :: rvdw(:, :)
      !> Van der Waals radii for true nuclear charges
      real(wp), allocatable :: rvdw_srb(:, :)
      !> Effective nuclear charges
      integer, allocatable :: zeff(:)
      !> Scaling factor for the counter-poise correction
      real(wp) :: sigma = 0.0_wp
      !> Exponential parameter
      real(wp) :: alpha = 0.0_wp
      !> Power parameter
      real(wp) :: beta = 0.0_wp
      !> Damping enabled
      logical :: damp = .false.
      !> Damping scaling factor
      real(wp) :: dmp_scal = 4.0_wp
      !> Damping exponent
      real(wp) :: dmp_exp = 6.0_wp
      !> Short-range bond correction
      logical :: srb = .false.
      !> Short-range bond correction radii scaling factor
      real(wp) :: rscal = 0.0_wp
      !> Short-range bond correction scaling factor
      real(wp) :: qscal = 0.0_wp
      !> Short-range bond correction for HF-3c
      logical :: base = .false.
    end type

   enum, bind(c)
      enumerator :: &
         p_unknown_bas, &
         p_sv_bas, &
         p_sv_p_bas, &
         p_svx_bas, &
         p_svp_bas, &
         p_svp_old_bas, &
         p_minis_bas, &
         p_631gd_bas, &
         p_tz_bas, &
         p_deftzvp_bas, &
         p_ccdz_bas, &
         p_accdz_bas, &
         p_pobtz_bas, &
         ! p_pobdzvp_bas, &
         p_minix_bas, &
         p_gcore_bas, &
         p_2g_bas, &
         p_dzp_bas, &
         ! p_hsv_bas, &
         p_dz_bas, &
         p_msvp_bas, &
         p_lanl2_bas, &
         p_pbeh3c_bas, &
         p_def2mtzvpp_bas, &
         p_def2mtzvp_bas, &
         p_vmb_bas
   end enum

   enum, bind(c)
      enumerator :: &
         p_unknown_method, &
         p_hf_method, &
         p_dft_method, &
         p_hyb_method, &
         p_gga_method, &
         p_b3lyp_method, &
         p_blyp_method, &
         p_pbe_method, &
         p_tpss_method, &
         p_pw6b95_method, &
         p_hf3c_method, &
         p_pbeh3c_method, &
         p_hse3c_method, &
         p_b973c_method, &
         p_b3pbe3c_method, &
         p_r2scan3c_method
   end enum

   real(wp), parameter :: emiss_hf_sv(*) = [ &
      & 0.009037_wp, 0.008843_wp, &  ! He,He
      & 0.204189_wp, 0.107747_wp, 0.049530_wp, 0.055482_wp, 0.072823_wp, 0.100847_wp, 0.134029_wp, 0.174222_wp, &  ! Li-Ne
      & 0.315616_wp, 0.261123_wp, 0.168568_wp, 0.152287_wp, 0.146909_wp, 0.168248_wp, 0.187882_wp, 0.211160_wp, &  !Na -Ar
      & 0.374252_wp, 0.460972_wp, &  ! K-Ca
      & 0.444886_wp, 0.404993_wp, 0.378406_wp, 0.373439_wp, 0.361245_wp, &
      & 0.360014_wp, 0.362928_wp, 0.243801_wp, 0.405299_wp, 0.396510_wp, &   ! 3d-TM
      & 0.362671_wp, 0.360457_wp, 0.363355_wp, 0.384170_wp, 0.399698_wp, 0.417307_wp] !Ga-Kr

   real(wp), parameter :: emiss_hf_minis(*) = [ &
      & 0.042400_wp, 0.028324_wp, &
      & 0.252661_wp, 0.197201_wp, 0.224237_wp, 0.279950_wp, 0.357906_wp, 0.479012_wp, 0.638518_wp, 0.832349_wp, &
      & 1.232920_wp, 1.343390_wp, 1.448280_wp, 1.613360_wp, 1.768140_wp, 1.992010_wp, 2.233110_wp, 2.493230_wp, &
      & 3.029640_wp, 3.389980_wp, &  ! H-Ca
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

   real(wp), parameter :: emiss_hf_631gd(*) = [ &! H-Ca + Br (no 3d)
      & 0.010083_wp, 0.008147_wp, &
      & 0.069260_wp, 0.030540_wp, 0.032736_wp, 0.021407_wp, 0.024248_wp, 0.036746_wp, 0.052733_wp, 0.075120_wp, &
      & 0.104255_wp, 0.070691_wp, 0.100260_wp, 0.072534_wp, 0.054099_wp, 0.056408_wp, 0.056025_wp, 0.057578_wp, &
      & 0.079198_wp, 0.161462_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.000000_wp, 0.381049_wp, 0.000000_wp]

   real(wp), parameter :: emiss_hf_svp_old(*) = [real(wp):: & ! Li,Na,Mg,K had wrong parameters
      & 0.008107,0.008045,&
      & 0.142957,0.028371,0.049369,0.055376,0.072785,0.100310,0.133273,0.173600,&
      & 0.191109,0.222839,0.167188,0.149843,0.145396,0.164308,0.182990,0.205668,&
      & 0.221189,0.299661,&
      & 0.325995,0.305488,0.291723,0.293801,0.29179,0.296729,0.304603,0.242041,0.354186,0.350715,&
      & 0.350021,0.345779,0.349532,0.367305,0.382008,0.399709]

   real(wp), parameter :: emiss_hf_svp(*) = [  & ! H-Kr
      & 0.008107,0.008045,&
      & 0.113583,0.028371,0.049369,0.055376,0.072785,0.100310,0.133273,0.173600,&
      & 0.181140,0.125558,0.167188,0.149843,0.145396,0.164308,0.182990,0.205668,&
      & 0.200956,0.299661, &
      & 0.325995,0.305488,0.291723,0.293801,0.29179,0.296729,0.304603,0.242041,0.354186,0.350715,&
      & 0.350021,0.345779,0.349532,0.367305,0.382008,0.399709]

   real(wp), parameter :: emiss_hf_sv_p(*) = [  & ! H-Kr
      & emiss_hf_sv(1), emiss_hf_svp(2:)]

   real(wp), parameter :: emiss_hf_svx(*) = [  & ! H-Kr
      & emiss_hf_sv(1), emiss_hf_svp(2:5), emiss_hf_sv(6), emiss_hf_svp(7:)]

   real(wp), parameter :: emiss_hf_tz(*) = [ &  ! H-Kr !def2-TZVP
      & 0.007577,0.003312,&
      & 0.086763,0.009962,0.013964,0.005997,0.004731,0.005687,0.006367,0.007511,&
      & 0.077721,0.050003,0.068317,0.041830,0.025796,0.025512,0.023345,0.022734,&
      & 0.097241,0.099167,&
      & 0.219194,0.189098,0.164378,0.147238,0.137298,0.12751,0.118589,0.0318653,0.120985,0.0568313, &
      & 0.090996,0.071820,0.063562,0.064241,0.061848,0.061021]
   
   real(wp), parameter :: emiss_hf_def2mtzvp(*) = [&  ! m def2-TZVP, no f for B-Ne
      & 0.007930,0.003310,&
      & 0.086760,0.009960,0.013960,0.006000,0.003760,0.004430,0.005380,0.006750,&
      & 0.077720,0.050000,0.068320,0.041830,0.025800,0.025510,0.023340,0.022730,&
      & 0.097240,0.099170,&
      & 0.219190,0.189100,0.164380,0.147240,0.137300,0.127510,0.118590,0.031870,0.120990,0.056830,&
      & 0.091000,0.071820,0.063560,0.064240,0.061850,0.061020]
   
   real(wp), parameter :: emiss_hf_vmb(*) = [ &
      & 0.042400_wp, 0.028324_wp, &
      & 0.252661_wp, 0.197201_wp, 0.156009_wp, 0.164586_wp, 0.169273_wp, 0.214704_wp, 0.729138_wp, 0.336072_wp, &
      & 0.262329_wp, 0.231722_wp, 0.251169_wp, 0.287795_wp, 0.213590_wp, 0.250524_wp, 0.728579_wp, 0.260658_wp, &
      & 0.0_wp, 0.0_wp,&
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
   
   real(wp), parameter :: emiss_hf_minisd(*) = [& !Al-Ar MINIS + Ahlrichs "P" funktions
      & 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 1.446950_wp, 1.610980_wp, 1.766610_wp, 1.988230_wp, 2.228450_wp, 2.487960_wp, &
      & 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

   real(wp), parameter :: emiss_hf_minix(*) = [ &
      & emiss_hf_minis(1:2), 0.177871_wp, 0.171596_wp, emiss_hf_minis(5:10), 1.114110_wp, 1.271150_wp, &
      & emiss_hf_minisd(13:18), emiss_hf_sv(19:30), emiss_hf_svp(31:36)]
   
   real(wp), parameter :: emiss_hf_lanl2(*) = [ & !  LANL2TZ+ vs LANL2DZ (ORCA), only Sc-Zn
      & emiss_hf_631gd(1:20), &
      & 0.102545_wp, 0.0719529_wp, 0.0491798_wp, 0.0362976_wp, 0.0266369_wp, &
      & 0.0235484_wp, 0.0171578_wp, 0.0438906_wp, 0.0100259_wp, 0.016208_wp, &
      & emiss_hf_631gd(31:)]
   
   real(wp), parameter :: emiss_hf_pobtz(*) = [ & ! H-Kr, no RG
      & 0.010077_wp, 0.000000_wp, &
      & 0.173239_wp, 0.101973_wp, 0.131181_wp, 0.032234_wp, 0.011630_wp, 0.008475_wp, 0.011673_wp, 0.000000_wp, &
      & 0.240653_wp, 0.661819_wp, 0.522306_wp, 0.141630_wp, 0.052456_wp, 0.184547_wp, 0.040837_wp, 0.000000_wp, &
      & 0.318136_wp, 0.564721_wp, &
      & 0.523194_wp, 0.767449_wp, 0.620122_wp, 0.390227_wp, 0.237571_wp, &
      & 0.247940_wp, 0.249589_wp, 0.117864_wp, 0.325725_wp, 0.197183_wp, &
      & 0.264489_wp, 0.180375_wp, 0.111262_wp, 0.147239_wp, 0.081747_wp, 0.000000_wp]

   real(wp), parameter :: emiss_hf_pobdzvp(*) = [ &  ! FIXME: https://github.com/grimme-lab/gcp/issues/22
      & 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
   
   real(wp), parameter :: emiss_hf_2gcore(*) = [real(wp)::  & ! only HCNOF yet
      & 0.000539,0.000000,&
      & 0.000000,0.000000,0.000000,0.173663,0.269952,0.364341,0.384923,0.000000,&
      & 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,&
      & 0.000000,0.000000,&
      & 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,&
      & 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000]
   
   real(wp), parameter :: emiss_hf_def1tzvp(*) = [real(wp):: &  !  org
      & 0.007577,0.003312,&
      & 0.136371,0.011163,0.017129,0.008140,0.005826,0.006777,0.007108,0.008132,&
      & 0.134992,0.147417,0.085253,0.054238,0.033790,0.032862,0.029038,0.026555,&
      & 0.141595,0.207980,&
      & 0.223252,0.193038,0.167892,0.148726,0.140473,0.130220,0.121166,0.113839,0.121855,0.107138,&
      & 0.105637,0.086639,0.075084,0.075089,0.070868,0.068706]
   
   real(wp), parameter :: emiss_hf_def2mtzvpp(*) = [real(wp):: &    !SG
      & 0.027000,0.000000,&
      & 0.000000,0.000000,0.200000,0.020000,0.180000,0.080000,0.070000,0.065000,&
      & 0.000000,0.000000,0.000000,0.200000,0.600000,0.600000,0.600000,0.300000,&
      & 0.000000,0.000000,&
      & 0.300000,0.300000,0.300000,0.300000,0.300000,0.300000,0.300000,0.300000,0.300000,0.300000,&
      & 0.300000,0.300000,0.300000,0.300000,0.300000,0.000000]
   
   real(wp), parameter :: emiss_hf_2g(*) = [real(wp):: & !no ne, ar ecp
      & 0.0539181,0.161846,&
      & 0.1581960,0.214318,0.808955,0.470398,0.724457,1.260960,2.014430,0.000000,&
      & 0.3072290,0.265975,0.354139,0.307818,0.356962,0.461661,0.588346,0.000000,&
      & 0.000000,0.000000,&
      & 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,&
      & 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000]
   
   real(wp), parameter :: emiss_hf_ccdz(*) = [real(wp):: &
      & 0.007907,0.008287,&
      & 0.047380,0.014240,0.022133,0.014999,0.018148,0.028240,0.042261,0.061485,&
      & 0.073185,0.056218,0.082660,0.052975,0.033874,0.034056,0.031433,0.030034,&
      & 0.000000,0.078016,& !no k cc-pVDZ Basis
      & 0.036885,0.038540,0.036474,0.036061,0.030289,0.027959,0.025177,0.022709,0.027386,0.015816,&
      & 0.135176,0.115515,0.102761,0.102967,0.097891,0.097331]
   
   real(wp), parameter :: emiss_hf_accdz(*) = [real(wp):: & !for li,be,na,mg,k-zn energy below def2-QZVPD reference
      & 0.001183,0.005948,&
      & 0.000000,0.000000,0.005269,0.006380,0.011700,0.021199,0.034160,0.051481,&
      & 0.000000,0.000000,0.016018,0.009268,0.010076,0.015153,0.016889,0.018563,&
      & 0.000000,0.000000,&
      & 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,&
      & 0.069963,0.065687,0.072944,0.077585,0.078777,0.080746]
   
   real(wp), parameter :: emiss_hf_dzp(*) = [real(wp):: &
      & 0.008107,0.008045,&
      & 0.136751,0.016929,0.026729,0.021682,0.027391,0.040841,0.058747,0.082680,&
      & 0.153286,0.162296,0.102704,0.073144,0.056217,0.061333,0.065045,0.071398,&
      & 0.145642,0.212865,&
      & 0.232821,0.204796,0.182933,0.169554,0.164701,0.160112,0.157723,0.158037,0.179104,0.169782,&
      & 0.159396,0.140611,0.129645,0.132664,0.132121,0.134081]
   
   real(wp), parameter :: emiss_hf_hsv(*) = [real(wp):: &
      & 0.030224,0.028324,&
      & 0.125379,0.064094,0.059751,0.079387,0.108929,0.167264,0.245786,0.347818,&
      & 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,&
      & 0.000000,0.000000,&
      & 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,&
      & 0.000000,0.000000,0.000000,0.000000,0.000000,0.000000]
   
   real(wp), parameter :: emiss_hf_dz(*) = [real(wp):: &
      & 0.009037,0.008843,&
      & 0.198254,0.000000,0.026921,0.021817,0.027458,0.041391,0.059495,0.083286,&
      & 0.268608,0.202374,0.104146,0.075686,0.057826,0.065300,0.069912,0.076845,&
      & 0.296046,0.370399,&
      & 0.349482,0.302284,0.267639,0.244306,0.232237,0.221488,0.214153,0.032694,0.226865,0.213902,&
      & 0.172296,0.155496,0.143646,0.149642,0.149871,0.151705]
   
   real(wp), parameter :: emiss_hf_msvp(*) = [real(wp):: &     !H-Kr modified Ahlrichs DZ, supplemented by def2-SV(P)
      & 0.000000_wp,0.000000_wp,& !RG,H set to zero,  F adjusted empirically, Be corrected due to ROHF problems
      & 0.107750_wp,0.020000_wp,0.026850_wp,0.021740_wp,0.027250_wp,0.039930_wp,0.030000_wp,0.000000_wp,&
      & 0.153290_wp,0.162300_wp,0.102700_wp,0.073140_wp,0.056220_wp,0.061330_wp,0.065040_wp,0.000000_wp,&
      & 0.200960_wp,0.299660_wp,&
      & 0.325990_wp,0.305490_wp,0.291720_wp,0.293800_wp,0.291790_wp,0.296730_wp,0.304600_wp,0.242040_wp,0.354190_wp,0.350720_wp,&
      & 0.350020_wp,0.345780_wp,0.349530_wp,0.367310_wp,0.382010_wp,0.000000_wp]

   real(wp), parameter :: emiss_hf_pbeh3c(*) = [ &
      & emiss_hf_msvp(1:18), emiss_hf_dzp(19:35), 0.0_wp]
   
   ! *********************
   ! * nr. of basis fkt  *
   ! *********************
   integer, parameter :: nbas_sv(*) = [ &
      & 2, 2, &
      & 3, 3, 9, 9, 9, 9, 9, 9, &
      & 7, 7, 13, 13, 13, 13, 13, 13, &
      & 11, 11, &
      & 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, &
      & 27, 27, 27, 27, 27, 27]
   integer, parameter :: nbas_minis(*) = [ &
      & 1, 1, &
      & 2, 2, 5, 5, 5, 5, 5, 5, &
      & 6, 6, 9, 9, 9, 9, 9, 9, &
      & 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0]
   integer, parameter :: nbas_631gd(*) = [ &
      & 2, 5, &
      & 14, 14, 14, 14, 14, 14, 14, 14, &
      & 18, 18, 18, 18, 18, 18, 18, 18, &
      & 22, 22, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 32, 0]
   integer, parameter :: nbas_svp(*) = [ &
      & 5, 5, &
      & 9, 9, 14, 14, 14, 14, 14, 14, &
      & 15, 18, 18, 18, 18, 18, 18, 18, &
      & 24, 24, &
      & 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, &
      & 32, 32, 32, 32, 32, 32]
   integer, parameter :: nbas_sv_p(*) = [ &
      & nbas_sv(1), nbas_svp(2:)]
   integer, parameter :: nbas_svx(*) = [ &
      & nbas_sv(1), nbas_svp(2:5), nbas_sv(6), nbas_svp(7:)]
   integer, parameter :: nbas_svp_old(*) = [ &
      & 5, 5, &
      & 6, 9, 14, 14, 14, 14, 14, 14, &
      & 10, 10, 18, 18, 18, 18, 18, 18, &
      & 14, 24, &
      & 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, &
      & 32, 32, 32, 32, 32, 32]
   integer, parameter :: nbas_tz(*) = [ &
      & 6, 6, &
      & 14, 19, 31, 31, 31, 31, 31, 31, &
      & 32, 32, 37, 37, 37, 37, 37, 37, &
      & 33, 36, &
      & 45, 45, 45, 45, 45, 45, 45, 45, 45, 48, &
      & 48, 48, 48, 48, 48, 48]
   integer, parameter :: nbas_def2mtzvp(*) = [ &
      & 3, 3, &
      & 8, 11, 19, 19, 19, 24, 19, 19, &
      & 14, 14, 22, 22, 22, 22, 22, 22, &
      & 18, 28, &
      & 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, &
      & 36, 36, 36, 36, 36, 36]   !def2-mTZVP
   integer, parameter :: nbas_def2mtzvpp(*) = [ &
      & 5, 5, &
      & 9, 11, 19, 19, 24, 24, 24, 24, &
      & 14, 14, 27, 27, 27, 27, 27, 27, &
      & 18, 28, &
      & 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, &
      & 36, 36, 36, 36, 36, 36]  !def2-mTZVPP
   integer, parameter :: nbas_vmb(*) = [ &
      & 1, 1, &
      & 2, 2, 4, 4, 4, 4, 4, 4, &
      & 1, 1, 4, 4, 4, 4, 4, 4, &
      & 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0] ! minimal basis set with ECPs
   integer, parameter :: nbas_minisd(*) = [ &
      & 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 14, 14, 14, 14, 14, 14, &
      & 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0]
   integer, parameter :: nbas_hsv(*) = [ &  ! FIXME: https://github.com/grimme-lab/gcp/issues/23
      & 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0]
   integer, parameter :: nbas_minix(*) = [ &
      & nbas_minis(1:2), 5, 5, nbas_minis(5:10), 9, 9, &
      & nbas_minisd(13:18), nbas_sv(19:30), nbas_svp(31:36)]
   integer, parameter :: nbas_lanl2(*) = [ &
      & nbas_631gd(1:20), &
      & 22, 22, 22, 22, 22, 22, 22, 22, 22, 18, &
      & nbas_631gd(31:)] ! Sc-Zn LANL2DZ
   integer, parameter :: nbas_pobtz(*) = [ & ! H-Kr no RG
      & 6, 0, &
      & 7, 7, 18, 18, 18, 18, 18, 0, &
      & 19, 19, 22, 22, 22, 22, 22, 0, &
      & 23, 23, &
      & 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, &
      & 41, 41, 41, 41, 41, 0 ]
   
   integer, parameter :: nbas_pobdzvp(*) = [ & ! H-KR, no RG
      & 5, 0, &
      & 6, 6, 14, 14, 14, 14, 14, 0, &
      & 15, 18, 18, 18, 18, 18, 18, 0, &
      & 19, 19, &
      & 31, 31, 31, 31, 31, 31, 31, 31, 31, 31,&
      & 32, 32, 32, 32, 32, 0]
   
   integer, parameter :: nbas_2gcore(*) = [ & ! Only HCNOF yet
      & 1, 0, &
      & 5, 5, 5, 5, 5, 5, 5, 5, &
      & 9, 9, 14, 14, 14, 14, 14, 14, &
      & 0, 0,&
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0]
   
   integer, parameter :: nbas_2g(*) = [ &
      & 1, 1, &
      & 5, 5, 5, 5, 5, 5, 5, 5, &
      & 9, 9, 14, 14, 14, 14, 14, 14, &
      & 0, 0, &
      & 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
      & 0, 0, 0, 0, 0, 0]
   
   integer, parameter :: nbas_def1tzvp(*) = [ &  ! FIXME: https://github.com/grimme-lab/gcp/issues/24
      & 6, 6, &
      & 8, 11, 19, 19, 19, 19, 19, 19, &
      & 14, 14, 22, 22, 22, 22, 22, 22, &
      & 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, &
      & 18, 28,&
      & 36, 36, 36, 36, 36, 36]
   
   integer, parameter :: nbas_ccdz(*) = [ &
      & 5, 5, &
      & 14, 14, 14, 14, 14, 14, 14, 14, &
      & 18, 18, 18, 18, 18, 18, 18, 18, &
      & 0,27,&
      & 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, &
      & 27, 27, 27, 27, 27, 27]
   
   integer, parameter :: nbas_accdz(*) = [ &
      & 9, 9,&
      & 23, 23, 23, 23, 23, 23, 23, 23,&
      & 27, 27, 27, 27, 27, 27, 27, 27,&
      & 0, 0,&
      & 59, 59, 59, 59, 59, 59, 59, 59, 59, 59,&
      & 36, 36, 36, 36, 36, 36]
   
   integer, parameter :: nbas_dzp(*) = [ &
      & 5, 5,&
      & 7, 10, 15, 15, 15, 15, 15, 15,&
      & 15, 15, 23, 23, 23, 23, 23, 23,&
      & 26, 36,&
      & 41, 41, 41, 41, 41, 41, 41, 41, 41, 41,&
      & 41, 41, 41, 41, 41, 41]
   
   integer, parameter :: nbas_dz(*) = [ &
      & 2, 2, &
      & 4, 0, 10, 10, 10, 10, 10, 10, &
      & 12, 12, 18, 18, 18, 18, 18, 18, &
      & 23, 23, &
      & 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, &
      & 36, 36, 36, 36, 36, 36]
   
   integer, parameter :: nbas_msvp(*) = [ &  ! modified Ahlrichs DZ, supplemented by def2-SV(P)
      & 2, 2,&
      & 10, 10, 15, 15, 15, 15, 15, 15,&
      & 15, 18, 18, 18, 18, 18, 18, 18,&
      & 24, 24, &
      & 31, 31, 31, 31, 31, 31, 31, 31, 31, 31,&
      & 32, 32, 32, 32, 32, 32]

   real(wp), parameter :: slater_s(*) = [real(wp) :: &
      & 1.2000,1.6469,0.6534,1.0365,1.3990,1.7210,2.0348,2.2399,2.5644,2.8812,&
      & 0.8675,1.1935,1.5143,1.7580,1.9860,2.1362,2.3617,2.5796,0.9362,1.2112,&
      & 1.2870,1.3416,1.3570,1.3804,1.4761,1.5465,1.5650,1.5532,1.5781,1.7778,&
      & 2.0675,2.2702,2.4546,2.5680,2.7523,2.9299]
   real(wp), parameter :: slater_p(*) = [real(wp) :: &
      & 0.0000,0.0000,0.5305,0.8994,1.2685,1.6105,1.9398,2.0477,2.4022,2.7421,&
      & 0.6148,0.8809,1.1660,1.4337,1.6755,1.7721,2.0176,2.2501,0.6914,0.9329,&
      & 0.9828,1.0104,0.9947,0.9784,1.0641,1.1114,1.1001,1.0594,1.0527,1.2448,&
      & 1.5073,1.7680,1.9819,2.0548,2.2652,2.4617]
   real(wp), parameter :: slater_d(*) = [real(wp) :: &
      & 0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,&
      & 0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,&
      & 2.4341,2.6439,2.7809,2.9775,3.2208,3.4537,3.6023,3.7017,3.8962,2.0477,&
      & 2.4022,2.7421,0.6148,0.8809,1.1660,1.4337]
   real(wp), parameter :: slater_exp(*) = [ &
      slater_s(1:2), (slater_s(3:20) + slater_p(3:20))/2, &
      (slater_s(21:30) + slater_d(21:30) + slater_d(21:30))/3, slater_s(31:36)]

contains

subroutine get_gcp_param(param, mol, method, basis, eta)
   type(gcp_param), intent(out) :: param
   type(structure_type), intent(in) :: mol
   character(len=*), intent(in), optional :: method
   character(len=*), intent(in), optional :: basis
   real(wp), intent(in), optional :: eta

   real(wp) :: eta_, eta_spec
   real(wp), allocatable :: nbas(:)
   integer :: isp, izp, jsp, basis_id, method_id
   integer, allocatable :: nel(:)
   logical :: valence_minimal_basis

   method_id = p_unknown_method
   if (present(method)) then
      method_id = get_method_id(method)
      basis_id = get_basis_id(method)
   else
      basis_id = p_unknown_bas
   end if
   if (present(basis)) basis_id = get_basis_id(basis)
   if (basis_id == p_unknown_bas .and. method_id == p_unknown_method) return

   eta_ = 0.0_wp
   if (present(eta)) eta_ = eta
   eta_spec = 0.0_wp

   allocate(param%zeff(mol%nid))
   do isp = 1, mol%nid
      izp = mol%num(isp)
      param%zeff(isp) = effective_atomic_number(izp)
   end do

   allocate(param%rvdw(mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         param%rvdw(isp, jsp) = get_vdw_rad(param%zeff(isp), param%zeff(jsp))
      end do
   end do

   allocate(param%rvdw_srb(mol%nid, mol%nid))
   do isp = 1, mol%nid
      do jsp = 1, mol%nid
         param%rvdw_srb(isp, jsp) = get_vdw_rad(mol%num(isp), mol%num(jsp))
      end do
   end do

   valence_minimal_basis = basis_id == p_vmb_bas

   allocate(nel(mol%nid))
   do isp = 1, mol%nid
      izp = param%zeff(isp)
      nel(isp) = number_of_electrons(izp, valence_minimal_basis)
   end do

   select case(basis_id)
   case(p_vmb_bas)
      param%emiss = emiss_hf_vmb(param%zeff)
      nbas = nbas_vmb(param%zeff)

   case(p_sv_bas)
      param%emiss = emiss_hf_sv(param%zeff)
      nbas = nbas_sv(param%zeff)

      select case(method_id)
      case (p_hf_method) ! RMS=0.3218975
         param%sigma = 0.1724_wp
         eta_ = 1.2804_wp
         param%alpha = 0.8568_wp
         param%beta = 1.2342_wp
      case (p_hyb_method, p_b3lyp_method, p_pw6b95_method) ! RMS= 0.557
         param%sigma = 0.4048_wp
         eta_ = 1.1626_wp
         param%alpha = 0.8652_wp
         param%beta = 1.2375_wp
      case (p_gga_method, p_tpss_method, p_blyp_method) ! RMS = 0.6652
         param%sigma = 0.2727_wp
         eta_ = 1.4022_wp
         param%alpha = 0.8055_wp
         param%beta = 1.3000_wp
      end select

   case(p_sv_p_bas)
      param%emiss = emiss_hf_sv_p(param%zeff)
      nbas = nbas_sv_p(param%zeff)

      if (method_id == p_hf_method) then !RMS=0.3502
         param%sigma = 0.1373_wp
         eta_ = 1.4271_wp
         param%alpha = 0.8141_wp
         param%beta = 1.2760_wp
      elseif (is_dft_method(method_id)) then ! RMS= 0.57 ! def2-SV(P)
         param%sigma = 0.2424_wp
         eta_ = 1.2371_wp
         param%alpha = 0.6076_wp
         param%beta = 1.4078_wp
      end if

   case(p_svx_bas) ! RMS=  ! def2-SV(P/h,c)  = SV at h,c
      param%emiss = emiss_hf_svx(param%zeff)
      nbas = nbas_svx(param%zeff)

      if (is_dft_method(method_id)) then ! RMS=  0.56 ! def2-SV(P/h,c)  = SV at h,c
         param%sigma = 0.1861_wp
         eta_ = 1.3200_wp
         param%alpha = 0.6171_wp
         param%beta = 1.4019_wp
      end if

   case(p_svp_bas)
      param%emiss = emiss_hf_svp(param%zeff)
      nbas = nbas_svp(param%zeff)

      if (method_id == p_hf_method) then ! RMS=0.4065
         param%sigma = 0.2054_wp
         eta_ = 1.3157_wp
         param%alpha = 0.8136_wp
         param%beta = 1.2572_wp
      elseif (method_id == p_tpss_method) then ! RMS=  0.618
         param%sigma = 0.6647_wp
         eta_ = 1.3306_wp
         param%alpha = 1.0792_wp
         param%beta = 1.1651_wp
      elseif (method_id == p_pw6b95_method) then  ! RMS = 0.58312
         param%sigma = 0.3098_wp
         eta_ = 1.2373_wp
         param%alpha = 0.6896_wp
         param%beta = 1.3347_wp
      elseif (is_hyb_method(method_id)) then ! RMS=0.6498
         param%sigma = 0.2990_wp
         eta_ = 1.2605_wp
         param%alpha = 0.6438_wp
         param%beta = 1.3694_wp
      elseif (is_gga_method(method_id)) then ! RMS=
         param%sigma = 0.6823_wp
         eta_ = 1.2491_wp
         param%alpha = 0.8225_wp
         param%beta = 1.2811_wp
      end if

   case(p_svp_old_bas)
      param%emiss = emiss_hf_svp_old(param%zeff)
      nbas = nbas_svp_old(param%zeff)

      if (method_id == p_hf_method) then ! RMS=0.4065
         param%sigma = 0.2054_wp
         eta_ = 1.3157_wp
         param%alpha = 0.8136_wp
         param%beta = 1.2572_wp
      elseif (is_dft_method(method_id)) then ! RMS=0.6498
         param%sigma = 0.2990_wp
         eta_ = 1.2605_wp
         param%alpha = 0.6438_wp
         param%beta = 1.3694_wp
      end if

   case(p_minis_bas)
      param%emiss = emiss_hf_minis(param%zeff)
      nbas = nbas_minis(param%zeff)

      if (method_id == p_hf_method) then ! RMS= 0.3040
         param%sigma = 0.1290_wp
         eta_ = 1.1526_wp
         param%alpha = 1.1549_wp
         param%beta = 1.1763_wp
      elseif (method_id == p_tpss_method) then ! RMS=
         param%sigma = 0.22982_wp
         eta_ = 1.35401_wp
         param%alpha = 1.47633_wp
         param%beta = 1.11300_wp
      elseif (method_id == p_pw6b95_method) then  ! RMS = 0.3279929
         param%sigma = 0.21054_wp
         eta_ = 1.25458_wp
         param%alpha = 1.35003_wp
         param%beta = 1.14061_wp
      elseif (is_gga_method(method_id)) then ! RMS= 0.3462
         param%sigma = 0.1566_wp
         eta_ = 1.0271_wp
         param%alpha = 1.0732_wp
         param%beta = 1.1968_wp
      elseif (is_hyb_method(method_id)) then ! RMS= 0.3400
         param%sigma = 0.2059_wp
         eta_ = 0.9722_wp
         param%alpha = 1.1961_wp
         param%beta = 1.1456_wp
      end if

   case(p_631gd_bas)
      param%emiss = emiss_hf_631gd(param%zeff)
      nbas = nbas_631gd(param%zeff)

      if (method_id == p_hf_method) then ! RMS= 0.40476
         param%sigma = 0.2048_wp
         eta_ = 1.5652_wp
         param%alpha = 0.9447_wp
         param%beta = 1.2100_wp
      elseif (is_dft_method(method_id)) then ! RMS=  0.47856
         param%sigma = 0.3405_wp
         eta_ = 1.6127_wp
         param%alpha = 0.8589_wp
         param%beta = 1.2830_wp
      end if

   case(p_tz_bas)
      param%emiss = emiss_hf_tz(param%zeff)
      nbas = nbas_tz(param%zeff)

      if (method_id == p_hf_method) then !  RMS= 0.1150
         param%sigma = 0.3127_wp
         eta_ = 1.9914_wp
         param%alpha = 1.0216_wp
         param%beta = 1.2833_wp
      elseif (is_hyb_method(method_id)) then ! RMS=0.19648
         param%sigma = 0.2905_wp
         eta_ = 2.2495_wp
         param%alpha = 0.8120_wp
         param%beta = 1.4412_wp
      elseif (is_gga_method(method_id)) then !RMS = 0.21408
         param%sigma = 0.1182_wp
         eta_ = 1.0631_wp
         param%alpha = 1.0510_wp
         param%beta = 1.1287_wp
      end if

   case(p_deftzvp_bas)
      param%emiss = emiss_hf_def1tzvp(param%zeff)
      nbas = nbas_def1tzvp(param%zeff)

      if (method_id == p_hf_method) then ! RMS=0.209
         param%sigma = 0.2600_wp
         eta_ = 2.2448_wp
         param%alpha = 0.7998_wp
         param%beta = 1.4381_wp
      elseif (is_dft_method(method_id)) then ! RMS=0.1817
         param%sigma = 0.2393_wp
         eta_ = 2.2247_wp
         param%alpha = 0.8185_wp
         param%beta = 1.4298_wp
      end if

   case(p_ccdz_bas)
      param%emiss = emiss_hf_ccdz(param%zeff)
      nbas = nbas_ccdz(param%zeff)

      if (method_id == p_hf_method) then ! RMS=0.4968
         param%sigma = 0.4416_wp
         eta_ = 1.5185_wp
         param%alpha = 0.6902_wp
         param%beta = 1.3713_wp
      elseif (is_dft_method(method_id)) then ! RMS=0.7610
         param%sigma = 0.5383_wp
         eta_ = 1.6482_wp
         param%alpha = 0.6230_wp
         param%beta = 1.4523_wp
      end if

   case(p_accdz_bas)
      param%emiss = emiss_hf_accdz(param%zeff)
      nbas = nbas_accdz(param%zeff)

      if (method_id == p_hf_method) then !RMS=0.2222
         param%sigma = 0.0748_wp
         eta_ = 0.0663_wp
         param%alpha = 0.3811_wp
         param%beta = 1.0155_wp
      elseif (is_dft_method(method_id)) then ! RMS=0.1840
         param%sigma = 0.1465_wp
         eta_ = 0.0500_wp
         param%alpha = 0.6003_wp
         param%beta = 0.8761_wp
      end if

   case(p_pobtz_bas)
      param%emiss = emiss_hf_pobtz(param%zeff)
      nbas = nbas_pobtz(param%zeff)

      if (is_dft_method(method_id)) then
         param%sigma = 0.1300_wp
         eta_ = 1.3743_wp
         param%alpha = 0.4792_wp
         param%beta = 1.3962_wp
      end if

   ! case(p_pobdzvp_bas)
   !    param%emiss = emiss_hf_pobdzvp(param%zeff)
   !    nbas = nbas_pobdzvp(param%zeff)

   case(p_minix_bas)
      param%emiss = emiss_hf_minix(param%zeff)
      nbas = nbas_minix(param%zeff)

      if (method_id == p_hf_method .or. method_id == p_hf3c_method) then
         param%sigma = 0.1290_wp
         eta_ = 1.1526_wp
         param%alpha = 1.1549_wp
         param%beta = 1.1763_wp
         param%base = method_id == p_hf3c_method
         if (param%base) then
            param%rscal = 0.7_wp
            param%qscal = 0.03_wp
         end if
      elseif (is_dft_method(method_id)) then
         param%sigma = 0.2059_wp
         eta_ = 0.9722_wp
         param%alpha = 1.1961_wp
         param%beta = 1.1456_wp
      end if

   case(p_gcore_bas)
      param%emiss = emiss_hf_2gcore(param%zeff)
      nbas = nbas_2gcore(param%zeff)

   case(p_2g_bas)
      param%emiss = emiss_hf_2g(param%zeff)
      nbas = nbas_2g(param%zeff)

      if (method_id == p_hf_method) then
         param%sigma = 0.2461_wp
         eta_ = 1.1616_wp
         param%alpha = 0.7335_wp
         param%beta = 1.4709_wp
      end if

   case(p_dzp_bas)
      param%emiss = emiss_hf_dzp(param%zeff)
      nbas = nbas_dzp(param%zeff)

      if (method_id == p_hf_method) then !RMS=0.4571
         param%sigma = 0.1443_wp
         eta_ = 1.4547_wp
         param%alpha = 0.3711_wp
         param%beta = 1.6300_wp
      elseif (is_dft_method(method_id)) then !RMS=0.7184
         param%sigma = 0.2687_wp
         eta_ = 1.4634_wp
         param%alpha = 0.3513_wp
         param%beta = 1.6880_wp
      end if

   ! case(p_hsv_bas)
   !    param%emiss = emiss_hf_hsv(param%zeff)
   !    nbas = nbas_hsv(param%zeff)

   case(p_dz_bas)
      param%emiss = emiss_hf_dz(param%zeff)
      nbas = nbas_dz(param%zeff)

      if (method_id == p_hf_method) then  !RMS=0.3754
         param%sigma = 0.1059_wp
         eta_ = 1.4554_wp
         param%alpha = 0.3711_wp
         param%beta = 1.6342_wp
      elseif (is_dft_method(method_id)) then
         param%sigma = 0.2687_wp
         eta_ = 1.4634_wp
         param%alpha = 0.3513_wp
         param%beta = 1.6880_wp
      end if

   case(p_msvp_bas)
      param%emiss = emiss_hf_msvp(param%zeff)
      nbas = nbas_msvp(param%zeff)

   case(p_def2mtzvp_bas)
      param%emiss = emiss_hf_def2mtzvp(param%zeff)
      nbas = nbas_def2mtzvp(param%zeff)

      if (method_id == p_b3pbe3c_method) then
         param%sigma = 1.0000_wp
         eta_ = 2.98561_wp
         param%alpha = 0.3011_wp
         param%beta = 2.4405_wp
      end if

   case(p_pbeh3c_bas)
      param%emiss = emiss_hf_pbeh3c(param%zeff)
      nbas = nbas_msvp(param%zeff)

      if (method_id == p_pbeh3c_method) then
         param%sigma = 1.00000_wp
         eta_ = 1.32492_wp
         param%alpha = 0.27649_wp
         param%beta = 1.95600_wp
         param%damp = .true.
      elseif (method_id == p_hse3c_method) then
         param%sigma = 1.00000_wp
         eta_ = 1.32378_wp
         param%alpha = 0.28314_wp
         param%beta = 1.94527_wp
         param%damp = .true.
      end if

   case(p_lanl2_bas)
      param%emiss = emiss_hf_lanl2(param%zeff)
      nbas = nbas_lanl2(param%zeff)

      if (is_dft_method(method_id)) then
         param%sigma = 0.3405_wp
         eta_ = 1.6127_wp
         param%alpha = 0.8589_wp
         param%beta = 1.2830_wp
      end if

   case(p_def2mtzvpp_bas)
      param%emiss = emiss_hf_def2mtzvpp(param%zeff)
      nbas = nbas_def2mtzvpp(param%zeff)

      if (method_id == p_r2scan3c_method) then
         param%sigma = 1.0000_wp
         eta_ = 1.3150_wp
         eta_spec = 1.15_wp
         param%alpha = 0.9410_wp
         param%beta = 1.4636_wp
         param%damp = .true.
      end if

   end select
   if (allocated(nbas)) then
      param%xv = nbas - 0.5_wp * nel
      if (basis_id == p_def2mtzvpp_bas) then
         param%xv(:) = 1.0_wp
         where(param%zeff == 6) param%xv = 3.0_wp
         where(param%zeff == 7) param%xv = 0.5_wp
         where(param%zeff == 8) param%xv = 0.5_wp
      end if
   end if

   if (eta_ > 0.0_wp) then
      param%slater = eta_ * slater_exp(param%zeff)
      if (eta_spec > 0.0_wp) then
         where(param%zeff > 10) param%slater = eta_spec * param%slater
      end if
   end if

   if (method_id == p_b973c_method) then
      param%srb=.true.
      param%rscal=10.00_wp
      param%qscal=0.08_wp
   end if
end subroutine get_gcp_param

pure function get_method_id(method) result(id)
   character(len=*), intent(in) :: method
   integer :: id
   select case(method)
   case default
      id = p_unknown_method
   case('hf')
      id = p_hf_method
   case('dft')
      id = p_dft_method
   case('b3lyp')
      id = p_b3lyp_method
   case('blyp')
      id = p_blyp_method
   case('gga')
      id = p_gga_method
   case('tpss')
      id = p_tpss_method
   case('pw6b95')
      id = p_pw6b95_method
   case('b973c')
      id = p_b973c_method
   case('r2scan3c')
      id = p_r2scan3c_method
   case('pbe')
      id = p_pbe_method
   case('pbeh3c')
      id = p_pbeh3c_method
   case('hse3c')
      id = p_hse3c_method
   case('hf3c')
      id = p_hf3c_method
   case('b3pbe3c')
      id = p_b3pbe3c_method
   end select
end function get_method_id

pure function is_dft_method(method_id) result(res)
   integer, intent(in) :: method_id
   logical :: res
   res = method_id == p_dft_method .or. method_id == p_b3lyp_method .or. &
         method_id == p_blyp_method .or. method_id == p_gga_method .or. &
         method_id == p_tpss_method .or. method_id == p_pw6b95_method
end function is_dft_method

pure function is_hyb_method(method_id) result(res)
   integer, intent(in) :: method_id
   logical :: res
   res = method_id == p_b3lyp_method .or. method_id == p_pw6b95_method
end function is_hyb_method

pure function is_gga_method(method_id) result(res)
   integer, intent(in) :: method_id
   logical :: res
   res = method_id == p_gga_method .or. method_id == p_tpss_method .or. &
         method_id == p_blyp_method
end function is_gga_method

pure function get_basis_id(basis) result(id)
   character(len=*), intent(in) :: basis
   integer :: id
   select case(basis)
   case default
      id = p_unknown_bas
   case('sv')
      id = p_sv_bas
   case('sv(p)','def2sv(p)', 'sv_p', 'def2sv_p')
      id = p_sv_p_bas
   case ('svx') ! RMS=  ! def2-SV(P/h,c)  = SV at h,c
      id = p_svx_bas
   case('svp')
      id = p_svp_bas
   case('minis')
      id = p_minis_bas
   case('631gd', '631gs')
      id = p_631gd_bas
   case('tz', 'def2tzvp')
      id = p_tz_bas
   case('deftzvp', 'def1tzvp')
      id = p_deftzvp_bas
   case('ccdz', 'ccpvdz')
      id = p_ccdz_bas
   case('accdz', 'augccpvdz', 'accpvdz')
      id = p_accdz_bas
   case('pobtz', 'pobtzvp')
      id = p_pobtz_bas
   ! case('pobdzvp')
   !    id = p_pobdzvp_bas
   case ('minix', 'hf3c')
      id = p_minix_bas
   case('gcore')
      id = p_gcore_bas
   case('2g', 'twog', 'fitg')
      id = p_2g_bas
   case('dzp')
      id = p_dzp_bas
   ! case('hsv')
   !    id = p_hsv_bas
   case('lanl')
      id = p_lanl2_bas
   case('dz')
      id = p_dz_bas
   case('msvp', 'def2msvp')
      id = p_msvp_bas
   case('pbeh3c', 'hse3c')
      id = p_pbeh3c_bas
   case('mtzvp', 'def2mtzvp')
      id = p_def2mtzvp_bas
   case('mtzvpp', 'def2mtzvpp', 'r2scan3c')
      id = p_def2mtzvpp_bas
   end select
end function get_basis_id

pure function effective_atomic_number(number) result(zeff)
   integer, intent(in) :: number
   integer :: zeff

   select case (number)
   case default
      zeff = number
   case(37:54)
      zeff = number - 18
   case(55:57)
      zeff = number - 18*2
   case(58:71, 90:94)
      zeff = 21
   case(72:89)
      zeff = number - 2*18 - 14
   end select
end function effective_atomic_number

pure function number_of_electrons(number, valence_minimal_basis) result(nel)
   integer, intent(in) :: number
   logical, intent(in) :: valence_minimal_basis
   integer :: nel

   if (valence_minimal_basis) then
      select case(number)
      case default
         nel = number
      case(5:10)
         nel = number - 2
      case(11:18)
         nel = number - 10
      end select
   else
      nel = number
   end if
end function number_of_electrons

endmodule dftd3_gcp_param