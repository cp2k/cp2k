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

!> @dir tblite/solvation/data
!> Provides solvent specific data for parametrized solvation models

!> @file tblite/solvation/data.f90
!> Provides element specific data

!> Element specific data for solvation models
module tblite_solvation_data
   use mctc_env, only : wp
   use mctc_io, only : to_number
   use mctc_io_convert, only : aatoau
   use tblite_utils, only : compact
   implicit none
   private

   public :: get_vdw_rad_d3, get_vdw_rad_cosmo, get_vdw_rad_bondi
   public :: get_solvent_data


   !> Solvent specific parameters
   type, public :: solvent_data
      real(wp) :: eps
      real(wp) :: n
      real(wp) :: alpha
      real(wp) :: beta
      real(wp) :: gamma
      real(wp) :: phi2
      real(wp) :: psi2
      character(len=:), allocatable :: solvent
   end type solvent_data

   interface solvent_data
      module procedure get_solvent_data
   end interface solvent_data


   !> Get D3 van-der-Waals radius for a species
   interface get_vdw_rad_d3
      module procedure :: get_vdw_rad_d3_sym
      module procedure :: get_vdw_rad_d3_num
   end interface get_vdw_rad_d3

   !> Get COSMO van-der-Waals radius for a species
   interface get_vdw_rad_cosmo
      module procedure :: get_vdw_rad_cosmo_sym
      module procedure :: get_vdw_rad_cosmo_num
   end interface get_vdw_rad_cosmo

   !> Get Bondi van-der-Waals radius for a species
   interface get_vdw_rad_bondi
      module procedure :: get_vdw_rad_bondi_sym
      module procedure :: get_vdw_rad_bondi_num
   end interface get_vdw_rad_bondi


   !> In case no van-der-Waals value is provided
   real(wp), parameter :: missing = -1.0_wp

   !> D3 pairwise van-der-Waals radii (only homoatomic pairs present here)
   real(wp), parameter :: vdw_rad_d3(1:94) = aatoau * [&
      & 1.09155_wp, 0.86735_wp, 1.74780_wp, 1.54910_wp, &  ! H-Be
      & 1.60800_wp, 1.45515_wp, 1.31125_wp, 1.24085_wp, &  ! B-O
      & 1.14980_wp, 1.06870_wp, 1.85410_wp, 1.74195_wp, &  ! F-Mg
      & 2.00530_wp, 1.89585_wp, 1.75085_wp, 1.65535_wp, &  ! Al-S
      & 1.55230_wp, 1.45740_wp, 2.12055_wp, 2.05175_wp, &  ! Cl-Ca
      & 1.94515_wp, 1.88210_wp, 1.86055_wp, 1.72070_wp, &  ! Sc-Cr
      & 1.77310_wp, 1.72105_wp, 1.71635_wp, 1.67310_wp, &  ! Mn-Ni
      & 1.65040_wp, 1.61545_wp, 1.97895_wp, 1.93095_wp, &  ! Cu-Ge
      & 1.83125_wp, 1.76340_wp, 1.68310_wp, 1.60480_wp, &  ! As-Kr
      & 2.30880_wp, 2.23820_wp, 2.10980_wp, 2.02985_wp, &  ! Rb-Zr
      & 1.92980_wp, 1.87715_wp, 1.78450_wp, 1.73115_wp, &  ! Nb-Ru
      & 1.69875_wp, 1.67625_wp, 1.66540_wp, 1.73100_wp, &  ! Rh-Cd
      & 2.13115_wp, 2.09370_wp, 2.00750_wp, 1.94505_wp, &  ! In-Te
      & 1.86900_wp, 1.79445_wp, 2.52835_wp, 2.59070_wp, &  ! I-Ba
      & 2.31305_wp, 2.31005_wp, 2.28510_wp, 2.26355_wp, &  ! La-Nd
      & 2.24480_wp, 2.22575_wp, 2.21170_wp, 2.06215_wp, &  ! Pm-Gd
      & 2.12135_wp, 2.07705_wp, 2.13970_wp, 2.12250_wp, &  ! Tb-Er
      & 2.11040_wp, 2.09930_wp, 2.00650_wp, 2.12250_wp, &  ! Tm-Hf
      & 2.04900_wp, 1.99275_wp, 1.94775_wp, 1.87450_wp, &  ! Ta-Os
      & 1.72280_wp, 1.67625_wp, 1.62820_wp, 1.67995_wp, &  ! Ir-Hg
      & 2.15635_wp, 2.13820_wp, 2.05875_wp, 2.00270_wp, &  ! Tl-Po
      & 1.93220_wp, 1.86080_wp, 2.53980_wp, 2.46470_wp, &  ! At-Ra
      & 2.35215_wp, 2.21260_wp, 2.22970_wp, 2.19785_wp, &  ! Ac-U
      & 2.17695_wp, 2.21705_wp]                            ! Np-Pu


   !> Default value for unoptimized van-der-Waals radii
   real(wp), parameter :: cosmostub = 2.223_wp

   !> COSMO optimized van-der-Waals radii
   real(wp), parameter :: vdw_rad_cosmo(1:94) = aatoau * [ &
      & 1.3000_wp, 1.6380_wp, 1.5700_wp, 1.0530_wp, &   ! H-Be
      & 2.0480_wp, 2.0000_wp, 1.8300_wp, 1.7200_wp, &   ! B-O
      & 1.7200_wp, 1.8018_wp, 1.8000_wp, 1.6380_wp, &   ! F-Mg
      & 2.1530_wp, 2.2000_wp, 2.1060_wp, 2.1600_wp, &   ! Al-S
      & 2.0500_wp, 2.2000_wp, 2.2230_wp, cosmostub, &   ! Cl-Ca
      & cosmostub, 2.2930_wp, cosmostub, cosmostub, &   ! Sc-Cr
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Mn-Ni
      & cosmostub, 1.6260_wp, cosmostub, 2.7000_wp, &   ! Cu-Ge
      & 2.3500_wp, 2.2000_wp, 2.1600_wp, 2.3630_wp, &   ! As-Kr
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Rb-Zr
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Nb-Ru
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Rh-Cd
      & 2.2580_wp, 2.5500_wp, 2.4100_wp, 2.4100_wp, &   ! In-Te
      & 2.3200_wp, 2.5270_wp, cosmostub, cosmostub, &   ! I-Ba
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! La-Nd
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Pm-Gd
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Tb-Er
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Tm-Hf
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ta-Os
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ir-Hg
      & cosmostub, 2.3600_wp, 2.4220_wp, 2.3050_wp, &   ! Tl-Po
      & 2.3630_wp, 2.5740_wp, cosmostub, cosmostub, &   ! At-Ra
      & cosmostub, cosmostub, cosmostub, cosmostub, &   ! Ac-U
      & cosmostub, cosmostub]                           ! Np-Pu

   !> Van-der-Waals radii from
   !> Manjeera Mantina, Adam C. Chamberlin, Rosendo Valero, Christopher J. Cramer,
   !> and Donald G. Truhlar, Consistent van der Waals Radii for the Whole Main Group,
   !> J. Phys. Chem. A 2009, 113, 19, 5806–5812. https://doi.org/10.1021/jp8111556
   real(wp), parameter :: vdw_rad_bondi(1:88) = aatoau * [ &
      & 1.10_wp, 1.40_wp, 1.81_wp, 1.53_wp, 1.92_wp, 1.70_wp, 1.55_wp, 1.52_wp, &  ! H-O
      & 1.47_wp, 1.54_wp, 2.27_wp, 1.73_wp, 1.84_wp, 2.10_wp, 1.80_wp, 1.80_wp, &  ! F-S
      & 1.75_wp, 1.88_wp, 2.75_wp, 2.31_wp, missing, missing, missing, missing, &  ! Cl-Cr
      & missing, missing, missing, missing, missing, missing, 1.87_wp, 2.11_wp, &  ! Mn-Ge
      & 1.85_wp, 1.90_wp, 1.83_wp, 2.02_wp, 3.03_wp, 2.49_wp, missing, missing, &  ! As-Zr
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! Nb-Cd
      & 1.93_wp, 2.17_wp, 2.06_wp, 2.06_wp, 1.98_wp, 2.16_wp, 3.43_wp, 2.68_wp, &  ! I-Ba
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! La-Gd
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! Tb-Hf
      & missing, missing, missing, missing, missing, missing, missing, missing, &  ! Ta-Hg
      & 1.96_wp, 2.02_wp, 2.07_wp, 1.97_wp, 2.02_wp, 2.20_wp, 3.48_wp, 2.83_wp]    ! Tl-Ra


contains


!> Get van-der-Waals radius for species with a given symbol
elemental function get_vdw_rad_d3_sym(symbol) result(radius)
   !> Element symbol
   character(len=*), intent(in) :: symbol
   !> van-der-Waals radius
   real(wp) :: radius

   radius = get_vdw_rad_d3(to_number(symbol))

end function get_vdw_rad_d3_sym


!> Get van-der-Waals radius for species with a given atomic number
elemental function get_vdw_rad_d3_num(number) result(radius)
   !> Atomic number
   integer, intent(in) :: number
   !> van-der-Waals radius
   real(wp) :: radius

   if (number > 0 .and. number <= size(vdw_rad_d3, dim=1)) then
      radius = vdw_rad_d3(number)
   else
      radius = missing
   end if

end function get_vdw_rad_d3_num


!> Get van-der-Waals radius for species with a given symbol
elemental function get_vdw_rad_cosmo_sym(symbol) result(radius)
   !> Element symbol
   character(len=*), intent(in) :: symbol
   !> van-der-Waals radius
   real(wp) :: radius

   radius = get_vdw_rad_cosmo(to_number(symbol))

end function get_vdw_rad_cosmo_sym


!> Get van-der-Waals radius for species with a given atomic number
elemental function get_vdw_rad_cosmo_num(number) result(radius)
   !> Atomic number
   integer, intent(in) :: number
   !> van-der-Waals radius
   real(wp) :: radius

   if (number > 0 .and. number <= size(vdw_rad_cosmo, dim=1)) then
      radius = vdw_rad_cosmo(number)
   else
      radius = missing
   end if

end function get_vdw_rad_cosmo_num


!> Get van-der-Waals radius for species with a given symbol
elemental function get_vdw_rad_bondi_sym(symbol) result(radius)
   !> Element symbol
   character(len=*), intent(in) :: symbol
   !> van-der-Waals radius
   real(wp) :: radius

   radius = get_vdw_rad_bondi(to_number(symbol))

end function get_vdw_rad_bondi_sym


!> Get van-der-Waals radius for species with a given atomic number
elemental function get_vdw_rad_bondi_num(number) result(radius)
   !> Atomic number
   integer, intent(in) :: number
   !> van-der-Waals radius
   real(wp) :: radius

   if (number > 0 .and. number <= size(vdw_rad_bondi, dim=1)) then
      radius = vdw_rad_bondi(number)
   else
      radius = missing
   end if

end function get_vdw_rad_bondi_num


!> Get solvent parameters, in case the solvent is not known the ideal conductor is assumed
pure function get_solvent_data(solvent) result(data)
   !> Name of solvent
   character(len=*), intent(in) :: solvent
   !> Solvent data
   type(solvent_data) :: data

   select case(compact(solvent))
   case default
      data = solvent_data(-1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, "none")
   case("inf")
      data = solvent_data(huge(1.0_wp), 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, "inf")
   case("2methylpyridine")
      data = solvent_data(9.9533_wp, 1.4957_wp, 0.0_wp, 0.58_wp, 47.5_wp, 0.509796_wp, 0.0_wp, "2methylpyridine")
   case("4methyl2pentanone")
      data = solvent_data(12.8871_wp, 1.3962_wp, 0.0_wp, 0.51_wp, 33.83_wp, 0.0_wp, 0.0_wp, "4methyl2pentanone")
   case("aceticacid")
      data = solvent_data(6.2528_wp, 1.372_wp, 0.61_wp, 0.44_wp, 39.01_wp, 0.0_wp, 0.0_wp, "aceticacid")
   case("acetone")
      data = solvent_data(20.493_wp, 1.3588_wp, 0.04_wp, 0.49_wp, 33.77_wp, 0.0_wp, 0.0_wp, "acetone") 
   case("acetonitrile","mecn")
      data = solvent_data(35.6881_wp, 1.3442_wp, 0.07_wp, 0.32_wp, 41.25_wp, 0.0_wp, 0.0_wp, "acetonitrile")
   case("acetophenone")
      data = solvent_data(17.44_wp, 1.5372_wp, 0.0_wp, 0.48_wp, 56.19_wp, 0.444889_wp, 0.0_wp, "acetophenone")
   case("aniline")
      data = solvent_data(6.8882_wp, 1.5863_wp, 0.26_wp, 0.41_wp, 60.62_wp, 0.734449_wp, 0.0_wp, "aniline")
   case("anisole")
      data = solvent_data(4.2247_wp, 1.5174_wp, 0.0_wp, 0.29_wp, 50.52_wp, 0.5625_wp, 0.0_wp, "anisole")
   case("benzaldehyde") ! no Minnisota database entry
      data = solvent_data(17.85_wp,1.5456_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, "benzaldehyde")
   case("benzene")
      data = solvent_data(2.2706_wp, 1.5011_wp, 0.0_wp, 0.14_wp, 40.62_wp, 1_wp, 0.0_wp, "benzene")
   case("benzonitrile")
      data = solvent_data(25.592_wp, 1.5289_wp, 0.0_wp, 0.33_wp, 55.83_wp, 0.5625_wp, 0.0_wp, "benzonitrile")
   case("benzylalcohol")
      data = solvent_data(12.4569_wp, 1.5396_wp, 0.33_wp, 0.56_wp, 52.96_wp, 0.5625_wp, 0.0_wp, "benzylalcohol")
   case("bromobenzene")
      data = solvent_data(5.3954_wp, 1.5597_wp, 0.0_wp, 0.09_wp, 50.72_wp, 0.734449_wp, 0.020449_wp, "bromobenzene")
   case("bromoethane")
      data = solvent_data(9.01_wp, 1.4239_wp, 0.0_wp, 0.12_wp, 34_wp, 0.0_wp, 0.110889_wp, "bromoethane")
   case("bromoform")
      data = solvent_data(4.2488_wp, 1.6005_wp, 0.15_wp, 0.06_wp, 64.58_wp, 0.0_wp, 0.5625_wp)
   case("bromooctane")
      data = solvent_data(5.0244_wp, 1.4524_wp, 0.0_wp, 0.12_wp, 41.28_wp, 0.0_wp, 0.0121_wp, "bromooctane")
   case("butanol")
      data = solvent_data(17.3323_wp, 1.3993_wp, 0.37_wp, 0.48_wp, 35.88_wp, 0.0_wp, 0.0_wp, "butanol")
   case("butanone")
      data = solvent_data(18.2457_wp, 1.3788_wp, 0.0_wp, 0.51_wp, 34.5_wp, 0.0_wp, 0.0_wp, "butanone")
   case("butylacetate")
      data = solvent_data(4.9941_wp, 1.3941_wp, 0.0_wp, 0.45_wp, 35.81_wp, 0.0_wp, 0.0_wp, "butylacetate")
   case("butylbenzene")
      data = solvent_data(2.36_wp, 1.4898_wp, 0.0_wp, 0.15_wp, 41.33_wp, 0.36_wp, 0.0_wp, "butylbenzene")
   case("carbondisulfide","cs2")
      data = solvent_data(2.6105_wp, 1.6319_wp, 0.0_wp, 0.07_wp, 45.45_wp, 0.0_wp, 0.0_wp, "carbondisulfide")
   case("carbontet")
      data = solvent_data(2.228_wp, 1.4601_wp, 0.0_wp, 0.0_wp, 38.04_wp, 0.0_wp, 0.64_wp, "carbontet")
   case("chlorobenzene")
      data = solvent_data(5.6968_wp, 1.5241_wp, 0.0_wp, 0.07_wp, 47.48_wp, 0.734449_wp, 0.020449_wp, "chlorobenzene")
   case("chloroform","chcl3")
      data = solvent_data(4.7113_wp, 1.4459_wp, 0.15_wp, 0.02_wp, 38.39_wp, 0.0_wp, 0.5625_wp, "chloroform")
   case("chlorohexane")
      data = solvent_data(5.9491_wp, 1.4199_wp, 0.0_wp, 0.1_wp, 37.03_wp, 0.0_wp, 0.020449_wp, "chlorohexane")
   case("cyclohexane")
      data = solvent_data(2.0165_wp, 1.4266_wp, 0.0_wp, 0.0_wp, 35.48_wp, 0.0_wp, 0.0_wp, "cyclohexane")
   case("cyclohexanone")
      data = solvent_data(15.6186_wp, 1.4507_wp, 0.0_wp, 0.56_wp, 49.76_wp, 0.0_wp, 0.0_wp, "cyclohexanone")
   case("decalin")
      data = solvent_data(2.196_wp, 1.4528_wp, 0.0_wp, 0.0_wp, 43.82_wp, 0.0_wp, 0.0_wp, "decalin")
   case("decane")
      data = solvent_data(1.9846_wp, 1.4102_wp, 0.0_wp, 0.0_wp, 33.64_wp, 0.0_wp, 0.0_wp, "decane")
   case("decanol")
      data = solvent_data(7.5305_wp, 1.4372_wp, 0.37_wp, 0.48_wp, 41.04_wp, 0.0_wp, 0.0_wp, "decanol")
   case("dibromoethane")
      data = solvent_data(4.9313_wp, 1.5387_wp, 0.1_wp, 0.17_wp, 56.93_wp, 0.0_wp, 0.25_wp, "dibromoethane")
   case("dibutylether")
      data = solvent_data(3.0473_wp, 1.3992_wp, 0.0_wp, 0.45_wp, 35.98_wp, 0.0_wp, 0.0_wp, "dibutylether")
   case("dichloroethane")
      data = solvent_data(10.125_wp, 1.4448_wp, 0.1_wp, 0.11_wp, 45.86_wp, 0.0_wp, 0.25_wp, "dichloroethane")
   case("diethylether","ether")
      data = solvent_data(4.24_wp, 1.3526_wp, 0.0_wp, 0.41_wp, 23.96_wp, 0.0_wp, 0.0_wp, "diethylether")
   case("diisopropylether")
      data = solvent_data(3.38_wp, 1.3679_wp, 0.0_wp, 0.41_wp, 24.86_wp, 0.0_wp, 0.0_wp, "diisopropylether")
   case("dimethylacetamide")
      data = solvent_data(37.7807_wp, 1.438_wp, 0.0_wp, 0.78_wp, 47.62_wp, 0.0_wp, 0.0_wp, "dimethylacetamide")
   case("dimethylformamide","dmf")
      data = solvent_data(37.219_wp, 1.4305_wp, 0.0_wp, 0.74_wp, 49.56_wp, 0.0_wp, 0.0_wp, "dimethylformamide")
   case("dimethylpyridine")
      data = solvent_data(7.1735_wp, 1.4953_wp, 0.0_wp, 0.63_wp, 44.64_wp, 0.390625_wp, 0.0_wp, "dimethylpyridine")
   case("dimethylsulfoxide","dmso")
      data = solvent_data(46.826_wp, 1.417_wp, 0.0_wp, 0.88_wp, 61.78_wp, 0.0_wp, 0.0_wp, "dimethylsulfoxide")
   case("dioxane") ! no Minnisota database entry
      data = solvent_data(2.2099_wp, 1.4224_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, "dioxane")
   case("dodecane")
      data = solvent_data(2.006_wp, 1.4216_wp, 0.0_wp, 0.0_wp, 35.85_wp, 0.0_wp, 0.0_wp, "dodecane")
   case("ethanol")
      data = solvent_data(24.852_wp, 1.3611_wp, 0.37_wp, 0.48_wp, 31.62_wp, 0.0_wp, 0.0_wp, "ethanol")
   case("ethoxybenzene")
      data = solvent_data(4.1797_wp, 1.5076_wp, 0.0_wp, 0.32_wp, 46.65_wp, 0.444889_wp, 0.0_wp, "ethoxybenzene")
   case("ethylacetate")
      data = solvent_data(5.9867_wp, 1.3723_wp, 0.0_wp, 0.45_wp, 33.67_wp, 0.0_wp, 0.0_wp, "ethylacetate")
   case("ethylbenzene")
      data = solvent_data(2.4339_wp, 1.4959_wp, 0.0_wp, 0.15_wp, 41.38_wp, 0.5625_wp, 0.0_wp, "ethylbenzene")
   case("fluorobenzene")
      data = solvent_data(5.42_wp, 1.4684_wp, 0.0_wp, 0.1_wp, 38.37_wp, 0.734449_wp, 0.020449_wp, "fluorobenzene")
   case("fluoroctane")
      data = solvent_data(3.89_wp, 1.3935_wp, 0.0_wp, 0.1_wp, 33.92_wp, 0.0_wp, 0.012321_wp, "fluoroctane")
   case("furan","furane") ! no Minnisota database entry
      data = solvent_data(2.953_wp, 1.4214_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, "furane")
   case("heptane")
      data = solvent_data(1.9113_wp, 1.3878_wp, 0.0_wp, 0.0_wp, 28.28_wp, 0.0_wp, 0.0_wp, "heptane")
   case("heptanol")
      data = solvent_data(11.321_wp, 1.4249_wp, 0.37_wp, 0.48_wp, 38.5_wp, 0.0_wp, 0.0_wp, "heptanol")
   case("hexadecane")
      data = solvent_data(2.0402_wp, 1.4345_wp, 0.0_wp, 0.0_wp, 38.93_wp, 0.0_wp, 0.0_wp, "hexadecane")
   case("hexadecyliodide")
      data = solvent_data(3.5338_wp, 1.4806_wp, 0.0_wp, 0.15_wp, 46.48_wp, 0.0_wp, 0.0_wp, "hexadecyliodide")
   case("hexane","nhexan","n-hexan","nhexane","n-hexane")
      data = solvent_data(1.8819_wp, 1.3749_wp, 0.0_wp, 0.0_wp, 25.7495_wp, 0.0_wp, 0.0_wp, "hexane")
   case("hexanol")
      data = solvent_data(12.5102_wp, 1.4178_wp, 0.37_wp, 0.48_wp, 37.15_wp, 0.0_wp, 0.0_wp, "hexanol")
   case("iodobenzene")
      data = solvent_data(4.547_wp, 1.62_wp, 0.0_wp, 0.12_wp, 55.72_wp, 0.734449_wp, 0.0_wp, "iodobenzene")
   case("isobutanol")
      data = solvent_data(16.7766_wp, 1.3955_wp, 0.37_wp, 0.48_wp, 32.38_wp, 0.0_wp, 0.0_wp, "isobutanol")
   case("isooctane")
      data = solvent_data(1.9358_wp, 1.3915_wp, 0.0_wp, 0.0_wp, 26.38_wp, 0.0_wp, 0.0_wp, "isooctane")
   case("isopropanol")
      data = solvent_data(19.2645_wp, 1.3776_wp, 0.33_wp, 0.56_wp, 30.13_wp, 0.0_wp, 0.0_wp, "isopropanol")
   case("isopropylbenzene")
      data = solvent_data(2.3712_wp, 1.4915_wp, 0.0_wp, 0.16_wp, 39.85_wp, 0.444889_wp, 0.0_wp, "isopropylbenzene")
   case("mcresol")
      data = solvent_data(12.44_wp, 1.5438_wp, 0.57_wp, 0.34_wp, 51.37_wp, 0.5625_wp, 0.0_wp, "mcresol")
   case("mesitylene")
      data = solvent_data(2.265_wp, 1.4994_wp, 0.0_wp, 0.19_wp, 39.65_wp, 0.444889_wp, 0.0_wp, "mesitylene")
   case("methanol")
      data = solvent_data(32.613_wp, 1.3288_wp, 0.43_wp, 0.47_wp, 31.77_wp, 0.0_wp, 0.0_wp, "methanol")
   case("methoxyethanol")
      data = solvent_data(17.2_wp, 1.4024_wp, 0.3_wp, 0.84_wp, 44.39_wp, 0.0_wp, 0.0_wp, "methoxyethanol")
   case("methylenechloride","dichloromethane","dichlormethane","ch2cl2")
      data = solvent_data(8.93_wp, 1.4242_wp, 0.1_wp, 0.05_wp, 39.15_wp, 0.0_wp, 0.444889_wp, "methylenechloride")
   case("methylformamide")
      data = solvent_data(181.5619_wp, 1.4319_wp, 0.4_wp, 0.55_wp, 55.4372_wp, 0.0_wp, 0.0_wp, "methylformamide")
   case("nitrobenzene")
      data = solvent_data(34.8091_wp, 1.5562_wp, 0.0_wp, 0.28_wp, 57.54_wp, 0.444889_wp, 0.0_wp, "nitrobenzene")
   case("nitroethane")
      data = solvent_data(28.2896_wp, 1.3917_wp, 0.02_wp, 0.33_wp, 46.25_wp, 0.0_wp, 0.0_wp, "nitroethane")
   case("nitromethane")
      data = solvent_data(36.5623_wp, 1.3817_wp, 0.06_wp, 0.31_wp, 52.58_wp, 0.0_wp, 0.0_wp, "nitromethane")
   case("nonane")
      data = solvent_data(1.9605_wp, 1.4054_wp, 0.0_wp, 0.0_wp, 32.21_wp, 0.0_wp, 0.0_wp, "nonane")
   case("nonanol")
      data = solvent_data(8.5991_wp, 1.4333_wp, 0.37_wp, 0.48_wp, 40.14_wp, 0.0_wp, 0.0_wp, "nonanol")
   case("octane")
      data = solvent_data(1.9406_wp, 1.3974_wp, 0.0_wp, 0.0_wp, 30.4273_wp, 0.0_wp, 0.0_wp, "octane")
   case("octanol")
      data = solvent_data(9.8629_wp, 1.4295_wp, 0.37_wp, 0.48_wp, 39.01_wp, 0.0_wp, 0.0_wp, "octanol")
   case("odichlorobenzene")
      data = solvent_data(9.9949_wp, 1.5515_wp, 0.0_wp, 0.04_wp, 52.72_wp, 0.5625_wp, 0.0625_wp, "odichlorobenzene")
   case("onitrotoluene")
      data = solvent_data(25.6692_wp, 1.545_wp, 0.0_wp, 0.27_wp, 59.12_wp, 0.36_wp, 0.0_wp, "onitrotoluene")
   case("pentadecane")
      data = solvent_data(2.0333_wp, 1.4315_wp, 0.0_wp, 0.0_wp, 38.34_wp, 0.0_wp, 0.0_wp, "pentadecane")
   case("pentane")
      data = solvent_data(1.8371_wp, 1.3575_wp, 0.0_wp, 0.0_wp, 22.3_wp, 0.0_wp, 0.0_wp, "pentane")
   case("pentanol")
      data = solvent_data(15.13_wp, 1.4101_wp, 0.37_wp, 0.48_wp, 36.5_wp, 0.0_wp, 0.0_wp, "pentanol")
   case("perfluorobenzene")
      data = solvent_data(2.029_wp, 1.3777_wp, 0.0_wp, 0.0_wp, 31.74_wp, 0.25_wp, 0.25_wp, "perfluorobenzene")
   case("phenol") ! no Minnisota database entry
      data = solvent_data(12.4_wp, 1.5408_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, "phenol")
   case("phenylether")
      data = solvent_data(3.73_wp, 1.5787_wp, 0.0_wp, 0.2_wp, 38.5_wp, 0.851929_wp, 0.0_wp, "phenylether")
   case("propanol")
      data = solvent_data(20.5237_wp, 1.385_wp, 0.37_wp, 0.48_wp, 33.57_wp, 0.0_wp, 0.0_wp, "propanol")
   case("pyridine")
      data = solvent_data(12.9776_wp, 1.5095_wp, 0.0_wp, 0.52_wp, 52.62_wp, 0.693889_wp, 0.0_wp, "pyridine")
   case("secbutanol")
      data = solvent_data(15.9436_wp, 1.3978_wp, 0.33_wp, 0.56_wp, 32.44_wp, 0.0_wp, 0.0_wp, "secbutanol")
   case("secbutylbenzene")
      data = solvent_data(2.3446_wp, 1.4895_wp, 0.0_wp, 0.16_wp, 40.35_wp, 0.36_wp, 0.0_wp, "secbutylbenzene")
   case("tbutylbenzene")
      data = solvent_data(2.3447_wp, 1.4927_wp, 0.0_wp, 0.16_wp, 39.78_wp, 0.36_wp, 0.0_wp, "tbutylbenzene")
   case("tetrachloroethene")
      data = solvent_data(2.268_wp, 1.5053_wp, 0.0_wp, 0.0_wp, 45.19_wp, 0.0_wp, 0.444889_wp, "tetrachloroethene")
   case("tetrahydrofuran","thf")
      data = solvent_data(7.4257_wp, 1.405_wp, 0.0_wp, 0.48_wp, 39.44_wp, 0.0_wp, 0.0_wp, "tetrahydrofuran")
   case("tetrahydrothiophenedioxide")
      data = solvent_data(43.9622_wp, 1.4833_wp, 0.0_wp, 0.88_wp, 87.49_wp, 0.0_wp, 0.0_wp, "tetrahydrothiophenedioxide")
   case("tetralin")
      data = solvent_data(2.771_wp, 1.5413_wp, 0.0_wp, 0.19_wp, 47.74_wp, 0.36_wp, 0.0_wp, "tetralin")
   case("toluene")
      data = solvent_data(2.3741_wp, 1.4961_wp, 0.0_wp, 0.14_wp, 40.2_wp, 0.734449_wp, 0.0_wp, "toluene")
   case("tributylphosphate")
      data = solvent_data(8.1781_wp, 1.4224_wp, 0.0_wp, 1.21_wp, 27.55_wp, 0.0_wp, 0.0_wp, "tributylphosphate")
   case("triethylamine")
      data = solvent_data(2.3832_wp, 1.401_wp, 0.0_wp, 0.79_wp, 29.1_wp, 0.0_wp, 0.0_wp, "triethylamine")
   case("trimethylbenzene")
      data = solvent_data(2.3653_wp, 1.5048_wp, 0.0_wp, 0.19_wp, 42.03_wp, 0.444889_wp, 0.0_wp, "trimethylbenzene")
   case("undecane")
      data = solvent_data(1.991_wp, 1.4398_wp, 0.0_wp, 0.0_wp, 34.85_wp, 0.0_wp, 0.0_wp, "undecane")
   case("water","h2o")
      data = solvent_data(78.36_wp, 1.3328_wp, 0.82_wp, 0.35_wp, 103.62_wp, 0.0_wp, 0.0_wp, "water")
   case("woctanol") ! no Minnisota database entry
      data = solvent_data(9.8629_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, "woctanol")
   case("xylene")
      data = solvent_data(2.3879_wp, 1.4995_wp, 0.0_wp, 0.16_wp, 41.38_wp, 0.5625_wp, 0.0_wp, "xylene")
   end select
end function get_solvent_data

end module tblite_solvation_data
