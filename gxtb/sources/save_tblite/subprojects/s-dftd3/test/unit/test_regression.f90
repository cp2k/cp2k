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

module test_regression
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use dftd3
   implicit none
   private

   public :: collect_regression

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_regression(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("dftbplus#871", test_dftbplus_871) &
      & ]

end subroutine collect_regression


!> Reported at https://github.com/dftbplus/dftbplus/issue/871
subroutine test_dftbplus_871(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param = rational_damping_param( &
      & s6=1.00_wp, s8=2.34_wp, a1=6.30_wp, a2=5.00_wp, s9=0.00_wp, alp=14.0_wp)
   type(d3_model) :: d3
   real(wp) :: energy, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)

   call structure_dftbplus_871(mol)
   allocate(gradient(3, mol%nat))
   call new_d3_model(d3, mol)
   call get_dispersion(mol, d3, param, realspace_cutoff(), energy, gradient, sigma)

   if (is_exceptional(norm2(gradient))) then
      call test_failed(error, "Exceptional value propagated to the gradient")
      return
   end if

   call check(error, energy, -8.9964956516666404E-5_wp, thr=thr)
   if (allocated(error)) return
   call check(error, norm2(gradient), 8.9844177273972976E-9_wp, thr=thr)

end subroutine test_dftbplus_871


!> Check whether we are dealing with an exceptional value, NaN or Inf
elemental function is_exceptional(val)
   use ieee_arithmetic, only : ieee_is_nan
   real(wp), intent(in) :: val
   logical :: is_exceptional
   is_exceptional = ieee_is_nan(val) .or. abs(val) > huge(val)
end function is_exceptional


subroutine structure_dftbplus_871(self)
   type(structure_type), intent(out) :: self
   integer, parameter :: nat = 90
   character(len=*), parameter :: sym(nat) = [character(len=4)::&
      & "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", &
      & "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", &
      & "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", &
      & "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", &
      & "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", &
      & "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", "Ta", &
      & "Ta", "Ta", "Ta", "Ta", "Ta", "Ta"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      &  2.19308288262168_wp,  8.11634393608512_wp, 32.35990015739237_wp, &
      & -2.30657594039162_wp, 13.98031135475578_wp, 39.25405756706034_wp, &
      &  6.69274265930930_wp,  2.25237699425163_wp, 25.46574084037577_wp, &
      &  9.46827252075950_wp, 24.66680713570060_wp, 30.44026438468729_wp, &
      & 14.28964169190208_wp,  4.14292021667899_wp, 36.22537485832011_wp, &
      & 13.96793110535422_wp, 18.80283924019279_wp, 23.54610506767070_wp, &
      &  3.99725945160666_wp, 23.81505390083732_wp, 31.29898325675761_wp, &
      &  8.81862957642356_wp,  3.29116745865287_wp, 37.08409754508769_wp, &
      &  8.49691803620139_wp, 17.95108600532951_wp, 24.40482584708964_wp, &
      &  1.31219442055503_wp, 23.86167331612052_wp, 25.80180041068827_wp, &
      &  6.13356430695335_wp,  3.33778639709892_wp, 31.58690897697245_wp, &
      &  1.63390524552146_wp,  9.20175357735099_wp, 38.48107020133769_wp, &
      &  2.75415952370445_wp, 13.39193530952873_wp, 35.64818827384746_wp, &
      & 11.75347740814964_wp,  1.66400094902458_wp, 21.85987345451152_wp, &
      &  7.25381882355491_wp,  7.52796836769523_wp, 28.75403086417948_wp, &
      & -2.44976942374429_wp,  7.58241315758171_wp, 29.78376070731913_wp, &
      & 14.90544636414329_wp, 13.44638057625236_wp, 36.12339655631816_wp, &
      &  2.04988963768760_wp,  1.71844550049247_wp, 22.88960520499980_wp, &
      &  6.54104073212425_wp, 17.28986164009513_wp, 38.50863901847637_wp, &
      & 15.54035885498802_wp,  5.56192680275382_wp, 24.72032419914042_wp, &
      & 11.04070027039329_wp, 11.42589374458732_wp, 31.61448160880839_wp, &
      & -2.47468643023690_wp, 19.96374698555412_wp, 22.94135157340800_wp, &
      & -6.97434584929665_wp, 25.82771488106193_wp, 29.83550898307596_wp, &
      & -2.15297584368905_wp,  5.30382820045890_wp, 35.62061945670879_wp, &
      & 16.45295460388938_wp, 12.58680244362296_wp, 30.45520464652812_wp, &
      & 11.95329697296897_wp, 18.45076938545646_wp, 37.34936014884746_wp, &
      & -0.90225934743008_wp,  6.72283454811515_wp, 24.11556879752909_wp, &
      & -1.97328333689889_wp, 10.88073058998527_wp, 34.69864908927715_wp, &
      & -2.29499404265603_wp, 25.54064746773185_wp, 22.01938502067362_wp, &
      &  2.52637584374229_wp,  5.01676269447746_wp, 27.80449549430644_wp, &
      &  0.18181331322471_wp, 19.54803272163810_wp, 28.73308054679667_wp, &
      & -4.31784574820718_wp, 25.41200061714592_wp, 35.62723795646464_wp, &
      &  4.68147261307517_wp, 13.68406577980461_wp, 21.83892504447734_wp, &
      &  4.28300650284568_wp, 15.69481941139640_wp, 27.27791468375956_wp, &
      & -0.21665261819085_wp, 21.55878826057853_wp, 34.17207590812480_wp, &
      &  8.78266556427756_wp,  9.83085246956291_wp, 20.38375918144023_wp, &
      &  8.48231060669700_wp, 20.29385753548088_wp, 33.51624933952128_wp, &
      & 17.48162968323509_wp,  8.56592269813957_wp, 19.72793452018535_wp, &
      & 12.98196919129173_wp, 14.42989059364738_wp, 26.62209192985332_wp, &
      &  8.83651383087913_wp, 16.24142166054191_wp, 29.78496424430643_wp, &
      &  4.33685476944724_wp, 22.10538860237541_wp, 36.67912165397440_wp, &
      & 13.33617336914817_wp, 10.37745376503410_wp, 22.89080683463847_wp, &
      &  6.37293036148826_wp, 23.27329631721916_wp, 22.19593302482401_wp, &
      & 11.19430000946799_wp,  2.74941035187187_wp, 27.98104540580546_wp, &
      &  6.69464094803611_wp,  8.61337848579826_wp, 34.87520090812480_wp, &
      &  3.88454086945335_wp, 17.70557399666252_wp, 32.71691004508769_wp, &
      & -0.61511843039712_wp, 23.56954189217033_wp, 39.61106745475566_wp, &
      &  8.38419945404807_wp, 11.84160705482902_wp, 25.82275072807109_wp, &
      &  6.24731667206565_wp, 20.88277431404533_wp, 27.85155550712382_wp, &
      &  1.74765737221519_wp, 26.74674030220451_wp, 34.74571100944316_wp, &
      & 10.74697525666038_wp, 15.01880546486320_wp, 20.95740000480448_wp, &
      &  0.06909485028068_wp, 13.43855472481193_wp, 30.15100352042948_wp, &
      & -4.43056433036050_wp, 19.30252262031975_wp, 37.04516093009746_wp, &
      &  4.56875355408469_wp,  7.57458730614128_wp, 23.25684801811015_wp, &
      & -6.80528609587869_wp, 23.02478023445549_wp, 21.21554247611796_wp, &
      & 19.87095768616478_wp,  2.50089355385246_wp, 26.44612948173319_wp, &
      & 15.37129910157005_wp,  8.36486144936027_wp, 33.34028879874980_wp, &
      &  8.79064209625999_wp, 14.35817428505363_wp, 35.06190745109355_wp, &
      & 17.78996021912376_wp,  2.63024016296806_wp, 21.27359453910624_wp, &
      & 13.29030163452903_wp,  8.49420782005729_wp, 28.16775004142558_wp, &
      & -2.18096260859689_wp, 16.37083621895255_wp, 33.59843317741190_wp, &
      &  6.81835539505760_wp,  4.64290042793693_wp, 19.81011835807597_wp, &
      &  2.31869657204429_wp, 10.50686832344474_wp, 26.70427576774394_wp, &
      & -0.27050064637383_wp, 15.14821906959953_wp, 24.77086703056132_wp, &
      & -4.77015964820107_wp, 21.01218601143303_wp, 31.66502253288066_wp, &
      & -9.26981894805154_wp, 26.87615390694084_wp, 38.55918375724590_wp, &
      & -4.72428839041909_wp, 22.89543338692130_wp, 26.38808123344218_wp, &
      &  0.09708173439780_wp,  2.37154754078330_wp, 32.17319361442363_wp, &
      & 17.45229656861106_wp,  8.23551507866325_wp, 38.51282755607401_wp, &
      &  5.53267700836936_wp, 12.26058574593009_wp, 31.01460520499980_wp, &
      &  1.03301806614677_wp, 18.12455364143791_wp, 37.90876261466777_wp, &
      & 10.03233559296409_wp,  6.39661928093376_wp, 24.12044970268046_wp, &
      & 12.31745560334007_wp, 20.31532855904045_wp, 28.59174792045390_wp, &
      &  7.81779701874534_wp, 26.17929645454826_wp, 35.48590342277323_wp, &
      & -5.03775899245461_wp, 14.45136257088127_wp, 22.25211207145487_wp, &
      &  6.03963692353050_wp, 26.37287898933830_wp, 26.75133768790996_wp, &
      & 10.86100704834739_wp,  5.84899230873527_wp, 32.53644816154276_wp, &
      & 10.53929646179954_wp, 20.50891109383048_wp, 19.85718027824199_wp, &
      & -3.75144259764870_wp, 11.07431026375236_wp, 25.96408335441386_wp, &
      & 13.60377247498313_wp, 16.93827815926017_wp, 32.30371729606425_wp, &
      &  9.10411293671409_wp, 22.80224605476799_wp, 39.19787852042948_wp, &
      & -4.41595666243753_wp, 16.95975109016838_wp, 27.93373934501444_wp, &
      & 12.93925888703147_wp, 22.82371707832755_wp, 34.27337519401347_wp, &
      &  0.08370233938972_wp, 11.09578319466057_wp, 21.03958193534648_wp, &
      &  3.03333575890342_wp, 19.12905498421134_wp, 23.54122606986796_wp, &
      & -1.46632324292382_wp, 24.99302287971916_wp, 30.43538347953593_wp, &
      &  3.35504658386985_wp,  4.46913548386039_wp, 36.22049204582011_wp, &
      & 11.01578266785423_wp, 23.80722804939689_wp, 24.77207056754862_wp, &
      & 15.83715374634544_wp,  3.28334208404960_wp, 30.55718294853007_wp, &
      & 11.33749325440208_wp,  9.14730926430168_wp, 37.45133845084941_wp],&
      & shape(xyz))
   real(wp), parameter :: lattice(3, 3) = reshape([&
      & 21.86190700000000_wp,  0.00000000000000_wp,  0.00000000000000_wp, &
      & -9.34616046617651_wp, 26.40226595548828_wp,  0.00000000000000_wp, &
      & -1.50349355759049_wp,  1.95765836507584_wp, 59.22338117110787_wp],&
      & shape(lattice))
   call new(self, sym, xyz, lattice=lattice)
end subroutine structure_dftbplus_871


end module test_regression
