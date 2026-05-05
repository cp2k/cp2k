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

module test_gfn1_xtb
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_context_type, only : context_type
   use tblite_lapack_solver, only : lapack_solver, lapack_algorithm
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_gfn1_xtb

   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 100_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 3.166808578545117e-06_wp

contains


!> Collect all exported unit tests
subroutine collect_gfn1_xtb(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-atom", test_e_pse), &
      new_unittest("energy-atom-cation", test_e_pse_cation), &
      new_unittest("energy-atom-anion", test_e_pse_anion), &
      new_unittest("energy-mol", test_e_mb01), &
      new_unittest("gradient-mol", test_g_mb02), &
      new_unittest("numgrad-mol", test_g_mb03), &
      !new_unittest("virial-mol", test_s_mb03), &
      new_unittest("error-uhf", test_error_mb01) &
      ]

end subroutine collect_gfn1_xtb


subroutine numdiff_grad(ctx, mol, calc, wfn, numgrad)
   type(context_type), intent(inout) :: ctx
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(in) :: wfn
   real(wp), intent(out) :: numgrad(:, :)

   integer :: iat, ic
   real(wp) :: er, el
   type(structure_type) :: moli
   type(wavefunction_type) :: wfni

   ! Step size is square root of SCF convergence threshold
   ! (SCF energy convergence threshold is 1e-6 * acc)
   real(wp), parameter :: step = 1.0e-4_wp

   do iat = 1, mol%nat
      do ic = 1, 3
         moli = mol
         wfni = wfn
         moli%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call xtb_singlepoint(ctx, moli, calc, wfni, acc, er, verbosity=0)

         moli = mol
         wfni = wfn
         moli%xyz(ic, iat) = mol%xyz(ic, iat) - step
         call xtb_singlepoint(ctx, moli, calc, wfni, acc, el, verbosity=0)

         numgrad(ic, iat) = 0.5_wp*(er - el)/step
      end do
   end do
end subroutine numdiff_grad


subroutine numdiff_sigma(ctx, mol, calc, wfn, numsigma)
   type(context_type), intent(inout) :: ctx
   type(structure_type), intent(in) :: mol
   type(xtb_calculator), intent(in) :: calc
   type(wavefunction_type), intent(in) :: wfn
   real(wp), intent(out) :: numsigma(:, :)

   integer :: ic, jc
   real(wp) :: er, el, eps(3, 3)
   type(structure_type) :: moli
   type(wavefunction_type) :: wfni
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))

   eps(:, :) = unity
   do ic = 1, 3
      do jc = 1, 3
         moli = mol
         wfni = wfn
         eps(jc, ic) = eps(jc, ic) + step
         moli%xyz(:, :) = matmul(eps, mol%xyz)
         if (any(mol%periodic)) moli%lattice(:, :) = matmul(eps, mol%lattice)
         call xtb_singlepoint(ctx, moli, calc, wfni, acc, er, verbosity=0)

         moli = mol
         wfni = wfn
         eps(jc, ic) = eps(jc, ic) - step
         moli%xyz(:, :) = matmul(eps, mol%xyz)
         if (any(mol%periodic)) moli%lattice(:, :) = matmul(eps, mol%lattice)
         call xtb_singlepoint(ctx, moli, calc, wfni, acc, el, verbosity=0)

         numsigma(jc, ic) = 0.5_wp*(er - el)/step
      end do
   end do
   numsigma = (numsigma + transpose(numsigma)) * 0.5_wp
end subroutine numdiff_sigma


subroutine test_e_pse(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp) :: energy
   real(wp), parameter :: xyz(3, 1) = 0.0_wp
   integer, parameter :: uhf(86) = [&
      & 1, 0, &
      & 1, 0, 1, 0, 1, 0, 1, 0, &
      & 1, 0, 1, 0, 1, 0, 1, 0, &
      & 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, &
      & 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, &
      & 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, &
      & 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
   real(wp), parameter :: ref(86) = [&
      &-0.4014294744618301_wp,-1.6258646654198721_wp,-0.2671714426384922_wp, &
      &-0.7012869049089436_wp,-1.1998696279038878_wp,-1.7411359557542052_wp, &
      &-2.8988862104065958_wp,-4.3526523408648030_wp,-4.9969448424630265_wp, &
      &-6.3190030086322286_wp,-0.1733674706866860_wp,-0.7328492086440124_wp, &
      &-1.0775966895776903_wp,-1.6252631326184159_wp,-2.4250620380497385_wp, &
      &-3.5359207433156339_wp,-4.1663539440166675_wp,-6.1469954466724239_wp, &
      &-0.2137179709637750_wp,-0.5864589388109813_wp,-0.9739522227703182_wp, &
      &-1.1487975500262675_wp,-1.4032495694507467_wp,-1.6258320924723875_wp, &
      &-2.0432134739535242_wp,-2.8490837602318742_wp,-3.4039476610517951_wp, &
      &-3.9783299398205854_wp,-4.3908896126848882_wp,-0.8278490035919727_wp, &
      &-1.1559605434617155_wp,-1.5172654693894492_wp,-2.2030246081688953_wp, &
      &-3.2227173551239399_wp,-3.8180395534595282_wp,-4.8421263290307968_wp, &
      &-0.2799176883500828_wp,-0.5027433177275922_wp,-0.8559774488167702_wp, &
      &-1.0015251619993044_wp,-1.7072155954305537_wp,-1.7786275624264398_wp, &
      &-2.2160358545342351_wp,-3.2181664696662833_wp,-3.9147305338504088_wp, &
      &-4.4068901473558819_wp,-3.7633756826766938_wp,-0.8892760127989151_wp, &
      &-1.3780364568041852_wp,-2.3412464596133091_wp,-2.3522450904243550_wp, &
      &-3.6750286938689936_wp,-3.8857572752278440_wp,-4.2321989780523870_wp, &
      &-0.2330413901515095_wp,-0.4742595783451311_wp,-0.6959185846380349_wp, &
      &-0.6485304192347468_wp,-0.6450212630347002_wp,-0.6415120701014334_wp, &
      &-0.6380028404220172_wp,-0.6344936474925504_wp,-0.6309844913125303_wp, &
      &-0.6274752983832106_wp,-0.6239660687045753_wp,-0.6204569125245851_wp, &
      &-0.6169477195952713_wp,-0.6134385634152814_wp,-0.6099292969873201_wp, &
      &-0.6064201408073303_wp,-0.6029109478780164_wp,-0.7596719486827229_wp, &
      &-1.6934927181631920_wp,-1.9903305781546048_wp,-2.0068018142203417_wp, &
      &-3.1552379430508148_wp,-3.9539731772481859_wp,-4.2489748260259601_wp, &
      &-3.9019226424065825_wp,-0.9152519783258882_wp,-1.1137676526000495_wp, &
      &-1.4989523256308206_wp,-2.0838915690368061_wp,-3.0676044864514225_wp, &
      &-3.4437503367483560_wp,-3.6081562853045592_wp]

   integer :: izp

   do izp = 1, 86
      if (izp == 25) cycle
      call new(mol, [izp], xyz, uhf=uhf(izp))
      call new_gfn1_calculator(calc, mol, error, accuracy=acc)
      if (allocated(error)) return
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
         & 1, calc%default_etemp * kt)

      energy = 0.0_wp
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)
      call check(error, .not.ctx%failed(), &
         & message="SCF does not converge for "//trim(mol%sym(1)))
      if (allocated(error)) exit
      call check(error, energy, ref(izp), thr=thr2, &
         & message="Atomic energy does not match for "//trim(mol%sym(1)))
      if (allocated(error)) exit
   end do

end subroutine test_e_pse


subroutine test_e_pse_cation(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp) :: energy
   real(wp), parameter :: xyz(3, 1) = 0.0_wp
   integer, parameter :: uhf(86) = [&
      & 0, 1, &
      & 0, 1, 0, 1, 0, 1, 0, 1, &
      & 0, 1, 0, 1, 0, 1, 0, 1, &
      & 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, &
      & 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, &
      & 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
      & 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
   real(wp), parameter :: ref(86) = [&
      & 0.2350495000000000_wp,-0.0422428129239876_wp, 0.1369166666666667_wp, &
      &-0.1836139857878052_wp,-0.7021610492209647_wp,-1.1060742863488979_wp, &
      &-2.1782414865973068_wp,-3.3929027848641802_wp,-3.9486100517085290_wp, &
      &-5.4089550104332629_wp, 0.1229540000000000_wp,-0.1526824376553396_wp, &
      &-0.7929165898783470_wp,-1.2064302987580611_wp,-1.7319786163576927_wp, &
      &-2.7725673355423974_wp,-3.5179743056411774_wp,-5.2547187063374068_wp, &
      & 0.0816097666666667_wp,-0.1889038694054907_wp,-0.5164286712433095_wp, &
      &-0.7054229863305306_wp,-0.8953744943379764_wp,-1.1145548860858483_wp, &
      &-1.6781450836464422_wp,-2.0768800363599724_wp,-2.8098681418659432_wp, &
      &-3.3918338006786093_wp,-3.9272193399045281_wp,-0.2625283351293197_wp, &
      &-0.8517810578775546_wp,-1.0787522669740286_wp,-1.5017772126895934_wp, &
      &-2.5547725893532198_wp,-3.0519213538372423_wp,-3.9953892543401519_wp, &
      & 0.0878675000000000_wp,-0.1466078255304628_wp,-0.4937299864801741_wp, &
      &-0.6513197157489261_wp,-1.2734962578300519_wp,-1.3411315800579504_wp, &
      &-1.7441817291584463_wp,-2.6949495714667360_wp,-3.3840208304086556_wp, &
      &-3.8665790293167479_wp,-3.2936887930479184_wp,-0.3338044397327909_wp, &
      &-1.1351439870281808_wp,-1.9118351274462178_wp,-1.8456782193519603_wp, &
      &-3.0738678691656975_wp,-3.2276435086854969_wp,-3.5173994843193732_wp, &
      & 0.0925550000000000_wp,-0.1182202891725656_wp,-0.3383862439154567_wp, &
      &-0.3380276085842653_wp,-0.3364513307550568_wp,-0.3350141427496482_wp, &
      &-0.3337162359090419_wp,-0.3325572206289337_wp,-0.3315373115184752_wp, &
      &-0.3306564922318164_wp,-0.3299149902654105_wp,-0.3293123596652018_wp, &
      &-0.3288488188887932_wp,-0.3285243679361848_wp,-0.3283389200405289_wp, &
      &-0.3282929953096204_wp,-0.3283858232519615_wp,-0.4538801172959733_wp, &
      &-1.2190319743563898_wp,-1.5399384468018016_wp,-1.5831499553557236_wp, &
      &-2.5862107373895995_wp,-3.4087534399601136_wp,-3.6246625551978040_wp, &
      &-3.2948988109803317_wp,-0.2990523224962774_wp,-0.8284450281457056_wp, &
      &-1.1284674720252441_wp,-1.6724837041091658_wp,-2.4123355716184571_wp, &
      &-2.6886076809857715_wp,-2.9576573558829384_wp]

   integer :: izp

   do izp = 1, 86
      if (izp == 79) cycle  ! SCF does not converge for gold
      call new(mol, [izp], xyz, uhf=uhf(izp), charge=1.0_wp)
      call new_gfn1_calculator(calc, mol, error, accuracy=acc)
      if (allocated(error)) return
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
         & 1, calc%default_etemp * kt)

      energy = 0.0_wp
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)
      call check(error, .not.ctx%failed(), &
         & message="SCF does not converge for "//trim(mol%sym(1)))
      if (allocated(error)) exit
      call check(error, energy, ref(izp), thr=thr2, &
         & message="Atomic energy does not match for "//trim(mol%sym(1)))
      if (allocated(error)) exit
   end do

end subroutine test_e_pse_cation


subroutine test_e_pse_anion(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp) :: energy
   real(wp), parameter :: xyz(3, 1) = 0.0_wp
   integer, parameter :: uhf(86) = [&
      & 0, 1, &
      & 0, 1, 0, 1, 0, 1, 0, 1, &
      & 0, 1, 0, 1, 0, 1, 0, 1, &
      & 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, &
      & 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, &
      & 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, &
      & 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
   real(wp), parameter :: ref(86) = [&
      &-0.5678094489236601_wp, 0.0000000000000000_wp,-0.4659175519436511_wp, &
      &-0.7944690166480025_wp,-1.3514075648859643_wp,-1.9170116002810331_wp, &
      &-3.1284233130873655_wp,-4.7053860329689856_wp,-5.3222970343111768_wp, &
      &-6.0135919541021741_wp,-0.3037809414057536_wp,-0.6889153966301176_wp, &
      &-1.1294568900412956_wp,-1.8636984933817820_wp,-2.5218710230846475_wp, &
      &-3.7615667932865060_wp,-4.5279482362587169_wp,-5.9532003542728642_wp, &
      &-0.3946877084813858_wp,-0.6601217419378675_wp,-1.1181158489120291_wp, &
      &-1.3177666092815929_wp,-1.2384479883964847_wp,-1.6487377642904479_wp, &
      &-0.5954698324087719_wp,-2.6284270918459773_wp,-3.5618386589525119_wp, &
      &-4.1031863005505338_wp,-4.5943678854652479_wp,-0.9409945569463229_wp, &
      &-1.3059077459899349_wp,-1.7645198843278482_wp,-2.1964765797004659_wp, &
      &-3.4039351935495046_wp,-3.9571135364820273_wp,-4.6240117892188310_wp, &
      &-0.5719678767001657_wp,-0.6583170517916040_wp,-1.0130429561259551_wp, &
      &-1.1770073221613631_wp,-1.9625368412303521_wp,-2.0338624073866263_wp, &
      &-2.1752623138381653_wp,-3.3802910855544024_wp,-4.0402591326658541_wp, &
      &-4.4158570177448571_wp,-3.7600845723053635_wp,-1.0180004822365232_wp, &
      &-1.4888020156580570_wp,-2.4537671300397492_wp,-2.3767986294370997_wp, &
      &-3.8685840856051446_wp,-4.0435926774613034_wp,-3.9095090580155762_wp, &
      &-0.4735277804823783_wp,-0.6639055123686469_wp,-0.8747426483751759_wp, &
      &-0.7930037109846064_wp,-0.7859031407166542_wp,-0.7774600374514302_wp, &
      &-0.7692826709429081_wp,-0.7615788166990638_wp,-0.7541577457962461_wp, &
      &-0.7468810501730819_wp,-0.7397437693038101_wp,-0.7327453598812527_wp, &
      &-0.7258859667911156_wp,-0.7191657370235338_wp,-0.7125844368143500_wp, &
      &-0.7061427332686532_wp,-0.6998397088976513_wp,-0.8340395816028910_wp, &
      &-1.9776974757621952_wp,-2.1930497428504525_wp,-2.0018467092139165_wp, &
      &-3.3776836100161272_wp,-4.0976008074736150_wp,-4.3290532643937203_wp, &
      &-4.0370809224168633_wp,-0.9642956138466999_wp,-1.2716085097941741_wp, &
      &-1.6601446321465445_wp,-2.2260645738865312_wp,-3.0614995810384213_wp, &
      &-3.3648843601107044_wp,-3.2192830086977700_wp]

   integer :: izp

   do izp = 1, 86
      if (izp == 2) cycle  ! Helium doesn't have enough orbitals for negative charge
      if (any(izp == [21, 22, 23, 25, 43, 57, 58, 59])) cycle  ! not converging
      call new(mol, [izp], xyz, uhf=uhf(izp), charge=-1.0_wp)
      call new_gfn1_calculator(calc, mol, error, accuracy=acc)
      if (allocated(error)) return
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
         & 1, calc%default_etemp * kt)

      energy = 0.0_wp
      call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)
      call check(error, .not.ctx%failed(), &
         & message="SCF does not converge for "//trim(mol%sym(1)))
      if (allocated(error)) exit
      call check(error, energy, ref(izp), thr=thr2, &
         & message="Atomic energy does not match for "//trim(mol%sym(1)))
      if (allocated(error)) exit
   end do

end subroutine test_e_pse_anion


subroutine test_e_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp) :: energy
   real(wp), parameter :: ref = -33.040345103604_wp  ! value calculated by xtb

   call get_structure(mol, "MB16-43", "01")

   energy = 0.0_wp

   call new_gfn1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

   call check(error, energy, ref, thr=1e-7_wp)

end subroutine test_e_mb01


subroutine test_g_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   real(wp), parameter :: ref(3, 16) = reshape([&  ! calculated by xtb
      & 3.0701068406331E-03_wp,  -2.8687993171440E-04_wp,  -2.7410991187148E-03_wp, &
      &-8.8629107012253E-03_wp,   1.0764495017149E-02_wp,  -2.2851740646278E-03_wp, &
      &-2.2557834233232E-03_wp,   9.7846226010149E-03_wp,   8.2851189762489E-03_wp, &
      &-7.7313039027486E-03_wp,   5.1272904824851E-03_wp,  -1.1338292273413E-02_wp, &
      &-5.6136126168979E-03_wp,  -1.4277687517607E-02_wp,  -8.7142306125955E-03_wp, &
      &-2.8860397448112E-04_wp,  -9.7634288764352E-03_wp,   7.6436790299120E-03_wp, &
      & 1.8264923404347E-03_wp,  -6.7638299477356E-04_wp,   1.4564513464460E-03_wp, &
      & 7.6140197746102E-03_wp,  -4.6761005573138E-03_wp,   1.7232054352028E-04_wp, &
      &-1.4988170602255E-02_wp,  -3.2937246100605E-03_wp,   2.8792223976237E-02_wp, &
      & 8.6150837651184E-03_wp,   1.5081696039773E-03_wp,   2.7847814325776E-03_wp, &
      & 1.1723867162304E-02_wp,  -1.2461660002080E-03_wp,  -6.5732067425589E-04_wp, &
      &-7.0967446809462E-03_wp,  -3.8359040408867E-03_wp,  -4.6752045418470E-03_wp, &
      & 1.0347555714522E-02_wp,   7.0999014051259E-03_wp,   7.5036479095391E-03_wp, &
      & 4.8130307607921E-03_wp,  -7.2967587464359E-03_wp,   2.1987757421243E-03_wp, &
      &-6.2043176813360E-03_wp,  -8.3018303879663E-03_wp,  -6.1290891860884E-03_wp, &
      & 5.0312912247987E-03_wp,   1.9370384553649E-02_wp,  -2.2296588485063E-02_wp],&
      & shape(ref))

   ctx%solver = lapack_solver(lapack_algorithm%gvr)
   call get_structure(mol, "MB16-43", "02")

   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_gfn1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, 0)

   if (any(abs(gradient - ref) > thr2)) then
      call test_failed(error, "Gradient of energy does not match")
      print'(3es21.14)', gradient
      print'("---")'
      print'(3es21.14)', ref
      print'("---")'
      print'(3es21.14)', gradient-ref
   end if

end subroutine test_g_mb02


subroutine test_g_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), numgrad(:, :), sigma(:, :)

   ctx%solver = lapack_solver(lapack_algorithm%gvr)
   call get_structure(mol, "MB16-43", "03")

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   numgrad(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_gfn1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, 0)

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)
   call numdiff_grad(ctx, mol, calc, wfn, numgrad)

   if (any(abs(gradient - numgrad) > thr2)) then
      call test_failed(error, "Gradient of energy does not match (numerical)")
      print'(3es21.14)', gradient
      print'("---")'
      print'(3es21.14)', numgrad
      print'("---")'
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_g_mb03


subroutine test_s_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp) :: energy
   real(wp), allocatable :: gradient(:, :), sigma(:, :), numsigma(:, :)

   call get_structure(mol, "MB16-43", "03")

   allocate(gradient(3, mol%nat), sigma(3, 3), numsigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_gfn1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, 0)

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)
   call numdiff_sigma(ctx, mol, calc, wfn, numsigma)

   if (any(abs(sigma - numsigma) > thr2)) then
      call test_failed(error, "Strain derivatives do not match")
      print'(3es21.14)', sigma
      print'("---")'
      print'(3es21.14)', numsigma
      print'("---")'
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_s_mb03


subroutine test_error_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp) :: energy
   real(wp), parameter :: ref = -33.040345103604_wp  ! value calculated by xtb

   call get_structure(mol, "MB16-43", "01")
   mol%uhf = mol%uhf + 1

   energy = 0.0_wp

   call new_gfn1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

   call check(error, ctx%failed(), .true.)
   if (allocated(error)) return
   block
      type(error_type), allocatable :: error2
      call ctx%get_error(error2)
      call check(error, allocated(error2))
      if (allocated(error)) return
      call check(error, error2%message, "Total number of electrons (52) and number unpaired electrons (1) is not compatible")
      if (allocated(error)) return
   end block
   call check(error, ctx%failed(), .false.)

end subroutine test_error_mb01


end module test_gfn1_xtb
