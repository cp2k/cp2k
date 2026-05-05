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

module test_ipea1_xtb
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_context_type, only : context_type
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_ipea1, only : new_ipea1_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_ipea1_xtb

   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 100_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 3.166808578545117e-06_wp

contains


!> Collect all exported unit tests
subroutine collect_ipea1_xtb(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-atom", test_e_pse), &
      new_unittest("energy-atom-cation", test_e_pse_cation), &
      new_unittest("energy-atom-anion", test_e_pse_anion), &
      new_unittest("energy-mol", test_e_mb01), &
      new_unittest("gradient-mol", test_g_mb02), &
      new_unittest("numgrad-mol", test_g_mb03) &
      !new_unittest("virial-mol", test_s_mb03) &
      ]

end subroutine collect_ipea1_xtb


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
      &-0.3866742167172414_wp,-1.6258646654198721_wp,-0.2491081210142224_wp, &
      &-0.6729203958693754_wp,-1.1779317882934104_wp,-1.8958661905978365_wp, &
      &-3.3610873801335304_wp,-4.3197653712023074_wp,-5.2322787376334388_wp, &
      &-6.3190030086322286_wp,-0.1993996626790595_wp,-0.7780120696208370_wp, &
      &-1.0775966895776903_wp,-1.6252631326184159_wp,-2.3375260307637005_wp, &
      &-3.5875767688126725_wp,-4.2732461508565125_wp,-6.1469954466724239_wp, &
      &-0.2137179709637750_wp,-0.5864589388109813_wp,-0.5052556906812852_wp, &
      & 0.7937419091568174_wp,-1.6449950369399147_wp,-1.7022631944075812_wp, &
      &-1.7263559458095221_wp,-2.4152287584704664_wp,-2.9396210675617334_wp, &
      &-3.4786563319156572_wp,-4.4169348599647522_wp,-0.7010875030780547_wp, &
      &-1.2381959018228419_wp,-1.9032223027699005_wp,-2.6475620316691373_wp, &
      &-3.3630683144828950_wp,-3.8180395534595282_wp,-4.8421263290307968_wp, &
      &-0.2799176883500828_wp,-0.5027433177275922_wp,-0.8559774488167702_wp, &
      &-0.8259470957003716_wp,-1.3376542753718252_wp,-1.6557219239736551_wp, &
      &-2.2160358545342351_wp,-3.2181664696662833_wp,-3.9147305338504088_wp, &
      &-4.0411480991493995_wp,-3.8776985063010745_wp,-0.6814552058238580_wp, &
      &-1.2557641283776890_wp,-2.2117635981910682_wp,-2.3199789636400441_wp, &
      &-3.5896871356541111_wp,-3.8857572752278440_wp,-4.2321989780523870_wp, &
      &-0.2330413901515095_wp,-0.4742595783451311_wp,-0.6959185846380349_wp, &
      &-0.6485304192347467_wp,-0.6450212630347002_wp,-0.6415120701014334_wp, &
      &-0.6380028404220172_wp,-0.6344936474925504_wp,-0.6309844913125303_wp, &
      &-0.6274752983832106_wp,-0.6239660687045753_wp,-0.6204569125245851_wp, &
      &-0.6169477195952713_wp,-0.6134385634152814_wp,-0.6099292969873201_wp, &
      &-0.6064201408073303_wp,-0.6029109478780164_wp,-0.5889837270645610_wp, &
      &-1.6680608750310100_wp,-1.7392411746034055_wp,-2.0215567377448149_wp, &
      &-2.9125503382127143_wp,-3.5331732272481489_wp,-4.3179020613576666_wp, &
      &-4.0419331894786650_wp,-0.8551180308258510_wp,-0.9551323610773361_wp, &
      &-1.1004919026223705_wp,-2.2682401786640844_wp,-3.0676044864514225_wp, &
      &-3.4437503367483560_wp,-3.6081562853045592_wp]

   integer :: izp

   do izp = 1, 86
      if (any(izp == [21, 22, 24, 26, 40, 73])) cycle ! SCF does not converge
      call new(mol, [izp], xyz, uhf=uhf(izp))
      call new_ipea1_calculator(calc, mol, error, accuracy=acc)
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
      & 0.2524195000000000_wp,-0.0422428129239876_wp, 0.1596375666666667_wp, &
      &-0.1288451979346878_wp,-0.6774537388560994_wp,-1.1736448210760286_wp, &
      &-2.6006847322776907_wp,-3.5331911641456291_wp,-4.2162790243959840_wp, &
      &-5.4089550104332629_wp, 0.1132375000000000_wp,-0.1858943681437519_wp, &
      &-0.7929165898783470_wp,-1.2064302987580611_wp,-1.7026143817446238_wp, &
      &-2.9245283793472860_wp,-3.5323797219917012_wp,-5.2547187063374068_wp, &
      & 0.0816097666666667_wp,-0.1889038694054907_wp,-0.1840771026585860_wp, &
      & 0.0000000000000000_wp,-1.1214099690902426_wp,-1.3598407694042540_wp, &
      & 0.0000000000000000_wp,-1.7792714340771387_wp,-2.5486583996246086_wp, &
      &-3.0005527001603487_wp,-3.9315697627987665_wp,-0.1115875848723608_wp, &
      &-0.8758637148878993_wp,-1.3299982534079502_wp,-2.0542518849438327_wp, &
      &-2.7306615327593713_wp,-3.0519213538372423_wp,-3.9953892543401519_wp, &
      & 0.0878675000000000_wp,-0.1466078255304628_wp,-0.4937299864801741_wp, &
      &-0.5519585048615329_wp,-1.0596813265940885_wp,-1.2582804562299617_wp, &
      &-1.7441817291584463_wp,-2.6949495714667360_wp,-3.3840208304086556_wp, &
      &-3.5206612682222165_wp,-3.3889629189385042_wp,-0.1292345362452623_wp, &
      &-0.9580171953907000_wp,-1.7400454800940868_wp,-1.8399584688273154_wp, &
      &-3.0318097395444616_wp,-3.2276435086854969_wp,-3.5173994843193732_wp, &
      & 0.0925550000000000_wp,-0.1182202891725656_wp,-0.3383862439154567_wp, &
      &-0.3380276085842653_wp,-0.3364513307550568_wp,-0.3350141427496482_wp, &
      &-0.3337162359090419_wp,-0.3325572206289337_wp,-0.3315373115184752_wp, &
      &-0.3306564922318164_wp,-0.3299149902654105_wp,-0.3293123596652018_wp, &
      &-0.3288488188887932_wp,-0.3285243679361848_wp,-0.3283389200405289_wp, &
      &-0.3282929953096204_wp,-0.3283858232519615_wp,-0.2841939827671290_wp, &
      &-1.3015199696124606_wp,-1.3422789879217938_wp,-1.6201970652257776_wp, &
      &-1.6774725026379260_wp,-3.1852214803626704_wp,-3.7295747828824615_wp, &
      & 0.0000000000000000_wp,-0.2609871820795922_wp,-0.6672190773979373_wp, &
      &-0.7062159969566204_wp,-1.7430127291334649_wp,-2.4123355716184571_wp, &
      &-2.6886076809857715_wp,-2.9576573558829384_wp]

   integer :: izp

   do izp = 1, 86
      if (any(izp == [22, 25, 79])) cycle  ! SCF does not converge
      call new(mol, [izp], xyz, uhf=uhf(izp), charge=1.0_wp)
      call new_ipea1_calculator(calc, mol, error, accuracy=acc)
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
      &-0.6209289334344827_wp, 0.0000000000000000_wp,-0.4616778086951114_wp, &
      &-0.8227993117072305_wp,-1.2321973146555785_wp,-1.9940631078868596_wp, &
      &-3.5943477724466160_wp,-4.5750846155937186_wp,-5.6277444764354101_wp, &
      &-6.0135919541021741_wp,-0.3655618253606222_wp,-0.7293187727421735_wp, &
      &-1.1294568900412956_wp,-1.8636984933817820_wp,-2.5355513957879081_wp, &
      &-3.8827514585949539_wp,-4.6781445359617297_wp,-5.9532003542728642_wp, &
      &-0.3946877084813858_wp,-0.6601217419378675_wp, 0.0000000000000000_wp, &
      & 0.0000000000000000_wp,-1.6404683240059497_wp,-1.7545829076187183_wp, &
      & 0.0000000000000000_wp,-2.6868691903680801_wp,-3.0365322826962613_wp, &
      &-3.5891575157336724_wp,-4.6694329571307351_wp,-0.7478911585526812_wp, &
      &-1.3856578375924125_wp,-2.1638075079276167_wp,-2.8550820492802025_wp, &
      &-3.6832851496456129_wp,-3.9571135364820273_wp,-4.6240117892188310_wp, &
      &-0.5719678767001657_wp,-0.6583170517916040_wp,-1.0130429561259551_wp, &
      & 0.0000000000000000_wp,-1.4288569875823489_wp,-1.8026564185949252_wp, &
      & 0.0000000000000000_wp,-3.3802910855544024_wp,-4.0402591326658541_wp, &
      &-4.0324827651799442_wp,-3.9263250936636416_wp,-0.7543768101126942_wp, &
      &-1.3714541989645710_wp,-2.4316641193407351_wp,-2.3609644225466564_wp, &
      &-3.8939401009678702_wp,-4.0435926774613034_wp,-3.9095090580155762_wp, &
      &-0.4735277804823783_wp,-0.6639055123686469_wp, 0.0000000000000000_wp, &
      & 0.0000000000000000_wp, 0.0000000000000000_wp,-0.7774600374514302_wp, &
      &-0.7692826709429081_wp,-0.7615788166990638_wp,-0.7541577457962461_wp, &
      &-0.7468810501730819_wp,-0.7397437693038101_wp,-0.7327453598812527_wp, &
      &-0.7258859667911156_wp,-0.7191657370235338_wp,-0.7125844368143500_wp, &
      &-0.7061427332686532_wp,-0.6998397088976513_wp,-0.5919814184334333_wp, &
      &-1.8228791156675879_wp,-1.7570231210779892_wp,-2.2234177404606594_wp, &
      &-2.9268153486138431_wp, 0.0000000000000000_wp,-4.3806893728555032_wp, &
      &-4.2787845242354283_wp,-0.9229174394105254_wp,-1.1068049530063353_wp, &
      & 0.0000000000000000_wp,-2.4396056714685725_wp,-3.0614995810384213_wp, &
      &-3.3648843601107044_wp,-3.2192830086977700_wp]

   integer :: izp

   do izp = 1, 86
      if (izp == 2) cycle  ! Helium doesn't have enough orbitals for negative charge
      if (any(izp == [21, 22, 25, 40, 43, 57, 58, 59, 77, 82])) cycle  ! not converging
      call new(mol, [izp], xyz, uhf=uhf(izp), charge=-1.0_wp)
      call new_ipea1_calculator(calc, mol, error, accuracy=acc)
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
   real(wp), parameter :: ref = -35.413681615685_wp  ! value calculated by xtb

   call get_structure(mol, "MB16-43", "01")

   energy = 0.0_wp

   call new_ipea1_calculator(calc, mol, error, accuracy=acc)
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
      &  5.7754024469751E-03_wp,  -2.6658478685737E-03_wp,  -9.9836761519953E-03_wp, &
      & -1.1531812851444E-02_wp,   9.1493770962829E-03_wp,  -1.0826562359077E-03_wp, &
      & -4.5105398889876E-03_wp,  -2.5940495132499E-03_wp,   8.5884774999636E-03_wp, &
      & -1.2583415176846E-02_wp,  -9.9164814583170E-04_wp,  -2.4865397207021E-02_wp, &
      & -6.2312269681391E-03_wp,  -4.0424657610046E-03_wp,  -1.5150696737034E-02_wp, &
      & -3.1829365241243E-03_wp,  -1.7354975024896E-02_wp,   1.4341675119837E-02_wp, &
      &  5.8816317748151E-03_wp,  -2.6319916485481E-03_wp,   1.4252771728896E-03_wp, &
      &  1.0550028517182E-02_wp,  -1.0434581209469E-03_wp,   2.8047985700563E-03_wp, &
      & -1.6519777094403E-02_wp,  -7.2806520709268E-05_wp,   2.1400870147015E-02_wp, &
      &  2.0722923477563E-02_wp,   7.7305758264574E-03_wp,   7.8951317863186E-03_wp, &
      & -1.8046425156298E-02_wp,   6.4969554747904E-03_wp,   9.2735137832130E-03_wp, &
      & -9.6336555428399E-03_wp,   1.0303753120328E-03_wp,  -4.2400659781200E-03_wp, &
      &  9.3274049575810E-03_wp,   8.3669130497955E-03_wp,   1.0175543467468E-02_wp, &
      &  1.2803185479064E-02_wp,  -1.3435513444702E-02_wp,   9.7585465719697E-03_wp, &
      &  1.3863921847308E-02_wp,  -6.9554030525201E-03_wp,  -7.9147209305042E-03_wp, &
      &  3.3152907025928E-03_wp,   1.9013962341624E-02_wp,  -2.2426620878149E-02_wp],&
      & shape(ref))

   call get_structure(mol, "MB16-43", "02")

   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_ipea1_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, gradient, sigma, 0)

   if (any(abs(gradient - ref) > 10*thr2)) then
      call test_failed(error, "Gradient of energy does not match (xTB)")
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

   call get_structure(mol, "MB16-43", "03")

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   numgrad(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_ipea1_calculator(calc, mol, error, accuracy=acc)
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
   real(wp), allocatable :: gradient(:, :)
   real(wp), allocatable :: sigma(:, :), numsigma(:, :)

   call get_structure(mol, "MB16-43", "03")

   allocate(gradient(3, mol%nat), sigma(3, 3), numsigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_ipea1_calculator(calc, mol, error, accuracy=acc)
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


end module test_ipea1_xtb
