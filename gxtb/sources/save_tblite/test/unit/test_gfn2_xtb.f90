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

module test_gfn2_xtb
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_context_type, only : context_type
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, eeq_guess
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_gfn2_xtb

   real(wp), parameter :: acc = 0.01_wp
   real(wp), parameter :: thr = 5e+6_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = 10*sqrt(epsilon(1.0_wp))
   real(wp), parameter :: kt = 3.166808578545117e-06_wp

contains


!> Collect all exported unit tests
subroutine collect_gfn2_xtb(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-atom", test_e_pse), &
      new_unittest("energy-atom-cation", test_e_pse_cation), &
      new_unittest("energy-atom-anion", test_e_pse_anion), &
      new_unittest("energy-mol", test_e_mb01), &
      new_unittest("gradient-mol", test_g_mb02), &
      new_unittest("convergence", test_convergence) &
      ]

end subroutine collect_gfn2_xtb


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
      &-0.3934827590437188_wp,-1.7431266329458670_wp,-0.1800716865751749_wp, &
      &-0.5691059816155888_wp,-0.9524366145568731_wp,-1.7951105194038208_wp, &
      &-2.6094524546320614_wp,-3.7694210954143372_wp,-4.6193399554659953_wp, &
      &-5.9322150527578561_wp,-0.1670967498216340_wp,-0.4659746637924263_wp, &
      &-0.9053286117929751_wp,-1.5714240856447492_wp,-2.3778070880856061_wp, &
      &-3.1482710158768512_wp,-4.4825251342921133_wp,-4.2790432675899392_wp, &
      &-0.1657522390614218_wp,-0.3716463524887721_wp,-0.8541832932455158_wp, &
      &-1.3670573060841023_wp,-1.7181251908054325_wp,-1.7474975140354929_wp, &
      &-2.6191360738389711_wp,-2.9467954760610584_wp,-3.5067898492173386_wp, &
      &-4.6683274353639224_wp,-3.7480061309852375_wp,-0.5275214022962790_wp, &
      &-1.0811118359513132_wp,-1.8099037845688559_wp,-2.2394259485957670_wp, &
      &-3.1204361961032472_wp,-4.0483393705697406_wp,-4.2718555408483345_wp, &
      &-0.1599989486753450_wp,-0.4624308530006268_wp,-1.1948528971312189_wp, &
      &-1.3106537576623940_wp,-1.7812008396036578_wp,-1.7846204189006885_wp, &
      &-2.4820150859709535_wp,-3.0202741319471276_wp,-3.8958154714151534_wp, &
      &-4.4098452999301987_wp,-3.8217382102713393_wp,-0.5330372553013433_wp, &
      &-1.1259377791295933_wp,-2.0128966169144547_wp,-2.1642287880347792_wp, &
      &-3.0090909635721594_wp,-3.7796302627467928_wp,-3.8835884981900253_wp, &
      &-0.1485299624614294_wp,-0.4336420207320540_wp,-1.2047930237499953_wp, &
      &-0.8995088859104427_wp,-0.8933692906260793_wp,-0.8872296591025184_wp, &
      &-0.8810900618954706_wp,-0.8749504289172132_wp,-0.8688108329556412_wp, &
      &-0.8626712738164258_wp,-0.8565316211158370_wp,-0.8503920620019818_wp, &
      &-0.8442524661402829_wp,-0.8381128335296505_wp,-0.8319732376684492_wp, &
      &-0.8258335683086163_wp,-0.8196940091967588_wp,-1.3119987657554162_wp, &
      &-1.9048135505602584_wp,-2.2150777726130775_wp,-3.0080349133279856_wp, &
      &-3.0705864555234426_wp,-3.6450039004282027_wp,-4.4374585721609048_wp, &
      &-3.8026194480681164_wp,-0.8480322467084417_wp,-1.4386851451418587_wp, &
      &-2.2048409808696365_wp,-2.2666534163614425_wp,-2.7347817366291736_wp, &
      &-3.0005430918001297_wp,-3.8578865356212368_wp]

   integer :: izp

   do izp = 1, 86
      call new(mol, [izp], xyz, uhf=uhf(izp))
      call new_gfn2_calculator(calc, mol, error, accuracy=acc)
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
      & 0.2295521666666667_wp,-0.4838821498062669_wp, 0.1659637000000000_wp, &
      & 0.0382879603144248_wp,-0.3028735733744940_wp,-0.8032323515803012_wp, &
      &-1.9127945803017519_wp,-2.9485380631567804_wp,-3.7083708531554338_wp, &
      &-4.7615133525389206_wp, 0.1954855666666667_wp, 0.0177291347704535_wp, &
      &-0.4909477993783444_wp,-1.1519234984400135_wp,-1.8635041144653153_wp, &
      &-2.5869831371080387_wp,-3.8074087383765951_wp,-3.5353072019821710_wp, &
      & 0.1915705000000000_wp, 0.0412624237556138_wp,-0.4876393461395803_wp, &
      &-1.0512553851526367_wp,-1.2276387992117728_wp,-1.3603969166047660_wp, &
      &-2.1531486571035261_wp,-2.6154080943715741_wp,-2.9630475933202085_wp, &
      &-4.1356689668589270_wp,-3.3754409391861020_wp, 0.0954118321851938_wp, &
      &-0.7770718988731139_wp,-1.3380776245766277_wp,-1.7274789356622700_wp, &
      &-2.5727065347449729_wp,-3.3672645851873964_wp,-3.6095838966360381_wp, &
      & 0.2024472666666667_wp,-0.0012154265003135_wp,-0.8205283144950449_wp, &
      &-0.9268352102582531_wp,-1.1831321107080399_wp,-1.3548034028896174_wp, &
      &-2.0946916383295942_wp,-2.4392774466851530_wp,-3.3116283806720501_wp, &
      &-3.8839746929916270_wp,-3.3861214800526964_wp,-0.0014053943173384_wp, &
      &-0.8313996799454058_wp,-1.5499915918867848_wp,-1.7004451789474953_wp, &
      &-2.3281239391233317_wp,-3.1287586760341632_wp,-3.3107276971589097_wp, &
      & 0.1984240000000000_wp,-0.0248826103660272_wp,-0.8043328582495408_wp, &
      &-0.5782798242579857_wp,-0.5658520874417623_wp,-0.5534913129723908_wp, &
      &-0.5411976806647678_wp,-0.5289709450687978_wp,-0.5168113516345733_wp, &
      &-0.5047187940458501_wp,-0.4926934202118626_wp,-0.4807349495874395_wp, &
      &-0.4688434826446651_wp,-0.4570190103215450_wp,-0.4452616792885216_wp, &
      &-0.4335712785474503_wp,-0.4219479871505273_wp,-0.8307280647116214_wp, &
      &-1.4308333086427949_wp,-1.5144142097252595_wp,-2.4712730319839222_wp, &
      &-2.6040096866434608_wp,-3.1254832081174122_wp,-3.8878107594341054_wp, &
      &-3.2726654766627710_wp,-0.1990876900208875_wp,-1.2218139762625606_wp, &
      &-1.8331228693386723_wp,-1.8241019780597743_wp,-2.2952188045868103_wp, &
      &-2.5146156635502148_wp,-3.3252215178719071_wp]

   integer :: izp

   do izp = 1, 86
      if (any(izp == [4, 5, 6])) cycle  ! not converging
      call new(mol, [izp], xyz, uhf=uhf(izp), charge=1.0_wp)
      call new_gfn2_calculator(calc, mol, error, accuracy=acc)
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
      & -0.6107466928548624_wp,-1.5242281622618639_wp,-0.2811010731503498_wp, &
      & -0.0225472225337581_wp,-0.8833252789733284_wp,-1.9209221941523817_wp, &
      & -2.8228473813671178_wp,-4.0689442489122936_wp,-4.9096355171858255_wp, &
      & -5.8525351484414445_wp,-0.2586230663099347_wp,-0.3777253521591805_wp, &
      & -0.9769314233128334_wp,-1.5740591001415243_wp,-2.6051960359906325_wp, &
      & -3.4046900452338020_wp,-4.7851339537609663_wp,-4.0730096623488325_wp, &
      & -0.2754729781228437_wp,-0.2991556954084400_wp,-0.8448287925548863_wp, &
      & -1.4204576650727849_wp,-1.6417599547371493_wp,-1.7349725353358738_wp, &
      & -2.7795518977731328_wp,-2.8451293928733801_wp,-3.7199033615164447_wp, &
      & -4.7200304836864646_wp,-3.9176023227843726_wp,-0.3029556119060911_wp, &
      & -1.1870387794933592_wp,-1.7866119117748076_wp,-2.3590107738630417_wp, &
      & -3.4032788583401246_wp,-4.3322311664719511_wp,-4.1134660808242316_wp, &
      & -0.3119641640173567_wp,-0.3850328498442683_wp,-1.0956747182030877_wp, &
      & -1.4355530563288270_wp,-1.9640589777716178_wp,-2.2204613546902952_wp, &
      & -2.4156522277327448_wp,-3.0555901876944822_wp,-4.0966362950940773_wp, &
      & -4.5493845820943051_wp,-3.9936149404906232_wp,-0.3539529575986469_wp, &
      & -1.2296302217900721_wp,-2.0597413560462923_wp,-2.2773377812690265_wp, &
      & -3.0548938936900618_wp,-4.0603555990360469_wp,-3.7399393131122793_wp, &
      & -0.2589149249228587_wp,-0.4142473193609475_wp,-1.3583939757463996_wp, &
      & -1.1939756095103720_wp,-1.1787219836159635_wp,-1.1634684731233553_wp, &
      & -1.1482147512203462_wp,-1.1329610795197471_wp,-1.1177073276791298_wp, &
      & -1.1024536912403247_wp,-1.0871998378275738_wp,-1.0719461382994526_wp, &
      & -1.0566921636869997_wp,-1.0414384012947371_wp,-1.0261844283637402_wp, &
      & -1.0109304558921322_wp,-0.9956765988223176_wp,-1.4562264589533362_wp, &
      & -2.0816650782323638_wp,-2.5504231469608891_wp,-3.1722176794464403_wp, &
      & -3.2020176312576272_wp,-3.6596592185939323_wp,-4.5670020773987314_wp, &
      & -3.8940033683322537_wp,-0.8351587238649311_wp,-1.5354913944230826_wp, &
      & -2.2747911870200199_wp,-2.3492048567342909_wp,-3.1176094241608174_wp, &
      & -3.2707609413543071_wp,-3.7962721880829253_wp]

   integer :: izp

   do izp = 1, 86
      if (izp == 24) cycle  ! not converging
      call new(mol, [izp], xyz, uhf=uhf(izp), charge=-1.0_wp)
      call new_gfn2_calculator(calc, mol, error, accuracy=acc)
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
   real(wp), parameter :: ref = -30.348902328491_wp  ! value calculated by xtb

   call get_structure(mol, "MB16-43", "01")

   energy = 0.0_wp

   call new_gfn2_calculator(calc, mol, error, accuracy=acc)
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
      & 3.8555117961272E-03_wp,  -2.2109558455227E-03_wp,  -9.1601598199788E-03_wp, &
      & 1.0970538791192E-03_wp,   1.3996304174584E-02_wp,  -5.3472436268014E-03_wp, &
      &-8.0586101184995E-03_wp,  -2.3489741441695E-03_wp,   1.9653192119336E-02_wp, &
      &-4.4938867765893E-03_wp,  -4.5671183892695E-04_wp,  -1.5746269201745E-02_wp, &
      & 4.6136808085153E-03_wp,  -8.3294420588366E-03_wp,  -6.7825658568501E-03_wp, &
      &-4.0646538342810E-05_wp,  -7.6816577111296E-03_wp,   7.7872936300304E-03_wp, &
      & 7.8804339745087E-04_wp,  -4.9765330722152E-03_wp,   5.9586244312432E-03_wp, &
      & 5.3161702509583E-04_wp,   5.1936992720862E-03_wp,   8.0564683826072E-03_wp, &
      &-9.2289677830967E-03_wp,  -1.5299699935381E-02_wp,  -2.1407939802147E-02_wp, &
      & 6.9445923135428E-03_wp,  -3.7433916593696E-04_wp,   5.8421624969437E-03_wp, &
      & 3.5614557208870E-03_wp,   7.6527551612558E-03_wp,   4.8040821101276E-03_wp, &
      &-1.1121294159972E-02_wp,  -6.4733243214965E-03_wp,   7.1240862931304E-04_wp, &
      & 5.0827933118851E-03_wp,   1.2142952146235E-02_wp,   9.1889595078325E-04_wp, &
      & 1.5094799594745E-03_wp,  -2.3321073116725E-03_wp,   4.8244303955350E-03_wp, &
      & 9.2444420193416E-03_wp,   1.2118997544732E-03_wp,  -1.2700147822202E-03_wp, &
      &-4.2852648549389E-03_wp,   1.0286134896654E-02_wp,   1.1566349438238E-03_wp],&
      & shape(ref))

   call get_structure(mol, "MB16-43", "02")

   allocate(gradient(3, mol%nat), sigma(3, 3))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call new_gfn2_calculator(calc, mol, error, accuracy=acc)
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


subroutine test_convergence(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   integer, parameter :: num(2) = [3, 8]
   real(wp), parameter :: xyz(3, 2) = reshape([&
      & +0.00000000000000_wp, +0.00000000000000_wp, +1.50105302628963_wp, &
      & +0.00000000000000_wp, +0.00000000000000_wp, -1.50105302628963_wp],&
      & shape(xyz))
   real(wp) :: energy
   real(wp), parameter :: ref = -4.228326553369_wp  ! value calculated by xtb

   call new(mol, num, xyz)

   energy = 0.0_wp

   call new_gfn2_calculator(calc, mol, error, accuracy=acc)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
      & 1, calc%default_etemp * kt)
   call eeq_guess(mol, calc, wfn, error)
   if (allocated(error)) return
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)

   call check(error, energy, ref, thr=1e-7_wp)

end subroutine test_convergence


end module test_gfn2_xtb
