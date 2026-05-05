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

module test_param
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use dftd3
   implicit none
   private

   public :: collect_param

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)


contains


!> Collect all exported unit tests
subroutine collect_param(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("DFT-D3(BJ)", test_d3bj_mb01), &
      & new_unittest("DFT-D3(0)", test_d3zero_mb09), &
      & new_unittest("DFT-D3(BJM)", test_d3bjm_mb02), &
      & new_unittest("DFT-D3(0M)", test_d3zerom_mb03), &
      & new_unittest("DFT-D3(op)", test_d3op_mb06), &
      & new_unittest("DFT-D3(BJ)-ATM", test_d3bjatm_mb17), &
      & new_unittest("DFT-D3(0)-ATM", test_d3zeroatm_mb25), &
      & new_unittest("DFT-D3(BJM)-ATM", test_d3bjmatm_mb04), &
      & new_unittest("DFT-D3(0M)-ATM", test_d3zeromatm_mb05), &
      & new_unittest("DFT-D3(op)-ATM", test_d3opatm_mb07), &
      & new_unittest("DFT-D3(CSO)", test_d3cso_mb01), &
      & new_unittest("DFT-D3(CSO)-ATM", test_d3csoatm_mb02), &
      & new_unittest("unknown-D3(BJ)", test_d3bj_unknown, should_fail=.true.), &
      & new_unittest("unknown-D3(0)", test_d3zero_unknown, should_fail=.true.), &
      & new_unittest("unknown-D3(BJM)", test_d3bjm_unknown, should_fail=.true.), &
      & new_unittest("unknown-D3(0M)", test_d3zerom_unknown, should_fail=.true.), &
      & new_unittest("unknown-D3(op)", test_d3op_unknown, should_fail=.true.), &
      & new_unittest("unknown-D3(CSO)", test_d3cso_unknown, should_fail=.true.) &
      & ]

end subroutine collect_param


subroutine test_dftd3_gen(error, mol, param, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Damping parameters
   class(damping_param), intent(in) :: param

   !> Expected dispersion energy
   real(wp), intent(in) :: ref

   type(d3_model) :: d3
   real(wp) :: energy

   call new_d3_model(d3, mol)
   call get_dispersion(mol, d3, param, realspace_cutoff(), energy)

   call check(error, energy, ref, thr=thr)
   if (allocated(error)) then
      print '(ES23.16)', energy
   end if

end subroutine test_dftd3_gen


subroutine test_d3bj_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp

   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      & "b1b95", "b2gpplyp", "b2plyp", "b3lyp", "b3lyp/631gd", "b3pw91", "b97d", &
      & "bhlyp", "blyp", "bmk", "bop", "bp", "bpbe", "camb3lyp", "dftb3", &
      & "dsdblyp", "dsdblypfc", "hcth120", "hf", "hf3c", "hf3cv", "hf/minis", &
      & "hf/mixed", "hf/sv", "hse06", "hsesol", "lcwpbe", "mpw1b95", "mpwlyp", &
      & "olyp", "opbe", "otpss", "pbe", "pbe0", "pbeh3c", "pbesol", "ptpss", &
      & "pw1pw", "pw6b95", "pwb6k", "pwgga", "pwpb95", "revpbe", "revpbe0", &
      & "revpbe38", "revssb", "rpbe", "rpw86pbe", "ssb", "tpss", "tpss0", "tpssh", &
      & "scan", "rscan", "r2scan", "b97m", "wb97m", "wb97x", &
      & "pbehpbe", "xlyp", "mpwpw", "hcth407", "revtpss", "tauhcth", "b3p", "b1p", &
      & "b1lyp", "mpw1pw", "mpw1kcis", "pbeh1pbe", "pbe1kcis", "x3lyp", "o3lyp", &
      & "b971", "b972", "b98", "hiss", "hse03", "revtpssh", "revtpss0", "tpss1kcis", &
      & "tauhcthhyb", "mn15", "lcwhpbe", "mpw2plyp", "m11", "sogga11x", "n12sx", "mn12sx", &
      & "r2scanh", "r2scan0", "r2scan50", "b973c", "b3lyp3", "b3lyp5", "b3lypg", &
      & "dm21", "dm21m", "dm21mc", "dm21mu", "dsdpbep86", "dsdpbeb95", "dsdpbe", &
      & "dodscan66", "revdsdblyp", "revdsdpbep86", "revdsdpbeb95", "revdsdpbe", &
      & "revdodblyp", "revdodpbep86", "revdodpbeb95", "revdodpbe", "pw91", "drpa75", &
      & "scsdrpa75", "optscsdrpa75", "dsdpbedrpa75", "dsdpbep86drpa75", &
      & "dsdsvwn5", "dsdsp86", "dsdslyp", "dsdspbe", "dsdbvwn5", "dsdblyp_2013", "dsdbpbe", &
      & "dsdbp86", "dsdbpw91", "dsdbb95", "dsdpbevwn5", "dsdpbelyp", "dsdpbepw91", &
      & "dsdpbehb95", "dsdpbehp86", "dsdmpwlyp", "dsdmpwpw91", "dsdmpwp86", "dsdmpwpbe", &
      & "dsdmpwb95", "dsdhsepbe", "dsdhsepw91", "dsdhsep86", "dsdhselyp", "dsdtpss", &
      & "dsdtpsstpss", "dsdtpssb95", "dsdolyp", "dsdxlyp", "dsdxb95", "dsdb98", "dsdbmk", &
      & "dsdthcth", "dsdhcth407", "dodsvwn5", "dodblyp", "dodpbe", "dodpbep86", &
      & "dodpbeb95", "dodhsep86", "dodpbehb95", "dsdpbep86_2011", "skala-1.0", &
      & "mpwkcis1k"]
   real(wp), parameter :: ref(*) = [&
      &-2.9551694676908012E-2_wp,-1.6638703086788331E-2_wp,-1.6725877716130381E-2_wp, &
      &-3.3014429592265318E-2_wp,-2.2051435219996540E-2_wp,-3.3481565825316001E-2_wp, &
      &-4.4319254094064855E-2_wp,-2.7768780192796199E-2_wp,-3.9702051070377456E-2_wp, &
      &-2.9455263494312305E-2_wp,-6.4253095971854787E-2_wp,-3.2374555740782837E-2_wp, &
      &-4.0184787340500967E-2_wp,-1.6566618761190383E-2_wp,-1.6283175554508783E-2_wp, &
      &-1.8049466820672877E-2_wp,-1.8968229257655511E-2_wp,-3.1793648125734875E-2_wp, &
      &-1.2496497728657682E-1_wp,-7.3593575280702983E-2_wp,-4.0603488727796531E-2_wp, &
      &-1.4163816509828264E-1_wp,-2.0131394596367987E-2_wp,-3.4278101154601628E-2_wp, &
      &-1.4158700987643606E-2_wp,-7.6707358124053241E-3_wp,-1.8974039990743818E-2_wp, &
      &-1.3757652907307207E-2_wp,-1.9765756364360013E-2_wp,-8.2952458571961249E-2_wp, &
      &-7.7457838043986788E-2_wp,-3.1481181733832270E-2_wp,-1.7882220155186028E-2_wp, &
      &-1.6523309049788101E-2_wp,-8.8508561096339745E-3_wp,-8.5090062012807676E-3_wp, &
      &-1.7255799803511284E-2_wp,-1.2406518080297352E-2_wp,-1.1629315923150214E-2_wp, &
      &-5.5085428675363066E-3_wp,-1.5678626937248866E-2_wp,-1.0788648445240084E-2_wp, &
      &-4.1947444551537844E-2_wp,-3.7430728490769900E-2_wp,-3.5154037836964684E-2_wp, &
      &-1.6039212958496831E-2_wp,-1.0199409751803663E-1_wp,-1.8473409930576331E-2_wp, &
      &-4.5113899201072810E-2_wp,-2.3465719708170546E-2_wp,-2.5011519099279463E-2_wp, &
      &-2.2089903358298878E-2_wp,-4.1432013834768475E-3_wp,-6.6123761503017447E-3_wp, &
      &-5.4507719557227875E-3_wp,-5.8854885361595541E-2_wp,-2.1655168402602998E-2_wp, &
      &-4.8891114763583890E-2_wp,-3.2173513754592320E-2_wp,-8.1794582931313600E-2_wp, &
      &-3.5666193218422140E-2_wp,-1.2744920440901275E-1_wp,-2.0878779644526745E-2_wp, &
      &-8.9701534664655344E-2_wp,-2.1917075091994892E-2_wp,-2.0696460275344126E-2_wp, &
      &-4.5547718946122910E-2_wp,-2.7939171637344121E-2_wp,-6.2045224493951778E-2_wp, &
      &-2.9161162131212684E-2_wp,-3.8203961549843249E-2_wp,-3.6103963511026425E-2_wp, &
      &-4.4759235392592758E-2_wp,-3.2738407689525253E-2_wp,-8.6635467284479059E-2_wp, &
      &-4.3851712866357175E-2_wp,-2.3700538936935606E-2_wp,-2.8136697973180295E-2_wp, &
      &-2.4751360640335753E-2_wp,-2.3986533737266124E-2_wp,-5.6113880355570932E-2_wp, &
      &-2.6879433016433214E-3_wp,-4.7068039876219803E-5_wp,-2.2830467677890126E-2_wp, &
      &-8.9046741008620024E-3_wp,-4.2882619531553470E-3_wp,-3.4182032740020921E-2_wp, &
      &-1.7851923991011053E-2_wp,-7.2749942776428288E-3_wp,-5.9043780290722733E-3_wp, &
      &-6.5370091760150479E-3_wp,-7.3564457249637944E-3_wp,-4.2958524322550894E-2_wp, &
      &-3.3014429592265318E-2_wp,-3.3014429592265318E-2_wp,-3.3014429592265318E-2_wp, &
      &-3.3014429592265318E-2_wp,-3.3014429592265318E-2_wp,-3.3014429592265318E-2_wp, &
      &-3.3014429592265318E-2_wp,-1.5694642957869546E-2_wp,-1.2938686000156000E-2_wp, &
      &-1.7758051042408327E-2_wp,-9.2323777698865989E-3_wp,-2.4042802546329414E-2_wp, &
      &-1.5411238160834008E-2_wp,-1.2978255394296128E-2_wp,-2.0231431225074755E-2_wp, &
      &-2.7103838130103503E-2_wp,-1.6794975103307797E-2_wp,-1.4460579192722246E-2_wp, &
      &-2.1361659109385409E-2_wp,-9.7460845051874183E-3_wp,-2.8391262582808205E-2_wp, &
      &-1.9116077511614916E-2_wp,-1.9252188823010906E-2_wp,-2.4371486479404614E-2_wp, &
      &-2.2775959440262707E-2_wp,-1.5040699501291652E-2_wp,-8.4730446664204769E-3_wp, &
      &-9.8091518486684674E-3_wp,-9.7798508159701748E-3_wp,-2.6905355995708924E-2_wp, &
      &-2.1623262868629957E-2_wp,-1.9597764351264244E-2_wp,-1.8581716550343332E-2_wp, &
      &-1.9615520517321777E-2_wp,-1.4301725581560211E-2_wp,-2.5702259085751418E-2_wp, &
      &-1.8966070619925963E-2_wp,-1.7848227739145565E-2_wp,-1.2302357180476200E-2_wp, &
      &-1.5040699501291648E-2_wp,-1.9629301479961281E-2_wp,-1.9089864590394103E-2_wp, &
      &-1.6663654510626937E-2_wp,-1.8981014496672229E-2_wp,-1.3172267842653015E-2_wp, &
      &-1.7985718363464843E-2_wp,-1.8092724009544817E-2_wp,-1.5040699501291648E-2_wp, &
      &-1.7642856390628806E-2_wp,-1.2388749800413754E-2_wp,-1.2388749800413754E-2_wp, &
      &-6.2394376475324361E-3_wp,-2.6266438465903479E-2_wp,-2.0856132822458857E-2_wp, &
      &-1.3803848373701565E-2_wp,-1.0144589478429857E-2_wp,-2.0904528949713987E-2_wp, &
      &-2.3395356334306699E-2_wp,-2.7235779991786224E-2_wp,-1.8637388512470085E-2_wp, &
      &-4.5692905041335857E-2_wp,-2.3906657015347409E-2_wp,-2.7313595202479946E-2_wp, &
      &-1.7359235198347055E-2_wp,-2.6175528735709948E-2_wp,-1.6381250116750040E-2_wp, &
      &-1.3173453514111508E-2_wp,-3.3014429592265318E-2_wp,-4.2066538266174014E-2_wp]

   call get_structure(mol, "MB16-43", "01")
   do ii = 1, size(func)
      call get_rational_damping(inp, trim(func(ii)), error)
      if (allocated(error)) exit
      call new_rational_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3bj_mb01


subroutine test_d3zero_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp
   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      & "slaterdiracexchange", "blyp", "bp", "b97d", "revpbe", "pbe", "pbesol", &
      & "rpw86pbe", "rpbe", "tpss", "b3lyp", "pbe0", "hse06", "revpbe38", &
      & "pw6b95", "tpss0", "b2plyp", "pwpb95", "b2gpplyp", "ptpss", "hf", &
      & "mpwlyp", "bpbe", "bhlyp", "tpssh", "pwb6k", "b1b95", "bop", "olyp", &
      & "opbe", "ssb", "revssb", "otpss", "b3pw91", "revpbe0", "pbe38", &
      & "mpw1b95", "mpwb1k", "bmk", "camb3lyp", "lcwpbe", "m05", "m052x", &
      & "m06l", "m06", "m062x", "m06hf", "hcth120", "scan", "wb97x", &
      & "mpw1lyp", "mpwkcis1k", "pkzb", "n12", "m08hx", "m11l", "mn15l", "pwp", &
      & "pw1pw", "pbehpbe", "xlyp", "mpwpw", "hcth407", "revtpss", "tauhcth", &
      & "b3p", "b1p", "b1lyp", "mpw1pw", "mpw1kcis", "pbeh1pbe", "pbe1kcis", &
      & "x3lyp", "o3lyp", "b971", "b972", "b98", "hiss", "hse03", "revtpssh", &
      & "revtpss0", "tpss1kcis", "tauhcthhyb", "mpw2plyp", "b973c", "cf22d"]
   real(wp), parameter :: ref(*) = [&
      & 1.4617000329030605E-1_wp,-1.4741267229767294E-2_wp,-1.3716369898239468E-2_wp, &
      &-2.0673049860038258E-2_wp,-1.8741296181572904E-2_wp,-6.7002000141365174E-3_wp, &
      &-4.9290963844481634E-3_wp,-7.6895496501738787E-3_wp,-2.0178760785797962E-2_wp, &
      &-9.7472759549454002E-3_wp,-1.2109364912853944E-2_wp,-7.2249645365211976E-3_wp, &
      &-4.4683641693542752E-3_wp,-1.2242138743812855E-2_wp,-5.7528633397095361E-3_wp, &
      &-9.4290214142670348E-3_wp,-6.6997968702894161E-3_wp,-4.6620737482248901E-3_wp, &
      &-4.8252459487184342E-3_wp,-5.6960358502059277E-3_wp,-1.3736367500197053E-2_wp, &
      &-8.7008716258934124E-3_wp,-1.7036681072722496E-2_wp,-9.7807895043091932E-3_wp, &
      &-9.6115127055308430E-3_wp,-3.6502600790935453E-3_wp,-1.1637668908954207E-2_wp, &
      &-2.4093376405435917E-2_wp,-3.5647727873373614E-2_wp,-3.3269256563457195E-2_wp, &
      &-6.3710140907560419E-3_wp,-5.6777733335476623E-3_wp,-1.2808040338504561E-2_wp, &
      &-1.3606439741693448E-2_wp,-1.5611919106945399E-2_wp,-7.3187556909606131E-3_wp, &
      &-7.1465584485645543E-3_wp,-6.7044103692842994E-3_wp,-1.3136092587759902E-2_wp, &
      &-8.3906362864872280E-3_wp,-8.8794769166744038E-3_wp,-4.6788817240305874E-3_wp, &
      &-9.1484925069627409E-4_wp,-4.7465277362242642E-4_wp,-1.3751381718529450E-3_wp, &
      &-4.0846404480368841E-4_wp,-8.1141559233198794E-4_wp,-9.5577071142980783E-3_wp, &
      &-1.3816180712373748E-3_wp,-4.6307006508886413E-3_wp,-1.1799350053812951E-2_wp, &
      &-1.1234976685244820E-2_wp,-6.9517313125531277E-2_wp,-1.5593189251623335E-2_wp, &
      &-3.9932890855235272E-4_wp,-6.7006180343081945E-3_wp,-2.6285355487601172E-7_wp, &
      &-5.3097811920292193E-3_wp,-7.7402900184866570E-3_wp,-8.9096143513293247E-3_wp, &
      &-1.6037859346177272E-2_wp,-1.2799753733644717E-2_wp,-1.6633282433869329E-2_wp, &
      &-9.4380662409846317E-3_wp,-1.5415561344423432E-2_wp,-9.9189752545548755E-3_wp, &
      &-9.5912576431646454E-3_wp,-1.2799753733644717E-2_wp,-1.0497300192168362E-2_wp, &
      &-1.4033156903040580E-2_wp,-7.3749889527684152E-3_wp,-1.0771401497706931E-2_wp, &
      &-9.7960459201707675E-3_wp,-1.1804087267700464E-2_wp,-9.8608424063707390E-3_wp, &
      &-1.0148393050164134E-2_wp,-1.1462890566584262E-2_wp,-5.8934014233026572E-3_wp, &
      &-7.1067312076729611E-3_wp,-8.9020818106658166E-3_wp,-8.0382680164769616E-3_wp, &
      &-1.2773548954558107E-2_wp,-1.0444034180767760E-2_wp,-4.8721698293720563E-3_wp, &
      &-1.4639383757783959E-2_wp,-5.8016654989627144E-4_wp]

   call get_structure(mol, "MB16-43", "09")
   do ii = 1, size(func)
      call get_zero_damping(inp, trim(func(ii)), error)
      if (allocated(error)) return
      call new_zero_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3zero_mb09


subroutine test_d3bjatm_mb17(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp
   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      & "b1b95", "b2gpplyp", "b2plyp", "b3lyp", "b3lyp/631gd", "b3pw91", "b97d", &
      & "bhlyp", "blyp", "bmk", "bop", "bp", "bpbe", "camb3lyp", "dftb3", &
      & "dsdblyp", "dsdblypfc", "hcth120", "hf", "hf3c", "hf3cv", "hf/minis", &
      & "hf/mixed", "hf/sv", "hse06", "hsesol", "lcwpbe", "mpw1b95", "mpwlyp", &
      & "olyp", "opbe", "otpss", "pbe", "pbe0", "pbeh3c", "pbesol", "ptpss", &
      & "pw1pw", "pw6b95", "pwb6k", "pwgga", "pwpb95", "revpbe", "revpbe0", &
      & "revpbe38", "revssb", "rpbe", "rpw86pbe", "ssb", "tpss", "tpss0", "tpssh", &
      & "scan", "rscan", "r2scan", "b97m", "wb97m", "wb97x", &
      & "r2scanh", "r2scan0", "r2scan50", "b973c", "dsdpbep86", "dsdpbeb95", "dsdpbe", &
      & "dodscan66", "revdsdblyp", "revdsdpbep86", "revdsdpbeb95", "revdsdpbe", &
      & "revdodblyp", "revdodpbep86", "revdodpbeb95", "revdodpbe", "skala-1.0"]
   real(wp), parameter :: ref(*) = [&
      &-2.3886024757749025E-2_wp,-1.2511386468651674E-2_wp,-1.4044660238061260E-2_wp, &
      &-2.8422909990177846E-2_wp,-1.9411341783127211E-2_wp,-2.9014208237027887E-2_wp, &
      &-3.9594831005456688E-2_wp,-2.3209274743657832E-2_wp,-3.4356534352042042E-2_wp, &
      &-2.3548743355403044E-2_wp,-5.6222757182437889E-2_wp,-2.7737517963949018E-2_wp, &
      &-3.4944171516576986E-2_wp,-1.4207748203259795E-2_wp,-1.4807550127492797E-2_wp, &
      &-1.3592755832923201E-2_wp,-1.4297170123439285E-2_wp,-2.7202041052423038E-2_wp, &
      &-1.0723400922969686E-1_wp,-6.4226909111939484E-2_wp,-3.4476703783611985E-2_wp, &
      &-1.1444511063258883E-1_wp,-1.7961031598334974E-2_wp,-2.9706520132830379E-2_wp, &
      &-1.2161446446433095E-2_wp,-6.6471542148879674E-3_wp,-1.6387020088052297E-2_wp, &
      &-1.1272763044393083E-2_wp,-1.7452571215069202E-2_wp,-7.3575769466178789E-2_wp, &
      &-6.8870148867061456E-2_wp,-2.7534795830810248E-2_wp,-1.5749517743615271E-2_wp, &
      &-1.4436296564231682E-2_wp,-8.0176628290002383E-3_wp,-7.3547355990133435E-3_wp, &
      &-1.3145371537682821E-2_wp,-1.0640115983678637E-2_wp,-9.6504448753380566E-3_wp, &
      &-4.4996101141431362E-3_wp,-1.2740992855028960E-2_wp,-8.3393303059628495E-3_wp, &
      &-3.7212354276045315E-2_wp,-3.2852410384024519E-2_wp,-3.0614546491936758E-2_wp, &
      &-1.4346949783307351E-2_wp,-8.2244399711126057E-2_wp,-1.6294715757294374E-2_wp, &
      &-3.9038610656489593E-2_wp,-2.0567213895689281E-2_wp,-2.1558424406190164E-2_wp, &
      &-1.9336394133874998E-2_wp,-3.7332351753141214E-3_wp,-5.8233052775258781E-3_wp, &
      &-4.8268739519799772E-3_wp,-4.5927344991332809E-2_wp,-1.9723672704756907E-2_wp, &
      &-3.7824098003146032E-2_wp,-5.1862651268351767E-3_wp,-5.7275846053051168E-3_wp, &
      &-6.4146049302302843E-3_wp,-3.6643019717926037E-2_wp,-1.2627827031284736E-2_wp, &
      &-1.0353655269088267E-2_wp,-1.4237858050055431E-2_wp,-7.3827880103746520E-3_wp, &
      &-1.9570555815441794E-2_wp,-1.2421141218402177E-2_wp,-1.0448640892029166E-2_wp, &
      &-1.6329032603068817E-2_wp,-2.2071530109002500E-2_wp,-1.3542982938032555E-2_wp, &
      &-1.1650410265984359E-2_wp,-1.7245346068721114E-2_wp,-2.8422909990177846E-2_wp]

   call get_structure(mol, "MB16-43", "17")
   do ii = 1, size(func)
      call get_rational_damping(inp, trim(func(ii)), error, s9=1.0_wp)
      if (allocated(error)) return
      call new_rational_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3bjatm_mb17


subroutine test_d3zeroatm_mb25(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(zero_damping_param) :: param
   type(d3_param) :: inp
   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      & "slaterdiracexchange", "blyp", "bp", "b97d", "revpbe", "pbe", "pbesol", &
      & "rpw86pbe", "rpbe", "tpss", "b3lyp", "pbe0", "hse06", "revpbe38", &
      & "pw6b95", "tpss0", "b2plyp", "pwpb95", "b2gpplyp", "ptpss", "hf", &
      & "mpwlyp", "bpbe", "bhlyp", "tpssh", "pwb6k", "b1b95", "bop", "olyp", &
      & "opbe", "ssb", "revssb", "otpss", "b3pw91", "revpbe0", "pbe38", &
      & "mpw1b95", "mpwb1k", "bmk", "camb3lyp", "lcwpbe", "m05", "m052x", &
      & "m06l", "m06", "m062x", "m06hf", "hcth120", "scan", "wb97x", &
      & "mpw1lyp", "mpwkcis1k", "pkzb", "n12", "m08hx", "m11l", "mn15l", "pwp", &
      & "pw1pw", "pbehpbe", "xlyp", "mpwpw", "hcth407", "revtpss", "tauhcth", &
      & "b3p", "b1p", "b1lyp", "mpw1pw", "mpw1kcis", "pbeh1pbe", "pbe1kcis", &
      & "x3lyp", "o3lyp", "b971", "b972", "b98", "hiss", "hse03", "revtpssh", &
      & "revtpss0", "tpss1kcis", "tauhcthhyb", "mpw2plyp", "b973c", "cf22d"]
   real(wp), parameter :: ref(*) = [&
      & 1.0613154373646663E-1_wp,-1.8876259384422459E-2_wp,-1.7576377305669005E-2_wp, &
      &-2.3748702681871504E-2_wp,-2.2303369775790879E-2_wp,-8.6007220574423112E-3_wp, &
      &-6.1997715812495803E-3_wp,-9.8451693951057510E-3_wp,-2.2252291588569017E-2_wp, &
      &-1.2523754839208579E-2_wp,-1.5397269759240936E-2_wp,-9.1754056201685599E-3_wp, &
      &-5.8634194145658200E-3_wp,-1.5486782768094242E-2_wp,-7.1334089761414159E-3_wp, &
      &-1.2012602963321843E-2_wp,-8.3954510274247326E-3_wp,-5.7681788549171359E-3_wp, &
      &-6.0023822786917569E-3_wp,-7.0862693280077070E-3_wp,-1.7586391714237665E-2_wp, &
      &-1.1106929705171003E-2_wp,-2.1778345895411035E-2_wp,-1.2322535067990021E-2_wp, &
      &-1.2281827130435118E-2_wp,-4.4791870771283089E-3_wp,-1.4591798101832332E-2_wp, &
      &-2.9238129103606399E-2_wp,-3.7369678688182618E-2_wp,-3.6796456367160557E-2_wp, &
      &-8.1862174778819484E-3_wp,-7.2999555919895575E-3_wp,-1.6431140091829909E-2_wp, &
      &-1.7402804487684997E-2_wp,-1.8922181622277653E-2_wp,-9.2405814515867588E-3_wp, &
      &-8.9000357024468183E-3_wp,-8.3515200684305940E-3_wp,-1.6554317187182983E-2_wp, &
      &-1.0553671984532416E-2_wp,-1.1195082572888794E-2_wp,-5.8550361165388594E-3_wp, &
      &-1.0495090973696125E-3_wp,-4.4388328691468787E-4_wp,-1.7181960835395168E-3_wp, &
      &-3.6397340658888454E-4_wp,-9.0020511932907923E-4_wp,-1.2216272189061111E-2_wp, &
      &-1.7275291070818054E-3_wp,-6.1782653468711764E-3_wp,-1.4873624091587857E-2_wp, &
      &-1.4090547469564809E-2_wp,-6.2753850597721081E-2_wp,-1.9707922360970302E-2_wp, &
      &-3.5327512784430563E-4_wp,-8.4306186265710459E-3_wp, 5.3592259491208568E-5_wp, &
      &-6.6548263945738987E-3_wp,-9.6587859311542661E-3_wp,-1.1131535816615383E-2_wp, &
      &-1.9247819733383506E-2_wp,-1.6145567164365680E-2_wp,-2.1022405295279017E-2_wp, &
      &-1.1908888764015774E-2_wp,-1.8314527978666263E-2_wp,-1.2714111247643373E-2_wp, &
      &-1.2308733757439341E-2_wp,-1.6145567164365680E-2_wp,-1.3319121973193983E-2_wp, &
      &-1.7650135255476308E-2_wp,-9.2723328966778720E-3_wp,-1.3594718067989876E-2_wp, &
      &-1.2231360411976182E-2_wp,-1.4855506321386058E-2_wp,-1.2440947033067303E-2_wp, &
      &-1.2723713762009774E-2_wp,-1.4469421175303927E-2_wp,-7.4336230287781746E-3_wp, &
      &-8.9126407786282984E-3_wp,-1.1258613133535007E-2_wp,-1.0204599751079367E-2_wp, &
      &-1.6065929734220732E-2_wp,-1.3083826111391977E-2_wp,-6.0511412375980119E-3_wp, &
      &-1.8691177475853686E-2_wp,-5.7926845557257388E-4_wp]

   call get_structure(mol, "MB16-43", "25")
   do ii = 1, size(func)
      call get_zero_damping(inp, trim(func(ii)), error, s9=1.0_wp)
      if (allocated(error)) return
      call new_zero_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3zeroatm_mb25


subroutine test_d3bjm_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp

   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      "b2plyp", "b3lyp", "b97d", "blyp", "bp", "pbe", "pbe0", "lcwpbe"]
   real(wp), parameter :: ref(*) = [&
      &-2.4496062549935024E-2_wp,-4.7306507711299510E-2_wp,-5.9503371130070551E-2_wp, &
      &-5.6673386960027376E-2_wp,-4.7342638767330296E-2_wp,-2.4952684302621722E-2_wp, &
      &-2.3521898974796757E-2_wp,-2.7484708182531545E-2_wp]

   call get_structure(mol, "MB16-43", "02")
   do ii = 1, size(func)
      call get_rational_damping(inp, trim(func(ii)), error)
      if (allocated(error)) exit
      call new_rational_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3bjm_mb02


subroutine test_d3zerom_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(mzero_damping_param) :: param
   type(d3_param) :: inp
   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      "b2plyp", "b3lyp", "b97d", "blyp", "bp", "pbe", "pbe0", "lcwpbe"]
   real(wp), parameter :: ref(*) = [&
      &-1.8298642060887996E-2_wp,-3.3517506681341125E-2_wp,-8.3081078831892127E-2_wp, &
      &-4.1380833998837407E-2_wp,-2.2212534670402604E-2_wp,-7.4973021425312925E-2_wp, &
      &-5.7366087952954614E-2_wp,-1.6369492826335569E-2_wp]


   call get_structure(mol, "MB16-43", "03")
   do ii = 1, size(func)
      call get_mzero_damping(inp, trim(func(ii)), error)
      if (allocated(error)) return
      call new_mzero_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3zerom_mb03


subroutine test_d3bjmatm_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(rational_damping_param) :: param
   type(d3_param) :: inp

   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      "b2plyp", "b3lyp", "b97d", "blyp", "bp", "pbe", "pbe0", "lcwpbe"]
   real(wp), parameter :: ref(*) = [&
      &-1.8210582708938286E-2_wp,-4.9517667234694238E-2_wp,-1.1799780932991713E-1_wp, &
      &-5.0481997336614820E-2_wp,-2.1296302814694204E-2_wp,-4.2996498228542929E-2_wp, &
      &-4.2021630904794396E-2_wp,-1.8994831113760144E-2_wp]


   call get_structure(mol, "MB16-43", "04")
   do ii = 1, size(func)
      call get_mrational_damping(inp, trim(func(ii)), error, s9=1.0_wp)
      if (allocated(error)) exit
      call new_rational_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3bjmatm_mb04


subroutine test_d3zeromatm_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(mzero_damping_param) :: param
   type(d3_param) :: inp
   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      & "b2plyp", "b3lyp", "b97d", "blyp", "bp", "pbe", "pbe0", "lcwpbe"]
   real(wp), parameter :: ref(*) = [&
      &-2.0123625206829979E-2_wp,-3.7721574372867127E-2_wp,-7.8081589749874739E-2_wp, &
      &-4.6609937658477547E-2_wp,-2.6471019009804603E-2_wp,-1.2081305036647522E-1_wp, &
      &-8.5996131518298685E-2_wp,-1.9517670998784716E-2_wp]

   call get_structure(mol, "MB16-43", "05")
   do ii = 1, size(func)
      call get_mzero_damping(inp, trim(func(ii)), error, s9=1.0_wp)
      if (allocated(error)) return
      call new_mzero_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3zeromatm_mb05


subroutine test_d3op_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(optimizedpower_damping_param) :: param
   type(d3_param) :: inp

   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      & "pbe", "pbe0", "revtpss", "revtpssh", "blyp", "b3lyp", "b97d", "b3lyp", "b971", &
      & "revpbe", "revpbe0", "tpss", "tpssh", "ms2", "ms2h"]
   real(wp), parameter :: ref(*) = [&
      &-8.5891879356226816E-3_wp,-1.1189585521489914E-2_wp,-6.6642513192637625E-3_wp, &
      &-6.6836852162675785E-3_wp,-2.8029636019647405E-2_wp,-1.6569059351798480E-2_wp, &
      &-5.0045344281734523E-2_wp,-1.6569059351798480E-2_wp,-2.0411405825821769E-2_wp, &
      &-4.9647291366122226E-2_wp,-3.1170803546237417E-2_wp,-8.6007550664757206E-3_wp, &
      &-8.2556637816396158E-3_wp,-3.9408761597590552E-3_wp,-7.0002355523013822E-3_wp]

   call get_structure(mol, "MB16-43", "06")
   do ii = 1, size(func)
      call get_optimizedpower_damping(inp, trim(func(ii)), error)
      if (allocated(error)) exit
      call new_optimizedpower_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3op_mb06


subroutine test_d3opatm_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(optimizedpower_damping_param) :: param
   type(d3_param) :: inp
   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      &"pbe", "pbe0", "revtpss", "revtpssh", "blyp", "b3lyp", "b97d", "b3lyp", "b971", &
      &"revpbe", "revpbe0", "tpss", "tpssh", "ms2", "ms2h"]
   real(wp), parameter :: ref(*) = [&
      &-2.2127517310697481E-2_wp,-2.8136226523866572E-2_wp,-1.2454661288657357E-2_wp, &
      &-1.3483016246834175E-2_wp,-6.4438026481408256E-2_wp,-4.2243577071262313E-2_wp, &
      &-9.2952134273516937E-2_wp,-4.2243577071262313E-2_wp,-4.6977205030535207E-2_wp, &
      &-9.2179068206361112E-2_wp,-5.6249707420366557E-2_wp,-1.8481081447338607E-2_wp, &
      &-1.7632416675380274E-2_wp,-7.8454173753530632E-3_wp,-1.4068302146268733E-2_wp]


   call get_structure(mol, "MB16-43", "05")
   do ii = 1, size(func)
      call get_optimizedpower_damping(inp, trim(func(ii)), error, s9=1.0_wp)
      if (allocated(error)) return
      call new_optimizedpower_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3opatm_mb07


subroutine test_d3bj_unknown(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(d3_param) :: inp

   call get_rational_damping(inp, "unknown", error)

end subroutine test_d3bj_unknown


subroutine test_d3zero_unknown(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(d3_param) :: inp

   call get_zero_damping(inp, "unknown", error)

end subroutine test_d3zero_unknown


subroutine test_d3bjm_unknown(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(d3_param) :: inp

   call get_mrational_damping(inp, "unknown", error)

end subroutine test_d3bjm_unknown


subroutine test_d3zerom_unknown(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(d3_param) :: inp

   call get_mzero_damping(inp, "unknown", error)

end subroutine test_d3zerom_unknown


subroutine test_d3op_unknown(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(d3_param) :: inp

   call get_optimizedpower_damping(inp, "unknown", error)

end subroutine test_d3op_unknown


subroutine test_d3cso_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(cso_damping_param) :: param
   type(d3_param) :: inp

   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      & "blyp", "bp86", "pbe", "tpss", "b3lyp", "pbe0", "pw6b95", "b2plyp"]
   real(wp), parameter :: ref(*) = [&
      &-4.6561914861742063E-2_wp,-4.1059723690412941E-2_wp,-2.5368289609215090E-2_wp, &
      &-3.5149962802689075E-2_wp,-3.8002950817452322E-2_wp,-2.4553150176425589E-2_wp, &
      &-1.7420680139517471E-2_wp,-1.9839377295845999E-2_wp]

   call get_structure(mol, "MB16-43", "01")
   do ii = 1, size(func)
      call get_cso_damping(inp, trim(func(ii)), error)
      if (allocated(error)) exit
      call new_cso_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3cso_mb01


subroutine test_d3csoatm_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(cso_damping_param) :: param
   type(d3_param) :: inp

   integer :: ii
   character(len=*), parameter :: func(*) = [character(len=20)::&
      & "blyp", "bp86", "pbe", "tpss", "b3lyp", "pbe0", "pw6b95", "b2plyp"]
   real(wp), parameter :: ref(*) = [&
      &-7.0301894479094934E-2_wp,-6.1975513588852130E-2_wp,-3.8229908827789323E-2_wp, &
      &-5.3032363743776519E-2_wp,-5.7349746427606110E-2_wp,-3.6996370918123717E-2_wp, &
      &-2.6202914208549720E-2_wp,-2.9882002966364205E-2_wp]

   call get_structure(mol, "MB16-43", "02")
   do ii = 1, size(func)
      call get_cso_damping(inp, trim(func(ii)), error, s9=1.0_wp)
      if (allocated(error)) return
      call new_cso_damping(param, inp)
      call test_dftd3_gen(error, mol, param, ref(ii))
      if (allocated(error)) exit
   end do

end subroutine test_d3csoatm_mb02


subroutine test_d3cso_unknown(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(d3_param) :: inp

   call get_cso_damping(inp, "unknown", error)

end subroutine test_d3cso_unknown


end module test_param
