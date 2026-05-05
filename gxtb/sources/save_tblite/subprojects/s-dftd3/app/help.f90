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

module dftd3_app_help
   use dftd3, only : get_dftd3_version
   implicit none
   private

   public :: prog_name, header, version
   public :: help_text, run_help_text, param_help_text, gcp_help_text


   character(len=*), parameter :: prog_name = "s-dftd3"

   character(len=*), parameter :: nl = new_line('a')

   character(len=*), parameter :: run_options_text = &
      "-i,--input <format>      Hint for the format of the input file"//nl//&
      "   --bj <method>         Use rational (Becke-Johnson) damping function"//nl//&
      "   --bj-param <list>     Specify parameters for rational damping,"//nl//&
      "                         expected order is s6, s8, a1, a2 (requires four arguments)"//nl//&
      "   --zero <method>       Use zero (Chai-Head-Gordon) damping function"//nl//&
      "   --zero-param <list>   Specify parameters for zero damping,"//nl//&
      "                         expected order is s6, s8, rs6 (requires three arguments)"//nl//&
      "   --bjm <method>        Use modified rational damping function"//nl//&
      "   --bjm-param <list>    Specify parameters for rational damping,"//nl//&
      "                         expected order is s6, s8, a1, a2 (requires four arguments)"//nl//&
      "   --zerom <method>      Use modified zero damping function"//nl//&
      "   --zerom-param <list>  Specify parameters for modified zero damping,"//nl//&
      "                         expected order is s6, s8, rs6, bet (requires four arguments)"//nl//&
      "   --op <method>         Use optimized power damping function"//nl//&
      "   --op-param <list>     Specify parameters for optimized power,"//nl//&
      "                         expected order is s6, s8, a1, a2, bet (requires five arguments)"//nl//&
      "   --cso <method>        Use CSO (C6-scaled only) damping function"//nl//&
      "   --cso-param <list>    Specify parameters for CSO damping,"//nl//&
      "                         expected order is s6, a1 (requires two arguments)"//nl//&
      "   --atm                 Use ATM three-body dispersion"//nl//&
      "   --atm-scale <s9>      Use scaled ATM three-body dispersion"//nl//&
      "   --gcp <basis>         Include geometric counter-poise correction for given basis set"//nl//&
      "   --db <file>           Load parameters from external data file"//nl//&
      "   --noedisp             Disable writing of dispersion energy to .EDISP file"//nl//&
      "   --json [file]         Dump results to JSON output (default: dftd3.json)"//nl//&
      "   --grad [file]         Request gradient evaluation,"//nl//&
      "                         write results to file (default: dftd3.txt),"//nl//&
      "                         attempts to add to Turbomole gradient and gradlatt files"//nl//&
      "   --property            Evaluate dispersion related properties"//nl//&
      "   --pair-resolved       Calculate pairwise representation of dispersion energy"//nl//&
      "   --citation [file]     Print citation information to file (default: dftd3.bib)"//nl//&
      "-v,--verbose             Show more, can be used multiple times"//nl//&
      "-s,--silent              Show less, use twice to supress all output"//nl//&
      "   --version             Print program version and exit"//nl//&
      "   --help                Show this help message"

   character(len=*), parameter :: run_help_text = &
      "Usage: "//prog_name//" [run] [options] <input>"//nl//&
      ""//nl//&
      "Takes an geometry input to calculate the D3 dispersion correction."//nl//&
      "Periodic calculations are performed automatically for periodic input formats."//nl//&
      "Specify the functional to select the correct parameters."//nl//&
      ""//nl//&
      run_options_text//nl//&
      ""

   character(len=*), parameter :: param_help_text = &
      "Usage: "//prog_name//" param [options] <input> [method] [damping]"//nl//&
      ""//nl//&
      "Takes a damping parameter data file and performs queries for damping"//nl//&
      "parameters if a method is provided, if no damping function is provided"//nl//&
      "the default damping functions as provided in the data file will be used."//nl//&
      "The data file is provided in TOML format."//nl//&
      ""//nl//&
      "Example:"//nl//&
      ""//nl//&
      "    [default]"//nl//&
      "    d3 = [""bj"", ""zero""]"//nl//&
      ""//nl//&
      "    [default.parameter]"//nl//&
      "    d3.bj = {s6=1.0, s9=0.0, alp=14.0, damping=""rational""}"//nl//&
      "    d3.zero = {s6=1.0, s9=0.0, rs8=1.0, alp=14.0, damping=""zero""}"//nl//&
      "    d3.bjm = {s6=1.0, s9=0.0, alp=14.0, damping=""rational""}"//nl//&
      "    d3.zerom = {s6=1.0, s9=0.0, rs8=1.0, alp=14.0, damping=""mzero""}"//nl//&
      "    d3.op = {s9=0.0, alp=14.0, damping=""optimizedpower""}"//nl//&
      "    d3.cso = {s6=1.0, s9=0.0, a2=2.5, rs6=0.0, rs8=6.25, alp=14.0, damping=""cso""}"//nl//&
      ""//nl//&
      "    [parameter.bp]"//nl//&
      "    d3.bj = {a1=0.3946, s8=3.2822, a2=4.8516}"//nl//&
      "    d3.zero = {rs6=1.139, s8=1.683}"//nl//&
      "    d3.bjm = {a1=0.821850, s8=3.140281, a2=2.728151}"//nl//&
      "    d3.zerom = {rs6=1.233460, s8=1.945174, bet=0.000000}"//nl//&
      ""//nl//&
      "    [parameter.blyp]"//nl//&
      "    d3.bj = {a1=0.4298, s8=2.6996, a2=4.2359}"//nl//&
      "    d3.zero = {rs6=1.094, s8=1.682}"//nl//&
      "    d3.bjm = {a1=0.448486, s8=1.875007, a2=3.610679}"//nl//&
      "    d3.zerom = {rs6=1.279637, s8=1.841686, bet=0.014370}"//nl//&
      "    d3.op = {s6=1.0, s8=1.31867, a1=0.425, a2=3.50, bet=2.0}"//nl//&
      ""//nl//&
      "    [parameter.revpbe]"//nl//&
      "    d3.bj = {a1=0.5238, s8=2.3550, a2=3.5016}"//nl//&
      "    d3.zero = {rs6=0.923, s8=1.010}"//nl//&
      "    d3.op = {s6=1.0, s8=1.44765, a1=0.600, a2=2.50, bet=0.0}"//nl//&
      ""

   character(len=*), parameter :: gcp_help_text = &
      "Usage: "//prog_name//" gcp [options] <input>"//nl//&
      ""//nl//&
      "Takes an geometry input to calculate the geometric counter-poise correction."//nl//&
      "Periodic calculations are performed automatically for periodic input formats."//nl//&
      "Specify the level of theory to select the correct parameters."//nl//&
      ""//nl//&
      "-i,--input <format>      Hint for the format of the input file"//nl//&
      "-l,--level <method>[/<basis>]"//nl//&
      "                         Level of theory, basis set is inferred for 3c methods"//nl//&
      "   --nocpc               Disable writing of dispersion energy to .CPC file"//nl//&
      "   --json [file]         Dump results to JSON output (default: gcp.json)"//nl//&
      "   --grad [file]         Request gradient evaluation,"//nl//&
      "                         write results to file (default: gcp.txt),"//nl//&
      "                         attempts to add to Turbomole gradient and gradlatt files"//nl//&
      "-v,--verbose             Show more, can be used multiple times"//nl//&
      "-s,--silent              Show less, use twice to supress all output"//nl//&
      "   --version             Print program version and exit"//nl//&
      "   --help                Show this help message"//nl//&
      ""

   character(len=*), parameter :: help_text = &
      "Usage: "//prog_name//" [run|param|gcp] [options] ..."//nl//&
      ""//nl//&
      "Commands"//nl//&
      ""//nl//&
      "  run       Evaluate dispersion correction on the provided input structure."//nl//&
      "            Periodic calculations are performed automatically for periodic inputs"//nl//&
      "            If no command is specified run is selected by default."//nl//&
      ""//nl//&
      "  param     Inspect and manipulate damping parameter data file."//nl//&
      ""//nl//&
      "Options"//nl//&
      ""//nl//&
      run_options_text//nl//&
      ""


contains


subroutine header(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_dftd3_version(string=version_string)
   write(unit, '(a)') &
      "-----------------------------------", &
      " s i m p l e   D F T - D 3  v"// version_string, &
      "-----------------------------------", ""

end subroutine header


subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_dftd3_version(string=version_string)
   write(unit, '(a, *(1x, a))') &
      & prog_name, "version", version_string

end subroutine version


end module dftd3_app_help
