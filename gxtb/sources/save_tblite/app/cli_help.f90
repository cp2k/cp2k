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

!> Module to collect all help printout for the tblite program
module tblite_cli_help
   implicit none
   private

   public :: prog_name, help_text, help_text_run, help_text_param, help_text_fit, &
      & help_text_tagdiff, help_text_guess

   character(len=*), parameter :: nl = new_line('a')

   !> Name of the main executable
   character(len=*), parameter :: prog_name = "tblite"

   !> General help text regarding response file support
   character(len=*), parameter :: help_text_response = &
      "Command line arguments can be read from an indirect file / response file"//nl//&
      "by specifying the file with @name in the command line. Each line in the file"//nl//&
      "is interpreted as command line argument, shell like escape sequences are not"//nl//&
      "available. The file can contain further @name inputs. If the file cannot be"//nl//&
      "found the argument is used literally."

   !> General help text for version and help commands
   character(len=*), parameter :: help_text_general = &
      "      --version           Print program version and exit"//nl//&
      "      --help              Show this help message"

   !> General help text regarding geometry input formats
   character(len=*), parameter :: help_text_geometry = &
      "Supported geometry input formats are:"//nl//&
      ""//nl//&
      "- Xmol/xyz files (xyz, log)"//nl//&
      "- Turbomole's coord, riper's periodic coord (tmol, coord)"//nl//&
      "- DFTB+ genFormat geometry inputs as cluster, supercell or fractional (gen)"//nl//&
      "- VASP's POSCAR/CONTCAR input files (vasp, poscar, contcar)"//nl//&
      "- Protein Database files, only single files (pdb)"//nl//&
      "- Connection table files, molfile (mol) and structure data format (sdf)"//nl//&
      "- Gaussian's external program input (ein)"

   !> Help text for main driver
   character(len=*), parameter :: help_text = &
      "Usage: "//prog_name//" <run|guess|param|fit|tagdiff> [options]"//nl//&
      ""//nl//&
      "Commands"//nl//&
      ""//nl//&
      "  run       Evaluate tight-binding module on the provided input structure."//nl//&
      "            If no command is specified run is selected by default."//nl//&
      ""//nl//&
      "  guess     Execute only the atomic charge or density guess, respectively."//nl//&
      "            Prints the initializing charges and dipole moments."//nl//&
      ""//nl//&
      "  param     Inspect and manipulate tight-binding parametrization data."//nl//&
      ""//nl//&
      "  fit       Optimize tight-binding parameters."//nl//&
      ""//nl//&
      "  tagdiff   Auxilary tool to compute differences between data outputs."//nl//&
      ""//nl//&
      "Options"//nl//&
      ""//nl//&
      "      --color string      Support colorful terminal output,"//nl//&
      "                          possible values are never (default), always, auto."//nl//&
      help_text_general//nl//&
      ""//nl//&
      help_text_response

   !> Help text for param command
   character(len=*), parameter :: help_text_param = &
      "Usage: "//prog_name//" param [options] <input>"//nl//&
      ""//nl//&
      "Import, export and manipulate tight-binding parameter files."//nl//&
      "An input parameter file must be present as positional argument or"//nl//&
      "an internal parametrization should be selected."//nl//&
      ""//nl//&
      "Options"//nl//&
      ""//nl//&
      "      --method <name>     Base parametrization of the xTB Hamiltonian to use"//nl//&
      "                          Available methods: gfn1, gfn2, gxtb, ipea1"//nl//&
      "      --output <file>     Output file for writing parametrization"//nl//&
      help_text_general//nl//&
      ""//nl//&
      help_text_response

   !> Help text for run command
   character(len=*), parameter :: help_text_run = &
      "Usage: "//prog_name//" [run] [options] <input>"//nl//&
      ""//nl//&
      "Evaluates the tight-binding model on the provided input structure."//nl//&
      "Periodic calculations are performed automatically for periodic input formats."//nl//&
      "Reads .CHRG and .UHF files if present in the same directory as the input."//nl//&
      ""//nl//&
      help_text_geometry//nl//&
      ""//nl//&
      "Options"//nl//&
      ""//nl//&
      "  -c, --charge <real>      Set charge to molecule, overwrites .CHRG file"//nl//&
      "      --spin <int>         Set number of unpaired electrons, overwrites .UHF file"//nl//&
      "      --method <name>      Parametrization of the xTB Hamiltonian to use"//nl//&
      "                           Available methods: gfn1, gfn2, gxtb, ipea1 (Default: gfn2)"//nl//&
      "      --param <file>       Parametrization file to use for calculation"//nl//&
      "      --acc <real>         Convergence criterion for SCF (Default: 1.0)"//nl//&
      "      --etemp <real>       Electronic temperature for calculation (Default: 300K for GFN, 0K for g-xTB)"//nl//&
      "      --guess <name>       Guess for the initial populations, possible options:"//nl//&
      "                           sad (default), eeq, eeqbc (default for g-xTB), ceh (Charge-Extended Hückel method)"//nl//&
      "      --restart [<file>]   Restart calculation from previous run"//nl//&
      "                           File name for restart file is optional,"//nl//&
      "                           (Default: tblite-restart.npz)"//nl//&
      "      --no-restart         Do not restart calculation from previous run"//nl//&
      "      --etemp-guess <real> Electronic temperature for ceh-guess (Default: 5000K)"//nl//&
      "      --iterations <int>   Maximum number of SCF iterations (Default: 250)"//nl//&
      "      --broyden-start <int>"//nl//&
      "                           Iteration to start the Broyden mixing (Default: 0)"//nl//&
      "      --diis-start <int>   Iteration to start the DIIS mixing (Default: 0)"//nl//&
      "      --solver <name>      Electronic solver for SCF, possible options:"//nl//&
      "                           gvd (default), and gvr"//nl//&
      "      --efield <real>,<real>,<real>"//nl//&
      "                           Homogeneous electric field in V/Å."//nl//&
      "      --alpb <name>        Use analytical linearized Poisson-Boltzmann (ALPB) solvation model."//nl//&
      "                           Solvent is specified by the solvent name."//nl//&
      "                           Uses parametrized ALPB with CDS and empirical shift."//nl//&
      "      --gbsa <name>        Use generalized Born solvation model (GBSA)."//nl//&
      "                           Solvent is specified by the solvent name."//nl//&
      "                           Uses parametrized GB with CDS and empirical shift."//nl//&
      "      --gbe <real>/<name>  Use generalized Born for finite epsilion (GBE) solvation model."//nl//&
      "                           Solvent is specified by dielectric constant or the solvent name."//nl//&
      "      --gb <real>/<name>   Use generalized Born solvation model (GB)."//nl//&
      "                           Solvent is specified by dielectric constant or solvent name."//nl//&
      "      --cpcm <real>/<name> Use polarizable continuum solvation model (CPCM)."//nl//&
      "                           Solvent is specified by dielectric constant or solvent name."//nl//&
      "      --born-kernel <name> Specify Born kernel to use with ALPB, GBSA or GB solvation model."//nl//&
      "                           Possible options are p16 (default for ALPB) and still (default for GB/GBSA)."//nl//&
      "      --solv-state <name>  Solution state correction: gsolv (default), bar1mol, reference."//nl//& 
      "      --spin-polarized     Use spin-polarized xTB Hamiltonian"//nl//&
      "      --post-processing <file>"//nl//&
      "                           Add post processing methods to the calculation"//nl//&
      "                           by using a toml file as input."//nl//&
      "      --post-processing <name>"//nl//&
      "                           Add post processing methods to the calculation,"//nl//&
      "                           Mayer-Wiberg bond orders are computed by default."//nl//&
      "                           Options: molmom (molecular multipole moments)"//nl//&
      "                           Options: xtbml (atomistic properties based on Mulliken partitioning)"//nl//& 
      "      --post-processing-output <file>"//nl//&
      "                           Output file for post processing results in npz format."//nl//&
      "                           (default: tblite-data.npz)"//nl//&
      "      --grad [file]        Evaluate molecular gradient and virial"//nl//&
      "                           Results are stored in file (default: tblite.txt)"//nl//&
      "      --json [file]        Dump results as JSON output (default: tblite.json)"//nl//&
      "  -i, --input <format>     Hint for the format of the input file"//nl//&
      "  -v, --verbose            Increase verbosity of printout"//nl//&
      "  -s, --silent             Reduce verbosity of printout"//nl//&
      help_text_general//nl//&
      ""//nl//&
      help_text_response

   !> Help text for guess command
   character(len=*), parameter :: help_text_guess = &
      "Usage: "//prog_name//" [guess] [options] <input>"//nl//&
      ""//nl//&
      "Execute only the atomic charge / density guess."//nl//&
      "Prints the initializing charges and dipole moments."//nl//&
      "Reads .CHRG and .UHF files if present in the same directory as the input."//nl//&
      ""//nl//&
      help_text_geometry//nl//&
      ""//nl//&
      "Options"//nl//&
      ""//nl//&
      "  -c, --charge <real>      Set charge to molecule, overwrites .CHRG file"//nl//&
      "      --spin <int>         Set number of unpaired electrons, overwrites .UHF file"//nl//&
      "      --method <name>      Guess for the initial populations, possible options:"//nl//&
      "                           sad, eeq, eeqbc, ceh (Charge-Extended Hückel method, default)"//nl//&
      "      --etemp-guess <real> Electronic temperature for ceh-guess (Default: 4000K)"//nl//&
      "      --solver <name>      Electronic solver for SCF, possible options:"//nl//&
      "                           gvd (default), and gvr"//nl//&
      "      --efield <real>,<real>,<real>"//nl//&
      "                           Homogeneous electric field in V/Å."//nl//&
      "      --grad               Evaluate analytic gradient of charges."//nl//&
      "      --json [file]        Dump results as JSON output (default: tblite.json)"//nl//&
      "  -i, --input <format>     Hint for the format of the input file"//nl//&
      "  -v, --verbose            Increase verbosity of printout"//nl//&
      "  -s, --silent             Reduce verbosity of printout"//nl//&
      help_text_general//nl//&
      ""//nl//&
      help_text_response

   !> Help text for fit command
   character(len=*), parameter :: help_text_fit = &
      "Usage: "//prog_name//" fit [options] <param-file> <input-file>"//nl//&
      ""//nl//&
      "Takes the name of the starting parameter file as first positional argument and"//nl//&
      "an input file for the settings of the run as second positional argument."//nl//&
      "A starting parameter file can be produced using the tblite-param(1) command."//nl//&
      ""//nl//&
      "Options"//nl//&
      ""//nl//&
      "      --dry-run           Do not run start the fitting procedure"//nl//&
      "      --copy <file>       Write the full representation of the input to <file>,"//nl//&
      "                          all defaults will be filled in and the mask expanded"//nl//&
      "  -v, --verbose           Increase verbosity of printout"//nl//&
      "  -s, --silent            Reduce verbosity of printout"//nl//&
      help_text_general//nl//&
      ""//nl//&
      help_text_response

   !> Help text for tagdiff command
   character(len=*), parameter :: help_text_tagdiff = &
      "Usage: "//prog_name//" tagdiff [options] <actual> <reference>"//nl//&
      ""//nl//&
      "Takes two positional arguments and computes the difference between the"//nl//&
      "entries in the reference and the actual data files. Only keys of the"//nl//&
      "reference data file are used from the actual one."//nl//&
      ""//nl//&
      "Options"//nl//&
      ""//nl//&
      "      --fit               Produce output for fit command"//nl//&
      help_text_general//nl//&
      ""//nl//&
      help_text_response

end module tblite_cli_help
