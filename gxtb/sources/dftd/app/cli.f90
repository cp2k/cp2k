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

!> Definition of the command line interface
module dftd4_cli
   use, intrinsic :: iso_fortran_env, only : output_unit
   use mctc_env, only : error_type, fatal_error, get_argument, wp
   use mctc_io, only : get_filetype
   use dftd4, only : param_type, get_dftd4_version, twobody_damping_function, &
      & threebody_damping_function
   use dftd4_argument, only : argument_list, len
   use dftd4_help, only : citation, header, help_text, help_text_param, &
      & help_text_run, license, prog_name, version

   implicit none
   private

   public :: cli_config, param_config, run_config, get_arguments


   !> Base command line configuration
   type, abstract :: cli_config
   end type cli_config

   !> Configuration data for running single point calculations
   type, extends(cli_config) :: run_config
      !> Name of input geometry file
      character(len=:), allocatable :: input
      !> Format of the geometry input
      integer, allocatable :: input_format
      !> Damping parameters
      type(param_type) :: inp
      !> Name of the functional/method used
      character(len=:), allocatable :: method
      !> Dispersion model name
      character(len=:), allocatable :: model
      !> Flag to write JSON output
      logical :: json = .false.
      !> Name of JSON output file
      character(len=:), allocatable :: json_output
      !> Flag to wrap coordinates to central cell
      logical :: wrap = .true.
      !> Flag to write the disperison energy to the .EDISP file
      logical :: tmer = .true.
      !> Flag to output dispersion related properties
      logical :: properties = .false.
      !> Flag to use ATM three-body dispersion
      logical :: atm = .true.
      !> Flag to change the default charge scaling parameters
      logical :: zeta = .false.
      !> Flag to activate gradient calculation
      logical :: grad = .false.
      !> Flag to activate numerical hessian calculation
      logical :: hessian = .false.
      !> Name of gradient output file
      character(len=:), allocatable :: grad_output
      !> Identifier for the used two-body damping function
      integer :: damping_2b = -1
      !> Identifier for the used three-body damping function
      integer :: damping_3b = -1
      !> Flag to indicate the presence of custom damping parameters
      logical :: has_param = .false.
      !> Verbosity level
      integer :: verbosity = 2
      !> System charge
      real(wp), allocatable :: charge
      !> Charge scaling height for D4
      real(wp) :: ga = 3.0_wp
      !> Charge scaling steepness for D4
      real(wp) :: gc = 2.0_wp
      !> Weighting factor for the CN interpolation in D4
      real(wp) :: wf = 6.0_wp
      !> Name of the used charge model
      character(len=:), allocatable :: charge_model
      !> Flag to calculate pairwise dispersion energies
      logical :: pair_resolved = .false.
   end type run_config
   type, extends(cli_config) :: param_config
      ! Flag to list all supported functionals
      logical :: list = .false.
   end type param_config

contains


subroutine get_arguments(config, error)

   !> Configuation data
   class(cli_config), allocatable, intent(out) :: config

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(argument_list) :: list
   integer :: iarg, narg
   character(len=:), allocatable :: arg

   iarg = 0
   list = argument_list()
   narg = len(list)
   do while(iarg < narg)
      iarg = iarg + 1
      call list%get(iarg, arg)
      select case(arg)
      case("-h", "--help")
         write(output_unit, '(a)') help_text
         stop
      case("--version")
         call version(output_unit)
         stop
      case("--citation")
         call citation(output_unit)
         stop
      case("--license")
         call license(output_unit)
         stop
      case default
         iarg = iarg - 1
         allocate(run_config :: config)
         exit
      case("run")
         allocate(run_config :: config)
         exit
      case("param")
         allocate(param_config :: config)
         exit
      end select
   end do
   if (allocated(error)) return

   if (.not.allocated(config)) then
      write(output_unit, '(a)') help_text
      call fatal_error(error, "Insufficient arguments provided")
      return
   end if

   select type(config)
   type is(run_config)
      call get_run_arguments(config, list, iarg, error)
   type is(param_config)
      call get_param_arguments(config, list, iarg, error)
   end select

end subroutine get_arguments


!> Read configuration for the single point driver
subroutine get_run_arguments(config, list, start, error)

   !> Configuation data
   type(run_config), intent(out) :: config

   !> List of command line arguments
   type(argument_list), intent(in) :: list

   !> First command line argument
   integer, intent(in) :: start

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   logical :: getopts
   character(len=:), allocatable :: arg

   iarg = start
   getopts = .true.
   narg = len(list)
   do while(iarg < narg)
      iarg = iarg + 1
      call list%get(iarg, arg)
      if (.not.getopts) then
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      end if
      select case(arg)
      case("--")
         getopts = .false.
      case("-h", "--help")
         write(output_unit, '(a)') help_text_run
         stop
      case("-v", "--verbose")
         config%verbosity = config%verbosity + 1
      case("-s", "--silent")
         config%verbosity = config%verbosity - 1
      case default
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         if (arg(1:1) == "-") then
            call fatal_error(error, "Unknown command option '"//arg//"'")
         else
            call fatal_error(error, "Too many positional arguments present")
         end if
         exit
      case("-m", "--model") 
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for model")
            exit
         end if
         call move_alloc(arg, config%model)
      case("-i", "--input")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         config%input_format = get_filetype("."//arg)
      case("-c", "--charge")
         iarg = iarg + 1
         allocate(config%charge)
         call get_argument_as_real(iarg, config%charge, error)
         if (allocated(error)) exit
      case("--json")
         config%json = .true.
         config%json_output = "dftd4.json"
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%json_output)
         end if
      case("--property")
         config%properties = .true.
      case("--pair-resolved")
         config%pair_resolved = .true.
      case("--noedisp")
         config%tmer = .false.
      case("--nowrap")
         config%wrap = .false.
      case("-g", "--grad")
         config%grad = .true.
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%grad_output)
         end if
      case("--hessian")
         config%grad = .true.
         config%hessian = .true.
      case("--noatm")
         config%atm = .false.
      case("--zeta")
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%ga, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%gc, error)
         if (allocated(error)) exit
         config%zeta = .true.
      case("--wfactor")
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%wf, error)
         if (allocated(error)) exit
      case("--qmodel")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for charge model")
            exit
         end if
         call move_alloc(arg, config%charge_model)
      case("-f", "--func")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
         call move_alloc(arg, config%method)
      case("-2b", "--damp-2b")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for two-body damping function")
            exit
         end if
         ! Select the two-body damping function
         select case(arg)
         case default
            call fatal_error(error, "Unknown two-body damping function '"//arg//"' specified")
            exit
         case("bj", "rational")
            config%damping_2b = twobody_damping_function%rational
         case("screened")
            config%damping_2b = twobody_damping_function%rational
         case("zero")
            config%damping_2b = twobody_damping_function%zero
         case("mzero")
            config%damping_2b = twobody_damping_function%mzero
         case("op")
            config%damping_2b = twobody_damping_function%optpower
         case("cso")
            config%damping_2b = twobody_damping_function%cso
         case("koide")
            config%damping_2b = twobody_damping_function%koide
         end select
      case("-3b", "--damp-3b")
         config%atm = .true.
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for three-body damping function")
            exit
         end if
         ! Select the three-body damping function
         select case(arg)
         case default
            call fatal_error(error, "Unknown three-body damping function '"//arg//"' specified")
            exit
         case("bj", "rational")
            config%damping_3b = threebody_damping_function%rational
         case("screened")
            config%damping_3b = threebody_damping_function%screened
         case("zero")
            config%damping_3b = threebody_damping_function%zero
         case("zero-avg")
            config%damping_3b = threebody_damping_function%zero_avg
         end select
      case("--s6")
         config%has_param = .true.
         iarg = iarg + 1
         allocate(config%inp%s6)
         call get_argument_as_real(iarg, config%inp%s6, error)
         if (allocated(error)) exit
      case("--s8")
         config%has_param = .true.
         iarg = iarg + 1
         allocate(config%inp%s8)
         call get_argument_as_real(iarg, config%inp%s8, error)
         if (allocated(error)) exit
      case("--s9")
         config%has_param = .true.
         iarg = iarg + 1
         allocate(config%inp%s9)
         call get_argument_as_real(iarg, config%inp%s9, error)
         if (allocated(error)) exit
      case("--a1")
         config%has_param = .true.
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%inp%a1, error)
         if (allocated(error)) exit
      case("--a2")
         config%has_param = .true.
         iarg = iarg + 1
         call get_argument_as_real(iarg, config%inp%a2, error)
         if (allocated(error)) exit
      case("--a3")
         config%has_param = .true.
         iarg = iarg + 1
         allocate(config%inp%a3)
         call get_argument_as_real(iarg, config%inp%a3, error)
         if (allocated(error)) exit
      case("--a4")
         config%has_param = .true.
         iarg = iarg + 1
         allocate(config%inp%a4)
         call get_argument_as_real(iarg, config%inp%a4, error)
         if (allocated(error)) exit
      case("--rs6")
         config%has_param = .true.
         iarg = iarg + 1
         allocate(config%inp%rs6)
         call get_argument_as_real(iarg, config%inp%rs6, error)
         if (allocated(error)) exit
      case("--rs8")
         config%has_param = .true.
         iarg = iarg + 1
         allocate(config%inp%rs8)
         call get_argument_as_real(iarg, config%inp%rs8, error)
         if (allocated(error)) exit
      case("--rs9")
         config%has_param = .true.
         iarg = iarg + 1
         allocate(config%inp%rs9)
         call get_argument_as_real(iarg, config%inp%rs9, error)
         if (allocated(error)) exit
      case("--alp")
         config%has_param = .true.
         iarg = iarg + 1
         allocate(config%inp%alp)
         call get_argument_as_real(iarg, config%inp%alp, error)
         if (allocated(error)) exit
      case("--bet")
         config%has_param = .true.
         iarg = iarg + 1
         allocate(config%inp%bet)
         call get_argument_as_real(iarg, config%inp%bet, error)
         if (allocated(error)) exit
      end select
   end do
   if (allocated(error)) return

   if (.not.config%has_param .and. .not.allocated(config%method)) then
      config%properties = .true.
   end if

   if (.not.config%atm .and. allocated(config%inp%s9)) then
      if (abs(config%inp%s9) > epsilon(1.0_wp)) then
         call fatal_error(error, "s9 parameter provided without ATM three-body dispersion")
         return
      end if
   end if

   if (.not.allocated(config%grad_output)) then
      config%grad_output = "dftd4.txt"
   end if

   if (.not.allocated(config%model)) then
      config%model = "d4"
   end if

   if (.not.allocated(config%input)) then
      if (.not.allocated(error)) then
         write(output_unit, '(a)') help_text
         error stop
      end if
   end if

end subroutine get_run_arguments


subroutine get_param_arguments(config, list, start, error)

   !> Configuation data
   type(param_config), intent(out) :: config

   !> List of command line arguments
   type(argument_list), intent(in) :: list

   !> First command line argument
   integer, intent(in) :: start

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   character(len=:), allocatable :: arg

   iarg = start
   narg = len(list)
   do while(iarg < narg)
      iarg = iarg + 1
      call list%get(iarg, arg)
      select case(arg)
      case("--help")
         write(output_unit, '(a)') help_text_param
         stop
      case("-l", "--list", "--funcs")
         config%list = .true.
      case default
         write(output_unit, '(a)') help_text_param
         stop
      end select
   end do
   if (allocated(error)) return

   ! check if anything was set
   if (.not.config%list) then
      write(output_unit, '(a)') help_text
      call fatal_error(error, "No arguments provided")
   end if

end subroutine get_param_arguments


subroutine get_argument_as_real(iarg, val, error)

   !> Index of command line argument, range [0:command_argument_count()]
   integer, intent(in) :: iarg

   !> Real value
   real(wp), intent(out) :: val

   !> Error handling
   type(error_type), allocatable :: error

   integer :: stat
   character(len=:), allocatable :: arg

   call get_argument(iarg, arg)
   if (.not.allocated(arg)) then
      call fatal_error(error, "Cannot read real value, argument missing")
      return
   end if
   read(arg, *, iostat=stat) val
   if (stat /= 0) then
      call fatal_error(error, "Cannot read real value from '"//arg//"'")
      return
   end if

end subroutine get_argument_as_real


end module dftd4_cli
