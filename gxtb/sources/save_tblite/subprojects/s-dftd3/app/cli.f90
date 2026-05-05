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

module dftd3_app_cli
   use, intrinsic :: iso_fortran_env, only : output_unit
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type, read_structure, filetype, get_filetype
   use dftd3, only : d3_param
   use dftd3_app_argument, only : argument_list, len
   use dftd3_app_help, only : prog_name, header, help_text, run_help_text, param_help_text, &
      & gcp_help_text, version
   implicit none
   private

   public :: app_config, run_config, param_config, gcp_config, get_arguments

   type, abstract :: app_config
   end type app_config

   type, extends(app_config) :: run_config
      !> Geometry input file
      character(len=:), allocatable :: input
      !> Format of the geometry input
      integer, allocatable :: input_format
      !> Method name
      character(len=:), allocatable :: method
      !> Basis name
      character(len=:), allocatable :: basis
      !> Damping paramaters
      type(d3_param) :: inp
      logical :: json = .false.
      character(len=:), allocatable :: json_output
      logical :: wrap = .true.
      logical :: tmer = .true.
      logical :: properties = .false.
      logical :: atm = .false.
      logical :: grad = .false.
      character(len=:), allocatable :: grad_output
      logical :: zero = .false.
      logical :: rational = .false.
      logical :: mzero = .false.
      logical :: mrational = .false.
      logical :: optimizedpower = .false.
      logical :: cso = .false.
      logical :: has_param = .false.
      integer :: verbosity = 2
      logical :: pair_resolved = .false.
      logical :: citation = .false.
      logical :: gcp = .false.
      character(len=:), allocatable :: citation_output
      !> Parameter data base
      character(len=:), allocatable :: db
   end type run_config

   type, extends(app_config) :: param_config
      !> Data base input file
      character(len=:), allocatable :: input
      !> Method name
      character(len=:), allocatable :: method
      !> Damping function
      character(len=:), allocatable :: damping
   end type param_config

   type, extends(app_config) :: gcp_config
      !> Geometry input file
      character(len=:), allocatable :: input
      !> Format of the geometry input
      integer, allocatable :: input_format
      !> Method name
      character(len=:), allocatable :: method
      !> Basis name
      character(len=:), allocatable :: basis
      logical :: json = .false.
      character(len=:), allocatable :: json_output
      logical :: wrap = .true.
      logical :: tmer = .true.
      logical :: grad = .false.
      character(len=:), allocatable :: grad_output
      integer :: verbosity = 2
      logical :: citation = .false.
      character(len=:), allocatable :: citation_output
   end type gcp_config

contains


subroutine get_argument_as_real(arg, val, error)
   !> Index of command line argument, range [0:command_argument_count()]
   character(len=:), intent(in), allocatable :: arg
   !> Real value
   real(wp), intent(out) :: val
   !> Error handling
   type(error_type), allocatable :: error

   integer :: stat

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


subroutine get_arguments(config, error)

   !> Configuation data
   class(app_config), allocatable, intent(out) :: config

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
      case("--help")
         call info_message(error, help_text)
         exit
      case("--version")
         call version(output_unit)
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
      case("gcp")
         allocate(gcp_config :: config)
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
   type is(gcp_config)
      call get_gcp_arguments(config, list, iarg, error)
   end select
end subroutine get_arguments


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
   logical :: read_args
   character(len=:), allocatable :: arg

   read_args = .true.
   iarg = start
   narg = len(list)
   do while(iarg < narg)
      iarg = iarg + 1
      call list%get(iarg, arg)
      if (.not.read_args) then
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      end if
      select case(arg)
      case("--")
         read_args = .false.
      case("--help")
         call info_message(error, run_help_text)
         exit
      case("--version")
         call version(output_unit)
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
            call fatal_error(error, "Unknown argument encountered: '"//arg//"'")
         else
            call fatal_error(error, "Too many positional arguments present")
         end if
         exit
      case("-i", "--input")
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         config%input_format = get_filetype("."//arg)
      case("--json")
         config%json = .true.
         config%json_output = "dftd3.json"
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%json_output)
         end if
      case("--citation")
         config%citation = .true.
         config%citation_output = "dftd3.bib"
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%citation_output)
         end if
      case("--property")
         config%properties = .true.
      case("--pair-resolved")
         config%pair_resolved = .true.
      case("--noedisp")
         config%tmer = .false.
      case("--nowrap")
         config%wrap = .false.
      case("--grad")
         config%grad = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%grad_output)
         end if
      case("--atm")
         config%inp%s9 = 1.0_wp
         config%atm = .true.
      case("--atm-scale")
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%s9, error)
         if (allocated(error)) exit
         config%atm = .true.
      case("--gcp")
         config%gcp = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%basis)
         end if
      case("--zero")
         config%zero = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
         call move_alloc(arg, config%method)
      case("--zerom")
         config%mzero = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
         call move_alloc(arg, config%method)
      case("--zero-param")
         config%zero = .true.
         config%has_param = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%s6, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%s8, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%rs6, error)
         if (allocated(error)) exit
      case("--zerom-param")
         config%mzero = .true.
         config%has_param = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%s6, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%s8, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%rs6, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%bet, error)
         if (allocated(error)) exit
      case("--bj")
         config%rational = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
         call move_alloc(arg, config%method)
      case("--bjm")
         config%mrational = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
         call move_alloc(arg, config%method)
      case("--bj-param", "--bjm-param")
         config%rational = .true.
         config%has_param = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%s6, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%s8, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%a1, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%a2, error)
         if (allocated(error)) exit
      case("--op")
         config%optimizedpower = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
         call move_alloc(arg, config%method)
      case("--op-param")
         config%optimizedpower = .true.
         config%has_param = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%s6, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%s8, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%a1, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%a2, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%bet, error)
         if (allocated(error)) exit
      case("--cso")
         config%cso = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
         call move_alloc(arg, config%method)
      case("--cso-param")
         config%cso = .true.
         config%has_param = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%s6, error)
         if (allocated(error)) exit
         iarg = iarg + 1
         call list%get(iarg, arg)
         call get_argument_as_real(arg, config%inp%a1, error)
         if (allocated(error)) exit
      case("--db")
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%db)
         end if
         if (.not.allocated(config%db)) then
            call fatal_error(error, "No argument provided for data base path")
            exit
         end if
      end select
   end do
   if (allocated(error)) return

   if (.not.config%has_param .and. .not.allocated(config%method)) then
      config%properties = .true.
   end if

   if (count([config%zero, config%rational, config%mzero, config%mrational, &
      & config%optimizedpower, config%cso]) > 1) then
      call fatal_error(error, "Can only select one damping scheme (zero, rational, mzero, mrational, optimizedpower, or cso)")
      return
   end if

   if (config%grad.and. .not.config%json .and. .not.allocated(config%grad_output)) then
      config%grad_output = "dftd3.txt"
   end if

   if (.not.allocated(config%input)) then
      if (.not.allocated(error)) then
         write(output_unit, '(a)') run_help_text
         call fatal_error(error, "Insufficient arguments provided")
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
         call info_message(error, param_help_text)
         exit
      case("--version")
         call version(output_unit)
         stop
      case default
         if (.not.allocated(config%input)) then
            call move_alloc(arg, config%input)
            cycle
         end if
         if (.not.allocated(config%method)) then
            call move_alloc(arg, config%method)
            cycle
         end if
         if (.not.allocated(config%damping)) then
            call move_alloc(arg, config%damping)
            cycle
         end if
         if (arg(1:1) == "-") then
            call fatal_error(error, "Unknown argument encountered: '"//arg//"'")
         else
            call fatal_error(error, "Too many positional arguments present")
         end if
         exit
      end select
   end do
   if (allocated(error)) return

   if (.not.allocated(config%input)) then
      if (.not.allocated(error)) then
         write(output_unit, '(a)') param_help_text
         call fatal_error(error, "Insufficient arguments provided")
      end if
   end if

end subroutine get_param_arguments

subroutine get_gcp_arguments(config, list, start, error)

   !> Configuation data
   type(gcp_config), intent(out) :: config

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
         call info_message(error, gcp_help_text)
         exit
      case("--version")
         call version(output_unit)
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
            call fatal_error(error, "Unknown argument encountered: '"//arg//"'")
         else
            call fatal_error(error, "Too many positional arguments present")
         end if
         exit
      case("-i", "--input")
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         config%input_format = get_filetype("."//arg)
      case("-l", "--level")
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         if (index(arg, "/") > 0) then
            config%method = arg(1:index(arg, "/")-1)
            config%basis = arg(index(arg, "/")+1:)
         else
            call move_alloc(arg, config%method)
         end if
      case("--json")
         config%json = .true.
         config%json_output = "dftd3.json"
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%json_output)
         end if
      case("--nocpc")
         config%tmer = .false.
      case("--nowrap")
         config%wrap = .false.
      case("--grad")
         config%grad = .true.
         iarg = iarg + 1
         call list%get(iarg, arg)
         if (allocated(arg)) then
            if (arg(1:1) == "-") then
               iarg = iarg - 1
               cycle
            end if
            call move_alloc(arg, config%grad_output)
         end if
      end select
   end do
   if (allocated(error)) return

   if (config%grad.and. .not.config%json .and. .not.allocated(config%grad_output)) then
      config%grad_output = "dftd3.txt"
   end if

   if (.not.allocated(config%input)) then
      if (.not.allocated(error)) then
         write(output_unit, '(a)') gcp_help_text
         call fatal_error(error, "Insufficient arguments provided")
      end if
   end if

end subroutine get_gcp_arguments

subroutine info_message(error, message)
   type(error_type), allocatable, intent(out) :: error
   character(len=*), intent(in) :: message

   allocate(error)
   error%stat = 0
   error%message = message
end subroutine info_message

end module dftd3_app_cli
