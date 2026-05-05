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

!> Implementation of the driver entry points for the singlepoint runner
module tblite_driver_guess
   use, intrinsic :: iso_fortran_env, only : error_unit, input_unit
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type, read_structure, filetype
   use mctc_io_constants, only : codata
   use mctc_io_convert, only : aatoau, ctoau
   use tblite_cli, only : guess_config
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_type
   use tblite_context, only : context_type, context_terminal, escape
   use tblite_external_field, only : electric_field
   use tblite_lapack_solver, only : lapack_solver
   use tblite_output_ascii
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, &
      & sad_guess, eeq_guess, eeqbc_guess, get_molecular_dipole_moment
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_ceh_singlepoint, only : ceh_singlepoint
   use tblite_ceh_ceh, only : new_ceh_calculator

   implicit none
   private

   public :: main

   interface main
      module procedure :: guess_main
   end interface

   real(wp), parameter :: kt = 3.166808578545117e-06_wp
   real(wp), parameter :: jtoau = 1.0_wp / (codata%me*codata%c**2*codata%alpha**2)
   !> Convert V/Å = J/(C·Å) to atomic units
   real(wp), parameter :: vatoau = jtoau / (ctoau * aatoau)

contains


   subroutine guess_main(config, error)
      type(guess_config), intent(in) :: config
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      character(len=:), allocatable :: method, filename
      integer :: unpaired, charge, unit
      real(wp) :: dpmom(3), etemp
      type(context_type) :: ctx
      type(xtb_calculator):: calc_ceh
      type(wavefunction_type) :: wfn_ceh
      real(wp), allocatable :: qat(:)
      real(wp), allocatable :: dpat(:,:)

      ctx%terminal = context_terminal(config%color)
      ctx%solver = lapack_solver(config%solver)

      if (config%input == "-") then
         if (allocated(config%input_format)) then
            call read_structure(mol, input_unit, config%input_format, error)
         else
            call read_structure(mol, input_unit, filetype%xyz, error)
         end if
      else
         call read_structure(mol, config%input, error, config%input_format)
      end if
      if (allocated(error)) return

      if (config%grad) then
         call fatal_error(error, "Charge gradient not yet implemented.")
      end if
      if (allocated(error)) return

      if (allocated(config%charge)) then
         mol%charge = config%charge
      else
         filename = join(dirname(config%input), ".CHRG")
         if (exists(filename)) then
            call read_file(filename, charge, error)
            if (allocated(error)) return
            if (config%verbosity > 0) &
            & call info(ctx, "Molecular charge read from '"//filename//"'")
            mol%charge = charge
         end if
      end if
      if (allocated(error)) return

      if (allocated(config%spin)) then
         mol%uhf = config%spin
      else
         filename = join(dirname(config%input), ".UHF")
         if (exists(filename)) then
            call read_file(filename, unpaired, error)
            if (allocated(error)) return
            if (config%verbosity > 0) &
            & call info(ctx, "Molecular spin read from '"//filename//"'")
            mol%uhf = unpaired
         end if
      end if
      if (allocated(error)) return

      method = "ceh"
      if (allocated(config%method)) method = config%method
      if (method == "ceh") then
         call new_ceh_calculator(calc_ceh, mol, error)
         if (allocated(error)) return
         if (allocated(config%etemp_guess)) then
            etemp = config%etemp_guess
         else
            etemp = calc_ceh%default_etemp
         end if
         call new_wavefunction(wfn_ceh, mol%nat, calc_ceh%bas%nsh, calc_ceh%bas%nao, 1, etemp * kt)
      end if

      if (allocated(config%efield) .and. config%method == "ceh") then
         block
            class(container_type), allocatable :: cont
            cont = electric_field(config%efield*vatoau)
            call calc_ceh%push_back(cont)
         end block
      end if

      if (config%verbosity > 0) then
         select case(method)
         case default
            call fatal_error(error, "Unknown method '"//method//"' requested")
         case("sad")
            call ctx%message("Superposition of atomic densities (SAD) guess")
            call ctx%message("")
         case("eeq")
            call ctx%message("Electronegativity equilibration (EEQ) guess")
            call ctx%message("")
         case("eeqbc")
            call ctx%message("Bond Capacity electronegativity equilibration (EEQ) guess")
            call ctx%message("")
         case("ceh")
            call ctx%message(calc_ceh%info(config%verbosity, " | "))
         end select
      end if

      allocate(qat(mol%nat), dpat(3, mol%nat), source=0.0_wp)
      select case(method)
      case default
         call fatal_error(error, "Unknown method '"//method//"' requested")
      case("sad")
         call sad_guess(mol, qat, dpat)
      case("eeq")
         call eeq_guess(mol, qat, dpat, error)
      case("eeqbc")
         call eeqbc_guess(mol, qat, dpat, error)
      case("ceh")
         call ceh_singlepoint(ctx, calc_ceh, mol, wfn_ceh, config%accuracy, config%verbosity)
         if (ctx%failed()) then
            call fatal(ctx, "CEH singlepoint calculation failed")
            do while(ctx%failed())
               call ctx%get_error(error)
               write(error_unit, '("->", 1x, a)') error%message
            end do
            error stop
         end if
         qat = wfn_ceh%qat(:, 1)
         dpat = wfn_ceh%dpat(:, :, 1)
      end select
      if (allocated(error)) return

      call get_molecular_dipole_moment(mol, qat, dpat, dpmom)
      call ascii_atomic_charges(ctx%unit, 1, mol, qat)
      call ascii_dipole_moments(ctx%unit, 1, mol, dpat, dpmom)

      if (config%json) then
         open(file=config%json_output, newunit=unit)
         call json_results(unit, "  ", charges=qat)
         close(unit)
         if (config%verbosity > 0) then
            call info(ctx, "JSON dump of results written to '"//config%json_output//"'")
         end if
      end if
   end subroutine guess_main


   subroutine info(ctx, message)
      type(context_type), intent(inout) :: ctx
      character(len=*), intent(in) :: message

      call ctx%message( &
      & escape(ctx%terminal%bold) // "[" // &
      & escape(ctx%terminal%bold_cyan) // "Info" // &
      & escape(ctx%terminal%bold) // "]" // &
      & escape(ctx%terminal%reset) // " " // &
      & message)
   end subroutine info


   subroutine warn(ctx, message)
      type(context_type), intent(inout) :: ctx
      character(len=*), intent(in) :: message

      call ctx%message( &
      & escape(ctx%terminal%bold) // "[" // &
      & escape(ctx%terminal%bold_yellow) // "Warn" // &
      & escape(ctx%terminal%bold) // "]" // &
      & escape(ctx%terminal%reset) // " " // &
      & message)
   end subroutine warn


   subroutine fatal(ctx, message)
      type(context_type), intent(inout) :: ctx
      character(len=*), intent(in) :: message

      call ctx%message( &
      & escape(ctx%terminal%bold) // "[" // &
      & escape(ctx%terminal%bold_red) // "Fatal" // &
      & escape(ctx%terminal%bold) // "]" // &
      & escape(ctx%terminal%reset) // " " // &
      & message)
   end subroutine fatal

   !> Extract dirname from path
   function dirname(filename)
      character(len=*), intent(in) :: filename
      character(len=:), allocatable :: dirname

      dirname = filename(1:scan(filename, "/\", back=.true.))
      if (len_trim(dirname) == 0) dirname = "."
   end function dirname


   !> Construct path by joining strings with os file separator
   function join(a1, a2) result(path)
      use mctc_env_system, only : is_windows
      character(len=*), intent(in) :: a1, a2
      character(len=:), allocatable :: path
      character :: filesep

      if (is_windows()) then
         filesep = '\'
      else
         filesep = '/'
      end if

      path = a1 // filesep // a2
   end function join


   !> test if pathname already exists
   function exists(filename)
      character(len=*), intent(in) :: filename
      logical :: exists
      inquire(file=filename, exist=exists)
   end function exists


   subroutine read_file(filename, val, error)
      use mctc_io_utils, only : next_line, read_next_token, io_error, token_type
      character(len=*), intent(in) :: filename
      integer, intent(out) :: val
      type(error_type), allocatable, intent(out) :: error

      integer :: io, stat, lnum, pos
      type(token_type) :: token
      character(len=:), allocatable :: line

      lnum = 0

      open(file=filename, newunit=io, status='old', iostat=stat)
      if (stat /= 0) then
         call fatal_error(error, "Error: Could not open file '"//filename//"'")
         return
      end if

      call next_line(io, line, pos, lnum, stat)
      if (stat == 0) &
         call read_next_token(line, pos, token, val, stat)
      if (stat /= 0) then
         call io_error(error, "Cannot read value from file", line, token, &
            filename, lnum, "expected integer value")
         return
      end if

      close(io, iostat=stat)

   end subroutine read_file


end module tblite_driver_guess
