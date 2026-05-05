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
module tblite_driver_run
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type, read_structure, filetype
   use mctc_io_constants, only : codata
   use mctc_io_convert, only : aatoau, ctoau
   use multicharge, only : get_charges
   use tblite_cli, only : run_config
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_type
   use tblite_context, only : context_type, context_terminal, escape
   use tblite_data_spin, only : get_spin_constant
   use tblite_external_field, only : electric_field
   use tblite_lapack_solver, only : lapack_solver
   use tblite_output_ascii
   use tblite_param, only : param_record
   use tblite_results, only : results_type
   use tblite_spin, only : spin_polarization, new_spin_polarization
   use tblite_solvation, only : new_solvation, new_solvation_cds, &
      & new_solvation_shift, solvation_type
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, &
      & sad_guess, eeq_guess, eeqbc_guess, shell_partition, & 
      & load_wavefunction, save_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator, new_xtb_calculator
   use tblite_xtb_gxtb, only : new_gxtb_calculator, export_gxtb_param
   use tblite_xtb_gfn2, only : new_gfn2_calculator, export_gfn2_param
   use tblite_xtb_gfn1, only : new_gfn1_calculator, export_gfn1_param
   use tblite_xtb_h0, only : get_number_electrons
   use tblite_xtb_ipea1, only : new_ipea1_calculator, export_ipea1_param
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   use tblite_ceh_singlepoint, only : ceh_singlepoint
   use tblite_ceh_ceh, only : new_ceh_calculator
   use tblite_post_processing_list, only : add_post_processing, post_processing_type, &
      & post_processing_list

   implicit none
   private

   public :: main

   interface main
      module procedure :: run_main
   end interface

   real(wp), parameter :: kt = 3.166808578545117e-06_wp
   real(wp), parameter :: jtoau = 1.0_wp / (codata%me*codata%c**2*codata%alpha**2)
   !> Convert V/Å = J/(C·Å) to atomic units
   real(wp), parameter :: vatoau = jtoau / (ctoau * aatoau)
   character(len=:), allocatable :: wbo_label, molmom_label

contains


subroutine run_main(config, error)
   type(run_config), intent(in) :: config
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   character(len=:), allocatable :: method, filename, guess_
   integer :: unpaired, charge, unit, nspin
   logical :: restart_exist, use_guess, defined_spin
   real(wp) :: energy, etemp, nel
   real(wp), allocatable :: dpmom(:), qpmom(:)
   real(wp), allocatable :: gradient(:, :), sigma(:, :)
   real(wp), allocatable :: wll(:, :, :)
   type(param_record) :: param
   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(xtb_calculator) :: calc_ceh
   type(wavefunction_type) :: wfn
   type(wavefunction_type), allocatable :: wfn_ceh, wfn_aux
   type(results_type) :: results
   class(post_processing_list), allocatable :: post_proc

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

   defined_spin = .false.
   if (allocated(config%spin)) then
      mol%uhf = config%spin
      defined_spin = .true.
   else
      filename = join(dirname(config%input), ".UHF")
      if (exists(filename)) then
         call read_file(filename, unpaired, error)
         if (allocated(error)) return
         if (config%verbosity > 0) &
            & call info(ctx, "Molecular spin read from '"//filename//"'")
         mol%uhf = unpaired
      end if
      defined_spin = .true.
   end if

   nspin = merge(2, 1, config%spin_polarized)

   if (config%grad) then
      allocate(gradient(3, mol%nat), sigma(3, 3))
   end if
   
   if (allocated(error)) return

   if (allocated(config%param)) then
      call param%load(config%param, error)
      if (.not. allocated(error)) then
         call new_xtb_calculator(calc, mol, param, error, config%accuracy)
      end if
   else
      method = "gfn2"
      if (allocated(config%method)) method = config%method
      select case(method)
      case default
         call fatal_error(error, "Unknown method '"//method//"' requested")
      case("gxtb")
         call new_gxtb_calculator(calc, mol, error, config%accuracy)
         ! For g-xTB, we always do a UHF calculation if unpaired electrons
         ! are defined or there is an uneven number of electrons
         if (defined_spin) then
            nspin = 2
         else
            call get_number_electrons(mol, calc%bas, calc%h0, nel)
            if (mod(nint(nel), 2) > 0) then
               nspin = 2
            end if
         end if
      case("gfn2")
         call new_gfn2_calculator(calc, mol, error, config%accuracy)
      case("gfn1")
         call new_gfn1_calculator(calc, mol, error, config%accuracy)
      case("ipea1")
         call new_ipea1_calculator(calc, mol, error, config%accuracy)
      end select
   end if
   if (allocated(error)) return

   if (allocated(config%max_iter)) calc%iterator%max_iter = config%max_iter
   if (allocated(config%broyden_start) .or. allocated(config%diis_start)) then
      call calc%iterator%set_mixer_start(broyden_start=config%broyden_start, &
         & diis_start=config%diis_start)
   end if

   ! Optionally overwrite the default electronic temperature
   if (allocated(config%etemp)) then
      etemp = config%etemp
   else 
      etemp = calc%default_etemp
   end if

   use_guess = .true.
   restart_exist = .false.
   if (config%restart) restart_exist = exists(config%restart_file)
   if (restart_exist) then
      call load_wavefunction(wfn, config%restart_file, error, &
         & nat=mol%nat, nsh=calc%bas%nsh, nao=calc%bas%nao, nspin=nspin)
      if (allocated(error)) then
         call warn(ctx, error%message)
         deallocate(error)
         call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, nspin, etemp * kt)
      else
         call info(ctx, "Restarting from wavefunction '"//config%restart_file//"'")
         use_guess = .false.
         ! Optionally overwrite the electronic temperature from restart file
         if (allocated(config%etemp)) then
            wfn%kt = config%etemp * kt
         end if
      end if      
   else
      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, nspin, etemp * kt)
   end if


   ! Create auxiliary wavefunction and atomic charges
   if (allocated(calc%charge_model)) then
      allocate(wfn_aux)
      call new_wavefunction(wfn_aux, mol%nat, calc%bas%nsh, 0, 1, 0.0_wp, config%grad)
      if (config%grad) then
         call get_charges(calc%charge_model, mol, error, wfn_aux%qat(:, 1), &
            & dqdr=wfn_aux%dqatdr(:, :, :, 1), dqdL=wfn_aux%dqatdL(:, :, :, 1))
      else
         call get_charges(calc%charge_model, mol, error, wfn_aux%qat(:, 1))
      end if
      if (allocated(error)) return
   end if

   if (allocated(config%efield)) then
      block
         class(container_type), allocatable :: cont
         cont = electric_field(config%efield*vatoau)
         call calc%push_back(cont)
      end block
   end if

   if (use_guess) then
      ! Default to calculator specific calculator specific guess if needed
      if (.not. allocated(config%guess)) then
         guess_ = calc%default_guess
      else
         guess_ = config%guess
      end if

      select case(guess_)
      case default
         call fatal_error(error, "Unknown starting guess '"//guess_//"' requested")
      case("sad")
         call sad_guess(mol, calc, wfn)
      case("eeq")
         if (.not. allocated(wfn_aux)) then
            call eeq_guess(mol, calc, wfn, error)
         else
            wfn%qat(:, 1) = wfn_aux%qat(:, 1)
            wfn%qsh(:, 1) = wfn_aux%qsh(:, 1)
         end if
      case("eeqbc")
         if (.not. allocated(wfn_aux)) then
            call eeqbc_guess(mol, calc, wfn, error)
         else
            wfn%qat(:, 1) = wfn_aux%qat(:, 1)
            wfn%qsh(:, 1) = wfn_aux%qsh(:, 1)
         end if
      case("ceh")
         ! Setup CEH calculator and wavefunction
         allocate(wfn_ceh)
         call new_ceh_calculator(calc_ceh, mol, error)
         if (allocated(error)) return
         if (allocated(config%etemp_guess)) then
            etemp = config%etemp_guess
         else
            etemp = calc_ceh%default_etemp
         end if
         call new_wavefunction(wfn_ceh, mol%nat, calc_ceh%bas%nsh, calc_ceh%bas%nao, 1,&
            & etemp * kt)
         ! Check if an electric field is present
         if (allocated(config%efield)) then
            block
            class(container_type), allocatable :: cont
            cont = electric_field(config%efield*vatoau)
            call calc_ceh%push_back(cont)
            end block
         end if
         ! Perform CEH calculation
         call ceh_singlepoint(ctx, calc_ceh, mol, wfn_ceh, config%accuracy, config%verbosity)
         if (ctx%failed()) then
            call fatal(ctx, "CEH singlepoint calculation failed")
            do while(ctx%failed())
               call ctx%get_error(error)
               write(error_unit, '("->", 1x, a)') error%message
            end do
            error stop
         end if
         wfn%qat(:, 1) = wfn_ceh%qat(:, 1)
         ! Remove CEH wavefunction as it is no longer needed
         ! if CEH is only used as a guess
         deallocate(wfn_ceh)
         call shell_partition(mol, calc, wfn)
      end select
      if (allocated(error)) return
   end if

   ! Optionally add spin-polarization to the calculator
   if (config%spin_polarized .and. .not. allocated(calc%spin_polarization)) then
      call get_spin_constants(mol, calc%bas, wll)
      allocate(calc%spin_polarization)
      call new_spin_polarization(calc%spin_polarization, mol, wll, calc%bas%nsh_id)
   end if

   if (allocated(config%solvation)) then
      method = "gfn2"
      if (allocated(config%method)) method = config%method
      block
         class(container_type), allocatable :: cont
         class(solvation_type), allocatable :: solv
         call new_solvation(solv, mol, config%solvation, error, method)
         if (allocated(error)) return
         call move_alloc(solv, cont)
         call calc%push_back(cont)
      end block
      if (allocated(config%solvation%cds)) then
         block
            class(container_type), allocatable :: cont
            class(solvation_type), allocatable :: cds
            call new_solvation_cds(cds, mol, config%solvation, error, method)
            if (allocated(error)) return
            call move_alloc(cds, cont)
            call calc%push_back(cont)
         end block
      end if
      if (allocated(config%solvation%shift)) then
         block
            class(container_type), allocatable :: cont
            class(solvation_type), allocatable :: shift
            call new_solvation_shift(shift, config%solvation, error, method)
            if (allocated(error)) return
            call move_alloc(shift, cont)
            call calc%push_back(cont)
         end block
      end if
   end if

   wbo_label = "bond-orders"
   allocate(post_proc)
   call add_post_processing(post_proc, wbo_label, error)
   if (allocated(error)) return

   if (config%verbosity > 2) then
      molmom_label = "molmom"
      call add_post_processing(post_proc, molmom_label, error)
      if (allocated(error)) return
   end if

   if (allocated(config%post_processing)) then
      call add_post_processing(post_proc, config%post_processing, error)
      if (allocated(error)) return
   end if

   if (allocated(param%post_proc)) then
      call add_post_processing(post_proc, param%post_proc)
      if (allocated(error)) return
   end if

   if (config%verbosity > 0) then
      call ctx%message(calc%info(config%verbosity, " | "))
      call ctx%message("")
   end if

   call xtb_singlepoint(ctx, mol, calc, wfn, config%accuracy, energy, gradient, sigma, &
      & config%verbosity, results, post_proc, wfn_aux)
   if (ctx%failed()) then
      call fatal(ctx, "Singlepoint calculation failed")
      do while(ctx%failed())
         call ctx%get_error(error)
         write(error_unit, '("->", 1x, a)') error%message
      end do
      error stop
   end if

   if (allocated(config%restart_file)) then
      call info(ctx, "Writing wavefunction information to '"//config%restart_file//"'")
      call save_wavefunction(wfn, config%restart_file, error)
      if (allocated(error)) return
   end if

   if (allocated(post_proc) .and. allocated(config%post_proc_output)) then
      call results%dict%dump(config%post_proc_output, error)
      if (allocated(error)) return
   end if

   if (config%verbosity > 2) then
      call ascii_levels(ctx%unit, config%verbosity, wfn%emo, wfn%focc, 7)
      call post_proc%dict%get_entry("molecular-dipole", dpmom)
      call post_proc%dict%get_entry("molecular-quadrupole", qpmom)
      
      call ascii_dipole_moments(ctx%unit, 1, mol, wfn%dpat(:, :, 1), dpmom)
      call ascii_quadrupole_moments(ctx%unit, 1, mol, wfn%qpat(:, :, 1), qpmom)
   end if

   deallocate(post_proc)

   if (allocated(config%grad_output)) then
      open(file=config%grad_output, newunit=unit)
      call tagged_result(unit, energy, gradient, sigma, energies=results%energies)
      close(unit)
      if (config%verbosity > 0) then
         call info(ctx, "Tight-binding results written to '"//config%grad_output//"'")
      end if
   end if

   if (config%json) then
      open(file=config%json_output, newunit=unit)
      call json_results(unit, "  ", energy=energy, gradient=gradient, sigma=sigma, &
         & energies=results%energies)
      close(unit)
      if (config%verbosity > 0) then
         call info(ctx, "JSON dump of results written to '"//config%json_output//"'")
      end if
   end if
end subroutine run_main


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


subroutine get_spin_constants(mol, bas, wll)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Spin-constants
   real(wp), allocatable, intent(out) :: wll(:, :, :)
   
   integer :: izp, ish, jsh, il, jl

   allocate(wll(bas%nsh, bas%nsh, mol%nid), source=0.0_wp)

   do izp = 1, mol%nid
      do ish = 1, bas%nsh_id(izp)
         il = bas%cgto(ish, izp)%raw%ang
         do jsh = 1, bas%nsh_id(izp)
            jl = bas%cgto(jsh, izp)%raw%ang
            wll(jsh, ish, izp) = get_spin_constant(jl, il, mol%num(izp))
         end do
      end do
   end do
end subroutine get_spin_constants


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


end module tblite_driver_run
