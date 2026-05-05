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

module dftd3_app_driver
   use, intrinsic :: iso_fortran_env, only : output_unit, input_unit
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type, read_structure, filetype, get_filetype
   use dftd3, only : damping_param, d3_param, d3_model, &
      & get_dispersion, get_zero_damping, zero_damping_param, new_zero_damping, &
      & get_rational_damping, rational_damping_param, new_rational_damping, &
      & get_mzero_damping, mzero_damping_param, new_mzero_damping, get_mrational_damping, &
      & get_optimizedpower_damping, optimizedpower_damping_param, &
      & new_optimizedpower_damping, &
      & get_cso_damping, cso_damping_param, new_cso_damping, &
      & new_d3_model, get_pairwise_dispersion, &
      & realspace_cutoff, get_lattice_points, get_coordination_number
   use dftd3_gcp, only : gcp_param, get_gcp_param, get_geometric_counterpoise
   use dftd3_output, only : ascii_damping_param, ascii_atomic_radii, &
      & ascii_atomic_references, ascii_system_properties, ascii_energy_atom, &
      & ascii_results, ascii_pairwise, tagged_result, json_results, &
      & turbomole_gradient, turbomole_gradlatt, ascii_gcp_param
   use dftd3_utils, only : wrap_to_central_cell
   use dftd3_citation, only : format_bibtex, is_citation_present, citation_type, &
      & get_citation, doi_dftd3_0, doi_dftd3_bj, doi_dftd3_m, doi_dftd3_op, &
      & doi_dftd3_cso, doi_joss, same_citation
   use dftd3_app_help, only : header
   use dftd3_app_cli, only : app_config, run_config, param_config, gcp_config, get_arguments
   use dftd3_app_toml, only : param_database
   implicit none
   private

   public :: app_driver

contains

subroutine app_driver(config, error)
   class(app_config), intent(in) :: config
   type(error_type), allocatable, intent(out) :: error

   select type(config)
   type is(run_config)
      call run_driver(config, error)
   type is(param_config)
      call param_driver(config, error)
   type is(gcp_config)
      call gcp_driver(config, error)
   end select
end subroutine app_driver

subroutine run_driver(config, error)
   type(run_config), intent(in) :: config
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(damping_param), allocatable :: param
   type(d3_param) :: inp
   type(d3_model) :: d3
   type(gcp_param) :: gcp
   real(wp), allocatable :: energies(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: pair_disp2(:, :), pair_disp3(:, :)
   real(wp), allocatable :: s9
   real(wp) :: energy
   character(len=:), allocatable :: output
   integer :: unit
   type(citation_type) :: citation, param_citation

   if (config%verbosity > 1) then
      call header(output_unit)
   end if

   if (config%input == "-") then
      if (.not.allocated(config%input_format)) then
         call read_structure(mol, input_unit, filetype%xyz, error)
      else
         call read_structure(mol, input_unit, config%input_format, error)
      end if
   else
      call read_structure(mol, config%input, error, config%input_format)
   end if
   if (allocated(error)) return
   if (config%wrap) then
      call wrap_to_central_cell(mol%xyz, mol%lattice, mol%periodic)
   end if

   if (config%has_param) inp = config%inp
   if (config%atm) s9 = config%inp%s9
   if (config%zero) then
      citation = get_citation(doi_dftd3_0)
      if (.not.config%has_param) then
         if (allocated(config%db)) then
            call from_db(param, config%db, config%method, "zero", error)
         else
            call get_zero_damping(inp, config%method, error, s9, param_citation)
         end if
         if (allocated(error)) return
      end if
      if (.not.allocated(param)) then
         block
            type(zero_damping_param), allocatable :: zparam
            allocate(zparam)
            call new_zero_damping(zparam, inp)
            call move_alloc(zparam, param)
         end block
      end if
   end if
   if (config%mzero) then
      citation = get_citation(doi_dftd3_m)
      if (.not.config%has_param) then
         if (allocated(config%db)) then
            call from_db(param, config%db, config%method, "zerom", error)
         else
            call get_mzero_damping(inp, config%method, error, s9, param_citation)
         end if
         if (allocated(error)) return
      end if
      if (.not.allocated(param)) then
         block
            type(mzero_damping_param), allocatable :: mparam
            allocate(mparam)
            call new_mzero_damping(mparam, inp)
            call move_alloc(mparam, param)
         end block
      end if
   end if
   if (config%rational .or. config%mrational) then
      if (config%rational) then
         citation = get_citation(doi_dftd3_bj)
      else
         citation = get_citation(doi_dftd3_m)
      end if
      if (.not.config%has_param) then
         if (config%mrational) then
            if (allocated(config%db)) then
               call from_db(param, config%db, config%method, "bjm", error)
            else
               call get_mrational_damping(inp, config%method, error, s9, param_citation)
            end if
         else
            if (allocated(config%db)) then
               call from_db(param, config%db, config%method, "bj", error)
            else
               call get_rational_damping(inp, config%method, error, s9, param_citation)
            end if
         end if
         if (allocated(error)) return
      end if
      if (.not.allocated(param)) then
         block
            type(rational_damping_param), allocatable :: rparam
            allocate(rparam)
            call new_rational_damping(rparam, inp)
            call move_alloc(rparam, param)
         end block
      end if
   end if
   if (config%optimizedpower) then
      citation = get_citation(doi_dftd3_op)
      if (.not.config%has_param) then
         if (allocated(config%db)) then
            call from_db(param, config%db, config%method, "op", error)
         else
            call get_optimizedpower_damping(inp, config%method, error, s9, param_citation)
         end if
         if (allocated(error)) return
      end if
      if (.not.allocated(param)) then
         block
            type(optimizedpower_damping_param), allocatable :: oparam
            allocate(oparam)
            call new_optimizedpower_damping(oparam, inp)
            call move_alloc(oparam, param)
         end block
      end if
   end if
   if (config%cso) then
      citation = get_citation(doi_dftd3_cso)
      if (.not.config%has_param) then
         if (allocated(config%db)) then
            call from_db(param, config%db, config%method, "cso", error)
         else
            call get_cso_damping(inp, config%method, error, s9, param_citation)
         end if
         if (allocated(error)) return
      else
         inp%a2 = 2.5_wp
         inp%rs6 = 0.0_wp
         inp%rs8 = 6.25_wp
      end if
      if (.not.allocated(param)) then
         block
            type(cso_damping_param), allocatable :: cparam
            allocate(cparam)
            call new_cso_damping(cparam, inp)
            call move_alloc(cparam, param)
         end block
      end if
   end if
   if (config%gcp) then
      call get_gcp_param(gcp, mol, config%method, config%basis)
   end if

   if (config%verbosity > 0) then
      if (allocated(param)) &
         call ascii_damping_param(output_unit, param, config%method)
      if (config%gcp) &
         call ascii_gcp_param(output_unit, mol, gcp)
   end if

   if (allocated(param)) then
      allocate(energies(mol%nat))
      if (config%grad) then
         allocate(gradient(3, mol%nat), sigma(3, 3))
      end if
   end if

   call new_d3_model(d3, mol)

   if (config%properties) then
      call property_calc(output_unit, mol, d3, config%verbosity)
   end if

   if (allocated(param)) then
      call get_dispersion(mol, d3, param, realspace_cutoff(), energies, &
         & gradient, sigma)
      if (config%gcp) then
         call get_geometric_counterpoise(mol, gcp, realspace_cutoff(), energies, gradient, sigma)
      end if
      energy = sum(energies)

      if (config%pair_resolved) then
         allocate(pair_disp2(mol%nat, mol%nat), pair_disp3(mol%nat, mol%nat))
         call get_pairwise_dispersion(mol, d3, param, realspace_cutoff(), pair_disp2, &
            & pair_disp3)
      end if
      if (config%verbosity > 0) then
         if (config%verbosity > 2) then
            call ascii_energy_atom(output_unit, mol, energies)
         end if
         call ascii_results(output_unit, mol, energy, gradient, sigma)
         if (config%pair_resolved) then
            call ascii_pairwise(output_unit, mol, pair_disp2, pair_disp3)
         end if
      end if
      if (config%tmer) then
         call tmer_writer(".EDISP", energy, "Dispersion", config%verbosity)
      end if
      if (config%grad) then
         if (allocated(config%grad_output)) then
            call results_writer(config%grad_output, energy, gradient, sigma, &
               & "Dispersion", config%verbosity)
         end if

         call turbomole_writer(mol, energy, gradient, sigma, config%verbosity, "Dispersion")
      end if

      if (config%json) then
         open(file=config%json_output, newunit=unit)
         call json_results(unit, "  ", energy=energy, gradient=gradient, sigma=sigma, &
            & pairwise_energy2=pair_disp2, pairwise_energy3=pair_disp3, param=param)
         close(unit)
         if (config%verbosity > 0) then
            write(output_unit, '(a)') &
               & "[Info] JSON dump of results written to '"//config%json_output//"'"
         end if
      end if

   end if

   if (config%citation) then
      open(file=config%citation_output, newunit=unit)
      call format_bibtex(output, get_citation(doi_joss))
      if (allocated(output)) write(unit, '(a)') output
      if (.not.same_citation(citation, param_citation)) then
         call format_bibtex(output, citation)
         if (allocated(output)) write(unit, '(a)') output
      end if
      call format_bibtex(output, param_citation)
      if (allocated(output)) write(unit, '(a)') output
      close(unit)
      if (config%verbosity > 0) then
         write(output_unit, '(a)') &
            & "[Info] Citation information written to '"//config%citation_output//"'"
      end if
   end if

end subroutine run_driver

subroutine property_calc(unit, mol, disp, verbosity)

   !> Unit for output
   integer, intent(in) :: unit

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Dispersion model
   class(d3_model), intent(in) :: disp

   !> Printout verbosity
   integer, intent(in) :: verbosity

   integer :: mref
   real(wp), allocatable :: cn(:), gwvec(:, :), c6(:, :), lattr(:, :)

   if (verbosity > 1) then
      call ascii_atomic_radii(unit, mol, disp)
      write(unit, '(a)')
      call ascii_atomic_references(unit, mol, disp)
      write(unit, '(a)')
   end if

   mref = maxval(disp%ref)
   allocate(cn(mol%nat), gwvec(mref, mol%nat), c6(mol%nat, mol%nat))
   call get_lattice_points(mol%periodic, mol%lattice, 30.0_wp, lattr)
   call get_coordination_number(mol, lattr, 30.0_wp, disp%rcov, cn)
   call disp%weight_references(mol, cn, gwvec)
   call disp%get_atomic_c6(mol, gwvec, c6=c6)

   if (verbosity > 0) then
      call ascii_system_properties(unit, mol, disp, cn, c6)
      write(unit, '(a)')
   end if

end subroutine property_calc

subroutine param_driver(config, error)
   type(param_config), intent(in) :: config
   type(error_type), allocatable, intent(out) :: error

   type(param_database) :: db
   class(damping_param), allocatable :: param

   call db%load(config%input, error)
   if (allocated(error)) return

   if (allocated(config%method)) then
      call db%get(param, config%method, config%damping)
      if (.not.allocated(param)) then
         call fatal_error(error, "No entry for '"//config%method//"' found in '"//config%input//"'")
         return
      end if
      call ascii_damping_param(output_unit, param, config%method)
   else
      write(output_unit, '(a, *(1x, g0))') "[Info] Found", size(db%records), &
         "damping parameters in '"//config%input//"'"
   end if

end subroutine param_driver

subroutine from_db(param, input, method, damping, error)
   class(damping_param), allocatable, intent(out) :: param
   character(len=*), intent(in) :: input
   character(len=*), intent(in) :: method
   character(len=*), intent(in) :: damping
   type(error_type), allocatable, intent(out) :: error

   type(param_database) :: db

   call db%load(input, error)
   if (allocated(error)) return

   call db%get(param, method, damping)
   if (.not.allocated(param)) then
      call fatal_error(error, "No entry for '"//method//"' found in '"//input//"'")
      return
   end if
end subroutine from_db

subroutine gcp_driver(config, error)
   type(gcp_config), intent(in) :: config
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(gcp_param) :: param
   integer :: unit
   real(wp) :: energy
   real(wp), allocatable :: energies(:), gradient(:, :), sigma(:, :)

   if (config%verbosity > 1) then
      call header(output_unit)
   end if

   if (config%input == "-") then
      if (.not.allocated(config%input_format)) then
         call read_structure(mol, input_unit, filetype%xyz, error)
      else
         call read_structure(mol, input_unit, config%input_format, error)
      end if
   else
      call read_structure(mol, config%input, error, config%input_format)
   end if
   if (allocated(error)) return
   if (config%wrap) then
      call wrap_to_central_cell(mol%xyz, mol%lattice, mol%periodic)
   end if

   call get_gcp_param(param, mol, config%method, config%basis)
   if (config%verbosity > 0) then
      call ascii_gcp_param(output_unit, mol, param)
   end if

   allocate(energies(mol%nat), source=0.0_wp)
   if (config%grad) then
      allocate(gradient(3, mol%nat), source=0.0_wp)
      allocate(sigma(3, 3), source=0.0_wp)
   end if

   call get_geometric_counterpoise(mol, param, realspace_cutoff(), energies, gradient, sigma)
   energy = sum(energies)

   if (config%verbosity > 0) then
      if (config%verbosity > 2) then
         call ascii_energy_atom(output_unit, mol, energies, label="counter-poise")
      end if
      call ascii_results(output_unit, mol, energy, gradient, sigma, label="Counter-poise")
   end if

   if (config%tmer) then
      call tmer_writer(".CPC", energy, "Counter-poise", config%verbosity)
   end if

   if (config%grad) then
      if (allocated(config%grad_output)) then
         call results_writer(config%grad_output, energy, gradient, sigma, &
            & "Counter-poise", config%verbosity)
      end if

      call turbomole_writer(mol, energy, gradient, sigma, config%verbosity, "Counter-poise")
   end if

   if (config%json) then
      open(file=config%json_output, newunit=unit)
      call json_results(unit, "  ", energy=energy, gradient=gradient, sigma=sigma)
      close(unit)
      if (config%verbosity > 0) then
         write(output_unit, '(a)') &
            & "[Info] JSON dump of results written to '"//config%json_output//"'"
      end if
   end if
end subroutine gcp_driver

!> Write the energy to a file
subroutine tmer_writer(filename, energy, label, verbosity)
   !> File name to write to
   character(len=*), intent(in) :: filename

   !> Energy to write
   real(wp), intent(in) :: energy

   !> Label for the output
   character(len=*), intent(in) :: label

   !> Printout verbosity
   integer, intent(in) :: verbosity

   integer :: unit

   if (verbosity > 0) then
      if (verbosity > 1) then
         write(output_unit, '(a)') "[Info] Writing "//label//" energy to '"//filename//"'"
      end if
      open(file=filename, newunit=unit)
      write(unit, '(f24.14)') energy
      close(unit)
   end if
end subroutine tmer_writer

!> Write the results to a file
subroutine results_writer(filename, energy, gradient, sigma, label, verbosity)
   !> File name to write to
   character(len=*), intent(in) :: filename

   !> Energy to write
   real(wp), intent(in) :: energy

   !> Gradient to write
   real(wp), intent(in) :: gradient(:, :)

   !> Virial tensor to write
   real(wp), intent(in) :: sigma(:, :)

   !> Label for the output
   character(len=*), intent(in) :: label

   !> Printout verbosity
   integer, intent(in) :: verbosity

   integer :: unit

   open(file=filename, newunit=unit)
   call tagged_result(unit, energy, gradient, sigma)
   close(unit)
   if (verbosity > 0) then
      write(output_unit, '(a)') "[Info] "//label//" results written to '"//filename//"'"
   end if
end subroutine results_writer

!> Write the gradient and virial tensor to Turbomole files
subroutine turbomole_writer(mol, energy, gradient, sigma, verbosity, label)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Energy of the system
   real(wp), intent(in) :: energy

   !> Gradient of the system
   real(wp), intent(in) :: gradient(:, :)

   !> Virial tensor of the system
   real(wp), intent(in) :: sigma(:, :)

   !> Printout verbosity
   integer, intent(in) :: verbosity

   !> Label for the output
   character(len=*), intent(in) :: label

   logical :: exist
   integer :: stat

   inquire(file="gradient", exist=exist)
   if (exist) then
      call turbomole_gradient(mol, "gradient", energy, gradient, stat)
      if (verbosity > 0) then
         if (stat == 0) then
            write(output_unit, '(a)') &
               & "[Info] "//label//" gradient added to Turbomole gradient file"
         else
            write(output_unit, '(a)') &
               & "[Warn] Could not add to Turbomole gradient file"
         end if
      end if
   end if
   inquire(file="gradlatt", exist=exist)
   if (exist) then
      call turbomole_gradlatt(mol, "gradlatt", energy, sigma, stat)
      if (verbosity > 0) then
         if (stat == 0) then
            write(output_unit, '(a)') &
               & "[Info] "//label//" virial added to Turbomole gradlatt file"
         else
            write(output_unit, '(a)') &
               & "[Warn] Could not add to Turbomole gradlatt file"
         end if
      end if
   end if
end subroutine turbomole_writer

end module dftd3_app_driver
