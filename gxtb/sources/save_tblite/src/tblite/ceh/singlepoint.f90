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

!> @file tblite/ceh/singlepoint.f90
!> Provides main entry point for performing a CEH calculation
!> follows in close analogy the xtb singlepoint

!> Implementation of the single point calculation for the CEH model
module tblite_ceh_singlepoint
   use mctc_env, only : error_type, wp
   use mctc_io, only: structure_type
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_cutoff, only : get_lattice_points
   use tblite_basis_cache, only : basis_cache
   use tblite_basis_type, only : get_cutoff
   use tblite_context, only : context_type
   use tblite_output_format, only: format_string
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_wavefunction, only : wavefunction_type, &
   & get_alpha_beta_occupation
   use tblite_wavefunction_mulliken, only: get_mulliken_shell_charges, &
   & get_mulliken_atomic_multipoles
   use tblite_scf_iterator, only: next_density, get_qat_from_qsh
   use tblite_scf, only: new_potential, potential_type 
   use tblite_container, only : container_cache
   use tblite_scf_potential, only: add_pot_to_h1
   use tblite_scf_solver, only : solver_type
   use tblite_blas, only : gemv
   use tblite_ceh_ceh, only : get_effective_qat
   use tblite_xtb_h0, only : get_hamiltonian, get_selfenergy, get_occupation
   use tblite_xtb_spec, only : tb_h0spec 
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_timer, only : timer_type, format_time
   implicit none
   private

   public :: ceh_singlepoint

   real(wp), parameter :: cn_cutoff = 25.0_wp


   character(len=*), parameter :: real_format = "(es20.13)"
   character(len=25), parameter :: &
      label_cutoff = "integral cutoff", &
      label_charges = "CEH atomic charges", &
      label_dipole = "CEH mol. dip. mom. / a.u."

contains


   !> Run the CEH calculation (equivalent to xtb_singlepoint)
   subroutine ceh_singlepoint(ctx, calc, mol, wfn, accuracy, verbosity)
      !> Calculation context
      type(context_type), intent(inout) :: ctx
      !> CEH calculator
      type(xtb_calculator), intent(inout) :: calc
      !> Molecular structure data
      type(structure_type), intent(in)  :: mol
      !> Wavefunction data
      type(wavefunction_type), intent(inout) :: wfn
      !> Accuracy for computation
      real(wp), intent(in) :: accuracy
      !> Verbosity level of output
      integer, intent(in), optional :: verbosity

      ! Molecular dipole moment
      real(wp) :: dipole(3)
      ! Integral container
      type(integral_type) :: ints
      ! Electronic solver
      class(solver_type), allocatable :: solver
      ! Adjacency list
      type(adjacency_list) :: list
      ! Potential type
      type(potential_type) :: pot
      ! Restart data for interaction containers and coulomb 
      type(container_cache) :: icache, ccache
      ! Timer
      type(timer_type) :: timer
      ! Error container
      type(error_type), allocatable :: error
      
      logical :: grad
      real(wp) :: elec_entropy
      real(wp) :: nel, cutoff
      real(wp), allocatable :: tmp(:)

      integer :: prlevel

      real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
      real(wp), allocatable :: cn_en(:), dcn_endr(:, :, :), dcn_endL(:, :, :)
      real(wp), allocatable :: aniso_dip(:, :)
      real(wp), allocatable :: selfenergy(:), dsedcn(:), dsedcn_en(:), lattr(:, :)
      type(basis_cache) :: bcache

      call timer%push("total CEH")

      if (present(verbosity)) then
         prlevel = verbosity
      else
         prlevel = ctx%verbosity
      end if

      if (prlevel > 1) then
         call ctx%message("CEH singlepoint")
      endif
      ! Gradient logical as future starting point (not implemented yet)
      ! Entry point could either be (i) modified wavefunction type (including derivatives),
      ! (iii) additional wavefunction derivative type (see old commits) or (ii) optional
      ! dqdR variable in this routine
      grad = .false.

      ! Define occupation
      call get_occupation(mol, calc%bas, calc%h0, wfn%nocc, wfn%n0at, wfn%n0sh)
      nel = sum(wfn%n0at) - mol%charge
      if (mod(mol%uhf, 2) == mod(nint(nel), 2)) then
         wfn%nuhf = mol%uhf
      else
         wfn%nuhf = mod(nint(nel), 2)
      end if
      call get_alpha_beta_occupation(wfn%nocc, wfn%nuhf, wfn%nel(1), wfn%nel(2))

      call timer%push("hamiltonian")
      ! calculate coordination number (CN) and the EN-weighted coordination number
      if (allocated(calc%ncoord)) then
         allocate(cn(mol%nat))
         if (grad) then
            allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
         end if
         call calc%ncoord%get_cn(mol, cn, dcndr, dcndL)
      end if
      if (allocated(calc%ncoord_en)) then
         allocate(cn_en(mol%nat))
         if (grad) then
            allocate(dcn_endr(3, mol%nat, mol%nat), dcn_endL(3, 3, mol%nat))
         end if
         call calc%ncoord_en%get_cn(mol, cn_en, dcn_endr, dcn_endL)
      end if

      ! calculate the scaled self energies
      allocate(selfenergy(calc%bas%nsh))
      if(grad) allocate(dsedcn(calc%bas%nsh), dsedcn_en(calc%bas%nsh))
      call get_selfenergy(calc%h0, mol%id, calc%bas%ish_at, calc%bas%nsh_id, cn=cn, &
         & cn_en=cn_en, selfenergy=selfenergy, dsedcn=dsedcn, dsedcn_en=dsedcn_en)

      allocate(aniso_dip(3, mol%nat), source=0.0_wp)

      cutoff = get_cutoff(calc%bas, accuracy)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      call new_adjacency_list(list, mol, lattr, cutoff)

      if (prlevel > 1) then
         call ctx%message(label_cutoff // format_string(cutoff, real_format) // " bohr")
         call ctx%message("")
      end if

      ! Update basis cache
      call calc%bas%update(mol, bcache, grad)

      ! Get Hamiltonian and integrals
      call new_integral(ints, calc%bas%nao)
      call get_hamiltonian(mol, lattr, list, calc%bas, bcache, calc%h0, selfenergy, &
         & aniso_dip, ints%overlap, ints%dipole, ints%quadrupole, ints%hamiltonian)
      call timer%pop

      ! Get initial potential for external fields and Coulomb
      call new_potential(pot, mol, calc%bas, wfn%nspin)
      ! Set potential to zero
      call pot%reset
      ! Add potential due to external field
      if (allocated(calc%interactions)) then
         call timer%push("interactions")
         call calc%interactions%update(mol, icache)
         call calc%interactions%get_potential(mol, icache, wfn, pot)
         call timer%pop
      endif

      ! Add potential due to Coulomb
      if (allocated(calc%coulomb)) then
         call timer%push("coulomb")
         ! Use electronegativity-weighted CN as 0th-order charge guess
         call get_effective_qat(mol, cn_en, wfn%qat)
      
         call calc%coulomb%update(mol, ccache)
         call calc%coulomb%get_potential(mol, ccache, wfn, pot)
         call timer%pop
      end if

      ! Add effective Hamiltonian to potential
      call add_pot_to_h1(calc%bas, ints, pot, wfn%coeff)

      call timer%push("diagonalization")
      ! Solve the effective Hamiltonian
      call ctx%new_solver(solver, ints%overlap, wfn%nel, wfn%kt)

      ! Get the density matrix
      call next_density(wfn, solver, ints, elec_entropy, error)
      if (allocated(error)) then
         call ctx%set_error(error)
      end if
      call timer%pop

      ! Get charges and dipole moment from density and integrals
      call get_mulliken_shell_charges(calc%bas, ints%overlap, wfn%density, wfn%n0sh, &
      & wfn%qsh)
      call get_qat_from_qsh(calc%bas, wfn%qsh, wfn%qat)
      call get_mulliken_atomic_multipoles(calc%bas, ints%dipole, wfn%density, &
      & wfn%dpat)
      allocate(tmp(3), source = 0.0_wp)
      call gemv(mol%xyz, wfn%qat(:, 1), tmp)
      dipole(:) = tmp + sum(wfn%dpat(:, :, 1), 2)

      call timer%pop
      
      block
         integer :: it
         real(wp) :: ttime, stime
         character(len=*), parameter :: label(*) = [character(len=20):: &
            & "coulomb", "hamiltonian", "diagonalization"]
         if (prlevel > 1) then
            ttime = timer%get("total CEH")
            call ctx%message(" total CEH:"//repeat(" ", 16)//format_time(ttime))
         end if
         if (prlevel > 2) then
            do it = 1, size(label)
               stime = timer%get(label(it))
               if (stime <= epsilon(0.0_wp)) cycle
               call ctx%message(" - "//label(it)//format_time(stime) &
                  & //" ("//format_string(int(stime/ttime*100), '(i3)')//"%)")
            end do
            call ctx%message("")
         end if
      end block

   end subroutine ceh_singlepoint

end module tblite_ceh_singlepoint
