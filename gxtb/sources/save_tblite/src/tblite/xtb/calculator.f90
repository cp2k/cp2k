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

!> @file tblite/xtb/calculator.f90
!> Provides the calculator type for holding xTB Hamiltonian parametrization.

!> Implementation of calculator type for the extended-tight binding Hamiltonian.
!> The #tblite_xtb_calculator::xtb_calculator collects the basic interactions
!> required to perform a tight-binding calculation.
module tblite_xtb_calculator
   use mctc_data_vdwrad, only : get_vdw_rad
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mctc_io_convert, only : aatoau
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count, get_cn_count_id
   use multicharge, only : mchrg_model_type, new_eeq2019_model, new_eeqbc2025_model
   use tblite_acp, only : acp_type, new_acp
   use tblite_basis_ortho, only : orthogonalize
   use tblite_basis_qvszp, only : new_qvszp_cgto, qvszp_cgto_type, & 
      & qvszp_basis_type, new_qvszp_basis
   use tblite_basis_slater, only : slater_to_gauss
   use tblite_basis_type, only : basis_type, new_basis, basis_set, cgto_type, new_cgto, &
      & cgto_container
   use tblite_classical_halogen, only : halogen_correction, new_halogen_correction
   use tblite_classical_increment, only : core_increment, new_core_increment
   use tblite_container, only : container_type, container_list
   use tblite_coulomb_charge, only : coulomb_kernel, new_gamma_coulomb, gamma_coulomb, &
      & new_effective_coulomb, effective_coulomb
   use tblite_coulomb_firstorder, only : new_onsite_firstorder
   use tblite_coulomb_fourthorder, only : new_onsite_fourthorder
   use tblite_coulomb_multipole, only : new_gfn2_multipole, gfn2_multipole, &
      & gxtb_multipole, new_gxtb_multipole, multipole_damping
   use tblite_coulomb_thirdorder, only : new_onsite_thirdorder, onsite_thirdorder, & 
      & new_twobody_thirdorder, twobody_thirdorder, thirdorder_kernel
   use tblite_data_onecxints, only : get_onecxints
   use tblite_data_spin, only : get_spin_constant
   use tblite_disp, only : dispersion_type, d4_dispersion, new_d4_dispersion, &
      & new_d4s_dispersion, new_d4srev_dispersion, d3_dispersion, new_d3_dispersion, &
      & get_damping_function_id
   use tblite_exchange, only : exchange_type, exchange_fock, new_exchange_fock
   use tblite_param, only : param_record
   use tblite_repulsion, only : repulsion_type, gfn_repulsion, gxtb_repulsion, &
      & new_repulsion_gfn, new_repulsion_gxtb, repulsion_kernel
   use tblite_scf_iterator, only : iterator_type, new_iterator
   use tblite_scf_info, only : scf_info, max, orbital_resolved
   use tblite_scf_mixer_broyden, only : broyden_input
   use tblite_scf_mixer_diis, only : diis_input
   use tblite_scf_mixer_input, only : mixer_input_container, mixer_mode
   use tblite_scf_mixer_simple, only : simple_input
   use tblite_spin, only : spin_polarization, new_spin_polarization
   use tblite_utils_average, only : average_type, new_average, average_id
   use tblite_utils_string, only : lowercase
   use tblite_xtb_coulomb, only : tb_coulomb, new_coulomb
   use tblite_xtb_h0, only : tb_hamiltonian, new_hamiltonian
   use tblite_xtb_spec, only : tb_h0spec
   implicit none
   private

   public :: new_xtb_calculator
   public :: param_h0spec

   !> Default maximum number of self-consistent iterations
   integer, parameter :: max_iter_default = 250

   !> Extended tight-binding calculator
   type, public :: xtb_calculator
      !> Basis set definition
      class(basis_type), allocatable :: bas
      !> Core Hamiltonian
      type(tb_hamiltonian) :: h0
      !> Coordination number for modifying the self-energies
      class(ncoord_type), allocatable :: ncoord
      !> Electronegativity-weighted coordination number for modifying the self-energies
      class(ncoord_type), allocatable :: ncoord_en
      !> Auxiliary charge model 
      class(mchrg_model_type), allocatable :: charge_model
      !> Atomic correction potential
      class(acp_type), allocatable :: acp
      !> Collection of all Coulombic interactions
      type(tb_coulomb), allocatable :: coulomb
      !> Exchange interaction
      class(exchange_type), allocatable :: exchange
      !> Repulsion energy interactions
      class(repulsion_type), allocatable :: repulsion
      !> Energy increment for missing core electrons
      type(core_increment), allocatable :: increment
      !> Halogen bonding correction
      type(halogen_correction), allocatable :: halogen
      !> London-dispersion interaction
      class(dispersion_type), allocatable :: dispersion
      !> Spin-polarization interaction
      type(spin_polarization), allocatable :: spin_polarization
      !> List of additional interaction containers
      type(container_list), allocatable :: interactions
      !> Iterator for self-consistent field procedure
      type(iterator_type), allocatable :: iterator
      !> Store calculated integral intermediates
      logical :: save_integrals = .false.
      !> string with method or "custom"
      character(len=:), allocatable :: method
      !> Default electronic temperature
      real(wp) :: default_etemp = 300.0_wp
      !> Default guess method for the SCF
      character(len=:), allocatable :: default_guess
   contains
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
      !> Information on calculator
      procedure :: info
      !> Add an interaction container
      procedure :: push_back
      !> Remove an interaction container
      procedure :: pop
   end type xtb_calculator


   !> Specification of the Hamiltonian
   type, extends(tb_h0spec) :: param_h0spec
      type(param_record), pointer :: param => null()
      integer, pointer :: irc(:) => null()
      logical, allocatable :: valence(:, :)
   contains
      !> Generator for the self energy / atomic levels of the Hamiltonian
      procedure :: get_selfenergy
      !> Generator for the coordination number dependent shift of the self energy
      procedure :: get_cnshift
      !> Generator for the enhancement factor to for scaling Hamiltonian elements
      procedure :: get_hscale
      !> Generator for the atomic radii used in the distant dependent scaling
      procedure :: get_rad
      !> Generator for the polynomial parameters for the distant dependent scaling
      procedure :: get_shpoly
      !> Generator for the polynomial parameters for the linear distant dependent scaling
      procedure :: get_shpoly2
      !> Generator for the polynomial parameters for the square distant dependent scaling
      procedure :: get_shpoly4
      !> Generator for the coefficients of the anisotropic Hamiltonian
      procedure :: get_anisotropy
      !> Generator for the reference occupation numbers of the atoms
      procedure :: get_reference_occ
      !> Generator for the scaled van der Waals radii
      procedure :: get_rvdw
      !> Generator for the diatomic frame scaling factors
      procedure :: get_diat_scale
   end type param_h0spec

   !> Constructor for Hamiltonian specification
   interface param_h0spec
      module procedure :: new_param_h0spec
   end interface param_h0spec


contains


!> Create new xTB Hamiltonian calculator from parametrization data
subroutine new_xtb_calculator(calc, mol, param, error, accuracy)
   !> Instance of the xTB calculator
   type(xtb_calculator), intent(out) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Optional accuracy specification
   real(wp), intent(in), optional :: accuracy

   integer :: isp
   integer, allocatable :: irc(:)

   allocate(irc(mol%nid))

   do isp = 1, mol%nid
      call param%get(mol%sym(isp), mol%num(isp), irc(isp))
      if (irc(isp) == 0) then
         call fatal_error(error, "No entry in parametrization for element "//mol%sym(isp))
         exit
      end if
   end do
   if (allocated(error)) return

   call add_charge_model(calc, mol, param, error)
   if (allocated(error)) return
   call add_basis(calc, mol, param, irc, error, accuracy)
   if (allocated(error)) return
   call add_acp(calc, mol, param, irc, accuracy)
   call add_ncoord(calc, mol, param, irc, error)
   if (allocated(error)) return
   call add_increment(calc, mol, param, irc)
   call add_hamiltonian(calc, mol, param, irc)
   call add_repulsion(calc, mol, param, irc, error)
   if(allocated(error)) return
   call add_halogen(calc, mol, param, irc)
   call add_dispersion(calc, mol, param, error)
   if (allocated(error)) return
   call add_coulomb(calc, mol, param, irc, error)
   if (allocated(error)) return
   call add_exchange(calc, mol, param, irc)
   call add_spin_polarization(calc, mol, param, irc)

   call add_iterator(calc, mol, param, accuracy)

   calc%default_etemp = param%hamiltonian%default_etemp
   calc%default_guess = param%hamiltonian%default_guess

   calc%method = get_method(param%name)

end subroutine new_xtb_calculator


subroutine add_charge_model(calc, mol, param, error)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (allocated(param%hamiltonian%charge_model)) then
      select case(lowercase(param%hamiltonian%charge_model))
      case ("eeq")
         call new_eeq2019_model(mol, calc%charge_model, error)
      case ("eeqbc")
         call new_eeqbc2025_model(mol, calc%charge_model, error)
      case ("none")
         return
      case default
         call fatal_error(error, "Unknown charge model: "//trim(param%hamiltonian%charge_model))
         return
      end select
   end if 
end subroutine add_charge_model


subroutine add_basis(calc, mol, param, irc, error, accuracy)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Optional accuracy specification
   real(wp), intent(in), optional :: accuracy

   integer :: isp, ir, ish, stat, ng, il
   integer, allocatable :: nsh_id(:)
   integer :: ang_idx(0:4), ortho(10)
   type(cgto_container), allocatable :: cgto(:, :)
   type(cgto_container), allocatable :: cgto_h0(:, :)
   type(basis_type), allocatable :: basis
   type(qvszp_basis_type), allocatable :: basis_qvszp
   type(qvszp_cgto_type), allocatable :: cgto_qvszp
   logical :: env_dep

   env_dep = .false.
   nsh_id = param%record(irc)%nsh
   
   allocate(cgto(maxval(nsh_id), mol%nid))
   if (param%hamiltonian%scaled_h0) allocate(cgto_h0(maxval(nsh_id), mol%nid))
   do isp = 1, mol%nid
      ang_idx = 0
      ortho = 0
      ir = irc(isp)
      do ish = 1, nsh_id(isp)
         il = param%record(ir)%record(ish)%lsh
         ng = param%record(ir)%record(ish)%ngauss
         if (ang_idx(il) > 0) then
            ortho(ish) = ang_idx(il)
         else
            ang_idx(il) = ish
         end if

         ! STO-nG basis set
         if (param%record(ir)%record(ish)%basis == basis_set%sto_ng) then
            allocate(cgto_type :: cgto(ish, isp)%raw)
            call slater_to_gauss(ng, param%record(ir)%record(ish)%pqn, il, &
               & param%record(ir)%record(ish)%slater, cgto(ish, isp)%raw, .true., stat)

         ! q-vSZP basis set
         else if (param%record(ir)%record(ish)%basis == basis_set%q_vszp) then
            allocate(cgto_qvszp)
            call new_qvszp_cgto(cgto_qvszp, ir, ish, .true., error, &
               & param%record(ir)%record(ish)%expos, &
               & param%record(ir)%record(ish)%coeffs, &
               & param%record(ir)%record(ish)%coeffs_env, &
               & param%record(ir)%k0, param%record(ir)%k1, &
               & param%record(ir)%k2, param%record(ir)%k3)
            call move_alloc(cgto_qvszp, cgto(ish, isp)%raw)
            env_dep = .true.

         ! Custom basis set
         else
            allocate(cgto_type :: cgto(ish, isp)%raw)
            call new_cgto(cgto(ish, isp)%raw, ng, il, param%record(ir)%record(ish)%expos, &
                & param%record(ir)%record(ish)%coeffs, .true.)

         end if

         ! Setup the scaled CGTOs for H0 construction
         if (param%hamiltonian%scaled_h0) then
            allocate(cgto_qvszp)
            call new_qvszp_cgto(cgto_qvszp, ir, ish, .true., error, &
               & param%record(ir)%record(ish)%expos * param%record(ir)%record(ish)%h0_exp_scale, &
               & param%record(ir)%record(ish)%coeffs, &
               & param%record(ir)%record(ish)%coeffs_env, &
               & param%record(ir)%k0 * param%record(ir)%k0_h0_scale, &
               & param%record(ir)%k1 * param%record(ir)%k1_h0_scale, &
               & param%record(ir)%k2 * param%record(ir)%k2_h0_scale, &
               & param%record(ir)%k3 * param%record(ir)%k3_h0_scale)
            if (allocated(error)) return
            call move_alloc(cgto_qvszp, cgto_h0(ish, isp)%raw)
         end if
      end do

      ! Orthonormalize CGTOs with the same angular momentum
      do ish = 1, nsh_id(isp)
         if (ortho(ish) > 0) then
            call orthogonalize(cgto(ortho(ish), isp)%raw, cgto(ish, isp)%raw)
            if (param%hamiltonian%scaled_h0) then
               call orthogonalize(cgto_h0(ortho(ish), isp)%raw, cgto_h0(ish, isp)%raw)
            end if
         end if
      end do
   end do

   ! Setup (charge-dependent or default) basis set 
   if (env_dep) then
      allocate(basis_qvszp)
      call new_qvszp_basis(basis_qvszp, mol, nsh_id, cgto, error, accuracy, cgto_h0)
      if (allocated(error)) return
      call move_alloc(basis_qvszp, calc%bas)
   else
      allocate(basis)
      call new_basis(basis, mol, nsh_id, cgto, accuracy, cgto_h0)
      call move_alloc(basis, calc%bas)
   end if

end subroutine add_basis


subroutine add_acp(calc, mol, param, irc, accuracy)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Optional accuracy specification
   real(wp), intent(in), optional :: accuracy

   integer :: isp, ir, iproj, lproj, mproj
   type(cgto_container), allocatable :: cgtp(:, :)
   real(wp) :: alpha(1), coeff(1)
   real(wp), allocatable :: levels(:, :)
   integer, allocatable :: nproj(:)

   if (allocated(param%acp)) then
      nproj = param%record(irc)%n_acp
      mproj = maxval(nproj)
      coeff(1) = 1.0_wp
      allocate(cgtp(mproj, mol%nid))
      allocate(levels(mproj, mol%nid), source=0.0_wp)
      do isp = 1, mol%nid
         ir = irc(isp)
         do iproj = 1, nproj(isp)
            lproj = param%record(ir)%l_acp(iproj)
            alpha(1) = param%record(ir)%acp_expos(iproj)
            levels(iproj, isp) = param%record(ir)%acp_levels(iproj)
            allocate(cgto_type :: cgtp(iproj, isp)%raw)
            call new_cgto(cgtp(iproj, isp)%raw, 1, lproj, alpha, coeff, .true.)
         end do
      end do

      allocate(calc%acp)
      call new_acp(calc%acp, mol, nproj, cgtp, levels, accuracy)
   end if

end subroutine add_acp


subroutine add_ncoord(calc, mol, param, irc, error)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: cn_count_type
   real(wp), allocatable :: cn_rcov(:)

   if (allocated(param%hamiltonian%cn)) then
      cn_count_type = get_cn_count_id(param%hamiltonian%cn)
      cn_rcov = param%record(irc)%cn_rcov
      call new_ncoord(calc%ncoord, mol, cn_count_type=cn_count_type, &
         & error=error, kcn=param%hamiltonian%cn_exp, rcov=cn_rcov)
   end if
end subroutine add_ncoord


subroutine add_increment(calc, mol, param, irc)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)

   real(wp), allocatable :: increment(:)

   if (allocated(param%increment)) then
      increment = param%record(irc)%increment
      allocate(calc%increment)
      call new_core_increment(calc%increment, mol, increment)
   end if
end subroutine add_increment


subroutine add_hamiltonian(calc, mol, param, irc)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)

   call new_hamiltonian(calc%h0, mol, calc%bas, new_param_h0spec(mol, param, irc))
end subroutine add_hamiltonian


subroutine add_dispersion(calc, mol, param, error)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: damping_2b_id, damping_3b_id
   type(d4_dispersion), allocatable :: d4
   type(d3_dispersion), allocatable :: d3

   if (.not.allocated(param%dispersion)) return
   associate(par => param%dispersion)
      if (par%d3) then
         allocate(d3)
         call new_d3_dispersion(d3, mol, s6=par%s6, s8=par%s8, &
            & a1=par%a1, a2=par%a2, s9=par%s9, error=error)
         call move_alloc(d3, calc%dispersion)
      else
         allocate(d4)
         call get_damping_function_id(error, par%damping_2b, par%damping_3b, &
            & damping_2b_id, damping_3b_id)
         if (par%smooth) then
            if (par%rev) then
               call new_d4srev_dispersion(d4, mol, damping_2b=damping_2b_id, &
                  & damping_3b=damping_3b_id, s6=par%s6, s8=par%s8, &
                  & a1=par%a1, a2=par%a2, a3=par%a3, a4=par%a4, s9=par%s9, &
                  & kcn=par%kcn, error=error)
            else 
               call new_d4s_dispersion(d4, mol, damping_2b=damping_2b_id, &
                  & damping_3b=damping_3b_id, s6=par%s6, s8=par%s8, &
                  & a1=par%a1, a2=par%a2, a3=par%a3, a4=par%a4, s9=par%s9, &
                  & error=error)
            end if 
         else
            call new_d4_dispersion(d4, mol, damping_2b=damping_2b_id, &
               & damping_3b=damping_3b_id, s6=par%s6, s8=par%s8, &
               & a1=par%a1, a2=par%a2, a3=par%a3, a4=par%a4, s9=par%s9, &
               & error=error)
         end if
         call move_alloc(d4, calc%dispersion)
      end if
   end associate
end subroutine add_dispersion

subroutine add_repulsion(calc, mol, param, irc, error)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(gfn_repulsion), allocatable :: gfn
   type(gxtb_repulsion), allocatable :: gxtb

   integer :: isp, jsp, ir, jr
   real(wp), allocatable :: alpha(:), zeff(:)
   real(wp), allocatable :: rep_cn(:), rep_q(:), rep_q2(:), rep_k1(:)
   real(wp), allocatable :: rep_roffset(:), cn_rcov(:), rvdw(:, :)
   type(average_type), allocatable :: average

   select case(param%repulsion%kernel)
   case(repulsion_kernel%gfn)
      block
         type(gfn_repulsion), allocatable :: repulsion
         zeff = param%record(irc)%zeff
         alpha = param%record(irc)%alpha
         allocate(repulsion)
         call new_repulsion_gfn(repulsion, mol, alpha=alpha, zeff=zeff, &
            & kexp=param%repulsion%kexp, kexp_light=param%repulsion%klight, &
            & rexp=param%repulsion%rexp) 
         call move_alloc(repulsion, calc%repulsion)
      end block
   case(repulsion_kernel%gxtb)
      block
         type(gxtb_repulsion), allocatable :: repulsion
         zeff = param%record(irc)%zeff
         alpha = param%record(irc)%alpha
         rep_cn = param%record(irc)%rep_cn
         rep_q = param%record(irc)%rep_q
         rep_roffset = param%record(irc)%rep_roffset
         rep_k1 = param%record(irc)%rep_k1
         cn_rcov = param%record(irc)%cn_rcov
         allocate(average)
         call new_average(average, average_id%arithmetic)
         allocate(rvdw(mol%nid, mol%nid))
         do isp = 1, mol%nid
            ir = irc(isp)
            do jsp = 1, isp
               jr = irc(jsp)
               rvdw(jsp, isp) = get_vdw_rad(jr, ir) &
                  & * average%value(param%record(ir)%rvdw_scale, &
                  & param%record(jr)%rvdw_scale)
               rvdw(isp, jsp) = rvdw(jsp, isp)
            end do
         end do
         allocate(repulsion)
         call new_repulsion_gxtb(repulsion, mol, alpha=alpha, zeff=zeff, &
            & kcn=rep_cn, kq=rep_q, roffset=rep_roffset, kexp=param%repulsion%kexp, &
            & k1=rep_k1, k2=param%repulsion%k2, k2_light=param%repulsion%k2light, &
            & k3=param%repulsion%k3, k4=param%repulsion%k4, &
            & kshort=param%repulsion%short, kshort_alpha=param%repulsion%short_alpha, &
            & kshort_exp=param%repulsion%short_exp, rvdw=rvdw, cn_rcov=cn_rcov, &
            & cn_exp=param%hamiltonian%cn_exp, error=error)
         if (allocated(error)) return
         call move_alloc(repulsion, calc%repulsion)
      end block
   end select

end subroutine add_repulsion


subroutine add_halogen(calc, mol, param, irc)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)

   real(wp), allocatable :: bond_strength(:)

   if (.not.allocated(param%halogen)) return
   allocate(calc%halogen)
   bond_strength = param%record(irc)%xbond
   associate(par => param%halogen)
      call new_halogen_correction(calc%halogen, mol, par%damping, par%rscale, &
         & bond_strength)
   end associate
end subroutine add_halogen

subroutine add_coulomb(calc, mol, param, irc, error)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), allocatable :: refqsh(:), ipea(:, :), hardness(:, :), &
      & hubbard_derivs(:, :), hubbard_second_derivs(:, :)
   real(wp), allocatable :: tb_cn_rcov(:), ipea_cn(:), hardness_cn(:)
   real(wp), allocatable :: dkernel(:), qkernel(:), rad(:), vcn(:)
   real(wp), allocatable :: rvdw(:, :), aes_dip_scale(:)
   
   integer :: cn_count_type, isp, jsp, ir, jr
   type(average_type), allocatable :: average_es2, average_es3, average_aes2

   if (.not.any([allocated(param%firstorder), allocated(param%charge), &
      & allocated(param%thirdorder), allocated(param%fourthorder), &
      & allocated(param%multipole)])) return

   tb_cn_rcov = param%record(irc)%cn_rcov

   ! Make sure that the GFN2 multipole CN has the right (default) covalent radii
   allocate(calc%coulomb)
   if (allocated(param%hamiltonian%cn)) then
      cn_count_type = get_cn_count_id(param%hamiltonian%cn)
   end if
   ! Collect the valence coordination number in case on CN was defined 
   if (allocated(param%multipole)) then
      vcn = param%record(irc)%mpvcn
   end if
   call new_coulomb(calc%coulomb, mol, error, cn_count_type=cn_count_type, &
      & cn_rcov=tb_cn_rcov, cn_exp=param%hamiltonian%cn_exp, valence_cn=vcn)

   if (allocated(param%firstorder)) then
      allocate(calc%coulomb%es1)
      call get_shell_ipea(mol, param, irc, ipea)
      ipea_cn = param%record(irc)%ipea_cn
      call get_reference_shell_charge(mol, param, irc, calc%bas%ish_at, refqsh)
      call new_onsite_firstorder(calc%coulomb%es1, mol, ipea, ipea_cn, &
         & param%firstorder%split_exp, param%firstorder%split_slope, &
         & param%firstorder%split_offset, calc%bas%nsh_id)
   end if

   if (allocated(param%charge)) then
      call get_shell_hardness(mol, param, irc, hardness)
      hardness_cn = param%record(irc)%gam_cn
      select case(param%charge%kernel)
      case(coulomb_kernel%effective)
         block
            type(effective_coulomb), allocatable :: es2
            allocate(es2)
            call get_average(param%charge%average, average_es2)
            call new_effective_coulomb(es2, mol, param%charge%gexp, hardness, &
               & average_es2, hardness_cn, param%charge%hubbard_exp, refqsh, &
               & calc%bas%nsh_id)
            call move_alloc(es2, calc%coulomb%es2)
         end block
      case(coulomb_kernel%dftb_gamma)
         block
            type(gamma_coulomb), allocatable :: es2
            allocate(es2)
            call new_gamma_coulomb(es2, mol, hardness, refqsh, calc%bas%nsh_id)
            call move_alloc(es2, calc%coulomb%es2)
         end block
      end select
   end if

   if (allocated(param%thirdorder)) then
      if (param%thirdorder%shell) then
         call get_hubbard_derivs(mol, param, irc, hubbard_derivs)
      else
         allocate(hubbard_derivs(1, mol%nid))
         hubbard_derivs(1, :) = param%record(irc)%gam3
      end if
      select case(param%thirdorder%kernel)
      case(thirdorder_kernel%onsite)
         block
            type(onsite_thirdorder), allocatable :: es3
            allocate(es3)
            if (param%thirdorder%shell) then
               call new_onsite_thirdorder(es3, mol, hubbard_derivs, calc%bas%nsh_id)
            else
               call new_onsite_thirdorder(es3, mol, hubbard_derivs)
            end if
            call move_alloc(es3, calc%coulomb%es3)
         end block
      case(thirdorder_kernel%twobody)  
         block
            type(twobody_thirdorder), allocatable :: es3
            if (.not.allocated(hardness)) then
               call get_shell_hardness(mol, param, irc, hardness)
            end if
            if (.not.allocated(hardness_cn)) then
               hardness_cn = param%record(irc)%gam_cn
            end if
            call get_average(param%thirdorder%average, average_es3)
            allocate(es3)
            if (param%thirdorder%shell) then
               call new_twobody_thirdorder(es3, mol, param%thirdorder%texp, &
                  & param%thirdorder%onsite_scale, param%thirdorder%offsite_scale, & 
                  & hardness, average_es3, hubbard_derivs, hardness_cn, calc%bas%nsh_id)
            else
               call new_twobody_thirdorder(es3, mol, param%thirdorder%texp, &
                  & param%thirdorder%onsite_scale, param%thirdorder%offsite_scale, &
                  & hardness, average_es3, hubbard_derivs, hardness_cn)
            end if
            call move_alloc(es3, calc%coulomb%es3)
         end block
      end select
   end if

   if (allocated(param%fourthorder)) then
      allocate(calc%coulomb%es4)
      if (param%fourthorder%shell) then
         call get_hubbard_second_derivs(mol, param, irc, hubbard_second_derivs)
         call new_onsite_fourthorder(calc%coulomb%es4, mol, hubbard_second_derivs, &
            & calc%bas%nsh_id)
      else
         allocate(hubbard_second_derivs(1, mol%nid))
         hubbard_second_derivs(1, :) = param%record(irc)%gam4
         call new_onsite_fourthorder(calc%coulomb%es4, mol, hubbard_second_derivs)
      end if
   end if

   if (allocated(param%multipole)) then
      call get_average(param%multipole%average, average_aes2)
      select case(param%multipole%damping)
      case (multipole_damping%gfn2)
         block
            type(gfn2_multipole), allocatable :: aes2
            dkernel = param%record(irc)%dkernel
            qkernel = param%record(irc)%qkernel
            rad = param%record(irc)%mprad
            vcn = param%record(irc)%mpvcn
            allocate(aes2)
            associate(par => param%multipole)
               call new_gfn2_multipole(aes2, mol, dkernel, qkernel, &
                  & par%shift, par%kradexp, par%rmax, rad, average_aes2, vcn, &
                  & par%kdmp3, par%kdmp5, par%kdmp7, par%kdmp9, &
                  & par%kexp3, par%kexp5, par%kexp7, par%kexp9)
            end associate
            call move_alloc(aes2, calc%coulomb%aes2)
         end block
      case (multipole_damping%gxtb)
         block
            type(gxtb_multipole), allocatable :: aes2
            allocate(rvdw(mol%nid, mol%nid))
            aes_dip_scale = param%record(irc)%aes_dip_scale
            do isp = 1, mol%nid
               ir = irc(isp)
               do jsp = 1, isp
                  jr = irc(jsp)
                  rvdw(jsp, isp) = get_vdw_rad(jr, ir) &
                     & * average_aes2%value(param%record(ir)%rvdw_scale, &
                     & param%record(jr)%rvdw_scale)
                  rvdw(isp, jsp) = rvdw(jsp, isp)
               end do
            end do
            allocate(aes2)
            associate(par => param%multipole)
               call new_gxtb_multipole(aes2, mol, rvdw, aes_dip_scale, &
                  & par%kdmp3,par%kdmp5, par%kdmp7, par%kdmp9, &
                  & par%kexp3, par%kexp5, par%kexp7, par%kexp9)
            end associate
            call move_alloc(aes2, calc%coulomb%aes2)
         end block
      end select
   end if

end subroutine add_coulomb


subroutine add_exchange(calc, mol, param, irc)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)

   real(wp), allocatable :: hardness(:, :), onecxints(:, :, :), offdiag(:, :)
   real(wp), allocatable :: avg_exponents(:, :), rad(:, :), kq(:, :)
   real(wp), allocatable :: cscale(:), crad(:)
   integer :: isp, jsp, ir, jr
   type(average_type), allocatable :: rvdw_average
   type(average_type), allocatable :: hubbard_average
   type(average_type), allocatable :: offdiag_average
   type(average_type), allocatable :: cscale_average
   type(average_type), allocatable :: crad_average

   if (.not.allocated(param%exchange)) return

   ! Obtain radii for radius-dependent hubbard scaling
   allocate(rad(mol%nid, mol%nid))
   allocate(rvdw_average)
   call new_average(rvdw_average, average_id%arithmetic)
   associate(record => param%record)
      do isp = 1, mol%nid
         ir = irc(isp)
         do jsp = 1, isp
            jr = irc(jsp)
            rad(jsp, isp) = get_vdw_rad(jr, ir) * aatoau &
               & * rvdw_average%value(record(ir)%rvdw_scale, record(jr)%rvdw_scale)
            rad(isp, jsp) = rad(jsp, isp)
         end do
      end do
   end associate

   call get_shell_hardness_exchange(mol, param, irc, hardness)
   call get_exchange_averaging_exponents(mol, param, irc, avg_exponents)
   call get_onecenter_exchange_ints(mol, param, irc, onecxints)
   call get_offdiagonal_exchange_scaling(mol, param, irc, offdiag)
   call get_charge_exchange_fraction_scaling(mol, param, irc, kq)
   cscale = param%record(irc)%cscale
   crad = param%record(irc)%crad
   associate(par => param%exchange)
      call get_average(par%hubbard_average, hubbard_average)
      call get_average(par%offdiag_average, offdiag_average)
      call get_average(par%cscale_average, cscale_average)
      call get_average(par%crad_average, crad_average)
      block 
         type(exchange_fock), allocatable :: tmp
         allocate(tmp)
         call new_exchange_fock(tmp, mol, calc%bas, hardness, hubbard_average, &
            & avg_exponents, par%ondiag, offdiag, offdiag_average, par%hubbard_exp, &
            & par%hubbard_exp_r0, rad, par%gexp, onecxints, kq, cscale, &
            & cscale_average, par%cexp, crad, crad_average, par%frscale, &
            & par%omega, par%lrscale)
         call move_alloc(tmp, calc%exchange)
      end block
   end associate

end subroutine add_exchange


subroutine add_spin_polarization(calc, mol, param, irc)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)

   real(wp), allocatable :: wll(:, :, :)

   if (.not. allocated(param%spin)) return

   call get_spin_constants(mol, param, irc, wll)
   allocate(calc%spin_polarization)
   call new_spin_polarization(calc%spin_polarization, mol, wll, calc%bas%nsh_id)

end subroutine add_spin_polarization


subroutine add_iterator(calc, mol, param, accuracy)
   !> Instance of the xTB evaluator
   type(xtb_calculator), intent(inout) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Optional accuracy specification
   real(wp), intent(in), optional :: accuracy

   type(scf_info) :: info
   type(mixer_input_container), allocatable :: mixer_inputs(:)

   info = calc%variable_info()

   if (info%density == orbital_resolved) then
      allocate(mixer_inputs(2))
      allocate(mixer_inputs(1)%raw, source=simple_input(mode=mixer_mode%potential, &
         & damp=0.2_wp))
      allocate(mixer_inputs(2)%raw, source=diis_input(mode=mixer_mode%potential, &
         & output_fraction=1.0_wp))
   else
      allocate(mixer_inputs(2))
      allocate(mixer_inputs(1)%raw, source=simple_input(mode=mixer_mode%density))
      allocate(mixer_inputs(2)%raw, source=broyden_input(memory=max_iter_default, &
         & mode=mixer_mode%density))
   end if

   allocate(calc%iterator)
   call new_iterator(calc%iterator, mol, calc%bas, mixer_inputs, &
      & max_iter=max_iter_default, accuracy=accuracy)

end subroutine add_iterator


subroutine get_average(average_name, averager)
   !> Name of the averaging scheme
   character(len=*), intent(in) :: average_name
   !> Averaging object
   type(average_type), allocatable, intent(out) :: averager
   integer :: id

   select case(average_name)
   case default
      id = -1 
   case("arithmetic")
      id = average_id%arithmetic
   case("geometric")
      id = average_id%geometric
   case("harmonic")
      id = average_id%harmonic
   case("general")
      id = average_id%general
   end select

   allocate(averager)
   call new_average(averager, id)
end subroutine get_average


subroutine get_shell_hardness(mol, param, irc, hardness)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Shell resolved hardness parameters
   real(wp), allocatable, intent(out) :: hardness(:, :)

   integer :: isp, ir, ish, il

   allocate(hardness(maxval(param%record(irc)%nsh), mol%nid))
   hardness(:, :) = 0.0_wp
   do isp = 1, mol%nid
      ir = irc(isp)
      do ish = 1, param%record(ir)%nsh
         il = param%record(ir)%record(ish)%lsh
         hardness(ish, isp) = param%record(ir)%gam &
            & * param%record(ir)%record(ish)%lgam
      end do
   end do
end subroutine get_shell_hardness


subroutine get_shell_hardness_exchange(mol, param, irc, hardness_fx)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Shell resolved hardness parameters for exchange
   real(wp), allocatable, intent(out) :: hardness_fx(:, :)

   integer :: isp, ir, ish

   allocate(hardness_fx(maxval(param%record(irc)%nsh), mol%nid))
   hardness_fx(:, :) = 0.0_wp
   do isp = 1, mol%nid
      ir = irc(isp)
      do ish = 1, param%record(ir)%nsh
         hardness_fx(ish, isp) = param%record(ir)%gam &
            & * param%record(ir)%record(ish)%lgam_fx
      end do
   end do


end subroutine get_shell_hardness_exchange


subroutine get_exchange_averaging_exponents(mol, param, irc, avg_exponents)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Shell resolved averaging exponents for exchange
   real(wp), allocatable, intent(out) :: avg_exponents(:, :)

   integer :: isp, ir, ish

   allocate(avg_exponents(maxval(param%record(irc)%nsh), mol%nid))
   avg_exponents(:, :) = 0.0_wp
   do isp = 1, mol%nid
      ir = irc(isp)
      do ish = 1, param%record(ir)%nsh
         avg_exponents(ish, isp) = param%record(ir)%record(ish)%avg_exp_fx
      end do
   end do

end subroutine get_exchange_averaging_exponents


subroutine get_onecenter_exchange_ints(mol, param, irc, onecxints)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Shell resolved averaging exponents for exchange
   real(wp), allocatable, intent(out) :: onecxints(:, :, :)

   integer :: isp, ir, ish, jsh, il, jl

   allocate(onecxints(maxval(param%record(irc)%nsh), &
      & maxval(param%record(irc)%nsh), mol%nid))
   onecxints(:, :, :) = 0.0_wp
   do isp = 1, mol%nid
      ir = irc(isp)
      do ish = 1, param%record(ir)%nsh
         il = param%record(ir)%record(ish)%lsh
         do jsh = 1, param%record(ir)%nsh
            jl = param%record(ir)%record(jsh)%lsh
            onecxints(jsh, ish, isp) = get_onecxints(jl, il, ir)
         end do
      end do
   end do

end subroutine get_onecenter_exchange_ints


subroutine get_offdiagonal_exchange_scaling(mol, param, irc, offdiag)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Shell resolved off-diagonal exchange scaling factors
   real(wp), allocatable, intent(out) :: offdiag(:, :)

   integer :: isp, ir, ish, il

   allocate(offdiag(maxval(param%record(irc)%nsh), mol%nid))
   offdiag(:, :) = 0.0_wp
   do isp = 1, mol%nid
      ir = irc(isp)
      do ish = 1, param%record(ir)%nsh
         il = param%record(ir)%record(ish)%lsh
         offdiag(ish, isp) = param%exchange%offdiag_l(il)
      end do
   end do

end subroutine get_offdiagonal_exchange_scaling


subroutine get_charge_exchange_fraction_scaling(mol, param, irc, kq)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Shell resolved charge scaling of the exchange fraction
   real(wp), allocatable, intent(out) :: kq(:, :)

   integer :: isp, ir, ish, il

   allocate(kq(maxval(param%record(irc)%nsh), mol%nid))
   kq(:, :) = 0.0_wp
   do isp = 1, mol%nid
      ir = irc(isp)
      do ish = 1, param%record(ir)%nsh
         il = param%record(ir)%record(ish)%lsh
         kq(ish, isp) = param%exchange%kq(il)
      end do
   end do

end subroutine get_charge_exchange_fraction_scaling


subroutine get_spin_constants(mol, param, irc, wll)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Scaled spin constants
   real(wp), allocatable, intent(out) :: wll(:, :, :)

   integer :: isp, ir, ish, jsh, il, jl

   allocate(wll(maxval(param%record(irc)%nsh), &
      & maxval(param%record(irc)%nsh), mol%nid))
   wll(:, :, :) = 0.0_wp

   do isp = 1, mol%nid
      ir = irc(isp)
      do ish = 1, param%record(ir)%nsh
         il = param%record(ir)%record(ish)%lsh
         do jsh = 1, param%record(ir)%nsh
            jl = param%record(ir)%record(jsh)%lsh
            wll(jsh, ish, isp) = param%record(ir)%wll_scale &
               & * get_spin_constant(jl, il, ir, wb97mv=param%spin%wb97mv)
         end do
      end do
   end do

end subroutine get_spin_constants


subroutine get_shell_ipea(mol, param, irc, ipea)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Shell resolved IP/EA chemical potential parameters
   real(wp), allocatable, intent(out) :: ipea(:, :)

   integer :: isp, ir, ish

   allocate(ipea(maxval(param%record(irc)%nsh), mol%nid))
   ipea(:, :) = 0.0_wp
   do isp = 1, mol%nid
      ir = irc(isp)
      do ish = 1, param%record(ir)%nsh
         ipea(ish, isp) = param%record(ir)%record(ish)%ipea
      end do
   end do
end subroutine get_shell_ipea


subroutine get_reference_shell_charge(mol, param, irc, ish_at, refqsh)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Index offset for each atom in the shell space
   integer, intent(in) :: ish_at(:)
   !> Reference occupation numbers
   real(wp), allocatable, intent(out) :: refqsh(:)

   integer :: iat, isp, ir, ish, ii, nsh

   nsh = 0
   do iat = 1, mol%nat
      isp = mol%id(iat)
      ir = irc(isp)
      nsh = nsh + param%record(ir)%nsh
   end do

   allocate(refqsh(nsh))
   refqsh(:) = 0.0_wp
   do iat = 1, mol%nat
      isp = mol%id(iat)
      ir = irc(isp)
      ii = ish_at(iat)
      do ish = 1, param%record(ir)%nsh
         refqsh(ii + ish) = param%record(ir)%record(ish)%refocc &
            & - param%record(ir)%record(ish)%zeffsh
      end do
   end do
end subroutine get_reference_shell_charge


subroutine get_hubbard_derivs(mol, param, irc, hubbard_derivs)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Shell resolved Hubbard derivatives
   real(wp), allocatable, intent(out) :: hubbard_derivs(:, :)

   integer :: isp, ir, ish, il

   allocate(hubbard_derivs(maxval(param%record(irc)%nsh), mol%nid))
   hubbard_derivs(:, :) = 0.0_wp
   do isp = 1, mol%nid
      ir = irc(isp)
      do ish = 1, param%record(ir)%nsh
         il = param%record(ir)%record(ish)%lsh
         hubbard_derivs(ish, isp) = param%record(ir)%gam3 * param%thirdorder%ksh(il)
      end do
   end do
end subroutine get_hubbard_derivs


subroutine get_hubbard_second_derivs(mol, param, irc, hubbard_second_derivs)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), intent(in) :: param
   !> Record identifiers
   integer, intent(in) :: irc(:)
   !> Shell resolved Hubbard second derivatives
   real(wp), allocatable, intent(out) :: hubbard_second_derivs(:, :)

   integer :: isp, ir, ish, il

   allocate(hubbard_second_derivs(maxval(param%record(irc)%nsh), mol%nid))
   hubbard_second_derivs(:, :) = 0.0_wp
   do isp = 1, mol%nid
      ir = irc(isp)
      do ish = 1, param%record(ir)%nsh
         il = param%record(ir)%record(ish)%lsh
         hubbard_second_derivs(ish, isp) = param%record(ir)%gam4 &
            & * param%fourthorder%ksh(il)
      end do
   end do
end subroutine get_hubbard_second_derivs


function new_param_h0spec(mol, param, irc) result(self)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Parametrization records
   type(param_record), target, intent(in) :: param
   !> Record identifiers
   integer, target, intent(in) :: irc(:)
   !> Instance of the Hamiltonian specification
   type(param_h0spec) :: self

   integer :: isp, il, ir, ish
   integer :: ang_idx(0:4)

   self%param => param
   self%irc => irc

   allocate(self%valence(maxval(param%record(irc)%nsh), mol%nid))
   do isp = 1, mol%nid
      ang_idx = 0
      ir = irc(isp)
      do ish = 1, param%record(ir)%nsh
         il = param%record(ir)%record(ish)%lsh
         self%valence(ish, isp) = ang_idx(il) == 0
         if (self%valence(ish, isp)) ang_idx(il) = ish
      end do
   end do
end function new_param_h0spec


!> Generator for the enhancement factor to for scaling Hamiltonian elements
subroutine get_hscale(self, mol, bas, hscale)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Scaling parameters for the Hamiltonian elements
   real(wp), intent(out) :: hscale(:, :, :, :)

   integer :: isp, jsp, ish, jsh, il, jl, ir, jr
   real(wp) :: zi, zj, zij, den, enp, km

   hscale(:, :, :, :) = 0.0_wp

   associate(par => self%param%hamiltonian, record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         ir = irc(isp)
         do jsp = 1, mol%nid
            jr = irc(jsp)
            den = (record(ir)%en - record(jr)%en)**2
            do ish = 1, bas%nsh_id(isp)
               il = bas%cgto(ish, isp)%raw%ang
               do jsh = 1, bas%nsh_id(jsp)
                  jl = bas%cgto(jsh, jsp)%raw%ang
                  ! Do Slater-exponent bast scaling only for STO-nG basis functions
                  if (record(ir)%record(ish)%basis == basis_set%sto_ng) then
                     zi = record(ir)%record(ish)%slater
                     zj = record(jr)%record(jsh)%slater
                     if (abs(par%wexp) < epsilon(par%wexp)) then
                        zij = 1.0_wp
                     else
                        zij = (2*sqrt(zi*zj)/(zi+zj))**par%wexp
                     end if
                  else
                     zij = 1.0_wp
                  end if
                  if (self%valence(ish, isp) .and. self%valence(jsh, jsp)) then
                     enp = 1.0_wp + par%enscale * den
                     km = par%kpair(jr, ir) * par%ksh(jl, il) * enp
                  else if (self%valence(ish, isp) .and. .not.self%valence(jsh, jsp)) then
                     km = 0.5_wp * (par%ksh(il, il) + par%kpol)
                  else if (.not.self%valence(ish, isp) .and. self%valence(jsh, jsp)) then
                     km = 0.5_wp * (par%ksh(jl, jl) + par%kpol)
                  else
                     km = par%kpol
                  end if
                  hscale(jsh, ish, jsp, isp) = zij * km
               end do
            end do
         end do
      end do
   end associate
end subroutine get_hscale


!> Generator for the self energy / atomic levels of the Hamiltonian
subroutine get_selfenergy(self, mol, bas, selfenergy)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Self energy / atomic levels
   real(wp), intent(out) :: selfenergy(:, :)

   integer :: isp, ir, ish

   selfenergy(:, :) = 0.0_wp

   associate(record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         ir = irc(isp)
         do ish = 1, bas%nsh_id(isp)
            selfenergy(ish, isp) = record(ir)%record(ish)%level
         end do
      end do
   end associate
end subroutine get_selfenergy


!> Generator of the coordination number dependent shift of the self energy
subroutine get_cnshift(self, mol, bas, kcn)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Coordination number dependent shift
   real(wp), intent(out) :: kcn(:, :)

   integer :: isp, ir, ish

   kcn(:, :) = 0.0_wp

   associate(record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         ir = irc(isp)
         do ish = 1, bas%nsh_id(isp)
            kcn(ish, isp) = record(ir)%record(ish)%kcn
         end do
      end do
   end associate
end subroutine get_cnshift


!> Generator for the atomic radii used in the distant dependent scaling
subroutine get_rad(self, mol, bas, rad)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Atomic radii
   real(wp), intent(out) :: rad(:)

   integer :: isp

   rad(:) = 0.0_wp

   associate(record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         rad(isp) = record(irc(isp))%h0_rad
      end do
   end associate
end subroutine get_rad


!> Generator for the polynomial parameters for the square-root distant dependent scaling
subroutine get_shpoly(self, mol, bas, shpoly)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Polynomial parameters for square-root distant dependent scaleing
   real(wp), intent(out) :: shpoly(:, :)

   integer :: isp, ir, ish

   shpoly(:, :) = 0.0_wp

   associate(record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         ir = irc(isp)
         do ish = 1, bas%nsh_id(isp)
            shpoly(ish, isp) = record(ir)%record(ish)%shpoly
         end do
      end do
   end associate
end subroutine get_shpoly


!> Generator for the polynomial parameters for the linear distant dependent scaling
subroutine get_shpoly2(self, mol, bas, shpoly2)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Polynomial parameters for linear distant dependent scaleing
   real(wp), intent(out) :: shpoly2(:, :)

   integer :: isp, ir, ish

   shpoly2(:, :) = 0.0_wp

   associate(record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         ir = irc(isp)
         do ish = 1, bas%nsh_id(isp)
            shpoly2(ish, isp) = record(ir)%record(ish)%shpoly2
         end do
      end do
   end associate
end subroutine get_shpoly2


!> Generator for the polynomial parameters for the square distant dependent scaling
subroutine get_shpoly4(self, mol, bas, shpoly4)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Polynomial parameters for square distant dependent scaleing
   real(wp), intent(out) :: shpoly4(:, :)

   integer :: isp, ir, ish

   shpoly4(:, :) = 0.0_wp

   associate(record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         ir = irc(isp)
         do ish = 1, bas%nsh_id(isp)
            shpoly4(ish, isp) = record(ir)%record(ish)%shpoly4
         end do
      end do
   end associate
end subroutine get_shpoly4


!> Generator for the coefficients of the anisotropic Hamiltonian
subroutine get_anisotropy(self, mol, bas, dip_scale, aniso_exp)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Parameters for dipole Hamiltonian correction
   real(wp), intent(out) :: dip_scale(:, :)
   !> Exponent for the anisotropic damping function
   real(wp), intent(out) :: aniso_exp

   type(average_type), allocatable :: average
   integer :: isp, jsp, ir, jr

   dip_scale(:, :) = 0.0_wp
   aniso_exp = self%param%hamiltonian%aniso_exp

   allocate(average)
   call new_average(average, average_id%arithmetic)

   associate(record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         ir = irc(isp)
         do jsp = 1, mol%nid
            jr = irc(jsp)
            dip_scale(jsp, isp) = average%value(record(ir)%h0_dip_scale, &
               & record(jr)%h0_dip_scale)
         end do
      end do
   end associate
end subroutine get_anisotropy


!> Generator for the scaled van der Waals radii
subroutine get_rvdw(self, mol, bas, rvdw)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Pairwise van der Waals radii
   real(wp), intent(out) :: rvdw(:, :)

   type(average_type), allocatable :: average
   integer :: isp, jsp, ir, jr

   allocate(average)
   call new_average(average, average_id%arithmetic)

   associate(record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         ir = irc(isp)
         do jsp = 1, isp
            jr = irc(jsp)
            rvdw(jsp, isp) = get_vdw_rad(jr, ir) &
               & * average%value(record(ir)%rvdw_scale, record(jr)%rvdw_scale)
            rvdw(isp, jsp) = rvdw(jsp, isp)
         end do
      end do
   end associate
end subroutine get_rvdw


!> Generator for the reference occupation numbers of the atoms
subroutine get_reference_occ(self, mol, bas, refocc)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Reference occupation numbers
   real(wp), intent(out) :: refocc(:, :)

   integer :: isp, ir, ish

   refocc(:, :) = 0.0_wp

   associate(record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         ir = irc(isp)
         do ish = 1, bas%nsh_id(isp)
            refocc(ish, isp) = merge(record(ir)%record(ish)%refocc, 0.0_wp, &
               & self%valence(ish, isp))
         end do
      end do
   end associate
end subroutine get_reference_occ


!> Generator for the diatomic frame scaling factors
subroutine get_diat_scale(self, mol, bas, ksig, kpi, kdel, do_diat_scale)
   !> Instance of the Hamiltonian specification
   class(param_h0spec), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Diatomic frame scaling of sigma bonding contribution
   real(wp), intent(out) :: ksig(:, :)
   !> Diatomic frame scaling of pi bonding contribution
   real(wp), intent(out) :: kpi(:, :)
   !> Diatomic frame scaling of delta bonding contribution
   real(wp), intent(out) :: kdel(:, :)
   !> Flag to indicate if diatomic frame scaling is used
   logical, intent(out) :: do_diat_scale

   integer :: isp, jsp, ir, jr

   do_diat_scale = self%param%hamiltonian%scaled_h0

   ksig(:, :) = 0.0_wp
   kpi(:, :) = 0.0_wp
   kdel(:, :) = 0.0_wp

   if (.not. do_diat_scale) return

   associate(record => self%param%record, irc => self%irc)
      do isp = 1, mol%nid
         ir = irc(isp)
         do jsp = 1, mol%nid
            jr = irc(jsp)
            ! mean in a GPT sense
            ksig(isp, jsp) = 2.0_wp/((1.0_wp / record(ir)%h0_diat_scale_sig) &
               & + (1.0_wp / record(jr)%h0_diat_scale_sig))
            kpi (isp, jsp) = 2.0_wp/((1.0_wp / record(ir)%h0_diat_scale_pi) &
               & + (1.0_wp / record(jr)%h0_diat_scale_pi))
            kdel(isp, jsp) = 2.0_wp/((1.0_wp / record(ir)%h0_diat_scale_del) &
               & + (1.0_wp / record(jr)%h0_diat_scale_del))
         end do
      end do
   end associate
end subroutine get_diat_scale


subroutine update(self, mol)
   class(xtb_calculator), intent(inout) :: self
   type(structure_type), intent(in) :: mol
end subroutine update


!> Add an interaction container
subroutine push_back(self, cont)
   !> Instance of the tight-binding calculator
   class(xtb_calculator), intent(inout) :: self
   !> Container to be added
   class(container_type), allocatable, intent(inout) :: cont

   if (.not.allocated(self%interactions)) allocate(self%interactions)

   call self%interactions%push_back(cont)
end subroutine push_back


!> Add a container
subroutine pop(self, cont, idx)
   !> Instance of the tight-binding calculator
   class(xtb_calculator), intent(inout) :: self
   !> Container to be removed
   class(container_type), allocatable, intent(out) :: cont
   !> Index to remove container from
   integer, intent(in), optional :: idx

   if (.not.allocated(self%interactions)) return

   call self%interactions%pop(cont, idx)
end subroutine pop


pure function variable_info(self) result(info)
   !> Instance of the electrostatic container
   class(xtb_calculator), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info()

   if (allocated(self%coulomb)) then
      info = max(info, self%coulomb%variable_info())
   end if

   if (allocated(self%dispersion)) then
      info = max(info, self%dispersion%variable_info())
   end if

   if (allocated(self%exchange)) then
      info = max(info, self%exchange%variable_info())
   end if

   if (allocated(self%interactions)) then
      info = max(info, self%interactions%variable_info())
   end if

end function variable_info


!> Information on container
pure function info(self, verbosity, indent) result(str)
   !> Instance of the interaction container
   class(xtb_calculator), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   character(len=*), parameter :: nl = new_line('a')

   str = "xTB calculator"

   if (allocated(self%repulsion)) then
      str = str // nl // indent // self%repulsion%info(verbosity, indent)
   end if

   if (allocated(self%coulomb)) then
      str = str // nl // indent // self%coulomb%info(verbosity, indent)
   end if

   if (allocated(self%dispersion)) then
      str = str // nl // indent // self%dispersion%info(verbosity, indent)
   end if

   if (allocated(self%exchange)) then
      str = str // nl // indent // self%exchange%info(verbosity, indent)
   end if

   if (allocated(self%interactions)) then
      str = str // nl // indent // self%interactions%info(verbosity, indent)
   end if
end function info

!> Find the base method the parametrization is based on
pure function get_method(method) result(str)
   !> Method name from parameterization
   character(len=*), intent(in) :: method
   !> Information on the container
   character(len=:), allocatable :: str

   if (method == "GFN2-xTB" .or. method == "gfn2") then
      str = "gfn2"
   else if (method == "GFN1-xTB" .or. method == "gfn1") then
      str = "gfn1"
   else if (method == "g-xTB" .or. method == "gxtb") then
      str = "gxtb"
   else
      str = "custom"
   end if

end function get_method


end module tblite_xtb_calculator
