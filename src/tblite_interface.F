!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief interface to tblite
!> \author JVP
!> \history creation 09.2024
! **************************************************************************************************

MODULE tblite_interface

#if defined(__TBLITE)
   USE mctc_env, ONLY: error_type
   USE mctc_io, ONLY: structure_type, new
   USE mctc_io_symbols, ONLY: symbol_to_number
   USE tblite_adjlist, ONLY: adjacency_list, new_adjacency_list
   USE tblite_basis_type, ONLY: get_cutoff
   USE tblite_container, ONLY: container_cache
   USE tblite_cutoff, ONLY: get_lattice_points
   USE tblite_integral_multipole, ONLY: multipole_cgto, multipole_grad_cgto, maxl, msao
   USE tblite_integral_type, ONLY: integral_type, new_integral
   USE tblite_scf, ONLY: get_mixer_dimension, new_mixer
   USE tblite_scf_info, ONLY: scf_info, atom_resolved, shell_resolved, &
                              orbital_resolved, not_used
   USE tblite_scf_potential, ONLY: potential_type, new_potential, add_pot_to_h1
   USE tblite_wavefunction_type, ONLY: wavefunction_type, new_wavefunction
   USE tblite_xtb_calculator, ONLY: xtb_calculator, new_xtb_calculator
   USE tblite_xtb_gfn1, ONLY: new_gfn1_calculator
   USE tblite_xtb_gfn2, ONLY: new_gfn2_calculator
   USE tblite_xtb_h0, ONLY: get_selfenergy, get_hamiltonian, get_occupation, &
                            get_hamiltonian_gradient, tb_hamiltonian
   USE tblite_xtb_ipea1, ONLY: new_ipea1_calculator
#endif
   USE ai_contraction, ONLY: block_add, &
                             contraction
   USE ai_overlap, ONLY: overlap_ab
   USE atomic_kind_types, ONLY: atomic_kind_type, get_atomic_kind, get_atomic_kind_set
   USE atprop_types, ONLY: atprop_type
   USE basis_set_types, ONLY: gto_basis_set_type, gto_basis_set_p_type, &
                              & allocate_gto_basis_set, write_gto_basis_set, process_gto_basis
   USE cell_types, ONLY: cell_type, get_cell
   USE cp_blacs_env, ONLY: cp_blacs_env_type
   USE cp_control_types, ONLY: dft_control_type, xtb_reference_cli_type
   USE cp_dbcsr_api, ONLY: dbcsr_type, dbcsr_p_type, dbcsr_create, dbcsr_add, &
                           dbcsr_get_block_p, dbcsr_finalize, &
                           dbcsr_iterator_type, dbcsr_iterator_blocks_left, &
                           dbcsr_iterator_start, dbcsr_iterator_stop, &
                           dbcsr_iterator_next_block
   USE cp_dbcsr_contrib, ONLY: dbcsr_print
   USE cp_dbcsr_cp2k_link, ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations, ONLY: dbcsr_allocate_matrix_set
   USE cp_log_handling, ONLY: cp_get_default_logger, &
                              cp_logger_type, cp_logger_get_default_io_unit
   USE cp_output_handling, ONLY: cp_print_key_should_output, &
                                 cp_print_key_unit_nr, cp_print_key_finished_output
   USE cp_units, ONLY: cp_unit_from_cp2k
   USE input_constants, ONLY: gfn1xtb, gfn2xtb, ipea1xtb, &
                              tblite_scc_mixer_auto, tblite_scc_mixer_cp2k, &
                              tblite_scc_mixer_none, tblite_scc_mixer_tblite
   USE input_section_types, ONLY: section_vals_val_get
   USE kinds, ONLY: default_path_length, default_string_length, dp
   USE kpoint_types, ONLY: get_kpoint_info, kpoint_type
   USE memory_utilities, ONLY: reallocate
   USE message_passing, ONLY: mp_para_env_type
   USE mulliken, ONLY: ao_charges
   USE orbital_pointers, ONLY: ncoset
   USE particle_types, ONLY: particle_type
   USE qs_charge_mixing, ONLY: charge_mixing
   USE qs_condnum, ONLY: overlap_condnum
   USE qs_energy_types, ONLY: qs_energy_type
   USE qs_environment_types, ONLY: get_qs_env, qs_environment_type
   USE qs_force_types, ONLY: qs_force_type, total_qs_force
   USE qs_integral_utils, ONLY: basis_set_list_setup, get_memory_usage
   USE qs_kind_types, ONLY: get_qs_kind, qs_kind_type, get_qs_kind_set
   USE qs_ks_types, ONLY: get_ks_env, qs_ks_env_type, set_ks_env
   USE qs_neighbor_list_types, ONLY: neighbor_list_iterator_create, neighbor_list_iterate, &
                                     get_iterator_info, neighbor_list_set_p_type, &
                                     neighbor_list_iterator_p_type, neighbor_list_iterator_release
   USE qs_overlap, ONLY: create_sab_matrix
   USE cp_dbcsr_output, ONLY: cp_dbcsr_write_sparse_matrix
   USE qs_rho_types, ONLY: qs_rho_get, qs_rho_type
   USE qs_scf_types, ONLY: qs_scf_env_type
   USE scf_control_types, ONLY: scf_control_type
   USE input_section_types, ONLY: section_vals_get_subs_vals, section_vals_type
   USE string_utilities, ONLY: integer_to_string
   USE tblite_types, ONLY: tblite_type, deallocate_tblite_type, allocate_tblite_type
   USE virial_types, ONLY: virial_type
   USE xtb_types, ONLY: get_xtb_atom_param, xtb_atom_type
   USE xtb_types, ONLY: xtb_atom_type

#include "./base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tblite_interface'

   INTEGER, PARAMETER                                 :: dip_n = 3
   INTEGER, PARAMETER                                 :: quad_n = 6
   REAL(KIND=dp), PARAMETER                           :: same_atom = 0.00001_dp
   REAL(KIND=dp), PARAMETER                           :: tblite_scc_pconv = 2.0E-5_dp

   PUBLIC :: tb_set_calculator, tb_init_geometry, tb_init_wf
   PUBLIC :: tb_get_basis, build_tblite_matrices
   PUBLIC :: tb_get_energy, tb_update_charges, tb_ham_add_coulomb
   PUBLIC :: tb_native_scc_mixer_active
   PUBLIC :: tb_scf_mixer_error
   PUBLIC :: tb_get_multipole
   PUBLIC :: tb_derive_dH_diag, tb_derive_dH_off
   PUBLIC :: tb_reference_cli_compare

CONTAINS

! **************************************************************************************************
!> \brief intialize geometry objects ...
!> \param qs_env ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_init_geometry(qs_env, tb)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb

#if defined(__TBLITE)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tblite_init_geometry'

      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER                                            :: iatom, natom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xyz
      INTEGER                                            :: handle, ikind
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL, DIMENSION(3)                              :: lperiod

      CALL timeset(routineN, handle)

      !get info from environment vaiarable
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell, qs_kind_set=qs_kind_set)

      !get information about particles
      natom = SIZE(particle_set)
      ALLOCATE (xyz(3, natom))
      CALL allocate_tblite_type(tb)
      ALLOCATE (tb%el_num(natom))
      tb%el_num = -9
      DO iatom = 1, natom
         xyz(:, iatom) = particle_set(iatom)%r(:)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), zatom=tb%el_num(iatom))
         IF (tb%el_num(iatom) < 1 .OR. tb%el_num(iatom) > 85) THEN
            CPABORT("only elements 1-85 are supported by tblite")
         END IF
      END DO

      !get information about cell / lattice
      CALL get_cell(cell=cell, periodic=periodic)
      lperiod(1) = periodic(1) == 1
      lperiod(2) = periodic(2) == 1
      lperiod(3) = periodic(3) == 1

      !prepare for the call to the dispersion function
      CALL new(tb%mol, tb%el_num, xyz, lattice=cell%hmat, periodic=lperiod)

      DEALLOCATE (xyz)

      CALL timestop(handle)

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_init_geometry

! **************************************************************************************************
!> \brief updating coordinates...
!> \param qs_env ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_update_geometry(qs_env, tb)

      TYPE(qs_environment_type)                          :: qs_env
      TYPE(tblite_type)                                  :: tb

#if defined(__TBLITE)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tblite_update_geometry'

      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER                                            :: iatom, natom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: xyz
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      !get info from environment vaiarable
      NULLIFY (cell)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell)

      !get information about particles
      natom = SIZE(particle_set)
      ALLOCATE (xyz(3, natom))
      DO iatom = 1, natom
         xyz(:, iatom) = particle_set(iatom)%r(:)
      END DO
      tb%mol%xyz(:, :) = xyz
      tb%mol%lattice(:, :) = cell%hmat

      DEALLOCATE (xyz)

      CALL timestop(handle)

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_update_geometry

! **************************************************************************************************
!> \brief initialize wavefunction ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_init_wf(tb)

      TYPE(tblite_type), POINTER                         :: tb

#if defined(__TBLITE)

      INTEGER, PARAMETER                                 :: nSpin = 1 !number of spin channels

      TYPE(scf_info)                                     :: info

      info = tb%calc%variable_info()
      IF (info%charge > shell_resolved) CPABORT("tblite: no support for orbital resolved charge")
      IF (info%dipole > atom_resolved) CPABORT("tblite: no support for shell resolved dipole moment")
      IF (info%quadrupole > atom_resolved) &
         CPABORT("tblite: no support shell resolved quadrupole moment")

      CALL new_wavefunction(tb%wfn, tb%mol%nat, tb%calc%bas%nsh, tb%calc%bas%nao, nSpin, 0.0_dp)
      CALL tb_reset_mixer(tb)

      CALL new_potential(tb%pot, tb%mol, tb%calc%bas, tb%wfn%nspin)

      !allocate quantities later required
      ALLOCATE (tb%e_hal(tb%mol%nat), tb%e_rep(tb%mol%nat), tb%e_disp(tb%mol%nat))
      ALLOCATE (tb%e_scd(tb%mol%nat), tb%e_es(tb%mol%nat))
      ALLOCATE (tb%selfenergy(tb%calc%bas%nsh))
      IF (ALLOCATED(tb%calc%ncoord)) ALLOCATE (tb%cn(tb%mol%nat))

#else
      MARK_USED(tb)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_init_wf

! **************************************************************************************************
!> \brief Reset tblite's internal SCC mixer for a new CP2K SCF cycle.
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_reset_mixer(tb)

      TYPE(tblite_type), POINTER                         :: tb

#if defined(__TBLITE)

      TYPE(scf_info)                                     :: info

      info = tb%calc%variable_info()
      IF (ALLOCATED(tb%mixer)) DEALLOCATE (tb%mixer)
      CALL new_mixer(tb%mixer, tb%calc%max_iter, &
                     tb%wfn%nspin*get_mixer_dimension(tb%mol, tb%calc%bas, info), &
                     tb%calc%mixer_damping)

#else
      MARK_USED(tb)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_reset_mixer

! **************************************************************************************************
!> \brief Configure tblite's internal SCC mixer from CP2K input.
!> \param tb ...
!> \param max_scf ...
!> \param damping ...
! **************************************************************************************************
   SUBROUTINE tb_configure_mixer(tb, max_scf, damping)

      TYPE(tblite_type), POINTER                         :: tb
      INTEGER, INTENT(IN)                                :: max_scf
      REAL(KIND=dp), INTENT(IN)                          :: damping

#if defined(__TBLITE)

      IF (max_scf < 1) CPABORT("tblite SCC mixer MAX_SCF must be positive")
      IF (damping <= 0.0_dp) CPABORT("tblite SCC mixer damping must be positive")

      tb%calc%max_iter = max_scf
      tb%calc%mixer_damping = damping

#else
      MARK_USED(tb)
      MARK_USED(max_scf)
      MARK_USED(damping)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_configure_mixer

! **************************************************************************************************
!> \brief Return whether the tblite native SCC mixer is active for this run.
!> \param dft_control ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_native_scc_mixer_active(dft_control) RESULT(use_native_mixer)

      TYPE(dft_control_type), POINTER                    :: dft_control
      LOGICAL                                            :: use_native_mixer

      use_native_mixer = .FALSE.
      IF (.NOT. ASSOCIATED(dft_control)) RETURN

      SELECT CASE (dft_control%qs_control%xtb_control%tblite_scc_mixer)
      CASE (tblite_scc_mixer_auto)
         use_native_mixer = .TRUE.
      CASE (tblite_scc_mixer_tblite)
         use_native_mixer = .TRUE.
      CASE (tblite_scc_mixer_cp2k, tblite_scc_mixer_none)
         use_native_mixer = .FALSE.
      CASE DEFAULT
         CPABORT("Unknown tblite SCC mixer")
      END SELECT

   END FUNCTION tb_native_scc_mixer_active

! **************************************************************************************************
!> \brief Return the native tblite SCC mixer residual for CP2K SCF convergence.
!> \param dft_control ...
!> \param tb ...
!> \param eps_scf ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_scf_mixer_error(dft_control, tb, eps_scf) RESULT(mixer_error)

      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(tblite_type), POINTER                         :: tb
      REAL(KIND=dp), INTENT(IN)                          :: eps_scf
      REAL(KIND=dp)                                      :: mixer_error

      REAL(KIND=dp)                                      :: raw_error

      mixer_error = 0.0_dp

#if defined(__TBLITE)
      IF (.NOT. ASSOCIATED(tb)) RETURN
      IF (.NOT. tb_native_scc_mixer_active(dft_control)) RETURN
      IF (ALLOCATED(tb%mixer)) THEN
         raw_error = REAL(tb%mixer%get_error(), KIND=dp)
         IF (eps_scf > 0.0_dp) THEN
            mixer_error = eps_scf*raw_error/tblite_scc_pconv
         ELSE
            mixer_error = raw_error
         END IF
      END IF
#else
      MARK_USED(dft_control)
      MARK_USED(tb)
      MARK_USED(eps_scf)
#endif

   END FUNCTION tb_scf_mixer_error

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param typ ...
! **************************************************************************************************
   SUBROUTINE tb_set_calculator(tb, typ)

      TYPE(tblite_type), POINTER                         :: tb
      INTEGER                                            :: typ

#if defined(__TBLITE)

      TYPE(error_type), ALLOCATABLE                      :: error

      SELECT CASE (typ)
      CASE default
         CPABORT("Unknown xtb type")
      CASE (gfn1xtb)
         CALL new_gfn1_calculator(tb%calc, tb%mol, error)
      CASE (gfn2xtb)
         CALL new_gfn2_calculator(tb%calc, tb%mol, error)
      CASE (ipea1xtb)
         CALL new_ipea1_calculator(tb%calc, tb%mol, error)
      END SELECT

#else
      MARK_USED(tb)
      MARK_USED(typ)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_set_calculator

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_init_ham(qs_env, tb, para_env)

      TYPE(qs_environment_type)                          :: qs_env
      TYPE(tblite_type)                                  :: tb
      TYPE(mp_para_env_type)                             :: para_env

#if defined(__TBLITE)

      TYPE(container_cache)                              :: hcache, rcache

      tb%e_hal = 0.0_dp
      tb%e_rep = 0.0_dp
      tb%e_disp = 0.0_dp
      IF (ALLOCATED(tb%grad)) THEN
         tb%grad = 0.0_dp
         CALL tb_zero_force(qs_env)
      END IF
      tb%sigma = 0.0_dp

      IF (ALLOCATED(tb%calc%halogen)) THEN
         CALL tb%calc%halogen%update(tb%mol, hcache)
         IF (ALLOCATED(tb%grad)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%halogen%get_engrad(tb%mol, hcache, tb%e_hal, &
            & tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_halogen", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 0)
         ELSE
            CALL tb%calc%halogen%get_engrad(tb%mol, hcache, tb%e_hal)
         END IF
      END IF

      IF (ALLOCATED(tb%calc%repulsion)) THEN
         CALL tb%calc%repulsion%update(tb%mol, rcache)
         IF (ALLOCATED(tb%grad)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%repulsion%get_engrad(tb%mol, rcache, tb%e_rep, &
            & tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_repulsion", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 1)
         ELSE
            CALL tb%calc%repulsion%get_engrad(tb%mol, rcache, tb%e_rep)
         END IF
      END IF

      IF (ALLOCATED(tb%calc%dispersion)) THEN
         CALL tb%calc%dispersion%update(tb%mol, tb%dcache)
         IF (ALLOCATED(tb%grad)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%dispersion%get_engrad(tb%mol, tb%dcache, tb%e_disp, &
            & tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_dispersion_static", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 2)
         ELSE
            CALL tb%calc%dispersion%get_engrad(tb%mol, tb%dcache, tb%e_disp)
         END IF
      END IF

      CALL new_potential(tb%pot, tb%mol, tb%calc%bas, tb%wfn%nspin)
      IF (ALLOCATED(tb%calc%coulomb)) THEN
         CALL tb%calc%coulomb%update(tb%mol, tb%cache)
      END IF

      IF (ALLOCATED(tb%grad)) THEN
         IF (ALLOCATED(tb%calc%ncoord)) THEN
            CALL tb%calc%ncoord%get_cn(tb%mol, tb%cn, tb%dcndr, tb%dcndL)
         END IF
         CALL get_selfenergy(tb%calc%h0, tb%mol%id, tb%calc%bas%ish_at, &
         &  tb%calc%bas%nsh_id, cn=tb%cn, selfenergy=tb%selfenergy, dsedcn=tb%dsedcn)
      ELSE
         IF (ALLOCATED(tb%calc%ncoord)) THEN
            CALL tb%calc%ncoord%get_cn(tb%mol, tb%cn)
         END IF
         CALL get_selfenergy(tb%calc%h0, tb%mol%id, tb%calc%bas%ish_at, &
         &  tb%calc%bas%nsh_id, cn=tb%cn, selfenergy=tb%selfenergy, dsedcn=tb%dsedcn)
      END IF

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      MARK_USED(para_env)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_init_ham

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
!> \param energy ...
! **************************************************************************************************
   SUBROUTINE tb_get_energy(qs_env, tb, energy)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb
      TYPE(qs_energy_type), POINTER                      :: energy

#if defined(__TBLITE)

      INTEGER                                            :: iounit
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: scf_section

      NULLIFY (scf_section, logger)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)
      scf_section => section_vals_get_subs_vals(qs_env%input, "DFT%SCF")

      energy%repulsive = SUM(tb%e_rep)
      energy%el_stat = SUM(tb%e_es)
      energy%dispersion = SUM(tb%e_disp)
      energy%dispersion_sc = SUM(tb%e_scd)
      energy%xtb_xb_inter = SUM(tb%e_hal)

      energy%total = energy%core + energy%repulsive + energy%el_stat + energy%dispersion &
                     + energy%dispersion_sc + energy%xtb_xb_inter + energy%kTS &
                     + energy%efield + energy%qmmm_el

      iounit = cp_print_key_unit_nr(logger, scf_section, "PRINT%DETAILED_ENERGY", &
                                    extension=".scfLog")
      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,(T9,A,T60,F20.10))") &
            "Repulsive pair potential energy:               ", energy%repulsive, &
            "Zeroth order Hamiltonian energy:               ", energy%core, &
            "Electrostatic energy:                          ", energy%el_stat, &
            "Self-consistent dispersion energy:             ", energy%dispersion_sc, &
            "Non-self consistent dispersion energy:         ", energy%dispersion
         WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
            "Correction for halogen bonding:                ", energy%xtb_xb_inter
         IF (ABS(energy%efield) > 1.e-12_dp) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "Electric field interaction energy:             ", energy%efield
         END IF
         IF (qs_env%qmmm) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "QM/MM Electrostatic energy:                    ", energy%qmmm_el
         END IF
      END IF
      CALL cp_print_key_finished_output(iounit, logger, scf_section, &
                                        "PRINT%DETAILED_ENERGY")

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      MARK_USED(energy)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_get_energy

! **************************************************************************************************
!> \brief ...
!> \param tb ...
!> \param gto_basis_set ...
!> \param element_symbol ...
!> \param param ...
!> \param occ ...
! **************************************************************************************************
   SUBROUTINE tb_get_basis(tb, gto_basis_set, element_symbol, param, occ)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(gto_basis_set_type), POINTER                  :: gto_basis_set
      CHARACTER(len=2), INTENT(IN)                       :: element_symbol
      TYPE(xtb_atom_type), POINTER                       :: param
      INTEGER, DIMENSION(5), INTENT(out)                 :: occ

#if defined(__TBLITE)

      REAL(KIND=dp)                                      :: docc
      CHARACTER(LEN=default_string_length)               :: sng
      INTEGER                                            :: ang, i_type, id_atom, ind_ao, ipgf, iSH, &
                                                            ishell, ityp, maxl, mprim, natorb, &
                                                            nset, nshell
      LOGICAL                                            :: do_ortho

      CALL allocate_gto_basis_set(gto_basis_set)

      !identifying element in the bas data
      CALL symbol_to_number(i_type, element_symbol)
      DO id_atom = 1, tb%mol%nat
         IF (i_type == tb%el_num(id_atom)) EXIT
      END DO
      param%z = i_type
      param%symbol = element_symbol
      param%defined = .TRUE.
      ityp = tb%mol%id(id_atom)

      !getting size information
      nset = tb%calc%bas%nsh_id(ityp)
      nshell = 1
      mprim = 0
      DO ishell = 1, nset
         mprim = MAX(mprim, tb%calc%bas%cgto(ishell, ityp)%nprim)
      END DO
      param%nshell = nset
      natorb = 0

      !write basis set information
      CALL integer_to_string(mprim, sng)
      gto_basis_set%name = element_symbol//"_STO-"//TRIM(sng)//"G"
      gto_basis_set%nset = nset
      CALL reallocate(gto_basis_set%lmax, 1, nset)
      CALL reallocate(gto_basis_set%lmin, 1, nset)
      CALL reallocate(gto_basis_set%npgf, 1, nset)
      CALL reallocate(gto_basis_set%nshell, 1, nset)
      CALL reallocate(gto_basis_set%n, 1, 1, 1, nset)
      CALL reallocate(gto_basis_set%l, 1, 1, 1, nset)
      CALL reallocate(gto_basis_set%zet, 1, mprim, 1, nset)
      CALL reallocate(gto_basis_set%gcc, 1, mprim, 1, 1, 1, nset)

      ind_ao = 0
      maxl = 0
      DO ishell = 1, nset
         ang = tb%calc%bas%cgto(ishell, ityp)%ang
         natorb = natorb + (2*ang + 1)
         param%lval(ishell) = ang
         maxl = MAX(ang, maxl)
         gto_basis_set%lmax(ishell) = ang
         gto_basis_set%lmin(ishell) = ang
         gto_basis_set%npgf(ishell) = tb%calc%bas%cgto(ishell, ityp)%nprim
         gto_basis_set%nshell(ishell) = nshell
         gto_basis_set%n(1, ishell) = ang + 1
         gto_basis_set%l(1, ishell) = ang
         DO ipgf = 1, gto_basis_set%npgf(ishell)
            gto_basis_set%gcc(ipgf, 1, ishell) = tb%calc%bas%cgto(ishell, ityp)%coeff(ipgf)
            gto_basis_set%zet(ipgf, ishell) = tb%calc%bas%cgto(ishell, ityp)%alpha(ipgf)
         END DO
         DO ipgf = 1, (2*ang + 1)
            ind_ao = ind_ao + 1
            param%lao(ind_ao) = ang
            param%nao(ind_ao) = ishell
         END DO
      END DO

      do_ortho = .FALSE.
      CALL process_gto_basis(gto_basis_set, do_ortho, nset, maxl)

      !setting additional values in parameter
      param%rcut = get_cutoff(tb%calc%bas)
      param%natorb = natorb
      param%lmax = maxl !max angular momentum

      !getting occupation
      occ = 0
      docc = 0.0_dp
      IF (tb%calc%bas%nsh_at(id_atom) > 5) CPABORT("too many shells in tblite")
      DO iSh = 1, tb%calc%bas%nsh_at(id_atom)
         occ(iSh) = NINT(tb%calc%h0%refocc(iSh, ityp) + docc)
         docc = docc + tb%calc%h0%refocc(iSh, ityp) - REAL(occ(iSh))
         param%occupation(iSh) = occ(iSh)
      END DO
      IF (ABS(docc) > 0.1_dp) CPABORT("Getting occupation numbers from tblite fails")
      param%zeff = SUM(occ) !effective core charge

      !set normalization process
      gto_basis_set%norm_type = 3

#else
      occ = 0
      MARK_USED(tb)
      MARK_USED(gto_basis_set)
      MARK_USED(element_symbol)
      MARK_USED(param)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_get_basis

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE build_tblite_matrices(qs_env, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

#if defined(__TBLITE)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_tblite_matrices'

      INTEGER                                            :: handle, maxder, nderivatives, nimg, img, nkind, i, j, &
                                                            ic, iw, iatom, jatom, ikind, jkind, iset, jset, n1, n2, icol, &
                                                            irow, ia, ib, sgfa, sgfb, atom_a, atom_b, ldsab, nseta, nsetb, &
                                                            natorb_a, natorb_b, sgfa0, mp_test_mode, i0, j0
      LOGICAL                                            :: found, norml1, norml2, use_arnoldi
      REAL(KIND=dp)                                      :: dr, rr
      INTEGER, DIMENSION(3)                              :: cell
      REAL(KIND=dp)                                      :: hij, shpoly
      REAL(KIND=dp), DIMENSION(2)                        :: condnum
      REAL(KIND=dp), DIMENSION(3)                        :: rij
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      INTEGER                                            :: debug_status
      CHARACTER(LEN=32)                                  :: debug_value
#endif
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind
      INTEGER, DIMENSION(:, :, :), POINTER                :: cell_to_index
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: owork, native_s, native_h, native_lattr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: oint, sint, hint, native_dip, native_quad
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min
      INTEGER, DIMENSION(:), POINTER                     :: npgfa, npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, zeta, zetb, scon_a, scon_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: sblock, fblock

      TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
      TYPE(adjacency_list)                                  :: native_list
      TYPE(atprop_type), POINTER                            :: atprop
      TYPE(cp_blacs_env_type), POINTER                      :: blacs_env
      TYPE(cp_logger_type), POINTER                         :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER          :: matrix_h, matrix_s, matrix_p, matrix_w
      TYPE(dft_control_type), POINTER                       :: dft_control
      TYPE(qs_force_type), DIMENSION(:), POINTER            :: force
      TYPE(gto_basis_set_type), POINTER                     :: basis_set_a, basis_set_b
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER     :: basis_set_list
      TYPE(kpoint_type), POINTER                            :: kpoints
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                              :: nl_iterator
      TYPE(mp_para_env_type), POINTER                       :: para_env
      TYPE(qs_energy_type), POINTER                         :: energy
      TYPE(qs_ks_env_type), POINTER                         :: ks_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER             :: qs_kind_set
      TYPE(qs_rho_type), POINTER                            :: rho
      TYPE(tblite_type), POINTER                            :: tb
      TYPE(tb_hamiltonian), POINTER                         :: h0
      TYPE(virial_type), POINTER                            :: virial

      CALL timeset(routineN, handle)

      mp_test_mode = 0
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_MULTIPOLE_TEST", debug_value, STATUS=debug_status)
      IF (debug_status == 0) READ (debug_value, *, IOSTAT=debug_status) mp_test_mode
#endif

      NULLIFY (ks_env, energy, atomic_kind_set, qs_kind_set)
      NULLIFY (matrix_h, matrix_s, atprop, dft_control)
      NULLIFY (sab_orb, rho, tb, kpoints, cell_to_index)

      CALL get_qs_env(qs_env=qs_env, &
                      ks_env=ks_env, para_env=para_env, &
                      energy=energy, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      matrix_h_kp=matrix_h, &
                      matrix_s_kp=matrix_s, &
                      atprop=atprop, &
                      dft_control=dft_control, &
                      sab_orb=sab_orb, &
                      rho=rho, tb_tblite=tb)
      h0 => tb%calc%h0

      !update geometry (required for debug / geometry optimization)
      CALL tb_update_geometry(qs_env, tb)

      nkind = SIZE(atomic_kind_set)
      nderivatives = 0
      IF (calculate_forces) THEN
         nderivatives = 1
         IF (ALLOCATED(tb%grad)) DEALLOCATE (tb%grad)
         ALLOCATE (tb%grad(3, tb%mol%nat))
         IF (ALLOCATED(tb%dsedcn)) DEALLOCATE (tb%dsedcn)
         ALLOCATE (tb%dsedcn(tb%calc%bas%nsh))
         IF (ALLOCATED(tb%calc%ncoord)) THEN
            IF (ALLOCATED(tb%dcndr)) DEALLOCATE (tb%dcndr)
            ALLOCATE (tb%dcndr(3, tb%mol%nat, tb%mol%nat))
            IF (ALLOCATED(tb%dcndL)) DEALLOCATE (tb%dcndL)
            ALLOCATE (tb%dcndL(3, 3, tb%mol%nat))
         END IF
      ELSE
         IF (ALLOCATED(tb%grad)) DEALLOCATE (tb%grad)
         IF (ALLOCATED(tb%dcndr)) DEALLOCATE (tb%dcndr)
         IF (ALLOCATED(tb%dcndL)) DEALLOCATE (tb%dcndL)
      END IF
      maxder = ncoset(nderivatives)
      nimg = dft_control%nimages
      IF (nimg > 1) THEN
         CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
      END IF

      !intialise hamiltonian
      CALL tb_init_ham(qs_env, tb, para_env)

      ! get density matrtix
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)

      ! set up matrices for force calculations
      IF (calculate_forces) THEN
         NULLIFY (force, matrix_w, virial)
         CALL get_qs_env(qs_env=qs_env, &
                         matrix_w_kp=matrix_w, &
                         virial=virial, force=force)

         IF (SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimg
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
               CALL dbcsr_add(matrix_w(1, img)%matrix, matrix_w(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            END DO
         END IF
         tb%use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      END IF

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind)

      ! set up basis set lists
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      ! allocate overlap matrix
      CALL dbcsr_allocate_matrix_set(matrix_s, maxder, nimg)
      CALL create_sab_matrix(ks_env, matrix_s, "OVERLAP MATRIX", basis_set_list, basis_set_list, &
                             sab_orb, .TRUE.)
      CALL set_ks_env(ks_env, matrix_s_kp=matrix_s)

      ! initialize H matrix
      CALL dbcsr_allocate_matrix_set(matrix_h, 1, nimg)
      DO img = 1, nimg
         ALLOCATE (matrix_h(1, img)%matrix)
         CALL dbcsr_create(matrix_h(1, img)%matrix, template=matrix_s(1, img)%matrix, &
                           name="HAMILTONIAN MATRIX")
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_h(1, img)%matrix, sab_orb)
      END DO
      CALL set_ks_env(ks_env, matrix_h_kp=matrix_h)

      IF (mp_test_mode == 81 .OR. mp_test_mode == 82 .OR. &
          mp_test_mode == 232 .OR. mp_test_mode == 233) THEN
         IF (nimg > 1) CPABORT("Native tblite matrix test with k-points not implemented")
         ALLOCATE (native_s(tb%calc%bas%nao, tb%calc%bas%nao))
         ALLOCATE (native_h(tb%calc%bas%nao, tb%calc%bas%nao))
         ALLOCATE (native_dip(dip_n, tb%calc%bas%nao, tb%calc%bas%nao))
         ALLOCATE (native_quad(quad_n, tb%calc%bas%nao, tb%calc%bas%nao))
         CALL get_lattice_points(tb%mol%periodic, tb%mol%lattice, get_cutoff(tb%calc%bas, 1.0_dp), native_lattr)
         CALL new_adjacency_list(native_list, tb%mol, native_lattr, get_cutoff(tb%calc%bas, 1.0_dp))
         CALL get_hamiltonian(tb%mol, native_lattr, native_list, tb%calc%bas, tb%calc%h0, tb%selfenergy, &
                              native_s, native_dip, native_quad, native_h)

         NULLIFY (nl_iterator)
         CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom)
            icol = MAX(iatom, jatom)
            irow = MIN(iatom, jatom)
            i0 = tb%calc%bas%iao_sh(tb%calc%bas%ish_at(irow) + 1)
            j0 = tb%calc%bas%iao_sh(tb%calc%bas%ish_at(icol) + 1)

            NULLIFY (sblock, fblock)
            CALL dbcsr_get_block_p(matrix=matrix_s(1, 1)%matrix, row=irow, col=icol, BLOCK=sblock, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=matrix_h(1, 1)%matrix, row=irow, col=icol, BLOCK=fblock, found=found)
            CPASSERT(found)

            DO j = 1, SIZE(sblock, 1)
               DO i = 1, SIZE(sblock, 2)
                  sblock(j, i) = native_s(i0 + j, j0 + i)
                  fblock(j, i) = native_h(i0 + j, j0 + i)
               END DO
            END DO
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

         DO img = 1, nimg
            DO i = 1, SIZE(matrix_s, 1)
               CALL dbcsr_finalize(matrix_s(i, img)%matrix)
            END DO
            DO i = 1, SIZE(matrix_h, 1)
               CALL dbcsr_finalize(matrix_h(i, img)%matrix)
            END DO
         END DO

         IF (dft_control%qs_control%xtb_control%tblite_method == gfn2xtb) &
            CALL tb_get_multipole(qs_env, tb)

         DEALLOCATE (native_s, native_h, native_dip, native_quad, native_lattr)
         DEALLOCATE (basis_set_list)
         CALL timestop(handle)
         RETURN
      END IF

      ! loop over all atom pairs with a non-zero overlap (sab_orb)
      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)

         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)
         IF (icol == jatom) THEN
            rij = -rij
            i = ikind
            ikind = jkind
            jkind = i
         END IF

         dr = NORM2(rij(:))

         IF (nimg == 1) THEN
            ic = 1
         ELSE
            ic = cell_to_index(cell(1), cell(2), cell(3))
            CPASSERT(ic > 0)
         END IF
         NULLIFY (sblock)
         CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=sblock, found=found)
         CPASSERT(found)
         NULLIFY (fblock)
         CALL dbcsr_get_block_p(matrix=matrix_h(1, ic)%matrix, &
                                row=irow, col=icol, BLOCK=fblock, found=found)
         CPASSERT(found)

         ! --------- Overlap
         !get basis information
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         atom_a = atom_of_kind(icol)
         atom_b = atom_of_kind(irow)
         ! basis a
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         la_min => basis_set_a%lmin
         npgfa => basis_set_a%npgf
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         rpgfa => basis_set_a%pgf_radius
         set_radius_a => basis_set_a%set_radius
         scon_a => basis_set_a%scon
         zeta => basis_set_a%zet
         ! basis b
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         lb_min => basis_set_b%lmin
         npgfb => basis_set_b%npgf
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set
         rpgfb => basis_set_b%pgf_radius
         set_radius_b => basis_set_b%set_radius
         scon_b => basis_set_b%scon
         zetb => basis_set_b%zet

         ldsab = get_memory_usage(qs_kind_set, "ORB", "ORB")
         ALLOCATE (oint(ldsab, ldsab, maxder), owork(ldsab, ldsab))
         natorb_a = 0
         DO iset = 1, nseta
            natorb_a = natorb_a + (2*basis_set_a%l(1, iset) + 1)
         END DO
         natorb_b = 0
         DO iset = 1, nsetb
            natorb_b = natorb_b + (2*basis_set_b%l(1, iset) + 1)
         END DO
         ALLOCATE (sint(natorb_a, natorb_b, maxder))
         sint = 0.0_dp
         ALLOCATE (hint(natorb_a, natorb_b, maxder))
         hint = 0.0_dp

         !----------------- overlap integrals
         DO iset = 1, nseta
            n1 = npgfa(iset)*(ncoset(la_max(iset)) - ncoset(la_min(iset) - 1))
            sgfa = first_sgfa(1, iset)
            DO jset = 1, nsetb
               IF (set_radius_a(iset) + set_radius_b(jset) < dr) CYCLE
               n2 = npgfb(jset)*(ncoset(lb_max(jset)) - ncoset(lb_min(jset) - 1))
               sgfb = first_sgfb(1, jset)
               IF (calculate_forces) THEN
                  CALL overlap_ab(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                  lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                                  rij, sab=oint(:, :, 1), dab=oint(:, :, 2:4))
               ELSE
                  CALL overlap_ab(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                  lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                                  rij, sab=oint(:, :, 1))
               END IF
               ! Contraction
               CALL contraction(oint(:, :, 1), owork, ca=scon_a(:, sgfa:), na=n1, ma=nsgfa(iset), &
                                cb=scon_b(:, sgfb:), nb=n2, mb=nsgfb(jset), fscale=1.0_dp, trans=.FALSE.)
               CALL block_add("IN", owork, nsgfa(iset), nsgfb(jset), sint(:, :, 1), sgfa, sgfb, trans=.FALSE.)
               IF (calculate_forces) THEN
                  DO i = 2, 4
                     CALL contraction(oint(:, :, i), owork, ca=scon_a(:, sgfa:), na=n1, ma=nsgfa(iset), &
                                      cb=scon_b(:, sgfb:), nb=n2, mb=nsgfb(jset), fscale=1.0_dp, trans=.FALSE.)
                     CALL block_add("IN", owork, nsgfa(iset), nsgfb(jset), sint(:, :, i), sgfa, sgfb, trans=.FALSE.)
                  END DO
               END IF
            END DO
         END DO

         ! update S matrix
         IF (icol <= irow) THEN
            sblock(:, :) = sblock(:, :) + sint(:, :, 1)
         ELSE
            sblock(:, :) = sblock(:, :) + TRANSPOSE(sint(:, :, 1))
         END IF

         ! --------- Hamiltonian
         IF (icol == irow .AND. dr < same_atom) THEN
            !get diagonal F matrix from selfenergy
            n1 = tb%calc%bas%ish_at(icol)
            DO iset = 1, nseta
               sgfa = first_sgfa(1, iset)
               hij = tb%selfenergy(n1 + iset)
               DO ia = sgfa, sgfa + nsgfa(iset) - 1
                  hint(ia, ia, 1) = hij
               END DO
            END DO
         ELSE
            !get off-diagonal F matrix
            rr = SQRT(dr/(h0%rad(jkind) + h0%rad(ikind)))
            n1 = tb%calc%bas%ish_at(icol)
            sgfa0 = 1
            DO iset = 1, nseta
               sgfa = first_sgfa(1, iset)
               sgfa0 = sgfa0 + tb%calc%bas%cgto(iset, tb%mol%id(icol))%nprim
               n2 = tb%calc%bas%ish_at(irow)
               DO jset = 1, nsetb
                  sgfb = first_sgfb(1, jset)
                  shpoly = (1.0_dp + h0%shpoly(iset, ikind)*rr) &
                           *(1.0_dp + h0%shpoly(jset, jkind)*rr)
                  hij = 0.5_dp*(tb%selfenergy(n1 + iset) + tb%selfenergy(n2 + jset)) &
                        *h0%hscale(iset, jset, ikind, jkind)*shpoly
                  DO ia = sgfa, sgfa + nsgfa(iset) - 1
                     DO ib = sgfb, sgfb + nsgfb(jset) - 1
                        hint(ia, ib, 1) = hij*sint(ia, ib, 1)
                     END DO
                  END DO
               END DO
            END DO
         END IF

         ! update F matrix
         IF (icol <= irow) THEN
            fblock(:, :) = fblock(:, :) + hint(:, :, 1)
         ELSE
            fblock(:, :) = fblock(:, :) + TRANSPOSE(hint(:, :, 1))
         END IF

         DEALLOCATE (oint, owork, sint, hint)

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO img = 1, nimg
         DO i = 1, SIZE(matrix_s, 1)
            CALL dbcsr_finalize(matrix_s(i, img)%matrix)
         END DO
         DO i = 1, SIZE(matrix_h, 1)
            CALL dbcsr_finalize(matrix_h(i, img)%matrix)
         END DO
      END DO

      !compute multipole moments for gfn2
      IF (dft_control%qs_control%xtb_control%tblite_method == gfn2xtb) &
         CALL tb_get_multipole(qs_env, tb)

      ! output overlap information
      NULLIFY (logger)
      logger => cp_get_default_logger()
      IF (.NOT. calculate_forces) THEN
         IF (cp_print_key_should_output(logger%iter_info, qs_env%input, &
                                        "DFT%PRINT%OVERLAP_CONDITION") /= 0) THEN
            iw = cp_print_key_unit_nr(logger, qs_env%input, "DFT%PRINT%OVERLAP_CONDITION", &
                                      extension=".Log")
            CALL section_vals_val_get(qs_env%input, "DFT%PRINT%OVERLAP_CONDITION%1-NORM", l_val=norml1)
            CALL section_vals_val_get(qs_env%input, "DFT%PRINT%OVERLAP_CONDITION%DIAGONALIZATION", l_val=norml2)
            CALL section_vals_val_get(qs_env%input, "DFT%PRINT%OVERLAP_CONDITION%ARNOLDI", l_val=use_arnoldi)
            CALL get_qs_env(qs_env=qs_env, blacs_env=blacs_env)
            CALL overlap_condnum(matrix_s, condnum, iw, norml1, norml2, use_arnoldi, blacs_env)
         END IF
      END IF

      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

#else
      MARK_USED(qs_env)
      MARK_USED(calculate_forces)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE build_tblite_matrices

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dft_control ...
!> \param tb ...
!> \param calculate_forces ...
!> \param use_rho ...
! **************************************************************************************************
   SUBROUTINE tb_update_charges(qs_env, dft_control, tb, calculate_forces, use_rho)

      TYPE(qs_environment_type), POINTER                                 :: qs_env
      TYPE(dft_control_type), POINTER                                    :: dft_control
      TYPE(tblite_type), POINTER                                         :: tb
      LOGICAL, INTENT(IN)                                                :: calculate_forces
      LOGICAL, INTENT(IN)                                                :: use_rho

#if defined(__TBLITE)

      INTEGER                                            :: iatom, ikind, is, ns, atom_a, ii, im, irow, icol, i, j
      INTEGER                                            :: iao, jao, ish, ispin, i0, j0, nspin
      INTEGER                                            :: nimg, nkind, nsgf, natorb, na, n_mix_cols, mix_offset
      INTEGER                                            :: n_atom, max_orb, max_shell
      INTEGER                                            :: mp_test_mode, raw_state_status, raw_state_unit, &
                                                            state_status, state_unit, last_mix_slot, ia_local, ns_mix
      LOGICAL                                            :: discard_mixed_output, do_combined_mixing, do_dipole, do_quadrupole, &
                                                            native_sign_mixing, phased_combined_mixing, skip_charge_mixing, &
                                                            reuse_native_input, skip_scf_dispersion, skip_scf_dispersion_energy, &
                                                            skip_scf_dispersion_gradient, skip_scf_dispersion_potential, &
                                                            use_native_mixer, use_no_mixer
      REAL(KIND=dp)                                      :: mix_alpha, norm, moment_alpha, new_charge, new_moment, pao
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      INTEGER                                            :: debug_status
      CHARACTER(LEN=32)                                  :: debug_value
#endif
      CHARACTER(LEN=default_path_length)                 :: raw_state_file, state_file
      INTEGER, DIMENSION(5)                              :: occ
      INTEGER, DIMENSION(25)                             :: lao
      INTEGER, DIMENSION(25)                             :: nao
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ch_atom, ch_shell, ch_ref, ch_orb, mix_vars, prev_vars
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aocg, ao_dip, ao_quad, native_lattr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: native_p

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(adjacency_list)                               :: native_list
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, matrix_p
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_matrix
      TYPE(dbcsr_type), POINTER                          :: s_matrix
      TYPE(error_type), ALLOCATABLE                      :: error
      TYPE(integral_type)                                :: native_ints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p_block
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      ! also compute multipoles needed by GFN2
      do_dipole = .FALSE.
      do_quadrupole = .FALSE.
      mp_test_mode = 0
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_MULTIPOLE_TEST", debug_value, STATUS=debug_status)
      IF (debug_status == 0) READ (debug_value, *, IOSTAT=debug_status) mp_test_mode
#endif
      skip_scf_dispersion = .FALSE.
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DEBUG_SKIP_SCF_DISPERSION", debug_value, STATUS=debug_status)
      IF (debug_status == 0) READ (debug_value, *, IOSTAT=debug_status) skip_scf_dispersion
#endif
      skip_scf_dispersion_energy = skip_scf_dispersion
      skip_scf_dispersion_gradient = skip_scf_dispersion
      skip_scf_dispersion_potential = skip_scf_dispersion
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DEBUG_SKIP_SCF_DISPERSION_ENERGY", debug_value, STATUS=debug_status)
      IF (debug_status == 0) READ (debug_value, *, IOSTAT=debug_status) skip_scf_dispersion_energy
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DEBUG_SKIP_SCF_DISPERSION_GRADIENT", debug_value, STATUS=debug_status)
      IF (debug_status == 0) READ (debug_value, *, IOSTAT=debug_status) skip_scf_dispersion_gradient
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DEBUG_SKIP_SCF_DISPERSION_POTENTIAL", debug_value, STATUS=debug_status)
      IF (debug_status == 0) READ (debug_value, *, IOSTAT=debug_status) skip_scf_dispersion_potential
#endif
      SELECT CASE (dft_control%qs_control%xtb_control%tblite_scc_mixer)
      CASE (tblite_scc_mixer_auto)
         use_native_mixer = tb_native_scc_mixer_active(dft_control)
         use_no_mixer = .FALSE.
      CASE (tblite_scc_mixer_tblite)
         use_native_mixer = .TRUE.
         use_no_mixer = .FALSE.
      CASE (tblite_scc_mixer_cp2k)
         use_native_mixer = .FALSE.
         use_no_mixer = .FALSE.
      CASE (tblite_scc_mixer_none)
         use_native_mixer = .FALSE.
         use_no_mixer = .TRUE.
      CASE DEFAULT
         CPABORT("Unknown tblite SCC mixer")
      END SELECT

      ! compute mulliken charges required for charge update
      NULLIFY (particle_set, qs_kind_set, atomic_kind_set, scf_control)
      CALL get_qs_env(qs_env=qs_env, scf_env=scf_env, particle_set=particle_set, qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, matrix_s_kp=matrix_s, rho=rho, para_env=para_env, &
                      scf_control=scf_control)
      IF (use_native_mixer .AND. use_rho .AND. (.NOT. calculate_forces) .AND. &
          scf_env%iter_count == 1) THEN
         CALL tb_configure_mixer(tb, scf_control%max_scf, &
                                 dft_control%qs_control%xtb_control%tblite_mixer_damping)
         CALL tb_reset_mixer(tb)
      END IF
      NULLIFY (matrix_p)
      IF (use_rho) THEN
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      ELSE
         matrix_p => scf_env%p_mix_new
      END IF
      IF (ASSOCIATED(tb%dipbra)) do_dipole = .TRUE.
      IF (ASSOCIATED(tb%quadbra)) do_quadrupole = .TRUE.
      reuse_native_input = .FALSE.
      IF (use_native_mixer .AND. scf_env%iter_count == 1) THEN
         reuse_native_input = ANY(ABS(tb%wfn%qsh) > 1.0E-14_dp)
         IF (do_dipole) reuse_native_input = reuse_native_input .OR. &
                                             ANY(ABS(tb%wfn%dpat) > 1.0E-14_dp)
         IF (do_quadrupole) reuse_native_input = reuse_native_input .OR. &
                                                 ANY(ABS(tb%wfn%qpat) > 1.0E-14_dp)
      END IF
      n_atom = SIZE(particle_set)
      nkind = SIZE(atomic_kind_set)
      nimg = dft_control%nimages
      CALL get_qs_kind_set(qs_kind_set, maxsgf=nsgf)
      ALLOCATE (aocg(nsgf, n_atom))
      aocg = 0.0_dp
      IF (do_dipole) ALLOCATE (ao_dip(n_atom, dip_n))
      IF (do_quadrupole) ALLOCATE (ao_quad(n_atom, quad_n))
      max_orb = 0
      max_shell = 0
      DO ikind = 1, nkind
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
         CALL get_xtb_atom_param(xtb_kind, natorb=natorb)
         max_orb = MAX(max_orb, natorb)
      END DO
      DO is = 1, n_atom
         max_shell = MAX(max_shell, tb%calc%bas%nsh_at(is))
      END DO
      ALLOCATE (ch_atom(n_atom, 1), ch_shell(n_atom, max_shell))
      ALLOCATE (ch_orb(max_orb, n_atom), ch_ref(max_orb, n_atom))
      ch_atom = 0.0_dp
      ch_shell = 0.0_dp
      ch_orb = 0.0_dp
      ch_ref = 0.0_dp
      IF (nimg > 1) THEN
         CALL ao_charges(matrix_p, matrix_s, aocg, para_env)
         IF (do_dipole) THEN
            DO im = 1, dip_n
               CALL contract_dens_kp(matrix_p, tb%dipbra, tb%dipket, im, dip_n, ao_dip(:, im), para_env)
            END DO
         END IF
         IF (do_quadrupole) THEN
            DO im = 1, quad_n
               CALL contract_dens_kp(matrix_p, tb%quadbra, tb%quadket, im, quad_n, ao_quad(:, im), para_env)
            END DO
         END IF
      ELSE
         NULLIFY (p_matrix, s_matrix)
         p_matrix => matrix_p(:, 1)
         s_matrix => matrix_s(1, 1)%matrix
         CALL ao_charges(p_matrix, s_matrix, aocg, para_env)
         IF (do_dipole) THEN
            DO im = 1, dip_n
               CALL contract_dens(p_matrix, tb%dipbra(im)%matrix, tb%dipket(im)%matrix, ao_dip(:, im), para_env)
            END DO
         END IF
         IF (do_quadrupole) THEN
            DO im = 1, quad_n
               CALL contract_dens(p_matrix, tb%quadbra(im)%matrix, tb%quadket(im)%matrix, ao_quad(:, im), para_env)
            END DO
         END IF
      END IF
      NULLIFY (xtb_kind)
      DO ikind = 1, nkind
         CALL get_atomic_kind(atomic_kind_set(ikind), natom=na)
         CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
         CALL get_xtb_atom_param(xtb_kind, natorb=natorb, lao=lao, nao=nao, occupation=occ)
         DO iatom = 1, na
            atom_a = atomic_kind_set(ikind)%atom_list(iatom)
            DO is = 1, natorb
               ns = lao(is) + 1
               norm = 2*lao(is) + 1
               ch_ref(is, atom_a) = tb%calc%h0%refocc(nao(is), ikind)/norm
               ch_orb(is, atom_a) = aocg(is, atom_a) - ch_ref(is, atom_a)
               ch_shell(atom_a, ns) = ch_orb(is, atom_a) + ch_shell(atom_a, ns)
            END DO
            ch_atom(atom_a, 1) = SUM(ch_orb(:, atom_a))
         END DO
      END DO
      DEALLOCATE (aocg)

      raw_state_status = 1
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_RAW_STATE_DUMP", raw_state_file, STATUS=raw_state_status)
#endif
      IF (raw_state_status == 0) THEN
         OPEN (NEWUNIT=raw_state_unit, FILE=TRIM(raw_state_file), STATUS="REPLACE", ACTION="WRITE")
         WRITE (raw_state_unit, *) "qat"
         DO iatom = 1, n_atom
            WRITE (raw_state_unit, "(I0,1X,ES24.16)") iatom, -ch_atom(iatom, 1)
         END DO
         WRITE (raw_state_unit, *) "qsh"
         DO iatom = 1, n_atom
            DO is = 1, tb%calc%bas%nsh_at(iatom)
               WRITE (raw_state_unit, "(I0,1X,ES24.16)") tb%calc%bas%ish_at(iatom) + is, -ch_shell(iatom, is)
            END DO
         END DO
         IF (do_dipole) THEN
            WRITE (raw_state_unit, *) "dpat"
            DO iatom = 1, n_atom
               WRITE (raw_state_unit, "(I0,3(1X,ES24.16))") iatom, -ao_dip(iatom, :)
            END DO
         END IF
         IF (do_quadrupole) THEN
            WRITE (raw_state_unit, *) "qpat"
            DO iatom = 1, n_atom
               WRITE (raw_state_unit, "(I0,6(1X,ES24.16))") iatom, -ao_quad(iatom, :)
            END DO
         END IF
         CLOSE (raw_state_unit)
      END IF

      IF (mp_test_mode == 150) THEN
         n_mix_cols = max_shell
         IF (do_dipole) n_mix_cols = n_mix_cols + dip_n
         IF (do_quadrupole) n_mix_cols = n_mix_cols + quad_n
         ALLOCATE (mix_vars(n_atom, n_mix_cols))
         mix_vars = 0.0_dp

         mix_vars(:, 1:max_shell) = -ch_shell(:, 1:max_shell)
         mix_offset = max_shell
         IF (do_dipole) THEN
            mix_vars(:, mix_offset + 1:mix_offset + dip_n) = -ao_dip(:, 1:dip_n)
            mix_offset = mix_offset + dip_n
         END IF
         IF (do_quadrupole) THEN
            mix_vars(:, mix_offset + 1:mix_offset + quad_n) = -ao_quad(:, 1:quad_n)
         END IF

         IF (use_rho .AND. (.NOT. calculate_forces) .AND. scf_env%mixing_store%ncall > 0) THEN
            mix_vars = 0.0_dp
            ns_mix = MIN(n_mix_cols, scf_env%mixing_store%max_shell)
            last_mix_slot = MOD(scf_env%mixing_store%ncall - 1, scf_env%mixing_store%nbuffer) + 1
            DO ia_local = 1, scf_env%mixing_store%nat_local
               iatom = scf_env%mixing_store%atlist(ia_local)
               mix_vars(iatom, 1:ns_mix) = scf_env%mixing_store%acharge(ia_local, 1:ns_mix, last_mix_slot)
            END DO
            CALL para_env%sum(mix_vars)
         ELSEIF ((.NOT. use_rho) .AND. (.NOT. calculate_forces)) THEN
            ALLOCATE (prev_vars(n_atom, n_mix_cols))
            prev_vars(:, :) = mix_vars
            IF (scf_env%mixing_store%ncall > 0) THEN
               prev_vars = 0.0_dp
               ns_mix = MIN(n_mix_cols, scf_env%mixing_store%max_shell)
               last_mix_slot = MOD(scf_env%mixing_store%ncall - 1, scf_env%mixing_store%nbuffer) + 1
               DO ia_local = 1, scf_env%mixing_store%nat_local
                  iatom = scf_env%mixing_store%atlist(ia_local)
                  prev_vars(iatom, 1:ns_mix) = scf_env%mixing_store%acharge(ia_local, 1:ns_mix, last_mix_slot)
               END DO
               CALL para_env%sum(prev_vars)
            END IF
            mix_alpha = scf_env%mixing_store%alpha
            prev_vars(:, :) = (1.0_dp - mix_alpha)*prev_vars + mix_alpha*mix_vars
            scf_env%mixing_store%ncall = scf_env%mixing_store%ncall + 1
            ns_mix = MIN(n_mix_cols, scf_env%mixing_store%max_shell)
            last_mix_slot = MOD(scf_env%mixing_store%ncall - 1, scf_env%mixing_store%nbuffer) + 1
            DO ia_local = 1, scf_env%mixing_store%nat_local
               iatom = scf_env%mixing_store%atlist(ia_local)
               scf_env%mixing_store%acharge(ia_local, 1:ns_mix, last_mix_slot) = prev_vars(iatom, 1:ns_mix)
            END DO
            DEALLOCATE (prev_vars)
         END IF

         DO iatom = 1, n_atom
            ii = tb%calc%bas%ish_at(iatom)
            DO is = 1, tb%calc%bas%nsh_at(iatom)
               tb%wfn%qsh(ii + is, 1) = mix_vars(iatom, is)
            END DO
            tb%wfn%qat(iatom, 1) = SUM(tb%wfn%qsh(ii + 1:ii + tb%calc%bas%nsh_at(iatom), 1))
         END DO
         mix_offset = max_shell
         IF (do_dipole) THEN
            DO iatom = 1, n_atom
               tb%wfn%dpat(:, iatom, 1) = mix_vars(iatom, mix_offset + 1:mix_offset + dip_n)
            END DO
            mix_offset = mix_offset + dip_n
            DEALLOCATE (ao_dip)
         END IF
         IF (do_quadrupole) THEN
            DO iatom = 1, n_atom
               tb%wfn%qpat(:, iatom, 1) = mix_vars(iatom, mix_offset + 1:mix_offset + quad_n)
            END DO
            DEALLOCATE (ao_quad)
         END IF
         DEALLOCATE (mix_vars)
      ELSEIF (use_native_mixer .OR. mp_test_mode == 120 .OR. mp_test_mode == 220 .OR. &
              mp_test_mode == 230 .OR. mp_test_mode == 231 .OR. &
              mp_test_mode == 232 .OR. mp_test_mode == 233) THEN
         IF (.NOT. ALLOCATED(tb%mixer)) CPABORT("tblite mixer not initialized")
         IF (use_rho .AND. (.NOT. calculate_forces) .AND. scf_env%iter_count > 1) THEN
            CALL tb%mixer%next(error)
            IF (ALLOCATED(error)) CPABORT("tblite native mixer failed")
            CALL tb%mixer%get(tb%wfn%qsh)
            tb%wfn%qat(:, 1) = 0.0_dp
            DO iatom = 1, n_atom
               ii = tb%calc%bas%ish_at(iatom)
               tb%wfn%qat(iatom, 1) = SUM(tb%wfn%qsh(ii + 1:ii + tb%calc%bas%nsh_at(iatom), 1))
            END DO
            IF (do_dipole) THEN
               CALL tb%mixer%get(tb%wfn%dpat)
               DEALLOCATE (ao_dip)
            END IF
            IF (do_quadrupole) THEN
               CALL tb%mixer%get(tb%wfn%qpat)
               DEALLOCATE (ao_quad)
            END IF
         ELSE
            IF ((use_native_mixer .OR. mp_test_mode == 230 .OR. mp_test_mode == 231 .OR. &
                 mp_test_mode == 232 .OR. mp_test_mode == 233) .AND. &
                use_rho .AND. (.NOT. calculate_forces)) THEN
               IF (.NOT. reuse_native_input) THEN
                  tb%wfn%qsh(:, :) = 0.0_dp
                  tb%wfn%qat(:, :) = 0.0_dp
                  IF (do_dipole) tb%wfn%dpat(:, :, :) = 0.0_dp
                  IF (do_quadrupole) tb%wfn%qpat(:, :, :) = 0.0_dp
               END IF
               IF (do_dipole) DEALLOCATE (ao_dip)
               IF (do_quadrupole) DEALLOCATE (ao_quad)
            ELSE
               IF ((.NOT. use_rho) .AND. (.NOT. calculate_forces)) THEN
                  CALL tb%mixer%set(tb%wfn%qsh)
                  IF (do_dipole) CALL tb%mixer%set(tb%wfn%dpat)
                  IF (do_quadrupole) CALL tb%mixer%set(tb%wfn%qpat)
               END IF
               DO iatom = 1, n_atom
                  ii = tb%calc%bas%ish_at(iatom)
                  DO is = 1, tb%calc%bas%nsh_at(iatom)
                     tb%wfn%qsh(ii + is, 1) = -ch_shell(iatom, is)
                  END DO
                  tb%wfn%qat(iatom, 1) = SUM(tb%wfn%qsh(ii + 1:ii + tb%calc%bas%nsh_at(iatom), 1))
               END DO
               IF (do_dipole) THEN
                  DO iatom = 1, n_atom
                     tb%wfn%dpat(:, iatom, 1) = -ao_dip(iatom, :)
                  END DO
                  DEALLOCATE (ao_dip)
               END IF
               IF (do_quadrupole) THEN
                  DO iatom = 1, n_atom
                     tb%wfn%qpat(:, iatom, 1) = -ao_quad(iatom, :)
                  END DO
                  DEALLOCATE (ao_quad)
               END IF
               IF ((.NOT. use_rho) .AND. (.NOT. calculate_forces)) THEN
                  CALL tb%mixer%diff(tb%wfn%qsh)
                  IF (do_dipole) CALL tb%mixer%diff(tb%wfn%dpat)
                  IF (do_quadrupole) CALL tb%mixer%diff(tb%wfn%qpat)
               END IF
            END IF
         END IF
      ELSE
         ! charge mixing
         IF (dft_control%qs_control%do_ls_scf) THEN
            !
         ELSE
            ! CP2K SCC mixing has to see all tblite SCC variables, not only shell charges.
            do_combined_mixing = do_dipole .OR. do_quadrupole .OR. &
                                 (mp_test_mode == 110 .OR. mp_test_mode == 111 .OR. &
                                  mp_test_mode == 112 .OR. mp_test_mode == 114 .OR. &
                                  mp_test_mode == 115 .OR. &
                                  mp_test_mode == 116 .OR. mp_test_mode == 117 .OR. &
                                  mp_test_mode == 118 .OR. mp_test_mode == 122 .OR. &
                                  mp_test_mode == 130 .OR. mp_test_mode == 131 .OR. &
                                  mp_test_mode == 140 .OR. mp_test_mode == 141 .OR. &
                                  mp_test_mode == 142 .OR. mp_test_mode == 183)
            native_sign_mixing = ((do_dipole .OR. do_quadrupole) .AND. mp_test_mode == 0) .OR. &
                                 (mp_test_mode == 116 .OR. mp_test_mode == 117 .OR. &
                                  mp_test_mode == 118)
            phased_combined_mixing = (mp_test_mode == 130 .OR. mp_test_mode == 131)
            discard_mixed_output = (mp_test_mode == 114 .AND. .NOT. use_rho) .OR. &
                                   (mp_test_mode == 115 .AND. (.NOT. use_rho .OR. calculate_forces)) .OR. &
                                   (mp_test_mode == 117 .AND. .NOT. use_rho) .OR. &
                                   (mp_test_mode == 118 .AND. (.NOT. use_rho .OR. calculate_forces)) .OR. &
                                   (mp_test_mode == 131 .AND. (.NOT. use_rho .OR. calculate_forces))
            skip_charge_mixing = use_no_mixer .OR. &
                                 (mp_test_mode == 111 .AND. use_rho) .OR. &
                                 (mp_test_mode == 112 .AND. .NOT. use_rho) .OR. &
                                 (phased_combined_mixing .AND. use_rho)
            IF (phased_combined_mixing .AND. use_rho .AND. scf_env%mixing_store%ncall > 0 .AND. &
                (mp_test_mode /= 131 .OR. .NOT. calculate_forces)) THEN
               n_mix_cols = max_shell
               IF (do_dipole) n_mix_cols = n_mix_cols + dip_n
               IF (do_quadrupole) n_mix_cols = n_mix_cols + quad_n
               ALLOCATE (mix_vars(n_atom, n_mix_cols))
               mix_vars = 0.0_dp
               ns_mix = MIN(n_mix_cols, scf_env%mixing_store%max_shell)
               last_mix_slot = MOD(scf_env%mixing_store%ncall - 1, scf_env%mixing_store%nbuffer) + 1
               DO ia_local = 1, scf_env%mixing_store%nat_local
                  iatom = scf_env%mixing_store%atlist(ia_local)
                  mix_vars(iatom, 1:ns_mix) = scf_env%mixing_store%acharge(ia_local, 1:ns_mix, last_mix_slot)
               END DO
               CALL para_env%sum(mix_vars)
               ch_shell(:, 1:max_shell) = mix_vars(:, 1:max_shell)
               mix_offset = max_shell
               IF (do_dipole) THEN
                  ao_dip(:, 1:dip_n) = mix_vars(:, mix_offset + 1:mix_offset + dip_n)
                  mix_offset = mix_offset + dip_n
               END IF
               IF (do_quadrupole) THEN
                  ao_quad(:, 1:quad_n) = mix_vars(:, mix_offset + 1:mix_offset + quad_n)
               END IF
               DEALLOCATE (mix_vars)
            END IF
            IF (skip_charge_mixing) THEN
               !
            ELSEIF (do_combined_mixing) THEN
               n_mix_cols = max_shell
               IF (do_dipole) n_mix_cols = n_mix_cols + dip_n
               IF (do_quadrupole) n_mix_cols = n_mix_cols + quad_n
               ALLOCATE (mix_vars(n_atom, n_mix_cols))
               IF (native_sign_mixing) THEN
                  mix_vars(:, 1:max_shell) = -ch_shell(:, 1:max_shell)
               ELSE
                  mix_vars(:, 1:max_shell) = ch_shell(:, 1:max_shell)
               END IF
               mix_offset = max_shell
               IF (do_dipole) THEN
                  IF (native_sign_mixing) THEN
                     mix_vars(:, mix_offset + 1:mix_offset + dip_n) = -ao_dip(:, 1:dip_n)
                  ELSE
                     mix_vars(:, mix_offset + 1:mix_offset + dip_n) = ao_dip(:, 1:dip_n)
                  END IF
                  mix_offset = mix_offset + dip_n
               END IF
               IF (do_quadrupole) THEN
                  IF (native_sign_mixing) THEN
                     mix_vars(:, mix_offset + 1:mix_offset + quad_n) = -ao_quad(:, 1:quad_n)
                  ELSE
                     mix_vars(:, mix_offset + 1:mix_offset + quad_n) = ao_quad(:, 1:quad_n)
                  END IF
               END IF
               CALL charge_mixing(scf_env%mixing_method, scf_env%mixing_store, &
                                  mix_vars, para_env, scf_env%iter_count)
               IF (.NOT. discard_mixed_output) THEN
                  IF (native_sign_mixing) THEN
                     ch_shell(:, 1:max_shell) = -mix_vars(:, 1:max_shell)
                  ELSE
                     ch_shell(:, 1:max_shell) = mix_vars(:, 1:max_shell)
                  END IF
                  mix_offset = max_shell
                  IF (do_dipole) THEN
                     IF (native_sign_mixing) THEN
                        ao_dip(:, 1:dip_n) = -mix_vars(:, mix_offset + 1:mix_offset + dip_n)
                     ELSE
                        ao_dip(:, 1:dip_n) = mix_vars(:, mix_offset + 1:mix_offset + dip_n)
                     END IF
                     mix_offset = mix_offset + dip_n
                  END IF
                  IF (do_quadrupole) THEN
                     IF (native_sign_mixing) THEN
                        ao_quad(:, 1:quad_n) = -mix_vars(:, mix_offset + 1:mix_offset + quad_n)
                     ELSE
                        ao_quad(:, 1:quad_n) = mix_vars(:, mix_offset + 1:mix_offset + quad_n)
                     END IF
                  END IF
               END IF
               DEALLOCATE (mix_vars)
            ELSE
               CALL charge_mixing(scf_env%mixing_method, scf_env%mixing_store, &
                                  ch_shell, para_env, scf_env%iter_count)
            END IF
         END IF

         !setting new wave function
         CALL tb%pot%reset
         tb%e_es = 0.0_dp
         tb%e_scd = 0.0_dp
         moment_alpha = scf_env%p_mix_alpha
         IF (mp_test_mode == 61 .OR. mp_test_mode == 91) moment_alpha = 0.05_dp
         IF (mp_test_mode == 62 .OR. mp_test_mode == 92) moment_alpha = 0.5_dp
         DO iatom = 1, n_atom
            ii = tb%calc%bas%ish_at(iatom)
            DO is = 1, tb%calc%bas%nsh_at(iatom)
               new_charge = -ch_shell(iatom, is)
               IF ((mp_test_mode == 90 .OR. mp_test_mode == 91 .OR. mp_test_mode == 92) .AND. &
                   scf_env%iter_count > 1) THEN
                  new_charge = moment_alpha*new_charge + (1.0_dp - moment_alpha)*tb%wfn%qsh(ii + is, 1)
               END IF
               tb%wfn%qsh(ii + is, 1) = new_charge
            END DO
            IF (native_sign_mixing .OR. mp_test_mode == 122) THEN
               tb%wfn%qat(iatom, 1) = SUM(tb%wfn%qsh(ii + 1:ii + tb%calc%bas%nsh_at(iatom), 1))
            ELSE
               new_charge = -ch_atom(iatom, 1)
               IF ((mp_test_mode == 90 .OR. mp_test_mode == 91 .OR. mp_test_mode == 92) .AND. &
                   scf_env%iter_count > 1) THEN
                  new_charge = moment_alpha*new_charge + (1.0_dp - moment_alpha)*tb%wfn%qat(iatom, 1)
               END IF
               tb%wfn%qat(iatom, 1) = new_charge
            END IF
         END DO

         IF (do_dipole) THEN
            DO iatom = 1, n_atom
               DO im = 1, dip_n
                  new_moment = -ao_dip(iatom, im)
                  IF ((mp_test_mode == 60 .OR. mp_test_mode == 61 .OR. mp_test_mode == 62 .OR. &
                       mp_test_mode == 90 .OR. mp_test_mode == 91 .OR. mp_test_mode == 92) .AND. &
                      scf_env%iter_count > 1) THEN
                     new_moment = moment_alpha*new_moment + (1.0_dp - moment_alpha)*tb%wfn%dpat(im, iatom, 1)
                  END IF
                  tb%wfn%dpat(im, iatom, 1) = new_moment
               END DO
            END DO
            DEALLOCATE (ao_dip)
         END IF
         IF (do_quadrupole) THEN
            DO iatom = 1, n_atom
               DO im = 1, quad_n
                  new_moment = -ao_quad(iatom, im)
                  IF ((mp_test_mode == 60 .OR. mp_test_mode == 61 .OR. mp_test_mode == 62 .OR. &
                       mp_test_mode == 90 .OR. mp_test_mode == 91 .OR. mp_test_mode == 92) .AND. &
                      scf_env%iter_count > 1) THEN
                     new_moment = moment_alpha*new_moment + (1.0_dp - moment_alpha)*tb%wfn%qpat(im, iatom, 1)
                  END IF
                  tb%wfn%qpat(im, iatom, 1) = new_moment
               END DO
            END DO
            DEALLOCATE (ao_quad)
         END IF
      END IF
      IF (mp_test_mode == 2 .OR. mp_test_mode == 7 .OR. mp_test_mode == 8 .OR. &
          mp_test_mode == 13 .OR. mp_test_mode == 14 .OR. mp_test_mode == 15 .OR. &
          mp_test_mode == 16 .OR. mp_test_mode == 17 .OR. mp_test_mode == 18 .OR. &
          mp_test_mode == 20 .OR. mp_test_mode == 21 .OR. &
          mp_test_mode == 32 .OR. mp_test_mode == 33 .OR. &
          mp_test_mode == 34 .OR. mp_test_mode == 35 .OR. &
          mp_test_mode == 36 .OR. mp_test_mode == 37 .OR. &
          mp_test_mode == 10 .OR. mp_test_mode == 11 .OR. mp_test_mode == 12) THEN
         tb%wfn%dpat(:, :, :) = -tb%wfn%dpat(:, :, :)
         tb%wfn%qpat(:, :, :) = -tb%wfn%qpat(:, :, :)
      ELSEIF (mp_test_mode == 4) THEN
         tb%wfn%dpat(:, :, :) = 0.0_dp
         tb%wfn%qpat(:, :, :) = 0.0_dp
      END IF

      IF (mp_test_mode == 100 .OR. mp_test_mode == 101 .OR. mp_test_mode == 102) THEN
         nspin = SIZE(p_matrix)
         CALL new_integral(native_ints, tb%calc%bas%nao)
         CALL get_lattice_points(tb%mol%periodic, tb%mol%lattice, get_cutoff(tb%calc%bas, 1.0_dp), native_lattr)
         CALL new_adjacency_list(native_list, tb%mol, native_lattr, get_cutoff(tb%calc%bas, 1.0_dp))
         CALL get_hamiltonian(tb%mol, native_lattr, native_list, tb%calc%bas, tb%calc%h0, tb%selfenergy, &
                              native_ints%overlap, native_ints%dipole, native_ints%quadrupole, &
                              native_ints%hamiltonian)
         ALLOCATE (native_p(tb%calc%bas%nao, tb%calc%bas%nao, nspin))
         native_p = 0.0_dp
         DO ispin = 1, nspin
            CALL dbcsr_iterator_start(iter, p_matrix(ispin)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               NULLIFY (p_block)
               CALL dbcsr_iterator_next_block(iter, irow, icol, p_block)
               i0 = tb%calc%bas%iao_sh(tb%calc%bas%ish_at(irow) + 1)
               j0 = tb%calc%bas%iao_sh(tb%calc%bas%ish_at(icol) + 1)
               DO j = 1, SIZE(p_block, 1)
                  DO i = 1, SIZE(p_block, 2)
                     IF (mp_test_mode == 102 .AND. irow /= icol) THEN
                        native_p(i0 + j, j0 + i, ispin) = 0.5_dp*p_block(j, i)
                        native_p(j0 + i, i0 + j, ispin) = 0.5_dp*p_block(j, i)
                     ELSE
                        native_p(i0 + j, j0 + i, ispin) = p_block(j, i)
                        IF (mp_test_mode == 100 .AND. irow /= icol) &
                           native_p(j0 + i, i0 + j, ispin) = p_block(j, i)
                     END IF
                  END DO
               END DO
            END DO
            CALL dbcsr_iterator_stop(iter)
         END DO

         tb%wfn%qsh(:, :) = 0.0_dp
         DO ispin = 1, nspin
            DO iao = 1, tb%calc%bas%nao
               pao = 0.0_dp
               DO jao = 1, tb%calc%bas%nao
                  pao = pao + native_p(jao, iao, ispin)*native_ints%overlap(jao, iao)
               END DO
               ish = tb%calc%bas%ao2sh(iao)
               tb%wfn%qsh(ish, ispin) = tb%wfn%qsh(ish, ispin) - pao
            END DO
         END DO
         tb%wfn%qsh(:, 1) = tb%wfn%qsh(:, 1) + tb%wfn%n0sh(:)

         tb%wfn%qat(:, :) = 0.0_dp
         DO ispin = 1, nspin
            DO ish = 1, tb%calc%bas%nsh
               tb%wfn%qat(tb%calc%bas%sh2at(ish), ispin) = &
                  tb%wfn%qat(tb%calc%bas%sh2at(ish), ispin) + tb%wfn%qsh(ish, ispin)
            END DO
         END DO

         IF (do_dipole) THEN
            tb%wfn%dpat(:, :, :) = 0.0_dp
            DO ispin = 1, nspin
               DO iao = 1, tb%calc%bas%nao
                  atom_a = tb%calc%bas%ao2at(iao)
                  DO jao = 1, tb%calc%bas%nao
                     DO im = 1, dip_n
                        tb%wfn%dpat(im, atom_a, ispin) = tb%wfn%dpat(im, atom_a, ispin) &
                                                         - native_p(jao, iao, ispin)*native_ints%dipole(im, jao, iao)
                     END DO
                  END DO
               END DO
            END DO
         END IF
         IF (do_quadrupole) THEN
            tb%wfn%qpat(:, :, :) = 0.0_dp
            DO ispin = 1, nspin
               DO iao = 1, tb%calc%bas%nao
                  atom_a = tb%calc%bas%ao2at(iao)
                  DO jao = 1, tb%calc%bas%nao
                     DO im = 1, quad_n
                        tb%wfn%qpat(im, atom_a, ispin) = tb%wfn%qpat(im, atom_a, ispin) &
                                                         - native_p(jao, iao, ispin)*native_ints%quadrupole(im, jao, iao)
                     END DO
                  END DO
               END DO
            END DO
         END IF
         DEALLOCATE (native_p, native_lattr)
      END IF

      CALL tb%pot%reset
      tb%e_es = 0.0_dp
      tb%e_scd = 0.0_dp
      IF (ALLOCATED(tb%calc%coulomb)) THEN
         CALL tb%calc%coulomb%get_potential(tb%mol, tb%cache, tb%wfn, tb%pot)
         CALL tb%calc%coulomb%get_energy(tb%mol, tb%cache, tb%wfn, tb%e_es)
      END IF
      IF (ALLOCATED(tb%calc%dispersion)) THEN
         IF (.NOT. skip_scf_dispersion_potential) &
            CALL tb%calc%dispersion%get_potential(tb%mol, tb%dcache, tb%wfn, tb%pot)
         IF (.NOT. skip_scf_dispersion_energy) &
            CALL tb%calc%dispersion%get_energy(tb%mol, tb%dcache, tb%wfn, tb%e_scd)
      END IF

      state_status = 1
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_STATE_DUMP", state_file, STATUS=state_status)
#endif
      IF (state_status == 0) THEN
         OPEN (NEWUNIT=state_unit, FILE=TRIM(state_file), STATUS="REPLACE", ACTION="WRITE")
         WRITE (state_unit, *) "qat"
         DO iatom = 1, n_atom
            WRITE (state_unit, "(I0,1X,ES24.16)") iatom, tb%wfn%qat(iatom, 1)
         END DO
         WRITE (state_unit, *) "qsh"
         DO iatom = 1, SIZE(tb%wfn%qsh, 1)
            WRITE (state_unit, "(I0,1X,ES24.16)") iatom, tb%wfn%qsh(iatom, 1)
         END DO
         IF (do_dipole) THEN
            WRITE (state_unit, *) "dpat"
            DO iatom = 1, n_atom
               WRITE (state_unit, "(I0,3(1X,ES24.16))") iatom, tb%wfn%dpat(:, iatom, 1)
            END DO
         END IF
         IF (do_quadrupole) THEN
            WRITE (state_unit, *) "qpat"
            DO iatom = 1, n_atom
               WRITE (state_unit, "(I0,6(1X,ES24.16))") iatom, tb%wfn%qpat(:, iatom, 1)
            END DO
         END IF
         CLOSE (state_unit)
      END IF

      IF (calculate_forces) THEN
         IF (ALLOCATED(tb%calc%coulomb)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%coulomb%get_gradient(tb%mol, tb%cache, tb%wfn, tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_coulomb", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 3)
         END IF

         IF (ALLOCATED(tb%calc%dispersion) .AND. .NOT. skip_scf_dispersion_gradient) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%dispersion%get_gradient(tb%mol, tb%dcache, tb%wfn, tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_dispersion_scf", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 2)
         END IF
      END IF

      DEALLOCATE (ch_atom, ch_shell, ch_orb, ch_ref)

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      MARK_USED(dft_control)
      MARK_USED(calculate_forces)
      MARK_USED(use_rho)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_update_charges

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
!> \param dft_control ...
! **************************************************************************************************
   SUBROUTINE tb_ham_add_coulomb(qs_env, tb, dft_control)

      TYPE(qs_environment_type), POINTER                       :: qs_env
      TYPE(tblite_type), POINTER                               :: tb
      TYPE(dft_control_type), POINTER                          :: dft_control

#if defined(__TBLITE)

      INTEGER                                            :: ikind, jkind, iatom, jatom, icol, irow
      INTEGER                                            :: ic, id1, id2, id3, iq1, iq2, iq3, iq4, iq5, iq6, &
                                                            is, nimg, ni, nj, i, j, i0, j0
      INTEGER                                            :: la, lb, za, zb
      INTEGER                                            :: mp_test_mode
      LOGICAL                                            :: found
      INTEGER, DIMENSION(3)                              :: cellind
      INTEGER, DIMENSION(25)                             :: naoa, naob
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      REAL(KIND=dp)                                      :: mpfac
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      INTEGER                                            :: debug_status
      CHARACTER(LEN=32)                                  :: debug_value
#endif
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of, sum_shell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ashift, bshift, native_lattr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: native_h1
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: ksblock, sblock
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dip_ket1, dip_ket2, dip_ket3, &
                                                            dip_bra1, dip_bra2, dip_bra3
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: quad_ket1, quad_ket2, quad_ket3, &
                                                            quad_ket4, quad_ket5, quad_ket6, &
                                                            quad_bra1, quad_bra2, quad_bra3, &
                                                            quad_bra4, quad_bra5, quad_bra6

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(adjacency_list)                               :: native_list
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(integral_type)                                :: native_ints
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(potential_type)                               :: native_pot
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b

      nimg = dft_control%nimages
      mp_test_mode = 0
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_MULTIPOLE_TEST", debug_value, STATUS=debug_status)
      IF (debug_status == 0) READ (debug_value, *, IOSTAT=debug_status) mp_test_mode
#endif
      mpfac = -0.5_dp
      IF (mp_test_mode == 3 .OR. mp_test_mode == 10 .OR. mp_test_mode == 11 .OR. &
          mp_test_mode == 12) mpfac = 0.5_dp
      IF (mp_test_mode == 5) mpfac = 0.0_dp

      NULLIFY (matrix_s, ks_matrix, n_list, qs_kind_set)
      CALL get_qs_env(qs_env=qs_env, sab_orb=n_list, matrix_s_kp=matrix_s, matrix_ks_kp=ks_matrix, qs_kind_set=qs_kind_set)

      IF (mp_test_mode == 83 .OR. mp_test_mode == 183 .OR. &
          mp_test_mode == 220 .OR. mp_test_mode == 231 .OR. &
          mp_test_mode == 233) THEN
         IF (nimg > 1) CPABORT("Native tblite Hamiltonian test with k-points not implemented")
         CALL new_integral(native_ints, tb%calc%bas%nao)
         CALL get_lattice_points(tb%mol%periodic, tb%mol%lattice, get_cutoff(tb%calc%bas, 1.0_dp), native_lattr)
         CALL new_adjacency_list(native_list, tb%mol, native_lattr, get_cutoff(tb%calc%bas, 1.0_dp))
         CALL get_hamiltonian(tb%mol, native_lattr, native_list, tb%calc%bas, tb%calc%h0, tb%selfenergy, &
                              native_ints%overlap, native_ints%dipole, native_ints%quadrupole, &
                              native_ints%hamiltonian)
         ALLOCATE (native_h1(tb%calc%bas%nao, tb%calc%bas%nao, 1))
         native_pot = tb%pot
         CALL add_pot_to_h1(tb%calc%bas, native_ints, native_pot, native_h1)

         CALL dbcsr_iterator_start(iter, ks_matrix(1, 1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, irow, icol, ksblock)
            i0 = tb%calc%bas%iao_sh(tb%calc%bas%ish_at(irow) + 1)
            j0 = tb%calc%bas%iao_sh(tb%calc%bas%ish_at(icol) + 1)
            DO j = 1, SIZE(ksblock, 1)
               DO i = 1, SIZE(ksblock, 2)
                  ksblock(j, i) = native_h1(i0 + j, j0 + i, 1)
               END DO
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)
         DEALLOCATE (native_h1, native_lattr)
         RETURN
      END IF

      !creating sum of shell lists
      ALLOCATE (sum_shell(tb%mol%nat))
      i = 0
      DO j = 1, tb%mol%nat
         sum_shell(j) = i
         i = i + tb%calc%bas%nsh_at(j)
      END DO

      IF (nimg == 1) THEN
         ! no k-points; all matrices have been transformed to periodic bsf
         CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  kind_of=kind_of)
         CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, irow, icol, sblock)

            ikind = kind_of(irow)
            jkind = kind_of(icol)

            ! atomic parameters
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_a, z=za, nao=naoa)
            CALL get_xtb_atom_param(xtb_atom_b, z=zb, nao=naob)

            ni = SIZE(sblock, 1)
            ALLOCATE (ashift(ni, ni))
            ashift = 0.0_dp
            DO i = 1, ni
               la = naoa(i) + sum_shell(irow)
               ashift(i, i) = tb%pot%vsh(la, 1)
            END DO

            nj = SIZE(sblock, 2)
            ALLOCATE (bshift(nj, nj))
            bshift = 0.0_dp
            DO j = 1, nj
               lb = naob(j) + sum_shell(icol)
               bshift(j, j) = tb%pot%vsh(lb, 1)
            END DO

            DO is = 1, SIZE(ks_matrix, 1)
               NULLIFY (ksblock)
               CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                      row=irow, col=icol, block=ksblock, found=found)
               CPASSERT(found)
               ksblock = ksblock - 0.5_dp*(MATMUL(ashift, sblock) &
                                           + MATMUL(sblock, bshift))
               ksblock = ksblock - 0.5_dp*(tb%pot%vat(irow, 1) &
                                           + tb%pot%vat(icol, 1))*sblock
            END DO
            DEALLOCATE (ashift, bshift)
         END DO
         CALL dbcsr_iterator_stop(iter)

         IF (ASSOCIATED(tb%dipbra)) THEN
            CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, sblock)

               NULLIFY (dip_bra1, dip_bra2, dip_bra3)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(1)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(2)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(3)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra3, found=found)
               CPASSERT(found)
               NULLIFY (dip_ket1, dip_ket2, dip_ket3)
               CALL dbcsr_get_block_p(matrix=tb%dipket(1)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipket(2)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipket(3)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket3, found=found)
               CPASSERT(found)

               DO is = 1, SIZE(ks_matrix, 1)
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock + mpfac*(dip_ket1*tb%pot%vdp(1, irow, 1) &
                                             + dip_ket2*tb%pot%vdp(2, irow, 1) &
                                             + dip_ket3*tb%pot%vdp(3, irow, 1) &
                                             + dip_bra1*tb%pot%vdp(1, icol, 1) &
                                             + dip_bra2*tb%pot%vdp(2, icol, 1) &
                                             + dip_bra3*tb%pot%vdp(3, icol, 1))
               END DO
            END DO
            CALL dbcsr_iterator_stop(iter)
         END IF

         IF (ASSOCIATED(tb%quadbra)) THEN
            CALL dbcsr_iterator_start(iter, matrix_s(1, 1)%matrix)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, irow, icol, sblock)

               NULLIFY (quad_bra1, quad_bra2, quad_bra3, quad_bra4, quad_bra5, quad_bra6)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(1)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(2)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(3)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra3, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(4)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra4, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(5)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra5, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(6)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra6, found=found)
               CPASSERT(found)

               NULLIFY (quad_ket1, quad_ket2, quad_ket3, quad_ket4, quad_ket5, quad_ket6)
               CALL dbcsr_get_block_p(matrix=tb%quadket(1)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(2)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(3)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket3, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(4)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket4, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(5)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket5, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(6)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket6, found=found)
               CPASSERT(found)

               DO is = 1, SIZE(ks_matrix, 1)
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)

                  ksblock = ksblock + mpfac*(quad_ket1*tb%pot%vqp(1, irow, 1) &
                                             + quad_ket2*tb%pot%vqp(2, irow, 1) &
                                             + quad_ket3*tb%pot%vqp(3, irow, 1) &
                                             + quad_ket4*tb%pot%vqp(4, irow, 1) &
                                             + quad_ket5*tb%pot%vqp(5, irow, 1) &
                                             + quad_ket6*tb%pot%vqp(6, irow, 1) &
                                             + quad_bra1*tb%pot%vqp(1, icol, 1) &
                                             + quad_bra2*tb%pot%vqp(2, icol, 1) &
                                             + quad_bra3*tb%pot%vqp(3, icol, 1) &
                                             + quad_bra4*tb%pot%vqp(4, icol, 1) &
                                             + quad_bra5*tb%pot%vqp(5, icol, 1) &
                                             + quad_bra6*tb%pot%vqp(6, icol, 1))
               END DO
            END DO
            CALL dbcsr_iterator_stop(iter)
         END IF

      ELSE
         NULLIFY (kpoints)
         CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
         NULLIFY (cell_to_index)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)

         NULLIFY (nl_iterator)
         CALL neighbor_list_iterator_create(nl_iterator, n_list)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                   iatom=iatom, jatom=jatom, r=rij, cell=cellind)

            icol = MAX(iatom, jatom)
            irow = MIN(iatom, jatom)

            IF (icol == jatom) THEN
               rij = -rij
               i = ikind
               ikind = jkind
               jkind = i
            END IF

            ic = cell_to_index(cellind(1), cellind(2), cellind(3))
            CPASSERT(ic > 0)

            NULLIFY (sblock)
            CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                   row=irow, col=icol, block=sblock, found=found)
            CPASSERT(found)

            ! atomic parameters
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_atom_a)
            CALL get_qs_kind(qs_kind_set(jkind), xtb_parameter=xtb_atom_b)
            CALL get_xtb_atom_param(xtb_atom_a, z=za, nao=naoa)
            CALL get_xtb_atom_param(xtb_atom_b, z=zb, nao=naob)

            ni = SIZE(sblock, 1)
            ALLOCATE (ashift(ni, ni))
            ashift = 0.0_dp
            DO i = 1, ni
               la = naoa(i) + sum_shell(irow)
               ashift(i, i) = tb%pot%vsh(la, 1)
            END DO

            nj = SIZE(sblock, 2)
            ALLOCATE (bshift(nj, nj))
            bshift = 0.0_dp
            DO j = 1, nj
               lb = naob(j) + sum_shell(icol)
               bshift(j, j) = tb%pot%vsh(lb, 1)
            END DO

            DO is = 1, SIZE(ks_matrix, 1)
               NULLIFY (ksblock)
               CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                      row=irow, col=icol, block=ksblock, found=found)
               CPASSERT(found)
               ksblock = ksblock - 0.5_dp*(MATMUL(ashift, sblock) &
                                           + MATMUL(sblock, bshift))
               ksblock = ksblock - 0.5_dp*(tb%pot%vat(irow, 1) &
                                           + tb%pot%vat(icol, 1))*sblock
            END DO
            DEALLOCATE (ashift, bshift)
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)

         IF (ASSOCIATED(tb%dipbra)) THEN
            NULLIFY (nl_iterator)
            CALL neighbor_list_iterator_create(nl_iterator, n_list)
            DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
               CALL get_iterator_info(nl_iterator, &
                                      iatom=iatom, jatom=jatom, cell=cellind)
               icol = MAX(iatom, jatom)
               irow = MIN(iatom, jatom)
               ic = cell_to_index(cellind(1), cellind(2), cellind(3))
               CPASSERT(ic > 0)
               id1 = 1 + dip_n*(ic - 1)
               id2 = 2 + dip_n*(ic - 1)
               id3 = 3 + dip_n*(ic - 1)

               NULLIFY (dip_bra1, dip_bra2, dip_bra3)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(id1)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(id2)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipbra(id3)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_bra3, found=found)
               CPASSERT(found)
               NULLIFY (dip_ket1, dip_ket2, dip_ket3)
               CALL dbcsr_get_block_p(matrix=tb%dipket(id1)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipket(id2)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%dipket(id3)%matrix, &
                                      row=irow, col=icol, BLOCK=dip_ket3, found=found)
               CPASSERT(found)

               DO is = 1, SIZE(ks_matrix, 1)
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock + mpfac*(dip_ket1*tb%pot%vdp(1, irow, 1) &
                                             + dip_ket2*tb%pot%vdp(2, irow, 1) &
                                             + dip_ket3*tb%pot%vdp(3, irow, 1) &
                                             + dip_bra1*tb%pot%vdp(1, icol, 1) &
                                             + dip_bra2*tb%pot%vdp(2, icol, 1) &
                                             + dip_bra3*tb%pot%vdp(3, icol, 1))
               END DO
            END DO
            CALL neighbor_list_iterator_release(nl_iterator)
         END IF

         IF (ASSOCIATED(tb%quadbra)) THEN
            NULLIFY (nl_iterator)
            CALL neighbor_list_iterator_create(nl_iterator, n_list)
            DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
               CALL get_iterator_info(nl_iterator, &
                                      iatom=iatom, jatom=jatom, cell=cellind)
               icol = MAX(iatom, jatom)
               irow = MIN(iatom, jatom)
               ic = cell_to_index(cellind(1), cellind(2), cellind(3))
               CPASSERT(ic > 0)
               iq1 = 1 + quad_n*(ic - 1)
               iq2 = 2 + quad_n*(ic - 1)
               iq3 = 3 + quad_n*(ic - 1)
               iq4 = 4 + quad_n*(ic - 1)
               iq5 = 5 + quad_n*(ic - 1)
               iq6 = 6 + quad_n*(ic - 1)

               NULLIFY (quad_bra1, quad_bra2, quad_bra3, quad_bra4, quad_bra5, quad_bra6)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(iq1)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(iq2)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(iq3)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra3, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(iq4)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra4, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(iq5)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra5, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadbra(iq6)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_bra6, found=found)
               CPASSERT(found)

               NULLIFY (quad_ket1, quad_ket2, quad_ket3, quad_ket4, quad_ket5, quad_ket6)
               CALL dbcsr_get_block_p(matrix=tb%quadket(iq1)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket1, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(iq2)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket2, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(iq3)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket3, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(iq4)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket4, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(iq5)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket5, found=found)
               CPASSERT(found)
               CALL dbcsr_get_block_p(matrix=tb%quadket(iq6)%matrix, &
                                      row=irow, col=icol, BLOCK=quad_ket6, found=found)
               CPASSERT(found)

               DO is = 1, SIZE(ks_matrix, 1)
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)

                  ksblock = ksblock + mpfac*(quad_ket1*tb%pot%vqp(1, irow, 1) &
                                             + quad_ket2*tb%pot%vqp(2, irow, 1) &
                                             + quad_ket3*tb%pot%vqp(3, irow, 1) &
                                             + quad_ket4*tb%pot%vqp(4, irow, 1) &
                                             + quad_ket5*tb%pot%vqp(5, irow, 1) &
                                             + quad_ket6*tb%pot%vqp(6, irow, 1) &
                                             + quad_bra1*tb%pot%vqp(1, icol, 1) &
                                             + quad_bra2*tb%pot%vqp(2, icol, 1) &
                                             + quad_bra3*tb%pot%vqp(3, icol, 1) &
                                             + quad_bra4*tb%pot%vqp(4, icol, 1) &
                                             + quad_bra5*tb%pot%vqp(5, icol, 1) &
                                             + quad_bra6*tb%pot%vqp(6, icol, 1))
               END DO
            END DO
            CALL neighbor_list_iterator_release(nl_iterator)
         END IF

      END IF

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      MARK_USED(dft_control)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_ham_add_coulomb

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_get_multipole(qs_env, tb)

      TYPE(qs_environment_type), POINTER                       :: qs_env
      TYPE(tblite_type), POINTER                               :: tb

#if defined(__TBLITE)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tb_get_multipole'

      INTEGER                                     :: ikind, jkind, iatom, jatom, icol, irow, iset, jset, ityp, jtyp
      INTEGER                                     :: ic, idx, id1, id2, id3, img, iq1, iq2, iq3, iq4, iq5, iq6
      INTEGER                                     :: nkind, natom, handle, nimg, i, j, inda, indb
      INTEGER                                     :: mp_test_mode
      INTEGER                                     :: atom_a, atom_b, nseta, nsetb, ia, ib, ij, i0, j0
      LOGICAL                                     :: found
      REAL(KIND=dp)                               :: r2
      INTEGER, DIMENSION(3)                       :: cell
      INTEGER, DIMENSION(:, :, :), POINTER        :: cell_to_index
      REAL(KIND=dp), DIMENSION(3)                 :: rij, mp_int_vec, mp_shift_vec
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      INTEGER                                     :: debug_status
      CHARACTER(LEN=32)                           :: debug_value
#endif
      INTEGER, DIMENSION(:), POINTER              :: la_max, lb_max
      INTEGER, DIMENSION(:), POINTER              :: nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER           :: first_sgfa, first_sgfb
      INTEGER, ALLOCATABLE                        :: atom_of_kind(:)
      REAL(KIND=dp), ALLOCATABLE                  :: stmp(:), native_s(:, :), native_h(:, :), native_lattr(:, :)
      REAL(KIND=dp), ALLOCATABLE                  :: dtmp(:, :), qtmp(:, :), dtmpj(:, :), qtmpj(:, :)
      REAL(KIND=dp), ALLOCATABLE                  :: native_dip(:, :, :), native_quad(:, :, :)
      REAL(KIND=dp), DIMENSION(:, :), POINTER     :: dip_ket1, dip_ket2, dip_ket3, &
                                                     dip_bra1, dip_bra2, dip_bra3
      REAL(KIND=dp), DIMENSION(:, :), POINTER     :: quad_ket1, quad_ket2, quad_ket3, &
                                                     quad_ket4, quad_ket5, quad_ket6, &
                                                     quad_bra1, quad_bra2, quad_bra3, &
                                                     quad_bra4, quad_bra5, quad_bra6

      TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                       :: dft_control
      TYPE(dbcsr_iterator_type)                             :: iter
      TYPE(gto_basis_set_type), POINTER                     :: basis_set_a, basis_set_b
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER     :: basis_set_list
      TYPE(adjacency_list)                                  :: native_list
      TYPE(kpoint_type), POINTER                            :: kpoints
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                              :: nl_iterator
      TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER             :: qs_kind_set

      CALL timeset(routineN, handle)

      !get info from environment vaiarable
      NULLIFY (atomic_kind_set, qs_kind_set, sab_orb, particle_set)
      NULLIFY (dft_control, matrix_s, kpoints, cell_to_index)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      sab_orb=sab_orb, &
                      particle_set=particle_set, &
                      dft_control=dft_control, &
                      kpoints=kpoints, &
                      matrix_s_kp=matrix_s)
      natom = SIZE(particle_set)
      nkind = SIZE(atomic_kind_set)
      nimg = dft_control%nimages
      IF (nimg > 1) THEN
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
      END IF
      mp_test_mode = 0
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_MULTIPOLE_TEST", debug_value, STATUS=debug_status)
      IF (debug_status == 0) READ (debug_value, *, IOSTAT=debug_status) mp_test_mode
#endif

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind)

      !set up basis set lists
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      ALLOCATE (stmp(msao(tb%calc%bas%maxl)**2))
      ALLOCATE (dtmp(dip_n, msao(tb%calc%bas%maxl)**2))
      ALLOCATE (qtmp(quad_n, msao(tb%calc%bas%maxl)**2))
      ALLOCATE (dtmpj(dip_n, msao(tb%calc%bas%maxl)**2))
      ALLOCATE (qtmpj(quad_n, msao(tb%calc%bas%maxl)**2))

      ! allocate dipole/quadrupole moment matrix elemnts
      CALL dbcsr_allocate_matrix_set(tb%dipbra, dip_n*nimg)
      CALL dbcsr_allocate_matrix_set(tb%dipket, dip_n*nimg)
      CALL dbcsr_allocate_matrix_set(tb%quadbra, quad_n*nimg)
      CALL dbcsr_allocate_matrix_set(tb%quadket, quad_n*nimg)
      DO img = 1, nimg
         DO i = 1, dip_n
            idx = i + dip_n*(img - 1)
            ALLOCATE (tb%dipbra(idx)%matrix)
            ALLOCATE (tb%dipket(idx)%matrix)
            CALL dbcsr_create(tb%dipbra(idx)%matrix, template=matrix_s(1, img)%matrix, &
                              name="DIPOLE BRAMATRIX")
            CALL dbcsr_create(tb%dipket(idx)%matrix, template=matrix_s(1, img)%matrix, &
                              name="DIPOLE KETMATRIX")
            CALL cp_dbcsr_alloc_block_from_nbl(tb%dipbra(idx)%matrix, sab_orb)
            CALL cp_dbcsr_alloc_block_from_nbl(tb%dipket(idx)%matrix, sab_orb)
         END DO
         DO i = 1, quad_n
            idx = i + quad_n*(img - 1)
            ALLOCATE (tb%quadbra(idx)%matrix)
            ALLOCATE (tb%quadket(idx)%matrix)
            CALL dbcsr_create(tb%quadbra(idx)%matrix, template=matrix_s(1, img)%matrix, &
                              name="QUADRUPOLE BRAMATRIX")
            CALL dbcsr_create(tb%quadket(idx)%matrix, template=matrix_s(1, img)%matrix, &
                              name="QUADRUPOLE KETMATRIX")
            CALL cp_dbcsr_alloc_block_from_nbl(tb%quadbra(idx)%matrix, sab_orb)
            CALL cp_dbcsr_alloc_block_from_nbl(tb%quadket(idx)%matrix, sab_orb)
         END DO
      END DO

      IF (mp_test_mode == 80 .OR. mp_test_mode == 82) THEN
         IF (nimg > 1) CPABORT("Native tblite multipole test with k-points not implemented")
         ALLOCATE (native_s(tb%calc%bas%nao, tb%calc%bas%nao))
         ALLOCATE (native_h(tb%calc%bas%nao, tb%calc%bas%nao))
         ALLOCATE (native_dip(dip_n, tb%calc%bas%nao, tb%calc%bas%nao))
         ALLOCATE (native_quad(quad_n, tb%calc%bas%nao, tb%calc%bas%nao))
         CALL get_lattice_points(tb%mol%periodic, tb%mol%lattice, get_cutoff(tb%calc%bas, 1.0_dp), native_lattr)
         CALL new_adjacency_list(native_list, tb%mol, native_lattr, get_cutoff(tb%calc%bas, 1.0_dp))
         CALL get_hamiltonian(tb%mol, native_lattr, native_list, tb%calc%bas, tb%calc%h0, tb%selfenergy, &
                              native_s, native_dip, native_quad, native_h)

         CALL dbcsr_iterator_start(iter, tb%dipbra(1)%matrix)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, irow, icol, dip_bra1)
            CALL dbcsr_get_block_p(matrix=tb%dipbra(2)%matrix, row=irow, col=icol, BLOCK=dip_bra2, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%dipbra(3)%matrix, row=irow, col=icol, BLOCK=dip_bra3, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%dipket(1)%matrix, row=irow, col=icol, BLOCK=dip_ket1, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%dipket(2)%matrix, row=irow, col=icol, BLOCK=dip_ket2, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%dipket(3)%matrix, row=irow, col=icol, BLOCK=dip_ket3, found=found)
            CPASSERT(found)

            CALL dbcsr_get_block_p(matrix=tb%quadbra(1)%matrix, row=irow, col=icol, BLOCK=quad_bra1, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%quadbra(2)%matrix, row=irow, col=icol, BLOCK=quad_bra2, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%quadbra(3)%matrix, row=irow, col=icol, BLOCK=quad_bra3, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%quadbra(4)%matrix, row=irow, col=icol, BLOCK=quad_bra4, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%quadbra(5)%matrix, row=irow, col=icol, BLOCK=quad_bra5, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%quadbra(6)%matrix, row=irow, col=icol, BLOCK=quad_bra6, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%quadket(1)%matrix, row=irow, col=icol, BLOCK=quad_ket1, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%quadket(2)%matrix, row=irow, col=icol, BLOCK=quad_ket2, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%quadket(3)%matrix, row=irow, col=icol, BLOCK=quad_ket3, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%quadket(4)%matrix, row=irow, col=icol, BLOCK=quad_ket4, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%quadket(5)%matrix, row=irow, col=icol, BLOCK=quad_ket5, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=tb%quadket(6)%matrix, row=irow, col=icol, BLOCK=quad_ket6, found=found)
            CPASSERT(found)

            i0 = tb%calc%bas%iao_sh(tb%calc%bas%ish_at(irow) + 1)
            j0 = tb%calc%bas%iao_sh(tb%calc%bas%ish_at(icol) + 1)
            DO j = 1, SIZE(dip_bra1, 1)
               DO i = 1, SIZE(dip_bra1, 2)
                  dip_bra1(j, i) = native_dip(1, i0 + j, j0 + i)
                  dip_bra2(j, i) = native_dip(2, i0 + j, j0 + i)
                  dip_bra3(j, i) = native_dip(3, i0 + j, j0 + i)
                  dip_ket1(j, i) = native_dip(1, j0 + i, i0 + j)
                  dip_ket2(j, i) = native_dip(2, j0 + i, i0 + j)
                  dip_ket3(j, i) = native_dip(3, j0 + i, i0 + j)

                  quad_bra1(j, i) = native_quad(1, i0 + j, j0 + i)
                  quad_bra2(j, i) = native_quad(2, i0 + j, j0 + i)
                  quad_bra3(j, i) = native_quad(3, i0 + j, j0 + i)
                  quad_bra4(j, i) = native_quad(4, i0 + j, j0 + i)
                  quad_bra5(j, i) = native_quad(5, i0 + j, j0 + i)
                  quad_bra6(j, i) = native_quad(6, i0 + j, j0 + i)
                  quad_ket1(j, i) = native_quad(1, j0 + i, i0 + j)
                  quad_ket2(j, i) = native_quad(2, j0 + i, i0 + j)
                  quad_ket3(j, i) = native_quad(3, j0 + i, i0 + j)
                  quad_ket4(j, i) = native_quad(4, j0 + i, i0 + j)
                  quad_ket5(j, i) = native_quad(5, j0 + i, i0 + j)
                  quad_ket6(j, i) = native_quad(6, j0 + i, i0 + j)
               END DO
            END DO
         END DO
         CALL dbcsr_iterator_stop(iter)

         DO i = 1, dip_n
            CALL dbcsr_finalize(tb%dipbra(i)%matrix)
            CALL dbcsr_finalize(tb%dipket(i)%matrix)
         END DO
         DO i = 1, quad_n
            CALL dbcsr_finalize(tb%quadbra(i)%matrix)
            CALL dbcsr_finalize(tb%quadket(i)%matrix)
         END DO

         DEALLOCATE (native_s, native_h, native_dip, native_quad, native_lattr, basis_set_list)
         CALL timestop(handle)
         RETURN
      END IF

      !loop over all atom pairs with a non-zero overlap (sab_orb)
      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)

         r2 = NORM2(rij(:))**2

         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)

         IF (icol == jatom) THEN
            rij = -rij
            i = ikind
            ikind = jkind
            jkind = i
         END IF

         IF (nimg == 1) THEN
            ic = 1
         ELSE
            ic = cell_to_index(cell(1), cell(2), cell(3))
            CPASSERT(ic > 0)
         END IF
         id1 = 1 + dip_n*(ic - 1)
         id2 = 2 + dip_n*(ic - 1)
         id3 = 3 + dip_n*(ic - 1)
         iq1 = 1 + quad_n*(ic - 1)
         iq2 = 2 + quad_n*(ic - 1)
         iq3 = 3 + quad_n*(ic - 1)
         iq4 = 4 + quad_n*(ic - 1)
         iq5 = 5 + quad_n*(ic - 1)
         iq6 = 6 + quad_n*(ic - 1)

         ityp = tb%mol%id(icol)
         jtyp = tb%mol%id(irow)

         NULLIFY (dip_bra1, dip_bra2, dip_bra3)
         CALL dbcsr_get_block_p(matrix=tb%dipbra(id1)%matrix, &
                                row=irow, col=icol, BLOCK=dip_bra1, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%dipbra(id2)%matrix, &
                                row=irow, col=icol, BLOCK=dip_bra2, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%dipbra(id3)%matrix, &
                                row=irow, col=icol, BLOCK=dip_bra3, found=found)
         CPASSERT(found)

         NULLIFY (dip_ket1, dip_ket2, dip_ket3)
         CALL dbcsr_get_block_p(matrix=tb%dipket(id1)%matrix, &
                                row=irow, col=icol, BLOCK=dip_ket1, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%dipket(id2)%matrix, &
                                row=irow, col=icol, BLOCK=dip_ket2, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%dipket(id3)%matrix, &
                                row=irow, col=icol, BLOCK=dip_ket3, found=found)
         CPASSERT(found)

         NULLIFY (quad_bra1, quad_bra2, quad_bra3, quad_bra4, quad_bra5, quad_bra6)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(iq1)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra1, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(iq2)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra2, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(iq3)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra3, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(iq4)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra4, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(iq5)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra5, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadbra(iq6)%matrix, &
                                row=irow, col=icol, BLOCK=quad_bra6, found=found)
         CPASSERT(found)

         NULLIFY (quad_ket1, quad_ket2, quad_ket3, quad_ket4, quad_ket5, quad_ket6)
         CALL dbcsr_get_block_p(matrix=tb%quadket(iq1)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket1, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(iq2)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket2, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(iq3)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket3, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(iq4)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket4, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(iq5)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket5, found=found)
         CPASSERT(found)
         CALL dbcsr_get_block_p(matrix=tb%quadket(iq6)%matrix, &
                                row=irow, col=icol, BLOCK=quad_ket6, found=found)
         CPASSERT(found)

         !get basis information
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         atom_a = atom_of_kind(icol)
         atom_b = atom_of_kind(irow)
         ! basis a
         first_sgfa => basis_set_a%first_sgf
         la_max => basis_set_a%lmax
         nseta = basis_set_a%nset
         nsgfa => basis_set_a%nsgf_set
         ! basis b
         first_sgfb => basis_set_b%first_sgf
         lb_max => basis_set_b%lmax
         nsetb = basis_set_b%nset
         nsgfb => basis_set_b%nsgf_set

         ! --------- Hamiltonian
         ! Periodic self-images are off-diagonal lattice contributions, not the on-site block.
         IF (icol == irow .AND. r2 < same_atom**2) THEN
            DO iset = 1, nseta
               DO jset = 1, nsetb
                  CALL multipole_cgto(tb%calc%bas%cgto(jset, ityp), tb%calc%bas%cgto(iset, ityp), &
                       & r2, -rij, tb%calc%bas%intcut, stmp, dtmp, qtmp)

                  DO inda = 1, nsgfa(iset)
                     ia = first_sgfa(1, iset) - first_sgfa(1, 1) + inda
                     DO indb = 1, nsgfb(jset)
                        ib = first_sgfb(1, jset) - first_sgfb(1, 1) + indb
                        ij = indb + nsgfb(jset)*(inda - 1)

                        dip_ket1(ib, ia) = dip_ket1(ib, ia) + dtmp(1, ij)
                        dip_ket2(ib, ia) = dip_ket2(ib, ia) + dtmp(2, ij)
                        dip_ket3(ib, ia) = dip_ket3(ib, ia) + dtmp(3, ij)

                        quad_ket1(ib, ia) = quad_ket1(ib, ia) + qtmp(1, ij)
                        quad_ket2(ib, ia) = quad_ket2(ib, ia) + qtmp(2, ij)
                        quad_ket3(ib, ia) = quad_ket3(ib, ia) + qtmp(3, ij)
                        quad_ket4(ib, ia) = quad_ket4(ib, ia) + qtmp(4, ij)
                        quad_ket5(ib, ia) = quad_ket5(ib, ia) + qtmp(5, ij)
                        quad_ket6(ib, ia) = quad_ket6(ib, ia) + qtmp(6, ij)

                        dip_bra1(ib, ia) = dip_bra1(ib, ia) + dtmp(1, ij)
                        dip_bra2(ib, ia) = dip_bra2(ib, ia) + dtmp(2, ij)
                        dip_bra3(ib, ia) = dip_bra3(ib, ia) + dtmp(3, ij)

                        quad_bra1(ib, ia) = quad_bra1(ib, ia) + qtmp(1, ij)
                        quad_bra2(ib, ia) = quad_bra2(ib, ia) + qtmp(2, ij)
                        quad_bra3(ib, ia) = quad_bra3(ib, ia) + qtmp(3, ij)
                        quad_bra4(ib, ia) = quad_bra4(ib, ia) + qtmp(4, ij)
                        quad_bra5(ib, ia) = quad_bra5(ib, ia) + qtmp(5, ij)
                        quad_bra6(ib, ia) = quad_bra6(ib, ia) + qtmp(6, ij)
                     END DO
                  END DO
               END DO
            END DO
         ELSE
            DO iset = 1, nseta
               DO jset = 1, nsetb
                  mp_int_vec = -rij
                  mp_shift_vec = -rij
                  IF (mp_test_mode == 14) THEN
                     mp_int_vec = rij
                     mp_shift_vec = rij
                  ELSEIF (mp_test_mode == 15) THEN
                     mp_shift_vec = rij
                  ELSEIF (mp_test_mode == 16) THEN
                     mp_int_vec = rij
                  END IF
                  CALL multipole_cgto(tb%calc%bas%cgto(jset, jtyp), tb%calc%bas%cgto(iset, ityp), &
                      & r2, mp_int_vec, tb%calc%bas%intcut, stmp, dtmp, qtmp)

                  DO inda = 1, nsgfa(iset)
                     ia = first_sgfa(1, iset) - first_sgfa(1, 1) + inda
                     DO indb = 1, nsgfb(jset)
                        ib = first_sgfb(1, jset) - first_sgfb(1, 1) + indb

                        ij = indb + nsgfb(jset)*(inda - 1)
                        CALL tb_shift_multipole(mp_shift_vec, stmp(ij), dtmp(:, ij), qtmp(:, ij), &
                                                dtmpj(:, ij), qtmpj(:, ij))
                        IF (icol == irow .AND. r2 >= same_atom**2) THEN
                           IF (nimg > 1) THEN
                              ! Match the k-point image block orientation used for self-image quadrupoles.
                              qtmp(:, ij) = -qtmp(:, ij)
                              qtmpj(:, ij) = -qtmpj(:, ij)
                           END IF
                           IF (mp_test_mode == 30 .OR. mp_test_mode == 32 .OR. mp_test_mode == 33) THEN
                              dtmpj(:, ij) = dtmp(:, ij)
                              qtmpj(:, ij) = qtmp(:, ij)
                           ELSEIF (mp_test_mode == 31) THEN
                              CALL tb_shift_multipole(-mp_shift_vec, stmp(ij), dtmp(:, ij), qtmp(:, ij), &
                                                      dtmpj(:, ij), qtmpj(:, ij))
                           END IF
                           IF (mp_test_mode == 40 .OR. mp_test_mode == 41 .OR. mp_test_mode == 44) THEN
                              dtmp(:, ij) = -dtmp(:, ij)
                           END IF
                           IF (mp_test_mode == 40 .OR. mp_test_mode == 42 .OR. mp_test_mode == 44) THEN
                              dtmpj(:, ij) = -dtmpj(:, ij)
                           END IF
                           IF (mp_test_mode == 40 .OR. mp_test_mode == 41 .OR. mp_test_mode == 45) THEN
                              qtmp(:, ij) = -qtmp(:, ij)
                           END IF
                           IF (mp_test_mode == 40 .OR. mp_test_mode == 42 .OR. mp_test_mode == 45) THEN
                              qtmpj(:, ij) = -qtmpj(:, ij)
                           END IF
                           IF (mp_test_mode == 43) THEN
                              dtmp(:, ij) = 0.0_dp
                              qtmp(:, ij) = 0.0_dp
                              dtmpj(:, ij) = 0.0_dp
                              qtmpj(:, ij) = 0.0_dp
                           END IF
                        END IF
                        IF (r2 >= same_atom**2) THEN
                           IF (mp_test_mode == 46 .OR. mp_test_mode == 48) THEN
                              dtmp(:, ij) = -dtmp(:, ij)
                              dtmpj(:, ij) = -dtmpj(:, ij)
                           END IF
                           IF (mp_test_mode == 46 .OR. mp_test_mode == 49) THEN
                              qtmp(:, ij) = -qtmp(:, ij)
                              qtmpj(:, ij) = -qtmpj(:, ij)
                           END IF
                           IF (mp_test_mode == 47) THEN
                              dtmp(:, ij) = 0.0_dp
                              qtmp(:, ij) = 0.0_dp
                              dtmpj(:, ij) = 0.0_dp
                              qtmpj(:, ij) = 0.0_dp
                           END IF
                           IF (mp_test_mode == 50) THEN
                              dtmp(:, ij) = 0.25_dp*dtmp(:, ij)
                              qtmp(:, ij) = 0.25_dp*qtmp(:, ij)
                              dtmpj(:, ij) = 0.25_dp*dtmpj(:, ij)
                              qtmpj(:, ij) = 0.25_dp*qtmpj(:, ij)
                           ELSEIF (mp_test_mode == 51) THEN
                              dtmp(:, ij) = 0.5_dp*dtmp(:, ij)
                              qtmp(:, ij) = 0.5_dp*qtmp(:, ij)
                              dtmpj(:, ij) = 0.5_dp*dtmpj(:, ij)
                              qtmpj(:, ij) = 0.5_dp*qtmpj(:, ij)
                           ELSEIF (mp_test_mode == 52) THEN
                              dtmp(:, ij) = 0.75_dp*dtmp(:, ij)
                              qtmp(:, ij) = 0.75_dp*qtmp(:, ij)
                              dtmpj(:, ij) = 0.75_dp*dtmpj(:, ij)
                              qtmpj(:, ij) = 0.75_dp*qtmpj(:, ij)
                           ELSEIF (mp_test_mode == 53) THEN
                              dtmp(:, ij) = -0.25_dp*dtmp(:, ij)
                              qtmp(:, ij) = -0.25_dp*qtmp(:, ij)
                              dtmpj(:, ij) = -0.25_dp*dtmpj(:, ij)
                              qtmpj(:, ij) = -0.25_dp*qtmpj(:, ij)
                           ELSEIF (mp_test_mode == 54) THEN
                              dtmp(:, ij) = -0.5_dp*dtmp(:, ij)
                              qtmp(:, ij) = -0.5_dp*qtmp(:, ij)
                              dtmpj(:, ij) = -0.5_dp*dtmpj(:, ij)
                              qtmpj(:, ij) = -0.5_dp*qtmpj(:, ij)
                           ELSEIF (mp_test_mode == 55) THEN
                              dtmp(:, ij) = -0.75_dp*dtmp(:, ij)
                              qtmp(:, ij) = -0.75_dp*qtmp(:, ij)
                              dtmpj(:, ij) = -0.75_dp*dtmpj(:, ij)
                              qtmpj(:, ij) = -0.75_dp*qtmpj(:, ij)
                           END IF
                        END IF

                        dip_bra1(ib, ia) = dip_bra1(ib, ia) + dtmp(1, ij)
                        dip_bra2(ib, ia) = dip_bra2(ib, ia) + dtmp(2, ij)
                        dip_bra3(ib, ia) = dip_bra3(ib, ia) + dtmp(3, ij)

                        quad_bra1(ib, ia) = quad_bra1(ib, ia) + qtmp(1, ij)
                        quad_bra2(ib, ia) = quad_bra2(ib, ia) + qtmp(2, ij)
                        quad_bra3(ib, ia) = quad_bra3(ib, ia) + qtmp(3, ij)
                        quad_bra4(ib, ia) = quad_bra4(ib, ia) + qtmp(4, ij)
                        quad_bra5(ib, ia) = quad_bra5(ib, ia) + qtmp(5, ij)
                        quad_bra6(ib, ia) = quad_bra6(ib, ia) + qtmp(6, ij)

                        dip_ket1(ib, ia) = dip_ket1(ib, ia) + dtmpj(1, ij)
                        dip_ket2(ib, ia) = dip_ket2(ib, ia) + dtmpj(2, ij)
                        dip_ket3(ib, ia) = dip_ket3(ib, ia) + dtmpj(3, ij)

                        quad_ket1(ib, ia) = quad_ket1(ib, ia) + qtmpj(1, ij)
                        quad_ket2(ib, ia) = quad_ket2(ib, ia) + qtmpj(2, ij)
                        quad_ket3(ib, ia) = quad_ket3(ib, ia) + qtmpj(3, ij)
                        quad_ket4(ib, ia) = quad_ket4(ib, ia) + qtmpj(4, ij)
                        quad_ket5(ib, ia) = quad_ket5(ib, ia) + qtmpj(5, ij)
                        quad_ket6(ib, ia) = quad_ket6(ib, ia) + qtmpj(6, ij)
                     END DO
                  END DO
               END DO
            END DO
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO i = 1, SIZE(tb%dipbra)
         CALL dbcsr_finalize(tb%dipbra(i)%matrix)
         CALL dbcsr_finalize(tb%dipket(i)%matrix)
      END DO
      DO i = 1, SIZE(tb%quadbra)
         CALL dbcsr_finalize(tb%quadbra(i)%matrix)
         CALL dbcsr_finalize(tb%quadket(i)%matrix)
      END DO

      DEALLOCATE (basis_set_list)

      CALL timestop(handle)

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_get_multipole

! **************************************************************************************************
!> \brief Shift a multipole operator from one center to the other.
!> \param vec displacement vector between the two centers
!> \param s overlap integral
!> \param di dipole integral on the original center
!> \param qi quadrupole integral on the original center
!> \param dj dipole integral on the shifted center
!> \param qj quadrupole integral on the shifted center
! **************************************************************************************************
   PURE SUBROUTINE tb_shift_multipole(vec, s, di, qi, dj, qj)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: vec
      REAL(KIND=dp), INTENT(IN)                          :: s
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: di, qi
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: dj, qj

      REAL(KIND=dp)                                      :: tr

      dj(1) = di(1) + vec(1)*s
      dj(2) = di(2) + vec(2)*s
      dj(3) = di(3) + vec(3)*s

      qj(1) = 2*vec(1)*di(1) + vec(1)**2*s
      qj(3) = 2*vec(2)*di(2) + vec(2)**2*s
      qj(6) = 2*vec(3)*di(3) + vec(3)**2*s
      qj(2) = vec(1)*di(2) + vec(2)*di(1) + vec(1)*vec(2)*s
      qj(4) = vec(1)*di(3) + vec(3)*di(1) + vec(1)*vec(3)*s
      qj(5) = vec(2)*di(3) + vec(3)*di(2) + vec(2)*vec(3)*s
      tr = 0.5_dp*(qj(1) + qj(3) + qj(6))

      qj(1) = qi(1) + 1.5_dp*qj(1) - tr
      qj(2) = qi(2) + 1.5_dp*qj(2)
      qj(3) = qi(3) + 1.5_dp*qj(3) - tr
      qj(4) = qi(4) + 1.5_dp*qj(4)
      qj(5) = qi(5) + 1.5_dp*qj(5)
      qj(6) = qi(6) + 1.5_dp*qj(6) - tr

   END SUBROUTINE tb_shift_multipole

! **************************************************************************************************
!> \brief compute the mulliken properties (AO resolved)
!> \param p_matrix ...
!> \param bra_mat ...
!> \param ket_mat ...
!> \param output ...
!> \param para_env ...
!> \par History
!>      adapted from ao_charges_2
!> \note
! **************************************************************************************************
   SUBROUTINE contract_dens(p_matrix, bra_mat, ket_mat, output, para_env)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_matrix
      TYPE(dbcsr_type), POINTER                          :: bra_mat, ket_mat
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: output
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'contract_dens'

      INTEGER                                            :: handle, i, iblock_col, &
                                                            iblock_row, ispin, j, mp_test_mode, &
                                                            nspin
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      INTEGER                                            :: debug_status
      CHARACTER(LEN=32)                                  :: debug_value
#endif
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: bra, ket, p_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      mp_test_mode = 0
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_MULTIPOLE_TEST", debug_value, STATUS=debug_status)
      IF (debug_status == 0) READ (debug_value, *, IOSTAT=debug_status) mp_test_mode
#endif

      nspin = SIZE(p_matrix)
      output = 0.0_dp
      DO ispin = 1, nspin
         CALL dbcsr_iterator_start(iter, bra_mat)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            NULLIFY (p_block, bra, ket)

            CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, bra)
            CALL dbcsr_get_block_p(matrix=p_matrix(ispin)%matrix, &
                                   row=iblock_row, col=iblock_col, BLOCK=p_block, found=found)
            IF (.NOT. found) CYCLE
            CALL dbcsr_get_block_p(matrix=ket_mat, &
                                   row=iblock_row, col=iblock_col, BLOCK=ket, found=found)
            IF (.NOT. found) CPABORT("missing block")

            IF (.NOT. (ASSOCIATED(bra) .AND. ASSOCIATED(p_block))) CYCLE
            IF (iblock_col == iblock_row) THEN
               IF (mp_test_mode == 70) THEN
                  DO j = 1, SIZE(p_block, 1)
                     DO i = 1, SIZE(p_block, 2)
                        output(iblock_row) = output(iblock_row) + p_block(j, i)*ket(j, i)
                     END DO
                  END DO
               ELSEIF (mp_test_mode == 71) THEN
                  DO j = 1, SIZE(p_block, 1)
                     DO i = 1, SIZE(p_block, 2)
                        output(iblock_row) = output(iblock_row) + 0.5_dp*p_block(j, i)*(bra(j, i) + ket(j, i))
                     END DO
                  END DO
               ELSE
                  DO j = 1, SIZE(p_block, 1)
                     DO i = 1, SIZE(p_block, 2)
                        output(iblock_row) = output(iblock_row) + p_block(j, i)*bra(j, i)
                     END DO
                  END DO
               END IF
            ELSE
               DO j = 1, SIZE(p_block, 1)
                  DO i = 1, SIZE(p_block, 2)
                     output(iblock_row) = output(iblock_row) + p_block(j, i)*ket(j, i)
                  END DO
               END DO
               DO j = 1, SIZE(p_block, 1)
                  DO i = 1, SIZE(p_block, 2)
                     output(iblock_col) = output(iblock_col) + p_block(j, i)*bra(j, i)
                  END DO
               END DO
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)
      END DO

      CALL para_env%sum(output)
      CALL timestop(handle)

   END SUBROUTINE contract_dens

! **************************************************************************************************
!> \brief compute the AO-resolved density contraction for real-space k-point image matrices
!> \param p_matrix ...
!> \param bra_mat ...
!> \param ket_mat ...
!> \param iop ...
!> \param nops ...
!> \param output ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE contract_dens_kp(p_matrix, bra_mat, ket_mat, iop, nops, output, para_env)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: p_matrix
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: bra_mat, ket_mat
      INTEGER, INTENT(IN)                                :: iop, nops
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: output
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: ic, idx, nimg
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: image_output
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_image

      nimg = SIZE(p_matrix, 2)
      output = 0.0_dp
      ALLOCATE (image_output(SIZE(output)))
      DO ic = 1, nimg
         idx = iop + nops*(ic - 1)
         CPASSERT(idx <= SIZE(bra_mat))
         CPASSERT(idx <= SIZE(ket_mat))
         NULLIFY (p_image)
         p_image => p_matrix(:, ic)
         image_output = 0.0_dp
         CALL contract_dens(p_image, bra_mat(idx)%matrix, ket_mat(idx)%matrix, image_output, para_env)
         output = output + image_output
      END DO
      DEALLOCATE (image_output)

   END SUBROUTINE contract_dens_kp

! **************************************************************************************************
!> \brief save gradient to force
!> \param qs_env ...
!> \param tb ...
!> \param para_env ...
!> \param ityp ...
!> \note
! **************************************************************************************************
   SUBROUTINE tb_grad2force(qs_env, tb, para_env, ityp)

      TYPE(qs_environment_type)                          :: qs_env
      TYPE(tblite_type)                                  :: tb
      TYPE(mp_para_env_type)                             :: para_env
      INTEGER                                            :: ityp

      CHARACTER(len=*), PARAMETER                        :: routineN = 'tb_grad2force'

      CHARACTER(LEN=default_path_length)                 :: dump_file
      INTEGER                                            :: atoma, dump_status, dump_unit, handle, &
                                                            iatom, ikind, natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force

      CALL timeset(routineN, handle)

      NULLIFY (force, atomic_kind_set)
      CALL get_qs_env(qs_env=qs_env, force=force, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of)

      natom = SIZE(particle_set)

      dump_status = 1
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_FORCE_DUMP", dump_file, STATUS=dump_status)
#endif
      IF (dump_status == 0) THEN
         OPEN (NEWUNIT=dump_unit, FILE=TRIM(dump_file), STATUS="UNKNOWN", &
               POSITION="APPEND", ACTION="WRITE")
         WRITE (dump_unit, "(A,1X,I0)") "component", ityp
         DO iatom = 1, natom
            WRITE (dump_unit, "(I0,3(1X,ES24.16))") iatom, tb%grad(:, iatom)/para_env%num_pe
         END DO
         CLOSE (dump_unit)
      END IF

      SELECT CASE (ityp)
      CASE DEFAULT
         CPABORT("unknown force type")
      CASE (0)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%all_potential(:, atoma) = &
               force(ikind)%all_potential(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (1)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%repulsive(:, atoma) = &
               force(ikind)%repulsive(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (2)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%dispersion(:, atoma) = &
               force(ikind)%dispersion(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (3)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%rho_elec(:, atoma) = &
               force(ikind)%rho_elec(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (4)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%overlap(:, atoma) = &
               force(ikind)%overlap(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      CASE (5)
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atoma = atom_of_kind(iatom)
            force(ikind)%efield(:, atoma) = &
               force(ikind)%efield(:, atoma) + tb%grad(:, iatom)/para_env%num_pe
         END DO
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE tb_grad2force

! **************************************************************************************************
!> \brief set gradient to zero
!> \param qs_env ...
!> \note
! **************************************************************************************************
   SUBROUTINE tb_zero_force(qs_env)

      TYPE(qs_environment_type)                          :: qs_env

      INTEGER                                            :: iatom, ikind, natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force

      NULLIFY (force, atomic_kind_set)
      CALL get_qs_env(qs_env=qs_env, force=force, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               kind_of=kind_of)

      natom = SIZE(particle_set)

      DO iatom = 1, natom
         ikind = kind_of(iatom)
         force(ikind)%all_potential = 0.0_dp
         force(ikind)%repulsive = 0.0_dp
         force(ikind)%dispersion = 0.0_dp
         force(ikind)%rho_elec = 0.0_dp
         force(ikind)%overlap = 0.0_dp
         force(ikind)%efield = 0.0_dp
      END DO

   END SUBROUTINE tb_zero_force

! **************************************************************************************************
!> \brief compute the derivative of the energy
!> \param qs_env ...
!> \param use_rho ...
!> \param nimg ...
! **************************************************************************************************
   SUBROUTINE tb_derive_dH_diag(qs_env, use_rho, nimg)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: use_rho
      INTEGER, INTENT(IN)                                :: nimg

#if defined(__TBLITE)
      INTEGER                                            :: i, iatom, ic, ikind, ind, is, ish, &
                                                            jatom, jkind
      INTEGER, DIMENSION(3)                              :: cellind
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dE
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pblock
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(tblite_type), POINTER                         :: tb

      ! compute mulliken charges required for charge update
      NULLIFY (scf_env, rho, tb, sab_orb, para_env, kpoints)
      CALL get_qs_env(qs_env=qs_env, scf_env=scf_env, rho=rho, tb_tblite=tb, &
                      sab_orb=sab_orb, para_env=para_env, kpoints=kpoints)

      NULLIFY (cell_to_index)
      IF (nimg > 1) THEN
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
      END IF

      NULLIFY (matrix_p)
      IF (use_rho) THEN
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      ELSE
         matrix_p => scf_env%p_mix_new
      END IF

      ALLOCATE (dE(tb%mol%nat))

      dE = 0.0_dp
      ! loop over all atom pairs with a non-zero overlap (sab_orb)
      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, cell=cellind)

         IF (iatom /= jatom) CYCLE

         IF (ikind /= jkind) CPABORT("Type wrong")

         is = tb%calc%bas%ish_at(iatom)

         IF (nimg == 1) THEN
            ic = 1
         ELSE
            ic = cell_to_index(cellind(1), cellind(2), cellind(3))
            CPASSERT(ic > 0)
         END IF

         IF (cellind(1) /= 0) CYCLE
         IF (cellind(2) /= 0) CYCLE
         IF (cellind(3) /= 0) CYCLE

         CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                row=iatom, col=jatom, BLOCK=pblock, found=found)

         IF (.NOT. found) CPABORT("block not found")

         ind = 0
         DO ish = 1, tb%calc%bas%nsh_id(ikind)
            DO i = 1, msao(tb%calc%bas%cgto(ish, ikind)%ang)
               ind = ind + 1
               dE(iatom) = dE(iatom) + tb%dsedcn(is + ish)*pblock(ind, ind)
            END DO
         END DO
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)
      CALL para_env%sum(dE)

      tb%grad = 0.0_dp
      CALL tb_add_grad(tb%grad, tb%dcndr, dE, tb%mol%nat)
      IF (tb%use_virial) CALL tb_add_sig(tb%sigma, tb%dcndL, dE, tb%mol%nat)
      CALL tb_grad2force(qs_env, tb, para_env, 4)
      DEALLOCATE (dE)

#else
      MARK_USED(qs_env)
      MARK_USED(use_rho)
      MARK_USED(nimg)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_derive_dH_diag

! **************************************************************************************************
!> \brief compute the derivative of the energy
!> \param qs_env ...
!> \param use_rho ...
!> \param nimg ...
! **************************************************************************************************
   SUBROUTINE tb_derive_dH_off(qs_env, use_rho, nimg)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: use_rho
      INTEGER, INTENT(IN)                                :: nimg

#if defined(__TBLITE)
      INTEGER                                            :: i, ij, iatom, ic, icol, ikind, ni, nj, nkind, nel, &
                                                            sgfi, sgfj, ityp, jatom, jkind, jrow, jtyp, iset, jset, &
                                                            nseti, nsetj, ia, ib, inda, indb, sampled_odd_axes
      INTEGER, DIMENSION(3)                              :: cellind, nkp_grid
      INTEGER, DIMENSION(:), POINTER                     :: nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found, image_pair, sampled_image_pair
      REAL(KIND=dp)                                      :: r2, dr, itemp, jtemp, rr, hij, shpoly, dshpoly, idHdc, jdHdc, &
                                                            scal, hp, i_a_shift, j_a_shift, ishift, jshift, &
                                                            dip_deriv_fac, quad_deriv_fac, dip_deriv_fac_pair, &
                                                            quad_deriv_fac_pair, overlap_deriv_fac, pot_deriv_scale, &
                                                            w_deriv_scale, h0_deriv_scale, mp_deriv_scale, &
                                                            pot_image_deriv_scale, w_image_deriv_scale, &
                                                            h0_image_deriv_scale, mp_image_deriv_scale, &
                                                            pot_pair_scale, w_pair_scale, h0_pair_scale, mp_pair_scale, &
                                                            pot_self_image_deriv_scale, w_self_image_deriv_scale, &
                                                            h0_self_image_deriv_scale, mp_self_image_deriv_scale, &
                                                            cn_image_scale
      INTEGER                                            :: mp_test_mode
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      INTEGER                                            :: debug_status
      CHARACTER(LEN=32)                                  :: debug_value, scale_value
#endif
      REAL(KIND=dp), DIMENSION(3)                        :: rij, dgrad
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hsigma
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: dE, t_ov, idip, jdip, iquad, jquad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: t_dip, t_quad, t_d_ov
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: t_i_dip, t_i_quad, t_j_dip, t_j_quad
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pblock, wblock

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_w
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER            :: basis_set_a, basis_set_b
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(tb_hamiltonian), POINTER                      :: h0
      TYPE(tblite_type), POINTER                         :: tb

      ! compute mulliken charges required for charge update
      NULLIFY (scf_env, rho, tb, sab_orb, para_env, kpoints, matrix_w)
      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      scf_env=scf_env, &
                      rho=rho, &
                      tb_tblite=tb, &
                      sab_orb=sab_orb, &
                      para_env=para_env, &
                      qs_kind_set=qs_kind_set, &
                      kpoints=kpoints, &
                      matrix_w_kp=matrix_w)

      NULLIFY (cell_to_index)
      IF (nimg > 1) THEN
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index, nkp_grid=nkp_grid)
      ELSE
         nkp_grid = 1
      END IF
      sampled_odd_axes = 0
      IF (nkp_grid(1) > 1 .AND. MODULO(nkp_grid(1), 2) == 1) sampled_odd_axes = sampled_odd_axes + 1
      IF (nkp_grid(2) > 1 .AND. MODULO(nkp_grid(2), 2) == 1) sampled_odd_axes = sampled_odd_axes + 1
      IF (nkp_grid(3) > 1 .AND. MODULO(nkp_grid(3), 2) == 1) sampled_odd_axes = sampled_odd_axes + 1

      h0 => tb%calc%h0
      mp_test_mode = 0
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_MULTIPOLE_TEST", debug_value, STATUS=debug_status)
      IF (debug_status == 0) READ (debug_value, *, IOSTAT=debug_status) mp_test_mode
#endif
      dip_deriv_fac = -1.0_dp
      quad_deriv_fac = -1.0_dp
      IF (mp_test_mode == 7 .OR. mp_test_mode == 11) THEN
         dip_deriv_fac = 1.0_dp
         quad_deriv_fac = 1.0_dp
      END IF
      IF (mp_test_mode == 8 .OR. mp_test_mode == 12) THEN
         dip_deriv_fac = 0.0_dp
         quad_deriv_fac = 0.0_dp
      END IF
      IF (mp_test_mode == 13 .OR. mp_test_mode == 14 .OR. mp_test_mode == 15 .OR. &
          mp_test_mode == 16 .OR. mp_test_mode == 21 .OR. mp_test_mode == 33 .OR. &
          mp_test_mode == 34 .OR. mp_test_mode == 37) THEN
         dip_deriv_fac = -0.5_dp
         quad_deriv_fac = -0.5_dp
      END IF
      IF (mp_test_mode == 17) THEN
         dip_deriv_fac = -1.0_dp
         quad_deriv_fac = 0.0_dp
      END IF
      IF (mp_test_mode == 18) THEN
         dip_deriv_fac = 0.0_dp
         quad_deriv_fac = -1.0_dp
      END IF
      pot_deriv_scale = 1.0_dp
      w_deriv_scale = 1.0_dp
      h0_deriv_scale = 1.0_dp
      mp_deriv_scale = 1.0_dp
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DH_POT_SCALE", scale_value, STATUS=debug_status)
      IF (debug_status == 0) READ (scale_value, *, IOSTAT=debug_status) pot_deriv_scale
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DH_W_SCALE", scale_value, STATUS=debug_status)
      IF (debug_status == 0) READ (scale_value, *, IOSTAT=debug_status) w_deriv_scale
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DH_H0_SCALE", scale_value, STATUS=debug_status)
      IF (debug_status == 0) READ (scale_value, *, IOSTAT=debug_status) h0_deriv_scale
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DH_MP_SCALE", scale_value, STATUS=debug_status)
      IF (debug_status == 0) READ (scale_value, *, IOSTAT=debug_status) mp_deriv_scale
#endif
      pot_image_deriv_scale = pot_deriv_scale
      w_image_deriv_scale = w_deriv_scale
      h0_image_deriv_scale = h0_deriv_scale
      mp_image_deriv_scale = mp_deriv_scale
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DH_POT_IMAGE_SCALE", scale_value, STATUS=debug_status)
      IF (debug_status == 0) READ (scale_value, *, IOSTAT=debug_status) pot_image_deriv_scale
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DH_W_IMAGE_SCALE", scale_value, STATUS=debug_status)
      IF (debug_status == 0) READ (scale_value, *, IOSTAT=debug_status) w_image_deriv_scale
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DH_H0_IMAGE_SCALE", scale_value, STATUS=debug_status)
      IF (debug_status == 0) READ (scale_value, *, IOSTAT=debug_status) h0_image_deriv_scale
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DH_MP_IMAGE_SCALE", scale_value, STATUS=debug_status)
      IF (debug_status == 0) READ (scale_value, *, IOSTAT=debug_status) mp_image_deriv_scale
#endif
      pot_self_image_deriv_scale = pot_image_deriv_scale
      w_self_image_deriv_scale = w_image_deriv_scale
      h0_self_image_deriv_scale = h0_image_deriv_scale
      mp_self_image_deriv_scale = mp_image_deriv_scale
      ! Periodic self-images are stored in the sorted on-site image block for
      ! k-point calculations. Their multipole derivatives need the same
      ! opposite quarter-weighting as the explicit self-image virial
      ! contribution below. The Gamma periodic bsf path keeps the native factor.
      IF (nimg > 1 .AND. ALL(tb%mol%periodic)) &
         mp_self_image_deriv_scale = -0.25_dp*mp_image_deriv_scale
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DH_POT_SELF_IMAGE_SCALE", scale_value, STATUS=debug_status)
      IF (debug_status == 0) READ (scale_value, *, IOSTAT=debug_status) pot_self_image_deriv_scale
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DH_W_SELF_IMAGE_SCALE", scale_value, STATUS=debug_status)
      IF (debug_status == 0) READ (scale_value, *, IOSTAT=debug_status) w_self_image_deriv_scale
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DH_H0_SELF_IMAGE_SCALE", scale_value, STATUS=debug_status)
      IF (debug_status == 0) READ (scale_value, *, IOSTAT=debug_status) h0_self_image_deriv_scale
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DH_MP_SELF_IMAGE_SCALE", scale_value, STATUS=debug_status)
      IF (debug_status == 0) READ (scale_value, *, IOSTAT=debug_status) mp_self_image_deriv_scale
#endif
      cn_image_scale = 1.0_dp
      IF (nimg > 1 .AND. ANY(tb%mol%periodic)) cn_image_scale = -0.25_dp
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DH_CN_IMAGE_SCALE", scale_value, STATUS=debug_status)
      IF (debug_status == 0) READ (scale_value, *, IOSTAT=debug_status) cn_image_scale
#endif

      NULLIFY (matrix_p)
      IF (use_rho) THEN
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      ELSE
         matrix_p => scf_env%p_mix_new
      END IF

      ! set up basis set lists
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      ALLOCATE (dE(tb%mol%nat))

      nel = msao(tb%calc%bas%maxl)**2
      ALLOCATE (t_ov(nel))
      ALLOCATE (t_d_ov(3, nel))
      ALLOCATE (t_dip(dip_n, nel))
      ALLOCATE (t_i_dip(3, dip_n, nel), t_j_dip(3, dip_n, nel))
      ALLOCATE (t_quad(quad_n, nel))
      ALLOCATE (t_i_quad(3, quad_n, nel), t_j_quad(3, quad_n, nel))

      ALLOCATE (idip(dip_n), jdip(dip_n))
      ALLOCATE (iquad(quad_n), jquad(quad_n))

      dE = 0.0_dp
      tb%grad = 0.0_dp
      hsigma = 0.0_dp
      ! loop over all atom pairs with a non-zero overlap (sab_orb)
      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cellind)

         icol = MAX(iatom, jatom)
         jrow = MIN(iatom, jatom)

         IF (icol == jatom) THEN
            rij = -rij
            i = ikind
            ikind = jkind
            jkind = i
         END IF

         ityp = tb%mol%id(icol)
         jtyp = tb%mol%id(jrow)

         r2 = DOT_PRODUCT(rij, rij)
         dr = SQRT(r2)
         IF (icol == jrow .AND. dr < same_atom) CYCLE
         IF (mp_test_mode == 140 .AND. icol == jrow) CYCLE
         rr = SQRT(dr/(h0%rad(ikind) + h0%rad(jkind)))

         !get basis information
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         first_sgfa => basis_set_a%first_sgf
         nsgfa => basis_set_a%nsgf_set
         nseti = basis_set_a%nset
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
         first_sgfb => basis_set_b%first_sgf
         nsgfb => basis_set_b%nsgf_set
         nsetj = basis_set_b%nset

         IF (nimg == 1) THEN
            ic = 1
         ELSE
            ic = cell_to_index(cellind(1), cellind(2), cellind(3))
            CPASSERT(ic > 0)
         END IF
         image_pair = (cellind(1) /= 0 .OR. cellind(2) /= 0 .OR. cellind(3) /= 0)
         sampled_image_pair = image_pair .AND. &
                              ((cellind(1) /= 0 .AND. nkp_grid(1) > 1) .OR. &
                               (cellind(2) /= 0 .AND. nkp_grid(2) > 1) .OR. &
                               (cellind(3) /= 0 .AND. nkp_grid(3) > 1))
         pot_pair_scale = MERGE(pot_image_deriv_scale, pot_deriv_scale, image_pair)
         w_pair_scale = MERGE(w_image_deriv_scale, w_deriv_scale, image_pair)
         h0_pair_scale = MERGE(h0_image_deriv_scale, h0_deriv_scale, image_pair)
         mp_pair_scale = MERGE(mp_image_deriv_scale, mp_deriv_scale, image_pair)
         IF (icol == jrow .AND. sampled_image_pair) THEN
            pot_pair_scale = pot_self_image_deriv_scale
            w_pair_scale = w_self_image_deriv_scale
            ! Odd Monkhorst-Pack meshes do not contain a Nyquist image along the
            ! sampled axis. Match the self-image W derivative to the reduced
            ! real-space representation used by the k-point density transform.
            w_pair_scale = w_pair_scale*(1.0_dp - 0.125_dp*REAL(sampled_odd_axes, dp))
            h0_pair_scale = h0_self_image_deriv_scale
            mp_pair_scale = mp_self_image_deriv_scale
         END IF

         NULLIFY (pblock)
         CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                row=jrow, col=icol, BLOCK=pblock, found=found)
         IF (.NOT. found) CPABORT("pblock not found")

         NULLIFY (wblock)
         CALL dbcsr_get_block_p(matrix=matrix_w(1, ic)%matrix, &
                                row=jrow, col=icol, block=wblock, found=found)
         CPASSERT(found)

         i_a_shift = tb%pot%vat(icol, 1)
         j_a_shift = tb%pot%vat(jrow, 1)
         idip(:) = tb%pot%vdp(:, icol, 1)
         jdip(:) = tb%pot%vdp(:, jrow, 1)
         iquad(:) = tb%pot%vqp(:, icol, 1)
         jquad(:) = tb%pot%vqp(:, jrow, 1)
         dip_deriv_fac_pair = dip_deriv_fac
         quad_deriv_fac_pair = quad_deriv_fac
         overlap_deriv_fac = 1.0_dp
         IF (icol == jrow) THEN
            IF (mp_test_mode == 141) THEN
               dip_deriv_fac_pair = 0.0_dp
               quad_deriv_fac_pair = 0.0_dp
            ELSEIF (mp_test_mode == 142) THEN
               overlap_deriv_fac = 0.0_dp
            END IF
         END IF

         ni = tb%calc%bas%ish_at(icol)
         DO iset = 1, nseti
            sgfi = first_sgfa(1, iset)
            ishift = i_a_shift + tb%pot%vsh(ni + iset, 1)
            nj = tb%calc%bas%ish_at(jrow)
            DO jset = 1, nsetj
               sgfj = first_sgfb(1, jset)
               jshift = j_a_shift + tb%pot%vsh(nj + jset, 1)

               !get integrals and derivatives
               IF (mp_test_mode == 20 .OR. mp_test_mode == 21) THEN
                  CALL multipole_grad_cgto(tb%calc%bas%cgto(jset, jtyp), tb%calc%bas%cgto(iset, ityp), &
                     & r2, rij, tb%calc%bas%intcut, t_ov, t_dip, t_quad, t_d_ov, t_j_dip, t_j_quad, &
                     & t_i_dip, t_i_quad)
               ELSEIF (mp_test_mode == 34 .OR. mp_test_mode == 35) THEN
                  CALL multipole_grad_cgto(tb%calc%bas%cgto(jset, jtyp), tb%calc%bas%cgto(iset, ityp), &
                     & r2, -rij, tb%calc%bas%intcut, t_ov, t_dip, t_quad, t_d_ov, t_j_dip, t_j_quad, &
                     & t_i_dip, t_i_quad)
               ELSEIF (mp_test_mode == 36 .OR. mp_test_mode == 37) THEN
                  CALL multipole_grad_cgto(tb%calc%bas%cgto(iset, ityp), tb%calc%bas%cgto(jset, jtyp), &
                     & r2, -rij, tb%calc%bas%intcut, t_ov, t_dip, t_quad, t_d_ov, t_i_dip, t_i_quad, &
                     & t_j_dip, t_j_quad)
               ELSE
                  CALL multipole_grad_cgto(tb%calc%bas%cgto(iset, ityp), tb%calc%bas%cgto(jset, jtyp), &
                     & r2, rij, tb%calc%bas%intcut, t_ov, t_dip, t_quad, t_d_ov, t_i_dip, t_i_quad, &
                     & t_j_dip, t_j_quad)
               END IF

               shpoly = (1.0_dp + h0%shpoly(iset, ikind)*rr) &
                        *(1.0_dp + h0%shpoly(jset, jkind)*rr)
               dshpoly = ((1.0_dp + h0%shpoly(iset, ikind)*rr)*h0%shpoly(jset, jkind)*rr &
                          + (1.0_dp + h0%shpoly(jset, jkind)*rr)*h0%shpoly(iset, ikind)*rr) &
                         *0.5_dp/r2
               scal = h0%hscale(iset, jset, ikind, jkind)
               hij = 0.5_dp*(tb%selfenergy(ni + iset) + tb%selfenergy(nj + jset))*scal

               idHdc = tb%dsedcn(ni + iset)*shpoly*scal
               jdHdc = tb%dsedcn(nj + jset)*shpoly*scal

               itemp = 0.0_dp
               jtemp = 0.0_dp
               dgrad = 0.0_dp
               DO inda = 1, nsgfa(iset)
                  ia = first_sgfa(1, iset) - first_sgfa(1, 1) + inda
                  DO indb = 1, nsgfb(jset)
                     ib = first_sgfb(1, jset) - first_sgfb(1, 1) + indb

                     IF (mp_test_mode == 20 .OR. mp_test_mode == 21 .OR. &
                         mp_test_mode == 34 .OR. mp_test_mode == 35) THEN
                        ij = indb + nsgfb(jset)*(inda - 1)
                     ELSE
                        ij = inda + nsgfa(iset)*(indb - 1)
                     END IF

                     itemp = itemp + overlap_deriv_fac*idHdc*pblock(ib, ia)*t_ov(ij)
                     jtemp = jtemp + overlap_deriv_fac*jdHdc*pblock(ib, ia)*t_ov(ij)

                     hp = 2*hij*pblock(ib, ia)
                     dgrad(:) = dgrad(:) &
                                + overlap_deriv_fac*(-pot_pair_scale*(ishift + jshift)*pblock(ib, ia)*t_d_ov(:, ij) &
                                                     - 2*w_pair_scale*wblock(ib, ia)*t_d_ov(:, ij) &
                                                     + h0_pair_scale*(hp*shpoly*t_d_ov(:, ij) &
                                                                      + hp*dshpoly*t_ov(ij)*rij)) &
                                + mp_pair_scale*pblock(ib, ia)*( &
                                dip_deriv_fac_pair*(MATMUL(t_i_dip(:, :, ij), idip) &
                                                    + MATMUL(t_j_dip(:, :, ij), jdip)) &
                                + quad_deriv_fac_pair*(MATMUL(t_i_quad(:, :, ij), iquad) &
                                                       + MATMUL(t_j_quad(:, :, ij), jquad)))

                  END DO
               END DO
               ! Periodic self-image pairs share one atom in the sorted CP2K block representation.
               IF (icol == jrow) THEN
                  IF (sampled_image_pair) THEN
                     dE(icol) = dE(icol) + cn_image_scale*0.5_dp*(itemp + jtemp)
                  ELSE
                     dE(icol) = dE(icol) + 0.5_dp*(itemp + jtemp)
                  END IF
               ELSE
                  dE(icol) = dE(icol) + itemp
                  dE(jrow) = dE(jrow) + jtemp
               END IF
               tb%grad(:, icol) = tb%grad(:, icol) - dgrad
               tb%grad(:, jrow) = tb%grad(:, jrow) + dgrad
               IF (tb%use_virial) THEN
                  IF (icol == jrow) THEN
                     DO ia = 1, 3
                        DO ib = 1, 3
                           IF (sampled_image_pair) THEN
                              hsigma(ia, ib) = hsigma(ia, ib) - 0.25_dp*(rij(ia)*dgrad(ib) + rij(ib)*dgrad(ia))
                           ELSE
                              hsigma(ia, ib) = hsigma(ia, ib) + 0.25_dp*(rij(ia)*dgrad(ib) + rij(ib)*dgrad(ia))
                           END IF
                        END DO
                     END DO
                  ELSE
                     DO ia = 1, 3
                        DO ib = 1, 3
                           hsigma(ia, ib) = hsigma(ia, ib) + 0.50_dp*(rij(ia)*dgrad(ib) + rij(ib)*dgrad(ia))
                        END DO
                     END DO
                  END IF
               END IF
            END DO
         END DO
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      CALL para_env%sum(hsigma)
      CALL para_env%sum(dE)
      CALL para_env%sum(tb%grad)

      CALL tb_add_grad(tb%grad, tb%dcndr, dE, tb%mol%nat)
      CALL tb_grad2force(qs_env, tb, para_env, 4)

      CALL tb_dump_sigma_component("before_h0_off", tb%sigma, para_env)
      tb%sigma = tb%sigma + hsigma
      CALL tb_dump_sigma_component("after_h0_off_direct", tb%sigma, para_env)
      IF (tb%use_virial) THEN
         CALL tb_add_sig(tb%sigma, tb%dcndL, dE, tb%mol%nat)
         CALL tb_dump_sigma_component("after_h0_off_cn", tb%sigma, para_env)
      END IF

      DEALLOCATE (dE)
      DEALLOCATE (basis_set_list)
      DEALLOCATE (t_ov, t_d_ov)
      DEALLOCATE (t_dip, t_i_dip, t_j_dip)
      DEALLOCATE (t_quad, t_i_quad, t_j_quad)
      DEALLOCATE (idip, jdip, iquad, jquad)

      IF (tb%use_virial) CALL tb_add_stress(qs_env, tb, para_env)

#else
      MARK_USED(qs_env)
      MARK_USED(use_rho)
      MARK_USED(nimg)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_derive_dH_off

! **************************************************************************************************
!> \brief Run native tblite CLI and compare against CP2K/tblite.
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE tb_reference_cli_compare(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tb_reference_cli_compare'

      CHARACTER(LEN=32)                                  :: acc_str, charge_str, iter_str, spin_str
      CHARACTER(LEN=8)                                   :: method
      CHARACTER(LEN=8*default_path_length)               :: command
      CHARACTER(LEN=default_path_length)                 :: file_base, gen_file, grad_file, &
                                                            json_file, log_file
      INTEGER                                            :: cmdstat, exitstat, handle, iounit, &
                                                            natom, nkp, spin
      LOGICAL                                            :: do_kpoints, have_energy, have_gradient, &
                                                            have_virial, too_large, &
                                                            unsupported_kpoints
      REAL(KIND=dp)                                      :: cli_energy, cp_energy, ediff, fmax, &
                                                            fsum, vmax, vsum
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: cli_gradient, cli_virial, cp_gradient
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xtb_reference_cli_type)                       :: ref

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, dft_control, energy, force, kpoints, logger, para_env, particle_set, &
               virial, xkp)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                      dft_control=dft_control, energy=energy, force=force, &
                      para_env=para_env, particle_set=particle_set, virial=virial, &
                      do_kpoints=do_kpoints, kpoints=kpoints)

      ref = dft_control%qs_control%xtb_control%reference_cli
      IF (.NOT. ref%enabled) THEN
         CALL timestop(handle)
         RETURN
      END IF
      IF (.NOT. para_env%is_source()) THEN
         CALL timestop(handle)
         RETURN
      END IF

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)
      unsupported_kpoints = .FALSE.
      IF (do_kpoints .AND. ASSOCIATED(kpoints)) THEN
         nkp = 0
         CALL get_kpoint_info(kpoint=kpoints, nkp=nkp, xkp=xkp)
         unsupported_kpoints = nkp > 1
         IF (nkp == 1 .AND. ASSOCIATED(xkp)) unsupported_kpoints = ANY(ABS(xkp(:, 1)) > 1.0E-12_dp)
      END IF
      IF (unsupported_kpoints) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") &
            "tblite reference CLI check skipped: CP2K KPOINTS are active."
         WRITE (UNIT=iounit, FMT="(T2,A)") &
            "The native tblite CLI reference path does not reproduce CP2K multi-k-point sampling."
         IF (ref%stop_on_error) CPABORT("tblite reference CLI cannot check CP2K k-point calculations")
         CALL timestop(handle)
         RETURN
      END IF
      natom = SIZE(particle_set)
      method = tb_reference_method_name(dft_control%qs_control%xtb_control%tblite_method)
      file_base = tb_join_path(ref%work_directory, ref%prefix)
      gen_file = TRIM(file_base)//".gen"
      grad_file = TRIM(file_base)//".grad"
      json_file = TRIM(file_base)//".json"
      log_file = TRIM(file_base)//".log"

      WRITE (charge_str, "(I0)") dft_control%charge
      spin = MAX(0, dft_control%multiplicity - 1)
      WRITE (spin_str, "(I0)") spin
      WRITE (acc_str, "(ES16.8)") ref%accuracy
      WRITE (iter_str, "(I0)") ref%iterations

      CALL tb_write_reference_gen(qs_env, TRIM(gen_file))

      command = TRIM(tb_shell_quote(ref%program_name))//" run --method "//TRIM(method)// &
                " --charge "//TRIM(ADJUSTL(charge_str))// &
                " --spin "//TRIM(ADJUSTL(spin_str))// &
                " --acc "//TRIM(ADJUSTL(acc_str))// &
                " --iterations "//TRIM(ADJUSTL(iter_str))// &
                " --etemp 0.0"// &
                " --no-restart --input gen --grad "//TRIM(tb_shell_quote(grad_file))// &
                " --json "//TRIM(tb_shell_quote(json_file))//" "//TRIM(tb_shell_quote(gen_file))// &
                " > "//TRIM(tb_shell_quote(log_file))//" 2>&1"

      cmdstat = 0
      exitstat = 0
      CALL execute_command_line(TRIM(command), exitstat=exitstat, cmdstat=cmdstat)
      IF (cmdstat /= 0 .OR. exitstat /= 0) THEN
         WRITE (UNIT=iounit, FMT="(/,T2,A)") "tblite reference CLI check failed to run."
         WRITE (UNIT=iounit, FMT="(T2,A,A)") "Command: ", TRIM(command)
         WRITE (UNIT=iounit, FMT="(T2,A,I0,T32,A,I0)") "cmdstat:", cmdstat, "exitstat:", exitstat
         IF (ref%stop_on_error) CPABORT("tblite reference CLI command failed")
         CALL tb_reference_cleanup(ref, gen_file, grad_file, json_file, log_file)
         CALL timestop(handle)
         RETURN
      END IF

      ALLOCATE (cli_gradient(3, natom), cli_virial(3, 3))
      CALL tb_read_reference_grad(TRIM(grad_file), natom, cli_energy, cli_gradient, cli_virial, &
                                  have_energy, have_gradient, have_virial)

      WRITE (UNIT=iounit, FMT="(/,T2,A)") "tblite reference CLI check"
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Executable: ", TRIM(ref%program_name)
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Method:     ", TRIM(method)
      WRITE (UNIT=iounit, FMT="(T2,A,A)") "Log file:   ", TRIM(log_file)

      too_large = .FALSE.
      IF (ref%check_energy) THEN
         IF (have_energy) THEN
            cp_energy = energy%total
            ediff = ABS(cp_energy - cli_energy)
            WRITE (UNIT=iounit, FMT="(T2,A,3ES22.12)") &
               "Energy CP2K/CLI/absdiff:", cp_energy, cli_energy, ediff
            too_large = too_large .OR. ediff > ref%error_limit
         ELSE
            WRITE (UNIT=iounit, FMT="(T2,A)") "Energy check skipped: no CLI energy found."
         END IF
      END IF

      IF (ref%check_forces) THEN
         IF (have_gradient .AND. ASSOCIATED(force)) THEN
            ALLOCATE (cp_gradient(3, natom))
            CALL total_qs_force(cp_gradient, force, atomic_kind_set)
            fsum = SUM(ABS(cp_gradient - cli_gradient))
            fmax = MAXVAL(ABS(cp_gradient - cli_gradient))
            WRITE (UNIT=iounit, FMT="(T2,A,2ES22.12)") "Gradient diff sum/max:", fsum, fmax
            too_large = too_large .OR. fmax > ref%error_limit
            DEALLOCATE (cp_gradient)
         ELSE
            WRITE (UNIT=iounit, FMT="(T2,A)") "Gradient check skipped: no CLI gradient or CP2K force found."
         END IF
      END IF

      IF (ref%check_virial) THEN
         IF (have_virial .AND. ASSOCIATED(virial)) THEN
            ! Native tblite prints the positive cell derivative; CP2K stores the PV virial
            ! with the opposite sign.
            vsum = SUM(ABS(-virial%pv_virial - cli_virial))
            vmax = MAXVAL(ABS(-virial%pv_virial - cli_virial))
            WRITE (UNIT=iounit, FMT="(T2,A,2ES22.12)") "Virial diff sum/max:", vsum, vmax
            too_large = too_large .OR. vmax > ref%error_limit
         ELSE
            WRITE (UNIT=iounit, FMT="(T2,A)") "Virial check skipped: no CLI virial or CP2K virial found."
         END IF
      END IF

      IF (too_large) THEN
         WRITE (UNIT=iounit, FMT="(T2,A,ES12.4)") &
            "tblite reference CLI deviation exceeded ERROR_LIMIT = ", ref%error_limit
         IF (ref%stop_on_error) CPABORT("tblite reference CLI deviation exceeded ERROR_LIMIT")
      END IF

      CALL tb_reference_cleanup(ref, gen_file, grad_file, json_file, log_file)
      DEALLOCATE (cli_gradient, cli_virial)

      CALL timestop(handle)

   END SUBROUTINE tb_reference_cli_compare

! **************************************************************************************************
!> \brief Map CP2K tblite method id to native tblite CLI method name.
!> \param method_id ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_reference_method_name(method_id) RESULT(method)
      INTEGER, INTENT(IN)                                :: method_id
      CHARACTER(LEN=8)                                   :: method

      SELECT CASE (method_id)
      CASE (gfn1xtb)
         method = "gfn1"
      CASE (gfn2xtb)
         method = "gfn2"
      CASE (ipea1xtb)
         method = "ipea1"
      CASE DEFAULT
         CPABORT("Unknown tblite reference CLI method")
      END SELECT

   END FUNCTION tb_reference_method_name

! **************************************************************************************************
!> \brief Join directory and filename.
!> \param directory ...
!> \param filename ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_join_path(directory, filename) RESULT(path)
      CHARACTER(LEN=*), INTENT(IN)                       :: directory, filename
      CHARACTER(LEN=default_path_length)                 :: path

      IF (LEN_TRIM(directory) == 0 .OR. TRIM(directory) == ".") THEN
         path = TRIM(filename)
      ELSE IF (directory(LEN_TRIM(directory):LEN_TRIM(directory)) == "/") THEN
         path = TRIM(directory)//TRIM(filename)
      ELSE
         path = TRIM(directory)//"/"//TRIM(filename)
      END IF

   END FUNCTION tb_join_path

! **************************************************************************************************
!> \brief Shell-quote a filename or executable path.
!> \param text ...
!> \return ...
! **************************************************************************************************
   FUNCTION tb_shell_quote(text) RESULT(quoted)
      CHARACTER(LEN=*), INTENT(IN)                       :: text
      CHARACTER(LEN=4*default_path_length)               :: quoted

      INTEGER                                            :: i

      quoted = "'"
      DO i = 1, LEN_TRIM(text)
         IF (text(i:i) == "'") THEN
            quoted = TRIM(quoted)//"'\\''"
         ELSE
            quoted = TRIM(quoted)//text(i:i)
         END IF
      END DO
      quoted = TRIM(quoted)//"'"

   END FUNCTION tb_shell_quote

! **************************************************************************************************
!> \brief Write current CP2K geometry as DFTB+ gen format for native tblite.
!> \param qs_env ...
!> \param filename ...
! **************************************************************************************************
   SUBROUTINE tb_write_reference_gen(qs_env, filename)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(LEN=*), INTENT(IN)                       :: filename

      CHARACTER(LEN=2), ALLOCATABLE, DIMENSION(:)        :: symbols, unique_symbols
      INTEGER                                            :: iatom, ikind, ios, natom, nuniq, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: species
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: to_angstrom
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      NULLIFY (cell, particle_set)
      CALL get_qs_env(qs_env=qs_env, cell=cell, particle_set=particle_set)

      natom = SIZE(particle_set)
      to_angstrom = cp_unit_from_cp2k(1.0_dp, "angstrom")
      ALLOCATE (symbols(natom), unique_symbols(natom), species(natom))
      nuniq = 0
      DO iatom = 1, natom
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, element_symbol=symbols(iatom))
         found = .FALSE.
         DO ikind = 1, nuniq
            IF (TRIM(unique_symbols(ikind)) == TRIM(symbols(iatom))) THEN
               found = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT. found) THEN
            nuniq = nuniq + 1
            unique_symbols(nuniq) = symbols(iatom)
            ikind = nuniq
         END IF
         species(iatom) = ikind
      END DO

      OPEN (NEWUNIT=unit_nr, FILE=TRIM(filename), STATUS="REPLACE", ACTION="WRITE", &
            FORM="FORMATTED", IOSTAT=ios)
      IF (ios /= 0) CPABORT("Could not open tblite reference CLI geometry file")

      CALL get_cell(cell=cell, periodic=periodic)
      IF (ANY(periodic == 1)) THEN
         WRITE (UNIT=unit_nr, FMT="(I0,1X,A)") natom, "S"
      ELSE
         WRITE (UNIT=unit_nr, FMT="(I0,1X,A)") natom, "C"
      END IF
      WRITE (UNIT=unit_nr, FMT="(*(A,1X))") (TRIM(unique_symbols(ikind)), ikind=1, nuniq)
      DO iatom = 1, natom
         WRITE (UNIT=unit_nr, FMT="(I0,1X,I0,3(1X,ES24.16))") &
            iatom, species(iatom), particle_set(iatom)%r(:)*to_angstrom
      END DO
      IF (ANY(periodic == 1)) THEN
         WRITE (UNIT=unit_nr, FMT="(3(1X,ES24.16))") 0.0_dp, 0.0_dp, 0.0_dp
         DO ikind = 1, 3
            WRITE (UNIT=unit_nr, FMT="(3(1X,ES24.16))") cell%hmat(:, ikind)*to_angstrom
         END DO
      END IF
      CLOSE (unit_nr)

      DEALLOCATE (symbols, unique_symbols, species)

   END SUBROUTINE tb_write_reference_gen

! **************************************************************************************************
!> \brief Read native tblite gradient file.
!> \param filename ...
!> \param natom ...
!> \param energy ...
!> \param gradient ...
!> \param virial ...
!> \param have_energy ...
!> \param have_gradient ...
!> \param have_virial ...
! **************************************************************************************************
   SUBROUTINE tb_read_reference_grad(filename, natom, energy, gradient, virial, &
                                     have_energy, have_gradient, have_virial)

      CHARACTER(LEN=*), INTENT(IN)                       :: filename
      INTEGER, INTENT(IN)                                :: natom
      REAL(KIND=dp), INTENT(OUT)                         :: energy
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: gradient, virial
      LOGICAL, INTENT(OUT)                               :: have_energy, have_gradient, have_virial

      CHARACTER(LEN=1024)                                :: line
      INTEGER                                            :: ios, nread, unit_nr
      LOGICAL                                            :: exists
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: values

      have_energy = .FALSE.
      have_gradient = .FALSE.
      have_virial = .FALSE.
      energy = 0.0_dp
      gradient = 0.0_dp
      virial = 0.0_dp

      INQUIRE (FILE=TRIM(filename), EXIST=exists)
      IF (.NOT. exists) RETURN

      OPEN (NEWUNIT=unit_nr, FILE=TRIM(filename), STATUS="OLD", ACTION="READ", &
            FORM="FORMATTED", IOSTAT=ios)
      IF (ios /= 0) RETURN

      DO
         READ (UNIT=unit_nr, FMT="(A)", IOSTAT=ios) line
         IF (ios /= 0) EXIT
         IF (INDEX(line, "energy             :real:0:") > 0) THEN
            READ (UNIT=unit_nr, FMT="(A)", IOSTAT=ios) line
            IF (ios == 0) THEN
               READ (line, *, IOSTAT=ios) energy
               have_energy = ios == 0
            END IF
         ELSE IF (INDEX(line, "gradient           :real:2:3,") > 0) THEN
            ALLOCATE (values(3*natom))
            CALL tb_read_real_values(unit_nr, values, nread)
            IF (nread == 3*natom) THEN
               CALL tb_values_to_matrix(values, gradient)
               have_gradient = .TRUE.
            END IF
            DEALLOCATE (values)
         ELSE IF (INDEX(line, "virial             :real:2:3,3") > 0) THEN
            ALLOCATE (values(9))
            CALL tb_read_real_values(unit_nr, values, nread)
            IF (nread == 9) THEN
               CALL tb_values_to_matrix(values, virial)
               have_virial = .TRUE.
            END IF
            DEALLOCATE (values)
         END IF
      END DO
      CLOSE (unit_nr)

   END SUBROUTINE tb_read_reference_grad

! **************************************************************************************************
!> \brief Read a fixed number of real values from following lines.
!> \param unit_nr ...
!> \param values ...
!> \param nread ...
! **************************************************************************************************
   SUBROUTINE tb_read_real_values(unit_nr, values, nread)

      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: values
      INTEGER, INTENT(OUT)                               :: nread

      CHARACTER(LEN=1024)                                :: line
      INTEGER                                            :: ios

      nread = 0
      DO WHILE (nread < SIZE(values))
         READ (UNIT=unit_nr, FMT="(A)", IOSTAT=ios) line
         IF (ios /= 0) EXIT
         CALL tb_parse_real_line(line, values, nread)
      END DO

   END SUBROUTINE tb_read_real_values

! **************************************************************************************************
!> \brief Parse real values from one text line.
!> \param line ...
!> \param values ...
!> \param nread ...
! **************************************************************************************************
   SUBROUTINE tb_parse_real_line(line, values, nread)

      CHARACTER(LEN=*), INTENT(IN)                       :: line
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: values
      INTEGER, INTENT(INOUT)                             :: nread

      CHARACTER(LEN=128)                                 :: token
      INTEGER                                            :: first, ios, last, pos

      pos = 1
      DO WHILE (pos <= LEN_TRIM(line) .AND. nread < SIZE(values))
         DO WHILE (pos <= LEN_TRIM(line) .AND. INDEX(" ,[]", line(pos:pos)) > 0)
            pos = pos + 1
         END DO
         IF (pos > LEN_TRIM(line)) EXIT
         first = pos
         DO WHILE (pos <= LEN_TRIM(line) .AND. INDEX(" ,[]", line(pos:pos)) == 0)
            pos = pos + 1
         END DO
         last = pos - 1
         token = line(first:last)
         READ (token, *, IOSTAT=ios) values(nread + 1)
         IF (ios == 0) nread = nread + 1
      END DO

   END SUBROUTINE tb_parse_real_line

! **************************************************************************************************
!> \brief Convert flat native tblite values to CP2K atom-major matrix layout.
!> \param values ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE tb_values_to_matrix(values, matrix)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: values
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: matrix

      INTEGER                                            :: i, j, n

      n = 0
      DO j = 1, SIZE(matrix, 2)
         DO i = 1, SIZE(matrix, 1)
            n = n + 1
            matrix(i, j) = values(n)
         END DO
      END DO

   END SUBROUTINE tb_values_to_matrix

! **************************************************************************************************
!> \brief Delete temporary files unless requested otherwise.
!> \param ref ...
!> \param gen_file ...
!> \param grad_file ...
!> \param json_file ...
!> \param log_file ...
! **************************************************************************************************
   SUBROUTINE tb_reference_cleanup(ref, gen_file, grad_file, json_file, log_file)

      TYPE(xtb_reference_cli_type), INTENT(IN)           :: ref
      CHARACTER(LEN=*), INTENT(IN)                       :: gen_file, grad_file, json_file, log_file

      IF (ref%keep_files) RETURN
      CALL tb_delete_file(gen_file)
      CALL tb_delete_file(grad_file)
      CALL tb_delete_file(json_file)
      CALL tb_delete_file(log_file)

   END SUBROUTINE tb_reference_cleanup

! **************************************************************************************************
!> \brief Delete a file if it exists.
!> \param filename ...
! **************************************************************************************************
   SUBROUTINE tb_delete_file(filename)

      CHARACTER(LEN=*), INTENT(IN)                       :: filename

      INTEGER                                            :: ios, unit_nr
      LOGICAL                                            :: exists

      INQUIRE (FILE=TRIM(filename), EXIST=exists)
      IF (.NOT. exists) RETURN
      OPEN (NEWUNIT=unit_nr, FILE=TRIM(filename), STATUS="OLD", IOSTAT=ios)
      IF (ios == 0) CLOSE (unit_nr, STATUS="DELETE")

   END SUBROUTINE tb_delete_file

! **************************************************************************************************
!> \brief Dump cumulative tblite virial pieces for local debugging.
!> \param label ...
!> \param sigma ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_dump_sigma_component(label, sigma, para_env)

      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sigma
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      CHARACTER(LEN=default_path_length)                 :: dump_file
      INTEGER                                            :: dump_status, dump_unit, i

      dump_status = 1
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_SIGMA_COMPONENT_DUMP", dump_file, STATUS=dump_status)
#endif
      IF (dump_status /= 0) RETURN

      OPEN (NEWUNIT=dump_unit, FILE=TRIM(dump_file), STATUS="UNKNOWN", &
            POSITION="APPEND", ACTION="WRITE")
      WRITE (dump_unit, "(A)") TRIM(label)
      DO i = 1, 3
         WRITE (dump_unit, "(3(1X,ES24.16))") sigma(i, :)/para_env%num_pe
      END DO
      CLOSE (dump_unit)

   END SUBROUTINE tb_dump_sigma_component

! **************************************************************************************************
!> \brief save stress tensor
!> \param qs_env ...
!> \param tb ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_add_stress(qs_env, tb, para_env)

      TYPE(qs_environment_type)                          :: qs_env
      TYPE(tblite_type)                                  :: tb
      TYPE(mp_para_env_type)                             :: para_env

      CHARACTER(LEN=default_path_length)                 :: dump_file
      INTEGER                                            :: dump_status, dump_unit, i
      INTEGER, DIMENSION(3)                              :: periodic
      TYPE(cell_type), POINTER                           :: cell
      TYPE(virial_type), POINTER                         :: virial

      NULLIFY (virial, cell)
      CALL get_qs_env(qs_env=qs_env, virial=virial, cell=cell)
      CALL get_cell(cell=cell, periodic=periodic)

      IF (ALL(periodic == 0)) THEN
         CALL cp_warn(__LOCATION__, &
                      "tblite stress tensor requested for an isolated system. "// &
                      "The reported virial is useful for finite-difference checks, "// &
                      "but it is not a physically meaningful bulk stress for an isolated molecule.")
      END IF

      dump_status = 1
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_VIRIAL_DUMP", dump_file, STATUS=dump_status)
#endif
      IF (dump_status == 0) THEN
         OPEN (NEWUNIT=dump_unit, FILE=TRIM(dump_file), STATUS="UNKNOWN", &
               POSITION="APPEND", ACTION="WRITE")
         WRITE (dump_unit, "(A)") "sigma"
         DO i = 1, 3
            WRITE (dump_unit, "(3(1X,ES24.16))") tb%sigma(i, :)/para_env%num_pe
         END DO
         CLOSE (dump_unit)
      END IF

      virial%pv_virial = virial%pv_virial - tb%sigma/para_env%num_pe

   END SUBROUTINE tb_add_stress

! **************************************************************************************************
!> \brief add contrib. to gradient
!> \param grad ...
!> \param deriv ...
!> \param dE ...
!> \param natom ...
! **************************************************************************************************
   SUBROUTINE tb_add_grad(grad, deriv, dE, natom)

      REAL(KIND=dp), DIMENSION(:, :)                     :: grad
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: deriv
      REAL(KIND=dp), DIMENSION(:)                        :: dE
      INTEGER                                            :: natom

      INTEGER                                            :: i, j

      DO i = 1, natom
         DO j = 1, natom
            grad(:, i) = grad(:, i) + deriv(:, i, j)*dE(j)
         END DO
      END DO

   END SUBROUTINE tb_add_grad

! **************************************************************************************************
!> \brief add contrib. to sigma
!> \param sig ...
!> \param deriv ...
!> \param dE ...
!> \param natom ...
! **************************************************************************************************
   SUBROUTINE tb_add_sig(sig, deriv, dE, natom)

      REAL(KIND=dp), DIMENSION(:, :)                     :: sig
      REAL(KIND=dp), DIMENSION(:, :, :)                  :: deriv
      REAL(KIND=dp), DIMENSION(:)                        :: dE
      INTEGER                                            :: natom

      INTEGER                                            :: i, j

      DO i = 1, 3
         DO j = 1, natom
            sig(:, i) = sig(:, i) + deriv(:, i, j)*dE(j)
         END DO
      END DO

   END SUBROUTINE tb_add_sig

END MODULE tblite_interface
