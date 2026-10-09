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
   USE mctc_io, ONLY: new
   USE mctc_io_symbols, ONLY: symbol_to_number
   USE tblite_basis_type, ONLY: get_cutoff
   USE tblite_container, ONLY: container_cache
   USE tblite_container_type, ONLY: container_type
   USE tblite_data_spin, ONLY: get_spin_constant
   USE tblite_integral_multipole, ONLY: multipole_cgto, multipole_grad_cgto, msao
   USE tblite_scf, ONLY: get_mixer_dimension
   USE tblite_scf_info, ONLY: scf_info, atom_resolved, shell_resolved
   USE tblite_scf_potential, ONLY: new_potential
   USE tblite_spin, ONLY: spin_polarization, new_spin_polarization
   USE tblite_wavefunction_type, ONLY: new_wavefunction
   USE tblite_xtb_calculator, ONLY: new_xtb_calculator
   USE tblite_xtb_gfn1, ONLY: new_gfn1_calculator
   USE tblite_xtb_gfn2, ONLY: new_gfn2_calculator
   USE tblite_xtb_h0, ONLY: get_selfenergy, get_occupation, tb_hamiltonian
   USE tblite_xtb_ipea1, ONLY: new_ipea1_calculator
#endif
   USE ai_contraction, ONLY: block_add, contraction
   USE ai_overlap, ONLY: overlap_ab
   USE atomic_kind_types, ONLY: atomic_kind_type, get_atomic_kind, get_atomic_kind_set
   USE atprop_types, ONLY: atprop_type
   USE basis_set_types, ONLY: gto_basis_set_type, gto_basis_set_p_type, allocate_gto_basis_set, &
                              process_gto_basis
   USE block_p_types, ONLY: block_p_type
   USE cell_types, ONLY: cell_type, get_cell
   USE cp_blacs_env, ONLY: cp_blacs_env_type
   USE cp_control_types, ONLY: dft_control_type
   USE cp_dbcsr_api, ONLY: dbcsr_type, dbcsr_p_type, dbcsr_create, dbcsr_add, dbcsr_copy, &
                           dbcsr_get_block_p, dbcsr_finalize, dbcsr_iterator_type, &
                           dbcsr_iterator_blocks_left, dbcsr_iterator_start, dbcsr_iterator_stop, &
                           dbcsr_iterator_next_block
   USE cp_dbcsr_contrib, ONLY: dbcsr_dot
   USE cp_dbcsr_cp2k_link, ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations, ONLY: dbcsr_allocate_matrix_set, dbcsr_deallocate_matrix_set
   USE cp_log_handling, ONLY: cp_get_default_logger, cp_logger_type, cp_logger_get_default_io_unit
   USE cp_output_handling, ONLY: cp_print_key_should_output, cp_print_key_unit_nr, &
                                 cp_print_key_finished_output
   USE input_constants, ONLY: gfn1xtb, gfn2xtb, ipea1xtb, tblite_scc_mixer_auto, tblite_scc_mixer_cp2k, &
                              tblite_scc_mixer_none, tblite_scc_mixer_tblite, tblite_solver_gvd, &
                              tblite_solver_gvr
   USE input_section_types, ONLY: section_vals_val_get, section_vals_get_subs_vals, section_vals_type
   USE kinds, ONLY: default_path_length, default_string_length, dp
   USE kpoint_types, ONLY: get_kpoint_info, kpoint_type
   USE memory_utilities, ONLY: reallocate
   USE message_passing, ONLY: mp_para_env_type
   USE orbital_pointers, ONLY: ncoset
   USE particle_types, ONLY: particle_type
   USE qs_charge_mixing, ONLY: charge_mixing, tblite_scc_error_on_cp2k_scale, tblite_scc_pconv
   USE qs_condnum, ONLY: overlap_condnum
   USE qs_density_mixing_types, ONLY: modified_broyden_mixing_nr
   USE qs_energy_types, ONLY: qs_energy_type
   USE qs_environment_types, ONLY: get_qs_env, qs_environment_type
   USE qs_force_types, ONLY: qs_force_type
   USE qs_integral_utils, ONLY: basis_set_list_setup, get_memory_usage
   USE qs_kind_types, ONLY: get_qs_kind, qs_kind_type, get_qs_kind_set
   USE qs_ks_types, ONLY: get_ks_env, qs_ks_env_type, set_ks_env
   USE qs_neighbor_list_types, ONLY: neighbor_list_iterator_create, neighbor_list_iterate, &
                                     get_iterator_info, neighbor_list_set_p_type, &
                                     neighbor_list_iterator_p_type, neighbor_list_iterator_release
   USE qs_overlap, ONLY: build_overlap_matrix, create_sab_matrix
   USE qs_rho_types, ONLY: qs_rho_get, qs_rho_type
   USE qs_scf_types, ONLY: qs_scf_env_type
   USE scf_control_types, ONLY: scf_control_type
   USE string_utilities, ONLY: integer_to_string
#if defined(__TBLITE)
   USE tblite_scc_mixer, ONLY: new_cp2k_tblite_mixer
#endif
   USE tblite_types, ONLY: tblite_type, allocate_tblite_type
   USE virial_types, ONLY: virial_type
   USE xtb_types, ONLY: get_xtb_atom_param, xtb_atom_type

!$ USE OMP_LIB, ONLY: omp_destroy_lock, omp_init_lock, omp_set_lock, &
!$                    omp_unset_lock, omp_lock_kind

#include "../base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tblite_interface'

   INTEGER, PARAMETER                                 :: dip_n = 3
   INTEGER, PARAMETER                                 :: quad_n = 6
   REAL(KIND=dp), PARAMETER                           :: same_atom = 0.00001_dp

   PUBLIC :: tb_set_calculator, tb_init_geometry, tb_init_wf
   PUBLIC :: tb_get_basis, build_tblite_matrices
   PUBLIC :: tb_get_energy, tb_update_charges, tb_ham_add_coulomb
   PUBLIC :: tb_native_scc_mixer_active
   PUBLIC :: tb_scf_mixer_error
   PUBLIC :: tb_get_multipole
   PUBLIC :: tb_derive_dH_off

CONTAINS

#if defined(__TBLITE)
! **************************************************************************************************
!> \brief Project a tblite charge/magnetization quantity to a CP2K spin channel.
!> \param values ...
!> \param ispin ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION tb_spin_project(values, ispin) RESULT(value)

      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: values
      INTEGER, INTENT(IN)                                :: ispin
      REAL(KIND=dp)                                      :: value

      value = values(1)
      IF (SIZE(values) > 1) THEN
         SELECT CASE (ispin)
         CASE (1)
            value = values(1) + values(2)
         CASE (2)
            value = values(1) - values(2)
         CASE DEFAULT
            value = values(1)
         END SELECT
      END IF

   END FUNCTION tb_spin_project

! **************************************************************************************************
!> \brief Preserve spin-resolved density before force assembly combines spin channels.
!> \param tb ...
!> \param matrix_p ...
! **************************************************************************************************
   SUBROUTINE tb_store_density_ref(tb, matrix_p)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p

      INTEGER                                            :: img, ispin, nimg, nspin

      nspin = SIZE(matrix_p, 1)
      nimg = SIZE(matrix_p, 2)
      IF (ASSOCIATED(tb%rho_ao_kp_ref)) THEN
         IF (SIZE(tb%rho_ao_kp_ref, 1) /= nspin .OR. SIZE(tb%rho_ao_kp_ref, 2) /= nimg) THEN
            CALL dbcsr_deallocate_matrix_set(tb%rho_ao_kp_ref)
         END IF
      END IF
      IF (.NOT. ASSOCIATED(tb%rho_ao_kp_ref)) THEN
         CALL dbcsr_allocate_matrix_set(tb%rho_ao_kp_ref, nspin, nimg)
         DO img = 1, nimg
            DO ispin = 1, nspin
               ALLOCATE (tb%rho_ao_kp_ref(ispin, img)%matrix)
            END DO
         END DO
      END IF
      DO img = 1, nimg
         DO ispin = 1, nspin
            CALL dbcsr_copy(tb%rho_ao_kp_ref(ispin, img)%matrix, matrix_p(ispin, img)%matrix)
         END DO
      END DO

   END SUBROUTINE tb_store_density_ref
#endif

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
      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      !get info from environment vaiarable
      NULLIFY (cell)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell)

      !get information about particles
      natom = SIZE(particle_set)
      DO iatom = 1, natom
         tb%mol%xyz(:, iatom) = particle_set(iatom)%r(:)
      END DO
      tb%mol%lattice(:, :) = cell%hmat

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
!> \param dft_control ...
! **************************************************************************************************
   SUBROUTINE tb_init_wf(tb, dft_control)

      TYPE(tblite_type), POINTER                         :: tb
      TYPE(dft_control_type), POINTER                    :: dft_control

#if defined(__TBLITE)

      INTEGER                                            :: nSpin

      TYPE(scf_info)                                     :: info

      nSpin = dft_control%nspins
      IF (nSpin /= 1 .AND. nSpin /= 2) CPABORT("tblite supports only one or two spin channels")

      tb%mol%charge = dft_control%charge
      tb%mol%uhf = MAX(0, dft_control%multiplicity - 1)
      IF (nSpin == 2) CALL tb_add_spin_polarization(tb)

      info = tb%calc%variable_info()
      IF (info%charge > shell_resolved) CPABORT("tblite: no support for orbital resolved charge")
      IF (info%dipole > atom_resolved) CPABORT("tblite: no support for shell resolved dipole moment")
      IF (info%quadrupole > atom_resolved) THEN
         CPABORT("tblite: no support shell resolved quadrupole moment")
      END IF

      CALL new_wavefunction(tb%wfn, tb%mol%nat, tb%calc%bas%nsh, tb%calc%bas%nao, nSpin, 0.0_dp)
      CALL get_occupation(tb%mol, tb%calc%bas, tb%calc%h0, tb%wfn%nocc, tb%wfn%n0at, tb%wfn%n0sh)
      CALL tb_reset_mixer(tb)

      CALL new_potential(tb%pot, tb%mol, tb%calc%bas, tb%wfn%nspin)

      !allocate quantities later required
      ALLOCATE (tb%e_hal(tb%mol%nat), tb%e_rep(tb%mol%nat), tb%e_disp(tb%mol%nat))
      ALLOCATE (tb%e_scd(tb%mol%nat), tb%e_es(tb%mol%nat), tb%e_int(tb%mol%nat))
      ALLOCATE (tb%selfenergy(tb%calc%bas%nsh))
      IF (ALLOCATED(tb%calc%ncoord)) ALLOCATE (tb%cn(tb%mol%nat))

#else
      MARK_USED(tb)
      MARK_USED(dft_control)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_init_wf

#if defined(__TBLITE)
! **************************************************************************************************
!> \brief Add tblite's on-site spin-polarization interaction.
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_add_spin_polarization(tb)

      TYPE(tblite_type), POINTER                         :: tb

      CLASS(container_type), ALLOCATABLE                 :: cont
      TYPE(spin_polarization), ALLOCATABLE               :: spin
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: wll

      ALLOCATE (spin)
      CALL tb_get_spin_constants(tb, wll)
      CALL new_spin_polarization(spin, tb%mol, wll, tb%calc%bas%nsh_id)
      CALL MOVE_ALLOC(spin, cont)
      CALL tb%calc%push_back(cont)

   END SUBROUTINE tb_add_spin_polarization

! **************************************************************************************************
!> \brief Build tblite spin constants for the current basis.
!> \param tb ...
!> \param wll ...
! **************************************************************************************************
   SUBROUTINE tb_get_spin_constants(tb, wll)

      TYPE(tblite_type), POINTER                         :: tb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: wll

      INTEGER                                            :: il, ish, izp, jl, jsh, max_shell

      max_shell = MAXVAL(tb%calc%bas%nsh_id)
      ALLOCATE (wll(max_shell, max_shell, tb%mol%nid))
      wll = 0.0_dp
      DO izp = 1, tb%mol%nid
         DO ish = 1, tb%calc%bas%nsh_id(izp)
            il = tb%calc%bas%cgto(ish, izp)%ang
            DO jsh = 1, tb%calc%bas%nsh_id(izp)
               jl = tb%calc%bas%cgto(jsh, izp)%ang
               wll(jsh, ish, izp) = get_spin_constant(jl, il, tb%mol%num(izp))
            END DO
         END DO
      END DO

   END SUBROUTINE tb_get_spin_constants
#endif

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
      CALL new_cp2k_tblite_mixer(tb%mixer, tb%mixer_memory, &
                                 tb%wfn%nspin*get_mixer_dimension(tb%mol, tb%calc%bas, info), &
                                 tb%mixer_damping, tb%mixer_omega0, tb%mixer_min_weight, &
                                 tb%mixer_max_weight, tb%mixer_weight_factor)

#else
      MARK_USED(tb)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_reset_mixer

! **************************************************************************************************
!> \brief Configure tblite's internal SCC mixer from CP2K input.
!> \param tb ...
!> \param iterations ...
!> \param memory ...
!> \param damping ...
!> \param omega0 ...
!> \param min_weight ...
!> \param max_weight ...
!> \param weight_factor ...
!> \param solver ...
! **************************************************************************************************
   SUBROUTINE tb_configure_mixer(tb, iterations, memory, damping, omega0, min_weight, max_weight, &
                                 weight_factor, solver)

      TYPE(tblite_type), POINTER                         :: tb
      INTEGER, INTENT(IN)                                :: iterations, memory, solver
      REAL(KIND=dp), INTENT(IN)                          :: damping, max_weight, min_weight, omega0, &
                                                            weight_factor

#if defined(__TBLITE)

      IF (iterations < 1) CPABORT("tblite SCC mixer ITERATIONS must be positive")
      IF (memory < 1) CPABORT("tblite SCC mixer MEMORY must be positive")
      IF (damping <= 0.0_dp) CPABORT("tblite SCC mixer damping must be positive")
      IF (omega0 <= 0.0_dp) CPABORT("tblite SCC mixer OMEGA0 must be positive")
      IF (min_weight <= 0.0_dp) CPABORT("tblite SCC mixer MIN_WEIGHT must be positive")
      IF (max_weight <= 0.0_dp) CPABORT("tblite SCC mixer MAX_WEIGHT must be positive")
      IF (max_weight < min_weight) THEN
         CPABORT("tblite SCC mixer MAX_WEIGHT must not be smaller than MIN_WEIGHT")
      END IF
      IF (weight_factor <= 0.0_dp) CPABORT("tblite SCC mixer WEIGHT_FACTOR must be positive")
      SELECT CASE (solver)
      CASE (tblite_solver_gvd, tblite_solver_gvr)
      CASE DEFAULT
         CPABORT("Unknown tblite SCC mixer SOLVER")
      END SELECT

      tb%calc%max_iter = iterations
      tb%mixer_memory = memory
      tb%mixer_damping = damping
      tb%calc%mixer_input%damping = damping
      tb%mixer_omega0 = omega0
      tb%mixer_min_weight = min_weight
      tb%mixer_max_weight = max_weight
      tb%mixer_weight_factor = weight_factor

#else
      MARK_USED(tb)
      MARK_USED(iterations)
      MARK_USED(memory)
      MARK_USED(damping)
      MARK_USED(omega0)
      MARK_USED(min_weight)
      MARK_USED(max_weight)
      MARK_USED(weight_factor)
      MARK_USED(solver)
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
      IF (dft_control%qs_control%do_ls_scf) RETURN

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
!> \brief Return the native tblite SCC mixer residual on the CP2K iter_delta scale.
!> \param dft_control ...
!> \param tb ...
!> \param eps_scf CP2K reporting scale for the native residual.
!> \return ...
! **************************************************************************************************
   FUNCTION tb_scf_mixer_error(dft_control, tb, eps_scf) RESULT(mixer_error)

      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(tblite_type), POINTER                         :: tb
      REAL(KIND=dp), INTENT(IN)                          :: eps_scf
      REAL(KIND=dp)                                      :: mixer_error

#if defined(__TBLITE)
      REAL(KIND=dp)                                      :: raw_error
#endif

      mixer_error = 0.0_dp

#if defined(__TBLITE)
      IF (.NOT. ASSOCIATED(tb)) RETURN
      IF (.NOT. tb_native_scc_mixer_active(dft_control)) RETURN
      IF (ALLOCATED(tb%mixer)) THEN
         raw_error = REAL(tb%mixer%get_error(), KIND=dp)
         mixer_error = tblite_scc_error_on_cp2k_scale( &
                       raw_error, eps_scf, &
                       tblite_scc_pconv*dft_control%qs_control%xtb_control%tblite_accuracy)
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
!> \param accuracy ...
!> \param param_file ...
! **************************************************************************************************
   SUBROUTINE tb_set_calculator(tb, typ, accuracy, param_file)

      TYPE(tblite_type), POINTER                         :: tb
      INTEGER                                            :: typ
      REAL(KIND=dp), INTENT(IN)                          :: accuracy
      CHARACTER(LEN=*), INTENT(IN)                       :: param_file

#if defined(__TBLITE)

      TYPE(error_type), ALLOCATABLE                      :: error

      IF (ALLOCATED(tb%param)) DEALLOCATE (tb%param)
      IF (LEN_TRIM(param_file) > 0) THEN
         ALLOCATE (tb%param)
         CALL tb%param%load(TRIM(param_file), error)
         IF (ALLOCATED(error)) CPABORT("Could not load tblite PARAM file: "//TRIM(param_file))
         CALL new_xtb_calculator(tb%calc, tb%mol, tb%param, error)
      ELSE
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
      END IF
      IF (ALLOCATED(error)) CPABORT("tblite calculator setup failed")

      tb%accuracy = accuracy

#else
      MARK_USED(tb)
      MARK_USED(typ)
      MARK_USED(accuracy)
      MARK_USED(param_file)
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
      tb%e_int = 0.0_dp
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

      IF (ALLOCATED(tb%calc%interactions)) THEN
         CALL tb%calc%interactions%update(tb%mol, tb%icache)
      END IF

      CALL new_potential(tb%pot, tb%mol, tb%calc%bas, tb%wfn%nspin)
      IF (ALLOCATED(tb%calc%coulomb)) THEN
         CALL tb%calc%coulomb%update(tb%mol, tb%cache)
      END IF

      IF (ALLOCATED(tb%grad)) THEN
         IF (ALLOCATED(tb%calc%ncoord)) THEN
            CALL tb%calc%ncoord%get_cn(tb%mol, tb%cn, tb%dcndr, tb%dcndL)
         END IF
      ELSE
         IF (ALLOCATED(tb%calc%ncoord)) THEN
            CALL tb%calc%ncoord%get_cn(tb%mol, tb%cn)
         END IF
      END IF

      CALL get_selfenergy(tb%calc%h0, tb%mol%id, tb%calc%bas%ish_at, &
                          tb%calc%bas%nsh_id, cn=tb%cn, selfenergy=tb%selfenergy, dsedcn=tb%dsedcn)

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
      REAL(KIND=dp)                                      :: xtb_inter
      NULLIFY (scf_section, logger)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)
      scf_section => section_vals_get_subs_vals(qs_env%input, "DFT%SCF")

      energy%repulsive = SUM(tb%e_rep)
      energy%el_stat = SUM(tb%e_es)
      energy%dispersion = SUM(tb%e_disp)
      energy%dispersion_sc = SUM(tb%e_scd)
      energy%xtb_xb_inter = SUM(tb%e_hal)
      xtb_inter = SUM(tb%e_int)

      energy%total = energy%core + energy%repulsive + energy%el_stat + energy%dispersion &
                     + energy%dispersion_sc + energy%xtb_xb_inter + xtb_inter &
                     + energy%kTS + energy%efield + energy%qmmm_el

      iounit = cp_print_key_unit_nr(logger, scf_section, "PRINT%DETAILED_ENERGY", &
                                    extension=".scfLog")
      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT="(/,(T9,A,T60,F20.10))") &
            "Repulsive pair potential energy:               ", energy%repulsive, &
            "Zeroth order Hamiltonian energy:               ", energy%core, &
            "Electrostatic energy:                          ", energy%el_stat, &
            "Self-consistent dispersion energy:             ", energy%dispersion_sc, &
            "Non-self consistent dispersion energy:         ", energy%dispersion
         IF (ABS(energy%xtb_xb_inter) > 1.e-9_dp) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "Correction for halogen bonding:                ", energy%xtb_xb_inter
         END IF
         IF (ABS(xtb_inter) > 1.e-9_dp) THEN
            WRITE (UNIT=iounit, FMT="(T9,A,T60,F20.10)") &
               "Additional interaction (e.g. spin):            ", xtb_inter
         END IF
         IF (ABS(energy%efield) > 1.e-9_dp) THEN
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
      param%rcut = get_cutoff(tb%calc%bas, tb%accuracy)
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

      INTEGER                                            :: handle, maxder, nderivatives, nimg, img, nkind, i, &
                                                            ic, iw, iatom, jatom, ikind, jkind, iset, jset, n1, n2, icol, &
                                                            irow, ia, ib, sgfa, sgfb, ldsab, nseta, nsetb, &
                                                            natorb_a, natorb_b, raw_iatom, raw_jatom, &
                                                            slot, idim, jdim
      LOGICAL                                            :: found, norml1, norml2, use_arnoldi
      REAL(KIND=dp)                                      :: dr, dshpoly, ff, hij_base, r2, rr
      REAL(KIND=dp)                                      :: native_dot_tmp, native_cn_icol, native_cn_irow, &
                                                            native_dot_weight
      INTEGER, DIMENSION(3)                              :: cell
      REAL(KIND=dp)                                      :: hij, shpoly
      REAL(KIND=dp), DIMENSION(2)                        :: condnum
      REAL(KIND=dp), DIMENSION(3)                        :: native_h0_overlap_force, native_radial_force, raw_rij, rij
      REAL(KIND=dp), DIMENSION(3, 3)                     :: native_radial_dot
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of
      INTEGER, DIMENSION(:, :, :), POINTER                :: cell_to_index
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: owork
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: native_cn_deriv_thread
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: native_radial_force_thread
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: oint, sint, hint
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: radial_hint
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min
      INTEGER, DIMENSION(:), POINTER                     :: npgfa, npgfb, nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgfa, rpgfb, zeta, zetb, scon_a, scon_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: fblock, pblock, sblock
      TYPE(block_p_type), DIMENSION(3, 3)                 :: radial_blocks
!$    INTEGER                                            :: hash, hash1, lock_num
!$    INTEGER(KIND=int_8)                                :: iatom8
!$    INTEGER(kind=omp_lock_kind), ALLOCATABLE, DIMENSION(:) :: locks
      INTEGER, PARAMETER                                :: nlock = 501

      TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
      TYPE(atprop_type), POINTER                            :: atprop
      TYPE(cp_blacs_env_type), POINTER                      :: blacs_env
      TYPE(cp_logger_type), POINTER                         :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER          :: matrix_dh_native, matrix_h, matrix_p, &
                                                               matrix_q_native, matrix_s, matrix_s_native, matrix_w
      TYPE(dft_control_type), POINTER                       :: dft_control
      TYPE(qs_force_type), DIMENSION(:), POINTER            :: force
      TYPE(gto_basis_set_type), POINTER                     :: basis_set_a, basis_set_b
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER     :: basis_set_list
      TYPE(kpoint_type), POINTER                            :: kpoints
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_kp
      TYPE(mp_para_env_type), POINTER                       :: para_env
      TYPE(qs_energy_type), POINTER                         :: energy
      TYPE(qs_ks_env_type), POINTER                         :: ks_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER             :: qs_kind_set
      TYPE(qs_rho_type), POINTER                            :: rho
      TYPE(tblite_type), POINTER                            :: tb
      TYPE(tb_hamiltonian), POINTER                         :: h0
      TYPE(virial_type), POINTER                            :: virial

      CALL timeset(routineN, handle)

      NULLIFY (ks_env, energy, atomic_kind_set, qs_kind_set)
      NULLIFY (matrix_dh_native, matrix_h, matrix_q_native, matrix_s, matrix_s_native, atprop, dft_control)
      NULLIFY (sab_orb, sab_kp, rho, tb, kpoints, cell_to_index)

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
                      sab_kp=sab_kp, &
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
         IF (.NOT. ASSOCIATED(sab_kp)) CPABORT("Missing k-point neighbor list for tblite")
         sab_orb => sab_kp
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

      IF (calculate_forces) THEN
         CALL dbcsr_allocate_matrix_set(matrix_q_native, 1, nimg)
         DO img = 1, nimg
            ALLOCATE (matrix_q_native(1, img)%matrix)
            CALL dbcsr_copy(matrix_q_native(1, img)%matrix, matrix_w(1, img)%matrix, &
                            name="TBLITE NATIVE OVERLAP FORCE MATRIX")
         END DO
      END IF

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
      IF (calculate_forces) THEN
         ALLOCATE (native_radial_force_thread(3, SIZE(atom_of_kind)))
         ALLOCATE (native_cn_deriv_thread(SIZE(atom_of_kind)))
         native_cn_deriv_thread = 0.0_dp
         native_radial_force_thread = 0.0_dp
      END IF

      ! set up basis set lists
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      ! allocate overlap matrix
      CALL dbcsr_allocate_matrix_set(matrix_s, 1, nimg)
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
      IF (calculate_forces .AND. nimg > 1) THEN
         CALL dbcsr_allocate_matrix_set(matrix_dh_native, 9, nimg)
         DO img = 1, nimg
            DO idim = 1, 9
               ALLOCATE (matrix_dh_native(idim, img)%matrix)
               CALL dbcsr_create(matrix_dh_native(idim, img)%matrix, template=matrix_s(1, img)%matrix, &
                                 name="TBLITE H0 STRAIN DERIVATIVE")
               CALL cp_dbcsr_alloc_block_from_nbl(matrix_dh_native(idim, img)%matrix, sab_orb)
            END DO
         END DO
      END IF
      ldsab = get_memory_usage(qs_kind_set, "ORB", "ORB")

      native_radial_dot = 0.0_dp

      ! loop over all atom pairs with a non-zero overlap (sab_orb)
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED (calculate_forces, basis_set_list, sab_orb, matrix_dh_native, matrix_s, matrix_h, matrix_p, &
!$OMP        cell_to_index, nimg, atom_of_kind, qs_kind_set, tb, h0, ldsab, maxder, locks, &
!$OMP        native_cn_deriv_thread, native_radial_force_thread, ncoset) &
!$OMP PRIVATE (slot, hash, hash1, iatom8, lock_num, ikind, jkind, iatom, jatom, cell, rij, ic, irow, &
!$OMP        icol, dr, dshpoly, ff, hij_base, r2, raw_iatom, raw_jatom, &
!$OMP        raw_rij, native_cn_icol, native_cn_irow, native_dot_weight, native_h0_overlap_force, &
!$OMP        native_radial_force, &
!$OMP        n1, n2, ia, ib, i, idim, jdim, iset, jset, sgfa, sgfb, &
!$OMP        nseta, nsetb, natorb_a, &
!$OMP        natorb_b, found, zeta, first_sgfa, la_max, la_min, npgfa, nsgfa, rpgfa, set_radius_a, &
!$OMP        scon_a, first_sgfb, lb_max, lb_min, npgfb, nsgfb, rpgfb, set_radius_b, scon_b, zetb, rr, hij, shpoly, &
!$OMP        owork, oint, sint, hint, radial_hint, basis_set_a, basis_set_b, sblock, fblock, &
!$OMP        pblock, radial_blocks) &
!$OMP REDUCTION (+ : native_radial_dot)

!$OMP SINGLE
!$    ALLOCATE (locks(nlock))
!$OMP END SINGLE
!$OMP DO
!$    DO lock_num = 1, nlock
!$       CALL omp_init_lock(locks(lock_num))
!$    END DO
!$OMP END DO

      ALLOCATE (oint(ldsab, ldsab, maxder), owork(ldsab, ldsab))

!$OMP DO SCHEDULE(GUIDED)
      DO slot = 1, sab_orb(1)%nl_size

         ikind = sab_orb(1)%nlist_task(slot)%ikind
         jkind = sab_orb(1)%nlist_task(slot)%jkind
         iatom = sab_orb(1)%nlist_task(slot)%iatom
         jatom = sab_orb(1)%nlist_task(slot)%jatom
         cell(:) = sab_orb(1)%nlist_task(slot)%cell(:)
         rij(1:3) = sab_orb(1)%nlist_task(slot)%r(1:3)
         raw_iatom = iatom
         raw_jatom = jatom
         raw_rij = rij
         native_cn_icol = 0.0_dp
         native_cn_irow = 0.0_dp
         native_h0_overlap_force = 0.0_dp
         native_radial_force = 0.0_dp

         ! canonicalize pair ordering as in current serial flow
         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)
         IF (iatom < jatom) THEN
            rij = -rij
            i = ikind
            ikind = jkind
            jkind = i
         END IF

         dr = NORM2(rij(:))
         r2 = dr*dr

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
         IF (calculate_forces) THEN
            NULLIFY (pblock)
            CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                   row=irow, col=icol, BLOCK=pblock, found=found)
            CPASSERT(found)
            IF (nimg > 1) THEN
               DO jdim = 1, 3
                  DO idim = 1, 3
                     NULLIFY (radial_blocks(idim, jdim)%block)
                     CALL dbcsr_get_block_p(matrix=matrix_dh_native(idim + 3*(jdim - 1), ic)%matrix, &
                                            row=irow, col=icol, BLOCK=radial_blocks(idim, jdim)%block, &
                                            found=found)
                     CPASSERT(found)
                  END DO
               END DO
            END IF
         END IF

         ! --------- Overlap
         !get basis information
         basis_set_a => basis_set_list(ikind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
         basis_set_b => basis_set_list(jkind)%gto_basis_set
         IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE
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
         IF (calculate_forces .AND. nimg > 1) THEN
            ALLOCATE (radial_hint(natorb_a, natorb_b, 3, 3))
            radial_hint = 0.0_dp
         END IF

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

!$       iatom8 = INT(iatom - 1, int_8)*INT(SIZE(atom_of_kind), int_8) + INT(jatom, int_8)
!$       hash1 = INT(MOD(iatom8, INT(nlock, int_8)) + 1)
         ! update S matrix
!$       hash = hash1
!$       CALL omp_set_lock(locks(hash))
         IF (icol <= irow) THEN
            sblock(:, :) = sblock(:, :) + sint(:, :, 1)
         ELSE
            sblock(:, :) = sblock(:, :) + TRANSPOSE(sint(:, :, 1))
         END IF
!$       CALL omp_unset_lock(locks(hash))

         ! --------- Hamiltonian
         IF (icol == irow .AND. dr < same_atom) THEN
            !get diagonal F matrix from selfenergy
            n1 = tb%calc%bas%ish_at(icol)
            DO iset = 1, nseta
               sgfa = first_sgfa(1, iset)
               hij = tb%selfenergy(n1 + iset)
               DO ia = sgfa, sgfa + nsgfa(iset) - 1
                  hint(ia, ia, 1) = hij
                  IF (calculate_forces) THEN
                     native_cn_icol = native_cn_icol + tb%dsedcn(n1 + iset)*pblock(ia, ia)
                     DO jdim = 1, 3
                        DO idim = 1, 3
                           IF (nimg == 1) THEN
                              native_radial_dot(idim, jdim) = native_radial_dot(idim, jdim) + &
                                                              tb%dsedcn(n1 + iset)*tb%dcndL(idim, jdim, icol)*pblock(ia, ia)
                           ELSE
                              radial_hint(ia, ia, idim, jdim) = radial_hint(ia, ia, idim, jdim) + &
                                                                tb%dsedcn(n1 + iset)*tb%dcndL(idim, jdim, icol)
                           END IF
                        END DO
                     END DO
                  END IF
               END DO
            END DO
         ELSE
            !get off-diagonal F matrix
            rr = SQRT(dr/(h0%rad(jkind) + h0%rad(ikind)))
            n1 = tb%calc%bas%ish_at(icol)
            DO iset = 1, nseta
               sgfa = first_sgfa(1, iset)
               n2 = tb%calc%bas%ish_at(irow)
               DO jset = 1, nsetb
                  sgfb = first_sgfb(1, jset)
                  shpoly = (1.0_dp + h0%shpoly(iset, ikind)*rr) &
                           *(1.0_dp + h0%shpoly(jset, jkind)*rr)
                  dshpoly = ((1.0_dp + h0%shpoly(iset, ikind)*rr)*h0%shpoly(jset, jkind)*rr &
                             + (1.0_dp + h0%shpoly(jset, jkind)*rr)*h0%shpoly(iset, ikind)*rr) &
                            *0.5_dp/r2
                  hij_base = 0.5_dp*(tb%selfenergy(n1 + iset) + tb%selfenergy(n2 + jset)) &
                             *h0%hscale(iset, jset, ikind, jkind)
                  hij = hij_base*shpoly
                  DO ia = sgfa, sgfa + nsgfa(iset) - 1
                     DO ib = sgfb, sgfb + nsgfb(jset) - 1
                        hint(ia, ib, 1) = hij*sint(ia, ib, 1)
                        IF (calculate_forces) THEN
                           native_dot_weight = 2.0_dp
                           IF (icol == irow) native_dot_weight = 1.0_dp
                           native_cn_icol = native_cn_icol + native_dot_weight*0.5_dp* &
                                            h0%hscale(iset, jset, ikind, jkind)*shpoly* &
                                            tb%dsedcn(n1 + iset)*pblock(ib, ia)*sint(ia, ib, 1)
                           native_cn_irow = native_cn_irow + native_dot_weight*0.5_dp* &
                                            h0%hscale(iset, jset, ikind, jkind)*shpoly* &
                                            tb%dsedcn(n2 + jset)*pblock(ib, ia)*sint(ia, ib, 1)
                           native_radial_force = native_radial_force + &
                                                 hij_base*dshpoly*pblock(ib, ia)*sint(ia, ib, 1)*raw_rij
                           native_h0_overlap_force = native_h0_overlap_force + &
                                                     hij*pblock(ib, ia)*sint(ia, ib, 2:4)
                           DO jdim = 1, 3
                              DO idim = 1, 3
                                 IF (nimg == 1) THEN
                                    native_radial_dot(idim, jdim) = native_radial_dot(idim, jdim) + &
                                                                    native_dot_weight*pblock(ib, ia)*( &
                                                                    (-hij)*sint(ia, ib, idim + 1)*rij(jdim) + &
                                                                    hij_base*dshpoly*sint(ia, ib, 1)*rij(idim)*rij(jdim) + &
                                                                0.5_dp*h0%hscale(iset, jset, ikind, jkind)*shpoly*sint(ia, ib, 1)* &
                                                                    (tb%dsedcn(n1 + iset)*tb%dcndL(idim, jdim, icol) + &
                                                                     tb%dsedcn(n2 + jset)*tb%dcndL(idim, jdim, irow)))
                                 ELSE
                                    radial_hint(ia, ib, idim, jdim) = radial_hint(ia, ib, idim, jdim) + &
                                                                      (-hij)*sint(ia, ib, idim + 1)*rij(jdim) + &
                                                                      hij_base*dshpoly*sint(ia, ib, 1)*rij(idim)*rij(jdim) + &
                                                                0.5_dp*h0%hscale(iset, jset, ikind, jkind)*shpoly*sint(ia, ib, 1)* &
                                                                      (tb%dsedcn(n1 + iset)*tb%dcndL(idim, jdim, icol) + &
                                                                       tb%dsedcn(n2 + jset)*tb%dcndL(idim, jdim, irow))
                                 END IF
                              END DO
                           END DO
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END IF

         ! update F matrix
!$       CALL omp_set_lock(locks(hash))
         IF (icol <= irow) THEN
            fblock(:, :) = fblock(:, :) + hint(:, :, 1)
            IF (calculate_forces .AND. nimg > 1) THEN
               DO jdim = 1, 3
                  DO idim = 1, 3
                     radial_blocks(idim, jdim)%block(:, :) = radial_blocks(idim, jdim)%block(:, :) + &
                                                             radial_hint(:, :, idim, jdim)
                  END DO
               END DO
            END IF
         ELSE
            fblock(:, :) = fblock(:, :) + TRANSPOSE(hint(:, :, 1))
            IF (calculate_forces .AND. nimg > 1) THEN
               DO jdim = 1, 3
                  DO idim = 1, 3
                     radial_blocks(idim, jdim)%block(:, :) = radial_blocks(idim, jdim)%block(:, :) + &
                                                             TRANSPOSE(radial_hint(:, :, idim, jdim))
                  END DO
               END DO
            END IF
         END IF
!$       CALL omp_unset_lock(locks(hash))

         IF (calculate_forces) THEN
            ff = 2.0_dp
!$OMP CRITICAL(tblite_native_radial_force)
            native_radial_force_thread(:, raw_iatom) = &
               native_radial_force_thread(:, raw_iatom) - ff*native_radial_force
            native_radial_force_thread(:, raw_jatom) = &
               native_radial_force_thread(:, raw_jatom) + ff*native_radial_force
            native_radial_force_thread(:, icol) = &
               native_radial_force_thread(:, icol) + ff*native_h0_overlap_force
            native_radial_force_thread(:, irow) = &
               native_radial_force_thread(:, irow) - ff*native_h0_overlap_force
            native_cn_deriv_thread(icol) = native_cn_deriv_thread(icol) + native_cn_icol
            native_cn_deriv_thread(irow) = native_cn_deriv_thread(irow) + native_cn_irow
!$OMP END CRITICAL(tblite_native_radial_force)
         END IF

         DEALLOCATE (sint, hint)
         IF (ALLOCATED(radial_hint)) DEALLOCATE (radial_hint)

      END DO
!$OMP END DO
      DEALLOCATE (oint, owork)
!$OMP DO
!$    DO lock_num = 1, nlock
!$       CALL omp_destroy_lock(locks(lock_num))
!$    END DO
!$OMP END DO
!$OMP SINGLE
!$    DEALLOCATE (locks)
!$OMP END SINGLE NOWAIT
!$OMP END PARALLEL

      IF (calculate_forces) THEN
         IF (nimg > 1) THEN
            native_radial_dot = 0.0_dp
            DO img = 1, nimg
               DO jdim = 1, 3
                  DO idim = 1, 3
                     CALL dbcsr_finalize(matrix_dh_native(idim + 3*(jdim - 1), img)%matrix)
                     CALL dbcsr_dot(matrix_dh_native(idim + 3*(jdim - 1), img)%matrix, &
                                    matrix_p(1, img)%matrix, native_dot_tmp)
                     native_radial_dot(idim, jdim) = native_radial_dot(idim, jdim) + native_dot_tmp
                  END DO
               END DO
            END DO
            CALL dbcsr_deallocate_matrix_set(matrix_dh_native)
         ELSE
            CALL para_env%sum(native_radial_dot)
         END IF
         CALL para_env%sum(native_cn_deriv_thread)
         tb%grad = 0.0_dp
         CALL tb_add_grad(tb%grad, tb%dcndr, native_cn_deriv_thread, tb%mol%nat)
         CALL tb_grad2force(qs_env, tb, para_env, 4)
         DO iatom = 1, SIZE(atom_of_kind)
            ikind = kind_of(iatom)
            force(ikind)%overlap(:, atom_of_kind(iatom)) = &
               force(ikind)%overlap(:, atom_of_kind(iatom)) + native_radial_force_thread(:, iatom)
         END DO
         IF (tb%use_virial) THEN
            virial%pv_overlap = virial%pv_overlap - native_radial_dot/para_env%num_pe
            virial%pv_virial = virial%pv_virial - native_radial_dot/para_env%num_pe
         END IF
      END IF

      DO img = 1, nimg
         DO i = 1, SIZE(matrix_s, 1)
            CALL dbcsr_finalize(matrix_s(i, img)%matrix)
         END DO
         DO i = 1, SIZE(matrix_h, 1)
            CALL dbcsr_finalize(matrix_h(i, img)%matrix)
         END DO
      END DO

      IF (calculate_forces) THEN
         CALL build_overlap_matrix(ks_env, matrixkp_s=matrix_s_native, &
                                   matrix_name="OVERLAP MATRIX", nderivative=nderivatives, &
                                   basis_type_a="ORB", basis_type_b="ORB", sab_nl=sab_orb, &
                                   calculate_forces=.TRUE., matrixkp_p=matrix_q_native)
         CALL dbcsr_deallocate_matrix_set(matrix_s_native)
         CALL dbcsr_deallocate_matrix_set(matrix_q_native)
      END IF

      !compute multipole moments for gfn2
      IF (dft_control%qs_control%xtb_control%tblite_method == gfn2xtb) THEN
         CALL tb_get_multipole(qs_env, tb)
      END IF

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
      IF (ALLOCATED(native_cn_deriv_thread)) DEALLOCATE (native_cn_deriv_thread)
      IF (ALLOCATED(native_radial_force_thread)) DEALLOCATE (native_radial_force_thread)

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

      INTEGER                                            :: iatom, ikind, is, ns, atom_a, im
      INTEGER                                            :: ispin, nspin
      INTEGER                                            :: nimg, nkind, nsgf, natorb, na
      INTEGER                                            :: n_atom, max_orb, max_shell
      INTEGER                                            :: raw_state_status, raw_state_unit
      LOGICAL                                            :: advance_native_mixer, do_dipole, do_quadrupole, &
                                                            reuse_native_input, skip_scf_dispersion, &
                                                            seed_native_from_rho, use_native_mixer, use_no_mixer
      REAL(KIND=dp)                                      :: native_seed_charge, norm, pao
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      INTEGER                                            :: debug_status
      CHARACTER(LEN=32)                                  :: debug_value
#endif
      CHARACTER(LEN=default_path_length)                 :: raw_state_file
      INTEGER, DIMENSION(5)                              :: occ
      INTEGER, DIMENSION(25)                             :: lao
      INTEGER, DIMENSION(25)                             :: nao
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ch_atom
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aocg
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: aocg_spin, ao_dip_spin, ao_quad_spin, &
                                                            ch_orb_spin, ch_shell_spin

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s, matrix_p
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: p_matrix
      TYPE(dbcsr_type), POINTER                          :: s_matrix
      TYPE(error_type), ALLOCATABLE                      :: error
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      ! compute mulliken charges required for charge update
      NULLIFY (particle_set, qs_kind_set, atomic_kind_set, scf_control, scf_env)
      CALL get_qs_env(qs_env=qs_env, scf_env=scf_env, particle_set=particle_set, qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set, matrix_s_kp=matrix_s, rho=rho, para_env=para_env, &
                      scf_control=scf_control)

      ! also compute multipoles needed by GFN2
      do_dipole = .FALSE.
      do_quadrupole = .FALSE.
      skip_scf_dispersion = .FALSE.
#if defined(__TBLITE_DEBUG_DIAGNOSTICS)
      CALL GET_ENVIRONMENT_VARIABLE("CP2K_TBLITE_DEBUG_SKIP_SCF_DISPERSION", debug_value, STATUS=debug_status)
      IF (debug_status == 0) THEN
         READ (debug_value, *, IOSTAT=debug_status) skip_scf_dispersion
         IF (debug_status /= 0) skip_scf_dispersion = .FALSE.
      END IF
#endif
      IF (dft_control%qs_control%do_ls_scf .OR. scf_control%use_ot) THEN
         use_native_mixer = .FALSE.
         use_no_mixer = .TRUE.
      ELSE
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
      END IF
      IF (use_native_mixer) THEN
         IF (.NOT. ASSOCIATED(scf_env)) CPABORT("tblite SCC mixer requires a QS SCF environment")
         IF (use_rho .AND. (.NOT. calculate_forces)) THEN
            IF (scf_env%iter_count > dft_control%qs_control%xtb_control%tblite_mixer_iterations) THEN
               CPABORT("tblite SCC mixer exceeded TBLITE_MIXER/ITERATIONS")
            END IF
            IF (scf_env%iter_count == 1) THEN
               CALL tb_configure_mixer(tb, dft_control%qs_control%xtb_control%tblite_mixer_iterations, &
                                       dft_control%qs_control%xtb_control%tblite_mixer_memory, &
                                       dft_control%qs_control%xtb_control%tblite_mixer_damping, &
                                       dft_control%qs_control%xtb_control%tblite_mixer_omega0, &
                                       dft_control%qs_control%xtb_control%tblite_mixer_min_weight, &
                                       dft_control%qs_control%xtb_control%tblite_mixer_max_weight, &
                                       dft_control%qs_control%xtb_control%tblite_mixer_weight_factor, &
                                       dft_control%qs_control%xtb_control%tblite_mixer_solver)
               CALL tb_reset_mixer(tb)
            END IF
         END IF
      END IF
      nspin = dft_control%nspins
      IF (nspin /= tb%wfn%nspin) CPABORT("CP2K/tblite spin channel mismatch")

      NULLIFY (matrix_p)
      IF (use_rho) THEN
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      ELSE IF (calculate_forces .AND. nspin > 1) THEN
         IF (.NOT. ASSOCIATED(tb%rho_ao_kp_ref)) THEN
            CPABORT("Missing converged tblite density for UKS/LSD forces")
         END IF
         matrix_p => tb%rho_ao_kp_ref
      ELSE
         matrix_p => scf_env%p_mix_new
      END IF
      IF (nspin > 1 .AND. (.NOT. calculate_forces)) CALL tb_store_density_ref(tb, matrix_p)
      IF (ASSOCIATED(tb%dipbra)) do_dipole = .TRUE.
      IF (ASSOCIATED(tb%quadbra)) do_quadrupole = .TRUE.
      reuse_native_input = .FALSE.
      IF (use_native_mixer) THEN
         IF (scf_env%iter_count == 1) THEN
            reuse_native_input = ANY(ABS(tb%wfn%qsh) > 1.0E-14_dp)
            IF (do_dipole) reuse_native_input = reuse_native_input .OR. &
                                                ANY(ABS(tb%wfn%dpat) > 1.0E-14_dp)
            IF (do_quadrupole) reuse_native_input = reuse_native_input .OR. &
                                                    ANY(ABS(tb%wfn%qpat) > 1.0E-14_dp)
         END IF
      END IF
      n_atom = SIZE(particle_set)
      nkind = SIZE(atomic_kind_set)
      nimg = dft_control%nimages
      CALL get_qs_kind_set(qs_kind_set, maxsgf=nsgf)
      ALLOCATE (aocg(nsgf, n_atom))
      ALLOCATE (aocg_spin(nsgf, n_atom, nspin))
      aocg = 0.0_dp
      aocg_spin = 0.0_dp
      IF (do_dipole) THEN
         ALLOCATE (ao_dip_spin(n_atom, dip_n, nspin))
         ao_dip_spin = 0.0_dp
      END IF
      IF (do_quadrupole) THEN
         ALLOCATE (ao_quad_spin(n_atom, quad_n, nspin))
         ao_quad_spin = 0.0_dp
      END IF
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
      ALLOCATE (ch_atom(n_atom, nspin))
      ALLOCATE (ch_orb_spin(max_orb, n_atom, nspin), ch_shell_spin(n_atom, max_shell, nspin))
      ch_atom = 0.0_dp
      ch_orb_spin = 0.0_dp
      ch_shell_spin = 0.0_dp
      IF (nimg > 1) THEN
         DO ispin = 1, nspin
            CALL tb_ao_charges_kp_spin(matrix_p, matrix_s, aocg_spin(:, :, ispin), ispin, para_env)
            IF (do_dipole) THEN
               DO im = 1, dip_n
                  CALL tb_contract_dens_kp_spin(matrix_p, tb%dipbra, tb%dipket, im, dip_n, &
                                                ao_dip_spin(:, im, ispin), ispin, para_env)
               END DO
            END IF
            IF (do_quadrupole) THEN
               DO im = 1, quad_n
                  CALL tb_contract_dens_kp_spin(matrix_p, tb%quadbra, tb%quadket, im, quad_n, &
                                                ao_quad_spin(:, im, ispin), ispin, para_env)
               END DO
            END IF
         END DO
      ELSE
         NULLIFY (p_matrix, s_matrix)
         p_matrix => matrix_p(:, 1)
         s_matrix => matrix_s(1, 1)%matrix
         DO ispin = 1, nspin
            CALL tb_ao_charges_matrix(matrix_p(ispin, 1)%matrix, s_matrix, aocg_spin(:, :, ispin), para_env)
            IF (do_dipole) THEN
               DO im = 1, dip_n
                  CALL tb_contract_dens_matrix(matrix_p(ispin, 1)%matrix, tb%dipbra(im)%matrix, &
                                               tb%dipket(im)%matrix, ao_dip_spin(:, im, ispin), para_env)
               END DO
            END IF
            IF (do_quadrupole) THEN
               DO im = 1, quad_n
                  CALL tb_contract_dens_matrix(matrix_p(ispin, 1)%matrix, tb%quadbra(im)%matrix, &
                                               tb%quadket(im)%matrix, ao_quad_spin(:, im, ispin), para_env)
               END DO
            END IF
         END DO
      END IF
      IF (nspin == 1) THEN
         aocg(:, :) = aocg_spin(:, :, 1)
      ELSE
         aocg(:, :) = aocg_spin(:, :, 1) + aocg_spin(:, :, 2)
         IF (do_dipole) THEN
            DO im = 1, dip_n
               DO iatom = 1, n_atom
                  pao = ao_dip_spin(iatom, im, 1)
                  ao_dip_spin(iatom, im, 1) = pao + ao_dip_spin(iatom, im, 2)
                  ao_dip_spin(iatom, im, 2) = pao - ao_dip_spin(iatom, im, 2)
               END DO
            END DO
         END IF
         IF (do_quadrupole) THEN
            DO im = 1, quad_n
               DO iatom = 1, n_atom
                  pao = ao_quad_spin(iatom, im, 1)
                  ao_quad_spin(iatom, im, 1) = pao + ao_quad_spin(iatom, im, 2)
                  ao_quad_spin(iatom, im, 2) = pao - ao_quad_spin(iatom, im, 2)
               END DO
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
               ! Preserve tblite's native shell order, including repeated angular momenta.
               ns = nao(is)
               norm = 2*lao(is) + 1
               ch_orb_spin(is, atom_a, 1) = aocg(is, atom_a) - tb%calc%h0%refocc(nao(is), ikind)/norm
               IF (nspin == 2) ch_orb_spin(is, atom_a, 2) = &
                  aocg_spin(is, atom_a, 1) - aocg_spin(is, atom_a, 2)
               DO ispin = 1, nspin
                  ch_shell_spin(atom_a, ns, ispin) = ch_orb_spin(is, atom_a, ispin) + &
                                                     ch_shell_spin(atom_a, ns, ispin)
               END DO
            END DO
            DO ispin = 1, nspin
               ch_atom(atom_a, ispin) = SUM(ch_orb_spin(:, atom_a, ispin))
            END DO
         END DO
      END DO
      native_seed_charge = -SUM(ch_atom(:, 1))
      seed_native_from_rho = SUM(ABS(aocg_spin)) > 1.0E-10_dp .AND. &
                             ABS(native_seed_charge - REAL(dft_control%charge, dp)) < 1.0E-5_dp
      DEALLOCATE (aocg, aocg_spin)

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
               WRITE (raw_state_unit, "(I0,1X,ES24.16)") tb%calc%bas%ish_at(iatom) + is, -ch_shell_spin(iatom, is, 1)
            END DO
         END DO
         IF (do_dipole) THEN
            WRITE (raw_state_unit, *) "dpat"
            DO iatom = 1, n_atom
               WRITE (raw_state_unit, "(I0,3(1X,ES24.16))") iatom, -ao_dip_spin(iatom, :, 1)
            END DO
         END IF
         IF (do_quadrupole) THEN
            WRITE (raw_state_unit, *) "qpat"
            DO iatom = 1, n_atom
               WRITE (raw_state_unit, "(I0,6(1X,ES24.16))") iatom, -ao_quad_spin(iatom, :, 1)
            END DO
         END IF
         CLOSE (raw_state_unit)
      END IF

      ! Population arrays use CP2K electron-excess signs. The spin axis already
      ! contains total charge and magnetization, as required by tblite.
      IF (use_native_mixer) THEN
         IF (.NOT. ALLOCATED(tb%mixer)) CPABORT("tblite mixer not initialized")
         advance_native_mixer = .FALSE.
         IF (use_rho .AND. (.NOT. calculate_forces)) THEN
            advance_native_mixer = scf_env%iter_count > 1
         END IF
         IF (advance_native_mixer) THEN
            ! Build the next Hamiltonian from the previous SCC residual.
            CALL tb%mixer%next(error)
            IF (ALLOCATED(error)) CPABORT("tblite native mixer failed: "//error%message)
            CALL tb%mixer%get(tb%wfn%qsh)
            CALL tb_sum_shell_charges(tb)
            IF (do_dipole) CALL tb%mixer%get(tb%wfn%dpat)
            IF (do_quadrupole) CALL tb%mixer%get(tb%wfn%qpat)
         ELSE IF (use_rho .AND. (.NOT. calculate_forces)) THEN
            IF (.NOT. reuse_native_input) THEN
               IF (seed_native_from_rho) THEN
                  ! SCF_GUESS/RESTART controls the initial SCC state.
                  CALL tb_set_population(tb, ch_shell_spin, ao_dip_spin, ao_quad_spin)
               ELSE
                  tb%wfn%qsh(:, :) = 0.0_dp
                  tb%wfn%qat(:, :) = 0.0_dp
                  IF (do_dipole) tb%wfn%dpat(:, :, :) = 0.0_dp
                  IF (do_quadrupole) tb%wfn%qpat(:, :, :) = 0.0_dp
               END IF
            END IF
         ELSE
            ! Record input/output pairs only after diagonalization, never for forces.
            IF (.NOT. calculate_forces) THEN
               CALL tb%mixer%set(tb%wfn%qsh)
               IF (do_dipole) CALL tb%mixer%set(tb%wfn%dpat)
               IF (do_quadrupole) CALL tb%mixer%set(tb%wfn%qpat)
            END IF
            CALL tb_set_population(tb, ch_shell_spin, ao_dip_spin, ao_quad_spin)
            IF (.NOT. calculate_forces) THEN
               CALL tb%mixer%diff(tb%wfn%qsh)
               IF (do_dipole) CALL tb%mixer%diff(tb%wfn%dpat)
               IF (do_quadrupole) CALL tb%mixer%diff(tb%wfn%qpat)
            END IF
         END IF
      ELSE
         IF (.NOT. (dft_control%qs_control%do_ls_scf .OR. scf_control%use_ot)) THEN
            IF (.NOT. ASSOCIATED(scf_env)) CPABORT("CP2K SCC mixer requires a QS SCF environment")
            IF (.NOT. (do_dipole .OR. do_quadrupole) .AND. &
                scf_env%mixing_method == modified_broyden_mixing_nr) THEN
               CPABORT("MODIFIED_BROYDEN_MIXING with SCC_MIXER CP2K requires GFN2")
            END IF
         END IF
         IF (.NOT. use_no_mixer) THEN
            CALL tb_mix_population(scf_env, para_env, ch_shell_spin, ao_dip_spin, ao_quad_spin)
         END IF
         CALL tb_set_population(tb, ch_shell_spin, ao_dip_spin, ao_quad_spin)
         IF (nspin == 1) THEN
            IF (.NOT. (do_dipole .OR. do_quadrupole) .OR. &
                dft_control%qs_control%do_ls_scf .OR. scf_control%use_ot) THEN
               ! Preserve the existing atom-charge convention for shell-only mixing
               ! and direct density optimization: qat comes from the unmixed density.
               tb%wfn%qat(:, 1) = -ch_atom(:, 1)
            END IF
         END IF
      END IF

      CALL tb_update_scf_potential(qs_env, tb, para_env, calculate_forces, skip_scf_dispersion)

      IF (ALLOCATED(ao_dip_spin)) DEALLOCATE (ao_dip_spin)
      IF (ALLOCATED(ao_quad_spin)) DEALLOCATE (ao_quad_spin)
      DEALLOCATE (ch_atom, ch_orb_spin, ch_shell_spin)

#else
      MARK_USED(qs_env)
      MARK_USED(tb)
      MARK_USED(dft_control)
      MARK_USED(calculate_forces)
      MARK_USED(use_rho)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_update_charges

#if defined(__TBLITE)
! **************************************************************************************************
!> \brief Convert CP2K shell/multipole populations to tblite charge/magnetization state.
!> \param tb ...
!> \param shell Electron excess per atom, shell and spin channel.
!> \param dipole Electron dipoles, allocated only when required.
!> \param quadrupole Electron quadrupoles, allocated only when required.
! **************************************************************************************************
   SUBROUTINE tb_set_population(tb, shell, dipole, quadrupole)
      TYPE(tblite_type), INTENT(INOUT)                   :: tb
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: shell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(IN)                                      :: dipole, quadrupole

      INTEGER                                            :: iatom, ii, is, ispin

      DO iatom = 1, tb%mol%nat
         ii = tb%calc%bas%ish_at(iatom)
         DO ispin = 1, tb%wfn%nspin
            DO is = 1, tb%calc%bas%nsh_at(iatom)
               tb%wfn%qsh(ii + is, ispin) = -shell(iatom, is, ispin)
            END DO
            IF (ALLOCATED(dipole)) tb%wfn%dpat(:, iatom, ispin) = -dipole(iatom, :, ispin)
            IF (ALLOCATED(quadrupole)) tb%wfn%qpat(:, iatom, ispin) = -quadrupole(iatom, :, ispin)
         END DO
      END DO
      CALL tb_sum_shell_charges(tb)

   END SUBROUTINE tb_set_population

! **************************************************************************************************
!> \brief Reconstruct tblite atomic charges from its current shell charges.
!> \param tb ...
! **************************************************************************************************
   SUBROUTINE tb_sum_shell_charges(tb)
      TYPE(tblite_type), INTENT(INOUT)                   :: tb

      INTEGER                                            :: iatom, ii, ispin

      DO iatom = 1, tb%mol%nat
         ii = tb%calc%bas%ish_at(iatom)
         DO ispin = 1, tb%wfn%nspin
            tb%wfn%qat(iatom, ispin) = SUM(tb%wfn%qsh(ii + 1:ii + tb%calc%bas%nsh_at(iatom), ispin))
         END DO
      END DO

   END SUBROUTINE tb_sum_shell_charges

! **************************************************************************************************
!> \brief Apply the CP2K mixer to shell and multipole populations for all spin channels.
!> \param scf_env ...
!> \param para_env ...
!> \param shell ...
!> \param dipole ...
!> \param quadrupole ...
! **************************************************************************************************
   SUBROUTINE tb_mix_population(scf_env, para_env, shell, dipole, quadrupole)
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: shell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: dipole, quadrupole

      INTEGER                                            :: ispin, ncols, nsh, nspin, offset
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mix_vars

      nspin = SIZE(shell, 3)
      IF (nspin == 1 .AND. .NOT. ALLOCATED(dipole) .AND. .NOT. ALLOCATED(quadrupole)) THEN
         ! Retain CP2K's electron-excess convention for single-spin shell-only mixing.
         CALL charge_mixing(scf_env%mixing_method, scf_env%mixing_store, &
                            shell(:, :, 1), para_env, scf_env%iter_count)
         RETURN
      END IF

      nsh = SIZE(shell, 2)
      ncols = nspin*nsh
      IF (ALLOCATED(dipole)) ncols = ncols + nspin*dip_n
      IF (ALLOCATED(quadrupole)) ncols = ncols + nspin*quad_n
      ALLOCATE (mix_vars(SIZE(shell, 1), ncols))

      ! Preserve the mixer history layout: all shell channels, then dipoles, then quadrupoles.
      offset = 0
      DO ispin = 1, nspin
         mix_vars(:, offset + 1:offset + nsh) = -shell(:, :, ispin)
         offset = offset + nsh
      END DO
      IF (ALLOCATED(dipole)) THEN
         DO ispin = 1, nspin
            mix_vars(:, offset + 1:offset + dip_n) = -dipole(:, :, ispin)
            offset = offset + dip_n
         END DO
      END IF
      IF (ALLOCATED(quadrupole)) THEN
         DO ispin = 1, nspin
            mix_vars(:, offset + 1:offset + quad_n) = -quadrupole(:, :, ispin)
            offset = offset + quad_n
         END DO
      END IF
      CALL charge_mixing(scf_env%mixing_method, scf_env%mixing_store, &
                         mix_vars, para_env, scf_env%iter_count)
      offset = 0
      DO ispin = 1, nspin
         shell(:, :, ispin) = -mix_vars(:, offset + 1:offset + nsh)
         offset = offset + nsh
      END DO
      IF (ALLOCATED(dipole)) THEN
         DO ispin = 1, nspin
            dipole(:, :, ispin) = -mix_vars(:, offset + 1:offset + dip_n)
            offset = offset + dip_n
         END DO
      END IF
      IF (ALLOCATED(quadrupole)) THEN
         DO ispin = 1, nspin
            quadrupole(:, :, ispin) = -mix_vars(:, offset + 1:offset + quad_n)
            offset = offset + quad_n
         END DO
      END IF

   END SUBROUTINE tb_mix_population

! **************************************************************************************************
!> \brief Evaluate SCC potentials, energies and optional gradients for the current tblite state.
!> \param qs_env ...
!> \param tb ...
!> \param para_env ...
!> \param calculate_forces ...
!> \param skip_scf_dispersion ...
! **************************************************************************************************
   SUBROUTINE tb_update_scf_potential(qs_env, tb, para_env, calculate_forces, skip_scf_dispersion)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tblite_type), POINTER                         :: tb
      TYPE(mp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: calculate_forces, skip_scf_dispersion

      ! These library results are replicated across MPI ranks. The force conversion
      ! retains the existing normalization before CP2K's global force reduction.
      CALL tb%pot%reset
      tb%e_es = 0.0_dp
      tb%e_scd = 0.0_dp
      tb%e_int = 0.0_dp
      IF (ALLOCATED(tb%calc%coulomb)) THEN
         CALL tb%calc%coulomb%get_potential(tb%mol, tb%cache, tb%wfn, tb%pot)
         CALL tb%calc%coulomb%get_energy(tb%mol, tb%cache, tb%wfn, tb%e_es)
      END IF
      IF (ALLOCATED(tb%calc%dispersion)) THEN
         IF (.NOT. skip_scf_dispersion) THEN
            CALL tb%calc%dispersion%get_potential(tb%mol, tb%dcache, tb%wfn, tb%pot)
            CALL tb%calc%dispersion%get_energy(tb%mol, tb%dcache, tb%wfn, tb%e_scd)
         END IF
      END IF
      IF (ALLOCATED(tb%calc%interactions)) THEN
         CALL tb%calc%interactions%get_potential(tb%mol, tb%icache, tb%wfn, tb%pot)
         CALL tb%calc%interactions%get_energy(tb%mol, tb%icache, tb%wfn, tb%e_int)
      END IF

      IF (calculate_forces) THEN
         IF (ALLOCATED(tb%calc%coulomb)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%coulomb%get_gradient(tb%mol, tb%cache, tb%wfn, tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_coulomb", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 3)
         END IF

         IF (ALLOCATED(tb%calc%dispersion) .AND. .NOT. skip_scf_dispersion) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%dispersion%get_gradient(tb%mol, tb%dcache, tb%wfn, tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_dispersion_scf", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 2)
         END IF

         IF (ALLOCATED(tb%calc%interactions)) THEN
            tb%grad = 0.0_dp
            CALL tb%calc%interactions%get_gradient(tb%mol, tb%icache, tb%wfn, tb%grad, tb%sigma)
            CALL tb_dump_sigma_component("after_interactions_scf", tb%sigma, para_env)
            CALL tb_grad2force(qs_env, tb, para_env, 3)
         END IF
      END IF

   END SUBROUTINE tb_update_scf_potential
#endif

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
                                                            is, nimg, ni, nj, i, j, nspin
      INTEGER                                            :: la, lb, za, zb
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
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ashift, bshift
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: ksblock, sblock
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dip_ket1, dip_ket2, dip_ket3, &
                                                            dip_bra1, dip_bra2, dip_bra3
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: quad_ket1, quad_ket2, quad_ket3, &
                                                            quad_ket4, quad_ket5, quad_ket6, &
                                                            quad_bra1, quad_bra2, quad_bra3, &
                                                            quad_bra4, quad_bra5, quad_bra6

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: kp_list
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(xtb_atom_type), POINTER                       :: xtb_atom_a, xtb_atom_b

      nimg = dft_control%nimages
      mpfac = -0.5_dp

      NULLIFY (matrix_s, ks_matrix, n_list, kp_list, qs_kind_set)
      CALL get_qs_env(qs_env=qs_env, sab_orb=n_list, sab_kp=kp_list, &
                      matrix_s_kp=matrix_s, matrix_ks_kp=ks_matrix, qs_kind_set=qs_kind_set)
      IF (nimg > 1) THEN
         IF (.NOT. ASSOCIATED(kp_list)) THEN
            CPABORT("Missing k-point neighbor list for tblite Hamiltonian")
         END IF
         n_list => kp_list
      END IF
      nspin = SIZE(ks_matrix, 1)

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

            nj = SIZE(sblock, 2)
            ALLOCATE (bshift(nj, nj))

            DO is = 1, nspin
               ashift = 0.0_dp
               DO i = 1, ni
                  la = naoa(i) + sum_shell(irow)
                  ashift(i, i) = tb_spin_project(tb%pot%vsh(la, :), is)
               END DO
               bshift = 0.0_dp
               DO j = 1, nj
                  lb = naob(j) + sum_shell(icol)
                  bshift(j, j) = tb_spin_project(tb%pot%vsh(lb, :), is)
               END DO
               NULLIFY (ksblock)
               CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                      row=irow, col=icol, block=ksblock, found=found)
               CPASSERT(found)
               ksblock = ksblock - 0.5_dp*(MATMUL(ashift, sblock) &
                                           + MATMUL(sblock, bshift))
               ksblock = ksblock - 0.5_dp*(tb_spin_project(tb%pot%vat(irow, :), is) &
                                           + tb_spin_project(tb%pot%vat(icol, :), is))*sblock
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

               DO is = 1, nspin
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock + mpfac*(dip_ket1*tb_spin_project(tb%pot%vdp(1, irow, :), is) &
                                             + dip_ket2*tb_spin_project(tb%pot%vdp(2, irow, :), is) &
                                             + dip_ket3*tb_spin_project(tb%pot%vdp(3, irow, :), is) &
                                             + dip_bra1*tb_spin_project(tb%pot%vdp(1, icol, :), is) &
                                             + dip_bra2*tb_spin_project(tb%pot%vdp(2, icol, :), is) &
                                             + dip_bra3*tb_spin_project(tb%pot%vdp(3, icol, :), is))
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

               DO is = 1, nspin
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, 1)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)

                  ksblock = ksblock + mpfac*(quad_ket1*tb_spin_project(tb%pot%vqp(1, irow, :), is) &
                                             + quad_ket2*tb_spin_project(tb%pot%vqp(2, irow, :), is) &
                                             + quad_ket3*tb_spin_project(tb%pot%vqp(3, irow, :), is) &
                                             + quad_ket4*tb_spin_project(tb%pot%vqp(4, irow, :), is) &
                                             + quad_ket5*tb_spin_project(tb%pot%vqp(5, irow, :), is) &
                                             + quad_ket6*tb_spin_project(tb%pot%vqp(6, irow, :), is) &
                                             + quad_bra1*tb_spin_project(tb%pot%vqp(1, icol, :), is) &
                                             + quad_bra2*tb_spin_project(tb%pot%vqp(2, icol, :), is) &
                                             + quad_bra3*tb_spin_project(tb%pot%vqp(3, icol, :), is) &
                                             + quad_bra4*tb_spin_project(tb%pot%vqp(4, icol, :), is) &
                                             + quad_bra5*tb_spin_project(tb%pot%vqp(5, icol, :), is) &
                                             + quad_bra6*tb_spin_project(tb%pot%vqp(6, icol, :), is))
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

            IF (iatom > jatom) THEN
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

            nj = SIZE(sblock, 2)
            ALLOCATE (bshift(nj, nj))

            DO is = 1, nspin
               ashift = 0.0_dp
               DO i = 1, ni
                  la = naoa(i) + sum_shell(irow)
                  ashift(i, i) = tb_spin_project(tb%pot%vsh(la, :), is)
               END DO
               bshift = 0.0_dp
               DO j = 1, nj
                  lb = naob(j) + sum_shell(icol)
                  bshift(j, j) = tb_spin_project(tb%pot%vsh(lb, :), is)
               END DO
               NULLIFY (ksblock)
               CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                      row=irow, col=icol, block=ksblock, found=found)
               CPASSERT(found)
               ksblock = ksblock - 0.5_dp*(MATMUL(ashift, sblock) &
                                           + MATMUL(sblock, bshift))
               ksblock = ksblock - 0.5_dp*(tb_spin_project(tb%pot%vat(irow, :), is) &
                                           + tb_spin_project(tb%pot%vat(icol, :), is))*sblock
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

               DO is = 1, nspin
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)
                  ksblock = ksblock + mpfac*(dip_ket1*tb_spin_project(tb%pot%vdp(1, irow, :), is) &
                                             + dip_ket2*tb_spin_project(tb%pot%vdp(2, irow, :), is) &
                                             + dip_ket3*tb_spin_project(tb%pot%vdp(3, irow, :), is) &
                                             + dip_bra1*tb_spin_project(tb%pot%vdp(1, icol, :), is) &
                                             + dip_bra2*tb_spin_project(tb%pot%vdp(2, icol, :), is) &
                                             + dip_bra3*tb_spin_project(tb%pot%vdp(3, icol, :), is))
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

               DO is = 1, nspin
                  NULLIFY (ksblock)
                  CALL dbcsr_get_block_p(matrix=ks_matrix(is, ic)%matrix, &
                                         row=irow, col=icol, block=ksblock, found=found)
                  CPASSERT(found)

                  ksblock = ksblock + mpfac*(quad_ket1*tb_spin_project(tb%pot%vqp(1, irow, :), is) &
                                             + quad_ket2*tb_spin_project(tb%pot%vqp(2, irow, :), is) &
                                             + quad_ket3*tb_spin_project(tb%pot%vqp(3, irow, :), is) &
                                             + quad_ket4*tb_spin_project(tb%pot%vqp(4, irow, :), is) &
                                             + quad_ket5*tb_spin_project(tb%pot%vqp(5, irow, :), is) &
                                             + quad_ket6*tb_spin_project(tb%pot%vqp(6, irow, :), is) &
                                             + quad_bra1*tb_spin_project(tb%pot%vqp(1, icol, :), is) &
                                             + quad_bra2*tb_spin_project(tb%pot%vqp(2, icol, :), is) &
                                             + quad_bra3*tb_spin_project(tb%pot%vqp(3, icol, :), is) &
                                             + quad_bra4*tb_spin_project(tb%pot%vqp(4, icol, :), is) &
                                             + quad_bra5*tb_spin_project(tb%pot%vqp(5, icol, :), is) &
                                             + quad_bra6*tb_spin_project(tb%pot%vqp(6, icol, :), is))
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
      INTEGER                                     :: nkind, natom, handle, nimg, i, inda, indb
      INTEGER                                     :: atom_a, atom_b, nseta, nsetb, ia, ib, ij
      LOGICAL                                     :: found
      REAL(KIND=dp)                               :: r2
      INTEGER, DIMENSION(3)                       :: cell
      INTEGER, DIMENSION(:, :, :), POINTER        :: cell_to_index
      REAL(KIND=dp), DIMENSION(3)                 :: rij
      INTEGER, DIMENSION(:), POINTER              :: la_max, lb_max
      INTEGER, DIMENSION(:), POINTER              :: nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER           :: first_sgfa, first_sgfb
      INTEGER, ALLOCATABLE                        :: atom_of_kind(:)
      REAL(KIND=dp), ALLOCATABLE                  :: stmp(:)
      REAL(KIND=dp), ALLOCATABLE                  :: dtmp(:, :), qtmp(:, :), dtmpj(:, :), qtmpj(:, :)
      REAL(KIND=dp), DIMENSION(:, :), POINTER     :: dip_ket1, dip_ket2, dip_ket3, &
                                                     dip_bra1, dip_bra2, dip_bra3
      REAL(KIND=dp), DIMENSION(:, :), POINTER     :: quad_ket1, quad_ket2, quad_ket3, &
                                                     quad_ket4, quad_ket5, quad_ket6, &
                                                     quad_bra1, quad_bra2, quad_bra3, &
                                                     quad_bra4, quad_bra5, quad_bra6

      TYPE(atomic_kind_type), DIMENSION(:), POINTER         :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                       :: dft_control
      TYPE(gto_basis_set_type), POINTER                     :: basis_set_a, basis_set_b
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER     :: basis_set_list
      TYPE(kpoint_type), POINTER                            :: kpoints
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_orb
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: sab_kp
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                              :: nl_iterator
      TYPE(particle_type), DIMENSION(:), POINTER            :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER             :: qs_kind_set

      CALL timeset(routineN, handle)

      !get info from environment vaiarable
      NULLIFY (atomic_kind_set, qs_kind_set, sab_orb, sab_kp, particle_set)
      NULLIFY (dft_control, matrix_s, kpoints, cell_to_index)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      sab_orb=sab_orb, &
                      sab_kp=sab_kp, &
                      particle_set=particle_set, &
                      dft_control=dft_control, &
                      kpoints=kpoints, &
                      matrix_s_kp=matrix_s)
      natom = SIZE(particle_set)
      nkind = SIZE(atomic_kind_set)
      nimg = dft_control%nimages
      IF (nimg > 1) THEN
         IF (.NOT. ASSOCIATED(sab_kp)) THEN
            CPABORT("Missing k-point neighbor list for tblite multipoles")
         END IF
         sab_orb => sab_kp
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
      END IF

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

      !loop over all atom pairs with a non-zero overlap (sab_orb)
      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)

         r2 = NORM2(rij(:))**2

         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)

         IF (iatom < jatom) THEN
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
                  CALL multipole_cgto(tb%calc%bas%cgto(jset, jtyp), tb%calc%bas%cgto(iset, ityp), &
                      & r2, -rij, tb%calc%bas%intcut, stmp, dtmp, qtmp)

                  DO inda = 1, nsgfa(iset)
                     ia = first_sgfa(1, iset) - first_sgfa(1, 1) + inda
                     DO indb = 1, nsgfb(jset)
                        ib = first_sgfb(1, jset) - first_sgfb(1, 1) + indb

                        ij = indb + nsgfb(jset)*(inda - 1)
                        CALL tb_shift_multipole(-rij, stmp(ij), dtmp(:, ij), qtmp(:, ij), &
                                                dtmpj(:, ij), qtmpj(:, ij))

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
!> \param p_mat ...
!> \param s_matrix ...
!> \param charges ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_ao_charges_matrix(p_mat, s_matrix, charges, para_env)
      TYPE(dbcsr_type), POINTER                          :: p_mat, s_matrix
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: charges
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: i, iblock_col, iblock_row, j
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: p_block, s_block
      TYPE(dbcsr_iterator_type)                          :: iter

      charges = 0.0_dp
      CALL dbcsr_iterator_start(iter, s_matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         NULLIFY (s_block, p_block)
         CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, s_block)
         CALL dbcsr_get_block_p(matrix=p_mat, row=iblock_row, col=iblock_col, BLOCK=p_block, found=found)
         IF (.NOT. found) CYCLE
         IF (.NOT. (ASSOCIATED(s_block) .AND. ASSOCIATED(p_block))) CYCLE

         DO j = 1, SIZE(p_block, 2)
            DO i = 1, SIZE(p_block, 1)
               charges(i, iblock_row) = charges(i, iblock_row) + p_block(i, j)*s_block(i, j)
            END DO
         END DO
         IF (iblock_col /= iblock_row) THEN
            DO j = 1, SIZE(p_block, 2)
               DO i = 1, SIZE(p_block, 1)
                  charges(j, iblock_col) = charges(j, iblock_col) + p_block(i, j)*s_block(i, j)
               END DO
            END DO
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)
      CALL para_env%sum(charges)

   END SUBROUTINE tb_ao_charges_matrix

! **************************************************************************************************
!> \brief compute the AO-resolved Mulliken charges for one k-point spin channel.
!> \param p_matrix_kp ...
!> \param s_matrix_kp ...
!> \param charges ...
!> \param ispin ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_ao_charges_kp_spin(p_matrix_kp, s_matrix_kp, charges, ispin, para_env)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: p_matrix_kp, s_matrix_kp
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: charges
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: ic
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: image_charges
      TYPE(dbcsr_type), POINTER                          :: p_mat, s_mat

      charges = 0.0_dp
      ALLOCATE (image_charges(SIZE(charges, 1), SIZE(charges, 2)))
      DO ic = 1, SIZE(s_matrix_kp, 2)
         NULLIFY (p_mat, s_mat)
         p_mat => p_matrix_kp(ispin, ic)%matrix
         s_mat => s_matrix_kp(1, ic)%matrix
         IF (ASSOCIATED(p_mat) .AND. ASSOCIATED(s_mat)) THEN
            image_charges = 0.0_dp
            CALL tb_ao_charges_matrix(p_mat, s_mat, image_charges, para_env)
            charges(:, :) = charges(:, :) + image_charges(:, :)
         END IF
      END DO
      DEALLOCATE (image_charges)

   END SUBROUTINE tb_ao_charges_kp_spin

! **************************************************************************************************
!> \brief compute the mulliken properties (AO resolved)
!> \param p_mat ...
!> \param bra_mat ...
!> \param ket_mat ...
!> \param output ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_contract_dens_matrix(p_mat, bra_mat, ket_mat, output, para_env)
      TYPE(dbcsr_type), POINTER                          :: p_mat, bra_mat, ket_mat
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: output
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: i, iblock_col, iblock_row, j
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: bra, ket, p_block
      TYPE(dbcsr_iterator_type)                          :: iter

      output = 0.0_dp
      CALL dbcsr_iterator_start(iter, bra_mat)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         NULLIFY (p_block, bra, ket)
         CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, bra)
         CALL dbcsr_get_block_p(matrix=p_mat, row=iblock_row, col=iblock_col, BLOCK=p_block, found=found)
         IF (.NOT. found) CYCLE
         CALL dbcsr_get_block_p(matrix=ket_mat, row=iblock_row, col=iblock_col, BLOCK=ket, found=found)
         IF (.NOT. found) CPABORT("missing block")

         IF (.NOT. (ASSOCIATED(bra) .AND. ASSOCIATED(p_block))) CYCLE
         IF (iblock_col == iblock_row) THEN
            DO j = 1, SIZE(p_block, 1)
               DO i = 1, SIZE(p_block, 2)
                  output(iblock_row) = output(iblock_row) + p_block(j, i)*bra(j, i)
               END DO
            END DO
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
      CALL para_env%sum(output)

   END SUBROUTINE tb_contract_dens_matrix

! **************************************************************************************************
!> \brief compute the AO-resolved density contraction for one k-point spin channel.
!> \param p_matrix ...
!> \param bra_mat ...
!> \param ket_mat ...
!> \param iop ...
!> \param nops ...
!> \param output ...
!> \param ispin ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE tb_contract_dens_kp_spin(p_matrix, bra_mat, ket_mat, iop, nops, output, ispin, para_env)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: p_matrix
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: bra_mat, ket_mat
      INTEGER, INTENT(IN)                                :: iop, nops
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: output
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: ic, idx, nimg
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: image_output
      TYPE(dbcsr_type), POINTER                          :: p_mat

      nimg = SIZE(p_matrix, 2)
      output = 0.0_dp
      ALLOCATE (image_output(SIZE(output)))
      DO ic = 1, nimg
         idx = iop + nops*(ic - 1)
         CPASSERT(idx <= SIZE(bra_mat))
         CPASSERT(idx <= SIZE(ket_mat))
         NULLIFY (p_mat)
         p_mat => p_matrix(ispin, ic)%matrix
         image_output = 0.0_dp
         CALL tb_contract_dens_matrix(p_mat, bra_mat(idx)%matrix, ket_mat(idx)%matrix, image_output, para_env)
         output = output + image_output
      END DO
      DEALLOCATE (image_output)

   END SUBROUTINE tb_contract_dens_kp_spin

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
!> \brief Add SCC-overlap and direct multipole Hamiltonian derivatives.
!> \param qs_env ...
!> \param use_rho ...
!> \param nimg ...
! **************************************************************************************************
   SUBROUTINE tb_derive_dH_off(qs_env, use_rho, nimg)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: use_rho
      INTEGER, INTENT(IN)                                :: nimg

#if defined(__TBLITE)
      INTEGER                                            :: i, idim, ij, iatom, ic, icol, ikind, img, ispin, &
                                                            jdim, ni, nj, nkind, nel, &
                                                            ityp, jatom, jkind, jrow, jtyp, iset, jset, nseti, nsetj, &
                                                            ia, ib, inda, indb, sampled_axes, sampled_even_axes, &
                                                            sampled_gamma_axes, nspin, ikp_axis
      INTEGER, DIMENSION(3)                              :: cellind, nkp_cellind, nkp_grid
      INTEGER, DIMENSION(:), POINTER                     :: nsgfa, nsgfb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: found, gamma_centered, gamma_sampled_image_pair, &
                                                            has_multipole_response, sampled_image_pair, &
                                                            use_matrix_scc_stress
      LOGICAL, DIMENSION(3)                              :: mesh_has_gamma
      REAL(KIND=dp)                                      :: r2, dr, i_a_shift, j_a_shift, i_a_shift_mag, j_a_shift_mag, &
                                                            ishift, jshift, ishift_mag, jshift_mag, pij_charge, &
                                                            pij_magnet, mp_pair_scale, kpoint_coordinate, native_dot_tmp
      REAL(KIND=dp), DIMENSION(3)                        :: kp_shift
      REAL(KIND=dp), DIMENSION(3)                        :: rij, dgrad, dhgrad_charge, dhgrad_magnet, &
                                                            mpgrad_charge, mpgrad_magnet
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hsigma
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: t_ov, idip, jdip, idip_mag, jdip_mag, &
                                                            iquad, jquad, iquad_mag, jquad_mag
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: t_dip, t_quad, t_d_ov
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: t_i_dip, t_i_quad, t_j_dip, t_j_quad
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :, :) :: scc_strain_hint
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: pblock, pblock_beta
      TYPE(block_p_type), DIMENSION(3, 3, 2)             :: scc_strain_blocks

      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_dh_scc, matrix_p
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list
      TYPE(gto_basis_set_type), POINTER            :: basis_set_a, basis_set_b
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_kp
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(tblite_type), POINTER                         :: tb
      TYPE(virial_type), POINTER                         :: virial

      ! compute mulliken charges required for charge update
      NULLIFY (scf_env, rho, tb, sab_orb, sab_kp, para_env, kpoints, matrix_dh_scc, virial)
      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      scf_env=scf_env, &
                      rho=rho, &
                      tb_tblite=tb, &
                      sab_orb=sab_orb, &
                      sab_kp=sab_kp, &
                      para_env=para_env, &
                      qs_kind_set=qs_kind_set, &
                      kpoints=kpoints)

      NULLIFY (cell_to_index)
      IF (nimg > 1) THEN
         IF (.NOT. ASSOCIATED(sab_kp)) CPABORT("Missing tblite k-point neighbor list")
         sab_orb => sab_kp
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index, nkp_grid=nkp_grid, &
                              kp_shift=kp_shift, gamma_centered=gamma_centered)
      ELSE
         nkp_grid = 1
         kp_shift = 0.0_dp
         gamma_centered = .TRUE.
      END IF
      mesh_has_gamma = .FALSE.
      DO i = 1, 3
         DO ikp_axis = 1, nkp_grid(i)
            IF (gamma_centered .AND. MODULO(nkp_grid(i), 2) == 0) THEN
               kpoint_coordinate = REAL(2*ikp_axis - nkp_grid(i), KIND=dp)/ &
                                   REAL(2*nkp_grid(i), KIND=dp) + kp_shift(i)
            ELSE
               kpoint_coordinate = REAL(2*ikp_axis - nkp_grid(i) - 1, KIND=dp)/ &
                                   REAL(2*nkp_grid(i), KIND=dp) + kp_shift(i)
            END IF
            IF (ABS(kpoint_coordinate - ANINT(kpoint_coordinate)) < 1.0e-12_dp) THEN
               mesh_has_gamma(i) = .TRUE.
            END IF
         END DO
      END DO
      has_multipole_response = ASSOCIATED(tb%dipbra) .OR. ASSOCIATED(tb%quadbra)

      NULLIFY (matrix_p)
      IF (use_rho) THEN
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      ELSE IF (ASSOCIATED(tb%rho_ao_kp_ref)) THEN
         matrix_p => tb%rho_ao_kp_ref
      ELSE
         matrix_p => scf_env%p_mix_new
      END IF
      nspin = SIZE(matrix_p, 1)
      use_matrix_scc_stress = tb%use_virial .AND. nimg > 1 .AND. ALL(tb%mol%periodic)

      IF (use_matrix_scc_stress) THEN
         CALL dbcsr_allocate_matrix_set(matrix_dh_scc, 9*nspin, nimg)
         DO img = 1, nimg
            DO ispin = 1, nspin
               DO jdim = 1, 3
                  DO idim = 1, 3
                     i = idim + 3*(jdim - 1) + 9*(ispin - 1)
                     ALLOCATE (matrix_dh_scc(i, img)%matrix)
                     CALL dbcsr_create(matrix_dh_scc(i, img)%matrix, template=matrix_p(ispin, img)%matrix, &
                                       name="TBLITE SCC STRAIN DERIVATIVE")
                     CALL cp_dbcsr_alloc_block_from_nbl(matrix_dh_scc(i, img)%matrix, sab_orb)
                  END DO
               END DO
            END DO
         END DO
      END IF

      ! set up basis set lists
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (basis_set_list(nkind))
      CALL basis_set_list_setup(basis_set_list, "ORB", qs_kind_set)

      nel = msao(tb%calc%bas%maxl)**2
      ALLOCATE (t_ov(nel))
      ALLOCATE (t_d_ov(3, nel))
      ALLOCATE (t_dip(dip_n, nel))
      ALLOCATE (t_i_dip(3, dip_n, nel), t_j_dip(3, dip_n, nel))
      ALLOCATE (t_quad(quad_n, nel))
      ALLOCATE (t_i_quad(3, quad_n, nel), t_j_quad(3, quad_n, nel))

      ALLOCATE (idip(dip_n), jdip(dip_n), idip_mag(dip_n), jdip_mag(dip_n))
      ALLOCATE (iquad(quad_n), jquad(quad_n), iquad_mag(quad_n), jquad_mag(quad_n))

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

         IF (iatom < jatom) THEN
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
         nkp_cellind = 0
         DO i = 1, 3
            IF (nkp_grid(i) > 1) THEN
               nkp_cellind(i) = MODULO(cellind(i), nkp_grid(i))
               IF (2*nkp_cellind(i) > nkp_grid(i)) nkp_cellind(i) = nkp_cellind(i) - nkp_grid(i)
            END IF
         END DO
         sampled_axes = 0
         IF (nkp_cellind(1) /= 0) sampled_axes = sampled_axes + 1
         IF (nkp_cellind(2) /= 0) sampled_axes = sampled_axes + 1
         IF (nkp_cellind(3) /= 0) sampled_axes = sampled_axes + 1
         sampled_even_axes = 0
         sampled_gamma_axes = 0
         DO i = 1, 3
            IF (nkp_grid(i) > 1 .AND. MODULO(nkp_grid(i), 2) == 0 .AND. &
                ABS(2*nkp_cellind(i)) == nkp_grid(i)) THEN
               sampled_even_axes = sampled_even_axes + 1
               IF (mesh_has_gamma(i)) sampled_gamma_axes = sampled_gamma_axes + 1
            END IF
         END DO
         sampled_image_pair = sampled_axes > 0 .AND. sampled_even_axes > 0
         gamma_sampled_image_pair = sampled_image_pair .AND. sampled_gamma_axes == sampled_even_axes
         mp_pair_scale = 1.0_dp
         IF (icol == jrow .AND. sampled_image_pair .AND. has_multipole_response .AND. &
             .NOT. gamma_sampled_image_pair) mp_pair_scale = -1.0_dp

         NULLIFY (pblock, pblock_beta)
         CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                row=jrow, col=icol, BLOCK=pblock, found=found)
         IF (.NOT. found) CPABORT("pblock not found")
         IF (nspin > 1) THEN
            CALL dbcsr_get_block_p(matrix=matrix_p(2, ic)%matrix, &
                                   row=jrow, col=icol, BLOCK=pblock_beta, found=found)
            IF (.NOT. found) CPABORT("pblock beta not found")
         END IF
         IF (use_matrix_scc_stress) THEN
            ALLOCATE (scc_strain_hint(SIZE(pblock, 2), SIZE(pblock, 1), 3, 3, nspin))
            scc_strain_hint = 0.0_dp
            DO ispin = 1, nspin
               DO jdim = 1, 3
                  DO idim = 1, 3
                     NULLIFY (scc_strain_blocks(idim, jdim, ispin)%block)
                     i = idim + 3*(jdim - 1) + 9*(ispin - 1)
                     CALL dbcsr_get_block_p(matrix=matrix_dh_scc(i, ic)%matrix, &
                                            row=jrow, col=icol, &
                                            BLOCK=scc_strain_blocks(idim, jdim, ispin)%block, found=found)
                     IF (.NOT. found) CPABORT("SCC strain derivative block not found")
                  END DO
               END DO
            END DO
         END IF
         i_a_shift = tb%pot%vat(icol, 1)
         j_a_shift = tb%pot%vat(jrow, 1)
         i_a_shift_mag = 0.0_dp
         j_a_shift_mag = 0.0_dp
         IF (SIZE(tb%pot%vat, 2) > 1) THEN
            i_a_shift_mag = tb%pot%vat(icol, 2)
            j_a_shift_mag = tb%pot%vat(jrow, 2)
         END IF
         idip(:) = tb%pot%vdp(:, icol, 1)
         jdip(:) = tb%pot%vdp(:, jrow, 1)
         idip_mag(:) = 0.0_dp
         jdip_mag(:) = 0.0_dp
         IF (SIZE(tb%pot%vdp, 3) > 1) THEN
            idip_mag(:) = tb%pot%vdp(:, icol, 2)
            jdip_mag(:) = tb%pot%vdp(:, jrow, 2)
         END IF
         iquad(:) = tb%pot%vqp(:, icol, 1)
         jquad(:) = tb%pot%vqp(:, jrow, 1)
         iquad_mag(:) = 0.0_dp
         jquad_mag(:) = 0.0_dp
         IF (SIZE(tb%pot%vqp, 3) > 1) THEN
            iquad_mag(:) = tb%pot%vqp(:, icol, 2)
            jquad_mag(:) = tb%pot%vqp(:, jrow, 2)
         END IF
         ni = tb%calc%bas%ish_at(icol)
         DO iset = 1, nseti
            ishift = i_a_shift + tb%pot%vsh(ni + iset, 1)
            ishift_mag = 0.0_dp
            IF (SIZE(tb%pot%vsh, 2) > 1) ishift_mag = i_a_shift_mag + tb%pot%vsh(ni + iset, 2)
            nj = tb%calc%bas%ish_at(jrow)
            DO jset = 1, nsetj
               jshift = j_a_shift + tb%pot%vsh(nj + jset, 1)
               jshift_mag = 0.0_dp
               IF (SIZE(tb%pot%vsh, 2) > 1) jshift_mag = j_a_shift_mag + tb%pot%vsh(nj + jset, 2)

               !get integrals and derivatives
               CALL multipole_grad_cgto(tb%calc%bas%cgto(iset, ityp), tb%calc%bas%cgto(jset, jtyp), &
                  & r2, rij, tb%calc%bas%intcut, t_ov, t_dip, t_quad, t_d_ov, t_i_dip, t_i_quad, &
                  & t_j_dip, t_j_quad)

               dgrad = 0.0_dp
               DO inda = 1, nsgfa(iset)
                  ia = first_sgfa(1, iset) - first_sgfa(1, 1) + inda
                  DO indb = 1, nsgfb(jset)
                     ib = first_sgfb(1, jset) - first_sgfb(1, 1) + indb

                     ij = inda + nsgfa(iset)*(indb - 1)

                     pij_charge = pblock(ib, ia)
                     pij_magnet = 0.0_dp
                     IF (nspin > 1) THEN
                        pij_charge = pij_charge + pblock_beta(ib, ia)
                        pij_magnet = pblock(ib, ia) - pblock_beta(ib, ia)
                     END IF
                     mpgrad_charge = MATMUL(t_i_dip(:, :, ij), idip) &
                                     + MATMUL(t_j_dip(:, :, ij), jdip) &
                                     + MATMUL(t_i_quad(:, :, ij), iquad) &
                                     + MATMUL(t_j_quad(:, :, ij), jquad)
                     mpgrad_magnet = MATMUL(t_i_dip(:, :, ij), idip_mag) &
                                     + MATMUL(t_j_dip(:, :, ij), jdip_mag) &
                                     + MATMUL(t_i_quad(:, :, ij), iquad_mag) &
                                     + MATMUL(t_j_quad(:, :, ij), jquad_mag)
                     dhgrad_charge = -(ishift + jshift)*t_d_ov(:, ij) - mp_pair_scale*mpgrad_charge
                     dhgrad_magnet = -(ishift_mag + jshift_mag)*t_d_ov(:, ij) - mp_pair_scale*mpgrad_magnet
                     dgrad(:) = dgrad(:) - &
                                ((ishift + jshift)*pij_charge + &
                                 (ishift_mag + jshift_mag)*pij_magnet)*t_d_ov(:, ij) - &
                                mp_pair_scale*(pij_charge*mpgrad_charge + pij_magnet*mpgrad_magnet)

                     IF (ALLOCATED(scc_strain_hint)) THEN
                        DO jdim = 1, 3
                           DO idim = 1, 3
                              scc_strain_hint(ia, ib, idim, jdim, 1) = &
                                 scc_strain_hint(ia, ib, idim, jdim, 1) &
                                 + (dhgrad_charge(idim) + MERGE(dhgrad_magnet(idim), 0.0_dp, nspin > 1))*rij(jdim)
                              IF (nspin > 1) THEN
                                 scc_strain_hint(ia, ib, idim, jdim, 2) = &
                                    scc_strain_hint(ia, ib, idim, jdim, 2) &
                                    + (dhgrad_charge(idim) - dhgrad_magnet(idim))*rij(jdim)
                              END IF
                           END DO
                        END DO
                     END IF

                  END DO
               END DO
               tb%grad(:, icol) = tb%grad(:, icol) - dgrad
               tb%grad(:, jrow) = tb%grad(:, jrow) + dgrad
               IF (tb%use_virial .AND. .NOT. use_matrix_scc_stress) THEN
                  IF (icol == jrow) THEN
                     DO ia = 1, 3
                        DO ib = 1, 3
                           IF (sampled_image_pair .AND. .NOT. gamma_sampled_image_pair) THEN
                              hsigma(ia, ib) = hsigma(ia, ib) - 0.25_dp* &
                                               (rij(ia)*dgrad(ib) + rij(ib)*dgrad(ia))
                           ELSE
                              hsigma(ia, ib) = hsigma(ia, ib) + 0.25_dp* &
                                               (rij(ia)*dgrad(ib) + rij(ib)*dgrad(ia))
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
         IF (ALLOCATED(scc_strain_hint)) THEN
            DO ispin = 1, nspin
               DO jdim = 1, 3
                  DO idim = 1, 3
                     IF (icol <= jrow) THEN
                        scc_strain_blocks(idim, jdim, ispin)%block(:, :) = &
                           scc_strain_blocks(idim, jdim, ispin)%block(:, :) &
                           + scc_strain_hint(:, :, idim, jdim, ispin)
                     ELSE
                        scc_strain_blocks(idim, jdim, ispin)%block(:, :) = &
                           scc_strain_blocks(idim, jdim, ispin)%block(:, :) &
                           + TRANSPOSE(scc_strain_hint(:, :, idim, jdim, ispin))
                     END IF
                  END DO
               END DO
            END DO
            DEALLOCATE (scc_strain_hint)
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      IF (use_matrix_scc_stress) THEN
         hsigma = 0.0_dp
         DO img = 1, nimg
            DO ispin = 1, nspin
               DO jdim = 1, 3
                  DO idim = 1, 3
                     i = idim + 3*(jdim - 1) + 9*(ispin - 1)
                     CALL dbcsr_finalize(matrix_dh_scc(i, img)%matrix)
                     CALL dbcsr_dot(matrix_dh_scc(i, img)%matrix, matrix_p(ispin, img)%matrix, native_dot_tmp)
                     ! Symmetric DBCSR storage represents both orientations of each assembled pair.
                     hsigma(idim, jdim) = hsigma(idim, jdim) + 0.5_dp*native_dot_tmp
                  END DO
               END DO
            END DO
         END DO
         CALL dbcsr_deallocate_matrix_set(matrix_dh_scc)
      ELSE
         CALL para_env%sum(hsigma)
      END IF
      CALL para_env%sum(tb%grad)
      CALL tb_grad2force(qs_env, tb, para_env, 4)

      IF (.NOT. use_matrix_scc_stress) tb%sigma = tb%sigma + hsigma

      DEALLOCATE (basis_set_list)
      DEALLOCATE (t_ov, t_d_ov)
      DEALLOCATE (t_dip, t_i_dip, t_j_dip)
      DEALLOCATE (t_quad, t_i_quad, t_j_quad)
      DEALLOCATE (idip, jdip, idip_mag, jdip_mag, iquad, jquad, iquad_mag, jquad_mag)

      IF (tb%use_virial) THEN
         CALL tb_add_stress(qs_env, tb, para_env)
         IF (use_matrix_scc_stress) THEN
            CALL get_qs_env(qs_env=qs_env, virial=virial)
            virial%pv_virial = virial%pv_virial - hsigma/para_env%num_pe
         END IF
      END IF

#else
      MARK_USED(qs_env)
      MARK_USED(use_rho)
      MARK_USED(nimg)
      CPABORT("Built without TBLITE")
#endif

   END SUBROUTINE tb_derive_dH_off

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

END MODULE tblite_interface
