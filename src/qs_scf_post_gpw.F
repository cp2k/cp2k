!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Does all kind of post scf calculations for GPW/GAPW
!> \par History
!>      Started as a copy from the relevant part of qs_scf
!>      Start to adapt for k-points [07.2015, JGH]
!> \author Joost VandeVondele (10.2003)
! **************************************************************************************************
MODULE qs_scf_post_gpw
   USE admm_types,                      ONLY: admm_type
   USE admm_utils,                      ONLY: admm_correct_for_eigenvalues,&
                                              admm_uncorrect_for_eigenvalues
   USE ai_onecenter,                    ONLY: sg_overlap
   USE atom_kind_orbitals,              ONLY: calculate_atomic_density
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_deallocate_matrix_set
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
   USE cp_ddapc_util,                   ONLY: get_ddapc
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_init_random,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE dct,                             ONLY: pw_shrink
   USE et_coupling_types,               ONLY: set_et_coupling_type
   USE hirshfeld_methods,               ONLY: comp_hirshfeld_charges,&
                                              comp_hirshfeld_i_charges,&
                                              create_shape_function,&
                                              write_hirshfeld_charges
   USE hirshfeld_types,                 ONLY: create_hirshfeld_type,&
                                              hirshfeld_type,&
                                              release_hirshfeld_type,&
                                              set_hirshfeld_info
   USE input_constants,                 ONLY: do_loc_both,&
                                              do_loc_homo,&
                                              do_loc_lumo,&
                                              ot_precond_full_all,&
                                              radius_covalent,&
                                              radius_user,&
                                              ref_charge_atomic,&
                                              ref_charge_mulliken
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_ivals,&
                                              section_get_lval,&
                                              section_get_rval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE lapack,                          ONLY: lapack_sgesv
   USE mao_wfn_analysis,                ONLY: mao_analysis
   USE mathconstants,                   ONLY: pi
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_sum
   USE minbas_wfn_analysis,             ONLY: minbas_analysis
   USE molden_utils,                    ONLY: write_mos_molden
   USE molecule_types,                  ONLY: molecule_type
   USE mulliken,                        ONLY: mulliken_charges
   USE orbital_pointers,                ONLY: indso
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE population_analyses,             ONLY: lowdin_population_analysis,&
                                              mulliken_population_analysis
   USE preconditioner_types,            ONLY: preconditioner_type
   USE ps_implicit_types,               ONLY: MIXED_BC,&
                                              MIXED_PERIODIC_BC,&
                                              NEUMANN_BC,&
                                              PERIODIC_BC
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grids,                        ONLY: get_pw_grid_info
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_derive,&
                                              pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_implicit,&
                                              pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_type
   USE qs_active_space_methods,         ONLY: active_space_main
   USE qs_collocate_density,            ONLY: calculate_wavefunction
   USE qs_commutators,                  ONLY: build_com_hr_matrix
   USE qs_core_energies,                ONLY: calculate_ptrace
   USE qs_dos,                          ONLY: calculate_dos,&
                                              calculate_dos_kp
   USE qs_electric_field_gradient,      ONLY: qs_efg_calc
   USE qs_elf_methods,                  ONLY: qs_elf_calc
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_energy_window,                ONLY: energy_windows
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_epr_hyp,                      ONLY: qs_epr_hyp_calc
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace,&
                                              qs_ks_update_qs_env
   USE qs_ks_types,                     ONLY: qs_ks_did_change
   USE qs_loc_dipole,                   ONLY: loc_dipole
   USE qs_loc_states,                   ONLY: get_localization_info
   USE qs_loc_types,                    ONLY: qs_loc_env_create,&
                                              qs_loc_env_destroy,&
                                              qs_loc_env_new_type
   USE qs_loc_utils,                    ONLY: loc_write_restart,&
                                              qs_loc_control_init,&
                                              qs_loc_env_init,&
                                              qs_loc_init,&
                                              retain_history
   USE qs_local_properties,             ONLY: qs_local_energy,&
                                              qs_local_stress
   USE qs_mo_io,                        ONLY: write_dm_binary_restart,&
                                              write_mo_set
   USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues,&
                                              make_mo_eig
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_moments,                      ONLY: qs_moment_berry_phase,&
                                              qs_moment_locop
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              get_neighbor_list_set_p,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_ot_eigensolver,               ONLY: ot_eigensolver
   USE qs_pdos,                         ONLY: calculate_projected_dos
   USE qs_resp,                         ONLY: resp_fit
   USE qs_rho0_types,                   ONLY: get_rho0_mpole,&
                                              mpole_rho_atom,&
                                              rho0_mpole_type
   USE qs_rho_atom_types,               ONLY: rho_atom_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_csr_write,                ONLY: write_ks_matrix_csr,&
                                              write_s_matrix_csr
   USE qs_scf_types,                    ONLY: ot_method_nr,&
                                              qs_scf_env_type
   USE qs_scf_wfn_mix,                  ONLY: wfn_mix
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
   USE qs_wannier90,                    ONLY: wannier90_interface
   USE s_square_methods,                ONLY: compute_s_square
   USE scf_control_types,               ONLY: scf_control_type
   USE stm_images,                      ONLY: th_stm_image
   USE transport,                       ONLY: qs_scf_post_transport
   USE virial_types,                    ONLY: virial_type
   USE xray_diffraction,                ONLY: calculate_rhotot_elec_gspace,&
                                              xray_diffraction_spectrum
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Global parameters
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_post_gpw'
   PUBLIC :: scf_post_calculation_gpw, &
             qs_scf_post_moments, &
             write_available_results, &
             write_mo_free_results
   PUBLIC :: make_lumo_gpw

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief collects possible post - scf calculations and prints info / computes properties.
!> \param qs_env the qs_env in which the qs_env lives
!> \param wf_type ...
!> \par History
!>      02.2003 created [fawzi]
!>      10.2004 moved here from qs_scf [Joost VandeVondele]
!>              started splitting out different subroutines
!>      10.2015 added header for wave-function correlated methods [Vladimir Rybkin]
!> \author fawzi
!> \note
!>      this function changes mo_eigenvectors and mo_eigenvalues, depending on the print keys.
!>      In particular, MO_CUBES causes the MOs to be rotated to make them eigenstates of the KS
!>      matrix, and mo_eigenvalues is updated accordingly. This can, for unconverged wavefunctions,
!>      change afterwards slightly the forces (hence small numerical differences between MD
!>      with and without the debug print level). Ideally this should not happen...
! **************************************************************************************************
   SUBROUTINE scf_post_calculation_gpw(qs_env, wf_type)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      CHARACTER(6), OPTIONAL                             :: wf_type

      CHARACTER(len=*), PARAMETER :: routineN = 'scf_post_calculation_gpw', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, homo, ispin, min_lumos, n_rep, &
                                                            nhomo, nlumo, nlumo_stm, nlumo_tddft, &
                                                            nlumos, nmo, output_unit, unit_nr
      INTEGER, DIMENSION(:, :, :), POINTER               :: marked_states
      LOGICAL :: check_write, compute_lumos, do_homo, do_kpoints, do_mo_cubes, do_stm, &
         do_wannier_cubes, has_homo, has_lumo, loc_explicit, loc_print_explicit, my_localized_wfn, &
         p_loc, p_loc_homo, p_loc_lumo
      REAL(dp)                                           :: e_kin
      REAL(KIND=dp)                                      :: gap, homo_lumo(2, 2)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: occupied_evals, unoccupied_evals, &
                                                            unoccupied_evals_stm
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: homo_localized, lumo_localized, lumo_ptr, &
         mo_loc_history, occupied_orbs, unoccupied_orbs, unoccupied_orbs_stm
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_rmpv, matrix_s, mo_derivs
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: kinetic_m, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(molecule_type), POINTER                       :: molecule_set(:)
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: wf_g, wf_r
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env_homo, qs_loc_env_lumo
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input, loc_print_section, &
                                                            localize_section, print_key, &
                                                            stm_section

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      ! Print out the type of wavefunction to distinguish between SCF and post-SCF
      IF (PRESENT(wf_type)) THEN
         IF (output_unit > 0) THEN
            WRITE (UNIT=output_unit, FMT='(/,(T1,A))') REPEAT("-", 40)
            WRITE (UNIT=output_unit, FMT='(/,(T3,A,T19,A,T25,A))') "Properties from ", wf_type, " density"
            WRITE (UNIT=output_unit, FMT='(/,(T1,A))') REPEAT("-", 40)
         ENDIF
      ENDIF

      ! Writes the data that is already available in qs_env
      CALL get_qs_env(qs_env, scf_env=scf_env)
      CALL write_available_results(qs_env, scf_env)

      my_localized_wfn = .FALSE.
      NULLIFY (admm_env, dft_control, pw_env, auxbas_pw_pool, pw_pools, mos, rho, &
               mo_coeff, ks_rmpv, matrix_s, qs_loc_env_homo, qs_loc_env_lumo, scf_control, &
               unoccupied_orbs, unoccupied_orbs_stm, mo_eigenvalues, unoccupied_evals, &
               unoccupied_evals_stm, molecule_set, mo_derivs, &
               subsys, particles, input, print_key, kinetic_m, marked_states)
      NULLIFY (homo_localized, lumo_localized, lumo_ptr, rho_ao)

      has_homo = .FALSE.
      has_lumo = .FALSE.
      p_loc = .FALSE.
      p_loc_homo = .FALSE.
      p_loc_lumo = .FALSE.

      CPASSERT(ASSOCIATED(scf_env))
      CPASSERT(scf_env%ref_count > 0)
      CPASSERT(ASSOCIATED(qs_env))
      ! Here we start with data that needs a postprocessing...
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      molecule_set=molecule_set, &
                      scf_control=scf_control, &
                      do_kpoints=do_kpoints, &
                      input=input, &
                      subsys=subsys, &
                      rho=rho, &
                      pw_env=pw_env, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set)
      CALL qs_subsys_get(subsys, particles=particles)

      CALL qs_rho_get(rho, rho_ao_kp=rho_ao)

      ! In MP2 case update the Hartree potential
      IF (ASSOCIATED(qs_env%mp2_env)) CALL update_hartree_with_mp2(rho, qs_env)

      !    **** the kinetic energy
      IF (cp_print_key_should_output(logger%iter_info, input, &
                                     "DFT%PRINT%KINETIC_ENERGY") /= 0) THEN
         CALL get_qs_env(qs_env, kinetic_kp=kinetic_m)
         CPASSERT(ASSOCIATED(kinetic_m))
         CPASSERT(ASSOCIATED(kinetic_m(1, 1)%matrix))
         CALL calculate_ptrace(kinetic_m, rho_ao, e_kin, dft_control%nspins)
         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%KINETIC_ENERGY", &
                                        extension=".Log")
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(T3,A,T55,F25.14)') "Electronic kinetic energy:", e_kin
         ENDIF
         CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                           "DFT%PRINT%KINETIC_ENERGY")
      END IF

      ! Atomic Charges that require further computation
      CALL qs_scf_post_charges(input, logger, qs_env)

      ! Moments of charge distribution
      CALL qs_scf_post_moments(input, logger, qs_env, output_unit)

      ! Determine if we need to computer properties using the localized centers
      dft_section => section_vals_get_subs_vals(input, "DFT")
      localize_section => section_vals_get_subs_vals(dft_section, "LOCALIZE")
      loc_print_section => section_vals_get_subs_vals(localize_section, "PRINT")
      CALL section_vals_get(localize_section, explicit=loc_explicit)
      CALL section_vals_get(loc_print_section, explicit=loc_print_explicit)

      ! Print_keys controlled by localization
      IF (loc_print_explicit) THEN
         print_key => section_vals_get_subs_vals(loc_print_section, "MOLECULAR_DIPOLES")
         p_loc = BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)
         print_key => section_vals_get_subs_vals(loc_print_section, "TOTAL_DIPOLE")
         p_loc = p_loc .OR. BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)
         print_key => section_vals_get_subs_vals(loc_print_section, "WANNIER_CENTERS")
         p_loc = p_loc .OR. BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)
         print_key => section_vals_get_subs_vals(loc_print_section, "WANNIER_SPREADS")
         p_loc = p_loc .OR. BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)
         print_key => section_vals_get_subs_vals(loc_print_section, "WANNIER_CUBES")
         p_loc = p_loc .OR. BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)
         print_key => section_vals_get_subs_vals(loc_print_section, "MOLECULAR_STATES")
         p_loc = p_loc .OR. BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)
         print_key => section_vals_get_subs_vals(loc_print_section, "MOLECULAR_MOMENTS")
         p_loc = p_loc .OR. BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)
      ELSE
         p_loc = .FALSE.
      END IF
      IF (loc_explicit) THEN
         p_loc_homo = (section_get_ival(localize_section, "STATES") == do_loc_homo .OR. &
                       section_get_ival(localize_section, "STATES") == do_loc_both) .AND. p_loc
         p_loc_lumo = (section_get_ival(localize_section, "STATES") == do_loc_lumo .OR. &
                       section_get_ival(localize_section, "STATES") == do_loc_both) .AND. p_loc
         CALL section_vals_val_get(localize_section, "LIST_UNOCCUPIED", n_rep_val=n_rep)
      ELSE
         p_loc_homo = .FALSE.
         p_loc_lumo = .FALSE.
         n_rep = 0
      END IF

      IF (n_rep == 0 .AND. p_loc_lumo) THEN
         CALL cp_abort(__LOCATION__, "No LIST_UNOCCUPIED was specified, "// &
                       "therefore localization of unoccupied states will be skipped!")
         p_loc_lumo = .FALSE.
      END IF
      print_key => section_vals_get_subs_vals(loc_print_section, "WANNIER_STATES")
      p_loc = p_loc .OR. BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)

      ! Control for STM
      stm_section => section_vals_get_subs_vals(input, "DFT%PRINT%STM")
      CALL section_vals_get(stm_section, explicit=do_stm)
      nlumo_stm = 0
      IF (do_stm) nlumo_stm = section_get_ival(stm_section, "NLUMO")

      ! check for CUBES (MOs and WANNIERS)
      do_mo_cubes = BTEST(cp_print_key_should_output(logger%iter_info, dft_section, "PRINT%MO_CUBES") &
                          , cp_p_file)
      IF (loc_print_explicit) THEN
         do_wannier_cubes = BTEST(cp_print_key_should_output(logger%iter_info, loc_print_section, &
                                                             "WANNIER_CUBES"), cp_p_file)
      ELSE
         do_wannier_cubes = .FALSE.
      END IF
      nlumo = section_get_ival(dft_section, "PRINT%MO_CUBES%NLUMO")
      nhomo = section_get_ival(dft_section, "PRINT%MO_CUBES%NHOMO")
      nlumo_tddft = 0
      IF (dft_control%do_tddfpt_calculation) THEN
         nlumo_tddft = section_get_ival(dft_section, "TDDFPT%NLUMO")
      END IF

      ! Setup the grids needed to compute a wavefunction given a vector..
      IF (((do_mo_cubes .OR. do_wannier_cubes) .AND. (nlumo /= 0 .OR. nhomo /= 0)) .OR. p_loc) THEN
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                         pw_pools=pw_pools)
         CALL pw_pool_create_pw(auxbas_pw_pool, wf_r%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, wf_g%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
      END IF

      ! Makes the MOs eigenstates, computes eigenvalues, write cubes
      IF (do_kpoints) THEN
         IF (do_mo_cubes) THEN
            CPWARN("Print MO Cubes not implemented for k-point calculations!!")
         END IF
      ELSE
         CALL get_qs_env(qs_env, &
                         mos=mos, &
                         matrix_ks=ks_rmpv)
         IF ((do_mo_cubes .AND. nhomo /= 0) .OR. do_stm .OR. dft_control%do_tddfpt_calculation) THEN
            IF (dft_control%restricted) THEN
               IF (output_unit > 0) WRITE (output_unit, *) &
                  " Unclear how we define MOs in the restricted case ... skipping"
            ELSE
               CALL get_qs_env(qs_env, mo_derivs=mo_derivs)
               IF (dft_control%do_admm) THEN
                  CALL get_qs_env(qs_env, admm_env=admm_env)
                  CALL make_mo_eig(mos, dft_control%nspins, ks_rmpv, scf_control, mo_derivs, admm_env=admm_env)
               ELSE
                  CALL make_mo_eig(mos, dft_control%nspins, ks_rmpv, scf_control, mo_derivs)
               END IF
            END IF
            DO ispin = 1, dft_control%nspins
               CALL get_mo_set(mo_set=mos(ispin)%mo_set, eigenvalues=mo_eigenvalues, homo=homo)
               homo_lumo(ispin, 1) = mo_eigenvalues(homo)
            END DO
            has_homo = .TRUE.
         END IF
         IF (do_mo_cubes .AND. nhomo /= 0) THEN
            DO ispin = 1, dft_control%nspins
               ! Prints the cube files of OCCUPIED ORBITALS
               CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, &
                               eigenvalues=mo_eigenvalues, homo=homo, nmo=nmo)
               CALL qs_scf_post_occ_cubes(input, dft_section, dft_control, logger, qs_env, &
                                          mo_coeff, wf_g, wf_r, particles, homo, ispin)
            END DO
         ENDIF
      ENDIF

      ! Initialize the localization environment, needed e.g. for wannier functions and molecular states
      ! Gets localization info for the occupied orbs
      !  - Possibly gets wannier functions
      !  - Possibly gets molecular states
      IF (p_loc_homo) THEN
         IF (do_kpoints) THEN
            CPWARN("Localization not implemented for k-point calculations!!")
         ELSEIF (dft_control%restricted) THEN
            IF (output_unit > 0) WRITE (output_unit, *) &
               " Unclear how we define MOs / localization in the restricted case ... skipping"
         ELSE
            ALLOCATE (occupied_orbs(dft_control%nspins))
            ALLOCATE (occupied_evals(dft_control%nspins))
            ALLOCATE (homo_localized(dft_control%nspins))
            DO ispin = 1, dft_control%nspins
               CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, &
                               eigenvalues=mo_eigenvalues)
               occupied_orbs(ispin)%matrix => mo_coeff
               occupied_evals(ispin)%array => mo_eigenvalues
               CALL cp_fm_create(homo_localized(ispin)%matrix, occupied_orbs(ispin)%matrix%matrix_struct)
               CALL cp_fm_to_fm(occupied_orbs(ispin)%matrix, homo_localized(ispin)%matrix)
            END DO

            CALL get_qs_env(qs_env, mo_loc_history=mo_loc_history)
            do_homo = .TRUE.

            CALL qs_loc_env_create(qs_loc_env_homo)
            CALL qs_loc_control_init(qs_loc_env_homo, localize_section, do_homo=do_homo)
            CALL qs_loc_init(qs_env, qs_loc_env_homo, localize_section, homo_localized, do_homo, &
                             do_mo_cubes, mo_loc_history=mo_loc_history)
            CALL get_localization_info(qs_env, qs_loc_env_homo, localize_section, homo_localized, &
                                       wf_r, wf_g, particles, occupied_orbs, occupied_evals, marked_states)

            !retain the homo_localized for future use
            IF (qs_loc_env_homo%localized_wfn_control%use_history) THEN
               CALL retain_history(mo_loc_history, homo_localized)
               CALL set_qs_env(qs_env, mo_loc_history=mo_loc_history)
            ENDIF

            !write restart for localization of occupied orbitals
            CALL loc_write_restart(qs_loc_env_homo, loc_print_section, mos, &
                                   homo_localized, do_homo)
            CALL cp_fm_vect_dealloc(homo_localized)
            DEALLOCATE (occupied_orbs)
            DEALLOCATE (occupied_evals)
            ! Print Total Dipole if the localization has been performed
            IF (qs_loc_env_homo%do_localize) THEN
               CALL loc_dipole(input, dft_control, qs_loc_env_homo, logger, qs_env)
            END IF
         END IF
      ENDIF

      ! Gets the lumos, and eigenvalues for the lumos, and localize them if requested
      IF (do_kpoints) THEN
         IF (do_mo_cubes .OR. p_loc_lumo) THEN
            ! nothing at the moment, not implemented
            CPWARN("Localization and MO related output not implemented for k-point calculations!")
         END IF
      ELSE
         IF (nlumo .GT. -1) THEN
            nlumo = MAX(nlumo, nlumo_tddft)
         END IF
         compute_lumos = (do_mo_cubes .OR. dft_control%do_tddfpt_calculation) .AND. nlumo .NE. 0
         compute_lumos = compute_lumos .OR. p_loc_lumo

         DO ispin = 1, dft_control%nspins
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, homo=homo, nmo=nmo)
            compute_lumos = compute_lumos .AND. homo == nmo
         END DO

         IF (do_mo_cubes .AND. .NOT. compute_lumos) THEN

            nlumo = section_get_ival(dft_section, "PRINT%MO_CUBES%NLUMO")
            DO ispin = 1, dft_control%nspins

               CALL get_mo_set(mo_set=mos(ispin)%mo_set, homo=homo, nmo=nmo, eigenvalues=mo_eigenvalues)
               IF (nlumo > nmo - homo) THEN
                  ! this case not yet implemented
               ELSE
                  IF (nlumo .EQ. -1) THEN
                     nlumo = nmo - homo
                  END IF
                  IF (output_unit > 0) WRITE (output_unit, *) " "
                  IF (output_unit > 0) WRITE (output_unit, *) " Lowest eigenvalues of the unoccupied subspace spin ", ispin
                  IF (output_unit > 0) WRITE (output_unit, *) "---------------------------------------------"
                  IF (output_unit > 0) WRITE (output_unit, '(4(1X,1F16.8))') mo_eigenvalues(homo + 1:homo + nlumo)

                  ! Prints the cube files of UNOCCUPIED ORBITALS
                  CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
                  CALL qs_scf_post_unocc_cubes(input, dft_section, dft_control, logger, qs_env, &
                                               mo_coeff, wf_g, wf_r, particles, nlumo, homo, ispin, lumo=homo + 1)
               END IF
            END DO

         END IF

         IF (compute_lumos) THEN
            check_write = .TRUE.
            min_lumos = nlumo
            IF (nlumo == 0) check_write = .FALSE.
            IF (p_loc_lumo) THEN
               do_homo = .FALSE.
               CALL qs_loc_env_create(qs_loc_env_lumo)
               CALL qs_loc_control_init(qs_loc_env_lumo, localize_section, do_homo=do_homo)
               min_lumos = MAX(MAXVAL(qs_loc_env_lumo%localized_wfn_control%loc_states(:, :)), nlumo)
            END IF

            ALLOCATE (unoccupied_orbs(dft_control%nspins))
            ALLOCATE (unoccupied_evals(dft_control%nspins))
            CALL make_lumo_gpw(qs_env, scf_env, unoccupied_orbs, unoccupied_evals, min_lumos, nlumos)
            lumo_ptr => unoccupied_orbs
            DO ispin = 1, dft_control%nspins
               has_lumo = .TRUE.
               homo_lumo(ispin, 2) = unoccupied_evals(ispin)%array(1)
               CALL get_mo_set(mo_set=mos(ispin)%mo_set, homo=homo)
               IF (check_write) THEN
                  IF (p_loc_lumo .AND. nlumo .NE. -1) nlumos = MIN(nlumo, nlumos)
                  ! Prints the cube files of UNOCCUPIED ORBITALS
                  CALL qs_scf_post_unocc_cubes(input, dft_section, dft_control, logger, qs_env, &
                                               unoccupied_orbs(ispin)%matrix, wf_g, wf_r, particles, nlumos, homo, ispin)
               END IF
            END DO

            ! Save the info for tddfpt calculation
            IF (dft_control%do_tddfpt_calculation) THEN
               ALLOCATE (dft_control%tddfpt_control%lumos_eigenvalues(nlumos, dft_control%nspins))
               DO ispin = 1, dft_control%nspins
                  dft_control%tddfpt_control%lumos_eigenvalues(1:nlumos, ispin) = &
                     unoccupied_evals(ispin)%array(1:nlumos)
               END DO
               dft_control%tddfpt_control%lumos => unoccupied_orbs
            END IF

            IF (p_loc_lumo) THEN
               ALLOCATE (lumo_localized(dft_control%nspins))
               DO ispin = 1, dft_control%nspins
                  CALL cp_fm_create(lumo_localized(ispin)%matrix, unoccupied_orbs(ispin)%matrix%matrix_struct)
                  CALL cp_fm_to_fm(unoccupied_orbs(ispin)%matrix, lumo_localized(ispin)%matrix)
               END DO
               CALL qs_loc_init(qs_env, qs_loc_env_lumo, localize_section, lumo_localized, do_homo, do_mo_cubes, &
                                evals=unoccupied_evals)
               CALL qs_loc_env_init(qs_loc_env_lumo, qs_loc_env_lumo%localized_wfn_control, qs_env, &
                                    loc_coeff=unoccupied_orbs)
               CALL get_localization_info(qs_env, qs_loc_env_lumo, localize_section, &
                                          lumo_localized, wf_r, wf_g, particles, &
                                          unoccupied_orbs, unoccupied_evals, marked_states)
               CALL loc_write_restart(qs_loc_env_lumo, loc_print_section, mos, homo_localized, do_homo, &
                                      evals=unoccupied_evals)
               lumo_ptr => lumo_localized
            END IF
         ENDIF

         IF (has_homo .AND. has_lumo) THEN
            IF (output_unit > 0) WRITE (output_unit, *) " "
            DO ispin = 1, dft_control%nspins
               IF (.NOT. scf_control%smear%do_smear) THEN
                  gap = homo_lumo(ispin, 2) - homo_lumo(ispin, 1)
                  IF (output_unit > 0) WRITE (output_unit, '(T2,A,F12.6)') &
                     "HOMO - LUMO gap [eV] :", gap*evolt
               END IF
            ENDDO
         ENDIF
      ENDIF

      ! Deallocate grids needed to compute wavefunctions
      IF (((do_mo_cubes .OR. do_wannier_cubes) .AND. (nlumo /= 0 .OR. nhomo /= 0)) .OR. p_loc) THEN
         CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_g%pw)
      END IF

      ! Destroy the localization environment
      IF (.NOT. do_kpoints) THEN
         IF (p_loc_homo) CALL qs_loc_env_destroy(qs_loc_env_homo)
         IF (p_loc_lumo) CALL qs_loc_env_destroy(qs_loc_env_lumo)
      END IF

      ! generate a mix of wfns, and write to a restart
      IF (do_kpoints) THEN
         ! nothing at the moment, not implemented
      ELSE
         CALL get_qs_env(qs_env, matrix_s=matrix_s)
         CALL wfn_mix(mos, particle_set, dft_section, qs_kind_set, &
                      lumo_ptr, scf_env, matrix_s, output_unit, marked_states)
         IF (p_loc_lumo) CALL cp_fm_vect_dealloc(lumo_localized)
      END IF
      IF (ASSOCIATED(marked_states)) THEN
         DEALLOCATE (marked_states)
      END IF

      ! This is just a deallocation for printing MO_CUBES or TDDFPT
      IF (.NOT. do_kpoints) THEN
         IF (compute_lumos) THEN
            DO ispin = 1, dft_control%nspins
               DEALLOCATE (unoccupied_evals(ispin)%array)
               IF (.NOT. dft_control%do_tddfpt_calculation) &
                  CALL cp_fm_release(unoccupied_orbs(ispin)%matrix)
            ENDDO
            DEALLOCATE (unoccupied_evals)
            IF (.NOT. dft_control%do_tddfpt_calculation) THEN
               DEALLOCATE (unoccupied_orbs)
            END IF
         ENDIF
      ENDIF

      !stm images
      IF (do_stm) THEN
         IF (do_kpoints) THEN
            CPWARN("STM not implemented for k-point calculations!")
         ELSE
            IF (nlumo_stm > 0) THEN
               ALLOCATE (unoccupied_orbs_stm(dft_control%nspins))
               ALLOCATE (unoccupied_evals_stm(dft_control%nspins))
               CALL make_lumo_gpw(qs_env, scf_env, unoccupied_orbs_stm, unoccupied_evals_stm, &
                                  nlumo_stm, nlumos)
            END IF

            CALL th_stm_image(qs_env, stm_section, particles, unoccupied_orbs_stm, &
                              unoccupied_evals_stm)

            IF (nlumo_stm > 0) THEN
               DO ispin = 1, dft_control%nspins
                  DEALLOCATE (unoccupied_evals_stm(ispin)%array)
                  CALL cp_fm_release(unoccupied_orbs_stm(ispin)%matrix)
               ENDDO
               DEALLOCATE (unoccupied_evals_stm)
               DEALLOCATE (unoccupied_orbs_stm)
            END IF
         END IF
      END IF

      ! Print coherent X-ray diffraction spectrum
      CALL qs_scf_post_xray(input, dft_section, logger, qs_env, output_unit)

      ! Calculation of Electric Field Gradients
      CALL qs_scf_post_efg(input, logger, qs_env)

      ! Calculation of ET
      CALL qs_scf_post_et(input, qs_env, dft_control)

      ! Calculation of EPR Hyperfine Coupling Tensors
      CALL qs_scf_post_epr(input, logger, qs_env)

      ! Calculation of properties needed for BASIS_MOLOPT optimizations
      CALL qs_scf_post_molopt(input, logger, qs_env)

      ! Calculate ELF
      CALL qs_scf_post_elf(input, logger, qs_env)

      ! Use Wannier90 interface
      CALL wannier90_interface(input, logger, qs_env)

      ! Do active space calculation
      CALL active_space_main(input, logger, qs_env)

      CALL timestop(handle)

   END SUBROUTINE scf_post_calculation_gpw

! **************************************************************************************************
!> \brief Gets the lumos, and eigenvalues for the lumos
!> \param qs_env ...
!> \param scf_env ...
!> \param unoccupied_orbs ...
!> \param unoccupied_evals ...
!> \param nlumo ...
!> \param nlumos ...
! **************************************************************************************************
   SUBROUTINE make_lumo_gpw(qs_env, scf_env, unoccupied_orbs, unoccupied_evals, nlumo, nlumos)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: unoccupied_orbs
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: unoccupied_evals
      INTEGER, INTENT(IN)                                :: nlumo
      INTEGER, INTENT(OUT)                               :: nlumos

      CHARACTER(len=*), PARAMETER :: routineN = 'make_lumo_gpw', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, homo, ispin, n, nao, nmo, &
                                                            output_unit
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_rmpv, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(preconditioner_type), POINTER                 :: local_preconditioner
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (mos, ks_rmpv, scf_control, dft_control, admm_env, para_env, blacs_env)
      CALL get_qs_env(qs_env, &
                      mos=mos, &
                      matrix_ks=ks_rmpv, &
                      scf_control=scf_control, &
                      dft_control=dft_control, &
                      matrix_s=matrix_s, &
                      admm_env=admm_env, &
                      para_env=para_env, &
                      blacs_env=blacs_env)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      DO ispin = 1, dft_control%nspins
         NULLIFY (unoccupied_orbs(ispin)%matrix)
         NULLIFY (unoccupied_evals(ispin)%array)
         ! Always write eigenvalues
         IF (output_unit > 0) WRITE (output_unit, *) " "
         IF (output_unit > 0) WRITE (output_unit, *) " Lowest Eigenvalues of the unoccupied subspace spin ", ispin
         IF (output_unit > 0) WRITE (output_unit, FMT='(1X,A)') "-----------------------------------------------------"
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, homo=homo, nao=nao, nmo=nmo)
         CALL cp_fm_get_info(mo_coeff, nrow_global=n)
         nlumos = MAX(1, MIN(nlumo, nao - nmo))
         IF (nlumo == -1) nlumos = nao - nmo
         ALLOCATE (unoccupied_evals(ispin)%array(nlumos))
         CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=blacs_env, &
                                  nrow_global=n, ncol_global=nlumos)
         CALL cp_fm_create(unoccupied_orbs(ispin)%matrix, fm_struct_tmp, name="lumos")
         CALL cp_fm_struct_release(fm_struct_tmp)
         CALL cp_fm_init_random(unoccupied_orbs(ispin)%matrix, nlumos)

         ! the full_all preconditioner makes not much sense for lumos search
         NULLIFY (local_preconditioner)
         IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
            local_preconditioner => scf_env%ot_preconditioner(1)%preconditioner
            ! this one can for sure not be right (as it has to match a given C0)
            IF (local_preconditioner%in_use == ot_precond_full_all) THEN
               NULLIFY (local_preconditioner)
            ENDIF
         ENDIF

         ! ** If we do ADMM, we add have to modify the kohn-sham matrix
         IF (dft_control%do_admm) THEN
            CALL admm_correct_for_eigenvalues(ispin, admm_env, ks_rmpv(ispin)%matrix)
         END IF

         CALL ot_eigensolver(matrix_h=ks_rmpv(ispin)%matrix, matrix_s=matrix_s(1)%matrix, &
                             matrix_c_fm=unoccupied_orbs(ispin)%matrix, &
                             matrix_orthogonal_space_fm=mo_coeff, &
                             eps_gradient=scf_control%eps_lumos, &
                             preconditioner=local_preconditioner, &
                             iter_max=scf_control%max_iter_lumos, &
                             size_ortho_space=nmo)

         CALL calculate_subspace_eigenvalues(unoccupied_orbs(ispin)%matrix, ks_rmpv(ispin)%matrix, &
                                             unoccupied_evals(ispin)%array, scr=output_unit, &
                                             ionode=output_unit > 0)

         ! ** If we do ADMM, we restore the original kohn-sham matrix
         IF (dft_control%do_admm) THEN
            CALL admm_uncorrect_for_eigenvalues(ispin, admm_env, ks_rmpv(ispin)%matrix)
         END IF

      END DO

      CALL timestop(handle)

   END SUBROUTINE make_lumo_gpw
! **************************************************************************************************
!> \brief Computes and Prints Atomic Charges with several methods
!> \param input ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
! **************************************************************************************************
   SUBROUTINE qs_scf_post_charges(input, logger, qs_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_charges', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, print_level, unit_nr
      LOGICAL                                            :: do_kpoints, print_it
      TYPE(section_vals_type), POINTER                   :: density_fit_section, print_key

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, do_kpoints=do_kpoints)

      ! Mulliken charges require no further computation and are printed from write_mo_free_results

      ! Compute the Lowdin charges
      print_key => section_vals_get_subs_vals(input, "DFT%PRINT%LOWDIN")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%LOWDIN", extension=".lowdin", &
                                        log_filename=.FALSE.)
         print_level = 1
         CALL section_vals_val_get(print_key, "PRINT_GOP", l_val=print_it)
         IF (print_it) print_level = 2
         CALL section_vals_val_get(print_key, "PRINT_ALL", l_val=print_it)
         IF (print_it) print_level = 3
         IF (do_kpoints) THEN
            CPWARN("Lowdin charges not implemented for k-point calculations!")
         ELSE
            CALL lowdin_population_analysis(qs_env, unit_nr, print_level)
         END IF
         CALL cp_print_key_finished_output(unit_nr, logger, input, "DFT%PRINT%LOWDIN")
      END IF

      ! Compute the RESP charges
      CALL resp_fit(qs_env)

      ! Compute the Density Derived Atomic Point charges with the Bloechl scheme
      print_key => section_vals_get_subs_vals(input, "PROPERTIES%FIT_CHARGE")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, input, "PROPERTIES%FIT_CHARGE", extension=".Fitcharge", &
                                        log_filename=.FALSE.)
         density_fit_section => section_vals_get_subs_vals(input, "DFT%DENSITY_FITTING")
         CALL get_ddapc(qs_env, .FALSE., density_fit_section, iwc=unit_nr)
         CALL cp_print_key_finished_output(unit_nr, logger, input, "PROPERTIES%FIT_CHARGE")
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_charges

! **************************************************************************************************
!> \brief Computes and prints the Cube Files for MO
!> \param input ...
!> \param dft_section ...
!> \param dft_control ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param mo_coeff ...
!> \param wf_g ...
!> \param wf_r ...
!> \param particles ...
!> \param homo ...
!> \param ispin ...
! **************************************************************************************************
   SUBROUTINE qs_scf_post_occ_cubes(input, dft_section, dft_control, logger, qs_env, &
                                    mo_coeff, wf_g, wf_r, particles, homo, ispin)
      TYPE(section_vals_type), POINTER                   :: input, dft_section
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(pw_p_type)                                    :: wf_g, wf_r
      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER, INTENT(IN)                                :: homo, ispin

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_occ_cubes', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename, my_pos_cube, title
      INTEGER                                            :: handle, i, ir, ivector, n_rep, nhomo, &
                                                            nlist, unit_nr
      INTEGER, DIMENSION(:), POINTER                     :: list, list_index
      LOGICAL                                            :: append_cube, mpi_io
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (list_index)

      IF (BTEST(cp_print_key_should_output(logger%iter_info, dft_section, "PRINT%MO_CUBES") &
                , cp_p_file) .AND. section_get_lval(dft_section, "PRINT%MO_CUBES%WRITE_CUBE")) THEN
         nhomo = section_get_ival(dft_section, "PRINT%MO_CUBES%NHOMO")
         append_cube = section_get_lval(dft_section, "PRINT%MO_CUBES%APPEND")
         my_pos_cube = "REWIND"
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         END IF
         CALL section_vals_val_get(dft_section, "PRINT%MO_CUBES%HOMO_LIST", n_rep_val=n_rep)
         IF (n_rep > 0) THEN ! write the cubes of the list
            nlist = 0
            DO ir = 1, n_rep
               NULLIFY (list)
               CALL section_vals_val_get(dft_section, "PRINT%MO_CUBES%HOMO_LIST", i_rep_val=ir, &
                                         i_vals=list)
               IF (ASSOCIATED(list)) THEN
                  CALL reallocate(list_index, 1, nlist + SIZE(list))
                  DO i = 1, SIZE(list)
                     list_index(i + nlist) = list(i)
                  END DO
                  nlist = nlist + SIZE(list)
               END IF
            END DO
         ELSE

            IF (nhomo == -1) nhomo = homo
            nlist = homo - MAX(1, homo - nhomo + 1) + 1
            ALLOCATE (list_index(nlist))
            DO i = 1, nlist
               list_index(i) = MAX(1, homo - nhomo + 1) + i - 1
            END DO
         END IF
         DO i = 1, nlist
            ivector = list_index(i)
            CALL get_qs_env(qs_env=qs_env, &
                            atomic_kind_set=atomic_kind_set, &
                            qs_kind_set=qs_kind_set, &
                            cell=cell, &
                            particle_set=particle_set, &
                            pw_env=pw_env)
            CALL calculate_wavefunction(mo_coeff, ivector, wf_r, wf_g, atomic_kind_set, qs_kind_set, &
                                        cell, dft_control, particle_set, pw_env)
            WRITE (filename, '(a4,I5.5,a1,I1.1)') "WFN_", ivector, "_", ispin
            mpi_io = .TRUE.
            unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%MO_CUBES", extension=".cube", &
                                           middle_name=TRIM(filename), file_position=my_pos_cube, log_filename=.FALSE., &
                                           mpi_io=mpi_io)
            WRITE (title, *) "WAVEFUNCTION ", ivector, " spin ", ispin, " i.e. HOMO - ", ivector - homo
            CALL cp_pw_to_cube(wf_r%pw, unit_nr, title, particles=particles, &
                               stride=section_get_ivals(dft_section, "PRINT%MO_CUBES%STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr, logger, input, "DFT%PRINT%MO_CUBES", mpi_io=mpi_io)
         ENDDO
         IF (ASSOCIATED(list_index)) DEALLOCATE (list_index)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_occ_cubes

! **************************************************************************************************
!> \brief Computes and prints the Cube Files for MO
!> \param input ...
!> \param dft_section ...
!> \param dft_control ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param unoccupied_orbs ...
!> \param wf_g ...
!> \param wf_r ...
!> \param particles ...
!> \param nlumos ...
!> \param homo ...
!> \param ispin ...
!> \param lumo ...
! **************************************************************************************************
   SUBROUTINE qs_scf_post_unocc_cubes(input, dft_section, dft_control, logger, qs_env, &
                                      unoccupied_orbs, wf_g, wf_r, particles, nlumos, homo, ispin, lumo)
      TYPE(section_vals_type), POINTER                   :: input, dft_section
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), POINTER                          :: unoccupied_orbs
      TYPE(pw_p_type)                                    :: wf_g, wf_r
      TYPE(particle_list_type), POINTER                  :: particles
      INTEGER, INTENT(IN)                                :: nlumos, homo, ispin
      INTEGER, INTENT(IN), OPTIONAL                      :: lumo

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_unocc_cubes', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename, my_pos_cube, title
      INTEGER                                            :: handle, ifirst, index_mo, ivector, &
                                                            unit_nr
      LOGICAL                                            :: append_cube, mpi_io
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      IF (BTEST(cp_print_key_should_output(logger%iter_info, dft_section, "PRINT%MO_CUBES"), cp_p_file) &
          .AND. section_get_lval(dft_section, "PRINT%MO_CUBES%WRITE_CUBE")) THEN
         NULLIFY (qs_kind_set, particle_set, pw_env, cell)
         append_cube = section_get_lval(dft_section, "PRINT%MO_CUBES%APPEND")
         my_pos_cube = "REWIND"
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         END IF
         ifirst = 1
         IF (PRESENT(lumo)) ifirst = lumo
         DO ivector = ifirst, ifirst + nlumos - 1
            CALL get_qs_env(qs_env=qs_env, &
                            atomic_kind_set=atomic_kind_set, &
                            qs_kind_set=qs_kind_set, &
                            cell=cell, &
                            particle_set=particle_set, &
                            pw_env=pw_env)
            CALL calculate_wavefunction(unoccupied_orbs, ivector, wf_r, wf_g, atomic_kind_set, &
                                        qs_kind_set, cell, dft_control, particle_set, pw_env)

            IF (ifirst == 1) THEN
               index_mo = homo + ivector
            ELSE
               index_mo = ivector
            END IF
            WRITE (filename, '(a4,I5.5,a1,I1.1)') "WFN_", index_mo, "_", ispin
            mpi_io = .TRUE.

            unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%MO_CUBES", extension=".cube", &
                                           middle_name=TRIM(filename), file_position=my_pos_cube, log_filename=.FALSE., &
                                           mpi_io=mpi_io)
            WRITE (title, *) "WAVEFUNCTION ", index_mo, " spin ", ispin, " i.e. LUMO + ", ifirst + ivector - 2
            CALL cp_pw_to_cube(wf_r%pw, unit_nr, title, particles=particles, &
                               stride=section_get_ivals(dft_section, "PRINT%MO_CUBES%STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr, logger, input, "DFT%PRINT%MO_CUBES", mpi_io=mpi_io)
         ENDDO
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_unocc_cubes

! **************************************************************************************************
!> \brief Computes and prints electric moments
!> \param input ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE qs_scf_post_moments(input, logger, qs_env, output_unit)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_moments', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: handle, maxmom, reference, unit_nr
      LOGICAL                                            :: com_nl, do_kpoints, magnetic, periodic, &
                                                            second_ref_point, vel_reprs
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      print_key => section_vals_get_subs_vals(section_vals=input, &
                                              subsection_name="DFT%PRINT%MOMENTS")

      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN

         maxmom = section_get_ival(section_vals=input, &
                                   keyword_name="DFT%PRINT%MOMENTS%MAX_MOMENT")
         periodic = section_get_lval(section_vals=input, &
                                     keyword_name="DFT%PRINT%MOMENTS%PERIODIC")
         reference = section_get_ival(section_vals=input, &
                                      keyword_name="DFT%PRINT%MOMENTS%REFERENCE")
         magnetic = section_get_lval(section_vals=input, &
                                     keyword_name="DFT%PRINT%MOMENTS%MAGNETIC")
         vel_reprs = section_get_lval(section_vals=input, &
                                      keyword_name="DFT%PRINT%MOMENTS%VEL_REPRS")
         com_nl = section_get_lval(section_vals=input, &
                                   keyword_name="DFT%PRINT%MOMENTS%COM_NL")
         second_ref_point = section_get_lval(section_vals=input, &
                                             keyword_name="DFT%PRINT%MOMENTS%SECOND_REFERENCE_POINT")

         NULLIFY (ref_point)
         CALL section_vals_val_get(input, "DFT%PRINT%MOMENTS%REF_POINT", r_vals=ref_point)
         unit_nr = cp_print_key_unit_nr(logger=logger, basis_section=input, &
                                        print_key_path="DFT%PRINT%MOMENTS", extension=".dat", &
                                        middle_name="moments", log_filename=.FALSE.)

         IF (output_unit > 0) THEN
            IF (unit_nr /= output_unit) THEN
               INQUIRE (UNIT=unit_nr, NAME=filename)
               WRITE (UNIT=output_unit, FMT="(/,T2,A,2(/,T3,A),/)") &
                  "MOMENTS", "The electric/magnetic moments are written to file:", &
                  TRIM(filename)
            ELSE
               WRITE (UNIT=output_unit, FMT="(/,T2,A)") "ELECTRIC/MAGNETIC MOMENTS"
            END IF
         END IF

         CALL get_qs_env(qs_env, do_kpoints=do_kpoints)

         IF (do_kpoints) THEN
            CPWARN("Moments not implemented for k-point calculations!")
         ELSE
            IF (periodic) THEN
               CALL qs_moment_berry_phase(qs_env, magnetic, maxmom, reference, ref_point, unit_nr)
            ELSE
               CALL qs_moment_locop(qs_env, magnetic, maxmom, reference, ref_point, unit_nr, vel_reprs, com_nl)
            END IF
         END IF

         CALL cp_print_key_finished_output(unit_nr=unit_nr, logger=logger, &
                                           basis_section=input, print_key_path="DFT%PRINT%MOMENTS")

         IF (second_ref_point) THEN
            reference = section_get_ival(section_vals=input, &
                                         keyword_name="DFT%PRINT%MOMENTS%REFERENCE_2")

            NULLIFY (ref_point)
            CALL section_vals_val_get(input, "DFT%PRINT%MOMENTS%REF_POINT_2", r_vals=ref_point)
            unit_nr = cp_print_key_unit_nr(logger=logger, basis_section=input, &
                                           print_key_path="DFT%PRINT%MOMENTS", extension=".dat", &
                                           middle_name="moments_refpoint_2", log_filename=.FALSE.)

            IF (output_unit > 0) THEN
               IF (unit_nr /= output_unit) THEN
                  INQUIRE (UNIT=unit_nr, NAME=filename)
                  WRITE (UNIT=output_unit, FMT="(/,T2,A,2(/,T3,A),/)") &
                     "MOMENTS", "The electric/magnetic moments for the second reference point are written to file:", &
                     TRIM(filename)
               ELSE
                  WRITE (UNIT=output_unit, FMT="(/,T2,A)") "ELECTRIC/MAGNETIC MOMENTS"
               END IF
            END IF
            IF (do_kpoints) THEN
               CPWARN("Moments not implemented for k-point calculations!")
            ELSE
               IF (periodic) THEN
                  CALL qs_moment_berry_phase(qs_env, magnetic, maxmom, reference, ref_point, unit_nr)
               ELSE
                  CALL qs_moment_locop(qs_env, magnetic, maxmom, reference, ref_point, unit_nr, vel_reprs, com_nl)
               END IF
            END IF
            CALL cp_print_key_finished_output(unit_nr=unit_nr, logger=logger, &
                                              basis_section=input, print_key_path="DFT%PRINT%MOMENTS")
         ENDIF

      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_moments

! **************************************************************************************************
!> \brief Computes and prints the X-ray diffraction spectrum.
!> \param input ...
!> \param dft_section ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param output_unit ...
! **************************************************************************************************
   SUBROUTINE qs_scf_post_xray(input, dft_section, logger, qs_env, output_unit)

      TYPE(section_vals_type), POINTER                   :: input, dft_section
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: output_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_post_xray', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: handle, unit_nr
      REAL(KIND=dp)                                      :: q_max
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      print_key => section_vals_get_subs_vals(section_vals=input, &
                                              subsection_name="DFT%PRINT%XRAY_DIFFRACTION_SPECTRUM")

      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         q_max = section_get_rval(section_vals=dft_section, &
                                  keyword_name="PRINT%XRAY_DIFFRACTION_SPECTRUM%Q_MAX")
         unit_nr = cp_print_key_unit_nr(logger=logger, &
                                        basis_section=input, &
                                        print_key_path="DFT%PRINT%XRAY_DIFFRACTION_SPECTRUM", &
                                        extension=".dat", &
                                        middle_name="xrd", &
                                        log_filename=.FALSE.)
         IF (output_unit > 0) THEN
            INQUIRE (UNIT=unit_nr, NAME=filename)
            WRITE (UNIT=output_unit, FMT="(/,/,T2,A)") &
               "X-RAY DIFFRACTION SPECTRUM"
            IF (unit_nr /= output_unit) THEN
               WRITE (UNIT=output_unit, FMT="(/,T3,A,/,/,T3,A,/)") &
                  "The coherent X-ray diffraction spectrum is written to the file:", &
                  TRIM(filename)
            END IF
         END IF
         CALL xray_diffraction_spectrum(qs_env=qs_env, &
                                        unit_number=unit_nr, &
                                        q_max=q_max)
         CALL cp_print_key_finished_output(unit_nr=unit_nr, &
                                           logger=logger, &
                                           basis_section=input, &
                                           print_key_path="DFT%PRINT%XRAY_DIFFRACTION_SPECTRUM")
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_xray

! **************************************************************************************************
!> \brief Computes and prints Electric Field Gradient
!> \param input ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
! **************************************************************************************************
   SUBROUTINE qs_scf_post_efg(input, logger, qs_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_efg', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      print_key => section_vals_get_subs_vals(section_vals=input, &
                                              subsection_name="DFT%PRINT%ELECTRIC_FIELD_GRADIENT")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), &
                cp_p_file)) THEN
         CALL qs_efg_calc(qs_env=qs_env)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_efg

! **************************************************************************************************
!> \brief Computes the Electron Transfer Coupling matrix element
!> \param input ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param dft_control ...
! **************************************************************************************************
   SUBROUTINE qs_scf_post_et(input, qs_env, dft_control)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dft_control_type), POINTER                    :: dft_control

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_et', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin
      LOGICAL                                            :: do_et
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: my_mos
      TYPE(section_vals_type), POINTER                   :: et_section

      CALL timeset(routineN, handle)

      do_et = .FALSE.
      et_section => section_vals_get_subs_vals(input, "PROPERTIES%ET_COUPLING")
      CALL section_vals_get(et_section, explicit=do_et)
      IF (do_et) THEN
         IF (qs_env%et_coupling%first_run) THEN
            NULLIFY (my_mos)
            ALLOCATE (my_mos(dft_control%nspins))
            ALLOCATE (qs_env%et_coupling%et_mo_coeff(dft_control%nspins))
            DO ispin = 1, dft_control%nspins
               NULLIFY (my_mos(ispin)%matrix)
               CALL cp_fm_create(matrix=my_mos(ispin)%matrix, &
                                 matrix_struct=qs_env%mos(ispin)%mo_set%mo_coeff%matrix_struct, &
                                 name="FIRST_RUN_COEFF"//TRIM(ADJUSTL(cp_to_string(ispin)))//"MATRIX")
               CALL cp_fm_to_fm(qs_env%mos(ispin)%mo_set%mo_coeff, &
                                my_mos(ispin)%matrix)
            END DO
            CALL set_et_coupling_type(qs_env%et_coupling, et_mo_coeff=my_mos)
            DEALLOCATE (my_mos)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_et

! **************************************************************************************************
!> \brief compute the electron localization function
!>
!> \param input ...
!> \param logger ...
!> \param qs_env ...
!> \par History
!>      2012-07 Created [MI]
! **************************************************************************************************
   SUBROUTINE qs_scf_post_elf(input, logger, qs_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_elf', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename, mpi_filename, my_pos_cube, &
                                                            title
      INTEGER                                            :: handle, ispin, output_unit, unit_nr
      LOGICAL                                            :: append_cube, gapw, mpi_io
      REAL(dp)                                           :: rho_cutoff
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: elf_r
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: elf_section

      CALL timeset(routineN, handle)
      output_unit = cp_logger_get_default_io_unit(logger)

      elf_section => section_vals_get_subs_vals(input, "DFT%PRINT%ELF_CUBE")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%PRINT%ELF_CUBE"), cp_p_file)) THEN

         NULLIFY (dft_control, pw_env, auxbas_pw_pool, pw_pools, particles, subsys, elf_r)
         CALL get_qs_env(qs_env, dft_control=dft_control, pw_env=pw_env, subsys=subsys)
         CALL qs_subsys_get(subsys, particles=particles)

         gapw = dft_control%qs_control%gapw
         IF (.NOT. gapw) THEN
            ! allocate
            ALLOCATE (elf_r(dft_control%nspins))
            CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                            pw_pools=pw_pools)
            DO ispin = 1, dft_control%nspins
               CALL pw_pool_create_pw(auxbas_pw_pool, elf_r(ispin)%pw, &
                                      use_data=REALDATA3D, &
                                      in_space=REALSPACE)
               CALL pw_zero(elf_r(ispin)%pw)
            END DO

            IF (output_unit > 0) THEN
               WRITE (UNIT=output_unit, FMT="(/,T15,A,/)") &
                  " ----- ELF is computed on the real space grid -----"
            END IF
            rho_cutoff = section_get_rval(elf_section, "density_cutoff")
            CALL qs_elf_calc(qs_env, elf_r, rho_cutoff)

            ! write ELF into cube file
            append_cube = section_get_lval(elf_section, "APPEND")
            my_pos_cube = "REWIND"
            IF (append_cube) THEN
               my_pos_cube = "APPEND"
            END IF

            DO ispin = 1, dft_control%nspins
               WRITE (filename, '(a5,I1.1)') "ELF_S", ispin
               WRITE (title, *) "ELF spin ", ispin
               mpi_io = .TRUE.
               unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%ELF_CUBE", extension=".cube", &
                                              middle_name=TRIM(filename), file_position=my_pos_cube, log_filename=.FALSE., &
                                              mpi_io=mpi_io, fout=mpi_filename)
               IF (output_unit > 0) THEN
                  IF (.NOT. mpi_io) THEN
                     INQUIRE (UNIT=unit_nr, NAME=filename)
                  ELSE
                     filename = mpi_filename
                  END IF
                  WRITE (UNIT=output_unit, FMT="(/,T2,A,/,/,T2,A)") &
                     "ELF is written in cube file format to the file:", &
                     TRIM(filename)
               END IF

               CALL cp_pw_to_cube(elf_r(ispin)%pw, unit_nr, title, particles=particles, &
                                  stride=section_get_ivals(elf_section, "STRIDE"), mpi_io=mpi_io)
               CALL cp_print_key_finished_output(unit_nr, logger, input, "DFT%PRINT%ELF_CUBE", mpi_io=mpi_io)

               CALL pw_pool_give_back_pw(auxbas_pw_pool, elf_r(ispin)%pw)
            END DO

            ! deallocate
            DEALLOCATE (elf_r)

         ELSE
            ! not implemented
            CPWARN("ELF not implemented for GAPW calculations!!")

         END IF

      END IF ! print key

      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_elf

! **************************************************************************************************
!> \brief computes the condition number of the overlap matrix and
!>      prints the value of the total energy. This is needed
!>      for BASIS_MOLOPT optimizations
!> \param input ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \par History
!>      2007-07 Created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE qs_scf_post_molopt(input, logger, qs_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_molopt', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nao, unit_nr
      REAL(KIND=dp)                                      :: S_cond_number
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_fmstruct
      TYPE(cp_fm_type), POINTER                          :: fm_s, fm_work, mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      print_key => section_vals_get_subs_vals(section_vals=input, &
                                              subsection_name="DFT%PRINT%BASIS_MOLOPT_QUANTITIES")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), &
                cp_p_file)) THEN

         CALL get_qs_env(qs_env, energy=energy, matrix_s=matrix_s, mos=mos)

         ! set up the two needed full matrices, using mo_coeff as a template
         CALL get_mo_set(mo_set=mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao)
         CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct, &
                                  nrow_global=nao, ncol_global=nao, &
                                  template_fmstruct=mo_coeff%matrix_struct)
         CALL cp_fm_create(fm_s, matrix_struct=ao_ao_fmstruct, &
                           name="fm_s")
         CALL cp_fm_create(fm_work, matrix_struct=ao_ao_fmstruct, &
                           name="fm_work")
         CALL cp_fm_struct_release(ao_ao_fmstruct)
         ALLOCATE (eigenvalues(nao))

         CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, fm_s)
         CALL choose_eigv_solver(fm_s, fm_work, eigenvalues)

         CALL cp_fm_release(fm_s)
         CALL cp_fm_release(fm_work)

         S_cond_number = MAXVAL(ABS(eigenvalues))/MAX(MINVAL(ABS(eigenvalues)), EPSILON(0.0_dp))

         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%BASIS_MOLOPT_QUANTITIES", &
                                        extension=".molopt")

         IF (unit_nr > 0) THEN
            ! please keep this format fixed, needs to be grepable for molopt
            ! optimizations
            WRITE (unit_nr, '(T2,A28,2A25)') "", "Tot. Ener.", "S Cond. Numb."
            WRITE (unit_nr, '(T2,A28,2E25.17)') "BASIS_MOLOPT_QUANTITIES", energy%total, S_cond_number
         ENDIF

         CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                           "DFT%PRINT%BASIS_MOLOPT_QUANTITIES")

      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_molopt

! **************************************************************************************************
!> \brief Dumps EPR
!> \param input ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
! **************************************************************************************************
   SUBROUTINE qs_scf_post_epr(input, logger, qs_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_epr', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      print_key => section_vals_get_subs_vals(section_vals=input, &
                                              subsection_name="DFT%PRINT%HYPERFINE_COUPLING_TENSOR")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), &
                cp_p_file)) THEN
         CALL qs_epr_hyp_calc(qs_env=qs_env)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_epr

! **************************************************************************************************
!> \brief Interface routine to trigger writing of results available from normal
!>        SCF. Can write MO-dependent and MO free results (needed for call from
!>        the linear scaling code)
!> \param qs_env the qs_env in which the qs_env lives
!> \param scf_env ...
! **************************************************************************************************
   SUBROUTINE write_available_results(qs_env, scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), OPTIONAL, POINTER           :: scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'write_available_results', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! those properties that require MOs (not suitable density matrix based methods)
      CALL write_mo_dependent_results(qs_env, scf_env)

      ! those that depend only on the density matrix, they should be linear scaling in their implementation
      CALL write_mo_free_results(qs_env)

      CALL timestop(handle)

   END SUBROUTINE write_available_results

! **************************************************************************************************
!> \brief Write QS results available if MO's are present (if switched on through the print_keys)
!>        Writes only MO dependent results. Split is necessary as ls_scf does not
!>        provide MO's
!> \param qs_env the qs_env in which the qs_env lives
!> \param scf_env ...
! **************************************************************************************************
   SUBROUTINE write_mo_dependent_results(qs_env, scf_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), OPTIONAL, POINTER           :: scf_env

      CHARACTER(len=*), PARAMETER :: routineN = 'write_mo_dependent_results', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, homo, ik, ispin, nmo, output_unit
      LOGICAL                                            :: all_equal, do_kpoints
      REAL(KIND=dp)                                      :: maxocc, s_square, s_square_ideal, &
                                                            total_abs_spin_dens
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues, occupation_numbers
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_rmpv, matrix_s
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_deriv
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(mo_set_p_type), DIMENSION(:, :), POINTER      :: mos_k
      TYPE(molecule_type), POINTER                       :: molecule_set(:)
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: wf_r
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input, sprint_section

      CALL timeset(routineN, handle)
      NULLIFY (cell, dft_control, pw_env, auxbas_pw_pool, pw_pools, mo_coeff, &
               mo_coeff_deriv, mo_eigenvalues, mos, atomic_kind_set, qs_kind_set, &
               particle_set, rho, ks_rmpv, matrix_s, scf_control, dft_section, &
               molecule_set, input, particles, subsys, rho_r)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      molecule_set=molecule_set, &
                      atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      admm_env=admm_env, &
                      scf_control=scf_control, &
                      input=input, &
                      cell=cell, &
                      subsys=subsys)
      CALL qs_subsys_get(subsys, particles=particles)
      CALL get_qs_env(qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_r=rho_r)

      ! kpoints
      CALL get_qs_env(qs_env, do_kpoints=do_kpoints)

      ! *** if the dft_section tells you to do so, write last wavefunction to screen
      dft_section => section_vals_get_subs_vals(input, "DFT")
      IF (.NOT. qs_env%run_rtp) THEN
         IF (.NOT. do_kpoints) THEN
            CALL get_qs_env(qs_env, mos=mos)
            IF (dft_control%nspins == 2) THEN
               CALL write_mo_set(mos(1)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                                 dft_section, spin="ALPHA", last=.TRUE.)
               CALL write_mo_set(mos(2)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                                 dft_section, spin="BETA", last=.TRUE.)
            ELSE
               CALL write_mo_set(mos(1)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                                 dft_section, last=.TRUE.)
            END IF
            CALL get_qs_env(qs_env, matrix_ks=ks_rmpv)
            CALL write_dm_binary_restart(mos, dft_section, ks_rmpv)
            sprint_section => section_vals_get_subs_vals(dft_section, "PRINT%MO_MOLDEN")
            CALL write_mos_molden(mos, qs_kind_set, particle_set, sprint_section)
         ELSE ! *** if we have k points, let's print the eigenvalues at each of them
            CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)

            IF (kpoints%nkp /= 0) THEN
               DO ik = 1, SIZE(kpoints%kp_env)
                  mos_k => kpoints%kp_env(ik)%kpoint_env%mos
                  CPASSERT(ASSOCIATED(mos_k))

                  IF (dft_control%nspins == 2) THEN
                     CALL write_mo_set(mos_k(1, 1)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                                       dft_section, spin="ALPHA", last=.TRUE., kpt=ik)
                     CALL write_mo_set(mos_k(1, 2)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                                       dft_section, spin="BETA", last=.TRUE., kpt=ik)
                  ELSE
                     CALL write_mo_set(mos_k(1, 1)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                                       dft_section, last=.TRUE., kpt=ik)
                  END IF
               END DO
            END IF
         END IF

         ! *** at the end of scf print out dos
         IF (BTEST(cp_print_key_should_output(logger%iter_info, dft_section, "PRINT%DOS") &
                   , cp_p_file)) THEN
            IF (do_kpoints) THEN
               CALL get_qs_env(qs_env=qs_env, kpoints=kpoints)
               CALL calculate_dos_kp(kpoints, qs_env, dft_section)
            ELSE
               CALL get_qs_env(qs_env, mos=mos)
               CALL calculate_dos(mos, qs_env, dft_section)
            ENDIF
         ENDIF

         ! *** at the end of scf print out the projected dos per kind
         IF (BTEST(cp_print_key_should_output(logger%iter_info, dft_section, "PRINT%PDOS") &
                   , cp_p_file)) THEN
            IF (do_kpoints) THEN
               CPWARN("Projected density of states not implemented for k-points.")
            ELSE
               CALL get_qs_env(qs_env, &
                               mos=mos, &
                               matrix_ks=ks_rmpv)
               DO ispin = 1, dft_control%nspins

                  ! ** If we do ADMM, we add have to modify the kohn-sham matrix
                  IF (dft_control%do_admm) THEN
                     CALL admm_correct_for_eigenvalues(ispin, admm_env, ks_rmpv(ispin)%matrix)
                  END IF
                  IF (PRESENT(scf_env)) THEN
                     IF (scf_env%method == ot_method_nr) THEN
                        CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff, &
                                        eigenvalues=mo_eigenvalues)
                        IF (ASSOCIATED(qs_env%mo_derivs)) THEN
                           mo_coeff_deriv => qs_env%mo_derivs(ispin)%matrix
                        ELSE
                           mo_coeff_deriv => NULL()
                        ENDIF

                        CALL calculate_subspace_eigenvalues(mo_coeff, ks_rmpv(ispin)%matrix, mo_eigenvalues, &
                                                            do_rotation=.TRUE., &
                                                            co_rotate_dbcsr=mo_coeff_deriv)
                        CALL set_mo_occupation(mo_set=mos(ispin)%mo_set)

                     END IF
                  END IF
                  IF (dft_control%nspins == 2) THEN
                     CALL calculate_projected_dos(mos(ispin)%mo_set, atomic_kind_set, &
                                                  qs_kind_set, particle_set, qs_env, dft_section, ispin=ispin)
                  ELSE
                     CALL calculate_projected_dos(mos(ispin)%mo_set, atomic_kind_set, &
                                                  qs_kind_set, particle_set, qs_env, dft_section)
                  END IF

                  ! ** If we do ADMM, we add have to modify the kohn-sham matrix
                  IF (dft_control%do_admm) THEN
                     CALL admm_uncorrect_for_eigenvalues(ispin, admm_env, ks_rmpv(ispin)%matrix)
                  END IF

               END DO
            ENDIF
         ENDIF
      END IF

      !   *** Integrated absolute spin density and spin contamination ***
      IF (dft_control%nspins .EQ. 2) THEN
         CALL get_qs_env(qs_env, mos=mos)
         CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                         pw_pools=pw_pools)
         CALL pw_pool_create_pw(auxbas_pw_pool, wf_r%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_copy(rho_r(1)%pw, wf_r%pw)
         CALL pw_axpy(rho_r(2)%pw, wf_r%pw, alpha=-1._dp)
         total_abs_spin_dens = pw_integrate_function(wf_r%pw, oprt="ABS")
         IF (output_unit > 0) WRITE (UNIT=output_unit, FMT='(/,(T3,A,T61,F20.10))') &
            "Integrated absolute spin density  : ", total_abs_spin_dens
         CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r%pw)
         !
         ! XXX Fix Me XXX
         ! should be extended to the case where added MOs are present
         ! should be extended to the k-point case
         !
         IF (do_kpoints) THEN
            CPWARN("Spin contamination estimate not implemented for k-points.")
         ELSE
            all_equal = .TRUE.
            DO ispin = 1, dft_control%nspins
               CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
                               occupation_numbers=occupation_numbers, &
                               homo=homo, &
                               nmo=nmo, &
                               maxocc=maxocc)
               IF (nmo > 0) THEN
                  all_equal = all_equal .AND. &
                              (ALL(occupation_numbers(1:homo) == maxocc) .AND. &
                               ALL(occupation_numbers(homo + 1:nmo) == 0.0_dp))
               END IF
            END DO
            IF (.NOT. all_equal) THEN
               IF (output_unit > 0) WRITE (UNIT=output_unit, FMT="(T3,A)") &
                  "WARNING: S**2 computation does not yet treat fractional occupied orbitals"
            ELSE
               CALL get_qs_env(qs_env=qs_env, &
                               matrix_s=matrix_s, &
                               energy=energy)
               CALL compute_s_square(mos=mos, matrix_s=matrix_s, s_square=s_square, &
                                     s_square_ideal=s_square_ideal)
               IF (output_unit > 0) WRITE (UNIT=output_unit, FMT='(T3,A,T51,2F15.6)') &
                  "Ideal and single determinant S**2 : ", s_square_ideal, s_square
               energy%s_square = s_square
            END IF
         ENDIF
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE write_mo_dependent_results

! **************************************************************************************************
!> \brief Write QS results always available (if switched on through the print_keys)
!>        Can be called from ls_scf
!> \param qs_env the qs_env in which the qs_env lives
! **************************************************************************************************
   SUBROUTINE write_mo_free_results(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'write_mo_free_results', &
         routineP = moduleN//':'//routineN
      CHARACTER(len=1), DIMENSION(3), PARAMETER          :: cdir = (/"x", "y", "z"/)

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=default_path_length)                 :: filename, mpi_filename, my_pos_cube
      CHARACTER(LEN=default_string_length)               :: name
      INTEGER                                            :: after, handle, i, iat, id, ikind, img, &
                                                            iso, ispin, iw, l, nd(3), ngto, niso, &
                                                            nkind, np, nr, output_unit, &
                                                            print_level, unit_nr
      LOGICAL :: append_cube, do_kpoints, mpi_io, omit_headers, print_it, print_total_density, &
         rho_r_valid, write_ks, write_xc, xrd_interface
      REAL(KIND=dp)                                      :: q_max, rho_hard, rho_soft, rho_total, &
                                                            rho_total_rspace, udvol, volume
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: bfun
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: aedens, ccdens, ppdens
      REAL(KIND=dp), DIMENSION(3)                        :: dr
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hr
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_rmpv, matrix_vxc, rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: aux_g, aux_r, rho_elec_gspace, &
                                                            rho_elec_rspace, wf_r
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(pw_p_type), POINTER                           :: rho0_s_gs, rho_core, vee
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type), POINTER                             :: v_hartree_rspace
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set
      TYPE(rho_atom_type), POINTER                       :: rho_atom
      TYPE(section_vals_type), POINTER                   :: dft_section, input, print_key, xc_section

      CALL timeset(routineN, handle)
      NULLIFY (cell, dft_control, pw_env, auxbas_pw_pool, pw_pools, &
               atomic_kind_set, qs_kind_set, particle_set, rho, ks_rmpv, rho_ao, rho_r, &
               dft_section, xc_section, input, particles, subsys, matrix_vxc, v_hartree_rspace, vee)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      cell=cell, &
                      para_env=para_env, &
                      dft_control=dft_control, &
                      input=input, &
                      do_kpoints=do_kpoints, &
                      subsys=subsys)
      dft_section => section_vals_get_subs_vals(input, "DFT")
      CALL qs_subsys_get(subsys, particles=particles)

      ! Print the total density (electronic + core charge)
      CALL get_qs_env(qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_r=rho_r)
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%PRINT%TOT_DENSITY_CUBE"), cp_p_file)) THEN
         NULLIFY (rho_core, rho0_s_gs)
         append_cube = section_get_lval(input, "DFT%PRINT%TOT_DENSITY_CUBE%APPEND")
         my_pos_cube = "REWIND"
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         END IF

         CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, rho_core=rho_core, &
                         rho0_s_gs=rho0_s_gs)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                         pw_pools=pw_pools)
         CALL pw_pool_create_pw(auxbas_pw_pool, wf_r%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         IF (dft_control%qs_control%gapw) THEN
            CALL pw_transfer(rho0_s_gs%pw, wf_r%pw)
            IF (dft_control%qs_control%gapw_control%nopaw_as_gpw) THEN
               CALL pw_axpy(rho_core%pw, wf_r%pw)
            END IF
         ELSE
            CALL pw_transfer(rho_core%pw, wf_r%pw)
         END IF
         DO ispin = 1, dft_control%nspins
            CALL pw_axpy(rho_r(ispin)%pw, wf_r%pw)
         END DO
         filename = "TOTAL_DENSITY"
         mpi_io = .TRUE.
         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%TOT_DENSITY_CUBE", &
                                        extension=".cube", middle_name=TRIM(filename), file_position=my_pos_cube, &
                                        log_filename=.FALSE., mpi_io=mpi_io)
         CALL cp_pw_to_cube(wf_r%pw, unit_nr, "TOTAL DENSITY", &
                            particles=particles, &
                            stride=section_get_ivals(dft_section, "PRINT%TOT_DENSITY_CUBE%STRIDE"), mpi_io=mpi_io)
         CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                           "DFT%PRINT%TOT_DENSITY_CUBE", mpi_io=mpi_io)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, wf_r%pw)
      END IF

      ! Write cube file with electron density
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%PRINT%E_DENSITY_CUBE"), cp_p_file)) THEN
         CALL section_vals_val_get(dft_section, &
                                   keyword_name="PRINT%E_DENSITY_CUBE%TOTAL_DENSITY", &
                                   l_val=print_total_density)
         append_cube = section_get_lval(input, "DFT%PRINT%E_DENSITY_CUBE%APPEND")
         my_pos_cube = "REWIND"
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         END IF
         ! Write the info on core densities for the interface between cp2k and the XRD code
         ! together with the valence density they are used to compute the form factor (Fourier transform)
         xrd_interface = section_get_lval(input, "DFT%PRINT%E_DENSITY_CUBE%XRD_INTERFACE")
         IF (xrd_interface) THEN
            !cube file only contains soft density (GAPW)
            IF (dft_control%qs_control%gapw) print_total_density = .FALSE.

            filename = "ELECTRON_DENSITY"
            unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%E_DENSITY_CUBE", &
                                           extension=".xrd", middle_name=TRIM(filename), &
                                           file_position=my_pos_cube, log_filename=.FALSE.)
            ngto = section_get_ival(input, "DFT%PRINT%E_DENSITY_CUBE%NGAUSS")
            IF (output_unit > 0) THEN
               INQUIRE (UNIT=unit_nr, NAME=filename)
               WRITE (UNIT=output_unit, FMT="(/,T2,A,/,/,T2,A)") &
                  "The electron density (atomic part) is written to the file:", &
                  TRIM(filename)
            END IF

            xc_section => section_vals_get_subs_vals(input, "DFT%XC")
            nkind = SIZE(atomic_kind_set)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *) "Atomic (core) densities"
               WRITE (unit_nr, *) "Unit cell"
               WRITE (unit_nr, FMT="(3F20.12)") cell%hmat(1, 1), cell%hmat(1, 2), cell%hmat(1, 3)
               WRITE (unit_nr, FMT="(3F20.12)") cell%hmat(2, 1), cell%hmat(2, 2), cell%hmat(2, 3)
               WRITE (unit_nr, FMT="(3F20.12)") cell%hmat(3, 1), cell%hmat(3, 2), cell%hmat(3, 3)
               WRITE (unit_nr, *) "Atomic types"
               WRITE (unit_nr, *) nkind
            END IF
            ! calculate atomic density and core density
            ALLOCATE (ppdens(ngto, 2, nkind), aedens(ngto, 2, nkind), ccdens(ngto, 2, nkind))
            DO ikind = 1, nkind
               atomic_kind => atomic_kind_set(ikind)
               qs_kind => qs_kind_set(ikind)
               CALL get_atomic_kind(atomic_kind, name=name, element_symbol=element_symbol)
               CALL calculate_atomic_density(ppdens(:, :, ikind), atomic_kind, qs_kind, ngto, &
                                             iunit=output_unit, confine=.TRUE.)
               CALL calculate_atomic_density(aedens(:, :, ikind), atomic_kind, qs_kind, ngto, &
                                             iunit=output_unit, allelectron=.TRUE., confine=.TRUE.)
               ccdens(:, 1, ikind) = aedens(:, 1, ikind)
               ccdens(:, 2, ikind) = 0._dp
               CALL project_function_a(ccdens(1:ngto, 2, ikind), ccdens(1:ngto, 1, ikind), &
                                       ppdens(1:ngto, 2, ikind), ppdens(1:ngto, 1, ikind), 0)
               ccdens(:, 2, ikind) = aedens(:, 2, ikind) - ccdens(:, 2, ikind)
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, FMT="(I6,A10,A20)") ikind, TRIM(element_symbol), TRIM(name)
                  WRITE (unit_nr, FMT="(I6)") ngto
                  WRITE (unit_nr, *) "   Total density"
                  WRITE (unit_nr, FMT="(2G24.12)") (aedens(i, 1, ikind), aedens(i, 2, ikind), i=1, ngto)
                  WRITE (unit_nr, *) "    Core density"
                  WRITE (unit_nr, FMT="(2G24.12)") (ccdens(i, 1, ikind), ccdens(i, 2, ikind), i=1, ngto)
               END IF
               NULLIFY (atomic_kind)
            END DO

            IF (dft_control%qs_control%gapw) THEN
               CALL get_qs_env(qs_env=qs_env, rho_atom_set=rho_atom_set)

               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, *) "Coordinates and GAPW density"
               END IF
               np = particles%n_els
               DO iat = 1, np
                  CALL get_atomic_kind(particles%els(iat)%atomic_kind, kind_number=ikind)
                  CALL get_qs_kind(qs_kind_set(ikind), grid_atom=grid_atom)
                  rho_atom => rho_atom_set(iat)
                  IF (ASSOCIATED(rho_atom%rho_rad_h(1)%r_coef)) THEN
                     nr = SIZE(rho_atom%rho_rad_h(1)%r_coef, 1)
                     niso = SIZE(rho_atom%rho_rad_h(1)%r_coef, 2)
                  ELSE
                     nr = 0
                     niso = 0
                  ENDIF
                  CALL mp_sum(nr, para_env%group)
                  CALL mp_sum(niso, para_env%group)

                  ALLOCATE (bfun(nr, niso))
                  bfun = 0._dp
                  DO ispin = 1, dft_control%nspins
                     IF (ASSOCIATED(rho_atom%rho_rad_h(1)%r_coef)) THEN
                        bfun(:, :) = bfun + rho_atom%rho_rad_h(ispin)%r_coef - rho_atom%rho_rad_s(ispin)%r_coef
                     END IF
                  END DO
                  CALL mp_sum(bfun, para_env%group)
                  ccdens(:, 1, ikind) = ppdens(:, 1, ikind)
                  ccdens(:, 2, ikind) = 0._dp
                  IF (unit_nr > 0) THEN
                     WRITE (unit_nr, '(I10,I5,3f12.6)') iat, ikind, particles%els(iat)%r
                  END IF
                  DO iso = 1, niso
                     l = indso(1, iso)
                     CALL project_function_b(ccdens(:, 2, ikind), ccdens(:, 1, ikind), bfun(:, iso), grid_atom, l)
                     IF (unit_nr > 0) THEN
                        WRITE (unit_nr, FMT="(3I6)") iso, l, ngto
                        WRITE (unit_nr, FMT="(2G24.12)") (ccdens(i, 1, ikind), ccdens(i, 2, ikind), i=1, ngto)
                     END IF
                  END DO
                  DEALLOCATE (bfun)
               END DO
            ELSE
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, *) "Coordinates"
                  np = particles%n_els
                  DO iat = 1, np
                     CALL get_atomic_kind(particles%els(iat)%atomic_kind, kind_number=ikind)
                     WRITE (unit_nr, '(I10,I5,3f12.6)') iat, ikind, particles%els(iat)%r
                  END DO
               END IF
            END IF

            DEALLOCATE (ppdens, aedens, ccdens)

            CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                              "DFT%PRINT%E_DENSITY_CUBE")

         END IF
         IF (dft_control%qs_control%gapw .AND. print_total_density) THEN
            ! total density in g-space not implemented for k-points
            CPASSERT(.NOT. do_kpoints)
            ! Print total electronic density
            CALL get_qs_env(qs_env=qs_env, &
                            pw_env=pw_env)
            CALL pw_env_get(pw_env=pw_env, &
                            auxbas_pw_pool=auxbas_pw_pool, &
                            pw_pools=pw_pools)
            CALL pw_pool_create_pw(pool=auxbas_pw_pool, &
                                   pw=rho_elec_rspace%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)
            CALL pw_zero(rho_elec_rspace%pw)
            CALL pw_pool_create_pw(pool=auxbas_pw_pool, &
                                   pw=rho_elec_gspace%pw, &
                                   use_data=COMPLEXDATA1D, &
                                   in_space=RECIPROCALSPACE)
            CALL pw_zero(rho_elec_gspace%pw)
            CALL get_pw_grid_info(pw_grid=rho_elec_gspace%pw%pw_grid, &
                                  dr=dr, &
                                  vol=volume)
            q_max = SQRT(SUM((pi/dr(:))**2))
            CALL calculate_rhotot_elec_gspace(qs_env=qs_env, &
                                              auxbas_pw_pool=auxbas_pw_pool, &
                                              rhotot_elec_gspace=rho_elec_gspace, &
                                              q_max=q_max, &
                                              rho_hard=rho_hard, &
                                              rho_soft=rho_soft)
            rho_total = rho_hard + rho_soft
            CALL get_pw_grid_info(pw_grid=rho_elec_gspace%pw%pw_grid, &
                                  vol=volume)
            CALL pw_transfer(rho_elec_gspace%pw, rho_elec_rspace%pw, debug=.FALSE.)
            rho_total_rspace = pw_integrate_function(rho_elec_rspace%pw, isign=-1)/volume
            filename = "TOTAL_ELECTRON_DENSITY"
            mpi_io = .TRUE.
            unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%E_DENSITY_CUBE", &
                                           extension=".cube", middle_name=TRIM(filename), &
                                           file_position=my_pos_cube, log_filename=.FALSE., mpi_io=mpi_io, &
                                           fout=mpi_filename)
            IF (output_unit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=unit_nr, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF
               WRITE (UNIT=output_unit, FMT="(/,T2,A,/,/,T2,A)") &
                  "The total electron density is written in cube file format to the file:", &
                  TRIM(filename)
               WRITE (UNIT=output_unit, FMT="(/,(T2,A,F20.10))") &
                  "q(max) [1/Angstrom]              :", q_max/angstrom, &
                  "Soft electronic charge (G-space) :", rho_soft, &
                  "Hard electronic charge (G-space) :", rho_hard, &
                  "Total electronic charge (G-space):", rho_total, &
                  "Total electronic charge (R-space):", rho_total_rspace
            END IF
            CALL cp_pw_to_cube(rho_elec_rspace%pw, unit_nr, "TOTAL ELECTRON DENSITY", &
                               particles=particles, &
                               stride=section_get_ivals(dft_section, "PRINT%E_DENSITY_CUBE%STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                              "DFT%PRINT%E_DENSITY_CUBE", mpi_io=mpi_io)
            ! Print total spin density for spin-polarized systems
            IF (dft_control%nspins > 1) THEN
               CALL pw_zero(rho_elec_gspace%pw)
               CALL pw_zero(rho_elec_rspace%pw)
               CALL calculate_rhotot_elec_gspace(qs_env=qs_env, &
                                                 auxbas_pw_pool=auxbas_pw_pool, &
                                                 rhotot_elec_gspace=rho_elec_gspace, &
                                                 q_max=q_max, &
                                                 rho_hard=rho_hard, &
                                                 rho_soft=rho_soft, &
                                                 fsign=-1.0_dp)
               rho_total = rho_hard + rho_soft
               CALL pw_transfer(rho_elec_gspace%pw, rho_elec_rspace%pw, debug=.FALSE.)
               rho_total_rspace = pw_integrate_function(rho_elec_rspace%pw, isign=-1)/volume
               filename = "TOTAL_SPIN_DENSITY"
               mpi_io = .TRUE.
               unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%E_DENSITY_CUBE", &
                                              extension=".cube", middle_name=TRIM(filename), &
                                              file_position=my_pos_cube, log_filename=.FALSE., mpi_io=mpi_io, &
                                              fout=mpi_filename)
               IF (output_unit > 0) THEN
                  IF (.NOT. mpi_io) THEN
                     INQUIRE (UNIT=unit_nr, NAME=filename)
                  ELSE
                     filename = mpi_filename
                  END IF
                  WRITE (UNIT=output_unit, FMT="(/,T2,A,/,/,T2,A)") &
                     "The total spin density is written in cube file format to the file:", &
                     TRIM(filename)
                  WRITE (UNIT=output_unit, FMT="(/,(T2,A,F20.10))") &
                     "q(max) [1/Angstrom]                    :", q_max/angstrom, &
                     "Soft part of the spin density (G-space):", rho_soft, &
                     "Hard part of the spin density (G-space):", rho_hard, &
                     "Total spin density (G-space)           :", rho_total, &
                     "Total spin density (R-space)           :", rho_total_rspace
               END IF
               CALL cp_pw_to_cube(rho_elec_rspace%pw, unit_nr, "TOTAL SPIN DENSITY", &
                                  particles=particles, &
                                  stride=section_get_ivals(dft_section, "PRINT%E_DENSITY_CUBE%STRIDE"), mpi_io=mpi_io)
               CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                                 "DFT%PRINT%E_DENSITY_CUBE", mpi_io=mpi_io)
            END IF
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_elec_gspace%pw)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_elec_rspace%pw)
         ELSE
            IF (dft_control%nspins > 1) THEN
               CALL get_qs_env(qs_env=qs_env, &
                               pw_env=pw_env)
               CALL pw_env_get(pw_env=pw_env, &
                               auxbas_pw_pool=auxbas_pw_pool, &
                               pw_pools=pw_pools)
               CALL pw_pool_create_pw(pool=auxbas_pw_pool, &
                                      pw=rho_elec_rspace%pw, &
                                      use_data=REALDATA3D, &
                                      in_space=REALSPACE)
               CALL pw_copy(rho_r(1)%pw, rho_elec_rspace%pw)
               CALL pw_axpy(rho_r(2)%pw, rho_elec_rspace%pw)
               filename = "ELECTRON_DENSITY"
               mpi_io = .TRUE.
               unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%E_DENSITY_CUBE", &
                                              extension=".cube", middle_name=TRIM(filename), &
                                              file_position=my_pos_cube, log_filename=.FALSE., mpi_io=mpi_io, &
                                              fout=mpi_filename)
               IF (output_unit > 0) THEN
                  IF (.NOT. mpi_io) THEN
                     INQUIRE (UNIT=unit_nr, NAME=filename)
                  ELSE
                     filename = mpi_filename
                  END IF
                  WRITE (UNIT=output_unit, FMT="(/,T2,A,/,/,T2,A)") &
                     "The sum of alpha and beta density is written in cube file format to the file:", &
                     TRIM(filename)
               END IF
               CALL cp_pw_to_cube(rho_elec_rspace%pw, unit_nr, "SUM OF ALPHA AND BETA DENSITY", &
                                  particles=particles, stride=section_get_ivals(dft_section, "PRINT%E_DENSITY_CUBE%STRIDE"), &
                                  mpi_io=mpi_io)
               CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                                 "DFT%PRINT%E_DENSITY_CUBE", mpi_io=mpi_io)
               CALL pw_copy(rho_r(1)%pw, rho_elec_rspace%pw)
               CALL pw_axpy(rho_r(2)%pw, rho_elec_rspace%pw, alpha=-1.0_dp)
               filename = "SPIN_DENSITY"
               mpi_io = .TRUE.
               unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%E_DENSITY_CUBE", &
                                              extension=".cube", middle_name=TRIM(filename), &
                                              file_position=my_pos_cube, log_filename=.FALSE., mpi_io=mpi_io, &
                                              fout=mpi_filename)
               IF (output_unit > 0) THEN
                  IF (.NOT. mpi_io) THEN
                     INQUIRE (UNIT=unit_nr, NAME=filename)
                  ELSE
                     filename = mpi_filename
                  END IF
                  WRITE (UNIT=output_unit, FMT="(/,T2,A,/,/,T2,A)") &
                     "The spin density is written in cube file format to the file:", &
                     TRIM(filename)
               END IF
               CALL cp_pw_to_cube(rho_elec_rspace%pw, unit_nr, "SPIN DENSITY", &
                                  particles=particles, &
                                  stride=section_get_ivals(dft_section, "PRINT%E_DENSITY_CUBE%STRIDE"), mpi_io=mpi_io)
               CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                                 "DFT%PRINT%E_DENSITY_CUBE", mpi_io=mpi_io)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_elec_rspace%pw)
            ELSE
               filename = "ELECTRON_DENSITY"
               mpi_io = .TRUE.
               unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%E_DENSITY_CUBE", &
                                              extension=".cube", middle_name=TRIM(filename), &
                                              file_position=my_pos_cube, log_filename=.FALSE., mpi_io=mpi_io, &
                                              fout=mpi_filename)
               IF (output_unit > 0) THEN
                  IF (.NOT. mpi_io) THEN
                     INQUIRE (UNIT=unit_nr, NAME=filename)
                  ELSE
                     filename = mpi_filename
                  END IF
                  WRITE (UNIT=output_unit, FMT="(/,T2,A,/,/,T2,A)") &
                     "The electron density is written in cube file format to the file:", &
                     TRIM(filename)
               END IF
               CALL cp_pw_to_cube(rho_r(1)%pw, unit_nr, "ELECTRON DENSITY", &
                                  particles=particles, &
                                  stride=section_get_ivals(dft_section, "PRINT%E_DENSITY_CUBE%STRIDE"), mpi_io=mpi_io)
               CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                                 "DFT%PRINT%E_DENSITY_CUBE", mpi_io=mpi_io)
            END IF ! nspins
         END IF ! total density for GAPW
      END IF ! print key

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           dft_section, "PRINT%ENERGY_WINDOWS"), cp_p_file) .AND. .NOT. do_kpoints) THEN
         CALL energy_windows(qs_env)
      ENDIF

      ! Print the hartree potential
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%PRINT%V_HARTREE_CUBE"), cp_p_file)) THEN

         CALL get_qs_env(qs_env=qs_env, &
                         pw_env=pw_env, &
                         v_hartree_rspace=v_hartree_rspace)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL pw_pool_create_pw(auxbas_pw_pool, aux_r%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)

         append_cube = section_get_lval(input, "DFT%PRINT%V_HARTREE_CUBE%APPEND")
         my_pos_cube = "REWIND"
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         END IF
         mpi_io = .TRUE.
         CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
         CALL pw_env_get(pw_env)
         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%V_HARTREE_CUBE", &
                                        extension=".cube", middle_name="v_hartree", file_position=my_pos_cube, mpi_io=mpi_io)
         udvol = 1.0_dp/v_hartree_rspace%pw_grid%dvol

         CALL pw_copy(v_hartree_rspace, aux_r%pw)
         CALL pw_scale(aux_r%pw, udvol)

         CALL cp_pw_to_cube(aux_r%pw, unit_nr, "HARTREE POTENTIAL", particles=particles, &
                            stride=section_get_ivals(dft_section, &
                                                     "PRINT%V_HARTREE_CUBE%STRIDE"), mpi_io=mpi_io)
         CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                           "DFT%PRINT%V_HARTREE_CUBE", mpi_io=mpi_io)

         CALL pw_pool_give_back_pw(auxbas_pw_pool, aux_r%pw)
      ENDIF

      ! Print the external potential
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%PRINT%EXTERNAL_POTENTIAL_CUBE"), cp_p_file)) THEN
         IF (dft_control%apply_external_potential) THEN
            CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, vee=vee)
            CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
            CALL pw_pool_create_pw(auxbas_pw_pool, aux_r%pw, &
                                   use_data=REALDATA3D, &
                                   in_space=REALSPACE)

            append_cube = section_get_lval(input, "DFT%PRINT%EXTERNAL_POTENTIAL_CUBE%APPEND")
            my_pos_cube = "REWIND"
            IF (append_cube) THEN
               my_pos_cube = "APPEND"
            END IF
            mpi_io = .TRUE.
            CALL pw_env_get(pw_env)
            unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%EXTERNAL_POTENTIAL_CUBE", &
                                           extension=".cube", middle_name="ext_pot", file_position=my_pos_cube, mpi_io=mpi_io)

            CALL pw_copy(vee%pw, aux_r%pw)

            CALL cp_pw_to_cube(aux_r%pw, unit_nr, "EXTERNAL POTENTIAL", particles=particles, &
                               stride=section_get_ivals(dft_section, &
                                                        "PRINT%EXTERNAL_POTENTIAL_CUBE%STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                              "DFT%PRINT%EXTERNAL_POTENTIAL_CUBE", mpi_io=mpi_io)

            CALL pw_pool_give_back_pw(auxbas_pw_pool, aux_r%pw)
         ENDIF
      ENDIF

      ! Print the Electrical Field Components
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%PRINT%EFIELD_CUBE"), cp_p_file)) THEN

         CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL pw_pool_create_pw(auxbas_pw_pool, aux_r%pw, &
                                use_data=REALDATA3D, &
                                in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, aux_g%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)

         append_cube = section_get_lval(input, "DFT%PRINT%EFIELD_CUBE%APPEND")
         my_pos_cube = "REWIND"
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         END IF
         CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, &
                         v_hartree_rspace=v_hartree_rspace)
         CALL pw_env_get(pw_env)
         udvol = 1.0_dp/v_hartree_rspace%pw_grid%dvol
         DO id = 1, 3
            mpi_io = .TRUE.
            unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%EFIELD_CUBE", &
                                           extension=".cube", middle_name="efield_"//cdir(id), file_position=my_pos_cube, &
                                           mpi_io=mpi_io)

            CALL pw_transfer(v_hartree_rspace, aux_g%pw)
            nd = 0
            nd(id) = 1
            CALL pw_derive(aux_g%pw, nd)
            CALL pw_transfer(aux_g%pw, aux_r%pw)
            CALL pw_scale(aux_r%pw, udvol)

            CALL cp_pw_to_cube(aux_r%pw, &
                               unit_nr, "ELECTRIC FIELD", particles=particles, &
                               stride=section_get_ivals(dft_section, &
                                                        "PRINT%EFIELD_CUBE%STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                              "DFT%PRINT%EFIELD_CUBE", mpi_io=mpi_io)
         END DO

         CALL pw_pool_give_back_pw(auxbas_pw_pool, aux_r%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, aux_g%pw)
      END IF

      ! Write cube files from the local energy
      CALL qs_scf_post_local_energy(input, logger, qs_env)

      ! Write cube files from the local stress tensor
      CALL qs_scf_post_local_stress(input, logger, qs_env)

      ! Write cube files from the implicit Poisson solver
      CALL qs_scf_post_ps_implicit(input, logger, qs_env)

      ! post SCF Transport
      CALL qs_scf_post_transport(qs_env)

      CALL section_vals_val_get(input, "DFT%PRINT%AO_MATRICES%OMIT_HEADERS", l_val=omit_headers)
      ! Write the density matrices
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%PRINT%AO_MATRICES/DENSITY"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, input, "DFT%PRINT%AO_MATRICES/DENSITY", &
                                   extension=".Log")
         CALL section_vals_val_get(input, "DFT%PRINT%AO_MATRICES%NDIGITS", i_val=after)
         CALL qs_rho_get(rho, rho_ao_kp=rho_ao)
         after = MIN(MAX(after, 1), 16)
         DO ispin = 1, dft_control%nspins
            DO img = 1, dft_control%nimages
               CALL cp_dbcsr_write_sparse_matrix(rho_ao(ispin, img)%matrix, 4, after, qs_env, &
                                                 para_env, output_unit=iw, omit_headers=omit_headers)
            END DO
         END DO
         CALL cp_print_key_finished_output(iw, logger, input, &
                                           "DFT%PRINT%AO_MATRICES/DENSITY")
      END IF

      ! Write the Kohn-Sham matrices
      write_ks = BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                  "DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX"), cp_p_file)
      write_xc = BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                  "DFT%PRINT%AO_MATRICES/MATRIX_VXC"), cp_p_file)
      ! we need to update stuff before writing, potentially computing the matrix_vxc
      IF (write_ks .OR. write_xc) THEN
         IF (write_xc) qs_env%requires_matrix_vxc = .TRUE.
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
         CALL qs_ks_update_qs_env(qs_env, calculate_forces=.FALSE., &
                                  just_energy=.FALSE.)
         IF (write_xc) qs_env%requires_matrix_vxc = .FALSE.
      END IF

      ! Write the Kohn-Sham matrices
      IF (write_ks) THEN
         iw = cp_print_key_unit_nr(logger, input, "DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX", &
                                   extension=".Log")
         CALL get_qs_env(qs_env=qs_env, matrix_ks_kp=ks_rmpv)
         CALL section_vals_val_get(input, "DFT%PRINT%AO_MATRICES%NDIGITS", i_val=after)
         after = MIN(MAX(after, 1), 16)
         DO ispin = 1, dft_control%nspins
            DO img = 1, dft_control%nimages
               CALL cp_dbcsr_write_sparse_matrix(ks_rmpv(ispin, img)%matrix, 4, after, qs_env, &
                                                 para_env, output_unit=iw, omit_headers=omit_headers)
            END DO
         END DO
         CALL cp_print_key_finished_output(iw, logger, input, &
                                           "DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX")
      END IF

      ! write csr matrices
      ! matrices in terms of the PAO basis will be taken care of in pao_post_scf.
      IF (.NOT. dft_control%qs_control%pao) THEN
         CALL write_ks_matrix_csr(qs_env, input)
         CALL write_s_matrix_csr(qs_env, input)
      END IF

      ! write adjacency matrix
      CALL write_adjacency_matrix(qs_env, input)

      ! Write the xc matrix
      IF (write_xc) THEN
         CALL get_qs_env(qs_env=qs_env, matrix_vxc_kp=matrix_vxc)
         CPASSERT(ASSOCIATED(matrix_vxc))
         iw = cp_print_key_unit_nr(logger, input, "DFT%PRINT%AO_MATRICES/MATRIX_VXC", &
                                   extension=".Log")
         CALL section_vals_val_get(input, "DFT%PRINT%AO_MATRICES%NDIGITS", i_val=after)
         after = MIN(MAX(after, 1), 16)
         DO ispin = 1, dft_control%nspins
            DO img = 1, dft_control%nimages
               CALL cp_dbcsr_write_sparse_matrix(matrix_vxc(ispin, img)%matrix, 4, after, qs_env, &
                                                 para_env, output_unit=iw, omit_headers=omit_headers)
            END DO
         END DO
         CALL cp_print_key_finished_output(iw, logger, input, &
                                           "DFT%PRINT%AO_MATRICES/MATRIX_VXC")
      END IF

      ! Write the [H,r] commutator matrices
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%PRINT%AO_MATRICES/COMMUTATOR_HR"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, input, "DFT%PRINT%AO_MATRICES/COMMUTATOR_HR", &
                                   extension=".Log")
         CALL section_vals_val_get(input, "DFT%PRINT%AO_MATRICES%NDIGITS", i_val=after)
         NULLIFY (matrix_hr)
         CALL build_com_hr_matrix(qs_env, matrix_hr)
         after = MIN(MAX(after, 1), 16)
         DO img = 1, 3
            CALL cp_dbcsr_write_sparse_matrix(matrix_hr(img)%matrix, 4, after, qs_env, &
                                              para_env, output_unit=iw, omit_headers=omit_headers)
         END DO
         CALL dbcsr_deallocate_matrix_set(matrix_hr)
         CALL cp_print_key_finished_output(iw, logger, input, &
                                           "DFT%PRINT%AO_MATRICES/COMMUTATOR_HR")
      END IF

      ! Compute the Mulliken charges
      print_key => section_vals_get_subs_vals(input, "DFT%PRINT%MULLIKEN")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%MULLIKEN", extension=".mulliken", log_filename=.FALSE.)
         print_level = 1
         CALL section_vals_val_get(print_key, "PRINT_GOP", l_val=print_it)
         IF (print_it) print_level = 2
         CALL section_vals_val_get(print_key, "PRINT_ALL", l_val=print_it)
         IF (print_it) print_level = 3
         CALL mulliken_population_analysis(qs_env, unit_nr, print_level)
         CALL cp_print_key_finished_output(unit_nr, logger, input, "DFT%PRINT%MULLIKEN")
      END IF

      ! Compute the Hirshfeld charges
      print_key => section_vals_get_subs_vals(input, "DFT%PRINT%HIRSHFELD")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         ! we check if real space density is available
         NULLIFY (rho)
         CALL get_qs_env(qs_env=qs_env, rho=rho)
         CALL qs_rho_get(rho, rho_r_valid=rho_r_valid)
         IF (rho_r_valid) THEN
            unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%HIRSHFELD", extension=".hirshfeld", log_filename=.FALSE.)
            CALL hirshfeld_charges(qs_env, print_key, unit_nr)
            CALL cp_print_key_finished_output(unit_nr, logger, input, "DFT%PRINT%HIRSHFELD")
         END IF
      END IF

      ! MAO analysis
      print_key => section_vals_get_subs_vals(input, "DFT%PRINT%MAO_ANALYSIS")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%MAO_ANALYSIS", extension=".mao", log_filename=.FALSE.)
         CALL mao_analysis(qs_env, print_key, unit_nr)
         CALL cp_print_key_finished_output(unit_nr, logger, input, "DFT%PRINT%MAO_ANALYSIS")
      END IF

      ! MINBAS analysis
      print_key => section_vals_get_subs_vals(input, "DFT%PRINT%MINBAS_ANALYSIS")
      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_key), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%MINBAS_ANALYSIS", extension=".mao", log_filename=.FALSE.)
         CALL minbas_analysis(qs_env, print_key, unit_nr)
         CALL cp_print_key_finished_output(unit_nr, logger, input, "DFT%PRINT%MINBAS_ANALYSIS")
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_mo_free_results

! **************************************************************************************************
!> \brief Calculates Hirshfeld charges
!> \param qs_env the qs_env where to calculate the charges
!> \param input_section the input section for Hirshfeld charges
!> \param unit_nr the output unit number
! **************************************************************************************************
   SUBROUTINE hirshfeld_charges(qs_env, input_section, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input_section
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=*), PARAMETER :: routineN = 'hirshfeld_charges', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iat, ikind, natom, nkind, nspin, &
                                                            radius_type, refc, shapef
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: do_radius, do_sc, paw_atom
      REAL(KIND=dp)                                      :: zeff
      REAL(KIND=dp), DIMENSION(:), POINTER               :: radii
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: charges
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hirshfeld_type), POINTER                      :: hirshfeld_env
      TYPE(mpole_rho_atom), DIMENSION(:), POINTER        :: mp_rho
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(rho0_mpole_type), POINTER                     :: rho0_mpole

      NULLIFY (hirshfeld_env)
      NULLIFY (radii)
      CALL create_hirshfeld_type(hirshfeld_env)
      !
      CALL get_qs_env(qs_env, nkind=nkind, natom=natom)
      ALLOCATE (hirshfeld_env%charges(natom))
      ! input options
      CALL section_vals_val_get(input_section, "SELF_CONSISTENT", l_val=do_sc)
      CALL section_vals_val_get(input_section, "USER_RADIUS", l_val=do_radius)
      CALL section_vals_val_get(input_section, "SHAPE_FUNCTION", i_val=shapef)
      CALL section_vals_val_get(input_section, "REFERENCE_CHARGE", i_val=refc)
      IF (do_radius) THEN
         radius_type = radius_user
         CALL section_vals_val_get(input_section, "ATOMIC_RADII", r_vals=radii)
         IF (.NOT. SIZE(radii) == nkind) &
            CALL cp_abort(__LOCATION__, &
                          "Length of keyword HIRSHFELD\ATOMIC_RADII does not "// &
                          "match number of atomic kinds in the input coordinate file.")
      ELSE
         radius_type = radius_covalent
      END IF
      CALL set_hirshfeld_info(hirshfeld_env, shape_function_type=shapef, &
                              iterative=do_sc, ref_charge=refc, &
                              radius_type=radius_type)
      ! shape function
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set)
      CALL create_shape_function(hirshfeld_env, qs_kind_set, atomic_kind_set, &
                                 radii_list=radii)
      ! reference charges
      CALL get_qs_env(qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      nspin = SIZE(matrix_p, 1)
      ALLOCATE (charges(natom, nspin))
      SELECT CASE (refc)
      CASE (ref_charge_atomic)
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
            atomic_kind => atomic_kind_set(ikind)
            CALL get_atomic_kind(atomic_kind, atom_list=atom_list)
            DO iat = 1, SIZE(atom_list)
               i = atom_list(iat)
               hirshfeld_env%charges(i) = zeff
            END DO
         END DO
      CASE (ref_charge_mulliken)
         CALL get_qs_env(qs_env, matrix_s_kp=matrix_s, para_env=para_env)
         CALL mulliken_charges(matrix_p, matrix_s, para_env, charges)
         DO iat = 1, natom
            hirshfeld_env%charges(iat) = SUM(charges(iat, :))
         END DO
      CASE DEFAULT
         CPABORT("Unknown type of reference charge for Hirshfeld partitioning.")
      END SELECT
      !
      charges = 0.0_dp
      IF (hirshfeld_env%iterative) THEN
         ! Hirshfeld-I charges
         CALL comp_hirshfeld_i_charges(qs_env, hirshfeld_env, charges, unit_nr)
      ELSE
         ! Hirshfeld charges
         CALL comp_hirshfeld_charges(qs_env, hirshfeld_env, charges)
      END IF
      CALL get_qs_env(qs_env, particle_set=particle_set, dft_control=dft_control)
      IF (dft_control%qs_control%gapw) THEN
         ! GAPW: add core charges (rho_hard - rho_soft)
         CALL get_qs_env(qs_env, rho0_mpole=rho0_mpole)
         CALL get_rho0_mpole(rho0_mpole, mp_rho=mp_rho)
         DO iat = 1, natom
            atomic_kind => particle_set(iat)%atomic_kind
            CALL get_atomic_kind(atomic_kind, kind_number=ikind)
            CALL get_qs_kind(qs_kind_set(ikind), paw_atom=paw_atom)
            IF (paw_atom) THEN
               charges(iat, 1:nspin) = charges(iat, 1:nspin) + mp_rho(iat)%q0(1:nspin)
            END IF
         END DO
      END IF
      !
      IF (unit_nr > 0) THEN
         CALL write_hirshfeld_charges(charges, hirshfeld_env, particle_set, &
                                      qs_kind_set, unit_nr)
      END IF
      !
      CALL release_hirshfeld_type(hirshfeld_env)
      DEALLOCATE (charges)

   END SUBROUTINE hirshfeld_charges

! **************************************************************************************************
!> \brief ...
!> \param ca ...
!> \param a ...
!> \param cb ...
!> \param b ...
!> \param l ...
! **************************************************************************************************
   SUBROUTINE project_function_a(ca, a, cb, b, l)
      ! project function cb on ca
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: ca
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: a, cb, b
      INTEGER, INTENT(IN)                                :: l

      CHARACTER(len=*), PARAMETER :: routineN = 'project_function_a', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: info, n
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ipiv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: smat, tmat, v

      n = SIZE(ca)
      ALLOCATE (smat(n, n), tmat(n, n), v(n, 1), ipiv(n))

      CALL sg_overlap(smat, l, a, a)
      CALL sg_overlap(tmat, l, a, b)
      v(:, 1) = MATMUL(tmat, cb)
      CALL lapack_sgesv(n, 1, smat, n, ipiv, v, n, info)
      CPASSERT(info == 0)
      ca(:) = v(:, 1)

      DEALLOCATE (smat, tmat, v, ipiv)

   END SUBROUTINE project_function_a

! **************************************************************************************************
!> \brief ...
!> \param ca ...
!> \param a ...
!> \param bfun ...
!> \param grid_atom ...
!> \param l ...
! **************************************************************************************************
   SUBROUTINE project_function_b(ca, a, bfun, grid_atom, l)
      ! project function f on ca
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: ca
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: a, bfun
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      INTEGER, INTENT(IN)                                :: l

      CHARACTER(len=*), PARAMETER :: routineN = 'project_function_b', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, info, n, nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ipiv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: afun
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: smat, v

      n = SIZE(ca)
      nr = grid_atom%nr
      ALLOCATE (smat(n, n), v(n, 1), ipiv(n), afun(nr))

      CALL sg_overlap(smat, l, a, a)
      DO i = 1, n
         afun(:) = grid_atom%rad(:)**l*EXP(-a(i)*grid_atom%rad2(:))
         v(i, 1) = SUM(afun(:)*bfun(:)*grid_atom%wr(:))
      END DO
      CALL lapack_sgesv(n, 1, smat, n, ipiv, v, n, info)
      CPASSERT(info == 0)
      ca(:) = v(:, 1)

      DEALLOCATE (smat, v, ipiv, afun)

   END SUBROUTINE project_function_b

! **************************************************************************************************
!> \brief Performs printing of cube files from local energy
!> \param input input
!> \param logger the logger
!> \param qs_env the qs_env in which the qs_env lives
!> \par History
!>      07.2019 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE qs_scf_post_local_energy(input, logger, qs_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_local_energy', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename, my_pos_cube
      INTEGER                                            :: handle, io_unit, natom, unit_nr
      LOGICAL                                            :: append_cube, gapw, gapw_xc, mpi_io
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: eden
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: dft_section

      CALL timeset(routineN, handle)
      io_unit = cp_logger_get_default_io_unit(logger)
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%PRINT%LOCAL_ENERGY_CUBE"), cp_p_file)) THEN
         dft_section => section_vals_get_subs_vals(input, "DFT")
         CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, natom=natom)
         gapw = dft_control%qs_control%gapw
         gapw_xc = dft_control%qs_control%gapw_xc
         CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, subsys=subsys)
         CALL qs_subsys_get(subsys, particles=particles)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL pw_pool_create_pw(auxbas_pw_pool, eden%pw, use_data=REALDATA3D, in_space=REALSPACE)
         !
         CALL qs_local_energy(qs_env, eden)
         !
         append_cube = section_get_lval(input, "DFT%PRINT%LOCAL_ENERGY_CUBE%APPEND")
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         ELSE
            my_pos_cube = "REWIND"
         END IF
         mpi_io = .TRUE.
         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%LOCAL_ENERGY_CUBE", &
                                        extension=".cube", middle_name="local_energy", &
                                        file_position=my_pos_cube, mpi_io=mpi_io)
         CALL cp_pw_to_cube(eden%pw, &
                            unit_nr, "LOCAL ENERGY", particles=particles, &
                            stride=section_get_ivals(dft_section, &
                                                     "PRINT%LOCAL_ENERGY_CUBE%STRIDE"), mpi_io=mpi_io)
         IF (io_unit > 0) THEN
            INQUIRE (UNIT=unit_nr, NAME=filename)
            IF (gapw .OR. gapw_xc) THEN
               WRITE (UNIT=io_unit, FMT="(/,T3,A,A)") &
                  "The soft part of the local energy is written to the file: ", TRIM(ADJUSTL(filename))
            ELSE
               WRITE (UNIT=io_unit, FMT="(/,T3,A,A)") &
                  "The local energy is written to the file: ", TRIM(ADJUSTL(filename))
            END IF
         END IF
         CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                           "DFT%PRINT%LOCAL_ENERGY_CUBE", mpi_io=mpi_io)
         !
         CALL pw_pool_give_back_pw(auxbas_pw_pool, eden%pw)
      END IF
      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_local_energy

! **************************************************************************************************
!> \brief Performs printing of cube files from local energy
!> \param input input
!> \param logger the logger
!> \param qs_env the qs_env in which the qs_env lives
!> \par History
!>      07.2019 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE qs_scf_post_local_stress(input, logger, qs_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_local_stress', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename, my_pos_cube
      INTEGER                                            :: handle, io_unit, natom, unit_nr
      LOGICAL                                            :: append_cube, gapw, gapw_xc, mpi_io
      REAL(KIND=dp)                                      :: beta
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: stress
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: dft_section

      CALL timeset(routineN, handle)
      io_unit = cp_logger_get_default_io_unit(logger)
      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                           "DFT%PRINT%LOCAL_STRESS_CUBE"), cp_p_file)) THEN
         dft_section => section_vals_get_subs_vals(input, "DFT")
         CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, natom=natom)
         gapw = dft_control%qs_control%gapw
         gapw_xc = dft_control%qs_control%gapw_xc
         CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, subsys=subsys)
         CALL qs_subsys_get(subsys, particles=particles)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL pw_pool_create_pw(auxbas_pw_pool, stress%pw, use_data=REALDATA3D, in_space=REALSPACE)
         !
         ! use beta=0: kinetic energy density in symmetric form
         beta = 0.0_dp
         CALL qs_local_stress(qs_env, pressure=stress, beta=beta)
         !
         append_cube = section_get_lval(input, "DFT%PRINT%LOCAL_STRESS_CUBE%APPEND")
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         ELSE
            my_pos_cube = "REWIND"
         END IF
         mpi_io = .TRUE.
         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%LOCAL_STRESS_CUBE", &
                                        extension=".cube", middle_name="local_stress", &
                                        file_position=my_pos_cube, mpi_io=mpi_io)
         CALL cp_pw_to_cube(stress%pw, &
                            unit_nr, "LOCAL STRESS", particles=particles, &
                            stride=section_get_ivals(dft_section, &
                                                     "PRINT%LOCAL_STRESS_CUBE%STRIDE"), mpi_io=mpi_io)
         IF (io_unit > 0) THEN
            INQUIRE (UNIT=unit_nr, NAME=filename)
            WRITE (UNIT=io_unit, FMT="(/,T3,A)") "Write 1/3*Tr(sigma) to cube file"
            IF (gapw .OR. gapw_xc) THEN
               WRITE (UNIT=io_unit, FMT="(T3,A,A)") &
                  "The soft part of the local stress is written to the file: ", TRIM(ADJUSTL(filename))
            ELSE
               WRITE (UNIT=io_unit, FMT="(T3,A,A)") &
                  "The local stress is written to the file: ", TRIM(ADJUSTL(filename))
            END IF
         END IF
         CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                           "DFT%PRINT%LOCAL_STRESS_CUBE", mpi_io=mpi_io)
         !
         CALL pw_pool_give_back_pw(auxbas_pw_pool, stress%pw)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_local_stress

! **************************************************************************************************
!> \brief Performs printing of cube files related to the implicit Poisson solver
!> \param input input
!> \param logger the logger
!> \param qs_env the qs_env in which the qs_env lives
!> \par History
!>      03.2016 refactored from write_mo_free_results [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE qs_scf_post_ps_implicit(input, logger, qs_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_ps_implicit', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename, my_pos_cube
      INTEGER                                            :: boundary_condition, handle, i, j, &
                                                            n_cstr, n_tiles, unit_nr
      LOGICAL :: append_cube, do_cstr_charge_cube, do_dielectric_cube, do_dirichlet_bc_cube, &
         has_dirichlet_bc, has_implicit_ps, mpi_io, tile_cubes
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: aux_r
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type), POINTER                             :: dirichlet_tile
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: dft_section

      CALL timeset(routineN, handle)

      NULLIFY (pw_env, auxbas_pw_pool, dft_section, particles)

      dft_section => section_vals_get_subs_vals(input, "DFT")
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, subsys=subsys)
      CALL qs_subsys_get(subsys, particles=particles)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      has_implicit_ps = .FALSE.
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      IF (pw_env%poisson_env%parameters%solver .EQ. pw_poisson_implicit) has_implicit_ps = .TRUE.

      ! Write the dielectric constant into a cube file
      do_dielectric_cube = BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                            "DFT%PRINT%IMPLICIT_PSOLVER%DIELECTRIC_CUBE"), cp_p_file)
      IF (has_implicit_ps .AND. do_dielectric_cube) THEN
         append_cube = section_get_lval(input, "DFT%PRINT%IMPLICIT_PSOLVER%DIELECTRIC_CUBE%APPEND")
         my_pos_cube = "REWIND"
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         END IF
         mpi_io = .TRUE.
         unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%IMPLICIT_PSOLVER%DIELECTRIC_CUBE", &
                                        extension=".cube", middle_name="DIELECTRIC_CONSTANT", file_position=my_pos_cube, &
                                        mpi_io=mpi_io)
         CALL pw_env_get(pw_env, poisson_env=poisson_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL pw_pool_create_pw(auxbas_pw_pool, aux_r%pw, use_data=REALDATA3D, in_space=REALSPACE)

         boundary_condition = pw_env%poisson_env%parameters%ps_implicit_params%boundary_condition
         SELECT CASE (boundary_condition)
         CASE (PERIODIC_BC, MIXED_PERIODIC_BC)
            CALL pw_copy(poisson_env%implicit_env%dielectric%eps, aux_r%pw)
         CASE (MIXED_BC, NEUMANN_BC)
            CALL pw_shrink(pw_env%poisson_env%parameters%ps_implicit_params%neumann_directions, &
                           pw_env%poisson_env%implicit_env%dct_env%dests_shrink, &
                           pw_env%poisson_env%implicit_env%dct_env%srcs_shrink, &
                           pw_env%poisson_env%implicit_env%dct_env%bounds_local_shftd, &
                           poisson_env%implicit_env%dielectric%eps, aux_r%pw)
         END SELECT

         CALL cp_pw_to_cube(aux_r%pw, unit_nr, "DIELECTRIC CONSTANT", particles=particles, &
                            stride=section_get_ivals(dft_section, "PRINT%IMPLICIT_PSOLVER%DIELECTRIC_CUBE%STRIDE"), &
                            mpi_io=mpi_io)
         CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                           "DFT%PRINT%IMPLICIT_PSOLVER%DIELECTRIC_CUBE", mpi_io=mpi_io)

         CALL pw_pool_give_back_pw(auxbas_pw_pool, aux_r%pw)
      ENDIF

      ! Write Dirichlet constraint charges into a cube file
      do_cstr_charge_cube = BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                             "DFT%PRINT%IMPLICIT_PSOLVER%DIRICHLET_CSTR_CHARGE_CUBE"), cp_p_file)

      has_dirichlet_bc = .FALSE.
      IF (has_implicit_ps) THEN
         boundary_condition = pw_env%poisson_env%parameters%ps_implicit_params%boundary_condition
         IF (boundary_condition .EQ. MIXED_PERIODIC_BC .OR. boundary_condition .EQ. MIXED_BC) THEN
            has_dirichlet_bc = .TRUE.
         END IF
      END IF

      IF (has_implicit_ps .AND. do_cstr_charge_cube .AND. has_dirichlet_bc) THEN
         append_cube = section_get_lval(input, &
                                        "DFT%PRINT%IMPLICIT_PSOLVER%DIRICHLET_CSTR_CHARGE_CUBE%APPEND")
         my_pos_cube = "REWIND"
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         END IF
         mpi_io = .TRUE.
         unit_nr = cp_print_key_unit_nr(logger, input, &
                                        "DFT%PRINT%IMPLICIT_PSOLVER%DIRICHLET_CSTR_CHARGE_CUBE", &
                                        extension=".cube", middle_name="dirichlet_cstr_charge", file_position=my_pos_cube, &
                                        mpi_io=mpi_io)
         CALL pw_env_get(pw_env, poisson_env=poisson_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL pw_pool_create_pw(auxbas_pw_pool, aux_r%pw, use_data=REALDATA3D, in_space=REALSPACE)

         boundary_condition = pw_env%poisson_env%parameters%ps_implicit_params%boundary_condition
         SELECT CASE (boundary_condition)
         CASE (MIXED_PERIODIC_BC)
            CALL pw_copy(poisson_env%implicit_env%cstr_charge, aux_r%pw)
         CASE (MIXED_BC)
            CALL pw_shrink(pw_env%poisson_env%parameters%ps_implicit_params%neumann_directions, &
                           pw_env%poisson_env%implicit_env%dct_env%dests_shrink, &
                           pw_env%poisson_env%implicit_env%dct_env%srcs_shrink, &
                           pw_env%poisson_env%implicit_env%dct_env%bounds_local_shftd, &
                           poisson_env%implicit_env%cstr_charge, aux_r%pw)
         END SELECT

         CALL cp_pw_to_cube(aux_r%pw, unit_nr, "DIRICHLET CONSTRAINT CHARGE", particles=particles, &
                            stride=section_get_ivals(dft_section, "PRINT%IMPLICIT_PSOLVER%DIRICHLET_CSTR_CHARGE_CUBE%STRIDE"), &
                            mpi_io=mpi_io)
         CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                           "DFT%PRINT%IMPLICIT_PSOLVER%DIRICHLET_CSTR_CHARGE_CUBE", mpi_io=mpi_io)

         CALL pw_pool_give_back_pw(auxbas_pw_pool, aux_r%pw)
      ENDIF

      ! Write Dirichlet type constranits into cube files
      do_dirichlet_bc_cube = BTEST(cp_print_key_should_output(logger%iter_info, input, &
                                                              "DFT%PRINT%IMPLICIT_PSOLVER%DIRICHLET_BC_CUBE"), cp_p_file)
      has_dirichlet_bc = .FALSE.
      IF (has_implicit_ps) THEN
         boundary_condition = pw_env%poisson_env%parameters%ps_implicit_params%boundary_condition
         IF (boundary_condition .EQ. MIXED_PERIODIC_BC .OR. boundary_condition .EQ. MIXED_BC) THEN
            has_dirichlet_bc = .TRUE.
         END IF
      END IF

      IF (has_implicit_ps .AND. has_dirichlet_bc .AND. do_dirichlet_bc_cube) THEN
         append_cube = section_get_lval(input, "DFT%PRINT%IMPLICIT_PSOLVER%DIRICHLET_BC_CUBE%APPEND")
         my_pos_cube = "REWIND"
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         END IF
         tile_cubes = section_get_lval(input, "DFT%PRINT%IMPLICIT_PSOLVER%DIRICHLET_BC_CUBE%TILE_CUBES")

         CALL pw_env_get(pw_env, poisson_env=poisson_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL pw_pool_create_pw(auxbas_pw_pool, aux_r%pw, use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_zero(aux_r%pw)

         IF (tile_cubes) THEN
            ! one cube file per tile
            n_cstr = SIZE(poisson_env%implicit_env%contacts)
            DO j = 1, n_cstr
               n_tiles = poisson_env%implicit_env%contacts(j)%dirichlet_bc%n_tiles
               DO i = 1, n_tiles
                  filename = "dirichlet_cstr_"//TRIM(ADJUSTL(cp_to_string(j)))// &
                             "_tile_"//TRIM(ADJUSTL(cp_to_string(i)))
                  mpi_io = .TRUE.
                  unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%IMPLICIT_PSOLVER%DIRICHLET_BC_CUBE", &
                                                 extension=".cube", middle_name=filename, file_position=my_pos_cube, &
                                                 mpi_io=mpi_io)

                  CALL pw_copy(poisson_env%implicit_env%contacts(j)%dirichlet_bc%tiles(i)%tile%tile_pw, aux_r%pw)

                  CALL cp_pw_to_cube(aux_r%pw, unit_nr, "DIRICHLET TYPE CONSTRAINT", particles=particles, &
                                     stride=section_get_ivals(dft_section, "PRINT%IMPLICIT_PSOLVER%DIRICHLET_BC_CUBE%STRIDE"), &
                                     mpi_io=mpi_io)
                  CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                                    "DFT%PRINT%IMPLICIT_PSOLVER%DIRICHLET_BC_CUBE", mpi_io=mpi_io)
               END DO
            END DO
         ELSE
            ! a single cube file
            NULLIFY (dirichlet_tile)
            CALL pw_pool_create_pw(auxbas_pw_pool, dirichlet_tile, use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_zero(dirichlet_tile)
            mpi_io = .TRUE.
            unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%IMPLICIT_PSOLVER%DIRICHLET_BC_CUBE", &
                                           extension=".cube", middle_name="DIRICHLET_CSTR", file_position=my_pos_cube, &
                                           mpi_io=mpi_io)

            n_cstr = SIZE(poisson_env%implicit_env%contacts)
            DO j = 1, n_cstr
               n_tiles = poisson_env%implicit_env%contacts(j)%dirichlet_bc%n_tiles
               DO i = 1, n_tiles
                  CALL pw_copy(poisson_env%implicit_env%contacts(j)%dirichlet_bc%tiles(i)%tile%tile_pw, dirichlet_tile)
                  CALL pw_axpy(dirichlet_tile, aux_r%pw)
               END DO
            END DO

            CALL cp_pw_to_cube(aux_r%pw, unit_nr, "DIRICHLET TYPE CONSTRAINT", particles=particles, &
                               stride=section_get_ivals(dft_section, "PRINT%IMPLICIT_PSOLVER%DIRICHLET_BC_CUBE%STRIDE"), &
                               mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr, logger, input, &
                                              "DFT%PRINT%IMPLICIT_PSOLVER%DIRICHLET_BC_CUBE", mpi_io=mpi_io)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, dirichlet_tile)
         END IF

         CALL pw_pool_give_back_pw(auxbas_pw_pool, aux_r%pw)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_post_ps_implicit

!**************************************************************************************************
!> \brief write an adjacency (interaction) matrix
!> \param qs_env qs environment
!> \param input the input
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE write_adjacency_matrix(qs_env, input)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input

      CHARACTER(len=*), PARAMETER :: routineN = 'write_adjacency_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: adjm_size, colind, handle, iatom, ikind, &
                                                            ind, jatom, jkind, k, natom, nkind, &
                                                            output_unit, rowind, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: interact_adjm
      LOGICAL                                            :: do_adjm_write, do_symmetric
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_list_a, basis_set_list_b
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: dft_section

      CALL timeset(routineN, handle)

      NULLIFY (dft_section)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      dft_section => section_vals_get_subs_vals(input, "DFT")
      do_adjm_write = BTEST(cp_print_key_should_output(logger%iter_info, dft_section, &
                                                       "PRINT%ADJMAT_WRITE"), cp_p_file)

      IF (do_adjm_write) THEN
         NULLIFY (qs_kind_set, nl_iterator)
         NULLIFY (basis_set_list_a, basis_set_list_b, basis_set_a, basis_set_b)

         CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, sab_orb=nl, natom=natom, para_env=para_env)

         nkind = SIZE(qs_kind_set)
         CPASSERT(SIZE(nl) .GT. 0)
         CALL get_neighbor_list_set_p(neighbor_list_sets=nl, symmetric=do_symmetric)
         CPASSERT(do_symmetric)
         ALLOCATE (basis_set_list_a(nkind), basis_set_list_b(nkind))
         CALL basis_set_list_setup(basis_set_list_a, "ORB", qs_kind_set)
         CALL basis_set_list_setup(basis_set_list_b, "ORB", qs_kind_set)

         adjm_size = ((natom + 1)*natom)/2
         ALLOCATE (interact_adjm(4*adjm_size))
         interact_adjm = 0

         NULLIFY (nl_iterator)
         CALL neighbor_list_iterator_create(nl_iterator, nl)
         DO WHILE (neighbor_list_iterate(nl_iterator) .EQ. 0)
            CALL get_iterator_info(nl_iterator, &
                                   ikind=ikind, jkind=jkind, &
                                   iatom=iatom, jatom=jatom)

            basis_set_a => basis_set_list_a(ikind)%gto_basis_set
            IF (.NOT. ASSOCIATED(basis_set_a)) CYCLE
            basis_set_b => basis_set_list_b(jkind)%gto_basis_set
            IF (.NOT. ASSOCIATED(basis_set_b)) CYCLE

            ! move everything to the upper triangular part
            IF (iatom .LE. jatom) THEN
               rowind = iatom
               colind = jatom
            ELSE
               rowind = jatom
               colind = iatom
               ! swap the kinds too
               ikind = ikind + jkind
               jkind = ikind - jkind
               ikind = ikind - jkind
            END IF

            ! indexing upper triangular matrix
            ind = adjm_size - (natom - rowind + 1)*((natom - rowind + 1) + 1)/2 + colind - rowind + 1
            ! convert the upper triangular matrix into a adjm_size x 4 matrix
            ! columns are: iatom, jatom, ikind, jkind
            interact_adjm((ind - 1)*4 + 1) = rowind
            interact_adjm((ind - 1)*4 + 2) = colind
            interact_adjm((ind - 1)*4 + 3) = ikind
            interact_adjm((ind - 1)*4 + 4) = jkind
         ENDDO

         CALL mp_sum(interact_adjm, para_env%group)

         unit_nr = cp_print_key_unit_nr(logger, dft_section, "PRINT%ADJMAT_WRITE", &
                                        extension=".adjmat", file_form="FORMATTED", &
                                        file_status="REPLACE")
         IF (unit_nr .GT. 0) THEN
            WRITE (unit_nr, "(1A,2X,1A,5X,1A,4X,A5,3X,A5)") "#", "iatom", "jatom", "ikind", "jkind"
            DO k = 1, 4*adjm_size, 4
               ! print only the interacting atoms
               IF (interact_adjm(k) .GT. 0 .AND. interact_adjm(k + 1) .GT. 0) THEN
                  WRITE (unit_nr, "(I8,2X,I8,3X,I6,2X,I6)") interact_adjm(k:k + 3)
               END IF
            END DO
         END IF

         CALL cp_print_key_finished_output(unit_nr, logger, dft_section, "PRINT%ADJMAT_WRITE")

         CALL neighbor_list_iterator_release(nl_iterator)
         DEALLOCATE (basis_set_list_a, basis_set_list_b)
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_adjacency_matrix

! **************************************************************************************************
!> \brief Updates Hartree potential with MP2 density. Important for REPEAT charges
!> \param rho ...
!> \param qs_env ...
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE update_hartree_with_mp2(rho, qs_env)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_environment_type), POINTER                 :: qs_env

      LOGICAL                                            :: use_virial
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), POINTER                           :: rho_core
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(virial_type), POINTER                         :: virial

      NULLIFY (auxbas_pw_pool, pw_env, poisson_env, energy, rho_core, v_hartree_rspace%pw, virial)
      CALL get_qs_env(qs_env, pw_env=pw_env, energy=energy, &
                      rho_core=rho_core, virial=virial, &
                      v_hartree_rspace=v_hartree_rspace%pw)

      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      IF (.NOT. use_virial) THEN

         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                         poisson_env=poisson_env)
         CALL pw_pool_create_pw(auxbas_pw_pool, &
                                v_hartree_gspace%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, &
                                rho_tot_gspace%pw, &
                                use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)

         CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)
         CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, energy%hartree, &
                               v_hartree_gspace%pw, rho_core=rho_core)

         CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
         CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
      ENDIF

   END SUBROUTINE update_hartree_with_mp2

END MODULE qs_scf_post_gpw
