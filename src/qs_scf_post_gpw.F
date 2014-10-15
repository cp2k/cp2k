!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Does all kind of post scf calculations for GPW/GAPW
!> \par History
!>      Started as a copy from the relevant part of qs_scf
!> \author Joost VandeVondele (10.2003)
! *****************************************************************************
MODULE qs_scf_post_gpw
  USE admm_types,                      ONLY: admm_type
  USE admm_utils,                      ONLY: admm_correct_for_eigenvalues,&
                                             admm_uncorrect_for_eigenvalues
  USE ai_onecenter,                    ONLY: sg_overlap
  USE atom_kind_orbitals,              ONLY: calculate_atomic_density
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_methods,                    ONLY: cell_type,&
                                             get_cell,&
                                             pbc
  USE cp_array_r_utils,                ONLY: cp_1d_r_p_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
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
  USE cp_output_handling,              ONLY: cp_iter_string,&
                                             cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
  USE cp_result_methods,               ONLY: cp_results_erase,&
                                             get_results,&
                                             put_results
  USE cp_result_types,                 ONLY: cp_result_type
  USE et_coupling_types,               ONLY: set_et_coupling_type
  USE f77_blas
  USE hfx_ri_gemopt,                   ONLY: geminal_optimize
  USE input_constants,                 ONLY: do_loc_both,&
                                             do_loc_homo,&
                                             do_loc_lumo,&
                                             ot_precond_full_all,&
                                             use_orb_basis_set
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
  USE lapack,                          ONLY: lapack_sgesv
  USE mathconstants,                   ONLY: pi,&
                                             twopi
  USE message_passing,                 ONLY: mp_sum
  USE molecular_states,                ONLY: construct_molecular_states
  USE molecule_types_new,              ONLY: molecule_type
  USE moments_utils,                   ONLY: get_reference_point
  USE orbital_pointers,                ONLY: indso
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom,&
                                             debye,&
                                             evolt
  USE population_analyses,             ONLY: lowdin_population_analysis,&
                                             mulliken_population_analysis
  USE preconditioner_types,            ONLY: preconditioner_type
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
  USE qs_collocate_density,            ONLY: calculate_rho_elec,&
                                             calculate_wavefunction
  USE qs_conductivity,                 ONLY: optical_conductivity
  USE qs_core_energies,                ONLY: calculate_ptrace
  USE qs_electric_field_gradient,      ONLY: qs_efg_calc
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_epr_hyp,                      ONLY: qs_epr_hyp_calc
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE qs_kind_types,                   ONLY: get_qs_kind,&
                                             qs_kind_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_loc_methods,                  ONLY: qs_loc_driver
  USE qs_loc_molecules,                ONLY: wfc_to_molecule
  USE qs_loc_types,                    ONLY: qs_loc_env_create,&
                                             qs_loc_env_destroy,&
                                             qs_loc_env_new_type
  USE qs_loc_utils,                    ONLY: loc_write_restart,&
                                             qs_loc_control_init,&
                                             qs_loc_env_init,&
                                             qs_loc_init,&
                                             retain_history
  USE qs_mo_io,                        ONLY: write_mo_set
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues,&
                                             make_mo_eig
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             set_mo_occupation
  USE qs_moments,                      ONLY: qs_moment_berry_phase,&
                                             qs_moment_locop
  USE qs_ot_eigensolver,               ONLY: ot_eigensolver
  USE qs_pdos,                         ONLY: calculate_projected_dos
  USE qs_resp,                         ONLY: resp_fit
  USE qs_rho_atom_types,               ONLY: rho_atom_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_types,                    ONLY: ot_method_nr,&
                                             qs_scf_env_type
  USE qs_scf_wfn_mix,                  ONLY: wfn_mix
  USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                             qs_subsys_type
  USE s_square_methods,                ONLY: compute_s_square
  USE scf_control_types,               ONLY: scf_control_type
  USE stm_images,                      ONLY: th_stm_image
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE wannier_states,                  ONLY: construct_wannier_states
  USE wannier_states_types,            ONLY: wannier_centres_type
  USE xray_diffraction,                ONLY: calculate_rhotot_elec_gspace,&
                                             xray_diffraction_spectrum
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  ! Global parameters
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_post_gpw'
  PUBLIC :: scf_post_calculation_gpw,&
            qs_scf_post_occ_cubes,&
            qs_scf_post_unocc_cubes,&
            qs_scf_post_moments,&
            write_available_results,&
            write_mo_free_results

CONTAINS

! *****************************************************************************
!> \brief collects possible post - scf calculations and prints info / computes properties.
!> \param dft_section ...
!> \param scf_env the scf_env whose info should be written out
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2003 created [fawzi]
!>      10.2004 moved here from qs_scf [Joost VandeVondele]
!>              started splitting out different subroutines
!> \author fawzi
!> \note
!>      this function changes mo_eigenvectors and mo_eigenvalues, depending on the print keys.
!>      In particular, MO_CUBES causes the MOs to be rotated to make them eigenstates of the KS
!>      matrix, and mo_eigenvalues is updated accordingly. This can, for unconverged wavefunctions,
!>      change afterwards slightly the forces (hence small numerical differences between MD
!>      with and without the debug print level). Ideally this should not happen...
! *****************************************************************************
  SUBROUTINE scf_post_calculation_gpw(dft_section, scf_env,qs_env, error)

    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_post_calculation_gpw', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, homo, ispin, istat, min_lumos, n_rep, nhomo, nlumo, &
      nlumo_stm, nlumo_tddft, nlumos, nmo, output_unit, unit_nr
    INTEGER, DIMENSION(:, :, :), POINTER     :: marked_states
    LOGICAL :: check_write, compute_lumos, do_homo, do_mo_cubes, do_stm, &
      do_wannier_cubes, failure, has_homo, has_lumo, loc_explicit, &
      loc_print_explicit, my_localized_wfn, p_loc, p_loc_homo, p_loc_lumo
    REAL(dp)                                 :: e_kin
    REAL(KIND=dp)                            :: gap, homo_lumo(2,2)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
      POINTER                                :: occupied_evals, &
                                                unoccupied_evals, &
                                                unoccupied_evals_stm
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: kinetic_m, ks_rmpv, matrix_s, &
                                                mo_derivs
    TYPE(cp_fm_p_type), DIMENSION(:), POINTER :: homo_localized, &
      lumo_localized, lumo_ptr, mo_loc_history, occupied_orbs, &
      unoccupied_orbs, unoccupied_orbs_stm
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: wf_g, wf_r
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env_homo, &
                                                qs_loc_env_lumo
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: input, loc_print_section, &
                                                localize_section, print_key, &
                                                stm_section

    CALL timeset(routineN,handle)

    ! Writes the data that is already available in qs_env
    CALL write_available_results(qs_env,scf_env,error)

    failure=.FALSE.
    my_localized_wfn = .FALSE.
    NULLIFY(admm_env, dft_control, pw_env, auxbas_pw_pool, pw_pools, mos, rho, &
         mo_coeff, ks_rmpv, matrix_s, qs_loc_env_homo,qs_loc_env_lumo, scf_control, &
         unoccupied_orbs,unoccupied_orbs_stm,  mo_eigenvalues, unoccupied_evals, unoccupied_evals_stm,  molecule_set, mo_derivs,&
         subsys, particles, input, print_key, kinetic_m,marked_states)
    NULLIFY(homo_localized, lumo_localized,lumo_ptr)

    has_homo=.FALSE.
    has_lumo=.FALSE.
    p_loc = .FALSE.
    p_loc_homo = .FALSE.
    p_loc_lumo= .FALSE.


    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    ! Here we start with data that needs a postprocessing...
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,dft_control=dft_control,molecule_set=molecule_set, &
            mos=mos,rho=rho,matrix_ks=ks_rmpv,scf_control=scf_control,matrix_s=matrix_s, &
            input=input, subsys=subsys, pw_env=pw_env,&
            particle_set=particle_set, atomic_kind_set=atomic_kind_set, &
            qs_kind_set=qs_kind_set, error=error)
       CALL qs_subsys_get(subsys,particles=particles,error=error)

       !    **** the kinetic energy
       IF (cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%KINETIC_ENERGY",error=error)/=0) THEN
          CALL get_qs_env(qs_env,kinetic=kinetic_m,error=error)
          CPPrecondition(ASSOCIATED(kinetic_m),cp_failure_level,routineP,error,failure)
          CPPrecondition(ASSOCIATED(kinetic_m(1)%matrix),cp_failure_level,routineP,error,failure)
          CALL calculate_ptrace(kinetic_m,rho%rho_ao,e_kin,dft_control%nspins,error)
          unit_nr = cp_print_key_unit_nr(logger,input,"DFT%PRINT%KINETIC_ENERGY",&
               extension=".Log",error=error)
          IF (unit_nr>0) THEN
             WRITE (unit_nr,'(T3,A,T55,F25.14)') "Electronic kinetic energy:",e_kin
          ENDIF
          CALL cp_print_key_finished_output(unit_nr,logger,input,&
               "DFT%PRINT%KINETIC_ENERGY", error=error)
       END IF

       ! Atomic Charges that require further computation
       CALL qs_scf_post_charges(input, logger, qs_env, rho, matrix_s, error=error)

       ! Moments of charge distribution
       CALL qs_scf_post_moments(input, logger, qs_env, output_unit, error=error)

       ! Determine if we need to computer properties using the localized centers
       localize_section     => section_vals_get_subs_vals(dft_section,"LOCALIZE",error=error)
       loc_print_section => section_vals_get_subs_vals(localize_section,"PRINT",error=error)
       CALL section_vals_get(localize_section, explicit=loc_explicit, error=error)
       CALL section_vals_get(loc_print_section, explicit=loc_print_explicit, error=error)

       ! Print_keys controlled by localization
       IF(loc_print_explicit) THEN
         print_key => section_vals_get_subs_vals(loc_print_section,"MOLECULAR_DIPOLES",error=error)
         p_loc=BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)
         print_key => section_vals_get_subs_vals(loc_print_section,"TOTAL_DIPOLE",error=error)
         p_loc=p_loc.OR.BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)
         print_key => section_vals_get_subs_vals(loc_print_section,"WANNIER_CENTERS",error=error)
         p_loc=p_loc.OR.BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)
         print_key => section_vals_get_subs_vals(loc_print_section,"WANNIER_SPREADS",error=error)
         p_loc=p_loc.OR.BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)
         print_key => section_vals_get_subs_vals(loc_print_section,"WANNIER_CUBES",error=error)
         p_loc=p_loc.OR.BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)
         print_key => section_vals_get_subs_vals(loc_print_section,"MOLECULAR_STATES",error=error)
         p_loc=p_loc.OR.BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)
       ELSE
         p_loc=.FALSE.
       END IF
       IF(loc_explicit) THEN
         p_loc_homo=(section_get_ival(localize_section,"STATES",error=error)==do_loc_homo.OR.&
                  section_get_ival(localize_section,"STATES",error=error)==do_loc_both).AND.p_loc
         p_loc_lumo=(section_get_ival(localize_section,"STATES",error=error)==do_loc_lumo.OR.&
                  section_get_ival(localize_section,"STATES",error=error)==do_loc_both).AND.p_loc
         CALL section_vals_val_get(localize_section,"LIST_UNOCCUPIED",  n_rep_val=n_rep,error=error)
       ELSE
         p_loc_homo=.FALSE.
         p_loc_lumo=.FALSE.
         n_rep=0
       END IF

       IF(n_rep==0.AND.p_loc_lumo)THEN
           CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
              "No LIST_UNOCCUPIED was specified, therefore localization of unoccupied states will be skipped!"//&
              CPSourceFileRef,&
              only_ionode=.TRUE.)
           p_loc_lumo=.FALSE.
       END IF
       print_key => section_vals_get_subs_vals(loc_print_section,"WANNIER_STATES",error=error)
       p_loc=p_loc.OR.BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)


       ! Control for STM
       stm_section => section_vals_get_subs_vals(input,"DFT%PRINT%STM", error=error)
       CALL section_vals_get(stm_section,explicit=do_stm,error=error)
! BTEST(cp_print_key_should_output(logger%iter_info,stm_section,error=error),cp_p_file)
       nlumo_stm = 0
       IF(do_stm) nlumo_stm = section_get_ival(stm_section,"NLUMO",error=error)

       ! check for CUBES (MOs and WANNIERS)
       do_mo_cubes=BTEST(cp_print_key_should_output(logger%iter_info,dft_section,"PRINT%MO_CUBES",&
                         error=error),cp_p_file)
       IF(loc_print_explicit) THEN
         do_wannier_cubes=BTEST(cp_print_key_should_output(logger%iter_info,loc_print_section,&
                        "WANNIER_CUBES",error=error),cp_p_file)
       ELSE
         do_wannier_cubes=.FALSE.
       END IF
       nlumo=section_get_ival(dft_section,"PRINT%MO_CUBES%NLUMO",error=error)
       nhomo=section_get_ival(dft_section,"PRINT%MO_CUBES%NHOMO",error=error)
       nlumo_tddft = 0
       IF(dft_control%do_tddfpt_calculation) THEN
         nlumo_tddft=section_get_ival(dft_section,"TDDFPT%NLUMO",error=error)
       END IF

       ! Setup the grids needed to compute a wavefunction given a vector..
       IF ( ( ( do_mo_cubes .OR. do_wannier_cubes ).AND. (nlumo /= 0 .OR. nhomo /= 0  )) .OR. p_loc ) THEN
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                          pw_pools=pw_pools,error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,&
                                  use_data = REALDATA3D,&
                                  in_space = REALSPACE, error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,wf_g%pw,&
                                  use_data = COMPLEXDATA1D,&
                                  in_space = RECIPROCALSPACE, error=error)
       END IF

       ! Makes the MOs eigenstates, computes eigenvalues, write cubes
       IF ((do_mo_cubes .AND. nhomo /= 0) .OR. do_stm .OR. dft_control%do_tddfpt_calculation ) THEN
          IF (dft_control%restricted) THEN
            IF (output_unit>0) WRITE(output_unit,*) &
              " Unclear how we define MOs in the restricted case ... skipping"
          ELSE
            CALL  get_qs_env(qs_env,mo_derivs=mo_derivs,error=error)
            IF(dft_control%do_admm) THEN
              CALL get_qs_env(qs_env,admm_env=admm_env, error=error)
              CALL make_mo_eig(mos,dft_control%nspins,ks_rmpv,scf_control,mo_derivs,admm_env=admm_env,&
                 error=error)
            ELSE
               CALL make_mo_eig(mos,dft_control%nspins,ks_rmpv,scf_control,mo_derivs,error=error)
            END IF
          END IF 
          DO ispin=1,dft_control%nspins
             CALL get_mo_set(mo_set=mos(ispin)%mo_set, eigenvalues=mo_eigenvalues,homo=homo)
             homo_lumo(ispin,1)=mo_eigenvalues(homo)
          END DO
          has_homo=.TRUE.
       END IF
       IF (do_mo_cubes .AND. nhomo /= 0) THEN
          DO ispin=1,dft_control%nspins
             ! Prints the cube files of OCCUPIED ORBITALS
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                  eigenvalues=mo_eigenvalues,homo=homo,nmo=nmo)
             CALL qs_scf_post_occ_cubes(input, dft_section, dft_control, logger, qs_env,&
                  mo_coeff, wf_g, wf_r, particles, homo, ispin, error=error)
          END DO
       ENDIF

       ! Initialize the localization environment, needed e.g. for wannier functions and molecular states
       ! Gets localization info for the occupied orbs
       !  - Possibly gets wannier functions
       !  - Possibly gets molecular states
       IF (p_loc_homo) THEN
          IF (dft_control%restricted) THEN
             IF (output_unit>0) WRITE(output_unit,*) &
                  " Unclear how we define MOs / localization in the restricted case ... skipping"
          ELSE
             ALLOCATE(occupied_orbs(dft_control%nspins),stat=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(occupied_evals(dft_control%nspins),stat=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             ALLOCATE(homo_localized(dft_control%nspins),stat=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             DO ispin=1,dft_control%nspins
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                     eigenvalues=mo_eigenvalues)
                occupied_orbs(ispin)%matrix=>mo_coeff
                occupied_evals(ispin)%array=>mo_eigenvalues
                CALL cp_fm_create(homo_localized(ispin)%matrix,occupied_orbs(ispin)%matrix%matrix_struct,error=error)
                CALL cp_fm_to_fm(occupied_orbs(ispin)%matrix,homo_localized(ispin)%matrix,error=error)
             END DO

             CALL get_qs_env(qs_env,mo_loc_history=mo_loc_history,error=error)
             do_homo=.TRUE.

             CALL qs_loc_env_create(qs_loc_env_homo,error=error)
             CALL qs_loc_control_init(qs_loc_env_homo,localize_section,do_homo=do_homo,error=error)
             CALL qs_loc_init(qs_env,qs_loc_env_homo,localize_section,homo_localized,do_homo,&
                  do_mo_cubes,mo_loc_history=mo_loc_history,error=error)
             CALL get_localization_info(qs_env,qs_loc_env_homo,localize_section,homo_localized,&
                  wf_r, wf_g,particles,occupied_orbs,occupied_evals,marked_states,error=error)

             !retain the homo_localized for future use
             IF (qs_loc_env_homo%localized_wfn_control%use_history) THEN
                CALL retain_history(mo_loc_history,homo_localized,error)
                CALL set_qs_env(qs_env,mo_loc_history=mo_loc_history,error=error)
             ENDIF

             !write restart for localization of occupied orbitals
             CALL loc_write_restart(qs_env,qs_loc_env_homo,loc_print_section,mos,&
                  homo_localized, do_homo, error=error)
             CALL cp_fm_vect_dealloc(homo_localized,error)
             DEALLOCATE(occupied_orbs,stat=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             DEALLOCATE(occupied_evals,stat=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             ! Print Total Dipole if the localization has been performed
             IF (qs_loc_env_homo%do_localize) THEN
                CALL qs_scf_post_loc_dip(input, dft_control, qs_loc_env_homo, logger, qs_env, error)
             END IF
          END IF
       ENDIF

       ! Gets the lumos, and eigenvalues for the lumos, and localize them if requested
       nlumo = MAX(nlumo, nlumo_tddft)
       compute_lumos=(do_mo_cubes .OR.dft_control%do_tddfpt_calculation ) .AND. nlumo .NE. 0
       compute_lumos=compute_lumos.OR.p_loc_lumo

       DO ispin = 1,dft_control%nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, homo=homo, nmo=nmo)
         compute_lumos=compute_lumos.AND.homo==nmo
       END DO

       IF(do_mo_cubes .AND. .NOT. compute_lumos) THEN

         nlumo=section_get_ival(dft_section,"PRINT%MO_CUBES%NLUMO",error=error)
         DO ispin = 1,dft_control%nspins

            CALL get_mo_set(mo_set=mos(ispin)%mo_set, homo=homo, nmo=nmo, eigenvalues=mo_eigenvalues)
            IF(nlumo > nmo-homo) THEN
         ! this case not yet implemented
            ELSE
              IF (output_unit>0) WRITE(output_unit,*) " "
              IF (output_unit>0) WRITE(output_unit,*) " Lowest eigenvalues of the unoccupied subspace spin ",ispin
              IF (output_unit>0) WRITE(output_unit,*) "---------------------------------------------"
              IF (output_unit>0) WRITE(output_unit,'(4(1X,1F16.8))') mo_eigenvalues(homo+1:homo+nlumo)

             ! Prints the cube files of UNOCCUPIED ORBITALS
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff) 
             CALL qs_scf_post_unocc_cubes(input, dft_section, dft_control, logger, qs_env,&
                  mo_coeff, wf_g, wf_r, particles, nlumo, homo, ispin, lumo=homo+1, error=error)
            END IF
         END DO

       END IF

       IF (compute_lumos) THEN
          check_write=.TRUE.
          min_lumos=nlumo
          IF(nlumo==0)check_write=.FALSE.
          IF(p_loc_lumo)THEN
             do_homo=.FALSE.
             CALL qs_loc_env_create(qs_loc_env_lumo,error=error)
             CALL qs_loc_control_init(qs_loc_env_lumo,localize_section,do_homo=do_homo,error=error)
             min_lumos=MAX(MAXVAL(qs_loc_env_lumo%localized_wfn_control%loc_states(:,:)),nlumo)
          END IF

          ALLOCATE(unoccupied_orbs(dft_control%nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(unoccupied_evals(dft_control%nspins),STAT=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          CALL  make_lumo(qs_env,scf_env,unoccupied_orbs,unoccupied_evals,min_lumos,nlumos,error)
          lumo_ptr=>unoccupied_orbs
          DO ispin=1,dft_control%nspins
             has_lumo=.TRUE.
             homo_lumo(ispin,2)=unoccupied_evals(ispin)%array(1)
             CALL get_mo_set(mo_set=mos(ispin)%mo_set, homo=homo)
             IF(check_write)THEN
                IF(p_loc_lumo.AND.nlumo.NE.-1)nlumos=MIN(nlumo,nlumos)
                ! Prints the cube files of UNOCCUPIED ORBITALS
                CALL qs_scf_post_unocc_cubes(input, dft_section, dft_control, logger, qs_env,&
                     unoccupied_orbs(ispin)%matrix, wf_g, wf_r, particles, nlumos, homo, ispin, error=error)
             END IF
          END DO

          ! Save the info for tddfpt calculation
          IF (dft_control%do_tddfpt_calculation) THEN
             ALLOCATE(dft_control%tddfpt_control%lumos_eigenvalues(nlumos,dft_control%nspins),stat=istat)
             DO ispin=1, dft_control%nspins
                dft_control%tddfpt_control%lumos_eigenvalues(1:nlumos,ispin) = &
                     unoccupied_evals(ispin)%array(1:nlumos)
             END DO
             dft_control%tddfpt_control%lumos => unoccupied_orbs
          END IF

          IF(p_loc_lumo)THEN
             ALLOCATE(lumo_localized(dft_control%nspins),STAT=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             DO ispin=1,dft_control%nspins
                CALL cp_fm_create(lumo_localized(ispin)%matrix,unoccupied_orbs(ispin)%matrix%matrix_struct,error=error)
                CALL cp_fm_to_fm(unoccupied_orbs(ispin)%matrix,lumo_localized(ispin)%matrix,error=error)
             END DO
             CALL qs_loc_init(qs_env,qs_loc_env_lumo,localize_section,lumo_localized,do_homo,do_mo_cubes,&
                  evals=unoccupied_evals, error=error)
             CALL qs_loc_env_init(qs_loc_env_lumo,qs_loc_env_lumo%localized_wfn_control,qs_env,&
                  loc_coeff=unoccupied_orbs,error=error)
             CALL get_localization_info(qs_env,qs_loc_env_lumo,localize_section, &
                  lumo_localized, wf_r, wf_g,particles,&
                  unoccupied_orbs,unoccupied_evals,marked_states,error=error)
             CALL loc_write_restart(qs_env,qs_loc_env_lumo,loc_print_section,mos,homo_localized, do_homo,&
                   evals=unoccupied_evals, error=error)
             lumo_ptr=>lumo_localized
          END IF
       ENDIF

       IF (has_homo .AND. has_lumo) THEN
         IF (output_unit>0) WRITE(output_unit,*) " "
         DO ispin=1,dft_control%nspins
             IF(.NOT. scf_control%smear%do_smear) THEN
                gap = homo_lumo(ispin,2)-homo_lumo(ispin,1)
                IF (output_unit>0) WRITE(output_unit,'(T2,A,F12.6)')&
                                   "HOMO - LUMO gap [eV] :",gap*evolt
             END IF
         ENDDO
       ENDIF

       ! Deallocate grids needed to compute wavefunctions
       IF ( ( ( do_mo_cubes .OR. do_wannier_cubes ).AND. (nlumo /= 0 .OR. nhomo /= 0 )) .OR. p_loc ) THEN
          CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_r%pw, error=error)
          CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_g%pw, error=error)
       END IF

       ! Destroy the localization environment
       IF (p_loc_homo)  CALL qs_loc_env_destroy(qs_loc_env_homo, error=error)
       IF (p_loc_lumo)  CALL qs_loc_env_destroy(qs_loc_env_lumo, error=error)



       ! generate a mix of wfns, and write to a restart
       CALL wfn_mix(mos, particle_set, dft_section, atomic_kind_set, qs_kind_set, &
                    lumo_ptr, scf_env, matrix_s, output_unit, marked_states,error)
       IF(p_loc_lumo)CALL cp_fm_vect_dealloc(lumo_localized,error)
       IF(ASSOCIATED(marked_states))THEN
          DEALLOCATE(marked_states,stat=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF

       IF(compute_lumos) THEN
       ! Compute the optical conductivity (needs homos and lumos)
         CALL qs_scf_post_optc(input, dft_section, dft_control, logger, qs_env,&
              unoccupied_orbs, unoccupied_evals, output_unit, error)
       END IF

       ! This is just a deallocation for printing MO_CUBES or TDDFPT
       IF(compute_lumos) THEN
          DO ispin=1,dft_control%nspins
             DEALLOCATE(unoccupied_evals(ispin)%array,stat=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             IF (.NOT.dft_control%do_tddfpt_calculation) &
                  CALL cp_fm_release(unoccupied_orbs(ispin)%matrix,error=error)
          ENDDO
          DEALLOCATE(unoccupied_evals,stat=istat)
          CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          IF (.NOT.dft_control%do_tddfpt_calculation) THEN
            DEALLOCATE(unoccupied_orbs,stat=istat)
            CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
          END IF
       ENDIF

       !stm images
       IF(do_stm)THEN
         ALLOCATE(unoccupied_orbs_stm(dft_control%nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         ALLOCATE(unoccupied_evals_stm(dft_control%nspins),STAT=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         CALL  make_lumo(qs_env,scf_env,unoccupied_orbs_stm,unoccupied_evals_stm,nlumo_stm,nlumos,error=error)

         CALL th_stm_image(qs_env, stm_section, particles, unoccupied_orbs_stm, &
                         unoccupied_evals_stm, error=error)

         DO ispin=1,dft_control%nspins
             DEALLOCATE(unoccupied_evals_stm(ispin)%array,stat=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             CALL cp_fm_release(unoccupied_orbs_stm(ispin)%matrix,error=error)
         ENDDO
         DEALLOCATE(unoccupied_evals_stm,stat=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
         DEALLOCATE(unoccupied_orbs_stm,stat=istat)
         CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       END IF


       ! Print coherent X-ray diffraction spectrum
       CALL qs_scf_post_xray(input, dft_section, logger, qs_env, output_unit, error)

       ! Calculation of Electric Field Gradients
       CALL qs_scf_post_efg(input, logger, qs_env, error)

       ! Calculation of ET
       CALL qs_scf_post_et(input, qs_env, dft_control, error)

       ! Calculation of EPR Hyperfine Coupling Tensors
       CALL qs_scf_post_epr(input, logger, qs_env, error)

       ! Calculation of properties needed for BASIS_MOLOPT optimizations
       CALL qs_scf_post_molopt(input, logger, qs_env, error)

       ! Optimize Geminals
       CALL qs_scf_post_gemopt(input, logger, qs_env, error)

       ! Calculate ELF
       CALL qs_scf_post_elf(input, logger, qs_env, error)
    END IF
    CALL timestop(handle)


  END SUBROUTINE scf_post_calculation_gpw

! *****************************************************************************
!> \brief Performs localization of the orbitals
!> \param qs_env ...
!> \param qs_loc_env ...
!> \param loc_section ...
!> \param mo_local ...
!> \param wf_r ...
!> \param wf_g ...
!> \param particles ...
!> \param coeff ...
!> \param evals ...
!> \param marked_states ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_localization_info(qs_env,qs_loc_env,loc_section,mo_local,&
                                   wf_r, wf_g,particles,coeff,evals,marked_states,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(section_vals_type), POINTER         :: loc_section
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_local
    TYPE(pw_p_type)                          :: wf_r, wf_g
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: coeff
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
      POINTER                                :: evals
    INTEGER, DIMENSION(:, :, :), POINTER     :: marked_states
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_localization_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, mystate, ns, &
                                                output_unit, stat
    INTEGER, DIMENSION(:), POINTER           :: lstates, marked_states_spin
    LOGICAL                                  :: do_homo, failure
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: scenter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv, matrix_s
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(section_vals_type), POINTER         :: loc_print_section
    TYPE(wannier_centres_type), &
      DIMENSION(:), POINTER                  :: wc

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(mos, ks_rmpv, dft_control, loc_print_section, marked_states_spin,&
            matrix_s, molecule_set, scenter,wc)
    CALL get_qs_env(qs_env,mos=mos,matrix_ks=ks_rmpv,molecule_set=molecule_set,&
         dft_control=dft_control,matrix_s=matrix_s,error=error)
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)
    loc_print_section => section_vals_get_subs_vals(loc_section,"PRINT",error=error)
    do_homo=qs_loc_env%localized_wfn_control%do_homo
    IF (BTEST(cp_print_key_should_output(logger%iter_info,loc_print_section,&
               "WANNIER_STATES",error=error),cp_p_file)) THEN
        CALL get_qs_env(qs_env=qs_env,WannierCentres=wc,error=error)
        IF (.NOT. ASSOCIATED(wc)) THEN
          ALLOCATE(wc(dft_control%nspins),stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
        ENDIF
    ENDIF
    DO ispin=1,dft_control%nspins
       !
       IF(do_homo)THEN
          qs_loc_env%tag_mo="HOMO"
       ELSE
          qs_loc_env%tag_mo="LUMO"
       END IF

       IF (qs_loc_env%do_localize) THEN
          ! Do the Real localization..
          IF (output_unit>0.AND.do_homo) WRITE(output_unit,"(/,T2,A,I3)")&
               "LOCALIZATION| Computing localization properties "//&
               "for OCCUPIED ORBITALS. Spin:",ispin
          IF (output_unit>0.AND.(.NOT.do_homo)) WRITE(output_unit,"(/,T2,A,I3)")&
               "LOCALIZATION| Computing localization properties "//&
               "for UNOCCUPIED ORBITALS. Spin:",ispin

          scenter => qs_loc_env%localized_wfn_control%centers_set(ispin)%array

          CALL qs_loc_driver(qs_env,qs_loc_env,loc_section,loc_print_section,&
               myspin=ispin,ext_mo_coeff=mo_local(ispin)%matrix,error=error)

          ! maps wfc to molecules, and compute the molecular dipoles if required
          IF (( BTEST(cp_print_key_should_output(logger%iter_info,loc_print_section,&
               "MOLECULAR_DIPOLES",error=error),cp_p_file) .OR. &
               BTEST(cp_print_key_should_output(logger%iter_info,loc_print_section,&
               "MOLECULAR_STATES",error=error),cp_p_file))) THEN
             CALL wfc_to_molecule(qs_env, qs_loc_env, loc_print_section, scenter,&
                  molecule_set, dft_control%nspins, error)
          END IF

            ! Compute the wannier states
          IF (BTEST(cp_print_key_should_output(logger%iter_info,loc_print_section,&
               "WANNIER_STATES",error=error),cp_p_file)) THEN
             ns=SIZE(qs_loc_env%localized_wfn_control%loc_states,1)
             IF (.NOT. ASSOCIATED(wc(ispin)%centres)) THEN
               ALLOCATE(wc(ispin)%WannierHamDiag(ns),stat=stat)
               CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
               ALLOCATE(wc(ispin)%centres(3,ns),stat=stat)
               CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             ENDIF

             wc(ispin)%centres(:,:)=scenter(1+(ispin-1)*3:ispin*3,:)
             lstates => qs_loc_env%localized_wfn_control%loc_states(:,ispin)
             CALL construct_wannier_states(molecule_set, mo_local(ispin)%matrix,&
              ks_rmpv(ispin)%matrix, qs_env, loc_print_section=loc_print_section,&
              WannierCentres=wc(ispin),ns=ns,states=lstates, error=error)
          ENDIF
          ! Compute the molecular states
          IF ( BTEST(cp_print_key_should_output(logger%iter_info,loc_print_section,&
               "MOLECULAR_STATES",error=error),cp_p_file)) THEN
             CALL construct_molecular_states(molecule_set, mo_local(ispin)%matrix, coeff(ispin)%matrix,&
                  evals(ispin)%array,ks_rmpv(ispin)%matrix, matrix_s(1)%matrix, qs_env, wf_r, wf_g, &
                  loc_print_section=loc_print_section, particles=particles, tag=TRIM(qs_loc_env%tag_mo),&
                  marked_states=marked_states_spin,error=error)
             IF(ASSOCIATED(marked_states_spin))THEN
                IF(.NOT.ASSOCIATED(marked_states))THEN
                   ALLOCATE(marked_states(SIZE(marked_states_spin),dft_control%nspins,2),stat=stat)
                   CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
                END IF
                mystate=1
                IF(qs_loc_env%tag_mo=="LUMO")mystate=2
                marked_states(:,ispin,mystate)=marked_states_spin(:)
                DEALLOCATE(marked_states_spin,stat=stat)
                CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          ENDIF
       END IF
    ENDDO
    IF (BTEST(cp_print_key_should_output(logger%iter_info,loc_print_section,&
          "WANNIER_STATES",error=error),cp_p_file)) THEN
      CALL set_qs_env(qs_env=qs_env,WannierCentres=wc,error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE get_localization_info

! *****************************************************************************
!> \brief Gets the lumos, and eigenvalues for the lumos
!> \param qs_env ...
!> \param scf_env ...
!> \param unoccupied_orbs ...
!> \param unoccupied_evals ...
!> \param nlumo ...
!> \param nlumos ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE make_lumo(qs_env,scf_env,unoccupied_orbs,unoccupied_evals,nlumo,nlumos,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_orbs
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_evals
    INTEGER                                  :: nlumo
    INTEGER, INTENT(OUT)                     :: nlumos
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_lumo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homo, ispin, istat, &
                                                n, nao, nmo, output_unit
    LOGICAL                                  :: failure
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv, matrix_s
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(preconditioner_type), POINTER       :: local_preconditioner
    TYPE(scf_control_type), POINTER          :: scf_control

    CALL timeset(routineN,handle)

    NULLIFY(mos,ks_rmpv,scf_control,dft_control,admm_env,para_env,blacs_env)
    CALL get_qs_env(qs_env,&
                    mos=mos,&
                    matrix_ks=ks_rmpv,&
                    scf_control=scf_control,&
                    dft_control=dft_control,&
                    matrix_s=matrix_s,&
                    admm_env=admm_env,&
                    para_env=para_env,&
                    blacs_env=blacs_env,&
                    error=error)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    DO ispin=1,dft_control%nspins
       NULLIFY(unoccupied_orbs(ispin)%matrix)
       NULLIFY(unoccupied_evals(ispin)%array)
       ! Always write eigenvalues
       IF (output_unit>0) WRITE(output_unit,*) " "
       IF (output_unit>0) WRITE(output_unit,*) " Lowest Eigenvalues of the unoccupied subspace spin ",ispin
       IF (output_unit>0) WRITE(output_unit,FMT='(1X,A)') "-----------------------------------------------------"
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,homo=homo,nao=nao,nmo=nmo)
       CALL cp_fm_get_info(mo_coeff, nrow_global=n,error=error)
       nlumos=MAX(1,MIN(nlumo,nao-nmo))
       IF (nlumo==-1) nlumos=nao-nmo
       ALLOCATE(unoccupied_evals(ispin)%array(nlumos),stat=istat)
       CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
       CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=blacs_env, &
            nrow_global=n,ncol_global=nlumos,error=error)
       CALL cp_fm_create(unoccupied_orbs(ispin)%matrix, fm_struct_tmp,name="lumos",error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)
       CALL cp_fm_init_random(unoccupied_orbs(ispin)%matrix,nlumos,error=error)

       ! the full_all preconditioner makes not much sense for lumos search
       NULLIFY(local_preconditioner)
       IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
          local_preconditioner=>scf_env%ot_preconditioner(1)%preconditioner
          ! this one can for sure not be right (as it has to match a given C0)
          IF (local_preconditioner%in_use == ot_precond_full_all) THEN
             NULLIFY(local_preconditioner)
          ENDIF
       ENDIF

       ! ** If we do ADMM, we add have to modify the kohn-sham matrix
       IF( dft_control%do_admm ) THEN
         CALL admm_correct_for_eigenvalues(ispin, admm_env, ks_rmpv(ispin)%matrix, &
                                           error)
       END IF

       CALL ot_eigensolver(matrix_h=ks_rmpv(ispin)%matrix,matrix_s=matrix_s(1)%matrix, &
            matrix_c_fm=unoccupied_orbs(ispin)%matrix, &
            matrix_orthogonal_space_fm=mo_coeff, &
            eps_gradient=scf_control%eps_lumos, &
            preconditioner=local_preconditioner, &
            iter_max=scf_control%max_iter_lumos,&
            size_ortho_space=nmo,error=error)

       CALL calculate_subspace_eigenvalues(unoccupied_orbs(ispin)%matrix,ks_rmpv(ispin)%matrix,&
            unoccupied_evals(ispin)%array, scr=output_unit, &
            ionode=output_unit>0,error=error)


       ! ** If we do ADMM, we restore the original kohn-sham matrix
       IF( dft_control%do_admm ) THEN
         CALL admm_uncorrect_for_eigenvalues(ispin, admm_env, ks_rmpv(ispin)%matrix, &
                                             error)
       END IF

    END DO

    CALL timestop(handle)

  END SUBROUTINE make_lumo
! *****************************************************************************
!> \brief Computes and Prints Atomic Charges with several methods
!> \param input ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param rho ...
!> \param matrix_s ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_charges(input, logger, qs_env, rho, matrix_s, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_charges', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, print_level, unit_nr
    LOGICAL                                  :: print_it
    TYPE(section_vals_type), POINTER         :: density_fit_section, print_key

    CALL timeset(routineN,handle)

    ! Mulliken charges require no further computation and are printed from write_mo_free_results

    ! Compute the Lowdin charges
    print_key => section_vals_get_subs_vals(input,"DFT%PRINT%LOWDIN", error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
       unit_nr = cp_print_key_unit_nr(logger,input,"DFT%PRINT%LOWDIN",extension=".lowdin",&
                                      log_filename=.FALSE.,error=error)
       print_level = 1
       CALL section_vals_val_get(print_key,"PRINT_GOP",l_val=print_it,error=error)
       IF (print_it) print_level = 2
       CALL section_vals_val_get(print_key,"PRINT_ALL",l_val=print_it,error=error)
       IF (print_it) print_level = 3
       CALL lowdin_population_analysis(qs_env,unit_nr,print_level,error)
       CALL cp_print_key_finished_output(unit_nr, logger,input,"DFT%PRINT%LOWDIN", error=error)
    END IF

    ! Compute the RESP charges
    CALL resp_fit(qs_env,error)

    ! Compute the Density Derived Atomic Point charges with the Bloechl scheme
    print_key => section_vals_get_subs_vals(input,"PROPERTIES%FIT_CHARGE", error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
       unit_nr=cp_print_key_unit_nr(logger,input,"PROPERTIES%FIT_CHARGE",extension=".Fitcharge",&
            log_filename=.FALSE.,error=error)
       density_fit_section =>  section_vals_get_subs_vals(input,"DFT%DENSITY_FITTING", error=error)
       CALL get_ddapc(qs_env,.FALSE.,density_fit_section,iwc=unit_nr,error=error)
       CALL cp_print_key_finished_output(unit_nr, logger,input,"PROPERTIES%FIT_CHARGE", error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_post_charges

! *****************************************************************************
!> \brief Computes and prints the Dipole (using localized charges)
!> \param input ...
!> \param dft_control ...
!> \param qs_loc_env ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_loc_dip(input, dft_control, qs_loc_env, logger, qs_env, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_loc_dip', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: description, &
                                                descriptionThisDip, iter
    COMPLEX(KIND=dp)                         :: zeta
    COMPLEX(KIND=dp), DIMENSION(3)           :: ggamma, zphase
    INTEGER                                  :: handle, i, ikind, ispins, j, &
                                                n_rep, reference, unit_nr
    LOGICAL                                  :: do_berry, failure, &
                                                first_time, ghost
    REAL(KIND=dp)                            :: charge_tot, theta, zeff, zwfc
    REAL(KIND=dp), DIMENSION(3)              :: ci, dipole, dipole_old, gvec, &
                                                rcc, ria
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ref_point
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_result_type), POINTER            :: results
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    failure = .FALSE.
    print_key => section_vals_get_subs_vals(input,"DFT%LOCALIZE%PRINT%TOTAL_DIPOLE", error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,first_time=first_time,&
         error=error),cp_p_file)) THEN
       NULLIFY(cell, particle_set, qs_kind_set, ref_point, results)
       CALL get_qs_env(qs_env=qs_env,&
                       cell=cell,&
                       particle_set=particle_set,&
                       qs_kind_set=qs_kind_set,&
                       results=results,&
                       error=error)

       reference = section_get_ival(print_key,keyword_name="REFERENCE",error=error)
       CALL section_vals_val_get(print_key,"REF_POINT",r_vals=ref_point,error=error)
       CALL section_vals_val_get(print_key,"PERIODIC",l_val=do_berry,error=error)
       description='[DIPOLE]'
       descriptionThisDip='[TOTAL_DIPOLE]'
       CALL get_reference_point(rcc,qs_env=qs_env,reference=reference,ref_point=ref_point,error=error)

       dipole=0.0_dp
       IF (do_berry) THEN
          rcc = pbc(rcc,cell)
          charge_tot = REAL(dft_control%charge,KIND=dp)
          ria  = twopi * MATMUL(cell%h_inv, rcc)
          zphase  = CMPLX(COS(ria), SIN(ria), KIND=dp)**charge_tot
          ggamma  = CMPLX(1.0_dp,0.0_dp,KIND=dp)

          ! Nuclear charges
          DO i=1,SIZE(particle_set)
             CALL get_atomic_kind(particle_set(i)%atomic_kind,kind_number=ikind)
             CALL get_qs_kind(qs_kind_set(ikind), ghost=ghost)
             IF (.NOT.ghost) THEN
                CALL get_qs_kind(qs_kind_set(ikind), core_charge=zeff)
                ria = pbc(particle_set(i)%r, cell)
                DO j = 1, 3
                   gvec  = twopi*cell%h_inv(j,:)
                   theta = SUM(ria(:)*gvec(:))
                   zeta  = CMPLX(COS(theta),SIN(theta),KIND=dp)**(zeff)
                   ggamma(j)  = ggamma(j) * zeta
                END DO
             END IF
          END DO

          ! Charges of the wfc involved
          ! Warning, this assumes the same occupation for all states
          zwfc = 3.0_dp-REAL(dft_control%nspins,dp)

          DO ispins=1,dft_control%nspins
             DO i=1,SIZE(qs_loc_env%localized_wfn_control%centers_set(ispins)%array,2)
                ria = pbc(qs_loc_env%localized_wfn_control%centers_set(ispins)%array(1:3,i),cell)
                DO j = 1, 3
                   gvec  = twopi*cell%h_inv(j,:)
                   theta = SUM(ria(:)*gvec(:))
                   zeta  = CMPLX(COS(theta),SIN(theta),KIND=dp)**(-zwfc)
                   ggamma(j)  = ggamma(j) * zeta
                END DO
             ENDDO
          ENDDO
          ggamma = ggamma * zphase
          ci     = AIMAG(LOG(ggamma))/twopi
          dipole = MATMUL ( cell%hmat, ci )
       ELSE
          ! Charges of the atoms involved
          DO i=1,SIZE(particle_set)
             CALL get_atomic_kind(particle_set(i)%atomic_kind, kind_number=ikind)
             CALL get_qs_kind(qs_kind_set(ikind), ghost=ghost)
             IF (.NOT.ghost) THEN
                CALL get_qs_kind(qs_kind_set(ikind), core_charge=zeff)
                ria = pbc(particle_set(i)%r, cell)
                dipole=dipole + zeff*(ria-rcc)
             END IF
          END DO

          ! Charges of the wfc involved
          ! Warning, this assumes the same occupation for all states
          zwfc = 3.0_dp-REAL(dft_control%nspins,dp)

          DO ispins=1,dft_control%nspins
             DO i=1,SIZE(qs_loc_env%localized_wfn_control%centers_set(ispins)%array,2)
                ria = pbc(qs_loc_env%localized_wfn_control%centers_set(ispins)%array(1:3,i),cell)
                dipole=dipole - zwfc * (ria-rcc)
             ENDDO
          ENDDO
       END IF

       ! Print and possibly store results
       unit_nr=cp_print_key_unit_nr(logger,print_key,extension=".Dipole",&
            middle_name="TOTAL_DIPOLE",error=error)
       IF (unit_nr>0) THEN
          IF (first_time) THEN
             WRITE(unit=unit_nr,fmt="(A,T31,A,T88,A,T136,A)")&
                  "# iter_level","dipole(x,y,z)[atomic units]",&
                                 "dipole(x,y,z)[debye]",&
                                 "delta_dipole(x,y,z)[atomic units]"
          END IF
          iter=cp_iter_string(logger%iter_info,error=error)
          CALL get_results(results,descriptionThisDip,n_rep=n_rep,error=error)
          IF(n_rep==0)THEN
             dipole_old=0._dp
          ELSE
             CALL get_results(results,descriptionThisDip,dipole_old,nval=n_rep,error=error)
          END IF
          IF (do_berry) THEN
             WRITE(unit=unit_nr,fmt="(a,9(es18.8))")&
                  iter(1:15), dipole, dipole*debye, pbc(dipole-dipole_old,cell)
          ELSE
             WRITE(unit=unit_nr,fmt="(a,9(es18.8))")&
                  iter(1:15), dipole, dipole*debye, (dipole-dipole_old)
          END IF
       END IF
       CALL cp_print_key_finished_output(unit_nr,logger,print_key,error=error)
       CALL cp_results_erase(results,description,error=error)
       CALL put_results(results,description,dipole,error)
       CALL cp_results_erase(results,descriptionThisDip,error=error)
       CALL put_results(results,descriptionThisDip,dipole,error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_post_loc_dip

! *****************************************************************************
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
!> \param basis_set_id ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_occ_cubes(input, dft_section, dft_control, logger, qs_env,&
       mo_coeff, wf_g, wf_r, particles, homo, ispin, basis_set_id, error)
    TYPE(section_vals_type), POINTER         :: input, dft_section
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(pw_p_type)                          :: wf_g, wf_r
    TYPE(particle_list_type), POINTER        :: particles
    INTEGER, INTENT(IN)                      :: homo, ispin
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_occ_cubes', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename, my_pos_cube, title
    INTEGER                                  :: handle, ivector, &
                                                my_basis_set_id, nhomo, &
                                                unit_nr
    LOGICAL                                  :: append_cube
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL timeset(routineN,handle)

    IF(PRESENT(basis_set_id)) THEN
      my_basis_set_id = basis_set_id
    ELSE
      my_basis_set_id = use_orb_basis_set
    END IF

    IF (BTEST(cp_print_key_should_output(logger%iter_info,dft_section,"PRINT%MO_CUBES",&
         error=error),cp_p_file) .AND. section_get_lval(dft_section,"PRINT%MO_CUBES%WRITE_CUBE",error=error)) THEN
       nhomo=section_get_ival(dft_section,"PRINT%MO_CUBES%NHOMO",error=error)
       append_cube = section_get_lval(dft_section,"PRINT%MO_CUBES%APPEND",error=error)
       my_pos_cube="REWIND"
       IF(append_cube) THEN
          my_pos_cube="APPEND"
       END IF
       IF (nhomo==-1) nhomo=homo
       DO ivector=MAX(1,homo-nhomo+1),homo
          CALL get_qs_env(qs_env=qs_env,&
                          atomic_kind_set=atomic_kind_set,&
                          qs_kind_set=qs_kind_set,&
                          cell=cell,&
                          particle_set=particle_set,&
                          pw_env=pw_env,&
                          error=error)
          CALL calculate_wavefunction(mo_coeff,ivector,wf_r,wf_g,atomic_kind_set,qs_kind_set,&
               cell,dft_control,particle_set,pw_env ,basis_set_id=my_basis_set_id, error=error)
          WRITE(filename,'(a4,I5.5,a1,I1.1)')"WFN_",ivector,"_",ispin
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%MO_CUBES",extension=".cube",&
               middle_name=TRIM(filename),file_position=my_pos_cube,log_filename=.FALSE.,error=error)
          WRITE(title,*) "WAVEFUNCTION ",ivector," spin ",ispin," i.e. HOMO - ",ivector-homo
          CALL cp_pw_to_cube(wf_r%pw,unit_nr,title,particles=particles,&
               stride=section_get_ivals(dft_section,"PRINT%MO_CUBES%STRIDE",error=error),&
               error=error)
          CALL cp_print_key_finished_output(unit_nr,logger,input,"DFT%PRINT%MO_CUBES",error=error)
       ENDDO
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_post_occ_cubes

! *****************************************************************************
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
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_unocc_cubes(input, dft_section, dft_control, logger, qs_env,&
       unoccupied_orbs, wf_g, wf_r, particles, nlumos, homo, ispin, lumo, error)
    TYPE(section_vals_type), POINTER         :: input, dft_section
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_type), POINTER                :: unoccupied_orbs
    TYPE(pw_p_type)                          :: wf_g, wf_r
    TYPE(particle_list_type), POINTER        :: particles
    INTEGER, INTENT(IN)                      :: nlumos, homo, ispin
    INTEGER, INTENT(IN), OPTIONAL            :: lumo
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_unocc_cubes', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename, my_pos_cube, title
    INTEGER                                  :: handle, ifirst, index_mo, &
                                                ivector, unit_nr
    LOGICAL                                  :: append_cube
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set

    CALL timeset(routineN,handle)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,dft_section,"PRINT%MO_CUBES",error=error),cp_p_file)&
         .AND.  section_get_lval(dft_section,"PRINT%MO_CUBES%WRITE_CUBE",error=error) ) THEN
       NULLIFY(qs_kind_set, particle_set, pw_env, cell)
       append_cube = section_get_lval(dft_section,"PRINT%MO_CUBES%APPEND",error=error)
       my_pos_cube="REWIND"
       IF(append_cube) THEN
          my_pos_cube="APPEND"
       END IF
       ifirst = 1
       IF(PRESENT(lumo)) ifirst = lumo
       DO ivector=ifirst, ifirst+nlumos-1
          CALL get_qs_env(qs_env=qs_env,&
                          atomic_kind_set=atomic_kind_set,&
                          qs_kind_set=qs_kind_set,&
                          cell=cell,&
                          particle_set=particle_set,&
                          pw_env=pw_env,&
                          error=error)
          CALL calculate_wavefunction(unoccupied_orbs, ivector, wf_r, wf_g, atomic_kind_set, &
                           qs_kind_set, cell, dft_control, particle_set, pw_env, error=error)

          IF(ifirst==1) THEN
             index_mo =  homo+ivector
          ELSE
             index_mo =  ivector
          END IF
          WRITE(filename,'(a4,I5.5,a1,I1.1)')"WFN_",index_mo,"_",ispin
             
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%MO_CUBES",extension=".cube",&
               middle_name=TRIM(filename),file_position=my_pos_cube,log_filename=.FALSE.,error=error)
          WRITE(title,*) "WAVEFUNCTION ",index_mo," spin ",ispin," i.e. LUMO + ",ifirst+ivector-2
          CALL cp_pw_to_cube(wf_r%pw, unit_nr, title, particles=particles,&
               stride=section_get_ivals(dft_section,"PRINT%MO_CUBES%STRIDE",error=error),&
               error=error)
          CALL cp_print_key_finished_output(unit_nr,logger,input,"DFT%PRINT%MO_CUBES",error=error)
       ENDDO
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_post_unocc_cubes

! *****************************************************************************
!> \brief Computes and prints the optical conductivity
!> \param input ...
!> \param dft_section ...
!> \param dft_control ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param unoccupied_orbs ...
!> \param unoccupied_evals ...
!> \param output_unit ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_optc(input, dft_section, dft_control, logger, qs_env,&
       unoccupied_orbs, unoccupied_evals, output_unit, error)
    TYPE(section_vals_type), POINTER         :: input, dft_section
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_orbs
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_evals
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_optc', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: handle, homo, ispin, nmo, &
                                                unit_nr
    LOGICAL                                  :: homoEQnmo
    REAL(KIND=dp)                            :: volume
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos

    CALL timeset(routineN,handle)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
         "DFT%PRINT%OPTICAL_CONDUCTIVITY",error=error),cp_p_file)) THEN
       NULLIFY(mos, cell, mo_coeff, mo_eigenvalues)
       CALL get_qs_env(qs_env, mos=mos, cell=cell, matrix_s=matrix_s,error=error)
       IF (output_unit>0) WRITE(output_unit,*) " "
       IF (output_unit>0) WRITE(output_unit,*) " Computing the optical conductivity "
       IF (output_unit>0) WRITE(output_unit,*) " Experimental version "
       IF (output_unit>0) WRITE(output_unit,*) " Check the code before believing results "
       homoEQnmo=.TRUE.
       DO ispin=1,dft_control%nspins
          CALL get_mo_set(mo_set=mos(ispin)%mo_set,homo=homo,nmo=nmo)
          IF (homo.NE.nmo) homoEQnmo=.FALSE.
       ENDDO
       IF (.NOT. homoEQnmo) THEN
          IF (output_unit>0) WRITE(output_unit,*) " homo.NE.nmo : skip optical conductivty "
       ELSE
          IF (.NOT.( BTEST(cp_print_key_should_output(logger%iter_info,dft_section,"PRINT%MO_CUBES",&
               error=error),cp_p_file) )) &
               CALL stop_program(routineN,moduleN,__LINE__,"Needs MO_CUBES to be activated")
          filename="CONDUCTIVITY"
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%OPTICAL_CONDUCTIVITY",&
               extension=".data",middle_name=TRIM(filename),log_filename=.FALSE.,error=error)
          CALL get_cell(cell,deth=volume)

          DO ispin=1,dft_control%nspins
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                             eigenvalues=mo_eigenvalues)
             CALL optical_conductivity(matrix_s, mo_coeff, mo_eigenvalues, unoccupied_orbs(ispin)%matrix,&
                                       unoccupied_evals(ispin)%array, volume, unit_nr, error=error)
          ENDDO
          CALL cp_print_key_finished_output(unit_nr,logger,input,"DFT%PRINT%OPTICAL_CONDUCTIVITY",&
               error=error)
          IF (output_unit>0) WRITE(output_unit,*) " Results written to file ",filename
       ENDIF
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_post_optc

! *****************************************************************************
!> \brief Computes and prints electric moments
!> \param input ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param output_unit ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_moments(input, logger, qs_env, output_unit, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_moments', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: handle, maxmom, reference, &
                                                unit_nr
    LOGICAL                                  :: magnetic, periodic
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ref_point
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    print_key => section_vals_get_subs_vals(section_vals=input,&
         subsection_name="DFT%PRINT%MOMENTS",error=error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN

       maxmom = section_get_ival(section_vals=input,&
             keyword_name="DFT%PRINT%MOMENTS%MAX_MOMENT",error=error)
       periodic = section_get_lval(section_vals=input,&
             keyword_name="DFT%PRINT%MOMENTS%PERIODIC",error=error)
       reference = section_get_ival(section_vals=input,&
             keyword_name="DFT%PRINT%MOMENTS%REFERENCE",error=error)
       magnetic = section_get_lval(section_vals=input,&
             keyword_name="DFT%PRINT%MOMENTS%MAGNETIC",error=error)
       NULLIFY ( ref_point )
       CALL section_vals_val_get(input,"DFT%PRINT%MOMENTS%REF_POINT",r_vals=ref_point,error=error)
       unit_nr = cp_print_key_unit_nr(logger=logger,basis_section=input,&
             print_key_path="DFT%PRINT%MOMENTS",extension=".dat",&
             middle_name="moments",log_filename=.FALSE.,error=error)

       IF (output_unit>0) THEN
          IF(unit_nr /= output_unit) THEN
            INQUIRE (UNIT=unit_nr,NAME=filename)
            WRITE (UNIT=output_unit,FMT="(/,T2,A,2(/,T3,A),/)")&
                 "MOMENTS","The electric/magnetic moments are written to file:",&
                 TRIM(filename)
          ELSE
            WRITE (UNIT=output_unit,FMT="(/,T2,A)") "ELECTRIC/MAGNETIC MOMENTS"
          END IF
       END IF

       IF (periodic) THEN
          CALL qs_moment_berry_phase(qs_env,magnetic,maxmom,reference,ref_point,unit_nr,error)
       ELSE
          CALL qs_moment_locop(qs_env,magnetic,maxmom,reference,ref_point,unit_nr,error)
       END IF

       CALL cp_print_key_finished_output(unit_nr=unit_nr,logger=logger,&
            basis_section=input,print_key_path="DFT%PRINT%MOMENTS",&
            error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_post_moments

! *****************************************************************************
!> \brief Computes and prints the X-ray diffraction spectrum.
!> \param input ...
!> \param dft_section ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param output_unit ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_xray(input,dft_section,logger,qs_env,output_unit,error)

    TYPE(section_vals_type), POINTER         :: input, dft_section
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_post_xray', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: handle, unit_nr
    REAL(KIND=dp)                            :: q_max
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    print_key => section_vals_get_subs_vals(section_vals=input,&
                                            subsection_name="DFT%PRINT%XRAY_DIFFRACTION_SPECTRUM",&
                                            error=error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
      q_max = section_get_rval(section_vals=dft_section,&
                               keyword_name="PRINT%XRAY_DIFFRACTION_SPECTRUM%Q_MAX",&
                               error=error)
      unit_nr = cp_print_key_unit_nr(logger=logger,&
                                     basis_section=input,&
                                     print_key_path="DFT%PRINT%XRAY_DIFFRACTION_SPECTRUM",&
                                     extension=".dat",&
                                     middle_name="xrd",&
                                     log_filename=.FALSE.,&
                                     error=error)
      IF (output_unit>0) THEN
        INQUIRE (UNIT=unit_nr,NAME=filename)
        WRITE (UNIT=output_unit,FMT="(/,/,T2,A)")&
          "X-RAY DIFFRACTION SPECTRUM"
        IF (unit_nr /= output_unit) THEN
          WRITE (UNIT=output_unit,FMT="(/,T3,A,/,/,T3,A,/)")&
            "The coherent X-ray diffraction spectrum is written to the file:",&
            TRIM(filename)
        END IF
      END IF
      CALL xray_diffraction_spectrum(qs_env=qs_env,&
                                     unit_number=unit_nr,&
                                     q_max=q_max,&
                                     error=error)
      CALL cp_print_key_finished_output(unit_nr=unit_nr,&
                                        logger=logger,&
                                        basis_section=input,&
                                        print_key_path="DFT%PRINT%XRAY_DIFFRACTION_SPECTRUM",&
                                        error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_post_xray

! *****************************************************************************
!> \brief Computes and prints Electric Field Gradient
!> \param input ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_efg(input, logger, qs_env, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_efg', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    print_key => section_vals_get_subs_vals(section_vals=input,&
                   subsection_name="DFT%PRINT%ELECTRIC_FIELD_GRADIENT",&
                   error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),&
              cp_p_file)) THEN
       CALL qs_efg_calc(qs_env=qs_env,error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_post_efg

! *****************************************************************************
!> \brief Computes the Electron Transfer Coupling matrix element
!> \param input ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param dft_control ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_et(input, qs_env, dft_control, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_et', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, stat
    LOGICAL                                  :: do_et, failure
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: my_mos
    TYPE(section_vals_type), POINTER         :: et_section

    CALL timeset(routineN,handle)
    failure=.FALSE.

    do_et=.FALSE.
    et_section =>  section_vals_get_subs_vals(input,"PROPERTIES%ET_COUPLING",&
                                                  error=error)
    CALL section_vals_get(et_section,explicit=do_et,error=error)
    IF(do_et)THEN
       IF(qs_env%et_coupling%first_run)THEN
          NULLIFY(my_mos)
          ALLOCATE(my_mos(dft_control%nspins),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(qs_env%et_coupling%et_mo_coeff(dft_control%nspins),STAT=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
          DO ispin =1,dft_control%nspins
             NULLIFY(my_mos(ispin)%matrix)
             CALL cp_fm_create(matrix=my_mos(ispin)%matrix,&
                  matrix_struct=qs_env%mos(ispin)%mo_set%mo_coeff%matrix_struct,&
                  name="FIRST_RUN_COEFF"//TRIM(ADJUSTL(cp_to_string(ispin)))//"MATRIX",&
                  error=error)
             CALL cp_fm_to_fm(qs_env%mos(ispin)%mo_set%mo_coeff,&
                              my_mos(ispin)%matrix,error=error)
          END DO
          CALL set_et_coupling_type(qs_env%et_coupling,et_mo_coeff=my_mos,error=error)
          DEALLOCATE(my_mos)
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_post_et



! *****************************************************************************
!> \brief compute the electron localization function 
!>
!> \param input ...
!> \param logger ...
!> \param qs_env ...
!> \param error ...
!> \par History
!>      2012-07 Created [MI]
! *****************************************************************************

  SUBROUTINE qs_scf_post_elf(input, logger, qs_env, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_elf', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename, my_pos_cube, title
    INTEGER                                  :: handle, ispin, output_unit, &
                                                stat, unit_nr
    LOGICAL                                  :: append_cube, do_elf, failure, &
                                                gapw
    REAL(dp)                                 :: rho_cutoff
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: elf_r
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(section_vals_type), POINTER         :: elf_section

    CALL timeset(routineN,handle)
    failure=.FALSE.
    output_unit= cp_logger_get_default_io_unit(logger)

    do_elf=.FALSE.
    elf_section =>  section_vals_get_subs_vals(input,"DFT%PRINT%ELF_CUBE",error=error)
    CALL section_vals_get(elf_section,explicit=do_elf,error=error)
    IF(do_elf)THEN
      NULLIFY(dft_control, pw_env, auxbas_pw_pool, pw_pools, particles, subsys, elf_r)
      CALL get_qs_env(qs_env,dft_control=dft_control,  pw_env=pw_env, subsys=subsys, error=error)
      CALL qs_subsys_get(subsys,particles=particles,error=error)

      gapw=dft_control%qs_control%gapw
      IF(.NOT. gapw) THEN
        ! allocate
        ALLOCATE(elf_r(dft_control%nspins),stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                            pw_pools=pw_pools,error=error)
        DO ispin = 1,dft_control%nspins
          CALL pw_pool_create_pw(auxbas_pw_pool,elf_r(ispin)%pw,&
                                    use_data = REALDATA3D,&
                                    in_space = REALSPACE, error=error)
          CALL pw_zero(elf_r(ispin)%pw, error=error)
        END DO
   
        IF (output_unit>0) THEN
          WRITE (UNIT=output_unit,FMT="(/,T15,A,/)")&
                     " ----- ELF is computed on the real space grid -----"
        END IF
        rho_cutoff=section_get_rval(elf_section,"density_cutoff",error=error)
        CALL qs_elf_calc(qs_env, elf_r, rho_cutoff, error=error)
   
        ! write ELF into cube file
        append_cube = section_get_lval(elf_section,"APPEND",error=error)
        my_pos_cube="REWIND"
        IF(append_cube) THEN
          my_pos_cube="APPEND"
        END IF

        DO ispin = 1,dft_control%nspins
          WRITE(filename,'(a5,I1.1)')"ELF_S",ispin
          WRITE(title,*) "ELF spin ", ispin
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%ELF_CUBE",extension=".cube",&
                 middle_name=TRIM(filename),file_position=my_pos_cube,log_filename=.FALSE.,error=error)
          IF (output_unit>0) THEN
              INQUIRE (UNIT=unit_nr,NAME=filename)
               WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A)")&
                     "ELF is written in cube file format to the file:",&
                     TRIM(filename)
          END IF

          CALL cp_pw_to_cube(elf_r(ispin)%pw,unit_nr,title,particles=particles,&
                 stride=section_get_ivals(elf_section,"STRIDE",error=error),&
                 error=error)
          CALL cp_print_key_finished_output(unit_nr,logger,input,"DFT%PRINT%ELF_CUBE",error=error)

          CALL pw_pool_give_back_pw(auxbas_pw_pool,elf_r(ispin)%pw, error=error)
        END DO

        ! deallocate
        DEALLOCATE(elf_r,stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      ELSE
         ! not implemented
         CALL cp_unimplemented_error(fromWhere=routineP, &
            message="ELF not implemented for GAPW calculations!!", &
            error=error, error_level=cp_warning_level)

      END IF

    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_post_elf


! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param elf_r ...
!> \param rho_cutoff ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE qs_elf_calc(qs_env, elf_r, rho_cutoff, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: elf_r
    REAL(kind=dp), INTENT(IN)                :: rho_cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_elf_calc', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER :: ELFCUT = 0.0001_dp, ELFEPS = 2.87E-05_dp, &
      f18 = (1.0_dp/8.0_dp), f23 = (2.0_dp/3.0_dp), f53 = (5.0_dp/3.0_dp)

    INTEGER                                  :: handle, i, idir, ispin, j, k, &
                                                nspin, stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    INTEGER, DIMENSION(3, 3)                 :: nd
    LOGICAL                                  :: deriv_pw, failure
    REAL(kind=dp)                            :: cfermi, dum, elf_kernel, &
                                                norm_drho, rho_53, udvol
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: drho_g, drho_r, rho_r, tau_g, &
                                                tau_r
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_type), POINTER                   :: tmp_g
    TYPE(qs_rho_type), POINTER               :: rho

! for spin restricted systems

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(rho,rho_r,  drho_g, drho_r, tau_r, tau_g, pw_env, auxbas_pw_pool, pw_pools,tmp_g )

    CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, rho=rho, error=error)

    CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                            pw_pools=pw_pools,error=error)
    nspin = SIZE(rho%rho_r)
    bo = rho%rho_r(1)%pw%pw_grid%bounds_local
    cfermi = (3.0_dp/10.0_dp)*(pi*pi*3.0_dp)**f23

    ALLOCATE(rho_r(nspin),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tau_r(nspin),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(tau_g(nspin),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(drho_r(3*nspin),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(drho_g(3*nspin),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    DO ispin = 1,nspin
       rho_r(ispin)%pw => rho%rho_r(ispin)%pw
       IF (rho%tau_r_valid) THEN
         tau_r(ispin)%pw => rho%tau_r(ispin)%pw
       ELSE
         CALL pw_pool_create_pw(auxbas_pw_pool,tau_r(ispin)%pw,&
                     use_data=REALDATA3D,in_space=REALSPACE,error=error)
         CALL pw_pool_create_pw(auxbas_pw_pool,tau_g(ispin)%pw,&
                     use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
         CALL calculate_rho_elec(matrix_p=rho%rho_ao(ispin)%matrix,&
               rho=tau_r(ispin),&
               rho_gspace=tau_g(ispin),&
               total_rho=dum, &  
               qs_env=qs_env, soft_valid=.FALSE., &
               compute_tau=.TRUE., error=error)
       END IF

       IF(rho%drho_r_valid) THEN
         DO idir = 1,3
           drho_r(3*(ispin-1)+idir)%pw => rho%drho_r(3*(ispin-1)+idir)%pw
         END DO
       ELSE
         deriv_pw = .FALSE.
 !        deriv_pw = .TRUE.
         IF(deriv_pw) THEN
           nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))
           CALL pw_pool_create_pw(auxbas_pw_pool,tmp_g,&
                      use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
           udvol = 1.0_dp/rho%rho_r(ispin)%pw%pw_grid%dvol
           DO idir = 1,3
             CALL pw_pool_create_pw(auxbas_pw_pool,drho_r(3*(ispin-1)+idir)%pw,&
                      use_data=REALDATA3D,in_space=REALSPACE,error=error)
             CALL pw_transfer ( rho%rho_r(ispin)%pw, tmp_g , error=error)
             CALL pw_derive ( tmp_g, nd(:,idir) , error=error)
             CALL pw_transfer (tmp_g, drho_r(3*(ispin-1)+idir)%pw , error=error)
!              CALL pw_scale(drho_r(3*(ispin-1)+idir)%pw,udvol,error=error)

           END DO
            
         ELSE
           DO idir = 1,3
             CALL pw_pool_create_pw(auxbas_pw_pool,drho_r(3*(ispin-1)+idir)%pw,&
                      use_data=REALDATA3D,in_space=REALSPACE,error=error)
             CALL pw_pool_create_pw(auxbas_pw_pool,drho_g(3*(ispin-1)+idir)%pw,&
                       use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,error=error)
             CALL calculate_rho_elec(matrix_p=rho%rho_ao(ispin)%matrix,&
                  rho=drho_r(3*(ispin-1)+idir),&
                  rho_gspace=drho_g(3*(ispin-1)+idir),&
                  total_rho=dum, &  
                  qs_env=qs_env, soft_valid=.FALSE., &
                  compute_tau=.FALSE., compute_grad=.TRUE., idir =idir, error=error)
   
           END DO
         END IF
       END IF

       ! Calculate elf_r
       !$omp parallel do default(none) shared(bo,elf_r, ispin, drho_r,rho_r, tau_r, cfermi, rho_cutoff)&
       !$omp             private(k,j,i, norm_drho, rho_53, elf_kernel)
       DO k = bo(1,3), bo(2,3)
          DO j = bo(1,2), bo(2,2)
             DO i = bo(1,1), bo(2,1)
               norm_drho = drho_r(3*(ispin-1)+1)%pw%cr3d(i,j,k)**2+&
                        drho_r(3*(ispin-1)+2)%pw%cr3d(i,j,k)**2+&
                        drho_r(3*(ispin-1)+3)%pw%cr3d(i,j,k)**2
               norm_drho = norm_drho/MAX( rho_r(ispin)%pw%cr3d(i,j,k),rho_cutoff)
               rho_53 = cfermi * MAX(rho_r(ispin)%pw%cr3d(i,j,k),rho_cutoff)**f53
               elf_kernel = (tau_r(ispin)%pw%cr3d(i,j,k) - f18*norm_drho) +2.87E-5_dp
               elf_kernel = (elf_kernel/rho_53)**2
               elf_r(ispin)%pw%cr3d(i,j,k) =  1.0_dp/(1.0_dp+elf_kernel)
               IF(elf_r(ispin)%pw%cr3d(i,j,k)<ELFCUT ) elf_r(ispin)%pw%cr3d(i,j,k) =0.0_dp
             END DO
          END DO
       END DO

      IF (.NOT. rho%tau_r_valid) THEN
        CALL pw_pool_give_back_pw(auxbas_pw_pool,tau_r(ispin)%pw, error=error)
        CALL pw_pool_give_back_pw(auxbas_pw_pool,tau_g(ispin)%pw, error=error)
      END IF
      IF (.NOT. rho%drho_r_valid) THEN
        IF(deriv_pw) THEN
          CALL pw_pool_give_back_pw(auxbas_pw_pool,tmp_g, error=error)
          DO idir = 1,3
            CALL pw_pool_give_back_pw(auxbas_pw_pool,drho_r(3*(ispin-1)+idir)%pw, error=error)
          END DO
        ELSE
          DO idir = 1,3
            CALL pw_pool_give_back_pw(auxbas_pw_pool,drho_r(3*(ispin-1)+idir)%pw, error=error)
            CALL pw_pool_give_back_pw(auxbas_pw_pool,drho_g(3*(ispin-1)+idir)%pw, error=error)
          END DO
        END IF
      END IF
    END DO  !ispin

    DEALLOCATE(rho_r,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tau_r,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tau_g,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(drho_r,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(drho_g,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE  qs_elf_calc

! *****************************************************************************
!> \brief computes the condition number of the overlap matrix and
!>      prints the value of the total energy. This is needed
!>      for BASIS_MOLOPT optimizations
!> \param input ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      2007-07 Created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE qs_scf_post_molopt(input, logger, qs_env, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_molopt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao, unit_nr
    REAL(KIND=dp)                            :: S_cond_number
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_fm_struct_type), POINTER         :: ao_ao_fmstruct
    TYPE(cp_fm_type), POINTER                :: fm_s, fm_work, mo_coeff
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    print_key => section_vals_get_subs_vals(section_vals=input,&
                   subsection_name="DFT%PRINT%BASIS_MOLOPT_QUANTITIES",&
                   error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),&
              cp_p_file)) THEN

       CALL get_qs_env(qs_env,energy=energy,matrix_s=matrix_s,mos=mos,error=error)

       ! set up the two needed full matrices, using mo_coeff as a template
       CALL get_mo_set(mo_set=mos(1)%mo_set,mo_coeff=mo_coeff,nao=nao)
       CALL cp_fm_struct_create(fmstruct=ao_ao_fmstruct,&
            nrow_global=nao, ncol_global=nao,&
            template_fmstruct=mo_coeff%matrix_struct, error=error)
       CALL cp_fm_create(fm_s, matrix_struct=ao_ao_fmstruct,&
                     name="fm_s", error=error)
       CALL cp_fm_create(fm_work, matrix_struct=ao_ao_fmstruct,&
                     name="fm_work", error=error)
       CALL cp_fm_struct_release(ao_ao_fmstruct,error=error)
       ALLOCATE(eigenvalues(nao))

       CALL copy_dbcsr_to_fm(matrix_s(1)%matrix,fm_s,error=error)
       CALL choose_eigv_solver(fm_s,fm_work,eigenvalues,error=error)

       CALL cp_fm_release(fm_s,error=error)
       CALL cp_fm_release(fm_work,error=error)

       S_cond_number=MAXVAL(ABS(eigenvalues))/MAX(MINVAL(ABS(eigenvalues)),EPSILON(0.0_dp))

       unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%BASIS_MOLOPT_QUANTITIES",&
          extension=".molopt",error=error)

       IF (unit_nr>0) THEN
           ! please keep this format fixed, needs to be grepable for molopt
           ! optimizations
           WRITE(unit_nr,'(T2,A28,2A25)') "","Tot. Ener.","S Cond. Numb."
           WRITE(unit_nr,'(T2,A28,2E25.17)') "BASIS_MOLOPT_QUANTITIES",energy%total,S_cond_number
       ENDIF

       CALL cp_print_key_finished_output(unit_nr,logger,input,&
            "DFT%PRINT%BASIS_MOLOPT_QUANTITIES",error=error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_post_molopt

! *****************************************************************************
!> \brief Dumps EPR
!> \param input ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_epr(input, logger, qs_env, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_epr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)

    print_key => section_vals_get_subs_vals(section_vals=input,&
                   subsection_name="DFT%PRINT%HYPERFINE_COUPLING_TENSOR",&
                   error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),&
              cp_p_file)) THEN
       CALL qs_epr_hyp_calc(qs_env=qs_env,error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_post_epr

! *****************************************************************************
!> \brief Geminal optimizer
!> \param input ...
!> \param logger ...
!> \param qs_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_gemopt(input, logger, qs_env, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_gemopt', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, unit_nr
    LOGICAL                                  :: explicit
    TYPE(section_vals_type), POINTER         :: print_key

    CALL timeset(routineN,handle)


    print_key => section_vals_get_subs_vals(section_vals=input,&
                   subsection_name="DFT%PRINT%OPTIMIZE_GEMINALS",&
                   error=error)
    unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%OPTIMIZE_GEMINALS",&
                 extension="",error=error)
    CALL section_vals_get(print_key,explicit=explicit, error=error)

    IF ( explicit) THEN
      IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),&
                cp_p_file)) THEN
         IF (unit_nr>0) THEN
            WRITE(unit_nr,"(T2,79('-'),/,T2,'-',T30,A,T80,'-',/,T2,79('-'))") " Optimize Geminals "
         ENDIF
         !
         CALL geminal_optimize(qs_env,print_key,unit_nr,error)
         !
         IF (unit_nr>0) WRITE(unit_nr,"(T2,79('-'))")
      END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_scf_post_gemopt

! *****************************************************************************
!> \brief Interface routine to trigger writing of results available from normal 
!>        SCF. Can write MO-dependent and MO free results (needed for call from
!>        the linear scaling code)
!> \param qs_env the qs_env in which the qs_env lives
!> \param scf_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************

  SUBROUTINE write_available_results(qs_env,scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), OPTIONAL, POINTER :: scf_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_available_results', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)
    
    ! those properties that require MOs (not suitable density matrix based methods)
    CALL write_mo_dependent_results(qs_env,scf_env,error)

    ! those that depend only on the density matrix, they should be linear scaling in their implementation
    CALL write_mo_free_results(qs_env,error)

    CALL timestop(handle)

  END SUBROUTINE write_available_results

! *****************************************************************************
!> \brief Write QS results available if MO's are present (if switched on through the print_keys)
!>        Writes only MO dependent results. Split is necessary as ls_scf does not
!>        provide MO's
!> \param qs_env the qs_env in which the qs_env lives
!> \param scf_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE write_mo_dependent_results(qs_env,scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), OPTIONAL, POINTER :: scf_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_mo_dependent_results', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, homo, ispin, nmo, &
                                                output_unit
    LOGICAL                                  :: all_equal, failure
    REAL(KIND=dp)                            :: maxocc, s_square, &
                                                s_square_ideal, &
                                                total_abs_spin_dens
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues, &
                                                occupation_numbers
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv, matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_deriv
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: wf_r
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: dft_section, input

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(cell, dft_control, pw_env, auxbas_pw_pool, pw_pools, mo_coeff, &
         mo_coeff_deriv, mo_eigenvalues, mos, atomic_kind_set, qs_kind_set,&
         particle_set, rho, ks_rmpv, matrix_s, scf_control, dft_section, &
         molecule_set, input, particles, subsys)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,&
                       dft_control=dft_control,&
                       molecule_set=molecule_set, &
                       mos=mos,&
                       atomic_kind_set=atomic_kind_set,&
                       particle_set=particle_set,&
                       qs_kind_set=qs_kind_set,&
                       rho=rho,matrix_ks=ks_rmpv,&
                       admm_env=admm_env,&
                       scf_control=scf_control,&
                       matrix_s=matrix_s,&
                       input=input,&
                       cell=cell,&
                       subsys=subsys,&
                       error=error)
       CALL qs_subsys_get(subsys,particles=particles,error=error)

       ! *** if the dft_section tells you to do so, write last wavefunction to screen
       dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
       IF(.NOT.qs_env%run_rtp)THEN
          IF (dft_control%nspins == 2) THEN
             CALL write_mo_set(mos(1)%mo_set,atomic_kind_set,qs_kind_set,particle_set,4,&
                  dft_section,spin="ALPHA",last=.TRUE.,error=error)
             CALL write_mo_set(mos(2)%mo_set,atomic_kind_set,qs_kind_set,particle_set,4,&
                  dft_section,spin="BETA",last=.TRUE.,error=error)
          ELSE
             CALL write_mo_set(mos(1)%mo_set,atomic_kind_set,qs_kind_set,particle_set,4,&
                  dft_section,last=.TRUE.,error=error)
          END IF

          ! *** at the end of scf print out the projected dos per kind
          IF (BTEST(cp_print_key_should_output(logger%iter_info,dft_section,"PRINT%PDOS",&
               error=error),cp_p_file) ) THEN
             DO ispin = 1,dft_control%nspins

                ! ** If we do ADMM, we add have to modify the kohn-sham matrix
                IF( dft_control%do_admm ) THEN
                  CALL admm_correct_for_eigenvalues(ispin, admm_env, ks_rmpv(ispin)%matrix, &
                                                     error)
                END IF
                IF(PRESENT(scf_env))THEN
                   IF (scf_env%method == ot_method_nr) THEN
                      CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                           eigenvalues=mo_eigenvalues)
                      IF (ASSOCIATED(qs_env%mo_derivs)) THEN
                         mo_coeff_deriv=>qs_env%mo_derivs(ispin)%matrix
                      ELSE
                         mo_coeff_deriv=>NULL()
                      ENDIF
   
                      CALL calculate_subspace_eigenvalues(mo_coeff,ks_rmpv(ispin)%matrix,mo_eigenvalues, &
                           do_rotation=.TRUE.,&
                           co_rotate_dbcsr=mo_coeff_deriv,error=error)
                      CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,error=error)
   
                   END IF
                END IF
                IF(dft_control%nspins==2) THEN
                   CALL calculate_projected_dos(mos(ispin)%mo_set,atomic_kind_set,&
                        qs_kind_set,particle_set,qs_env, dft_section,ispin=ispin,error=error)
                ELSE
                   CALL calculate_projected_dos(mos(ispin)%mo_set,atomic_kind_set,&
                        qs_kind_set,particle_set,qs_env,dft_section,error=error)
                END IF

                ! ** If we do ADMM, we add have to modify the kohn-sham matrix
                IF( dft_control%do_admm ) THEN
                  CALL admm_uncorrect_for_eigenvalues(ispin, admm_env, ks_rmpv(ispin)%matrix, &
                                                       error)
                END IF

             END DO
          ENDIF
       END IF
       
       !   *** Integrated absolute spin density and spin contamination ***
       IF (dft_control%nspins.eq.2) THEN
          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
               pw_pools=pw_pools,error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          CALL pw_copy(rho%rho_r(1)%pw,wf_r%pw, error=error)
          CALL pw_axpy(rho%rho_r(2)%pw,wf_r%pw,alpha=-1._dp, error=error)
          total_abs_spin_dens=pw_integrate_function(wf_r%pw,oprt="ABS", error=error)
          IF (output_unit > 0) WRITE(UNIT=output_unit,FMT='(/,(T3,A,T61,F20.10))')&
               "Integrated absolute spin density  : ",total_abs_spin_dens
          CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_r%pw, error=error)
          !
          ! XXX Fix Me XXX
          ! should be extended to the case where added MOs are present
          !
          all_equal = .TRUE.
          DO ispin=1,dft_control%nspins
            CALL get_mo_set(mo_set=mos(ispin)%mo_set,&
                            occupation_numbers=occupation_numbers,&
                            homo=homo,&
                            nmo=nmo,&
                            maxocc=maxocc)
            IF (nmo > 0) THEN
              all_equal = all_equal.AND.&
                          (ALL(occupation_numbers(1:homo) == maxocc).AND.&
                           ALL(occupation_numbers(homo+1:nmo) == 0.0_dp))
             END IF
          END DO
          IF (.NOT.all_equal) THEN
            IF (output_unit>0) WRITE(UNIT=output_unit,FMT="(T3,A)")&
              "WARNING: S**2 computation does not yet treat fractional occupied orbitals"
          ELSE
             CALL get_qs_env(qs_env=qs_env,&
                              energy=energy,error=error)
             CALL compute_s_square(mos=mos, matrix_s=matrix_s, s_square=s_square,&
                  s_square_ideal=s_square_ideal,error=error)
             IF (output_unit > 0) WRITE (UNIT=output_unit,FMT='(T3,A,T51,2F15.6)')&
                  "Ideal and single determinant S**2 : ",s_square_ideal,s_square
             energy%s_square=s_square
          END IF
       ENDIF

    END IF

    CALL timestop(handle)

  END SUBROUTINE write_mo_dependent_results

! *****************************************************************************
!> \brief Write QS results always available (if switched on through the print_keys)
!>        Can be called from ls_scf
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE write_mo_free_results(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_mo_free_results', &
      routineP = moduleN//':'//routineN
    CHARACTER(len=1), DIMENSION(3), &
      PARAMETER                              :: cdir = (/"x","y","z"/)

    CHARACTER(LEN=2)                         :: element_symbol
    CHARACTER(LEN=default_path_length)       :: filename, my_pos_cube
    CHARACTER(LEN=default_string_length)     :: name
    INTEGER :: handle, i, iat, id, ikind, iso, ispin, istat, iw, l, nd(3), &
      ngto, niso, nkind, np, nr, output_unit, print_level, unit_nr
    LOGICAL                                  :: append_cube, failure, &
                                                print_it, &
                                                print_total_density, &
                                                write_ks, write_xc, &
                                                xrd_interface
    REAL(KIND=dp)                            :: q_max, rho_hard, rho_soft, &
                                                rho_total, rho_total_rspace, &
                                                udvol, volume
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: bfun
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: aedens, ccdens, ppdens
    REAL(KIND=dp), DIMENSION(3)              :: dr
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv, matrix_s, matrix_vxc
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(grid_atom_type), POINTER            :: grid_atom
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: aux_g, aux_r, &
                                                rho_elec_gspace, &
                                                rho_elec_rspace, wf_r
    TYPE(pw_p_type), POINTER                 :: rho0_s_gs, rho_core
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(pw_type), POINTER                   :: v_hartree_rspace
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_kind_type), POINTER              :: qs_kind
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(rho_atom_type), DIMENSION(:), &
      POINTER                                :: rho_atom_set
    TYPE(rho_atom_type), POINTER             :: rho_atom
    TYPE(section_vals_type), POINTER         :: dft_section, input, &
                                                print_key, xc_section

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(cell, dft_control, pw_env, auxbas_pw_pool, pw_pools, &
         atomic_kind_set, qs_kind_set, particle_set, rho, ks_rmpv, matrix_s, &
         dft_section, xc_section, input, particles, subsys, matrix_vxc,v_hartree_rspace)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,&
                       dft_control=dft_control, &
                       atomic_kind_set=atomic_kind_set,&
                       qs_kind_set=qs_kind_set,&
                       particle_set=particle_set,&
                       rho=rho,matrix_ks=ks_rmpv,&
                       matrix_s=matrix_s,&
                       matrix_vxc=matrix_vxc,&
                       input=input,&
                       cell=cell,&
                       subsys=subsys,&
                       para_env=para_env,&
                       v_hartree_rspace=v_hartree_rspace,&
                       error=error)
       dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
       CALL qs_subsys_get(subsys,particles=particles,error=error)

       ! Print the total density (electronic + core charge)

       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%TOT_DENSITY_CUBE", error=error),cp_p_file)) THEN
          NULLIFY(rho_core,rho0_s_gs)
          append_cube = section_get_lval(input,"DFT%PRINT%TOT_DENSITY_CUBE%APPEND",error=error)
          my_pos_cube="REWIND"
          IF(append_cube) THEN
            my_pos_cube="APPEND"
          END IF

          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,rho_core=rho_core,&
               rho0_s_gs=rho0_s_gs,error=error)
!          CALL pw_env_get(pw_env, error=error)
!          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
               pw_pools=pw_pools,error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          IF (dft_control%qs_control%gapw) THEN
             CALL pw_transfer(rho0_s_gs%pw,wf_r%pw,error=error)
             IF(dft_control%qs_control%gapw_control%nopaw_as_gpw) THEN
                CALL pw_axpy(rho_core%pw,wf_r%pw, error=error)
             END IF
          ELSE
             CALL pw_transfer(rho_core%pw,wf_r%pw,error=error)
          END IF
          DO ispin=1,dft_control%nspins
             CALL pw_axpy(rho%rho_r(ispin)%pw,wf_r%pw, error=error)
          END DO
          filename = "TOTAL_DENSITY"
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%TOT_DENSITY_CUBE",&
               extension=".cube",middle_name=TRIM(filename),file_position=my_pos_cube,&
               log_filename=.FALSE.,error=error)
          CALL cp_pw_to_cube(wf_r%pw,unit_nr,"TOTAL DENSITY",&
               particles=particles,&
               stride=section_get_ivals(dft_section,"PRINT%TOT_DENSITY_CUBE%STRIDE",error=error),&
               error=error)
          CALL cp_print_key_finished_output(unit_nr,logger,input,&
               "DFT%PRINT%TOT_DENSITY_CUBE",error=error)
          CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_r%pw, error=error)
       END IF

       ! Write cube file with electron density
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
                 "DFT%PRINT%E_DENSITY_CUBE",error=error),cp_p_file)) THEN
          CALL section_vals_val_get(dft_section,&
                                    keyword_name="PRINT%E_DENSITY_CUBE%TOTAL_DENSITY",&
                                    l_val=print_total_density,&
                                    error=error)
          append_cube = section_get_lval(input,"DFT%PRINT%E_DENSITY_CUBE%APPEND",error=error)
          my_pos_cube="REWIND"
          IF(append_cube) THEN
            my_pos_cube="APPEND"
          END IF
          ! Write the info on core densities for the interface between cp2k and the XRD code
          ! together with the valence density they are used to compute the form factor (Fourier transform)
          xrd_interface = section_get_lval(input,"DFT%PRINT%E_DENSITY_CUBE%XRD_INTERFACE",error=error)
          IF(xrd_interface) THEN
             !cube file only contains soft density (GAPW)
             IF(dft_control%qs_control%gapw) print_total_density = .FALSE.

             filename = "ELECTRON_DENSITY"
             unit_nr = cp_print_key_unit_nr(logger,input,"DFT%PRINT%E_DENSITY_CUBE",&
                                            extension=".xrd",middle_name=TRIM(filename),&
                                            file_position=my_pos_cube,log_filename=.FALSE.,error=error)
             ngto = section_get_ival(input,"DFT%PRINT%E_DENSITY_CUBE%NGAUSS",error=error)
             IF (output_unit>0) THEN
                 INQUIRE (UNIT=unit_nr,NAME=filename)
                 WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A)")&
                     "The electron density (atomic part) is written to the file:",&
                     TRIM(filename)
             END IF

             xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
             nkind = SIZE(atomic_kind_set)
             IF (unit_nr>0) THEN
                WRITE(unit_nr,*) "Atomic (core) densities"
                WRITE(unit_nr,*) "Unit cell"
                WRITE(unit_nr,FMT="(3F20.12)") cell%hmat(1,1),cell%hmat(1,2),cell%hmat(1,3)
                WRITE(unit_nr,FMT="(3F20.12)") cell%hmat(2,1),cell%hmat(2,2),cell%hmat(2,3)
                WRITE(unit_nr,FMT="(3F20.12)") cell%hmat(3,1),cell%hmat(3,2),cell%hmat(3,3)
                WRITE(unit_nr,*) "Atomic types"
                WRITE(unit_nr,*) nkind
             END IF
             ! calculate atomic density and core density
             ALLOCATE(ppdens(ngto,2,nkind),aedens(ngto,2,nkind),ccdens(ngto,2,nkind),stat=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
             DO ikind=1,nkind
               atomic_kind => atomic_kind_set(ikind)
               qs_kind     => qs_kind_set(ikind)
               CALL get_atomic_kind(atomic_kind,name=name,element_symbol=element_symbol)
               CALL calculate_atomic_density(ppdens(:,:,ikind),atomic_kind,qs_kind, ngto,&
                                             iunit=output_unit, confine=.TRUE.,error=error)
               CALL calculate_atomic_density(aedens(:,:,ikind),atomic_kind,qs_kind, ngto,&
                           iunit=output_unit,allelectron=.TRUE.,confine=.TRUE.,error=error)
               ccdens(:,1,ikind) = aedens(:,1,ikind)
               ccdens(:,2,ikind) = 0._dp
               CALL project_function_a(ccdens(1:ngto,2,ikind),ccdens(1:ngto,1,ikind),&
                    ppdens(1:ngto,2,ikind),ppdens(1:ngto,1,ikind),0,error)
               ccdens(:,2,ikind) = aedens(:,2,ikind) - ccdens(:,2,ikind)
               IF (unit_nr>0) THEN
                  WRITE(unit_nr,FMT="(I6,A10,A20)") ikind,TRIM(element_symbol),TRIM(name)
                  WRITE(unit_nr,FMT="(I6)") ngto
                  WRITE(unit_nr,*) "   Total density"
                  WRITE(unit_nr,FMT="(2G24.12)") (aedens(i,1,ikind),aedens(i,2,ikind),i=1,ngto)
                  WRITE(unit_nr,*) "    Core density"
                  WRITE(unit_nr,FMT="(2G24.12)") (ccdens(i,1,ikind),ccdens(i,2,ikind),i=1,ngto)
               END IF
               NULLIFY(atomic_kind)
             END DO

             IF (dft_control%qs_control%gapw) THEN
                CALL get_qs_env(qs_env=qs_env,rho_atom_set=rho_atom_set,error=error)

                IF (unit_nr>0) THEN
                   WRITE(unit_nr,*) "Coordinates and GAPW density"
                END IF
                np = particles%n_els
                DO iat=1,np
                   CALL get_atomic_kind(particles%els(iat)%atomic_kind,kind_number=ikind)
                   CALL get_qs_kind(qs_kind_set(ikind), grid_atom=grid_atom)
                   rho_atom => rho_atom_set(iat)
                   IF(ASSOCIATED(rho_atom%rho_rad_h(1)%r_coef)) THEN
                      nr = SIZE(rho_atom%rho_rad_h(1)%r_coef,1)
                      niso = SIZE(rho_atom%rho_rad_h(1)%r_coef,2)
                   ELSE
                      nr = 0
                      niso = 0
                   ENDIF
                   CALL mp_sum(nr,para_env%group)
                   CALL mp_sum(niso,para_env%group)

                   ALLOCATE(bfun(nr,niso),stat=istat)
                   CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                   bfun = 0._dp
                   DO ispin = 1,dft_control%nspins
                      IF(ASSOCIATED(rho_atom%rho_rad_h(1)%r_coef)) THEN
                         bfun(:,:) = bfun + rho_atom%rho_rad_h(ispin)%r_coef - rho_atom%rho_rad_s(ispin)%r_coef
                      END IF
                   END DO
                   CALL mp_sum(bfun,para_env%group)
                   ccdens(:,1,ikind) = ppdens(:,1,ikind)
                   ccdens(:,2,ikind) = 0._dp
                   IF (unit_nr>0) THEN
                      WRITE(unit_nr,'(I10,I5,3f12.6)') iat,ikind,particles%els(iat)%r
                   END IF
                   DO iso=1,niso
                      l = indso(1,iso)
                      CALL project_function_b(ccdens(:,2,ikind),ccdens(:,1,ikind),bfun(:,iso),grid_atom,l,error)
                      IF (unit_nr>0) THEN
                         WRITE(unit_nr,FMT="(3I6)") iso,l,ngto
                         WRITE(unit_nr,FMT="(2G24.12)") (ccdens(i,1,ikind),ccdens(i,2,ikind),i=1,ngto)
                      END IF
                   END DO
                   DEALLOCATE(bfun,stat=istat)
                   CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)
                END DO
             ELSE
                IF (unit_nr>0) THEN
                   WRITE(unit_nr,*) "Coordinates"
                   np = particles%n_els
                   DO iat=1,np
                      CALL get_atomic_kind(particles%els(iat)%atomic_kind,kind_number=ikind)
                      WRITE(unit_nr,'(I10,I5,3f12.6)') iat,ikind,particles%els(iat)%r
                   END DO
                END IF
             END IF

             DEALLOCATE(ppdens,aedens,ccdens,stat=istat)
             CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

             CALL cp_print_key_finished_output(unit_nr,logger,input,&
                                               "DFT%PRINT%E_DENSITY_CUBE",&
                                               error=error)

          END IF
          IF (dft_control%qs_control%gapw.AND.print_total_density) THEN
             ! Print total electronic density
             CALL get_qs_env(qs_env=qs_env,&
                             pw_env=pw_env,&
                             error=error)
             CALL pw_env_get(pw_env=pw_env,&
                             auxbas_pw_pool=auxbas_pw_pool,&
                             pw_pools=pw_pools,&
                             error=error)
             CALL pw_pool_create_pw(pool=auxbas_pw_pool,&
                                    pw=rho_elec_rspace%pw,&
                                    use_data=REALDATA3D,&
                                    in_space=REALSPACE,&
                                    error=error)
             CALL pw_zero(rho_elec_rspace%pw,error=error)
             CALL pw_pool_create_pw(pool=auxbas_pw_pool,&
                                    pw=rho_elec_gspace%pw,&
                                    use_data=COMPLEXDATA1D,&
                                    in_space=RECIPROCALSPACE,&
                                    error=error)
             CALL pw_zero(rho_elec_gspace%pw,error=error)
             CALL get_pw_grid_info(pw_grid=rho_elec_gspace%pw%pw_grid,&
                                   dr=dr,&
                                   vol=volume,&
                                   error=error)
             q_max = SQRT(SUM((pi/dr(:))**2))
             CALL calculate_rhotot_elec_gspace(qs_env=qs_env,&
                                               auxbas_pw_pool=auxbas_pw_pool,&
                                               rhotot_elec_gspace=rho_elec_gspace,&
                                               q_max=q_max,&
                                               rho_hard=rho_hard,&
                                               rho_soft=rho_soft,&
                                               error=error)
             rho_total = rho_hard + rho_soft
             CALL get_pw_grid_info(pw_grid=rho_elec_gspace%pw%pw_grid,&
                                   vol=volume,&
                                   error=error)
             CALL pw_transfer(rho_elec_gspace%pw,rho_elec_rspace%pw,debug=.FALSE.,error=error)
             rho_total_rspace = pw_integrate_function(rho_elec_rspace%pw,isign=-1,error=error)/volume
             filename = "TOTAL_ELECTRON_DENSITY"
             unit_nr = cp_print_key_unit_nr(logger,input,"DFT%PRINT%E_DENSITY_CUBE",&
                                            extension=".cube",middle_name=TRIM(filename),&
                                            file_position=my_pos_cube,log_filename=.FALSE.,error=error)
             IF (output_unit>0) THEN
                INQUIRE (UNIT=unit_nr,NAME=filename)
                WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A)")&
                  "The total electron density is written in cube file format to the file:",&
                  TRIM(filename)
                WRITE (UNIT=output_unit,FMT="(/,(T2,A,F20.10))")&
                  "q(max) [1/Angstrom]              :",q_max/angstrom,&
                  "Soft electronic charge (G-space) :",rho_soft,&
                  "Hard electronic charge (G-space) :",rho_hard,&
                  "Total electronic charge (G-space):",rho_total,&
                  "Total electronic charge (R-space):",rho_total_rspace
              END IF
              CALL cp_pw_to_cube(rho_elec_rspace%pw,unit_nr,"TOTAL ELECTRON DENSITY",&
                   particles=particles,&
                   stride=section_get_ivals(dft_section,"PRINT%E_DENSITY_CUBE%STRIDE",error=error),&
                   error=error)
             CALL cp_print_key_finished_output(unit_nr,logger,input,&
                                               "DFT%PRINT%E_DENSITY_CUBE",&
                                               error=error)
             ! Print total spin density for spin-polarized systems
             IF (dft_control%nspins > 1) THEN
                CALL pw_zero(rho_elec_gspace%pw,error=error)
                CALL pw_zero(rho_elec_rspace%pw,error=error)
                CALL calculate_rhotot_elec_gspace(qs_env=qs_env,&
                                                  auxbas_pw_pool=auxbas_pw_pool,&
                                                  rhotot_elec_gspace=rho_elec_gspace,&
                                                  q_max=q_max,&
                                                  rho_hard=rho_hard,&
                                                  rho_soft=rho_soft,&
                                                  fsign=-1.0_dp,&
                                                  error=error)
                rho_total = rho_hard + rho_soft
                CALL pw_transfer(rho_elec_gspace%pw,rho_elec_rspace%pw,debug=.FALSE.,error=error)
                rho_total_rspace = pw_integrate_function(rho_elec_rspace%pw,isign=-1,error=error)/volume
                filename = "TOTAL_SPIN_DENSITY"
                unit_nr = cp_print_key_unit_nr(logger,input,"DFT%PRINT%E_DENSITY_CUBE",&
                                               extension=".cube",middle_name=TRIM(filename),&
                                               file_position=my_pos_cube,log_filename=.FALSE.,error=error)
                IF (output_unit>0) THEN
                   INQUIRE (UNIT=unit_nr,NAME=filename)
                   WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A)")&
                     "The total spin density is written in cube file format to the file:",&
                     TRIM(filename)
                   WRITE (UNIT=output_unit,FMT="(/,(T2,A,F20.10))")&
                     "q(max) [1/Angstrom]                    :",q_max/angstrom,&
                     "Soft part of the spin density (G-space):",rho_soft,&
                     "Hard part of the spin density (G-space):",rho_hard,&
                     "Total spin density (G-space)           :",rho_total,&
                     "Total spin density (R-space)           :",rho_total_rspace
                END IF
                CALL cp_pw_to_cube(rho_elec_rspace%pw,unit_nr,"TOTAL SPIN DENSITY",&
                     particles=particles,&
                     stride=section_get_ivals(dft_section,"PRINT%E_DENSITY_CUBE%STRIDE",error=error),&
                     error=error)
             END IF
             CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_elec_gspace%pw,error=error)
             CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_elec_rspace%pw,error=error)
          ELSE
             IF (dft_control%nspins > 1) THEN
                CALL get_qs_env(qs_env=qs_env,&
                                pw_env=pw_env,&
                                error=error)
                CALL pw_env_get(pw_env=pw_env,&
                                auxbas_pw_pool=auxbas_pw_pool,&
                                pw_pools=pw_pools,&
                                error=error)
                CALL pw_pool_create_pw(pool=auxbas_pw_pool,&
                                       pw=rho_elec_rspace%pw,&
                                       use_data=REALDATA3D,&
                                       in_space=REALSPACE,&
                                       error=error)
                CALL pw_copy(rho%rho_r(1)%pw,rho_elec_rspace%pw, error=error)
                CALL pw_axpy(rho%rho_r(2)%pw,rho_elec_rspace%pw, error=error)
                filename = "ELECTRON_DENSITY"
                unit_nr = cp_print_key_unit_nr(logger,input,"DFT%PRINT%E_DENSITY_CUBE",&
                                               extension=".cube",middle_name=TRIM(filename),&
                                               file_position=my_pos_cube,log_filename=.FALSE.,error=error)
                IF (output_unit>0) THEN
                   INQUIRE (UNIT=unit_nr,NAME=filename)
                   WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A)")&
                     "The sum of alpha and beta density is written in cube file format to the file:",&
                     TRIM(filename)
                END IF
                CALL cp_pw_to_cube(rho_elec_rspace%pw,unit_nr,"SUM OF ALPHA AND BETA DENSITY",&
                     particles=particles,stride=section_get_ivals(dft_section,"PRINT%E_DENSITY_CUBE%STRIDE",&
                     error=error),error=error)
                CALL cp_print_key_finished_output(unit_nr,logger,input,&
                                         "DFT%PRINT%E_DENSITY_CUBE",error=error)
                CALL pw_copy(rho%rho_r(1)%pw,rho_elec_rspace%pw, error=error)
                CALL pw_axpy(rho%rho_r(2)%pw,rho_elec_rspace%pw,alpha=-1.0_dp, error=error)
                filename = "SPIN_DENSITY"
                unit_nr = cp_print_key_unit_nr(logger,input,"DFT%PRINT%E_DENSITY_CUBE",&
                                               extension=".cube",middle_name=TRIM(filename),&
                                               file_position=my_pos_cube,log_filename=.FALSE.,error=error)
                IF (output_unit>0) THEN
                   INQUIRE (UNIT=unit_nr,NAME=filename)
                   WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A)")&
                     "The spin density is written in cube file format to the file:",&
                     TRIM(filename)
                END IF
                CALL cp_pw_to_cube(rho_elec_rspace%pw,unit_nr,"SPIN DENSITY",&
                     particles=particles,&
                     stride=section_get_ivals(dft_section,"PRINT%E_DENSITY_CUBE%STRIDE",error=error),&
                     error=error)
                CALL cp_print_key_finished_output(unit_nr,logger,input,&
                                                  "DFT%PRINT%E_DENSITY_CUBE",&
                                                  error=error)
                CALL pw_pool_give_back_pw(auxbas_pw_pool,rho_elec_rspace%pw,error=error)
             ELSE
                filename = "ELECTRON_DENSITY"
                unit_nr = cp_print_key_unit_nr(logger,input,"DFT%PRINT%E_DENSITY_CUBE",&
                                               extension=".cube",middle_name=TRIM(filename),&
                                               file_position=my_pos_cube,log_filename=.FALSE.,error=error)
                IF (output_unit>0) THEN
                   INQUIRE (UNIT=unit_nr,NAME=filename)
                   WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A)")&
                     "The electron density is written in cube file format to the file:",&
                     TRIM(filename)
                END IF
                CALL cp_pw_to_cube(rho%rho_r(1)%pw,unit_nr,"ELECTRON DENSITY",&
                     particles=particles,&
                     stride=section_get_ivals(dft_section,"PRINT%E_DENSITY_CUBE%STRIDE",error=error),&
                     error=error)
                CALL cp_print_key_finished_output(unit_nr,logger,input,&
                                                  "DFT%PRINT%E_DENSITY_CUBE",&
                                                  error=error)
             END IF ! nspins
          END IF ! total density for GAPW
       END IF ! print key


       ! Print the hartree potential
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%V_HARTREE_CUBE",error=error),cp_p_file)) THEN

          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,aux_r%pw,&
                                  use_data = REALDATA3D,&
                                  in_space = REALSPACE, error=error)

          append_cube = section_get_lval(input,"DFT%PRINT%V_HARTREE_CUBE%APPEND",error=error)
          my_pos_cube="REWIND"
          IF(append_cube) THEN
            my_pos_cube="APPEND"
          END IF
          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
          CALL pw_env_get(pw_env, error=error)
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%V_HARTREE_CUBE",&
               extension=".cube",middle_name="v_hartree",file_position=my_pos_cube,error=error)
          udvol = 1.0_dp/v_hartree_rspace%pw_grid%dvol

          CALL pw_copy(v_hartree_rspace,aux_r%pw, error=error)
          CALL pw_scale(aux_r%pw,udvol,error=error)

          CALL cp_pw_to_cube(aux_r%pw,unit_nr,"HARTREE POTENTIAL",particles=particles,&
               stride=section_get_ivals(dft_section,&
               "PRINT%V_HARTREE_CUBE%STRIDE",error=error),&
               error=error)
          CALL cp_print_key_finished_output(unit_nr,logger,input,&
               "DFT%PRINT%V_HARTREE_CUBE",error=error)

          CALL pw_pool_give_back_pw(auxbas_pw_pool,aux_r%pw, error=error)
       ENDIF

       ! Print the Electrical Field Components
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%EFIELD_CUBE",error=error),cp_p_file)) THEN

          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,aux_r%pw,&
                                  use_data = REALDATA3D,&
                                  in_space = REALSPACE, error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,aux_g%pw,&
                                  use_data = COMPLEXDATA1D,&
                                  in_space = RECIPROCALSPACE, error=error)

          append_cube = section_get_lval(input,"DFT%PRINT%EFIELD_CUBE%APPEND",error=error)
          my_pos_cube="REWIND"
          IF(append_cube) THEN
            my_pos_cube="APPEND"
          END IF
          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
          CALL pw_env_get(pw_env, error=error)
          udvol = 1.0_dp/v_hartree_rspace%pw_grid%dvol
          DO id=1,3
            unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%EFIELD_CUBE",&
              extension=".cube",middle_name="efield_"//cdir(id),file_position=my_pos_cube,error=error)

            CALL pw_transfer(v_hartree_rspace,aux_g%pw, error=error)
            nd=0
            nd(id)=1
            CALL pw_derive(aux_g%pw,nd,error=error)
            CALL pw_transfer(aux_g%pw,aux_r%pw, error=error)
            CALL pw_scale(aux_r%pw,udvol,error=error)

            CALL cp_pw_to_cube(aux_r%pw,&
               unit_nr,"ELECTRIC FIELD",particles=particles,&
               stride=section_get_ivals(dft_section,&
               "PRINT%EFIELD_CUBE%STRIDE",error=error),&
               error=error)
            CALL cp_print_key_finished_output(unit_nr,logger,input,&
                 "DFT%PRINT%EFIELD_CUBE",error=error)
          END DO

          CALL pw_pool_give_back_pw(auxbas_pw_pool,aux_r%pw, error=error)
          CALL pw_pool_give_back_pw(auxbas_pw_pool,aux_g%pw, error=error)
       END IF

       ! Write the density matrices
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
                 "DFT%PRINT%AO_MATRICES/DENSITY",error=error),cp_p_file)) THEN
          iw = cp_print_key_unit_nr(logger,input,"DFT%PRINT%AO_MATRICES/DENSITY",&
                                    extension=".Log",error=error)
          DO ispin=1,dft_control%nspins
             CALL cp_dbcsr_write_sparse_matrix(rho%rho_ao(ispin)%matrix,4,6,qs_env,&
                                            para_env,output_unit=iw,error=error)
          END DO
          CALL cp_print_key_finished_output(iw,logger,input,&
                                            "DFT%PRINT%AO_MATRICES/DENSITY",&
                                            error=error)
       END IF

       ! Write the Kohn-Sham matrices
       write_ks=BTEST(cp_print_key_should_output(logger%iter_info,input,&
                 "DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",error=error),cp_p_file)
       write_xc=BTEST(cp_print_key_should_output(logger%iter_info,input,&
                 "DFT%PRINT%AO_MATRICES/MATRIX_VXC",error=error),cp_p_file)
       ! we need to update stuff before writing, potentially computing the matrix_vxc
       IF (write_ks .OR. write_xc) THEN
          IF (write_xc) qs_env%requires_matrix_vxc=.TRUE.
          CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.,error=error)
          CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,error=error,&
                                    calculate_forces=.FALSE.,just_energy=.FALSE.)
          IF (write_xc) qs_env%requires_matrix_vxc=.FALSE.
       END IF

       ! Write the Kohn-Sham matrices
       IF (write_ks) THEN
          iw = cp_print_key_unit_nr(logger,input,"DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",&
                                    extension=".Log",error=error)
          DO ispin=1,dft_control%nspins
             CALL cp_dbcsr_write_sparse_matrix(ks_rmpv(ispin)%matrix,4,6,qs_env,&
                                            para_env,output_unit=iw,error=error)
          END DO
          CALL cp_print_key_finished_output(iw,logger,input,&
                                            "DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",&
                                            error=error)
       END IF

       ! Write the xc matrix
       IF (write_xc) THEN
          iw = cp_print_key_unit_nr(logger,input,"DFT%PRINT%AO_MATRICES/MATRIX_VXC",&
                                    extension=".Log",error=error)
          DO ispin=1,dft_control%nspins
             CALL cp_dbcsr_write_sparse_matrix(matrix_vxc(ispin)%matrix,4,6,qs_env,&
                                            para_env,output_unit=iw,error=error)
          END DO
          CALL cp_print_key_finished_output(iw,logger,input,&
                                            "DFT%PRINT%AO_MATRICES/MATRIX_VXC",&
                                            error=error)
       END IF

       ! Compute the Mulliken charges
       print_key => section_vals_get_subs_vals(input,"DFT%PRINT%MULLIKEN", error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%MULLIKEN",extension=".mulliken",&
                                       middle_name="",log_filename=.FALSE.,error=error)
          print_level = 1
          CALL section_vals_val_get(print_key,"PRINT_GOP",l_val=print_it,error=error)
          IF (print_it) print_level = 2
          CALL section_vals_val_get(print_key,"PRINT_ALL",l_val=print_it,error=error)
          IF (print_it) print_level = 3
          CALL mulliken_population_analysis(qs_env,unit_nr,print_level,error)
          CALL cp_print_key_finished_output(unit_nr, logger,input,"DFT%PRINT%MULLIKEN",error=error)
       END IF

    END IF

    CALL timestop(handle)

  END SUBROUTINE write_mo_free_results

! *****************************************************************************
!> \brief ...
!> \param ca ...
!> \param a ...
!> \param cb ...
!> \param b ...
!> \param l ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE project_function_a(ca,a,cb,b,l,error)
    ! project function cb on ca
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: ca
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a, cb, b
    INTEGER, INTENT(IN)                      :: l
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'project_function_a', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: info, istat, n
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipiv
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: smat, tmat, v

    failure = .FALSE.

    n = SIZE(ca)
    ALLOCATE(smat(n,n),tmat(n,n),v(n,1),ipiv(n),stat=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL sg_overlap ( smat, l, a, a )
    CALL sg_overlap ( tmat, l, a, b )
    v(:,1) = MATMUL(tmat,cb)
    CALL lapack_sgesv ( n, 1, smat, n, ipiv, v, n, info )
    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)
    ca(:) = v(:,1)
    
    DEALLOCATE(smat,tmat,v,ipiv,stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE project_function_a

! *****************************************************************************
!> \brief ...
!> \param ca ...
!> \param a ...
!> \param bfun ...
!> \param grid_atom ...
!> \param l ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE project_function_b(ca,a,bfun,grid_atom,l,error)
    ! project function f on ca
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: ca
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: a, bfun
    TYPE(grid_atom_type), POINTER            :: grid_atom
    INTEGER, INTENT(IN)                      :: l
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'project_function_b', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, info, istat, n, nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipiv
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: afun
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: smat, v

    failure=.FALSE.

    n = SIZE(ca)
    nr = grid_atom%nr
    ALLOCATE(smat(n,n),v(n,1),ipiv(n),afun(nr),stat=istat)
    CPPrecondition(istat==0,cp_failure_level,routineP,error,failure)

    CALL sg_overlap ( smat, l, a, a )
    DO i=1,n
       afun(:) = grid_atom%rad(:)**l * EXP(-a(i)*grid_atom%rad2(:))
       v(i,1) = SUM(afun(:)*bfun(:)*grid_atom%wr(:))
    END DO
    CALL lapack_sgesv ( n, 1, smat, n, ipiv, v, n, info )
    CPPostcondition(info==0,cp_failure_level,routineP,error,failure)
    ca(:) = v(:,1)

    DEALLOCATE(smat,v,ipiv,afun,stat=istat)
    CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE project_function_b

END MODULE qs_scf_post_gpw
