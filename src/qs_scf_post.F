!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Does all kind of post scf calculations
!> \par History
!>      Started as a copy from the relevant part of qs_scf
!> \author Joost VandeVondele (10.2003)
! *****************************************************************************
MODULE qs_scf_post
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type,&
                                             get_cell,&
                                             pbc
  USE cp_array_r_utils,                ONLY: cp_1d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             tddfpt_control_type
  USE cp_ddapc_util,                   ONLY: get_ddapc
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
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
  USE cp_output_handling,              ONLY: cp_iter_string,&
                                             cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_types,                 ONLY: erase_result,&
                                             get_results,&
                                             put_results
  USE cp_sm_fm_interactions,           ONLY: copy_sm_to_fm,&
                                             cp_sm_plus_fm_fm_t
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE et_coupling_types,               ONLY: set_et_coupling_type
  USE external_potential_types,        ONLY: gth_potential_type
  USE f77_blas
  USE input_constants,                 ONLY: ot_precond_full_all,&
                                             state_loc_all,&
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
  USE lowdin,                          ONLY: lowdin_charges
  USE mathconstants,                   ONLY: pi
  USE molecular_states,                ONLY: construct_molecular_states
  USE molecule_types_new,              ONLY: molecule_type
  USE mulliken,                        ONLY: mulliken_charges
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom,&
                                             debye
  USE preconditioner,                  ONLY: preconditioner_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grids,                        ONLY: get_pw_grid_info
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
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
                                             pw_p_type
  USE qs_collocate_density,            ONLY: calculate_rho_elec,&
                                             calculate_wavefunction
  USE qs_conductivity,                 ONLY: optical_conductivity
  USE qs_electric_field_gradient,      ONLY: qs_efg_calc
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_epr_hyp,                      ONLY: qs_epr_hyp_calc
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_loc_control,                  ONLY: localized_wfn_control_create,&
                                             read_loc_section,&
                                             set_loc_centers,&
                                             set_loc_wfn_lists
  USE qs_loc_methods,                  ONLY: qs_loc_driver
  USE qs_loc_molecules,                ONLY: wfc_to_molecule
  USE qs_loc_types,                    ONLY: qs_loc_env_create,&
                                             qs_loc_env_destroy,&
                                             qs_loc_env_new_type
  USE qs_loc_utils,                    ONLY: qs_loc_env_init
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             write_mo_set
  USE qs_moments,                      ONLY: qs_moment_berry_phase,&
                                             qs_moment_locop
  USE qs_ot_eigensolver,               ONLY: ot_eigensolver
  USE qs_pdos,                         ONLY: calculate_projected_dos
  USE qs_resp,                         ONLY: resp_fit
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_types,                    ONLY: ot_method_nr,&
                                             qs_scf_env_type
  USE realspace_grid_cube,             ONLY: pw_to_cube
  USE s_square_methods,                ONLY: compute_s_square
  USE scf_control_types,               ONLY: scf_control_type
  USE sparse_matrix_output,            ONLY: write_sparse_matrix
  USE sparse_matrix_types,             ONLY: deallocate_matrix,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure,&
                                             set_matrix
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xray_diffraction,                ONLY: calculate_rhotot_elec_gspace,&
                                             xray_diffraction_spectrum
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  ! Local type
! *****************************************************************************

  ! Global parameters
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_post'
  PUBLIC :: scf_post_calculation

CONTAINS

! *****************************************************************************
!> \brief collects possible post - scf calculations and prints info / computes properties.
!> \param scf_env the scf_env whose info should be written out
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      this function changes mo_eigenvectors and mo_eigenvalues, depending on the print keys.
!>      In particular, MO_CUBES causes the MOs to be rotated to make them eigenstates of the KS
!>      matrix, and mo_eigenvalues is updated accordingly. This can, for unconverged wavefunctions,
!>      change afterwards slightly the forces (hence small numerical differences between MD
!>      with and without the debug print level). Ideally this should not happen...
!> \par History
!>      02.2003 created [fawzi]
!>      10.2004 moved here from qs_scf [Joost VandeVondele]
!>              started splitting out different subroutines
!> \author fawzi
! *****************************************************************************
  SUBROUTINE scf_post_calculation(dft_section, scf_env,qs_env,error)

    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_post_calculation', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, homo, ispin, n, n_mo(2), nao, nelectron, nhomo, nlumo, &
      nlumos, nmo, nmoloc(2), nmos, output_unit
    LOGICAL :: do_mo_cubes, do_stm, do_wannier_cubes, failure, ionode, &
      my_localized_wfn, p_loc, p_molstates, uniform_occ
    REAL(KIND=dp)                            :: maxocc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues, stm_biases
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: scenter
    TYPE(cp_1d_r_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: evals
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_evals
    TYPE(cp_fm_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: mo_arrays
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_orbs
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_deriv, &
                                                mo_localized
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(preconditioner_type), POINTER       :: preconditioner
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: wf_g, wf_r
    TYPE(pw_p_type), POINTER                 :: rho_core
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv, matrix_s
    TYPE(real_matrix_type), POINTER          :: stm_density_ao
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: input, loc_print_key, &
                                                print_key
    TYPE(tddfpt_control_type), POINTER       :: tddfpt_control

    CALL timeset(routineN,"I",'',handle)

    ! Writes the data that is already available in qs_env
    para_env=>qs_env%para_env
    CALL write_available_results(qs_env,scf_env,error)

    failure=.FALSE.
    my_localized_wfn = .FALSE.
    NULLIFY(dft_control, pw_env, auxbas_pw_pool, pw_pools, mos, rho, &
         mo_coeff, ks_rmpv, matrix_s, qs_loc_env, scf_control, &
         unoccupied_orbs, mo_eigenvalues, unoccupied_evals, fm_struct_tmp, &
         molecule_set, rho_core, tddfpt_control, subsys, particles, input, &
         print_key)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    output_unit = -1
    IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    ! Here we start with data that needs a postprocessing...
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,dft_control=dft_control,molecule_set=molecule_set, &
            mos=mos,rho=rho,matrix_ks=ks_rmpv,scf_control=scf_control,matrix_s=matrix_s, &
            input=input, subsys=subsys,rho_core=rho_core,tddfpt_control=tddfpt_control,&
            mos_aux_fit=mos_aux_fit, error=error)
       CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
       CALL cp_subsys_get(subsys,particles=particles,error=error)

       ! Compute Atomic Charges
       CALL qs_scf_post_charges(input, logger, qs_env, rho, matrix_s, mos, para_env, error)

       ! Moments of charge distribution
       CALL qs_scf_post_moments(input, logger, qs_env, ionode, output_unit, error)

       ! Print Localized orbitals
       print_key => section_vals_get_subs_vals(dft_section,"PRINT%LOCALIZATION",error=error)
       p_loc=BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)
       ! Print Molecular States
       print_key => section_vals_get_subs_vals(dft_section,"PRINT%MOLECULAR_STATES",error=error)
       p_molstates=BTEST(cp_print_key_should_output(logger%iter_info,print_key,&
            error=error),cp_p_file)

       CALL section_vals_val_get(dft_section,"PRINT%MO_CUBES%STM_BIAS",r_vals=stm_biases,error=error)
       do_stm=ASSOCIATED(stm_biases)

       do_mo_cubes=BTEST(cp_print_key_should_output(logger%iter_info,dft_section,"PRINT%MO_CUBES", error=error),cp_p_file)
       do_wannier_cubes=BTEST(cp_print_key_should_output(logger%iter_info,dft_section,&
                        "PRINT%LOCALIZATION%LOCALIZE%PRINT%WANNIER_CUBES",error=error),cp_p_file)
       nlumo=section_get_ival(dft_section,"PRINT%MO_CUBES%NLUMO",error=error)
       nhomo=section_get_ival(dft_section,"PRINT%MO_CUBES%NHOMO",error=error)

       ! Setup the grids needed to compute a wavefunction given a vector..
       IF ( ( ( do_mo_cubes .OR. do_wannier_cubes ).AND. (nlumo /= 0 .OR. nhomo /= 0 .OR. do_stm )) .OR. p_loc ) THEN
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                          pw_pools=pw_pools,error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,&
                                  use_data = REALDATA3D,&
                                  in_space = REALSPACE, error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,wf_g%pw,&
                                  use_data = COMPLEXDATA1D,&
                                  in_space = RECIPROCALSPACE, error=error)
          IF (do_stm) THEN
             NULLIFY(stm_density_ao)
             CALL replicate_matrix_structure(source=rho%rho_ao(1)%matrix,&
                                       TARGET=stm_density_ao,&
                                       target_name="stm_density_ao",error=error)
          ENDIF
       END IF


       ! Makes the MOs eigenstates, computes eigenvalues, write cubes
       IF (do_mo_cubes .AND. (nhomo /= 0 .OR. do_stm)) THEN
          IF (qs_env%dft_control%restricted) THEN
             !
             ! The issue is that one state is different from the others, 
             ! so we certainly can't rotate all N+1 alpha states among themselves
             ! furthermore, alpha and beta orbitals do have the same shape, 
             ! but not the same v_xc, therefore these seems to be no reason that
             ! they should have the same eigenvalues.
             ! So it looks like I need a convincing definition of an eigenstate in this case :-)
             ! 
             IF (ionode) WRITE(output_unit,*) &
                " Unclear how we define MOs in the restricted case ... skipping"
          ELSE
             DO ispin=1,dft_control%nspins
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                                eigenvalues=mo_eigenvalues,homo=homo,nmo=nmo)
                IF (ionode) WRITE(output_unit,*) " "
                IF (ionode) WRITE(output_unit,*) " Eigenvalues of the occupied subspace spin ",ispin
                IF (homo .NE. nmo) THEN
                   IF (ionode) WRITE(output_unit,*)" and ",nmo-homo," added MO eigenvalues"
                ENDIF
                IF (ionode) WRITE(output_unit,*) "---------------------------------------------"
                ! Also rotate the OT derivs, since they are needed for force calculations
                IF (ASSOCIATED(qs_env%mo_derivs)) THEN
                   mo_coeff_deriv=>qs_env%mo_derivs(ispin)%matrix
                ELSE
                   mo_coeff_deriv=>NULL()
                ENDIF
                CALL calculate_subspace_eigenvalues(mo_coeff,ks_rmpv(ispin)%matrix,mo_eigenvalues, &
                  para_env=qs_env%para_env,scr=output_unit, ionode=ionode, do_rotation=.TRUE.,&
                  co_rotate=mo_coeff_deriv,error=error)

                ! Prints the cube files of OCCUPIED ORBITALS
                CALL qs_scf_post_occ_cubes(input, dft_section, dft_control, logger, qs_env,&
                     mo_coeff, wf_g, wf_r, particles, homo, ispin, ionode, error=error)
             ENDDO
          ENDIF
       ENDIF

       ! do_stm for the occupied states (requires eigenstates at this point)
       IF (do_mo_cubes .AND. ANY(stm_biases<0.0_dp)) THEN
          ALLOCATE(mo_arrays(dft_control%nspins))
          ALLOCATE(evals(dft_control%nspins))
          DO ispin=1,dft_control%nspins
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                  eigenvalues=mo_eigenvalues)
             mo_arrays(ispin)%matrix=>mo_coeff
             evals(ispin)%array=>mo_eigenvalues
          ENDDO
          CALL stm_cubes(qs_env, dft_section, stm_density_ao, wf_r, wf_g, mo_arrays, evals, &
                       stm_biases, .TRUE. , particles, ionode, output_unit,  error)
          DEALLOCATE(mo_arrays)
          DEALLOCATE(evals)
       ENDIF

       ! Initialize the localization environment, needed e.g. for wannier functions and molecular states
       IF (p_loc) THEN
          CALL qs_loc_env_create(qs_loc_env,error=error)
          CALL localized_wfn_control_create(qs_loc_env%localized_wfn_control,error=error)
          loc_print_key => section_vals_get_subs_vals(input,"DFT%PRINT%LOCALIZATION%LOCALIZE",error=error)
          CALL read_loc_section(qs_loc_env%localized_wfn_control,loc_print_key,qs_loc_env%do_localize,&
               error=error)
          IF (qs_loc_env%do_localize) THEN
             ! Some setup for MOs to be localized
             nmoloc = 0
             DO ispin = 1,dft_control%nspins
                CALL get_mo_set(mos(ispin)%mo_set,nmo=n_mo(ispin),nelectron=nelectron,maxocc=maxocc)
                IF(qs_loc_env%localized_wfn_control%set_of_states == state_loc_all) THEN
                   nmoloc(ispin) = NINT(nelectron/maxocc)
                ELSE
                   nmoloc(ispin) = MIN(qs_loc_env%localized_wfn_control%nloc_states(1) ,n_mo(ispin))
                ENDIF
             ENDDO  ! ispin
             CALL set_loc_wfn_lists(qs_loc_env%localized_wfn_control,nmoloc,n_mo,dft_control%nspins,error=error)
             CALL set_loc_centers(qs_loc_env%localized_wfn_control,nmoloc,dft_control%nspins,error=error)
             CALL qs_loc_env_init(qs_loc_env,qs_loc_env%localized_wfn_control,qs_env,error=error)
          ELSE
             ! Let's inform in case the section is not present in the input
             CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                  "User requested the calculation of the localized wavefunction but the section "//&
                  "LOCALIZE (located in: DFT%PRINT%LOCALIZATION%LOCALIZE) was not specified. Localization "//&
                  "will not be performed!"//&
CPSourceFileRef,&
                  only_ionode=.TRUE.)
          END IF
       END IF

       ! Gets localization info for the occupied orbs
       !  - Possibly gets wannier functions
       !  - Possibly gets molecular states
       IF (p_loc) THEN
          IF (qs_env%dft_control%restricted) THEN
             IF (ionode) WRITE(output_unit,*) &
              " Unclear how we define MOs / localization in the restricted case ... skipping"
          ELSE
             nmos=0
             DO ispin=1,dft_control%nspins
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,nmo=nmo,uniform_occupation=uniform_occ)
                CPPrecondition(uniform_occ,cp_failure_level,routineP,error,failure)
                nmos=MAX(nmos,nmo)
             END DO
             qs_loc_env%tag_mo="HOMO"
             DO ispin=1,dft_control%nspins
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                                eigenvalues=mo_eigenvalues,homo=homo,nmo=nmo)
                ! Get eigenstates
                IF (section_get_ival(dft_section,"PRINT%MO_CUBES%NHOMO",error=error)==0) THEN
                   IF (ionode) WRITE(output_unit,*) " "
                   IF (ionode) WRITE(output_unit,*) " Eigenvalues of the occupied subspace spin ",ispin
                   IF (homo .NE. nmo) THEN
                      IF (ionode) WRITE(output_unit,*) " and ",nmo-homo," added MO eigenvalues"
                   ENDIF
                   IF (ionode) WRITE(output_unit,*) "---------------------------------------------"
                   ! Also rotate the OT derivs, since they are needed for force calculations
                   IF (ASSOCIATED(qs_env%mo_derivs)) THEN
                      mo_coeff_deriv=>qs_env%mo_derivs(ispin)%matrix
                   ELSE
                      mo_coeff_deriv=>NULL()
                   ENDIF
                   CALL calculate_subspace_eigenvalues(mo_coeff,ks_rmpv(ispin)%matrix, &
                        evals_arg=mo_eigenvalues, para_env=qs_env%para_env,scr=output_unit, &
                        ionode=ionode, do_rotation=.TRUE.,error=error)
                END IF
                NULLIFY(mo_localized)
                CALL cp_fm_create(mo_localized,mo_coeff%matrix_struct,error=error)
                CALL cp_fm_to_fm(mo_coeff,mo_localized,error=error)
                IF (qs_loc_env%do_localize) THEN
                   ! Do the Real localization..
                   IF (ionode) WRITE(output_unit,"(/,T2,A,I3)") "LOCALIZATION| Computing localization properties "//&
                        "for OCCUPIED ORBITALS. Spin:",ispin

                   loc_print_key => section_vals_get_subs_vals(input,"DFT%PRINT%LOCALIZATION%LOCALIZE",&
                        error=error) 
                   scenter => qs_loc_env%localized_wfn_control%centers_set(ispin)%array
                
                   CALL qs_loc_driver(qs_env,qs_loc_env,loc_print_key,myspin=ispin,&
                        ext_mo_coeff=mo_localized,error=error)

                   ! maps wfc to molecules, and compute the molecular dipoles if required
                   IF (BTEST(cp_print_key_should_output(logger%iter_info,dft_section,"PRINT%LOCALIZATION/MOLECULAR_DIPOLES",&
                       error=error),cp_p_file) .OR. p_molstates) THEN
                      loc_print_key => section_vals_get_subs_vals(input,"DFT%PRINT%LOCALIZATION",error=error) 
                      CALL wfc_to_molecule(qs_loc_env, loc_print_key, scenter, molecule_set, error)
                   END IF

                   ! Compute the molecular states
                   IF ( p_molstates ) THEN
                      CALL construct_molecular_states(molecule_set, mo_localized, mo_coeff, &
                           mo_eigenvalues, ks_rmpv(ispin)%matrix, matrix_s(1)%matrix, qs_env, &
                           wf_r, wf_g, force_env_section=qs_env%input, particles=particles, error=error)
                   ENDIF
                END IF
                CALL cp_fm_release(mo_localized,error=error)
             ENDDO
             ! Print Total Dipole if the localization has been performed
             IF (qs_loc_env%do_localize) THEN
                CALL qs_scf_post_loc_dip(input, dft_control, qs_loc_env, logger, qs_env, error)
             END IF
          ENDIF
       ELSE
          ! Print a warning on screen if molecular states have been required but the LOCALIZATION
          ! section is not present
          CALL cp_assert(.NOT.p_molstates,cp_warning_level,cp_assertion_failed,routineP,&
               "Molecular States have been required but the LOCALIZATION keyword is not present! "//&
               "Skipping the calculation of the Molecular States. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
       ENDIF

       ! Gets the lumos, and eigenvalues for the lumos
       IF ((do_mo_cubes .OR.dft_control%do_tddfpt_calculation ) .AND. nlumo .NE. 0 ) THEN 
          ALLOCATE(unoccupied_orbs(dft_control%nspins))
          ALLOCATE(unoccupied_evals(dft_control%nspins))
          DO ispin=1,dft_control%nspins
             NULLIFY(unoccupied_orbs(ispin)%matrix)
             NULLIFY(unoccupied_evals(ispin)%array)
             ! Always write eigenvalues
             IF (ionode) WRITE(output_unit,*) " "
             IF (ionode) WRITE(output_unit,*) " Lowest Eigenvalues of the unoccupied subspace spin ",ispin
             IF (ionode) WRITE(output_unit,FMT='(1X,A)') "-----------------------------------------------------"
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,homo=homo,nao=nao,nmo=nmo)
             CALL cp_fm_get_info(mo_coeff, nrow_global=n,error=error)
             nlumos=MAX(1,MIN(nlumo,nao-homo))
             IF (nlumo==-1) nlumos=nao-homo
             ALLOCATE(unoccupied_evals(ispin)%array(nlumos))
             CALL cp_fm_struct_create(fm_struct_tmp,para_env=qs_env%para_env,context=qs_env%blacs_env, &
                                                    nrow_global=n,ncol_global=nlumos,error=error)
             CALL cp_fm_create(unoccupied_orbs(ispin)%matrix, fm_struct_tmp,name="lumos",error=error)
             CALL cp_fm_struct_release(fm_struct_tmp,error=error)
             CALL cp_fm_init_random(unoccupied_orbs(ispin)%matrix,nlumos,error=error)

             ! the full_all preconditioner makes not much sense for lumos search
             NULLIFY(preconditioner)
             IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
                preconditioner=>scf_env%ot_preconditioner(1)%preconditioner
                ! this one can for sure not be right (as it has to match a given C0)
                IF (preconditioner%in_use == ot_precond_full_all) THEN
                    NULLIFY(preconditioner)
                ENDIF
             ENDIF

             CALL ot_eigensolver(matrix_h=ks_rmpv(ispin)%matrix,matrix_s=matrix_s(1)%matrix, &
                  matrix_c=unoccupied_orbs(ispin)%matrix, &
                  matrix_orthogonal_space=mo_coeff, &
                  eps_gradient=scf_control%eps_lumos, &
                  preconditioner=preconditioner, &
                  iter_max=scf_control%max_iter_lumos,&
                  size_ortho_space=homo,error=error)

             CALL calculate_subspace_eigenvalues(unoccupied_orbs(ispin)%matrix,ks_rmpv(ispin)%matrix,&
                  unoccupied_evals(ispin)%array, para_env=qs_env%para_env,scr=output_unit, &
                  ionode=ionode,error=error)

             ! Prints the cube files of UNOCCUPIED ORBITALS
             CALL qs_scf_post_unocc_cubes(input, dft_section, dft_control, logger, qs_env,&
                  unoccupied_orbs, wf_g, wf_r, particles, nlumos, homo, ispin, ionode, error)
          ENDDO

          ! Save the info for tddfpt calculation
          IF (dft_control%do_tddfpt_calculation) THEN
             ALLOCATE(tddfpt_control%lumos_eigenvalues(nlumos,dft_control%nspins))
             DO ispin=1, dft_control%nspins
                tddfpt_control%lumos_eigenvalues(1:nlumos,ispin) = &
                     unoccupied_evals(ispin)%array(1:nlumos)
             END DO
             tddfpt_control%lumos => unoccupied_orbs
          END IF

       ENDIF

       ! do_stm for the unoccupied states (requires eigenstates at this point)
       IF (do_mo_cubes .AND. ANY(stm_biases>0.0_dp)) THEN
          CALL stm_cubes(qs_env, dft_section, stm_density_ao, wf_r, wf_g, unoccupied_orbs, unoccupied_evals, &
                       stm_biases, .FALSE. , particles, ionode, output_unit,  error)
       ENDIF


       ! Deallocate grids needed to compute wavefunctions
       IF ( ( ( do_mo_cubes .OR. do_wannier_cubes ).AND. (nlumo /= 0 .OR. nhomo /= 0 .OR. do_stm )) .OR. p_loc ) THEN
          CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_r%pw, error=error)
          CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_g%pw, error=error)
          IF (do_stm) THEN
             CALL deallocate_matrix(stm_density_ao,error=error)
          ENDIF
       END IF

       ! Destroy the localization environment
       IF (p_loc) THEN
          CALL qs_loc_env_destroy(qs_loc_env, error=error)
       ENDIF

       ! Compute the optical conductivity (needs homos and lumos)
       CALL qs_scf_post_optc(input, dft_section, dft_control, logger, qs_env,&
            unoccupied_orbs, unoccupied_evals, ionode, output_unit, error)

       ! This is just a deallocation for printing MO_CUBES or TDDFPT
       IF (( do_mo_cubes .OR. dft_control%do_tddfpt_calculation ) .AND. nlumo /= 0 ) THEN
          DO ispin=1,dft_control%nspins
             DEALLOCATE(unoccupied_evals(ispin)%array)
             IF (.NOT.dft_control%do_tddfpt_calculation) &
                  CALL cp_fm_release(unoccupied_orbs(ispin)%matrix,error=error)
          ENDDO
          DEALLOCATE(unoccupied_evals)
          IF (.NOT.dft_control%do_tddfpt_calculation) DEALLOCATE(unoccupied_orbs)
       ENDIF
       
       ! Print coherent X-ray diffraction spectrum
       CALL qs_scf_post_xray(input, dft_section, logger, qs_env, ionode, output_unit, error)
       
       ! Calculation of Electric Field Gradients
       CALL qs_scf_post_efg(input, logger, qs_env, dft_control, error)
       
       ! Calculation of EPR Hyperfine Coupling Tensors
       CALL qs_scf_post_epr(input, logger, qs_env, error)

       ! Calculation of properties needed for BASIS_MOLOPT optimizations
       CALL qs_scf_post_molopt(input, logger, qs_env, error)

    END IF
    
    CALL timestop(0.0_dp,handle)

  END SUBROUTINE scf_post_calculation

! *****************************************************************************
!> \brief computes a simple approximation to the tunneling current for STM
!> \param 
!> \note
!>      requires the MOs that are passed to be eigenstates, and energy ordered
!> \par History
!>      7.2008 Created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE stm_cubes(qs_env, dft_section, stm_density_ao, wf_r, wf_g, mo_arrays, evals, &
                       stm_biases, do_negative, particles, ionode, output_unit,  error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(real_matrix_type), POINTER          :: stm_density_ao
    TYPE(pw_p_type)                          :: wf_r, wf_g
    TYPE(cp_fm_p_type), DIMENSION(:), &
      INTENT(IN)                             :: mo_arrays
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
      INTENT(IN)                             :: evals
    REAL(KIND=dp), DIMENSION(:), POINTER     :: stm_biases
    LOGICAL                                  :: do_negative
    TYPE(particle_list_type), POINTER        :: particles
    LOGICAL                                  :: ionode
    INTEGER                                  :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'stm_cubes', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=80)                        :: filename, title
    INTEGER                                  :: handle, homo, ibias, ispin, &
                                                istates, nstates, &
                                                state_start, unit_nr
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha, e_end, e_start, &
                                                total_rho
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_v
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,"I",'',handle)

    failure=.FALSE.
    logger => cp_error_get_logger(error)
    NULLIFY(fm_struct_tmp)

    DO ibias=1,SIZE(stm_biases)
       IF (do_negative) THEN
          IF (stm_biases(ibias)>0.0_dp) CYCLE
       ELSE
          IF (stm_biases(ibias)<0.0_dp) CYCLE
       ENDIF

       IF (ionode) WRITE(output_unit,'(T2,A)')       ""
       IF (ionode) WRITE(output_unit,'(T2,A,F16.6)') "Preparing for STM image at bias [a.u.] ",stm_biases(ibias)

       ! a shared starting energy for alpha and beta
       IF (do_negative) THEN
          e_start=-HUGE(0.0_dp)
       ELSE
          e_start=+HUGE(0.0_dp)
       ENDIF
       DO ispin=1,SIZE(mo_arrays)
          IF (do_negative) THEN
             e_start=MAX(e_start,MAXVAL(evals(ispin)%array)) ! HOMO
          ELSE
             e_start=MIN(e_start,MINVAL(evals(ispin)%array)) ! LUMO
          ENDIF
       ENDDO

       istates=0
       DO ispin=1,SIZE(mo_arrays)
          IF (do_negative) THEN
             e_end =MINVAL(evals(ispin)%array)
             nstates=COUNT(evals(ispin)%array>e_start+stm_biases(ibias))
             IF (ionode .AND. e_start+stm_biases(ibias) < e_end) &
                 WRITE(output_unit,'(T4,A)') "Warning: HOMO+bias below lowest computed occupied MO"
          ELSE
             e_end =MAXVAL(evals(ispin)%array)
             nstates=COUNT(evals(ispin)%array<e_start+stm_biases(ibias))
             IF (ionode .AND. e_start+stm_biases(ibias) > e_end) &
                 WRITE(output_unit,'(T4,A)') "Warning: LUMO+bias above highest computed unoccupied MO"
          ENDIF
          istates=istates+nstates
       ENDDO
       IF (ionode) WRITE(output_unit,'(T4,A,I0,A)') "Using a total of ",istates," states"

       CALL cp_fm_struct_create(fm_struct_tmp, ncol_global=istates, &
                       template_fmstruct=mo_arrays(1)%matrix%matrix_struct, error=error)
       CALL cp_fm_create(matrix_v,fm_struct_tmp, name="matrix_v",error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)

       ! we sum both alpha and beta electrons together for this density of states
       istates=0
       DO ispin=1,SIZE(mo_arrays)
          IF (do_negative) THEN
             nstates=COUNT(evals(ispin)%array>e_start+stm_biases(ibias))
             state_start=SIZE(evals(ispin)%array)-nstates+1
          ELSE
             nstates=COUNT(evals(ispin)%array<e_start+stm_biases(ibias))
             state_start=1
          ENDIF
          CALL cp_fm_to_fm(mo_arrays(ispin)%matrix,matrix_v,nstates,state_start,istates+1)
          istates=istates+nstates
       ENDDO

       ! assume integer occupation, and construct the density matrix
       alpha=1.0_dp
       IF (SIZE(mo_arrays)==1) alpha=2.0_dp
       CALL set_matrix(stm_density_ao,0.0_dp)
       CALL cp_sm_plus_fm_fm_t(stm_density_ao,matrix_v=matrix_v,ncol=istates,&
                               alpha=alpha,error=error)
       CALL calculate_rho_elec(matrix_p=stm_density_ao,&
               rho=wf_r,rho_gspace=wf_g, total_rho=total_rho,&
               qs_env=qs_env, error=error)

       CALL cp_fm_release(matrix_v, error=error)

       WRITE(filename,'(a4,I5.5)')"STM_",ibias
       unit_nr=cp_print_key_unit_nr(logger,dft_section,"PRINT%MO_CUBES",extension=".cube",&
            middle_name=TRIM(filename),log_filename=.FALSE.,error=error)
       WRITE(title,*) "STM cube ",ibias," at bias ",stm_biases(ibias)
       CALL pw_to_cube(wf_r%pw,unit_nr,ionode,title,particles=particles,&
            stride=section_get_ivals(dft_section,"PRINT%MO_CUBES%STRIDE",error=error),&
            error=error)
       CALL cp_print_key_finished_output(unit_nr,logger,dft_section,"PRINT%MO_CUBES",error=error)

    ENDDO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE stm_cubes

! *****************************************************************************
!> \brief Computes and Prints Atomic Charges with several methods
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_charges(input, logger, qs_env, rho, matrix_s, mos, para_env, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_charges', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: lowdin_title, mulliken_title
    INTEGER                                  :: unit_nr
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: density_fit_section, print_key

    NULLIFY(particle_set)
    CALL get_qs_env(qs_env=qs_env,&
                    particle_set=particle_set,&
                    error=error)
    
    ! Compute the mulliken charges
    print_key => section_vals_get_subs_vals(input,"DFT%PRINT%MULLIKEN", error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
       unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%MULLIKEN",extension=".mulliken",&
            middle_name="",log_filename=.FALSE.,error=error)
       mulliken_title="  Mulliken charges : iteration level: "//cp_iter_string(logger%iter_info, print_key,error=error)
       CALL mulliken_charges(rho%rho_ao,matrix_s(1)%matrix,para_env,particle_set,unit_nr,mulliken_title)

       CALL cp_print_key_finished_output(unit_nr, logger,input,"DFT%PRINT%MULLIKEN",error=error)
    END IF

    ! Compute the Loewdin charges
    print_key => section_vals_get_subs_vals(input,"DFT%PRINT%LOWDIN", error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
       unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%LOWDIN",extension=".lowdin",&
            log_filename=.FALSE.,error=error)
       lowdin_title="  Lowdin charges : iteration level: "//cp_iter_string(logger%iter_info, print_key,error=error)
       CALL lowdin_charges(matrix_s(1)%matrix, mos, particle_set, unit_nr, lowdin_title, error=error)
       CALL cp_print_key_finished_output(unit_nr, logger,input,"DFT%PRINT%LOWDIN", error=error)
    END IF
    
    ! Compute the RESP charges
    CALL resp_fit(qs_env,error)

    ! Compute the Density Derived Atomic Point charges with the Bloechl scheme
    print_key => section_vals_get_subs_vals(input,"PROPERTIES%FIT_CHARGE", error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
       unit_nr=cp_print_key_unit_nr(logger,input,"PROPERTIES%FIT_CHARGE",extension=".Fitcharge",&
            log_filename=.FALSE.,error=error)
       density_fit_section => section_vals_get_subs_vals(input,"DFT%DENSITY_FITTING",error=error)
       CALL get_ddapc(qs_env,.FALSE.,density_fit_section,iwc=unit_nr,error=error)
       CALL cp_print_key_finished_output(unit_nr, logger,input,"PROPERTIES%FIT_CHARGE", error=error)
    END IF

  END SUBROUTINE qs_scf_post_charges

! *****************************************************************************
!> \brief Computes and prints the Dipole (using localized charges)
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

    CHARACTER(LEN=80)                        :: description
    CHARACTER(LEN=default_string_length)     :: iter
    INTEGER                                  :: i, ispins, n_rep, unit_nr
    LOGICAL                                  :: failure, first_time, ghost
    REAL(KIND=dp)                            :: zeff, zwfc
    REAL(KIND=dp), DIMENSION(3)              :: dipole, dipole_old, r_pbc
    TYPE(cell_type), POINTER                 :: cell
    TYPE(gth_potential_type), POINTER        :: gth_potential
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

    failure = .FALSE.
    IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
         "DFT%PRINT%LOCALIZATION/TOTAL_DIPOLE",first_time=first_time,error=error),cp_p_file)) THEN
       NULLIFY(gth_potential, cell, particle_set, gth_potential)
       CALL get_qs_env(qs_env=qs_env,&
                       cell=cell,&
                       particle_set=particle_set,&
                       error=error)
       
       description='[DIPOLE]'
       dipole=0.0_dp
       ! Charges of the atoms involved
       DO i=1,SIZE(particle_set)
          CALL get_atomic_kind(particle_set(i)%atomic_kind,ghost=ghost,&
               gth_potential=gth_potential)
          zeff = 0.0_dp
          IF (.NOT.ghost) THEN
             CALL get_atomic_kind(particle_set(i)%atomic_kind,core_charge=zeff)
             CPPrecondition(cell%orthorhombic,cp_failure_level,routineP,error,failure)
             r_pbc = pbc(particle_set(i)%r, cell)
             dipole=dipole + zeff*r_pbc
          END IF
       END DO

       ! Charges of the wfc involved
       ! Warning, this assumes the same occupation for all states
       zwfc = 3.0_dp-REAL(dft_control%nspins,dp)
       
       DO ispins=1,dft_control%nspins
          DO i=1,SIZE(qs_loc_env%localized_wfn_control%centers_set(ispins)%array,2)
             r_pbc = pbc(qs_loc_env%localized_wfn_control%centers_set(ispins)%array(1:3,i),cell)
             dipole=dipole - zwfc * r_pbc
          ENDDO
       ENDDO
       
       unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%LOCALIZATION",&
            extension=".Dipole",log_filename=.FALSE.,on_file=.TRUE.,error=error)
       IF (unit_nr>0) THEN
          IF (first_time) THEN
             WRITE(unit=unit_nr,fmt="(A,T31,A,T88,A,T136,A)")&
                  "# iter_level","dipole(x,y,z)[atomic units]",&
                                 "dipole(x,y,z)[debye]",&
                                 "delta_dipole(x,y,z)[atomic units]"
          END IF
          iter=cp_iter_string(logger%iter_info,error=error)
          CALL get_results(qs_env%results,description,n_rep=n_rep,error=error)
          IF(n_rep==0)THEN
             dipole_old=0._dp
          ELSE
             CALL get_results(qs_env%results,description,dipole_old,nval=n_rep,error=error)
          END IF
          WRITE(unit=unit_nr,fmt="(a,9(es18.8))")&
               iter(1:15), dipole, dipole*debye, (dipole-dipole_old)
       END IF
       CALL cp_print_key_finished_output(unit_nr,logger,input,&
            "DFT%PRINT%LOCALIZATION",on_file=.TRUE.,error=error)
       CALL erase_result(qs_env%results,description,error=error)
       CALL put_results(qs_env%results,description,dipole,error)
    END IF

  END SUBROUTINE qs_scf_post_loc_dip

! *****************************************************************************
!> \brief Computes and prints the Cube Files for MO
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_occ_cubes(input, dft_section, dft_control, logger, qs_env,&
       mo_coeff, wf_g, wf_r, particles, homo, ispin, ionode, basis_set_id, error)
    TYPE(section_vals_type), POINTER         :: input, dft_section
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(pw_p_type)                          :: wf_g, wf_r
    TYPE(particle_list_type), POINTER        :: particles
    INTEGER, INTENT(IN)                      :: homo, ispin
    LOGICAL, INTENT(IN)                      :: ionode
    INTEGER, INTENT(IN), OPTIONAL            :: basis_set_id
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_occ_cubes', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename, title
    INTEGER                                  :: ivector, my_basis_set_id, &
                                                nhomo, unit_nr
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env

    IF(PRESENT(basis_set_id)) THEN
      my_basis_set_id = basis_set_id
    ELSE
      my_basis_set_id = use_orb_basis_set
    END IF

    IF (BTEST(cp_print_key_should_output(logger%iter_info,dft_section,"PRINT%MO_CUBES",&
         error=error),cp_p_file) .AND. section_get_lval(dft_section,"PRINT%MO_CUBES%WRITE_CUBE",error=error)) THEN
       nhomo=section_get_ival(dft_section,"PRINT%MO_CUBES%NHOMO",error=error)
       IF (nhomo==-1) nhomo=homo
       DO ivector=MAX(1,homo-nhomo+1),homo
          CALL get_qs_env(qs_env=qs_env,&
                          atomic_kind_set=atomic_kind_set,&
                          cell=cell,&
                          particle_set=particle_set,&
                          pw_env=pw_env,&
                          error=error)
          CALL calculate_wavefunction(mo_coeff,ivector,wf_r,wf_g,atomic_kind_set,&
               cell,dft_control,particle_set,pw_env ,basis_set_id=my_basis_set_id, error=error)
          WRITE(filename,'(a4,I5.5,a1,I1.1)')"WFN_",ivector,"_",ispin
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%MO_CUBES",extension=".cube",&
               middle_name=TRIM(filename),log_filename=.FALSE.,error=error)
          WRITE(title,*) "WAVEFUNCTION ",ivector," spin ",ispin," i.e. HOMO - ",ivector-homo
          CALL pw_to_cube(wf_r%pw,unit_nr,ionode,title,particles=particles,&
               stride=section_get_ivals(dft_section,"PRINT%MO_CUBES%STRIDE",error=error),&
               error=error)
          CALL cp_print_key_finished_output(unit_nr,logger,input,"DFT%PRINT%MO_CUBES",error=error)
       ENDDO
    END IF

  END SUBROUTINE qs_scf_post_occ_cubes

! *****************************************************************************
!> \brief Computes and prints the Cube Files for MO
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_unocc_cubes(input, dft_section, dft_control, logger, qs_env,&
       unoccupied_orbs, wf_g, wf_r, particles, nlumos, homo, ispin, ionode, error)
    TYPE(section_vals_type), POINTER         :: input, dft_section
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_orbs
    TYPE(pw_p_type)                          :: wf_g, wf_r
    TYPE(particle_list_type), POINTER        :: particles
    INTEGER, INTENT(IN)                      :: nlumos, homo, ispin
    LOGICAL, INTENT(IN)                      :: ionode
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_unocc_cubes', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename, title
    INTEGER                                  :: ivector, unit_nr
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env

    IF (BTEST(cp_print_key_should_output(logger%iter_info,dft_section,"PRINT%MO_CUBES",error=error),cp_p_file)&
         .AND.  section_get_lval(dft_section,"PRINT%MO_CUBES%WRITE_CUBE",error=error) ) THEN
       NULLIFY(atomic_kind_set, particle_set, pw_env, cell)
       DO ivector=1,nlumos
          CALL get_qs_env(qs_env=qs_env,&
                          atomic_kind_set=atomic_kind_set,&
                          cell=cell,&
                          particle_set=particle_set,&
                          pw_env=pw_env,&
                          error=error)
          CALL calculate_wavefunction(unoccupied_orbs(ispin)%matrix, ivector, wf_r, wf_g, atomic_kind_set,&
               cell, dft_control, particle_set, pw_env, error=error)
          WRITE(filename,'(a4,I5.5,a1,I1.1)')"WFN_",homo+ivector,"_",ispin
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%MO_CUBES",extension=".cube",&
               middle_name=TRIM(filename),log_filename=.FALSE.,error=error)
          WRITE(title,*) "WAVEFUNCTION ",homo+ivector," spin ",ispin," i.e. LUMO + ",ivector-1
          CALL pw_to_cube(wf_r%pw, unit_nr, ionode, title, particles=particles,&
               stride=section_get_ivals(dft_section,"PRINT%MO_CUBES%STRIDE",error=error),&
               error=error)
          CALL cp_print_key_finished_output(unit_nr,logger,input,"DFT%PRINT%MO_CUBES",error=error)
       ENDDO
    ENDIF

  END SUBROUTINE qs_scf_post_unocc_cubes

! *****************************************************************************
!> \brief Computes and prints the optical conductivity
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_optc(input, dft_section, dft_control, logger, qs_env,&
       unoccupied_orbs, unoccupied_evals, ionode, output_unit, error)
    TYPE(section_vals_type), POINTER         :: input, dft_section
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_orbs
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_evals
    LOGICAL, INTENT(IN)                      :: ionode
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_optc', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: homo, ispin, nmo, unit_nr
    LOGICAL                                  :: homoEQnmo
    REAL(KIND=dp)                            :: volume
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s

    IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
         "DFT%PRINT%OPTICAL_CONDUCTIVITY",error=error),cp_p_file)) THEN
       NULLIFY(mos, cell, mo_coeff, mo_eigenvalues)
       CALL get_qs_env(qs_env, mos=mos, cell=cell, matrix_s=matrix_s,error=error)
       IF (ionode) WRITE(output_unit,*) " "
       IF (ionode) WRITE(output_unit,*) " Computing the optical conductivity "
       IF (ionode) WRITE(output_unit,*) " Experimental version "
       IF (ionode) WRITE(output_unit,*) " Check the code before believing results "
       homoEQnmo=.TRUE.
       DO ispin=1,dft_control%nspins
          CALL get_mo_set(mo_set=mos(ispin)%mo_set,homo=homo,nmo=nmo)
          IF (homo.NE.nmo) homoEQnmo=.FALSE.
       ENDDO
       IF (.NOT. homoEQnmo) THEN
          IF (ionode) WRITE(output_unit,*) " homo.NE.nmo : skip optical conductivty "
       ELSE
          IF (.NOT.( BTEST(cp_print_key_should_output(logger%iter_info,dft_section,"PRINT%MO_CUBES",&
               error=error),cp_p_file) )) &
               CALL stop_program("optical_conductivity","Needs MO_CUBES to be activated")
          filename="CONDUCTIVITY"
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%OPTICAL_CONDUCTIVITY",&
               extension=".data",middle_name=TRIM(filename),log_filename=.FALSE.,error=error)
          CALL get_cell(cell,deth=volume)

          DO ispin=1,dft_control%nspins
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                             eigenvalues=mo_eigenvalues)
             CALL optical_conductivity(matrix_s, mo_coeff, mo_eigenvalues, unoccupied_orbs(ispin)%matrix,&
                                       unoccupied_evals(ispin)%array, volume, ionode, unit_nr, error=error)
          ENDDO
          CALL cp_print_key_finished_output(unit_nr,logger,input,"DFT%PRINT%OPTICAL_CONDUCTIVITY",&
               error=error)
          IF (ionode) WRITE(output_unit,*) " Results written to file ",filename
       ENDIF
    ENDIF
  END SUBROUTINE qs_scf_post_optc

! *****************************************************************************
!> \brief Computes and prints electric moments
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_moments(input, logger, qs_env, ionode, output_unit, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: ionode
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_moments', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: maxmom, reference, unit_nr
    LOGICAL                                  :: periodic
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ref_point
    TYPE(section_vals_type), POINTER         :: print_key

    print_key => section_vals_get_subs_vals(section_vals=input,&
         subsection_name="DFT%PRINT%MOMENTS",error=error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN

       maxmom = section_get_ival(section_vals=input,&
             keyword_name="DFT%PRINT%MOMENTS%MAX_MOMENT",error=error)
       periodic = section_get_lval(section_vals=input,&
             keyword_name="DFT%PRINT%MOMENTS%PERIODIC",error=error)
       reference = section_get_ival(section_vals=input,&
             keyword_name="DFT%PRINT%MOMENTS%REFERENCE",error=error)
       NULLIFY ( ref_point )
       CALL section_vals_val_get(input,"DFT%PRINT%MOMENTS%REF_POINT",r_vals=ref_point,error=error)
       unit_nr = cp_print_key_unit_nr(logger=logger,basis_section=input,&
             print_key_path="DFT%PRINT%MOMENTS",extension=".dat",&
             middle_name="moments",log_filename=.FALSE.,error=error)

       IF (ionode) THEN
          IF(unit_nr /= output_unit) THEN
            INQUIRE (UNIT=unit_nr,NAME=filename)
            WRITE (UNIT=output_unit,FMT="(/,T2,A,2(/,T3,A),/)")&
                 "MOMENTS","The electric moments are written to file:",&
                 TRIM(filename)
          ELSE
            WRITE (UNIT=output_unit,FMT="(/,T2,A)") "ELECTRIC MOMENTS"
          END IF
       END IF

       IF (periodic) THEN
          CALL qs_moment_berry_phase(qs_env,maxmom,reference,ref_point,unit_nr,error)
       ELSE
          CALL qs_moment_locop(qs_env,maxmom,reference,ref_point,unit_nr,error)
       END IF

       CALL cp_print_key_finished_output(unit_nr=unit_nr,logger=logger,&
            basis_section=input,print_key_path="DFT%PRINT%MOMENTS",&
            error=error)
    END IF

  END SUBROUTINE qs_scf_post_moments

! *****************************************************************************
!> \brief Computes and prints the X-ray diffraction spectrum.
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_xray(input,dft_section,logger,qs_env,ionode,&
                              output_unit,error)

    TYPE(section_vals_type), POINTER         :: input, dft_section
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN)                      :: ionode
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_post_xray', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: unit_nr
    REAL(KIND=dp)                            :: q_max
    TYPE(section_vals_type), POINTER         :: print_key

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
      IF (ionode) THEN
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

  END SUBROUTINE qs_scf_post_xray

! *****************************************************************************
!> \brief Computes and prints Electric Field Gradient
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_efg(input, logger, qs_env, dft_control, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_efg', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, stat
    LOGICAL                                  :: do_et
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: my_mos
    TYPE(section_vals_type), POINTER         :: et_section, print_key

    print_key => section_vals_get_subs_vals(section_vals=input,&
                   subsection_name="DFT%PRINT%ELECTRIC_FIELD_GRADIENT",&
                   error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),&
              cp_p_file)) THEN
       CALL qs_efg_calc(qs_env=qs_env,error=error)
    END IF

    do_et=.FALSE.
    et_section =>  section_vals_get_subs_vals(input,"PROPERTIES%ET_COUPLING",&
                                                  error=error)
    CALL section_vals_get(et_section,explicit=do_et,error=error)
    IF(do_et)THEN
       IF(qs_env%et_coupling%first_run)THEN
          NULLIFY(my_mos)
          ALLOCATE(my_mos(dft_control%nspins),STAT=stat)
          ALLOCATE(qs_env%et_coupling%et_mo_coeff(dft_control%nspins),STAT=stat)          
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

  END SUBROUTINE qs_scf_post_efg

! *****************************************************************************
!> \brief computes the condition number of the overlap matrix and 
!>      prints the value of the total energy. This is needed
!>      for BASIS_MOLOPT optimizations
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

    INTEGER                                  :: nao, unit_nr
    REAL(KIND=dp)                            :: S_cond_number
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: eigenvalues
    TYPE(cp_fm_struct_type), POINTER         :: ao_ao_fmstruct
    TYPE(cp_fm_type), POINTER                :: fm_s, fm_work, mo_coeff
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(section_vals_type), POINTER         :: print_key

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

       CALL copy_sm_to_fm(matrix_s(1)%matrix,fm_s,error=error)
       CALL cp_fm_syevd(fm_s,fm_work,eigenvalues,error=error)
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

  END SUBROUTINE qs_scf_post_molopt

! *****************************************************************************
!> \brief Dumps EPR 
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

    TYPE(section_vals_type), POINTER         :: print_key

    print_key => section_vals_get_subs_vals(section_vals=input,&
                   subsection_name="DFT%PRINT%HYPERFINE_COUPLING_TENSOR",&
                   error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),&
              cp_p_file)) THEN
       CALL qs_epr_hyp_calc(qs_env=qs_env,error=error)
    END IF

  END SUBROUTINE qs_scf_post_epr

! *****************************************************************************
!> \brief Write QS results always available (if switched on through the print_keys)
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE write_available_results(qs_env,scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_available_results', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    INTEGER                                  :: handle, homo, ispin, iw, nmo, &
                                                output_unit, unit_nr
    LOGICAL                                  :: all_equal, failure, ionode, &
                                                print_total_density
    REAL(KIND=dp) :: maxocc, q_max, rho_hard, rho_soft, rho_total, &
      rho_total_rspace, s_square, s_square_ideal, total_abs_spin_dens, udvol, &
      volume
    REAL(KIND=dp), DIMENSION(3)              :: dr
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues, &
                                                occupation_numbers
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_deriv
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rho_elec_gspace, &
                                                rho_elec_rspace, wf_r
    TYPE(pw_p_type), POINTER                 :: rho0_s_gs, rho_core
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: dft_section, input

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(cell,dft_control,pw_env,auxbas_pw_pool,pw_pools,&
            mo_coeff, mo_coeff_deriv,mo_eigenvalues,mos,atomic_kind_set,particle_set,rho,&
            ks_rmpv,matrix_s,scf_control,dft_section,&
            molecule_set,input,particles,subsys)
    para_env=>qs_env%para_env
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,dft_control=dft_control,molecule_set=molecule_set, &
            mos=mos,atomic_kind_set=atomic_kind_set,particle_set=particle_set,&
            rho=rho,matrix_ks=ks_rmpv,scf_control=scf_control,matrix_s=matrix_s,&
            input=input,cell=cell,subsys=subsys,error=error)
       CALL cp_subsys_get(subsys,particles=particles,error=error)

       ! *** if the dft_section tells you to do so, write last wavefunction to screen
       dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
       IF (dft_control%nspins == 2) THEN
         CALL write_mo_set(mos(1)%mo_set,atomic_kind_set,particle_set,4,6,&
                           dft_section,spin="ALPHA",last=.TRUE.,error=error)
         CALL write_mo_set(mos(2)%mo_set,atomic_kind_set,particle_set,4,6,&
                           dft_section,spin="BETA",last=.TRUE.,error=error)
       ELSE
         CALL write_mo_set(mos(1)%mo_set,atomic_kind_set,particle_set,4,6,&
                           dft_section,last=.TRUE.,error=error)
       END IF

       ! *** at the end of scf print out the projected dos per kind
       IF (BTEST(cp_print_key_should_output(logger%iter_info,dft_section,"PRINT%PDOS",&
           error=error),cp_p_file) ) THEN
         DO ispin = 1,dft_control%nspins
            IF (scf_env%method == ot_method_nr) THEN
              CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                              eigenvalues=mo_eigenvalues)
              IF (ASSOCIATED(qs_env%mo_derivs)) THEN
                     mo_coeff_deriv=>qs_env%mo_derivs(ispin)%matrix
              ELSE
                     mo_coeff_deriv=>NULL()
              ENDIF
  
              CALL calculate_subspace_eigenvalues(mo_coeff,ks_rmpv(ispin)%matrix,mo_eigenvalues, &
                   para_env=qs_env%para_env, do_rotation=.TRUE.,&
                   co_rotate=mo_coeff_deriv,error=error)
            END IF
            IF(dft_control%nspins==2) THEN
              CALL calculate_projected_dos(mos(ispin)%mo_set,atomic_kind_set,&
                   particle_set,qs_env, scf_env, dft_section,ispin=ispin,error=error)
            ELSE
              CALL calculate_projected_dos(mos(ispin)%mo_set,atomic_kind_set,&
                   particle_set,qs_env,scf_env,dft_section,error=error)
            END IF
         END DO
       ENDIF

       !   *** Integrated absolute spin density and spin contamination ***
       IF (dft_control%nspins.eq.2) THEN
          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
               pw_pools=pw_pools,error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          CALL pw_copy(rho%rho_r(1)%pw,wf_r%pw)
          CALL pw_axpy(rho%rho_r(2)%pw,wf_r%pw,alpha=-1._dp)
          total_abs_spin_dens=pw_integrate_function(wf_r%pw,OPERATOR="ABS")
          IF (ionode) WRITE(UNIT=output_unit,fmt='(/,(T3,A,T60,F20.10))') &
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
            IF (ionode) WRITE(UNIT=output_unit,FMT="(T3,A)")&
              "WARNING: S**2 computation does not yet treat fractional occupied orbitals"
          END IF
          CALL compute_s_square(mos=mos, matrix_s=matrix_s, s_square=s_square,&
               s_square_ideal=s_square_ideal,error=error)
          IF (ionode) WRITE(UNIT=output_unit,fmt='(T3,A,T60,2F10.5)') &
               "Ideal and single determinant S**2 : ",s_square_ideal,s_square
       ENDIF

       ! Print the total density (electronic + core charge)

       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%TOT_DENSITY_CUBE", error=error),cp_p_file)) THEN
          NULLIFY(rho_core,rho0_s_gs)
          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,rho_core=rho_core,&
               rho0_s_gs=rho0_s_gs,error=error)
          CALL pw_env_get(pw_env, error=error)
          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
               pw_pools=pw_pools,error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          IF (dft_control%qs_control%gapw) THEN
             CALL pw_transfer(rho0_s_gs%pw,wf_r%pw)
          ELSE
             CALL pw_transfer(rho_core%pw,wf_r%pw)
          END IF
          DO ispin=1,dft_control%nspins
             CALL pw_axpy(rho%rho_r(ispin)%pw,wf_r%pw)
          END DO
          filename = "TOTAL_DENSITY"
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%TOT_DENSITY_CUBE",&
               extension=".cube",middle_name=TRIM(filename),log_filename=.FALSE.,error=error)
          CALL pw_to_cube(wf_r%pw,unit_nr,ionode,"TOTAL DENSITY",&
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
             CALL pw_zero(rho_elec_rspace%pw)
             CALL pw_pool_create_pw(pool=auxbas_pw_pool,&
                                    pw=rho_elec_gspace%pw,&
                                    use_data=COMPLEXDATA1D,&
                                    in_space=RECIPROCALSPACE,&
                                    error=error)
             CALL pw_zero(rho_elec_gspace%pw)
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
             CALL pw_transfer(rho_elec_gspace%pw,rho_elec_rspace%pw,debug=.FALSE.)
             rho_total_rspace = pw_integrate_function(rho_elec_rspace%pw,isign=-1)/volume
             filename = "TOTAL_ELECTRON_DENSITY"
             unit_nr = cp_print_key_unit_nr(logger,input,"DFT%PRINT%E_DENSITY_CUBE",&
                                            extension=".cube",middle_name=TRIM(filename),&
                                            log_filename=.FALSE.,error=error)
             IF (ionode) THEN
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
              CALL pw_to_cube(rho_elec_rspace%pw,unit_nr,ionode,"TOTAL ELECTRON DENSITY",&
                   particles=particles,&
                   stride=section_get_ivals(dft_section,"PRINT%E_DENSITY_CUBE%STRIDE",error=error),&
                   error=error)
             CALL cp_print_key_finished_output(unit_nr,logger,input,&
                                               "DFT%PRINT%E_DENSITY_CUBE",&
                                               error=error)
             ! Print total spin density for spin-polarized systems
             IF (dft_control%nspins > 1) THEN
                CALL pw_zero(rho_elec_gspace%pw)
                CALL pw_zero(rho_elec_rspace%pw)
                CALL calculate_rhotot_elec_gspace(qs_env=qs_env,&
                                                  auxbas_pw_pool=auxbas_pw_pool,&
                                                  rhotot_elec_gspace=rho_elec_gspace,&
                                                  q_max=q_max,&
                                                  rho_hard=rho_hard,&
                                                  rho_soft=rho_soft,&
                                                  fsign=-1.0_dp,&
                                                  error=error)
                rho_total = rho_hard + rho_soft
                CALL pw_transfer(rho_elec_gspace%pw,rho_elec_rspace%pw,debug=.FALSE.)
                rho_total_rspace = pw_integrate_function(rho_elec_rspace%pw,isign=-1)/volume
                filename = "TOTAL_SPIN_DENSITY"
                unit_nr = cp_print_key_unit_nr(logger,input,"DFT%PRINT%E_DENSITY_CUBE",&
                                               extension=".cube",middle_name=TRIM(filename),&
                                               log_filename=.FALSE.,error=error)
                IF (ionode) THEN
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
                CALL pw_to_cube(rho_elec_rspace%pw,unit_nr,ionode,"TOTAL SPIN DENSITY",&
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
                CALL pw_copy(rho%rho_r(1)%pw,rho_elec_rspace%pw)
                CALL pw_axpy(rho%rho_r(2)%pw,rho_elec_rspace%pw)
                filename = "ELECTRON_DENSITY"
                unit_nr = cp_print_key_unit_nr(logger,input,"DFT%PRINT%E_DENSITY_CUBE",&
                                               extension=".cube",middle_name=TRIM(filename),&
                                               log_filename=.FALSE.,error=error)
                IF (ionode) THEN
                   INQUIRE (UNIT=unit_nr,NAME=filename)
                   WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A)")&
                     "The sum of alpha and beta density is written in cube file format to the file:",&
                     TRIM(filename)
                END IF
                CALL pw_to_cube(rho_elec_rspace%pw,unit_nr,ionode,"SUM OF ALPHA AND BETA DENSITY",&
                     particles=particles,stride=section_get_ivals(dft_section,"PRINT%E_DENSITY_CUBE%STRIDE",error=error),&
                     error=error)
                CALL cp_print_key_finished_output(unit_nr,logger,input,&
                                         "DFT%PRINT%E_DENSITY_CUBE",error=error)
                CALL pw_copy(rho%rho_r(1)%pw,rho_elec_rspace%pw)
                CALL pw_axpy(rho%rho_r(2)%pw,rho_elec_rspace%pw,alpha=-1.0_dp)
                filename = "SPIN_DENSITY"
                unit_nr = cp_print_key_unit_nr(logger,input,"DFT%PRINT%E_DENSITY_CUBE",&
                                               extension=".cube",middle_name=TRIM(filename),&
                                               log_filename=.FALSE.,error=error)
                IF (ionode) THEN
                   INQUIRE (UNIT=unit_nr,NAME=filename)
                   WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A)")&
                     "The spin density is written in cube file format to the file:",&
                     TRIM(filename)
                END IF
                CALL pw_to_cube(rho_elec_rspace%pw,unit_nr,ionode,"SPIN DENSITY",&
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
                                               log_filename=.FALSE.,error=error)
                IF (ionode) THEN
                   INQUIRE (UNIT=unit_nr,NAME=filename)
                   WRITE (UNIT=output_unit,FMT="(/,T2,A,/,/,T2,A)")&
                     "The electron density is written in cube file format to the file:",&
                     TRIM(filename)
                END IF
                CALL pw_to_cube(rho%rho_r(1)%pw,unit_nr,ionode,"ELECTRON DENSITY",&
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
          CALL pw_env_get(pw_env, error=error)
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%V_HARTREE_CUBE",&
               extension=".cube",middle_name="v_hartree",error=error)
          udvol = 1.0_dp/qs_env%ks_env%v_hartree_rspace%pw%pw_grid%dvol
          CALL pw_scale(qs_env%ks_env%v_hartree_rspace%pw,udvol)

          CALL pw_to_cube(qs_env%ks_env%v_hartree_rspace%pw,&
               unit_nr,unit_nr>0,"HARTREE POTENTIAL",particles=particles,&
               stride=section_get_ivals(dft_section,&
               "PRINT%V_HARTREE_CUBE%STRIDE",error=error),&
               error=error)
          CALL cp_print_key_finished_output(unit_nr,logger,input,&
               "DFT%PRINT%V_HARTREE_CUBE",error=error)
       ENDIF

       !    *** write the density matrix ***
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%AO_MATRICES/DENSITY",error=error),cp_p_file)) THEN
          iw = cp_print_key_unit_nr(logger,input,"DFT%PRINT%AO_MATRICES/DENSITY",&
               extension=".Log",error=error)
          DO ispin=1,dft_control%nspins
             CALL write_sparse_matrix(rho%rho_ao(ispin)%matrix,4,6,qs_env,&
                  para_env,output_unit=iw,error=error)
          END DO
          CALL cp_print_key_finished_output(iw,logger,input,&
               "DFT%PRINT%AO_MATRICES/DENSITY", error=error)
       END IF

       !    **** the Kohn-Sham matrix itself
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",error=error),cp_p_file)) THEN
          CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,&
               error=error,&
               calculate_forces=.FALSE.,just_energy=.FALSE.)
          CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.,error=error)
          iw = cp_print_key_unit_nr(logger,input,"DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",&
               extension=".Log",error=error)
          CALL Write_sparse_matrix(ks_rmpv(1)%matrix,4,6,Qs_env,para_env,output_unit=iw,error=error)
          CALL cp_print_key_finished_output(iw,logger,input,&
               "DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX", error=error)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE write_available_results

END MODULE
