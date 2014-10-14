!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Does all kind of post scf calculations for semi-empirical
!> \par History
!>      Started printing preliminary stuff for MO_CUBES and MO requires some
!>      more work to complete all other functionalities
!> \author Teodoro Laino (07.2008)
! *****************************************************************************
MODULE qs_scf_post_se
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cell_types,                      ONLY: cell_type
  USE cp_array_r_utils,                ONLY: cp_1d_r_p_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
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
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE input_constants,                 ONLY: ot_precond_full_all,&
                                             state_loc_all
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE molecular_states,                ONLY: construct_molecular_states
  USE molecule_types_new,              ONLY: molecule_type
  USE mulliken,                        ONLY: multipoles
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: evolt
  USE population_analyses,             ONLY: lowdin_population_analysis,&
                                             mulliken_population_analysis
  USE preconditioner_types,            ONLY: preconditioner_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_kind_types,                   ONLY: qs_kind_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_loc_methods,                  ONLY: qs_loc_driver
  USE qs_loc_molecules,                ONLY: wfc_to_molecule
  USE qs_loc_types,                    ONLY: localized_wfn_control_create,&
                                             qs_loc_env_create,&
                                             qs_loc_env_destroy,&
                                             qs_loc_env_new_type
  USE qs_loc_utils,                    ONLY: qs_loc_env_init,&
                                             read_loc_section,&
                                             set_loc_centers,&
                                             set_loc_wfn_lists
  USE qs_mo_io,                        ONLY: write_mo_set
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             set_mo_occupation
  USE qs_ot_eigensolver,               ONLY: ot_eigensolver
  USE qs_resp,                         ONLY: resp_fit
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_post_gpw,                 ONLY: qs_scf_post_occ_cubes,&
                                             qs_scf_post_unocc_cubes
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                             qs_subsys_type
  USE scf_control_types,               ONLY: scf_control_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  ! Global parameters
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_post_se'
  PUBLIC :: scf_post_calculation_se

CONTAINS

! *****************************************************************************
!> \brief collects possible post - scf calculations and prints info / computes properties.
!>        specific for Semi-empirical calculations
!> \param dft_section ...
!> \param scf_env the scf_env whose info should be written out
!> \param qs_env the qs_env in which the qs_env lives
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      07.2008 created [tlaino] - Splitted from qs_scf_post (general)
!> \author tlaino
!> \note
!>      this function changes mo_eigenvectors and mo_eigenvalues, depending on the print keys.
!>      In particular, MO_CUBES causes the MOs to be rotated to make them eigenstates of the KS
!>      matrix, and mo_eigenvalues is updated accordingly. This can, for unconverged wavefunctions,
!>      change afterwards slightly the forces (hence small numerical differences between MD
!>      with and without the debug print level). Ideally this should not happen...
! *****************************************************************************
  SUBROUTINE scf_post_calculation_se(dft_section, scf_env,qs_env,error)

    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_post_calculation_se', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, homo, ispin, n, n_mo(2), nao, nelectron, nhomo, nlumo, &
      nlumos, nmo, nmoloc(2), nmos, output_unit
    INTEGER, DIMENSION(:), POINTER           :: marked_states
    LOGICAL :: do_mo_cubes, do_wannier_cubes, explicit, failure, &
      my_localized_wfn, p_loc, p_molstates, p_wanstates, uniform_occ
    REAL(KIND=dp)                            :: maxocc
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: scenter
    TYPE(cp_1d_r_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_evals
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv, matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_deriv
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_orbs
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_localized
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(preconditioner_type), POINTER       :: local_preconditioner
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: wf_g, wf_r
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_loc_env_new_type), POINTER       :: qs_loc_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: input, loc_print_section, &
                                                localize_section, print_key

    CALL timeset(routineN,handle)

    ! Writes the data that is already available in qs_env
    para_env=>qs_env%para_env
    CALL write_available_results(qs_env,scf_env,error)

    failure=.FALSE.
    my_localized_wfn = .FALSE.
    NULLIFY(dft_control, pw_env, auxbas_pw_pool, pw_pools, mos, rho, &
         mo_coeff, ks_rmpv, matrix_s, qs_loc_env, scf_control, &
         unoccupied_orbs, mo_eigenvalues, unoccupied_evals, fm_struct_tmp, &
         molecule_set, subsys, particles, input, print_key)

    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    ! Here we start with data that needs a postprocessing...
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,dft_control=dft_control,molecule_set=molecule_set,&
            mos=mos,rho=rho,matrix_ks=ks_rmpv,scf_control=scf_control,matrix_s=matrix_s,&
            input=input, subsys=subsys, pw_env=pw_env, error=error)
       CALL qs_subsys_get(subsys,particles=particles,error=error)

       ! Compute Atomic Charges
       CALL qs_scf_post_charges(input, logger, qs_env, rho, matrix_s, mos, para_env, error)

       ! Moments of charge distribution
       CALL qs_scf_post_moments(input, logger, qs_env, output_unit, error)

       ! Determine if we need to computer properties using the localized centers
       localize_section     => section_vals_get_subs_vals(dft_section,"LOCALIZE",error=error)
       loc_print_section => section_vals_get_subs_vals(localize_section,"PRINT",error=error)
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
       CALL section_vals_val_get(localize_section,"_SECTION_PARAMETERS_",l_val=explicit,error=error)
       p_loc=p_loc.OR.explicit

       ! Print Molecular States
       print_key => section_vals_get_subs_vals(loc_print_section,"MOLECULAR_STATES",error=error)
       p_molstates=BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)
       p_loc=p_loc.OR.p_molstates
       ! Print Wannier States
       print_key => section_vals_get_subs_vals(loc_print_section,"WANNIER_STATES",error=error)
       p_wanstates=BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)
       p_loc=p_loc.OR.p_wanstates

       ! Check for cubes MOs and WANNIERS
       do_mo_cubes=BTEST(cp_print_key_should_output(logger%iter_info,dft_section,"PRINT%MO_CUBES",&
                   error=error),cp_p_file)
       do_wannier_cubes=BTEST(cp_print_key_should_output(logger%iter_info,loc_print_section,&
                        "WANNIER_CUBES",error=error),cp_p_file)
       nlumo=section_get_ival(dft_section,"PRINT%MO_CUBES%NLUMO",error=error)
       nhomo=section_get_ival(dft_section,"PRINT%MO_CUBES%NHOMO",error=error)

       ! Setup the grids needed to compute a wavefunction given a vector..
       IF (((do_mo_cubes.OR.do_wannier_cubes).AND.(nlumo/=0.OR.nhomo/=0.)).OR.p_loc) THEN
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, pw_pools=pw_pools,error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,wf_r%pw, use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,wf_g%pw, use_data = COMPLEXDATA1D,&
               in_space = RECIPROCALSPACE, error=error)
       END IF


       ! Makes the MOs eigenstates, computes eigenvalues, write cubes
       IF (do_mo_cubes .AND. (nhomo/=0)) THEN
          IF (dft_control%restricted) THEN
             !
             ! The issue is that one state is different from the others,
             ! so we certainly can't rotate all N+1 alpha states among themselves
             ! furthermore, alpha and beta orbitals do have the same shape,
             ! but not the same v_xc, therefore these seems to be no reason that
             ! they should have the same eigenvalues.
             ! So it looks like I need a convincing definition of an eigenstate in this case :-)
             !
             IF (output_unit>0) WRITE(output_unit,*) &
                " Unclear how we define MOs in the restricted case ... skipping"
          ELSE
             DO ispin=1,dft_control%nspins
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                                eigenvalues=mo_eigenvalues,homo=homo,nmo=nmo)
                IF (output_unit>0) WRITE(output_unit,*) " "
                IF (output_unit>0) WRITE(output_unit,*) " Eigenvalues of the occupied subspace spin ",ispin
                IF (homo .NE. nmo) THEN
                   IF (output_unit>0) WRITE(output_unit,*)" and ",nmo-homo," added MO eigenvalues"
                ENDIF
                IF (output_unit>0) WRITE(output_unit,*) "---------------------------------------------"
                ! Also rotate the OT derivs, since they are needed for force calculations
                IF (ASSOCIATED(qs_env%mo_derivs)) THEN
                   mo_coeff_deriv=>qs_env%mo_derivs(ispin)%matrix
                ELSE
                   mo_coeff_deriv=>NULL()
                ENDIF
                CALL calculate_subspace_eigenvalues(mo_coeff,ks_rmpv(ispin)%matrix,mo_eigenvalues, &
                  scr=output_unit, ionode=output_unit>0, do_rotation=.TRUE.,&
                  co_rotate_dbcsr=mo_coeff_deriv,error=error)

                CALL set_mo_occupation(mo_set=mos(ispin)%mo_set,smear=scf_control%smear,error=error)
                IF (output_unit>0) WRITE(output_unit,'(T2,A,F12.6)') "Fermi Energy [eV] :",mos(ispin)%mo_set%mu*evolt
                ! Prints the cube files of OCCUPIED ORBITALS
                CALL qs_scf_post_occ_cubes(input, dft_section, dft_control, logger, qs_env,&
                     mo_coeff, wf_g, wf_r, particles, homo, ispin, error=error)
             ENDDO
          ENDIF
       ENDIF

       ! Initialize the localization environment, needed e.g. for wannier functions and molecular states
       IF (p_loc) THEN
          CALL qs_loc_env_create(qs_loc_env,error=error)
          CALL localized_wfn_control_create(qs_loc_env%localized_wfn_control,error=error)
          CALL read_loc_section(qs_loc_env%localized_wfn_control,localize_section,&
               qs_loc_env%do_localize, error=error)
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
                  "LOCALIZE (located in: DFT%LOCALIZE) was not specified. Localization "//&
                  "will not be performed!"//&
CPSourceFileRef,&
                  only_ionode=.TRUE.)
          END IF
       END IF

       ! Gets localization info for the occupied orbs
       !  - Possibly gets wannier functions
       !  - Possibly gets molecular states
       IF (p_loc) THEN
          IF (dft_control%restricted) THEN
             IF (output_unit>0) WRITE(output_unit,*) &
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
                   IF (output_unit>0) WRITE(output_unit,*) " "
                   IF (output_unit>0) WRITE(output_unit,*) " Eigenvalues of the occupied subspace spin ",ispin
                   IF (homo .NE. nmo) THEN
                      IF (output_unit>0) WRITE(output_unit,*) " and ",nmo-homo," added MO eigenvalues"
                   ENDIF
                   IF (output_unit>0) WRITE(output_unit,*) "---------------------------------------------"
                   ! Also rotate the OT derivs, since they are needed for force calculations
                   IF (ASSOCIATED(qs_env%mo_derivs)) THEN
                      mo_coeff_deriv=>qs_env%mo_derivs(ispin)%matrix
                   ELSE
                      mo_coeff_deriv=>NULL()
                   ENDIF
                   CALL calculate_subspace_eigenvalues(mo_coeff,ks_rmpv(ispin)%matrix, &
                        evals_arg=mo_eigenvalues, scr=output_unit, &
                        ionode=output_unit>0, do_rotation=.TRUE.,error=error)
                END IF
                NULLIFY(mo_localized)
                CALL cp_fm_create(mo_localized,mo_coeff%matrix_struct,error=error)
                CALL cp_fm_to_fm(mo_coeff,mo_localized,error=error)
                IF (qs_loc_env%do_localize) THEN
                   ! Do the Real localization..
                   IF (output_unit>0) WRITE(output_unit,"(/,T2,A,I3)")&
                        "LOCALIZATION| Computing localization properties "//&
                        "for OCCUPIED ORBITALS. Spin:",ispin

                   scenter => qs_loc_env%localized_wfn_control%centers_set(ispin)%array
                   CALL qs_loc_driver(qs_env,qs_loc_env,localize_section,loc_print_section,&
                        myspin=ispin,ext_mo_coeff=mo_localized,error=error)

                   ! maps wfc to molecules, and compute the molecular dipoles if required
                   IF (BTEST(cp_print_key_should_output(logger%iter_info,loc_print_section,&
                       "MOLECULAR_DIPOLES",error=error),cp_p_file) .OR. p_molstates) THEN
                      CALL wfc_to_molecule(qs_env, qs_loc_env, loc_print_section, scenter, molecule_set,&
                           dft_control%nspins, error)
                   END IF

                   ! Compute the molecular states
                   IF ( p_molstates ) THEN
                      CALL construct_molecular_states(molecule_set, mo_localized, mo_coeff, mo_eigenvalues,&
                           ks_rmpv(ispin)%matrix, matrix_s(1)%matrix, qs_env,  wf_r, wf_g,&
                           loc_print_section=loc_print_section, particles=particles, tag="HOMO",&
                           marked_states=marked_states, error=error)
                      IF(ASSOCIATED(marked_states)) DEALLOCATE(marked_states)
                   ENDIF
                   IF ( p_wanstates ) THEN
                      IF(ASSOCIATED(marked_states)) DEALLOCATE(marked_states)
                   ENDIF
                END IF
                CALL cp_fm_release(mo_localized,error=error)
             ENDDO
             ! Print Total Dipole if the localization has been performed
             IF (qs_loc_env%do_localize) THEN
                CALL qs_scf_post_loc_dip(input, dft_control, qs_loc_env, logger, qs_env, error)
             END IF
          ENDIF
       ENDIF

       ! Gets the lumos, and eigenvalues for the lumos
       IF (do_mo_cubes .AND. nlumo .NE. 0 ) THEN
          ALLOCATE(unoccupied_orbs(dft_control%nspins))
          ALLOCATE(unoccupied_evals(dft_control%nspins))
          DO ispin=1,dft_control%nspins
             NULLIFY(unoccupied_orbs(ispin)%matrix)
             NULLIFY(unoccupied_evals(ispin)%array)
             ! Always write eigenvalues
             IF (output_unit>0) WRITE(output_unit,*) " "
             IF (output_unit>0) WRITE(output_unit,*) " Lowest Eigenvalues of the unoccupied subspace spin ",ispin
             IF (output_unit>0) WRITE(output_unit,FMT='(1X,A)')&
                                      "-----------------------------------------------------"
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
             NULLIFY(local_preconditioner)
             IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
                local_preconditioner=>scf_env%ot_preconditioner(1)%preconditioner
                ! this one can for sure not be right (as it has to match a given C0)
                IF (local_preconditioner%in_use == ot_precond_full_all) THEN
                    NULLIFY(local_preconditioner)
                ENDIF
             ENDIF

             CALL ot_eigensolver(matrix_h=ks_rmpv(ispin)%matrix,&
                  matrix_s=matrix_s(1)%matrix, &
                  matrix_c_fm=unoccupied_orbs(ispin)%matrix, &
                  matrix_orthogonal_space_fm=mo_coeff, &
                  eps_gradient=scf_control%eps_lumos, &
                  preconditioner=local_preconditioner, &
                  iter_max=scf_control%max_iter_lumos,&
                  size_ortho_space=homo,error=error)

             CALL calculate_subspace_eigenvalues(unoccupied_orbs(ispin)%matrix,ks_rmpv(ispin)%matrix,&
                  unoccupied_evals(ispin)%array, scr=output_unit, &
                  ionode=output_unit>0,error=error)

             ! Prints the cube files of UNOCCUPIED ORBITALS
             CALL qs_scf_post_unocc_cubes(input, dft_section, dft_control, logger, qs_env,&
                  unoccupied_orbs(ispin)%matrix, wf_g, wf_r, particles, nlumos, homo, ispin, error=error)
          ENDDO

       ENDIF


       ! Deallocate grids needed to compute wavefunctions
       IF ( ( ( do_mo_cubes .OR. do_wannier_cubes ).AND. (nlumo /= 0 .OR. nhomo /= 0 )) .OR. p_loc ) THEN
          CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_r%pw, error=error)
          CALL pw_pool_give_back_pw(auxbas_pw_pool,wf_g%pw, error=error)
       END IF

       ! Destroy the localization environment
       IF (p_loc) THEN
          CALL qs_loc_env_destroy(qs_loc_env, error=error)
       ENDIF

       ! Compute the optical conductivity (needs homos and lumos)
       CALL qs_scf_post_optc(input, dft_section, dft_control, logger, qs_env,&
            unoccupied_orbs, unoccupied_evals, output_unit, error)

       ! This is just a deallocation for printing MO_CUBES
       IF (do_mo_cubes .AND. nlumo /= 0 ) THEN
          DO ispin=1,dft_control%nspins
             DEALLOCATE(unoccupied_evals(ispin)%array)
             CALL cp_fm_release(unoccupied_orbs(ispin)%matrix,error=error)
          ENDDO
          DEALLOCATE(unoccupied_evals)
          DEALLOCATE(unoccupied_orbs)
       ENDIF

       ! Print coherent X-ray diffraction spectrum
       CALL qs_scf_post_xray(input, dft_section, logger, qs_env, output_unit, error)

       ! Calculation of Electric Field Gradients
       CALL qs_scf_post_efg(input, logger, error)

       ! Calculation of ET
       CALL qs_scf_post_et(input, qs_env, dft_control, error)

       ! Calculation of EPR Hyperfine Coupling Tensors
       CALL qs_scf_post_epr(input, logger, qs_env, error)

       ! Calculation of properties needed for BASIS_MOLOPT optimizations
       CALL qs_scf_post_molopt(input, logger, qs_env, error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE scf_post_calculation_se

! *****************************************************************************
!> \brief Computes and prints Electric Field Gradient
!> \param input ...
!> \param logger ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_efg(input, logger, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_efg', &
      routineP = moduleN//':'//routineN

    TYPE(section_vals_type), POINTER         :: print_key

    print_key => section_vals_get_subs_vals(section_vals=input,&
                   subsection_name="DFT%PRINT%ELECTRIC_FIELD_GRADIENT",&
                   error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),&
              cp_p_file)) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="ELECTRIC_FIELD_GRADIENT not implemented for Semi-Empirical calculations!!", &
            error=error, error_level=cp_warning_level)
    END IF

  END SUBROUTINE qs_scf_post_efg

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

    TYPE(section_vals_type), POINTER         :: print_key

    print_key => section_vals_get_subs_vals(section_vals=input,&
         subsection_name="DFT%PRINT%MOMENTS",error=error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="ELECTRIC MOMENTS not implemented for Semi-Empirical calculations!!", &
            error=error, error_level=cp_warning_level)
    END IF

  END SUBROUTINE qs_scf_post_moments

! *****************************************************************************
!> \brief Computes and Prints Atomic Charges with several methods
!> \param input ...
!> \param logger ...
!> \param qs_env the qs_env in which the qs_env lives
!> \param rho ...
!> \param matrix_s ...
!> \param mos ...
!> \param para_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE qs_scf_post_charges(input, logger, qs_env, rho, matrix_s, mos, para_env, error)
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_scf_post_charges', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: print_level, unit_nr
    LOGICAL                                  :: print_it, scp_nddo
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_type), POINTER             :: pscp
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(section_vals_type), POINTER         :: print_key

    NULLIFY(particle_set)
    CALL get_qs_env(qs_env=qs_env,&
                    atomic_kind_set=atomic_kind_set,&
                    qs_kind_set=qs_kind_set,&
                    local_particles=local_particles,&
                    particle_set=particle_set,&
                    scp_env=scp_env,&
                    dft_control=dft_control,&
                    error=error)
    scp_nddo = dft_control%qs_control%se_control%scp
    IF (scp_nddo ) &
    CALL get_scp_env(scp_env=scp_env, pscp=pscp, error=error)
    ! Compute the mulliken charges
    print_key => section_vals_get_subs_vals(input,"DFT%PRINT%MULLIKEN", error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="MULLIKEN charges not validated for Semi-Empirical calculations!!", &
            error=error, error_level=cp_warning_level)
       unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%MULLIKEN",extension=".mulliken",&
                                    middle_name="",log_filename=.FALSE.,error=error)
       print_level = 1
       CALL section_vals_val_get(print_key,"PRINT_GOP",l_val=print_it,error=error)
       IF (print_it) print_level = 2
       CALL section_vals_val_get(print_key,"PRINT_ALL",l_val=print_it,error=error)
       IF (print_it) print_level = 3
       CALL mulliken_population_analysis(qs_env,unit_nr,print_level,error)
       IF ( scp_nddo ) &
       CALL multipoles(atomic_kind_set,qs_kind_set,particle_set,&
         rho%rho_ao,pscp,matrix_s(1)%matrix,para_env,unit_nr,"DFT%PRINT%MULLIKEN",error=error)

       CALL cp_print_key_finished_output(unit_nr, logger,input,"DFT%PRINT%MULLIKEN",error=error)
    END IF

    ! Compute the Lowdin charges
    print_key => section_vals_get_subs_vals(input,"DFT%PRINT%LOWDIN", error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
       CALL cp_unimplemented_error(fromWhere=routineP,&
                                   message="Lowdin charges not validated for semi-empirical calculations!", &
                                   error=error, error_level=cp_warning_level)
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
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="FIT CHARGES not implemented for Semi-Empirical calculations!!", &
            error=error, error_level=cp_warning_level)
    END IF

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

    LOGICAL                                  :: failure, first_time

    failure = .FALSE.
    IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
         "DFT%PRINT%LOCALIZATION%TOTAL_DIPOLE",first_time=first_time,error=error),cp_p_file)) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="TOTAL DIPOLE  not implemented for Semi-Empirical calculations!!", &
            error=error, error_level=cp_warning_level)
    END IF

  END SUBROUTINE qs_scf_post_loc_dip

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

    IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
         "DFT%PRINT%OPTICAL_CONDUCTIVITY",error=error),cp_p_file)) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="OPTICAL_CONDUCTIVITY  not implemented for Semi-Empirical calculations!!", &
            error=error, error_level=cp_warning_level)
    ENDIF
  END SUBROUTINE qs_scf_post_optc

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
  SUBROUTINE qs_scf_post_xray(input,dft_section,logger,qs_env,&
                              output_unit,error)

    TYPE(section_vals_type), POINTER         :: input, dft_section
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_environment_type), POINTER       :: qs_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_post_xray', &
      routineP = moduleN//':'//routineN

    TYPE(section_vals_type), POINTER         :: print_key

    print_key => section_vals_get_subs_vals(section_vals=input,&
                                            subsection_name="DFT%PRINT%XRAY_DIFFRACTION_SPECTRUM",&
                                            error=error)

    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),cp_p_file)) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="XRAY_DIFFRACTION_SPECTRUM  not implemented for Semi-Empirical calculations!!", &
            error=error, error_level=cp_warning_level)
    END IF

  END SUBROUTINE qs_scf_post_xray

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

    TYPE(section_vals_type), POINTER         :: print_key

    print_key => section_vals_get_subs_vals(section_vals=input,&
                   subsection_name="DFT%PRINT%HYPERFINE_COUPLING_TENSOR",&
                   error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),&
              cp_p_file)) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="HYPERFINE_COUPLING_TENSOR  not implemented for Semi-Empirical calculations!!", &
            error=error, error_level=cp_warning_level)
    END IF

  END SUBROUTINE qs_scf_post_epr

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

    TYPE(section_vals_type), POINTER         :: print_key

    print_key => section_vals_get_subs_vals(section_vals=input,&
                   subsection_name="DFT%PRINT%BASIS_MOLOPT_QUANTITIES",&
                   error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,error=error),&
              cp_p_file)) THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="BASIS_MOLOPT_QUANTITIES  not implemented for Semi-Empirical calculations!!", &
            error=error, error_level=cp_warning_level)
    END IF

  END SUBROUTINE qs_scf_post_molopt

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

    LOGICAL                                  :: do_et
    TYPE(section_vals_type), POINTER         :: et_section

    do_et=.FALSE.
    et_section =>  section_vals_get_subs_vals(input,"PROPERTIES%ET_COUPLING",&
                                                  error=error)
    CALL section_vals_get(et_section,explicit=do_et,error=error)
    IF(do_et)THEN
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="ET_COUPLING  not implemented for Semi-Empirical calculations!!", &
            error=error, error_level=cp_warning_level)
    END IF

  END SUBROUTINE qs_scf_post_et

! *****************************************************************************
!> \brief Write QS results always available (if switched on through the print_keys)
!> \param qs_env the qs_env in which the qs_env lives
!> \param scf_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE write_available_results(qs_env,scf_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_available_results', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, iw, output_unit
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_kind_type), DIMENSION(:), &
      POINTER                                :: qs_kind_set
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(qs_subsys_type), POINTER            :: subsys
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(section_vals_type), POINTER         :: dft_section, input

    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(cell, dft_control, mos, atomic_kind_set, particle_set, rho, &
         ks_rmpv, scf_control, dft_section, molecule_set, input, &
         particles, subsys)
    para_env=>qs_env%para_env
    logger => cp_error_get_logger(error)
    output_unit= cp_logger_get_default_io_unit(logger)

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,dft_control=dft_control,molecule_set=molecule_set, &
            mos=mos,atomic_kind_set=atomic_kind_set,&
            qs_kind_set=qs_kind_set,particle_set=particle_set,&
            rho=rho,matrix_ks=ks_rmpv,scf_control=scf_control,&
            input=input,cell=cell,subsys=subsys,error=error)
       CALL qs_subsys_get(subsys,particles=particles,error=error)

       ! *** if the dft_section tells you to do so, write last wavefunction to screen
       dft_section => section_vals_get_subs_vals(input,"DFT",error=error)
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
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="PDOS  not implemented for Semi-Empirical calculations!!", &
               error=error, error_level=cp_warning_level)
       ENDIF

       ! Print the total density (electronic + core charge)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%TOT_DENSITY_CUBE", error=error),cp_p_file)) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="TOT_DENSITY_CUBE  not implemented for Semi-Empirical calculations!!", &
               error=error, error_level=cp_warning_level)
       END IF

       ! Write cube file with electron density
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
                 "DFT%PRINT%E_DENSITY_CUBE",error=error),cp_p_file)) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="E_DENSITY_CUBE  not implemented for Semi-Empirical calculations!!", &
               error=error, error_level=cp_warning_level)
       END IF ! print key

       ! Print the hartree potential
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%V_HARTREE_CUBE",error=error),cp_p_file)) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="V_HARTREE_CUBE  not implemented for Semi-Empirical calculations!!", &
               error=error, error_level=cp_warning_level)
       ENDIF

       !    *** write the density matrix ***
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%AO_MATRICES/DENSITY",error=error),cp_p_file)) THEN
          iw = cp_print_key_unit_nr(logger,input,"DFT%PRINT%AO_MATRICES/DENSITY",&
               extension=".Log",error=error)
          DO ispin=1,dft_control%nspins
             CALL cp_dbcsr_write_sparse_matrix(rho%rho_ao(ispin)%matrix,4,6,qs_env,&
                  para_env,output_unit=iw,error=error)
          END DO
          CALL cp_print_key_finished_output(iw,logger,input,&
               "DFT%PRINT%AO_MATRICES/DENSITY", error=error)
       END IF

       !    **** the Kohn-Sham matrix itself
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",error=error),cp_p_file)) THEN
          CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,&
               error=error,calculate_forces=.FALSE.,just_energy=.FALSE.)
          CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.,error=error)
          iw = cp_print_key_unit_nr(logger,input,"DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",&
               extension=".Log",error=error)
          CALL cp_dbcsr_write_sparse_matrix(ks_rmpv(1)%matrix,4,6,Qs_env,para_env,output_unit=iw,error=error)
          CALL cp_print_key_finished_output(iw,logger,input,&
               "DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX", error=error)
       END IF
    END IF
    CALL timestop(handle)
  END SUBROUTINE write_available_results

END MODULE qs_scf_post_se
