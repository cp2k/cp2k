!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2003  CP2K developers group
!-----------------------------------------------------------------------------!


!!****** cp2k/qs_scf_post [1.0] *
!!
!!   NAME
!!     qs_scf_post
!!
!!   FUNCTION
!!     Does all kind of post scf calculations
!!
!!   AUTHOR
!!     Joost VandeVondele (10.2003)
!!
!!   MODIFICATION HISTORY
!!     Started as a copy from the relevant part of qs_scf
!!
!!   SOURCE
!******************************************************************************

MODULE qs_scf_post

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE coefficient_types,               ONLY: coeff_copy,&
                                             coeff_reduce,&
                                             coeff_sumup,&
                                             coeff_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             tddfpt_control_type
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
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
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_iter_string,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr,&
                                             cp_p_file
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE cp_rs_pool_types,                ONLY: cp_rs_pool_type,&
                                             rs_pool_create_rs,&
                                             rs_pool_give_back_rs
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_type
  USE global_types,                    ONLY: WANNIER_NONE,&
                                             global_environment_type
  USE input_constants,                 ONLY: ot_precond_full_all
  USE input_section_types,             ONLY: section_vals_type, &
                                             section_vals_val_get,&
                                             section_vals_get_subs_vals, &
                                             section_vals_get_subs_vals2, &
                                             section_get_rval, &
                                             section_get_ival, &
                                             section_vals_get
  USE cp_ddapc,                        ONLY: get_ddapc
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE level_shift,                     ONLY: level_shift_testing
  USE message_passing,                 ONLY: mp_sum
  USE molecular_states,                ONLY: construct_molecular_states
  USE molecule_types_new,              ONLY: molecule_type
  USE mulliken,                        ONLY: mulliken_charges
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom
  USE preconditioner,                  ONLY: preconditioner_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_coeff,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_sumup,&
                                             pw_transfer,&
                                             pw_type
  USE qs_collocate_density,            ONLY: calculate_total_abs_rho,&
                                             calculate_wavefunction
  USE qs_conductivity,                 ONLY: optical_conductivity
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_loc_control,                  ONLY: localized_wfn_control_type,&
                                             localized_wfn_control_release,&
                                             localized_wfn_control_retain
  USE qs_localization,                 ONLY: qs_loc_compute_properties,&
                                             qs_loc_destroy,&
                                             qs_loc_env_type,&
                                             qs_loc_init
  USE qs_loc_methods,                  ONLY: qs_loc_driver
  USE qs_loc_types       ,             ONLY: qs_loc_env_create,&
                                             qs_loc_env_destroy,& 
                                             qs_loc_env_new_type
  USE qs_loc_utils,                    ONLY: qs_loc_env_init
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             write_mo_set
  USE qs_ot_eigensolver,               ONLY: ot_eigensolver
  USE qs_overlap,                      ONLY: write_sparse_matrix
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE realspace_grid_types,            ONLY: rs_pw_to_cube, & 
                                             realspace_grid_type
  USE scf_control_types,               ONLY: scf_control_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE sparse_matrix_types,             ONLY: real_matrix_p_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xray_diffraction,                ONLY: xray_diffraction_spectrum
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "qs_scf_post"

  PUBLIC :: scf_post_calculation

! *****************************************************************************

CONTAINS

!!****f* qs_scf/scf_post_calculation [1.0] *
!!
!!   NAME
!!     scf_post_calculation
!!
!!   FUNCTION
!!     collects possible post - scf calculations and prints info / computes properties.
!!
!!   NOTES
!!     globenv should be eliminated
!!
!!   ARGUMENTS
!!     - scf_env: the scf_env whose info should be written out
!!     - qs_env: the qs_env in which the qs_env lives
!!     - globenv: the global environement (should be removed)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!     10.2004 moved here from qs_scf [Joost VandeVondele]
!!             started splitting out different subroutines
!!
!!*** **********************************************************************
  SUBROUTINE scf_post_calculation(dft_section, scf_env,qs_env,globenv,error)

    TYPE(section_vals_type), POINTER         :: dft_section
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type),POINTER    :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_post_calc', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, ionode
    LOGICAL                                  :: p_loc,p_molstates,p_cubes
    LOGICAL                                  :: my_localized_wfn, was_present
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type)                         :: wf_g, wf_r
    TYPE(coeff_type), POINTER                :: rho_core
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: unoccupied_orbs
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_localized
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(localized_wfn_control_type), POINTER :: my_localized_wfn_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(particle_list_type), pointer :: particles
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(qs_loc_env_new_type), POINTER           :: qs_loc_env_new
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(tddfpt_control_type), POINTER       :: tddfpt_control

    TYPE evals_p_type
         REAL(KIND = dp), DIMENSION(:), POINTER :: evals
    END TYPE
    TYPE (evals_p_type), DIMENSION(:), POINTER :: unoccupied_evals
    INTEGER ::  handle, ispin, istat, n, nmo, output_unit
    CHARACTER(LEN=200) :: file_name
    INTEGER            :: file_unit,ivector,nao,nlumos, homo, unit_nr
    LOGICAL            :: homoEQnmo
    CHARACTER(LEN=80)  :: filename
    CHARACTER(LEN=80) :: title
    CHARACTER(LEN=default_string_length) :: mulliken_title
    REAL(KIND = dp)           :: q_max,volume
    REAL(KIND = dp), POINTER, DIMENSION(:,:) :: center
    TYPE(section_vals_type), POINTER :: print_key, input
    TYPE(cp_rs_pool_type), POINTER           :: auxbas_rs_pool
    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(preconditioner_type), POINTER       :: preconditioner
    TYPE(cp_subsystem_type), pointer :: subsys
 
    CALL timeset(routineN,"I",'',handle)

! writes the data that is already available in qs_env 
    CALL write_available_results(qs_env,globenv,error)

    failure=.FALSE.
    my_localized_wfn = .FALSE.
    NULLIFY(cell,dft_control,pw_env,auxbas_pw_pool,pw_pools,&
            mos,atomic_kind_set,particle_set, rho,mo_coeff,&
            ks_rmpv,matrix_s,qs_loc_env,scf_control,unoccupied_orbs,&
            mo_eigenvalues,unoccupied_evals,fm_struct_tmp,molecule_set,&
            center,rho_core,tddfpt_control,subsys,particles,input,print_key)
    
    NULLIFY(qs_loc_env_new,my_localized_wfn_control)

    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) THEN
       output_unit= cp_logger_get_default_unit_nr(logger)
    ELSE
       output_unit = -1
    END IF

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,dft_control=dft_control,molecule_set=molecule_set, &
            mos=mos,atomic_kind_set=atomic_kind_set,particle_set=particle_set,&
            rho=rho,matrix_ks=ks_rmpv,scf_control=scf_control,matrix_s=matrix_s, &
            input=input, subsys=subsys,&
            cell=cell,rho_core=rho_core,tddfpt_control=tddfpt_control,error=error)
       CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
       CALL cp_subsys_get(subsys,particles=particles)

!    *** compute the mulliken charges ***
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%MULLIKEN",error=error),cp_p_file)) THEN
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%MULLIKEN",&
               extension=".mulliken",middle_name="",log_filename=.FALSE.,error=error)
          mulliken_title="  Mulliken charges : iteration level: " &
                                //cp_iter_string(logger%iter_info, print_key)
          CALL mulliken_charges(rho%rho_ao,matrix_s(1)%matrix,globenv%para_env,particle_set, &
                                  unit_nr,unit_nr>0,mulliken_title)

          CALL cp_print_key_finished_output(unit_nr, logger,input,&
               "DFT%PRINT%MULLIKEN", error=error)
       END IF

!    *** Fit charges with the Bloechl scheme
       print_key => section_vals_get_subs_vals(input,"PROPERTIES%FIT_CHARGE", error=error)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,&
           error=error),cp_p_file)) THEN
          CALL get_ddapc(qs_env,.FALSE.,.FALSE.,print_key,just_charge=.TRUE.,error=error)
       END IF

!    ***  Print Localized orbitals
       print_key => section_vals_get_subs_vals(dft_section,"PRINT%LOCALIZATION",error=error)
       p_loc=BTEST(cp_print_key_should_output(logger%iter_info,print_key,&
           error=error),cp_p_file)
!    ***  Print Molecular States
       print_key => section_vals_get_subs_vals(dft_section,"PRINT%MOLECULAR_STATES",error=error)
       p_molstates=BTEST(cp_print_key_should_output(logger%iter_info,print_key,&
           error=error),cp_p_file)
!    *** Print Cubes     
       print_key => section_vals_get_subs_vals(dft_section,"PRINT%MO_CUBES",error=error)
       p_cubes=BTEST(cp_print_key_should_output(logger%iter_info,print_key,&
            error=error),cp_p_file)

!    *** setup the grids needed to compute a wavefunction given a vector
       IF ( (p_cubes .AND. &
             ( section_get_ival(dft_section,"PRINT%MO_CUBES%NLUMO") /= 0 .OR. &
               section_get_ival(dft_section,"PRINT%MO_CUBES%NHOMO") /=0 ) ).OR. &
            p_loc ) THEN
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
               pw_pools=pw_pools,auxbas_rs_pool=auxbas_rs_pool)
          CALL rs_pool_create_rs(auxbas_rs_pool,rs)
          CALL pw_pool_init_coeff(auxbas_pw_pool,wf_r,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          CALL pw_pool_init_coeff(auxbas_pw_pool,wf_g,&
               use_data = COMPLEXDATA1D,&
               in_space = RECIPROCALSPACE, error=error)
       END IF

! *** makes the MOs eigenstates, computes eigenvalues, write cubes 
       IF (p_cubes .AND. section_get_ival(dft_section,"PRINT%MO_CUBES%NHOMO") /=0) THEN
          DO ispin=1,dft_control%nspins
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                             eigenvalues=mo_eigenvalues,homo=homo,nmo=nmo)
             IF (ionode) WRITE(output_unit,*) " "
             IF (ionode) WRITE(output_unit,*) " Eigenvalues of the occupied subspace spin ",ispin
             IF (homo .NE. nmo) THEN
                IF (ionode) WRITE(output_unit,*) &
                    " and ",nmo-homo," added MO eigenvalues" 
             ENDIF
             IF (ionode) WRITE(output_unit,*) " ---------------------------------------------"
             CALL calculate_subspace_eigenvalues(mo_coeff,&
                  ks_rmpv(ispin)%matrix,mo_eigenvalues, &
                  para_env=qs_env%para_env,scr=output_unit, &
                  ionode=ionode, do_rotation=.TRUE.)

             IF (p_cubes) THEN
                DO ivector=MAX(1,&
                   homo-section_get_ival(dft_section,"PRINT%MO_CUBES%NHOMO") +1),homo
                   CALL calculate_wavefunction(mo_coeff,ivector,wf_r, &
                        wf_g,&
                        qs_env,error=error)

                   WRITE(filename,'(a13,I5.5,a1,I1.1,a5)') "WAVEFUNCTION_",ivector,"_",ispin,".cube"
                   IF (ionode) OPEN(UNIT=21,FILE=filename)
                   WRITE(title,*) "WAVEFUNCTION ",ivector," spin ",ispin," i.e. HOMO - ",ivector-homo
                   CALL rs_pw_to_cube(wf_r%pw,21,ionode,title,rs,particles=particles,&
                        stride=section_get_ival(dft_section,"PRINT%MO_CUBES%STRIDE"))
                   IF (ionode) CLOSE(21)
                ENDDO
             END IF
          ENDDO
       ENDIF

! *** init the localization environment, needed e.g. for wannier functions and molecular states
       IF(dft_control%qs_control%localized_wfn) THEN

          my_localized_wfn_control =>  dft_control%qs_control%localized_wfn_control
          CALL localized_wfn_control_retain(my_localized_wfn_control,error=error)
          CPPrecondition(ASSOCIATED(my_localized_wfn_control),cp_failure_level,routineP,error,failure)
!     this should be done only in the at selected steps of the iter_level we are in (MD or..) 
!     for this I need to know the niter and the iter_level to refer to
          my_localized_wfn = (.NOT. failure .AND. .TRUE.)
          IF(my_localized_wfn) THEN
             CALL qs_loc_env_create(qs_loc_env_new,error=error)
             CALL qs_loc_env_init(qs_loc_env_new,my_localized_wfn_control,qs_env,error=error)
          END IF
          CALL localized_wfn_control_release(my_localized_wfn_control,error=error)

       ELSEIF (p_loc)THEN
          NULLIFY(qs_loc_env)
          ALLOCATE(qs_loc_env)
          CALL qs_loc_init(qs_loc_env,qs_env)
          CALL section_vals_val_get(input,"DFT%PRINT%LOCALIZATION%LOC_TYPE",&
               i_val=qs_loc_env%localization_method,error=error)
          qs_loc_env%wannier = .NOT. (qs_loc_env % localization_method .EQ. WANNIER_NONE)
          IF (p_molstates) qs_loc_env%molecular_states = .TRUE.
       ENDIF

! *** gets localization info for the occupied orbs
! *** possibly gets wannier functions
! *** possibly gets molecular states

       IF (my_localized_wfn) THEN

          CALL qs_loc_driver(qs_env,qs_loc_env_new,ionode,error=error)
 
       ELSEIF (p_loc) THEN
          DO ispin=1,dft_control%nspins
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                                              eigenvalues=mo_eigenvalues,homo=homo,nmo=nmo)
             ! just get eigenstates
             IF (section_get_ival(dft_section,"PRINT%MO_CUBES%NHOMO")==0) THEN
                IF (ionode) WRITE(output_unit,*) " "
                IF (ionode) WRITE(output_unit,*) " Eigenvalues of the occupied subspace spin ",ispin
                IF (homo .NE. nmo) THEN
                   IF (ionode) WRITE(output_unit,*) &
                       " and ",nmo-homo," added MO eigenvalues" 
                ENDIF
                IF (ionode) WRITE(output_unit,*) " ---------------------------------------------"
                CALL calculate_subspace_eigenvalues(mo_coeff,ks_rmpv(ispin)%matrix, &
                     evals_arg=mo_eigenvalues, para_env=qs_env%para_env,scr=output_unit, &
                     ionode=ionode, do_rotation=.TRUE.)
             ENDIF
             NULLIFY(mo_localized)
             CALL cp_fm_create(mo_localized,mo_coeff%matrix_struct)
             ALLOCATE(center(mo_localized%matrix_struct%ncol_global,6)) 
             CALL cp_fm_to_fm(mo_coeff,mo_localized)
             ! do localization computation
             WRITE(filename,'(a13,I1.1)') "HOMO_CENTERS_",ispin
             IF (ionode) WRITE(output_unit,*) ""
             IF (ionode) WRITE(output_unit,*) "Computing localization properties"
             IF (ionode) OPEN(UNIT=21,FILE=filename)
             CALL qs_loc_compute_properties(qs_loc_env,molecule_set,mo_localized,center,ionode,21)
             IF (ionode) CLOSE(21)
             IF (ionode) WRITE(output_unit,*) "Results in ",filename
             ! write the wannier functions to cubes
             IF (p_cubes.AND. .NOT. p_molstates ) THEN
                ! obviously we need a smarter scheme to select them. 
                DO ivector=1,mo_localized%matrix_struct%ncol_global   
                   CALL calculate_wavefunction(mo_localized,ivector,wf_r, &
                        wf_g,&
                        qs_env,error=error)

                   WRITE(filename,'(a7,I5.5,a1,I1.1,a5)') "WANNIER",ivector,"_",ispin,".cube"
                   IF (ionode) OPEN(UNIT=21,FILE=filename)
                   ! title name
                   WRITE(title,'(A7,I5.5,A1,I1.1,A1,6G10.4)') "WANNIER",ivector,"_",ispin," ",&
                     center(ivector,1:3)*angstrom, &
                     center(ivector,4:6)*angstrom
                   CALL rs_pw_to_cube(wf_r%pw,21,ionode,title,rs,&
                        particles=particles,&
                        stride=section_get_ival(dft_section,"PRINT%MO_CUBES%STRIDE"))
                   IF (ionode) CLOSE(21)
                ENDDO
             END IF
             ! compute the molecular states
             IF ( p_molstates ) THEN
                  CALL construct_molecular_states ( molecule_set, mo_localized, &
                                                    mo_coeff, mo_eigenvalues, &
                                                    ks_rmpv ( ispin ) % matrix, &
                                                    matrix_s(1)%matrix, &
                                                    globenv, qs_env,    &
                                                    wf_r, wf_g, error)
             ENDIF
             CALL cp_fm_release(mo_localized)

             IF (.FALSE.) &
             CALL level_shift_testing(atomic_kind_set, particle_set, cell, center, &
                                      pw_env, rho, rho_core)

             DEALLOCATE(center)
          ENDDO
       ENDIF

! gets the lumos, and eigenvalues for the lumos 
       IF ((p_cubes .OR.dft_control%do_tddfpt_calculation ) .AND. &
           section_get_ival(dft_section,"PRINT%MO_CUBES%NLUMO") /= 0) THEN
          ALLOCATE(unoccupied_orbs(dft_control%nspins))
          ALLOCATE(unoccupied_evals(dft_control%nspins))
          DO ispin=1,dft_control%nspins
             NULLIFY(unoccupied_orbs(ispin)%matrix)
             NULLIFY(unoccupied_evals(ispin)%evals)
! write eigenvalues and a set of cubes
             IF (ionode) WRITE(output_unit,*) " "
             IF (ionode) WRITE(output_unit,*) " Lowest Eigenvalues of the unoccupied subspace spin ",ispin
             IF (ionode) WRITE(output_unit,*) " -----------------------------------------------------"
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,homo=homo,nao=nao,nmo=nmo)
             CALL cp_fm_get_info(mo_coeff, nrow_global=n)
             nlumos=MAX(1,MIN(section_get_ival(dft_section,"PRINT%MO_CUBES%NLUMO"),nao-homo))
             ALLOCATE(unoccupied_evals(ispin)%evals(nlumos))
             CALL cp_fm_struct_create(fm_struct_tmp,para_env=qs_env%para_env,context=qs_env%blacs_env, &
                                                    nrow_global=n,ncol_global=nlumos)
             CALL cp_fm_create(unoccupied_orbs(ispin)%matrix, fm_struct_tmp,name="lumos")
             CALL cp_fm_struct_release(fm_struct_tmp)
             CALL cp_fm_init_random(unoccupied_orbs(ispin)%matrix,nlumos)

             ! the full_all preconditioner makes not much sense for lumos search
             NULLIFY(preconditioner)
             IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
                preconditioner=>scf_env%ot_preconditioner(1)%preconditioner
                IF (scf_env%qs_ot_env(1)%settings%preconditioner_type == ot_precond_full_all) THEN
                    NULLIFY(preconditioner)
                ENDIF
             ENDIF 

             CALL ot_eigensolver(matrix_h=ks_rmpv(ispin)%matrix,matrix_s=matrix_s(1)%matrix, &
                  matrix_c=unoccupied_orbs(ispin)%matrix, &
                  matrix_orthogonal_space=mo_coeff, &
                  eps_gradient=scf_control%eps_lumos, &
                  preconditioner=preconditioner, &
                  iter_max=scf_control%max_iter_lumos,globenv=globenv,&
                  size_ortho_space=homo)

             CALL calculate_subspace_eigenvalues(unoccupied_orbs(ispin)%matrix,ks_rmpv(ispin)%matrix,&
                  unoccupied_evals(ispin)%evals, para_env=qs_env%para_env,scr=output_unit, &
                  ionode=ionode)

             IF (p_loc) THEN
                IF (ionode) WRITE(output_unit,*) ""
                IF (ionode) WRITE(output_unit,*) "Computing localization properties"
                WRITE(filename,'(a13,I1.1)') "LUMO_CENTERS_",ispin
                IF (ionode) OPEN(UNIT=21,FILE=filename)
                ALLOCATE(center(unoccupied_orbs(ispin)%matrix%matrix_struct%ncol_global,6))
                CALL qs_loc_compute_properties(qs_loc_env,molecule_set,unoccupied_orbs(ispin)%matrix,&
                                center,ionode,21)
                DEALLOCATE(center)
                IF (ionode) CLOSE(21)
                IF (ionode) WRITE(output_unit,*) "Results in ",filename
             ENDIF

             IF (p_cubes) THEN
                DO ivector=1,nlumos
                   CALL calculate_wavefunction(mo_vectors=unoccupied_orbs(ispin)%matrix,&
                        ivector=ivector,rho=wf_r, &
                        rho_gspace=wf_g,&
                        qs_env=qs_env,error=error)

                   WRITE(filename,'(a13,I5.5,a1,I1.1,a5)') "WAVEFUNCTION_",homo+ivector,"_",ispin,".cube"
                   IF (ionode) OPEN(UNIT=21,FILE=filename)
                   WRITE(title,*) "WAVEFUNCTION ",homo+ivector," spin ",ispin," i.e. LUMO + ",ivector-1
                   CALL rs_pw_to_cube(wf_r%pw,21,ionode,title,rs,&
                        particles=particles,&
                        stride=section_get_ival(dft_section,"PRINT%MO_CUBES%STRIDE"))
                   IF (ionode) CLOSE(21)
                ENDDO
             ENDIF
          ENDDO

          !--------------------------------------!
          ! save the info for tddfpt calculation !
          !--------------------------------------!
          IF (dft_control%do_tddfpt_calculation) THEN
             ALLOCATE(tddfpt_control%lumos_eigenvalues(nlumos,dft_control%nspins))
             DO ispin=1, dft_control%nspins
                tddfpt_control%lumos_eigenvalues(1:nlumos,ispin) = &
                     unoccupied_evals(ispin)%evals(1:nlumos)
             END DO
             tddfpt_control%lumos => unoccupied_orbs
          END IF

       ENDIF

! *** deallocate grids needed to compute wavefunctions
       IF ((p_cubes .AND. (section_get_ival(dft_section,"PRINT%MO_CUBES%NLUMO") /= 0 .OR. &
            section_get_ival(dft_section,"PRINT%MO_CUBES%NHOMO") /= 0 )).OR. &
            p_loc ) THEN
          CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_r,&
               error=error)
          CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_g,&
               error=error)
          CALL rs_pool_give_back_rs(auxbas_rs_pool,rs, error=error)
       END IF

! *** destroy the localization environment
       IF(dft_control%qs_control%localized_wfn .AND. my_localized_wfn) THEN
          CALL qs_loc_env_destroy(qs_loc_env_new, error=error)
          DEALLOCATE(qs_loc_env_new,STAT=istat)
          CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ELSEIF (p_loc) THEN
          CALL qs_loc_destroy(qs_loc_env)
          DEALLOCATE(qs_loc_env) 
       ENDIF

! *** compute the optical conductivity (needs homos and lumos)
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%OPTICAL_CONDUCTIVITY",error=error),cp_p_file)) THEN
           IF (ionode) WRITE(output_unit,*) " "
           IF (ionode) WRITE(output_unit,*) " Computing the optical conductivity "
           IF (ionode) WRITE(output_unit,*) " Experimental version " 
           IF (ionode) WRITE(output_unit,*) " Check the code before believing results " 
           homoEQnmo=.true.
           DO ispin=1,dft_control%nspins
              CALL get_mo_set(mo_set=mos(ispin)%mo_set,homo=homo,nmo=nmo)
              if (homo.NE.nmo) homoEQnmo=.false. 
           ENDDO
           IF (.NOT. homoEQnmo) THEN
              IF (ionode) WRITE(output_unit,*) " homo.NE.nmo : skip optical conductivty "
           ELSE
             IF (.NOT.(p_cubes .AND. section_get_ival(dft_section,"PRINT%MO_CUBES%NLUMO") /= 0)) &
                CALL stop_program("optical_conductivity","put LUMOS XX in the print section")
             IF (.NOT. (p_cubes .AND. section_get_ival(dft_section,"PRINT%MO_CUBES%NHOMO") /= 0)) &
                CALL stop_program("optical_conductivity","put HOMOS XX in the print section")
             IF (ionode) THEN
                File_name="CONDUCTIVITY"
                CALL Open_file(File_name=File_name,&
                    File_action="Write",&
                    File_form="Formatted",&
                    File_status="Replace",&
                    Unit_number=File_unit)
             ENDIF
             CALL get_cell(cell,deth=volume)
             DO ispin=1,dft_control%nspins
                CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                                                       eigenvalues=mo_eigenvalues)
                CALL optical_conductivity(matrix_s,mo_coeff,mo_eigenvalues, &
                                            unoccupied_orbs(ispin)%matrix,unoccupied_evals(ispin)%evals, &
                                            volume,ionode,File_unit) 
             ENDDO
             IF (ionode) CALL close_file(Unit_number=File_unit)
             IF (ionode) WRITE(output_unit,*) " Results written to file ",File_name
           ENDIF
       ENDIF

       IF ((p_cubes .OR.dft_control%do_tddfpt_calculation ) &
          .AND. section_get_ival(dft_section,"PRINT%MO_CUBES%NLUMO") /= 0) THEN
        DO ispin=1,dft_control%nspins
           DEALLOCATE(unoccupied_evals(ispin)%evals)
           IF (.NOT.dft_control%do_tddfpt_calculation) &
                CALL cp_fm_release(unoccupied_orbs(ispin)%matrix)
        ENDDO
        DEALLOCATE(unoccupied_evals)
        IF (.NOT.dft_control%do_tddfpt_calculation) DEALLOCATE(unoccupied_orbs)
       ENDIF

    END IF

    ! Print coherent X-ray diffraction spectrum

    print_key => section_vals_get_subs_vals(&
                   section_vals=input,&
                   subsection_name="DFT%PRINT%XRAY_DIFFRACTION_SPECTRUM",&
                   error=error)
    IF (BTEST(cp_print_key_should_output(logger%iter_info,print_key,&
            error=error),cp_p_file)) THEN
      q_max = section_get_rval(section_vals=dft_section,&
                               keyword_name="PRINT%XRAY_DIFFRACTION_SPECTRUM%Q_MAX",&
                               error=error)
      file_name = "x-ray.dat"
      CALL xray_diffraction_spectrum(qs_env=qs_env,&
                                     file_name=file_name,&
                                     q_max=q_max,&
                                     error=error)
      IF (ionode) THEN
        WRITE (UNIT=output_unit,FMT="(/,/,T2,A,2(/,/,T3,A),/)")&
          "X-RAY DIFFRACTION SPECTRUM",&
          "The coherent X-ray diffraction spectrum was written to the file:",&
          TRIM(file_name)
      END IF
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE scf_post_calculation

  ! ***************************************************************************

  SUBROUTINE write_available_results(qs_env,globenv,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      pointer                     :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'write_available_results', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, ionode
    REAL(KIND=dp)                            :: total_abs_spin_dens
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type)                         :: wf_r
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(particle_list_type), pointer :: particles
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(cp_subsystem_type), pointer :: subsys

    TYPE evals_p_type
         REAL(KIND = dp), DIMENSION(:), POINTER :: evals
    END TYPE
    INTEGER :: ispin, output_unit
    REAL(KIND=dp)      :: s_square_ideal,s_square
    LOGICAL            :: all_equal
    REAL(KIND=dp), DIMENSION(:), POINTER :: occupation_numbers
    INTEGER            :: nmo, unit_nr
    TYPE(section_vals_type), POINTER :: dft_section, input
    TYPE(cp_rs_pool_type), POINTER           :: auxbas_rs_pool
    TYPE(realspace_grid_type), POINTER       :: rs
    INTEGER            :: handle


    CALL timeset(routineN,handle)
    failure=.FALSE.
    NULLIFY(cell,dft_control,pw_env,auxbas_pw_pool,pw_pools,&
            mos,atomic_kind_set,particle_set,rho,&
            ks_rmpv,matrix_s,scf_control,dft_section,&
            molecule_set,input,particles,subsys)
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

    dft_section => section_vals_get_subs_vals(input,"DFT",error=error)

!   *** Write last wavefunction to screen ***
       DO ispin=1,dft_control%nspins
          CALL write_mo_set(mos(ispin)%mo_set,atomic_kind_set,particle_set,4,6,dft_section,globenv)
       ENDDO

!   *** Integrated absolute spin density and spin contamination ***
       IF (dft_control%nspins.eq.2) THEN
          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
               pw_pools=pw_pools,auxbas_rs_pool=auxbas_rs_pool)
          CALL pw_pool_init_coeff(auxbas_pw_pool,wf_r,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          CALL coeff_copy(rho%rho_r(1),wf_r)
          CALL coeff_reduce(rho%rho_r(2),wf_r)
          total_abs_spin_dens=calculate_total_abs_rho(wf_r)
          IF (globenv%ionode) WRITE(UNIT=output_unit,fmt='(/,(T3,A,T60,F20.10))') &
               "Integrated absolute spin density  : ",total_abs_spin_dens
          CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_r, error=error)
          !
          ! XXX Fix Me XXX
          ! should be extended to the case where added MOs are present 
          !
          all_equal=.TRUE.
          DO ispin=1,dft_control%nspins
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,occupation_numbers=occupation_numbers,nmo=nmo)
             IF (nmo.GT.0) THEN
                all_equal=all_equal .AND. ALL(occupation_numbers.EQ.occupation_numbers(1))
             ENDIF
          ENDDO
          IF (.NOT. all_equal) THEN
             IF (ionode) WRITE(UNIT=output_unit,fmt='(T3,A)') &
                   "Warning: S**2 computation does not yet treat empty/fractional orbitals"
          ENDIF
          CALL compute_s_square(mos=mos, matrix_s=matrix_s, s_square=s_square, s_square_ideal=s_square_ideal)
          IF (ionode) WRITE(UNIT=output_unit,fmt='(T3,A,T60,2F10.5)') &
               "Ideal and single determinant S**2 : ",s_square_ideal,s_square
       ENDIF

!   *** Cube file of the electron density ***
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%E_DENSITY_CUBE", error=error),cp_p_file)) THEN ! write the electron densoity
          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
          CALL pw_env_get(pw_env, auxbas_rs_pool=auxbas_rs_pool)
          CALL rs_pool_create_rs(auxbas_rs_pool,rs)
          IF (dft_control%nspins.ne.1) THEN
             CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
             CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                  pw_pools=pw_pools)
             CALL pw_pool_init_coeff(auxbas_pw_pool,wf_r,&
                  use_data = REALDATA3D,&
                  in_space = REALSPACE, error=error)
             CALL coeff_copy(rho%rho_r(1),wf_r)
             CALL coeff_sumup(rho%rho_r(2),wf_r)
             IF (ionode) OPEN(UNIT=21,FILE="ELECTRON_DENSITY.cube")
             CALL rs_pw_to_cube(wf_r%pw,21,ionode,"SPIN 1 + SPIN 2 ELECTRON DENSITY",&
                  rs,particles=particles,&
                  stride=section_get_ival(dft_section,"PRINT%E_DENSITY_CUBE%STRIDE"))
             IF (ionode) CLOSE(21)
             CALL coeff_copy(rho%rho_r(1),wf_r)
             CALL coeff_reduce(rho%rho_r(2),wf_r)
             IF (ionode) OPEN(UNIT=21,FILE="SPIN_DENSITY.cube")
             CALL rs_pw_to_cube(wf_r%pw,21,ionode,"SPIN 1 DENSITY - SPIN 2 DENSITY",rs,&
                  particles=particles,&
                  stride=section_get_ival(dft_section,"PRINT%E_DENSITY_CUBE%STRIDE"))
             IF (ionode) CLOSE(21)
             CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_r, error=error)
          ELSE
             IF (ionode) OPEN(UNIT=21,FILE="ELECTRON_DENSITY.cube")
             CALL rs_pw_to_cube(rho%rho_r(1)%pw,21,ionode,"ELECTRON DENSITY",&
                  rs,particles=particles,&
                  stride=section_get_ival(dft_section,"PRINT%E_DENSITY_CUBE%STRIDE"))
             IF (ionode) CLOSE(21)
          ENDIF
          CALL rs_pool_give_back_rs(auxbas_rs_pool,rs,error=error)
       ENDIF

!   *** write the hartree potential ***
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%V_HARTREE_CUBE",error=error),cp_p_file)) THEN
          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
          CALL pw_env_get(pw_env, auxbas_rs_pool=auxbas_rs_pool)
          CALL rs_pool_create_rs(auxbas_rs_pool,rs)
          unit_nr=cp_print_key_unit_nr(logger,input,"DFT%PRINT%V_HARTREE_CUBE",&
               extension=".cube",middle_name="v_hartree",error=error)
          CALL rs_pw_to_cube(qs_env%ks_env%v_hartree_rspace%pw,&
               unit_nr,unit_nr>0,"HARTREE POTENTIAL",rs,particles=particles,&
               stride=section_get_ival(dft_section,"PRINT%V_HARTREE_CUBE%STRIDE"))
          CALL cp_print_key_finished_output(unit_nr,logger,input,&
               "DFT%PRINT%V_HARTREE_CUBE",error=error)
          CALL rs_pool_give_back_rs(auxbas_rs_pool,rs,error=error)
       ENDIF

!    *** stupid (non periodic) calculation of the dipole moment ***
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,"DFT%PRINT%DIPOLE",&
            error=error),cp_p_file)) THEN
          unit_nr=cp_print_key_unit_nr(logger,input, "DFT%PRINT%DIPOLE",&
               extension=".dipole", log_filename=.FALSE., error=error)
          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env,error=error)
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
               error=error)
          CALL calc_stupid_dipole(qs_env%subsys,rho%rho_r,qs_env%rho_core%pw,&
               aux_bas_pool=auxbas_pw_pool,cell=cell,globenv=globenv,&
               unit_nr=unit_nr,&
               print_header=.FALSE.,&
               error=error)
          CALL cp_print_key_finished_output(unit_nr, logger,input,&
               "DFT%PRINT%DIPOLE",error=error)
       END IF
 
!    *** write the density matrix ***
       IF (BTEST(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%AO_MATRICES/DENSITY",error=error),cp_p_file)) THEN
          DO ispin=1,dft_control%nspins
             CALL write_sparse_matrix(rho%rho_ao(ispin)%matrix,4,6,qs_env,globenv)
          END DO
       END IF

!    **** the Kohn-Sham matrix itself
       IF (btest(cp_print_key_should_output(logger%iter_info,input,&
            "DFT%PRINT%AO_MATRICES/KOHN_SHAM_MATRIX",error=error),cp_p_file)) THEN
          CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,&
               error=error,&
               calculate_forces=.FALSE.,just_energy=.FALSE.)
          CALL qs_ks_did_change(qs_env%ks_env,rho_changed=.TRUE.)
          CALL Write_sparse_matrix(ks_rmpv(1)%matrix,4,6,Qs_env,Globenv)
       END IF
    END IF

    CALL timestop(handle)

  END SUBROUTINE write_available_results

!!****f* qs_scf_post/compute_s_square *
!!
!!   NAME
!!      compute_s_square
!!
!!   FUNCTION
!!      compute s square of the single determinant defined by the 
!!      by the spin up and spin down orbitals
!!
!!   NOTES
!!      see E.q. 2.271 in Modern Quantum Chemistry
!!      szabo and ostlund
!!
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     07.2004 created [ Joost VandeVondele ]
!!
!!   SOURCE
!!*** **********************************************************************
   SUBROUTINE compute_s_square(mos, matrix_s, s_square, s_square_ideal)
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    REAL(KIND=dp)                            :: s_square, s_square_ideal

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_s_square', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, na, nb, ncol_local, &
                                                nrow, nrow_local, handle
    REAL(KIND=dp)                            :: tmp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data
    TYPE(cp_blacs_env_type), POINTER         :: context
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: c_alpha, c_beta, &
                                                matrix_overlap, matrix_sc
    TYPE(cp_para_env_type), POINTER          :: para_env

    CALL timeset(routineN,handle)

    NULLIFY(fm_struct_tmp,matrix_sc,matrix_overlap,para_env,context,local_data)

    SELECT CASE (SIZE(mos))
    CASE (1)
          s_square=0.0_dp
          s_square_ideal=0.0_dp
    CASE (2)
          CALL get_mo_set(mo_set=mos(1)%mo_set,mo_coeff=c_alpha)
          CALL get_mo_set(mo_set=mos(2)%mo_set,mo_coeff=c_beta)
          CALL cp_fm_get_info(c_alpha,ncol_global=na)
          CALL cp_fm_get_info(c_beta,ncol_global=nb)
          s_square_ideal=(na-nb)*(na-nb+2)/4.0_dp
          ! create overlap matrix
          CALL cp_fm_get_info(c_alpha,para_env=para_env,context=context)
          CALL cp_fm_struct_create(fm_struct_tmp,para_env=para_env,context=context, &
                                   nrow_global=na,ncol_global=nb)
          CALL cp_fm_create(matrix_overlap, fm_struct_tmp,name="matrix_overlap")
          CALL cp_fm_struct_release(fm_struct_tmp)
          ! create S C_beta and compute overlap
          CALL cp_fm_get_info(c_beta, matrix_struct=fm_struct_tmp,nrow_global=nrow)
          CALL cp_fm_create(matrix_sc, fm_struct_tmp,name="matrix_sc")
          CALL cp_sm_fm_multiply(matrix_s(1)%matrix,c_beta,matrix_sc,nb)
          CALL cp_fm_gemm('T','N',na,nb,nrow,1.0_dp,c_alpha,matrix_sc,0.0_dp,matrix_overlap)
          ! invoke formula 2.271
          CALL cp_fm_get_info(matrix_overlap,local_data=local_data, nrow_local=nrow_local, ncol_local=ncol_local)
          tmp=0.0_dp
          DO j=1,ncol_local
          DO i=1,nrow_local
             tmp=tmp+local_data(i,j)**2
          ENDDO
          ENDDO
          CALL mp_sum(tmp,para_env%group)
          s_square=s_square_ideal+nb-tmp
          CALL cp_fm_release(matrix_sc)
          CALL cp_fm_release(matrix_overlap)
    CASE DEFAULT
          CALL stop_program("compute_s_square","alpha, beta, what else ?")
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE

!!****f* qs_scf_post/calc_stupid_dipole *
!!
!!   NAME
!!     calc_stupid_dipole
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - subsystem: the atoms to calculate the center of mass
!!     - rho_r: the electronic density in real space
!!     - rho_core_g: the charge of the core in the reciprocal space
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
SUBROUTINE calc_stupid_dipole(subsystem, rho_r, rho_core_g,aux_bas_pool,&
     cell,globenv,unit_nr,print_header,error)
    TYPE(cp_subsystem_type), POINTER         :: subsystem
    TYPE(coeff_type), DIMENSION(:), &
      INTENT(in)                             :: rho_r
    TYPE(pw_type), POINTER                   :: rho_core_g
    TYPE(pw_pool_type), POINTER              :: aux_bas_pool
    TYPE(cell_type), POINTER                 :: cell
    TYPE(global_environment_type), &
      pointer                             :: globenv
    INTEGER, INTENT(in)                      :: unit_nr
    LOGICAL, INTENT(in), OPTIONAL            :: print_header
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_stupid_dipole', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_string_length)     :: iter
    INTEGER                                  :: i, ip, j, k
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure, my_print_header
    REAL(kind=dp)                            :: mass, tot_m
    REAL(kind=dp), DIMENSION(3)              :: com, dipole_center, &
                                                dipole_com, dr, r, r_com,&
                                                shift
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(pw_type), POINTER                   :: total_charge

  failure=.FALSE.
  my_print_header=.FALSE.
  IF (PRESENT(print_header)) my_print_header=print_header 

  logger => cp_error_get_logger(error)
  com=0._dp
  tot_m=0._dp
  CALL cp_subsys_get(subsystem,particles=particles,error=error)
  DO ip=1,particles%n_els
     CALL get_atomic_kind(particles%els(ip)%atomic_kind,mass=mass)
     tot_m=tot_m+mass
     r=particles%els(ip)%r
     com=com+mass*r
     !com(1)=com(1)+mass*(MODULO(r(1),cell%hmat(1,1))-cell%hmat(1,1)/2._dp)
     !com(2)=com(2)+mass*(MODULO(r(2),cell%hmat(2,2))-cell%hmat(2,2)/2._dp)
     !com(3)=com(3)+mass*(MODULO(r(3),cell%hmat(3,3))-cell%hmat(3,3)/2._dp)
  END DO
  com=com/tot_m
  com(1)=MODULO(com(1),cell%hmat(1,1))
  com(2)=MODULO(com(2),cell%hmat(2,2))
  com(3)=MODULO(com(3),cell%hmat(3,3))

  
  CALL pw_pool_create_pw(aux_bas_pool,total_charge,use_data=REALDATA3D,&
       in_space=REALSPACE,error=error)
  CALL pw_transfer(rho_core_g,total_charge)
  DO i=1,SIZE(rho_r)
     CALL pw_sumup(rho_r(i)%pw,total_charge)
  END DO

  dipole_center=0._dp
  dipole_com=0._dp
  bo=total_charge%pw_grid%bounds_local
  dr=total_charge%pw_grid%dr
  shift=REAL(MODULO(total_charge%pw_grid%npts,2),dp)*dr
  DO k=bo(1,3),bo(2,3)
     r(3)=k*dr(3)+shift(3)
     r_com(3)=MODULO(r(3)-com(3),cell%hmat(3,3))-cell%hmat(3,3)/2._dp
     DO j=bo(1,2),bo(2,2)
        r(2)=j*dr(2)+shift(2)
        r_com(2)=MODULO(r(2)-com(2),cell%hmat(2,2))-cell%hmat(2,2)/2._dp
        r(1)=bo(1,1)*dr(1)+shift(1)
        r_com(1)=MODULO(r(1)-com(1),cell%hmat(1,1))-cell%hmat(1,1)/2._dp
        DO i=bo(1,1),bo(2,1)
           dipole_com=dipole_com+total_charge%cr3d(i,j,k)*r_com
           dipole_center=dipole_center+total_charge%cr3d(i,j,k)*r
           r(1)=r(1)+dr(1)
           r_com(1)=MODULO(r(1)-com(1),cell%hmat(1,1))-cell%hmat(1,1)/2._dp
        END DO
     END DO
  END DO

  dipole_center=dipole_center*total_charge%pw_grid%dvol
  dipole_com=dipole_com*total_charge%pw_grid%dvol

  CALL mp_sum(dipole_center,total_charge%pw_grid%para%group)
  CALL mp_sum(dipole_com,total_charge%pw_grid%para%group)

  IF (unit_nr>0) THEN
     IF (my_print_header) THEN
        WRITE(unit=unit_nr,fmt="('# iteration com(x,y,z)[Angstrom]"//&
             " dipole_com(x,y,z)[au] dipole_cell_center(x,y,z)[au]')")
     END IF
     iter=cp_iter_string(logger%iter_info,error=error)
     WRITE(unit=unit_nr,fmt="(a,9(es18.8))")&
          iter(1:15),&
          com*angstrom,-dipole_com,-dipole_center
  END IF
  CALL pw_pool_give_back_pw(aux_bas_pool,total_charge,error=error)

END SUBROUTINE calc_stupid_dipole
!***************************************************************************

END MODULE
