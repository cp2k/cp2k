!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2003  CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

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

  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE coefficient_types,               ONLY: coeff_copy,&
                                             coeff_reduce,&
                                             coeff_sumup,&
                                             coeff_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             tddfpt_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
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
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: WANNIER_NONE,&
                                             global_environment_type
  USE kinds,                           ONLY: dp
  USE level_shift,                     ONLY: level_shift_testing
  USE message_passing,                 ONLY: mp_sum
  USE molecular_states,                ONLY: construct_molecular_states
  USE molecule_types_new,              ONLY: molecule_type
  USE mulliken,                        ONLY: mulliken_charges
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE
  USE qs_collocate_density,            ONLY: calculate_total_abs_rho,&
                                             calculate_wavefunction
  USE qs_conductivity,                 ONLY: optical_conductivity
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_methods,                   ONLY: qs_ks_update_qs_env
  USE qs_localization,                 ONLY: qs_loc_compute_properties,&
                                             qs_loc_destroy,&
                                             qs_loc_env_type,&
                                             qs_loc_init
  USE qs_mo_methods,                   ONLY: calculate_density_matrix,&
                                             calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             write_mo_set
  USE qs_ot_eigensolver,               ONLY: ot_eigensolver
  USE qs_overlap,                      ONLY: write_sparse_matrix
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE realspace_grid_types,            ONLY: rs_pw_to_cube
  USE scf_control_types,               ONLY: scf_control_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE sparse_matrix_types,             ONLY: add_all_real_matrix_blocks,&
                                             deallocate_matrix,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix,&
                                             write_blocks_maxabsval
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

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
  SUBROUTINE scf_post_calculation(scf_env,qs_env,globenv,error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scf_env_post_calc', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, ionode
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
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_loc_env_type), POINTER           :: qs_loc_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv, matrix_s, matrix_p
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(tddfpt_control_type), POINTER       :: tddfpt_control

    TYPE evals_p_type
         REAL(KIND = dp), DIMENSION(:), POINTER :: evals
    END TYPE
    TYPE (evals_p_type), DIMENSION(:), POINTER :: unoccupied_evals
    INTEGER :: ispin,n, output_unit, nmo,handle
    CHARACTER(LEN=200) :: file_name
    INTEGER            :: file_unit,ivector,nao,nlumos, homo
    CHARACTER(LEN=80)  :: filename
    CHARACTER(LEN=80) :: title
    REAL(KIND = dp)           :: volume
    REAL(KIND = dp), POINTER, DIMENSION(:,:) :: center
 
    CALL timeset(routineN,"I",'',handle)

! writes the data that is already available in qs_env 
    CALL write_available_results(qs_env,globenv,error)

    failure=.FALSE.
    NULLIFY(cell,dft_control,pw_env,auxbas_pw_pool,pw_pools,&
            mos,atomic_kind_set,particle_set,rho,mo_coeff,&
            ks_rmpv,matrix_s,qs_loc_env,scf_control,unoccupied_orbs,&
            mo_eigenvalues,unoccupied_evals,fm_struct_tmp,molecule_set,&
            center,rho_core,tddfpt_control)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)

    CPPrecondition(ASSOCIATED(scf_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(scf_env%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,dft_control=dft_control,molecule_set=molecule_set, &
            mos=mos,atomic_kind_set=atomic_kind_set,particle_set=particle_set,&
            rho=rho,matrix_ks=ks_rmpv,scf_control=scf_control,matrix_s=matrix_s, &
            cell=cell,rho_core=rho_core,tddfpt_control=tddfpt_control,error=error)
       CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)

!    *** compute the mulliken charges ***
       IF (globenv%print%mulliken_analysis)  &
            CALL mulliken_charges(rho%rho_ao,matrix_s(1)%matrix,globenv%para_env,particle_set, &
                                  globenv%scr,globenv%ionode)

!    *** setup the grids needed to compute a wavefunction given a vector
       IF ((globenv%print%homos .OR. globenv%print%lumos .OR. globenv%print%localization) &
                   .AND. globenv%print%cubes ) THEN
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
               pw_pools=pw_pools)
          CALL pw_pool_init_coeff(auxbas_pw_pool,wf_r,&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          CALL pw_pool_init_coeff(auxbas_pw_pool,wf_g,&
               use_data = COMPLEXDATA1D,&
               in_space = RECIPROCALSPACE, error=error)
       END IF

! *** makes the MOs eigenstates, computes eigenvalues, write cubes 
       IF (globenv%print%homos) THEN
          DO ispin=1,dft_control%nspins
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                                              eigenvalues=mo_eigenvalues,homo=homo,nmo=nmo)
             IF (globenv%ionode) WRITE(globenv%scr,*) " "
             IF (globenv%ionode) WRITE(globenv%scr,*) " Eigenvalues of the occupied subspace spin ",ispin
             IF (homo .NE. nmo) THEN
                IF (globenv%ionode) WRITE(globenv%scr,*) &
                    " and ",nmo-homo," added MO eigenvalues" 
             ENDIF
             IF (globenv%ionode) WRITE(globenv%scr,*) " ---------------------------------------------"
             CALL calculate_subspace_eigenvalues(mo_coeff,&
                  ks_rmpv(ispin)%matrix,mo_eigenvalues, &
                  para_env=qs_env%para_env,scr=globenv%scr, &
                  ionode=globenv%ionode, do_rotation=.TRUE.)

             IF (globenv%print%cubes) THEN
                DO ivector=MAX(1,homo-globenv%print%nhomos+1),homo
                   CALL calculate_wavefunction(mo_coeff,ivector,wf_r, &
                        wf_g,&
                        qs_env,error=error)

                   WRITE(filename,'(a13,I5.5,a1,I1.1,a5)') "WAVEFUNCTION_",ivector,"_",ispin,".cube"
                   IF (globenv%ionode) OPEN(UNIT=21,FILE=filename)
                   WRITE(title,*) "WAVEFUNCTION ",ivector," spin ",ispin," i.e. HOMO - ",ivector-homo
                   CALL rs_pw_to_cube(wf_r%pw,21,globenv%ionode,title)
                   IF (globenv%ionode) CLOSE(21)
                ENDDO
             END IF
          ENDDO
       ENDIF

! *** init the localization environment, needed e.g. for wannier functions and molecular states
       IF (globenv%print%localization)THEN
          NULLIFY(qs_loc_env)
          ALLOCATE(qs_loc_env)
          CALL qs_loc_init(qs_loc_env,qs_env)
          qs_loc_env%localization_method = globenv % print % localization_method
          qs_loc_env%wannier = .NOT. (qs_loc_env % localization_method .EQ. WANNIER_NONE)
          IF (globenv%print%molecular_states) qs_loc_env%molecular_states = .TRUE.
       ENDIF

! *** gets localization info for the occupied orbs
! *** possibly gets wannier functions
! *** possibly gets molecular states
       IF (globenv%print%localization) THEN
          DO ispin=1,dft_control%nspins
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff, &
                                              eigenvalues=mo_eigenvalues,homo=homo,nmo=nmo)
             ! just get eigenstates
             IF (.NOT.globenv%print%homos) THEN
                IF (globenv%ionode) WRITE(globenv%scr,*) " "
                IF (globenv%ionode) WRITE(globenv%scr,*) " Eigenvalues of the occupied subspace spin ",ispin
                IF (homo .NE. nmo) THEN
                   IF (globenv%ionode) WRITE(globenv%scr,*) &
                       " and ",nmo-homo," added MO eigenvalues" 
                ENDIF
                IF (globenv%ionode) WRITE(globenv%scr,*) " ---------------------------------------------"
                CALL calculate_subspace_eigenvalues(mo_coeff,ks_rmpv(ispin)%matrix, &
                     mo_eigenvalues, para_env=qs_env%para_env,scr=globenv%scr, &
                     ionode=globenv%ionode, do_rotation=.TRUE.)
             ENDIF
             NULLIFY(mo_localized)
             CALL cp_fm_create(mo_localized,mo_coeff%matrix_struct)
             ALLOCATE(center(mo_localized%matrix_struct%ncol_global,6)) 
             CALL cp_fm_to_fm(mo_coeff,mo_localized)
             ! do localization computation
             WRITE(filename,'(a13,I1.1)') "HOMO_CENTERS_",ispin
             IF (globenv%ionode) WRITE(globenv%scr,*) ""
             IF (globenv%ionode) WRITE(globenv%scr,*) "Computing localization properties"
             IF (globenv%ionode) OPEN(UNIT=21,FILE=filename)
             CALL qs_loc_compute_properties(qs_loc_env,molecule_set,mo_localized,center,globenv%ionode,21)
             IF (globenv%ionode) CLOSE(21)
             IF (globenv%ionode) WRITE(globenv%scr,*) "Results in ",filename
             ! write the wannier functions to cubes
             IF (globenv%print%cubes.AND. .NOT. globenv % print % molecular_states ) THEN
                ! obviously we need a smarter scheme to select them. 
                DO ivector=1,mo_localized%matrix_struct%ncol_global   
                   CALL calculate_wavefunction(mo_localized,ivector,wf_r, &
                        wf_g,&
                        qs_env,error=error)

                   WRITE(filename,'(a7,I5.5,a1,I1.1,a5)') "WANNIER",ivector,"_",ispin,".cube"
                   IF (globenv%ionode) OPEN(UNIT=21,FILE=filename)
                   ! title name
                   WRITE(title,'(A7,I5.5,A1,I1.1,A1,6G10.4)') "WANNIER",ivector,"_",ispin," ",&
                                             center(ivector,1:3)*angstrom, &
                                             center(ivector,4:6)*angstrom
                   CALL rs_pw_to_cube(wf_r%pw,21,globenv%ionode,title)
                   IF (globenv%ionode) CLOSE(21)
                ENDDO
             END IF
             ! compute the molecular states
             IF ( globenv % print % molecular_states ) THEN
                  CALL construct_molecular_states ( molecule_set, mo_localized, &
                                                    mo_coeff, mo_eigenvalues, &
                                                    ks_rmpv ( ispin ) % matrix, &
                                                    matrix_s(1)%matrix, &
                                                    globenv, qs_env,    &
                                                    wf_r, wf_g)
             ENDIF
             CALL cp_fm_release(mo_localized)

             IF (.FALSE.) &
             CALL level_shift_testing(atomic_kind_set, particle_set, cell, center, &
                                      pw_env, rho, rho_core)

             DEALLOCATE(center)
          ENDDO
       ENDIF

! gets the lumos, and eigenvalues for the lumos 
       IF (globenv%print%lumos) THEN
          ALLOCATE(unoccupied_orbs(dft_control%nspins))
          ALLOCATE(unoccupied_evals(dft_control%nspins))
          DO ispin=1,dft_control%nspins
             NULLIFY(unoccupied_orbs(ispin)%matrix)
             NULLIFY(unoccupied_evals(ispin)%evals)
! write eigenvalues and a set of cubes
             IF (globenv%ionode) WRITE(globenv%scr,*) " "
             IF (globenv%ionode) WRITE(globenv%scr,*) " Lowest Eigenvalues of the unoccupied subspace spin ",ispin
             IF (globenv%ionode) WRITE(globenv%scr,*) " -----------------------------------------------------"
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff,homo=homo,nao=nao,nmo=nmo)
             CALL cp_fm_get_info(mo_coeff, nrow_global=n)
             nlumos=MAX(1,MIN(globenv%print%nlumos,nao-homo))
             ALLOCATE(unoccupied_evals(ispin)%evals(nlumos))
             CALL cp_fm_struct_create(fm_struct_tmp,para_env=qs_env%para_env,context=qs_env%blacs_env, &
                                                    nrow_global=n,ncol_global=nlumos)
             CALL cp_fm_create(unoccupied_orbs(ispin)%matrix, fm_struct_tmp,name="lumos")
             CALL cp_fm_struct_release(fm_struct_tmp)
             CALL cp_fm_init_random(unoccupied_orbs(ispin)%matrix,nlumos)

             CALL ot_eigensolver(matrix_h=ks_rmpv(ispin)%matrix,matrix_s=matrix_s(1)%matrix, &
                  matrix_c=unoccupied_orbs(ispin)%matrix, &
                  matrix_orthogonal_space=mo_coeff, &
                  eps_gradient=scf_control%eps_lumos, &
                  preconditioner=scf_env%ot_preconditioner, &
                  iter_max=scf_control%max_iter_lumos,globenv=globenv,&
                  size_ortho_space=homo)

             CALL calculate_subspace_eigenvalues(unoccupied_orbs(ispin)%matrix,ks_rmpv(ispin)%matrix,&
                  unoccupied_evals(ispin)%evals, para_env=qs_env%para_env,scr=globenv%scr, &
                  ionode=globenv%ionode)

             IF (globenv%print%localization) THEN
                IF (globenv%ionode) WRITE(globenv%scr,*) ""
                IF (globenv%ionode) WRITE(globenv%scr,*) "Computing localization properties"
                WRITE(filename,'(a13,I1.1)') "LUMO_CENTERS_",ispin
                IF (globenv%ionode) OPEN(UNIT=21,FILE=filename)
                ALLOCATE(center(unoccupied_orbs(ispin)%matrix%matrix_struct%ncol_global,6))
                CALL qs_loc_compute_properties(qs_loc_env,molecule_set,unoccupied_orbs(ispin)%matrix,&
                                center,globenv%ionode,21)
                DEALLOCATE(center)
                IF (globenv%ionode) CLOSE(21)
                IF (globenv%ionode) WRITE(globenv%scr,*) "Results in ",filename
             ENDIF

             IF (globenv%print%cubes) THEN
                DO ivector=1,nlumos
                   CALL calculate_wavefunction(mo_vectors=unoccupied_orbs(ispin)%matrix,&
                        ivector=ivector,rho=wf_r, &
                        rho_gspace=wf_g,&
                        qs_env=qs_env,error=error)

                   WRITE(filename,'(a13,I5.5,a1,I1.1,a5)') "WAVEFUNCTION_",homo+ivector,"_",ispin,".cube"
                   IF (globenv%ionode) OPEN(UNIT=21,FILE=filename)
                   WRITE(title,*) "WAVEFUNCTION ",homo+ivector," spin ",ispin," i.e. LUMO + ",ivector-1
                   CALL rs_pw_to_cube(wf_r%pw,21,globenv%ionode,title)
                   IF (globenv%ionode) CLOSE(21)
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
       IF ((globenv%print%lumos .OR. globenv%print%homos .OR. globenv%print%localization) &
                          .AND. globenv%print%cubes) THEN
          CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_r,&
               error=error)
          CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_g,&
               error=error)
       END IF

! *** destroy the localization environment
       IF (globenv%print%localization) THEN
          CALL qs_loc_destroy(qs_loc_env)
          DEALLOCATE(qs_loc_env) 
       ENDIF

! *** compute the optical conductivity (needs homos and lumos)
       IF (globenv%print%optical_conductivity) THEN
           IF (globenv%ionode) WRITE(globenv%scr,*) " "
           IF (globenv%ionode) WRITE(globenv%scr,*) " Computing the optical conductivity "
           IF (globenv%ionode) WRITE(globenv%scr,*) " Experimental version " 
           IF (globenv%ionode) WRITE(globenv%scr,*) " Check the code before believing results " 
           IF (.NOT. globenv%print%lumos) &
              CALL stop_program("optical_conductivity","put LUMOS XX in the print section")
           IF (.NOT. globenv%print%homos) &
              CALL stop_program("optical_conductivity","put HOMOS XX in the print section")
           IF (globenv%ionode) THEN
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
                                          volume,globenv%ionode,File_unit) 
           ENDDO
           IF (globenv%ionode) CALL Close_file(Unit_number=File_unit)
           IF (globenv%ionode) WRITE(globenv%scr,*) " Results written to file ",File_name
       ENDIF

       IF (globenv%print%lumos) THEN
        DO ispin=1,dft_control%nspins
           DEALLOCATE(unoccupied_evals(ispin)%evals)
           IF (.NOT.dft_control%do_tddfpt_calculation) &
                CALL cp_fm_release(unoccupied_orbs(ispin)%matrix)
        ENDDO
        DEALLOCATE(unoccupied_evals)
        IF (.NOT.dft_control%do_tddfpt_calculation) DEALLOCATE(unoccupied_orbs)
       ENDIF

    END IF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE scf_post_calculation

  SUBROUTINE write_available_results(qs_env,globenv,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(global_environment_type), &
      INTENT(in), TARGET                     :: globenv
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
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_rmpv, matrix_s
    TYPE(scf_control_type), POINTER          :: scf_control

    TYPE evals_p_type
         REAL(KIND = dp), DIMENSION(:), POINTER :: evals
    END TYPE
    TYPE(real_matrix_type), POINTER :: pfull
    INTEGER :: ispin, output_unit
    CHARACTER(LEN=6)   :: extension
    CHARACTER(LEN=200) :: file_name
    INTEGER            :: file_unit
    REAL(KIND=dp)      :: s_square_ideal,s_square
    LOGICAL            :: all_equal
    REAL(KIND=dp), DIMENSION(:), POINTER :: occupation_numbers

    failure=.FALSE.
    NULLIFY(cell,dft_control,pw_env,auxbas_pw_pool,pw_pools,&
            mos,atomic_kind_set,particle_set,rho,&
            ks_rmpv,matrix_s,scf_control,pfull,&
            molecule_set)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source
    IF (ionode) output_unit= cp_logger_get_default_unit_nr(logger)

    CPPrecondition(ASSOCIATED(qs_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL get_qs_env(qs_env,dft_control=dft_control,molecule_set=molecule_set, &
            mos=mos,atomic_kind_set=atomic_kind_set,particle_set=particle_set,&
            rho=rho,matrix_ks=ks_rmpv,scf_control=scf_control,matrix_s=matrix_s,cell=cell,error=error)

!   *** Write last wavefunction to screen ***
       DO ispin=1,dft_control%nspins
          CALL write_mo_set(mos(ispin)%mo_set,atomic_kind_set,particle_set,4,6,globenv)
       ENDDO

!   *** Integrated absolute spin density and spin contamination ***
       IF (dft_control%nspins.eq.2) THEN
          CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
          CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
               pw_pools=pw_pools)
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
             CALL get_mo_set(mo_set=mos(ispin)%mo_set,occupation_numbers=occupation_numbers)
             all_equal=all_equal .AND. ALL(occupation_numbers.EQ.occupation_numbers(1))
          ENDDO
          IF (.NOT. all_equal) THEN
             IF (globenv%ionode) WRITE(UNIT=output_unit,fmt='(T3,A)') &
                   "Warning: S**2 computation does not yet treat empty/fractional orbitals"
          ENDIF
          CALL compute_s_square(mos=mos, matrix_s=matrix_s, s_square=s_square, s_square_ideal=s_square_ideal)
          IF (globenv%ionode) WRITE(UNIT=output_unit,fmt='(T3,A,T60,2F10.5)') &
               "Ideal and single determinant S**2 : ",s_square_ideal,s_square
       ENDIF

!   *** Cube file of the electron density ***
       IF (globenv%print%e_density_cube) THEN ! write the electron density
          IF (dft_control%nspins.ne.1) THEN
             CALL get_qs_env(qs_env=qs_env,pw_env=pw_env)
             CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
                  pw_pools=pw_pools)
             CALL pw_pool_init_coeff(auxbas_pw_pool,wf_r,&
                  use_data = REALDATA3D,&
                  in_space = REALSPACE, error=error)
             CALL coeff_copy(rho%rho_r(1),wf_r)
             CALL coeff_sumup(rho%rho_r(2),wf_r)
             IF (globenv%ionode) OPEN(UNIT=21,FILE="ELECTRON_DENSITY.cube")
             CALL rs_pw_to_cube(wf_r%pw,21,globenv%ionode,"SPIN 1 + SPIN 2 ELECTRON DENSITY")
             IF (globenv%ionode) CLOSE(21)
             CALL coeff_copy(rho%rho_r(1),wf_r)
             CALL coeff_reduce(rho%rho_r(2),wf_r)
             IF (globenv%ionode) OPEN(UNIT=21,FILE="SPIN_DENSITY.cube")
             CALL rs_pw_to_cube(wf_r%pw,21,globenv%ionode,"SPIN 1 DENSITY - SPIN 2 DENSITY")
             IF (globenv%ionode) CLOSE(21)
             CALL pw_pool_give_back_coeff(auxbas_pw_pool,wf_r, error=error)
          ELSE
             IF (globenv%ionode) OPEN(UNIT=21,FILE="ELECTRON_DENSITY.cube")
             CALL rs_pw_to_cube(rho%rho_r(1)%pw,21,globenv%ionode,"ELECTRON DENSITY")
             IF (globenv%ionode) CLOSE(21)
          ENDIF
       ENDIF

!   *** write the hartree potential ***
       IF (globenv%print%v_hartree_cube) THEN
          IF (globenv%ionode) OPEN(UNIT=21,FILE="V_HARTREE.cube")
          CALL rs_pw_to_cube(qs_env%ks_env%v_hartree_rspace%pw,&
               21,globenv%ionode,"HARTREE POTENTIAL")
          IF (globenv%ionode) CLOSE(21)
       ENDIF

!    *** write the density matrix ***
       IF (globenv%print%density_matrix) THEN
         DO ispin=1,dft_control%nspins
           CALL write_sparse_matrix(rho%rho_ao(ispin)%matrix,4,6,qs_env,globenv)
         END DO
       END IF

!    *** a coarse grained version of the density matrix ***
       IF (Globenv%print%density_matrix_magnitude) THEN
          IF (globenv%num_pe > 1) THEN
             WRITE (unit=extension,fmt="(I6)") globenv%mepos
             file_name = "P_magnitude"//"."//ADJUSTL(extension)
          ELSE
             file_name = "P_magnitude"
          END IF
          CALL open_file(file_name=file_name,&
               file_action="Write",&
               file_form="Formatted",&
               file_status="Replace",&
               unit_number=file_unit)

          NULLIFY(pfull) 
          CALL replicate_matrix(matrix_s(1)%matrix,Pfull,"Full Density Matrix",.FALSE.)
          CALL add_all_real_matrix_blocks(Pfull,Globenv)
          CALL calculate_density_matrix(mos(1)%mo_set,Pfull,error=error)
          CALL write_blocks_maxabsval(Pfull,File_unit)
          CALL deallocate_matrix(Pfull)
          CALL close_file(Unit_number=File_unit)
       ENDIF

!    **** the Kohn-Sham matrix itself
       IF (globenv%print%kohn_sham_matrix) THEN
          CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,&
               error=error,&
               calculate_forces=.FALSE.,just_energy=.FALSE.)
          CALL Write_sparse_matrix(ks_rmpv(1)%matrix,4,6,Qs_env,Globenv)
       END IF
    END IF
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

    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: matrix_sc,matrix_overlap
    TYPE(cp_fm_type), POINTER                :: c_alpha,c_beta
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_blacs_env_type), POINTER         :: context 
    REAL(KIND=dp), DIMENSION(:,:), POINTER   :: local_data
    INTEGER :: na, nb, i, j, nrow
    REAL(KIND=dp) :: tmp

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
          CALL cp_fm_get_info(matrix_overlap,local_data=local_data)
          tmp=SUM(local_data**2)
          CALL mp_sum(tmp,para_env%group)
          s_square=s_square_ideal+nb-tmp
          CALL cp_fm_release(matrix_sc)
          CALL cp_fm_release(matrix_overlap)
    CASE DEFAULT
          CALL stop_program("compute_s_square","alpha, beta, what else ?")
    END SELECT

  END SUBROUTINE

END MODULE
