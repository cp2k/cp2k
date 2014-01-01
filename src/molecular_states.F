!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for the calculation of molecular states
!> \author CJM
! *****************************************************************************
MODULE molecular_states
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE bibliography,                    ONLY: Hunt2003,&
                                             cite_reference
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_diag,                      ONLY: choose_eigv_solver
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_section_types,             ONLY: section_get_ivals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_maxloc
  USE molecule_types_new,              ONLY: molecule_type
  USE particle_list_types,             ONLY: particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_type
  USE pw_types,                        ONLY: pw_p_type
  USE qs_collocate_density,            ONLY: calculate_wavefunction
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE realspace_grid_cube,             ONLY: pw_to_cube
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'molecular_states'

  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

! *** Public subroutines ***

  PUBLIC :: construct_molecular_states

CONTAINS

! *****************************************************************************
!> \brief constructs molecular states. mo_localized gets overwritten!
! *****************************************************************************
  SUBROUTINE construct_molecular_states(molecule_set, mo_localized, &
       mo_coeff, mo_eigenvalues, Hks, matrix_S, qs_env, wf_r, wf_g,&
       loc_print_section, particles, tag, marked_states,error)

    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(cp_fm_type), POINTER                :: mo_localized, mo_coeff
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(cp_dbcsr_type), POINTER             :: Hks, matrix_S
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type), INTENT(INOUT)           :: wf_r, wf_g
    TYPE(section_vals_type), POINTER         :: loc_print_section
    TYPE(particle_list_type), POINTER        :: particles
    CHARACTER(LEN=*), INTENT(IN)             :: tag
    INTEGER, DIMENSION(:), POINTER           :: marked_states
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'construct_molecular_states', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename
    CHARACTER(LEN=default_string_length)     :: title
    INTEGER :: handle, i, imol, iproc, isos, k, n_rep, ncol_global, nproc, &
      nrow_global, ns, output_unit, unit_nr, unit_report
    INTEGER, DIMENSION(:), POINTER           :: ind, mark_list
    INTEGER, DIMENSION(:, :), POINTER        :: mark_states
    INTEGER, POINTER                         :: nstates(:), states(:)
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: tmp
    REAL(KIND=dp), ALLOCATABLE               :: evals( : )
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eval_range
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: b, c, d, D_igk, e_vectors, &
                                                rot_e_vectors, smo, storage
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env

    CALL cite_reference(Hunt2003)

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(logger,mark_states,mark_list)
    logger => cp_error_get_logger(error)
    !-----------------------------------------------------------------------------
    ! 1.
    !-----------------------------------------------------------------------------
    para_env => qs_env % para_env
    nproc = para_env%num_pe
    output_unit = cp_logger_get_default_io_unit(logger)
    CALL section_vals_val_get(loc_print_section,"MOLECULAR_STATES%CUBE_EVAL_RANGE",&
               explicit=explicit,error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(loc_print_section,"MOLECULAR_STATES%CUBE_EVAL_RANGE",&
                  r_vals=eval_range,error=error)
    ELSE
       ALLOCATE(eval_range(2))
       eval_range(1)=-HUGE(0.0_dp)
       eval_range(2)=+HUGE(0.0_dp)
    ENDIF
    CALL section_vals_val_get(loc_print_section,"MOLECULAR_STATES%MARK_STATES", &
             n_rep_val=n_rep, error=error)
    IF(n_rep.GT.0)THEN
       ALLOCATE (mark_states(2,n_rep),STAT=isos)
       IF (isos /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "mark_states",int_size*2*n_rep)
       IF (.NOT.ASSOCIATED(marked_states)) THEN
          ALLOCATE (marked_states(n_rep),stat=isos)
          IF (isos /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                          "mark_states",int_size*n_rep)
       END IF
       DO i=1,n_rep
          CALL section_vals_val_get(loc_print_section,"MOLECULAR_STATES%MARK_STATES",&
                i_rep_val=i,i_vals=mark_list, error=error)
          mark_states(:,i)=mark_list(:)
       END DO
    ELSE
       NULLIFY(marked_states)
    END IF

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 2.
    !-----------------------------------------------------------------------------
    unit_report=cp_print_key_unit_nr(logger,loc_print_section,"MOLECULAR_STATES",&
         extension=".data",middle_name="Molecular_DOS",log_filename=.FALSE.,error=error)
    IF (unit_report>0) THEN
       WRITE(unit_report,*) SIZE(mo_eigenvalues)," number of states "
       DO i=1,SIZE(mo_eigenvalues)
          WRITE(unit_report,*) mo_eigenvalues(i)
       ENDDO
    ENDIF

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 3.
    !-----------------------------------------------------------------------------
    CALL cp_fm_get_info(mo_localized, &
                        ncol_global=ncol_global, &
                        nrow_global=nrow_global ,error=error)
    NULLIFY(smo)
    CALL cp_fm_create(smo,mo_coeff%matrix_struct,error=error)
    CALL cp_dbcsr_sm_fm_multiply(matrix_S,mo_coeff,smo,ncol_global,error=error)

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    ! 4.
    !-----------------------------------------------------------------------------
    ALLOCATE (nstates(2),STAT=isos)
    IF (isos /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "nstates",int_size*2)

    CALL cp_fm_create(storage,mo_localized%matrix_struct,name='storage',error=error)

    DO imol = 1, SIZE(molecule_set)
       IF (ASSOCIATED(molecule_set(imol)%lmi)) THEN
          nstates(1) = molecule_set(imol)%lmi%nstates
       ELSE
          nstates(1) = 0
       ENDIF
       nstates(2) = para_env%mepos

       CALL mp_maxloc(nstates,para_env%group)

       IF (nstates(1)==0) CYCLE
       NULLIFY(states)
       ALLOCATE(states(nstates(1)),STAT=isos)
       IF (isos /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "states",int_size*nstates(1))
       states(:) = 0

       iproc=nstates(2)
       IF(iproc == para_env%mepos) THEN
          states(:) = molecule_set(imol)%lmi%states(:)
       END IF
       !!BCAST from here root = iproc
       CALL mp_bcast(states,iproc,para_env%group)

       ns = nstates(1)
       ind => states( : )
       ALLOCATE ( evals ( ns ), STAT = isos )
       IF (isos /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                       "evals",dp_size*ns)

       NULLIFY(b,c,d,fm_struct_tmp,e_vectors,rot_e_vectors)

       CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nrow_global, &
            ncol_global=ns, &
            para_env=mo_localized%matrix_struct%para_env, &
            context=mo_localized%matrix_struct%context,error=error)

       CALL cp_fm_create(b,fm_struct_tmp, name="b",error=error)
       CALL cp_fm_create(c,fm_struct_tmp, name="c",error=error)
       CALL cp_fm_create(rot_e_vectors,fm_struct_tmp, name="rot_e_vectors",error=error)

       CALL cp_fm_struct_release(fm_struct_tmp,error=error)

       CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ns,  ncol_global=ns, &
            para_env=mo_localized%matrix_struct%para_env, &
            context=mo_localized%matrix_struct%context,error=error)

       CALL cp_fm_create(d,fm_struct_tmp, name="d",error=error)
       CALL cp_fm_create(e_vectors,fm_struct_tmp, name="e_vectors",error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)

       DO i=1,ns
          CALL cp_fm_to_fm ( mo_localized, b, 1, ind ( i ), i)
       END DO

       CALL cp_dbcsr_sm_fm_multiply(Hks,b,c,ns,error=error)

       CALL cp_fm_gemm('T','N',ns,ns,nrow_global,1.0_dp, &
            b,c,0.0_dp,d,error=error)

       CALL choose_eigv_solver( d, e_vectors, evals , error=error)

       IF (output_unit>0) WRITE(output_unit,*)""
       IF (output_unit>0) WRITE(output_unit,*)"MOLECULE ",imol
       IF (output_unit>0) WRITE(output_unit,*)"NUMBER OF STATES  ", ns
       IF (output_unit>0) WRITE(output_unit,*)"EIGENVALUES"
       IF (output_unit>0) WRITE(output_unit,*)""
       IF (output_unit>0) WRITE(output_unit,*)"ENERGY      original MO-index"

       DO k=1,ns
          IF(ASSOCIATED(mark_states))THEN
             DO i=1,n_rep
                IF(imol==mark_states(1,i).AND.k==mark_states(2,i))marked_states(i)=ind(k)
             END DO
          END IF
          IF (output_unit>0) WRITE(output_unit,*) evals(k), ind(k)
       END DO
       IF (unit_report>0) THEN
          WRITE(unit_report,*) imol,ns, " imol, number of states"
          DO k=1,ns
             WRITE(unit_report,*) evals(k)
          ENDDO
       ENDIF

       CALL cp_fm_gemm('N','N',nrow_global,ns,ns,1.0_dp, &
            b,e_vectors,0.0_dp,rot_e_vectors,error=error)

       DO i=1,ns
          CALL cp_fm_to_fm ( rot_e_vectors, storage, 1, i, ind ( i ))
       END DO

       IF (.FALSE.) THEN ! this is too much data for large systems
         ! compute Eq. 6 from P. Hunt et al. (CPL 376, p. 68-74)
         CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ns, &
              ncol_global=ncol_global, &
              para_env=mo_localized%matrix_struct%para_env, &
              context=mo_localized%matrix_struct%context,error=error)
         CALL cp_fm_create(D_igk,fm_struct_tmp,error=error)
         CALL cp_fm_struct_release(fm_struct_tmp,error=error)
         CALL cp_fm_gemm('T','N',ns,ncol_global,nrow_global,1.0_dp, &
              rot_e_vectors,smo,0.0_dp,D_igk,error=error)
            DO i=1,ns
               DO k=1,ncol_global
                  CALL cp_fm_get_element(D_igk,i,k,tmp)
                  IF (unit_report>0) THEN
                     WRITE(unit_report,*) tmp**2
                  ENDIF
               ENDDO
            ENDDO
         CALL cp_fm_release(D_igk,error=error)
       ENDIF

       IF ( BTEST(cp_print_key_should_output(logger%iter_info,loc_print_section,&
            "MOLECULAR_STATES%CUBES",error=error),cp_p_file) ) THEN

          CALL get_qs_env(qs_env=qs_env,&
               atomic_kind_set=atomic_kind_set,&
               cell=cell,&
               dft_control=dft_control,&
               particle_set=particle_set,&
               pw_env=pw_env,error=error)

          DO i=1,ns
             IF (evals(i)<eval_range(1) .OR. evals(i)>eval_range(2)) CYCLE

             CALL calculate_wavefunction(rot_e_vectors,i,wf_r, &
                  wf_g, atomic_kind_set,cell,dft_control,particle_set, &
                  pw_env, error = error )

             WRITE(filename,'(a9,I4.4,a1,I5.5,a4)')"MOLECULE_",imol,"_",i,tag
             WRITE(title,'(A,I0,A,I0,A,F14.6,A,I0)') "Mol. Eigenstate ",i," of ",ns," E [a.u.] = ",&
                                                    evals(i)," Orig. index ",ind(i)
             unit_nr=cp_print_key_unit_nr(logger,loc_print_section,"MOLECULAR_STATES%CUBES",&
                  extension=".cube",middle_name=TRIM(filename),log_filename=.FALSE.,error=error)
             CALL pw_to_cube(wf_r%pw,unit_nr,particles=particles,title=title,&
                  stride=section_get_ivals(loc_print_section,&
                  "MOLECULAR_STATES%CUBES%STRIDE",error=error), error=error)
             CALL cp_print_key_finished_output(unit_nr,logger,loc_print_section,&
                  "MOLECULAR_STATES%CUBES",error=error)
          END DO
       ENDIF

       DEALLOCATE (evals,STAT=isos)
       IF (isos /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"eval")
       CALL cp_fm_release ( b ,error=error)
       CALL cp_fm_release ( c ,error=error)
       CALL cp_fm_release ( d ,error=error)
       CALL cp_fm_release ( e_vectors ,error=error)
       CALL cp_fm_release ( rot_e_vectors ,error=error)

       DEALLOCATE (states,STAT=isos)
       IF (isos /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"states")

    END DO
    CALL cp_fm_release(smo,error=error)
    CALL cp_fm_to_fm(storage,mo_localized,error)
    CALL cp_fm_release(storage,error)
    IF (ASSOCIATED(mark_states))THEN
       DEALLOCATE (mark_states,stat=isos)
       IF (isos /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"mark_states")
    END IF
    DEALLOCATE (nstates,STAT=isos)
    IF (isos /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"nstates")
    CALL cp_print_key_finished_output(unit_report,logger,loc_print_section,&
                                      "MOLECULAR_STATES",error=error)
    !------------------------------------------------------------------------------

    IF (.NOT.explicit) THEN
      DEALLOCATE(eval_range)
    END IF

    CALL timestop(handle)

  END SUBROUTINE construct_molecular_states

END MODULE molecular_states
