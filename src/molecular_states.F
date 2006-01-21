!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!


!!****** cp2k/molecular_states [1.0] *
!!
!!   NAME
!!     molecular_states
!!
!!   FUNCTION
!!     Routines for the calculation of molecular states
!!
!!   AUTHOR
!!     Manu Manu (transcribed by CJM) 
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE molecular_states
! *****************************************************************************
  USE coefficient_types,               ONLY: coeff_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_sm_fm_interactions,           ONLY: cp_sm_fm_multiply
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_maxloc
  USE molecule_types_new,              ONLY: molecule_type
  USE pw_env_types,                    ONLY: pw_env_type
  USE qs_collocate_density,            ONLY: calculate_wavefunction
  USE qs_environment_types,            ONLY: qs_environment_type
  USE realspace_grid_types,            ONLY: rs_pw_to_cube
  USE sparse_matrix_types,             ONLY: real_matrix_type
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr, &
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_p_file
  USE input_section_types,             ONLY: section_vals_type, &
                                             section_vals_val_get,&
                                             section_vals_get_subs_vals, &
                                             section_vals_get_subs_vals2, &
                                             section_get_rval,&
                                             section_get_ival,&
                                             section_get_cval,&
                                             section_vals_get
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "molecular_states"

  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

! *** Public subroutines ***

  PUBLIC :: construct_molecular_states

! *****************************************************************************

CONTAINS

! *****************************************************************************

SUBROUTINE construct_molecular_states(molecule_set, mo_localized, &
                                   mo_coeff, mo_eigenvalues, &
                                   Hks,matrix_S,globenv,qs_env,wf_r,wf_g,&
                                   force_env_section, error)


    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(cp_fm_type), POINTER                :: mo_localized, mo_coeff
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mo_eigenvalues
    TYPE(real_matrix_type), POINTER          :: Hks, matrix_S
    TYPE(global_environment_type), &
      pointer                                :: globenv
    TYPE(section_vals_type), POINTER :: force_env_section
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(coeff_type), INTENT(INOUT)          :: wf_r, wf_g

    CHARACTER(LEN=default_string_length)     :: filename
    INTEGER                                  :: i, imol, iproc, isos, j, k, &
                                                ncol_global, nproc, &
                                                nrow_global, ns, unit_report
    INTEGER, DIMENSION(:), POINTER           :: ind
    INTEGER, POINTER                         :: nstates(:), states(:)
    REAL(KIND=dp)                            :: tmp
    REAL(KIND=dp), ALLOCATABLE               :: evals( : )
    TYPE(cp_error_type), OPTIONAL, INTENT(INOUT) :: error
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: b, c, d, D_igk, e_vectors, &
                                                rot_e_vectors, smo
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_logger_type), POINTER            :: logger
    LOGICAL :: failure
    INTEGER :: unit_nr

    failure = .FALSE.
    NULLIFY(logger)
    logger => cp_error_get_logger(error)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
! 1. 
!-----------------------------------------------------------------------------

  para_env => qs_env % para_env
  nproc = para_env%num_pe
 

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 2. 
  !-----------------------------------------------------------------------------
  IF (globenv%ionode) THEN 
     unit_report=171
     OPEN(unit_report,FILE="Molecular_DOS")
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
                      nrow_global=nrow_global )
  NULLIFY(smo)
  CALL cp_fm_create(smo,mo_coeff%matrix_struct)
  CALL cp_sm_fm_multiply(matrix_S,mo_coeff,smo,ncol_global)


  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! 4. 
  !-----------------------------------------------------------------------------
  ALLOCATE(nstates(2),STAT=isos)
  IF(isos/=0) CALL stop_memory('construct_molecular_states','nstates')

  DO imol = 1, SIZE(molecule_set)
    nstates(1) = molecule_set(imol)%lmi%nstates
    nstates(2) = para_env%mepos

    CALL mp_maxloc(nstates,para_env%group)

    NULLIFY(states)
    ALLOCATE(states(nstates(1)),STAT=isos)
    IF(isos/=0) CALL stop_memory('construct_molecular_states','states')
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
    IF ( isos /= 0 ) CALL stop_memory ( "construct_molecular_states", &
                                        "eval", ns )

    NULLIFY(b,c,d,fm_struct_tmp,e_vectors,rot_e_vectors)

          
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nrow_global, &
                             ncol_global=ns, &
                             para_env=mo_localized%matrix_struct%para_env, &
                             context=mo_localized%matrix_struct%context)

    CALL cp_fm_create(b,fm_struct_tmp, name="b")
    CALL cp_fm_create(c,fm_struct_tmp, name="c")
    CALL cp_fm_create(rot_e_vectors,fm_struct_tmp, name="rot_e_vectors")

    CALL cp_fm_struct_release(fm_struct_tmp)

    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ns,  ncol_global=ns, &
                             para_env=mo_localized%matrix_struct%para_env, &
                             context=mo_localized%matrix_struct%context)

    CALL cp_fm_create(d,fm_struct_tmp, name="d")
    CALL cp_fm_create(e_vectors,fm_struct_tmp, name="e_vectors")
    CALL cp_fm_struct_release(fm_struct_tmp)
                
    DO i=1,ns
      CALL cp_fm_to_fm ( mo_localized, b, 1, ind ( i ), i)
    END DO

    CALL cp_sm_fm_multiply(Hks,b,c,ns)

    CALL cp_fm_gemm('T','N',ns,ns,nrow_global,1.0_dp, &             
                         b,c,0.0_dp,d)    

    CALL cp_fm_syevd ( d, e_vectors, evals )

    IF (globenv%ionode) WRITE(globenv%scr,*)"" 
    IF (globenv%ionode) WRITE(globenv%scr,*)"MOLECULE ",imol
    IF (globenv%ionode) WRITE(globenv%scr,*)"NUMBER OF STATES  ", ns
    IF (globenv%ionode) WRITE(globenv%scr,*)"EIGENVALUES"
    IF (globenv%ionode) WRITE(globenv%scr,*)"" 

    DO k=1,ns
      IF (globenv%ionode) WRITE(globenv%scr,*) evals(k)
    END DO
    IF (globenv%ionode) THEN
       WRITE(unit_report,*) imol,ns, " imol, number of states"
       DO k=1,ns
          WRITE(unit_report,*) evals(k)
       ENDDO
    ENDIF
                  
    CALL cp_fm_gemm('N','N',nrow_global,ns,ns,1.0_dp, &
                    b,e_vectors,0.0_dp,rot_e_vectors)

    DEALLOCATE ( evals,  STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( "construct_molecular_states", &
                                        "eval", ns )
    
    ! compute Eq. 6 from P. Hunt et al. (CPL 376, p. 68-74)
    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ns, &
                             ncol_global=ncol_global, &
                             para_env=mo_localized%matrix_struct%para_env, &
                             context=mo_localized%matrix_struct%context)
    CALL cp_fm_create(D_igk,fm_struct_tmp)
    CALL cp_fm_struct_release(fm_struct_tmp) 
    CALL cp_fm_gemm('T','N',ns,ncol_global,nrow_global,1.0_dp, &             
                         rot_e_vectors,smo,0.0_dp,D_igk)    
    DO i=1,ns
       DO k=1,ncol_global
          CALL cp_fm_get_element(D_igk,i,k,tmp)
          IF (globenv%ionode) THEN
             write(unit_report,*) tmp**2
          ENDIF
       ENDDO
    ENDDO
    CALL cp_fm_release(D_igk)

    IF ( BTEST(cp_print_key_should_output(logger%iter_info,force_env_section,&
         "DFT%PRINT%MOLECULAR_STATES%CUBES",error=error),cp_p_file) ) THEN
       DO i=1,ns
          j=ind(i)
          CALL calculate_wavefunction(rot_e_vectors,i,wf_r, &
                                      wf_g, qs_env, error = error )

          WRITE(filename,'(a9,I2.2,a1,I3.3)')"MOLECULE_",imol,"_",j          
          unit_nr=cp_print_key_unit_nr(logger,force_env_section,&
               "DFT%PRINT%MOLECULAR_STATES%CUBES",&
               extension=".cube",middle_name=TRIM(filename),log_filename=.FALSE.,error=error)
          CALL rs_pw_to_cube(wf_r%pw,unit_nr,globenv%ionode,filename,&
               stride=section_get_ival(force_env_section,&
               "DFT%PRINT%MOLECULAR_STATES%CUBES%STRIDE"),force_env_section=qs_env%input,&
               error=error)
          CALL cp_print_key_finished_output(unit_nr,logger,force_env_section,&
               "DFT%PRINT%MOLECULAR_STATES%CUBES",error=error)
          IF (globenv%ionode) WRITE(globenv%scr,*) filename
       END DO
    ENDIF

    CALL cp_fm_release ( b )
    CALL cp_fm_release ( c )
    CALL cp_fm_release ( d )
    CALL cp_fm_release ( e_vectors )
    CALL cp_fm_release ( rot_e_vectors )

    DEALLOCATE(states,STAT=isos)
    IF(isos/=0) CALL stop_memory('construct_molecular_states','states')

  END DO 
  CALL cp_fm_release(smo)
  DEALLOCATE(nstates,STAT=isos)
  IF(isos/=0) CALL stop_memory('construct_molecular_states','nstates')
  IF (globenv%ionode) CLOSE(unit_report)
  !CALL stop_program("molecular_states","checkpoint")
!------------------------------------------------------------------------------
END SUBROUTINE construct_molecular_states

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

END MODULE molecular_states
