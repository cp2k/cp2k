!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for the calculation of wannier states
!> \author Alin M Elena
! *****************************************************************************
MODULE wannier_states
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_element,&
                                             cp_fm_get_info,&
                                             cp_fm_release,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE molecule_types_new,              ONLY: molecule_type
  USE qs_environment_types,            ONLY: qs_environment_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE wannier_states_types,            ONLY: wannier_centres_type
!!!! this ones are needed to mapping
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'wannier_states'

  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.TRUE.

! *** Public subroutines ***

  PUBLIC :: construct_wannier_states

CONTAINS

! *****************************************************************************
!> \brief constructs wannier states. mo_localized should not be overwritten!
! *****************************************************************************
  SUBROUTINE construct_wannier_states(molecule_set, mo_localized, &
       Hks, qs_env,loc_print_section,WannierCentres,ns,states,error)

    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(cp_fm_type), INTENT(in), POINTER    :: mo_localized
    TYPE(cp_dbcsr_type), POINTER             :: Hks
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: loc_print_section
    TYPE(wannier_centres_type), &
      INTENT(INOUT)                          :: WannierCentres
    INTEGER, INTENT(IN)                      :: ns
    INTEGER, INTENT(IN), POINTER             :: states(:)
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'construct_wannier_states', &
      routineP = moduleN//':'//routineN

    CHARACTER(default_string_length)         :: unit_str
    INTEGER                                  :: handle, i, iproc, &
                                                ncol_global, nproc, &
                                                nrow_global, nstates(2), &
                                                output_unit, unit_mat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: unit_conv
    TYPE(cp_error_type)                      :: suberror
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: b, c, d
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: print_key

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

    CALL timeset(routineN,handle)
    failure = .FALSE.

    para_env => qs_env % para_env
    nproc = para_env%num_pe

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    output_unit = cp_logger_get_default_io_unit(logger)
    CALL cp_fm_get_info(mo_localized, &
                        ncol_global=ncol_global, &
                        nrow_global=nrow_global ,error=error)

    nstates(1)= ns
    nstates(2) = para_env%mepos
    iproc=nstates(2)
    NULLIFY(b,c,d,fm_struct_tmp,print_key)

    print_key => section_vals_get_subs_vals(loc_print_section,"WANNIER_CENTERS",error=error)
    CALL section_vals_val_get(print_key,"UNIT",c_val=unit_str,error=error)
    unit_conv = cp_unit_from_cp2k(1.0_dp,TRIM(unit_str),error=error)

    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=nrow_global, &
      ncol_global=1, &
      para_env=mo_localized%matrix_struct%para_env, &
      context=mo_localized%matrix_struct%context,error=error)

    CALL cp_fm_create(b,fm_struct_tmp, name="b",error=error)
    CALL cp_fm_create(c,fm_struct_tmp, name="c",error=error)

    CALL cp_fm_struct_release(fm_struct_tmp,error=error)


    CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=1,  ncol_global=1, &
      para_env=mo_localized%matrix_struct%para_env, &
      context=mo_localized%matrix_struct%context,error=error)

    CALL cp_fm_create(d,fm_struct_tmp, name="d",error=error)
    CALL cp_fm_struct_release(fm_struct_tmp,error=error)

    WannierCentres%WannierHamDiag=0.0_dp
      ! try to print the matrix

    CALL cp_error_init(suberror,template_error=error)
       unit_mat=cp_print_key_unit_nr(logger,loc_print_section,&
       "WANNIER_STATES",extension=".whks",&
       ignore_should_output=.FALSE.,error=error)
    IF (unit_mat > 0) THEN
      WRITE(unit_mat,'(a16,1(i0,1x))')"Wannier states: ",ns
      WRITE(unit_mat,'(a16)')"#No x y z energy "
    ENDIF
    DO i=1,ns
      CALL cp_fm_to_fm ( mo_localized, b, 1, states ( i ), 1)
      CALL cp_dbcsr_sm_fm_multiply(Hks,b,c,1,error=error)
      CALL cp_fm_gemm('T','N',1,1,nrow_global,1.0_dp, &
          b,c,0.0_dp,d,error=error)
      CALL cp_fm_get_element(d,1,1,WannierCentres%WannierHamDiag(i))
  !               if (iproc==para_env%mepos) WRITE(unit_mat,'(f16.8,2x)', advance='no')WannierCentres%WannierHamDiag(i)
      IF (unit_mat > 0)WRITE(unit_mat,'(i0,1x,4(f16.8,2x))')states(i),&
            WannierCentres%centres(1:3,states(i))*unit_conv,WannierCentres%WannierHamDiag(states(i))
    END DO

    IF (unit_mat > 0) WRITE(unit_mat,*)
    CALL cp_error_reset(suberror)
    CALL cp_print_key_finished_output(unit_mat,logger,loc_print_section,&
                "WANNIER_STATES",error=error)
    CALL cp_error_dealloc_ref(suberror)
    IF (output_unit>0) THEN
      WRITE(output_unit,*)""
      WRITE(output_unit,*)"NUMBER OF Wannier STATES  ", ns
      WRITE(output_unit,*)"ENERGY      original MO-index"
      DO i=1,ns
        WRITE(output_unit,'(f16.8,2x,i0)')WannierCentres%WannierHamDiag(i),states(i)
      ENDDO
    ENDIF
    CALL cp_fm_release ( b ,error=error)
    CALL cp_fm_release ( c ,error=error)
    CALL cp_fm_release ( d ,error=error)
    CALL timestop(handle)
  END SUBROUTINE construct_wannier_states

END MODULE wannier_states

