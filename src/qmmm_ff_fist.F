!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/qmmm_ff_fist [1.0] *
!!
!!   NAME
!!     qmmm_ff_fist
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     teo
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************


MODULE qmmm_ff_fist
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp,&
                                             default_string_length
  IMPLICIT NONE
  PRIVATE
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "qmmm_ff_fist"

  PUBLIC :: qmmm_ff_precond_only_qm

CONTAINS
  
  !!****f* qmmm_ff_fist/qmmm_ff_precond_only_qm [1.0] *
  !!
  !!   NAME
  !!     qmmm_ff_precond_only_qm
  !!
  !!   FUNCTION
  !!
  !!     This function handles the atom names and modifies
  !!     the "_QM_" prefix, in order to find the parameters
  !!     and return .TRUE. if all input labels belong to QM atoms
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     11.2004 created [tlaino]
  !!
  !!*************************************************************************
  FUNCTION qmmm_ff_precond_only_qm(id1, id2, id3, id4) RESULT(only_qm)
    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=default_string_length), INTENT(INOUT) :: id1
    CHARACTER(LEN=default_string_length), INTENT(INOUT), OPTIONAL :: id2, id3, id4
    LOGICAL :: only_qm
    ! Local Variables
    CHARACTER(LEN=default_string_length) :: tmp
    INTEGER :: j

    only_qm = .FALSE.
    IF (INDEX(id1,"_QM_") /= 0 ) THEN
       only_qm = .TRUE.
       tmp=TRIM(id1(5:))
       clean_string_1: DO j = 1, default_string_length
          id1(j:j) = " "
       END DO clean_string_1
       id1 = TRIM(tmp)
    END IF
    
    IF (PRESENT(id2)) THEN
       IF (INDEX(id2,"_QM_") /= 0 ) THEN
          tmp=TRIM(id2(5:))
          clean_string_2: DO j = 1, default_string_length
             id2(j:j) = " "
          END DO clean_string_2
          id2 = TRIM(tmp)
       ELSE
          only_qm = .FALSE.
       END IF
    END IF

    IF (PRESENT(id3)) THEN
       IF (INDEX(id3,"_QM_") /= 0 ) THEN
          tmp=TRIM(id3(5:))
          clean_string_3: DO j = 1, default_string_length
             id3(j:j) = " "
          END DO clean_string_3
          id3 = TRIM(tmp)
       ELSE
          only_qm = .FALSE.
       END IF
    END IF

    IF (PRESENT(id4)) THEN
       IF (INDEX(id4,"_QM_") /= 0 ) THEN
          tmp=TRIM(id4(5:))
          clean_string_4: DO j = 1, default_string_length
             id4(j:j) = " "
          END DO clean_string_4
          id4 = TRIM(tmp)
       ELSE
          only_qm = .FALSE.
       END IF
    END IF

   END FUNCTION qmmm_ff_precond_only_qm

END MODULE qmmm_ff_fist
