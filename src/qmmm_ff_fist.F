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
  FUNCTION qmmm_ff_precond_only_qm(id1, id2, id3, id4, is_link) RESULT(only_qm)
    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=default_string_length), INTENT(INOUT) :: id1
    CHARACTER(LEN=default_string_length), INTENT(INOUT), OPTIONAL :: id2, id3, id4
    LOGICAL, OPTIONAL  :: is_link
    LOGICAL :: only_qm, my_link
    ! Local Variables
    CHARACTER(LEN=default_string_length) :: tmp
    INTEGER :: j, my_index

    only_qm = .FALSE.
    my_link = .FALSE.
    IF (INDEX(id1,"_QM_") /= 0 ) THEN
       my_link = qmmm_ff_precond_only_link(id1).OR.my_link
       my_index=INDEX(id1,"_QM_")+LEN_TRIM("_QM_")
       only_qm = .TRUE.
       tmp=TRIM(id1(my_index:))
       clean_string_1: DO j = 1, default_string_length
          id1(j:j) = " "
       END DO clean_string_1
       id1 = TRIM(tmp)
    END IF
    
    IF (PRESENT(id2)) THEN
       IF (INDEX(id2,"_QM_") /= 0 ) THEN
          my_link = qmmm_ff_precond_only_link(id2).OR.my_link
          my_index=INDEX(id2,"_QM_")+LEN_TRIM("_QM_")
          tmp=TRIM(id2(my_index:))
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
          my_link = qmmm_ff_precond_only_link(id3).OR.my_link
          my_index=INDEX(id3,"_QM_")+LEN_TRIM("_QM_")
          tmp=TRIM(id3(my_index:))
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
          my_link = qmmm_ff_precond_only_link(id4).OR.my_link
          my_index=INDEX(id4,"_QM_")+LEN_TRIM("_QM_")
          tmp=TRIM(id4(my_index:))
          clean_string_4: DO j = 1, default_string_length
             id4(j:j) = " "
          END DO clean_string_4
          id4 = TRIM(tmp)
       ELSE
          only_qm = .FALSE.
       END IF
    END IF

    IF (PRESENT(is_link)) is_link = my_link

   END FUNCTION qmmm_ff_precond_only_qm

   
   FUNCTION qmmm_ff_precond_only_link(id1) RESULT(is_link)
     IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=default_string_length), INTENT(INOUT) :: id1
    LOGICAL :: is_link
    ! Local Variables
    CHARACTER(LEN=default_string_length) :: tmp
    INTEGER :: j, my_index
    
    is_link = .FALSE.
    IF (INDEX(id1,"_LNK") /= 0 ) THEN
       my_index=INDEX(id1,"_LNK")+LEN_TRIM("_LNK")
       is_link = .TRUE.
       tmp=TRIM(id1(my_index:))
       clean_string_1: DO j = 1, default_string_length
          id1(j:j) = " "
       END DO clean_string_1
       id1 = TRIM(tmp)
    END IF

   END FUNCTION qmmm_ff_precond_only_link


END MODULE qmmm_ff_fist
