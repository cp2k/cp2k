!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_fm_vect [1.0] *
!!
!!   NAME
!!     cp_fm_vect
!!
!!   FUNCTION
!!     routine to handle vectors of full matrixes
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_fm_vect
  USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                             cp_blacs_env_from_globenv,&
                                             cp_blacs_env_release,&
                                             cp_blacs_env_retain,&
                                             finish_blacs,&
                                             start_blacs
  USE cp_error_handling,               ONLY: cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_fm_struct,                    ONLY: cp_f_matrix_struct_type,&
                                             cp_fmstruct_create,&
                                             cp_fmstruct_equivalent,&
                                             cp_fmstruct_release,&
                                             cp_fmstruct_retain,&
                                             optimal_blacs_col_block_size,&
                                             optimal_blacs_row_block_size
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_create2,&
                                             cp_fm_get_info,&
                                             cp_fm_get_submatrix,&
                                             cp_fm_init_random,&
                                             cp_fm_release,&
                                             cp_fm_retain,&
                                             cp_fm_set_all,&
                                             cp_fm_set_element,&
                                             cp_fm_set_submatrix,&
                                             cp_full_matrix_p_type,&
                                             cp_full_matrix_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_env,                     ONLY: cp_para_env_from_globenv,&
                                             cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: wp=>dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_fm_vect'

  PUBLIC :: cp_fm_vect_create2, cp_fm_vect_dealloc, cp_fm_vect_copy
!***
!****************************************************************************
CONTAINS

!!****f* cp_fm_vect/cp_fm_vect_create2 [1.0] *
!!
!!   NAME
!!     cp_fm_vect_create2
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_vect_create2(matrixes, n_matrixes, nrow_global,&
!!         ncol_global, globenv, name, nrow_block, ncol_block, error)
!!       Type(cp_full_matrix_p_type), Dimension(:), Pointer:: matrixes
!!       Type(global_environment_type), Intent (IN):: globenv
!!       Character(Len=*), Intent (IN):: name
!!       Integer, Intent (IN):: n_matrixes, nrow_global, ncol_global
!!       Integer, Intent (IN), Optional:: nrow_block, ncol_block
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fm_vect_create2
!!
!!   FUNCTION
!!     allocates an array of identical blacs matrixes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrixes: a pointer that will point to the array of matrixes
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_fm_vect_create2(matrixes,n_matrixes,&
       nrow_global,ncol_global,globenv,name,&
       nrow_block,ncol_block,error)
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: matrixes
    TYPE(global_environment_type), INTENT(in):: globenv
    CHARACTER(len=*), INTENT(in) :: name
    INTEGER, INTENT(in) :: n_matrixes, nrow_global, ncol_global
    INTEGER, INTENT(in), OPTIONAL :: nrow_block, ncol_block
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_vect_create2',&
         routineP=moduleN//':'//routineN
    INTEGER :: i, stat, row_bl, col_bl
    TYPE(cp_f_matrix_struct_type), POINTER :: matrix_struct
    TYPE(cp_para_env_type), POINTER :: para_env
    TYPE(cp_blacs_env_type), POINTER :: blacs_env
    failure=.FALSE.
    para_env => cp_para_env_from_globenv(globenv)
    blacs_env => cp_blacs_env_from_globenv(globenv)
    row_bl=optimal_blacs_row_block_size
    col_bl=optimal_blacs_col_block_size
    IF (PRESENT(nrow_block)) row_bl=nrow_block
    IF (PRESENT(ncol_block)) col_bl=ncol_block

    ALLOCATE(matrixes(n_matrixes),stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL cp_fmstruct_create(matrix_struct,nrow_global=nrow_global,&
            ncol_global=ncol_global,nrow_block=row_bl,&
            ncol_block=col_bl, context=blacs_env,&
            para_env=para_env)
       DO i=1,n_matrixes
          CALL cp_fm_create(matrixes(i)%matrix,&
               matrix_struct=matrix_struct, name=name//cp_to_string(i),&
               error=error)
       END DO
       CALL cp_fmstruct_release(matrix_struct)
    END IF
  END SUBROUTINE cp_fm_vect_create2
!***************************************************************************

!!****f* cp_fm_vect/cp_fm_vect_dealloc [1.0] *
!!
!!   NAME
!!     cp_fm_vect_dealloc
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_vect_dealloc(matrixes, error)
!!       Type(cp_full_matrix_p_type), Dimension(:), Pointer:: matrixes
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine cp_fm_vect_dealloc
!!
!!   FUNCTION
!!     deallocate an array of pointers to blacs matrixes
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrixes: the array of matrixes to deallocate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE cp_fm_vect_dealloc(matrixes, error)
    TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: matrixes
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_fm_vect_dealloc',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat, i
    failure=.FALSE.

    IF (ASSOCIATED(matrixes)) THEN
       DO i=1,SIZE(matrixes)
          CALL cp_fm_release(matrixes(i)%matrix,error=error)
       END DO
       DEALLOCATE(matrixes,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
  END SUBROUTINE cp_fm_vect_dealloc
!***************************************************************************

!!****f* cp_fm_vect/cp_fm_vect_copy[1.0] *
!!
!!   NAME
!!     cp_fm_vect_copy
!!
!!   SYNOPSIS
!!     Subroutine cp_fm_vect_copy(matrixes, copy, error)
!!       Type(cp_full_matrix_p_type), Dimension(:), Intent (IN):: matrixes
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!       Type(cp_full_matrix_p_type), Dimension(:), Pointer:: copy
!!     End Subroutine cp_fm_vect_copy
!!
!!   FUNCTION
!!     Does a shallow copy of an array of full matrices (i.e. just retains 
!!     the matrices)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - matrixes: the matrixes to copy
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE cp_fm_vect_copy(matrixes, copy, error)
  TYPE(cp_full_matrix_p_type), DIMENSION(:), INTENT(in) :: matrixes
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  TYPE(cp_full_matrix_p_type), DIMENSION(:), POINTER :: copy
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='cp_fm_vect_copy',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat,i
  failure=.FALSE.
  
  ALLOCATE(copy(SIZE(matrixes)),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i=1,SIZE(matrixes)
        copy(i)%matrix => matrixes(i)%matrix
        CALL cp_fm_retain(matrixes(i)%matrix,error=error)
     END DO
  END IF
END SUBROUTINE cp_fm_vect_copy
!***************************************************************************

END MODULE cp_fm_vect
