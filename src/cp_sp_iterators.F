!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/cp_sp_iterators [1.0] *
!!
!!   NAME
!!     cp_sp_iterators
!!
!!   FUNCTION
!!     iterators on the sparse block matrix
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE cp_sp_iterators
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_sp_iterators'

  PUBLIC :: cp_sp_local_b_iter_type, cp_sp_owned_b_iter_type,&
       cp_sp_cached_b_iter_type, cp_sp_g_block_iterator
  PUBLIC :: cp_init, cp_dealloc_ref, cp_next
  PUBLIC :: cp_sp_local_b_iter_init, cp_sp_local_b_iter_dealloc_ref,&
       cp_sp_local_b_iter_next, cp_sp_owned_b_iter_init,&
       cp_sp_owned_b_iter_dealloc_ref, cp_sp_cached_b_iter_init,&
       cp_sp_owned_b_iter_dealloc_ref

  INTERFACE cp_init
     MODULE PROCEDURE cp_sp_local_b_iter_init
  END INTERFACE
  INTERFACE cp_dealloc_ref
     MODULE PROCEDURE cp_sp_local_b_iter_dealloc_ref
  END INTERFACE
  INTERFACE cp_next
     MODULE PROCEDURE cp_sp_local_b_iter_next
  END INTERFACE
  INTERFACE cp_init
     MODULE PROCEDURE cp_sp_owned_b_iter_init
  END INTERFACE
  INTERFACE cp_dealloc_ref
     MODULE PROCEDURE cp_sp_owned_b_iter_dealloc_ref
  END INTERFACE
  INTERFACE cp_init
     MODULE PROCEDURE cp_sp_cached_b_iter_init
  END INTERFACE
  INTERFACE cp_dealloc_ref
     MODULE PROCEDURE cp_sp_owned_b_iter_dealloc_ref
  END INTERFACE

  !***
  !****************************************************************************
  !!****s* cp_block_matrix_type/cp_sp_local_b_iter_type [1.0] *
  !!
  !!   NAME
  !!     cp_sp_local_b_iter_type
  !!
  !!   FUNCTION
  !!     iterates on the local blocks
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     @author Fawzi Mohamed
  !!     @version 3.2002
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE cp_sp_local_b_iter_type
     LOGICAL :: initialized, all_blocks
     TYPE(cp_block_matrix_type), POINTER :: matrix
     TYPE(cp_matrix_block_iterator) :: iterator
  END TYPE cp_sp_local_b_iter_type
  !!***
  !****************************************************************************

  !!****s* cp_block_matrix_type/cp_sp_owned_b_iter_type [1.0] *
  !!
  !!   NAME
  !!     cp_sp_owned_b_iter_type
  !!
  !!   FUNCTION
  !!     iterates on the owned blocks
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     @author Fawzi Mohamed
  !!     @version 3.2002
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE cp_sp_owned_b_iter_type
     TYPE(cp_sp_local_b_iter_type), POINTER :: local_iterator
  END TYPE cp_sp_owned_b_iter_type
  !!***
  !****************************************************************************

  !!****s* cp_block_matrix_type/cp_sp_cached_b_iter_type [1.0] *
  !!
  !!   NAME
  !!     cp_sp_cached_b_iter_type
  !!
  !!   FUNCTION
  !!     iterates on the cached blocks
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     @author Fawzi Mohamed
  !!     @version 3.2002
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE cp_sp_cached_b_iter_type
     TYPE(cp_sp_local_b_iter_type), POINTER :: local_iterator
  END TYPE cp_sp_cached_b_iter_type
  !!***
  !****************************************************************************

  !!****s* cp_block_matrix_type/cp_sp_g_block_iterator [1.0] *
  !!
  !!   NAME
  !!     cp_sp_g_block_iterator
  !!
  !!   FUNCTION
  !!     iterates on all the non zero blocks of the matrix
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     -
  !!
  !!   AUTHOR
  !!     @author Fawzi Mohamed
  !!     @version 3.2002
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE cp_sp_g_block_iterator
     LOGICAL :: initialized
     TYPE(cp_block_matrix_type) :: matrix
     INTEGER :: i,j
  END TYPE cp_sp_g_block_iterator
  !!***
  !****************************************************************************

CONTAINS
  !=========== block iterators ============

  !==== local block iterator ====

  !!****f* cp_sp_iterators/cp_sp_local_b_iter_init [1.0] *
  !!
  !!   NAME
  !!     cp_sp_local_b_iter_init
  !!
  !!   FUNCTION
  !!     an iteratotor that goes trough all the non zero blocks local
  !!     blocks, both owned and cached.
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     iterator: the iterator to be initialized
  !!     matrix: the matrix to iterate on
  !!     error: variable to control error logging, stopping,... 
  !!            see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     @author Fawzi Mohamed
  !!     @version 3.2002
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*** **********************************************************************
  SUBROUTINE cp_sp_local_b_iter_init(iterator,matrix,error)
    TYPE(cp_sp_owned_b_iter_type), INTENT(out) :: iterator
    TYPE(cp_block_matrix_type), INTENT(in), TARGET :: matrix
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_sp_local_b_iter_init',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    IF (.NOT. failure) THEN

    END IF
  END SUBROUTINE cp_sp_local_b_iter_init
  !***************************************************************************

  !!****f* cp_sp_local_b_iter_dealloc_ref [1.0] *
  !!
  !!   NAME
  !!      cp_sp_local_b_iter_dealloc_ref
  !!
  !!   FUNCTION
  !!     dealloc the memory used by the iterator
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     iterator: the iterator to be released
  !!     error: variable to control error logging, stopping,... 
  !!            see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     @author Fawzi Mohamed
  !!     @version 3.2002
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*** **********************************************************************
  SUBROUTINE cp_sp_local_b_iter_dealloc_ref(iterator, error)
    TYPE(cp_sp_owned_b_iter_dealloc_ref), INTENT(inout) :: iterator
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_sp_local_b_iter_dealloc_ref',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    IF (.NOT. failure) THEN

    END IF
  END SUBROUTINE cp_sp_local_b_iter_dealloc_ref
  !***************************************************************************

  !!****f* cp_sp_iterators/cp_sp_local_b_iter_next [1.0] *
  !!
  !!   NAME
  !!     cp_sp_local_b_iter_next
  !!
  !!   FUNCTION
  !!     goes to the next local block, returns true if not at end
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     iterator: the iterator that sould move to the next block
  !!     block_handler: the handler to the actual block (to do some changes 
  !!            to the values you MUST release the handler with
  !!            call cp_release_block (block_handler,changed_values=.true.)
  !!     block_row: the actual block row
  !!     block_col: the actual block col
  !!     block_val: the actual block value (becomes invalid after
  !!                another call to next)
  !!     trans_to_do: the transposition still to do [remove ?]
  !!     local_scratch: if it is a local scratch block
  !!     owned: if it is an owned block
  !!     local_copy: if it is a local copy
  !!     error: variable to control error logging, stopping,... 
  !!            see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     @author Fawzi Mohamed
  !!     @version 3.2002
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*** **********************************************************************
  FUNCTION cp_sp_local_b_iter_next(iterator, block_row, block_col, block_val,&
       trans_to_do, local_scatch, owned, local_copy,error)RESULT(res)
    TYPE(cp_sp_local_b_iter_next), INTENT(inout) :: iterator
    INTEGER, INTENT(out), OPTIONAL :: block_row, block_col
    REAL(KIND = dp), DIMENSION(:,:), POINTER, OPTIONAL :: block_val
    CHARACTER(len=1), INTENT(out), OPTIONAL :: trans_to_do
    LOGICAL, INTENT(out), OPTIONAL :: local_scratch, owned, local_copy
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='cp_sp_local_b_iter_next',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    IF (.NOT. failure) THEN

    END IF
  END FUNCTION cp_sp_local_b_iter_next
  !***************************************************************************

  !=== owned blocks iterator ===

  !!****f* cp_sp_iterators/cp_sp_owned_b_iter_init [1.0] *
  !!
  !!   NAME
  !!     cp_sp_owned_b_iter_init
  !!
  !!   FUNCTION
  !!     an iteratotor that goes trough the non zero owned blocks.
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     iterator: the iterator to be initialized
  !!     matrix: the matrix to iterate on
  !!     error: variable to control error logging, stopping,... 
  !!            see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     @author Fawzi Mohamed
  !!     @version 3.2002
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*** **********************************************************************
  SUBROUTINE cp_sp_owned_b_iter_init(iterator,matrix,error)
    TYPE(cp_sp_owned_b_iter_type), INTENT(out) :: iterator
    TYPE(cp_block_matrix_type), INTENT(in), TARGET :: matrix
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN='cp_sp_owned_b_iter_init',&
         routineP=moduleN//':'//routineN

    CALL cp_sp_local_b_iter_init(iterator%local_iterator,&
         matrix=matrix,error=error)
  END SUBROUTINE cp_sp_owned_b_iter_init
  !***************************************************************************

  !!****f* cp_sp_owned_b_iter_dealloc_ref [1.0] *
  !!
  !!   NAME
  !!      cp_sp_owned_b_iter_dealloc_ref
  !!
  !!   FUNCTION
  !!     dealloc the memory used by the iterator
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     iterator: the iterator to be released
  !!     error: variable to control error logging, stopping,... 
  !!            see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     @author Fawzi Mohamed
  !!     @version 3.2002
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*** **********************************************************************
  SUBROUTINE cp_sp_owned_b_iter_dealloc_ref(iterator, error)
    TYPE(cp_sp_owned_b_iter_dealloc_ref), INTENT(inout) :: iterator
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN='cp_sp_owned_b_iter_dealloc_ref',&
         routineP=moduleN//':'//routineN

    CALL cp_sp_local_b_iter_dealloc_ref(iterator%local_iterator,&
         error=error)
  END SUBROUTINE cp_sp_owned_b_iter_dealloc_ref
  !***************************************************************************


  !=== cached blocks iterator ===

  !!****f* cp_sp_iterators/cp_sp_cached_b_iter_init [1.0] *
  !!
  !!   NAME
  !!     cp_sp_cached_b_iter_init
  !!
  !!   FUNCTION
  !!     an iteratotor that goes trough the non zero cached blocks.
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     iterator: the iterator to be initialized
  !!     matrix: the matrix to iterate on
  !!     error: variable to control error logging, stopping,... 
  !!            see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     @author Fawzi Mohamed
  !!     @version 3.2002
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*** **********************************************************************
  SUBROUTINE cp_sp_cached_b_iter_init(iterator,matrix,error)
    TYPE(cp_sp_cached_b_iter_type), INTENT(out) :: iterator
    TYPE(cp_block_matrix_type), INTENT(in), TARGET :: matrix
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN='cp_sp_cached_b_iter_init',&
         routineP=moduleN//':'//routineN

    CALL cp_sp_local_b_iter_init(iterator%local_iterator,&
         matrix=matrix,error=error)
  END SUBROUTINE cp_sp_cached_b_iter_init
  !***************************************************************************


  !!****f* cp_sp_owned_b_iter_dealloc_ref [1.0] *
  !!
  !!   NAME
  !!      cp_sp_owned_b_iter_dealloc_ref
  !!
  !!   FUNCTION
  !!     dealloc the memory used by the iterator
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     iterator: the iterator to be released
  !!     error: variable to control error logging, stopping,... 
  !!            see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     @author Fawzi Mohamed
  !!     @version 3.2002
  !!
  !!   MODIFICATION HISTORY
  !!     none
  !!
  !!*** **********************************************************************
  SUBROUTINE cp_sp_owned_b_iter_dealloc_ref(iterator, error)
    TYPE(cp_sp_owned_b_iter_dealloc_ref), INTENT(inout) :: iterator
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    CHARACTER(len=*), PARAMETER :: routineN='cp_sp_owned_b_iter_dealloc_ref',&
         routineP=moduleN//':'//routineN

    CALL cp_sp_local_b_iter_dealloc_ref(iterator%local_iterator,&
         error=error)
  END SUBROUTINE cp_sp_owned_b_iter_dealloc_ref
  !***************************************************************************

END MODULE cp_sp_iterators
