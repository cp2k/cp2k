!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2003 CP2K developers group
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* distribution_2d_types [1.0] *
!!
!!   NAME
!!     distribution_2d_types
!!
!!   FUNCTION
!!     stores a mapping of 2D info (e.g. matrix) on a 
!!     2D processor distribution (i.e. blacs grid)
!!     where cpus in the same blacs row own the same rows of the 2D info 
!!     (and similar for the cols)
!!
!!   NOTES
!!
!!   AUTHOR
!!     Joost VandeVondele (2003-08)
!!
MODULE distribution_2d_types
  
  USE kinds,                           ONLY: int_size,&
                                             dp,&
                                             dp_size

  USE cp_para_types,                   ONLY: cp_blacs_env_type
  USE cp_blacs_env,                    ONLY: cp_blacs_env_retain, &
                                             cp_blacs_env_release,&
                                             cp_blacs_env_write
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
                                             cp_logger_generate_filename,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_array_i_utils, ONLY: cp_1d_i_p_type, cp_1d_i_write
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE machine, only: m_flush
  USE message_passing, ONLY: mp_max, mp_sum
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='cp_distribution_2d_type'
  INTEGER, SAVE, private :: last_distribution_2d_id=0

  PUBLIC :: distribution_2d_type

  PUBLIC :: distribution_2d_create, &
            distribution_2d_release, &
            distribution_2d_retain,&
            distribution_2d_build_gdist,&
            distribution_2d_destroy_gdist,&
            distribution_2d_write,&
            distribution_2d_get

!!*** **********************************************************************

  !!****s* distribution_2d_types/distribution_2d_type [1.0] *
  !!
  !!   NAME
  !!     distribution_2d_type
  !!
  !!   FUNCTION
  !!     distributes pairs on a 2d grid of processors
  !!
  !!   NOTES
  !!     -
  !!
  !!   ATTRIBUTES
  !!     - row_distribution(i): processor row that owns the row i
  !!     - col_distribution(i): processor col that owns the col i
  !!     - n_row_distribution: nuber of global rows
  !!     - n_col_distribution: number of global cols
  !!     - n_local_rows(ikind): number of local rows of kind ikind
  !!     - n_local_cols(ikind): number of local cols of kind ikind
  !!     - local_cols(ikind)%array: ordered global indexes of the local cols
  !!       of kind ikind (might be oversized)
  !!     - local_rows(ikind)%array: ordered global indexes of the local 
  !!       rows of kind ikind (might be oversized)
  !!     - flat_local_rows: ordered global indexes of the local rows
  !!       (allocated on request, might be oversized)
  !!     - flat_local_cols: ordered global indexes of the local cols
  !!       (allocated on request, might be oversized)
  !!     - blacs_env: parallel environment in which the pairs are distributed
  !!     - ref_count: reference count (see doc/ReferenceCounting.html)
  !!     - id_nr: identification number (unique)
  !!
  !!   AUTHOR
  !!     Joost & Fawzi
  !!
  !!   MODIFICATION HISTORY
  !!     08.2003 created [joost]
  !!     09.2003 kind separation, minor cleanup [fawzi]
  !!
  !!   SOURCE
  !***************************************************************************
  TYPE distribution_2d_type
    INTEGER, DIMENSION(:), POINTER       :: row_distribution
    INTEGER, DIMENSION(:), POINTER       :: col_distribution
    INTEGER                              :: n_row_distribution
    INTEGER                              :: n_col_distribution
    INTEGER, DIMENSION(:), POINTER       :: n_local_rows
    INTEGER, DIMENSION(:), POINTER       :: n_local_cols
    TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER :: local_rows
    TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER :: local_cols
    INTEGER, DIMENSION(:), POINTER       :: flat_local_rows
    INTEGER, DIMENSION(:), POINTER       :: flat_local_cols
    TYPE(cp_blacs_env_type), POINTER     :: blacs_env
    INTEGER                              :: ref_count
    INTEGER                              :: id_nr
  END TYPE distribution_2d_type
  !!*** ***********************************************************************

CONTAINS

!!****f* distribution_2d_types/distribution_2d_create [1.0] *
!!
!!   NAME
!!     distribution_2d_create
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!     initializes the distribution_2d
!!
!!   NOTES
!!     the row and col_distribution are not allocated if not given
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!     09.2003 rewamped [fawzi]
!!
!!*** *********************************************************************
SUBROUTINE distribution_2d_create(distribution_2d, blacs_env,&
     local_rows_ptr, n_local_rows, &
     local_cols_ptr, row_distribution_ptr, col_distribution_ptr,&
     n_local_cols, n_row_distribution, n_col_distribution, error)
   TYPE(distribution_2d_type), POINTER, optional :: distribution_2d
   TYPE(cp_blacs_env_type), POINTER :: blacs_env
   TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER, OPTIONAL :: local_cols_ptr, &
        local_rows_ptr
   INTEGER, DIMENSION(:), POINTER, optional :: row_distribution_ptr,&
        col_distribution_ptr
   INTEGER, DIMENSION(:), INTENT(in), optional :: n_local_rows, n_local_cols
   INTEGER, INTENT(in), optional :: n_row_distribution, n_col_distribution
   TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

   CHARACTER(len=*), PARAMETER :: routineN='distribution_2d_create',&
                                  routineP=moduleN//':'//routineN
   LOGICAL :: failure
   INTEGER :: stat, i, local_row_size, local_col_size

   failure=.FALSE.
   CPPrecondition(ASSOCIATED(blacs_env),cp_failure_level,routineP,error,failure)
   CPPrecondition(.NOT.ASSOCIATED(distribution_2d),cp_failure_level,routineP,error,failure)

   IF (.NOT. failure) THEN
      ALLOCATE(distribution_2d,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
   END IF
   IF (.NOT. failure) THEN
      distribution_2d%ref_count = 1
      last_distribution_2d_id=last_distribution_2d_id+1
      distribution_2d%id_nr=last_distribution_2d_id

      NULLIFY(distribution_2d%col_distribution, distribution_2d%row_distribution,&
           distribution_2d%local_rows,distribution_2d%local_cols,&
           distribution_2d%blacs_env,distribution_2d%n_local_cols, &
           distribution_2d%n_local_rows, distribution_2d%flat_local_rows,&
           distribution_2d%flat_local_cols)

      distribution_2d%n_col_distribution=-huge(0)
      IF (PRESENT(col_distribution_ptr)) THEN
         distribution_2d%col_distribution => col_distribution_ptr
         distribution_2d%n_col_distribution =SIZE(distribution_2d%col_distribution)
      END IF
      IF(PRESENT(n_col_distribution)) THEN
         IF (ASSOCIATED(distribution_2d%col_distribution)) THEN
            CPPrecondition(n_col_distribution<=distribution_2d%n_col_distribution,cp_failure_level,routineP,error,failure)
            ! else alloc col_distribution?
         END IF
         distribution_2d%n_col_distribution=n_col_distribution
      END IF
      distribution_2d%n_row_distribution=-huge(0)
      IF (PRESENT(row_distribution_ptr)) THEN
         distribution_2d%row_distribution => row_distribution_ptr
         distribution_2d%n_row_distribution=size(distribution_2d%row_distribution)
      END IF
      IF(PRESENT(n_row_distribution)) THEN
         IF (ASSOCIATED(distribution_2d%row_distribution)) THEN
            CPPrecondition(n_row_distribution<=distribution_2d%n_row_distribution,cp_failure_level,routineP,error,failure)
            ! else alloc row_distribution?
         END IF
         distribution_2d%n_row_distribution=n_row_distribution
      END IF
         
      IF (PRESENT(local_rows_ptr)) &
           distribution_2d%local_rows => local_rows_ptr
      IF (.NOT.ASSOCIATED(distribution_2d%local_rows)) THEN
         CPPrecondition(PRESENT(n_local_rows),cp_failure_level,routineP,error,failure)
         ALLOCATE(distribution_2d%local_rows(SIZE(n_local_rows)),stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         DO i=1,SIZE(distribution_2d%local_rows)
            ALLOCATE(distribution_2d%local_rows(i)%array(n_local_rows(i)),&
                 stat=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            distribution_2d%local_rows(i)%array=-huge(0)
         END DO
      END IF
      ALLOCATE(distribution_2d%n_local_rows(SIZE(distribution_2d%local_rows)),&
           stat=stat)
      IF (PRESENT(n_local_rows)) THEN
         CPPrecondition(SIZE(distribution_2d%n_local_rows)==SIZE(n_local_rows),cp_failure_level,routineP,error,failure)
         DO i=1,SIZE(distribution_2d%n_local_rows)
            CPPrecondition(SIZE(distribution_2d%local_rows(i)%array)>=n_local_rows(i),cp_failure_level,routineP,error,failure)
            distribution_2d%n_local_rows(i) = n_local_rows(i)
         END DO
      ELSE
         DO i=1,SIZE(distribution_2d%n_local_rows)
            distribution_2d%n_local_rows(i) = &
                 SIZE(distribution_2d%local_rows(i)%array)
         END DO
      END IF

      IF (PRESENT(local_cols_ptr)) &
           distribution_2d%local_cols => local_cols_ptr
      IF (.NOT.ASSOCIATED(distribution_2d%local_cols)) THEN
         CPPrecondition(PRESENT(n_local_cols),cp_failure_level,routineP,error,failure)
         ALLOCATE(distribution_2d%local_cols(SIZE(n_local_cols)),stat=stat)
         CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         DO i=1,SIZE(distribution_2d%local_cols)
            ALLOCATE(distribution_2d%local_cols(i)%array(n_local_cols(i)),&
                 stat=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
            distribution_2d%local_cols(i)%array=-huge(0)
         END DO
      END IF
      ALLOCATE(distribution_2d%n_local_cols(SIZE(distribution_2d%local_cols)),&
           stat=stat)
      IF (PRESENT(n_local_cols)) THEN
         CPPrecondition(SIZE(distribution_2d%n_local_cols)==SIZE(n_local_cols),cp_failure_level,routineP,error,failure)
         DO i=1,SIZE(distribution_2d%n_local_cols)
            CPPrecondition(SIZE(distribution_2d%local_cols(i)%array)>=n_local_cols(i),cp_failure_level,routineP,error,failure)
            distribution_2d%n_local_cols(i) = n_local_cols(i)
         END DO
      ELSE
         DO i=1,SIZE(distribution_2d%n_local_cols)
            distribution_2d%n_local_cols(i) = &
                 SIZE(distribution_2d%local_cols(i)%array)
         END DO
      END IF

      distribution_2d%blacs_env => blacs_env
      CALL cp_blacs_env_retain(distribution_2d%blacs_env)
      
   END IF
 END SUBROUTINE distribution_2d_create
!***************************************************************************

!!****f* distribution_2d_types/distribution_2d_retain [1.0] *
!!
!!   NAME
!!     distribution_2d_retain
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
SUBROUTINE distribution_2d_retain(distribution_2d,error)
  TYPE(distribution_2d_type), POINTER :: distribution_2d
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  CHARACTER(len=*), PARAMETER :: routineN='distribution_2d_retain',&
       routineP=moduleN//':'//routineN
  LOGICAL :: failure
  
  failure=.false.
  CPPrecondition(ASSOCIATED(distribution_2d),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPreconditionNoFail(distribution_2d%ref_count>0,cp_failure_level,routineP,error)
     distribution_2d%ref_count=distribution_2d%ref_count+1
  ENDIF
END SUBROUTINE distribution_2d_retain

!!****f* distribution_2d_types/distribution_2d_release [1.0] *
!!
!!   NAME
!!     distribution_2d_release
!!
!!   SYNOPSIS
!!
!!   FUNCTION
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!
!!*** *********************************************************************
SUBROUTINE distribution_2d_release(distribution_2d,error)
   TYPE(distribution_2d_type), POINTER :: distribution_2d
   TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

   CHARACTER(len=*), PARAMETER :: routineN='distribution_2d_release',&
                                  routineP=moduleN//':'//routineN
   LOGICAL :: failure
   INTEGER :: istat,i

   failure=.FALSE.
   IF (ASSOCIATED(distribution_2d)) THEN
      CPPrecondition(distribution_2d%ref_count>0,cp_failure_level,routineP,error,failure)
      IF (.NOT. failure ) THEN
         distribution_2d%ref_count=distribution_2d%ref_count-1
         IF (distribution_2d%ref_count == 0 ) THEN
            CALL cp_blacs_env_release(distribution_2d%blacs_env)
            IF (ASSOCIATED(distribution_2d%col_distribution)) THEN
               DEALLOCATE(distribution_2d%col_distribution,stat=istat)
               CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
            END IF
            IF (ASSOCIATED(distribution_2d%row_distribution)) THEN
               DEALLOCATE(distribution_2d%row_distribution,stat=istat)
               CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
            END IF
            DO i=1,SIZE(distribution_2d%local_rows)
               DEALLOCATE(distribution_2d%local_rows(i)%array,stat=istat)
               CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
            END DO
            DEALLOCATE(distribution_2d%local_rows,stat=istat)
            CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
            DO i=1,SIZE(distribution_2d%local_cols)
               DEALLOCATE(distribution_2d%local_cols(i)%array,stat=istat)
               CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
            END DO
            DEALLOCATE(distribution_2d%local_cols,stat=istat)
            CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
            IF (ASSOCIATED(distribution_2d%flat_local_rows)) THEN
               DEALLOCATE(distribution_2d%flat_local_rows,stat=istat)
               CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
            END IF
            IF (ASSOCIATED(distribution_2d%flat_local_cols)) THEN
               DEALLOCATE(distribution_2d%flat_local_cols,stat=istat)
               CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
            END IF
            DEALLOCATE(distribution_2d,stat=istat)
            CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
         ENDIF
      ENDIF
   ENDIF
   NULLIFY(distribution_2d)
END SUBROUTINE distribution_2d_release

!!****f* cp_local_lists_types/distribution_2d_build_gdist [1.0] *
!!
!!   NAME
!!     distribution_2d_build_gdist
!!
!!   FUNCTION
!!     build the arrays with the prow/pcol of the processor that owns that index
!!     (-1 if no processor owns it)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - distribution_2d: the distribution for which you want to build the global
!!       distribution
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE distribution_2d_build_gdist(distribution_2d,error)
  TYPE(distribution_2d_type), POINTER :: distribution_2d
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='distribution_2d_build_gdist',&
       routineP=moduleN//':'//routineN
  INTEGER :: ilist, stat, max_el, max_loc, max_el1, i, my_id, nprow, npcol
  INTEGER, DIMENSION(:), POINTER :: g_dist

  failure=.FALSE.
  NULLIFY(g_dist)

  CPPrecondition(ASSOCIATED(distribution_2d),cp_failure_level,routineP,error,failure)
  CPPrecondition(distribution_2d%ref_count>0,cp_failure_level,routineP,error,failure)
  nprow=distribution_2d%blacs_env%num_pe(1)
  npcol=distribution_2d%blacs_env%num_pe(2)
  IF (.NOT. failure) THEN
     ! ** rows **
     max_el=0
     DO ilist=1,SIZE(distribution_2d%local_rows)
        max_el1=0
        !$omp parallel do private(i,max_loc) reduction(max:max_el1)
        DO i=1,distribution_2d%n_local_rows(ilist)
           max_el1=MAX(max_el1,distribution_2d%local_rows(ilist)%array(i))
        END DO
        IF (max_el < max_el1) max_el=max_el1
     END DO
     CALL mp_max(max_el, distribution_2d%blacs_env%para_env%group)
     if (distribution_2d%n_row_distribution<0) &
          distribution_2d%n_row_distribution=max_el
     CPPrecondition(distribution_2d%n_row_distribution>=max_el,cp_failure_level,routineP,error,failure)
     IF (ASSOCIATED(distribution_2d%row_distribution)) THEN
        IF (SIZE(distribution_2d%row_distribution)<&
             distribution_2d%n_row_distribution) THEN
           DEALLOCATE(distribution_2d%row_distribution,&
                stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
     END IF
     IF (.NOT.ASSOCIATED(distribution_2d%row_distribution)) THEN
        ALLOCATE(distribution_2d%row_distribution &
             (distribution_2d%n_row_distribution),&
             stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF

     g_dist => distribution_2d%row_distribution
     g_dist=0
     my_id=distribution_2d%blacs_env%para_env%mepos+1

     DO ilist=1,SIZE(distribution_2d%local_rows)
        !$omp parallel do
        DO i=1,distribution_2d%n_local_rows(ilist)
           ! use += to catch multiple int on same proc (= would be faster)
           g_dist(distribution_2d%local_rows(ilist)%array(i))=&
                g_dist(distribution_2d%local_rows(ilist)%array(i))+my_id
        END DO
     END DO
#ifdef __SCALAPACK
     CALL IGSUM2D(distribution_2d%blacs_env%group,"C"," ",&
          distribution_2d%n_row_distribution,1,g_dist(1),&
          distribution_2d%n_row_distribution,-1,-1)
     ! check if different cols have differend distributions?
     !$omp parallel do
     DO i=1,SIZE(g_dist)
        g_dist(i)=g_dist(i)-1
     END DO
#else
     CALL mp_sum(g_dist,distribution_2d%blacs_env%para_env%group)
     !$omp parallel do
     DO i=1,SIZE(g_dist)
        g_dist(i)=g_dist(i)/npcol-1
     END DO
#endif

     my_id=my_id-1
     DO ilist=1,SIZE(distribution_2d%local_rows)
        !$omp parallel do
        DO i=1,distribution_2d%n_local_rows(ilist)
           IF (g_dist(distribution_2d%local_rows(ilist)%array(i))/=my_id) THEN
              CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                   routineP,"invalid distribution, two processors have "//&
                   "the same number (list="//cp_to_string(ilist)//", i="//&
                   cp_to_string(i)//",g_index="//&
                   cp_to_string(distribution_2d%local_rows(ilist)%array(i))//")"//&
                   CPSourceFileRef,&
                   error,failure)
           END IF
        END DO
     END DO

     ! ** cols **
     max_el=0
     DO ilist=1,SIZE(distribution_2d%local_cols)
        max_el1=0
        !$omp parallel do private(i,max_loc) reduction(max:max_el1)
        DO i=1,distribution_2d%n_local_cols(ilist)
           max_el1=MAX(max_el1,distribution_2d%local_cols(ilist)%array(i))
        END DO
        IF (max_el < max_el1) max_el=max_el1
     END DO
     CALL mp_max(max_el, distribution_2d%blacs_env%para_env%group)
     if (distribution_2d%n_col_distribution<0) &
          distribution_2d%n_col_distribution=max_el
     CPPrecondition(distribution_2d%n_col_distribution>=max_el,cp_failure_level,routineP,error,failure)
     IF (ASSOCIATED(distribution_2d%col_distribution)) THEN
        IF (SIZE(distribution_2d%col_distribution)<&
             distribution_2d%n_col_distribution) THEN
           DEALLOCATE(distribution_2d%col_distribution,stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
     END IF
     IF (.NOT.ASSOCIATED(distribution_2d%col_distribution)) THEN
        ALLOCATE(distribution_2d%col_distribution&
             (distribution_2d%n_col_distribution),&
             stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
     END IF

     g_dist => distribution_2d%col_distribution
     g_dist=0
     my_id=distribution_2d%blacs_env%para_env%mepos+1

     DO ilist=1,SIZE(distribution_2d%local_cols)
        !$omp parallel do
        DO i=1,distribution_2d%n_local_cols(ilist)
           ! use += to catch multiple int on same proc (= would be faster)
           g_dist(distribution_2d%local_cols(ilist)%array(i))=&
                g_dist(distribution_2d%local_cols(ilist)%array(i))+my_id
        END DO
     END DO
#ifdef __SCALAPACK
     CALL IGSUM2D(distribution_2d%blacs_env%group,"R"," ",&
          distribution_2d%n_row_distribution,1,g_dist(1),&
          distribution_2d%n_row_distribution,-1,-1)
     ! check if different rows have different distributions?
     !$omp parallel do
     DO i=1,SIZE(g_dist)
        g_dist(i)=g_dist(i)-1
     END DO
#else
     CALL mp_sum(g_dist,distribution_2d%blacs_env%para_env%group)
     !$omp parallel do
     DO i=1,SIZE(g_dist)
        g_dist(i)=g_dist(i)/nprow-1
     END DO
#endif

     my_id=my_id-1
     DO ilist=1,SIZE(distribution_2d%local_cols)
        !$omp parallel do
        DO i=1,distribution_2d%n_local_cols(ilist)
           IF (g_dist(distribution_2d%local_cols(ilist)%array(i))/=my_id) THEN
              CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                   routineP,"invalid distribution, two processors have "//&
                   "the same number (list="//cp_to_string(ilist)//", i="//&
                   cp_to_string(i)//",g_index="//&
                   cp_to_string(distribution_2d%local_cols(ilist)%array(i))//")"//&
                   CPSourceFileRef,&
                   error,failure)
           END IF
        END DO
     END DO

  END IF
END SUBROUTINE distribution_2d_build_gdist
!***************************************************************************

!!****f* distribution_2d_types/distribution_2d_destroy_gdist [1.0] *
!!
!!   NAME
!!     distribution_2d_destroy_gdist
!!
!!   FUNCTION
!!     deletes the global distribution
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - distribution_2d: the distribution from which you want to remove
!!       the global distribution
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE distribution_2d_destroy_gdist(distribution_2d,error)
  TYPE(distribution_2d_type), pointer :: distribution_2d
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='distribution_2d_destroy_gdist',&
        routineP=moduleN//':'//routineN
  integer :: stat

  failure=.false.

  CPPrecondition(ASSOCIATED(distribution_2d),cp_failure_level,routineP,error,failure)
  CPPrecondition(distribution_2d%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (ASSOCIATED(distribution_2d%row_distribution)) THEN
        DEALLOCATE(distribution_2d%row_distribution, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
     IF (ASSOCIATED(distribution_2d%col_distribution)) THEN
        DEALLOCATE(distribution_2d%col_distribution, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
END SUBROUTINE distribution_2d_destroy_gdist
!***************************************************************************

!!****f* distribution_2d_types/distribution_2d_write [1.0] *
!!
!!   NAME
!!     distribution_2d_write
!!
!!   FUNCTION
!!     writes out the given distribution
!!
!!   NOTES
!!     to clean up, make safer wrt. grabage in distribution_2d%n_*
!!
!!   INPUTS
!!     - distribution_2d: the distribution to write out
!!     - unit_nr: the unit to write to
!!     - local: if the unit is local to to each processor (otherwise
!!       only the processor with logger%para_env%source==
!!       logger%para_env%mepos writes), defaults to false.
!!     - long_description: if a long description should be given,
!!       defaults to false
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2003 adapted qs_distribution_2d_create write done by Matthias[fawzi]
!!
!!*** **********************************************************************
SUBROUTINE distribution_2d_write(distribution_2d, unit_nr, local,&
     long_description, error)
  TYPE(distribution_2d_type), pointer :: distribution_2d
  INTEGER, INTENT(in) :: unit_nr
  LOGICAL, INTENT(in), OPTIONAL :: local, long_description
  type(cp_error_type), optional, intent(inout) :: error
  
  LOGICAL :: failure, my_local, my_long_description
  character(len=*), parameter :: routineN='distribution_2d',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat, i
  INTEGER, DIMENSION(:), ALLOCATABLE :: nparticle
  TYPE(cp_logger_type), POINTER :: logger

  failure=.FALSE.
  logger => cp_error_get_logger(error)
  
  my_long_description=.FALSE.
  IF (PRESENT(long_description)) my_long_description=long_description
  my_local=.FALSE.
  IF (PRESENT(local)) my_local=local
  IF (.NOT.my_local) my_local=(logger%para_env%source==logger%para_env%mepos)
  
  IF (ASSOCIATED(distribution_2d)) THEN
     IF (my_local) THEN
        WRITE (unit=unit_nr,&
             fmt="(' <distribution_2d>{ id_nr=',i10,'ref_count=',i10,',')")&
             distribution_2d%id_nr, distribution_2d%ref_count

        WRITE (unit=unit_nr,fmt="(' n_row_distribution=',i15,',')")distribution_2d%n_row_distribution
        IF (ASSOCIATED(distribution_2d%row_distribution)) THEN
           IF (my_long_description) THEN
              WRITE (unit=unit_nr,fmt="(' row_distribution= (')", advance="no")
              DO i=1,SIZE(distribution_2d%row_distribution)
                 WRITE(unit=unit_nr,fmt="(i6,',')", advance="no")distribution_2d%row_distribution(i)
              END DO
              WRITE (unit=unit_nr,fmt="('),')")
           ELSE
              WRITE (unit=unit_nr,fmt="(' row_distribution= array(',i6,':',i6,'),')")&
                   LBOUND(distribution_2d%row_distribution),&
                   UBOUND(distribution_2d%row_distribution)
           END IF
        ELSE
           WRITE (unit=unit_nr,fmt="(' row_distribution=*null*,')")
        END IF

        WRITE (unit=unit_nr,fmt="(' n_col_distribution=',i15,',')")distribution_2d%n_col_distribution
        IF (ASSOCIATED(distribution_2d%col_distribution)) THEN
           IF (my_long_description) THEN
              WRITE (unit=unit_nr,fmt="(' col_distribution= (')", advance="no")
              DO i=1,SIZE(distribution_2d%col_distribution)
                 WRITE(unit=unit_nr,fmt="(i6,',')", advance="no")distribution_2d%col_distribution(i)
              END DO
              WRITE (unit=unit_nr,fmt="('),')")
           ELSE
              WRITE (unit=unit_nr,fmt="(' col_distribution= array(',i6,':',i6,'),')")&
                   LBOUND(distribution_2d%col_distribution),&
                   UBOUND(distribution_2d%col_distribution)
           END IF
        ELSE
           WRITE (unit=unit_nr,fmt="(' col_distribution=*null*,')")
        END IF

        IF (ASSOCIATED(distribution_2d%n_local_rows)) THEN
           IF (my_long_description) THEN
              WRITE (unit=unit_nr,fmt="(' n_local_rows= (')", advance="no")
              DO i=1,SIZE(distribution_2d%n_local_rows)
                 WRITE(unit=unit_nr,fmt="(i6,',')", advance="no")distribution_2d%n_local_rows(i)
              END DO
              WRITE (unit=unit_nr,fmt="('),')")
           ELSE
              WRITE (unit=unit_nr,fmt="(' n_local_rows= array(',i6,':',i6,'),')")&
                   LBOUND(distribution_2d%n_local_rows),&
                   UBOUND(distribution_2d%n_local_rows)
           END IF
        ELSE
           WRITE (unit=unit_nr,fmt="(' n_local_rows=*null*,')")
        END IF

        IF (ASSOCIATED(distribution_2d%local_rows)) THEN
           WRITE (unit=unit_nr,fmt="(' local_rows=(')")
           DO i=1,SIZE(distribution_2d%local_rows)
              IF (ASSOCIATED(distribution_2d%local_rows(i)%array)) THEN
                 IF (my_long_description) THEN
                    CALL cp_1d_i_write(array=distribution_2d%local_rows(i)%array, &
                         unit_nr=unit_nr, error=error)
                 ELSE
                    WRITE (unit=unit_nr,fmt="(' array(',i6,':',i6,'),')")&
                         LBOUND(distribution_2d%local_rows(i)%array),&
                         UBOUND(distribution_2d%local_rows(i)%array)
                 END IF
              ELSE
                 WRITE (unit=unit_nr,fmt="('*null*')")
              END IF
           END DO
           WRITE (unit=unit_nr,fmt="(' ),')")
        ELSE
           WRITE (unit=unit_nr,fmt="(' local_rows=*null*,')")           
        END IF

        IF (ASSOCIATED(distribution_2d%n_local_cols)) THEN
           IF (my_long_description) THEN
              WRITE (unit=unit_nr,fmt="(' n_local_cols= (')", advance="no")
              DO i=1,SIZE(distribution_2d%n_local_cols)
                 WRITE(unit=unit_nr,fmt="(i6,',')", advance="no")distribution_2d%n_local_cols(i)
              END DO
              WRITE (unit=unit_nr,fmt="('),')")
           ELSE
              WRITE (unit=unit_nr,fmt="(' n_local_cols= array(',i6,':',i6,'),')")&
                   LBOUND(distribution_2d%n_local_cols),&
                   UBOUND(distribution_2d%n_local_cols)
           END IF
        ELSE
           WRITE (unit=unit_nr,fmt="(' n_local_cols=*null*,')")
        END IF

        IF (ASSOCIATED(distribution_2d%local_cols)) THEN
           WRITE (unit=unit_nr,fmt="(' local_cols=(')")
           DO i=1,SIZE(distribution_2d%local_cols)
              IF (ASSOCIATED(distribution_2d%local_cols(i)%array)) THEN
                 IF (my_long_description) THEN
                    CALL cp_1d_i_write(array=distribution_2d%local_cols(i)%array, &
                         unit_nr=unit_nr, error=error)
                 ELSE
                    WRITE (unit=unit_nr,fmt="(' array(',i6,':',i6,'),')")&
                         LBOUND(distribution_2d%local_cols(i)%array),&
                         UBOUND(distribution_2d%local_cols(i)%array)
                 END IF
              ELSE
                 WRITE (unit=unit_nr,fmt="('*null*')")
              END IF
           END DO
           WRITE (unit=unit_nr,fmt="(' ),')")
        ELSE
           WRITE (unit=unit_nr,fmt="(' local_cols=*null*,')")           
        END IF

        IF (ASSOCIATED(distribution_2d%blacs_env)) THEN
           IF (my_long_description) THEN
              WRITE (unit=unit_nr,fmt="(' blacs_env=')",advance="no")
              CALL cp_blacs_env_write(distribution_2d%blacs_env,unit_nr=unit_nr,&
                   error=error)
           ELSE
              WRITE (unit=unit_nr,fmt="(' blacs_env=<blacs_env id=',i6,'>')")&
                   distribution_2d%blacs_env%group
           END IF
        ELSE
           WRITE (unit=unit_nr,fmt="(' blacs_env=*null*')")
        END IF

        WRITE (unit=unit_nr,fmt="(' }')")
     END IF
     
  ELSE IF (my_local) THEN
     WRITE (unit=unit_nr,&
          fmt="(' <distribution_2d *null*>')")
  END IF
  
  CALL m_flush(unit_nr)
  
END SUBROUTINE distribution_2d_write
!***************************************************************************

!!****f* distribution_2d_types/distribution_2d_get [1.0] *
!!
!!   NAME
!!     distribution_2d_get
!!
!!   FUNCTION
!!     returns various attributes about the distribution_2d
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - distribution_2d: the object you want info about
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!     see distribution_2d_type attributes for a description of the other
!!     arguments
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     09.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE distribution_2d_get(distribution_2d, row_distribution, &
     col_distribution, n_row_distribution, n_col_distribution,&
     n_local_rows, n_local_cols, local_rows, local_cols,&
     flat_local_rows, flat_local_cols, n_flat_local_rows, n_flat_local_cols,&
     blacs_env, id_nr,&
     error)
  TYPE(distribution_2d_type), POINTER :: distribution_2d
  INTEGER, DIMENSION(:), POINTER, OPTIONAL :: row_distribution, col_distribution
  INTEGER, INTENT(out), OPTIONAL :: n_row_distribution,n_col_distribution
  INTEGER, DIMENSION(:), POINTER, OPTIONAL :: n_local_rows, n_local_cols
  TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER, OPTIONAL :: local_rows
  TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER, OPTIONAL :: local_cols
  INTEGER, DIMENSION(:), POINTER, OPTIONAL :: flat_local_rows, flat_local_cols
  INTEGER, INTENT(out), OPTIONAL :: n_flat_local_rows, n_flat_local_cols
  TYPE(cp_blacs_env_type), POINTER, OPTIONAL :: blacs_env
  INTEGER, INTENT(out), optional :: id_nr
  type(cp_error_type), optional, intent(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='distribution_2d_get',&
       routineP=moduleN//':'//routineN
  INTEGER :: iblock_atomic, iblock_min, ikind_min,ikind, stat
  INTEGER, DIMENSION(:), allocatable :: multiindex
  
  failure=.FALSE.
  CPPrecondition(ASSOCIATED(distribution_2d),cp_failure_level,routineP,error,failure)
  CPPrecondition(distribution_2d%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(row_distribution)) row_distribution => distribution_2d%row_distribution
     IF (PRESENT(col_distribution)) col_distribution => distribution_2d%col_distribution
     IF (PRESENT(n_row_distribution)) n_row_distribution=distribution_2d%n_row_distribution
     IF (PRESENT(n_col_distribution)) n_col_distribution=distribution_2d%n_col_distribution
     IF (PRESENT(n_local_rows)) n_local_rows => distribution_2d%n_local_rows
     IF (PRESENT(n_local_cols)) n_local_cols => distribution_2d%n_local_cols
     IF (PRESENT(local_rows)) local_rows => distribution_2d%local_rows
     IF (PRESENT(local_cols)) local_cols => distribution_2d%local_cols
     IF (PRESENT(flat_local_rows)) THEN
        IF (.NOT.ASSOCIATED(distribution_2d%flat_local_rows)) THEN
           ALLOCATE(multiindex(SIZE(distribution_2d%local_rows)),&
                distribution_2d%flat_local_rows(SUM(distribution_2d%n_local_rows)),&
                stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           multiindex=1
           DO iblock_atomic=1,SIZE(distribution_2d%flat_local_rows)
              iblock_min=HUGE(0)
              ikind_min=-HUGE(0)
              DO ikind=1,SIZE(distribution_2d%local_rows)
                 IF (multiindex(ikind)<=distribution_2d%n_local_rows(ikind)) THEN
                    IF (distribution_2d%local_rows(ikind)%array(multiindex(ikind))<&
                         iblock_min) THEN
                       iblock_min=distribution_2d%local_rows(ikind)%array(multiindex(ikind))
                       ikind_min=ikind
                    END IF
                 END IF
              END DO
              CPPostcondition(ikind_min>0,cp_failure_level,routineP,error,failure)
              distribution_2d%flat_local_rows(iblock_atomic)=&
                   distribution_2d%local_rows(ikind_min)%array(multiindex(ikind_min))
              multiindex(ikind_min)=multiindex(ikind_min)+1
           END DO
           DEALLOCATE(multiindex, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        flat_local_rows => distribution_2d%flat_local_rows
     END IF
     IF (PRESENT(flat_local_cols)) THEN
        IF (.NOT.ASSOCIATED(distribution_2d%flat_local_cols)) THEN
           ALLOCATE(multiindex(SIZE(distribution_2d%local_cols)),&
                distribution_2d%flat_local_cols(SUM(distribution_2d%n_local_cols)),&
                stat=stat)
           CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
           multiindex=1
           DO iblock_atomic=1,SIZE(distribution_2d%flat_local_cols)
              iblock_min=HUGE(0)
              ikind_min=-HUGE(0)
              DO ikind=1,SIZE(distribution_2d%local_cols)
                 IF (multiindex(ikind)<=distribution_2d%n_local_cols(ikind)) THEN
                    IF (distribution_2d%local_cols(ikind)%array(multiindex(ikind))<&
                         iblock_min) THEN
                       iblock_min=distribution_2d%local_cols(ikind)%array(multiindex(ikind))
                       ikind_min=ikind
                    END IF
                 END IF
              END DO
              CPPostcondition(ikind_min>0,cp_failure_level,routineP,error,failure)
              distribution_2d%flat_local_cols(iblock_atomic)=&
                   distribution_2d%local_cols(ikind_min)%array(multiindex(ikind_min))
              multiindex(ikind_min)=multiindex(ikind_min)+1
           END DO
           DEALLOCATE(multiindex, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
        flat_local_cols => distribution_2d%flat_local_cols
     END IF
     IF (PRESENT(n_flat_local_rows)) n_flat_local_rows=SUM(distribution_2d%n_local_rows)
     IF (PRESENT(n_flat_local_cols)) n_flat_local_cols=SUM(distribution_2d%n_local_cols)
     IF (PRESENT(blacs_env)) blacs_env => distribution_2d%blacs_env
     IF (PRESENT(id_nr)) id_nr = distribution_2d%id_nr
  END IF
END SUBROUTINE distribution_2d_get
!***************************************************************************

END MODULE distribution_2d_types
