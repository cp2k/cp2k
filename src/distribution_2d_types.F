!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief stores a mapping of 2D info (e.g. matrix) on a
!>      2D processor distribution (i.e. blacs grid)
!>      where cpus in the same blacs row own the same rows of the 2D info
!>      (and similar for the cols)
!> \author Joost VandeVondele (2003-08)
! *****************************************************************************
MODULE distribution_2d_types

  USE cp_array_i_utils,                ONLY: cp_1d_i_p_type,&
                                             cp_1d_i_write
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_retain,&
                                             cp_blacs_env_type,&
                                             cp_blacs_env_write
  USE machine,                         ONLY: m_flush
#include "common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'distribution_2d_types'
  INTEGER, SAVE, PRIVATE :: last_distribution_2d_id=0

  PUBLIC :: distribution_2d_type

  PUBLIC :: distribution_2d_create,&
            distribution_2d_release,&
            distribution_2d_retain,&
            distribution_2d_write,&
            distribution_2d_get

! *****************************************************************************
!> \brief distributes pairs on a 2d grid of processors
!> \param row_distribution (i): processor row that owns the row i
!> \param col_distribution (i): processor col that owns the col i
!> \param n_row_distribution nuber of global rows
!> \param n_col_distribution number of global cols
!> \param n_local_rows (ikind): number of local rows of kind ikind
!> \param n_local_cols (ikind): number of local cols of kind ikind
!> \param local_cols (ikind)%array: ordered global indexes of the local cols
!>        of kind ikind (might be oversized)
!> \param local_rows (ikind)%array: ordered global indexes of the local
!>        rows of kind ikind (might be oversized)
!> \param flat_local_rows ordered global indexes of the local rows
!>        (allocated on request, might be oversized)
!> \param flat_local_cols ordered global indexes of the local cols
!>        (allocated on request, might be oversized)
!> \param blacs_env parallel environment in which the pairs are distributed
!> \param ref_count reference count (see doc/ReferenceCounting.html)
!> \param id_nr identification number (unique)
!> \par History
!>      08.2003 created [joost]
!>      09.2003 kind separation, minor cleanup [fawzi]
!> \author Joost & Fawzi
! *****************************************************************************
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

CONTAINS

! *****************************************************************************
!> \brief initializes the distribution_2d
!> \param distribution_2d ...
!> \param blacs_env ...
!> \param local_rows_ptr ...
!> \param n_local_rows ...
!> \param local_cols_ptr ...
!> \param row_distribution_ptr ...
!> \param col_distribution_ptr ...
!> \param n_local_cols ...
!> \param n_row_distribution ...
!> \param n_col_distribution ...
!> \param error ...
!> \par History
!>      09.2003 rewamped [fawzi]
!> \author Joost VandeVondele
!> \note
!>      the row and col_distribution are not allocated if not given
! *****************************************************************************
SUBROUTINE distribution_2d_create(distribution_2d, blacs_env,&
     local_rows_ptr, n_local_rows, &
     local_cols_ptr, row_distribution_ptr, col_distribution_ptr,&
     n_local_cols, n_row_distribution, n_col_distribution, error)
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_1d_i_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: local_rows_ptr
    INTEGER, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: n_local_rows
    TYPE(cp_1d_i_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: local_cols_ptr
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: row_distribution_ptr, &
                                                col_distribution_ptr
    INTEGER, DIMENSION(:), INTENT(in), &
      OPTIONAL                               :: n_local_cols
    INTEGER, INTENT(in), OPTIONAL            :: n_row_distribution, &
                                                n_col_distribution
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distribution_2d_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

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

      distribution_2d%n_col_distribution=-HUGE(0)
      IF (PRESENT(col_distribution_ptr)) THEN
         distribution_2d%col_distribution => col_distribution_ptr
         distribution_2d%n_col_distribution =SIZE(distribution_2d%col_distribution)
      END IF
      IF(PRESENT(n_col_distribution)) THEN
         IF (ASSOCIATED(distribution_2d%col_distribution)) THEN
            CALL cp_assert(n_col_distribution<=distribution_2d%n_col_distribution,&
                 cp_assertion_failed,cp_failure_level,routineP,&
                 "n_col_distribution<=distribution_2d%n_col_distribution",error,failure)
            ! else alloc col_distribution?
         END IF
         distribution_2d%n_col_distribution=n_col_distribution
      END IF
      distribution_2d%n_row_distribution=-HUGE(0)
      IF (PRESENT(row_distribution_ptr)) THEN
         distribution_2d%row_distribution => row_distribution_ptr
         distribution_2d%n_row_distribution=SIZE(distribution_2d%row_distribution)
      END IF
      IF(PRESENT(n_row_distribution)) THEN
         IF (ASSOCIATED(distribution_2d%row_distribution)) THEN
            CALL cp_assert(n_row_distribution<=distribution_2d%n_row_distribution,&
                 cp_failure_level,cp_assertion_failed,routineP,&
                 "n_row_distribution<=distribution_2d%n_row_distribution",error,failure)
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
            distribution_2d%local_rows(i)%array=-HUGE(0)
         END DO
      END IF
      ALLOCATE(distribution_2d%n_local_rows(SIZE(distribution_2d%local_rows)),&
           stat=stat)
      IF (PRESENT(n_local_rows)) THEN
         CALL cp_assert(SIZE(distribution_2d%n_local_rows)==SIZE(n_local_rows),&
              cp_failure_level,cp_assertion_failed,routineP,&
              "SIZE(distribution_2d%n_local_rows)==SIZE(n_local_rows)",error,failure)
         DO i=1,SIZE(distribution_2d%n_local_rows)
            CALL cp_assert(SIZE(distribution_2d%local_rows(i)%array)>=n_local_rows(i),&
                 cp_failure_level,cp_assertion_failed,routineP,&
                 "SIZE(distribution_2d%local_rows(i)%array)>=n_local_rows(i)",error,failure)
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
            distribution_2d%local_cols(i)%array=-HUGE(0)
         END DO
      END IF
      ALLOCATE(distribution_2d%n_local_cols(SIZE(distribution_2d%local_cols)),&
           stat=stat)
      IF (PRESENT(n_local_cols)) THEN
         CALL cp_assert(SIZE(distribution_2d%n_local_cols)==SIZE(n_local_cols),&
              cp_failure_level,cp_assertion_failed,routineP,&
              "SIZE(distribution_2d%n_local_cols)==SIZE(n_local_cols)",error,failure)
         DO i=1,SIZE(distribution_2d%n_local_cols)
            CALL cp_assert(SIZE(distribution_2d%local_cols(i)%array)>=n_local_cols(i),&
                 cp_failure_level,cp_assertion_failed,routineP,&
                 "SIZE(distribution_2d%local_cols(i)%array)>=n_local_cols(i)",error,failure)
            distribution_2d%n_local_cols(i) = n_local_cols(i)
         END DO
      ELSE
         DO i=1,SIZE(distribution_2d%n_local_cols)
            distribution_2d%n_local_cols(i) = &
                 SIZE(distribution_2d%local_cols(i)%array)
         END DO
      END IF

      distribution_2d%blacs_env => blacs_env
      CALL cp_blacs_env_retain(distribution_2d%blacs_env,error=error)

   END IF
 END SUBROUTINE distribution_2d_create

! *****************************************************************************
!> \brief ...
!> \param distribution_2d ...
!> \param error ...
!> \author Joost VandeVondele
! *****************************************************************************
SUBROUTINE distribution_2d_retain(distribution_2d,error)
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distribution_2d_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  CPPrecondition(ASSOCIATED(distribution_2d),cp_failure_level,routineP,error,failure)
  IF (.NOT.failure) THEN
     CPPreconditionNoFail(distribution_2d%ref_count>0,cp_failure_level,routineP,error)
     distribution_2d%ref_count=distribution_2d%ref_count+1
  ENDIF
END SUBROUTINE distribution_2d_retain

! *****************************************************************************
!> \brief ...
!> \param distribution_2d ...
!> \param error ...
! *****************************************************************************
SUBROUTINE distribution_2d_release(distribution_2d,error)
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distribution_2d_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istat
    LOGICAL                                  :: failure

   failure=.FALSE.
   IF (ASSOCIATED(distribution_2d)) THEN
      CPPrecondition(distribution_2d%ref_count>0,cp_failure_level,routineP,error,failure)
      IF (.NOT. failure ) THEN
         distribution_2d%ref_count=distribution_2d%ref_count-1
         IF (distribution_2d%ref_count == 0 ) THEN
            CALL cp_blacs_env_release(distribution_2d%blacs_env,error=error)
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
            IF (ASSOCIATED(distribution_2d%n_local_rows)) THEN
               DEALLOCATE(distribution_2d%n_local_rows,stat=istat)
               CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
            END IF
            IF (ASSOCIATED(distribution_2d%n_local_cols)) THEN
               DEALLOCATE(distribution_2d%n_local_cols,stat=istat)
               CPPostconditionNoFail(istat==0,cp_warning_level,routineP,error)
            END IF
            DEALLOCATE(distribution_2d,stat=istat)
            CPPostcondition(istat==0,cp_warning_level,routineP,error,failure)
         ENDIF
      ENDIF
   ENDIF
   NULLIFY(distribution_2d)
END SUBROUTINE distribution_2d_release

! *****************************************************************************
!> \brief writes out the given distribution
!> \param distribution_2d the distribution to write out
!> \param unit_nr the unit to write to
!> \param local if the unit is local to to each processor (otherwise
!>        only the processor with logger%para_env%source==
!>        logger%para_env%mepos writes), defaults to false.
!> \param long_description if a long description should be given,
!>        defaults to false
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      08.2003 adapted qs_distribution_2d_create write done by Matthias[fawzi]
!> \author Fawzi Mohamed
!> \note
!>      to clean up, make safer wrt. grabage in distribution_2d%n_*
! *****************************************************************************
SUBROUTINE distribution_2d_write(distribution_2d, unit_nr, local,&
     long_description, error)
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    INTEGER, INTENT(in)                      :: unit_nr
    LOGICAL, INTENT(in), OPTIONAL            :: local, long_description
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distribution_2d_write', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure, my_local, &
                                                my_long_description
    TYPE(cp_logger_type), POINTER            :: logger

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
             fmt="(/,' <distribution_2d> {      id_nr=',i10,'      ref_count=',i10,',')")&
             distribution_2d%id_nr, distribution_2d%ref_count

        WRITE (unit=unit_nr,fmt="('    n_row_distribution=',i15,',')")&
             distribution_2d%n_row_distribution
        IF (ASSOCIATED(distribution_2d%row_distribution)) THEN
           IF (my_long_description) THEN
              WRITE (unit=unit_nr,fmt="('      row_distribution= (')", advance="no")
              DO i=1,SIZE(distribution_2d%row_distribution)
                 WRITE(unit=unit_nr,fmt="(i6,',')", advance="no")distribution_2d%row_distribution(i)
                 ! keep lines finite, so that we can open outputs in vi
                 IF (MODULO(i,8)==0 .AND. i.NE.SIZE(distribution_2d%row_distribution)) &
                    WRITE(unit=unit_nr,fmt='()')
              END DO
              WRITE (unit=unit_nr,fmt="('),')")
           ELSE
              WRITE (unit=unit_nr,fmt="('      row_distribution= array(',i6,':',i6,'),')")&
                   LBOUND(distribution_2d%row_distribution),&
                   UBOUND(distribution_2d%row_distribution)
           END IF
        ELSE
           WRITE (unit=unit_nr,fmt="('       row_distribution=*null*,')")
        END IF

        WRITE (unit=unit_nr,fmt="('    n_col_distribution=',i15,',')")&
             distribution_2d%n_col_distribution
        IF (ASSOCIATED(distribution_2d%col_distribution)) THEN
           IF (my_long_description) THEN
              WRITE (unit=unit_nr,fmt="('      col_distribution= (')", advance="no")
              DO i=1,SIZE(distribution_2d%col_distribution)
                 WRITE(unit=unit_nr,fmt="(i6,',')", advance="no")distribution_2d%col_distribution(i)
                 ! keep lines finite, so that we can open outputs in vi
                 IF (MODULO(i,8)==0 .AND. i.NE.SIZE(distribution_2d%col_distribution)) &
                    WRITE(unit=unit_nr,fmt='()')
              END DO
              WRITE (unit=unit_nr,fmt="('),')")
           ELSE
              WRITE (unit=unit_nr,fmt="('      col_distribution= array(',i6,':',i6,'),')")&
                   LBOUND(distribution_2d%col_distribution),&
                   UBOUND(distribution_2d%col_distribution)
           END IF
        ELSE
           WRITE (unit=unit_nr,fmt="('       col_distribution=*null*,')")
        END IF

        IF (ASSOCIATED(distribution_2d%n_local_rows)) THEN
           IF (my_long_description) THEN
              WRITE (unit=unit_nr,fmt="('    n_local_rows= (')", advance="no")
              DO i=1,SIZE(distribution_2d%n_local_rows)
                 WRITE(unit=unit_nr,fmt="(i6,',')", advance="no")distribution_2d%n_local_rows(i)
                 ! keep lines finite, so that we can open outputs in vi
                 IF (MODULO(i,10)==0 .AND. i.NE.SIZE(distribution_2d%n_local_rows)) &
                    WRITE(unit=unit_nr,fmt='()')
              END DO
              WRITE (unit=unit_nr,fmt="('),')")
           ELSE
              WRITE (unit=unit_nr,fmt="('    n_local_rows= array(',i6,':',i6,'),')")&
                   LBOUND(distribution_2d%n_local_rows),&
                   UBOUND(distribution_2d%n_local_rows)
           END IF
        ELSE
           WRITE (unit=unit_nr,fmt="('    n_local_rows=*null*,')")
        END IF

        IF (ASSOCIATED(distribution_2d%local_rows)) THEN
           WRITE (unit=unit_nr,fmt="('      local_rows=(')")
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
           WRITE (unit=unit_nr,fmt="('      local_rows=*null*,')")
        END IF

        IF (ASSOCIATED(distribution_2d%n_local_cols)) THEN
           IF (my_long_description) THEN
              WRITE (unit=unit_nr,fmt="('    n_local_cols= (')", advance="no")
              DO i=1,SIZE(distribution_2d%n_local_cols)
                 WRITE(unit=unit_nr,fmt="(i6,',')", advance="no")distribution_2d%n_local_cols(i)
                 ! keep lines finite, so that we can open outputs in vi
                 IF (MODULO(i,10)==0 .AND. i.NE.SIZE(distribution_2d%n_local_cols)) &
                    WRITE(unit=unit_nr,fmt='()')
              END DO
              WRITE (unit=unit_nr,fmt="('),')")
           ELSE
              WRITE (unit=unit_nr,fmt="('    n_local_cols= array(',i6,':',i6,'),')")&
                   LBOUND(distribution_2d%n_local_cols),&
                   UBOUND(distribution_2d%n_local_cols)
           END IF
        ELSE
           WRITE (unit=unit_nr,fmt="('    n_local_cols=*null*,')")
        END IF

        IF (ASSOCIATED(distribution_2d%local_cols)) THEN
           WRITE (unit=unit_nr,fmt="('      local_cols=(')")
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
           WRITE (unit=unit_nr,fmt="('      local_cols=*null*,')")
        END IF

        IF (ASSOCIATED(distribution_2d%blacs_env)) THEN
           IF (my_long_description) THEN
              WRITE (unit=unit_nr,fmt="('    blacs_env=')",advance="no")
              CALL cp_blacs_env_write(distribution_2d%blacs_env,unit_nr=unit_nr,&
                   error=error)
           ELSE
              WRITE (unit=unit_nr,fmt="('    blacs_env=<blacs_env id=',i6,'>')")&
                   distribution_2d%blacs_env%group
           END IF
        ELSE
           WRITE (unit=unit_nr,fmt="('    blacs_env=*null*')")
        END IF

        WRITE (unit=unit_nr,fmt="(' }')")
     END IF

  ELSE IF (my_local) THEN
     WRITE (unit=unit_nr,&
          fmt="(' <distribution_2d *null*>')")
  END IF

  CALL m_flush(unit_nr)

END SUBROUTINE distribution_2d_write

! *****************************************************************************
!> \brief returns various attributes about the distribution_2d
!> \param distribution_2d the object you want info about
!> \param row_distribution ...
!> \param col_distribution ...
!> \param n_row_distribution ...
!> \param n_col_distribution ...
!> \param n_local_rows ...
!> \param n_local_cols ...
!> \param local_rows ...
!> \param local_cols ...
!> \param flat_local_rows ...
!> \param flat_local_cols ...
!> \param n_flat_local_rows ...
!> \param n_flat_local_cols ...
!> \param blacs_env ...
!> \param id_nr ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!>      see distribution_2d_type attributes for a description of the other
!>      arguments
!> \par History
!>      09.2003 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE distribution_2d_get(distribution_2d, row_distribution, &
     col_distribution, n_row_distribution, n_col_distribution,&
     n_local_rows, n_local_cols, local_rows, local_cols,&
     flat_local_rows, flat_local_cols, n_flat_local_rows, n_flat_local_cols,&
     blacs_env, id_nr,&
     error)
    TYPE(distribution_2d_type), POINTER      :: distribution_2d
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: row_distribution, &
                                                col_distribution
    INTEGER, INTENT(out), OPTIONAL           :: n_row_distribution, &
                                                n_col_distribution
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: n_local_rows, n_local_cols
    TYPE(cp_1d_i_p_type), DIMENSION(:), &
      OPTIONAL, POINTER                      :: local_rows, local_cols
    INTEGER, DIMENSION(:), OPTIONAL, POINTER :: flat_local_rows, &
                                                flat_local_cols
    INTEGER, INTENT(out), OPTIONAL           :: n_flat_local_rows, &
                                                n_flat_local_cols
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env
    INTEGER, INTENT(out), OPTIONAL           :: id_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'distribution_2d_get', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iblock_atomic, iblock_min, &
                                                ikind, ikind_min, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: multiindex
    LOGICAL                                  :: failure

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

END MODULE distribution_2d_types
