!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \note
!>      If parallel mode is distributed certain combination of
!>      "in_use" and "in_space" can not be used.
!>      For performance reasons it would be better to have the loops
!>      over g-vectros in the gather/scatter routines in new subprograms
!>      with the actual arrays (also the adressing) in the parameter list
!> \par History
!>      JGH (29-Dec-2000) : Changes for parallel use
!>      JGH (13-Mar-2001) : added timing calls
!>      JGH (26-Feb-2003) : OpenMP enabled
!>      JGH (17-Nov-2007) : Removed mass arrays
!>      JGH (01-Dec-2007) : Removed and renamed routines
!>      03.2008 [tlaino] : Splitting pw_types into pw_types and pw_methods
!> \author apsi
! *****************************************************************************
MODULE pw_types
  
  USE kinds,                           ONLY: dp
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_release,&
                                             pw_grid_retain
  USE timings,                         ONLY: print_stack,&
                                             timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pw_type, pw_p_type
  PUBLIC :: pw_retain, pw_release, pw_create

! *****************************************************************************
  TYPE pw_type
     REAL (KIND=dp), DIMENSION ( : ), POINTER :: cr
     REAL (KIND=dp), DIMENSION ( :, :, : ), POINTER :: cr3d
     COMPLEX (KIND=dp), DIMENSION ( : ), POINTER :: cc
     COMPLEX (KIND=dp), DIMENSION ( :, :, : ), POINTER :: cc3d

     INTEGER :: in_use              ! Which data is used [r1d/c1d/r3d/c3d]
     INTEGER :: in_space            ! Real/Reciprocal space
     INTEGER :: id_nr               ! unique identifier
     INTEGER :: ref_count           ! reference count

     TYPE ( pw_grid_type ), POINTER :: pw_grid
  END TYPE pw_type

! *****************************************************************************
  TYPE pw_p_type
     TYPE(pw_type), POINTER :: pw
  END TYPE pw_p_type

  ! Flags for the structure member 'in_use'
  INTEGER, PARAMETER, PUBLIC :: REALDATA1D = 301, COMPLEXDATA1D = 302
  INTEGER, PARAMETER, PUBLIC :: REALDATA3D = 303, COMPLEXDATA3D = 304, NODATA = 305

  ! Flags for the structure member 'in_space'
  INTEGER, PARAMETER, PUBLIC :: NOSPACE = 371, REALSPACE = 372, RECIPROCALSPACE = 373
  INTEGER, PUBLIC, PARAMETER :: SQUARE = 391, SQUAREROOT = 392

  ! to generate unique id_nr
  INTEGER, SAVE, PRIVATE :: last_pw_id_nr=0
  INTEGER, SAVE, PRIVATE :: allocated_pw_count=0

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_types'
  LOGICAL, PARAMETER, PRIVATE :: debug_this_module=.FALSE.

CONTAINS

! *****************************************************************************
!> \brief retains a pw type
!> \param pw the pw to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      03.2003 created [fawzi]
!> \author fawzi
!> \note
!>      see doc/ReferenceCounting.html
! *****************************************************************************
  SUBROUTINE pw_retain(pw, error)
    TYPE(pw_type), POINTER                   :: pw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pw),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(pw%ref_count>0,cp_failure_level,routineP,error)
       pw%ref_count=pw%ref_count+1
    END IF
  END SUBROUTINE pw_retain

! *****************************************************************************
!> \brief releases the given pw
!> \param pw the pw to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
!> \note
!>      see doc/ReferenceCounting.html
! *****************************************************************************
SUBROUTINE pw_release(pw, error)
    TYPE(pw_type), POINTER                   :: pw
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(pw)) THEN
     CPPreconditionNoFail(pw%ref_count>0,cp_failure_level,routineP,error)
     pw%ref_count=pw%ref_count-1
     IF (pw%ref_count==0) THEN
        pw%ref_count=1

        allocated_pw_count = allocated_pw_count - 1
        SELECT CASE(pw % in_use)
        CASE (REALDATA1D)
           DEALLOCATE ( pw % cr, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        CASE(COMPLEXDATA1D)
           DEALLOCATE ( pw % cc, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        CASE(REALDATA3D)
           IF (ASSOCIATED(pw%cr3d)) THEN
              !FM optimizations of pools might have removed the 3d field to cache it
              DEALLOCATE ( pw % cr3d, STAT = stat )
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
        CASE(COMPLEXDATA3D)
           DEALLOCATE ( pw % cc3d, STAT = stat )
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        CASE(NODATA)
        CASE default
           CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                "unknown data type "//cp_to_string(pw%in_use),error,failure)
        END SELECT
        CALL pw_grid_release(pw%pw_grid, error=error)
        pw%ref_count=0
        DEALLOCATE(pw, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(pw)
END SUBROUTINE pw_release

! *****************************************************************************
!> \brief allocates and initializes pw_type
!> \param pw the type that will bw allocated and initialized
!> \param pw_grid ...
!> \param use_data which kind of data will be used
!> \param in_space in which space the pw is (real or reciprocal)
!> \param cr3d_ptr pointer with the cr3d data (make sense only if
!>        use_data==REALDATA3D)
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE pw_create(pw, pw_grid, use_data, in_space, cr3d_ptr, error)
    TYPE(pw_type), POINTER                   :: pw
    TYPE(pw_grid_type), POINTER              :: pw_grid
    INTEGER, INTENT(in)                      :: use_data
    INTEGER, INTENT(in), OPTIONAL            :: in_space
    REAL(KIND=dp), DIMENSION(:, :, :), &
      OPTIONAL, POINTER                      :: cr3d_ptr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    INTEGER, DIMENSION(:, :), POINTER        :: bounds
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

  failure=.FALSE.

  CALL timeset(routineN,handle)
  CPPrecondition(.NOT.ASSOCIATED(pw),cp_failure_level,routineP,error,failure)
  ALLOCATE(pw,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
  logger => cp_error_get_logger(error)
  IF (debug_this_module) THEN
     WRITE (cp_logger_get_default_unit_nr(logger),"('*** allocated pw ***')")
     IF (PRESENT(cr3d_ptr)) THEN
        IF (ASSOCIATED(cr3d_ptr)) THEN
           WRITE (cp_logger_get_default_unit_nr(logger),"('*** cr3d associated ***')")
        END IF
     END IF
     CALL print_stack(cp_logger_get_default_unit_nr(logger))
  END IF

  IF (PRESENT(cr3d_ptr)) THEN
     IF (ASSOCIATED(cr3d_ptr)) THEN
        CPAssertNoFail(use_data==REALDATA3D,cp_failure_level,routineP,error)
     END IF
  END IF

  last_pw_id_nr=last_pw_id_nr+1
  pw % id_nr = last_pw_id_nr
  pw % ref_count = 1
  NULLIFY ( pw % pw_grid )
  pw % in_use = use_data
  pw % pw_grid => pw_grid
  CALL pw_grid_retain(pw%pw_grid, error=error)
  pw % in_space = NOSPACE
  bounds => pw % pw_grid % bounds_local

  allocated_pw_count = allocated_pw_count + 1

  NULLIFY ( pw % cr, pw % cc, pw % cr3d, pw % cc3d )

  SELECT CASE(use_data)
  CASE(REALDATA1D)
     ALLOCATE ( pw % cr ( pw % pw_grid % ngpts_cut_local ), STAT = stat )
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

  CASE(COMPLEXDATA1D)
     ALLOCATE ( pw % cc ( pw % pw_grid % ngpts_cut_local ), STAT = stat )
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)

  CASE (REALDATA3D)
     IF (PRESENT(cr3d_ptr)) THEN
        IF (ASSOCIATED(cr3d_ptr)) THEN
           IF (ALL(bounds(1,:) <= bounds(2,:))) THEN
             CPPreconditionNoFail(ALL(LBOUND(cr3d_ptr)==bounds(1,:)),cp_failure_level,routineP,error)
             CPPreconditionNoFail(ALL(UBOUND(cr3d_ptr)==bounds(2,:)),cp_failure_level,routineP,error)
           END IF
           pw%cr3d => cr3d_ptr
        END IF
     END IF
     IF (.NOT.ASSOCIATED(pw%cr3d)) THEN
        ALLOCATE ( pw % cr3d ( &
             bounds ( 1, 1 ) : bounds ( 2, 1 ), &
             bounds ( 1, 2 ) : bounds ( 2, 2 ), &
             bounds ( 1, 3 ) : bounds ( 2, 3 ) ), STAT = stat )
        CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
     END IF

  CASE(COMPLEXDATA3D)
     ALLOCATE ( pw % cc3d ( &
          bounds ( 1, 1 ) : bounds ( 2, 1 ), &
          bounds ( 1, 2 ) : bounds ( 2, 2 ), &
          bounds ( 1, 3 ) : bounds ( 2, 3 ) ), STAT = stat )
     CPPostcondition(stat==0,cp_fatal_level,routineP,error,failure)
  CASE(NODATA)
  CASE default
     CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
          "unknown data type",error,failure)
  END SELECT
  IF (PRESENT(in_space)) pw%in_space=in_space
  END IF
  CALL timestop(handle)
END SUBROUTINE pw_create

END MODULE pw_types
