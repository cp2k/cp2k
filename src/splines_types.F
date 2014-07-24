!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines for handling splines_types
!> \par History
!>      2001-09-21-HAF added this doc entry and changed formatting
!> \author various
! *****************************************************************************
MODULE splines_types

  USE f77_blas
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'splines_types'
  INTEGER, PRIVATE, SAVE :: last_spline_env_id_nr=0
  INTEGER, PRIVATE, SAVE :: last_spline_data_id_nr=0

  PUBLIC :: spline_env_release, spline_environment_type
  PUBLIC :: spline_env_create, spline_data_p_type
  PUBLIC :: spline_data_create, spline_data_p_copy
  PUBLIC :: spline_data_retain, spline_data_p_retain
  PUBLIC :: spline_data_release, spline_data_p_release
  PUBLIC :: spline_factor_copy, spline_factor_create, spline_factor_release
  PUBLIC :: spline_data_type  ! the data structure for spline table
  PUBLIC :: spline_factor_type  ! the multiplicative factors for splines

! *****************************************************************************
!> \brief Data-structure that holds all needed information about
!>      a specific spline interpolation.
!> \par History
!>      2001-09-19-HAF added this doc entry and changed formatting
!> \author unknown
! *****************************************************************************
  TYPE spline_data_type
     INTEGER :: ref_count, id_nr
     REAL (KIND=dp), POINTER :: y  ( : )    ! the function values y(x)
     REAL (KIND=dp), POINTER :: y2 ( : )    ! the 2nd derivative via interpolation
     INTEGER                 :: n           ! dimension of above arrays
     ! not used if uniform increments
     REAL (KIND=dp)          :: h           ! uniform increment of x if applicable
     REAL (KIND=dp)          :: invh        ! inverse of h
     REAL (KIND=dp)          :: h26         ! 1/6 * h**2 if uniform increments
     ! 1/6 otherwise
     REAL (KIND=dp)          :: x1          ! starting x value if uniform incr.
     REAL (KIND=dp)          :: xn          ! end x value if uniform incr.
  END TYPE spline_data_type

! *****************************************************************************
  TYPE spline_data_p_type
     TYPE ( spline_data_type ), POINTER :: spline_data
  END TYPE spline_data_p_type

! *****************************************************************************
  TYPE spline_data_pp_type
     TYPE ( spline_data_p_type ), POINTER, DIMENSION(:)  :: spl_p
  END TYPE spline_data_pp_type

! *****************************************************************************
  TYPE spline_environment_type
     INTEGER :: ref_count, id_nr
     TYPE ( spline_data_pp_type ), POINTER, DIMENSION ( : ) :: spl_pp
     INTEGER, POINTER, DIMENSION ( :, : ) :: spltab
  END TYPE spline_environment_type

! *****************************************************************************
  TYPE spline_factor_type
     REAL(KIND=dp) :: rcutsq_f, cutoff
     REAL(KIND=dp), DIMENSION(:), POINTER :: rscale
     REAL(KIND=dp), DIMENSION(:), POINTER :: fscale
     REAL(KIND=dp), DIMENSION(:), POINTER :: dscale
  END TYPE spline_factor_type

CONTAINS

! *****************************************************************************
!> \brief releases spline_env
!> \param spline_env ...
!> \param error ...
!> \author unknown
! *****************************************************************************
  SUBROUTINE spline_env_release(spline_env,error)
    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    TYPE(spline_data_p_type), DIMENSION(:), &
      POINTER                                :: spl_p

    failure=.FALSE.

    IF (ASSOCIATED(spline_env)) THEN
       CPPrecondition(spline_env%ref_count>0,cp_failure_level,routineP,error,failure)
       spline_env%ref_count=spline_env%ref_count-1
       IF (spline_env%ref_count<1) THEN
          DEALLOCATE (  spline_env % spltab, stat = stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          DO i = 1, SIZE (  spline_env % spl_pp )
             spl_p => spline_env % spl_pp ( i ) % spl_p
             CALL spline_data_p_release ( spl_p, error )
          END DO
          DEALLOCATE (  spline_env % spl_pp, stat = stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          NULLIFY ( spline_env % spltab )
          NULLIFY ( spline_env % spl_pp )
          DEALLOCATE (  spline_env , stat = stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF

  END SUBROUTINE spline_env_release

! *****************************************************************************
!> \brief releases spline_data
!> \param spline_data ...
!> \param error ...
!> \author CJM
! *****************************************************************************
  SUBROUTINE spline_data_release (spline_data,error)
    TYPE(spline_data_type), POINTER          :: spline_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_data_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(spline_data)) THEN
       CPPrecondition(spline_data%ref_count>0,cp_failure_level,routineP,error,failure)
       spline_data%ref_count=spline_data%ref_count-1
       IF (spline_data%ref_count<1) THEN
          IF ( ASSOCIATED ( spline_data%y)) THEN
             DEALLOCATE (spline_data%y,STAT=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             NULLIFY (spline_data%y)
          END IF
          IF ( ASSOCIATED ( spline_data%y2)) THEN
             DEALLOCATE (spline_data%y2,STAT=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             NULLIFY (spline_data%y2)
          END IF
          DEALLOCATE (  spline_data, stat = stat )
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
  END SUBROUTINE spline_data_release

! *****************************************************************************
!> \brief releases spline_data_p
!> \param spl_p ...
!> \param error ...
!> \author CJM
! *****************************************************************************
  SUBROUTINE spline_data_p_release (spl_p,error)
    TYPE(spline_data_p_type), DIMENSION(:), &
      POINTER                                :: spl_p
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_data_p_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure, release_kind

    failure=.FALSE.

    IF (ASSOCIATED(spl_p)) THEN
       release_kind = .TRUE.
       DO i = 1, SIZE(spl_p)
          CALL spline_data_release(spl_p(i)%spline_data,error)
          release_kind = release_kind .AND.(.NOT.ASSOCIATED(spl_p(i)%spline_data))
       END DO
       IF (release_kind) THEN
          DEALLOCATE(spl_p, stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF

  END SUBROUTINE spline_data_p_release

! *****************************************************************************
!> \brief retains spline_env
!> \param spline_data ...
!> \param error ...
!> \author CJM
! *****************************************************************************
  SUBROUTINE spline_data_retain(spline_data,error)
    TYPE(spline_data_type), POINTER          :: spline_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_data_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(spline_data),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(spline_data%ref_count>0,cp_failure_level,routineP,error,failure)
       spline_data%ref_count=spline_data%ref_count+1
    END IF
  END SUBROUTINE spline_data_retain

! *****************************************************************************
!> \brief retains spline_data_p_type
!> \param spl_p ...
!> \param error ...
!> \author CJM
! *****************************************************************************
  SUBROUTINE spline_data_p_retain(spl_p,error)
    TYPE(spline_data_p_type), DIMENSION(:), &
      POINTER                                :: spl_p
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_data_p_retain', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(spl_p),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       DO i = 1, SIZE(spl_p)
          CALL spline_data_retain(spl_p(i)%spline_data, error)
       END DO
    END IF
  END SUBROUTINE spline_data_p_retain

! *****************************************************************************
!> \brief retains spline_env
!> \param spline_env ...
!> \param error ...
!> \par History
!>      2001-09-19-HAF added this doc entry and changed formatting
!> \author unknown
! *****************************************************************************
  SUBROUTINE spline_env_retain(spline_env,error)
    TYPE(spline_environment_type), POINTER   :: spline_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(spline_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(spline_env%ref_count>0,cp_failure_level,routineP,error,failure)
       spline_env%ref_count=spline_env%ref_count+1
    END IF
  END SUBROUTINE spline_env_retain

! *****************************************************************************
!> \brief Data-structure that holds all needed information about
!>      a specific spline interpolation.
!> \param spline_env ...
!> \param ntype ...
!> \param ntab_in ...
!> \param error ...
!> \par History
!>      2001-09-19-HAF added this doc entry and changed formatting
!> \author unknown
! *****************************************************************************
  SUBROUTINE spline_env_create ( spline_env, ntype, ntab_in, error )
    TYPE(spline_environment_type), POINTER   :: spline_env
    INTEGER, INTENT(IN)                      :: ntype
    INTEGER, INTENT(IN), OPTIONAL            :: ntab_in
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, isize, j, ntab, &
                                                stat
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)

    failure=.FALSE.

    ALLOCATE(spline_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY ( spline_env % spl_pp )
       NULLIFY ( spline_env % spltab )
       spline_env%ref_count=1
       last_spline_env_id_nr=last_spline_env_id_nr+1
       spline_env%id_nr=last_spline_env_id_nr
       ! Allocate the number of spline data tables (upper triangular)
       IF (PRESENT(ntab_in)) THEN
          ntab = ntab_in
       ELSE
          ntab = ( ntype * ntype + ntype ) / 2
       END IF
       ALLOCATE ( spline_env % spl_pp( ntab ),  stat=stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ALLOCATE ( spline_env % spltab ( ntype, ntype ), stat=stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       DO i = 1, ntab
          NULLIFY (spline_env%spl_pp( i )%spl_p)
          isize = 1
          ALLOCATE(spline_env%spl_pp(i)%spl_p(isize), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO j = 1, SIZE(spline_env%spl_pp(i)%spl_p)
             CALL spline_data_create (spline_env%spl_pp(i)%spl_p(j)%spline_data,error=error)
          END DO
       END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE spline_env_create

! *****************************************************************************
!> \brief Copy Data-structure of spline_data_p_type
!> \param spl_p_source ...
!> \param spl_p_dest ...
!> \param error ...
!> \author teo 06.2007
! *****************************************************************************
  SUBROUTINE spline_data_p_copy ( spl_p_source, spl_p_dest, error )
    TYPE(spline_data_p_type), DIMENSION(:), &
      POINTER                                :: spl_p_source, spl_p_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_data_p_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nsized, nsizes, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPostcondition(ASSOCIATED(spl_p_source),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       nsizes = SIZE(spl_p_source)
       IF (.NOT.ASSOCIATED(spl_p_dest)) THEN
          ALLOCATE(spl_p_dest(nsizes), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i = 1, nsizes
             NULLIFY(spl_p_dest(i)%spline_data)
          END DO
       ELSE
          nsized = SIZE(spl_p_dest)
          CPPostcondition(nsizes==nsized,cp_failure_level,routineP,error,failure)
          DO i = 1, nsizes
             CALL spline_data_release(spl_p_dest(i)%spline_data, error=error)
          END DO
       END IF
       DO i = 1, nsizes
          CALL spline_data_copy (spl_p_source(i)%spline_data, spl_p_dest(i)%spline_data, error )
       END DO
    END IF
  END SUBROUTINE spline_data_p_copy

! *****************************************************************************
!> \brief Copy Data-structure that constains spline table
!> \param spline_data_source ...
!> \param spline_data_dest ...
!> \param error ...
!> \author teo 11.2005
! *****************************************************************************
  SUBROUTINE spline_data_copy ( spline_data_source, spline_data_dest, error )
    TYPE(spline_data_type), POINTER          :: spline_data_source, &
                                                spline_data_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_data_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPostcondition(ASSOCIATED(spline_data_source),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (.NOT.ASSOCIATED(spline_data_dest)) CALL spline_data_create(spline_data_dest,error)

       spline_data_dest%ref_count = spline_data_source%ref_count
       spline_data_dest%id_nr     = spline_data_source%id_nr
       spline_data_dest%n         = spline_data_source%n
       spline_data_dest%h         = spline_data_source%h
       spline_data_dest%invh      = spline_data_source%invh
       spline_data_dest%h26       = spline_data_source%h26
       spline_data_dest%x1        = spline_data_source%x1
       spline_data_dest%xn        = spline_data_source%xn
       IF (ASSOCIATED(spline_data_source%y)) THEN
          ALLOCATE(spline_data_dest%y(SIZE(spline_data_source%y)), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          spline_data_dest%y = spline_data_source%y
       END IF
       IF (ASSOCIATED(spline_data_source%y2)) THEN
          ALLOCATE(spline_data_dest%y2(SIZE(spline_data_source%y2)), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          spline_data_dest%y2 = spline_data_source%y2
       END IF
    END IF
  END SUBROUTINE spline_data_copy

! *****************************************************************************
!> \brief Data-structure that constains spline table
!> \param spline_data ...
!> \param error ...
!> \author unknown
! *****************************************************************************
  SUBROUTINE spline_data_create ( spline_data, error )
    TYPE(spline_data_type), POINTER          :: spline_data
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_data_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(spline_data, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       spline_data%ref_count=1
       last_spline_data_id_nr=last_spline_data_id_nr+1
       spline_data%id_nr = last_spline_data_id_nr
       NULLIFY ( spline_data % y )
       NULLIFY ( spline_data % y2 )
    END IF
  END SUBROUTINE spline_data_create

! *****************************************************************************
!> \brief releases spline_factor
!> \param spline_factor ...
!> \param error ...
!> \author teo
! *****************************************************************************
  SUBROUTINE spline_factor_release (spline_factor,error)
    TYPE(spline_factor_type), POINTER        :: spline_factor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_factor_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(spline_factor)) THEN
       IF (ASSOCIATED(spline_factor%rscale)) THEN
          DEALLOCATE(spline_factor%rscale,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(spline_factor%fscale)) THEN
          DEALLOCATE(spline_factor%fscale,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(spline_factor%dscale)) THEN
          DEALLOCATE(spline_factor%dscale,stat=stat)
          CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       DEALLOCATE (spline_factor, stat=stat )
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
  END SUBROUTINE spline_factor_release

! *****************************************************************************
!> \brief releases spline_factor
!> \param spline_factor ...
!> \param error ...
!> \author teo
! *****************************************************************************
  SUBROUTINE spline_factor_create (spline_factor, error)
    TYPE(spline_factor_type), POINTER        :: spline_factor
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_factor_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(spline_factor),cp_failure_level,routineP,error,failure)
    ALLOCATE(spline_factor, stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(spline_factor%rscale(1),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(spline_factor%fscale(1),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(spline_factor%dscale(1),stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    spline_factor%rscale = 1.0_dp
    spline_factor%fscale = 1.0_dp
    spline_factor%dscale = 1.0_dp
    spline_factor%rcutsq_f = 1.0_dp
    spline_factor%cutoff   = 0.0_dp
  END SUBROUTINE spline_factor_create

! *****************************************************************************
!> \brief releases spline_factor
!> \param spline_factor_source ...
!> \param spline_factor_dest ...
!> \param error ...
!> \author teo
! *****************************************************************************
  SUBROUTINE spline_factor_copy (spline_factor_source,spline_factor_dest, error)
    TYPE(spline_factor_type), POINTER        :: spline_factor_source, &
                                                spline_factor_dest
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'spline_factor_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: isize, jsize, ksize
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(spline_factor_dest)) CALL spline_factor_release(spline_factor_dest, error)
    IF (ASSOCIATED(spline_factor_source)) THEN
       isize = SIZE(spline_factor_source%rscale)
       jsize = SIZE(spline_factor_source%fscale)
       ksize = SIZE(spline_factor_source%dscale)
       CPPrecondition(isize==jsize,cp_failure_level,routineP,error,failure)
       CPPrecondition(isize==ksize,cp_failure_level,routineP,error,failure)
       CALL spline_factor_create(spline_factor_dest, error)
       spline_factor_dest%rscale = spline_factor_source%rscale
       spline_factor_dest%fscale = spline_factor_source%fscale
       spline_factor_dest%dscale = spline_factor_source%dscale
       spline_factor_dest%rcutsq_f = spline_factor_source%rcutsq_f
       spline_factor_dest%cutoff   = spline_factor_source%cutoff
    END IF
  END SUBROUTINE spline_factor_copy

END MODULE splines_types
