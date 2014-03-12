!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief This module returns additional info on PW grids
!> \par History
!>      JGH (09-06-2007) : Created from pw_grids
!> \author JGH
! *****************************************************************************
MODULE pw_grid_info
  USE f77_blas
  USE fft_tools,                       ONLY: FFT_RADIX_NEXT,&
                                             FFT_RADIX_NEXT_ODD,&
                                             fft_radix_operations
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE pw_grid_types,                   ONLY: pw_grid_type
#include "../cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pw_find_cutoff, pw_grid_init_setup, pw_grid_bounds_from_n

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_grid_info'

CONTAINS

! *****************************************************************************
  FUNCTION pw_grid_init_setup ( hmat, cutoff, spherical, odd, fft_usage, ncommensurate,&
       icommensurate, ref_grid, n_orig, error) RESULT(n)

    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: hmat
    REAL(KIND=dp), INTENT(IN)                :: cutoff
    LOGICAL, INTENT(IN)                      :: spherical, odd, fft_usage
    INTEGER, INTENT(IN)                      :: ncommensurate, icommensurate
    TYPE(pw_grid_type), INTENT(IN), OPTIONAL :: ref_grid
    INTEGER, INTENT(IN), OPTIONAL            :: n_orig(3)
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, DIMENSION(3)                    :: n

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_grid_init_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: my_icommensurate
    LOGICAL                                  :: failure

    IF (ncommensurate>0) THEN
       my_icommensurate=icommensurate
       CPPrecondition(icommensurate>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(icommensurate<=ncommensurate,cp_failure_level,routineP,error,failure)
    ELSE
       my_icommensurate=0
    END IF

    IF (my_icommensurate>1) THEN
       CPPostcondition(PRESENT(ref_grid),cp_failure_level,routineP,error,failure)
       n=ref_grid%npts/2**(my_icommensurate-1)
       CPPostcondition(ALL(ref_grid%npts==n*2**(my_icommensurate-1)),cp_failure_level,routineP,error,failure)
       CPPostcondition(ALL(pw_grid_n_for_fft(n,error=error)==n),cp_failure_level,routineP,error,failure)
    ELSE
       n=pw_grid_find_n(hmat,cutoff=cutoff, fft_usage=fft_usage, ncommensurate=ncommensurate,&
            spherical=spherical, odd=odd, n_orig=n_orig, error=error)
    END IF

  END FUNCTION pw_grid_init_setup

! *****************************************************************************
!> \brief returns the n needed for the grid with all the given constraints
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  FUNCTION pw_grid_find_n ( hmat, cutoff, fft_usage, spherical, odd, ncommensurate,&
       n_orig, error) RESULT(n)

    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: hmat
    REAL(KIND=dp), INTENT(IN)                :: cutoff
    LOGICAL, INTENT(IN)                      :: fft_usage, spherical, odd
    INTEGER, INTENT(IN)                      :: ncommensurate
    INTEGER, INTENT(IN), OPTIONAL            :: n_orig(3)
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, DIMENSION(3)                    :: n

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_grid_find_n', &
      routineP = moduleN//':'//routineN

    INTEGER :: idir, my_icommensurate, my_ncommensurate, nsubgrid, &
      nsubgrid_new, ntest(3), t_icommensurate
    LOGICAL                                  :: failure, ftest, subgrid_is_OK

    failure=.FALSE.
    ! ncommensurate is the number of commensurate grids
    ! in order to have non-commensurate grids ncommensurate must be 0
    ! icommensurte  is the level number of communensurate grids
    ! this implies that the number of grid points in each direction
    ! is k*2**(ncommensurate-icommensurate)
    my_ncommensurate=ncommensurate
    IF (my_ncommensurate > 0 ) THEN
       my_icommensurate=1
    ELSE
       my_icommensurate=0
    ENDIF
    CPPrecondition(my_icommensurate <= my_ncommensurate,cp_failure_level,routineP,error,failure)
    CPPrecondition((my_icommensurate > 0 .OR. my_ncommensurate <= 0),cp_failure_level,routineP,error,failure)
    CPPrecondition(my_ncommensurate >= 0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(n_orig)) THEN
       n = n_orig
    ELSE
       CPPrecondition(cutoff>0.0_dp,cp_failure_level,routineP,error,failure)
       n= pw_grid_n_from_cutoff (hmat, cutoff, error=error)
    END IF

    IF (fft_usage) THEN
       n=pw_grid_n_for_fft(n,odd=odd,error=error)

       IF (.NOT.spherical) THEN
          ntest = n

          IF ( my_ncommensurate>0 ) THEN
             DO idir=1,3
                DO
                   ! find valid radix >= ntest
                   CALL fft_radix_operations ( ntest(idir), n(idir), FFT_RADIX_NEXT )
                   ! check every subgrid of n
                   subgrid_is_OK = .TRUE.
                   DO t_icommensurate=1, my_ncommensurate-1
                      nsubgrid=n(idir)/2**(my_ncommensurate-t_icommensurate)
                      CALL fft_radix_operations ( nsubgrid,nsubgrid_new, FFT_RADIX_NEXT )
                      subgrid_is_OK = (nsubgrid == nsubgrid_new) .AND. &
                                      (MODULO(n(idir),2**(my_ncommensurate-t_icommensurate)) == 0)
                      IF (.NOT. subgrid_is_OK) EXIT
                   END DO
                   IF (subgrid_is_OK) THEN
                      EXIT
                   ELSE
                      ! subgrid wasn't OK, increment ntest and try again
                      ntest(idir)=n(idir)+1
                   ENDIF
                END DO
             END DO
          END IF
       END IF
    ELSE
       ! without a cutoff and HALFSPACE we have to be sure that there is
       ! a negative counterpart to every g vector (-> odd number of grid points)
       IF ( odd ) n = n + MOD ( n + 1, 2 )

    END IF

    ! final check if all went fine ...
    IF (my_ncommensurate>0) THEN
       DO my_icommensurate=1,my_ncommensurate
          ftest = ANY( MODULO(n,2**(my_ncommensurate-my_icommensurate)).NE.0 )
          CPPostcondition(.NOT. ftest,cp_failure_level,routineP,error,failure)
       END DO
    ENDIF

  END FUNCTION pw_grid_find_n

! *****************************************************************************
!> \brief returns the closest number of points >= n, on which you can perform
!>      ffts
!> \param n the minimum number of points you want
!> \param odd if the number has to be odd
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \note
!>      result<=n
!> \author fawzi
! *****************************************************************************
  FUNCTION pw_grid_n_for_fft(n,odd,error) RESULT(nout)
    INTEGER, DIMENSION(3), INTENT(in)        :: n
    LOGICAL, INTENT(in), OPTIONAL            :: odd
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, DIMENSION(3)                    :: nout

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_grid_n_for_fft', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, my_odd

    failure=.FALSE.
    my_odd=.FALSE.
    IF (PRESENT(odd)) my_odd=odd
    CPPrecondition(ALL(n>=0),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       IF (my_odd) THEN
          CALL fft_radix_operations ( n(1), nout(1), FFT_RADIX_NEXT_ODD )
          CALL fft_radix_operations ( n(2), nout(2), FFT_RADIX_NEXT_ODD )
          CALL fft_radix_operations ( n(3), nout(3), FFT_RADIX_NEXT_ODD )
       ELSE
          CALL fft_radix_operations ( n(1), nout(1), FFT_RADIX_NEXT )
          CALL fft_radix_operations ( n(2), nout(2), FFT_RADIX_NEXT )
          CALL fft_radix_operations ( n(3), nout(3), FFT_RADIX_NEXT )
       END IF
    END IF

  END FUNCTION pw_grid_n_for_fft

! *****************************************************************************
!> \brief Find the number of points that give at least the requested cutoff
!> \par History
!>      JGH (21-12-2000) : Simplify parameter list, bounds will be global
!>      JGH ( 8-01-2001) : Add check to FFT allowd grids (this now depends
!>                         on the FFT library.
!>                         Should the pw_grid_type have a reference to the FFT
!>                         library ?
!>      JGH (28-02-2001) : Only do conditional check for FFT
!>      JGH (21-05-2002) : Optimise code, remove orthorhombic special case
!> \author apsi
!>      Christopher Mundy
! *****************************************************************************
  FUNCTION pw_grid_n_from_cutoff ( hmat, cutoff, error ) RESULT(n)

    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: hmat
    REAL(KIND=dp), INTENT(IN)                :: cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, DIMENSION(3)                    :: n

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_grid_n_from_cutoff', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alat( 3 )

    failure=.FALSE.

    DO i=1,3
       alat(i) = SUM ( hmat(:,i)**2 )
    ENDDO
    CPPostcondition(ALL(alat/=0._dp),cp_failure_level,routineP,error,failure)
    IF ( failure ) THEN
       n = -HUGE(0)
    ELSE
       n = 2*FLOOR ( SQRT ( 2.0_dp * cutoff * alat ) / twopi ) + 1
    END IF

  END FUNCTION pw_grid_n_from_cutoff

! *****************************************************************************
!> \brief returns the bounds that distribute n points evenly around 0
!> \param npts the number of points in each direction
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fawzi
! *****************************************************************************
  FUNCTION pw_grid_bounds_from_n(npts,error) RESULT(bounds)
    INTEGER, DIMENSION(3), INTENT(in)        :: npts
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, DIMENSION(2, 3)                 :: bounds

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_grid_bounds_from_n', &
      routineP = moduleN//':'//routineN

    bounds(1,:)=-npts/2
    bounds(2,:)=bounds(1,:)+npts-1

  END FUNCTION pw_grid_bounds_from_n

! *****************************************************************************
!> \brief Given a grid and a box, calculate the corresponding cutoff
!>      *** This routine calculates the cutoff in MOMENTUM UNITS! ***
!> \note
!>      This routine is local. It works independent from the distribution
!>      of PW on processors.
!>      npts is the grid size for the full box.
!> \par History
!>      JGH (20-12-2000) : Deleted some strange comments
!> \author apsi
!>      Christopher Mundy
! *****************************************************************************
  FUNCTION pw_find_cutoff ( npts, h_inv, error ) RESULT(cutoff)

    INTEGER, DIMENSION(:), INTENT(IN)        :: npts
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: h_inv
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    REAL(KIND=dp)                            :: cutoff

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_find_cutoff', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: gcut, gdum( 3 ), length

! compute 2*pi*h_inv^t*g  where g = (nmax[1],0,0)

    gdum(:) = twopi * h_inv(1,:) * REAL((npts(1)-1)/2,KIND=dp)
    length = SQRT ( gdum(1)**2 + gdum(2)**2 + gdum(3)**2 )
    gcut = length

    ! compute 2*pi*h_inv^t*g  where g = (0,nmax[2],0)
    gdum(:) = twopi * h_inv(2,:) * REAL((npts(2)-1)/2,KIND=dp)
    length = SQRT ( gdum(1)**2 + gdum(2)**2 + gdum(3)**2 )
    gcut = MIN ( gcut, length )

    ! compute 2*pi*h_inv^t*g  where g = (0,0,nmax[3])
    gdum(:) = twopi * h_inv(3,:) * REAL((npts(3)-1)/2,KIND=dp)
    length = SQRT ( gdum(1)**2 + gdum(2)**2 + gdum(3)**2 )
    gcut = MIN ( gcut, length )

    cutoff = gcut - 1.e-8_dp

  END FUNCTION pw_find_cutoff

END MODULE pw_grid_info

