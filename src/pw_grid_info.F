!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!
!!****** cp2k/pw_grid_info [1.1] *
!!
!!   NAME
!!     pw_grid_info
!!
!!   FUNCTION
!!     This module returns additional info on PW grids
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (09-06-2007) : Created from pw_grids
!!
!!***
!******************************************************************************

MODULE pw_grid_info
  USE cell_types,                      ONLY: cell_type
  USE fft_tools,                       ONLY: FFT_RADIX_NEXT,&
                                             FFT_RADIX_NEXT_ODD,&
                                             fft_radix_operations
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE mathlib,                         ONLY: dotprod_3d,&
                                             inv_3x3,&
                                             matvec_3x3
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             HALFSPACE
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: pw_find_cutoff
  PUBLIC :: pw_grid_find_n, pw_grid_n_from_cutoff,&
       pw_grid_n_for_fft, pw_grid_bounds_from_n, pw_cutoff_from_n

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_grid_info'

  INTERFACE pw_find_cutoff
     MODULE PROCEDURE pw_find_cutoff_qs, pw_find_cutoff_fist
  END INTERFACE
!******************************************************************************

CONTAINS

!******************************************************************************
  !!****f* pw_grids/pw_grid_find_n *
  !!
  !!   NAME
  !!     pw_grid_find_n
  !!
  !!   FUNCTION
  !!     returns the n needed for the grid with all the given constraints
  !!
  !!   NOTES
  !!     -
  !!
  !!   ARGUMENTS
  !!     - error: variable to control error logging, stopping,...
  !!       see module cp_error_handling
  !!
  !!   AUTHOR
  !!     fawzi
  !!
  !!*** *********************************************************************
FUNCTION pw_grid_find_n ( cell, cutoff, fft_usage, grid_span,&
     symm_usage, spherical, ncommensurate, error) RESULT(n)

    TYPE(cell_type), INTENT(IN)              :: cell
    REAL(KIND=dp), INTENT(IN)                :: cutoff
    LOGICAL, INTENT(IN), OPTIONAL            :: fft_usage
    INTEGER, INTENT(IN), OPTIONAL            :: grid_span
    LOGICAL, INTENT(IN), OPTIONAL            :: symm_usage, spherical
    INTEGER, INTENT(IN), OPTIONAL            :: ncommensurate
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, DIMENSION(3)                    :: n

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_grid_find_n', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, i, idir, my_grid_span, my_icommensurate, &
      my_ncommensurate, nlowest, nlowest_new, ntest(3)
    LOGICAL                                  :: failure, fft, my_spherical, &
                                                symmetry, ftest

! cutoff could be inout, and set to the effective cutoff of the given grid
! in atomic units
! parent group for this grid
! output unit for information on grid
! has the grid size to be
! compatible with the FFT
! has the grid size to be
! symmetric (g <-> -g)
! block or ray-distribution
!------------------------------------------------------------------------------

  CALL timeset("pw_grid_find_n","I","",handle)
  failure=.FALSE.

  my_grid_span=FULLSPACE
  IF (PRESENT(grid_span)) my_grid_span=grid_span
  ! default is to use only fft compatible grids
  IF ( PRESENT ( fft_usage ) ) THEN
    fft = fft_usage
  ELSE
    fft = .TRUE.
  END IF

  ! default is to use only symmetric grids
  IF ( PRESENT ( symm_usage ) ) THEN
    symmetry = symm_usage
  ELSE
    symmetry = .TRUE.
  END IF

  ! ncommensurate is the number of commensurate grids
  ! in order to have non-commensurate grids ncommensurate must be 0
  ! icommensurte  is the level number of communensurate grids
  ! this implies that the number of grid points in each direction
  ! is k*2**(ncommensurate-icommensurate)
  IF ( PRESENT ( ncommensurate ) ) THEN
     my_ncommensurate=ncommensurate
     IF (my_ncommensurate > 0 ) THEN
         my_icommensurate=1
     ELSE
         my_icommensurate=0
     ENDIF
     CPPrecondition(my_icommensurate <= my_ncommensurate,cp_failure_level,routineP,error,failure)
     CPPrecondition((my_icommensurate > 0 .OR. my_ncommensurate <= 0),cp_failure_level,routineP,error,failure)
     CPPrecondition(my_ncommensurate >= 0,cp_failure_level,routineP,error,failure)
  ELSE
     my_ncommensurate=0
     my_icommensurate=0
  ENDIF

  CPPrecondition(cutoff>0.0_dp,cp_failure_level,routineP,error,failure)
  n= pw_grid_n_from_cutoff ( cell % hmat, cutoff,error=error)
  my_spherical = .FALSE.
  IF (PRESENT(spherical)) my_spherical=spherical

  IF (fft) THEN
     n=pw_grid_n_for_fft(n,odd=(my_grid_span == HALFSPACE .AND. symmetry),error=error)

     IF (.NOT.spherical) THEN
        ntest = n

        IF ( my_ncommensurate>0 ) THEN
           DO idir=1,3
              DO
                 CALL fft_radix_operations ( ntest(idir), n(idir), FFT_RADIX_NEXT )
                 ! is also the lowest grid allowed (e.g could be 17, which is too large, but might be 5)
                 nlowest=n(idir)/2**(my_ncommensurate-my_icommensurate)
                 CALL fft_radix_operations ( nlowest,nlowest_new, FFT_RADIX_NEXT )
                 IF (nlowest==nlowest_new .AND. MODULO(n(idir),2**(my_ncommensurate-my_icommensurate)).EQ.0) THEN
                    EXIT
                 ELSE
                    ntest(idir)=n(idir)+1
                 ENDIF
              END DO
           END DO
        END IF
     END IF
  ELSE
     ! without a cutoff and HALFSPACE we have to be sure that there is
     ! a negative counterpart to every g vector (-> odd number of grid points)
     IF ( my_grid_span == HALFSPACE .AND. symmetry ) &
          n = n + MOD ( n + 1, 2 )

  END IF

  ! final check if all went fine ...
  IF (my_ncommensurate>0) THEN
     DO my_icommensurate=1,my_ncommensurate
        ftest = ANY( MODULO(n,2**(my_ncommensurate-my_icommensurate)).NE.0 )
        CPPrecondition(.NOT. ftest,cp_failure_level,routineP,error,failure)
     END DO
  ENDIF

  CALL timestop(0.0_dp,handle)

END FUNCTION pw_grid_find_n
!****************************************************************************

!!****f* pw_grids/pw_grid_n_for_fft *
!!
!!   NAME
!!     pw_grid_n_for_fft
!!
!!   FUNCTION
!!     returns the closest number of points >= n, on which you can perform
!!     ffts
!!
!!   NOTES
!!     result<=n
!!
!!   ARGUMENTS
!!     - n: the minimum number of points you want
!!     - odd: if the number has to be odd
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
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
!***************************************************************************

!!****** pw_grids/pw_grid_n_from_cutoff [1.1] *
!!
!!   NAME
!!     pw_grid_n_from_cutoff
!!
!!   FUNCTION
!!     Find the number of points that give at least the requested cutoff
!!
!!   AUTHOR
!!     apsi
!!     Christopher Mundy
!!
!!   MODIFICATION HISTORY
!!     JGH (21-12-2000) : Simplify parameter list, bounds will be global
!!     JGH ( 8-01-2001) : Add check to FFT allowd grids (this now depends
!!                        on the FFT library.
!!                        Should the pw_grid_type have a reference to the FFT
!!                        library ?
!!     JGH (28-02-2001) : Only do conditional check for FFT
!!     JGH (21-05-2002) : Optimise code, remove orthorhombic special case
!!
!!*** *************************************************************************
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

!------------------------------------------------------------------------------

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
!****************************************************************************

!!****** pw_grids/pw_cutoff_from_n [1.1] *
!!
!!   NAME
!!     pw_cutoff_from_n
!!
!!   FUNCTION
!!     Given a n and a box, calculate the corresponding maximum cutoff in a.u.
!!     If npts is odd returns the energy of the next biggest g vector, because
!!     if you take the points with energy smaller than that they are still
!!     representable with that cutoff. Whereas if npts is even the energy of
!!     the last (unpaired) g vector is returned. This implies that you always
!!     have that 2*i+1, and 2*i+2 have the same cutoff.
!!
!!   AUTHOR
!!     fawzi (taken from pw_find_cutoff)
!!
!!*** *************************************************************************

FUNCTION pw_cutoff_from_n ( npts, box, error ) RESULT(cutoff)
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts
    TYPE(cell_type), INTENT(IN)              :: box
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp)                            :: cutoff

    INTEGER                                  :: idir
    REAL(KIND=dp)                            :: cut_att, gdum( 3 )

  cutoff=0._dp
  DO idir=1,3
     gdum ( : ) = twopi * box % h_inv ( idir, : ) &
          * REAL ( ( npts ( idir ) + 1 )/2,KIND=dp)
     cut_att = dotprod_3d ( gdum, gdum )/2._dp
     IF (idir==1) cutoff=cut_att
     cutoff = MIN ( cutoff, cut_att )
  END DO

END FUNCTION pw_cutoff_from_n

!******************************************************************************

!!****f* pw_grids/pw_grid_bounds_from_n *
!!
!!   NAME
!!     pw_grid_bounds_from_n
!!
!!   FUNCTION
!!     returns the bounds that distribute n points evenly around 0
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - npts: the number of points in each direction
!!     - error: variable to control error logging, stopping,...
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     fawzi
!!
!!*** **********************************************************************
FUNCTION pw_grid_bounds_from_n(npts,error) RESULT(bounds)
    INTEGER, DIMENSION(3), INTENT(in)        :: npts
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER, DIMENSION(2, 3)                 :: bounds

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_grid_bounds_from_n', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  bounds(1,:)=-npts/2
  bounds(2,:)=bounds(1,:)+npts-1
END FUNCTION pw_grid_bounds_from_n
!***************************************************************************

!!****** pw_grids/pw_grid_find_bounds [1.1] *
!!
!!   NAME
!!     pw_grid_find_bounds
!!
!!   FUNCTION
!!     Find bounds of grid array
!!
!!   AUTHOR
!!     apsi
!!     Christopher Mundy
!!
!!   MODIFICATION HISTORY
!!     JGH (21-12-2000) : Simplify parameter list, bounds will be global
!!     JGH ( 8-01-2001) : Add check to FFT allowd grids (this now depends
!!                        on the FFT library.
!!                        Should the pw_grid_type have a reference to the FFT
!!                        library ?
!!     JGH (28-02-2001) : Only do conditional check for FFT
!!     JGH (21-05-2002) : Optimise code, remove orthorhombic special case
!!
!!*** *************************************************************************

SUBROUTINE pw_grid_find_bounds ( bounds, hmat, cutoff, fft, error)


    INTEGER, DIMENSION(2, 3), INTENT(OUT)    :: bounds
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(IN)                             :: hmat
    REAL(KIND=dp), INTENT(IN)                :: cutoff
    LOGICAL, INTENT(IN)                      :: fft
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_grid_find_bounds', &
      routineP = moduleN//':'//routineN

    INTEGER, PARAMETER                       :: lwork = 20

    INTEGER                                  :: handle, n(3)

!------------------------------------------------------------------------------

  CALL timeset("pw_grid_find_bounds","I","",handle)

  n=pw_grid_n_from_cutoff(hmat,cutoff,error=error)
  IF (fft) n=pw_grid_n_for_fft(n,error=error)
  bounds = pw_grid_bounds_from_n(n,error=error)

  CALL timestop(0.0_dp,handle)

END SUBROUTINE pw_grid_find_bounds

!****************************************************************************
!!****** pw_grids/pw_find_cutoff_qs [1.1] *
!!
!!   NAME
!!     pw_find_cutoff
!!
!!   FUNCTION
!!     Given a grid and a box, calculate the corresponding cutoff
!!     *** This routine calculates the cutoff in MOMENTUM UNITS! ***
!!
!!   AUTHOR
!!     apsi
!!     Christopher Mundy
!!
!!   MODIFICATION HISTORY
!!     JGH (20-12-2000) : Deleted some strange comments
!!
!!   NOTES
!!     This routine is local. It works independent from the distribution
!!     of PW on processors.
!!     npts is the grid size for the full box.
!!
!!*** *************************************************************************

SUBROUTINE pw_find_cutoff_qs ( bounds, h_inv, cutoff, error )


    INTEGER, DIMENSION(:, :), INTENT(IN)     :: bounds
    REAL(KIND=dp), DIMENSION(:, :)           :: h_inv
    REAL(KIND=dp), INTENT(OUT)               :: cutoff
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER, PARAMETER                       :: lwork = 20

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_find_cutoff_qs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: info
    REAL(KIND=dp)                            :: eig( 3 ), gcut, length, &
                                                vec(3), vec1(3), work( lwork )
    REAL(KIND=dp), DIMENSION(3, 3)           :: gmat, gmat_inv
    LOGICAL                                  :: failure

!------------------------------------------------------------------------------

   failure = .FALSE.
   gmat = MATMUL ( h_inv, TRANSPOSE ( h_inv ) )
   CALL dsyev("V","U",3,gmat,3,eig,work,lwork,info)
   CPPrecondition(info==0,cp_failure_level,routineP,error,failure)
   eig = twopi * SQRT ( eig )
   gmat_inv = inv_3x3 ( gmat )
   vec1(:)=REAL ( bounds(2,:), dp )
   CALL MATVEC_3x3 (vec, gmat_inv (:,:), vec1( : ) )
   length = vec ( 1 ) * eig ( 1 )
   length = CEILING ( length * length * 0.5_dp )
   gcut = length
   length =  vec ( 2 ) * eig ( 2 )
   length = CEILING ( length * length * 0.5_dp )
   gcut = MIN ( gcut, length )
   length =  vec ( 3 ) * eig ( 3 )
   length = CEILING ( length * length * 0.5_dp )
   gcut = MIN ( gcut, length )
   cutoff = gcut

END SUBROUTINE pw_find_cutoff_qs

!****************************************************************************
!!****** pw_grids/pw_find_cutoff_fist [1.1] *
!!
!!   NAME
!!     pw_find_cutoff
!!
!!   FUNCTION
!!     Given a grid and a box, calculate the corresponding cutoff
!!     *** This routine calculates the cutoff in MOMENTUM UNITS! ***
!!
!!   AUTHOR
!!     apsi
!!     Christopher Mundy
!!
!!   MODIFICATION HISTORY
!!     JGH (20-12-2000) : Deleted some strange comments
!!
!!   NOTES
!!     This routine is local. It works independent from the distribution
!!     of PW on processors.
!!     npts is the grid size for the full box.
!!
!!*** *************************************************************************

SUBROUTINE pw_find_cutoff_fist ( npts, box, cutoff, error )


    INTEGER, DIMENSION(:), INTENT(IN)        :: npts
    TYPE(cell_type), INTENT(IN)              :: box
    REAL(KIND=dp), INTENT(OUT)               :: cutoff
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_find_cutoff_fist', &
      routineP = moduleN//':'//routineN

    REAL(KIND=dp)                            :: gcut, gdum( 3 ), length

!------------------------------------------------------------------------------
! compute 2*pi*h_inv^t*g  where g = (bounds[1],0,0)

  gdum ( : ) = twopi *  box % h_inv ( 1, : ) &
       * REAL ( ( npts ( 1 ) - 1 ) / 2 ,KIND=dp )
  length = SQRT ( DOTPROD_3D ( gdum, gdum ) )
  gcut = length

! compute 2*pi*h_inv^t*g  where g = (0,bounds[2],0)
  gdum ( : ) = twopi *  box % h_inv ( 2, : ) &
       * REAL ( ( npts ( 2 ) - 1 ) / 2 ,KIND=dp )
  length = SQRT ( DOTPROD_3D ( gdum, gdum ) )
  gcut = MIN ( gcut, length )

! compute 2*pi*h_inv^t*g  where g = (0,0,bounds[3])
  gdum ( : ) = twopi *  box % h_inv ( 3, : ) &
       * REAL ( ( npts ( 3 ) - 1 ) / 2 ,KIND=dp )
  length = SQRT ( DOTPROD_3D ( gdum, gdum ) )
  gcut = MIN ( gcut, length )
  cutoff = gcut - 1.e-5_dp

END SUBROUTINE pw_find_cutoff_fist

!******************************************************************************

END MODULE pw_grid_info

!******************************************************************************
