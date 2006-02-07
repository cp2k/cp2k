!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/structure_factors [1.0] *
!!
!!   NAME
!!     structure_factors
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE structure_factors

  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: twopi
  USE structure_factor_types,          ONLY: structure_factor_type
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: structure_factor_evaluate, structure_factor_allocate
  PUBLIC :: structure_factor_deallocate

!!*****
!*****************************************************************************

CONTAINS

!*****************************************************************************

SUBROUTINE structure_factor_deallocate ( exp_igr )


    TYPE(structure_factor_type), &
      INTENT(INOUT)                          :: exp_igr

    INTEGER                                  :: isos

!------------------------------------------------------------------------------

  DEALLOCATE ( exp_igr % ex, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ('ewald_polynomial', 'ex')
  DEALLOCATE ( exp_igr % ey, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ('ewald_polynomial', 'ey')
  DEALLOCATE ( exp_igr % ez, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ('ewald_polynomial', 'ez')

  IF ( ASSOCIATED ( exp_igr % centre ) ) THEN
     DEALLOCATE ( exp_igr % centre, STAT = isos )
     IF ( isos /= 0 ) CALL stop_memory ('ewald_polynomial', 'centre' )
  END IF

END SUBROUTINE structure_factor_deallocate

!*****************************************************************************

SUBROUTINE structure_factor_allocate ( bds, nparts, exp_igr, &
                                      allocate_centre)

    INTEGER, DIMENSION(:, :), INTENT(IN)     :: bds
    INTEGER, INTENT(IN)                      :: nparts
    TYPE(structure_factor_type), INTENT(OUT) :: exp_igr
    LOGICAL, INTENT(IN), OPTIONAL            :: allocate_centre

    INTEGER                                  :: isos

!------------------------------------------------------------------------------

  ALLOCATE ( exp_igr % ex ( bds ( 1, 1 ):bds ( 2, 1 ) + 1, nparts ), &
       STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'structure_factor_allocate', 'ex', 0 )
  ALLOCATE ( exp_igr % ey ( bds ( 1, 2 ):bds ( 2, 2 ) + 1, nparts ), &
       STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ('structure_factor_allocate', 'ey', 0 )
  ALLOCATE ( exp_igr % ez ( bds ( 1, 3 ):bds ( 2, 3 ) + 1, nparts ), &
       STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ('structure_factor_allocate', 'ez', 0 )
  NULLIFY ( exp_igr % centre )

  exp_igr % lb(1) = LBOUND ( exp_igr % ex, 1 )
  exp_igr % lb(2) = LBOUND ( exp_igr % ey, 1 )
  exp_igr % lb(3) = LBOUND ( exp_igr % ez, 1 )

  IF ( PRESENT ( allocate_centre ) ) THEN
    IF ( allocate_centre ) THEN
      ALLOCATE ( exp_igr % centre ( 3, nparts ), STAT = isos )
      IF ( isos /= 0 ) &
         CALL stop_memory ( 'structure_factor_allocate', 'centre', 3 * nparts)
    END IF
  END IF

END SUBROUTINE structure_factor_allocate

!*****************************************************************************

SUBROUTINE structure_factor_evaluate ( delta, npts, lb, ex, ey, ez )


    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: delta
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts
    INTEGER, DIMENSION(3), INTENT(IN)        :: lb
    COMPLEX(KIND=dp), DIMENSION(lb(1):), &
      INTENT(out)                            :: ex
    COMPLEX(KIND=dp), DIMENSION(lb(2):), &
      INTENT(out)                            :: ey
    COMPLEX(KIND=dp), DIMENSION(lb(3):), &
      INTENT(out)                            :: ez

    COMPLEX(KIND=dp)                         :: fm, fp
    INTEGER                                  :: j, l0, l1, m0, m1, n0, n1
    REAL(KIND=dp)                            :: vec( 3 )

!------------------------------------------------------------------------------

  l0 = LBOUND ( ex, 1 )
  l1 = UBOUND ( ex, 1 )
  m0 = LBOUND ( ey, 1 )
  m1 = UBOUND ( ey, 1 )
  n0 = LBOUND ( ez, 1 )
  n1 = UBOUND ( ez, 1 )

! delta is in scaled coordinates
  vec ( : ) = twopi * ( delta ( : ) + 0.5_dp  )

  ex ( l0 ) = 1.0_dp
  ey ( m0 ) = 1.0_dp
  ez ( n0 ) = 1.0_dp
  ex ( l1 ) = 1.0_dp
  ey ( m1 ) = 1.0_dp
  ez ( n1 ) = 1.0_dp

  fp = CMPLX ( COS ( vec ( 1 ) ), -SIN ( vec ( 1 ) ),KIND=dp)
  fm = CONJG ( fp )
  DO j = 1, -l0
     ex (  j + l0 ) = ex (  j + l0 - 1 ) * fp
     ex ( -j + l1 ) = ex ( -j + l1 + 1 ) * fm
  END DO

  fp = CMPLX ( COS ( vec ( 2 ) ), -SIN ( vec ( 2 ) ),KIND=dp)
  fm = CONJG ( fp )
  DO j = 1, -m0
     ey (  j + m0 ) = ey (  j + m0 - 1 ) * fp
     ey ( -j + m1 ) = ey ( -j + m1 + 1 ) * fm
  END DO

  fp = CMPLX ( COS ( vec ( 3 ) ), -SIN ( vec ( 3 ) ),KIND=dp)
  fm = CONJG ( fp )
  DO j = 1, -n0
     ez (  j + n0 ) = ez (  j + n0 - 1 ) * fp
     ez ( -j + n1 ) = ez ( -j + n1 + 1 ) * fm
  END DO

END SUBROUTINE structure_factor_evaluate

!*****************************************************************************

END MODULE structure_factors
