!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/radial_operators [1.0] *
!!
!!   NAME
!!     radial_operators
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH 17-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE radial_operators

  USE kinds,                           ONLY: dp
  USE radial_grids,                    ONLY: radial_grid_type
  USE radial_util,                     ONLY: radial_diff,&
                                             radial_integrate
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: rad_kinetic, rad_potential

!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
!!****** radial_operators/rad_kinetic [1.0] *
!!
!!   NAME
!!     Kinetic Energy operator
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH 17-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!-----------------------------------------------------------------------------!

  SUBROUTINE rad_kinetic ( p1, p2, l, rg, result )

! Kinetic energy integral < f1 T f2 > where p1 = r*f1(r) and p2 = r*f2(r)
    REAL(KIND=dp), DIMENSION(0:), INTENT(IN)     :: p1, p2
    INTEGER, INTENT(IN)                      :: l
    TYPE(radial_grid_type), INTENT(IN)       :: rg
    REAL(KIND=dp), INTENT(OUT)                   :: result

    INTEGER                                  :: i, ierr, n
    REAL(KIND=dp)                                :: int
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)     :: dp1, dp2, f

!------------------------------------------------------------------------------

    n = SIZE ( p1 ) - 2
    ALLOCATE ( dp1 ( 0:n+1 ), STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "rad_kinetic", "dp1", n+2 )
    ALLOCATE ( dp2 ( 0:n+1 ), STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "rad_kinetic", "dp2", n+2 )
    ALLOCATE ( f ( 0:n+1 ), STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "rad_kinetic", "f", n+2 )

    CALL radial_diff ( n+2, p1, dp1, 1, rg, ierr )
    IF ( ierr /= 0 ) call stop_program ( "rad_kinetic", "radial_diff: dp1" )
    CALL radial_diff ( n+2, p2, dp2, 1, rg, ierr )
    IF ( ierr /= 0 ) call stop_program ( "rad_kinetic", "radial_diff: dp2" )

    f = dp1 * dp2
    CALL radial_integrate ( n+2, f, rg, int, ierr )
    IF ( ierr /= 0 ) call stop_program ( "rad_kinetic", "radial_integrate: f" )
    result = 0.5_dp * int
    IF ( l > 0 ) THEN
      f ( 0 ) = 0.0_dp
      DO i = 1, n+1
        f ( i ) = p1 ( i ) * p2 ( i ) / rg % r ( i ) ** 2
      END DO
      CALL radial_integrate ( n+2, f, rg, int, ierr )
      IF ( ierr /= 0 ) call stop_program ( "rad_kinetic", "radial_integrate: l" )
      result = result + 0.5_dp * REAL ( l*(l+1),KIND=dp) * int
    END IF

    DEALLOCATE ( dp1, STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "rad_kinetic", "dp1" )
    DEALLOCATE ( dp2, STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "rad_kinetic", "dp2" )
    DEALLOCATE ( f, STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "rad_kinetic", "f" )

  END SUBROUTINE rad_kinetic

!!*****
!******************************************************************************
!!****** radial_operators/rad_potential [1.0] *
!!
!!   NAME
!!     Potential operator
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH 17-NOV-2000
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!-----------------------------------------------------------------------------!

  SUBROUTINE rad_potential ( p1, p2, v, rg, result )

! Calculates the integral < f1 V f2 >
! It is assumed that the s=0 term has been included correctly into V(s=0)
! This is only a problem if V(s->0) = s^2 or worse
    REAL(KIND=dp), DIMENSION(0:), INTENT(IN)     :: p1, p2, v
    TYPE(radial_grid_type), INTENT(IN)       :: rg
    REAL(KIND=dp), INTENT(OUT)                   :: result

    INTEGER                                  :: ierr, n
    REAL(KIND=dp)                                :: int
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)     :: f

!------------------------------------------------------------------------------

    n = SIZE ( p1 ) - 2
    ALLOCATE ( f ( 0:n+1 ), STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "rad_potential", "f", n+2 )

    f = p1 * p2 * V
    CALL radial_integrate ( n+2, f, rg, int, ierr )
    IF ( ierr /= 0 ) call stop_program ( "rad_potential", "radial_integrate" )
    result = int

    DEALLOCATE ( f, STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "rad_potential", "f" )

  END SUBROUTINE rad_potential

!!*****
!******************************************************************************

END MODULE radial_operators

!******************************************************************************
