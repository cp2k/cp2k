!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/transformations [1.0] *
!!
!!   NAME
!!     transformations
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     HAF
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE transformations
  
!  USE kinds, ONLY: dp
!  USE termination, ONLY : stop_program, stop_memory
!
!  USE structure_types, ONLY : structure_type
!  USE md, ONLY : simulation_parameters_type
!  
!  IMPLICIT NONE
!  
!  PRIVATE
!  PUBLIC :: transformation_eigenvalues, transform, back_transform
!  PUBLIC :: force_transform
!  PUBLIC :: r2u, u2r, fr2fu
!  PUBLIC :: v2ud, ud2v
!
!!!*****
!!******************************************************************************
!
!CONTAINS
!
!!******************************************************************************
!!!****** pimd/transformation_eigenvalues [1.0] *
!!!
!!!   NAME
!!!     transformation_eigenvalues
!!!
!!!   FUNCTION
!!!
!!!   AUTHOR
!!!     HAF
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   NOTES
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUTINE transformation_eigenvalues( lamb, tran )
!  
!  IMPLICIT NONE
!  
!! Arguments
!  REAL (KIND=dp), DIMENSION ( : ), INTENT(OUT) :: lamb
!  CHARACTER ( LEN = 20 ), INTENT(IN) :: tran
!  
!! Locals
!
!  INTEGER :: i, beads
!  
!!------------------------------------------------------------------------------
!
!  beads = SIZE(lamb)
!  SELECT CASE ( tran )
!  CASE DEFAULT
!    CALL stop_program('transformation_eigenvalues', &
!       'That transformation not supported atm.')
!  CASE ( 'STAGE' )
!    lamb(1) = 0.0_dp
!    DO i=2, beads
!       lamb(i) = i*(1.0_dp/(i-1))
!    END DO
!  END SELECT
!  
!END SUBROUTINE transformation_eigenvalues
!
!!!*****
!!******************************************************************************
!
!!******************************************************************************
!!!****** pimd/transform [1.0] *
!!!
!!!   NAME
!!!     transform
!!!
!!!   FUNCTION
!!!
!!!   AUTHOR
!!!     HAF
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   NOTES
!!!     transform(tran,n,p,r,u) transforms a 2D array ofKIND=dpvia 
!!!     the transformation 'tran'. 'n' is the number of sets to transform
!!!     'p' is the number of elements per set.
!!!     'r' is the input (real coordinates) of dimension (1:n,1:p+1)
!!!     'u' is the output of the transformed coordinates, if given.
!!!     If 'u' is ommitted, the output will be in 'r'.
!!!
!!!     Periodic boundary conditions in the sets have to be put in
!!!     by duplicating r(:,p+1)=r(:,1). The idea is to be able to
!!!     use these routines for partial (staging) transformations also.
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUTINE transform( tran, n, p, r, u )
!  
!  IMPLICIT NONE
!  
!! Arguments
!  CHARACTER ( LEN = 20 ), INTENT(IN) :: tran    ! which transform
!  INTEGER, INTENT(IN) :: n                      ! how many coordinate sets
!  INTEGER, INTENT(IN) :: p                      ! how many coordinates per set
!  REAL (KIND=dp), DIMENSION ( : , : ), INTENT(INOUT) :: r
!  REAL (KIND=dp), DIMENSION ( : , : ), INTENT(INOUT), OPTIONAL :: u
!  
!! Locals
!
!  INTEGER :: i
!  
!!------------------------------------------------------------------------------
!
!  SELECT CASE ( tran )
!  CASE DEFAULT
!    CALL stop_program('transform', 'That transformation not supported atm.')
!  CASE ( 'STAGE' )
!     IF ( PRESENT ( u ) ) THEN
!       u(1:n,1) = r(1:n,1)
!       DO i = 2, p
!         u(1:n,i) = r(1:n,i) - ((i-1) * r(1:n,i+1) + r(1:n,1))/i
!       END DO
!     ELSE
!       DO i = 2, p
!         r(1:n,i) = r(1:n,i) - ((i-1) * r(1:n,i+1) + r(1:n,1))/i
!       END DO
!     END IF
!  END SELECT
!  
!END SUBROUTINE transform
!
!!!*****
!!******************************************************************************
!
!!******************************************************************************
!!!****** pimd/back_transform [1.0] *
!!!
!!!   NAME
!!!     back_transform
!!!
!!!   FUNCTION
!!!
!!!   AUTHOR
!!!     HAF
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   NOTES
!!!     back_transform(tran,n,p,u,r) un-does  transform(tran,n,p,r,u)
!!!     'u' is the input ('normal' coordinates) of dimension (1:n,1:p+1)
!!!     'r' is the output of the real coordinates, if given.
!!!     If 'r' is ommitted, the output will be in 'u'.
!!!
!!!     Periodic boundary conditions in the sets have to be put in
!!!     by duplicating u(:,p+1)=u(:,1). The idea is to be able to
!!!     use these routines for partial (staging) transformations also.
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUTINE back_transform( tran, n, p, u, r )
!  
!  IMPLICIT NONE
!  
!! Arguments
!  CHARACTER ( LEN = 20 ), INTENT(IN) :: tran    ! which transform
!  INTEGER, INTENT(IN) :: n                      ! how many coordinate sets
!  INTEGER, INTENT(IN) :: p                      ! how many coordinates per set
!  REAL (KIND=dp), DIMENSION ( : , : ), INTENT(INOUT) :: u
!  REAL (KIND=dp), DIMENSION ( : , : ), INTENT(INOUT), OPTIONAL :: r
!  
!! Locals
!
!  INTEGER :: i
!  
!!------------------------------------------------------------------------------
!
!  SELECT CASE ( tran )
!  CASE DEFAULT
!    CALL stop_program('back_transform', 'That transformation not supported atm.')
!  CASE ( 'STAGE' )
!     IF ( PRESENT ( r ) ) THEN
!       r(1:n,1) = u(1:n,1)
!       DO i = p , 2, -1
!          r(1:n,i) = u(1:n,i) + ((i-1) * r(1:n,i+1) + r(1:n,1))/i
!       END DO
!     ELSE
!       DO i = p , 2, -1
!          u(1:n,i) = u(1:n,i) + ((i-1) * u(1:n,i+1) + u(1:n,1))/i
!       END DO
!     END IF
!  END SELECT
!  
!END SUBROUTINE back_transform
!
!!!*****
!!******************************************************************************
!
!!******************************************************************************
!!!****** pimd/force_transform [1.0] *
!!!
!!!   NAME
!!!     force_transform
!!!
!!!   FUNCTION
!!!
!!!   AUTHOR
!!!     HAF
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   NOTES
!!!     force_transform(tran,n,p,fr,fu) transforms the forces calculated
!!!     in the real coordinate frame into the frame of the of 'normal'
!!!     coordinates according to the transformation 'tran'.
!!!     'n' is the number of set to transform.
!!!     'p' is the number of elements per set.
!!!     'fr' is the input ('normal' coordinates forces) of dimension (1:n,1:p)
!!!     'fu' is the output of the transformed forces, if given.
!!!     If 'fu' is ommitted, the output will be in 'fr'.
!!!
!!!     Periodic boundary conditions in the sets are ASSUMED! You have to
!!!     post-process the "fu" in the following way to get non-periodic results
!!!     fu(:,1) = fu(:,1) + (p-1)*(fu(:,0)-fu(:,p))/p
!!!     where fu(:,0) is fu(:,p) of the previous segment
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUTINE force_transform( tran, n, p, fr, fu )
!  
!  IMPLICIT NONE
!  
!! Arguments
!  CHARACTER ( LEN = 20 ), INTENT(IN) :: tran    ! which transform
!  INTEGER, INTENT(IN) :: n                      ! how many coordinate sets
!  INTEGER, INTENT(IN) :: p                      ! how many coordinates per set
!  REAL (KIND=dp), DIMENSION ( : , : ), INTENT(INOUT) :: fr
!  REAL (KIND=dp), DIMENSION ( : , : ), INTENT(INOUT), OPTIONAL :: fu
!  
!! Locals
!
!  INTEGER :: i, j
!  REAL(KIND=dp) :: w
!  
!!------------------------------------------------------------------------------
!
!  SELECT CASE ( tran )
!  CASE DEFAULT
!    CALL stop_program('force_transform', 'That transformation not supported atm.')
!  CASE ( 'STAGE' )
!
!     IF ( PRESENT ( fu ) ) THEN
!        fu(1:n,1) = 0.0_dp
!        DO i = 1, p
!           fu(1:n,1) = fu(1:n,1) + fr(1:n,i) 
!        END DO
!
!        DO i = 2, p
!           w = (i - 2.0_dp) / (i - 1.0_dp)
!           fu(1:n,i) = w*fu(1:n,i-1) + fr(1:n,i)
!        END DO
!     ELSE
!        DO i = 2, p
!           fr(1:n,1) = fr(1:n,1) + fr(1:n,i)
!           w = (i - 2.0_dp) / (i - 1.0_dp)
!           fr(1:n,i) = fr(1:n,i) + w*fr(1:n,i-1)
!        END DO
!     END IF
!  END SELECT
!
!END SUBROUTINE force_transform
!
!!!*****
!!******************************************************************************
!
!!******************************************************************************
!!!****** pimd/r2u [1.0] *
!!!
!!!   NAME
!!!     r2u
!!!
!!!   FUNCTION
!!!
!!!   AUTHOR
!!!     HAF
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   NOTES
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUTINE r2u( struc, simpar )
!
!  TYPE ( structure_type ), DIMENSION ( : ), INTENT(INOUT) :: struc
!  TYPE ( simulation_parameters_type ), INTENT(IN) :: simpar
!
!  REAL(KIND=dp), ALLOCATABLE, DIMENSION ( : , : ) :: rtmp
!  INTEGER :: beads, natoms, ios
!  INTEGER :: i,j,k
!
!  beads = SIZE(struc)
!  natoms = SIZE(struc(1) % part)
!
!  ALLOCATE(rtmp(3*natoms,beads+1), STAT=ios)
!  IF (ios /= 0) CALL stop_memory ( 'r2u', 'rtmp', 3*natoms*(beads+1))
!
!  DO j=1, beads
!     DO k=1, natoms
!        DO i= 1, 3
!           rtmp(3*k+i-3,j)= struc(j)%part(k)%r(i)
!        END DO
!     END DO
!  END DO
!  rtmp(:,beads+1)=rtmp(:,1)
!
!  CALL transform( simpar%pimd_params%transformation, 3*natoms, beads, rtmp)
!
!  DO j= 1, beads
!     DO k=1, natoms
!        DO i= 1, 3
!           struc(j)%pimd_env%u(i,k)=rtmp(3*k+i-3,j)
!        END DO
!     END DO
!  END DO
!
!  DEALLOCATE( rtmp , STAT = ios)
!  IF (ios /= 0) CALL stop_memory ( 'r2u', 'rtmp' )
!
!END SUBROUTINE r2u
!
!!!*****
!!******************************************************************************
!
!!******************************************************************************
!!!****** pimd/u2r [1.0] *
!!!
!!!   NAME
!!!     u2r
!!!
!!!   FUNCTION
!!!
!!!   AUTHOR
!!!     HAF
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   NOTES
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUTINE u2r( struc, simpar )
!
!  TYPE ( structure_type ), DIMENSION ( : ), INTENT(INOUT) :: struc
!  TYPE ( simulation_parameters_type ), INTENT(IN) :: simpar
!
!  REAL(KIND=dp), ALLOCATABLE, DIMENSION ( : , : ) :: utmp
!  INTEGER :: beads, natoms, ios
!  INTEGER :: i,j,k
!
!  beads = SIZE(struc)
!  natoms = SIZE(struc(1) % part)
!
!  ALLOCATE(utmp(3*natoms,beads+1), STAT=ios)
!  IF (ios /= 0) CALL stop_memory ( 'u2r', 'utmp', 3*natoms*(beads+1))
!
!  DO j= 1, beads
!     DO k=1, natoms
!        DO i= 1, 3
!           utmp(3*k+i-3,j) = struc(j)%pimd_env%u(i,k)
!        END DO
!     END DO
!  END DO
!  utmp(:,beads+1)=utmp(:,1)
!
!  CALL back_transform( simpar%pimd_params%transformation,3*natoms, beads, utmp)
!
!  DO j=1, beads
!     DO k=1, natoms
!        DO i= 1, 3
!           struc(j)%part(k)%r(i) = utmp(3*k+i-3,j)
!        END DO
!     END DO
!  END DO
!
!  DEALLOCATE( utmp , STAT = ios)
!  IF (ios /= 0) CALL stop_memory ( 'u2r', 'utmp' )
!
!END SUBROUTINE u2r
!
!!!*****
!!******************************************************************************
!
!!******************************************************************************
!!!****** pimd/fr2fu [1.0] *
!!!
!!!   NAME
!!!     fr2fu
!!!
!!!   FUNCTION
!!!
!!!   AUTHOR
!!!     HAF
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   NOTES
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUTINE fr2fu( struc, simpar )
!
!  TYPE ( structure_type ), DIMENSION ( : ), INTENT(INOUT) :: struc
!  TYPE ( simulation_parameters_type ), INTENT(IN) :: simpar
!
!  REAL(KIND=dp), ALLOCATABLE, DIMENSION ( : , : ) :: ftmp
!  INTEGER :: beads, natoms, ios
!  INTEGER :: i,j,k
!
!  beads = SIZE(struc)
!  natoms = SIZE(struc(1) % part)
!
!  ALLOCATE(ftmp(3*natoms,beads+1), STAT=ios)
!  IF (ios /= 0) CALL stop_memory ( 'fr2fu', 'ftmp', 3*natoms*(beads+1))
!
!  DO j=1, beads
!     DO k=1, natoms
!        DO i= 1, 3
!           ftmp(3*k+i-3,j) = struc(j)%part(k)%f(i)
!        END DO
!     END DO
!  END DO
! 
!  CALL force_transform( simpar%pimd_params%transformation,3*natoms,beads, ftmp)
!
!  DO j=1, beads
!     DO k=1, natoms
!        DO i= 1, 3
!           struc(j)%part(k)%f(i) = ftmp(3*k+i-3,j) 
!        END DO
!     END DO
!  END DO
! 
!  DEALLOCATE( ftmp , STAT = ios)
!  IF (ios /= 0) CALL stop_memory ( 'fr2fu', 'ftmp' )
!
!END SUBROUTINE fr2fu
!
!!!*****
!!******************************************************************************
!
!!******************************************************************************
!!!****** pimd/v2ud [1.0] *
!!!
!!!   NAME
!!!     v2ud
!!!
!!!   FUNCTION
!!!
!!!   AUTHOR
!!!     HAF
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   NOTES
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUTINE v2ud( struc, simpar )
!
!  TYPE ( structure_type ), DIMENSION ( : ), INTENT(INOUT) :: struc
!  TYPE ( simulation_parameters_type ), INTENT(IN) :: simpar
!
!  REAL(KIND=dp), ALLOCATABLE, DIMENSION ( : , : ) :: rtmp
!  REAL(KIND=dp) :: factor
!  INTEGER :: beads, natoms, ios
!  INTEGER :: i,j,k
!
!  beads = SIZE(struc)
!  natoms = SIZE(struc(1) % part)
!
!  ALLOCATE(rtmp(3*natoms,beads+1), STAT=ios)
!  IF (ios /= 0) CALL stop_memory ( 'v2ud', 'rtmp', 3*natoms*(beads+1))
!
!  DO j=1, beads
!     DO k=1, natoms
!        DO i= 1, 3
!           rtmp(3*k+i-3,j)= struc(j)%part(k)%v(i)
!        END DO
!     END DO
!  END DO
!  rtmp(:,beads+1)=rtmp(:,1)
!
!  CALL transform( simpar%pimd_params%transformation, 3*natoms, beads, rtmp)
!
!  DO j= 1, beads
!     factor = 1.0_dp/struc(j)%pimd_env%ismeff
!     DO k=1, natoms
!        DO i= 1, 3
!           struc(j)%part(k)%v(i)=factor*rtmp(3*k+i-3,j)
!        END DO
!     END DO
!  END DO
!
!  DEALLOCATE( rtmp , STAT = ios)
!  IF (ios /= 0) CALL stop_memory ( 'v2ud', 'rtmp' )
!
!END SUBROUTINE v2ud
!
!!!*****
!!******************************************************************************
!
!!******************************************************************************
!!!****** pimd/ud2v [1.0] *
!!!
!!!   NAME
!!!     ud2v
!!!
!!!   FUNCTION
!!!
!!!   AUTHOR
!!!     HAF
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   NOTES
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUTINE ud2v( struc, simpar )
!
!  TYPE ( structure_type ), DIMENSION ( : ), INTENT(INOUT) :: struc
!  TYPE ( simulation_parameters_type ), INTENT(IN) :: simpar
!
!  REAL(KIND=dp), ALLOCATABLE, DIMENSION ( : , : ) :: utmp
!  INTEGER :: beads, natoms, ios
!  INTEGER :: i,j,k
!
!  beads = SIZE(struc)
!  natoms = SIZE(struc(1) % part)
!
!  ALLOCATE(utmp(3*natoms,beads+1), STAT=ios)
!  IF (ios /= 0) CALL stop_memory ( 'ud2v', 'utmp', 3*natoms*(beads+1))
!
!  DO j= 1, beads
!     DO k=1, natoms
!        DO i= 1, 3
!           utmp(3*k+i-3,j) = struc(j)%part(k)%v(i)*struc(j)%pimd_env%ismeff
!        END DO
!     END DO
!  END DO
!  utmp(:,beads+1)=utmp(:,1)
!
!  CALL back_transform( simpar%pimd_params%transformation,3*natoms, beads, utmp)
!
!  DO j=1, beads
!     DO k=1, natoms
!        DO i= 1, 3
!           struc(j)%part(k)%v(i) = utmp(3*k+i-3,j)
!        END DO
!     END DO
!  END DO
!
!  DEALLOCATE( utmp , STAT = ios)
!  IF (ios /= 0) CALL stop_memory ( 'ud2v', 'utmp' )
!
!END SUBROUTINE ud2v
!
!!!*****
!!******************************************************************************
!
END MODULE transformations

!******************************************************************************
