!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/initialize_pimd_types [1.0] *
!!
!!   NAME
!!     initialize_pimd_types
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

MODULE initialize_pimd_types
  
!  USE termination, ONLY : stop_memory
!  USE kinds, ONLY : dbl
!  USE md, ONLY : simulation_parameters_type
!  USE transformations, ONLY : r2u, transformation_eigenvalues
!  IMPLICIT NONE
  
!  PRIVATE
!  PUBLIC :: initialize_pimd
  
!!*****
!******************************************************************************

!CONTAINS

!******************************************************************************
!!****** md/initialize_pimd [1.0] *
!!
!!   NAME
!!     initialize_pimd
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     HAF
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!
!!   SOURCE
!******************************************************************************

!SUBROUTINE initialize_pimd ( struc, simpar )
  
!  IMPLICIT NONE
  
! Arguments
!  TYPE ( structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: struc
!  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar

! Locals:
!  INTEGER :: ios,i,j,k,beads,natoms
!  REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: lambdas
!  REAL ( dbl ) :: fac, edt, phi

!  beads = SIZE(struc)
!  natoms = SIZE(struc(1) % part)

!  ALLOCATE(lambdas(beads), STAT = ios)
!  IF (ios /= 0) CALL stop_memory ( 'initialize_pimd', 'lambdas', beads )

!  CALL transformation_eigenvalues( lambdas, &
!    simpar % pimd_params % transformation)

!  edt = simpar % dt / simpar % pimd_params % ndelta
!  DO i= 1, beads
!     ALLOCATE(struc(i)%pimd_env, STAT = ios)
!     IF (ios /= 0) CALL stop_memory ( 'initialize_pimd', 'pimd_env', i )
!     struc(i) % pimd_env % lbd  = lambdas(i)
!     struc(i) % pimd_env % slbd  = SQRT(lambdas(i))
!
!+++ here is the place for facstage and wmass
!     IF (lambdas(i) == 0.0_dbl) THEN
!       struc(i) % pimd_env % meff = 1.0
!       struc(i) % pimd_env % ismeff = 1.0/SQRT(struc(i)%pimd_env%meff)
!       struc(i) % pimd_env % harmonic_prop (1,1)=1.0
!       struc(i) % pimd_env % harmonic_prop(1,2)=edt*struc(i)%pimd_env%ismeff
!       struc(i) % pimd_env % harmonic_prop(2,1)=0.0
!       struc(i) % pimd_env % harmonic_prop(2,2)=1.0
!     ELSE
!       struc(i) % pimd_env % meff = lambdas(i)
!       struc(i) % pimd_env % ismeff = 1.0/SQRT(struc(i)%pimd_env%meff)
!       phi = simpar % pimd_params % wp * edt * &
!          struc(i) % pimd_env % slbd * struc(i) % pimd_env % ismeff 
!       fac = struc(i) % pimd_env % slbd * simpar % pimd_params % wp
!       struc(i) % pimd_env % harmonic_prop(1,1) = cos( phi )
!       struc(i) % pimd_env % harmonic_prop(1,2) = sin( phi ) / fac
!       struc(i) % pimd_env % harmonic_prop(2,1) = - sin( phi ) * fac
!       struc(i) % pimd_env % harmonic_prop(2,2) = cos( phi )
!     END IF
!  END DO

!  DEALLOCATE( lambdas , STAT = ios)
!  IF (ios /= 0) CALL stop_memory ( 'initialize_pimd', 'lambdas' )

!  DO j= 1, beads
!     ALLOCATE(struc(j)%pimd_env%u(3,natoms), STAT=ios)
!     IF (ios /= 0) CALL stop_memory ( 'initialize_pimd', 'u', 3*natoms)
!  END DO

!  CALL r2u( struc, simpar)

!END SUBROUTINE initialize_pimd

!!*****
!******************************************************************************

END MODULE initialize_pimd_types

!******************************************************************************
