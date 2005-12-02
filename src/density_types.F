!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/density_types [1.0] *
!!
!!   NAME
!!     density_types
!!
!!   FUNCTION
!!     definition of the electronic density type
!!
!!   AUTHOR
!!     JGH (12-Feb-2001); based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE density_types
  
  USE coefficient_types,               ONLY: coeff_allocate,&
                                             coeff_deallocate,&
                                             coeff_type
  USE kinds,                           ONLY: dp
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: density_type, density_allocate, density_deallocate, density_init
  PUBLIC :: D_NOT_USED, D_IN_USE, D_UP_TO_DATE

  TYPE density_type

     INTEGER :: npart 
     INTEGER :: nspin 
     INTEGER :: active
     TYPE ( coeff_type ), DIMENSION ( :, : ), POINTER :: dpart
     INTEGER, DIMENSION ( : ), POINTER :: status
     REAL (KIND=dp), DIMENSION ( :, : ), POINTER :: charge

  END TYPE density_type

  INTEGER, PARAMETER :: D_NOT_USED      = 1
  INTEGER, PARAMETER :: D_IN_USE        = 2
  INTEGER, PARAMETER :: D_UP_TO_DATE    = 3

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** density_types/density_allocate [1.0] *
!!
!!   NAME
!!     density_allocate
!!
!!   FUNCTION
!!     allocate a density_type variable
!!
!!   AUTHOR
!!     JGH (12-Feb-2001); based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE density_allocate ( rho, nspin, npart )
  
  
    TYPE(density_type), INTENT(OUT)          :: rho
    INTEGER, INTENT(IN), OPTIONAL            :: nspin, npart

    INTEGER                                  :: ierr, np, ns

!------------------------------------------------------------------------------

  IF ( PRESENT ( nspin ) ) THEN
    rho % nspin = nspin
  ELSE
    rho % nspin = 1
  END IF

  IF ( PRESENT ( npart ) ) THEN
    rho % npart = npart
  ELSE
    rho % npart = 1
  END IF

  rho % active = 1

  ns = rho % nspin
  np = rho % npart

  ALLOCATE ( rho % charge ( ns, np ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "density_allocate", &
     "rho % charge", ns*np )
  rho % charge = 0.0_dp

  ALLOCATE ( rho % status ( np ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "density_allocate", &
     "rho % status", np )
  rho % status = D_NOT_USED

  ALLOCATE ( rho % dpart ( ns, np ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "density_allocate", &
     "rho % dpart", ns*np )
  
!******************************************************************************

END SUBROUTINE density_allocate
  
!!*****
!******************************************************************************
!!****** density_types/density_init [1.0] *
!!
!!   NAME
!!     density_init
!!
!!   FUNCTION
!!     initialise a realisation of a density_type
!!
!!   AUTHOR
!!     JGH (13-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE density_init ( rho, part, pw_grid, use_data, c_model, use_basis )


    TYPE(density_type), INTENT(INOUT)        :: rho
    INTEGER, INTENT(IN), OPTIONAL            :: part
    TYPE(pw_grid_type),POINTER, OPTIONAL     :: pw_grid
    INTEGER, INTENT(IN), OPTIONAL            :: use_data
    TYPE(coeff_type), OPTIONAL               :: c_model
    INTEGER, INTENT(IN), OPTIONAL            :: use_basis

    INTEGER                                  :: is, np
    TYPE(coeff_type), DIMENSION(:), POINTER  :: dloc

!------------------------------------------------------------------------------

  IF ( PRESENT ( part ) ) THEN
    np = part
  ELSE
    np = 1
  END IF

  IF ( PRESENT ( c_model ) ) THEN

    rho % status ( np ) = D_IN_USE
    dloc => rho % dpart ( :, np )
    DO is = 1, rho % nspin
      CALL coeff_allocate ( dloc ( np ), c_model, use_basis )
    END DO

  ELSEIF ( PRESENT ( pw_grid ) ) THEN

    IF ( PRESENT ( use_data ) ) THEN
      rho % status ( np ) = D_IN_USE
      dloc => rho % dpart ( :, np )
      DO is = 1, rho % nspin
        CALL coeff_allocate ( dloc ( is ), pw_grid, use_data )
      END DO
    ELSE
      CALL stop_program ( "density_init", "use_data is needed with pw_grid" )
    END IF

  END IF

!******************************************************************************

END SUBROUTINE density_init

!!*****
!******************************************************************************
!!****** density_types/density_deallocate [1.0] *
!!
!!   NAME
!!     density_deallocate
!!
!!   FUNCTION
!!     deallocate a density_type variable
!!
!!   AUTHOR
!!     JGH (12-Feb-2001); based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE density_deallocate ( rho )


    TYPE(density_type), INTENT(OUT)          :: rho

    INTEGER                                  :: ierr, part, spin

!------------------------------------------------------------------------------

  DO part = 1, rho % npart
    DO spin = 1, rho % nspin
      CALL coeff_deallocate ( rho % dpart ( spin, part ) )
    END DO
  END DO

  DEALLOCATE ( rho % charge, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "density_deallocate", "rho % charge" )

  DEALLOCATE ( rho % status, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "density_deallocate", "rho % status" )

  rho % npart = 0
  rho % nspin = 0

!******************************************************************************

END SUBROUTINE density_deallocate

!!*****
!******************************************************************************

END MODULE density_types

!******************************************************************************
