!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/potential_types [1.0] *
!!
!!   NAME
!!     potential_types
!!
!!   FUNCTION
!!     definition of the local potential type
!!
!!   AUTHOR
!!     JGH (13-Feb-2001); based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE potential_types
  
  USE coefficient_types,               ONLY: coeff_allocate,&
                                             coeff_deallocate,&
                                             coeff_type
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: potential_type, potential_allocate, potential_deallocate, &
            potential_init, P_NOT_USED, P_IN_USE, P_UP_TO_DATE

  TYPE potential_type

     INTEGER :: npart 
     INTEGER :: nspin 
     INTEGER :: active
     TYPE ( coeff_type ), DIMENSION ( :, : ), POINTER :: pot_part
     INTEGER, DIMENSION ( : ), POINTER :: status

  END TYPE potential_type

  INTEGER, PARAMETER :: P_NOT_USED      = 1
  INTEGER, PARAMETER :: P_IN_USE        = 2
  INTEGER, PARAMETER :: P_UP_TO_DATE    = 3

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** potential_types/potential_allocate [1.0] *
!!
!!   NAME
!!     potential_allocate
!!
!!   FUNCTION
!!     allocate a potential_type variable
!!
!!   AUTHOR
!!     JGH (13-Feb-2001); based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE potential_allocate ( potential, nspin, npart )
  
  
    TYPE(potential_type), INTENT(OUT)        :: potential
    INTEGER, INTENT(IN), OPTIONAL            :: nspin, npart

    INTEGER                                  :: ierr, np, ns

!------------------------------------------------------------------------------

  IF ( PRESENT ( nspin ) ) THEN
    potential % nspin = nspin
  ELSE
    potential % nspin = 1
  END IF

  IF ( PRESENT ( npart ) ) THEN
    potential % npart = npart
  ELSE
    potential % npart = 1
  END IF

  ns = potential % nspin
  np = potential % npart

  ALLOCATE ( potential % status ( np ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "potential_allocate", &
     "potential % status", np )
  potential % status = P_NOT_USED

  ALLOCATE ( potential % pot_part ( ns, np ), STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "potential_allocate", &
     "potential % pot_part", ns*np )
  
!******************************************************************************

END SUBROUTINE potential_allocate
  
!!*****
!******************************************************************************
!!****** potential_types/potential_init [1.0] *
!!
!!   NAME
!!     potential_init
!!
!!   FUNCTION
!!     initialise a realisation of a potential_type
!!
!!   AUTHOR
!!     JGH (13-Feb-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE potential_init ( potential, part, pw_grid, use_data, &
                            c_model, use_basis )


    TYPE(potential_type), INTENT(INOUT)      :: potential
    INTEGER, INTENT(IN), OPTIONAL            :: part
    TYPE(pw_grid_type), OPTIONAL             :: pw_grid
    INTEGER, INTENT(IN), OPTIONAL            :: use_data
    TYPE(coeff_type), OPTIONAL               :: c_model
    INTEGER, INTENT(IN), OPTIONAL            :: use_basis

    INTEGER                                  :: is, np
    TYPE(coeff_type), DIMENSION(:), POINTER  :: potloc

!------------------------------------------------------------------------------

  IF ( PRESENT ( part ) ) THEN
    np = part
  ELSE
    np = 1
  END IF

  IF ( PRESENT ( c_model ) ) THEN

    potential % status ( np ) = P_IN_USE
    potloc => potential % pot_part ( :, np )
    DO is = 1, potential % nspin
      CALL coeff_allocate ( potloc ( np ), c_model, use_basis )
    END DO

  ELSEIF ( PRESENT ( pw_grid ) ) THEN

    IF ( PRESENT ( use_data ) ) THEN
      potential % status ( np ) = P_IN_USE
      potloc => potential % pot_part ( :, np )
      DO is = 1, potential % nspin
        CALL coeff_allocate ( potloc ( is ), pw_grid, use_data )
      END DO
    ELSE
      CALL stop_program ( "potential_init", "use_data is needed with pw_grid" )
    END IF

  END IF

!******************************************************************************

END SUBROUTINE potential_init

!!*****
!******************************************************************************
!!****** potential_types/potential_deallocate [1.0] *
!!
!!   NAME
!!     potential_deallocate
!!
!!   FUNCTION
!!     deallocate a potential_type variable
!!
!!   AUTHOR
!!     JGH (13-Feb-2001); based on earlier version of apsi
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE potential_deallocate ( potential )


    TYPE(potential_type), INTENT(OUT)        :: potential

    INTEGER                                  :: ierr, part, spin

!------------------------------------------------------------------------------

  DO part = 1, potential % npart
    DO spin = 1, potential % nspin
      CALL coeff_deallocate ( potential % pot_part ( spin, part ) )
    END DO
  END DO

  DEALLOCATE ( potential % status, STAT = ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "potential_deallocate", &
    "potential % status" )

  potential % npart = 0
  potential % nspin = 0

!******************************************************************************

END SUBROUTINE potential_deallocate

!!*****
!******************************************************************************

END MODULE potential_types

!******************************************************************************
