!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE coefficient_types
  
  USE kinds, ONLY : dbl
  USE pw_grid_types, ONLY : pw_grid_type
  USE pw_types, ONLY : pw_type, pw_zero, pw_allocate, pw_deallocate, &
       REALDATA1D, COMPLEXDATA1D, REALDATA3D, COMPLEXDATA3D
  USE stop_program, ONLY : stop_prg
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: coeff_type, coeff_allocate, coeff_deallocate, coeff_zero
  PUBLIC :: PLANEWAVES, GAUSSIAN
  PUBLIC :: PW_REALDATA1D, PW_COMPLEXDATA1D, PW_REALDATA3D, PW_COMPLEXDATA3D
  
  !! 'group_id' to grid?
  TYPE coeff_type
     TYPE ( pw_type ) :: pw
!     TYPE ( gaussian_type ) :: gaussian
     
     INTEGER :: in_use
  END TYPE coeff_type
  
  !! Flags for the structure member 'in_use'
  INTEGER, PARAMETER :: PLANEWAVES = 401, GAUSSIAN = 402
  
  INTERFACE coeff_allocate
     MODULE PROCEDURE &
          coeff_allocate_from_coeff, coeff_allocate_from_grid
  END INTERFACE
  
  !! Flags for 'integral'
  PUBLIC :: SQUARE, SQUAREROOT
  INTEGER, PARAMETER :: SQUARE = 391, SQUAREROOT = 392
  
  INTEGER, PARAMETER :: PW_REALDATA1D = REALDATA1D
  INTEGER, PARAMETER :: PW_COMPLEXDATA1D = COMPLEXDATA1D
  INTEGER, PARAMETER :: PW_REALDATA3D = REALDATA3D
  INTEGER, PARAMETER :: PW_COMPLEXDATA3D = COMPLEXDATA3D
  
CONTAINS

!******************************************************************************

SUBROUTINE coeff_allocate_from_coeff ( c_out, c_model, use_basis )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( IN ), TARGET :: c_model
  TYPE ( coeff_type ), INTENT ( OUT ) :: c_out
  INTEGER, INTENT ( IN ), OPTIONAL :: use_basis
  
! Locals
  INTEGER :: use_data
  
!------------------------------------------------------------------------------
  
  IF ( PRESENT ( use_basis ) ) THEN
     IF ( use_basis == PLANEWAVES ) THEN
        c_out % in_use = PLANEWAVES
        CALL pw_allocate ( c_out % pw, c_model % pw % pw_grid, &
             use_data = c_model % pw % in_use )
        
     ELSE IF ( use_basis == PW_REALDATA1D .or. use_basis == PW_COMPLEXDATA1D &
          .or. use_basis == PW_REALDATA3D .or. use_basis == PW_COMPLEXDATA3D &
          ) THEN
        
        c_out % in_use = PLANEWAVES
        
        CALL pw_allocate ( c_out % pw, c_model % pw % pw_grid, &
             use_data = use_basis )
        
     ELSE
        CALL stop_prg ( "coeff_allocate_from_coeff", &
             "no suitable data with use_basis" )
        
     END IF
     
  ELSE
     use_data = c_model % pw % in_use
     IF ( use_data == PW_REALDATA1D .or. use_data == PW_COMPLEXDATA1D .or. &
          use_data == PW_REALDATA3D .or. use_data == PW_COMPLEXDATA3D ) THEN
        
        c_out % in_use = PLANEWAVES
        CALL pw_allocate ( c_out % pw, c_model % pw % pw_grid, &
             use_data = c_model % pw % in_use )
        
     ELSE
        CALL stop_prg ( "coeff_allocate_from_coeff", &
             "no suitable data without use_basis" )
     END IF
  END IF
  
END SUBROUTINE coeff_allocate_from_coeff

!******************************************************************************

SUBROUTINE coeff_allocate_from_grid ( c_out, pw_grid, use_data )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( pw_grid_type ), INTENT ( IN ), OPTIONAL :: pw_grid
  INTEGER, INTENT ( IN ) :: use_data
  TYPE ( coeff_type ), INTENT ( OUT ) :: c_out
  
! Locals
  
!------------------------------------------------------------------------------
  
  IF ( PRESENT ( pw_grid ) ) THEN
     c_out % in_use = PLANEWAVES
     CALL pw_allocate ( c_out % pw, pw_grid, use_data )
  ELSE
     CALL stop_prg ( "coeff_allocate_from_grid", "no suitable grid" )
  END IF
  
END SUBROUTINE coeff_allocate_from_grid

!******************************************************************************

SUBROUTINE coeff_deallocate ( coeff )
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( INOUT ) :: coeff
  
!------------------------------------------------------------------------------
  
  IF ( coeff % in_use == PLANEWAVES ) THEN
     CALL pw_deallocate ( coeff % pw )
  ELSE
     CALL stop_prg ( "coeff_deallocate", "no valid data type" )
  END IF
  
END SUBROUTINE coeff_deallocate

!******************************************************************************

SUBROUTINE coeff_zero ( coeff )
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( INOUT ) :: coeff
  
!------------------------------------------------------------------------------
  
  IF ( coeff % in_use == PLANEWAVES ) THEN
     CALL pw_zero ( coeff % pw )
  ELSE
     CALL stop_prg ( "coeff_zero", "no valid data type" )
  END IF
  
END SUBROUTINE coeff_zero

!******************************************************************************

END MODULE coefficient_types
