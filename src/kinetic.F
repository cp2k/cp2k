!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/kinetic [1.0] *
!!
!!   NAME
!!     kinetic
!!
!!   FUNCTION
!!     Wrapper to call the correponding kinetic energy/Hamiltonian operator/
!!        stress routines
!!
!!   AUTHOR
!!     Christopher J. Mundy
!!
!!   MODIFICATION HISTORY
!!     APS: Added the possibility of stress
!!
!!   SOURCE
!******************************************************************************

MODULE kinetic
  
  USE coefficient_types,               ONLY: PLANEWAVES,&
                                             coeff_type
  USE kinds,                           ONLY: dbl
  USE pws,                             ONLY: pw_kinetic_evaluate
  USE termination,                     ONLY: stop_program

  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: kinetic_evaluate
  
CONTAINS

!!*****
!******************************************************************************

SUBROUTINE kinetic_evaluate ( psi_iks, x_kpt, e_kinetic, kinetic_oper_psi, &
     stress_kinetic )
  
  
    TYPE(coeff_type), INTENT(IN)             :: psi_iks
    REAL(dbl), INTENT(IN)                    :: x_kpt( 3 )
    REAL(dbl), INTENT(INOUT), OPTIONAL       :: e_kinetic
    TYPE(coeff_type), INTENT(INOUT), &
      OPTIONAL                               :: kinetic_oper_psi
    REAL(dbl), INTENT(INOUT), OPTIONAL       :: stress_kinetic( 3, 3 )

!------------------------------------------------------------------------------

  IF ( psi_iks % in_use == PLANEWAVES ) THEN
     CALL pw_kinetic_evaluate ( psi_iks, x_kpt, e_kinetic, kinetic_oper_psi, &
          stress_kinetic )
     
!  ELSE IF ( psi_iks % in_use == GAUSSIAN ) THEN
!     stop "kinetic, gaussians not yet implemented"
     
  ELSE
     CALL stop_program ( "kinetic", "other grids not yet implemented" )
  END IF
  
END SUBROUTINE kinetic_evaluate

!******************************************************************************

END MODULE kinetic
