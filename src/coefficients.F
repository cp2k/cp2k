!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/coefficients [1.0] *
!!
!!   NAME
!!     coefficients
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

MODULE coefficients

  USE coefficient_types, ONLY : coeff_type, PLANEWAVES
  USE kinds, ONLY : dbl
  USE pws, ONLY : pw_integral_a, pw_integral_ab

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: integral

  INTERFACE integral
     MODULE PROCEDURE integral_a, integral_ab
  END INTERFACE

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

FUNCTION integral_a ( coeff, flag ) RESULT ( integral_value )

  TYPE ( coeff_type ), INTENT ( IN ) :: coeff
  INTEGER, INTENT ( IN ), OPTIONAL :: flag

  REAL ( dbl ) :: integral_value

!------------------------------------------------------------------------------

  IF ( coeff % in_use == PLANEWAVES ) THEN
     integral_value = pw_integral_a ( coeff % pw, flag )
  ELSE
     stop "integral_a, not defined"
  END IF

END FUNCTION integral_a

!******************************************************************************

FUNCTION integral_ab ( coeff1, coeff2, flag ) RESULT ( integral_value )

  TYPE ( coeff_type ), INTENT ( IN ) :: coeff1, coeff2
  INTEGER, INTENT ( IN ), OPTIONAL :: flag

  REAL ( dbl ) :: integral_value

!------------------------------------------------------------------------------

  IF ( coeff1 % in_use == PLANEWAVES .AND. coeff2 % in_use == PLANEWAVES ) THEN
     integral_value = pw_integral_ab ( coeff1 % pw, coeff2 % pw, flag )
  ELSE
     stop "integral_ab, not defined"
  END IF

END FUNCTION integral_ab

!******************************************************************************

END MODULE coefficients
