!!
!!
!!@Author Christopher Mundy

MODULE hartree
  
  USE kinds, ONLY : dbl
  USE coefficient_types, ONLY : coeff_type, PLANEWAVES
  USE coefficients, ONLY : integral
  USE pws, ONLY : pw_poisson_solver
  
  IMPLICIT NONE
  
CONTAINS

!******************************************************************************

SUBROUTINE calculate_hartree ( n, ehartree, vhartree, dvhartree )
  
! Arguments
  TYPE ( coeff_type ), INTENT ( IN ) :: n   ! Total density
  
  REAL ( dbl ), INTENT ( OUT ) :: ehartree
  !! vxc is optional in case you only want the energy
  TYPE ( coeff_type ), INTENT ( OUT ), OPTIONAL :: vhartree
  TYPE ( coeff_type ), DIMENSION ( 3 ), INTENT ( OUT ), &
       OPTIONAL :: dvhartree
  
! Local
  REAL ( dbl ) :: eps0 ! TMPTMPTMP
  
!------------------------------------------------------------------------------
  
  IF ( PRESENT ( dvhartree ) ) THEN
     CALL poisson_solver ( n, vhartree, dvhartree )
  ELSE
     CALL poisson_solver ( n, vhartree )
  END IF
  
eps0 = 1._dbl !TMPTMPTMP
  ehartree = 0.5_dbl * integral ( n, vhartree ) / eps0
  
END SUBROUTINE calculate_hartree

!******************************************************************************

SUBROUTINE poisson_solver ( n, vhartree, dvhartree )
  
!! obtains electrostatic potential and its derivatives with respect to r
!!    from the density
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( IN ) :: n
  TYPE ( coeff_type ), INTENT ( OUT ), OPTIONAL :: vhartree
  TYPE ( coeff_type ), DIMENSION ( 3 ), INTENT ( OUT ), &
       OPTIONAL :: dvhartree
  
!------------------------------------------------------------------------------
  
  IF ( n % in_use == PLANEWAVES ) THEN
     IF ( PRESENT ( dvhartree ) ) THEN
        CALL pw_poisson_solver ( n, vhartree, dvhartree )
     ELSE
        CALL pw_poisson_solver ( n, vhartree )
     END IF
  ELSE
     stop "poisson_solver, no grid type implemented"
  END IF
  
END SUBROUTINE poisson_solver

!******************************************************************************

END MODULE hartree
