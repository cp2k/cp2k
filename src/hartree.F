!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/hartree [1.0] *
!!
!!   NAME
!!     hartree
!!
!!   FUNCTION
!!     Calculate Hartree energy and potential
!!
!!   AUTHOR
!!     CJM 
!!
!!   MODIFICATION HISTORY
!!     JGH (12-Mar-2001) : 
!!
!!   SOURCE
!******************************************************************************

MODULE hartree
  
  USE kinds, ONLY : dbl
  USE coefficient_types, ONLY : coeff_type, coeff_allocate, coeff_deallocate, &
      PLANEWAVES
  USE pws, ONLY : pw_poisson_solver
  USE termination, ONLY: stop_memory
  USE timings, ONLY: timeset, timestop
  
  IMPLICIT NONE

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** hartree/calculate_hartree [1.0] *
!!
!!   NAME
!!     calculate_hartree
!!
!!   FUNCTION
!!     
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (12-Mar-2001) : stress
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE calculate_hartree ( density, ehartree, vhartree, dvhartree, h_stress )
  
! Arguments
  TYPE ( coeff_type ), INTENT ( IN ) :: density
  REAL ( dbl ), INTENT ( OUT ), OPTIONAL :: ehartree
  TYPE ( coeff_type ), TARGET, INTENT ( INOUT ), OPTIONAL :: vhartree
  TYPE ( coeff_type ), DIMENSION ( 3 ), INTENT ( INOUT ), OPTIONAL :: dvhartree
  REAL ( dbl ), DIMENSION ( 3, 3 ), INTENT ( OUT ), OPTIONAL :: h_stress
  
! Local
  INTEGER :: handle, ierr
  TYPE ( coeff_type ), POINTER :: vh
  
!------------------------------------------------------------------------------
  
  CALL timeset ( 'hartree', 'I', ' ', handle )
  NULLIFY (vh)

  IF ( PRESENT ( vhartree ) ) THEN
    vh => vhartree
  ELSE
    ALLOCATE ( vh, STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "calculate_hartree", "vh", 0 )
    CALL coeff_allocate ( vh, density )
  END IF

  CALL poisson_solver ( density, ehartree, vh, dvhartree, h_stress )
  
  IF ( .NOT. PRESENT ( vhartree ) ) THEN
    CALL coeff_deallocate ( vh )
    DEALLOCATE ( vh, STAT = ierr )
    IF ( ierr /= 0 ) call stop_memory ( "calculate_hartree", "vh" )
  END IF

  CALL timestop ( 0.0_dbl, handle )
  
END SUBROUTINE calculate_hartree

!!*****
!******************************************************************************
!!****** hartree/poisson_solver [1.0] *
!!
!!   NAME
!!     poisson_solver
!!
!!   FUNCTION
!!
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     JGH (12-Mar-2001) : stress
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE poisson_solver ( density, ehartree, vhartree, dvhartree, h_stress )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( coeff_type ), INTENT ( IN ) :: density
  REAL ( dbl ), INTENT ( OUT ), OPTIONAL :: ehartree
  TYPE ( coeff_type ), INTENT ( INOUT ), OPTIONAL :: vhartree
  TYPE ( coeff_type ), DIMENSION ( 3 ), INTENT ( INOUT ), OPTIONAL :: dvhartree
  REAL ( dbl ), DIMENSION ( 3, 3 ), INTENT ( OUT ), OPTIONAL :: h_stress
  
!------------------------------------------------------------------------------
  
  IF ( density % in_use == PLANEWAVES ) THEN
    CALL pw_poisson_solver ( density, ehartree, vhartree, dvhartree, h_stress )
  ELSE
     stop "poisson_solver, no grid type implemented"
  END IF
  
END SUBROUTINE poisson_solver

!!*****
!******************************************************************************

END MODULE hartree

!******************************************************************************
