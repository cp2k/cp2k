!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!***** cp2k/ao_types [1.0] *
!!
!!   NAME
!!     ao_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!   NOTES
!!
!!
!!   SOURCE
!******************************************************************************

MODULE ao_types
  
  USE atomic_kinds, ONLY :  kind_info_type
  USE kinds, ONLY : dbl
  USE termination, ONLY : stop_memory, stop_program
  USE timings, ONLY : timeset, timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ao_type
  PUBLIC :: ao_allocate, ao_deallocate
  PUBLIC :: ao_zero

  TYPE ao_type
     REAL ( dbl ), DIMENSION ( : ), POINTER :: cr
     REAL ( dbl ), DIMENSION ( : ), POINTER :: mass
     INTEGER, DIMENSION ( : ), POINTER :: coef_to_basis_set 
     INTEGER, DIMENSION ( : ), POINTER :: coef_to_part 
     INTEGER, DIMENSION ( : ), POINTER :: grid_id
     TYPE ( kind_info_type ), DIMENSION ( : ),  POINTER :: kind_info
  END TYPE ao_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** ao_types/ao_allocate [1.0] *
!!
!!   NAME
!!     ao_allocate
!!
!!   FUNCTION
!!     Allocate a atomic oribital type
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE ao_allocate ( ao, ncoefs, mass_flag, grid_flag )

  IMPLICIT NONE

! Arguments
  TYPE ( ao_type ), INTENT ( INOUT ), TARGET :: ao
  INTEGER, INTENT ( IN ) :: ncoefs
  LOGICAL, INTENT ( IN ), OPTIONAL :: mass_flag, grid_flag

! Locals
  INTEGER :: allocstat
  LOGICAL :: local_mass_flag, local_grid_flag

!------------------------------------------------------------------------------

  IF ( PRESENT ( mass_flag ) ) THEN
    local_mass_flag=mass_flag
  ELSE
    local_mass_flag=.FALSE.
  ENDIF

  IF ( PRESENT ( grid_flag ) ) THEN
    local_grid_flag=grid_flag
  ELSE
    local_grid_flag=.FALSE.
  ENDIF
   
  NULLIFY ( ao % cr, ao % coef_to_basis_set, ao % coef_to_part )
  NULLIFY ( ao % mass )
  NULLIFY ( ao % grid_id )

  ALLOCATE ( ao % cr ( ncoefs ), STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "ao_allocate", &
     "ao % cr",  ncoefs )
  ALLOCATE ( ao % coef_to_basis_set ( ncoefs ), STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "ao_allocate", &
     "ao % coef_to_basis_set",  ncoefs )
  ALLOCATE ( ao % coef_to_part ( ncoefs ), STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "ao_allocate", &
     "ao % coef_to_part",  ncoefs )
  IF ( local_mass_flag ) THEN  
    ALLOCATE ( ao % mass ( ncoefs ), STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_allocate", &
       "ao % mass",  ncoefs )
  ENDIF
  IF ( local_grid_flag ) THEN  
    ALLOCATE ( ao % grid_id ( ncoefs ), STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_allocate", &
       "ao % grid_id",  ncoefs )
  ENDIF


END SUBROUTINE ao_allocate

!!*****
!******************************************************************************
!!****** ao_types/ao_deallocate [1.0] *
!!
!!   NAME
!!     ao_deallocate
!!
!!   FUNCTION
!!     Deallocate a atomic orbital type variable
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE ao_deallocate ( ao )

  IMPLICIT NONE

! Arguments
  TYPE ( ao_type ), INTENT ( INOUT ) :: ao

! Locals
  INTEGER :: allocstat

!------------------------------------------------------------------------------

  DEALLOCATE ( ao % cr, STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "ao_deallocate", &
     "ao % cr",  SIZE ( ao % cr ) )

  DEALLOCATE ( ao % coef_to_basis_set, STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "ao_deallocate", &
     "ao % coef_to_basis_set",  SIZE ( ao % coef_to_basis_set ) )

  DEALLOCATE ( ao % coef_to_part, STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "ao_deallocate", &
     "ao % coef_to_part",  SIZE ( ao % coef_to_part ) )

  IF ( ASSOCIATED ( ao % mass )) THEN
    DEALLOCATE ( ao % mass, STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_deallocate", &
       "ao % mass",  SIZE ( ao % mass ) )
  ENDIF

  IF ( ASSOCIATED ( ao % grid_id )) THEN
    DEALLOCATE ( ao % mass, STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_deallocate", &
       "ao % grid_id",  SIZE ( ao % grid_id ) )
  ENDIF

END SUBROUTINE ao_deallocate

!!*****
!******************************************************************************
!!****** ao_types/ao_zero [1.0] *
!!
!!   NAME
!!     ao_zero
!!
!!   FUNCTION
!!     Set values of a ao type to zero
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE ao_zero ( ao )

  IMPLICIT NONE

! Arguments
  TYPE ( ao_type ), INTENT ( INOUT ) :: ao

! Local
  INTEGER :: ns, handle
  REAL ( dbl ) :: zr

!------------------------------------------------------------------------------

  CALL timeset ( 'AO_ZERO', 'I', 'Mzero', handle )

  ns = SIZE ( ao % cr)
  ao % cr = 0.0_dbl
  ao % coef_to_basis_set = 0
  IF ( ASSOCIATED ( ao % mass ) ) ao % mass = 0.0_dbl
  IF ( ASSOCIATED ( ao % grid_id ) ) ao % grid_id = 0
  zr = REAL ( ns, dbl ) * 1.e-6_dbl
  CALL timestop ( zr, handle )

END SUBROUTINE ao_zero

!!*****
!******************************************************************************

END MODULE ao_types

!******************************************************************************
