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
  
  USE kinds,                           ONLY: dbl
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ao_type
  PUBLIC :: ao_allocate, ao_deallocate, ao_structure_factor_type
  PUBLIC :: ao_zero

  TYPE ao_structure_factor_type
     COMPLEX ( dbl ), DIMENSION ( : ), POINTER :: ex, ey, ez
     INTEGER, DIMENSION ( : ), POINTER :: centre
  END TYPE ao_structure_factor_type

  TYPE ao_type
     REAL ( dbl ), DIMENSION ( : ), POINTER :: cr
     REAL ( dbl ), DIMENSION ( : ), POINTER :: norm
     REAL ( dbl ), DIMENSION ( : ), POINTER :: mass
     INTEGER, DIMENSION ( : ), POINTER :: coef_to_basis 
     INTEGER, DIMENSION ( : ), POINTER :: coef_to_set 
     INTEGER, DIMENSION ( : ), POINTER :: coef_to_part 
     INTEGER, DIMENSION ( : ), POINTER :: grid_type_id
     INTEGER, DIMENSION ( : ), POINTER :: grid_id
     TYPE ( ao_structure_factor_type ), DIMENSION ( : ),  POINTER :: eigr
  END TYPE ao_type

!!***
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
!!*** *************************************************************************

SUBROUTINE ao_allocate ( ao, ncoefs, mass_flag, grid_flag, eigr_flag )


    TYPE(ao_type), INTENT(INOUT), TARGET     :: ao
    INTEGER, INTENT(IN)                      :: ncoefs
    LOGICAL, INTENT(IN), OPTIONAL            :: mass_flag, grid_flag, &
                                                eigr_flag

    INTEGER                                  :: allocstat
    LOGICAL                                  :: local_eigr_flag, &
                                                local_grid_flag, &
                                                local_mass_flag

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

  IF ( PRESENT ( eigr_flag ) ) THEN
    local_eigr_flag=eigr_flag
  ELSE
    local_eigr_flag=.FALSE.
  ENDIF
   
  NULLIFY ( ao % cr )
  NULLIFY ( ao % coef_to_basis )
  NULLIFY ( ao % coef_to_set )
  NULLIFY ( ao % coef_to_part )
  NULLIFY ( ao % mass )
  NULLIFY ( ao % norm )
  NULLIFY ( ao % grid_id )
  NULLIFY ( ao % grid_type_id )
  NULLIFY ( ao % eigr )

  ALLOCATE ( ao % cr ( ncoefs ), STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "ao_allocate", &
     "ao % cr",  ncoefs )
  IF ( local_mass_flag ) THEN  
    ALLOCATE ( ao % mass ( ncoefs ), STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_allocate", &
       "ao % mass",  ncoefs )
  ENDIF

  IF ( local_eigr_flag ) THEN  
    ALLOCATE ( ao % eigr ( ncoefs ), STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_allocate", &
     "ao % eigr",  ncoefs )
    ALLOCATE ( ao % norm ( ncoefs ), STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_allocate", &
     "ao % norm",  ncoefs )
  ENDIF

  IF ( local_grid_flag ) THEN  
    ALLOCATE ( ao % grid_id ( ncoefs ), STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_allocate", &
       "ao % grid_id",  ncoefs )
    ALLOCATE ( ao % grid_type_id ( ncoefs ), STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_allocate", &
       "ao % grid_type_id",  ncoefs )
    ALLOCATE ( ao % coef_to_basis ( ncoefs ), STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_allocate", &
       "ao % coef_to_basis",  ncoefs )
    ALLOCATE ( ao % coef_to_set ( ncoefs ), STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_allocate", &
       "ao % coef_to_set",  ncoefs )
    ALLOCATE ( ao % coef_to_part ( ncoefs ), STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_allocate", &
       "ao % coef_to_part",  ncoefs )
  ENDIF


END SUBROUTINE ao_allocate

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
!!*** *************************************************************************

SUBROUTINE ao_deallocate ( ao )


    TYPE(ao_type), INTENT(INOUT)             :: ao

    INTEGER                                  :: allocstat

!------------------------------------------------------------------------------

  DEALLOCATE ( ao % cr, STAT = allocstat )
  IF ( allocstat /= 0 ) CALL stop_memory ( "ao_deallocate", &
     "ao % cr",  SIZE ( ao % cr ) )

  IF ( ASSOCIATED ( ao % coef_to_basis )) THEN
    DEALLOCATE ( ao % coef_to_basis, STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_deallocate", &
       "ao % coef_to_basis",  SIZE ( ao % coef_to_basis ) )
  ENDIF

  IF ( ASSOCIATED ( ao % coef_to_set )) THEN
    DEALLOCATE ( ao % coef_to_set, STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_deallocate", &
       "ao % coef_to_set",  SIZE ( ao % coef_to_set ) )
  ENDIF

  IF ( ASSOCIATED ( ao % coef_to_part )) THEN
    DEALLOCATE ( ao % coef_to_part, STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_deallocate", &
       "ao % coef_to_part",  SIZE ( ao % coef_to_part ) )
  ENDIF

  IF ( ASSOCIATED ( ao % mass )) THEN
    DEALLOCATE ( ao % mass, STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_deallocate", &
       "ao % mass",  SIZE ( ao % mass ) )
  ENDIF

  IF ( ASSOCIATED ( ao % grid_id )) THEN
    DEALLOCATE ( ao % grid_id, STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_deallocate", &
       "ao % grid_id",  SIZE ( ao % grid_id ) )
  ENDIF

  IF ( ASSOCIATED ( ao % grid_type_id )) THEN
    DEALLOCATE ( ao % grid_type_id, STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_deallocate", &
       "ao % grid_type_id",  SIZE ( ao % grid_type_id ) )
  ENDIF
  IF ( ASSOCIATED ( ao % eigr )) THEN
    DEALLOCATE ( ao % eigr, STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_deallocate", &
       "ao % eigr",  SIZE ( ao % eigr ) )
  ENDIF
  IF ( ASSOCIATED ( ao % norm )) THEN
    DEALLOCATE ( ao % norm, STAT = allocstat )
    IF ( allocstat /= 0 ) CALL stop_memory ( "ao_deallocate", &
       "ao % norm",  SIZE ( ao % norm ) )
  ENDIF

END SUBROUTINE ao_deallocate

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
!!*** *************************************************************************

SUBROUTINE ao_zero ( ao )


    TYPE(ao_type), INTENT(INOUT)             :: ao

    INTEGER                                  :: handle, ns
    REAL(dbl)                                :: zr

!------------------------------------------------------------------------------

  CALL timeset ( 'AO_ZERO', 'I', 'Mzero', handle )

  ns = SIZE ( ao % cr)
  ao % cr = 0.0_dbl
  IF ( ASSOCIATED ( ao % coef_to_basis ) ) ao % coef_to_basis = 0
  IF ( ASSOCIATED ( ao % coef_to_set ) ) ao % coef_to_set = 0
  IF ( ASSOCIATED ( ao % coef_to_part ) ) ao % coef_to_part = 0
  IF ( ASSOCIATED ( ao % mass ) ) ao % mass = 0.0_dbl
  IF ( ASSOCIATED ( ao % norm ) ) ao % norm = 0.0_dbl
  IF ( ASSOCIATED ( ao % grid_id ) ) ao % grid_id = 0
  IF ( ASSOCIATED ( ao % grid_type_id ) ) ao % grid_type_id = 0
  zr = REAL ( ns, dbl ) * 1.e-6_dbl
  CALL timestop ( zr, handle )

END SUBROUTINE ao_zero

!******************************************************************************

END MODULE ao_types

!******************************************************************************
