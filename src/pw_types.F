!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****s* cp2k/pw_types [1.0] *
!!
!!   NAME
!!     pw_types
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

MODULE pw_types
  
  USE kinds, ONLY : dbl
  USE pw_grid_types, ONLY : pw_grid_type
  USE stop_program, ONLY : stop_prg, stop_memory
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: pw_type, pw_allocate, pw_deallocate, pw_zero
  PUBLIC :: REALDATA1D, COMPLEXDATA1D, REALDATA3D, COMPLEXDATA3D
  PUBLIC :: NOSPACE, REALSPACE, RECIPROCALSPACE
  
  TYPE pw_type
     REAL ( dbl ), DIMENSION ( : ), POINTER :: cr
     REAL ( dbl ), DIMENSION ( :, :, : ), POINTER :: cr3d
     COMPLEX ( dbl ), DIMENSION ( : ), POINTER :: cc
     COMPLEX ( dbl ), DIMENSION ( :, :, : ), POINTER :: cc3d
     
     TYPE ( pw_grid_type ), POINTER :: pw_grid
     
     INTEGER :: in_use              ! Which data is used [r1d/c1d/r3d/c3d]
     INTEGER :: in_space            ! Real/Reciprocal space
     
     INTEGER :: group_id            ! In which parallel group am I
     LOGICAL :: group_head          ! Whether I am the head of this group
     INTEGER :: group_head_id       ! The group id for the head of this group
  END TYPE pw_type
  
  !! Flags for the structure member 'in_use'
  INTEGER, PARAMETER :: REALDATA1D = 301, COMPLEXDATA1D = 302
  INTEGER, PARAMETER :: REALDATA3D = 303, COMPLEXDATA3D = 304
  
  !! Flags for the structure member 'in_space'
  INTEGER, PARAMETER :: NOSPACE = 371, REALSPACE = 372, RECIPROCALSPACE = 373
  
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE pw_allocate ( pw, pw_grid, use_data )
  
  IMPLICIT NONE
  
! Arguments
  INTEGER, INTENT ( IN ) :: use_data
  TYPE ( pw_grid_type ), INTENT ( IN ), TARGET :: pw_grid
  TYPE ( pw_type ), INTENT ( INOUT ), TARGET :: pw
  
! Locals
  INTEGER :: allocstat
  INTEGER, DIMENSION ( :, : ), POINTER :: bounds
  
!------------------------------------------------------------------------------
  
  pw % in_use = use_data
  pw % pw_grid => pw_grid
  pw % in_space = NOSPACE
  bounds => pw % pw_grid % bounds
  
  NULLIFY ( pw % cr, pw % cc, pw % cr3d, pw % cc3d )
  
  IF ( use_data == REALDATA1D ) THEN
     ALLOCATE ( pw % cr ( pw % pw_grid % ngpts_cut ), STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_allocate", "cr", pw % pw_grid % ngpts_cut )
     END IF
     
  ELSE IF ( use_data == COMPLEXDATA1D ) THEN
     ALLOCATE ( pw % cc ( pw % pw_grid % ngpts_cut ), STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_allocate", "cc", pw % pw_grid % ngpts_cut )
     END IF
     
  ELSE IF ( use_data == REALDATA3D ) THEN
     ALLOCATE ( pw % cr3d ( &
          bounds ( 1, 1 ) : bounds ( 2, 1 ), &
          bounds ( 1, 2 ) : bounds ( 2, 2 ), &
          bounds ( 1, 3 ) : bounds ( 2, 3 ) ), STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_allocate", "cr3d", 0 )
     END IF
     
  ELSE IF ( use_data == COMPLEXDATA3D ) THEN
     ALLOCATE ( pw % cc3d ( &
          bounds ( 1, 1 ) : bounds ( 2, 1 ), &
          bounds ( 1, 2 ) : bounds ( 2, 2 ), &
          bounds ( 1, 3 ) : bounds ( 2, 3 ) ), STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_allocate", "cc3d", 0 )
     END IF
     
  ELSE
     CALL stop_prg ( "pw_allocate", "no possible data field" )
  END IF
  
END SUBROUTINE pw_allocate

!******************************************************************************

SUBROUTINE pw_deallocate ( pw )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( pw_type ), INTENT ( INOUT ) :: pw
  
! Locals
  INTEGER :: allocstat
  
!------------------------------------------------------------------------------
  
  IF ( pw % in_use == REALDATA1D ) THEN
     DEALLOCATE ( pw % cr, STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_deallocate", "cr" )
     END IF
  ELSE IF ( pw % in_use == COMPLEXDATA1D ) THEN
     DEALLOCATE ( pw % cc, STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_deallocate", "cc" )
     END IF
  ELSE IF ( pw % in_use == REALDATA3D ) THEN
     DEALLOCATE ( pw % cr3d, STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_deallocate", "cr3d" )
     END IF
  ELSE IF ( pw % in_use == COMPLEXDATA3D ) THEN
     DEALLOCATE ( pw % cc3d, STAT = allocstat )
     IF ( allocstat /= 0 ) THEN
        CALL stop_memory ( "pw_deallocate", "cc3d" )
     END IF
  ELSE
     CALL stop_prg ( "pw_deallocate", "no possible data field" )
  END IF
  
END SUBROUTINE pw_deallocate

!******************************************************************************

SUBROUTINE pw_zero ( pw )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( pw_type ), INTENT ( INOUT ) :: pw
  
!------------------------------------------------------------------------------
  
  IF ( pw % in_use == REALDATA1D ) THEN
     pw % cr = 0.0_dbl
  ELSE IF ( pw % in_use == COMPLEXDATA1D ) THEN
     pw % cc = 0.0_dbl
  ELSE IF ( pw % in_use == REALDATA3D ) THEN
     pw % cr3d = 0.0_dbl
  ELSE IF ( pw % in_use == COMPLEXDATA3D ) THEN
     pw % cc3d = 0.0_dbl
  ELSE
     CALL stop_prg ( "pw_zero", "no possible data field" )
  END IF
  
END SUBROUTINE pw_zero

!******************************************************************************



!******************************************************************************
  
END MODULE pw_types
