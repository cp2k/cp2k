!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/coefficient_list [1.0] *
!!
!!   NAME
!!     coefficient_lists
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     gloria
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE coefficient_lists

  USE ao_types,                        ONLY: ao_type
  USE kinds,                           ONLY: dbl
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: coef_list_type, coef_list_grid      

  TYPE coef_list_type
     INTEGER :: coef_index
     INTEGER :: particle_index
     INTEGER :: grid_type_id
     INTEGER :: grid_id
     REAL ( dbl ), DIMENSION ( : ), POINTER :: r
     REAL ( dbl ), POINTER :: c
     TYPE ( coef_list_type ), POINTER :: next
  END TYPE coef_list_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE coef_list_grid ( ao, clist, part, ncoef, grid_type_tag, &
                             grid_tag)

    TYPE(ao_type), INTENT(IN)                :: ao
    TYPE(coef_list_type), POINTER            :: clist
    TYPE(particle_type), DIMENSION(:), &
      intent(IN), TARGET                     :: PART
    INTEGER, INTENT(OUT)                     :: ncoef
    INTEGER, INTENT(IN)                      :: grid_type_tag, grid_tag

    INTEGER                                  :: ic, ipart, isos, nc
    TYPE(coef_list_type), POINTER            :: clist_loc

!------------------------------------------------------------------------------

  nc = SIZE ( ao%cr )

  ! 'plist' will always point to the beginning of the list, 'plist_loc'
  !    to the current position
  NULLIFY ( clist )
  NULLIFY ( clist_loc )

  ncoef = 0
  DO ic = 1, nc
     IF ( ao%grid_type_id ( ic ) == grid_type_tag ) THEN

        IF ( .NOT. ASSOCIATED ( clist_loc ) ) THEN
           ALLOCATE ( clist, STAT = isos )
           IF ( isos /= 0 ) &
                CALL stop_memory ( "get_atom_list", "plist", 0 )
           clist_loc => clist
           NULLIFY ( clist_loc % next )
        ELSE
           ALLOCATE ( clist_loc % next, STAT = isos )
           IF ( isos /= 0 ) &
                CALL stop_memory ( "get_atom_list", "next", 0 )

           clist_loc => clist_loc % next

           NULLIFY ( clist_loc % next )
        END IF

        ncoef = ncoef + 1
        clist_loc % coef_index = ic 
        ipart = ao % coef_to_part ( ic )  
        clist_loc % particle_index =  ipart
        clist_loc % grid_type_id = grid_type_tag
        clist_loc % grid_id = grid_tag
        clist_loc % r => part ( ipart ) % r
        clist_loc % c => ao  % cr (ic)              

     END IF
  END DO

END SUBROUTINE coef_list_grid

!******************************************************************************

END MODULE coefficient_lists
