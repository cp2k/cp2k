!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

MODULE particle_lists
  
  USE kinds, ONLY : dbl
  USE molecule_types, ONLY : particle_node_type
  USE stop_program, ONLY : stop_memory
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: particle_list_type, particle_list_select
  
  TYPE particle_list_type
     INTEGER :: particle_index
     REAL ( dbl ), DIMENSION ( : ), POINTER :: r
     REAL ( dbl ), POINTER :: charge
     TYPE ( particle_list_type ), POINTER :: next
  END TYPE particle_list_type
  
CONTAINS

!******************************************************************************

SUBROUTINE particle_list_select ( pnode, plist, nparts)
  IMPLICIT NONE
  
! Arguments
  TYPE ( particle_node_type ), DIMENSION ( : ), INTENT ( IN ), TARGET :: pnode
  TYPE ( particle_list_type ), POINTER :: plist
  INTEGER , INTENT ( OUT )  :: nparts
  
! Locals
  INTEGER :: node, nnodes, isos
  TYPE ( particle_list_type ), POINTER :: plist_loc
  
!------------------------------------------------------------------------------
  
  nnodes = SIZE ( pnode )
  
  ! 'plist' will always point to the beginning of the list, 'plist_loc'
  !    to the current position
  NULLIFY ( plist )
  NULLIFY ( plist_loc )
  
  nparts = 0
  DO node = 1, nnodes
     IF ( pnode ( node ) % p % prop % charge /= 0.0_dbl ) THEN
        
        IF ( .NOT. ASSOCIATED ( plist_loc ) ) THEN
           ALLOCATE ( plist, STAT = isos )
           IF ( isos /= 0 ) &
                CALL stop_memory ( "get_atom_list", "plist", 0 )
           plist_loc => plist
           NULLIFY ( plist_loc % next )
        ELSE
           ALLOCATE ( plist_loc % next, STAT = isos )
           IF ( isos /= 0 ) &
                CALL stop_memory ( "get_atom_list", "next", 0 )
           
           plist_loc => plist_loc % next
           
           NULLIFY ( plist_loc % next )
        END IF
        
        nparts = nparts + 1
        plist_loc % particle_index = node
        plist_loc % r => pnode ( node ) % p % r
        plist_loc % charge => pnode ( node ) % p % prop % charge
        
     END IF
  END DO
  
END SUBROUTINE particle_list_select

!******************************************************************************
  
END MODULE particle_lists
