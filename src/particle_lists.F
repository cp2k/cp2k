!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/particle_lists [1.0] *
!!
!!   NAME
!!     particle_lists
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE particle_lists

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE kinds,                           ONLY: dbl
  USE molecule_types,                  ONLY: particle_node_type
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: particle_list_type, particle_list_select, particle_list_grid 

  TYPE particle_list_type
     INTEGER :: particle_index
     INTEGER :: grid_id
     REAL ( dbl ), DIMENSION ( : ), POINTER :: r
     REAL ( dbl ), POINTER :: charge
     TYPE ( particle_list_type ), POINTER :: next
  END TYPE particle_list_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE particle_list_select ( pnode, plist, nparts)

    TYPE(particle_node_type), DIMENSION(:), &
      INTENT(IN), TARGET                     :: pnode
    TYPE(particle_list_type), POINTER        :: plist
    INTEGER, INTENT(OUT)                     :: nparts

    INTEGER                                  :: isos, nnodes, node
    REAL(dbl), TARGET                        :: charge
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(particle_list_type), POINTER        :: plist_loc

!------------------------------------------------------------------------------

  nnodes = SIZE ( pnode )

  ! 'plist' will always point to the beginning of the list, 'plist_loc'
  !    to the current position
  NULLIFY ( plist )
  NULLIFY ( plist_loc )

  nparts = 0
  DO node = 1, nnodes
     atomic_kind => pnode ( node ) % p % atomic_kind
     CALL get_atomic_kind (atomic_kind=atomic_kind, qeff=charge)
     IF ( charge /= 0.0_dbl ) THEN

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
        plist_loc % charge => charge

     END IF
  END DO

END SUBROUTINE particle_list_select

!****************************************************************************

SUBROUTINE particle_list_grid ( part, plist, atom_list, nparts, grid_tag )

    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT), TARGET                  :: part
    TYPE(particle_list_type), POINTER        :: plist
    INTEGER, DIMENSION(:), INTENT(IN)        :: ATOM_LIST
    INTEGER, INTENT(OUT)                     :: nparts
    INTEGER, INTENT(IN)                      :: grid_tag

    INTEGER                                  :: isos, npart, part_index
    REAL(dbl), TARGET                        :: charge
    TYPE(particle_list_type), POINTER        :: plist_loc

!------------------------------------------------------------------------------
! 'plist' will always point to the beginning of the list, 'plist_loc'
!    to the current position

  NULLIFY ( plist )
  NULLIFY ( plist_loc )
 
  nparts = 0
  charge = 1.0_dbl
  DO npart = 1, SIZE ( atom_list )
 
        IF ( .NOT. ASSOCIATED ( plist_loc ) ) THEN
           ALLOCATE ( plist, STAT = isos )
           IF ( isos /= 0 ) &
                CALL stop_memory ( "build_particle_list", "plist", 0 )
           plist_loc => plist
           NULLIFY ( plist_loc % next )
        ELSE
           ALLOCATE ( plist_loc % next, STAT = isos )
           IF ( isos /= 0 ) &
                CALL stop_memory ( "build_particle_list", "next", 0 )
 
           plist_loc => plist_loc % next
 
           NULLIFY ( plist_loc % next )
        END IF
 
        nparts = nparts + 1
        part_index = atom_list (npart)
        part ( part_index ) %  grid_id = grid_tag
        plist_loc % particle_index = atom_list (npart)
        plist_loc % grid_id = grid_tag
        plist_loc % r => part ( part_index )  % r
        plist_loc % charge =>  charge 
 
  END DO
 
END SUBROUTINE particle_list_grid
                                                           

!******************************************************************************

END MODULE particle_lists
