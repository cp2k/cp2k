!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/linklist_utilities [1.0] *
!!
!!   NAME
!!     linklist_utilities
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!    gloria  8 may 2001 : linklist_verlet_list : added argument list_type,
!!    adapt the routine to be used as well for the overlaps lists
!!    cjm & haf 27 july 2001: fixed problem with cutoffs larger
!!    than half the boxsize and absolute coordinates far of the original box
!!   SOURCE
!******************************************************************************

MODULE linklist_utilities

  IMPLICIT NONE

  PUBLIC :: get_match
  PRIVATE

!!*****
!------------------------------------------------------------------------------!

  CONTAINS

!------------------------------------------------------------------------------!
  FUNCTION get_match(jatom,elist,nexcl) RESULT (match)
    INTEGER, INTENT(IN)                      :: jatom, elist( : ), nexcl
    LOGICAL                                  :: match

    INTEGER                                  :: iexclude

    match = .FALSE.
    DO iexclude = 1, nexcl
      IF (elist(iexclude)>jatom) EXIT
      IF (elist(iexclude)==jatom) THEN
        match = .TRUE.
        EXIT
      END IF
    END DO
  END FUNCTION get_match
!------------------------------------------------------------------------------!
END MODULE linklist_utilities
