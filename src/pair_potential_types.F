!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/pair_potential_types [1.0] *
!!
!!   NAME
!!     pair_potential_types
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

MODULE pair_potential_types

  USE kinds,                           ONLY: dbl,&
                                             default_string_length
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: module_name = "pair_potential_types"

  PRIVATE
  PUBLIC :: pair_potential_type
  PUBLIC :: lj_pot_type, williams_pot_type, goodwin_pot_type

  PUBLIC :: reallocate_potparm

  TYPE lj_pot_type
    REAL ( dbl ) :: epsilon
    REAL ( dbl ) :: sigma6
    REAL ( dbl ) :: sigma12
  END TYPE lj_pot_type

  TYPE williams_pot_type
    REAL ( dbl ) :: a
    REAL ( dbl ) :: b
    REAL ( dbl ) :: c
  END TYPE williams_pot_type

  TYPE goodwin_pot_type
    REAL ( dbl ) :: vr0
    REAL ( dbl ) :: m, mc
    REAL ( dbl ) :: d, dc
  END TYPE goodwin_pot_type

  TYPE pair_potential_type
    CHARACTER ( LEN = default_string_length ) :: type
    CHARACTER ( LEN = default_string_length ) :: at1
    CHARACTER ( LEN = default_string_length ) :: at2
    TYPE (lj_pot_type ) :: lj
    TYPE (williams_pot_type ) :: willis
    TYPE (goodwin_pot_type ) :: goodwin
    REAL ( dbl ) :: rcutsq
    REAL ( dbl ) :: energy_cutoff
    REAL ( dbl ) :: e_cutoff_coul
  END TYPE pair_potential_type

!!*****
!******************************************************************************

CONTAINS

  SUBROUTINE reallocate_potparm(p,lb1_new,ub1_new)

!   ***************************************************************************

    TYPE (pair_potential_type), POINTER       :: p(:)
    INTEGER, INTENT(IN)                :: lb1_new, ub1_new

    CHARACTER(LEN=*), PARAMETER        :: routine_name = "reallocate_potparm"

    INTEGER                            :: istat, lb1, lb1_old, &
                                          new_size, old_size, ub1, &
                                          ub1_old
    TYPE (pair_potential_type), ALLOCATABLE :: work(:)

!   ---------------------------------------------------------------------------

    IF (ASSOCIATED(p)) THEN
      lb1_old = LBOUND(p,1)
      ub1_old = UBOUND(p,1)
      ALLOCATE (work(lb1_old:ub1_old),STAT=istat)
      IF (istat /= 0) THEN
      END IF
      work(:) = p(:)
      DEALLOCATE (p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine_name,module_name,__LINE__,"p")
    END IF

    ALLOCATE (p(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
      new_size = (ub1_new - lb1_new + 1)
      CALL stop_memory(routine_name,module_name,__LINE__,"p",new_size)
    END IF

    p ( : ) % type = 'NULL'
    p ( : ) % lj % epsilon = 0._dbl
    p ( : ) % lj % sigma6 = 0._dbl
    p ( : ) % lj % sigma12 = 0._dbl
    p ( : ) % willis % a = 0._dbl
    p ( : ) % willis % b = 0._dbl
    p ( : ) % willis % c = 0._dbl
    p ( : ) % goodwin % vr0 = 0._dbl
    p ( : ) % goodwin % d = 0._dbl
    p ( : ) % goodwin % dc = 0._dbl
    p ( : ) % goodwin % m = 0._dbl
    p ( : ) % goodwin % mc = 0._dbl

    IF (ALLOCATED(work)) THEN
      lb1 = MAX(lb1_new,lb1_old)
      ub1 = MIN(ub1_new,ub1_old)
      p(lb1:ub1) = work(lb1:ub1)
      DEALLOCATE (work,STAT=istat)
      IF (istat /= 0) THEN
        CALL stop_memory(routine_name,module_name,__LINE__,"work")
      END IF
    END IF

  END SUBROUTINE reallocate_potparm

END MODULE pair_potential_types
