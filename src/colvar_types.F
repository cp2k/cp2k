!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/colvar_types [1.0] *
!!
!!   NAME
!!     colvar_types
!!
!!   FUNCTION
!!     Initialize the collective variables types
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     5.2004 created [fawzi and alessandro]
!!
!!   SOURCE
!****************************************************************************
MODULE colvar_types
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE 

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'colvar_types'

  PUBLIC :: colvar_type,&
            colvar_init,&
            colvar_dealloc_ref,&
            dist_colvar_type,&
            coord_colvar_type,&
            torsion_colvar_type,&
            angle_colvar_type,&
            plane_distance_colvar_type

  INTEGER, PARAMETER, PUBLIC :: no_colvar_id=-2,&
                                dist_colvar_id=1,&
                                coord_colvar_id=2,&
                                torsion_colvar_id=3,&
                                angle_colvar_id=4,&
                                plane_distance_colvar_id=5

!!****s* colvar/dist_colvar_param [1.0] *
!!
!!   NAME
!!     dist_colvar_param
!!
!!   FUNCTION
!!     parameters for the distance collective variable
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - i_at,j_at: indexes of the two atoms between which you calculate
!!       the distance
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!   SOURCE
!***************************************************************************
  TYPE dist_colvar_type
     INTEGER :: i_at,j_at
  END TYPE dist_colvar_type

  TYPE coord_colvar_type
     INTEGER :: n_atoms_to,n_atoms_from,nncrd,ndcrd
     INTEGER, POINTER, DIMENSION (:) :: i_at_from,i_at_to
     REAL(dp) :: r_0
  END TYPE coord_colvar_type

  TYPE torsion_colvar_type
     INTEGER, DIMENSION (4) :: i_at_tors
  END TYPE torsion_colvar_type

  TYPE plane_distance_colvar_type
     INTEGER, DIMENSION(3) :: plane
     INTEGER               :: point
  END TYPE plane_distance_colvar_type

  TYPE angle_colvar_type
     INTEGER, DIMENSION(3) :: i_at_angle
  END TYPE angle_colvar_type

!!****s* colvar_types/dist_colvar_param [1.0] *
!!
!!   NAME
!!     dist_colvar_param
!!
!!   FUNCTION
!!     parameters for a collective variable
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!   SOURCE
  !***************************************************************************
  TYPE colvar_type
     INTEGER :: type_id
     REAL(kind=dp)                           :: ss   ! value of the colvar
     REAL(kind=dp), DIMENSION (:,:), POINTER :: dsdr ! derivative of colvar 
     INTEGER, DIMENSION (:), POINTER :: i_atom       ! mapping of dsdr
     INTEGER                         :: n_atom_s
     TYPE(dist_colvar_type)    :: dist_param
     TYPE(coord_colvar_type)   :: coord_param
     TYPE(torsion_colvar_type) :: torsion_param
     TYPE(angle_colvar_type)         :: angle_param
     TYPE(plane_distance_colvar_type) :: plane_distance_param
  END TYPE colvar_type
!!***
  !****************************************************************************

CONTAINS

!!****f* colvar_types/colvar_init *
!!
!!   NAME
!!     colvar_init
!!
!!   FUNCTION
!!     initializes a colvar_param type
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - colvar: the colvat to initialize
!!     - index: the index of the colvar
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!*** **********************************************************************
SUBROUTINE colvar_init(colvar,error)
    TYPE(colvar_type), INTENT(inout)         :: colvar
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_init', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  colvar%type_id=no_colvar_id
  colvar%dist_param%i_at=-HUGE(0)
  colvar%dist_param%j_at=-HUGE(0)
  colvar%torsion_param%i_at_tors(1)=-HUGE(0)
  colvar%torsion_param%i_at_tors(2)=-HUGE(0)
  colvar%torsion_param%i_at_tors(3)=-HUGE(0)
  colvar%torsion_param%i_at_tors(4)=-HUGE(0)
  colvar%coord_param%n_atoms_from=-HUGE(0)
  colvar%coord_param%n_atoms_to=-HUGE(0)
  colvar%coord_param%nncrd=6
  colvar%coord_param%ndcrd=12
  colvar%coord_param%R_0=3.0_dp
  colvar%angle_param%i_at_angle=-HUGE(0)
  colvar%plane_distance_param%plane=-HUGE(0)
  colvar%plane_distance_param%point=-HUGE(0)
  NULLIFY(colvar%coord_param%i_at_from,colvar%coord_param%i_at_to)
  NULLIFY(colvar%i_atom,colvar%dsdr)
END SUBROUTINE colvar_init
!***************************************************************************

!!****f* colvar_types/colvar_dealloc_ref *
!!
!!   NAME
!!     colvar_dealloc_ref
!!
!!   FUNCTION
!!     releases the memory that might have been allocated by the colvar
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - colvar: the colvar to deallocate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!*** **********************************************************************
SUBROUTINE colvar_dealloc_ref(colvar,error)
    TYPE(colvar_type), INTENT(inout)         :: colvar
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_dealloc_ref', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  SELECT CASE(colvar%type_id)
  CASE (dist_colvar_id)
     DEALLOCATE(colvar%dsdr, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
     DEALLOCATE(colvar%i_atom, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)

     colvar%dist_param%i_at=-HUGE(0)
     colvar%dist_param%j_at=-HUGE(0)
  CASE (coord_colvar_id)
     DEALLOCATE(colvar%dsdr, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
     DEALLOCATE(colvar%i_atom, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)

     colvar%coord_param%n_atoms_from=-HUGE(0)
     colvar%coord_param%n_atoms_to=-HUGE(0)
     DEALLOCATE(colvar%coord_param%i_at_from,&
          colvar%coord_param%i_at_to,stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  CASE (torsion_colvar_id)
     DEALLOCATE(colvar%dsdr, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
     DEALLOCATE(colvar%i_atom, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)

     colvar%torsion_param%i_at_tors(1)=-HUGE(0)
     colvar%torsion_param%i_at_tors(2)=-HUGE(0)
     colvar%torsion_param%i_at_tors(3)=-HUGE(0)
     colvar%torsion_param%i_at_tors(4)=-HUGE(0)
  CASE(angle_colvar_id)
     DEALLOCATE(colvar%dsdr, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
     DEALLOCATE(colvar%i_atom, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)

     colvar%angle_param%i_at_angle=-HUGE(0)
  CASE(plane_distance_colvar_id)
     DEALLOCATE(colvar%dsdr, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
     DEALLOCATE(colvar%i_atom, stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)

     colvar%plane_distance_param%plane=-HUGE(0)
     colvar%plane_distance_param%point=-HUGE(0)
  CASE(no_colvar_id)
  CASE default
     CPAssert(.FALSE.,cp_failure_level,routinep,error,failure)
  END SELECT
  colvar%type_id=no_colvar_id
END SUBROUTINE colvar_dealloc_ref

END MODULE colvar_types
