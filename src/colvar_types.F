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
            colvar_p_type,&
            colvar_create,&
            colvar_setup,&
            colvar_release,&
            dist_colvar_type,&
            coord_colvar_type,&
            torsion_colvar_type,&
            angle_colvar_type,&
            plane_distance_colvar_type

  INTEGER, PARAMETER, PUBLIC :: no_colvar_id              =-2,&
                                dist_colvar_id            =1,&
                                coord_colvar_id           =2,&
                                torsion_colvar_id         =3,&
                                angle_colvar_id           =4,&
                                plane_distance_colvar_id  =5

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
     INTEGER       :: i_at, j_at
  END TYPE dist_colvar_type

  TYPE coord_colvar_type
     INTEGER                         :: n_atoms_to,&
                                        n_atoms_from,&
                                        nncrd,&
                                        ndcrd
     INTEGER, POINTER, DIMENSION (:) :: i_at_from,&
                                        i_at_to
     REAL(KIND=dp)                   :: r_0
  END TYPE coord_colvar_type

  TYPE torsion_colvar_type
     INTEGER, DIMENSION (4)          :: i_at_tors
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
     REAL(kind=dp)                            :: ss               ! Value of the colvar
     REAL(kind=dp), DIMENSION (:,:), POINTER  :: dsdr             ! Derivative of colvar 
     INTEGER, DIMENSION (:), POINTER          :: i_atom           ! Mapping of dsdr
     INTEGER                                  :: n_atom_s
     TYPE(dist_colvar_type),POINTER           :: dist_param
     TYPE(coord_colvar_type),POINTER          :: coord_param
     TYPE(torsion_colvar_type),POINTER        :: torsion_param
     TYPE(angle_colvar_type),POINTER          :: angle_param
     TYPE(plane_distance_colvar_type),POINTER :: plane_distance_param
  END TYPE colvar_type

  TYPE colvar_p_type
     TYPE(colvar_type), POINTER :: colvar
  END TYPE colvar_p_type
  !!***
  !****************************************************************************

CONTAINS

  !!****f* colvar_types/colvar_create *
  !!
  !!   NAME
  !!     colvar_create
  !!
  !!   FUNCTION
  !!     initializes a colvar_param type
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - colvar: the colvat to initialize
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     alessandro laio and fawzi mohamed
  !!
  !!*** **********************************************************************
  SUBROUTINE colvar_create(colvar,colvar_id,error)
    TYPE(colvar_type), POINTER               :: colvar
    INTEGER, INTENT(in)                      :: colvar_id
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    INTEGER :: stat

    failure=.FALSE.
    CPPostcondition(.NOT.ASSOCIATED(colvar),cp_warning_level,routineP,error,failure)
    ALLOCATE(colvar, stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    !
    NULLIFY(colvar%dist_param)
    NULLIFY(colvar%angle_param)
    NULLIFY(colvar%torsion_param)
    NULLIFY(colvar%plane_distance_param)
    NULLIFY(colvar%coord_param)
    NULLIFY(colvar%i_atom)
    NULLIFY(colvar%dsdr)
    IF (.NOT.failure) THEN
       colvar%type_id=colvar_id
       SELECT CASE (colvar_id)
       CASE(dist_colvar_id)
          ALLOCATE(colvar%dist_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(coord_colvar_id)
          ALLOCATE(colvar%coord_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure) 
       CASE(angle_colvar_id)
          ALLOCATE(colvar%angle_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)           
       CASE(torsion_colvar_id)
          ALLOCATE(colvar%torsion_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)           
       CASE(plane_distance_colvar_id)
          ALLOCATE(colvar%plane_distance_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)           
       CASE(no_colvar_id)
          ! Do nothing
       CASE DEFAULT
          CPPostcondition(.FALSE.,cp_warning_level,routineP,error,failure)          
       END SELECT
    END IF

  END SUBROUTINE colvar_create
  !***************************************************************************
  
  !!****f* colvar_types/colvar_setup *
  !!
  !!   NAME
  !!     colvar_setup
  !!
  !!   FUNCTION
  !!     Finalize the setup of the collective variable
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - colvar: the colvat to initialize
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino, [teo] 09.03.2006
  !!
  !!*** **********************************************************************
  SUBROUTINE colvar_setup(colvar,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_setup', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    INTEGER :: stat, ii

    CPPostcondition(ASSOCIATED(colvar),cp_warning_level,routineP,error,failure)
    failure=.FALSE.
    IF (.NOT.failure) THEN

       SELECT CASE(colvar%type_id)
       CASE(dist_colvar_id)
          colvar%n_atom_s=2
       CASE(coord_colvar_id)
          colvar%n_atom_s=colvar%coord_param%n_atoms_from+colvar%coord_param%n_atoms_to
       CASE(angle_colvar_id)
          colvar%n_atom_s=3
       CASE(torsion_colvar_id)
          colvar%n_atom_s=4
       CASE(plane_distance_colvar_id)
          colvar%n_atom_s=4
       END SELECT

       ALLOCATE(colvar%dsdr(3,colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(colvar%i_atom(colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)

       SELECT CASE(colvar%type_id)
       CASE(dist_colvar_id)
          colvar%i_atom(1)=colvar%dist_param%i_at
          colvar%i_atom(2)=colvar%dist_param%j_at
       CASE(coord_colvar_id)
          DO ii=1,colvar%coord_param%n_atoms_from
             colvar%i_atom(ii)=colvar%coord_param%i_at_from(ii)
          ENDDO
          DO ii=1,colvar%coord_param%n_atoms_to
             colvar%i_atom(colvar%coord_param%n_atoms_from+ii)=colvar%coord_param%i_at_to(ii)
          ENDDO
       CASE(angle_colvar_id)
          DO ii=1,3
             colvar%i_atom(ii)=colvar%angle_param%i_at_angle(ii)
          ENDDO
       CASE(torsion_colvar_id)
          DO ii=1,4
             colvar%i_atom(ii)=colvar%torsion_param%i_at_tors(ii)
          ENDDO
       CASE(plane_distance_colvar_id)
          DO ii=1,3
             colvar%i_atom(ii)=colvar%plane_distance_param%plane(ii)
          ENDDO
          colvar%i_atom(4)=colvar%plane_distance_param%point
       END SELECT
       
    END IF    
  END SUBROUTINE colvar_setup
  !!****f* colvar_types/colvar_release *
  !!
  !!   NAME
  !!     colvar_release
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
  SUBROUTINE colvar_release(colvar,error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPostcondition(ASSOCIATED(colvar),cp_warning_level,routineP,error,failure)
    IF (ASSOCIATED(colvar%dsdr)) THEN
       DEALLOCATE(colvar%dsdr, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    END IF
    IF (ASSOCIATED(colvar%i_atom)) THEN
       DEALLOCATE(colvar%i_atom, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    END IF

    SELECT CASE(colvar%type_id)
    CASE (dist_colvar_id)
       DEALLOCATE(colvar%dist_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (coord_colvar_id)
       DEALLOCATE(colvar%coord_param%i_at_from, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%coord_param%i_at_to,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       DEALLOCATE(colvar%coord_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE(angle_colvar_id)
       DEALLOCATE(colvar%angle_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)     
    CASE (torsion_colvar_id)
       DEALLOCATE(colvar%torsion_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE(plane_distance_colvar_id)
       DEALLOCATE(colvar%plane_distance_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE(no_colvar_id)
       ! Do nothing
    CASE default
       CPAssert(.FALSE.,cp_failure_level,routinep,error,failure)
    END SELECT
    DEALLOCATE(colvar, stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)

  END SUBROUTINE colvar_release

END MODULE colvar_types
