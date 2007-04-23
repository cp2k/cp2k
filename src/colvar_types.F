!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
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
  USE input_constants,                 ONLY: do_clv_fix_point,&
                                             do_clv_geo_center
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE 

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'colvar_types'

  PUBLIC :: colvar_type,&
            colvar_p_type,&
            colvar_p_reallocate,&
            colvar_create,&
            colvar_clone,&
            colvar_setup,&
            colvar_release,&
            colvar_counters,&
            dist_colvar_type,&
            coord_colvar_type,&
            torsion_colvar_type,&
            angle_colvar_type,&
            plane_distance_colvar_type,&
            rotation_colvar_type,&
            ddiff_colvar_type,&
            qparm_colvar_type,&
            point_type,&
            colvar_check_points,&
            eval_point_der,&
            eval_point_pos

  INTEGER, PARAMETER, PUBLIC :: no_colvar_id              =-2,&
                                dist_colvar_id            = 1,&
                                coord_colvar_id           = 2,&
                                torsion_colvar_id         = 3,&
                                angle_colvar_id           = 4,&
                                plane_distance_colvar_id  = 5,&
                                rotation_colvar_id        = 6,&
                                ddiff_colvar_id           = 7,&
                                qparm_colvar_id           = 8

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
     REAL(KIND=dp)                   :: o0
     INTEGER, DIMENSION (4)          :: i_at_tors
  END TYPE torsion_colvar_type

  TYPE plane_distance_colvar_type
     INTEGER, DIMENSION(3) :: plane
     INTEGER               :: point     
  END TYPE plane_distance_colvar_type

  TYPE angle_colvar_type
     INTEGER, DIMENSION(3) :: i_at_angle
  END TYPE angle_colvar_type

  TYPE rotation_colvar_type
     INTEGER :: i_at1_bond1,&
                i_at2_bond1,&
                i_at1_bond2,&
                i_at2_bond2
  END TYPE rotation_colvar_type

  TYPE ddiff_colvar_type
     INTEGER, DIMENSION(4)         :: i_at_ddiff
  END TYPE ddiff_colvar_type

  TYPE qparm_colvar_type
     INTEGER                         :: l
     INTEGER                         :: n_atoms_to,&
                                        n_atoms_from
     INTEGER, POINTER, DIMENSION (:) :: i_at_from,&
                                        i_at_to
     REAL(KIND=dp)                   :: rcut, alpha
  END TYPE qparm_colvar_type

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
     LOGICAL :: use_points
     REAL(kind=dp)                            :: ss               ! Value of the colvar
     REAL(kind=dp), DIMENSION (:,:), POINTER  :: dsdr             ! Derivative of colvar (3,:)
     INTEGER, DIMENSION (:), POINTER          :: i_atom           ! Mapping of dsdr
     INTEGER                                  :: n_atom_s
     TYPE(dist_colvar_type),POINTER           :: dist_param
     TYPE(coord_colvar_type),POINTER          :: coord_param
     TYPE(torsion_colvar_type),POINTER        :: torsion_param
     TYPE(angle_colvar_type),POINTER          :: angle_param
     TYPE(plane_distance_colvar_type),POINTER :: plane_distance_param
     TYPE(rotation_colvar_type), POINTER      :: rotation_param
     TYPE(ddiff_colvar_type), POINTER         :: ddiff_param
     TYPE(qparm_colvar_type), POINTER         :: qparm_param
     TYPE(point_type), DIMENSION(:), POINTER  :: points
  END TYPE colvar_type

  TYPE point_type
     INTEGER :: type_id
     INTEGER, DIMENSION (:), POINTER          :: atoms
     REAL(KIND=dp), DIMENSION(3)              :: r
  END TYPE point_type

  TYPE colvar_p_type
     TYPE(colvar_type), POINTER :: colvar
  END TYPE colvar_p_type

  TYPE colvar_counters
     INTEGER :: ndist
     INTEGER :: nangle
     INTEGER :: ntorsion
     INTEGER :: ncoord
     INTEGER :: nplane_dist
     INTEGER :: nrot
     INTEGER :: nddiff
     INTEGER :: nqparm
     INTEGER :: nrestraint
     !
     INTEGER :: ntot
  END TYPE colvar_counters
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
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPostcondition(.NOT.ASSOCIATED(colvar),cp_warning_level,routineP,error,failure)
    ALLOCATE(colvar, stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    !
    NULLIFY(colvar%dist_param)
    NULLIFY(colvar%angle_param)
    NULLIFY(colvar%torsion_param)
    NULLIFY(colvar%plane_distance_param)
    NULLIFY(colvar%rotation_param)
    NULLIFY(colvar%ddiff_param)
    NULLIFY(colvar%qparm_param)
    NULLIFY(colvar%coord_param)
    NULLIFY(colvar%i_atom)
    NULLIFY(colvar%dsdr)
    NULLIFY(colvar%points)
    IF (.NOT.failure) THEN
       colvar%type_id    = colvar_id
       colvar%use_points = .FALSE.
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
       CASE(rotation_colvar_id)
          ALLOCATE(colvar%rotation_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(ddiff_colvar_id)
          ALLOCATE(colvar%ddiff_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(qparm_colvar_id)
          ALLOCATE(colvar%qparm_param, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       CASE(no_colvar_id)
          ! Do nothing
       CASE DEFAULT
          CPPostcondition(.FALSE.,cp_warning_level,routineP,error,failure)          
       END SELECT
    END IF

  END SUBROUTINE colvar_create
  !***************************************************************************
  
  !!****f* colvar_types/colvar_check_points *
  !!
  !!   NAME
  !!     colvar_check_points
  !!
  !!   FUNCTION
  !!     Check and setup about the use of geometrical points instead of atoms
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
  !!     Teodoro Laino, [teo] 03.2007
  !!
  !!*** **********************************************************************
  SUBROUTINE  colvar_check_points(colvar, section, error)
    TYPE(colvar_type), POINTER               :: colvar
    TYPE(section_vals_type), POINTER         :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_check_points', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, natoms, npoints, nrep, stat
    INTEGER, DIMENSION(:), POINTER           :: atoms
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: r
    TYPE(section_vals_type), POINTER         :: point_sections

    failure=.FALSE.
    NULLIFY(point_sections)
    CPPostcondition(ASSOCIATED(colvar),cp_warning_level,routineP,error,failure)
    point_sections => section_vals_get_subs_vals(section,"POINT",error=error)
    CALL section_vals_get(point_sections, explicit=explicit, error=error)
    IF (explicit) THEN
       colvar%use_points = .TRUE.
       CALL section_vals_get(point_sections,n_repetition=npoints,error=error)
       ALLOCATE(colvar%points(npoints),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ! Read points definition
       DO i = 1, npoints
          NULLIFY(colvar%points(i)%atoms)
          CALL section_vals_val_get(point_sections,"TYPE",i_rep_section=i,i_val=colvar%points(i)%type_id,error=error)
          CALL section_vals_val_get(point_sections,"ATOMS",i_rep_section=i,n_rep_val=nrep,error=error)
          natoms = 0
          IF (nrep==1) THEN
             CALL section_vals_val_get(point_sections,"ATOMS",i_rep_section=i,i_vals=atoms,error=error)
             natoms = SIZE(atoms)
             ALLOCATE(colvar%points(i)%atoms(natoms),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             colvar%points(i)%atoms = atoms
          END IF
          CALL section_vals_val_get(point_sections,"XYZ",i_rep_section=i,n_rep_val=nrep,error=error)
          IF (nrep==1) THEN
             CALL section_vals_val_get(point_sections,"XYZ",i_rep_section=i,r_vals=r,error=error)
             colvar%points(i)%r = r
             IF (natoms /=0 ) THEN
                CALL stop_program(routineP,"COLVAR section: Specify only ATOMS or XYZ. Not allowed to specify both!")
             END IF
          END IF
       END DO
    END IF
    
  END SUBROUTINE colvar_check_points

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
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, idum, iend, ii, istart, j, &
                                                np, stat
    INTEGER, DIMENSION(:), POINTER           :: list
    LOGICAL                                  :: failure

    CPPostcondition(ASSOCIATED(colvar),cp_warning_level,routineP,error,failure)
    failure=.FALSE.
    IF (.NOT.failure) THEN
       SELECT CASE(colvar%type_id)
       CASE(dist_colvar_id)
          np = 2
          i = colvar%dist_param%i_at
          j = colvar%dist_param%j_at
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s=COLV_SIZE(colvar,i)+&
                          COLV_SIZE(colvar,j)
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          list(1) = colvar%dist_param%i_at
          list(2) = colvar%dist_param%j_at
       CASE(coord_colvar_id)
          np = colvar%coord_param%n_atoms_from+colvar%coord_param%n_atoms_to
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s=0
          DO ii=1,colvar%coord_param%n_atoms_from
             i = colvar%coord_param%i_at_from(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          DO ii=1,colvar%coord_param%n_atoms_to
             i = colvar%coord_param%i_at_to(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1,colvar%coord_param%n_atoms_from
             idum = idum + 1
             i = colvar%coord_param%i_at_from(ii)
             list(idum) = i
          ENDDO
          DO ii=1,colvar%coord_param%n_atoms_to
             idum = idum + 1
             i = colvar%coord_param%i_at_to(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(angle_colvar_id)
          np = 3
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          DO ii=1,3
             i = colvar%angle_param%i_at_angle(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1, 3
             idum = idum + 1
             i = colvar%angle_param%i_at_angle(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(torsion_colvar_id)
          np = 4
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          DO ii=1,4
             i = colvar%torsion_param%i_at_tors(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1, 4
             idum = idum + 1
             i = colvar%torsion_param%i_at_tors(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(plane_distance_colvar_id)
          np = 4
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          DO ii=1,3
             i = colvar%plane_distance_param%plane(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          i = colvar%plane_distance_param%point
          colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1, 3
             idum = idum + 1
             i = colvar%plane_distance_param%plane(ii)
             list(idum) = i
          ENDDO
          i = colvar%plane_distance_param%point
          list(4) = i
          idum = idum + 1
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)            
       CASE(ddiff_colvar_id)
          np = 4
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          DO ii=1,4
             i = colvar%ddiff_param%i_at_ddiff(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1, 4
             idum = idum + 1
             i = colvar%ddiff_param%i_at_ddiff(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       CASE(rotation_colvar_id)
          np = 4
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s= 0
          i = colvar%rotation_param%i_at1_bond1
          colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          i = colvar%rotation_param%i_at2_bond1
          colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          i = colvar%rotation_param%i_at1_bond2
          colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          i = colvar%rotation_param%i_at2_bond2
          colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          i = colvar%rotation_param%i_at1_bond1
          list(1) = i
          i = colvar%rotation_param%i_at2_bond1
          list(2) = i
          i = colvar%rotation_param%i_at1_bond2
          list(3) = i
          i = colvar%rotation_param%i_at2_bond2
          list(4) = i
       CASE(qparm_colvar_id)
          np = colvar%qparm_param%n_atoms_from+colvar%qparm_param%n_atoms_to
          ! Number of real atoms involved in the colvar
          colvar%n_atom_s=0
          DO ii=1,colvar%qparm_param%n_atoms_from
             i = colvar%qparm_param%i_at_from(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          DO ii=1,colvar%qparm_param%n_atoms_to
             i = colvar%qparm_param%i_at_to(ii)
             colvar%n_atom_s=colvar%n_atom_s+COLV_SIZE(colvar,i)
          ENDDO
          ! Create a List of points...
          ALLOCATE(list(np),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          idum = 0
          DO ii=1,colvar%qparm_param%n_atoms_from
             idum = idum + 1
             i = colvar%qparm_param%i_at_from(ii)
             list(idum) = i
          ENDDO
          DO ii=1,colvar%qparm_param%n_atoms_to
             idum = idum + 1
             i = colvar%qparm_param%i_at_to(ii)
             list(idum) = i
          ENDDO
          CPPostcondition(idum==np,cp_failure_level,routinep,error,failure)
       END SELECT
       
       ALLOCATE(colvar%dsdr(3,colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(colvar%i_atom(colvar%n_atom_s), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ! And now map real atoms
       istart = 0
       iend   = 0
       DO i = 1, SIZE(list)
          IF (.NOT.colvar%use_points) THEN
             ! No point centers
             colvar%i_atom(i) = list(i)
             iend = iend + 1
          ELSE
             IF (ASSOCIATED(colvar%points(list(i))%atoms)) THEN
                iend = istart + SIZE(colvar%points(list(i))%atoms)
                colvar%i_atom(istart+1:iend) = colvar%points(list(i))%atoms
                istart = iend
             END IF
          END IF
       END DO
       CPPostcondition(iend==colvar%n_atom_s,cp_failure_level,routinep,error,failure)
       DEALLOCATE(list, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)

    END IF    
  END SUBROUTINE colvar_setup

  !!****f* colvar_types/colv_size *
  !!
  !!   NAME
  !!     colv_size
  !!
  !!   FUNCTION
  !!     Gives back the size of an array of integer. If not associated gives back 1
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!
  !!   AUTHOR
  !!     Teodoro Laino - 03.2007
  !!
  !!*** **********************************************************************
  FUNCTION colv_size(colvar,i) RESULT(my_size)
    TYPE(colvar_type), POINTER               :: colvar
    INTEGER                                  :: i, my_size

    my_size = 1
    IF (ASSOCIATED(colvar%points)) THEN
       IF (ASSOCIATED(colvar%points(i)%atoms)) THEN
          my_size = SIZE(colvar%points(i)%atoms)
       ELSE
          my_size = 0
       END IF
    END IF
  END FUNCTION colv_size

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
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
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
    IF (ASSOCIATED(colvar%points)) THEN
       DO i = 1, SIZE(colvar%points)
          IF (ASSOCIATED(colvar%points(i)%atoms)) THEN
             DEALLOCATE(colvar%points(i)%atoms, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          END IF
       END DO
       DEALLOCATE(colvar%points, stat=stat)
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
    CASE(ddiff_colvar_id)
       DEALLOCATE(colvar%ddiff_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (rotation_colvar_id)
       DEALLOCATE(colvar%rotation_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE (qparm_colvar_id)
       DEALLOCATE(colvar%qparm_param%i_at_from, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%qparm_param%i_at_to,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       DEALLOCATE(colvar%qparm_param, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CASE(no_colvar_id)
       ! Do nothing
    CASE default
       CPAssert(.FALSE.,cp_failure_level,routinep,error,failure)
    END SELECT
    DEALLOCATE(colvar, stat=stat)
    CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)

  END SUBROUTINE colvar_release

  !!****f* colvar_types/colvar_clone *
  !!
  !!   NAME
  !!     colvar_clone
  !!
  !!   FUNCTION
  !!     Clone a colvar type
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
  !!     Teodoro Laino [tlaino] 04.2006
  !!
  !!*** **********************************************************************
  SUBROUTINE colvar_clone(colvar_out, colvar_in, i_atom_offset, error)
    TYPE(colvar_type), POINTER               :: colvar_out, colvar_in
    INTEGER, INTENT(IN), OPTIONAL            :: i_atom_offset
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_clone', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: my_offset, ndim, stat
    LOGICAL                                  :: failure

    my_offset = 0   
    failure=.FALSE.
    CPPostcondition(ASSOCIATED(colvar_in),cp_warning_level,routineP,error,failure)
    CPPostcondition(.NOT.ASSOCIATED(colvar_out),cp_warning_level,routineP,error,failure)
    IF (PRESENT(i_atom_offset)) my_offset = i_atom_offset
    CALL colvar_create(colvar_out, colvar_in%type_id, error)
    CALL colvar_clone_points(colvar_out, colvar_in, my_offset, error)
    IF (colvar_in%use_points) my_offset = 0
    SELECT CASE(colvar_out%type_id)
    CASE(dist_colvar_id)
       colvar_out%dist_param%i_at = colvar_in%dist_param%i_at + my_offset
       colvar_out%dist_param%j_at = colvar_in%dist_param%j_at + my_offset
    CASE(coord_colvar_id)
       colvar_out%coord_param%n_atoms_to   = colvar_in%coord_param%n_atoms_to  
       colvar_out%coord_param%n_atoms_from = colvar_in%coord_param%n_atoms_from
       colvar_out%coord_param%nncrd        = colvar_in%coord_param%nncrd       
       colvar_out%coord_param%ndcrd        = colvar_in%coord_param%ndcrd       
       colvar_out%coord_param%r_0          = colvar_in%coord_param%r_0         
       ndim = SIZE(colvar_out%coord_param%i_at_from)
       ALLOCATE(colvar_out%coord_param%i_at_from(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ndim = SIZE(colvar_out%coord_param%i_at_to)
       ALLOCATE(colvar_out%coord_param%i_at_from(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       colvar_out%coord_param%i_at_from    = colvar_in%coord_param%i_at_from + my_offset
       colvar_out%coord_param%i_at_to      = colvar_in%coord_param%i_at_from + my_offset
    CASE(angle_colvar_id)
       colvar_out%angle_param%i_at_angle   = colvar_in%angle_param%i_at_angle + my_offset
    CASE(torsion_colvar_id)
       colvar_out%torsion_param%i_at_tors     = colvar_in%torsion_param%i_at_tors + my_offset
       colvar_out%torsion_param%o0            = colvar_in%torsion_param%o0
    CASE(plane_distance_colvar_id)
       colvar_out%plane_distance_param%plane = colvar_in%plane_distance_param%plane + my_offset
       colvar_out%plane_distance_param%point = colvar_in%plane_distance_param%point + my_offset
    CASE(rotation_colvar_id)
       colvar_out%rotation_param%i_at1_bond1  = colvar_in%rotation_param%i_at1_bond1 + my_offset
       colvar_out%rotation_param%i_at2_bond1  = colvar_in%rotation_param%i_at2_bond1 + my_offset
       colvar_out%rotation_param%i_at1_bond2  = colvar_in%rotation_param%i_at1_bond2 + my_offset
       colvar_out%rotation_param%i_at2_bond2  = colvar_in%rotation_param%i_at2_bond2 + my_offset
    CASE(ddiff_colvar_id)
       colvar_out%ddiff_param%i_at_ddiff   = colvar_in%ddiff_param%i_at_ddiff + my_offset
    CASE(qparm_colvar_id)
       colvar_out%qparm_param%n_atoms_to   = colvar_in%qparm_param%n_atoms_to  
       colvar_out%qparm_param%n_atoms_from = colvar_in%qparm_param%n_atoms_from
       colvar_out%qparm_param%rcut         = colvar_in%qparm_param%rcut
       colvar_out%qparm_param%l            = colvar_in%qparm_param%l
       colvar_out%qparm_param%alpha        = colvar_in%qparm_param%alpha
       ndim = SIZE(colvar_out%qparm_param%i_at_from)
       ALLOCATE(colvar_out%qparm_param%i_at_from(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ndim = SIZE(colvar_out%qparm_param%i_at_to)
       ALLOCATE(colvar_out%qparm_param%i_at_from(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       colvar_out%qparm_param%i_at_from    = colvar_in%qparm_param%i_at_from + my_offset
       colvar_out%qparm_param%i_at_to      = colvar_in%qparm_param%i_at_from + my_offset
    END SELECT
    CALL colvar_setup(colvar_out, error) 
  END SUBROUTINE colvar_clone

  !!****f* colvar_types/colvar_clone_points *
  !!
  !!   NAME
  !!     colvar_clone_points
  !!
  !!   FUNCTION
  !!     Clone points type of a colvar type
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
  !!     Teodoro Laino [tlaino] 03.2007
  !!
  !!*** **********************************************************************
  SUBROUTINE colvar_clone_points(colvar_out, colvar_in, offset, error)
    TYPE(colvar_type), POINTER               :: colvar_out, colvar_in
    INTEGER, INTENT(IN)                      :: offset
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_clone_points', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, natoms, npoints, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    colvar_out%use_points = colvar_in%use_points
    IF (colvar_in%use_points) THEN
       CPPostcondition(ASSOCIATED(colvar_in%points),cp_failure_level,routineP,error,failure)
       npoints = SIZE(colvar_in%points)
       ALLOCATE(colvar_out%points(npoints),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       DO i = 1, npoints
          IF (ASSOCIATED(colvar_in%points(i)%atoms)) THEN
             natoms = SIZE(colvar_in%points(i)%atoms)
             ALLOCATE(colvar_out%points(i)%atoms(natoms),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             colvar_out%points(i)%atoms = colvar_in%points(i)%atoms+offset
          ELSE
             NULLIFY(colvar_out%points(i)%atoms)
          END IF
          colvar_out%points(i)%type_id = colvar_in%points(i)%type_id
          colvar_out%points(i)%r       = colvar_in%points(i)%r
       END DO
    ELSE
       NULLIFY(colvar_out%points)
    END IF
    
  END SUBROUTINE colvar_clone_points

  !!****f* colvar_types/colvar_p_reallocate *
  !!
  !!   NAME
  !!     colvar_p_reallocate
  !!
  !!   FUNCTION
  !!     Change the dimension of a colvar_p_type
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
  !!     Teodoro Laino [tlaino] 04.2006
  !!
  !!*** **********************************************************************

  SUBROUTINE colvar_p_reallocate(colvar_set,lb1_new,ub1_new, error)
    TYPE(colvar_p_type), DIMENSION(:), &
      POINTER                                :: colvar_set
    INTEGER, INTENT(IN)                      :: lb1_new, ub1_new
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'colvar_p_reallocate'

    INTEGER                                  :: istat, j, lb1, lb1_old, ub1, &
                                                ub1_old
    TYPE(colvar_p_type), DIMENSION(:), &
      POINTER                                :: work

    NULLIFY(work)
    IF (ASSOCIATED(colvar_set)) THEN
       lb1_old = LBOUND(colvar_set,1)
       ub1_old = UBOUND(colvar_set,1)
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       ALLOCATE (work(lb1:ub1),STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"work")
       END IF
       DO j = lb1, ub1
          NULLIFY(work(j)%colvar)
          CALL colvar_clone(work(j)%colvar, colvar_set(j)%colvar, error=error)
       END DO
       DO j = lb1, ub1
          CALL colvar_release(colvar_set(j)%colvar, error)
          NULLIFY(colvar_set(j)%colvar)
       END DO
       DEALLOCATE(colvar_set, STAT=istat)
       IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"colvar_set")
    END IF

    ALLOCATE (colvar_set(lb1_new:ub1_new),STAT=istat)
    IF (istat /= 0) THEN
       CALL stop_memory(routineN,moduleN,__LINE__,"colvar_set")
    END IF
    DO j = lb1_new, ub1_new
       NULLIFY(colvar_set(j)%colvar)
    END DO
    IF (ASSOCIATED(work)) THEN
       lb1 = MAX(lb1_new,lb1_old)
       ub1 = MIN(ub1_new,ub1_old)
       DO j = lb1, ub1
          CALL colvar_clone(colvar_set(j)%colvar, work(j)%colvar, error=error)
       END DO
       DO j = lb1, ub1
          CALL colvar_release(work(j)%colvar, error)
          NULLIFY(work(j)%colvar)
       END DO
       DEALLOCATE (work,STAT=istat)
       IF (istat /= 0) THEN
          CALL stop_memory(routineN,moduleN,__LINE__,"work")
       END IF
    END IF
  END SUBROUTINE colvar_p_reallocate

  !!****f* colvar_types/eval_point_pos *
  !!
  !!   NAME
  !!     eval_point_pos
  !!
  !!   FUNCTION
  !!     Evaluate the position of the geometrical point
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!
  !!   AUTHOR
  !!     Teodoro Laino - 03.2007
  !!
  !!*** **********************************************************************
  SUBROUTINE eval_point_pos(point, particles, r)
    TYPE(point_type)                         :: point
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: r

    INTEGER                                  :: i

    SELECT CASE(point%type_id)
    CASE(do_clv_geo_center)
       r = 0.0_dp
       DO i = 1,SIZE(point%atoms)
          r = r + particles(point%atoms(i))%r
       END DO
       r = r/REAL(SIZE(point%atoms),KIND=dp)
    CASE(do_clv_fix_point )
       r = point%r
    END SELECT

  END SUBROUTINE eval_point_pos

  !!****f* colvar_types/eval_point_pos *
  !!
  !!   NAME
  !!     eval_point_pos
  !!
  !!   FUNCTION
  !!     Evaluate the position of the geometrical point
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!
  !!   AUTHOR
  !!     Teodoro Laino - 03.2007
  !!
  !!*** **********************************************************************
  SUBROUTINE eval_point_der(points, i, dsdr, f)
    TYPE(point_type), DIMENSION(:), POINTER  :: points
    INTEGER, INTENT(IN)                      :: i
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: dsdr
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: f

    INTEGER                                  :: ind, j
    REAL(KIND=dp)                            :: fac

    SELECT CASE(points(i)%type_id)
    CASE(do_clv_geo_center)
       fac = REAL(SIZE(points(i)%atoms),KIND=dp)
       ind = 0
       DO j = 1, i-1
          IF (ASSOCIATED(points(j)%atoms)) THEN
             ind = ind + SIZE(points(j)%atoms)
          END IF
       END DO
       DO j = 1,SIZE(points(i)%atoms)
          dsdr(:,ind+j) = dsdr(:,ind+j) + f / fac 
       END DO
    CASE(do_clv_fix_point )
       ! Do nothing if it's a fixed point in space
    END SELECT

  END SUBROUTINE eval_point_der

END MODULE colvar_types
