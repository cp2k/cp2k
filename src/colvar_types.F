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
  USE termination,                     ONLY: stop_memory
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
            rotation_colvar_type

  INTEGER, PARAMETER, PUBLIC :: no_colvar_id              =-2,&
                                dist_colvar_id            =1,&
                                coord_colvar_id           =2,&
                                torsion_colvar_id         =3,&
                                angle_colvar_id           =4,&
                                plane_distance_colvar_id  =5,&
                                rotation_colvar_id        =6

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

  TYPE rotation_colvar_type
     INTEGER, DIMENSION(:),POINTER :: i_at1_bond1,&
                                      i_at2_bond1,&
                                      i_at1_bond2,&
                                      i_at2_bond2
  END TYPE rotation_colvar_type

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
     REAL(kind=dp), DIMENSION (:,:), POINTER  :: dsdr             ! Derivative of colvar (3,:)
     INTEGER, DIMENSION (:), POINTER          :: i_atom           ! Mapping of dsdr
     INTEGER                                  :: n_atom_s
     TYPE(dist_colvar_type),POINTER           :: dist_param
     TYPE(coord_colvar_type),POINTER          :: coord_param
     TYPE(torsion_colvar_type),POINTER        :: torsion_param
     TYPE(angle_colvar_type),POINTER          :: angle_param
     TYPE(plane_distance_colvar_type),POINTER :: plane_distance_param
     TYPE(rotation_colvar_type), POINTER      :: rotation_param
  END TYPE colvar_type

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
       CASE(rotation_colvar_id)
          ALLOCATE(colvar%rotation_param, stat=stat)
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

    INTEGER                                  :: ii, stat, idum
    LOGICAL                                  :: failure

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
       CASE(rotation_colvar_id)
          colvar%n_atom_s=SIZE(colvar%rotation_param%i_at1_bond1)+&
                          SIZE(colvar%rotation_param%i_at2_bond1)+&
                          SIZE(colvar%rotation_param%i_at1_bond2)+&
                          SIZE(colvar%rotation_param%i_at2_bond2)
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
       CASE(rotation_colvar_id)
          idum = 0
          DO ii=1, SIZE(colvar%rotation_param%i_at1_bond1)
             idum = idum+1
             colvar%i_atom(idum) = colvar%rotation_param%i_at1_bond1(ii)
          END DO
          DO ii=1, SIZE(colvar%rotation_param%i_at2_bond1)
             idum = idum+1
             colvar%i_atom(idum) = colvar%rotation_param%i_at2_bond1(ii)
          END DO
          DO ii=1, SIZE(colvar%rotation_param%i_at1_bond2)
             idum = idum+1
             colvar%i_atom(idum) = colvar%rotation_param%i_at1_bond2(ii)
          END DO
          DO ii=1, SIZE(colvar%rotation_param%i_at2_bond2)
             idum = idum+1
             colvar%i_atom(idum) = colvar%rotation_param%i_at2_bond2(ii)
          END DO
          CPPostcondition(idum==colvar%n_atom_s,cp_failure_level,routinep,error,failure)
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
    CASE (rotation_colvar_id)
       DEALLOCATE(colvar%rotation_param%i_at1_bond1, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%rotation_param%i_at2_bond1, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%rotation_param%i_at1_bond2, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%rotation_param%i_at2_bond2, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(colvar%rotation_param, stat=stat)
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
    IMPLICIT NONE
    ! Arguments
    TYPE(colvar_type), POINTER :: colvar_out, colvar_in
    INTEGER, INTENT(IN), OPTIONAL :: i_atom_offset
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'colvar_clone', &
         routineP = moduleN//':'//routineN
    ! Local Variables
    INTEGER                                  :: stat,ndim, my_offset
    LOGICAL                                  :: failure
    
    my_offset = 0   
    failure=.FALSE.
    CPPostcondition(ASSOCIATED(colvar_in),cp_warning_level,routineP,error,failure)
    CPPostcondition(.NOT.ASSOCIATED(colvar_out),cp_warning_level,routineP,error,failure)
    IF (PRESENT(i_atom_offset)) my_offset = i_atom_offset
    CALL colvar_create(colvar_out, colvar_in%type_id, error)
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
    CASE(plane_distance_colvar_id)
       colvar_out%plane_distance_param%plane = colvar_in%plane_distance_param%plane + my_offset
       colvar_out%plane_distance_param%point = colvar_in%plane_distance_param%point + my_offset
    CASE(rotation_colvar_id)
       ndim = SIZE(colvar_out%rotation_param%i_at1_bond1)
       ALLOCATE(colvar_out%rotation_param%i_at1_bond1(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ndim = SIZE(colvar_out%rotation_param%i_at2_bond1)
       ALLOCATE(colvar_out%rotation_param%i_at2_bond1(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ndim = SIZE(colvar_out%rotation_param%i_at1_bond2)
       ALLOCATE(colvar_out%rotation_param%i_at1_bond2(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       ndim = SIZE(colvar_out%rotation_param%i_at2_bond2)
       ALLOCATE(colvar_out%rotation_param%i_at2_bond2(ndim),stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       colvar_out%rotation_param%i_at1_bond1  = colvar_in%rotation_param%i_at1_bond1 + my_offset
       colvar_out%rotation_param%i_at2_bond1  = colvar_in%rotation_param%i_at2_bond1 + my_offset
       colvar_out%rotation_param%i_at1_bond2  = colvar_in%rotation_param%i_at1_bond2 + my_offset
       colvar_out%rotation_param%i_at2_bond2  = colvar_in%rotation_param%i_at2_bond2 + my_offset
    END SELECT    
    CALL colvar_setup(colvar_out, error) 
  END SUBROUTINE colvar_clone

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
    ! Arguments
    TYPE(colvar_p_type), POINTER, DIMENSION(:) :: colvar_set
    INTEGER, INTENT(IN)                        :: lb1_new, ub1_new
    TYPE(cp_error_type), INTENT(inout), &
         OPTIONAL                              :: error

    CHARACTER(LEN=*), PARAMETER                :: routineN = 'colvar_p_reallocate'
    INTEGER                                    :: istat, lb1, lb1_old, &
                                                  ub1, ub1_old, j
    TYPE(colvar_p_type), POINTER, DIMENSION(:) :: work

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

END MODULE colvar_types
