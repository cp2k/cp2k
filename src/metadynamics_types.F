#include "cp_prep_globals.h"

MODULE metadynamics_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp

  PRIVATE 

  CHARACTER(len=*), PARAMETER, PRIVATE :: modulen='metadynamics_types'

   PUBLIC :: colvar_type, dist_colvar_type,colvar_init, colvar_dealloc_ref
   PUBLIC :: hills_env_type, meta_env_type, meta_env_retain, meta_env_release

  INTEGER, PARAMETER, PUBLIC :: no_colvar_id=-2,dist_colvar_id=1

!!****s* metadynamics/dist_colvar_param [1.0] *
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
!!     - i_subsys,j_subsys: subsystem to which the indexes refer
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!   SOURCE
  !***************************************************************************
  TYPE dist_colvar_type
     INTEGER :: i_at,i_subsys,j_at,j_subsys
  END TYPE dist_colvar_type
!!***
  !****************************************************************************

!!****s* metadynamics/dist_colvar_param [1.0] *
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
!!     - index: index of this the colvar
!!     - type_id: the type of this colvar
!!     - mass_s mass of the collective variable
!!     - lambda: lambda of the gaussian in the direction of this colvar
!!     - dist_param: parameters for the distance collective variable
!!       (meaningful only if type_id==dist_colvar_id)
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!   SOURCE
  !***************************************************************************
  TYPE colvar_type
     INTEGER :: index,type_id
     REAL(kind=dp) :: mass,lambda,ss, ss0, ff_s,vv, vvp, &
                      cv_energy,delta_s,ff_hills
     TYPE(dist_colvar_type) :: dist_param
  END TYPE colvar_type
!!***
  !****************************************************************************

   TYPE hills_env_type
     REAL(kind=dp) :: ww
     INTEGER :: n_hills, nt_hills
     REAL(kind=dp), DIMENSION(:,:), POINTER :: ss_history
     LOGICAL :: restart_hills
  END TYPE hills_env_type

  TYPE meta_env_type
     INTEGER :: ref_count, id_nr
     LOGICAL :: do_hills
     INTEGER :: n_colvar
     REAL(kind=dp) :: ekin_s,cv_energy, dt
     INTEGER :: n_steps
     TYPE(hills_env_type) :: hills_env
     TYPE(colvar_type), DIMENSION(:), POINTER :: colvars
  END TYPE meta_env_type

CONTAINS
  

!!****f* metadynamics/colvar_init *
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
SUBROUTINE colvar_init(colvar,index,error)
    TYPE(colvar_type), INTENT(inout)         :: colvar
    INTEGER, INTENT(in)                      :: index
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routinen = 'colvar_init', &
      routinep = modulen//':'//routinen

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  colvar%index=index
  colvar%type_id=no_colvar_id
  colvar%mass=-HUGE(0.0_dp)
  colvar%lambda=-HUGE(0.0_dp)
  colvar%dist_param%i_at=-HUGE(0)
  colvar%dist_param%j_at=-HUGE(0)
  colvar%dist_param%i_subsys=-HUGE(0)
  colvar%dist_param%j_subsys=-HUGE(0)

  colvar%ss=0.0_dp
  colvar%ss0=0.0_dp
  colvar%ff_s=0.0_dp
  colvar%vv=0.0_dp
  colvar%vvp=0.0_dp
END SUBROUTINE colvar_init
!***************************************************************************

!!****f* metadynamics/colvar_dealloc_ref *
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

    CHARACTER(len=*), PARAMETER :: routinen = 'colvar_dealloc_ref', &
      routinep = modulen//':'//routinen

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  SELECT CASE(colvar%type_id)
  CASE (dist_colvar_id)
  CASE(no_colvar_id)
  CASE default
     CPAssert(.FALSE.,cp_failure_level,routinep,error,failure)
  END SELECT
  colvar%dist_param%i_at=-HUGE(0)
  colvar%dist_param%j_at=-HUGE(0)
  colvar%index=-HUGE(0)
  colvar%type_id=no_colvar_id
  colvar%mass=-HUGE(0.0_dp)
  colvar%lambda=-HUGE(0.0_dp)
  colvar%delta_s=-HUGE(0.0_dp)
END SUBROUTINE colvar_dealloc_ref
!***************************************************************************

!!****f* metadynamics/meta_env_retain *
!!
!!   NAME
!!     meta_env_retain
!!
!!   FUNCTION
!!     retains the meta_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!*** **********************************************************************
SUBROUTINE meta_env_retain(meta_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routinen = 'meta_env_retain', &
      routinep = modulen//':'//routinen

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPreconditionNoFail(ASSOCIATED(meta_env),cp_failure_level,routineP,error)
  CPPreconditionNoFail(meta_env%ref_count>0,cp_failure_level,routineP,error)
  meta_env%ref_count=meta_env%ref_count+1
END SUBROUTINE meta_env_retain
!***************************************************************************

!!****f* metadynamics/meta_env_release *
!!
!!   NAME
!!     meta_env_release
!!
!!   FUNCTION
!!     releases the meta_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!
!!*** **********************************************************************
SUBROUTINE meta_env_release(meta_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routinen = 'meta_env_release', &
      routinep = modulen//':'//routinen

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(meta_env)) THEN
     CPPreconditionNoFail(meta_env%ref_count>0,cp_failure_level,routineP,error)
     meta_env%ref_count=meta_env%ref_count-1
     IF (meta_env%ref_count==0) THEN
        IF (ASSOCIATED(meta_env%colvars)) THEN
           DO i=1,SIZE(meta_env%colvars)
              CALL colvar_dealloc_ref(meta_env%colvars(i),error=error)
           END DO
           DEALLOCATE(meta_env%colvars,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        IF (ASSOCIATED(meta_env%hills_env%ss_history)) THEN
           DEALLOCATE(meta_env%hills_env%ss_history,stat=stat)
           CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
        END IF
        DEALLOCATE(meta_env, stat=stat)
        CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
     END IF
  END IF
  NULLIFY(meta_env)
END SUBROUTINE meta_env_release
!***************************************************************************

END MODULE metadynamics_types
