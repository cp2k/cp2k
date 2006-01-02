MODULE metadynamics_types
  USE colvar_types,                    ONLY: colvar_type, coord_colvar_type, &
                                             dist_colvar_type,colvar_init,&
                                             colvar_dealloc_ref,torsion_colvar_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_para_env,                     ONLY: cp_para_env_release
  USE kinds,                           ONLY: default_path_length,&
                                             dp
#include "cp_common_uses.h"
  IMPLICIT NONE

  PRIVATE 

  CHARACTER(len=*), PARAMETER, PRIVATE :: modulen='metadynamics_types'

  PUBLIC :: hills_env_type, meta_env_type, meta_env_retain, meta_env_release

  TYPE hills_env_type
     REAL(kind=dp) :: ww
     INTEGER :: n_hills, nt_hills
     CHARACTER(len=default_path_length) :: hills_file_name
     REAL(kind=dp), DIMENSION(:,:), POINTER :: ss_history
     LOGICAL :: restart_hills
  END TYPE hills_env_type

  TYPE meta_env_type
     INTEGER :: ref_count, id_nr
     LOGICAL :: do_hills
     LOGICAL :: extended_lagrange
     INTEGER :: n_colvar
     REAL(kind=dp) :: ekin_s,cv_energy, dt
     LOGICAL :: tempcontrol
     REAL(kind=dp) :: temp_wanted,toll_temp
     INTEGER :: n_steps, colvar_file_unit_nr
     CHARACTER(len=default_path_length) :: colvar_file_name     
     TYPE(hills_env_type) :: hills_env
     TYPE(colvar_type), DIMENSION(:), POINTER :: colvars
     TYPE(cp_para_env_type), POINTER :: para_env
  END TYPE meta_env_type

CONTAINS
  
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
        CALL cp_para_env_release(meta_env%para_env)
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
