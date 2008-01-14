!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/dimer_methods [1.0] *
!!
!!   NAME
!!     dimer_methods
!!
!!   FUNCTION
!!     Contains types used for a Dimer Method calculations
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!*****************************************************************************
MODULE dimer_methods
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dimer_types,                     ONLY: dimer_env_type
  USE f77_blas
  USE f77_interface,                   ONLY: f_env_add_defaults,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE gopt_f_types,                    ONLY: gopt_f_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE replica_types,                   ONLY: rep_env_calc_e_f,&
                                             rep_env_create,&
                                             rep_env_release,&
                                             replica_env_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dimer_methods'

  PUBLIC :: cp_eval_at_ts

CONTAINS
  
!!****f* dimer_methods/dimer_env_create [1.0] *
!!
!!   NAME
!!     dimer_env_create
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     -
!!
!!   AUTHOR
!!     Luca Bellucci and Teodoro Laino - created [tlaino] - 01.2008
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  SUBROUTINE  cp_eval_at_ts (gopt_env, x, f, gradient, calc_force, error)
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: x
    REAL(KIND=dp), INTENT(OUT)               :: f
    REAL(KIND=dp), DIMENSION(:), POINTER     :: gradient
    LOGICAL, INTENT(IN)                      :: calc_force
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_eval_at_ts', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ierr, j, nrep, &
                                                prep, stat
    LOGICAL                                  :: failure, use_replica
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: dgrad
    REAL(KIND=dp), DIMENSION(2)              :: energy_loc
    TYPE(cp_error_type)                      :: new_error
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dimer_env_type), POINTER            :: dimer_env
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(replica_env_type), POINTER          :: rep_env
    TYPE(section_vals_type), POINTER         :: dimer_section, root_section

    failure=.FALSE.
    NULLIFY(para_env,dimer_env,f_env,rep_env,dimer_section, root_section)
    CALL timeset(routineN,"I","",handle)
    CPPostcondition(ASSOCIATED(gopt_env),cp_failure_level,routineP,error,failure)
    dimer_env => gopt_env%dimer_env
    para_env  => gopt_env%force_env%para_env
    root_section => gopt_env%force_env%root_section
    CPPostcondition(ASSOCIATED(dimer_env),cp_failure_level,routineP,error,failure)
    ALLOCATE(dgrad(SIZE(x)), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    dimer_section => section_vals_get_subs_vals(gopt_env%motion_section,&
         "GEO_OPT%TRANSITION_STATE%DIMER",error=error)
    CALL section_vals_val_get(dimer_section,"USE_REPLICA",l_val=use_replica,error=error)
    nrep = 2
    IF (.NOT.use_replica) nrep = 1
    prep = MAX(1,para_env%num_pe/nrep)
    ! Create Replica Environments
    CALL rep_env_create(rep_env, para_env=para_env, input=root_section, nrep=nrep,&
         prep=prep, error=error)
    IF (ASSOCIATED(rep_env)) THEN
       CALL f_env_add_defaults(f_env_id=rep_env%f_env_id,f_env=f_env,&
            new_error=new_error, failure=failure)
       gradient = 0.0_dp
       f = 0.0_dp
       DO i = 1, 2, nrep
          DO j = 1, nrep
             CPPostcondition(SIZE(rep_env%r,1)==SIZE(x),cp_failure_level,routineP,error,failure)
             rep_env%r(:,j) = x + (-1)**(j+i) * gopt_env%dimer_env%nvec * gopt_env%dimer_env%dr
          END DO

          ! Computes energy and forces with replica
          CALL rep_env_calc_e_f(rep_env,calc_f=calc_force,error=error)

          DO j = 1, nrep
             dgrad     = dgrad + (-1)**(j+i) * rep_env%f(1:SIZE(x),j)
             gradient  = gradient + rep_env%f(1:SIZE(x),j)
             energy_loc(j+i-1) = rep_env%f(SIZE(x)+1,j)
          END DO
       END DO
       gradient = - gradient * 0.5_dp
       f = SUM(energy_loc)*0.5_dp + gopt_env%dimer_env%dr*0.25_dp*DOT_PRODUCT(dgrad,gopt_env%dimer_env%nvec)

       CALL f_env_rm_defaults(f_env,new_error,ierr)
    END IF
    CALL rep_env_release(rep_env,error=error)
    DEALLOCATE(dgrad, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE cp_eval_at_ts
!***************************************************************************

END MODULE dimer_methods
