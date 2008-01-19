!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****f* gopt_f77_methods/cp_eval_at [1.0] *
!!
!!   NAME
!!     cp_eval_at
!!
!!   FUNCTION
!!     evaluete the potential energy and its gradients using an array
!!     with same dimension as the particle_set
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     gopt_env: the geometry optimization environment
!!     x: the position where the function should be evaluated
!!     f: the function value
!!     gradient: the value of its gradient
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino [tlaino] - University of Zurich - 01.2008
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************
  RECURSIVE SUBROUTINE cp_eval_at(gopt_env, x, f, gradient, master, para_env, error)
    USE cell_types,                      ONLY: cell_type
    USE cp_para_types,                   ONLY: cp_para_env_type
    USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                               cp_subsys_type
    USE dimer_methods,                   ONLY: cp_eval_at_ts
    USE f77_blas
    USE force_env_methods,               ONLY: force_env_calc_energy_force
    USE force_env_types,                 ONLY: force_env_get,&
                                               force_env_type
    USE gopt_f_types,                    ONLY: gopt_f_type
    USE gopt_param_types,                ONLY: gopt_param_type
    USE input_constants,                 ONLY: default_minimization_method_id,&
                                               default_ts_method_id
    USE kinds,                           ONLY: dp
    USE particle_list_types,             ONLY: particle_list_type
    USE particle_types,                  ONLY: write_structure_data
    USE timings,                         ONLY: timeset,&
                                               timestop
#include "cp_common_uses.h"
    TYPE(gopt_f_type), POINTER               :: gopt_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: x
    REAL(KIND=dp), INTENT(out), OPTIONAL     :: f
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: gradient
    INTEGER, INTENT(in)                      :: master
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_eval_at', moduleN='gopt_f77_methods',&
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, idg, idir, ip, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: f_ts
    REAL(KIND=dp), DIMENSION(:), POINTER     :: gradient_ts
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(particle_list_type), POINTER        :: particles

    failure=.FALSE.
    NULLIFY(particles, subsys, cell, gradient_ts)

    CALL timeset(routineN,"I","",handle)
    CALL force_env_get(gopt_env%force_env, subsys=subsys, cell=cell,error=error)
    idg=0
    CALL cp_subsys_get(subsys, particles=particles, error=error)
    DO ip=1,particles%n_els
       DO idir=1,3
          idg=idg+1
          particles%els(ip)%r(idir)=x(idg)
       END DO
    END DO
    CALL write_structure_data(particles%els,cell,gopt_env%motion_section,error)

    SELECT CASE(gopt_env%type_id)
    CASE(default_minimization_method_id)
       ! Geometry Minimization
       CALL force_env_calc_energy_force(gopt_env%force_env,&
            calc_force=PRESENT(gradient),error=error)

       ! Possibly take the potential energy
       IF (PRESENT(f)) THEN
          CALL force_env_get(gopt_env%force_env,potential_energy=f,error=error)
       END IF
       
       ! Possibly take the gradients
       IF (PRESENT(gradient)) THEN
          IF (master==para_env%mepos) THEN ! we are on the master
             CPPrecondition(ASSOCIATED(gradient),cp_failure_level,routineP,error,failure)
             CALL force_env_get(gopt_env%force_env, subsys=subsys,error=error)
             idg=0
             CALL cp_subsys_get(subsys, particles=particles, error=error)
             DO ip=1,particles%n_els
                DO idir=1,3
                   idg=idg+1
                   CPInvariant(SIZE(gradient)>=idg,cp_failure_level,routineP,error,failure)
                   gradient(idg)=-particles%els(ip)%f(idir)
                END DO
             END DO
          END IF
       END IF
    CASE(default_ts_method_id)
       ! Transition State Optimization
       ALLOCATE(gradient_ts(particles%n_els*3),stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

       ! Real calculation of energy and forces for transition state optimization
       CALL cp_eval_at_ts (gopt_env, x, f_ts, gradient_ts, calc_force=PRESENT(gradient), error=error)

       ! Possibly take the potential energy
       IF (PRESENT(f)) f = f_ts
       
       ! Possibly take the gradients
       IF (PRESENT(gradient)) THEN
          IF (master==para_env%mepos) THEN ! we are on the master
             CPPrecondition(ASSOCIATED(gradient),cp_failure_level,routineP,error,failure)
             gradient = gradient_ts
          END IF
       END IF
       DEALLOCATE(gradient_ts, stat=stat)
       CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    END SELECT

    ! This call is necessary for QM/MM if a Translation is applied
    ! this makes the geometry optimizer consistent
    idg=0
    CALL cp_subsys_get(subsys, particles=particles,error=error)
    DO ip=1,particles%n_els
       DO idir=1,3
          idg=idg+1
          particles%els(ip)%r(idir)=x(idg)
       END DO
    END DO
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE cp_eval_at
