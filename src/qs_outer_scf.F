!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2006  CP2K developers group
!-----------------------------------------------------------------------------!

!!****** cp2k/qs_outer_scf [1.0] *
!!
!!   NAME
!!     qs_outer_scf
!!
!!   FUNCTION
!!     Routines for performing an outer scf loop
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   MODIFICATION HISTORY
!!     Created [2006.03]
!!
!!   SOURCE
!******************************************************************************

MODULE qs_outer_scf
  USE cp_control_types,                ONLY: ddapc_restraint_type,&
                                             s2_restraint_type,&
                                             dft_control_type
  USE input_constants,                 ONLY: do_ddapc_constraint,&
                                             do_s2_constraint,&
                                             outer_scf_ddapc_constraint,&
                                             outer_scf_s2_constraint,&
                                             outer_scf_optimizer_diis,&
                                             outer_scf_optimizer_sd
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: diamat_all
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_scf_types,                    ONLY: qs_scf_env_type
  USE scf_control_types,               ONLY: scf_control_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_outer_scf'

! *** Public subroutines ***

  PUBLIC :: outer_loop_gradient,outer_loop_optimize,outer_loop_update_qs_env, &
            outer_loop_variables_count,outer_loop_extrapolate

CONTAINS
!!****f* qs_outer_scf/outer_loop_variables_count ***
!!
!!   NAME
!!      outer_loop_variables_count
!!
!!   FUNCTION
!!      returns the number of variables that is employed in the outer loop
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
FUNCTION outer_loop_variables_count(scf_control,error) RESULT(res)
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error
    INTEGER                                  :: res

    SELECT CASE(scf_control%outer_scf%type)
    CASE (outer_scf_ddapc_constraint)
         res=1
    CASE (outer_scf_s2_constraint)
         res=1
    CASE DEFAULT
         res=0
    END SELECT

END FUNCTION outer_loop_variables_count


!!****f* qs_outer_scf/outer_loop_gradient ***
!!
!!   NAME
!!      outer_loop_gradient
!!
!!   FUNCTION
!!      computes the gradient wrt to the outer loop variables
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE outer_loop_gradient(qs_env, scf_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'outer_loop_gradient', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ihistory
    LOGICAL                                  :: failure, is_constraint
    TYPE(ddapc_restraint_type), POINTER      :: ddapc_restraint_control
    TYPE(s2_restraint_type), POINTER         :: s2_restraint_control
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(scf_control_type), POINTER          :: scf_control

    CALL timeset(routineN,handle)
    failure = .FALSE.

    CALL get_qs_env(qs_env=qs_env,scf_control=scf_control, dft_control=dft_control,energy=energy)
    CPPrecondition(scf_control%outer_scf%have_scf,cp_failure_level,routineP,error,failure)

    ihistory=scf_env%outer_scf%iter_count
    CPPrecondition(ihistory<=SIZE(scf_env%outer_scf%energy,1),cp_failure_level,routineP,error,failure)
    
    scf_env%outer_scf%energy(ihistory)=energy%total

    SELECT CASE(scf_control%outer_scf%type)
    CASE (outer_scf_ddapc_constraint)
       CPPrecondition(dft_control%qs_control%ddapc_restraint,cp_failure_level,routineP,error,failure)
       ddapc_restraint_control=>dft_control%qs_control%ddapc_restraint_control
       is_constraint=(ddapc_restraint_control%functional_form==do_ddapc_constraint)
       CPPrecondition(is_constraint,cp_failure_level,routineP,error,failure)

       scf_env%outer_scf%variables(:,ihistory)=ddapc_restraint_control%strength
       scf_env%outer_scf%gradient(:,ihistory) =ddapc_restraint_control%ddapc_order_p- &
                                      ddapc_restraint_control%target
    CASE (outer_scf_s2_constraint)
       CPPrecondition(dft_control%qs_control%s2_restraint,cp_failure_level,routineP,error,failure)
       s2_restraint_control=>dft_control%qs_control%s2_restraint_control
       is_constraint=(s2_restraint_control%functional_form==do_s2_constraint)
       CPPrecondition(is_constraint,cp_failure_level,routineP,error,failure)

       scf_env%outer_scf%variables(:,ihistory)=s2_restraint_control%strength
       scf_env%outer_scf%gradient(:,ihistory) =s2_restraint_control%s2_order_p- &
                                      s2_restraint_control%target
    CASE DEFAULT
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL timestop(handle)

END SUBROUTINE

!!****f* qs_outer_scf/outer_loop_optimize ***
!!
!!   NAME
!!      outer_loop_optimize
!!
!!   FUNCTION
!!      optimizes the parameters of the outer_scf
!!
!!   NOTES
!!      ought to be general, and independent of the actual kind of variables
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE outer_loop_optimize(scf_env, scf_control, error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'outer_loop_optimize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ibuf, ihistory, j, &
                                                jbuf, nb
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipivot
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ev
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: a, b

    CALL timeset(routineN,handle)

    failure=.FALSE.
    ihistory=scf_env%outer_scf%iter_count

    SELECT CASE(scf_control%outer_scf%optimizer)
    CASE(outer_scf_optimizer_sd) 
      ! Notice that we are just trying to find a stationary point
      ! e.g. the ddpac_constraint, one maximizes the function, so the stepsize might have
      ! to be negative 
      scf_env%outer_scf%variables(:,ihistory+1)=scf_env%outer_scf%variables(:,ihistory) - &
                scf_control%outer_scf%step_size*scf_env%outer_scf%gradient(:,ihistory)
    CASE(outer_scf_optimizer_diis)
      CPPrecondition(scf_control%outer_scf%diis_buffer_length>0,cp_failure_level,routineP,error,failure)
      ! set up DIIS matrix
      nb=MIN(ihistory,scf_control%outer_scf%diis_buffer_length)
      IF (nb<2) THEN
         scf_env%outer_scf%variables(:,ihistory+1)=scf_env%outer_scf%variables(:,ihistory) - &
                        scf_control%outer_scf%step_size*scf_env%outer_scf%gradient(:,ihistory)

      ELSE
         ALLOCATE(b(nb+1,nb+1),a(nb+1,nb+1),ev(nb+1),ipivot(nb+1))
         DO I=1,nb
          DO J=I,nb
           ibuf=ihistory-nb+i
           jbuf=ihistory-nb+j
           b(I,J)= DOT_PRODUCT(scf_env%outer_scf%gradient(:,ibuf), &
                               scf_env%outer_scf%gradient(:,jbuf))
           b(J,I)=b(I,J)
          ENDDO
         ENDDO
         b(nb+1,:)=-1.0_dp
         b(:,nb+1)=-1.0_dp
         b(nb+1,nb+1)=0.0_dp

         CALL diamat_all(b,ev)
         a=b
         DO I=1,nb+1
            IF (ABS(ev(I)).LT.1.0E-12_dp) THEN
               a(:,I)=0.0_dp
            ELSE 
               a(:,I)=a(:,I)/ev(I)             
            ENDIF
         END DO
         ev=-MATMUL(a,b(nb+1,:))

         scf_env%outer_scf%variables(:,ihistory+1)=0.0_dp
         DO i=1,nb
            ibuf=ihistory-nb+i
            scf_env%outer_scf%variables(:,ihistory+1)=scf_env%outer_scf%variables(:,ihistory+1)+ &
                                                      ev(i)*scf_env%outer_scf%variables(:,ibuf)
         ENDDO
         DEALLOCATE(b,ev)
       ENDIF
    CASE DEFAULT
      CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT 

    CALL timestop(handle)

END SUBROUTINE outer_loop_optimize

!!****f* qs_outer_scf/outer_loop_update_qs_env ***
!!
!!   NAME
!!      outer_loop_set
!!
!!   FUNCTION
!!      propagates the updated variables to wherever they need to be set in 
!!      qs_env
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE outer_loop_update_qs_env(qs_env, scf_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'outer_loop_update_qs_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ihistory
    LOGICAL                                  :: failure
    TYPE(ddapc_restraint_type), POINTER      :: ddapc_restraint_control
    TYPE(s2_restraint_type), POINTER         :: s2_restraint_control
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(scf_control_type), POINTER          :: scf_control

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env,scf_control=scf_control, dft_control=dft_control)
    ihistory=scf_env%outer_scf%iter_count

    SELECT CASE(scf_control%outer_scf%type)
    CASE (outer_scf_ddapc_constraint)
       ddapc_restraint_control=>dft_control%qs_control%ddapc_restraint_control
       ddapc_restraint_control%strength=scf_env%outer_scf%variables(1,ihistory+1)
    CASE (outer_scf_s2_constraint)
       s2_restraint_control=>dft_control%qs_control%s2_restraint_control
       s2_restraint_control%strength=scf_env%outer_scf%variables(1,ihistory+1)
    CASE DEFAULT
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL timestop(handle)

END SUBROUTINE outer_loop_update_qs_env
!!****f* qs_outer_scf/outer_loop_extrapolate ***
!!
!!   NAME
!!      outer_loop_extrapolate
!!
!!   FUNCTION
!!      uses the outer_scf_history to extrapolate new values for the variables
!!      and updates their value in qs_env accordingly
!!
!!   NOTES
!!      it assumes that the current value of qs_env still needs to be added to the history
!!      simple multilinear extrapolation is employed
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     03.2006 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE outer_loop_extrapolate(qs_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'outer_loop_extrapolate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ihis, ivec, nhistory, &
                                                nvariables, nvec, &
                                                outer_scf_ihistory, stat
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: alpha
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: extrapolation
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: outer_scf_history
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(scf_control_type), POINTER          :: scf_control

    CALL timeset(routineN,handle)

    failure = .FALSE.
 
    CALL get_qs_env(qs_env,outer_scf_history=outer_scf_history, &
                           outer_scf_ihistory=outer_scf_ihistory, &
                           scf_control=scf_control, dft_control=dft_control)

    nvariables=SIZE(outer_scf_history,1)
    nhistory=SIZE(outer_scf_history,2)
    ALLOCATE(extrapolation(nvariables),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    CPPrecondition(nhistory>0,cp_failure_level,routineP,error,failure)

    ! add the current version of qs_env to the history
    outer_scf_ihistory=outer_scf_ihistory+1
    ivec=1+MODULO(outer_scf_ihistory-1,nhistory)
    SELECT CASE(scf_control%outer_scf%type)
    CASE (outer_scf_ddapc_constraint)
         outer_scf_history(1,ivec)= &
           dft_control%qs_control%ddapc_restraint_control%strength
    CASE (outer_scf_s2_constraint)
         outer_scf_history(1,ivec)= &
           dft_control%qs_control%s2_restraint_control%strength
    CASE DEFAULT
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT
    CALL set_qs_env(qs_env,outer_scf_ihistory=outer_scf_ihistory)     

    ! multilinear extrapolation
    nvec=MIN(nhistory,outer_scf_ihistory)
    alpha=nvec
    ivec=1+MODULO(outer_scf_ihistory-1,nhistory)
    extrapolation=alpha*outer_scf_history(:,ivec)
    DO ihis=2,nvec
       alpha = -1.0_dp * alpha * REAL(nvec - ihis + 1 , dp ) / REAL(ihis , dp )
       ivec=1+MODULO(outer_scf_ihistory-ihis,nhistory)
       extrapolation=extrapolation+alpha*outer_scf_history(:,ivec)
    ENDDO

    ! update qs_env to use this extrapolation        
    SELECT CASE(scf_control%outer_scf%type)
    CASE (outer_scf_ddapc_constraint)
           dft_control%qs_control%ddapc_restraint_control%strength=extrapolation(1)
    CASE (outer_scf_s2_constraint)
           dft_control%qs_control%s2_restraint_control%strength=extrapolation(1)
    CASE DEFAULT
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    DEALLOCATE(extrapolation,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

END SUBROUTINE outer_loop_extrapolate

END MODULE qs_outer_scf
