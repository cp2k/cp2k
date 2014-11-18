!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for performing an outer scf loop
!> \par History
!>      Created [2006.03]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE qs_outer_scf
  USE cp_control_types,                ONLY: becke_restraint_type,&
                                             ddapc_restraint_type,&
                                             dft_control_type,&
                                             s2_restraint_type
  USE input_constants,                 ONLY: &
       do_ddapc_constraint, do_s2_constraint, outer_scf_becke_constraint, &
       outer_scf_ddapc_constraint, outer_scf_none, &
       outer_scf_optimizer_bisect, outer_scf_optimizer_diis, &
       outer_scf_optimizer_none, outer_scf_optimizer_sd, &
       outer_scf_s2_constraint, outer_scf_scp
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
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_outer_scf'

! *** Public subroutines ***

  PUBLIC :: outer_loop_gradient,outer_loop_optimize,outer_loop_update_qs_env, &
            outer_loop_variables_count,outer_loop_extrapolate

CONTAINS

! *****************************************************************************
!> \brief returns the number of variables that is employed in the outer loop
!> \param scf_control ...
!> \param error ...
!> \retval res ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  FUNCTION outer_loop_variables_count(scf_control,error) RESULT(res)
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(cp_error_type), INTENT(inout)       :: error
    INTEGER                                  :: res

    SELECT CASE(scf_control%outer_scf%type)
    CASE (outer_scf_ddapc_constraint)
       res=1
    CASE (outer_scf_s2_constraint)
       res=1
    CASE (outer_scf_becke_constraint)
       res=1
    CASE (outer_scf_none) ! just needed to communicate the gradient criterium
       res=1
    CASE (outer_scf_scp)
       res=1
    CASE DEFAULT
       res=0
    END SELECT

  END FUNCTION outer_loop_variables_count

! *****************************************************************************
!> \brief computes the gradient wrt to the outer loop variables
!> \param qs_env ...
!> \param scf_env ...
!> \param error ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE outer_loop_gradient(qs_env, scf_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'outer_loop_gradient', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ihistory, n
    LOGICAL                                  :: failure, is_constraint
    TYPE(becke_restraint_type), POINTER      :: becke_control
    TYPE(ddapc_restraint_type), POINTER      :: ddapc_restraint_control
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(s2_restraint_type), POINTER         :: s2_restraint_control
    TYPE(scf_control_type), POINTER          :: scf_control

    CALL timeset(routineN,handle)
    failure = .FALSE.

    CALL get_qs_env(qs_env=qs_env,scf_control=scf_control, &
         dft_control=dft_control,energy=energy,error=error)
    CPPrecondition(scf_control%outer_scf%have_scf,cp_failure_level,routineP,error,failure)

    ihistory=scf_env%outer_scf%iter_count
    CPPrecondition(ihistory<=SIZE(scf_env%outer_scf%energy,1),cp_failure_level,routineP,error,failure)

    scf_env%outer_scf%energy(ihistory)=energy%total

    SELECT CASE(scf_control%outer_scf%type)
    CASE (outer_scf_none)
       ! just pass the inner loop scf criterium to the outer loop one
       scf_env%outer_scf%variables(1,ihistory)=scf_env%iter_delta
       scf_env%outer_scf%gradient(1,ihistory)=scf_env%iter_delta
    CASE (outer_scf_ddapc_constraint)
       CPPrecondition(dft_control%qs_control%ddapc_restraint,cp_failure_level,routineP,error,failure)
       DO n=1,SIZE(dft_control%qs_control%ddapc_restraint_control)
          NULLIFY(ddapc_restraint_control)
          ddapc_restraint_control=>dft_control%qs_control%ddapc_restraint_control(n)%ddapc_restraint_control
          is_constraint=(ddapc_restraint_control%functional_form==do_ddapc_constraint)
          IF(is_constraint)EXIT
       END DO
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
    CASE (outer_scf_becke_constraint)
       CPPrecondition(dft_control%qs_control%becke_restraint,cp_failure_level,routineP,error,failure)
       becke_control=>dft_control%qs_control%becke_control
       scf_env%outer_scf%variables(:,ihistory)=becke_control%strength
       scf_env%outer_scf%gradient(:,ihistory) =becke_control%becke_order_p- &
            becke_control%target
    CASE (outer_scf_scp)
       scf_env%outer_scf%variables(1,ihistory)=scf_env%outer_scf%energy(ihistory)

       IF ( ihistory == 1 ) THEN
          scf_env%outer_scf%gradient(1,ihistory)=scf_env%outer_scf%energy(ihistory)
       ELSE
          scf_env%outer_scf%gradient(1,ihistory)=scf_env%outer_scf%energy(ihistory)- &
               scf_env%outer_scf%energy(ihistory-1)
       ENDIF
    CASE DEFAULT
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE outer_loop_gradient

! *****************************************************************************
!> \brief optimizes the parameters of the outer_scf
!> \param scf_env ...
!> \param scf_control ...
!> \param error ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
!> \note
!>       ought to be general, and independent of the actual kind of variables
! *****************************************************************************
  SUBROUTINE outer_loop_optimize(scf_env, scf_control, error)
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(scf_control_type), POINTER          :: scf_control
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'outer_loop_optimize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ibuf, ihigh, &
                                                ihistory, ilow, j, jbuf, nb, &
                                                optimizer_type
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: ipivot
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: interval, tmp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: ev
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: a, b

    CALL timeset(routineN,handle)

    failure=.FALSE.
    ihistory=scf_env%outer_scf%iter_count
    optimizer_type=scf_control%outer_scf%optimizer

100 CONTINUE ! if we need a different run type we'll restart here

    SELECT CASE(optimizer_type)
    CASE(outer_scf_optimizer_bisect) ! bisection on the gradient, needs to be 1D
       CPPrecondition(SIZE(scf_env%outer_scf%gradient(:,1))==1,cp_failure_level,routineP,error,failure)
       ! find the pair of points that bracket a zero of the gradient, with the smallest interval possible
       ilow=-1
       ihigh=-1
       interval=HUGE(interval)
       DO i=1,ihistory
          DO j=i+1,ihistory
             ! distrust often used points
             IF (scf_env%outer_scf%count(i).GT.scf_control%outer_scf%bisect_trust_count) CYCLE
             IF (scf_env%outer_scf%count(j).GT.scf_control%outer_scf%bisect_trust_count) CYCLE

             ! if they bracket a zero use them
             IF (scf_env%outer_scf%gradient(1,i) * &
                  scf_env%outer_scf%gradient(1,j) < 0.0_dp) THEN
                tmp = ABS(scf_env%outer_scf%variables(1,i)-scf_env%outer_scf%variables(1,j))
                IF (tmp < interval) THEN
                   ilow = i
                   ihigh = j
                   interval = tmp
                ENDIF
             ENDIF
          ENDDO
       ENDDO
       IF (ilow == -1) THEN ! we didn't bracket a minimum yet, try something else
          optimizer_type=outer_scf_optimizer_diis
          GOTO 100
       ENDIF
       scf_env%outer_scf%count(ilow)=scf_env%outer_scf%count(ilow)+1
       scf_env%outer_scf%count(ihigh)=scf_env%outer_scf%count(ihigh)+1
       scf_env%outer_scf%variables(:,ihistory+1)=0.5_dp*(scf_env%outer_scf%variables(:,ilow)+&
            scf_env%outer_scf%variables(:,ihigh))
    CASE(outer_scf_optimizer_none)
       scf_env%outer_scf%variables(:,ihistory+1)=scf_env%outer_scf%variables(:,ihistory)
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
          optimizer_type = outer_scf_optimizer_sd
          GOTO 100
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

          CALL diamat_all(b,ev,error=error)
          a(:,:)=b
          DO I=1,nb+1
             IF (ABS(ev(I)).LT.1.0E-12_dp) THEN
                a(:,I)=0.0_dp
             ELSE
                a(:,I)=a(:,I)/ev(I)
             ENDIF
          END DO
          ev(:)=-MATMUL(a,b(nb+1,:))

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

! *****************************************************************************
!> \brief propagates the updated variables to wherever they need to be set in
!>       qs_env
!> \param qs_env ...
!> \param scf_env ...
!> \param error ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE outer_loop_update_qs_env(qs_env, scf_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_scf_env_type), POINTER           :: scf_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'outer_loop_update_qs_env', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ihistory, n
    LOGICAL                                  :: failure, is_constraint
    TYPE(becke_restraint_type), POINTER      :: becke_control
    TYPE(ddapc_restraint_type), POINTER      :: ddapc_restraint_control
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(s2_restraint_type), POINTER         :: s2_restraint_control
    TYPE(scf_control_type), POINTER          :: scf_control

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env=qs_env,scf_control=scf_control, dft_control=dft_control,error=error)
    ihistory=scf_env%outer_scf%iter_count

    SELECT CASE(scf_control%outer_scf%type)
    CASE (outer_scf_none)
       ! do nothing
    CASE (outer_scf_scp)
       ! do nothing
    CASE (outer_scf_ddapc_constraint)
       DO n=1,SIZE(dft_control%qs_control%ddapc_restraint_control)
          NULLIFY(ddapc_restraint_control)
          ddapc_restraint_control=>dft_control%qs_control%ddapc_restraint_control(n)%ddapc_restraint_control
          is_constraint=(ddapc_restraint_control%functional_form==do_ddapc_constraint)
          IF(is_constraint)EXIT
       END DO
       ddapc_restraint_control%strength=scf_env%outer_scf%variables(1,ihistory+1)
    CASE (outer_scf_s2_constraint)
       s2_restraint_control=>dft_control%qs_control%s2_restraint_control
       s2_restraint_control%strength=scf_env%outer_scf%variables(1,ihistory+1)
    CASE (outer_scf_becke_constraint)
       becke_control=>dft_control%qs_control%becke_control
       becke_control%strength=scf_env%outer_scf%variables(1,ihistory+1)
    CASE DEFAULT
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE outer_loop_update_qs_env

! *****************************************************************************
!> \brief uses the outer_scf_history to extrapolate new values for the variables
!>       and updates their value in qs_env accordingly
!> \param qs_env ...
!> \param error ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
!> \note
!>       it assumes that the current value of qs_env still needs to be added to the history
!>       simple multilinear extrapolation is employed
! *****************************************************************************
  SUBROUTINE outer_loop_extrapolate(qs_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'outer_loop_extrapolate', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ihis, ivec, n, &
                                                nhistory, nvariables, nvec, &
                                                outer_scf_ihistory, stat
    LOGICAL                                  :: failure, is_constraint
    REAL(kind=dp)                            :: alpha
    REAL(kind=dp), ALLOCATABLE, DIMENSION(:) :: extrapolation
    REAL(kind=dp), DIMENSION(:, :), POINTER  :: outer_scf_history
    TYPE(ddapc_restraint_type), POINTER      :: ddapc_restraint_control
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(scf_control_type), POINTER          :: scf_control

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CALL get_qs_env(qs_env,outer_scf_history=outer_scf_history, &
         outer_scf_ihistory=outer_scf_ihistory, &
         scf_control=scf_control, dft_control=dft_control,error=error)

    nvariables=SIZE(outer_scf_history,1)
    nhistory=SIZE(outer_scf_history,2)
    ALLOCATE(extrapolation(nvariables),STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    CPPrecondition(nhistory>0,cp_failure_level,routineP,error,failure)

    ! add the current version of qs_env to the history
    outer_scf_ihistory=outer_scf_ihistory+1
    ivec=1+MODULO(outer_scf_ihistory-1,nhistory)
    SELECT CASE(scf_control%outer_scf%type)
    CASE (outer_scf_none)
       outer_scf_history(1,ivec)=0.0_dp
    CASE (outer_scf_scp)
       outer_scf_history(1,ivec)=0.0_dp
    CASE (outer_scf_ddapc_constraint)
       DO n=1,SIZE(dft_control%qs_control%ddapc_restraint_control)
          NULLIFY(ddapc_restraint_control)
          ddapc_restraint_control=>dft_control%qs_control%ddapc_restraint_control(n)%ddapc_restraint_control
          is_constraint=(ddapc_restraint_control%functional_form==do_ddapc_constraint)
          IF(is_constraint)EXIT
       END DO
       outer_scf_history(1,ivec)= &
            ddapc_restraint_control%strength
    CASE (outer_scf_s2_constraint)
       outer_scf_history(1,ivec)= &
            dft_control%qs_control%s2_restraint_control%strength
    CASE (outer_scf_becke_constraint)
       outer_scf_history(1,ivec)= &
            dft_control%qs_control%becke_control%strength
    CASE DEFAULT
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT
    CALL set_qs_env(qs_env,outer_scf_ihistory=outer_scf_ihistory,error=error)

    ! multilinear extrapolation
    nvec=MIN(nhistory,outer_scf_ihistory)
    alpha=nvec
    ivec=1+MODULO(outer_scf_ihistory-1,nhistory)
    extrapolation(:)=alpha*outer_scf_history(:,ivec)
    DO ihis=2,nvec
       alpha = -1.0_dp * alpha * REAL(nvec - ihis + 1 , dp ) / REAL(ihis , dp )
       ivec=1+MODULO(outer_scf_ihistory-ihis,nhistory)
       extrapolation(:)=extrapolation+alpha*outer_scf_history(:,ivec)
    ENDDO

    ! update qs_env to use this extrapolation
    SELECT CASE(scf_control%outer_scf%type)
    CASE (outer_scf_none)
       ! nothing
    CASE (outer_scf_scp)
       ! nothing
    CASE (outer_scf_ddapc_constraint)
       ddapc_restraint_control%strength=extrapolation(1)
    CASE (outer_scf_s2_constraint)
       dft_control%qs_control%s2_restraint_control%strength=extrapolation(1)
    CASE (outer_scf_becke_constraint)
       dft_control%qs_control%becke_control%strength=extrapolation(1)
    CASE DEFAULT
       CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    DEALLOCATE(extrapolation,STAT=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE outer_loop_extrapolate

END MODULE qs_outer_scf
