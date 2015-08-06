!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief computes preconditioners, and implements methods to apply them
!>      currently used in qs_ot
!> \par History
!>      - [UB] 2009-05-13 Adding stable approximate inverse (full and sparse)
!> \author Joost VandeVondele (09.2002)
! *****************************************************************************
MODULE preconditioner
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: ot_precond_full_all,&
                                             ot_precond_full_kinetic,&
                                             ot_precond_full_single,&
                                             ot_precond_full_single_inverse,&
                                             ot_precond_none,&
                                             ot_precond_s_inverse,&
                                             ot_precond_solver_update
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE preconditioner_apply,            ONLY: apply_preconditioner_dbcsr,&
                                             apply_preconditioner_fm
  USE preconditioner_makes,            ONLY: make_preconditioner_matrix
  USE preconditioner_solvers,          ONLY: solve_preconditioner,&
                                             transfer_dbcsr_to_fm,&
                                             transfer_fm_to_dbcsr
  USE preconditioner_types,            ONLY: destroy_preconditioner,&
                                             init_preconditioner,&
                                             preconditioner_p_type,&
                                             preconditioner_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'preconditioner'

  PUBLIC :: make_preconditioner, restart_preconditioner
  PUBLIC :: apply_preconditioner, prepare_preconditioner

! The public interface for apply preconditioner, the routines can be found in preconditioner_apply.F
  INTERFACE apply_preconditioner
     MODULE PROCEDURE apply_preconditioner_dbcsr
     MODULE PROCEDURE apply_preconditioner_fm
  END INTERFACE

! *****************************************************************************

CONTAINS

! *****************************************************************************

! creates a preconditioner for the system (H-energy_homo S)
! this preconditioner is (must be) symmetric positive definite.
! currently uses a atom-block-diagonal form
! each block will be  ....
! might overwrite matrix_h, matrix_t

! *****************************************************************************
!> \brief ...
!> \param preconditioner_env ...
!> \param precon_type ...
!> \param solver_type ...
!> \param matrix_h ...
!> \param matrix_s ...
!> \param matrix_t ...
!> \param mo_set ...
!> \param energy_gap ...
!> \param mixed_precision ...
!> \param convert_precond_to_dbcsr ...
!> \param error ...
!> \par History
!>      09.2014 removed some unused or unfinished methods
!>              removed sparse preconditioners and the 
!>              sparse approximate inverse at rev 14341 [Florian Schiffmann]
! *****************************************************************************
SUBROUTINE make_preconditioner(preconditioner_env, precon_type, solver_type, matrix_h, matrix_s, &
     matrix_t, mo_set, energy_gap, mixed_precision, convert_precond_to_dbcsr, error)

    TYPE(preconditioner_type)                :: preconditioner_env
    INTEGER, INTENT(IN)                      :: precon_type, solver_type
    TYPE(cp_dbcsr_type), POINTER             :: matrix_h
    TYPE(cp_dbcsr_type), OPTIONAL, POINTER   :: matrix_s, matrix_t
    TYPE(mo_set_type), POINTER               :: mo_set
    REAL(KIND=dp)                            :: energy_gap
    LOGICAL, INTENT(IN), OPTIONAL            :: mixed_precision, &
                                                convert_precond_to_dbcsr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'make_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, istat, k, &
                                                my_solver_type
    LOGICAL :: failure, my_convert_precond_to_dbcsr, my_mixed_precision, &
      needs_full_spectrum, needs_homo, use_mo_coeff_b
    REAL(KIND=dp)                            :: energy_homo
    REAL(KIND=dp), DIMENSION(:), POINTER     :: eigenvalues_ot
    TYPE(cp_dbcsr_type), POINTER             :: mo_coeff_b
    TYPE(cp_fm_type), POINTER                :: mo_coeff

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, mo_coeff_b=mo_coeff_b)
    use_mo_coeff_b = mo_set%use_mo_coeff_b
    CALL cp_fm_get_info(mo_coeff,ncol_global=k,error=error)

    ! Starting some matrix mess, check where to store the result in preconditioner_env, fm or dbcsr_matrix
    my_convert_precond_to_dbcsr = .FALSE.
    IF(PRESENT(convert_precond_to_dbcsr)) my_convert_precond_to_dbcsr = convert_precond_to_dbcsr

    ! Check whether we perform everything in single precision
    my_mixed_precision = .FALSE.
    IF(PRESENT(mixed_precision)) my_mixed_precision = mixed_precision
    IF(my_mixed_precision) THEN
       SELECT CASE(precon_type)
       CASE(ot_precond_full_kinetic,ot_precond_full_single_inverse)
          !supported
       CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "This precond with mixed precision is not supported yet")
       END SELECT
    ENDIF

    ! Thanks to the mess with the matrices we need to make sure in this case that the 
    ! Previous inverse is properly stored as a sparse matrix, fm gets deallocated here 
    ! if it wasn't anyway   
    IF(preconditioner_env%solver==ot_precond_solver_update)&
       CALL transfer_fm_to_dbcsr(preconditioner_env%fm,preconditioner_env%dbcsr_matrix,&
                                 matrix_h,my_mixed_precision,error)

    needs_full_spectrum=.FALSE.
    needs_homo=.FALSE.

    SELECT CASE(precon_type)
    CASE (ot_precond_full_all)
       needs_full_spectrum=.TRUE.
       ! both of them need the coefficients as fm's, more matrix mess
       IF(use_mo_coeff_b) THEN
          CALL copy_dbcsr_to_fm(mo_coeff_b,mo_coeff,error=error)
       ENDIF
    CASE (ot_precond_full_single)
       needs_homo=.TRUE.
       ! XXXX to be removed if homo estimate only is implemented
       needs_full_spectrum=.TRUE.
    CASE (ot_precond_full_kinetic,ot_precond_s_inverse,ot_precond_full_single_inverse)
       ! these should be happy without an estimate for the homo energy
       ! preconditioning can  not depend on an absolute eigenvalue, only on eigenvalue differences
    CASE DEFAULT
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "The preconditioner is unknown ...")
    END SELECT

    energy_homo=0.0_dp
    IF (needs_full_spectrum) THEN
       ALLOCATE(eigenvalues_ot(k),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       ! XXXXXXXXXXXXXXXX do not touch the initial MOs, could be harmful for either
       !                  the case of non-equivalent MOs but also for the derivate
       ! we could already have all eigenvalues e.g. full_all and we could skip this
       ! to be optimised later.
       ! one flaw is that not all SCF methods (i.e. that go over mo_derivs directly)
       ! have a 'valid' matrix_h... (we even don't know what evals are in that case)
       IF(use_mo_coeff_b) THEN
          CALL calculate_subspace_eigenvalues(mo_coeff_b,matrix_h,&
               eigenvalues_ot, do_rotation = .FALSE.,&
               para_env=mo_coeff%matrix_struct%para_env,&
               blacs_env=mo_coeff%matrix_struct%context,error=error)
       ELSE
          CALL calculate_subspace_eigenvalues(mo_coeff,matrix_h,&
               eigenvalues_ot, do_rotation = .FALSE.,error=error)
       ENDIF
       IF (k>0) energy_homo=eigenvalues_ot(k)
    ELSE
       IF (needs_homo) THEN
          CALL stop_program(routineN,moduleN,__LINE__,&
                            "Not yet implemented")
       ENDIF
    ENDIF

    ! After all bits and pieces of checking and initialization, here comes the 
    ! part where the preconditioner matrix gets created and solved. 
    ! This will give the matrices for later use
    my_solver_type = solver_type
    preconditioner_env%in_use=precon_type
    CALL make_preconditioner_matrix(preconditioner_env, matrix_h, matrix_s, matrix_t, mo_coeff,&
                       energy_homo, eigenvalues_ot, energy_gap, &
                       my_solver_type, error)

    CALL solve_preconditioner(my_solver_type,preconditioner_env,matrix_s,matrix_h,my_mixed_precision,error)
  
    ! Here comes more matrix mess, make sure to output the correct matrix format, 
    ! A bit pointless to convert the cholesky factorized version as it doesn't work in
    ! dbcsr form and will crash later,...
    IF(my_convert_precond_to_dbcsr) THEN
       CALL transfer_fm_to_dbcsr(preconditioner_env%fm,preconditioner_env%dbcsr_matrix,&
                                 matrix_h,my_mixed_precision,error)
    ELSE
       CALL transfer_dbcsr_to_fm(preconditioner_env%dbcsr_matrix,preconditioner_env%fm, &
                                 preconditioner_env%para_env, preconditioner_env%ctxt,&
                                 my_mixed_precision,error)
    ENDIF


    IF (needs_full_spectrum) THEN
      DEALLOCATE(eigenvalues_ot,STAT=istat)
      CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE make_preconditioner

! *****************************************************************************
!> \brief Allows for a restart of the preconditioner
!>        depending on the method it purges all arrays or keeps them
!> \param qs_env ...
!> \param preconditioner ...
!> \param prec_type ...
!> \param nspins ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE restart_preconditioner(qs_env,preconditioner,prec_type,nspins,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(preconditioner_p_type), &
      DIMENSION(:), POINTER                  :: preconditioner
    INTEGER, INTENT(IN)                      :: prec_type, nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'restart_preconditioner', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, stat
    LOGICAL                                  :: failure
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_para_env_type), POINTER          :: para_env

    NULLIFY(para_env, blacs_env)
    CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env, error=error)

    IF (ASSOCIATED(preconditioner)) THEN
       SELECT CASE(prec_type)
       CASE(ot_precond_full_all,ot_precond_full_single) ! these depend on the ks matrix
         DO ispin=1,SIZE(preconditioner)
            CALL destroy_preconditioner(preconditioner(ispin)%preconditioner,error=error)
            DEALLOCATE(preconditioner(ispin)%preconditioner)
         ENDDO
         DEALLOCATE(preconditioner,stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       CASE(ot_precond_none,ot_precond_full_kinetic,ot_precond_s_inverse,ot_precond_full_single_inverse) ! these are 'independent'
         ! do nothing
       CASE DEFAULT
           CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
    END IF

    ! add an OT preconditioner if none is present
    IF (.NOT.ASSOCIATED(preconditioner)) THEN
         SELECT CASE(prec_type)
         CASE(ot_precond_full_all,ot_precond_full_single_inverse)
            ALLOCATE(preconditioner(nspins), stat=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         CASE DEFAULT
            ALLOCATE(preconditioner(1), stat=stat)
            CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
         END SELECT
         DO ispin=1,SIZE(preconditioner)
            ALLOCATE(preconditioner(ispin)%preconditioner)
            CALL init_preconditioner(preconditioner(ispin)%preconditioner,&
                                     para_env=para_env,&
                                     blacs_env=blacs_env,error=error)
         ENDDO
    END IF

  END SUBROUTINE restart_preconditioner

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mos ...
!> \param matrix_ks ...
!> \param matrix_s ...
!> \param ot_preconditioner ...
!> \param prec_type ...
!> \param solver_type ...
!> \param energy_gap ...
!> \param nspins ...
!> \param has_unit_metric ...
!> \param mixed_precision ...
!> \param convert_to_dbcsr ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE prepare_preconditioner(qs_env,mos,matrix_ks,matrix_s,&
                                    ot_preconditioner,prec_type,solver_type,&
                                    energy_gap,nspins,has_unit_metric,mixed_precision,&
                                    convert_to_dbcsr,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(preconditioner_p_type), &
      DIMENSION(:), POINTER                  :: ot_preconditioner
    INTEGER, INTENT(IN)                      :: prec_type, solver_type
    REAL(dp), INTENT(IN)                     :: energy_gap
    INTEGER, INTENT(IN)                      :: nspins
    LOGICAL, INTENT(IN), OPTIONAL            :: has_unit_metric, &
                                                mixed_precision, &
                                                convert_to_dbcsr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'prepare_preconditioner', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_string_length)     :: msg
    INTEGER                                  :: handle, icall, ispin, n_loops
    LOGICAL :: do_co_rotate, failure, my_convert_to_dbcsr, &
      my_has_unit_metric, my_mixed_precision, use_mo_coeff_b
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: kinetic
    TYPE(cp_dbcsr_type), POINTER             :: matrix_t, mo_coeff_b
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control

    CALL timeset(routineN,handle)
    NULLIFY(matrix_t, mo_coeff_b, mo_coeff, kinetic, dft_control, para_env, blacs_env)
    failure = .FALSE.
    my_has_unit_metric = .FALSE.
    IF(PRESENT(has_unit_metric)) my_has_unit_metric = has_unit_metric
    my_mixed_precision = .FALSE.
    IF(PRESENT(mixed_precision)) my_mixed_precision = mixed_precision
    my_convert_to_dbcsr = .TRUE.
    IF(PRESENT(convert_to_dbcsr))my_convert_to_dbcsr=convert_to_dbcsr

    CALL get_qs_env(qs_env,&
                    dft_control=dft_control,&
                    para_env=para_env,&
                    blacs_env=blacs_env,&
                    error=error)

    IF(dft_control%qs_control%semi_empirical .OR. dft_control%qs_control%dftb) THEN
      IF(prec_type==ot_precond_full_kinetic) THEN
          msg="Full_kinetic not available for semi-empirical methods"
          CPErrorMessage(cp_failure_level,routineP,TRIM(msg),error)
      END IF
      matrix_t => matrix_s(1)%matrix
    ELSE
      CPPrecondition(.NOT. my_has_unit_metric,cp_failure_level,routineP,error,failure)
      CALL get_qs_env(qs_env, kinetic=kinetic, error=error)
      matrix_t => kinetic(1)%matrix
    END IF

    !determins how often make preconditioner is called, spin dependent methods have to be called twice
    n_loops=1
    IF(prec_type==ot_precond_full_single_inverse)n_loops=nspins
    ! check whether we need the ev and rotate the MOs
    SELECT CASE(prec_type)
    CASE(ot_precond_full_all)
       ! if one of these preconditioners is used every spin needs to call make_preconditioner
       n_loops=nspins

       do_co_rotate = ASSOCIATED(qs_env%mo_derivs)
       DO ispin=1,nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff_b=mo_coeff_b,mo_coeff=mo_coeff)
         use_mo_coeff_b =mos(ispin)%mo_set%use_mo_coeff_b
         IF(use_mo_coeff_b .AND. do_co_rotate) THEN
           CALL calculate_subspace_eigenvalues(mo_coeff_b,matrix_ks(ispin)%matrix,&
                do_rotation = .TRUE., &
                co_rotate=qs_env%mo_derivs(ispin)%matrix,&
                para_env=para_env,&
                blacs_env=blacs_env,error=error)
         ELSEIF(use_mo_coeff_b) THEN
           CALL calculate_subspace_eigenvalues(mo_coeff_b,matrix_ks(ispin)%matrix,&
                do_rotation = .TRUE., &
                para_env=para_env,&
                blacs_env=blacs_env,error=error)
         ELSE
           CALL calculate_subspace_eigenvalues(mo_coeff,matrix_ks(ispin)%matrix,&
                do_rotation = .TRUE., error=error)
         END IF
       END DO
    CASE DEFAULT
      ! No need to rotate the MOs
    END SELECT

    ! check whether we have a preconditioner
    SELECT CASE(prec_type)
       CASE(ot_precond_none)
          DO ispin = 1,SIZE(ot_preconditioner)
            ot_preconditioner(ispin)%preconditioner%in_use=0
          END DO
       CASE DEFAULT
          DO icall=1,n_loops
             IF(my_has_unit_metric) THEN
               CALL make_preconditioner(ot_preconditioner(icall)%preconditioner, &
                        prec_type, &
                        solver_type, &
                        matrix_h=matrix_ks(icall)%matrix,&
                        mo_set=mos(icall)%mo_set,&
                        energy_gap=energy_gap,&
                        convert_precond_to_dbcsr=my_convert_to_dbcsr,&
                        error=error)
             ELSE
               CALL make_preconditioner(ot_preconditioner(icall)%preconditioner, &
                        prec_type, &
                        solver_type, &
                        matrix_h=matrix_ks(icall)%matrix,&
                        matrix_s=matrix_s(1)%matrix,&
                        matrix_t=matrix_t, &
                        mo_set=mos(icall)%mo_set,&
                        energy_gap=energy_gap,&
                        mixed_precision=my_mixed_precision,&
                        convert_precond_to_dbcsr=my_convert_to_dbcsr,&
                        error=error)
             END IF
          END DO
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE prepare_preconditioner

END MODULE preconditioner

