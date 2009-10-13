!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author CJM
! *****************************************************************************
MODULE scp_environment_types
  USE cp_blacs_env,                    ONLY: cp_blacs_env_release,&
                                             cp_blacs_env_retain
  USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_blacs_env_type,&
                                             cp_para_env_type
  USE dbcsr_types,                     ONLY: dbcsr_obj
  USE distribution_2d_types,           ONLY: distribution_2d_release,&
                                             distribution_2d_retain,&
                                             distribution_2d_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE qs_ks_scp_types,                 ONLY: qs_ks_scp_env_type,&
                                             qs_ks_scp_release,&
                                             qs_ks_scp_retain
  USE scp_coeff_types,                 ONLY: aux_coeff_set_release,&
                                             aux_coeff_set_retain,&
                                             aux_coeff_set_type
  USE scp_energy_types,                ONLY: deallocate_scp_energy,&
                                             scp_energy_type
  USE scp_force_types,                 ONLY: deallocate_scp_force,&
                                             scp_force_type
  USE scp_rho_types,                   ONLY: scp_rho_release,&
                                             scp_rho_retain,&
                                             scp_rho_type
  USE scp_rspw_types,                  ONLY: scp_rspw_release,&
                                             scp_rspw_retain,&
                                             scp_rspw_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_environment_types'
  INTEGER, PRIVATE, SAVE :: last_scp_env_id_nr=0

! *****************************************************************************
  TYPE scp_environment_type
    INTEGER :: id_nr, ref_count
    LOGICAL :: use_ref_cell
    REAL ( dp )                                           :: eps_gvg_rspace, &
                                                             eps_rho_rspace
! FOR NDDO
    TYPE ( dbcsr_obj ), POINTER                    :: pscp, ks_scp
    REAL ( dp ), POINTER                                  :: ascp_diag ( : )
! FOR DFT
    TYPE ( aux_coeff_set_type ), POINTER                  :: aux_coeff_set
    TYPE ( scp_energy_type ), POINTER                     :: energy
    TYPE ( scp_force_type ), POINTER                      :: force ( : )
    TYPE ( scp_rho_type ),  POINTER                       :: rho_scp
    TYPE ( scp_rspw_type ),  POINTER                      :: rspw
    TYPE ( qs_ks_scp_env_type ),  POINTER                 :: ks_scp_env
    TYPE ( section_vals_type ), POINTER                   :: input
    TYPE ( cp_para_env_type ),  POINTER                   :: para_env
    TYPE ( cp_blacs_env_type ), POINTER                   :: blacs_env
    TYPE ( distribution_2d_type ), POINTER                :: distribution_2d
  END TYPE scp_environment_type

! *** Public data types ***

  PUBLIC :: scp_environment_type

! *** Public subroutines ***

  PUBLIC :: get_scp_env,&
            init_scp_env,&
            set_scp_env, &
            scp_env_create, &
            scp_env_retain, &
            scp_env_release

CONTAINS

! *****************************************************************************
  SUBROUTINE get_scp_env( scp_env, aux_coeff_set, energy, force, input,&
                          rho_scp, rspw, para_env, blacs_env, ks_scp_env, & 
                          distribution_2d, eps_gvg_rspace, eps_rho_rspace, &
                          pscp, ks_scp, ascp_diag, error)

    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(aux_coeff_set_type), OPTIONAL, &
      POINTER                                :: aux_coeff_set
    TYPE(scp_energy_type), OPTIONAL, POINTER :: energy
    TYPE(scp_force_type), OPTIONAL, POINTER  :: force( : )
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    TYPE(scp_rho_type), OPTIONAL, POINTER    :: rho_scp
    TYPE(scp_rspw_type), OPTIONAL, POINTER   :: rspw
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env
    TYPE(qs_ks_scp_env_type), OPTIONAL, &
      POINTER                                :: ks_scp_env
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: distribution_2d
    REAL(dp), OPTIONAL                       :: eps_gvg_rspace, eps_rho_rspace
    TYPE(dbcsr_obj), OPTIONAL, POINTER       :: pscp, ks_scp
    REAL(dp), OPTIONAL, POINTER              :: ascp_diag( : )
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_scp_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    CPPrecondition(scp_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(blacs_env)) blacs_env => scp_env%blacs_env
    IF (PRESENT(aux_coeff_set)) aux_coeff_set => scp_env % aux_coeff_set
    IF (PRESENT(energy)) energy => scp_env%energy
    IF (PRESENT(force)) force => scp_env%force
    IF (PRESENT(input)) input=>scp_env%input
    IF (PRESENT(para_env)) para_env => scp_env%para_env
    IF (PRESENT(rho_scp)) rho_scp => scp_env%rho_scp
    IF (PRESENT(rspw)) rspw => scp_env%rspw
    IF (PRESENT(ks_scp_env)) ks_scp_env => scp_env%ks_scp_env
    IF (PRESENT(distribution_2d)) distribution_2d => scp_env%distribution_2d
    IF (PRESENT(eps_gvg_rspace )) eps_gvg_rspace = scp_env%eps_gvg_rspace
    IF (PRESENT(eps_rho_rspace )) eps_rho_rspace = scp_env%eps_rho_rspace
    IF (PRESENT( pscp )) pscp => scp_env%pscp
    IF (PRESENT( ks_scp )) ks_scp => scp_env%ks_scp
    IF (PRESENT( ascp_diag )) ascp_diag => scp_env%ascp_diag

  END SUBROUTINE get_scp_env

! *****************************************************************************
!> \brief   Initialise the SCP environment
!> \version 1.0
! *****************************************************************************
  SUBROUTINE init_scp_env ( scp_env, error )

    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    NULLIFY (scp_env%blacs_env)
    NULLIFY (scp_env%energy)
    NULLIFY (scp_env%force)
    NULLIFY (scp_env%input)
    NULLIFY (scp_env%rho_scp)
    NULLIFY (scp_env%rspw)
    NULLIFY (scp_env%ks_scp_env)
    NULLIFY (scp_env%para_env)
    NULLIFY (scp_env%blacs_env)
    NULLIFY (scp_env%distribution_2d)
    NULLIFY (scp_env%aux_coeff_set)
    NULLIFY (scp_env%pscp)
    NULLIFY (scp_env%ks_scp)
    NULLIFY (scp_env%ascp_diag)

    scp_env%ref_count=1
    last_scp_env_id_nr=last_scp_env_id_nr+1
    scp_env%id_nr=last_scp_env_id_nr

  END SUBROUTINE init_scp_env

! *****************************************************************************
!> \brief  set the SCP environment. 
!> \version 1.0
! *****************************************************************************
  SUBROUTINE set_scp_env( scp_env, aux_coeff_set, energy, force, ks_scp_env,&
                          input, rho_scp, rspw, para_env, blacs_env, &
                          distribution_2d, eps_gvg_rspace, eps_rho_rspace, &
                          ascp_diag, pscp, ks_scp, error)

    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(aux_coeff_set_type), OPTIONAL, &
      POINTER                                :: aux_coeff_set
    TYPE(scp_energy_type), OPTIONAL, POINTER :: energy
    TYPE(scp_force_type), OPTIONAL, POINTER  :: force( : )
    TYPE(qs_ks_scp_env_type), OPTIONAL, &
      POINTER                                :: ks_scp_env
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: input
    TYPE(scp_rho_type), OPTIONAL, POINTER    :: rho_scp
    TYPE(scp_rspw_type), OPTIONAL, POINTER   :: rspw
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_blacs_env_type), OPTIONAL, &
      POINTER                                :: blacs_env
    TYPE(distribution_2d_type), OPTIONAL, &
      POINTER                                :: distribution_2d
    REAL(dp), OPTIONAL                       :: eps_gvg_rspace, eps_rho_rspace
    REAL(dp), OPTIONAL, POINTER              :: ascp_diag( : )
    TYPE(dbcsr_obj), OPTIONAL, POINTER       :: pscp, ks_scp
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_scp_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure = .FALSE. 
    CPPrecondition(scp_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (PRESENT(aux_coeff_set)) THEN
       CALL aux_coeff_set_retain(aux_coeff_set, error=error)
       CALL aux_coeff_set_release(scp_env%aux_coeff_set, error=error)
       scp_env%aux_coeff_set => aux_coeff_set
    END IF
    IF (PRESENT(energy)) scp_env%energy => energy
    IF (PRESENT(force)) scp_env%force => force
    IF (PRESENT(rho_scp)) THEN ! accepts also null pointers !
      IF (ASSOCIATED(rho_scp)) CALL scp_rho_retain(rho_scp,error=error)
      CALL scp_rho_release(scp_env%rho_scp,error=error)
      scp_env%rho_scp => rho_scp
    END IF
    IF ( PRESENT ( rspw ) ) THEN ! accepts also null pointers !
      IF ( ASSOCIATED ( rspw ) ) CALL scp_rspw_retain ( rspw, error = error )
      CALL scp_rspw_release(scp_env%rspw,error=error)
      scp_env%rspw => rspw
    END IF
    IF (PRESENT(input)) THEN
       CALL section_vals_retain(input,error=error)
       CALL section_vals_release(scp_env%input,error=error)
       scp_env%input => input
    END IF
    IF (PRESENT(para_env)) THEN
       CALL cp_para_env_retain(para_env,error=error)
       CALL cp_para_env_release(scp_env%para_env,error=error)
       scp_env%para_env => para_env
    END IF
    IF (PRESENT(blacs_env)) THEN
       CALL cp_blacs_env_retain(blacs_env,error=error)
       CALL cp_blacs_env_release(scp_env%blacs_env,error=error)
       scp_env%blacs_env => blacs_env
    END IF
    IF (PRESENT(ks_scp_env)) THEN ! accept also null pointers?
       CALL qs_ks_scp_retain(ks_scp_env, error=error)
       CALL qs_ks_scp_release(scp_env%ks_scp_env, error=error)
       scp_env%ks_scp_env => ks_scp_env
    END IF
    IF (PRESENT(distribution_2d)) THEN
       CALL distribution_2d_retain(distribution_2d,error=error)
       CALL distribution_2d_release(scp_env%distribution_2d,error=error)
       scp_env%distribution_2d => distribution_2d
    ENDIF
    IF (PRESENT(eps_gvg_rspace)) THEN
       scp_env%eps_gvg_rspace = eps_gvg_rspace
    ENDIF
    IF (PRESENT(eps_rho_rspace)) THEN
       scp_env%eps_rho_rspace = eps_rho_rspace
    ENDIF
    IF (PRESENT(pscp)) THEN
       scp_env%pscp => pscp
    ENDIF
    IF (PRESENT(ks_scp)) THEN
       scp_env%ks_scp => ks_scp
    ENDIF
    IF (PRESENT(ascp_diag)) THEN
       scp_env%ascp_diag => ascp_diag
    ENDIF

  END SUBROUTINE set_scp_env

! *****************************************************************************
!> \brief allocates and intitializes a scp_env
!> \param scp_env the object to create
!> \param para_env the parallel environement for the scp_env
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
SUBROUTINE scp_env_create(scp_env,error)
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  ALLOCATE(scp_env, stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CALL init_scp_env(scp_env,error=error)
  END IF
END SUBROUTINE scp_env_create

! *****************************************************************************
!> \brief retains the given scp_env (see doc/ReferenceCounting.html)
!> \param scp_env the object to retain
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
SUBROUTINE scp_env_retain(scp_env,error)
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(scp_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(scp_env%ref_count>0,cp_failure_level,routineP,error,failure)
     scp_env%ref_count=scp_env%ref_count+1
  END IF
END SUBROUTINE scp_env_retain

! *****************************************************************************
!> \brief releases the given scp_env (see doc/ReferenceCounting.html)
!> \param scp_env the object to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
SUBROUTINE scp_env_release(scp_env,error)
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  IF (ASSOCIATED(scp_env)) THEN
     CPPrecondition(scp_env%ref_count>0,cp_failure_level,routineP,error,failure)
     scp_env%ref_count=scp_env%ref_count-1
     IF (scp_env%ref_count<1) THEN
        IF (ASSOCIATED(scp_env%energy)) &
             CALL deallocate_scp_energy(scp_env%energy)
        IF (ASSOCIATED(scp_env%force)) &
             CALL deallocate_scp_force(scp_env%force)
        IF (ASSOCIATED(scp_env%distribution_2d)) THEN
           CALL distribution_2d_release(scp_env%distribution_2d,error=error)
        ENDIF
        IF (ASSOCIATED(scp_env%pscp)) THEN
           CALL dbcsr_deallocate_matrix ( scp_env%pscp, error )
        ENDIF
        IF (ASSOCIATED(scp_env%ks_scp)) THEN
           CALL dbcsr_deallocate_matrix ( scp_env%ks_scp, error )
        ENDIF
        IF (ASSOCIATED(scp_env%ascp_diag)) THEN
           DEALLOCATE ( scp_env%ascp_diag, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        ENDIF
        CALL scp_rho_release(scp_env%rho_scp, error=error)
        CALL scp_rspw_release(scp_env%rspw, error=error)
        CALL qs_ks_scp_release(scp_env%ks_scp_env, error=error)
        CALL aux_coeff_set_release(scp_env%aux_coeff_set,error=error)
        CALL section_vals_release(scp_env%input,error=error)
        CALL cp_para_env_release(scp_env%para_env,error=error)
        CALL cp_blacs_env_release(scp_env%blacs_env, error=error)
        DEALLOCATE(scp_env, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(scp_env)
END SUBROUTINE scp_env_release

END MODULE scp_environment_types
