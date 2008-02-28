!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author CJM
! *****************************************************************************
MODULE qs_ks_scp_methods
  USE cell_types,                      ONLY: cell_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_retain,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_fft_wrap,&
                                             pw_integral_ab,&
                                             pw_p_type,&
                                             pw_zero
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_ks_scp_types,                 ONLY: qs_ks_scp_env_type
  USE scp_energy_types,                ONLY: scp_energy_type
  USE scp_environment_types,           ONLY: get_scp_env,&
                                             scp_environment_type
  USE scp_rho_types,                   ONLY: scp_rho_get,&
                                             scp_rho_type
  USE scp_rspw_types,                  ONLY: scp_rspw_get,&
                                             scp_rspw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_scp_methods'
  INTEGER, SAVE, PRIVATE :: last_ks_scp_id_nr=0

  PUBLIC :: qs_ks_scp_create,  scp_calculate_gpw_energy, scp_calculate_qmmm_energy,&
            scp_modify_hartree_pot, qs_ks_scp_did_change, qs_ks_scp_update

CONTAINS

! *****************************************************************************
!> \brief allocates and initializes the given ks_scp_env.
!> \param ks_scp_env the ks_scp env to be initialized
!> \param qs_env the qs environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      04.2007 
!> \author CJM
! *****************************************************************************
  SUBROUTINE qs_ks_scp_create ( ks_scp_env, qs_env, error )
    TYPE(qs_ks_scp_env_type), POINTER        :: ks_scp_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_scp_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: auxbas_grid, handle, stat
    LOGICAL                                  :: failure
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    CALL timeset ( routineN, 'I', '', handle )

    failure=.FALSE.
    CPPreconditionNoFail(.NOT.ASSOCIATED(ks_scp_env),cp_failure_level,routineP,error)
    ALLOCATE(ks_scp_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    IF (.not.failure) THEN
       NULLIFY( pw_env, auxbas_pw_pool,&
            pw_env, cell,pw_pools)
       CALL get_qs_env(qs_env=qs_env,&
            pw_env=pw_env, cell=cell,error=error)
       CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,pw_pools=pw_pools,&
            auxbas_grid=auxbas_grid,error=error)

       ks_scp_env%rho_changed=.TRUE.
       ks_scp_env%potential_changed=.TRUE.
       ks_scp_env%full_reset=.TRUE.
       ks_scp_env%n_evals=0
       ks_scp_env%ref_count=1
       last_ks_scp_id_nr=last_ks_scp_id_nr+1
       ks_scp_env%id_nr=last_ks_scp_id_nr

    END IF
    IF (.NOT.failure) THEN
       ks_scp_env%auxbas_pw_pool=>auxbas_pw_pool
       CALL pw_pool_retain(ks_scp_env%auxbas_pw_pool,error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool,ks_scp_env%v_scp_rspace%pw,&
            use_data=REALDATA3D, in_space=REALSPACE,error=error)
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_ks_scp_create

! *****************************************************************************
!> \brief Computes the contribution to the hartree energy due to the SCP
!>      electrostatic coupling
!> \par History
!>      12.2007 created
!> \author CJM
! *****************************************************************************
  SUBROUTINE scp_calculate_qmmm_energy(scp_env, v_qmmm, error)
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(pw_p_type), INTENT(IN)              :: v_qmmm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_calculate_qmmm_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(pw_p_type), POINTER                 :: rhop_r
    TYPE(scp_energy_type), POINTER           :: energy
    TYPE(scp_rho_type), POINTER              :: rho_scp

    CALL timeset(routineN,'I','',handle)
    NULLIFY ( rho_scp, rhop_r, energy )
   
    CALL get_scp_env(scp_env=scp_env,&
                     rho_scp = rho_scp, & 
                     energy = energy, &
                     error=error)
    CALL scp_rho_get ( rho_scp, rhop_r = rhop_r, error = error )

    energy % e_scp_qmmm = pw_integral_ab ( rhop_r%pw, v_qmmm%pw )
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE scp_calculate_qmmm_energy
! *****************************************************************************
  SUBROUTINE scp_calculate_gpw_energy(scp_env, v, error)
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(pw_p_type), INTENT(IN)              :: v
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_calculate_gpw_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(pw_p_type), POINTER                 :: rhop_g
    TYPE(scp_energy_type), POINTER           :: energy
    TYPE(scp_rho_type), POINTER              :: rho_scp

    CALL timeset(routineN,'I','',handle)
    NULLIFY ( rho_scp, rhop_g, energy )
   
    CALL get_scp_env(scp_env=scp_env,&
                     rho_scp = rho_scp, & 
                     energy = energy, &
                     error=error)
    CALL scp_rho_get ( rho_scp, rhop_g = rhop_g, error = error )
     
    energy % e_scp_ks = pw_integral_ab ( rhop_g % pw, v%pw )

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE scp_calculate_gpw_energy

! *****************************************************************************
!> \brief Modify the hartree potential in order to include the SCP correction
!> \par History
!>      04.2007 
!> \author CJM
! *****************************************************************************
  SUBROUTINE scp_modify_hartree_pot(v_hartree,v_scp, error)
    TYPE(pw_p_type), INTENT(INOUT)           :: v_hartree
    TYPE(pw_p_type), INTENT(IN)              :: v_scp
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_modify_hartree_pot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,'I','',handle)
    failure=.FALSE.

    v_hartree%pw%cr3d = v_hartree%pw%cr3d + &
                        v_scp%pw%pw_grid%dvol*v_scp%pw%cr3d

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE scp_modify_hartree_pot

! *****************************************************************************
!> \brief tells that some of the things relevant to the ks_scp calculation
!>      did change. Has to be called when changes happen otherwise
!>      the calculation will give wrong results.
!> \param ks_scp_env the environement that is informed about the changes
!> \param rho_changed if true it means that the density has changed
!> \param potential_changed if true it means that the potential has changed
!> \param full_reset if true everything has changed
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      4.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE qs_ks_scp_did_change(qs_env,rho_changed,potential_changed,full_reset,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(in), OPTIONAL            :: rho_changed, &
                                                potential_changed, full_reset
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_scp_did_change', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(qs_ks_scp_env_type), POINTER        :: ks_scp_env
    TYPE(scp_environment_type), POINTER      :: scp_env

  failure=.FALSE.

  CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, error = error ) 
  CALL get_scp_env ( scp_env = scp_env, ks_scp_env = ks_scp_env, error = error ) 

  CPPrecondition(ASSOCIATED(ks_scp_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ks_scp_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(rho_changed)) THEN
        IF (rho_changed) ks_scp_env%rho_changed=.TRUE.
     END IF
     IF (PRESENT(potential_changed)) THEN
        IF (potential_changed) ks_scp_env%potential_changed=.TRUE.
     END IF
     IF (PRESENT(full_reset)) THEN
        IF (full_reset) THEN
           ks_scp_env%potential_changed=.TRUE.
        END IF
     END IF
  END IF
END SUBROUTINE qs_ks_scp_did_change

! *****************************************************************************
!> \brief Updates the qs_ks_scp_env if necessary
!> \author CJM
! *****************************************************************************
SUBROUTINE qs_ks_scp_update ( qs_env, just_energy, error )

    TYPE(qs_environment_type), POINTER       :: qs_env
    LOGICAL, INTENT(IN), OPTIONAL            :: just_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_scp_update', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, stat, output_unit
    LOGICAL                                  :: energy_only, failure, ionode
    REAL(KIND=dp)                            :: t1
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), POINTER                 :: rhop_g, vhscp_gspace
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_ks_scp_env_type), POINTER        :: ks_scp_env
    TYPE(scp_energy_type), POINTER           :: energy
    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(scp_rho_type), POINTER              :: rho_scp
    TYPE(scp_rspw_type), POINTER             :: rspw

  CALL timeset ( 'qs_ks_scp_update','I',' ',handle )

  failure = .FALSE.
  energy_only = .FALSE.
  IF ( PRESENT ( just_energy ) ) energy_only = just_energy
  NULLIFY(scp_env)
  NULLIFY(poisson_env,auxbas_pw_pool)
  NULLIFY(pw_env, rspw)
  NULLIFY(ks_scp_env)
  para_env=>qs_env%para_env
  logger => cp_error_get_logger(error)
  ionode=para_env%ionode
  output_unit=cp_logger_get_default_unit_nr(logger)
  stat = 0
  t1 = m_walltime()

  CALL get_qs_env ( qs_env = qs_env, scp_env = scp_env, error=error )
  CALL get_scp_env ( scp_env = scp_env, rspw = rspw, ks_scp_env = ks_scp_env, error = error )
  CALL scp_rspw_get ( scp_rspw=rspw, qs_pw_env = pw_env, error = error ) 

! Could also get pw_env direct from qs_env. But, in general,
! SCP should have its own pw_env (for example with using SCP
! for classical polarizable models in conjunction with FIST)

  CPPrecondition(ASSOCIATED(ks_scp_env),cp_failure_level,routineP,error,failure)
  CPPrecondition(ks_scp_env%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure ) THEN
    IF ( ks_scp_env % rho_changed .OR. ks_scp_env % potential_changed .OR. &
         ks_scp_env % full_reset ) THEN
! Proceed with update
      CALL pw_env_get ( pw_env = pw_env, &
                        auxbas_pw_pool = auxbas_pw_pool, &
                        poisson_env = poisson_env, &
                        error = error )
       ALLOCATE(vhscp_gspace,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pw_pool_create_pw(auxbas_pw_pool, vhscp_gspace%pw, &
                               use_data = COMPLEXDATA1D,&
                               in_space = RECIPROCALSPACE, error=error)
       CALL pw_zero(vhscp_gspace%pw) 

! Getting the energy and density structures
       CALL get_scp_env ( scp_env=scp_env,&
                          rho_scp=rho_scp, energy=energy, &
                          error=error)
       CALL scp_rho_get ( rho_scp = rho_scp, rhop_g=rhop_g, error=error)

! calculate V_HSCP on the realspace grid for later use with the KS matrix
       CALL pw_poisson_solve(poisson_env,rhop_g%pw, energy % e_scp_hartree, &
                             vhscp_gspace%pw,error=error)

! We know and updated KS_SCP_ENV and have V_HSCP in real-space
! Volume factors will get added in SCP_MODIFY_HARTEE_POTENTIAL
       CALL pw_fft_wrap(vhscp_gspace%pw,ks_scp_env%v_scp_rspace%pw)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,vhscp_gspace%pw,error=error)
       DEALLOCATE(vhscp_gspace,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

       IF(.NOT.energy_only) THEN
          ks_scp_env%rho_changed=.FALSE.
          ks_scp_env%potential_changed=.FALSE.
          ks_scp_env%full_reset=.FALSE.
       END IF
    END IF   
  END IF
  t1 = m_walltime()
  CALL timestop(0.0_dp,handle)
END SUBROUTINE qs_ks_scp_update

END MODULE qs_ks_scp_methods
