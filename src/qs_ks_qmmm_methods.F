!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author T.Laino
! *****************************************************************************
MODULE qs_ks_qmmm_methods
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cube_utils,                      ONLY: cube_info_type,&
                                             init_cube_info
  USE f77_blas
  USE input_constants,                 ONLY: do_qmmm_gauss,&
                                             do_qmmm_swave
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_retain
  USE pw_methods,                      ONLY: pw_integral_ab
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_p_type
  USE qmmm_types,                      ONLY: qmmm_env_qm_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type,&
                                             qs_ks_qmmm_release
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ks_qmmm_methods'
  INTEGER, SAVE, PRIVATE :: last_ks_qmmm_nr=0

  PUBLIC :: ks_qmmm_env_rebuild, qs_ks_qmmm_create,qmmm_calculate_energy,&
            qmmm_modify_hartree_pot

CONTAINS

! *****************************************************************************
!> \brief Initialize the ks_qmmm_env
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE ks_qmmm_env_rebuild(qs_env,qmmm_env,error)
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env

    NULLIFY(ks_qmmm_env)
    CALL get_qs_env(qs_env=qs_env,&
         ks_qmmm_env=ks_qmmm_env,error=error)

    !   *** allocate the ks_qmmm env if not allocated yet!**
    IF (.NOT.ASSOCIATED(ks_qmmm_env)) THEN
       CALL qs_ks_qmmm_create(ks_qmmm_env=ks_qmmm_env,qs_env=qs_env,&
            qmmm_env=qmmm_env,error=error)
       CALL set_qs_env(qs_env=qs_env, ks_qmmm_env=ks_qmmm_env,error=error)
       CALL qs_ks_qmmm_release(ks_qmmm_env=ks_qmmm_env, error=error)
    END IF
  END SUBROUTINE ks_qmmm_env_rebuild

! *****************************************************************************
!> \brief allocates and initializes the given ks_qmmm_env.
!> \param ks_qmmm_env the ks_qmmm env to be initialized
!> \param qs_env the qs environment
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qs_ks_qmmm_create(ks_qmmm_env, qs_env, qmmm_env, error)
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_qmmm_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, igrid, stat
    LOGICAL                                  :: failure
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    CALL timeset(routineN,handle)

    failure=.FALSE.
    CPPreconditionNoFail(.NOT.ASSOCIATED(ks_qmmm_env),cp_failure_level,routineP,error)
    ALLOCATE(ks_qmmm_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    NULLIFY( ks_qmmm_env%pw_env,&
             ks_qmmm_env%cube_info)
    IF (.not.failure) THEN
       NULLIFY( auxbas_pw_pool)
       CALL get_qs_env(qs_env=qs_env,&
            pw_env=ks_qmmm_env%pw_env,error=error)
       CALL pw_env_get(ks_qmmm_env%pw_env,auxbas_pw_pool=auxbas_pw_pool,error=error)
       CALL pw_env_retain(ks_qmmm_env%pw_env,error=error)

       ks_qmmm_env%n_evals=0
       ks_qmmm_env%ref_count=1
       last_ks_qmmm_nr=last_ks_qmmm_nr+1
       ks_qmmm_env%id_nr=last_ks_qmmm_nr

       CALL pw_pool_create_pw(auxbas_pw_pool,ks_qmmm_env%v_qmmm_rspace%pw,&
            use_data=REALDATA3D, in_space=REALSPACE,error=error)

       IF ((qmmm_env%qmmm_coupl_type==do_qmmm_gauss).OR.(qmmm_env%qmmm_coupl_type==do_qmmm_swave)) THEN
          CALL pw_env_get(ks_qmmm_env%pw_env,pw_pools=pools,error=error)
          ALLOCATE(cube_info(SIZE(pools)), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO igrid=1,SIZE(pools)
             CALL init_cube_info(cube_info(igrid),&
               pools(igrid)%pool%pw_grid%dr(:),&
               pools(igrid)%pool%pw_grid%dh(:,:),&
               pools(igrid)%pool%pw_grid%dh_inv(:,:),&
               pools(igrid)%pool%pw_grid%orthorhombic,&
               qmmm_env%maxRadius(igrid))
          END DO
          ks_qmmm_env%cube_info => cube_info
       END IF
       NULLIFY(ks_qmmm_env%matrix_h)
       !
    END IF

    CALL timestop(handle)

  END SUBROUTINE qs_ks_qmmm_create

! *****************************************************************************
!> \brief Computes the contribution to the total energy of the QM/MM
!>      electrostatic coupling
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_calculate_energy(qs_env, rho, v_qmmm, qmmm_energy, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho
    TYPE(pw_p_type), INTENT(IN)              :: v_qmmm
    REAL(KIND=dp), INTENT(INOUT)             :: qmmm_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_calculate_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, output_unit
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_p_type), POINTER                 :: rho0_s_rs
    TYPE(section_vals_type), POINTER         :: input

    CALL timeset(routineN,handle)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    NULLIFY(dft_control,input,logger)
    logger => cp_error_get_logger(error)

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    input=input,&
                    error=error)

    output_unit=cp_print_key_unit_nr(logger,input,"QMMM%PRINT%PROGRAM_RUN_INFO",&
               extension=".qmmmLog",error=error)
    IF (output_unit>0)  WRITE (UNIT=output_unit,FMT="(T3,A)")&
               "Adding QM/MM electrostatic potential to the Kohn-Sham potential."


    qmmm_energy = 0.0_dp
    DO ispin=1, dft_control%nspins
       qmmm_energy = qmmm_energy + pw_integral_ab ( rho(ispin)%pw, v_qmmm%pw, error=error)
    END DO
    IF(dft_control%qs_control%gapw) THEN
      CALL get_qs_env(qs_env=qs_env,&
                    rho0_s_rs=rho0_s_rs,&
                    error=error)
       CPPrecondition(ASSOCIATED(rho0_s_rs),cp_failure_level,routineP,error,failure)
       qmmm_energy = qmmm_energy + pw_integral_ab (rho0_s_rs%pw, v_qmmm%pw, error=error)
    END IF

    CALL cp_print_key_finished_output(output_unit,logger,input,&
               "QMMM%PRINT%PROGRAM_RUN_INFO", error=error)

    CALL timestop(handle)
  END SUBROUTINE qmmm_calculate_energy

! *****************************************************************************
!> \brief Modify the hartree potential in order to include the QM/MM correction
!> \par History
!>      05.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_modify_hartree_pot(v_hartree,v_qmmm,scale,error)
    TYPE(pw_p_type), INTENT(INOUT)           :: v_hartree
    TYPE(pw_p_type), INTENT(IN)              :: v_qmmm
    REAL(KIND=dp), INTENT(IN)                :: scale
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_modify_hartree_pot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: fac

    CALL timeset(routineN,handle)
    failure=.FALSE.

    fac = v_qmmm%pw%pw_grid%dvol*scale
    v_hartree%pw%cr3d = v_hartree%pw%cr3d + fac*v_qmmm%pw%cr3d

    CALL timestop(handle)
  END SUBROUTINE qmmm_modify_hartree_pot

END MODULE qs_ks_qmmm_methods
