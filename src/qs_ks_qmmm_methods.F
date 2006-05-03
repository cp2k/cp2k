
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2003  CP2K developers group                            !
!-----------------------------------------------------------------------------!
!!****** cp2k/qs_ks_qmmm_methods [1.0] *
!!
!!   NAME
!!     qs_ks_qmmm_methods
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     T.Laino
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE qs_ks_qmmm_methods
  USE coefficient_types,               ONLY: coeff_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cube_utils,                      ONLY: cube_info_type,&
                                             init_cube_info
  USE input_constants,                 ONLY: do_qmmm_gauss
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_retain
  USE pw_pool_types,                   ONLY: pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: REALDATA3D,&
                                             REALSPACE,&
                                             pw_integral_ab
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
  !!****f* force_env_methods/ks_qmmm_env_rebuild [1.0] *
  !!
  !!   NAME
  !!     ks_qmmm_env_rebuild
  !!
  !!   FUNCTION
  !!     Initialize the ks_qmmm_env
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Teodoro Laino
  !!
  !!   MODIFICATION HISTORY
  !!     05.2004 created [tlaino]
  !!
  !!*** **********************************************************************
  SUBROUTINE ks_qmmm_env_rebuild(qs_env,qmmm_env,error)
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env

    NULLIFY(ks_qmmm_env)
    CALL get_qs_env(qs_env=qs_env,&
         ks_qmmm_env=ks_qmmm_env)

    !   *** allocate the ks_qmmm env if not allocated yet!**
    IF (.NOT.ASSOCIATED(ks_qmmm_env)) THEN
       CALL qs_ks_qmmm_create(ks_qmmm_env=ks_qmmm_env,qs_env=qs_env,&
            qmmm_env=qmmm_env,error=error)
       CALL set_qs_env(qs_env=qs_env, ks_qmmm_env=ks_qmmm_env)
       CALL qs_ks_qmmm_release(ks_qmmm_env=ks_qmmm_env,force_env_section=qs_env%input,&
            error=error)
    END IF
  END SUBROUTINE ks_qmmm_env_rebuild
  !***************************************************************************

!!****f* qs_ks_methods/qs_ks_qmmm_create [1.0] *
!!
!!   NAME
!!     qs_ks_qmmm_create
!!
!!   FUNCTION
!!     allocates and initializes the given ks_qmmm_env.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_qmmm_env: the ks_qmmm env to be initialized
!!     - qs_env: the qs environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE qs_ks_qmmm_create(ks_qmmm_env, qs_env, qmmm_env, error)
    TYPE(qs_ks_qmmm_env_type), POINTER       :: ks_qmmm_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qs_ks_qmmm_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: cmax, handle, igrid, stat
    LOGICAL                                  :: failure
    TYPE(cube_info_type), DIMENSION(:), &
      POINTER                                :: cube_info
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    CPPreconditionNoFail(.NOT.ASSOCIATED(ks_qmmm_env),cp_failure_level,routineP,error)
    ALLOCATE(ks_qmmm_env, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    NULLIFY( ks_qmmm_env%pw_env,&
             ks_qmmm_env%cube_info)
    IF (.not.failure) THEN
       NULLIFY( auxbas_pw_pool)
       CALL get_qs_env(qs_env=qs_env,&
            pw_env=ks_qmmm_env%pw_env)
       CALL pw_env_get(ks_qmmm_env%pw_env,auxbas_pw_pool=auxbas_pw_pool)
       CALL pw_env_retain(ks_qmmm_env%pw_env,error=error)

       ks_qmmm_env%n_evals=0
       ks_qmmm_env%ref_count=1
       last_ks_qmmm_nr=last_ks_qmmm_nr+1
       ks_qmmm_env%id_nr=last_ks_qmmm_nr

       CALL pw_pool_init_coeff(auxbas_pw_pool,ks_qmmm_env%v_qmmm_rspace,&
            use_data=REALDATA3D, in_space=REALSPACE)
  
       IF (qmmm_env%qmmm_coupl_type.EQ.do_qmmm_gauss) THEN
          CALL pw_env_get(ks_qmmm_env%pw_env,pw_pools=pools,error=error)
          ALLOCATE(cube_info(SIZE(pools)), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO igrid=1,SIZE(pools)
             CALL init_cube_info(cube_info(igrid),&
               pools(igrid)%pool%pw_grid%dr(:),&
               pools(igrid)%pool%pw_grid%dh(:,:),&
               pools(igrid)%pool%pw_grid%orthorhombic,&
               qmmm_env%maxRadius(igrid),cmax)
          END DO
          ks_qmmm_env%cube_info => cube_info
       END IF
       !
    END IF

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE qs_ks_qmmm_create
!***************************************************************************
!!****f* qmmm_methods/qmmm_calculate_energy [1.0] *
!!
!!   NAME
!!     qmmm_calculate_energy
!!
!!   FUNCTION
!!     Computes the contribution to the total energy of the QM/MM 
!!     electrostatic coupling
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE qmmm_calculate_energy(qs_env, rho, v_qmmm, qmmm_energy, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(coeff_type), DIMENSION(:), POINTER  :: rho
    TYPE(coeff_type), INTENT(IN)             :: v_qmmm
    REAL(KIND=dp), INTENT(INOUT)             :: qmmm_energy
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_calculate_energy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: failure
    TYPE(dft_control_type), POINTER          :: dft_control

    CALL timeset(routineN,"I","",handle)    
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    NULLIFY(dft_control)

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    error=error)

    qmmm_energy = 0.0_dp
    DO ispin=1, dft_control%nspins
       qmmm_energy = qmmm_energy + pw_integral_ab ( rho(ispin)%pw, v_qmmm%pw )
    END DO

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_calculate_energy

!!****f* qmmm_methods/qmmm_modify_hartree_pot [1.0] *
!!
!!   NAME
!!     qmmm_modify_hartree_pot
!!
!!   FUNCTION
!!     Modify the hartree potential in order to include the QM/MM correction
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     05.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE qmmm_modify_hartree_pot(v_hartree,v_qmmm, error)
    TYPE(coeff_type), INTENT(INOUT)          :: v_hartree
    TYPE(coeff_type), INTENT(IN)             :: v_qmmm
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_modify_hartree_pot', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,"I","",handle)    
    failure=.FALSE.

    v_hartree%pw%cr3d = v_hartree%pw%cr3d + &
                        v_qmmm%pw%pw_grid%dvol*v_qmmm%pw%cr3d

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE qmmm_modify_hartree_pot

END MODULE qs_ks_qmmm_methods
