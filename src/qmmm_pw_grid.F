!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                          !
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"
!!****s* cp2k/qmmm_pw_grid [1.1] *
!!
!!   NAME
!!     qmmm_pw_grid
!!
!!   FUNCTION
!!     sets variables for the qmmm pool of pw_types
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
MODULE qmmm_pw_grid
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE global_types,                    ONLY: global_environment_type
  USE particle_list_types,             ONLY: particle_list_create,&
                                             particle_list_release,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE memory_utilities,                ONLY: reallocate
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE coefficient_types,               ONLY: coeff_add,&
                                             coeff_init,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE qmmm_types,                      ONLY: qmmm_env_qm_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_retain,&
                                             pw_pool_type,&
                                             pw_pools_dealloc,&
                                             pw_pool_create
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: pw_grid_type, &
                                             PW_MODE_LOCAL
  USE pw_grids,                        ONLY: pw_grid_create,&
                                             pw_grid_release,&
                                             pw_grid_create_copy_no_pbc
  USE qs_ks_qmmm_types,                ONLY: qs_ks_qmmm_env_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE mathconstants,                   ONLY: pi
  USE simulation_cell,                 ONLY: cell_type
  USE input_constants, ONLY: do_qmmm_grid, do_qmmm_gauss, do_qmmm_none
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: qmmm_pw_grid_init
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qmmm_pw_grid'

CONTAINS

!!****f* qmmm_pw_grid/qmmm_pw_grid_init [1.0] *
!!
!!   NAME
!!     qmmm_pw_grid_init
!!
!!   FUNCTION
!!     Initialize the qmmm pool of pw_type.
!!     Then Main difference w.r.t. QS pw_type pools is that this pool
!!     has [0,L] as boundaries.
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
!!     08.2004 created [tlaino]
!!
!!*** **********************************************************************
  SUBROUTINE qmmm_pw_grid_init(qmmm_env, pw_env, cell, para_env, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(qmmm_env_qm_type), INTENT(inout), &
         TARGET                               :: qmmm_env  
    TYPE(pw_env_type), POINTER                :: pw_env  
    TYPE(cp_para_env_type), POINTER           :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                                :: error
    ! Local Variables
    TYPE(pw_pool_p_type), DIMENSION(:), &
         POINTER                              :: pw_pools
    TYPE(pw_pool_type), POINTER               :: pool
    TYPE(cell_type), POINTER                  :: cell
    LOGICAL :: failure
    INTEGER :: auxbas_grid, stat
    INTEGER :: Ilevel
    TYPE (pw_grid_type), POINTER :: el_struct
    REAL(KIND=dp), DIMENSION(3) :: Maxdr, Mindr
    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_el_coupling', &
         routineP = moduleN//':'//routineN

    NULLIFY(el_struct)
    failure = .false.
    Maxdr   = TINY(0.0_dp)
    Mindr   = HUGE(0.0_dp)
    IF (qmmm_env%qmmm_coupl_type == do_qmmm_gauss) THEN 
       CALL pw_env_get(pw_env=pw_env,&
                       pw_pools=pw_pools,&
                       auxbas_grid=auxbas_grid,&
                       error=error)
       !
       IF (ASSOCIATED(qmmm_env%aug_pools)) THEN
          CALL pw_pools_dealloc(qmmm_env%aug_pools,error=error)
       END IF
       ALLOCATE(qmmm_env%aug_pools(SIZE(pw_pools)),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       !
       DO Ilevel = 1, SIZE(pw_pools)
          NULLIFY(pool, qmmm_env%aug_pools(Ilevel)%pool)
          pool      => pw_pools(Ilevel)%pool
          nullify(el_struct)
          CALL pw_grid_create_copy_no_pbc(pool%pw_grid,el_struct, cell=cell,&
               error=error)

          CALL pw_pool_create(qmmm_env%aug_pools(Ilevel)%pool,&
               pw_grid=el_struct,error=error)

          Maxdr = MAX(Maxdr,el_struct%dr)
          Mindr = MIN(Mindr,el_struct%dr)
          IF (ALL(Maxdr.EQ.el_struct%dr)) qmmm_env%gridlevel_info%coarser_grid = Ilevel
          IF (ALL(Mindr.EQ.el_struct%dr)) qmmm_env%gridlevel_info%auxbas_grid  = Ilevel

          CALL pw_grid_release(el_struct,error=error)
          
       END DO
    END IF

  END SUBROUTINE qmmm_pw_grid_init

END MODULE qmmm_pw_grid
