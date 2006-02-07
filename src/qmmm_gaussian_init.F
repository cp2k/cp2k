!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2004  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/qmmm_gaussian_init [1.0] *
!!
!!   NAME
!!     qmmm_gaussian_init
!!
!!   FUNCTION
!!     Initialize the use of the gaussians to treat the QMMM 
!!     coupling potential
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Teodoro Laino
!!
!!   MODIFICATION HISTORY
!!     6.2004 created [tlaino]
!!
!!   SOURCE
!****************************************************************************
MODULE qmmm_gaussian_init
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_should_output
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: do_qmmm_gauss
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_p_type
  USE qmmm_gaussian_data,              ONLY: max_geep_lib_gauss,&
                                             min_geep_lib_gauss
  USE qmmm_gaussian_input,             ONLY: read_mm_potential,&
                                             set_mm_potential
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
  USE qs_util,                         ONLY: exp_radius
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_gaussian_init'
  PUBLIC :: qmmm_gaussian_initialize

CONTAINS

!!****f* qmmm_gaussian_init/qmmm_gaussian_initialize [1.0] *
!!
!!   NAME
!!     qmmm_gaussian_initialize
!!
!!   FUNCTION
!!     Initialize the Gaussian QMMM Environment
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
!!     06.2004 created [tlaino]
!!
!!*** **********************************************************************  
  SUBROUTINE qmmm_gaussian_initialize(qmmm_gaussian_fns, globenv, pw_env,&
       qm_cell_small, mm_cell, mm_el_pot_radius, mm_el_pot_radius_corr,  &
       qmmm_coupl_type, eps_mm_rspace, maxradius, maxchrg, compatibility,&
       print_section, qmmm_section, error)
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: qmmm_gaussian_fns
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(cell_type), POINTER                 :: qm_cell_small, mm_cell
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_el_pot_radius, &
                                                mm_el_pot_radius_corr
    INTEGER, INTENT(IN)                      :: qmmm_coupl_type
    REAL(KIND=dp), INTENT(IN)                :: eps_mm_rspace
    REAL(KIND=dp), DIMENSION(:), POINTER     :: maxradius
    REAL(KIND=dp), INTENT(IN)                :: maxchrg
    LOGICAL, INTENT(IN)                      :: compatibility
    TYPE(section_vals_type), POINTER         :: print_section, qmmm_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_gaussian_initialize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ilevel, j, Ndim, &
                                                num_geep_gauss, stat
    LOGICAL                                  :: Failure, Found, use_geep_lib
    REAL(KIND=dp)                            :: alpha, mymaxradius, Prefactor
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c_radius, radius
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(qmmm_gaussian_type), POINTER        :: mypgf
    TYPE(section_vals_type), POINTER         :: print_gen

! Statements

    NULLIFY(mypgf, gridlevel_info, radius, c_radius, print_gen, logger)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       logger => cp_error_get_logger(error)
       print_gen => section_vals_get_subs_vals(print_section,"PROGRAM_RUN_INFO",error=error)
       CALL section_vals_val_get(qmmm_section,"USE_GEEP_LIB",i_val=num_geep_gauss,error=error)
       IF (num_geep_gauss == 0) THEN
          use_geep_lib = .FALSE.
       ELSE
          use_geep_lib = .TRUE.
          CPPostcondition(num_geep_gauss>=min_geep_lib_gauss,cp_failure_level,routineP,error,failure)
          CPPostcondition(num_geep_gauss<=max_geep_lib_gauss,cp_failure_level,routineP,error,failure)          
       END IF
       SELECT CASE(qmmm_coupl_type)
       CASE(do_qmmm_gauss)
          !
          ! Preprocessing...
          !
          ALLOCATE(radius(1),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
          ALLOCATE(c_radius(1),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
          Ndim = SIZE(radius)
          Loop_on_all_values: DO I = 1, SIZE(mm_el_pot_radius)
             Found=.FALSE.
             Loop_on_found_values: DO J = 1, SIZE(radius)-1
                IF (mm_el_pot_radius(i).EQ.radius(j)) THEN
                   Found=.TRUE.
                   EXIT Loop_on_found_values
                END IF
             END DO Loop_on_found_values
             IF (.NOT.Found) THEN
                Ndim = SIZE(radius)
                radius(Ndim)   = mm_el_pot_radius(i)
                c_radius(Ndim) = mm_el_pot_radius_corr(i)
                Ndim = Ndim + 1
                CALL REALLOCATE(radius,1,Ndim)
                CALL REALLOCATE(c_radius,1,Ndim)
             END IF
          END DO Loop_on_all_values
          !
          IF (Ndim-1 > 0) THEN
             CALL REALLOCATE(radius,1,Ndim-1)
             CALL REALLOCATE(c_radius,1,Ndim-1)
          ELSE IF (Ndim-1 == 0) THEN
             DEALLOCATE(radius,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
             DEALLOCATE(c_radius,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
          ELSE
             CPPostcondition(.FALSE.,cp_failure_level,routineP,error,Failure)      
          END IF
          !
          ALLOCATE(qmmm_gaussian_fns(Ndim-1),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)       
          DO I = 1, Ndim-1
             NULLIFY(qmmm_gaussian_fns(I)%pgf)
             ALLOCATE(qmmm_gaussian_fns(I)%pgf, stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)
             NULLIFY(qmmm_gaussian_fns(I)%pgf%Ak)
             NULLIFY(qmmm_gaussian_fns(I)%pgf%Gk)
             NULLIFY(qmmm_gaussian_fns(I)%pgf%grid_level)
             ! 
             ! Default Values
             !
             qmmm_gaussian_fns(I)%pgf%Elp_Radius          = radius(I)
             qmmm_gaussian_fns(I)%pgf%Elp_Radius_corr     = c_radius(I)
          END DO
          IF (ASSOCIATED(radius)) THEN 
             DEALLOCATE(radius,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)      
          END IF
          IF (ASSOCIATED(c_radius)) THEN 
             DEALLOCATE(c_radius,stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,Failure)      
          END IF
          !
          IF (use_geep_lib) THEN
             CALL  set_mm_potential(globenv, qmmm_gaussian_fns, qm_cell_small,&
                  compatibility, num_geep_gauss, error)
          ELSE
             CALL read_mm_potential(globenv, qmmm_gaussian_fns, qm_cell_small,&
                  compatibility, error)
          END IF
          !
          CALL pw_env_get(pw_env,pw_pools=pools, gridlevel_info=gridlevel_info,&
               error=error)
          ALLOCATE(maxradius(SIZE(pools)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          maxradius = 0.0_dp
          DO J = 1, SIZE(qmmm_gaussian_fns)
             mypgf => qmmm_gaussian_fns(J)%pgf
             ALLOCATE(mypgf%grid_level(SIZE(mypgf%Ak)),stat=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             mypgf%grid_level=0
             mymaxradius = 0.0_dp
             DO I = 1, mypgf%number_of_gaussians
                IF ( mypgf%Gk(I) /= 0.0_dp) THEN
                   alpha = 1.0_dp / mypgf%Gk(I)
                   alpha = alpha*alpha
                   ilevel=gaussian_gridlevel(gridlevel_info,alpha)
                   Prefactor = mypgf%Ak(I)*maxchrg
                   mymaxradius=exp_radius( 0, alpha, eps_mm_rspace, Prefactor)
                   maxradius(ilevel)=MAX(maxradius(ilevel),mymaxradius)
                   mypgf%grid_level(i)=ilevel
                END IF
             END DO
          END DO
          !
          ! End of gaussian initialization...
       CASE DEFAULT
          IF (BTEST(cp_print_key_should_output(logger%iter_info,print_gen,error=error),cp_p_file).AND.&
               logger%para_env%mepos==logger%para_env%source) &
          WRITE(*,'(A)')" QMMM Gaussian Data Not Initialized!"
       END SELECT
    END IF
  END SUBROUTINE qmmm_gaussian_initialize


END MODULE qmmm_gaussian_init
