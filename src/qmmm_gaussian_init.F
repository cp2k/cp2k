!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Initialize the use of the gaussians to treat the QMMM
!>      coupling potential
!> \par History
!>      6.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
MODULE qmmm_gaussian_init
  USE ao_util,                         ONLY: exp_radius
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                             gridlevel_info_type
  USE input_constants,                 ONLY: do_qmmm_gauss,&
                                             do_qmmm_swave
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_p_type
  USE qmmm_gaussian_data,              ONLY: max_geep_lib_gauss,&
                                             min_geep_lib_gauss
  USE qmmm_gaussian_input,             ONLY: read_mm_potential,&
                                             set_mm_potential_erf,&
                                             set_mm_potential_swave
  USE qmmm_gaussian_types,             ONLY: qmmm_gaussian_p_type,&
                                             qmmm_gaussian_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_gaussian_init'
  PUBLIC :: qmmm_gaussian_initialize

CONTAINS

! *****************************************************************************
!> \brief Initialize the Gaussian QMMM Environment
!> \param qmmm_gaussian_fns ...
!> \param para_env ...
!> \param pw_env ...
!> \param mm_el_pot_radius ...
!> \param mm_el_pot_radius_corr ...
!> \param qmmm_coupl_type ...
!> \param eps_mm_rspace ...
!> \param maxradius ...
!> \param maxchrg ...
!> \param compatibility ...
!> \param print_section ...
!> \param qmmm_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      06.2004 created [tlaino]
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE qmmm_gaussian_initialize(qmmm_gaussian_fns, para_env, pw_env,&
       mm_el_pot_radius, mm_el_pot_radius_corr,  &
       qmmm_coupl_type, eps_mm_rspace, maxradius, maxchrg, compatibility,&
       print_section, qmmm_section, error)
    TYPE(qmmm_gaussian_p_type), &
      DIMENSION(:), POINTER                  :: qmmm_gaussian_fns
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(pw_env_type), POINTER               :: pw_env
    REAL(KIND=dp), DIMENSION(:), POINTER     :: mm_el_pot_radius, &
                                                mm_el_pot_radius_corr
    INTEGER, INTENT(IN)                      :: qmmm_coupl_type
    REAL(KIND=dp), INTENT(IN)                :: eps_mm_rspace
    REAL(KIND=dp), DIMENSION(:), POINTER     :: maxradius
    REAL(KIND=dp), INTENT(IN)                :: maxchrg
    LOGICAL, INTENT(IN)                      :: compatibility
    TYPE(section_vals_type), POINTER         :: print_section, qmmm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_gaussian_initialize', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ilevel, j, Ndim, &
                                                num_geep_gauss, output_unit, &
                                                stat
    LOGICAL                                  :: Failure, Found, use_geep_lib
    REAL(KIND=dp)                            :: alpha, mymaxradius, Prefactor
    REAL(KIND=dp), DIMENSION(:), POINTER     :: c_radius, radius
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(gridlevel_info_type), POINTER       :: gridlevel_info
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(qmmm_gaussian_type), POINTER        :: mypgf

! Statements

    NULLIFY(mypgf, gridlevel_info, radius, c_radius, logger)
    failure = .FALSE.
    IF (.NOT.failure) THEN
       logger => cp_error_get_logger(error)
       CALL section_vals_val_get(qmmm_section,"USE_GEEP_LIB",i_val=num_geep_gauss,error=error)
       IF (num_geep_gauss == 0) THEN
          use_geep_lib = .FALSE.
       ELSE
          use_geep_lib = .TRUE.
          CPPostcondition(num_geep_gauss>=min_geep_lib_gauss,cp_failure_level,routineP,error,failure)
          CPPostcondition(num_geep_gauss<=max_geep_lib_gauss,cp_failure_level,routineP,error,failure)
       END IF
       SELECT CASE(qmmm_coupl_type)
       CASE(do_qmmm_gauss,do_qmmm_swave)
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
             IF     (qmmm_coupl_type==do_qmmm_gauss) THEN
                CALL  set_mm_potential_erf(qmmm_gaussian_fns,&
                     compatibility, num_geep_gauss, error)
             ELSEIF (qmmm_coupl_type==do_qmmm_swave) THEN
                CALL  set_mm_potential_swave(qmmm_gaussian_fns,&
                     num_geep_gauss, error)
             END IF
          ELSE
             CALL read_mm_potential(para_env, qmmm_gaussian_fns,&
                  (compatibility.AND.(qmmm_coupl_type==do_qmmm_gauss)), qmmm_section, error)
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
          output_unit =cp_print_key_unit_nr(logger,print_section,"PROGRAM_RUN_INFO",&
          extension=".qmmmLog",error=error)
          IF (output_unit>0) WRITE(output_unit,'(A)')" QMMM Gaussian Data Not Initialized!"
          CALL cp_print_key_finished_output(output_unit,logger,print_section,"PROGRAM_RUN_INFO",error=error)
       END SELECT
    END IF
  END SUBROUTINE qmmm_gaussian_initialize

END MODULE qmmm_gaussian_init
