!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief defines types for metadynamics calculation
!> \par History
!>      01.2005 created [fawzi and ale]
! *****************************************************************************
MODULE metadynamics_types
  USE cp_para_env,                     ONLY: cp_para_env_release,&
                                             cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE parallel_rng_types,              ONLY: delete_rng_stream,&
                                             rng_stream_p_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  INTEGER, SAVE, PRIVATE               :: last_meta_env_id=0
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'metadynamics_types'

  PUBLIC :: metadyn_create,&
            hills_env_type,&
            meta_env_type,&
            meta_env_retain,&
            set_meta_env,&
            meta_env_release,&
            metavar_type,&
            multiple_walkers_type

! *****************************************************************************
!> \brief defines types for HILLS
! *****************************************************************************
  TYPE hills_env_type
     LOGICAL                                   :: restart, slow_growth, wtcontrol
     !RG Adaptive hills
     REAL(KIND=dp)                             :: ww, min_disp, energy
     INTEGER                                   :: n_hills, nt_hills, min_nt_hills
     INTEGER                                   :: old_hill_number, old_hill_step
     !RG Adaptive hills
     REAL(KIND=dp), DIMENSION(:,:), POINTER    :: ss_history
     REAL(KIND=dp), DIMENSION(:,:), POINTER    :: delta_s_history
     REAL(KIND=dp), DIMENSION(:),   POINTER    :: ww_history
     REAL(KIND=dp), DIMENSION(:),   POINTER    :: invdt_history
     !Hills tail damping
     REAL(KIND=dp)                             :: tail_cutoff
     INTEGER                                   :: p_exp, q_exp
  END TYPE hills_env_type

! *****************************************************************************
!> \brief defines types for WALLS
! *****************************************************************************
  TYPE wall_type
     INTEGER                                   :: id_type, id_direction
     REAL(KIND=dp)                             :: pos, pos0
     REAL(KIND=dp)                             :: k_quadratic, k_quartic, ww_gauss, sigma_gauss
  END TYPE wall_type

! *****************************************************************************
!> \brief defines types for COLVAR used in the metadynamics
! *****************************************************************************
  TYPE metavar_type
     INTEGER                                   :: icolvar
     LOGICAL                                   :: do_wall, periodic
     REAL(KIND=dp)                             :: mass, lambda, vvp
     REAL(KIND=dp)                             :: gamma
     REAL(KIND=dp)                             :: epot_s, delta_s, epot_walls
     REAL(KIND=dp)                             :: ss, ss0, ff_s, ff_hills, ff_walls
     TYPE(wall_type), DIMENSION(:), POINTER    :: walls
  END TYPE metavar_type

! *****************************************************************************
!> \brief defines types for multiple walkers run
! *****************************************************************************
  TYPE multiple_walkers_type
     INTEGER                                   :: n_hills_local
     INTEGER                                   :: walker_id
     INTEGER                                   :: walkers_tot_nr
     INTEGER                                   :: walkers_freq_comm
     INTEGER, DIMENSION(:), POINTER            :: walkers_status
     CHARACTER(LEN=default_path_length),&
          DIMENSION(:), POINTER                :: walkers_file_name
  END TYPE multiple_walkers_type

! *****************************************************************************
!> \brief defines meta_env type
! *****************************************************************************
  TYPE meta_env_type
     INTEGER                                   :: ref_count, id_nr
     LOGICAL                                   :: do_hills, do_multiple_walkers
     LOGICAL                                   :: extended_lagrange
     LOGICAL                                   :: well_tempered
     LOGICAL                                   :: langevin
     LOGICAL                                   :: use_plumed
     CHARACTER(LEN=default_path_length)        :: plumed_input_file
     INTEGER                                   :: n_colvar
     REAL(KIND=dp)                             :: ekin_s, epot_s, dt, &
                                                  avg_temp, epot_walls
     LOGICAL                                   :: tempcontrol, restart
     REAL(KIND=dp)                             :: temp_wanted,toll_temp
     REAL(KIND=dp)                             :: delta_t, invdt, &
                                                  wtgamma, wttemperature
     INTEGER                                   :: n_steps
     ! time pointer should never be allocated itself.. that's a quite bad
     ! choice derived from the md_env.. So whenever the meta_env type is
     ! used the pointer time should be pointer to the high level time counter
     REAL(KIND=dp), POINTER                    :: time
     TYPE(hills_env_type), POINTER             :: hills_env
     TYPE(metavar_type), POINTER, DIMENSION(:) :: metavar
     TYPE(multiple_walkers_type), POINTER      :: multiple_walkers
     TYPE(cp_para_env_type), POINTER           :: para_env
     TYPE(section_vals_type), POINTER          :: metadyn_section
     TYPE(rng_stream_p_type), DIMENSION(:),&
                                       POINTER :: rng
     INTEGER                                   :: TAMCSteps
     REAL(kind=dp)                             :: zdt
  END TYPE meta_env_type

CONTAINS

! *****************************************************************************
!> \brief allocates a metadynamic environment (performs only minimal
!>      initialization)
!> \param meta_env the meta env_ that will be allocated
!> \param n_colvar number of collectiva variables
!> \param dt ...
!> \param para_env ...
!> \param metadyn_section ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \param do_hills if the hills should be added UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \par History
!>      04.2004 created
!>      02.2006 Reorganized the structure of the restart for Metadynamics (teo)
!>              cleaned the metadynamic type
!> \author - alessandro laio and fawzi mohamed
!>         - Teodoro Laino [tlaino] - University of Zurich. 11.2007
!>         - Teodoro Laino [tlaino] - University of Zurich. 10.2008
!>           Major rewriting and addition of multiple walkers
! *****************************************************************************
  SUBROUTINE metadyn_create(meta_env, n_colvar, dt, para_env, metadyn_section, error)
    TYPE(meta_env_type), POINTER             :: meta_env
    INTEGER, INTENT(in)                      :: n_colvar
    REAL(dp), INTENT(in)                     :: dt
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: metadyn_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: do_langevin, failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(meta_env),cp_failure_level,routinep,error,failure)
    ALLOCATE(meta_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(meta_env%multiple_walkers,&
               meta_env%metadyn_section,&
               meta_env%time,&
               meta_env%rng,&
               meta_env%hills_env)

       meta_env%use_plumed = .FALSE.        
       meta_env%plumed_input_file = ""
       meta_env%metadyn_section => metadyn_section
       meta_env%ref_count = 1
       meta_env%restart   = .TRUE.
       last_meta_env_id   =  last_meta_env_id+1
       meta_env%id_nr     =  last_meta_env_id
       meta_env%n_colvar  =  n_colvar
       meta_env%para_env  => para_env
       CALL cp_para_env_retain(para_env,error=error)

       meta_env%ekin_s      = 0.0_dp
       meta_env%epot_s      = 0.0_dp
       meta_env%epot_walls  = 0.0_dp
       meta_env%n_steps     = 0
       meta_env%dt          = dt
       meta_env%tempcontrol = .FALSE.

       ! Hills_env
       ALLOCATE(meta_env%hills_env, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(meta_env%hills_env%ss_history(n_colvar,0), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(meta_env%hills_env%delta_s_history(n_colvar,0), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(meta_env%hills_env%ww_history(0), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(meta_env%hills_env%invdt_history(0), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       meta_env%hills_env%n_hills = 0
       meta_env%hills_env%energy = 0.0_dp
       meta_env%hills_env%restart = .TRUE.

       ! Colvar
       ALLOCATE(meta_env%metavar(n_colvar), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DO i = 1, n_colvar
          NULLIFY(meta_env%metavar(i)%walls)
          meta_env%metavar(i)%mass       = -HUGE(0.0_dp)
          meta_env%metavar(i)%lambda     = -HUGE(0.0_dp)
          meta_env%metavar(i)%gamma      = 0.0_dp
          meta_env%metavar(i)%ss         = 0.0_dp
          meta_env%metavar(i)%ss0        = 0.0_dp
          meta_env%metavar(i)%ff_s       = 0.0_dp
          meta_env%metavar(i)%vvp        = 0.0_dp
          meta_env%metavar(i)%epot_s     = 0.0_dp
          meta_env%metavar(i)%epot_walls = 0.0_dp
          meta_env%metavar(i)%delta_s    = 0.0_dp
          meta_env%metavar(i)%ff_hills   = 0.0_dp
          meta_env%metavar(i)%ff_walls   = 0.0_dp
          meta_env%metavar(i)%do_wall    = .FALSE.
          meta_env%metavar(i)%periodic   = .FALSE.
          meta_env%metavar(i)%icolvar    = 0
       END DO

       ! Multiple Walkers
       CALL section_vals_val_get(metadyn_section,"MULTIPLE_WALKERS%_SECTION_PARAMETERS_",&
            l_val=meta_env%do_multiple_walkers,error=error)
       IF (meta_env%do_multiple_walkers) THEN
          ALLOCATE(meta_env%multiple_walkers,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)

          ! Walkers status and Walkers file name
          NULLIFY(meta_env%multiple_walkers%walkers_status,&
                  meta_env%multiple_walkers%walkers_file_name)
          meta_env%multiple_walkers%n_hills_local = 0
       END IF

       CALL section_vals_val_get(metadyn_section,"LANGEVIN",l_val=do_langevin,error=error)
       IF (do_langevin) THEN
          ALLOCATE (meta_env%rng(meta_env%n_colvar),STAT=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          DO i=1,meta_env%n_colvar
             NULLIFY(meta_env%rng(meta_env%n_colvar)%stream)
          END DO
       ENDIF
    END IF
  END SUBROUTINE metadyn_create

! *****************************************************************************
!> \brief sets the meta_env
!> \param meta_env ...
!> \param time ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author alessandro laio and fawzi mohamed
! *****************************************************************************
  SUBROUTINE set_meta_env(meta_env, time, error)
    TYPE(meta_env_type), POINTER             :: meta_env
    REAL(KIND=dp), OPTIONAL, POINTER         :: time
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'set_meta_env', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(meta_env)) THEN
       IF (PRESENT(time)) THEN
          NULLIFY(meta_env%time)
          meta_env%time => time
       END IF
    END IF
  END SUBROUTINE set_meta_env

! *****************************************************************************
!> \brief retains the meta_env
!> \param meta_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author alessandro laio and fawzi mohamed
! *****************************************************************************
  SUBROUTINE meta_env_retain(meta_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPreconditionNoFail(ASSOCIATED(meta_env),cp_failure_level,routineP,error)
    CPPreconditionNoFail(meta_env%ref_count>0,cp_failure_level,routineP,error)
    meta_env%ref_count=meta_env%ref_count+1
  END SUBROUTINE meta_env_retain

! *****************************************************************************
!> \brief releases the meta_env
!> \param meta_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author alessandro laio and fawzi mohamed
! *****************************************************************************
  SUBROUTINE meta_env_release(meta_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (ASSOCIATED(meta_env)) THEN
       CPPreconditionNoFail(meta_env%ref_count>0,cp_failure_level,routineP,error)
       meta_env%ref_count=meta_env%ref_count-1
       IF (meta_env%ref_count==0) THEN
          CALL cp_para_env_release(meta_env%para_env,error=error)
          IF (ASSOCIATED(meta_env%metavar)) THEN
             DO i = 1, SIZE(meta_env%metavar)
                IF (ASSOCIATED(meta_env%metavar(i)%walls)) THEN
                   DEALLOCATE(meta_env%metavar(i)%walls,stat=stat)
                   CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
                END IF
             END DO
             DEALLOCATE(meta_env%metavar,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF

          ! Hills env
          CALL hills_env_release(meta_env%hills_env,error)
          ! Walkers type
          IF (ASSOCIATED(meta_env%multiple_walkers)) THEN
             IF (ASSOCIATED(meta_env%multiple_walkers%walkers_status)) THEN
                DEALLOCATE(meta_env%multiple_walkers%walkers_status,stat=stat)
                CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
             END IF
             IF (ASSOCIATED(meta_env%multiple_walkers%walkers_file_name)) THEN
                DEALLOCATE(meta_env%multiple_walkers%walkers_file_name,stat=stat)
                CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
             END IF
             DEALLOCATE(meta_env%multiple_walkers,stat=stat)
             CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          END IF

          ! Langevin on COLVARS
          IF (meta_env%langevin) THEN
             DO i = 1, SIZE(meta_env%rng)
                IF (ASSOCIATED(meta_env%rng(i)%stream)) THEN
                   CALL delete_rng_stream(meta_env%rng(i)%stream,error=error)
                END IF
             END DO
             DEALLOCATE (meta_env%rng,STAT=stat)
             CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ENDIF
          NULLIFY(meta_env%time)
          NULLIFY(meta_env%metadyn_section)
          DEALLOCATE(meta_env, stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
  END SUBROUTINE meta_env_release

! *****************************************************************************
!> \brief releases the hills_env
!> \param hills_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - University of Zurich 10.2008
! *****************************************************************************
  SUBROUTINE hills_env_release(hills_env, error)
    TYPE(hills_env_type), POINTER            :: hills_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hills_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure = .FALSE.
    IF (ASSOCIATED(hills_env)) THEN
       IF (ASSOCIATED(hills_env%ss_history)) THEN
          DEALLOCATE(hills_env%ss_history,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(hills_env%delta_s_history)) THEN
          DEALLOCATE(hills_env%delta_s_history,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(hills_env%ww_history)) THEN
          DEALLOCATE(hills_env%ww_history,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED(hills_env%invdt_history)) THEN
          DEALLOCATE(hills_env%invdt_history,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
       DEALLOCATE(hills_env,stat=stat)
       CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
    END IF
  END SUBROUTINE hills_env_release

END MODULE metadynamics_types
