!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Performs the metadynamics calculation
!> \par History
!>      01.2005 created [fawzi and ale]
!>      11.2007 Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
MODULE metadynamics
  USE colvar_methods,                  ONLY: colvar_eval_glob_f
  USE colvar_types,                    ONLY: colvar_p_type,&
                                             torsion_colvar_id
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_para_env,                     ONLY: cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_constants,                 ONLY: do_fe_meta,&
                                             do_wall_gaussian,&
                                             do_wall_m,&
                                             do_wall_none,&
                                             do_wall_p,&
                                             do_wall_quadratic,&
                                             do_wall_reflective
  USE input_cp2k_free_energy,          ONLY: create_metavar_section
  USE input_enumeration_types,         ONLY: enum_i2c,&
                                             enumeration_type
  USE input_keyword_types,             ONLY: keyword_get,&
                                             keyword_type
  USE input_section_types,             ONLY: section_get_keyword,&
                                             section_get_subsection,&
                                             section_release,&
                                             section_type,&
                                             section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE machine,                         ONLY: m_flush,&
                                             m_mov
  USE metadynamics_types,              ONLY: hills_env_type,&
                                             meta_env_type,&
                                             metavar_type,&
                                             multiple_walkers_type
  USE parallel_rng_types,              ONLY: next_random_number
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: femtoseconds
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE 

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'metadynamics'

  INTEGER, SAVE, PRIVATE :: last_meta_env_id=0

  PUBLIC :: metadyn_create, metadyn_read, metadyn_forces

CONTAINS
  
! *****************************************************************************
!> \brief allocates a metadynamic environment (performs only minimal 
!>      initialization)
!> \param meta_env the meta env_ that will be allocated
!> \param n_colvar number of collectiva variables
!> \param do_hills if the hills should be added
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
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
    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(meta_env),cp_failure_level,routinep,error,failure)
    ALLOCATE(meta_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(meta_env%multiple_walkers,&
               meta_env%metadyn_section,&
               meta_env%time,&
               meta_env%hills_env)

       meta_env%metadyn_section => metadyn_section
       meta_env%ref_count = 1
       meta_env%restart   = .TRUE.
       last_meta_env_id   =  last_meta_env_id+1
       meta_env%id_nr     =  last_meta_env_id
       meta_env%n_colvar  =  n_colvar
       meta_env%para_env  => para_env
       CALL cp_para_env_retain(para_env,error=error)

       meta_env%ekin_s      = 0.0_dp
       meta_env%cv_energy   = 0.0_dp
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
       meta_env%hills_env%n_hills = 0
       meta_env%hills_env%restart = .TRUE.

       ! Colvar
       ALLOCATE(meta_env%metavar(n_colvar), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DO i = 1, n_colvar
          NULLIFY(meta_env%metavar(i)%walls)
          meta_env%metavar(i)%mass      = -HUGE(0.0_dp)
          meta_env%metavar(i)%lambda    = -HUGE(0.0_dp)
          meta_env%metavar(i)%ss        = 0.0_dp
          meta_env%metavar(i)%ss0       = 0.0_dp
          meta_env%metavar(i)%ff_s      = 0.0_dp
          meta_env%metavar(i)%vvp       = 0.0_dp
          meta_env%metavar(i)%cv_energy = 0.0_dp
          meta_env%metavar(i)%delta_s   = 0.0_dp
          meta_env%metavar(i)%ff_hills  = 0.0_dp
          meta_env%metavar(i)%do_wall   = .FALSE.
          meta_env%metavar(i)%periodic  = .FALSE.
          meta_env%metavar(i)%icolvar   = 0
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
    END IF
  END SUBROUTINE metadyn_create

! *****************************************************************************
!> \brief reads metadynamics section
!> \par History
!>      04.2004 created
!> \author Teodoro Laino [tlaino] - University of Zurich. 11.2007 
! *****************************************************************************
  SUBROUTINE metadyn_read(meta_env,root_section,para_env,fe_section,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: fe_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_read', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: walkers_file_name
    INTEGER                                  :: i, id_method, n_colvar, &
                                                n_rep, stat
    INTEGER, DIMENSION(:), POINTER           :: walkers_status
    LOGICAL                                  :: check, explicit, failure
    REAL(kind=dp)                            :: dt
    TYPE(section_vals_type), POINTER         :: md_section, metadyn_section, &
                                                metavar_section, &
                                                walkers_section

    failure = .FALSE.
    CALL section_vals_get(fe_section,explicit=explicit,error=error)
    IF (explicit) THEN
       CALL section_vals_val_get(fe_section,"METHOD",i_val=id_method,error=error)
       IF (id_method/=do_fe_meta) RETURN
       metadyn_section => section_vals_get_subs_vals(fe_section,"METADYN",error=error)
       CPPreconditionNoFail(.NOT.ASSOCIATED(meta_env),cp_failure_level,routineP,error)  

       md_section => section_vals_get_subs_vals(root_section,"MOTION%MD",error=error)
       CALL section_vals_val_get(md_section,"TIMESTEP",r_val=dt,error=error)

       metavar_section => section_vals_get_subs_vals(metadyn_section,"METAVAR",error=error)
       CALL section_vals_get(metavar_section,n_repetition=n_colvar,error=error)
       CALL metadyn_create(meta_env,n_colvar=n_colvar,&
            dt=dt,para_env=para_env,metadyn_section=metadyn_section,error=error)

       CALL section_vals_val_get(metadyn_section,"DO_HILLS",l_val=meta_env%do_hills,error=error)
       CALL section_vals_val_get(metadyn_section,"LAGRANGE",l_val=meta_env%extended_lagrange,error=error)
       CALL section_vals_val_get(metadyn_section,"WW",r_val=meta_env%hills_env%ww,error=error)
       CALL section_vals_val_get(metadyn_section,"NT_HILLS",i_val=meta_env%hills_env%nt_hills,error=error)
       CALL section_vals_val_get(metadyn_section,"STEP_START_VAL",i_val=meta_env%n_steps,error=error)
       CPPreconditionNoFail(meta_env%n_steps>=0,cp_failure_level,routineP,error)
       CALL section_vals_val_get(metadyn_section,"NHILLS_START_VAL",&
            i_val=meta_env%hills_env%n_hills,error=error)
       CALL section_vals_val_get(metadyn_section,"TEMPERATURE",r_val=meta_env%temp_wanted,error=error)
       CALL section_vals_val_get(metadyn_section,"TEMP_TOL",explicit=meta_env%tempcontrol,&
            r_val=meta_env%toll_temp,error=error)

       ! Parsing Metavar Section
       DO i= 1, n_colvar
          CALL metavar_read(meta_env%metavar(i), meta_env%extended_lagrange,i, metavar_section,error=error)
       END DO

       ! Parsing the Multiple Walkers Info
       IF (meta_env%do_multiple_walkers) THEN
          NULLIFY(walkers_status)
          walkers_section => section_vals_get_subs_vals(metadyn_section,"MULTIPLE_WALKERS",error=error)

          ! General setup for walkers
          CALL section_vals_val_get(walkers_section,"WALKER_ID",&
               i_val=meta_env%multiple_walkers%walker_id,error=error)
          CALL section_vals_val_get(walkers_section,"NUMBER_OF_WALKERS",&
               i_val=meta_env%multiple_walkers%walkers_tot_nr,error=error)
          CALL section_vals_val_get(walkers_section,"WALKER_COMM_FREQUENCY",&
               i_val=meta_env%multiple_walkers%walkers_freq_comm,error=error)
          
          ! Handle status and file names
          ALLOCATE(meta_env%multiple_walkers%walkers_status(meta_env%multiple_walkers%walkers_tot_nr),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          ALLOCATE(meta_env%multiple_walkers%walkers_file_name(meta_env%multiple_walkers%walkers_tot_nr),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          CALL section_vals_val_get(walkers_section,"WALKERS_STATUS",explicit=explicit,error=error)
          IF (explicit) THEN
             CALL section_vals_val_get(walkers_section,"WALKERS_STATUS",i_vals=walkers_status,&
                  error=error)
             check = (SIZE(walkers_status)==meta_env%multiple_walkers%walkers_tot_nr)
             CALL cp_assert(check,cp_failure_level,cp_assertion_failed,routineP,&
                  "Number of Walkers specified in the input does not match with the "//&
                  "size of the WALKERS_STATUS. Please check your input and in case "//&
                  "this is a restart run consider the possibility to switch off the "//&
                  "RESTART_WALKERS in the EXT_RESTART section! "//&
CPSourceFileRef,&
                  only_ionode=.TRUE.)
             meta_env%multiple_walkers%walkers_status = walkers_status
          ELSE
             meta_env%multiple_walkers%walkers_status = 0
          END IF
          meta_env%multiple_walkers%n_hills_local = &
               meta_env%multiple_walkers%walkers_status(meta_env%multiple_walkers%walker_id)

          CALL section_vals_val_get(walkers_section,"WALKERS_FILE_NAME%_DEFAULT_KEYWORD_",&
               n_rep_val=n_rep, error=error)
          check = (n_rep==meta_env%multiple_walkers%walkers_tot_nr)
          CALL cp_assert(check,cp_failure_level,cp_assertion_failed,routineP,&
               "Number of Walkers specified in the input does not match with the "//&
               "number of Walkes File names provided. Please check your input and in case "//&
               "this is a restart run consider the possibility to switch off the "//&
               "RESTART_WALKERS in the EXT_RESTART section! "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          DO i = 1, n_rep
             CALL section_vals_val_get(walkers_section,"WALKERS_FILE_NAME%_DEFAULT_KEYWORD_",&
                  i_rep_val=i, c_val=walkers_file_name, error=error)
             meta_env%multiple_walkers%walkers_file_name(i) = walkers_file_name
          END DO
       END IF

       ! Print Metadynamics Info
       CALL print_metadyn_info(meta_env, n_colvar, metadyn_section, error)
    END IF

  END SUBROUTINE metadyn_read

! *****************************************************************************
!> \brief prints information on the metadynamics run
!> \author Teodoro Laino [tlaino] - University of Zurich. 10.2008
! *****************************************************************************
  SUBROUTINE print_metadyn_info(meta_env, n_colvar, metadyn_section, error)
    TYPE(meta_env_type), POINTER             :: meta_env
    INTEGER, INTENT(IN)                      :: n_colvar
    TYPE(section_vals_type), POINTER         :: metadyn_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'print_metadyn_info', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=10)                        :: my_id, my_tag
    INTEGER                                  :: i, iw, j
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section, wall_section, &
                                                work_section

    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,metadyn_section,&
         "PRINT%PROGRAM_RUN_INFO",extension=".metadynLog",error=error)
    NULLIFY(section,enum,keyword)
    CALL create_metavar_section(section,error=error)
    wall_section => section_get_subsection(section,"WALL",error)
    IF (iw>0) THEN
       WRITE (iw,'( /A )')          ' METADYN| Meta Dynamics Protocol '
       WRITE (iw,'( A,T71,I10)')    ' METADYN| Number of interval time steps to spawn hills',&
            meta_env%hills_env%nt_hills
       WRITE (iw,'( A,T71,I10)')    ' METADYN| Number of previously spawned hills',&
            meta_env%hills_env%n_hills
       IF (meta_env%extended_lagrange) THEN
          WRITE (iw,'( A )')        ' METADYN| Extended Lagrangian Scheme '
          IF (meta_env%tempcontrol) WRITE (iw,'( A,T71,F10.2)')&
               ' METADYN| Collective Variables Temperature control',meta_env%toll_temp
       END IF

       IF (meta_env%do_hills) THEN 
          WRITE (iw,'( A )')' METADYN| Spanwing the Hills '
          WRITE (iw,'( A,T71,F10.2)')  ' METADYN| Width of the Spawned Gaussian',meta_env%hills_env%ww
       END IF

       IF (meta_env%do_multiple_walkers) THEN
          WRITE (iw,'( A,T71,A10)')    ' METADYN| Multiple Walkers','   ENABLED'
          WRITE (iw,'( A,T71,I10)')    ' METADYN| Number of Multiple Walkers',&
               meta_env%multiple_walkers%walkers_tot_nr
          WRITE (iw,'( A,T71,I10)')    ' METADYN| Local Walker ID',&
               meta_env%multiple_walkers%walker_id
          WRITE (iw,'( A,T71,I10)')    ' METADYN| Walker Communication Frequency',&
               meta_env%multiple_walkers%walkers_freq_comm
          DO i = 1, meta_env%multiple_walkers%walkers_tot_nr
             my_tag = ""
             IF (i==meta_env%multiple_walkers%walker_id) my_tag = " ( Local )"
             my_id = '( '//TRIM(ADJUSTL(cp_to_string(i)))//' )'
             WRITE (iw,'(/,A,T71,A10)')    ' WALKERS| Walker ID'//TRIM(my_tag),ADJUSTR(my_id)
             WRITE (iw,'(  A,T71,I10)')    ' WALKERS| Number of Hills communicated',&
                  meta_env%multiple_walkers%walkers_status(i)
             WRITE (iw,'(  A,T24,A57)')    ' WALKERS| Base Filename',&
                  ADJUSTR(meta_env%multiple_walkers%walkers_file_name(i)(1:57))
          END DO
          WRITE (iw,'(/)')
       END IF

       WRITE (iw,'( A,T71,I10)')    ' METADYN| Number of collective variables',meta_env%n_colvar
       DO i=1,n_colvar
          WRITE (iw,'( A )')'          '//'----------------------------------------------------------------------'
          WRITE (iw,'( A,T71,I10)')      ' METAVARS| Collective Variable Number',meta_env%metavar(i)%icolvar
          IF (meta_env%extended_lagrange) THEN
             WRITE (iw,'( A,T71,F10.6)')    ' METAVARS| Lambda Parameter',meta_env%metavar(i)%lambda
             WRITE (iw,'( A,T66,F15.6)')    ' METAVARS| Collective Variable Mass',meta_env%metavar(i)%mass
          END IF
          WRITE (iw,'( A,T71,F10.6)')    ' METAVARS| Scaling factor',meta_env%metavar(i)%delta_s
          IF (meta_env%metavar(i)%do_wall) THEN
             WRITE (iw,'( A,T71,I10)')    ' METAVARS| Number of Walls present',SIZE(meta_env%metavar(i)%walls)
             DO j = 1, SIZE(meta_env%metavar(i)%walls)
                keyword => section_get_keyword(wall_section,"TYPE",error=error)
                CALL keyword_get(keyword,enum=enum,error=error)
                WRITE (iw,'(/,A,5X,I10,T50,A,T70,A11)')    ' METAVARS| Wall Number:',j,'Type of Wall:',&
                     ADJUSTR(TRIM(enum_i2c(enum,meta_env%metavar(i)%walls(j)%id_type,error=error)))
                ! Type of wall IO
                SELECT CASE(meta_env%metavar(i)%walls(j)%id_type)
                CASE(do_wall_none)
                   ! Do Nothing
                   CYCLE
                CASE(do_wall_reflective)
                   work_section => section_get_subsection(wall_section,"REFLECTIVE",error)
                   keyword => section_get_keyword(work_section,"DIRECTION",error=error)
                   CALL keyword_get(keyword,enum=enum,error=error)
                   WRITE (iw,'(A,T70,A11)')    ' METAVARS| Wall direction',&
                        ADJUSTR(TRIM(enum_i2c(enum,meta_env%metavar(i)%walls(j)%id_direction,error=error)))
                CASE(do_wall_quadratic)
                   work_section => section_get_subsection(wall_section,"QUADRATIC",error)
                   keyword => section_get_keyword(work_section,"DIRECTION",error=error)
                   CALL keyword_get(keyword,enum=enum,error=error)
                   WRITE (iw,'(A,T70,A11)')    ' METAVARS| Wall direction',&
                        ADJUSTR(TRIM(enum_i2c(enum,meta_env%metavar(i)%walls(j)%id_direction,error=error)))
                   WRITE (iw,'(A,T70,F11.6)')    ' METAVARS| Constant K of the quadratic potential',&
                        meta_env%metavar(i)%walls(j)%k_quadratic
                CASE(do_wall_gaussian)
                   WRITE (iw,'(A,T70,F11.6)')    ' METAVARS| Height of the Wall Gaussian',&
                        meta_env%metavar(i)%walls(j)%ww_gauss
                   WRITE (iw,'(A,T70,F11.6)')    ' METAVARS| Amplitude of the Wall Gaussian',&
                        meta_env%metavar(i)%walls(j)%sigma_gauss 
                END SELECT
                WRITE (iw,'(A,T70,F11.6)')    ' METAVARS| Wall location',&
                     meta_env%metavar(i)%walls(j)%pos
             END DO
          END IF
          WRITE (iw,'( A )')'          '//'----------------------------------------------------------------------'
       ENDDO
    END IF
    CALL section_release(section,error=error)
    CALL cp_print_key_finished_output(iw,logger,metadyn_section,"PRINT%PROGRAM_RUN_INFO", error=error)

  END SUBROUTINE print_metadyn_info

! *****************************************************************************
!> \brief reads metavar section
!> \par History
!>      04.2004 created
!> \author alessandro laio and fawzi mohamed
!>      Teodoro Laino [tlaino] - University of Zurich. 11.2007
! *****************************************************************************
  SUBROUTINE metavar_read(metavar,extended_lagrange,icol,metavar_section,error)
    TYPE(metavar_type), INTENT(INOUT)        :: metavar
    LOGICAL                                  :: extended_lagrange
    INTEGER, INTENT(IN)                      :: icol
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: metavar_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metavar_read', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n_walls, stat
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: wall_section, work_section

    failure=.FALSE.
    CALL section_vals_val_get(metavar_section,"COLVAR",i_rep_section=icol,i_val=metavar%icolvar,error=error)
    CALL section_vals_val_get(metavar_section,"SCALE",i_rep_section=icol,r_val=metavar%delta_s,error=error)
    ! Walls
    wall_section => section_vals_get_subs_vals(metavar_section,"WALL",i_rep_section=icol,error=error)
    CALL section_vals_get(wall_section,n_repetition=n_walls,error=error)
    IF (n_walls/=0) THEN
       metavar%do_wall = .TRUE.
       ALLOCATE(metavar%walls(n_walls),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DO i = 1, n_walls
          CALL section_vals_val_get(wall_section,"TYPE",i_rep_section=i,i_val=metavar%walls(i)%id_type,error=error)
          SELECT CASE(metavar%walls(i)%id_type)
          CASE(do_wall_none)
             ! Just cycle..
             CYCLE
          CASE(do_wall_reflective)
             work_section => section_vals_get_subs_vals(wall_section,"REFLECTIVE",i_rep_section=i,error=error)
             CALL section_vals_val_get(work_section,"DIRECTION",i_val=metavar%walls(i)%id_direction,error=error)
          CASE(do_wall_quadratic)
             work_section => section_vals_get_subs_vals(wall_section,"QUADRATIC",i_rep_section=i,error=error)
             CALL section_vals_val_get(work_section,"DIRECTION",i_val=metavar%walls(i)%id_direction,error=error)
             CALL section_vals_val_get(work_section,"K",r_val=metavar%walls(i)%k_quadratic,error=error)
          CASE(do_wall_gaussian)
             work_section => section_vals_get_subs_vals(wall_section,"GAUSSIAN",i_rep_section=i,error=error)
             CALL section_vals_val_get(work_section,"WW",r_val=metavar%walls(i)%ww_gauss,error=error)
             CALL section_vals_val_get(work_section,"SIGMA",r_val=metavar%walls(i)%sigma_gauss,error=error)
          END SELECT
          CALL section_vals_val_get(wall_section,"POSITION",i_rep_section=i,r_val=metavar%walls(i)%pos,error=error)
       END DO
    END IF
    ! Setup few more parameters for extended lagrangian
    IF(extended_lagrange)THEN
       CALL section_vals_val_get(metavar_section,"MASS",i_rep_section=icol,r_val=metavar%mass,error=error)
       CALL section_vals_val_get(metavar_section,"LAMBDA",i_rep_section=icol,r_val=metavar%lambda,error=error)
    ENDIF

  END SUBROUTINE metavar_read

! *****************************************************************************
!> \brief add forces to the subsys due to the metadynamics run
!>      possibly modifies the velocites (if reflective walls are applied)
!> \par History
!>      04.2004 created
! *****************************************************************************
  SUBROUTINE metadyn_forces(force_env,vel,error)
    TYPE(force_env_type), POINTER            :: force_env
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: vel
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, i_c, icolvar, ii, &
                                                iw, iwall
    LOGICAL                                  :: explicit, failure
    REAL(kind=dp)                            :: check_val, diff_ss, dt, &
                                                ekin_w, fac_t, fft, norm, &
                                                rval, scal, scalf, ss0_test, &
                                                tol_ekin
    TYPE(colvar_p_type), DIMENSION(:), &
      POINTER                                :: colvar_p
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(metavar_type), POINTER              :: cv
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: ss0_section, vvp_section

    failure=.FALSE.
    NULLIFY(logger,meta_env)
    meta_env => force_env%meta_env
    IF (.NOT.ASSOCIATED(meta_env)) RETURN

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    NULLIFY(colvar_p,subsys,cv,ss0_section, vvp_section)
    CALL force_env_get(force_env, subsys=subsys, error=error)

    dt = meta_env%dt
    IF (.NOT.meta_env%restart) meta_env%n_steps=meta_env%n_steps+1

    ! Initialize velocity
    IF (meta_env%restart.AND.meta_env%extended_lagrange) THEN
       meta_env%ekin_s = 0.0_dp
       DO i_c=1,meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          cv%vvp = next_random_number(force_env%globenv%gaussian_rng_stream,error=error)
          meta_env%ekin_s = meta_env%ekin_s + 0.5E0_dp*cv%mass*cv%vvp**2
       END DO
       ekin_w = 0.5_dp*meta_env%temp_wanted*REAL(meta_env%n_colvar,KIND=dp)
       fac_t  = SQRT(ekin_w/MAX(meta_env%ekin_s,1.0E-8_dp))
       DO i_c = 1, meta_env%n_colvar
          cv     => meta_env%metavar(i_c)
          cv%vvp = cv%vvp*fac_t
       ENDDO
       meta_env%ekin_s = 0.0_dp
    END IF

    ! compute ss and the derivative of ss with respect to the atomic positions
    DO i_c=1,meta_env%n_colvar
       cv => meta_env%metavar(i_c)
       icolvar = cv%icolvar
       CALL colvar_eval_glob_f(icolvar,force_env,error=error)
       cv%ss       = subsys%colvar_p(icolvar)%colvar%ss

       ! Setup the periodic flag if the COLVAR is (-pi,pi] periodic
       cv%periodic = (subsys%colvar_p(icolvar)%colvar%type_id==torsion_colvar_id)

       ! Restart for Extended Lagrangian Metadynamics
       IF (meta_env%restart) THEN
          ! Initialize the position of the collective variable in the extended lagrange
          ss0_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_SS0",error=error)
          CALL section_vals_get(ss0_section, explicit=explicit, error=error)
          IF (explicit) THEN
             CALL section_vals_val_get(ss0_section,"_DEFAULT_KEYWORD_",&
                  i_rep_val=i_c, r_val=rval, error=error)
             cv%ss0 = rval
          ELSE
             cv%ss0 = cv%ss
          END IF
          vvp_section => section_vals_get_subs_vals(meta_env%metadyn_section,"EXT_LAGRANGE_VVP",error=error)
          CALL section_vals_get(vvp_section, explicit=explicit, error=error)
          IF (explicit) THEN
             CALL section_vals_val_get(vvp_section,"_DEFAULT_KEYWORD_",&
                  i_rep_val=i_c, r_val=rval, error=error)
             cv%vvp = rval
          END IF
       END IF
       !
       IF(.NOT.meta_env%extended_lagrange) THEN 
          cv%ss0 = cv%ss
          cv%vvp = 0.0_dp
       END IF
    ENDDO
    ! History dependent forces (evaluated at s0)
    IF(meta_env%do_hills) CALL  hills(meta_env,error)
    
    ! Apply walls to the colvars
    CALL walls(meta_env, error)

    meta_env%restart = .FALSE.
    IF(.NOT. meta_env%extended_lagrange)THEN
       meta_env%ekin_s    = 0.0_dp
       meta_env%cv_energy = 0.0_dp
       DO i_c = 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          cv%cv_energy = 0.0_dp
          cv%ff_s      = 0.0_dp
          icolvar=cv%icolvar
          NULLIFY(particles)
          CALL cp_subsys_get(subsys, colvar_p=colvar_p, &
               particles=particles,error=error)
          DO ii=1,colvar_p(icolvar)%colvar%n_atom_s
             i=colvar_p(icolvar)%colvar%i_atom(ii)
             particles%els(i)%f=particles%els(i)%f+cv%ff_hills*colvar_p(icolvar)%colvar%dsdr(:,ii)
          ENDDO
       ENDDO
    ELSE
       meta_env%ekin_s = 0.0_dp
       DO i_c= 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          diff_ss = cv%ss-cv%ss0
          IF (cv%periodic) THEN
             ! The difference of a periodic COLVAR is always within [-pi,pi]
             diff_ss = SIGN(1.0_dp,ASIN(SIN(diff_ss)))*ACOS(COS(diff_ss))
          END IF
          cv%cv_energy = 0.5_dp*cv%lambda*(diff_ss)**2.0_dp
          cv%ff_s      = cv%lambda*(diff_ss)
          icolvar      = cv%icolvar
          ! forces on the atoms
          NULLIFY(particles)
          CALL cp_subsys_get(subsys, colvar_p=colvar_p, &
               particles=particles,error=error)
          DO ii=1,colvar_p(icolvar)%colvar%n_atom_s
             i=colvar_p(icolvar)%colvar%i_atom(ii)
             particles%els(i)%f=particles%els(i)%f- cv%ff_s*colvar_p(icolvar)%colvar%dsdr(:,ii)
          ENDDO
          !  velocity verlet on the s0
          fft    = cv%ff_s+cv%ff_hills
          cv%vvp = cv%vvp+dt*fft/cv%mass
          meta_env%ekin_s = meta_env%ekin_s + 0.5E0_dp*cv%mass*cv%vvp**2
       ENDDO
       !  velocity rescaling on the s0
       IF(meta_env%tempcontrol)THEN
          ekin_w   = 0.5_dp*meta_env%temp_wanted*REAL(meta_env%n_colvar,KIND=dp)
          tol_ekin = 0.5_dp*meta_env%toll_temp*REAL(meta_env%n_colvar,KIND=dp)
          IF(ABS(ekin_w-meta_env%ekin_s)> tol_ekin)THEN
             fac_t=SQRT(ekin_w/MAX(meta_env%ekin_s,1.0E-8_dp))
             DO i_c = 1, meta_env%n_colvar
                cv     => meta_env%metavar(i_c)
                cv%vvp = cv%vvp*fac_t
             ENDDO
             meta_env%ekin_s = ekin_w
          ENDIF
       ENDIF
       ! Reflective Wall only for s0
       DO i_c = 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          IF (cv%do_wall) THEN
             DO iwall = 1, SIZE(cv%walls)
                SELECT CASE(cv%walls(iwall)%id_type)
                CASE(do_wall_reflective)
                   ss0_test=cv%ss0+dt*cv%vvp
                   IF (cv%periodic) THEN
                      ! A periodic COLVAR is always within [-pi,pi]
                      ss0_test = SIGN(1.0_dp,ASIN(SIN(ss0_test)))*ACOS(COS(ss0_test))
                   END IF
                   SELECT CASE(cv%walls(iwall)%id_direction)
                   CASE(do_wall_p)
                      IF ((ss0_test > cv%walls(iwall)%pos).AND.(cv%vvp>0)) cv%vvp = -cv%vvp
                   CASE(do_wall_m)
                      IF ((ss0_test < cv%walls(iwall)%pos).AND.(cv%vvp<0)) cv%vvp = -cv%vvp
                   END SELECT
                END SELECT
             END DO
          ENDIF
       ENDDO
       ! Update of ss0
       meta_env%cv_energy=0.0_dp
       DO i_c = 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          cv%ss0 = cv%ss0 + dt*cv%vvp
          IF (cv%periodic) THEN
             ! A periodic COLVAR is always within [-pi,pi]
             cv%ss0 = SIGN(1.0_dp,ASIN(SIN(cv%ss0)))*ACOS(COS(cv%ss0))
          END IF
          meta_env%cv_energy = meta_env%cv_energy + cv%cv_energy
       ENDDO
    ENDIF

    ! Reflective Wall only for ss
    DO i_c = 1, meta_env%n_colvar
       cv => meta_env%metavar(i_c)
       IF (cv%do_wall) THEN
          DO iwall = 1, SIZE(cv%walls)
             SELECT CASE(cv%walls(iwall)%id_type)
             CASE(do_wall_reflective)
                SELECT CASE(cv%walls(iwall)%id_direction)
                CASE(do_wall_p)
                   IF (cv%ss < cv%walls(iwall)%pos) CYCLE
                   check_val = -1.0_dp
                CASE(do_wall_m)
                   IF (cv%ss > cv%walls(iwall)%pos) CYCLE
                   check_val =  1.0_dp
                END SELECT
                NULLIFY(particles)
                icolvar=cv%icolvar
                CALL cp_subsys_get(subsys, colvar_p=colvar_p, particles=particles,error=error)
                scal=0.0_dp
                scalf=0.0_dp
                norm=0.0_dp
                DO ii=1,colvar_p(icolvar)%colvar%n_atom_s
                   i=colvar_p(icolvar)%colvar%i_atom(ii)
                   IF (PRESENT(vel)) THEN
                      scal=scal+vel(1,i)*colvar_p(icolvar)%colvar%dsdr(1,ii)
                      scal=scal+vel(2,i)*colvar_p(icolvar)%colvar%dsdr(2,ii)
                      scal=scal+vel(3,i)*colvar_p(icolvar)%colvar%dsdr(3,ii)                      
                   ELSE
                      scal=scal+particles%els(i)%v(1)*colvar_p(icolvar)%colvar%dsdr(1,ii)
                      scal=scal+particles%els(i)%v(2)*colvar_p(icolvar)%colvar%dsdr(2,ii)
                      scal=scal+particles%els(i)%v(3)*colvar_p(icolvar)%colvar%dsdr(3,ii)
                   END IF
                   scalf=scalf+particles%els(i)%f(1)*colvar_p(icolvar)%colvar%dsdr(1,ii)
                   scalf=scalf+particles%els(i)%f(2)*colvar_p(icolvar)%colvar%dsdr(2,ii)
                   scalf=scalf+particles%els(i)%f(3)*colvar_p(icolvar)%colvar%dsdr(3,ii)
                   norm=norm+colvar_p(icolvar)%colvar%dsdr(1,ii)**2
                   norm=norm+colvar_p(icolvar)%colvar%dsdr(2,ii)**2
                   norm=norm+colvar_p(icolvar)%colvar%dsdr(3,ii)**2
                ENDDO
                IF (norm /= 0.0_dp) scal = scal/norm
                IF (norm /= 0.0_dp) scalf= scalf/norm

                IF (scal*check_val>0.0_dp) CYCLE
                DO ii=1,colvar_p(icolvar)%colvar%n_atom_s
                   i=colvar_p(icolvar)%colvar%i_atom(ii)
                   IF (PRESENT(vel)) THEN
                      vel(:,i)=vel(:,i)-2.0_dp*colvar_p(icolvar)%colvar%dsdr(:,ii)*scal
                   ELSE
                      particles%els(i)%v(:)=particles%els(i)%v(:)-2.0_dp*colvar_p(icolvar)%colvar%dsdr(:,ii)*scal
                   END IF
                   ! Nullify forces along the colvar (this avoids the weird behaviors of the reflective wall)
                    particles%els(i)%f(:)=particles%els(i)%f(:)-colvar_p(icolvar)%colvar%dsdr(:,ii)*scalf
                ENDDO
             END SELECT
          END DO
       END IF
    ENDDO

    ! write COLVAR file
    iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
         "PRINT%COLVAR",extension=".metadynLog",error=error)
    IF (iw>0) THEN
       WRITE(iw,'(f12.1,40f13.5)')meta_env%time*femtoseconds, &
            (meta_env%metavar(i)%ss0,i=1,meta_env%n_colvar), &
            (meta_env%metavar(i)%ss,i=1,meta_env%n_colvar), &
            (meta_env%metavar(i)%ff_s,i=1,meta_env%n_colvar), &
            (meta_env%metavar(i)%ff_hills,i=1,meta_env%n_colvar)
       CALL m_flush(iw)
    END IF
    CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
         "PRINT%COLVAR", error=error)
    CALL timestop(handle)
  END SUBROUTINE metadyn_forces

! *****************************************************************************
!> \brief Major driver for adding hills and computing forces due to the history
!>        dependent term
!> \par History
!>      04.2004 created
!>      10.2008 Teodoro Laino [tlaino] - University of Zurich
!>              Major rewriting and addition of multiple walkers
! *****************************************************************************
  SUBROUTINE hills(meta_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hills', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, i_hills, ih, iter_nr, iw, &
                                                n_colvar, n_hills_start, &
                                                n_step, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: dfunc, diff_ss, dp2
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ddp
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(hills_env_type), POINTER            :: hills_env
    TYPE(metavar_type), DIMENSION(:), &
      POINTER                                :: colvars
    TYPE(multiple_walkers_type), POINTER     :: multiple_walkers

    failure   = .FALSE.
    NULLIFY(hills_env, multiple_walkers, logger, colvars, ddp)
    hills_env => meta_env%hills_env
    logger    => cp_error_get_logger(error)
    colvars   => meta_env%metavar
    n_colvar  =  meta_env%n_colvar
    n_step    =  meta_env%n_steps

    ! Create a temporary logger level specific for metadynamics
    CALL cp_add_iter_level(logger%iter_info,"METADYNAMICS",error=error)
    CALL get_meta_iter_level(meta_env, iter_nr, error)
    CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=iter_nr,error=error)   
 
    ! Set-up restart if any
    IF (meta_env%hills_env%restart) THEN
       meta_env%hills_env%restart = .FALSE.
       CALL restart_hills(hills_env%ss_history, hills_env%delta_s_history, hills_env%ww_history,&
            hills_env%ww, hills_env%n_hills, n_colvar, colvars, meta_env%metadyn_section, error)
    END IF

    ! Proceed with normal calculation
    IF ((MODULO(n_step,hills_env%nt_hills)==0).AND.(.NOT.meta_env%restart)) THEN
       IF (meta_env%do_multiple_walkers) multiple_walkers => meta_env%multiple_walkers

       n_hills_start = hills_env%n_hills
       ! Add the hill corresponding to this location
       CALL add_hill_single(hills_env, colvars, hills_env%ww, hills_env%n_hills, n_colvar, error)
       ! Update local n_hills counter
       IF (meta_env%do_multiple_walkers) multiple_walkers%n_hills_local = multiple_walkers%n_hills_local+1

       ! Update iteration level for printing
       CALL get_meta_iter_level(meta_env, iter_nr, error)
       CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=iter_nr,error=error)

       ! Print just program_run_info
       iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
            "PRINT%PROGRAM_RUN_INFO",extension=".metadynLog",error=error)
       IF (iw>0) THEN
          IF (meta_env%do_multiple_walkers) THEN
             WRITE(iw,'(/,1X,"METADYN|",A,I0,A,I0,A,/)')&
                  ' Global/Local Hills number (',hills_env%n_hills,'/',multiple_walkers%n_hills_local,&
                  ') added.'
          ELSE
             WRITE(iw,'(/,1X,"METADYN|",A,I0,A,/)')' Hills number  (',hills_env%n_hills,') added.'
          END IF
       END IF
       CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)

       ! Handle Multiple Walkers
       IF (meta_env%do_multiple_walkers) THEN
          ! Print Local Hills file if requested
          iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
               "PRINT%HILLS",middle_name="LOCAL",extension=".metadynLog",error=error)
          IF (iw>0) THEN
             WRITE(iw,'(f12.1,30f13.5)')meta_env%time*femtoseconds,&
                  (hills_env%ss_history(ih,hills_env%n_hills),ih=1,n_colvar),&
                  (hills_env%delta_s_history(ih,hills_env%n_hills),ih=1,n_colvar),&
                   hills_env%ww_history(hills_env%n_hills)
          END IF
          CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
               "PRINT%HILLS", error=error)

          ! Check the communicationa buffer of the other walkers
          CALL synchronize_multiple_walkers(multiple_walkers, hills_env, colvars,&
               n_colvar, meta_env%metadyn_section, error)
       END IF

       ! Print Hills file if requested (for multiple walkers this file includes
       ! the Hills coming from all the walkers).
       iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
            "PRINT%HILLS",extension=".metadynLog",error=error)
       IF (iw>0) THEN
          DO i_hills = n_hills_start+1, hills_env%n_hills
             WRITE(iw,'(f12.1,30f13.5)')meta_env%time*femtoseconds,&
                  (hills_env%ss_history(ih,i_hills),ih=1,n_colvar),&
                  (hills_env%delta_s_history(ih,i_hills),ih=1,n_colvar),&
                   hills_env%ww_history(i_hills)
          END DO
       END IF
       CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
            "PRINT%HILLS", error=error)
    END IF

    ! Computes forces due to the hills: history dependent term
    ALLOCATE(ddp(meta_env%n_colvar), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    DO ih=1,n_colvar
       colvars(ih)%ff_hills = 0.0_dp
    ENDDO
    DO ih=1,hills_env%n_hills
       dp2 = 0._dp
       DO i = 1, n_colvar
          diff_ss = colvars(i)%ss0-hills_env%ss_history(i,ih)
          IF (colvars(i)%periodic) THEN
             ! The difference of a periodic COLVAR is always within [-pi,pi]
             diff_ss = SIGN(1.0_dp,ASIN(SIN(diff_ss)))*ACOS(COS(diff_ss))
          END IF
          ddp(i) = (diff_ss)/hills_env%delta_s_history(i,ih)
          dp2 = dp2 + ddp(i)**2
       ENDDO
       ! ff_hills contains the "force" due to the hills
       dfunc=-hills_env%ww_history(ih)*EXP(-0.5_dp*dp2)
       DO i=1,n_colvar
          colvars(i)%ff_hills=colvars(i)%ff_hills-ddp(i)/hills_env%delta_s_history(i,ih)*dfunc
       ENDDO
    ENDDO
    DEALLOCATE(ddp, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CALL cp_rm_iter_level(logger%iter_info,"METADYNAMICS",error=error)    
  END SUBROUTINE hills

! *****************************************************************************
!> \brief  Synchronize with the rest of the walkers
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
! *****************************************************************************
  SUBROUTINE synchronize_multiple_walkers(multiple_walkers, hills_env, colvars,&
       n_colvar, metadyn_section, error)
    TYPE(multiple_walkers_type), POINTER     :: multiple_walkers
    TYPE(hills_env_type), POINTER            :: hills_env
    TYPE(metavar_type), DIMENSION(:), &
      POINTER                                :: colvars
    INTEGER, INTENT(IN)                      :: n_colvar
    TYPE(section_vals_type), POINTER         :: metadyn_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'synchronize_multiple_walkers', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: filename, tmpname
    INTEGER                                  :: delta_hills, i, i_hills, ih, &
                                                iw, stat, unit_nr
    LOGICAL                                  :: exist, failure
    REAL(KIND=dp)                            :: ww
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: delta_s_save, ss0_save
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    logger  => cp_error_get_logger(error)

    ! Locally dump information on file..
    IF (logger%para_env%mepos==logger%para_env%mepos) THEN
       ! Generate file name for the specific Hill
       i = multiple_walkers%walker_id
       filename = TRIM(multiple_walkers%walkers_file_name(i))//"_"//&
                  TRIM(ADJUSTL(cp_to_string(multiple_walkers%n_hills_local)))
       tmpname  = TRIM(filename)//".tmp"
       CALL open_file(file_name=tmpname,file_status="UNKNOWN",&
            file_form="FORMATTED",file_action="WRITE",&
            file_position="APPEND",unit_number=unit_nr)
       WRITE(unit_nr,*)hills_env%ww_history(hills_env%n_hills)
       DO ih = 1, n_colvar
          WRITE(unit_nr,*)hills_env%ss_history(ih,hills_env%n_hills)
          WRITE(unit_nr,*)hills_env%delta_s_history(ih,hills_env%n_hills)
       END DO
       CALL close_file(unit_nr)
       CALL m_mov(tmpname, filename)
    END IF
    
    IF (MODULO(multiple_walkers%n_hills_local,multiple_walkers%walkers_freq_comm)==0) THEN
       ! Store colvars information
       ALLOCATE(ss0_save(n_colvar),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(delta_s_save(n_colvar),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DO i = 1, n_colvar
          ss0_save(i)     = colvars(i)%ss0
          delta_s_save(i) = colvars(i)%delta_s
       END DO

       ! Watch for other walkers's file and update
       DO i = 1, multiple_walkers%walkers_tot_nr
          IF (i==multiple_walkers%walker_id) THEN
             ! Update local counter
             multiple_walkers%walkers_status(i) = multiple_walkers%n_hills_local
             CYCLE
          END IF
          
          i_hills  = multiple_walkers%walkers_status(i) + 1
          filename = TRIM(multiple_walkers%walkers_file_name(i))//"_"//&
                     TRIM(ADJUSTL(cp_to_string(i_hills)))
          INQUIRE(FILE=TRIM(filename),EXIST=exist)
          DO WHILE (exist)
             ! Read information from the walker's file
             ! We shouldn't care too much about the concurrency of these I/O instructions..
             ! In case, they can be fixed in the future..
             CALL open_file(file_name=filename,file_status="OLD",&
                  file_form="FORMATTED",file_action="READ",&
                  file_position="REWIND",unit_number=unit_nr)
             READ(unit_nr,*)ww
             DO ih = 1, n_colvar
                READ(unit_nr,*)colvars(ih)%ss0
                READ(unit_nr,*)colvars(ih)%delta_s
             END DO
             CALL close_file(unit_nr)

             ! Add this hill to the history dependent terms
             CALL add_hill_single(hills_env, colvars, ww, hills_env%n_hills, n_colvar, error)

             i_hills  = i_hills + 1
             filename = TRIM(multiple_walkers%walkers_file_name(i))//"_"//&
                        TRIM(ADJUSTL(cp_to_string(i_hills)))
             INQUIRE(FILE=TRIM(filename),EXIST=exist)
          END DO

          delta_hills = i_hills - 1 - multiple_walkers%walkers_status(i)
          multiple_walkers%walkers_status(i) = i_hills - 1
          iw = cp_print_key_unit_nr(logger,metadyn_section,"PRINT%PROGRAM_RUN_INFO",&
               extension=".metadynLog",error=error)
          IF (iw>0) THEN
             WRITE(iw,'(T2,A,I0,A,I0,A,I0,A)')'WALKERS| Walker #',i,'. Reading [',delta_hills,&
                  '] Hills. Total number of Hills acquired [',multiple_walkers%walkers_status(i),']'
          END IF
          CALL cp_print_key_finished_output(iw,logger,metadyn_section,&
               "PRINT%PROGRAM_RUN_INFO", error=error)
       END DO

       ! Restore colvars information
       DO i = 1, n_colvar
          colvars(i)%ss0     = ss0_save(i)
          colvars(i)%delta_s = delta_s_save(i)
       END DO
       DEALLOCATE(ss0_save,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(delta_s_save,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    END IF
  END SUBROUTINE synchronize_multiple_walkers

! *****************************************************************************
!> \brief Add a single Hill
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
! *****************************************************************************
  SUBROUTINE add_hill_single(hills_env, colvars, ww, n_hills, n_colvar, error)
    TYPE(hills_env_type), POINTER            :: hills_env
    TYPE(metavar_type), DIMENSION(:), &
      POINTER                                :: colvars
    REAL(KIND=dp), INTENT(IN)                :: ww
    INTEGER, INTENT(INOUT)                   :: n_hills
    INTEGER, INTENT(IN)                      :: n_colvar
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'add_hill_single', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tnp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: tmp

    failure = .FALSE.
    NULLIFY(tmp, tnp)
    IF(SIZE(hills_env%ss_history,2)< n_hills+1)THEN
       ALLOCATE(tmp(n_colvar,n_hills+100), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       tmp(:,:n_hills)=hills_env%ss_history
       tmp(:,n_hills+1:)=0.0_dp
       DEALLOCATE(hills_env%ss_history,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       hills_env%ss_history => tmp
       NULLIFY(tmp)
    ENDIF
    IF(SIZE(hills_env%delta_s_history,2)< n_hills+1)THEN
       ALLOCATE(tmp(n_colvar,n_hills+100), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       tmp(:,:n_hills)=hills_env%delta_s_history
       tmp(:,n_hills+1:)=0.0_dp
       DEALLOCATE(hills_env%delta_s_history,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       hills_env%delta_s_history => tmp
       NULLIFY(tmp)
    ENDIF
    IF(SIZE(hills_env%ww_history)< n_hills+1)THEN
       ALLOCATE(tnp(n_hills+100), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       tnp(1:n_hills) = hills_env%ww_history
       tnp(n_hills+1:)= 0.0_dp
       DEALLOCATE(hills_env%ww_history,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       hills_env%ww_history => tnp
       NULLIFY(tnp)
    ENDIF
    n_hills = n_hills + 1
    ! Now add the hill
    DO i=1,n_colvar
       hills_env%ss_history(i,n_hills)      = colvars(i)%ss0
       hills_env%delta_s_history(i,n_hills) = colvars(i)%delta_s
       hills_env%ww_history(n_hills)        = ww
    ENDDO
  END SUBROUTINE add_hill_single

! *****************************************************************************
!> \brief Restart Hills Information 
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
! *****************************************************************************
  SUBROUTINE restart_hills(ss_history, delta_s_history, ww_history, ww, n_hills,&
       n_colvar, colvars, metadyn_section, error)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ss_history, delta_s_history
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ww_history
    REAL(KIND=dp)                            :: ww
    INTEGER, INTENT(IN)                      :: n_hills, n_colvar
    TYPE(metavar_type), DIMENSION(:), &
      POINTER                                :: colvars
    TYPE(section_vals_type), POINTER         :: metadyn_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'restart_hills', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, ndum, stat
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp)                            :: rval
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rvals
    TYPE(section_vals_type), POINTER         :: hills_history

    failure = .FALSE.
    NULLIFY(rvals)
    hills_history => section_vals_get_subs_vals(metadyn_section,"SPAWNED_HILLS_POS",error=error)
    CALL section_vals_get(hills_history, explicit=explicit, error=error)
    IF (explicit) THEN
       CALL section_vals_val_get( hills_history,"_DEFAULT_KEYWORD_",n_rep_val=ndum, error=error)
       ! ss_history, delta_s_history, ww_history : deallocate and reallocate with the proper size
       DEALLOCATE(ss_history,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(delta_s_history,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(ww_history,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       !
       CPPostcondition(n_hills==ndum,cp_failure_level,routinep,error,failure)
       ALLOCATE(ss_history(n_colvar,n_hills), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(delta_s_history(n_colvar,n_hills), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(ww_history(n_hills), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       !
       DO i = 1, n_hills
          CALL section_vals_val_get(hills_history,"_DEFAULT_KEYWORD_",&
               i_rep_val=i, r_vals=rvals, error=error)
          CPPostcondition(SIZE(rvals)==n_colvar,cp_failure_level,routinep,error,failure)
          ss_history(1:n_colvar,i)=rvals
       END DO
       !
       hills_history => section_vals_get_subs_vals(metadyn_section,"SPAWNED_HILLS_AMPL",error=error)
       CALL section_vals_get(hills_history, explicit=explicit, error=error)
       IF (explicit) THEN
          ! delta_s_history
          CALL section_vals_val_get( hills_history,"_DEFAULT_KEYWORD_",n_rep_val=ndum, error=error)
          CPPostcondition(n_hills==ndum,cp_failure_level,routinep,error,failure)
          DO i = 1, n_hills
             CALL section_vals_val_get(hills_history,"_DEFAULT_KEYWORD_",&
                  i_rep_val=i, r_vals=rvals, error=error)
             CPPostcondition(SIZE(rvals)==n_colvar,cp_failure_level,routinep,error,failure)
             delta_s_history(1:n_colvar,i)=rvals
          END DO
       ELSE
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "Section SPAWNED_HILLS_AMPL is not present! Setting the amplitudes of the"//&
               "restarted hills according the parameters specified in the input file. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          DO i = 1, n_hills
             DO j = 1, n_colvar
                delta_s_history(j,i) = colvars(i)%delta_s
             END DO
          END DO
       END IF
       !
       hills_history => section_vals_get_subs_vals(metadyn_section,"SPAWNED_HILLS_HEIGHT",error=error)
       CALL section_vals_get(hills_history, explicit=explicit, error=error)
       IF (explicit) THEN
          ! ww_history
          CALL section_vals_val_get( hills_history,"_DEFAULT_KEYWORD_",&
               n_rep_val=ndum, error=error)
          CPPostcondition(n_hills==ndum,cp_failure_level,routinep,error,failure)
          DO i = 1, n_hills
             CALL section_vals_val_get(hills_history,"_DEFAULT_KEYWORD_",&
                  i_rep_val=i, r_val=rval, error=error)
             CPPostcondition(SIZE(rvals)==n_colvar,cp_failure_level,routinep,error,failure)
             ww_history(i)=rval
          END DO
       ELSE
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
               "Section SPAWNED_HILLS_HEIGHT is not present! Setting the height of the"//&
               "restarted hills according the parameters specified in the input file. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          ww_history = ww
       END IF
    ENDIF

  END SUBROUTINE restart_hills
  
! *****************************************************************************
!> \brief Retrieves the iteration level for the metadynamics loop
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
! *****************************************************************************
  SUBROUTINE get_meta_iter_level(meta_env, iter_nr, error)
    TYPE(meta_env_type), POINTER             :: meta_env
    INTEGER, INTENT(OUT)                     :: iter_nr
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_meta_iter_level', &
      routineP = moduleN//':'//routineN

    IF (meta_env%do_multiple_walkers) THEN
       iter_nr = meta_env%multiple_walkers%n_hills_local
    ELSE
       iter_nr = meta_env%hills_env%n_hills
    END IF

  END SUBROUTINE get_meta_iter_level

! *****************************************************************************
!> \par History
!>      11.2007 [created] [tlaino]
!> \author Teodoro Laino - University of Zurich - 11.2007
! *****************************************************************************
  SUBROUTINE walls(meta_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'walls', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ih, iwall
    LOGICAL                                  :: failure
    REAL(dp)                                 :: ddp, delta_s, dfunc, diff_ss, &
                                                dp2, ww
    TYPE(metavar_type), DIMENSION(:), &
      POINTER                                :: colvars

    failure = .FALSE.
    colvars => meta_env%metavar
    ! Forces from the Walls
    DO ih=1,SIZE(colvars)
       IF (colvars(ih)%do_wall) THEN 
          DO iwall = 1, SIZE(colvars(ih)%walls)
             SELECT CASE(colvars(ih)%walls(iwall)%id_type)
             CASE(do_wall_reflective,do_wall_none)
                ! Do Nothing.. treated in the main metadyn function
                CYCLE
             CASE(do_wall_quadratic)
                diff_ss = colvars(ih)%ss0-colvars(ih)%walls(iwall)%pos
                IF (colvars(ih)%periodic) THEN
                   ! The difference of a periodic COLVAR is always within [-pi,pi]
                   diff_ss = SIGN(1.0_dp,ASIN(SIN(diff_ss)))*ACOS(COS(diff_ss))
                END IF
                dfunc = 2.0_dp*colvars(ih)%walls(iwall)%k_quadratic*diff_ss
                SELECT CASE (colvars(ih)%walls(iwall)%id_direction)
                CASE(do_wall_p)
                   IF (diff_ss>0.0_dp) colvars(ih)%ff_hills=colvars(ih)%ff_hills-dfunc
                CASE(do_wall_m)
                   IF (diff_ss<0.0_dp) colvars(ih)%ff_hills=colvars(ih)%ff_hills-dfunc
                END SELECT
             CASE(do_wall_gaussian)
                diff_ss = colvars(ih)%ss0-colvars(ih)%walls(iwall)%pos
                IF (colvars(ih)%periodic) THEN
                   ! The difference of a periodic COLVAR is always within [-pi,pi]
                   diff_ss = SIGN(1.0_dp,ASIN(SIN(diff_ss)))*ACOS(COS(diff_ss))
                END IF
                ww      = colvars(ih)%walls(iwall)%ww_gauss
                delta_s = colvars(ih)%walls(iwall)%sigma_gauss
                ddp = (diff_ss)/delta_s
                dp2 = ddp**2
                dfunc=-ww*EXP(-0.5_dp*dp2)
                colvars(ih)%ff_hills=colvars(ih)%ff_hills-ddp/delta_s*dfunc
             END SELECT
          END DO
       END IF
    ENDDO
  END SUBROUTINE walls

END MODULE metadynamics
