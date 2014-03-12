!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Performs the metadynamics calculation
!> \par History
!>      01.2005 created [fawzi and ale]
!>      11.2007 Teodoro Laino [tlaino] - University of Zurich
! *****************************************************************************
MODULE metadynamics_utils
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE force_env_types,                 ONLY: force_env_type
  USE input_constants,                 ONLY: do_fe_meta,&
                                             do_wall_gaussian,&
                                             do_wall_m,&
                                             do_wall_none,&
                                             do_wall_p,&
                                             do_wall_quadratic,&
                                             do_wall_quartic,&
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
  USE machine,                         ONLY: m_mov
  USE message_passing,                 ONLY: mp_bcast
  USE metadynamics_types,              ONLY: hills_env_type,&
                                             meta_env_type,&
                                             metadyn_create,&
                                             metavar_type,&
                                             multiple_walkers_type
  USE physcon,                         ONLY: kelvin
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'metadynamics_utils'

  PUBLIC :: metadyn_read,&
            synchronize_multiple_walkers,&
            add_hill_single,&
            restart_hills,&
            get_meta_iter_level,&
            meta_walls

CONTAINS

! *****************************************************************************
!> \brief reads metadynamics section
!> \par History
!>      04.2004 created
!> \author Teodoro Laino [tlaino] - University of Zurich. 11.2007
! *****************************************************************************
  SUBROUTINE metadyn_read(meta_env,force_env,root_section,para_env,fe_section,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: fe_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_read', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=default_path_length)       :: walkers_file_name
    INTEGER                                  :: handle, i, id_method, &
                                                n_colvar, n_rep, &
                                                number_allocated_colvars, stat
    INTEGER, DIMENSION(:), POINTER           :: walkers_status
    LOGICAL                                  :: check, explicit, failure
    REAL(kind=dp)                            :: dt
    TYPE(section_vals_type), POINTER         :: md_section, metadyn_section, &
                                                metavar_section, &
                                                walkers_section

    CALL timeset(routineN,handle)

    failure = .FALSE.
    CALL section_vals_get(fe_section,explicit=explicit,error=error)
    IF (explicit) THEN
       number_allocated_colvars = 0
       IF (ASSOCIATED(force_env%subsys%colvar_p)) THEN
          number_allocated_colvars = SIZE(force_env%subsys%colvar_p)
       END IF
       CALL section_vals_val_get(fe_section,"METHOD",i_val=id_method,error=error)
       IF (id_method/=do_fe_meta) THEN
          CALL timestop(handle)
          RETURN
       ENDIF
       metadyn_section => section_vals_get_subs_vals(fe_section,"METADYN",error=error)
       CPPreconditionNoFail(.NOT.ASSOCIATED(meta_env),cp_failure_level,routineP,error)

       md_section => section_vals_get_subs_vals(root_section,"MOTION%MD",error=error)
       CALL section_vals_val_get(md_section,"TIMESTEP",r_val=dt,error=error)

       metavar_section => section_vals_get_subs_vals(metadyn_section,"METAVAR",error=error)
       CALL section_vals_get(metavar_section,n_repetition=n_colvar,error=error)
       CALL metadyn_create(meta_env,n_colvar=n_colvar,&
            dt=dt,para_env=para_env,metadyn_section=metadyn_section,error=error)

       !Check if using plumed. If so, only get the file name and read nothing else
       CALL section_vals_val_get(metadyn_section,"USE_PLUMED",l_val=meta_env%use_plumed,error=error)
       IF (meta_env%use_plumed .EQV. .TRUE.) THEN
           CALL section_vals_val_get(metadyn_section,"PLUMED_INPUT_FILE",c_val=meta_env%plumed_input_file,error=error)
           meta_env%plumed_input_file=TRIM(meta_env%plumed_input_file)//CHAR(0)
           meta_env%langevin=.FALSE.
           CALL timestop(handle)
           RETURN
       END IF

       CALL section_vals_val_get(metadyn_section,"DO_HILLS",l_val=meta_env%do_hills,error=error)
       CALL section_vals_val_get(metadyn_section,"LAGRANGE",l_val=meta_env%extended_lagrange,error=error)
       CALL section_vals_val_get(metadyn_section,"TAMCSteps",i_val=meta_env%TAMCSteps,error=error)
       IF (meta_env%TAMCSteps<0) THEN
          CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,routineP,&
             "TAMCSteps must be positive!",only_ionode=.TRUE.)
       ENDIF
       CALL section_vals_val_get(metadyn_section,"Timestep",r_val=meta_env%zdt,error=error)
       IF (meta_env%zdt<=0.0_dp) THEN
          CALL cp_assert(.FALSE.,cp_fatal_level,cp_assertion_failed,routineP,&
             "Timestep must be positive!",only_ionode=.TRUE.)
       ENDIF
       CALL section_vals_val_get(metadyn_section,"WW",r_val=meta_env%hills_env%ww,error=error)
       CALL section_vals_val_get(metadyn_section,"NT_HILLS",i_val=meta_env%hills_env%nt_hills,error=error)
       CALL section_vals_val_get(metadyn_section,"MIN_NT_HILLS",i_val=meta_env%hills_env%min_nt_hills,error=error)
       IF (meta_env%hills_env%nt_hills<=0) THEN
          meta_env%hills_env%min_nt_hills=meta_env%hills_env%nt_hills
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                "NT_HILLS has a value <= 0; Setting MIN_NT_HILLS to the same value! Overriding input specification!",&
                only_ionode=.TRUE.)
       END IF
       check=meta_env%hills_env%nt_hills>=meta_env%hills_env%min_nt_hills
       CALL cp_assert(check,cp_fatal_level,cp_assertion_failed,routineP,&
             "MIN_NT_HILLS must have a value smaller or equal to NT_HILLS! Cross check with the input reference!",&
             only_ionode=.TRUE.)
       !RG Adaptive hills
       CALL section_vals_val_get(metadyn_section,"MIN_DISP",r_val=meta_env%hills_env%min_disp,error=error)
       CALL section_vals_val_get(metadyn_section,"OLD_HILL_NUMBER",i_val=meta_env%hills_env%old_hill_number,&
                           error=error)
       CALL section_vals_val_get(metadyn_section,"OLD_HILL_STEP",i_val=meta_env%hills_env%old_hill_step,&
                           error=error)

       !Hills tail damping
       CALL section_vals_val_get(metadyn_section,"HILL_TAIL_CUTOFF",r_val=meta_env%hills_env%tail_cutoff,error=error)
       CALL section_vals_val_get(metadyn_section,"P_EXPONENT",i_val=meta_env%hills_env%p_exp, error=error)
       CALL section_vals_val_get(metadyn_section,"Q_EXPONENT",i_val=meta_env%hills_env%q_exp, error=error)

       CALL section_vals_val_get(metadyn_section,"SLOW_GROWTH",l_val=meta_env%hills_env%slow_growth,error=error)

       !RG Adaptive hills
       CALL section_vals_val_get(metadyn_section,"STEP_START_VAL",i_val=meta_env%n_steps,error=error)
       CPPreconditionNoFail(meta_env%n_steps>=0,cp_failure_level,routineP,error)
       CALL section_vals_val_get(metadyn_section,"NHILLS_START_VAL",&
            i_val=meta_env%hills_env%n_hills,error=error)
       CALL section_vals_val_get(metadyn_section,"TEMPERATURE",r_val=meta_env%temp_wanted,error=error)
       CALL section_vals_val_get(metadyn_section,"LANGEVIN",l_val=meta_env%langevin,error=error)
       CALL section_vals_val_get(metadyn_section,"TEMP_TOL",explicit=meta_env%tempcontrol,&
            r_val=meta_env%toll_temp,error=error)
       CALL section_vals_val_get(metadyn_section,"WELL_TEMPERED",l_val=meta_env%well_tempered,error=error)
       CALL section_vals_val_get(metadyn_section,"DELTA_T",explicit=meta_env%hills_env%wtcontrol,&
&           r_val=meta_env%delta_t,error=error)
       CALL section_vals_val_get(metadyn_section,"WTGAMMA",explicit=check,&
&           r_val=meta_env%wtgamma,error=error)
       IF (meta_env%well_tempered) THEN
         meta_env%hills_env%wtcontrol=meta_env%hills_env%wtcontrol.OR.check
         check=meta_env%hills_env%wtcontrol
         CALL cp_assert(check,cp_fatal_level,cp_assertion_failed,routineP,&
               "When using Well-Tempered metadynamics, DELTA_T (or WTGAMMA) should be explicitly specified.",&
               only_ionode=.TRUE.)
         CALL cp_assert(.NOT.meta_env%extended_lagrange,cp_fatal_level,cp_assertion_failed,routineP,&
               "Well-Tempered metadynamics not possible with extended-lagrangian formulation.",&
               only_ionode=.TRUE.)
         CALL cp_assert(.NOT.(meta_env%hills_env%min_disp > 0.0_dp),cp_fatal_level,cp_assertion_failed,routineP,&
               "Well-Tempered metadynamics not possible with Adaptive hills.",&
               only_ionode=.TRUE.)
       ENDIF

       CALL section_vals_val_get(metadyn_section,"COLVAR_AVG_TEMPERATURE_RESTART",&
            r_val=meta_env%avg_temp,error=error)
       ! Parsing Metavar Section
       DO i= 1, n_colvar
          CALL metavar_read(meta_env%metavar(i), meta_env%extended_lagrange, &
               meta_env%langevin, i, metavar_section,error=error)
          check=(meta_env%metavar(i)%icolvar<=number_allocated_colvars)
          CALL cp_assert(check,cp_fatal_level,cp_assertion_failed,routineP,&
               "An error occured in the specification of COLVAR for METAVAR. "//&
               "Specified COLVAR #("//TRIM(ADJUSTL(cp_to_string(meta_env%metavar(i)%icolvar)))//") "//&
               "is larger than the maximum number of COLVARS defined in the SUBSYS ("//&
               TRIM(ADJUSTL(cp_to_string(number_allocated_colvars)))//") !"//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
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
               "number of Walkers File names provided. Please check your input and in case "//&
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

    CALL timestop(handle)

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
    INTEGER                                  :: handle, i, iw, j
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section, wall_section, &
                                                work_section

    CALL timeset(routineN,handle)

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
          IF (meta_env%langevin) THEN
               WRITE( iw, '(A,T71)') ' METADYN| Langevin Thermostat in use for COLVAR '
               WRITE( iw, '(A,T71,F10.4)') ' METADYN| Langevin Thermostat. Target Temperature = ',&
                    meta_env%temp_wanted*kelvin
          ENDIF
          WRITE( iw, '(A,T71,F10.4)') ' METADYN| COLVARS restarted average temperature ',&
               meta_env%avg_temp
       END IF
       IF (meta_env%do_hills) THEN
          WRITE (iw,'( A )')' METADYN| Spawning the Hills '
          WRITE (iw,'( A,T71,F10.3)')  ' METADYN| Height of the Spawned Gaussian',meta_env%hills_env%ww
          !RG Adaptive hills
          IF (meta_env%hills_env%min_disp.gt.0.0_dp) THEN
             WRITE (iw,'(A)')    ' METADYN| Adapative meta time step is activated'
             WRITE( iw, '(A,T71,F10.4)') ' METADYN| Minimum displacement for next hill',&
             meta_env%hills_env%min_disp
          END IF
          !RG Adaptive hills
       END IF

       IF (meta_env%well_tempered) THEN
          WRITE (iw,'( A )')' METADYN| Well-Tempered metadynamics '
          IF (meta_env%delta_t>EPSILON(1._dp)) THEN
             WRITE (iw,'( A,T71,F10.3)')  ' METADYN| Temperature parameter (Delta T) [K]',meta_env%delta_t*kelvin
          ELSE
             WRITE (iw,'( A,T71,F10.3)')  ' METADYN| Temperature parameter (gamma)',meta_env%wtgamma
          END IF
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
          IF (meta_env%langevin)  WRITE ( iw, '( A,T71,F10.6)' )  ' METAVARS| Friction for Langevin Thermostat',&
                                  meta_env%metavar(i)%gamma
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
                   WRITE (iw,'(A,T70,F11.6)')    ' METAVARS| Scale of the Wall Gaussian',&
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

    CALL timestop(handle)

  END SUBROUTINE print_metadyn_info

! *****************************************************************************
!> \brief reads metavar section
!> \par History
!>      04.2004 created
!> \author alessandro laio and fawzi mohamed
!>      Teodoro Laino [tlaino] - University of Zurich. 11.2007
! *****************************************************************************
  SUBROUTINE metavar_read(metavar,extended_lagrange,langevin,icol,metavar_section,error)
    TYPE(metavar_type), INTENT(INOUT)        :: metavar
    LOGICAL, INTENT(IN)                      :: extended_lagrange, langevin
    INTEGER, INTENT(IN)                      :: icol
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: metavar_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metavar_read', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, n_walls, stat
    LOGICAL                                  :: failure
    TYPE(section_vals_type), POINTER         :: wall_section, work_section

    CALL timeset(routineN,handle)

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
          CALL section_vals_val_get(wall_section,"POSITION",i_rep_section=i,r_val=metavar%walls(i)%pos,error=error)
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
          CASE(do_wall_quartic)
             work_section => section_vals_get_subs_vals(wall_section,"QUARTIC",i_rep_section=i,error=error)
             CALL section_vals_val_get(work_section,"DIRECTION",i_val=metavar%walls(i)%id_direction,error=error)
             CALL section_vals_val_get(work_section,"K",r_val=metavar%walls(i)%k_quartic,error=error)
             SELECT CASE(metavar%walls(i)%id_direction)
             CASE(do_wall_m)
               metavar%walls(i)%pos0= metavar%walls(i)%pos+(0.05_dp/metavar%walls(i)%k_quartic**(1/4))
             CASE(do_wall_p)
               metavar%walls(i)%pos0= metavar%walls(i)%pos-(0.05_dp/metavar%walls(i)%k_quartic**(1/4))
             END SELECT
          CASE(do_wall_gaussian)
             work_section => section_vals_get_subs_vals(wall_section,"GAUSSIAN",i_rep_section=i,error=error)
             CALL section_vals_val_get(work_section,"WW",r_val=metavar%walls(i)%ww_gauss,error=error)
             CALL section_vals_val_get(work_section,"SIGMA",r_val=metavar%walls(i)%sigma_gauss,error=error)
          END SELECT
       END DO
    END IF
    ! Setup few more parameters for extended lagrangian
    IF(extended_lagrange)THEN
       CALL section_vals_val_get(metavar_section,"MASS",i_rep_section=icol,r_val=metavar%mass,error=error)
       CALL section_vals_val_get(metavar_section,"LAMBDA",i_rep_section=icol,r_val=metavar%lambda,error=error)
       IF (langevin) THEN
          CALL section_vals_val_get(metavar_section,"GAMMA",i_rep_section=icol,r_val=metavar%gamma,error=error)
       END IF
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE metavar_read

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
    INTEGER                                  :: delta_hills, handle, i, &
                                                i_hills, ih, iw, stat, unit_nr
    LOGICAL                                  :: exist, failure
    REAL(KIND=dp)                            :: invdt, ww
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: delta_s_save, ss0_save
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: delta_s_ss0_buf
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env

    CALL timeset(routineN,handle)

    failure = .FALSE.
    logger  => cp_error_get_logger(error)
    para_env => logger%para_env

    ! Locally dump information on file..
    IF (para_env%mepos==para_env%source) THEN
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
       IF (hills_env%wtcontrol) WRITE(unit_nr,*)hills_env%invdt_history(hills_env%n_hills)
       CALL close_file(unit_nr)
       CALL m_mov(tmpname, filename)
    END IF

    IF (MODULO(multiple_walkers%n_hills_local,multiple_walkers%walkers_freq_comm)==0) THEN
       ! Store colvars information
       ALLOCATE(ss0_save(n_colvar),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(delta_s_save(n_colvar),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(delta_s_ss0_buf(2,0:n_colvar),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       delta_s_ss0_buf=0
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

          IF (para_env%mepos==para_env%source) THEN
             INQUIRE(FILE=TRIM(filename),EXIST=exist)
          ENDIF
          CALL mp_bcast(exist,para_env%source,para_env%group)
          DO WHILE (exist)
             ! Read information from the walker's file
             ! We shouldn't care too much about the concurrency of these I/O instructions..
             ! In case, they can be fixed in the future..
             IF (para_env%mepos==para_env%source) THEN
                CALL open_file(file_name=filename,file_status="OLD",&
                     file_form="FORMATTED",file_action="READ",&
                     file_position="REWIND",unit_number=unit_nr)
                READ(unit_nr,*) delta_s_ss0_buf(1,0)
                DO ih = 1, n_colvar
                   READ(unit_nr,*) delta_s_ss0_buf(1,ih)
                   READ(unit_nr,*) delta_s_ss0_buf(2,ih)
                END DO
                IF (hills_env%wtcontrol) READ(unit_nr,*) delta_s_ss0_buf(2,0)
                CALL close_file(unit_nr)
             ENDIF
             CALL mp_bcast(delta_s_ss0_buf,para_env%source,para_env%group)
             ww=delta_s_ss0_buf(1,0)
             IF (hills_env%wtcontrol) invdt=delta_s_ss0_buf(2,0)
             DO ih = 1, n_colvar
                colvars(ih)%ss0=delta_s_ss0_buf(1,ih)
                colvars(ih)%delta_s=delta_s_ss0_buf(2,ih)
             ENDDO

             ! Add this hill to the history dependent terms
             IF (hills_env%wtcontrol) THEN
                CALL add_hill_single(hills_env, colvars, ww, hills_env%n_hills, n_colvar, error, invdt=invdt)
             ELSE
                CALL add_hill_single(hills_env, colvars, ww, hills_env%n_hills, n_colvar, error)
             END IF

             i_hills  = i_hills + 1
             filename = TRIM(multiple_walkers%walkers_file_name(i))//"_"//&
                        TRIM(ADJUSTL(cp_to_string(i_hills)))
             IF (para_env%mepos==para_env%source) THEN
                INQUIRE(FILE=TRIM(filename),EXIST=exist)
             ENDIF
             CALL mp_bcast(exist,para_env%source,para_env%group)
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

    CALL timestop(handle)

  END SUBROUTINE synchronize_multiple_walkers

! *****************************************************************************
!> \brief Add a single Hill
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
! *****************************************************************************
  SUBROUTINE add_hill_single(hills_env, colvars, ww, n_hills, n_colvar, error, invdt)
    TYPE(hills_env_type), POINTER            :: hills_env
    TYPE(metavar_type), DIMENSION(:), &
      POINTER                                :: colvars
    REAL(KIND=dp), INTENT(IN)                :: ww
    INTEGER, INTENT(INOUT)                   :: n_hills
    INTEGER, INTENT(IN)                      :: n_colvar
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: invdt

    CHARACTER(len=*), PARAMETER :: routineN = 'add_hill_single', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, stat
    LOGICAL                                  :: failure, wtcontrol
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tnp
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: tmp

    CALL timeset(routineN,handle)

    failure = .FALSE.; wtcontrol=PRESENT(invdt)
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
    IF(wtcontrol)THEN
       IF(SIZE(hills_env%invdt_history)< n_hills+1)THEN
          ALLOCATE(tnp(n_hills+100), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          tnp(1:n_hills) = hills_env%invdt_history
          tnp(n_hills+1:)= 0.0_dp
          DEALLOCATE(hills_env%invdt_history,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          hills_env%invdt_history => tnp
          NULLIFY(tnp)
       ENDIF
    ENDIF
    n_hills = n_hills + 1
    ! Now add the hill
    DO i=1,n_colvar
       hills_env%ss_history(i,n_hills)      = colvars(i)%ss0
       hills_env%delta_s_history(i,n_hills) = colvars(i)%delta_s
    ENDDO
    hills_env%ww_history(n_hills)           = ww
    IF (wtcontrol) hills_env%invdt_history(n_hills) = invdt

    CALL timestop(handle)

  END SUBROUTINE add_hill_single

! *****************************************************************************
!> \brief Restart Hills Information
!> \author Teodoro Laino [tlaino] - University of Zurich - 10.2008
! *****************************************************************************
  SUBROUTINE restart_hills(ss_history, delta_s_history, ww_history, ww,&
       n_hills, n_colvar, colvars, metadyn_section, error, invdt_history)
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: ss_history, delta_s_history
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ww_history
    REAL(KIND=dp)                            :: ww
    INTEGER, INTENT(IN)                      :: n_hills, n_colvar
    TYPE(metavar_type), DIMENSION(:), &
      POINTER                                :: colvars
    TYPE(section_vals_type), POINTER         :: metadyn_section
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: invdt_history

    CHARACTER(len=*), PARAMETER :: routineN = 'restart_hills', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, ndum, stat
    LOGICAL                                  :: explicit, failure, wtcontrol
    REAL(KIND=dp)                            :: rval
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rvals
    TYPE(section_vals_type), POINTER         :: hills_history

    CALL timeset(routineN,handle)

    failure = .FALSE.; wtcontrol=PRESENT(invdt_history)
    NULLIFY(rvals)
    hills_history => section_vals_get_subs_vals(metadyn_section,"SPAWNED_HILLS_POS",error=error)
    CALL section_vals_get(hills_history, explicit=explicit, error=error)
    IF (explicit) THEN
       CALL section_vals_val_get( hills_history,"_DEFAULT_KEYWORD_",n_rep_val=ndum, error=error)
       ! ss_history, delta_s_history, ww_history, invdt_history : deallocate and reallocate with the proper size
       DEALLOCATE(ss_history,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(delta_s_history,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DEALLOCATE(ww_history,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       IF (wtcontrol) THEN
          DEALLOCATE(invdt_history,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ENDIF
       !
       CPPostcondition(n_hills==ndum,cp_failure_level,routinep,error,failure)
       ALLOCATE(ss_history(n_colvar,n_hills), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(delta_s_history(n_colvar,n_hills), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(ww_history(n_hills), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       IF (wtcontrol) THEN
          ALLOCATE(invdt_history(n_hills), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ENDIF
       !
       DO i = 1, n_hills
          CALL section_vals_val_get(hills_history,"_DEFAULT_KEYWORD_",&
               i_rep_val=i, r_vals=rvals, error=error)
          CPPostcondition(SIZE(rvals)==n_colvar,cp_failure_level,routinep,error,failure)
          ss_history(1:n_colvar,i)=rvals
       END DO
       !
       hills_history => section_vals_get_subs_vals(metadyn_section,"SPAWNED_HILLS_SCALE",error=error)
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
               "Section SPAWNED_HILLS_SCALE is not present! Setting the scales of the"//&
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
               " restarted hills according the parameters specified in the input file. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
          ww_history = ww
       END IF
       !
       hills_history => section_vals_get_subs_vals(metadyn_section,"SPAWNED_HILLS_INVDT",error=error)
       CALL section_vals_get(hills_history, explicit=explicit, error=error)
       IF (wtcontrol) THEN
          IF (explicit) THEN
             ! invdt_history
             CALL section_vals_val_get( hills_history,"_DEFAULT_KEYWORD_",&
                  n_rep_val=ndum, error=error)
             CPPostcondition(n_hills==ndum,cp_failure_level,routinep,error,failure)
             DO i = 1, n_hills
                CALL section_vals_val_get(hills_history,"_DEFAULT_KEYWORD_",&
                     i_rep_val=i, r_val=rval, error=error)
                CPPostcondition(SIZE(rvals)==n_colvar,cp_failure_level,routinep,error,failure)
                invdt_history(i)=rval
             END DO
          ELSE
             CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                  "Section SPAWNED_HILLS_INVDT is not present! Restarting from standard"//&
                  " metadynamics run i.e. setting 1/(Delta T) equal to zero. "//&
CPSourceFileRef,&
               only_ionode=.TRUE.)
             invdt_history = 0._dp
          END IF
       ELSE
          IF (explicit) THEN
             CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                     "Found section SPAWNED_HILLS_INVDT while restarting a standard metadynamics run..."//&
                     " Cannot restart metadynamics from well-tempered MetaD runs. "//&
CPSourceFileRef,&
                  only_ionode=.TRUE.)
          END IF
       END IF
    ENDIF

    CALL timestop(handle)

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
  SUBROUTINE meta_walls(meta_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'meta_walls', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ih, iwall
    LOGICAL                                  :: failure
    REAL(dp)                                 :: ddp, delta_s, dfunc, diff_ss, &
                                                dp2, efunc, ww
    TYPE(metavar_type), DIMENSION(:), &
      POINTER                                :: colvars

    failure = .FALSE.
    colvars => meta_env%metavar
    ! Forces from the Walls
    DO ih=1,SIZE(colvars)
       IF (colvars(ih)%do_wall) THEN
          colvars(ih)%epot_walls=0.0_dp
          colvars(ih)%ff_walls=0.0_dp
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
                efunc = colvars(ih)%walls(iwall)%k_quadratic*diff_ss**2
                dfunc = 2.0_dp*colvars(ih)%walls(iwall)%k_quadratic*diff_ss
                SELECT CASE (colvars(ih)%walls(iwall)%id_direction)
                CASE(do_wall_p)
                   IF (diff_ss>0.0_dp) THEN
                     colvars(ih)%ff_walls=colvars(ih)%ff_walls-dfunc
                     colvars(ih)%epot_walls=colvars(ih)%epot_walls+efunc
                   END IF
                CASE(do_wall_m)
                   IF (diff_ss<0.0_dp) THEN
                     colvars(ih)%ff_walls=colvars(ih)%ff_walls-dfunc
                     colvars(ih)%epot_walls=colvars(ih)%epot_walls+efunc
                   END IF
                END SELECT
             CASE(do_wall_quartic)
                diff_ss = colvars(ih)%ss0-colvars(ih)%walls(iwall)%pos0
                IF (colvars(ih)%periodic) THEN
                   ! The difference of a periodic COLVAR is always within [-pi,pi]
                   diff_ss = SIGN(1.0_dp,ASIN(SIN(diff_ss)))*ACOS(COS(diff_ss))
                END IF
                efunc = colvars(ih)%walls(iwall)%k_quartic*diff_ss*diff_ss**4
                dfunc = 4.0_dp*colvars(ih)%walls(iwall)%k_quartic*diff_ss**3
                SELECT CASE (colvars(ih)%walls(iwall)%id_direction)
                CASE(do_wall_p)
                   IF (diff_ss>0.0_dp) THEN
                     colvars(ih)%ff_walls=colvars(ih)%ff_walls-dfunc
                     colvars(ih)%epot_walls=colvars(ih)%epot_walls+efunc
                   END IF
                CASE(do_wall_m)
                   IF (diff_ss<0.0_dp) THEN
                     colvars(ih)%ff_walls=colvars(ih)%ff_walls-dfunc
                     colvars(ih)%epot_walls=colvars(ih)%epot_walls+efunc
                   END IF
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
                efunc = ww*EXP(-0.5_dp*dp2)
                dfunc = -efunc*ddp/delta_s
                colvars(ih)%ff_walls=colvars(ih)%ff_walls-dfunc
                colvars(ih)%epot_walls=colvars(ih)%epot_walls+efunc
             END SELECT
          END DO
       END IF
    ENDDO
  END SUBROUTINE meta_walls

END MODULE metadynamics_utils
