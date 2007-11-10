!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/metadynamics *
!!
!!   NAME
!!     metadynamics
!!
!!   FUNCTION
!!     Performs the metadynamics calculation
!!
!!   NOTES
!!
!!
!!   AUTHOR
!!     
!!
!!   MODIFICATION HISTORY
!!     01.2005 created [fawzi and ale]
!!     11.2007 Teodoro Laino [tlaino] - University of Zurich
!!
!!   SOURCE
!!
!****************************************************************************
MODULE metadynamics
  USE colvar_methods,                  ONLY: colvar_eval_glob_f
  USE colvar_types,                    ONLY: colvar_p_type,&
                                             torsion_colvar_id
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             cp_rm_iter_level
  USE cp_para_env,                     ONLY: cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_constants,                 ONLY: do_fe_meta,&
                                             do_wall_gaussian,&
                                             do_wall_m,&
                                             do_wall_none,&
                                             do_wall_p,&
                                             do_wall_quadratic,&
                                             do_wall_reflective
  USE input_cp2k_motion,               ONLY: create_metavar_section
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
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE metadynamics_types,              ONLY: hills_env_type,&
                                             meta_env_type,&
                                             metavar_type
  USE parallel_rng_types,              ONLY: next_random_number
  USE particle_list_types,             ONLY: particle_list_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE 

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'metadynamics'

  INTEGER, SAVE, PRIVATE :: last_meta_env_id=0

  PUBLIC :: metadyn_create, metadyn_read, metadyn_forces

CONTAINS
  
!!****f* metadynamics/metadyn_create *
!!
!!   NAME
!!     metadyn_create
!!
!!   FUNCTION
!!     allocates a metadynamic environment (performs only minimal 
!!     initialization)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - meta_env: the meta env_ that will be allocated
!!     - n_colvar: number of collectiva variables
!!     - do_hills: if the hills should be added
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!     Teodoro Laino [tlaino] - University of Zurich. 11.2007
!!
!!   MODIFICATION HISTORY
!!     04.2004 created
!!     02.2006 Reorganized the structure of the restart for Metadynamics (teo)
!!             cleaned the metadynamic type
!!
!!*** **********************************************************************
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
    IF (.NOT. failure) THEN
       NULLIFY(meta_env%metadyn_section)
       meta_env%metadyn_section => metadyn_section
       meta_env%ref_count=1
       meta_env%restart = .TRUE.
       last_meta_env_id=last_meta_env_id+1
       meta_env%id_nr=last_meta_env_id
       meta_env%n_colvar=n_colvar
       meta_env%para_env => para_env
       CALL cp_para_env_retain(para_env,error=error)

       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       meta_env%ekin_s=0.0_dp
       meta_env%cv_energy=0.0_dp

       meta_env%n_steps=0
       meta_env%dt = dt
       meta_env%tempcontrol=.FALSE.

       ! hills_env
       NULLIFY(meta_env%hills_env%ss_history,&
               meta_env%hills_env%delta_s_history,&
               meta_env%hills_env%ww_history)
       ALLOCATE(meta_env%hills_env%ss_history(n_colvar,0), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(meta_env%hills_env%delta_s_history(n_colvar,0), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       ALLOCATE(meta_env%hills_env%ww_history(0), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       meta_env%hills_env%n_hills=0
       meta_env%hills_env%restart = .TRUE.

       ! colvar
       ALLOCATE(meta_env%metavar(n_colvar), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
       DO i=1,n_colvar
          NULLIFY(meta_env%metavar(i)%walls)
          meta_env%metavar(i)%mass=-HUGE(0.0_dp)
          meta_env%metavar(i)%lambda=-HUGE(0.0_dp)
          meta_env%metavar(i)%ss=0.0_dp
          meta_env%metavar(i)%ss0=0.0_dp
          meta_env%metavar(i)%ff_s=0.0_dp
          meta_env%metavar(i)%vvp=0.0_dp
          meta_env%metavar(i)%cv_energy=0.0_dp
          meta_env%metavar(i)%delta_s=0.0_dp
          meta_env%metavar(i)%ff_hills=0.0_dp
          meta_env%metavar(i)%do_wall=.FALSE.
          meta_env%metavar(i)%periodic=.FALSE.
          meta_env%metavar(i)%icolvar=0
       END DO
    END IF
  END SUBROUTINE metadyn_create
!***************************************************************************

!!****f* metadynamics/metadyn_read *
!!
!!   NAME
!!     metadyn_read
!!
!!   FUNCTION
!!     reads metadynamics section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!    Teodoro Laino [tlaino] - University of Zurich. 11.2007 
!!
!!   MODIFICATION HISTORY
!!     04.2004 created
!!
!!*** **********************************************************************
  SUBROUTINE metadyn_read(meta_env,root_section,para_env,fe_section,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: fe_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_read', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, id_method, iw, j, &
                                                n_colvar, n_rep
    LOGICAL                                  :: any_reflective, explicit
    REAL(kind=dp)                            :: dt
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(enumeration_type), POINTER          :: enum
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: section, wall_section, &
                                                work_section
    TYPE(section_vals_type), POINTER         :: md_section, metadyn_section, &
                                                metavar_section

    logger => cp_error_get_logger(error)
    any_reflective = .FALSE.
    CALL section_vals_get(fe_section,explicit=explicit,error=error)
    IF (.NOT.explicit) RETURN
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
    CALL section_vals_val_get(metadyn_section,"TEMP_TOL",n_rep_val=n_rep,error=error)
    IF ( n_rep /= 0 ) THEN
       meta_env%tempcontrol=.TRUE.
       CALL section_vals_val_get(metadyn_section,"TEMP_TOL",r_val=meta_env%toll_temp,error=error)
    END IF

    DO i= 1, n_colvar
       CALL metavar_read(meta_env%metavar(i), meta_env%extended_lagrange,i, metavar_section,error=error)
    END DO

    iw = cp_print_key_unit_nr(logger,metadyn_section,&
         "PRINT%PROGRAM_RUN_INFO",extension=".metadynLog",error=error)
    NULLIFY(section,enum,keyword)
    CALL create_metavar_section(section,error=error)
    wall_section => section_get_subsection(section,"WALL",error)
    IF (iw>0) THEN
       WRITE ( iw, '( /A )' )          ' METADYN| Meta Dynamics Protocol '
       WRITE ( iw, '( A,T71,I10)' )    ' METADYN| Number of interval time steps to spawn hills',&
            meta_env%hills_env%nt_hills
       WRITE ( iw, '( A,T71,I10)' )    ' METADYN| Number of previously spawned hills',&
            meta_env%hills_env%n_hills
       IF (meta_env%extended_lagrange) THEN
          WRITE ( iw, '( A )' )        ' METADYN| Extended Lagrangian Scheme '
          IF (meta_env%tempcontrol) WRITE ( iw, '( A,T71,F10.2)' )&
               ' METADYN| Collective Variables Temperature control',meta_env%toll_temp
       END IF
       IF (meta_env%do_hills) &
            WRITE ( iw, '( A )' )      ' METADYN| Spanwing the Hills '
       WRITE ( iw, '( A,T71,F10.2)' )  ' METADYN| Width of the Spawned Gaussian',meta_env%hills_env%ww
       WRITE ( iw, '( A,T71,I10)' )    ' METADYN| Number of collective variables',meta_env%n_colvar

       DO i=1,n_colvar
          WRITE ( iw, '( A )')'          ----------------------------------------------------------------------'
          WRITE ( iw, '( A,T71,I10)' )      ' METAVARS| Collective Variable Number',meta_env%metavar(i)%icolvar
          IF (meta_env%extended_lagrange) THEN
             WRITE ( iw, '( A,T71,F10.6)' )    ' METAVARS| Lambda Parameter',meta_env%metavar(i)%lambda
             WRITE ( iw, '( A,T66,F15.6)' )    ' METAVARS| Collective Variable Mass',meta_env%metavar(i)%mass
          END IF
          WRITE ( iw, '( A,T71,F10.6)' )    ' METAVARS| Scaling factor',meta_env%metavar(i)%delta_s
          IF (meta_env%metavar(i)%do_wall) THEN
             WRITE ( iw, '( A,T71,I10)' )    ' METAVARS| Number of Walls present',SIZE(meta_env%metavar(i)%walls)
             DO j = 1, SIZE(meta_env%metavar(i)%walls)
                keyword => section_get_keyword(wall_section,"TYPE",error=error)
                CALL keyword_get(keyword,enum=enum,error=error)
                WRITE ( iw, '(/,A,5X,I10,T50,A,T70,A11)' )    ' METAVARS| Wall Number:',j,'Type of Wall:',&
                     ADJUSTR(TRIM(enum_i2c(enum,meta_env%metavar(i)%walls(j)%id_type,error=error)))
                ! Type of wall IO
                SELECT CASE(meta_env%metavar(i)%walls(j)%id_type)
                CASE(do_wall_none)
                   ! Do Nothing
                   CYCLE
                CASE(do_wall_reflective)
                   any_reflective = .TRUE.
                   work_section => section_get_subsection(wall_section,"REFLECTIVE",error)
                   keyword => section_get_keyword(work_section,"DIRECTION",error=error)
                   CALL keyword_get(keyword,enum=enum,error=error)
                   WRITE ( iw, '(A,T70,A11)' )    ' METAVARS| Wall direction',&
                        ADJUSTR(TRIM(enum_i2c(enum,meta_env%metavar(i)%walls(j)%id_direction,error=error)))
                CASE(do_wall_quadratic)
                   work_section => section_get_subsection(wall_section,"QUADRATIC",error)
                   keyword => section_get_keyword(work_section,"DIRECTION",error=error)
                   CALL keyword_get(keyword,enum=enum,error=error)
                   WRITE ( iw, '(A,T70,A11)' )    ' METAVARS| Wall direction',&
                        ADJUSTR(TRIM(enum_i2c(enum,meta_env%metavar(i)%walls(j)%id_direction,error=error)))
                   WRITE ( iw, '(A,T70,F11.6)' )    ' METAVARS| Constant K of the quadratic potential',&
                        meta_env%metavar(i)%walls(j)%k_quadratic
                CASE(do_wall_gaussian)
                   WRITE ( iw, '(A,T70,F11.6)' )    ' METAVARS| Height of the Wall Gaussian',&
                        meta_env%metavar(i)%walls(j)%ww_gauss
                   WRITE ( iw, '(A,T70,F11.6)' )    ' METAVARS| Amplitude of the Wall Gaussian',&
                        meta_env%metavar(i)%walls(j)%sigma_gauss 
                END SELECT
                WRITE ( iw, '(A,T70,F11.6)' )    ' METAVARS| Wall location',&
                     meta_env%metavar(i)%walls(j)%pos
             END DO
          END IF
          WRITE ( iw, '( A )')'          ----------------------------------------------------------------------'
       ENDDO
    END IF
    CALL section_release(section,error=error)
    CALL cp_print_key_finished_output(iw,logger,metadyn_section,"PRINT%PROGRAM_RUN_INFO", error=error)

  END SUBROUTINE metadyn_read

!!****f* metadynamics/metavar_read *
!!
!!   NAME
!!     metavar_read
!!
!!   FUNCTION
!!     reads metavar section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!     alessandro laio and fawzi mohamed
!!     Teodoro Laino [tlaino] - University of Zurich. 11.2007
!!
!!   MODIFICATION HISTORY
!!     04.2004 created
!!
!!*** **********************************************************************  
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

!!****f* metadynamics/metadyn_forces *
!!
!!   NAME
!!     metadyn_forces
!!
!!   FUNCTION
!!     add forces to the subsys due to the metadynamics run
!!     possibly modifies the velocites (if reflective walls are applied)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!     04.2004 created
!!
!!*** **********************************************************************  
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
    REAL(kind=dp)                            :: dt, ekin_w, fac_t, fft, norm, &
                                                rval, scal, ss0_test, tol_ekin
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
          cv%cv_energy = 0.5_dp*cv%lambda*(cv%ss-cv%ss0)**2.0_dp
          cv%ff_s      = cv%lambda*(cv%ss-cv%ss0)
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
                   SELECT CASE(cv%walls(iwall)%id_direction)
                   CASE(do_wall_p)
                      IF (ss0_test > cv%walls(iwall)%pos) cv%vvp = -cv%vvp
                   CASE(do_wall_m)
                      IF (ss0_test < cv%walls(iwall)%pos) cv%vvp = -cv%vvp
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
                CASE(do_wall_m)
                   IF (cv%ss > cv%walls(iwall)%pos) CYCLE
                END SELECT
                NULLIFY(particles)
                icolvar=cv%icolvar
                CALL cp_subsys_get(subsys, colvar_p=colvar_p, particles=particles,error=error)
                scal=0.0_dp
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
                   norm=norm+colvar_p(icolvar)%colvar%dsdr(1,ii)**2
                   norm=norm+colvar_p(icolvar)%colvar%dsdr(2,ii)**2
                   norm=norm+colvar_p(icolvar)%colvar%dsdr(3,ii)**2
                ENDDO
                IF (norm /= 0.0_dp) THEN
                  scal = scal/norm
                ELSE 
                  scal = 0.0_dp
                END IF
                DO ii=1,colvar_p(icolvar)%colvar%n_atom_s
                   i=colvar_p(icolvar)%colvar%i_atom(ii)
                   IF (PRESENT(vel)) THEN
                      vel(:,i)=vel(:,i)-2.0_dp*colvar_p(icolvar)%colvar%dsdr(:,ii)*scal
                   ELSE
                      particles%els(i)%v(:)=particles%els(i)%v(:)-2.0_dp*colvar_p(icolvar)%colvar%dsdr(:,ii)*scal
                   END IF
                ENDDO
             END SELECT
          END DO
       END IF
    ENDDO

    ! write COLVAR file
    iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
         "PRINT%COLVAR",extension=".metadynLog",error=error)
    IF (iw>0) THEN
       WRITE(iw,'(f12.1,40f13.5)')dt*meta_env%n_steps, &
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

!!****f* metadynamics/hills *
!!
!!   NAME
!!     hills
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!
!!
!!   MODIFICATION HISTORY
!!     04.2004 created
!!
!!*** **********************************************************************  
  SUBROUTINE hills(meta_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hills', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ih, iter_level, iw, j, &
                                                n_colvar, n_step, ndum, stat
    INTEGER, POINTER                         :: n_hills
    LOGICAL                                  :: explicit, failure
    REAL(dp)                                 :: dfunc, diff_ss, dp2, rval, &
                                                tmp_ss
    REAL(dp), DIMENSION(:), POINTER          :: ddp
    REAL(dp), DIMENSION(:, :), POINTER       :: delta_s_history, ss_history, &
                                                tmp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rvals, tnp, ww_history
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(hills_env_type), POINTER            :: hills_env
    TYPE(metavar_type), DIMENSION(:), &
      POINTER                                :: colvars
    TYPE(section_vals_type), POINTER         :: hills_history

    failure=.FALSE.
    NULLIFY(ddp)
    hills_env       => meta_env%hills_env
    logger          => cp_error_get_logger(error)
    ss_history      => hills_env%ss_history
    delta_s_history => hills_env%delta_s_history
    ww_history      => hills_env%ww_history
    n_hills         => hills_env%n_hills
    colvars         => meta_env%metavar
    ! Create a temporary logger level specific for metadynamics
    CALL cp_add_iter_level(logger%iter_info,"METADYNAMICS",n_rlevel_new=iter_level,&
         error=error)
    CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=n_hills,error=error)   
    !
    ALLOCATE(ddp(meta_env%n_colvar), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    n_colvar = meta_env%n_colvar
    n_step   = meta_env%n_steps
    NULLIFY(tmp, tnp, hills_history, rvals)

    ! Set-up restart if any
    IF (meta_env%hills_env%restart) THEN
       meta_env%hills_env%restart = .FALSE.
       hills_history => section_vals_get_subs_vals(meta_env%metadyn_section,"SPAWNED_HILLS_POS",error=error)
       CALL section_vals_get(hills_history, explicit=explicit, error=error)
       IF (explicit) THEN
          CALL section_vals_val_get( hills_history,"_DEFAULT_KEYWORD_",n_rep_val=ndum, error=error)
          ! ss_history
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
          hills_env%ss_history=>ss_history
          !
          hills_history => section_vals_get_subs_vals(meta_env%metadyn_section,"SPAWNED_HILLS_AMPL",error=error)
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
             DO i = 1, n_hills
                DO j = 1, n_colvar
                   delta_s_history(j,i) = colvars(i)%delta_s
                END DO
             END DO
          END IF
          hills_env%delta_s_history=>delta_s_history
          !
          hills_history => section_vals_get_subs_vals(meta_env%metadyn_section,"SPAWNED_HILLS_HEIGHT",error=error)
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
             ww_history=hills_env%ww
          END IF
          hills_env%ww_history=>ww_history
       ENDIF
    END IF

    ! Proceed with normal calculation
    IF ((MODULO(n_step,hills_env%nt_hills)==0).AND.(.NOT.meta_env%restart)) THEN
       IF(SIZE(ss_history,2)< n_hills+1)THEN
          ALLOCATE(tmp(n_colvar,n_hills+100), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          tmp(:,:n_hills)=ss_history
          tmp(:,n_hills+1:)=0.0_dp
          DEALLOCATE(ss_history,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          hills_env%ss_history => tmp
          ss_history           => tmp
          NULLIFY(tmp)
       ENDIF
       IF(SIZE(delta_s_history,2)< n_hills+1)THEN
          ALLOCATE(tmp(n_colvar,n_hills+100), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          tmp(:,:n_hills)=delta_s_history
          tmp(:,n_hills+1:)=0.0_dp
          DEALLOCATE(delta_s_history,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          hills_env%delta_s_history => tmp
          delta_s_history           => tmp
          NULLIFY(tmp)
       ENDIF
       IF(SIZE(ww_history)< n_hills+1)THEN
          ALLOCATE(tnp(n_hills+100), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          tnp(1:n_hills) = ww_history
          tnp(n_hills+1:)= 0.0_dp
          DEALLOCATE(ww_history,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          hills_env%ww_history => tnp
          ww_history           => tnp
          NULLIFY(tnp)
       ENDIF
       n_hills = n_hills + 1
       CALL cp_iterate(logger%iter_info,last=.FALSE.,iter_nr=n_hills,error=error)      

       ! Now add the hill
       DO i=1,n_colvar
          ss_history(i,n_hills)      = colvars(i)%ss0
          delta_s_history(i,n_hills) = colvars(i)%delta_s
          ww_history(n_hills)        = hills_env%ww
       ENDDO
       ! Print just program_run_info
       iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
            "PRINT%PROGRAM_RUN_INFO",extension=".metadynLog",error=error)
       IF (iw>0) THEN
          WRITE(iw,'(/,1X,"METADYN|",A,I7,A,/)')' Hills number ',n_hills,' added.'
       END IF
       CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)
       ! Print Hills file
       iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
            "PRINT%HILLS",extension=".metadynLog",error=error)
       IF (iw>0) THEN
          WRITE(iw,'(f12.1,30f13.5)')&
               meta_env%dt*meta_env%n_steps,  &
               (hills_env%ss_history(ih,n_hills),ih=1,n_colvar), &
               (hills_env%delta_s_history(ih,n_hills),ih=1,n_colvar),&
               hills_env%ww_history(n_hills)
       END IF
       CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
            "PRINT%HILLS", error=error)
    END IF

    ! Forces from the hills
    DO ih=1,n_colvar
       colvars(ih)%ff_hills=0.0_dp
    ENDDO
    DO ih=1,n_hills
       dp2=0._dp
       DO i=1,n_colvar
          diff_ss = colvars(i)%ss0-ss_history(i,ih)
          IF (colvars(i)%periodic) THEN
             ! The difference of a periodic COLVAR is always within [0,pi]
             tmp_ss  = COS(diff_ss)
             diff_ss = ACOS(tmp_ss)
          END IF
          ddp(i)=(diff_ss)/delta_s_history(i,ih)
          dp2=dp2+ddp(i)**2
       ENDDO
       dfunc=-ww_history(ih)*EXP(-0.5_dp*dp2)
       DO i=1,n_colvar
          colvars(i)%ff_hills=colvars(i)%ff_hills-ddp(i)/delta_s_history(i,ih)*dfunc
       ENDDO
    ENDDO
    DEALLOCATE(ddp, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    CALL cp_rm_iter_level(logger%iter_info,n_rlevel_att=iter_level,&
         error=error)    
  END SUBROUTINE hills

!!****f* metadynamics/walls *
!!
!!   NAME
!!     walls
!!
!!   FUNCTION
!!     
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!
!!
!!   AUTHOR
!!     Teodoro Laino - University of Zurich - 11.2007
!!
!!   MODIFICATION HISTORY
!!     11.2007 [created] [tlaino]
!!
!!*** **********************************************************************  
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
                dfunc = 2.0_dp*colvars(ih)%walls(iwall)%k_quadratic*diff_ss
                SELECT CASE (colvars(ih)%walls(iwall)%id_direction)
                CASE(do_wall_p)
                   IF (diff_ss>0.0_dp) colvars(ih)%ff_hills=colvars(ih)%ff_hills-dfunc
                CASE(do_wall_m)
                   IF (diff_ss<0.0_dp) colvars(ih)%ff_hills=colvars(ih)%ff_hills-dfunc
                END SELECT
             CASE(do_wall_gaussian)
                diff_ss = colvars(ih)%ss0-colvars(ih)%walls(iwall)%pos
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
