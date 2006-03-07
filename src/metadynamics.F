!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
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
!!
!!   SOURCE
!!
!****************************************************************************
MODULE metadynamics
  USE colvar_methods,                  ONLY: colvar_eval_f
  USE colvar_types,                    ONLY: colvar_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_env,                     ONLY: cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_flush
  USE metadynamics_types,              ONLY: hills_env_type,&
                                             meta_env_type,&
                                             metavar_type
  USE particle_list_types,             ONLY: particle_list_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE 

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'metadynamics'

  INTEGER, SAVE, PRIVATE :: last_meta_env_id=0

  PUBLIC :: metadyn_create, metadyn_read, metadyn_add_forces

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
       meta_env%ekin_s=HUGE(0.0_dp)
       meta_env%cv_energy=-HUGE(0.0_dp)
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
          meta_env%metavar(i)%ss=0.0_dp
          meta_env%metavar(i)%ss0=0.0_dp
          meta_env%metavar(i)%ff_s=0.0_dp
          meta_env%metavar(i)%vv=0.0_dp
          meta_env%metavar(i)%vvp=0.0_dp
          meta_env%metavar(i)%wall_p=0.0_dp
          meta_env%metavar(i)%wall_m=0.0_dp
          meta_env%metavar(i)%mass=-HUGE(0.0_dp)
          meta_env%metavar(i)%lambda=-HUGE(0.0_dp)
          meta_env%metavar(i)%icolvar=0
          meta_env%metavar(i)%isubsys=0
          meta_env%metavar(i)%do_wall_p=.FALSE.
          meta_env%metavar(i)%do_wall_m=.FALSE.
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
!!     
!!
!!   MODIFICATION HISTORY
!!     04.2004 created
!!
!!*** **********************************************************************
  SUBROUTINE metadyn_read(meta_env,globenv,metadyn_section,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(global_environment_type)            :: globenv
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: metadyn_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_read', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iw, n_colvar, n_rep
    LOGICAL                                  :: explicit
    REAL(kind=dp)                            :: dt
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: md_section, metavar_section

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    CALL section_vals_get(metadyn_section,explicit=explicit,error=error)
    IF (.NOT.explicit) RETURN
    CPPreconditionNoFail(.NOT.ASSOCIATED(meta_env),cp_failure_level,routineP,error)  

    md_section => section_vals_get_subs_vals(globenv%input_file,"MOTION%MD",error=error)
    CALL section_vals_val_get(md_section,"TIMESTEP",r_val=dt,error=error)

    metavar_section => section_vals_get_subs_vals(metadyn_section,"METAVAR",error=error)
    CALL section_vals_get(metavar_section,n_repetition=n_colvar,error=error)
    CALL metadyn_create(meta_env,n_colvar=n_colvar,&
         dt=dt,para_env=globenv%para_env,metadyn_section=metadyn_section,error=error)

    CALL section_vals_val_get(metadyn_section,"DO_HILLS",l_val=meta_env%do_hills,error=error)
    CALL section_vals_val_get(metadyn_section,"LAGRANGE",l_val=meta_env%extended_lagrange,error=error)
    CALL section_vals_val_get(metadyn_section,"WW",r_val=meta_env%hills_env%ww,error=error)
    CALL section_vals_val_get(metadyn_section,"NT_HILLS",i_val=meta_env%hills_env%nt_hills,error=error)
    CALL section_vals_val_get(metadyn_section,"STEP_START_VAL",i_val=meta_env%n_steps,error=error)
    CALL section_vals_val_get(metadyn_section,"NHILLS_START_VAL",i_val=meta_env%hills_env%n_hills,error=error)
    CALL section_vals_val_get(metadyn_section,"TEMPERATURE",n_rep_val=n_rep,error=error)
    IF ( n_rep /= 0 ) THEN
       meta_env%tempcontrol=.TRUE.
       CALL section_vals_val_get(metadyn_section,"TEMPERATURE",r_val=meta_env%temp_wanted,error=error)
       meta_env%toll_temp=meta_env%temp_wanted/SQRT(REAL(meta_env%n_colvar,dp))
    END IF

    DO i= 1, n_colvar
       CALL metavar_read(meta_env%metavar(i), meta_env%extended_lagrange,i, metavar_section,error=error)
    END DO

    iw = cp_print_key_unit_nr(logger,metadyn_section,&
         "PRINT%PROGRAM_RUN_INFO",extension=".metadynLog",error=error)
    IF (iw>0) THEN
       WRITE ( iw, '( /A )' )          ' METADYN| Meta Dynamics Protocol '
       WRITE ( iw, '( A,T71,I10)' )    ' METADYN| Number of interval time steps to spawn hills',&
            meta_env%hills_env%nt_hills
       WRITE ( iw, '( A,T71,I10)' )    ' METADYN| Number of previously spawned hills',&
            meta_env%hills_env%n_hills
       IF (meta_env%hills_env%n_hills /=0) THEN
          WRITE ( iw, '( A)' )    ' METADYN|  -  The last hill in the restart file will be re-evaluated'
          WRITE ( iw, '( A)' )    ' METADYN|  -  in order to avoid double counting.'
          WRITE ( iw, '( A,T71,I10)' )    ' METADYN| Actual number of previously spawned hills',&
               meta_env%hills_env%n_hills-1
       END IF
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
          WRITE ( iw, '( A,T71,I10)' )      ' METAVARS| Subsys',meta_env%metavar(i)%isubsys
          IF (meta_env%extended_lagrange) THEN
             WRITE ( iw, '( A,T71,F10.6)' )    ' METAVARS| Lambda Parameter',meta_env%metavar(i)%lambda
             WRITE ( iw, '( A,T66,F15.6)' )    ' METAVARS| Collective Variable Mass',meta_env%metavar(i)%mass
          END IF
          WRITE ( iw, '( A,T71,F10.6)' )    ' METAVARS| Scaling factor',meta_env%metavar(i)%delta_s
          IF (meta_env%metavar(i)%do_wall_p) &
               WRITE ( iw, '( A,T61,F12.6)' )    ' METAVARS| Upper Wall Height',meta_env%metavar(i)%wall_p
          IF (meta_env%metavar(i)%do_wall_m) &
               WRITE ( iw, '( A,T61,F12.6)' )    ' METAVARS| Lower Wall Height',meta_env%metavar(i)%wall_m
          WRITE ( iw, '( A )')'          ----------------------------------------------------------------------'
       ENDDO
    END IF
    CALL cp_print_key_finished_output(iw,logger,metadyn_section,&
         "PRINT%PROGRAM_RUN_INFO", error=error)

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metavar_read', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n_rep
    LOGICAL                                  :: failure

    failure=.FALSE.
    CALL section_vals_val_get(metavar_section,"SUBSYS",i_rep_section=icol,i_val=metavar%isubsys,error=error)
    CALL section_vals_val_get(metavar_section,"COLVAR",i_rep_section=icol,i_val=metavar%icolvar,error=error)
    !
    CALL section_vals_val_get(metavar_section,"WALL_PLUS",i_rep_section=icol,n_rep_val=n_rep,error=error)
    IF (n_rep /= 0) THEN
       metavar%do_wall_p=.TRUE.
       CALL section_vals_val_get(metavar_section,"WALL_PLUS",i_rep_section=icol,r_val=metavar%wall_p,error=error)
    END IF
    CALL section_vals_val_get(metavar_section,"WALL_MINUS",i_rep_section=icol,n_rep_val=n_rep,error=error)
    IF (n_rep /= 0) THEN
       metavar%do_wall_m=.TRUE.
       CALL section_vals_val_get(metavar_section,"WALL_MINUS",i_rep_section=icol,r_val=metavar%wall_m,error=error)
    END IF
    CALL section_vals_val_get(metavar_section,"SCALE",i_rep_section=icol,r_val=metavar%delta_s,error=error)
    !
    ! Check Input Inconsistencies..
    !
    IF(extended_lagrange)THEN
       CALL section_vals_val_get(metavar_section,"MASS",i_rep_section=icol,n_rep_val=n_rep, error=error)
       CPPostcondition(n_rep /= 0,cp_failure_level,routinep,error,failure)
       CALL section_vals_val_get(metavar_section,"MASS",i_rep_section=icol,r_val=metavar%mass,error=error)
       !
       CALL section_vals_val_get(metavar_section,"LAMBDA",i_rep_section=icol, n_rep_val=n_rep, error=error)
       CPPostcondition(n_rep /= 0,cp_failure_level,routinep,error,failure)
       CALL section_vals_val_get(metavar_section,"LAMBDA",i_rep_section=icol,r_val=metavar%lambda,error=error)
    ENDIF

  END SUBROUTINE metavar_read

!!****f* metadynamics/metadyn_add_forces *
!!
!!   NAME
!!     metadyn_add_forces
!!
!!   FUNCTION
!!     add forces to the subsys due to the metadynamics run
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
  SUBROUTINE metadyn_add_forces(meta_env,force_env,error)
!FM not safe with multiple non overlapping para_env for the subsys (sounds likely)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'metadyn_add_forces', &
      routineP = moduleN//':'//routineN
    REAL(kind=dp), PARAMETER                 :: FACTEM = 315795.260 

    INTEGER                                  :: i, i_c, icolvar, ii, isubsys, &
                                                iw
    LOGICAL                                  :: explicit, failure
    REAL(kind=dp)                            :: dt, ekin_w, fac_t, fft, norm, &
                                                rval, scal, ss0_test, tol_ekin
    TYPE(colvar_type), DIMENSION(:), POINTER :: colvar
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(metavar_type), POINTER              :: cv
    TYPE(particle_list_type), POINTER        :: particles
    TYPE(section_vals_type), POINTER         :: ss0_section, vvp_section

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    NULLIFY(colvar,subsys,cv,ss0_section, vvp_section)
    failure=.FALSE.
    CPPrecondition(ASSOCIATED(meta_env),cp_failure_level,routineP,error,failure)
    CALL force_env_get(force_env, subsys=subsys, error=error)

    dt = meta_env%dt
    meta_env%n_steps=meta_env%n_steps+1

    ! compute ss and the derivative of ss with respect to the atomic positions
    DO i_c=1,meta_env%n_colvar
       cv => meta_env%metavar(i_c)
       isubsys = cv%isubsys
       icolvar = cv%icolvar
       CALL colvar_eval_f(isubsys,icolvar,force_env,error=error)
       cv%ss = subsys(isubsys)%subsys%colvar(icolvar)%ss
       !
       ! Restart for Extended Lagrangian Metadynamics
       !
       IF (meta_env%restart) THEN
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
       IF(.NOT.meta_env%extended_lagrange) cv%ss0 =  cv%ss
    ENDDO
    meta_env%restart = .FALSE.

    ! history dependent forces (evaluated at s0)
    IF(meta_env%do_hills)THEN
       CALL  hills(meta_env,error)
    ENDIF

    IF(.NOT. meta_env%extended_lagrange)THEN
       meta_env%ekin_s    = 0.0_dp
       meta_env%cv_energy = 0.0_dp
       DO i_c = 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          cv%cv_energy = 0.0_dp
          cv%ff_s      = 0.0_dp
          isubsys=cv%isubsys
          icolvar=cv%icolvar
          NULLIFY(particles)
          CALL cp_subsys_get(subsys(isubsys)%subsys, colvar=colvar, &
               particles=particles,error=error)
          DO ii=1,colvar(icolvar)%n_atom_s
             i=colvar(icolvar)%i_atom(ii)
             particles%els(i)%f=particles%els(i)%f+cv%ff_hills*colvar(icolvar)%dsdr(:,ii)
          ENDDO
       ENDDO
    ELSE
       meta_env%ekin_s = 0.0_dp
       DO i_c= 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          cv%cv_energy = 0.5_dp*cv%lambda*(cv%ss-cv%ss0)**2.0_dp
          cv%ff_s      = cv%lambda*(cv%ss-cv%ss0)
          isubsys      = cv%isubsys
          icolvar      = cv%icolvar
          ! forces on the atoms
          NULLIFY(particles)
          CALL cp_subsys_get(subsys(isubsys)%subsys, colvar=colvar, &
               particles=particles,error=error)
          DO ii=1,colvar(icolvar)%n_atom_s
             i=colvar(icolvar)%i_atom(ii)
             particles%els(i)%f=particles%els(i)%f- cv%ff_s*colvar(icolvar)%dsdr(:,ii)
          ENDDO
          !  velocity verlet on the s0
          fft    = cv%ff_s+cv%ff_hills
          cv%vv  = cv%vvp+0.5E0_dp*dt*fft/cv%mass
          cv%vvp = cv%vv+0.5E0_dp*dt*fft/cv%mass
          meta_env%ekin_s = meta_env%ekin_s + 0.5E0_dp*cv%mass*cv%vvp**2
       ENDDO
       !  velocity rescaling on the s0
       IF(meta_env%tempcontrol)THEN
          ekin_w   = 0.5_dp*meta_env%temp_wanted/FACTEM*REAL(meta_env%n_colvar,KIND=dp)
          tol_ekin = 0.5_dp*meta_env%toll_temp/FACTEM*REAL(meta_env%n_colvar,KIND=dp)
          IF(ABS(ekin_w-meta_env%ekin_s)> tol_ekin)THEN
             fac_t=SQRT(ekin_w/MAX(meta_env%ekin_s,1.0E-8_dp))
             DO i_c = 1, meta_env%n_colvar
                cv     => meta_env%metavar(i_c)
                cv%vvp = cv%vvp*fac_t
             ENDDO
             meta_env%ekin_s = ekin_w
          ENDIF
       ENDIF
       ! reflecting boundary conditions on the s0  
       DO i_c = 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          ss0_test=cv%ss0+dt*cv%vvp
          IF((cv%do_wall_p .AND. ss0_test > cv%wall_p) .OR. &
               (cv%do_wall_m .AND. ss0_test < cv%wall_m))  THEN    
             cv%vvp = -cv%vvp
          ENDIF
       ENDDO
       ! update of ss0
       meta_env%cv_energy=0.0_dp
       DO i_c = 1, meta_env%n_colvar
          cv => meta_env%metavar(i_c)
          cv%ss0 = cv%ss0 + dt*cv%vvp
          meta_env%cv_energy = meta_env%cv_energy + cv%cv_energy
       ENDDO
    ENDIF
    ! reflecting boundary conditions on the ss  
    DO i_c = 1, meta_env%n_colvar
       cv => meta_env%metavar(i_c)
       IF((cv%do_wall_p .AND. cv%ss  > cv%wall_p) .OR. &
            (cv%do_wall_m .AND. cv%ss  < cv%wall_m))  THEN    
          NULLIFY(particles)
          isubsys=cv%isubsys
          icolvar=cv%icolvar
          CALL cp_subsys_get(subsys(isubsys)%subsys, colvar=colvar, &
               particles=particles,error=error)
          scal=0.0_dp
          norm=0.0_dp
          DO ii=1,colvar(icolvar)%n_atom_s
             i=colvar(icolvar)%i_atom(ii)
             scal=scal+particles%els(i)%v(1)*colvar(icolvar)%dsdr(1,ii)
             scal=scal+particles%els(i)%v(2)*colvar(icolvar)%dsdr(2,ii)
             scal=scal+particles%els(i)%v(3)*colvar(icolvar)%dsdr(3,ii)
             norm=norm+colvar(icolvar)%dsdr(1,ii)**2
             norm=norm+colvar(icolvar)%dsdr(2,ii)**2
             norm=norm+colvar(icolvar)%dsdr(3,ii)**2
          ENDDO
          DO ii=1,colvar(icolvar)%n_atom_s
             i=colvar(icolvar)%i_atom(ii)
             particles%els(i)%v(:)=particles%els(i)%v(:)-2.0_dp*colvar(icolvar)%dsdr(:,ii)*scal/norm
          ENDDO
       ENDIF
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

  END SUBROUTINE metadyn_add_forces

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
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'hills', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ih, iw, j, n_colvar, &
                                                n_step, ndum, stat
    INTEGER, POINTER                         :: n_hills
    LOGICAL                                  :: explicit, failure
    REAL(dp)                                 :: dfunc, dp2, rval
    REAL(dp), DIMENSION(:, :), POINTER       :: delta_s_history, ss_history, &
                                                tmp
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rvals, tnp, ww_history
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(hills_env_type), POINTER            :: hills_env
    REAL(dp)                                 :: ddp(meta_env%n_colvar)
    TYPE(metavar_type), DIMENSION(:), &
      POINTER                                :: colvars
    TYPE(section_vals_type), POINTER         :: hills_history

    failure=.FALSE.
    hills_env       => meta_env%hills_env
    logger          => cp_error_get_logger(error)
    ss_history      => hills_env%ss_history
    delta_s_history => hills_env%delta_s_history
    ww_history      => hills_env%ww_history
    n_hills         => hills_env%n_hills
    colvars         => meta_env%metavar

    n_colvar=meta_env%n_colvar
    n_step=meta_env%n_steps
    NULLIFY(tmp, tnp, hills_history, rvals)
    !
    ! Set-up restart if any
    !
    IF (meta_env%hills_env%restart) THEN
       meta_env%hills_env%restart = .FALSE.
       hills_history => section_vals_get_subs_vals(meta_env%metadyn_section,"SPAWNED_HILLS_POS",error=error)
       CALL section_vals_get(hills_history, explicit=explicit, error=error)
       IF (explicit) THEN
          n_hills = n_hills - 1
          CALL section_vals_val_get( hills_history,"_DEFAULT_KEYWORD_",n_rep_val=ndum, error=error)
          !
          ! ss_history
          !
          DEALLOCATE(ss_history,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          DEALLOCATE(delta_s_history,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          DEALLOCATE(ww_history,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          !
          CPPostcondition(n_hills==ndum-1,cp_failure_level,routinep,error,failure)
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
             !
             ! delta_s_history
             !
             CALL section_vals_val_get( hills_history,"_DEFAULT_KEYWORD_",n_rep_val=ndum, error=error)
             CPPostcondition(n_hills==ndum-1,cp_failure_level,routinep,error,failure)
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
             !
             ! ww_history
             !
             CALL section_vals_val_get( hills_history,"_DEFAULT_KEYWORD_",&
                  n_rep_val=ndum, error=error)
             CPPostcondition(n_hills==ndum-1,cp_failure_level,routinep,error,failure)
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
    !
    ! Proceed with normal calculation
    !
    IF(MODULO(n_step,hills_env%nt_hills)==0)THEN
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
       !
       ! Now add the hill
       !
       DO i=1,n_colvar
          ss_history(i,n_hills)      = colvars(i)%ss0
          delta_s_history(i,n_hills) = colvars(i)%delta_s
          ww_history(n_hills)        = hills_env%ww
       ENDDO

       iw = cp_print_key_unit_nr(logger,meta_env%metadyn_section,&
            "PRINT%PROGRAM_RUN_INFO",extension=".metadynLog",error=error)
       IF (iw>0) THEN
          WRITE(iw,'(/,1X,"METADYN|",A,I7,A,/)')' Hills number ',n_hills,' added.'
       END IF
       CALL cp_print_key_finished_output(iw,logger,meta_env%metadyn_section,&
            "PRINT%PROGRAM_RUN_INFO", error=error)
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
    ENDIF
    !
    ! Forces from the hills
    !
    DO ih=1,n_colvar
       colvars(ih)%ff_hills=0.0_dp
    ENDDO
    DO ih=1,n_hills
       dp2=0._dp
       DO i=1,n_colvar
          ddp(i)=(colvars(i)%ss0-ss_history(i,ih))/delta_s_history(i,ih)
          dp2=dp2+ddp(i)**2
       ENDDO
       dfunc=-ww_history(ih)*EXP(-0.5*dp2)
       DO i=1,n_colvar
          colvars(i)%ff_hills=colvars(i)%ff_hills-ddp(i)/delta_s_history(i,ih)*dfunc
       ENDDO
    ENDDO
  END SUBROUTINE hills

END MODULE metadynamics
