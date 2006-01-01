
MODULE metadynamics
  USE colvar,                          ONLY: colvar_eval_f,&
                                             colvar_read
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE cp_files,                        ONLY: open_file
  USE cp_para_env,                     ONLY: cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type,cp_subsys_get
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE message_passing,                 ONLY: mp_bcast
  USE metadynamics_types,              ONLY: colvar_init,&
                                             colvar_type,&
                                             hills_env_type,&
                                             meta_env_type
  USE particle_list_types,             ONLY: particle_list_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_get
  USE metadynamics_types,              ONLY: colvar_type,&
                                             coord_colvar_id,&
                                             dist_colvar_id,&
                                             no_colvar_id,&
                                             torsion_colvar_id
#include "cp_common_uses.h"
  
  IMPLICIT NONE
  PRIVATE 

  CHARACTER(len=*), PARAMETER, PRIVATE :: modulen='metadynamics'

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
!!
!!*** **********************************************************************
SUBROUTINE metadyn_create(meta_env,n_colvar,dt,para_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    INTEGER, INTENT(in)                      :: n_colvar
    REAL(dp), INTENT(in)                     :: dt
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routinen = 'metadyn_create', &
      routinep = modulen//':'//routinen

    INTEGER                                  :: stat,i
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(.NOT.ASSOCIATED(meta_env),cp_failure_level,routinep,error,failure)
  ALLOCATE(meta_env,stat=stat)
  IF (.NOT. failure) THEN
     meta_env%ref_count=1
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
    meta_env%colvar_file_unit_nr=-1
    meta_env%colvar_file_name="COLVAR"
    meta_env%tempcontrol=.FALSE.

    ! hills_env
    ALLOCATE(meta_env%hills_env%ss_history(n_colvar,0), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)   
    meta_env%hills_env%n_hills=0
    meta_env%hills_env%restart_hills=.FALSE.
    meta_env%hills_env%hills_file_name='HILLS'

    ! colvar
    ALLOCATE(meta_env%colvars(n_colvar), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    DO i=1,n_colvar
       CALL colvar_init(meta_env%colvars(i),index=i,error=error)
    END DO

 END IF
END SUBROUTINE metadyn_create
!***************************************************************************

SUBROUTINE metadyn_read(meta_env,globenv,metadyn_section,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(global_environment_type)            :: globenv
    TYPE(section_vals_type), POINTER, OPTIONAL :: metadyn_section    
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routinen = 'read_metadyn', &
      routinep = modulen//':'//routinen

    CHARACTER(len=3)                         :: test_result
    CHARACTER(len=default_string_length)     :: keyword, time_units
  INTEGER                                  :: n_colvar, n_colvar_inp, i, n_rep, iw, kk
  LOGICAL                                  :: failure, explicit
  REAL(kind=dp)                            :: dt
  TYPE(section_vals_type), POINTER         :: md_section, colvar_section


  CALL section_vals_get(metadyn_section,explicit=explicit,error=error)
  IF (.NOT.explicit) RETURN
  CPPreconditionNoFail(.NOT.ASSOCIATED(meta_env),cp_failure_level,routineP,error)  

  md_section => section_vals_get_subs_vals(globenv%input_file,"MOTION%MD",error=error)
  CALL section_vals_val_get(md_section,"TIMESTEP",r_val=dt,error=error)

  CALL section_vals_val_get(metadyn_section,"NCOLVAR",i_val=n_colvar,error=error)
  CALL metadyn_create(meta_env,n_colvar=n_colvar,&
       dt=dt,para_env=globenv%para_env,error=error)

  CALL section_vals_val_get(metadyn_section,"DO_HILLS",l_val=meta_env%do_hills,error=error)
  CALL section_vals_val_get(metadyn_section,"LAGRANGE",l_val=meta_env%extended_lagrange,error=error)
  CALL section_vals_val_get(metadyn_section,"WW",r_val=meta_env%hills_env%ww,error=error)
  CALL section_vals_val_get(metadyn_section,"NT_HILLS",i_val=meta_env%hills_env%nt_hills,error=error)
  CALL section_vals_val_get(metadyn_section,"RESTART_HILLS",l_val=meta_env%hills_env%restart_hills,error=error)
  CALL section_vals_val_get(metadyn_section,"TEMPERATURE",n_rep_val=n_rep,error=error)
  IF ( n_rep /= 0 ) THEN
     meta_env%tempcontrol=.TRUE.
     CALL section_vals_val_get(metadyn_section,"TEMPERATURE",r_val=meta_env%temp_wanted,error=error)
     meta_env%toll_temp=meta_env%temp_wanted/SQRT(REAL(meta_env%n_colvar,dp))
  END IF

  IF (globenv%project_name/="") THEN
     meta_env%colvar_file_name='COLVAR-'//globenv%project_name
     meta_env%hills_env%hills_file_name='HILLS-'//globenv%project_name
  ELSE
     meta_env%colvar_file_name='COLVAR'
     meta_env%hills_env%hills_file_name="HILLS"
  END IF
  colvar_section => section_vals_get_subs_vals(metadyn_section,"COLVAR",error=error)
  CALL section_vals_get(colvar_section,n_repetition=n_colvar_inp,error=error)
  CPPreconditionNoFail((n_colvar == n_colvar_inp),cp_failure_level,routineP,error)
  DO i= 1, n_colvar
     CALL colvar_read(meta_env%colvars(i),meta_env%extended_lagrange,i,colvar_section,error=error)
  END DO

  IF (globenv % ionode .AND. globenv % print_level>=0) THEN
     iw = globenv % scr
     WRITE ( iw, '( /A )' )           ' METADYN| Meta Dynamics Protocol '
     WRITE ( iw, '( A,T71,I10)' )    ' METADYN| Number of interval time steps to spawn hills',meta_env%hills_env%nt_hills
     IF (meta_env%extended_lagrange) THEN
        WRITE ( iw, '( A )' )        ' METADYN| Extended Lagrangian Scheme '
        IF (meta_env%tempcontrol) WRITE ( iw, '( A,T71,F10.2)' )&
                                     ' METADYN| Collective Variables Temperature control',meta_env%toll_temp
     END IF
     IF (meta_env%hills_env%restart_hills) &
        WRITE ( iw, '( A )' )        ' METADYN| Restarting HILLS '
     IF (meta_env%do_hills) &
          WRITE ( iw, '( A )' )      ' METADYN| Spanwing the Hills '
     WRITE ( iw, '( A,T71,F10.2)' )  ' METADYN| Width of the Spawned Gaussian',meta_env%hills_env%ww
     WRITE ( iw, '( A,T71,I10)' )    ' METADYN| Number of collective variables',meta_env%n_colvar
     
     
     DO i=1,n_colvar
        WRITE ( iw, '( A )')'          ----------------------------------------------------------------------'
        WRITE ( iw, '( A,T71,I10)' )      ' COLVARS| Collective Variable Number',i
        WRITE ( iw, '( A,T71,I10)' )      ' COLVARS| Type of collective variable',meta_env%colvars(i)%type_id
        WRITE ( iw, '( A,T71,I10)' )      ' COLVARS| Subsys',meta_env%colvars(i)%i_subsys
        WRITE ( iw, '( A,T71,F10.6)' )    ' COLVARS| Lambda Parameter',meta_env%colvars(i)%lambda
        WRITE ( iw, '( A,T66,F15.6)' )    ' COLVARS| Collective Variable Mass',meta_env%colvars(i)%mass
        WRITE ( iw, '( A,T71,F10.6)' )    ' COLVARS| Scaling factor',meta_env%colvars(i)%delta_s
        IF (meta_env%colvars(i)%do_wall_p) &
        WRITE ( iw, '( A,T61,F12.6)' )    ' COLVARS| Upper Wall Height',meta_env%colvars(i)%wall_p
        IF (meta_env%colvars(i)%do_wall_m) &
        WRITE ( iw, '( A,T61,F12.6)' )    ' COLVARS| Lower Wall Height',meta_env%colvars(i)%wall_m
        SELECT CASE(meta_env%colvars(i)%type_id)
        CASE(torsion_colvar_id)
           WRITE ( iw, '( A,T61,4I5)' )   ' COLVARS| TORSION       >>> ATOMS:',meta_env%colvars(i)%torsion_param%i_at_tors
        CASE(dist_colvar_id)
           WRITE ( iw, '( A,T71,2I5)' )   ' COLVARS| BOND          >>> ATOMS:',meta_env%colvars(i)% dist_param%i_at,&                
                meta_env%colvars(i)% dist_param%j_at
        CASE(coord_colvar_id)
           WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| COORDINATION  >>> FROM ATOM:',meta_env%colvars(i)%coord_param%i_at_from(kk),&
                kk=1,SIZE(meta_env%colvars(i)%coord_param%i_at_from))
           WRITE ( iw, '( A,T71,I10)' )  (' COLVARS| COORDINATION  >>>   TO ATOM:',meta_env%colvars(i)%coord_param%i_at_to(kk),&
                kk=1,SIZE(meta_env%colvars(i)%coord_param%i_at_to))
           WRITE ( iw, '( A,T71,F10.5)' ) ' COLVARS| R0',meta_env%colvars(i)%coord_param%r_0
           WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| NN',meta_env%colvars(i)%coord_param%nncrd
           WRITE ( iw, '( A,T71,I10)' )   ' COLVARS| ND',meta_env%colvars(i)%coord_param%ndcrd
           WRITE ( iw, '( A )')'          ----------------------------------------------------------------------'
        END SELECT
     END DO
  END IF

END SUBROUTINE metadyn_read


SUBROUTINE metadyn_add_forces(meta_env,force_env,error)
  !FM not safe with multiple non overlapping para_env for the subsys
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routinen = 'metadyn_add_forces', &
      routinep = modulen//':'//routinen

    INTEGER                                  :: i_c, np, stat, ii, i
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: dt, fft, ekin_w,tol_ekin,fac_t, ss0_test
    TYPE(colvar_type), DIMENSION(:), POINTER :: colvars
    TYPE(colvar_type), POINTER               ::  CV
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    REAL(kind=dp) , PARAMETER                :: FACTEM=315795.260
    REAL(kind=dp)                            :: scal,norm
    TYPE(particle_list_type), POINTER        :: particles
  NULLIFY(colvars,subsys)
  failure=.FALSE.
  CPPrecondition(ASSOCIATED(meta_env),cp_failure_level,routineP,error,failure)
  CALL force_env_get(force_env, subsys=subsys, error=error)

  colvars => meta_env%colvars
  dt = meta_env%dt
  meta_env%n_steps=meta_env%n_steps+1

  ! compute ss and the derivative of ss with respect to the atomic positions
  DO i_c=1,meta_env%n_colvar
     CV=>colvars(i_c)
     CALL colvar_eval_f(CV,force_env,error=error)
     IF(.not.meta_env%extended_lagrange) CV%ss0=CV%ss
     IF(meta_env%n_steps==1)CV%ss0=CV%ss
  ENDDO
       

  ! history dependent forces (evaluated at s0)
  IF(meta_env%do_hills)THEN
     CALL  hills(meta_env,error)
  ENDIF

  IF(.not. meta_env%extended_lagrange)THEN
     meta_env%ekin_s=0.0_dp
     meta_env%cv_energy=0.0_dp
     DO i_c=1,meta_env%n_colvar
        CV=>colvars(i_c)
        CV%cv_energy=0.0_dp
        CV%ff_s=0.0_dp
        NULLIFY(particles)
        CALL cp_subsys_get(subsys(CV%i_subsys)%subsys, particles=particles,error=error)
        do ii=1,CV%n_atom_s
          i=CV%i_atom(ii)
          particles%els(i)%f=particles%els(i)%f+ CV%ff_hills*CV%dsdr(:,ii)
        enddo
     ENDDO
  ELSE
     meta_env%ekin_s=0.0_dp
     DO i_c=1,meta_env%n_colvar
        CV=>colvars(i_c)
        CV%cv_energy=0.5_dp*CV%lambda*(CV%ss-CV%ss0)**2
        CV%ff_s=CV%lambda* (CV%ss-CV%ss0)
        ! forces on the atoms
        NULLIFY(particles)
        CALL cp_subsys_get(subsys(CV%i_subsys)%subsys,particles=particles,error=error)
        do ii=1,CV%n_atom_s
          i=CV%i_atom(ii)
          particles%els(i)%f=particles%els(i)%f- CV%ff_s*CV%dsdr(:,ii)
        enddo
        !  velocity verlet on the s0
        fft=CV%ff_s+CV%ff_hills
        CV%vv=CV%vvp+0.5E0_dp*dt*fft/CV%mass
        CV%vvp=CV%vv+0.5E0_dp*dt*fft/CV%mass
        meta_env%ekin_s=meta_env%ekin_s + 0.5E0_dp*CV%mass*CV%vvp**2
     ENDDO
     !  velocity rescaling on the s0
     if(meta_env%tempcontrol)then
       ekin_w=0.5_dp*meta_env%temp_wanted/FACTEM*REAL(meta_env%n_colvar,dp)
       tol_ekin=0.5_dp*meta_env%toll_temp/FACTEM*REAL(meta_env%n_colvar,dp)
       if(abs(ekin_w-meta_env%ekin_s)> tol_ekin)then
         fac_t=sqrt(ekin_w/max(meta_env%ekin_s,1.0E-8_dp))
         DO i_c=1,meta_env%n_colvar
           CV=>colvars(i_c)
           CV%vvp=CV%vvp*fac_t
         ENDDO
         meta_env%ekin_s=ekin_w
       endif
     endif
     ! reflecting boundary conditions on the s0  
     DO i_c=1,meta_env%n_colvar
       CV=>colvars(i_c)
       ss0_test=CV%ss0+dt*CV%vvp
       if((CV%do_wall_p .and. ss0_test > CV%wall_p) .or. &
          (CV%do_wall_m .and. ss0_test < CV%wall_m))  then    
         CV%vvp=-CV%vvp
       endif
     ENDDO
     ! update of ss0
     meta_env%cv_energy=0.0_dp
     DO i_c=1,meta_env%n_colvar
       CV=>colvars(i_c)
       CV%ss0=CV%ss0+dt*CV%vvp
       meta_env%cv_energy=meta_env%cv_energy + CV%cv_energy
     ENDDO
  ENDIF
  ! reflecting boundary conditions on the ss  
  DO i_c=1,meta_env%n_colvar
    CV=>colvars(i_c)
    if((CV%do_wall_p .and. CV%ss  > CV%wall_p) .or. &
       (CV%do_wall_m .and. CV%ss  < CV%wall_m))  then    
      NULLIFY(particles)
      CALL cp_subsys_get(subsys(CV%i_subsys)%subsys, particles=particles,error=error)
      scal=0.0_dp
      norm=0.0_dp
      do ii=1,CV%n_atom_s
        i=CV%i_atom(ii)
        scal=scal+particles%els(i)%v(1)*CV%dsdr(1,ii)
        scal=scal+particles%els(i)%v(2)*CV%dsdr(2,ii)
        scal=scal+particles%els(i)%v(3)*CV%dsdr(3,ii)
        norm=norm+CV%dsdr(1,ii)**2
        norm=norm+CV%dsdr(2,ii)**2
        norm=norm+CV%dsdr(3,ii)**2
      enddo
      do ii=1,CV%n_atom_s
        i=CV%i_atom(ii)
        particles%els(i)%v(:)=particles%els(i)%v(:)-2.0_dp*CV%dsdr(:,ii)*scal/norm
      enddo
    endif
  ENDDO

  ! write COLVAR file
  IF (meta_env%para_env%source==meta_env%para_env%mepos) THEN
    IF(meta_env%n_steps==1)THEN
       IF(meta_env%hills_env%restart_hills)THEN
          CALL open_file(trim(meta_env%colvar_file_name),file_status="OLD",&
             file_action="READ",unit_number=meta_env%colvar_file_unit_nr)
          DO i=1,100000
             READ(meta_env%colvar_file_unit_nr,*,END=144)
          ENDDO
144       CONTINUE
          CLOSE(meta_env%colvar_file_unit_nr)
          meta_env%n_steps=i
          CALL open_file(trim(meta_env%colvar_file_name),file_status='OLD',&
             file_position='APPEND', file_action="WRITE",&
             unit_number=meta_env%colvar_file_unit_nr)
       ELSE
          CALL open_file(trim(meta_env%colvar_file_name),file_status='REPLACE',&
             file_action="WRITE",unit_number=meta_env%colvar_file_unit_nr)
       END IF
    ENDIF
    WRITE(meta_env%colvar_file_unit_nr,'(f12.1,40f13.5)')dt*meta_env%n_steps, &
          (colvars(i)%ss0,i=1,meta_env%n_colvar), &
          (colvars(i)%ss,i=1,meta_env%n_colvar), &
          (colvars(i)%ff_s,i=1,meta_env%n_colvar), &
          (colvars(i)%ff_hills,i=1,meta_env%n_colvar)
     CALL m_flush(meta_env%colvar_file_unit_nr)
  ENDIF
END SUBROUTINE metadyn_add_forces

  SUBROUTINE hills(meta_env,error)
    TYPE(meta_env_type), POINTER  :: meta_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routinen = 'hills', &
      routinep = modulen//':'//routinen

    INTEGER                                  :: i, ih, n_colvar, n_step, &
                                                stat, unit_number
    INTEGER, POINTER                         :: n_hills
    LOGICAL                                  :: failure
    REAL(dp)                                 :: dfunc, dp2
    REAL(dp), DIMENSION(:, :), POINTER       :: ss_history, tmp
    TYPE(colvar_type), DIMENSION(:), POINTER :: colvars
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(hills_env_type), POINTER            :: hills_env
    REAL(dp)                                 :: ddp(meta_env%n_colvar)
    REAL(dp)                                 :: dum
    LOGICAL                                  :: MASTER

    failure=.FALSE.
    hills_env => meta_env%hills_env
    logger => cp_error_get_logger(error)
    ss_history => hills_env%ss_history
    colvars => meta_env%colvars
    n_hills => hills_env%n_hills
    n_colvar=meta_env%n_colvar
    n_step=meta_env%n_steps
    NULLIFY(tmp)
    MASTER=meta_env%para_env%source==meta_env%para_env%mepos
    IF (n_step==1) THEN
       IF(hills_env%restart_hills)THEN
          DEALLOCATE(ss_history,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          IF (MASTER) THEN
            CALL open_file(hills_env%hills_file_name,file_status="old",&
               file_action="READ",unit_number=unit_number)
            DO i=1,100000
               READ(unit_number,*,END=144)
            ENDDO
144         CONTINUE
            REWIND(unit_number)
            n_hills=i-1
            ALLOCATE(ss_history(n_colvar,n_hills), stat=stat)
            CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
            DO i=1,n_hills
               READ(unit_number,*)dum,(ss_history(ih,i),ih=1,n_colvar)
            ENDDO
            CLOSE(unit_number)
          ENDIF
          CALL MP_BCAST(n_hills,meta_env%para_env%source,meta_env%para_env%group)
          IF (.not. MASTER) THEN
            ALLOCATE(ss_history(n_colvar,n_hills), stat=stat)
            CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          ENDIF
          CALL MP_BCAST(ss_history,meta_env%para_env%source,meta_env%para_env%group)
          hills_env%ss_history=>ss_history
       ELSE IF(MASTER) then
          CALL open_file(hills_env%hills_file_name,file_status='replace',&
               file_action="WRITE",unit_number=unit_number)
          CLOSE(unit_number)
       ENDIF
    ENDIF
    IF(MODULO(n_step,hills_env%nt_hills)==0)THEN
       IF(SIZE(ss_history,2)< n_hills+1)THEN
          ALLOCATE(tmp(n_colvar,n_hills+100), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          tmp(:,:n_hills)=ss_history
          tmp(:,n_hills+1:)=0.0_dp
          DEALLOCATE(ss_history,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          hills_env%ss_history => tmp
          ss_history => tmp
          NULLIFY(tmp)
       ENDIF
       n_hills=n_hills +1
       ! now add the hill
       DO i=1,n_colvar
          ss_history(i,n_hills)=colvars(i)%ss0
       ENDDO
       IF(MASTER) then
         WRITE(cp_logger_get_default_unit_nr(logger,local=.TRUE.),*)' hills number ',n_hills,' added'
         CALL open_file(hills_env%hills_file_name,file_status='unknown',&
            file_action="WRITE",file_position='append', unit_number=unit_number)
         WRITE(unit_number,'(f12.1,30f13.5)')      &
                    meta_env%dt*meta_env%n_steps,  &
                   (hills_env%ss_history(ih,n_hills),ih=1,n_colvar), &
                   (colvars(ih)%delta_s,ih=1,n_colvar),hills_env%ww
         CLOSE(unit_number)
       ENDIF
    ENDIF

    ! forces from the hills
    DO ih=1,n_colvar
      colvars(ih)%ff_hills=0.0_dp
    ENDDO
    DO ih=1,n_hills
       dp2=0._dp
       DO i=1,n_colvar
          ddp(i)=(colvars(i)%ss0-ss_history(i,ih))/colvars(i)%delta_s
          dp2=dp2+ddp(i)**2
       ENDDO
       dfunc=-hills_env%ww*EXP(-0.5*dp2)
       DO i=1,n_colvar
          colvars(i)%ff_hills=colvars(i)%ff_hills- &
           ddp(i)/colvars(i)%delta_s*dfunc
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE hills

END MODULE metadynamics
