#include "cp_prep_globals.h"

MODULE metadynamics
  USE colvar,                          ONLY: colvar_eval_f,&
                                             colvar_read
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE force_env_types,                 ONLY: force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE machine,                         ONLY: m_flush
  USE metadynamics_types,              ONLY: colvar_init,&
                                             colvar_type,&
                                             hills_env_type,&
                                             meta_env_type
  USE qs_parser,                       ONLY: finish_parser,&
                                             read_object,&
                                             search_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE convert_units, only: convert_to_cp2k_units
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
SUBROUTINE metadyn_create(meta_env,n_colvar,dt,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    INTEGER, INTENT(in)                      :: n_colvar
    REAL(dp), INTENT(in)                     :: dt
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routinen = 'metadyn_create', &
      routinep = modulen//':'//routinen

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(.NOT.ASSOCIATED(meta_env),cp_failure_level,routinep,error,failure)
  ALLOCATE(meta_env,stat=stat)
  IF (.NOT. failure) THEN
     meta_env%ref_count=1
     last_meta_env_id=last_meta_env_id+1
     meta_env%id_nr=last_meta_env_id
     meta_env%n_colvar=n_colvar
    
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    meta_env%ekin_s=HUGE(0.0_dp)
    meta_env%cv_energy=-HUGE(0.0_dp)
    meta_env%n_steps=0
    meta_env%dt = dt

    ! hills_env
    ALLOCATE(meta_env%hills_env%ss_history(n_colvar,0), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)   
    meta_env%hills_env%n_hills=0
    meta_env%hills_env%restart_hills=.FALSE.

    ! colvar
    ALLOCATE(meta_env%colvars(n_colvar), stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
    DO i=1,n_colvar
       CALL colvar_init(meta_env%colvars(i),index=i,error=error)
    END DO

 END IF
END SUBROUTINE metadyn_create
!***************************************************************************

SUBROUTINE metadyn_read(meta_env,globenv,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(global_environment_type)            :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routinen = 'read_metadyn', &
      routinep = modulen//':'//routinen

    CHARACTER(len=3)                         :: test_result
    CHARACTER(len=default_string_length)     :: keyword, time_units
    INTEGER                                  :: n_colvar
    LOGICAL                                  :: do_hills, failure, &
                                                section_found
    REAL(kind=dp)                            :: dt, ww

  CPPreconditionNoFail(.NOT.ASSOCIATED(meta_env),cp_failure_level,routineP,error)
  CALL start_parser(file_name=globenv%input_file_name,&
       globenv=globenv,&
       start_section_label="metadynamics",&
       section_found=section_found)
  IF (.NOT.section_found) THEN
     CALL finish_parser()
     RETURN
  END IF
  DO
     test_result = test_object(newline=.TRUE.)
     IF (test_result == "EOS") THEN
        EXIT
     ELSE IF (test_result == "EOF") THEN
        CALL stop_parser(routinen,modulen,__LINE__,"EOF")
     END IF

     CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
     do_hills=.FALSE.
     SELECT CASE(keyword)
     CASE("NCOLVAR","N_COLVAR")
        CALL read_object(n_colvar)
     END SELECT
  END DO
  CALL finish_parser()

  CALL start_parser(file_name=globenv%input_file_name,&
       globenv=globenv,&
       start_section_label="MD",&
       section_found=section_found)
  CPPreconditionNoFail(section_found,cp_failure_level,routineP,error)
  DO
     test_result = test_object(newline=.TRUE.)
     IF (test_result == "EOS") THEN
        EXIT
     ELSE IF (test_result == "EOF") THEN
        CALL stop_parser(routinen,modulen,__LINE__,"EOF")
     END IF

     CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
     do_hills=.FALSE.
     SELECT CASE(keyword)
     CASE("TIMESTEP","TIMESTEP:")
        time_units='FS'
        CALL read_object(dt)
        IF (test_object(newline=.FALSE.)=='STR') THEN
           CALL read_object(time_units,lower_to_upper=.TRUE.,newline=.false.)
        end if
        CALL convert_to_cp2k_units (time_units,time=dt)
     END SELECT
  END DO
  CALL finish_parser()

  CALL metadyn_create(meta_env,n_colvar=n_colvar,&
       dt=dt,error=error)
  
  
  CALL start_parser(file_name=globenv%input_file_name,&
       globenv=globenv,&
       start_section_label="METADYNAMICS",&
       section_found=section_found)
  DO
     test_result = test_object(newline=.TRUE.)
     IF (test_result == "EOS") THEN
        EXIT
     ELSE IF (test_result == "EOF") THEN
        CALL stop_parser(routinen,modulen,__LINE__,"EOF")
     END IF

     CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
     do_hills=.FALSE.
     SELECT CASE(keyword)
     CASE("NCOLVAR","N_COLVAR")
     CASE("DO_HILLS")
        meta_env%do_hills=.TRUE.
     CASE("WW")
        CALL read_object(meta_env%hills_env%ww)
     CASE("NT_HILLS")
        CALL read_object(meta_env%hills_env%nt_hills)
     CASE("RESTART_HILLS")
     END SELECT
  END DO
  CALL finish_parser()

  CALL start_parser(file_name=globenv%input_file_name,&
       globenv=globenv,&
       start_section_label="",&
       section_found=section_found)
  DO i=1,n_colvar
     CALL colvar_read(meta_env%colvars(i),error=error)
  END DO
  CALL search_object("&COLVAR",ignore_case=.TRUE.,&
       found=section_found)
  CALL cp_assert(.NOT.section_found,cp_failure_level,cp_assertion_failed,&
       routinep,"too many colvar sections",error,failure)
  CALL finish_parser()

END SUBROUTINE metadyn_read

  SUBROUTINE metadyn_add_forces(meta_env,force_env,error)
    TYPE(meta_env_type), POINTER             :: meta_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routinen = 'metadyn_add_forces', &
      routinep = modulen//':'//routinen

    INTEGER                                  :: i_c
    LOGICAL                                  :: failure
    REAL(kind=dp)                            :: dt, fft
    TYPE(colvar_type), DIMENSION(:), POINTER :: colvars

    failure=.FALSE.
    CPPrecondition(ASSOCIATED(meta_env),cp_failure_level,routineP,error,failure)
    colvars => meta_env%colvars
    dt = meta_env%dt
    meta_env%n_steps=meta_env%n_steps+1

    IF(meta_env%n_steps==1)THEN
       IF(meta_env%hills_env%restart_hills)THEN
          OPEN(1012,file='COLVAR',status='OLD',position='APPEND')
       ELSE
          OPEN(1012,file='COLVAR',status='REPLACE')
       END IF
    ENDIF
    DO i_c=1,meta_env%n_colvar
       CALL colvar_eval_f(meta_env%colvars(i_c),force_env,meta_env%n_steps, &
                          error=error)
    ENDDO

    IF(meta_env%do_hills)THEN
       CALL  hills(meta_env,error)
    ENDIF
    DO i_c=1,meta_env%n_colvar
       fft=colvars(i_c)%ff_s+colvars(i_c)%ff_hills
       colvars(i_c)%vv=colvars(i_c)%vvp+0.5d0*dt*fft/colvars(i_c)%mass
       colvars(i_c)%vvp=colvars(i_c)%vv+0.5d0*dt*fft/colvars(i_c)%mass
    ENDDO
    meta_env%ekin_s=0.0_dp
    meta_env%cv_energy=0.0_dp
    DO i_c=1,meta_env%n_colvar
       colvars(i_c)%ss0=colvars(i_c)%ss0+dt*colvars(i_c)%vvp
       meta_env%ekin_s=meta_env%ekin_s + &
                       0.5d0*colvars(i_c)%mass*colvars(i_c)%vvp**2
       meta_env%cv_energy=meta_env%cv_energy + &
                       colvars(i_c)%cv_energy
    ENDDO
    WRITE(1012,'(f12.1,40f13.5)')dt*meta_env%n_steps, &
                   (colvars(i)%ss0,i=1,meta_env%n_colvar), &
                   (colvars(i)%ss,i=1,meta_env%n_colvar), &
                   (colvars(i)%ff_s,i=1,meta_env%n_colvar), &
                   (colvars(i)%ff_hills,i=1,meta_env%n_colvar)
    CALL m_flush(1012)
  END SUBROUTINE metadyn_add_forces

  SUBROUTINE hills(meta_env,error)
    TYPE(meta_env_type), INTENT(inout), &
      TARGET                                 :: meta_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routinen = 'hills', &
      routinep = modulen//':'//routinen

    INTEGER                                  :: i, ih, n_colvar, n_step, stat
    INTEGER, POINTER                         :: n_hills
    LOGICAL                                  :: failure
    REAL(dp)                                 :: dfunc, dp2
    REAL(dp), DIMENSION(:, :), POINTER       :: ss_history, tmp
    TYPE(colvar_type), DIMENSION(:), POINTER :: colvars
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(hills_env_type), POINTER            :: hills_env
    REAL(dp)                                 :: ddp(meta_env%n_colvar)

    failure=.FALSE.
    hills_env => meta_env%hills_env
    logger => cp_error_get_logger(error)
    ss_history => hills_env%ss_history
    colvars => meta_env%colvars
    n_hills => hills_env%n_hills
    n_colvar=meta_env%n_colvar
    n_step=meta_env%n_steps
    NULLIFY(tmp)
    IF (n_step==1) THEN
       IF(hills_env%restart_hills)THEN
          OPEN(667,file='HILLS',status='old')
          DO i=1,100000
             READ(667,*,END=144)
          ENDDO
144       CONTINUE
          n_hills=i-1
          ALLOCATE(ss_history(n_colvar,n_hills), stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routinep,error,failure)
          REWIND(667)
          DO i=1,n_hills
             READ(667,*)ih,(ss_history(ih,i),ih=1,n_colvar)
          ENDDO
          CLOSE(667)
       ELSE
          OPEN(667,file='HILLS',status='replace')
          CLOSE(667)
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
       WRITE(cp_logger_get_default_unit_nr(logger,local=.TRUE.),*)' hills number ',n_hills,' added'
       OPEN(667,file='HILLS',status='unknown',position='append')
       WRITE(667,'(i10,30f13.5)')&
            n_step,(hills_env%ss_history(ih,n_hills),ih=1,n_colvar), &
                   (colvars(ih)%delta_s,ih=1,n_colvar),hills_env%ww
       CLOSE(667)
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
