!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2004  CP2K developers group
!-----------------------------------------------------------------------------!
#include "cp_prep_globals.h"

MODULE farming_methods
  USE farming_types,                   ONLY: farming_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: default_string_length
  USE cp_log_handling,                 ONLY: cp_logger_type, &
                                             cp_failure_level, &
                                             cp_to_string
  USE cp_error_handling,               ONLY: cp_error_type, &
                                             cp_a_l, &
                                             cp_assert, &
                                             cp_error_init, &
                                             cp_error_dealloc_ref, &
                                             cp_assertion_failed 
  USE message_passing,                 ONLY: mp_bcast
  USE qs_parser,                       ONLY: finish_parser,&
                                             read_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object,&
                                             get_unit_number
  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: farming_parse_input

  CHARACTER(LEN=*), PARAMETER, PRIVATE :: moduleN = "farming_runs"

CONTAINS

  SUBROUTINE farming_parse_input(farming_env,globenv,error)
    TYPE(farming_env_type), POINTER          :: farming_env
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'farming_parse_input', &
      routineP = moduleN//':'//routineN
    CHARACTER(LEN=3)                         :: test_result
    INTEGER :: stat,i,iunit
    LOGICAL :: found,failure
    CHARACTER(LEN=default_string_length) :: dum
    CHARACTER(LEN=60)                        :: keyword
    CHARACTER(LEN=3)                         :: text
    failure=.FALSE.

!
!  &FARMING section sets options for the farming 
!    group_size  ! wished number of MPI tasks per group
!    ngroups     ! wished number of groups 
!  &END FARMING 
!
   farming_env%group_size_wish_set=.FALSE.
   farming_env%ngroup_wish_set    =.FALSE.
   CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="FARMING",section_found=found)
   IF (found) THEN 
   DO
        test_result = test_object()
        IF (test_result /= "EOL") THEN
          CALL read_object(keyword)
          CALL stop_parser(routine=routineN,message="INVALID_KEYWORD")
        END IF
        test_result = test_object(newline=.TRUE.)
        IF (test_result == "EOS") THEN
          EXIT
        ELSE IF (test_result == "EOF") THEN
          CALL stop_parser(routine=routineN,message="EOF")
        ELSE
          CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
          SELECT CASE (TRIM(keyword))
          CASE("GROUP_SIZE")
            IF (farming_env%ngroup_wish_set) & 
               CALL stop_parser(routine=routineN,message="specify one of GROUP_SIZE or NGROUPS ")
            CALL read_object(farming_env%group_size_wish)
            farming_env%group_size_wish_set=.TRUE.
          CASE("NGROUPS","NGROUP") 
            IF (farming_env%group_size_wish_set) & 
               CALL stop_parser(routine=routineN,message="specify one of GROUP_SIZE or NGROUPS ")
            CALL read_object(farming_env%ngroup_wish)
            farming_env%ngroup_wish_set=.TRUE.
          !
          ! allows strictly setting the size of the different partitions, up to the user to do the
          ! right thing for the allocated number of CPUs
          !
          CASE ("GROUP_PARTITION")
            IF (.NOT. farming_env%ngroup_wish_set) &
               CALL stop_parser(routine=routineN,message="specify NGROUPS before GROUP_PARTITION")
            ALLOCATE(farming_env%group_partition(0:farming_env%ngroup_wish-1))
            DO I=0,farming_env%ngroup_wish-1
               CALL read_object(farming_env%group_partition(I)) 
            ENDDO
          CASE ("RESTART")
            farming_env%restart=.TRUE.
          CASE ("MAX_STEP","MAX_STEPS")
            CALL read_object(farming_env%max_steps)
          CASE ("CYCLE")
            farming_env%cycle=.TRUE.
          CASE DEFAULT
            CALL stop_parser(routine=routineN,message="INVALID_KEYWORD "&
                 //keyword)
          END SELECT
        ENDIF
   END DO
   ENDIF
   CALL finish_parser()
   ! set a default if nothing given
   IF (.NOT. farming_env%ngroup_wish_set .AND. .NOT. farming_env%group_size_wish_set) THEN
      farming_env%group_size_wish_set=.TRUE.
      farming_env%group_size_wish    = 8 ! groups of 8 CPUs are beautiful (on HPCx ;-)
   ENDIF
!
!   &JOBS section gives the list of jobs
!      running_dir1 input_file_name3
!      running_dir2 input_file_name3
!      running_dir3 input_file_name3
!   &END 
!
    ! count the number of jobs first
    farming_env % njobs = 0
    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="JOBS",section_found=found)

    IF (found) THEN
       DO WHILE (test_object(newline=.TRUE.) /= "EOS")
          farming_env % njobs = farming_env % njobs + 1
          CALL read_object(dum,newline=.TRUE.)
       ENDDO
    ELSE
       failure=.true.
       CALL stop_parser(routine=routineP,message="section &JOBS not found in "//globenv%input_file_name)      
    ENDIF
    CALL finish_parser()

    ! get the job info
    ALLOCATE(farming_env%Job_cwd(farming_env % njobs),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE(farming_env%Job_input(farming_env % njobs),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
        CALL start_parser(file_name=globenv%input_file_name,&
                          globenv=globenv,&
                          start_section_label="JOBS",section_found=found)
        i=0
        DO WHILE (test_object(newline=.TRUE.) /= "EOS")
           i=i+1
           IF (i<=farming_env % njobs) THEN
              CALL read_object(farming_env%Job_cwd(i),newline=.TRUE.)           
              CALL read_object(farming_env%Job_input(i))           
           ELSE
              failure=.true.
              CALL stop_parser(routine=routineP,message="inconsistent reading in the JOBS section")
           ENDIF
        ENDDO
        CALL finish_parser()
    ENDIF


    IF (globenv%ionode .AND. .NOT. failure) THEN
        write(globenv%scr,FMT="(T2,A,T71,I10)") "FARMING| Number of jobs found",farming_env%njobs
        IF (farming_env%ngroup_wish_set) THEN
           write(globenv%scr,FMT="(T2,A,T71,I10)") "FARMING| Ngroup wish:",farming_env%ngroup_wish
           IF (ASSOCIATED(farming_env%group_partition)) THEN
              write(globenv%scr,FMT="(T2,A)",ADVANCE="NO") "FARMING| User partition:"
              DO i=0,SIZE(farming_env%group_partition)-1
                 IF (MODULO(i,4)==0) write(globenv%scr,*)
                 write(globenv%scr,FMT='(I4)',ADVANCE="NO") farming_env%group_partition(i)
              END DO
              write(globenv%scr,*)
              IF (SUM(farming_env%group_partition).NE.globenv%num_pe) THEN
                 write(globenv%scr,FMT="(T2,A,T61,I10,T71,I10)") &
                     "FARMING| WARNING : partition CPUs not equal actual number ", &
                     globenv%num_pe, SUM(farming_env%group_partition)
                 write(globenv%scr,FMT="(T2,A)") "FARMING|          partition data ignored" ! any better idea ??
                 DEALLOCATE(farming_env%group_partition)
              ENDIF
           ENDIF
        ENDIF
        IF (farming_env%group_size_wish_set) THEN
           write(globenv%scr,FMT="(T2,A,T71,I10)") "FARMING| Group size wish:", &
                 farming_env%group_size_wish
        ENDIF
        write(globenv%scr,FMT="(T2,A,T71,I10)") "FARMING| Max steps      :",farming_env%max_steps
        IF (farming_env%cycle) THEN
            text="YES"
        ELSE
            text=" NO"
        ENDIF
        write(globenv%scr,FMT="(T2,A,T78,A3)") "FARMING| Cyclic jobs execution:",text
        IF (farming_env%restart) THEN
            text="YES"
        ELSE
            text=" NO"
        ENDIF
        write(globenv%scr,FMT="(T2,A,T78,A3)") "FARMING| Restarting farm:",text
        farming_env%restart_n = 1
        IF (farming_env%restart) THEN
           iunit=get_unit_number()
           OPEN(UNIT=iunit,FILE="RESTARTFARM",IOSTAT=stat)
           IF (stat==0) THEN
             READ(UNIT=iunit,FMT=*,IOSTAT=stat) farming_env%restart_n
             IF (stat/=0) THEN
                write(globenv%scr,"(T2,A)") &
                     "FARMING| ---- WARNING ---- failed to read from RESTARTFARM, starting at 1"
             ELSE
                write(globenv%scr,"(T2,A,T71,I10)") &
                     "FARMING| restarting from RESTARTFARM, starting at ",farming_env%restart_n
             ENDIF
           ELSE
             write(globenv%scr,"(T2,A)") &
                     "FARMING| ---- WARNING ---- failed to open RESTARTFARM, starting at 1"
           ENDIF
           CLOSE(iunit,IOSTAT=stat)
        ENDIF
    ENDIF
    CALL mp_bcast(farming_env%restart_n,globenv%source,globenv%group)

  END SUBROUTINE

END MODULE farming_methods
!!*****
