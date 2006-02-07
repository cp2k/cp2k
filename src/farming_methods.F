!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2004  CP2K developers group
!-----------------------------------------------------------------------------!

MODULE farming_methods
  USE cp_files,                        ONLY: get_unit_number
  USE farming_types,                   ONLY: farming_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE message_passing,                 ONLY: mp_bcast
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: farming_parse_input

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'farming_methods'

CONTAINS

  SUBROUTINE farming_parse_input(farming_env,globenv,error)
    TYPE(farming_env_type), POINTER          :: farming_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'farming_parse_input', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=3)                         :: text
    INTEGER                                  :: i, iunit, n_rep_val, stat
    INTEGER, DIMENSION(:), POINTER           :: i_vals
    LOGICAL                                  :: failure, found
    TYPE(section_vals_type), POINTER         :: farming_section, jobs_section

    failure=.FALSE.

   farming_env%group_size_wish_set=.FALSE.
   farming_env%ngroup_wish_set    =.FALSE.
   farming_section => section_vals_get_subs_vals(globenv%input_file,"FARMING",&
        error=error)

   IF (ASSOCIATED(farming_env%group_partition)) THEN
      DEALLOCATE(farming_env%group_partition,stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
   END IF
   CALL section_vals_val_get(farming_section,"group_partition",&
        n_rep_val=n_rep_val,error=error)
   IF (n_rep_val>0) THEN
      CALL section_vals_val_get(farming_section,"group_partition",&
           i_vals=i_vals,error=error)
      ALLOCATE(farming_env%group_partition(0:SIZE(i_vals)-1))
      farming_env%group_partition(:)=i_vals
      farming_env%ngroup_wish_set=.TRUE.
      farming_env%ngroup_wish=SIZE(i_vals)
   ELSE
      CALL section_vals_val_get(farming_section,"GROUP_SIZE",&
           n_rep_val=n_rep_val,error=error)
      IF (n_rep_val>0) THEN
         CALL section_vals_val_get(farming_section,"GROUP_SIZE",&
              i_val=farming_env%ngroup_wish,error=error)
         farming_env%ngroup_wish_set=.TRUE.
      ELSE
         CALL section_vals_val_get(farming_section,"GROUP_SIZE",&
              i_val=farming_env%group_size_wish,error=error)
         farming_env%group_size_wish_set=.FALSE.
      END IF
   END IF
   CALL section_vals_val_get(farming_section,"restart",&
        l_val=farming_env%restart,error=error)      
   CALL section_vals_val_get(farming_section,"max_steps",&
        i_val=farming_env%max_steps,error=error)      
   CALL section_vals_val_get(farming_section,"cycle",&
        l_val=farming_env%cycle,error=error)      

   jobs_section => section_vals_get_subs_vals(farming_section,"JOB",error=error)
   CALL section_vals_get(jobs_section,n_repetition=farming_env % njobs,error=error)
   
   ALLOCATE(farming_env%Job_cwd(farming_env % njobs),STAT=stat)
   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
   ALLOCATE(farming_env%Job_input(farming_env % njobs),STAT=stat)
   CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
   IF (.NOT. failure) THEN
      DO i=1,farming_env%njobs
         CALL section_vals_val_get(jobs_section,i_rep_section=i,&
              keyword_name="directory",c_val=farming_env%Job_cwd(i),error=error)
         CALL section_vals_val_get(jobs_section,i_rep_section=i,&
              keyword_name="input_file_name",c_val=farming_env%Job_input(i),error=error)
      END DO
   END IF

   IF (globenv%ionode .AND. .NOT. failure) THEN
        WRITE(globenv%scr,FMT="(T2,A,T71,I10)") "FARMING| Number of jobs found",farming_env%njobs
        IF (farming_env%ngroup_wish_set) THEN
           WRITE(globenv%scr,FMT="(T2,A,T71,I10)") "FARMING| Ngroup wish:",farming_env%ngroup_wish
           IF (ASSOCIATED(farming_env%group_partition)) THEN
              WRITE(globenv%scr,FMT="(T2,A)",ADVANCE="NO") "FARMING| User partition:"
              DO i=0,SIZE(farming_env%group_partition)-1
                 IF (MODULO(i,4)==0) WRITE(globenv%scr,*)
                 WRITE(globenv%scr,FMT='(I4)',ADVANCE="NO") farming_env%group_partition(i)
              END DO
              WRITE(globenv%scr,*)
              IF (SUM(farming_env%group_partition).NE.globenv%num_pe) THEN
                 WRITE(globenv%scr,FMT="(T2,A,T61,I10,T71,I10)") &
                     "FARMING| WARNING : partition CPUs not equal actual number ", &
                     globenv%num_pe, SUM(farming_env%group_partition)
                 WRITE(globenv%scr,FMT="(T2,A)") "FARMING|          partition data ignored" ! any better idea ??
                 DEALLOCATE(farming_env%group_partition)
              ENDIF
           ENDIF
        ENDIF
        IF (farming_env%group_size_wish_set) THEN
           WRITE(globenv%scr,FMT="(T2,A,T71,I10)") "FARMING| Group size wish:", &
                 farming_env%group_size_wish
        ENDIF
        WRITE(globenv%scr,FMT="(T2,A,T71,I10)") "FARMING| Max steps      :",farming_env%max_steps
        IF (farming_env%cycle) THEN
            text="YES"
        ELSE
            text=" NO"
        ENDIF
        WRITE(globenv%scr,FMT="(T2,A,T78,A3)") "FARMING| Cyclic jobs execution:",text
        IF (farming_env%restart) THEN
            text="YES"
        ELSE
            text=" NO"
        ENDIF
        WRITE(globenv%scr,FMT="(T2,A,T78,A3)") "FARMING| Restarting farm:",text
        farming_env%restart_n = 1
        IF (farming_env%restart) THEN
           iunit=get_unit_number()
           OPEN(UNIT=iunit,FILE="RESTARTFARM",IOSTAT=stat)
           IF (stat==0) THEN
             READ(UNIT=iunit,FMT=*,IOSTAT=stat) farming_env%restart_n
             IF (stat/=0) THEN
                WRITE(globenv%scr,"(T2,A)") &
                     "FARMING| ---- WARNING ---- failed to read from RESTARTFARM, starting at 1"
             ELSE
                WRITE(globenv%scr,"(T2,A,T71,I10)") &
                     "FARMING| restarting from RESTARTFARM, starting at ",farming_env%restart_n
             ENDIF
           ELSE
             WRITE(globenv%scr,"(T2,A)") &
                     "FARMING| ---- WARNING ---- failed to open RESTARTFARM, starting at 1"
           ENDIF
           CLOSE(iunit,IOSTAT=stat)
        ENDIF
    ENDIF
    CALL mp_bcast(farming_env%restart_n,globenv%source,globenv%group)

  END SUBROUTINE

END MODULE farming_methods
!!*****
