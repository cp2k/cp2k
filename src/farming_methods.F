!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2006  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE farming_methods
  USE cp_files,                        ONLY: get_unit_number
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE farming_types,                   ONLY: farming_env_type
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

  SUBROUTINE farming_parse_input(farming_env,root_section,para_env,error)
    TYPE(farming_env_type), POINTER          :: farming_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)  &
                                             :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'farming_parse_input', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=3)                         :: text
    INTEGER                                  :: i, iunit, n_rep_val, &
                                                output_unit, stat
    INTEGER, DIMENSION(:), POINTER           :: i_vals
    LOGICAL                                  :: failure, found
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: farming_section, jobs_section

    failure=.FALSE.

   farming_env%group_size_wish_set=.FALSE.
   farming_env%ngroup_wish_set    =.FALSE.
   farming_section => section_vals_get_subs_vals(root_section,"FARMING",&
        error=error)

   IF (ASSOCIATED(farming_env%group_partition)) THEN
      DEALLOCATE(farming_env%group_partition,stat=stat)
      CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
   END IF

   ! The following input order is used
   ! 1) GROUP_PARTITION
   ! 2) NGROUP
   ! 3) GROUP_SIZE (default 8)
   CALL section_vals_val_get(farming_section,"GROUP_PARTITION",&
        n_rep_val=n_rep_val,error=error)
   IF (n_rep_val>0) THEN
      CALL section_vals_val_get(farming_section,"GROUP_PARTITION",&
           i_vals=i_vals,error=error)
      ALLOCATE(farming_env%group_partition(0:SIZE(i_vals)-1))
      farming_env%group_partition(:)=i_vals
      farming_env%ngroup_wish_set=.TRUE.
      farming_env%ngroup_wish=SIZE(i_vals)
   ELSE
      CALL section_vals_val_get(farming_section,"NGROUP",&
           n_rep_val=n_rep_val,error=error)
      IF (n_rep_val>0) THEN
         CALL section_vals_val_get(farming_section,"NGROUP",&
              i_val=farming_env%ngroup_wish,error=error)
         farming_env%ngroup_wish_set=.TRUE.
      ELSE
         CALL section_vals_val_get(farming_section,"GROUP_SIZE",&
              i_val=farming_env%group_size_wish,error=error)
         farming_env%group_size_wish_set=.TRUE.
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

   IF (.NOT. failure) THEN
     
      logger => cp_error_get_logger(error)
      output_unit=cp_print_key_unit_nr(logger,root_section,"GLOBAL%PROGRAM_RUN_INFO",&
                                         extension=".log",error=error)
      IF (output_unit > 0) THEN
        WRITE(output_unit,FMT="(T2,A,T71,I10)") "FARMING| Number of jobs found",farming_env%njobs
        IF (farming_env%ngroup_wish_set) THEN
           WRITE(output_unit,FMT="(T2,A,T71,I10)") "FARMING| Ngroup wish:",farming_env%ngroup_wish
           IF (ASSOCIATED(farming_env%group_partition)) THEN
              WRITE(output_unit,FMT="(T2,A)",ADVANCE="NO") "FARMING| User partition:"
              DO i=0,SIZE(farming_env%group_partition)-1
                 IF (MODULO(i,4)==0) WRITE(output_unit,*)
                 WRITE(output_unit,FMT='(I4)',ADVANCE="NO") farming_env%group_partition(i)
              END DO
              WRITE(output_unit,*)
              IF (SUM(farming_env%group_partition).NE.para_env%num_pe) THEN
                 WRITE(output_unit,FMT="(T2,A,T61,I10,T71,I10)") &
                     "FARMING| WARNING : partition CPUs not equal actual number ", &
                     para_env%num_pe, SUM(farming_env%group_partition)
                 WRITE(output_unit,FMT="(T2,A)") "FARMING|          partition data ignored" ! any better idea ??
                 DEALLOCATE(farming_env%group_partition)
              ENDIF
           ENDIF
        ENDIF
        IF (farming_env%group_size_wish_set) THEN
           WRITE(output_unit,FMT="(T2,A,T71,I10)") "FARMING| Group size wish:", &
                 farming_env%group_size_wish
        ENDIF
        WRITE(output_unit,FMT="(T2,A,T71,I10)") "FARMING| Max steps      :",farming_env%max_steps
        IF (farming_env%cycle) THEN
            text="YES"
        ELSE
            text=" NO"
        ENDIF
        WRITE(output_unit,FMT="(T2,A,T78,A3)") "FARMING| Cyclic jobs execution:",text
        IF (farming_env%restart) THEN
            text="YES"
        ELSE
            text=" NO"
        ENDIF
        WRITE(output_unit,FMT="(T2,A,T78,A3)") "FARMING| Restarting farm:",text
        farming_env%restart_n = 1
        IF (farming_env%restart) THEN
           iunit=get_unit_number()
           OPEN(UNIT=iunit,FILE="RESTARTFARM",IOSTAT=stat)
           IF (stat==0) THEN
             READ(UNIT=iunit,FMT=*,IOSTAT=stat) farming_env%restart_n
             IF (stat/=0) THEN
                WRITE(output_unit,"(T2,A)") &
                     "FARMING| ---- WARNING ---- failed to read from RESTARTFARM, starting at 1"
             ELSE
                WRITE(output_unit,"(T2,A,T71,I10)") &
                     "FARMING| restarting from RESTARTFARM, starting at ",farming_env%restart_n
             ENDIF
           ELSE
             WRITE(output_unit,"(T2,A)") &
                     "FARMING| ---- WARNING ---- failed to open RESTARTFARM, starting at 1"
           ENDIF
           CLOSE(iunit,IOSTAT=stat)
        ENDIF

        CALL cp_print_key_finished_output(output_unit,logger,root_section,&
               "GLOBAL%PROGRAM_RUN_INFO", error=error)
      ENDIF
    ENDIF
    CALL mp_bcast(farming_env%restart_n,para_env%source,para_env%group)

  END SUBROUTINE

END MODULE farming_methods
!!*****
