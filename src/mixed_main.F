!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief perform biased molecular dynamics (H= k H1 + (1-k) H2 [linear or general mixing)
!> \author fschiff 11.06
! *****************************************************************************
MODULE mixed_main
  USE cp_files,                        ONLY: open_file
  USE cp_output_handling,              ONLY: cp_print_key_unit_nr
  USE cp_para_env,                     ONLY: cp_para_env_create
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_path_length
  USE message_passing,                 ONLY: mp_comm_split
  USE mixed_environment,               ONLY: mixed_init
  USE mixed_environment_types,         ONLY: mixed_env_create,&
                                             mixed_environment_type,&
                                             set_mixed_env
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

! *** Global parameters ***
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mixed_main'

  PUBLIC :: mixed_create_force_env

CONTAINS
! *****************************************************************************
!> \brief Controls program flow for mixed calculations
!> \param mixed_env ...
!> \param root_section ...
!> \param para_env ...
!> \param force_env_section ...
!> \param n_subforce_eval ...
!> \param use_motion_section ...
!> \author fschiff
! *****************************************************************************
  SUBROUTINE mixed_create_force_env ( mixed_env, root_section, para_env,&
       force_env_section, n_subforce_eval, use_motion_section)

    TYPE(mixed_environment_type), POINTER    :: mixed_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(section_vals_type), POINTER         :: force_env_section
    INTEGER, INTENT(IN)                      :: n_subforce_eval
    LOGICAL, INTENT(IN)                      :: use_motion_section

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mixed_create_force_env', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: c_val, input_file_path, &
                                                output_file_path
    INTEGER                                  :: group_size_wish, handle, i, &
                                                lp, n_rep_val, ngroup_wish, &
                                                output_unit, unit_nr
    INTEGER, DIMENSION(:), POINTER           :: group_partition, i_vals
    LOGICAL                                  :: failure, group_size_wish_set, &
                                                ngroup_wish_set
    TYPE(cp_logger_type), POINTER            :: logger

    failure = .FALSE.
    CALL timeset(routineN,handle)
    logger => cp_get_default_logger()
    output_unit=cp_print_key_unit_nr(logger,force_env_section,"MIXED%PRINT%PROGRAM_RUN_INFO",&
                                         extension=".log")

    CALL mixed_env_create( mixed_env, para_env = para_env)
    ! Setup the new parallel env
    NULLIFY(group_partition)
    group_size_wish_set=.FALSE.
    ngroup_wish_set    =.FALSE.
    CALL section_vals_val_get(force_env_section,"MIXED%GROUP_PARTITION",n_rep_val=n_rep_val)
    IF (n_rep_val>0) THEN
       CALL section_vals_val_get(force_env_section,"MIXED%GROUP_PARTITION",i_vals=i_vals)
       ALLOCATE(group_partition(0:SIZE(i_vals)-1))
       group_partition(:)=i_vals
       ngroup_wish_set=.TRUE.
       ngroup_wish=SIZE(i_vals)
    ELSE
       CALL section_vals_val_get(force_env_section,"MIXED%NGROUPS",n_rep_val=n_rep_val)
       IF (n_rep_val>0) THEN
          CALL section_vals_val_get(force_env_section,"MIXED%NGROUPS",i_val=ngroup_wish)
       ELSE
          ngroup_wish = n_subforce_eval
       END IF
       group_size_wish_set=.TRUE.
       group_size_wish = MAX(1,para_env%num_pe/ ngroup_wish)
    END IF

    ! Split the current communicator
    ALLOCATE(mixed_env%group_distribution(0:para_env%num_pe-1))
    IF (group_size_wish_set) THEN
       CALL mp_comm_split(para_env%group,mixed_env%new_group,mixed_env%ngroups,mixed_env%group_distribution,&
            subgroup_min_size=group_size_wish)
    ELSE IF (ngroup_wish_set) THEN
       CALL mp_comm_split(para_env%group,mixed_env%new_group,mixed_env%ngroups,mixed_env%group_distribution,&
            n_subgroups=ngroup_wish,&
            group_partition=group_partition)
    ENDIF
    IF (output_unit>0) THEN
       WRITE(output_unit,FMT="(T2,A,T71,I10)") "MIXED_ENV| Number of created MPI groups:",mixed_env%ngroups
       WRITE(output_unit,FMT="(T2,A)",ADVANCE="NO") "MIXED_ENV| Task to group correspondence:"
       DO i=0,para_env%num_pe-1
          IF (MODULO(i,4)==0) WRITE(output_unit,*)
          WRITE(output_unit,FMT='(A3,I4,A3,I4,A1)',ADVANCE="NO")&
               "  (",i," : ",mixed_env%group_distribution(i),")"
       END DO
       WRITE(output_unit,*)
    ENDIF
    IF (ASSOCIATED(group_partition)) THEN
       DEALLOCATE(group_partition)
    END IF
    ! Allocate para_env and handle the several loggers
    ALLOCATE(mixed_env%sub_para_env(mixed_env%ngroups))
    ALLOCATE(mixed_env%sub_logger(mixed_env%ngroups))
    ALLOCATE(mixed_env%energies(n_subforce_eval))
    !
    DO i = 1, mixed_env%ngroups
       NULLIFY(mixed_env%sub_para_env(i)%para_env,logger)
       IF (MODULO(i-1,mixed_env%ngroups)==mixed_env%group_distribution(para_env%mepos)) THEN
          ! Create sub_para_env
          CALL cp_para_env_create(mixed_env%sub_para_env(i)%para_env,&
               group=mixed_env%new_group,&
               owns_group=.TRUE.)
          ! Create sub_logger
          IF (mixed_env%sub_para_env(i)%para_env%mepos==mixed_env%sub_para_env(i)%para_env%source) THEN
             ! Redirecting output of subforce_eval to file..
             CALL section_vals_val_get(root_section,"GLOBAL%PROJECT_NAME",&
                  c_val=input_file_path)
             lp=LEN_TRIM(input_file_path)
             input_file_path(lp+1:LEN(input_file_path))="-r-"//&
                  ADJUSTL(cp_to_string(i))
             lp=LEN_TRIM(input_file_path)
             output_file_path=input_file_path(1:lp)//".out"
             CALL open_file(file_name=output_file_path,file_status="UNKNOWN",&
                  file_action="WRITE", file_position="APPEND",&
                  unit_number=unit_nr)
          ELSE
             unit_nr=-1
          END IF
          CALL cp_logger_create(mixed_env%sub_logger(i)%p,&
               para_env=mixed_env%sub_para_env(i)%para_env,&
               default_global_unit_nr=unit_nr, &
               close_global_unit_on_dealloc=.FALSE.)
          ! Try to use better names for the local log if it is not too late
          CALL section_vals_val_get(root_section,"GLOBAL%OUTPUT_FILE_NAME",&
               c_val=c_val)
          IF (c_val/="") THEN
             CALL cp_logger_set(mixed_env%sub_logger(i)%p,&
                  local_filename=TRIM(c_val)//"_localLog")
          END IF
          CALL section_vals_val_get(root_section,"GLOBAL%PROJECT",c_val=c_val)
          IF (c_val/="") THEN
             CALL cp_logger_set(mixed_env%sub_logger(i)%p,&
                  local_filename=TRIM(c_val)//"_localLog")
          END IF
          mixed_env%sub_logger(i)%p%iter_info%project_name=c_val
          CALL section_vals_val_get(root_section,"GLOBAL%PRINT_LEVEL",&
             i_val=mixed_env%sub_logger(i)%p%iter_info%print_level)

       END IF
    END DO

    ! *** initializations for the setup of the MIXED environment ***
    CALL set_mixed_env (mixed_env)
    CALL mixed_init ( mixed_env, root_section, para_env, force_env_section, &
                      use_motion_section)
    CALL timestop(handle)

  END SUBROUTINE mixed_create_force_env

END MODULE mixed_main
