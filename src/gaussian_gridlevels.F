  ! Joost VandeVondele (27.02.02)
  ! Code to return a gridlevel associated with a given gaussian exponent


  MODULE gaussian_gridlevels
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level,&
                                             cp_logger_get_default_unit_nr
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE termination,                     ONLY: stop_memory

  IMPLICIT NONE

  TYPE gridlevel_info_type
       INTEGER                         :: ngrid_levels
       REAL(KIND = dp), POINTER, DIMENSION(:) :: cutoff
       INTEGER , POINTER, DIMENSION(:) :: count
       INTEGER                         :: total_count
       REAL(KIND = dp)                        :: rel_cutoff
  END TYPE

  PUBLIC :: gridlevel_info_type
  PUBLIC :: gaussian_gridlevel
  PUBLIC :: init_gaussian_gridlevel
  PUBLIC :: destroy_gaussian_gridlevel

CONTAINS

  SUBROUTINE init_gaussian_gridlevel(gridlevel_info,ngrid_levels,cutoff,rel_cutoff)
    TYPE(gridlevel_info_type)                :: gridlevel_info
    INTEGER                                  :: ngrid_levels
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cutoff
    REAL(KIND=dp)                            :: rel_cutoff

    INTEGER                                  :: i

       ALLOCATE(gridlevel_info%cutoff(ngrid_levels))
       ALLOCATE(gridlevel_info%count(ngrid_levels))
       gridlevel_info%ngrid_levels=ngrid_levels
       gridlevel_info%rel_cutoff=rel_cutoff
       gridlevel_info%total_count=0
       DO i=1,ngrid_levels
          gridlevel_info%cutoff(i)=cutoff(i)
          gridlevel_info%count(i)=0
       ENDDO
  END SUBROUTINE
 
SUBROUTINE destroy_gaussian_gridlevel(gridlevel_info,globenv,error)
    TYPE(gridlevel_info_type)                :: gridlevel_info
    TYPE(global_environment_type), &
      INTENT(in), OPTIONAL, TARGET           :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE destroy_gaussian_gridlevel (MODULE gaussian_gridlevels)"

    INTEGER                                  :: group, i, istat,scr
    LOGICAL                                  :: ionode, pw_grid_information
    TYPE(cp_logger_type), POINTER            :: logger

  NULLIFY(logger)
  IF (PRESENT(globenv)) THEN
     group=globenv%group
     ionode=globenv%ionode
     pw_grid_information=globenv%print%pw_grid_information
     scr=globenv%scr
  ELSE
     logger => cp_error_get_logger(error)
     group=logger%para_env%group
     pw_grid_information=logger%print_keys%pw_grid_information
     ionode=logger%para_env%source==logger%para_env%mepos
     IF (ionode) scr=cp_logger_get_default_unit_nr(logger)
  END IF

  DEALLOCATE (gridlevel_info%cutoff,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"gridlevel_info%cutoff")
  CALL mp_sum(gridlevel_info%total_count,group)
  CALL mp_sum(gridlevel_info%count,group)
  IF (ionode.AND.pw_grid_information) THEN
    WRITE(scr,'(/,T3,A)')         "---------------------------------------"
    WRITE(scr,'(T3,A)')           "--           MULTIGRID INFO          --"
    WRITE(scr,'(T3,A)')           "---------------------------------------"
    WRITE(scr,'(T3,A,I14)')        "total gridlevel count  : ",&
      gridlevel_info%total_count
    DO i=1,gridlevel_info%ngrid_levels
      WRITE(scr,'(T3,A,I4,A,I14)') "count for grid     ",i,": ",&
       gridlevel_info%count(i)
    END DO
  END IF
  DEALLOCATE (gridlevel_info%count,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"gridlevel_info%count")
  END SUBROUTINE

! ***************************************************************
  FUNCTION gaussian_gridlevel(gridlevel_info,exponent) RESULT(gridlevel)
    TYPE(gridlevel_info_type), INTENT(INOUT) :: gridlevel_info
    REAL(KIND=dp)                            :: exponent
    INTEGER                                  :: gridlevel

    INTEGER                                  :: i
    REAL(KIND=dp)                            :: needed_cutoff

  gridlevel=1 
  needed_cutoff=ABS(exponent)*gridlevel_info%rel_cutoff
  DO i=1,gridlevel_info%ngrid_levels
     IF ((gridlevel_info%cutoff(i)+1E-6_dp).GE.needed_cutoff) THEN
        gridlevel=i
     ENDIF
  ENDDO
  gridlevel_info%total_count=gridlevel_info%total_count+1
  gridlevel_info%count(gridlevel)= gridlevel_info%count(gridlevel)+1

  END FUNCTION gaussian_gridlevel
! ***************************************************************

END MODULE
