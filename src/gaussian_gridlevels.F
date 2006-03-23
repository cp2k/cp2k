!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!


!!****h* cp2k/gaussian_gridlevels [1.0] *
!!
!!   NAME
!!     gaussian_gridlevels
!!
!!   AUTHOR
!!    Joost VandeVondele (27.02.02)
!!
!!   MODIFICATION HISTORY
!!     Code to return a gridlevel associated with a given gaussian exponent
!!
!!   SOURCE
!****************************************************************************
MODULE gaussian_gridlevels
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type, section_vals_retain, section_vals_release
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE termination,                     ONLY: stop_memory
#include "cp_common_uses.h"

  IMPLICIT NONE

  TYPE gridlevel_info_type
       INTEGER                         :: ngrid_levels
       REAL(KIND = dp), POINTER, DIMENSION(:) :: cutoff
       INTEGER , POINTER, DIMENSION(:) :: count
       INTEGER                         :: total_count
       REAL(KIND = dp)                  :: rel_cutoff
       TYPE(section_vals_type), POINTER :: print_section
  END TYPE

  PUBLIC :: gridlevel_info_type
  PUBLIC :: gaussian_gridlevel
  PUBLIC :: init_gaussian_gridlevel
  PUBLIC :: destroy_gaussian_gridlevel

CONTAINS

  SUBROUTINE init_gaussian_gridlevel(gridlevel_info,ngrid_levels,cutoff,rel_cutoff, print_section)
    TYPE(gridlevel_info_type)                :: gridlevel_info
    INTEGER                                  :: ngrid_levels
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cutoff
    REAL(KIND=dp)                            :: rel_cutoff
    TYPE(section_vals_type), POINTER         :: print_section

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
       gridlevel_info%print_section=>print_section
       CALL section_vals_retain(print_section)
  END SUBROUTINE
 
SUBROUTINE destroy_gaussian_gridlevel(gridlevel_info,globenv,error)
    TYPE(gridlevel_info_type)                :: gridlevel_info
    TYPE(global_environment_type), &
      OPTIONAL, POINTER                      :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routine = &
      "SUBROUTINE destroy_gaussian_gridlevel (MODULE gaussian_gridlevels)"

    INTEGER                                  :: group, i, istat, output_unit
    TYPE(cp_logger_type), POINTER            :: logger

  NULLIFY(logger)
  logger => cp_error_get_logger(error)
  IF (PRESENT(globenv)) THEN
     group=globenv%group
  ELSE
     group=logger%para_env%group
  END IF
  output_unit = cp_print_key_unit_nr(logger,gridlevel_info%print_section, &
                   "", extension=".Log",error=error)

  DEALLOCATE (gridlevel_info%cutoff,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"gridlevel_info%cutoff")
  CALL mp_sum(gridlevel_info%total_count,group)
  CALL mp_sum(gridlevel_info%count,group)
  IF (output_unit>0) THEN
    WRITE(output_unit,'(/,T3,A)')         "---------------------------------------"
    WRITE(output_unit,'(T3,A)')           "--           MULTIGRID INFO          --"
    WRITE(output_unit,'(T3,A)')           "---------------------------------------"
    WRITE(output_unit,'(T3,A,I14)')        "total gridlevel count  : ",&
      gridlevel_info%total_count
    DO i=1,gridlevel_info%ngrid_levels
      WRITE(output_unit,'(T3,A,I4,A,I14)') "count for grid     ",i,": ",&
       gridlevel_info%count(i)
    END DO
  END IF

  CALL cp_print_key_finished_output(output_unit,logger,gridlevel_info%print_section,&
                                    "", error=error)  

  CALL section_vals_release(gridlevel_info%print_section)

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
