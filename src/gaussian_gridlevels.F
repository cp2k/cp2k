!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      Code to return a gridlevel associated with a given gaussian exponent
!> \author Joost VandeVondele (27.02.02)
! *****************************************************************************
MODULE gaussian_gridlevels
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_section_types,             ONLY: section_vals_release,&
                                             section_vals_retain,&
                                             section_vals_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE message_passing,                 ONLY: mp_sum
  USE termination,                     ONLY: stop_memory
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gaussian_gridlevels'

! *****************************************************************************
  TYPE gridlevel_info_type
     INTEGER                         :: ngrid_levels
     REAL(KIND = dp), POINTER, DIMENSION(:) :: cutoff
     INTEGER(KIND=int_8) , POINTER, DIMENSION(:) :: count
     INTEGER(KIND=int_8)                         :: total_count
     REAL(KIND = dp)                  :: rel_cutoff
     TYPE(section_vals_type), POINTER :: print_section
  END TYPE gridlevel_info_type

  PUBLIC :: gridlevel_info_type
  PUBLIC :: gaussian_gridlevel
  PUBLIC :: init_gaussian_gridlevel
  PUBLIC :: destroy_gaussian_gridlevel

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param gridlevel_info ...
!> \param ngrid_levels ...
!> \param cutoff ...
!> \param rel_cutoff ...
!> \param print_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE init_gaussian_gridlevel(gridlevel_info,ngrid_levels,cutoff,rel_cutoff, print_section, error)
    TYPE(gridlevel_info_type)                :: gridlevel_info
    INTEGER                                  :: ngrid_levels
    REAL(KIND=dp), DIMENSION(:), POINTER     :: cutoff
    REAL(KIND=dp)                            :: rel_cutoff
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(cp_error_type), INTENT(inout)       :: error

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
    CALL section_vals_retain(print_section,error=error)
  END SUBROUTINE init_gaussian_gridlevel

! *****************************************************************************
!> \brief ...
!> \param gridlevel_info ...
!> \param para_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE destroy_gaussian_gridlevel(gridlevel_info,para_env,error)
    TYPE(gridlevel_info_type)                :: gridlevel_info
    TYPE(cp_para_env_type), OPTIONAL, &
      POINTER                                :: para_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'destroy_gaussian_gridlevel', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group, i, istat, output_unit
    TYPE(cp_logger_type), POINTER            :: logger

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    IF (PRESENT(para_env)) THEN
       group=para_env%group
    ELSE
       group=logger%para_env%group
    END IF
    output_unit = cp_print_key_unit_nr(logger,gridlevel_info%print_section, &
         "", extension=".Log",error=error)

    CALL mp_sum(gridlevel_info%total_count,group)
    CALL mp_sum(gridlevel_info%count,group)
    IF (output_unit>0) THEN
       WRITE(output_unit,'(/,T2,A,A)') "----------------------------------------",&
            "---------------------------------------"
       WRITE(output_unit,'(T2,A,T35,A,T77,A)') "----","MULTIGRID INFO","----"
       WRITE(output_unit,'(T2,A,A)') "----------------------------------------",&
            "---------------------------------------"
       IF ( gridlevel_info%ngrid_levels > 1 ) THEN
          DO i=1,gridlevel_info%ngrid_levels
             WRITE(output_unit,'(T2,A,I4,A,I14,9x,A,F12.2)') "count for grid     ",i,": ",&
                  gridlevel_info%count(i)," cutoff [a.u.]    ",gridlevel_info%cutoff(i)
          END DO
          WRITE(output_unit,'(T2,A,I14)')        "total gridlevel count  : ",&
               gridlevel_info%total_count
       ELSE
          WRITE(output_unit,'(T2,A,I14,T51,A,F12.2)') "total grid count     :",&
               gridlevel_info%count(1)," cutoff [a.u.]    ",gridlevel_info%cutoff(1)
       END IF
    END IF

    DEALLOCATE (gridlevel_info%cutoff,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "gridlevel_info%cutoff")

    CALL cp_print_key_finished_output(output_unit,logger,gridlevel_info%print_section,&
         "", error=error)

    CALL section_vals_release(gridlevel_info%print_section,error=error)

    DEALLOCATE (gridlevel_info%count,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "gridlevel_info%count")

  END SUBROUTINE destroy_gaussian_gridlevel

! *****************************************************************************
!> \brief ...
!> \param gridlevel_info ...
!> \param exponent ...
!> \retval gridlevel ...
! *****************************************************************************
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

END MODULE gaussian_gridlevels
