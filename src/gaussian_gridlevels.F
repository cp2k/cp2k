  ! Joost VandeVondele (27.02.02)
  ! Code to return a gridlevel associated with a given gaussian exponent



  MODULE gaussian_gridlevels

  USE kinds, ONLY: wp => dp
  USE global_types, ONLY: global_environment_type
  USE message_passing, ONLY: mp_sum
  USE termination, ONLY: stop_memory,stop_program

  IMPLICIT NONE

  TYPE gridlevel_info_type
       integer                         :: ngrid_levels
       real(wp), pointer, dimension(:) :: cutoff
       integer , pointer, dimension(:) :: count
       integer                         :: total_count
       real(wp)                        :: rel_cutoff
  END TYPE

  PUBLIC :: gridlevel_info_type
  PUBLIC :: gaussian_gridlevel
  PUBLIC :: init_gaussian_gridlevel
  PUBLIC :: destroy_gaussian_gridlevel

CONTAINS

  SUBROUTINE init_gaussian_gridlevel(gridlevel_info,ngrid_levels,cutoff,rel_cutoff)
       integer                :: ngrid_levels,i
       real(wp), dimension(:), POINTER :: cutoff
       real(wp)               :: rel_cutoff
       TYPE(gridlevel_info_type) :: gridlevel_info
       allocate(gridlevel_info%cutoff(ngrid_levels))
       allocate(gridlevel_info%count(ngrid_levels))
       gridlevel_info%ngrid_levels=ngrid_levels
       gridlevel_info%rel_cutoff=rel_cutoff
       gridlevel_info%total_count=0
       do i=1,ngrid_levels
          gridlevel_info%cutoff(i)=cutoff(i)
          gridlevel_info%count(i)=0
       enddo
  END SUBROUTINE
 
  SUBROUTINE destroy_gaussian_gridlevel(gridlevel_info,globenv)
  TYPE(global_environment_type), INTENT(IN) :: globenv
  TYPE(gridlevel_info_type) :: gridlevel_info
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE destroy_gaussian_gridlevel (MODULE gaussian_gridlevels)"
  INTEGER :: i,istat
  DEALLOCATE (gridlevel_info%cutoff,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"gridlevel_info%cutoff")
  CALL mp_sum(gridlevel_info%total_count,globenv%group)
  CALL mp_sum(gridlevel_info%count,globenv%group)
  IF (globenv%ionode.AND.globenv%print%pw_grid_information) THEN
    write(6,'(/,T3,A)')           "---------------------------------"
    write(6,'(T3,A)')           "--        MULTIGRID INFO       --"
    write(6,'(T3,A)')           "---------------------------------"
    write(6,'(T3,A,I8)')        "total gridlevel count  : ",&
      gridlevel_info%total_count
    DO i=1,gridlevel_info%ngrid_levels
      write(6,'(T3,A,I3,A,I8)') "count for grid     ",i,": ",&
       gridlevel_info%count(i)
    END DO
  END IF
  DEALLOCATE (gridlevel_info%count,STAT=istat)
  IF (istat /= 0) CALL stop_memory(routine,"gridlevel_info%count")
  END SUBROUTINE

! ***************************************************************
  FUNCTION gaussian_gridlevel(gridlevel_info,exponent) result(gridlevel)
  real(wp) :: exponent
  integer  :: gridlevel,i
  TYPE(gridlevel_info_type) :: gridlevel_info
  real :: needed_cutoff

  gridlevel=1 
  needed_cutoff=abs(exponent)*gridlevel_info%rel_cutoff
  do i=1,gridlevel_info%ngrid_levels
     if ((gridlevel_info%cutoff(i)+1D-6).ge.needed_cutoff) then
        gridlevel=i
     endif
  enddo
  gridlevel_info%total_count=gridlevel_info%total_count+1
  gridlevel_info%count(gridlevel)= gridlevel_info%count(gridlevel)+1

  END FUNCTION gaussian_gridlevel
! ***************************************************************

END MODULE
