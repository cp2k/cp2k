  ! Joost VandeVondele (27.02.02)
  ! Code to return a gridlevel associated with a given gaussian exponent



  MODULE gaussian_gridlevels

  USE kinds, ONLY: wp => dp
  USE termination, ONLY: stop_program

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
  IMPLICIT NONE
       integer                :: ngrid_levels,i
       real(wp), dimension(:) :: cutoff
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
 
  SUBROUTINE destroy_gaussian_gridlevel(gridlevel_info)
  IMPLICIT NONE
  TYPE(gridlevel_info_type) :: gridlevel_info
  integer i
  deallocate(gridlevel_info%cutoff)
  write(6,'(T3,A)')           " "
  write(6,'(T3,A)')           "---------------------------------"
  write(6,'(T3,A)')           "--        MULTIGRID INFO       --"
  write(6,'(T3,A)')           "---------------------------------"
  write(6,'(T3,A,I8)')        "total gridlevel count  : ",gridlevel_info%total_count
  DO i=1,gridlevel_info%ngrid_levels
    write(6,'(T3,A,I3,A,I8)') "count for grid     ",i," : ",gridlevel_info%count(i)
  ENDDO
  deallocate(gridlevel_info%count)
  END SUBROUTINE

! ***************************************************************
  FUNCTION gaussian_gridlevel(gridlevel_info,exponent) result(gridlevel)
  IMPLICIT NONE
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

  RETURN
  END FUNCTION gaussian_gridlevel
! ***************************************************************

END MODULE
