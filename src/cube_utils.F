! for a given dr() this will provide the bounds to be used if 
! one wants to go over a sphere-subregion of given radius

MODULE cube_utils

  USE kinds, ONLY: dp
  USE timings,               ONLY: timeset,timestop
  USE termination,           ONLY: stop_memory,stop_program

  IMPLICIT NONE

  public :: cube_info_type
  public :: cube_ptr
  public :: init_cube_info
  public :: destroy_cube_info
  public :: return_cube

  TYPE :: cube_ptr
       integer, pointer, dimension(:) :: p
  END TYPE cube_ptr

  TYPE :: cube_info_type
     integer                      :: max_radius
     real(KIND = dp)                     :: dr(3),drmin,drsafe
     integer, pointer             :: lb_cube(:,:)
     integer, pointer             :: ub_cube(:,:)
     TYPE(cube_ptr), pointer, dimension(:)  :: sphere_bounds
     integer, pointer             :: sphere_bounds_count(:)
  END TYPE cube_info_type

CONTAINS
  SUBROUTINE return_cube(info,radius,lb_cube,ub_cube,sphere_bounds)
  IMPLICIT NONE
    TYPE (cube_info_type)    :: info
    real(KIND = dp)                 :: radius
    integer                  :: lb_cube(3),ub_cube(3)
    integer,dimension(:),pointer :: sphere_bounds
    integer                  :: imr

    imr=MAX(1,CEILING((radius+info%drsafe)/info%drmin))
    if (imr.gt.info%max_radius) then
        CALL stop_program("return_cube","called with radius too large")
    endif
    lb_cube(:)       = info%lb_cube(:,imr)
    ub_cube(:)       = info%ub_cube(:,imr)
    sphere_bounds => info%sphere_bounds(imr)%p

  RETURN
  END SUBROUTINE return_cube

  SUBROUTINE destroy_cube_info(info)
    TYPE (cube_info_type) :: info
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE destroy_cube_info (MODULE cube_utils)"
    INTEGER :: i,istat
    DEALLOCATE (info%lb_cube,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"info%lb_cube")
    DEALLOCATE (info%ub_cube,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"info%ub_cube")
    DEALLOCATE (info%sphere_bounds_count,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"info%sphere_bounds_count")
    DO i=1,info%max_radius
      DEALLOCATE (info%sphere_bounds(i)%p,STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"info%sphere_bounds(i)%p")
    END DO
    DEALLOCATE (info%sphere_bounds,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routine,"info%sphere_bounds")
  END SUBROUTINE

  SUBROUTINE init_cube_info(info,dr,max_radius,cmax)
  IMPLICIT NONE
    TYPE (cube_info_type) :: info 
    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE init_cube_info (MODULE cube_utils)"
    real(KIND = dp)       :: dr(3)
    real(KIND = dp)       :: max_radius
    integer        :: cmax
  
    real(KIND = dp)       :: drmin,drsafe,dxi,dyi,dzi,dz2,dy2
    real(KIND = dp)       :: radius,radius2
    integer        :: imr,i,istat,k,kgmin,kgmax,kg2,kg,jgmin
    integer        :: jgmax,jg2,jg,igmin,igmax,handle
    CALL timeset("init_cube_info","I","",handle) 
    info%dr = dr
    drmin = MINVAL(dr)
    info%drmin=drmin
    drsafe=0
    info%drsafe=drsafe
    imr=CEILING((max_radius+drsafe)/drmin)
    info%max_radius=imr
    dzi=1.0_dp/dr(3)
    dyi=1.0_dp/dr(2)
    dxi=1.0_dp/dr(1)
    dz2=(dr(3))**2
    dy2=(dr(2))**2

    NULLIFY(info%lb_cube,info%ub_cube,&
            info%sphere_bounds_count,info%sphere_bounds)
    allocate(info%lb_cube(3,imr),info%ub_cube(3,imr), &
             info%sphere_bounds_count(imr),info%sphere_bounds(imr))
    cmax=0
!   count and allocate
    DO i=1,imr
      k=1
      radius=i*drmin 
      radius2=radius**2
      kgmin=CEILING(-dzi*radius)
      k=k+1
      DO kg=kgmin,0
         kg2=kg*kg
         jgmin=CEILING(-dyi*sqrt(max(radius2-kg2*dz2,0.0_dp)))
         k=k+1
         DO jg=jgmin,0
          jg2=jg*jg
          igmin=CEILING(-dxi*sqrt(max(radius2-kg2*dz2-jg2*dy2,0.0_dp)))
          k=k+1
         ENDDO
      ENDDO
      info%sphere_bounds_count(i)=k-1
      ALLOCATE (info%sphere_bounds(i)%p(info%sphere_bounds_count(i)),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routine,"info%sphere_bounds(i)%p",0)
    ENDDO

!   init sphere_bounds array
    ! notice : as many points in lb_cube..0 as 1..ub_cube
    do i=1,imr 
      k=1
      radius=i*drmin
      info%lb_cube(:,i)=CEILING(-radius/dr(:))-1 ! take care in case of rounding
      info%ub_cube(:,i)=1-info%lb_cube(:,i)
      radius2=radius**2
      kgmin=CEILING(-dzi*radius)
      info%sphere_bounds(i)%p(k)=kgmin
      k=k+1
      DO kg=kgmin,0
         kg2=kg*kg
         jgmin=CEILING(-dyi*sqrt(max(radius2-kg2*dz2,0.0_dp)))
         info%sphere_bounds(i)%p(k)=jgmin
         k=k+1
         DO jg=jgmin,0
          jg2=jg*jg
          igmin=CEILING(-dxi*sqrt(max(radius2-kg2*dz2-jg2*dy2,0.0_dp)))
          info%sphere_bounds(i)%p(k)=igmin
          k=k+1
         ENDDO
      ENDDO
    enddo
    cmax=MAXVAL(info%ub_cube)
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE
END MODULE
! ***************************************************************
