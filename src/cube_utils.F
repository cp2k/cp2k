! for a given dr() this will provide the bounds to be used if 
! one wants to go over a sphere-subregion of given radius

! the computation of the exact sphere radius is sensitive to roundoff (e.g.
! compiler optimization level) and hence this small roundoff can result in
! energy difference of about EPS_DEFAULT in QS energies (one gridpoint more or
! less in the density mapping)

MODULE cube_utils

  USE kinds,                           ONLY: dp
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop

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
    TYPE(cube_info_type)                     :: info
    real(KIND=dp)                            :: radius
    integer                                  :: lb_cube(3), ub_cube(3)
    integer, dimension(:), pointer           :: sphere_bounds

    integer                                  :: imr

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
    TYPE(cube_info_type)                     :: info

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE destroy_cube_info (MODULE cube_utils)"

    INTEGER                                  :: i, istat

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
    TYPE(cube_info_type)                     :: info
    real(KIND=dp)                            :: dr(3), max_radius
    integer                                  :: cmax

    CHARACTER(LEN=*), PARAMETER :: &
      routine = "SUBROUTINE init_cube_info (MODULE cube_utils)"

    integer                                  :: check_1, check_2, handle, i, &
                                                igmin, imr, istat, jg, jg2, &
                                                jgmin, k, kg, kg2, kgmin
    real(KIND=dp)                            :: drmin, drsafe, dxi, dy2, dyi, &
                                                dz2, dzi, radius, radius2

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
    check_1=0
    check_2=0
!   count and allocate
    DO i=1,imr
      k=1
      radius=i*drmin 
      radius2=radius**2
      kgmin=do_and_hide_it_1(dzi,i,drmin,0.0_dp,0.0_dp,0,0)
      k=k+1
      DO kg=kgmin,0
         kg2=kg*kg
         jgmin=do_and_hide_it_1(dyi,i,drmin,dz2,0.0_dp,kg2,0)
         k=k+1
         DO jg=jgmin,0
          jg2=jg*jg
          igmin=do_and_hide_it_1(dxi,i,drmin,dz2,dy2,kg2,jg2)
          check_1=MODULO((kgmin*97+jgmin*37+igmin*113)*check_1,9343)
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
      info%lb_cube(:,i)=-1
      radius2=radius**2
      kgmin=do_and_hide_it_1(dzi,i,drmin,0.0_dp,0.0_dp,0,0)
      info%lb_cube(3,i)=MIN(kgmin,info%lb_cube(3,i))
      info%sphere_bounds(i)%p(k)=kgmin
      k=k+1
      DO kg=kgmin,0
         kg2=kg*kg
         jgmin=do_and_hide_it_1(dyi,i,drmin,dz2,0.0_dp,kg2,0)
         info%lb_cube(2,i)=MIN(jgmin,info%lb_cube(2,i))
         info%sphere_bounds(i)%p(k)=jgmin
         k=k+1
         DO jg=jgmin,0
          jg2=jg*jg
          igmin=do_and_hide_it_1(dxi,i,drmin,dz2,dy2,kg2,jg2)
          check_2=MODULO((kgmin*97+jgmin*37+igmin*113)*check_2,9343)
          info%lb_cube(1,i)=MIN(jgmin,info%lb_cube(1,i))
          info%sphere_bounds(i)%p(k)=igmin
          k=k+1
         ENDDO
      ENDDO
      info%ub_cube(:,i)=1-info%lb_cube(:,i)
    enddo
    cmax=MAXVAL(info%ub_cube)
    IF (check_1.NE.check_2) THEN 
        CALL stop_program("init_cube_info"," :-/ irreproducable fp math caused memory corruption")
    ENDIF
    CALL timestop(0.0_dp,handle)
  END SUBROUTINE
  ! try to hide things from the optimizer, so that we get the same numbers,
  ! always (this solves the optimisation problems with the intel and nag compiler
  ! in which the counting loops and execution loops above are executed a different
  ! number of times, even at -O1
  FUNCTION do_and_hide_it_1(prefactor,i,drmin,dz2,dy2,kg2,jg2) RESULT(res)
    REAL(KIND=dp)                            :: prefactor
    INTEGER                                  :: i
    REAL(KIND=dp)                            :: drmin, dz2, dy2
    INTEGER                                  :: kg2, jg2, res

    REAL(KIND=dp), DIMENSION(:), POINTER     :: buf

    ALLOCATE(buf(4))
    buf(1)=prefactor
    buf(2)=drmin
    buf(3)=dz2
    buf(4)=dy2
    res=do_and_hide_it_2(buf,i,jg2,kg2)
    DEALLOCATE(buf)
  END FUNCTION do_and_hide_it_1
  FUNCTION do_and_hide_it_2(buf,i,jg2,kg2) RESULT(res)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: buf
    INTEGER                                  :: i, jg2, kg2, res

    buf(2)=(i*buf(2))**2
    res = CEILING(-0.1D-7-buf(1)*sqrt(max(buf(2)-kg2*buf(3)-jg2*buf(4),0.0_dp)))
  END FUNCTION do_and_hide_it_2
END MODULE
! ***************************************************************
