!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief for a given dr()/dh(r) this will provide the bounds to be used if
!>      one wants to go over a sphere-subregion of given radius
!> \note
!>      the computation of the exact sphere radius is sensitive to roundoff (e.g.
!>      compiler optimization level) and hence this small roundoff can result in
!>      energy difference of about EPS_DEFAULT in QS energies (one gridpoint more or
!>      less in the density mapping)
!> \author Joost VandeVondele
! *****************************************************************************
MODULE cube_utils

  
  USE kinds,                           ONLY: dp,&
                                             int_size
  USE mathlib,                         ONLY: matvec_3x3
  USE realspace_grid_types,            ONLY: realspace_grid_desc_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cube_utils'

  PUBLIC :: cube_info_type

  PUBLIC :: init_cube_info, destroy_cube_info,&
            return_cube, return_cube_max_iradius, return_cube_nonortho,&
            compute_cube_center

  TYPE :: cube_ptr
       INTEGER, POINTER, DIMENSION(:) :: p
  END TYPE cube_ptr

  TYPE :: cube_info_type
     INTEGER                      :: max_radius
     REAL(KIND = dp)              :: dr(3),drmin
     REAL(KIND = dp)              :: dh(3,3)
     REAL(KIND = dp)              :: dh_inv(3,3)
     LOGICAL                      :: orthorhombic
     INTEGER, POINTER             :: lb_cube(:,:)
     INTEGER, POINTER             :: ub_cube(:,:)
     TYPE(cube_ptr), POINTER, DIMENSION(:)  :: sphere_bounds
     INTEGER, POINTER             :: sphere_bounds_count(:)
     REAL(KIND = dp)              :: max_rad_ga
  END TYPE cube_info_type

CONTAINS
! *****************************************************************************
!> \brief unifies the computation of the cube center, so that differences in
!>        implementation, and thus roundoff and numerics can not lead to
!>        off-by-one errors (which can lead to out-of-bounds access with distributed grids).
!>        in principle, something similar would be needed for the computation of the cube bounds
!>
!> \param cube_center ...
!> \param rs_desc ...
!> \param zeta ...
!> \param zetb ...
!> \param ra ...
!> \param rab ...
!> \par History
!>      11.2008 created [Joost VandeVondele]
! *****************************************************************************
  SUBROUTINE compute_cube_center(cube_center,rs_desc,zeta,zetb,ra,rab)

    INTEGER, DIMENSION(3), INTENT(OUT)       :: cube_center
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    REAL(KIND=dp), INTENT(IN)                :: zeta, zetb, ra(3), rab(3)

    REAL(KIND=dp)                            :: zetp
    REAL(KIND=dp), DIMENSION(3)              :: rp

    zetp = zeta + zetb
    rp(:) = ra(:) + zetb/zetp*rab(:)
    cube_center(:) = FLOOR(MATMUL(rs_desc%dh_inv,rp))

  END SUBROUTINE compute_cube_center

! *****************************************************************************
!> \brief ...
!> \param info ...
!> \param radius ...
!> \param lb ...
!> \param ub ...
!> \param rp ...
! *****************************************************************************
  SUBROUTINE return_cube_nonortho(info,radius,lb,ub,rp)

    TYPE(cube_info_type)                     :: info
    REAL(KIND=dp)                            :: radius
    INTEGER                                  :: lb(3), ub(3)
    REAL(KIND=dp)                            :: rp(3)

    CHARACTER(LEN=*), PARAMETER :: routineN = 'return_cube_nonortho', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k
    REAL(KIND=dp)                            :: point(3), res(3)

    IF (radius>info%max_rad_ga) THEN
            !
            ! This is an important check. If the required radius for mapping the density is different
            ! from the actual computed one, (significant) errors can occur.
            ! This error can invariably be fixed by improving the computation of maxradius
            ! in the call to init_cube_info
            !
            WRITE(*,*) info%max_rad_ga,radius
            CALL stop_program(routineN,moduleN,__LINE__,"Called with too large radius.")
    ENDIF

    ! get the min max indices that contain at least the cube that contains a sphere around rp of radius radius
    ! if the cell is very non-orthogonal this implies that many useless points are included
    ! this estimate can be improved (i.e. not box but sphere should be used)
    lb= HUGE(lb)
    ub=-HUGE(ub)
    DO i=-1,1
    DO j=-1,1
    DO k=-1,1
       point(1)=rp(1)+i*radius
       point(2)=rp(2)+j*radius
       point(3)=rp(3)+k*radius
       CALL matvec_3x3(res,info%dh_inv,point)
       lb=MIN(lb,FLOOR(res))
       ub=MAX(ub,CEILING(res))
    ENDDO
    ENDDO
    ENDDO

  END SUBROUTINE return_cube_nonortho

! *****************************************************************************
!> \brief ...
!> \param info ...
!> \param radius ...
!> \param lb_cube ...
!> \param ub_cube ...
!> \param sphere_bounds ...
! *****************************************************************************
  SUBROUTINE return_cube(info,radius,lb_cube,ub_cube,sphere_bounds)

    TYPE(cube_info_type)                     :: info
    REAL(KIND=dp)                            :: radius
    INTEGER                                  :: lb_cube(3), ub_cube(3)
    INTEGER, DIMENSION(:), POINTER           :: sphere_bounds

    CHARACTER(LEN=*), PARAMETER :: routineN = 'return_cube', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: imr

    IF (info%orthorhombic) THEN
        imr=MAX(1,CEILING(radius/info%drmin))
        IF (imr.gt.info%max_radius) THEN
            !
            ! This is an important check. If the required radius for mapping the density is different
            ! from the actual computed one, (significant) errors can occur.
            ! This error can invariably be fixed by improving the computation of maxradius
            ! in the call to init_cube_info
            !
            CALL stop_program(routineN,moduleN,__LINE__,"Called with too large radius.")
        ENDIF
        lb_cube(:)       = info%lb_cube(:,imr)
        ub_cube(:)       = info%ub_cube(:,imr)
        sphere_bounds => info%sphere_bounds(imr)%p
    ELSE
        ! nothing yet, we should check the radius
    ENDIF

  END SUBROUTINE return_cube

  ! this is the integer max radius of the cube
! *****************************************************************************
!> \brief ...
!> \param info ...
! *****************************************************************************
  INTEGER FUNCTION return_cube_max_iradius(info)
    TYPE(cube_info_type)                     :: info
    return_cube_max_iradius=info%max_radius
  END FUNCTION return_cube_max_iradius

! *****************************************************************************
!> \brief ...
!> \param info ...
! *****************************************************************************
  SUBROUTINE destroy_cube_info(info)
    TYPE(cube_info_type)                     :: info

    CHARACTER(LEN=*), PARAMETER :: routineN = 'destroy_cube_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, istat

    IF (info%orthorhombic) THEN
        DEALLOCATE (info%lb_cube,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "info%lb_cube")
        DEALLOCATE (info%ub_cube,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "info%ub_cube")
        DEALLOCATE (info%sphere_bounds_count,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "info%sphere_bounds_count")
        DO i=1,info%max_radius
          DEALLOCATE (info%sphere_bounds(i)%p,STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "info%sphere_bounds(i)%p")
        END DO
        DEALLOCATE (info%sphere_bounds,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "info%sphere_bounds")
    ELSE
        ! no info to be deallocated
    ENDIF
  END SUBROUTINE

! *****************************************************************************
!> \brief ...
!> \param info ...
!> \param dr ...
!> \param dh ...
!> \param dh_inv ...
!> \param ortho ...
!> \param max_radius ...
! *****************************************************************************
  SUBROUTINE init_cube_info(info,dr,dh,dh_inv,ortho,max_radius)
    TYPE(cube_info_type), INTENT(OUT)        :: info
    REAL(KIND=dp), INTENT(IN)                :: dr(3), dh(3,3), dh_inv(3,3)
    LOGICAL, INTENT(IN)                      :: ortho
    REAL(KIND=dp), INTENT(IN)                :: max_radius

    CHARACTER(LEN=*), PARAMETER :: routineN = 'init_cube_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: check_1, check_2, handle, i, &
                                                igmin, imr, istat, jg, jg2, &
                                                jgmin, k, kg, kg2, kgmin, &
                                                lb(3), ub(3)
    REAL(KIND=dp)                            :: drmin, dxi, dy2, dyi, dz2, &
                                                dzi, radius, radius2, rp(3)

    CALL timeset(routineN,handle)
    info%dr = dr
    info%dh = dh
    info%dh_inv = dh_inv
    info%orthorhombic = ortho
    info%max_rad_ga = max_radius
    drmin = MINVAL(dr)
    info%drmin=drmin

    NULLIFY(info%lb_cube,info%ub_cube,&
            info%sphere_bounds_count,info%sphere_bounds)

    IF ( .NOT. info%orthorhombic ) THEN

        rp=0.0_dp
        !
        ! could still be wrong (maybe needs an additional +1 to account for off-gridpoint rp's)
        !
        CALL return_cube_nonortho(info,max_radius,lb,ub,rp)
        info%max_radius=MAX(MAXVAL(ABS(lb)),MAXVAL(ABS(ub)))

    ELSE

        ! this info is specialized to orthogonal grids
        imr=CEILING((max_radius)/drmin)
        info%max_radius=imr
        dzi=1.0_dp/dr(3)
        dyi=1.0_dp/dr(2)
        dxi=1.0_dp/dr(1)
        dz2=(dr(3))**2
        dy2=(dr(2))**2

        ALLOCATE(info%lb_cube(3,imr),info%ub_cube(3,imr), &
                 info%sphere_bounds_count(imr),info%sphere_bounds(imr))
        check_1=0
        check_2=0
!       count and allocate

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
              check_1=MODULO((kgmin*97+jgmin*37+igmin*113)*check_1+1277,9343)
              k=k+1
             ENDDO
          ENDDO
          info%sphere_bounds_count(i)=k-1
          ALLOCATE (info%sphere_bounds(i)%p(info%sphere_bounds_count(i)),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "info%sphere_bounds(i)%p",&
                                           int_size*info%sphere_bounds_count(i))
        ENDDO

!       init sphere_bounds array
        ! notice : as many points in lb_cube..0 as 1..ub_cube
        DO i=1,imr
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
              check_2=MODULO((kgmin*97+jgmin*37+igmin*113)*check_2+1277,9343)
              info%lb_cube(1,i)=MIN(igmin,info%lb_cube(1,i))
              info%sphere_bounds(i)%p(k)=igmin
              k=k+1
             ENDDO
          ENDDO
          info%ub_cube(:,i)=1-info%lb_cube(:,i)
        ENDDO
        IF (check_1.NE.check_2) THEN
            CALL stop_program(routineN,moduleN,__LINE__,&
                              "Irreproducible fp math caused memory corruption")
        END IF

    END IF

    CALL timestop(handle)

  END SUBROUTINE

  ! try to hide things from the optimizer, so that we get the same numbers,
  ! always (this solves the optimisation problems with the intel and nag compiler
  ! in which the counting loops and execution loops above are executed a different
  ! number of times, even at -O1
! *****************************************************************************
!> \brief ...
!> \param prefactor ...
!> \param i ...
!> \param drmin ...
!> \param dz2 ...
!> \param dy2 ...
!> \param kg2 ...
!> \param jg2 ...
!> \retval res ...
! *****************************************************************************
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

! *****************************************************************************
!> \brief ...
!> \param buf ...
!> \param i ...
!> \param jg2 ...
!> \param kg2 ...
!> \retval res ...
! *****************************************************************************
  FUNCTION do_and_hide_it_2(buf,i,jg2,kg2) RESULT(res)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: buf
    INTEGER                                  :: i, jg2, kg2, res

    buf(2)=(i*buf(2))**2
    res = CEILING(-0.1E-7_dp-buf(1)*SQRT(MAX(buf(2)-kg2*buf(3)-jg2*buf(4),0.0_dp)))
  END FUNCTION do_and_hide_it_2

END MODULE cube_utils

