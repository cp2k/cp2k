!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!>   \brief
!>     Routines to efficently collocate and integrate gaussians on a grid
!>     These use most of Joost's tricks and a couple more...
!>     result is *speed* and genericity
!>   \author Fawzi Mohamed, 2007
!>   \notes original available with BSD style license
! *****************************************************************************
MODULE gauss_colloc
  USE d3_poly
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE lgrid_types,                     ONLY: lgrid_type

  !$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
#include "./common/cp_common_uses.f90"

IMPLICIT NONE
PRIVATE

  PUBLIC :: collocGauss,&
            integrateGaussFull

#ifdef FD_DEBUG
#define IF_CHECK(x,y) x
#else
#define IF_CHECK(x,y) y
#endif

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gauss_colloc'

    REAL(dp), PARAMETER :: small=TINY(1.0_dp)

  ! keep prettify happy (does not see the include)
  INTEGER(KIND=int_8), PARAMETER           :: unused_import_of_int_8=1


CONTAINS

!**** radius calculation

! *****************************************************************************
!> \brief finds square radius of a gaussian time a polynomial with a gaussian
!> in poly_shift (wrt. the reference system of the polynomial)
!>
!> utility method
!> \param poly ...
!> \param poly_shift ...
!> \param alpha ...
!> \param epsilon ...
!> \param error ...
!> \param scale ...
!> \retval res ...
! *****************************************************************************
FUNCTION calc_max_r2(poly,poly_shift,alpha,epsilon,error,scale) RESULT(res)
    REAL(dp), DIMENSION(0:), INTENT(in)      :: poly
    REAL(dp), DIMENSION(3)                   :: poly_shift
    REAL(dp), INTENT(in)                     :: alpha, epsilon
    TYPE(cp_error_type), INTENT(inout)       :: error
    REAL(dp), INTENT(in), OPTIONAL           :: scale
    REAL(dp)                                 :: res

    CHARACTER(len=*), PARAMETER :: routineN = 'calc_max_r2', &
      routineP = moduleN//':'//routineN
    REAL(dp), DIMENSION(3, 3), PARAMETER :: &
      identity_m = RESHAPE((/ 1,0,0, 0,1,0, 0,0,1 /),(/ 3,3 /))

    INTEGER                                  :: grad, stat
    REAL(dp)                                 :: eps
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: poly_ijk

    eps=epsilon
    IF (PRESENT(scale)) THEN
        IF (scale==0.0_dp) THEN
            res=0.0_dp
            RETURN
        END IF
        eps=eps/scale
    END IF
    grad=grad_size3(SIZE(poly))
    ALLOCATE(poly_ijk(0:poly_size3(grad)-1),stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
    CALL poly_affine_t3(poly,identity_m,poly_shift,poly_ijk,error=error)
    res=calcRadius2(poly_ijk,alpha,eps)
    DEALLOCATE(poly_ijk,stat=stat)
    CPPostconditionNoFail(stat==0,cp_fatal_level,routineP,error)
END FUNCTION

! *****************************************************************************
!> \brief calculates the radius of polynomial times a gaussian in 0,0,0
!>
!> performace of bounding
!> random poly (grad 0-6), random h (det(h)/(a*b*c)>0.2) at least 1 grid point per unit, alpha (0.5-2.5)
!>
!> actual choice: igrad-norm, pade(1,2), 1 step iterative solver
!>  Nr values:  15996
!>  InnerBound too Low: 96 0.600150037509 % min value: 9.25569425556e-12 mean 0.000749140381379
!>  OuterBound too High: 0 0.0 % max value: 9.99925035941e-11 mean 3.18900144269e-11
!>
!> runner up: 2-norm, pade(1,2), 1 step iterative solver
!>  Nr values:  15996
!>  InnerBound too Low: 197 1.23155788947 % min value: 7.29606466146e-12 mean 0.000749166756885
!>  OuterBound too High: 0 0.0 % max value: 9.99925035941e-11 mean 2.72429081183e-11
!>
!> \param poly ...
!> \param alpha ...
!> \param epsilon ...
!> \retval res ...
!> \notes can bound either the value or the (absolute) surface integral
! *****************************************************************************
FUNCTION calcRadius2(poly,alpha,epsilon) RESULT(res)
    REAL(dp), DIMENSION(0:), INTENT(in)      :: poly
    REAL(dp), INTENT(in)                     :: alpha, epsilon
    REAL(dp)                                 :: res

    INTEGER                                  :: from_i, grad, gradSize, i, &
                                                igrad, p_size, to_i
    REAL(dp)                                 :: coeffAtt, rAtt, residual, t
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: rAtts

    res=0.0_dp
    grad=grad_size3(SIZE(poly))
    ALLOCATE(rAtts(0:grad))
    from_i=0
    gradSize=1
    p_size=SIZE(poly)
    DO igrad=0,grad
        to_i=MIN(from_i+gradSize,p_size)
        !! igrad-norm
        IF (igrad==0) THEN
            coeffAtt=ABS(poly(from_i))
        ELSE
            coeffAtt=0.0_dp
            DO i=from_i,to_i-1
                coeffAtt=coeffAtt+ABS(poly(i))**igrad
            END DO
            coeffAtt=coeffAtt**(1.0_dp/REAL(igrad))
        END IF
        !! infinity-norm
        !coeffAtt=maxval(abs(poly(from_i:to_i-1)))
        !! 2-norm (euclidean)
        !coeffAtt=sqrt(sum(poly(from_i:to_i-1)**2))
        rAtt=sphRad2(coeffAtt,igrad,alpha,epsilon) ! bounds value
        !rAtt=sphRad2(4*pi*coeffAtt,igrad+2,alpha,epsilon) ! bounds surface integral
        rAtts(igrad)=rAtt
        res=MAX(res,rAtts(igrad))
        gradSize=gradSize+(igrad+2)
        from_i=to_i
    END DO
    residual=0.0_dp
    DO igrad=0,grad
        ! pade expansion 1,2 of exp(-x) at x=0
        t=alpha*(res-rAtts(igrad))
        residual=residual+1.0_dp/(1.0_dp+t+0.5_dp*t**2)
    END DO
    res=res+LOG(residual)/alpha
    DEALLOCATE(rAtts)
END FUNCTION

! *****************************************************************************
!> \brief finds the largest real solution for r^2 of the equation coeff*r^l*exp(-alpha*r^2)
!> \param coeff ...
!> \param l ...
!> \param alpha ...
!> \param epsilon ...
!> \retval res ...
! *****************************************************************************
FUNCTION sphRad2(coeff,l,alpha,epsilon) RESULT(res)
    REAL(dp), INTENT(in)                     :: coeff
    INTEGER, INTENT(in)                      :: l
    REAL(dp), INTENT(in)                     :: alpha, epsilon
    REAL(dp)                                 :: res

    REAL(dp)                                 :: r1, t, x, xx

    IF (ABS(coeff)<small) THEN
        res=0.0_dp
        RETURN
    END IF
    xx=MAX(0.0_dp,-LOG(epsilon/coeff)/alpha)
    x=MAX(1.0_dp,xx)
    t=0.5_dp*REAL(l,dp)/alpha

    ! rational functions that approximate the solution of r=x+t*log(r)
    r1=MAX(x,(5.8679922317012645_dp - 13.433907351015852_dp*t + t**2 + 2.7189706320638596_dp*x)&
        /(2.1224230654286007_dp + 1.5855089862522989_dp*t + 0.05338142797029377_dp*t**2 &
         + 1.4761929491905326_dp*x - 0.059575647341275406_dp*t*x + 0.031067088387632627_dp*x**2)&
        +(-0.43804767398212585_dp - 0.49900615723405245_dp*t + 0.20629903250691906_dp*t**2&
         + 0.9017510689333377_dp*x - 0.07046283868925775_dp*t*x - 0.1743645447869195_dp*x**2)&
        /(-3.617044072726269_dp - 0.7662306349036399_dp*t - 0.7614179211580372_dp*t**2&
         + 0.17934761179419514_dp*t**3 + 6.788069591421069_dp*x + 0.5260819714537482_dp*t*x&
         + 1.4780624656229489_dp*t**2*x - 3.05963073688407_dp*x**2 - 1.5521580811226825_dp*t*x**2&
         + x**3) - 0.6840379938652372_dp*MAX(0.0_dp,0.4_dp - t + (1.0_dp - x)/3.0_dp)&
        - 0.16927886252742727_dp*MAX(0.0_dp,1.0_dp - t + (1.0_dp - x)/3.0_dp) &
        +(-2.1054595688424427_dp + 3.0262780685014503_dp*t + 0.1353476536688748_dp*t**2 + x&
         + 0.03244405654399495_dp*t*x + 0.0017806182354747218_dp*x**2)&
        /( 0.9552112747891883_dp + 0.016045994433298932_dp*t + 0.0018347038147891919_dp*x&
         - 1.8296261802177613e-6_dp*MIN(1000.0_dp,t)**2))
    res=MAX(0.0_dp,xx+t*LOG(r1)) ! one step of the iterative solver, relative error less than 0.004
END FUNCTION

! *****************************************************************************
!> \brief finds a box on the grid that contains a sphere of the given radius.
!> If guarantee_nearest is true at least the nearest points on the grid are
!> included by enlarging the radius
!> \param h ...
!> \param h_inv ...
!> \param posi ...
!> \param max_r2 ...
!> \param periodic ...
!> \param gdim ...
!> \param error ...
!> \param guarantee_nearest ...
!> \retval res ...
! *****************************************************************************
FUNCTION calcBox(h,h_inv,posi,max_r2,&
        periodic,gdim,error,guarantee_nearest) RESULT(res)
    REAL(dp), DIMENSION(0:2, 0:2), &
      INTENT(in)                             :: h, h_inv
    REAL(dp), DIMENSION(0:2), INTENT(in)     :: posi
    REAL(dp), INTENT(in)                     :: max_r2
    INTEGER, DIMENSION(0:2), INTENT(in)      :: periodic, gdim
    TYPE(cp_error_type), INTENT(inout)       :: error
    LOGICAL, INTENT(in), OPTIONAL            :: guarantee_nearest
    INTEGER                                  :: res(2,0:2)

    CHARACTER(len=*), PARAMETER :: routineN = 'calcBox', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, imax, imin, j, jmax, jmin, &
                                                k, kmax, kmin
    INTEGER, DIMENSION(0:2)                  :: cellShift, fullShift, ndim, &
                                                shiftPos
    LOGICAL                                  :: g_nearest
    REAL(dp) :: cci0, cci1, cci2, ccj0, ccj1, ccj2, cck0, cck1, cck2, &
      delta_i, delta_j, delta_k, m(0:2,0:2), maxr2, r_0, scaled_h(0:2,0:2), &
      sqDi, sqDj, sqDk
    REAL(dp), DIMENSION(0:2)                 :: l, normD, resPos, riPos, &
                                                rpos, wrPos

    g_nearest=.TRUE.
    IF (PRESENT(guarantee_nearest)) g_nearest=guarantee_nearest
    ndim=gdim
    rPos=0.0_dp
    DO j=0,2
        DO i=0,2
            rPos(i)=rPos(i)+h_inv(i,j)*posi(j)
        END DO
    END DO
    cellShift=FLOOR(rPos)*periodic
    wrPos=rPos-REAL(cellShift,dp)
    riPos=wrPos*ndim
    shiftPos=FLOOR(riPos+0.5_dp)
    fullShift=shiftPos+ndim*cellShift
    resPos=riPos-shiftPos
    normD=1.0_dp/REAL(ndim,dp)

    IF (g_nearest) THEN
        maxr2=0.0_dp
        DO j=0,2
            DO i=0,2
                maxr2=maxr2+(h(i,j)*normD(j))**2
            END DO
        END DO
        maxr2=MAX(max_r2,maxr2)
    ELSE
        maxr2=max_r2
    END IF

    scaled_h=0.0_dp
    DO j=0,2
        DO i=0,2
            scaled_h(i,2-j)=scaled_h(i,2-j)+h(i,j)*normD(j)
        END DO
    END DO

    ! build up quadratic form (ellipsoid)
    m=0.0_dp
    DO j=0,2
        DO i=0,2
            DO k=0,2
                m(i,j)=m(i,j)+scaled_h(k,i)*scaled_h(k,j)
            END DO
        END DO
    END DO

    l=0.0_dp
    DO j=0,2
        DO i=0,2
            l(j)=l(j)-2.0*resPos(2-i)*m(i,j)
        END DO
    END DO

    r_0=0.0_dp
    DO i=0,2
        r_0=r_0-0.5*resPos(2-i)*l(i)
    END DO

    ! i box bounds
    cci2 = (m(2,2) * m(0,0) * m(1,1) - m(2,2) * m(0,1) ** 2 - m(1,1) * m(0,2) ** 2 &
        + 2.0_dp * m(0,2) * m(0,1) * m(1,2) - m(0,0) * m(1,2) ** 2) / (m(2,2) * m(1,1) - m(1,2) ** 2)
    cci1 = -(-m(2,2) * l(0) * m(1,1) + m(2,2) * m(0,1) * l(1) + l(2) * m(0,2) * m(1,1) &
        + l(0) * m(1,2) ** 2 - l(2) * m(0,1) * m(1,2) - m(1,2) * l(1) * m(0,2)) / (m(2,2) * m(1,1) - m(1,2) ** 2)
    cci0 = -((-4.0_dp * m(2,2) * r_0 * m(1,1) + m(2,2) * l(1) ** 2 + l(2) ** 2 * m(1,1) &
        - 2.0_dp * l(1) * m(1,2) * l(2) + 4.0_dp * r_0 * m(1,2) ** 2) &
        / (m(2,2) * m(1,1) - m(1,2) ** 2)) / 4.0_dp-maxr2
    delta_i=cci1*cci1-4.0_dp*cci2*cci0
    IF (delta_i<0.0_dp) THEN
        imin=fullShift(2)+CEILING((-cci1)/(2.0_dp*cci2))
        imax=imin-1
    ELSE
        sqDi=SQRT(delta_i)
        imin=fullShift(2)+CEILING((-cci1-sqDi)/(2.0_dp*cci2))
        imax=fullShift(2)+FLOOR((-cci1+sqDi)/(2.0_dp*cci2))
    END IF
    IF (periodic(2)==0) THEN
        imin=MAX(0,imin)
        imax=MIN(imax,ndim(2)-1)
    END IF
    res(1,2)=imin
    res(2,2)=imax

    ! j box bounds
    ccj2 = (m(0,0) * m(2,2) * m(1,1) - m(0,0) * m(1,2) ** 2 - m(0,1) ** 2 * m(2,2) &
        + 2.0_dp * m(0,1) * m(0,2) * m(1,2) - m(1,1) * m(0,2) ** 2) &
        / (m(0,0) * m(2,2) - m(0,2) ** 2)
    ccj1 = -(-m(0,0) * l(1) * m(2,2) + m(0,0) * l(2) * m(1,2) + l(0) * m(0,1) * m(2,2) &
        - m(0,2) * l(2) * m(0,1) - l(0) * m(0,2) * m(1,2) + l(1) * m(0,2) ** 2) &
        / (m(0,0) * m(2,2) - m(0,2) ** 2)
    ccj0 = (4.0_dp * m(0,0) * m(2,2) * r_0 - m(0,0) * l(2) ** 2 - m(2,2) * l(0) ** 2 &
        + 2.0_dp * m(0,2) * l(2) * l(0) - 4.0_dp * r_0 * m(0,2) ** 2) &
        / (m(0,0) * m(2,2) - m(0,2) ** 2) / 4.0_dp-maxr2
    delta_j=ccj1*ccj1-4.0_dp*ccj2*ccj0
    IF (delta_j<0.0_dp) THEN
        jmin=fullShift(1)+CEILING((-ccj1)/(2.0_dp*ccj2))
        jmax=jmin-1
    ELSE
        sqDj=SQRT(delta_j)
        jmin=fullShift(1)+CEILING((-ccj1-sqDj)/(2.0_dp*ccj2))
        jmax=fullShift(1)+FLOOR((-ccj1+sqDj)/(2.0_dp*ccj2))
    END IF
    IF (periodic(1)==0) THEN
        jmin=MAX(0,jmin)
        jmax=MIN(jmax,ndim(1)-1)
    END IF
    res(1,1)=jmin
    res(2,1)=jmax

    ! k box bounds
    cck2 = (m(0,0) * m(2,2) * m(1,1) - m(0,0) * m(1,2) ** 2 - m(0,1) ** 2 * m(2,2) &
        + 2.0_dp * m(0,1) * m(0,2) * m(1,2) - m(1,1) * m(0,2) ** 2) / (m(0,0) * m(1,1) - m(0,1) ** 2)
    cck1 = (m(0,0) * l(2) * m(1,1) - m(0,0) * m(1,2) * l(1) - m(0,2) * l(0) * m(1,1) &
        + l(0) * m(0,1) * m(1,2) - l(2) * m(0,1) ** 2 + m(0,1) * l(1) * m(0,2)) / (m(0,0) * m(1,1) - m(0,1) ** 2)
    cck0 = (4.0_dp * m(0,0) * m(1,1) * r_0 - m(0,0) * l(1) ** 2 - m(1,1) * l(0) ** 2 &
        + 2.0_dp * m(0,1) * l(1) * l(0) - 4.0_dp * r_0 * m(0,1) ** 2) &
        / (m(0,0) * m(1,1) - m(0,1) ** 2) / 4.0_dp-maxr2
    delta_k=cck1*cck1-4.0_dp*cck2*cck0
    IF (delta_k<0.0_dp) THEN
        kmin=fullShift(0)+CEILING((-cck1)/(2.0_dp*cck2))
        kmax=kmin-1
    ELSE
        sqDk=SQRT(delta_k)
        kmin=fullShift(0)+CEILING((-cck1-sqDk)/(2.0_dp*cck2))
        kmax=fullShift(0)+FLOOR((-cck1+sqDk)/(2.0_dp*cck2))
    END IF
    IF (periodic(0)==0) THEN
        kmin=MAX(0,kmin)
        kmax=MIN(kmax,ndim(0)-1)
    END IF
    res(1,0)=kmin
    res(2,0)=kmax

END FUNCTION

! *****************************************************************************
!> \brief collocate a periodically repeated gaussian on a non orthormbic grid
!>
!> this routine has been tested and works well with cells with
!> det(h)/sqrt(tr(dot(h^T,h)))>0.2 (2 angles bigger than 24 deg or one angle
!> bigger than 11 deg).
!> Because of its numerics it might fail badly (infinity or NaN) with
!> with more deformed cells. Avoiding this would be bossible only using
!> IEEE numerics controls, which would also make everything slower and
!> less supported.
!> Still the actual numeric has been carefully tuned, and in normal cases
!> and most anormal it should work.
!> With det(h)/sqrt(tr(dot(h^T,h)))>0.2 I could not find any failure.
!>
!> \param h cell matrix
!> \param h_inv inverse of the cell matrix
!> \param grid the grid
!> \param poly polynomial (d3_poly format)
!> \param alphai exponential coeff
!> \param posi position of the gaussian
!> \param max_r2 maximum radius of collocation squared
!> \param periodic array of 0 or 1 that says which dimensions have pbc (1=pbc)
!> \param gdim dimension of the grid (grid might be a subset)
!> \param local_bounds local bounds of the grid piece that is kept locally
!>   (i.e. of grid) the global grid is assumed to atart at 0,0,0
!> \param local_shift start indexes of the local slice (i.e. of grid)
!> \param poly_shift position of posi in the polynomial reference system.
!>  Set it to posi to use the global reference system.
!> \param scale a global scale factor
!> \param lgrid ...
!> \param error type to control the error handling
! *****************************************************************************
SUBROUTINE collocGauss(h,h_inv,grid,poly,alphai,posi,max_r2,&
        periodic,gdim,local_bounds,local_shift,poly_shift,scale,lgrid,error)
    REAL(dp), DIMENSION(0:2, 0:2), &
      INTENT(in)                             :: h, h_inv
    REAL(dp), DIMENSION(0:, 0:, 0:), &
      INTENT(inout)                          :: grid
    REAL(dp), DIMENSION(:), INTENT(in)       :: poly
    REAL(dp), INTENT(in)                     :: alphai
    REAL(dp), DIMENSION(0:2), INTENT(in)     :: posi
    REAL(dp), INTENT(in)                     :: max_r2
    INTEGER, DIMENSION(0:2), INTENT(in)      :: periodic
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: gdim
    INTEGER, DIMENSION(2, 0:2), INTENT(in), &
      OPTIONAL                               :: local_bounds
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: local_shift
    REAL(dp), DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: poly_shift
    REAL(dp), INTENT(in), OPTIONAL           :: scale
    TYPE(lgrid_type), INTENT(inout), &
      OPTIONAL                               :: lgrid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'collocGauss', &
      routineP = moduleN//':'//routineN


#include "colloc_int_body.f90"

END SUBROUTINE

! *****************************************************************************
!> \brief collocate a periodically repeated gaussian on a non orthormbic grid
!>
!> like collocGauss, but takes a flattened grid as input.
!> A little bit less optimized, mostly to connect to other languages.
!> \param h ...
!> \param h_inv ...
!> \param grid ...
!> \param ngpts ...
!> \param ldim ...
!> \param poly ...
!> \param alphai ...
!> \param posi ...
!> \param max_r2 ...
!> \param periodic ...
!> \param gdim ...
!> \param local_bounds ...
!> \param local_shift ...
!> \param poly_shift ...
!> \param scale ...
!> \param error ...
! *****************************************************************************
SUBROUTINE collocGaussFlat(h,h_inv,grid,ngpts,ldim,poly,alphai,posi,max_r2,&
        periodic,gdim,local_bounds,local_shift,poly_shift,scale,error)
    REAL(dp), DIMENSION(0:2, 0:2), &
      INTENT(in)                             :: h, h_inv
    INTEGER, INTENT(in) :: ngpts,ldim(0:2)
    REAL(dp), DIMENSION(IF_CHECK(ngpts,*)), &
      INTENT(inout)                          :: grid
    REAL(dp), DIMENSION(:), INTENT(in)       :: poly
    REAL(dp), INTENT(in)                     :: alphai
    REAL(dp), DIMENSION(0:2), INTENT(in)     :: posi
    REAL(dp), INTENT(in)                     :: max_r2
    INTEGER, DIMENSION(0:2), INTENT(in)      :: periodic
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: gdim
    INTEGER, DIMENSION(2, 0:2), INTENT(in), &
      OPTIONAL                               :: local_bounds
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: local_shift
    REAL(dp), DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: poly_shift
    REAL(dp), INTENT(in), OPTIONAL           :: scale
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'collocGaussFlat', &
        routineP=moduleN//':'//routineN
#define FM_FLAT_GRID
#include "colloc_int_body.f90"
#undef FM_FLAT_GRID
END SUBROUTINE

! *****************************************************************************
!> \brief integrates a gaussian times a polynomial.
!>
!> Most things are the same as for collocGauss (see its comments).
!> Unlike collocGauss this function can receive more than one
!> polynomial to integrate at once.
!> How many polynomial are passed is controlled with npoly.
!> The polynomials must be in the d3_poly format.
!> res will then contain the corresponding integrals
!> \param h ...
!> \param h_inv ...
!> \param grid ...
!> \param poly ...
!> \param alphai ...
!> \param posi ...
!> \param max_r2 ...
!> \param periodic ...
!> \param npoly ...
!> \param res ...
!> \param gdim ...
!> \param local_bounds ...
!> \param local_shift ...
!> \param poly_shift ...
!> \param scale ...
!> \param error ...
! *****************************************************************************
SUBROUTINE integrateGauss(h,h_inv,grid,poly,alphai,posi,max_r2,&
        periodic,npoly,res,gdim,local_bounds,local_shift,poly_shift,scale,error)
    INTEGER, INTENT(in) :: npoly
    REAL(dp), DIMENSION(npoly), INTENT(out) :: res
    REAL(dp), DIMENSION(0:2, 0:2), &
      INTENT(in)                             :: h, h_inv
    REAL(dp), DIMENSION(0:, 0:, 0:), &
      INTENT(inout)                          :: grid
    REAL(dp), DIMENSION(:), INTENT(in)       :: poly
    REAL(dp), INTENT(in)                     :: alphai
    REAL(dp), DIMENSION(0:2), INTENT(in)     :: posi
    REAL(dp), INTENT(in)                     :: max_r2
    INTEGER, DIMENSION(0:2), INTENT(in)      :: periodic
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: gdim
    INTEGER, DIMENSION(2, 0:2), INTENT(in), &
      OPTIONAL                               :: local_bounds
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: local_shift
    REAL(dp), DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: poly_shift
    REAL(dp), INTENT(in), OPTIONAL           :: scale
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'integrateGauss', &
        routineP=moduleN//':'//routineN

#define FMG_INTEGRATE
#include "colloc_int_body.f90"
#undef FMG_INTEGRATE
END SUBROUTINE

! *****************************************************************************
!> \brief integrates a gaussian times any polynomial up to a give order.
!>
!> Most things are the same as for collocGauss (see its comments).
!> Returns the integrals of all the monomials in d3 format into poly
!> \param h ...
!> \param h_inv ...
!> \param grid ...
!> \param poly ...
!> \param alphai ...
!> \param posi ...
!> \param max_r2 ...
!> \param periodic ...
!> \param gdim ...
!> \param local_bounds ...
!> \param local_shift ...
!> \param poly_shift ...
!> \param scale ...
!> \param error ...
! *****************************************************************************
SUBROUTINE integrateGaussFull(h,h_inv,grid,poly,alphai,posi,max_r2,&
        periodic,gdim,local_bounds,local_shift,poly_shift,scale,error)
    REAL(dp), DIMENSION(0:2, 0:2), &
      INTENT(in)                             :: h, h_inv
    REAL(dp), DIMENSION(0:, 0:, 0:), &
      INTENT(inout)                          :: grid
    REAL(dp), DIMENSION(:), INTENT(out)      :: poly
    REAL(dp), INTENT(in)                     :: alphai
    REAL(dp), DIMENSION(0:2), INTENT(in)     :: posi
    REAL(dp), INTENT(in)                     :: max_r2
    INTEGER, DIMENSION(0:2), INTENT(in)      :: periodic
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: gdim
    INTEGER, DIMENSION(2, 0:2), INTENT(in), &
      OPTIONAL                               :: local_bounds
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: local_shift
    REAL(dp), DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: poly_shift
    REAL(dp), INTENT(in), OPTIONAL           :: scale
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'integrateGaussFull', &
        routineP=moduleN//':'//routineN

#define FMG_INTEGRATE_FULL
#include "colloc_int_body.f90"
#undef FMG_INTEGRATE_FULL
END SUBROUTINE

! *****************************************************************************
!> \brief integrates a gaussian times any polynomial up to a give order.
!>
!> like integrateGaussFull, but with a flat grid.
!> a little bit less optimized, mostly to connect to other languages
!> \param h ...
!> \param h_inv ...
!> \param grid ...
!> \param ngpts ...
!> \param ldim ...
!> \param poly ...
!> \param alphai ...
!> \param posi ...
!> \param max_r2 ...
!> \param periodic ...
!> \param gdim ...
!> \param local_bounds ...
!> \param local_shift ...
!> \param poly_shift ...
!> \param scale ...
!> \param error ...
! *****************************************************************************
SUBROUTINE integrateGaussFullFlat(h,h_inv,grid,ngpts,ldim,poly,alphai,posi,max_r2,&
        periodic,gdim,local_bounds,local_shift,poly_shift,scale,error)
    REAL(dp), DIMENSION(0:2, 0:2), &
      INTENT(in)                             :: h, h_inv
    INTEGER, INTENT(in) :: ngpts,ldim(0:2)
    REAL(dp), DIMENSION(IF_CHECK(ngpts,*)), &
      INTENT(inout)                          :: grid
    REAL(dp), DIMENSION(:), INTENT(out)      :: poly
    REAL(dp), INTENT(in)                     :: alphai
    REAL(dp), DIMENSION(0:2), INTENT(in)     :: posi
    REAL(dp), INTENT(in)                     :: max_r2
    INTEGER, DIMENSION(0:2), INTENT(in)      :: periodic
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: gdim
    INTEGER, DIMENSION(2, 0:2), INTENT(in), &
      OPTIONAL                               :: local_bounds
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: local_shift
    REAL(dp), DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: poly_shift
    REAL(dp), INTENT(in), OPTIONAL           :: scale
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'integrateGaussFullFlat', &
        routineP=moduleN//':'//routineN

#define FM_FLAT_GRID
#define FMG_INTEGRATE_FULL
#include "colloc_int_body.f90"
#undef FMG_INTEGRATE_FULL
#undef FM_FLAT_GRID
END SUBROUTINE

! *****************************************************************************
!> \brief collocate a periodically repeated gaussian on a non orthormbic grid (reference function)
!> \param h ...
!> \param h_inv ...
!> \param grid ...
!> \param poly ...
!> \param alphai ...
!> \param posi ...
!> \param max_r2 ...
!> \param periodic ...
!> \param gdim ...
!> \param local_bounds ...
!> \param local_shift ...
!> \param error ...
! *****************************************************************************
SUBROUTINE collocGauss_safe(h,h_inv,grid,poly,alphai,posi,max_r2,&
        periodic,gdim,local_bounds,local_shift,error)
    REAL(dp), DIMENSION(0:2, 0:2), &
      INTENT(in)                             :: h, h_inv
    REAL(dp), DIMENSION(0:, 0:, 0:), &
      INTENT(inout)                          :: grid
    REAL(dp), DIMENSION(:), INTENT(in)       :: poly
    REAL(dp), INTENT(in)                     :: alphai
    REAL(dp), DIMENSION(0:2), INTENT(in)     :: posi
    REAL(dp), INTENT(in)                     :: max_r2
    INTEGER, DIMENSION(0:2), INTENT(in)      :: periodic
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: gdim
    INTEGER, DIMENSION(2, 0:2), INTENT(in), &
      OPTIONAL                               :: local_bounds
    INTEGER, DIMENSION(0:2), INTENT(in), &
      OPTIONAL                               :: local_shift
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'collocGauss_safe', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bounds(0:1,0:2), i, idir, ii, &
                                                ij, ik, j, k, &
                                                l_bounds(0:1,0:2)
    INTEGER, DIMENSION(0:2)                  :: cellShift, l_shift, lb, ndim, &
                                                shiftPos, ub
    REAL(dp)                                 :: f, max_miss, maxr2, maxr2_0, &
                                                maxr2_2, min_miss, p_v(1), &
                                                r_2, scaled_h(0:2,0:2)
    REAL(dp), DIMENSION(0:2)                 :: normD, r, resPosReal, ri, &
                                                riPos, rj, rPos, wrPos

    IF (PRESENT(gdim)) THEN
        ndim=gdim
    ELSE
        ndim(0)=SIZE(grid,1)
        ndim(1)=SIZE(grid,2)
        ndim(2)=SIZE(grid,3)
    END IF
    IF (PRESENT(local_bounds)) THEN
        l_bounds=local_bounds
    ELSE
        l_bounds(0,:)=0
        l_bounds(1,:)=ndim-1
    END IF
    IF (PRESENT(local_shift)) THEN
        l_shift=local_shift
    ELSE
        l_shift=0 ! l_bounds(0,:)
    END IF
    lb=l_bounds(0,:)
    ub=l_bounds(1,:)-l_bounds(0,:)

    rPos=0
    DO j=0,2
        DO i=0,2
            rPos(i)=rPos(i)+h_inv(i,j)*posi(j)
        END DO
    END DO
    cellShift=FLOOR(rPos)*periodic
    wrPos=rPos-cellShift
    riPos=wrPos*ndim
    shiftPos=FLOOR(riPos+0.5_dp)
    normD=1.0_dp/REAL(ndim,dp)
    resPosReal=0.0_dp
    DO j=0,2
        DO i=0,2
            resPosReal(i)=resPosReal(i)+h(i,j)*(wrPos(j)-normD(j)*REAL(shiftPos(j),dp))
        END DO
    END DO

    IF (ALL(poly==0.0_dp)) RETURN ! remove if derivs

    maxr2=0.0_dp
    DO j=0,2
        DO i=0,2
            maxr2=maxr2+(h(i,j)*normD(j))**2
        END DO
    END DO
    maxr2_2=maxr2+0.2_dp/alphai
    maxr2_0=0.2_dp/alphai
    maxr2=MAX(max_r2,maxr2)
    maxr2_2=maxr2+maxr2_2+2.0_dp*SQRT(maxr2_2*maxr2)
    maxr2_0=maxr2+maxr2_0-2.0_dp*SQRT(maxr2_0*maxr2)

    DO j=0,2
        DO i=0,2
            scaled_h(i,j)=h(i,j)*normD(j)
        END DO
    END DO

    ! sqSize=1+floor(ndim*sqrt(maxr2/minimum.reduce(numpy.linalg.eig(dot(transpose(h),h))[0]))).astype(int)
    bounds=calcBox(h,h_inv,posi,max_r2,periodic,gdim,error)
    bounds(0,:)=bounds(0,:)-2
    bounds(1,:)=bounds(1,:)+2
    DO i=0,2
        IF (periodic(i)==0) THEN
            bounds(0,i)=MAX(l_bounds(0,i),bounds(0,i))
            bounds(1,i)=MIN(l_bounds(1,i),bounds(1,i))
        END IF
    END DO
    p_v=0
    max_miss=-1.0_dp
    min_miss=-1.0_dp
    DO i=bounds(0,0),bounds(1,0)
        ii=MODULO(i-lb(0),ndim(0))
        IF (ii>ub(0)) CYCLE
        ii=ii+l_shift(0)
        DO idir=0,2
            ri(idir)=scaled_h(idir,0)*i-posi(idir)
        END DO
        DO j=bounds(0,1),bounds(1,1)
            ij=MODULO(j-lb(1),ndim(1))
            IF (ij>ub(1)) CYCLE
            ij=ij+l_shift(1)
            DO idir=0,2
                rj(idir)=scaled_h(idir,1)*j+ri(idir)
            END DO
            DO k=bounds(0,2),bounds(1,2)
                ik=MODULO(k-lb(2),ndim(2))
                IF (ik>ub(2)) CYCLE
                ik=ik+l_shift(2)
                r_2=0.0_dp
                DO idir=0,2
                    r(idir)=scaled_h(idir,2)*k+rj(idir)
                    r_2=r_2+r(idir)**2
                END DO
                IF (r_2<maxr2_2) THEN
                    f=EXP(-alphai*r_2)
                    CALL poly_eval3(poly,r(0)+posi(0),r(1)+posi(1),r(2)+posi(2),p_v,error=error)
                    IF (r_2<=maxr2) THEN
                        IF (r_2>maxr2_0) min_miss=MAX(min_miss,ABS(p_v(1)*f))
                        grid(ii,ij,ik)=grid(ii,ij,ik)+p_v(1)*f
                    ELSE
                        max_miss=MAX(max_miss,ABS(p_v(1)*f))
                    END IF
                END IF
            END DO
        END DO
    END DO
    PRINT *,'miss2',min_miss,max_miss
END SUBROUTINE

END MODULE gauss_colloc
