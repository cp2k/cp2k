!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  numerical operations on real-space grid
!> \par History
!>       12.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE rs_methods

  USE kinds,                           ONLY: dp
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_methods,                      ONLY: pw_integrate_function,&
                                             pw_scale,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_type
  USE realspace_grid_types,            ONLY: &
       pw2rs, realspace_grid_desc_type, realspace_grid_p_type, &
       realspace_grid_type, rs2pw, rs_grid_create, rs_grid_release, &
       rs_grid_zero, rs_pw_transfer
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rs_methods'

  PUBLIC derive_fdm_cd3,&
         derive_fdm_cd5,&
         derive_fdm_cd7,&
         setup_grid_axes,&
         pw_mollifier

  REAL(dp), PARAMETER, PRIVATE         :: small_value = 1.0E-8_dp

CONTAINS 

! *****************************************************************************
!> \brief    2nd order finite difference derivative of a function on realspace grid
!> \param f  input function
!> \param df derivative of f
!> \param rs_grid real-space grid
!> \par History:
!>      - Creation (15.11.2013,MK)
!>      - Refactored and moved here from qs_sccs.F (12.2014, Hossein Bani-Hashemian)
!> \author     Matthias Krack (MK)
!> \version    1.0
! *****************************************************************************
  SUBROUTINE derive_fdm_cd3(f,df,rs_grid)

    TYPE(pw_type), POINTER                   :: f
    TYPE(pw_p_type), DIMENSION(3), &
      INTENT(OUT)                            :: df
    TYPE(realspace_grid_type), POINTER       :: rs_grid

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derive_fdm_cd3', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k
    INTEGER, DIMENSION(3)                    :: lb, ub
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: h
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: drdx, drdy, drdz, r
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    TYPE(realspace_grid_p_type), &
      DIMENSION(3)                           :: drs_grid

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPASSERT(ASSOCIATED(f))

    ! Setup
    rs_desc => rs_grid%desc
    CALL rs_pw_transfer(rs_grid,f,pw2rs)
    DO i=1,3
       CALL rs_grid_create(drs_grid(i)%rs_grid,rs_desc)
       CALL rs_grid_zero(drs_grid(i)%rs_grid)
    END DO

    lb(1:3) = rs_grid%lb_real(1:3)
    ub(1:3) = rs_grid%ub_real(1:3)
    r => rs_grid%r
    drdx => drs_grid(1)%rs_grid%r
    drdy => drs_grid(2)%rs_grid%r
    drdz => drs_grid(3)%rs_grid%r

    ! 3-point stencil central differences
    h(1:3) = 2.0_dp*f%pw_grid%dr(1:3)
    !$omp parallel do default(none) &
    !$omp             private(i,j,k) &
    !$omp             shared(drdx,drdy,drdz,h,lb,r,ub)
    DO k=lb(3),ub(3)
      DO j=lb(2),ub(2)
        DO i=lb(1),ub(1)
          drdx(i,j,k) = (r(i+1,j,k) - r(i-1,j,k))/h(1)
          drdy(i,j,k) = (r(i,j+1,k) - r(i,j-1,k))/h(2)
          drdz(i,j,k) = (r(i,j,k+1) - r(i,j,k-1))/h(3)
        END DO
      END DO
    END DO
    !$omp end parallel do

    ! Cleanup
    DO i=1,3
       CALL rs_pw_transfer(drs_grid(i)%rs_grid, df(i)%pw,rs2pw)
       CALL rs_grid_release(drs_grid(i)%rs_grid)
    END DO

    CALL timestop(handle)

  END SUBROUTINE derive_fdm_cd3

! *****************************************************************************
!> \brief    4th order finite difference derivative of a function on realspace grid
!> \param f  input function
!> \param df derivative of f
!> \param rs_grid real-space grid
!> \par History:
!>      - Creation (15.11.2013,MK)
!>      - Refactored and moved here from qs_sccs.F (12.2014, Hossein Bani-Hashemian)
!> \author     Matthias Krack (MK)
!> \version    1.0
! *****************************************************************************
  SUBROUTINE derive_fdm_cd5(f,df,rs_grid)

    TYPE(pw_type), POINTER                   :: f
    TYPE(pw_p_type), DIMENSION(3), &
      INTENT(INOUT)                          :: df
    TYPE(realspace_grid_type), POINTER       :: rs_grid

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derive_fdm_cd5', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k
    INTEGER, DIMENSION(3)                    :: lb, ub
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: h
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: drdx, drdy, drdz, r
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    TYPE(realspace_grid_p_type), &
      DIMENSION(3)                           :: drs_grid

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPASSERT(ASSOCIATED(f))

    ! Setup
    rs_desc => rs_grid%desc
    CALL rs_pw_transfer(rs_grid,f,pw2rs)
    DO i=1,3
       CALL rs_grid_create(drs_grid(i)%rs_grid,rs_desc)
       CALL rs_grid_zero(drs_grid(i)%rs_grid)
    END DO

    lb(1:3) = rs_grid%lb_real(1:3)
    ub(1:3) = rs_grid%ub_real(1:3)
    r => rs_grid%r
    drdx => drs_grid(1)%rs_grid%r
    drdy => drs_grid(2)%rs_grid%r
    drdz => drs_grid(3)%rs_grid%r

    ! 5-point stencil central differences
    h(1:3) = 12.0_dp*f%pw_grid%dr(1:3)
    !$omp parallel do default(none) &
    !$omp             private(i,j,k) &
    !$omp             shared(drdx,drdy,drdz,h,lb,r,ub)
    DO k=lb(3),ub(3)
      DO j=lb(2),ub(2)
        DO i=lb(1),ub(1)
          drdx(i,j,k) = (r(i-2,j,k) - r(i+2,j,k) + 8.0_dp*(r(i+1,j,k) - r(i-1,j,k)))/h(1)
          drdy(i,j,k) = (r(i,j-2,k) - r(i,j+2,k) + 8.0_dp*(r(i,j+1,k) - r(i,j-1,k)))/h(2)
          drdz(i,j,k) = (r(i,j,k-2) - r(i,j,k+2) + 8.0_dp*(r(i,j,k+1) - r(i,j,k-1)))/h(3)
        END DO
      END DO
    END DO
    !$omp end parallel do

    ! Cleanup
    DO i=1,3
       CALL rs_pw_transfer(drs_grid(i)%rs_grid, df(i)%pw,rs2pw)
       CALL rs_grid_release(drs_grid(i)%rs_grid)
    END DO

    CALL timestop(handle)

  END SUBROUTINE derive_fdm_cd5

! *****************************************************************************
!> \brief    6th order finite difference derivative of a function on realspace grid
!> \param f  input function
!> \param df derivative of f
!> \param rs_grid real-space grid
!> \par History:
!>      - Creation (15.11.2013,MK)
!>      - Refactored and moved here from qs_sccs.F (12.2014, Hossein Bani-Hashemian)
!> \author     Matthias Krack (MK)
!> \version    1.0
! *****************************************************************************
  SUBROUTINE derive_fdm_cd7(f,df,rs_grid)

    TYPE(pw_type), POINTER                   :: f
    TYPE(pw_p_type), DIMENSION(3), &
      INTENT(OUT)                            :: df
    TYPE(realspace_grid_type), POINTER       :: rs_grid

    CHARACTER(LEN=*), PARAMETER :: routineN = 'derive_fdm_cd7', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k
    INTEGER, DIMENSION(3)                    :: lb, ub
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(3)              :: h
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: drdx, drdy, drdz, r
    TYPE(realspace_grid_desc_type), POINTER  :: rs_desc
    TYPE(realspace_grid_p_type), &
      DIMENSION(3)                           :: drs_grid

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPASSERT(ASSOCIATED(f))

    ! Setup
    rs_desc => rs_grid%desc
    CALL rs_pw_transfer(rs_grid,f,pw2rs)
    DO i=1,3
       CALL rs_grid_create(drs_grid(i)%rs_grid,rs_desc)
       CALL rs_grid_zero(drs_grid(i)%rs_grid)
    END DO

    lb(1:3) = rs_grid%lb_real(1:3)
    ub(1:3) = rs_grid%ub_real(1:3)
    r => rs_grid%r
    drdx => drs_grid(1)%rs_grid%r
    drdy => drs_grid(2)%rs_grid%r
    drdz => drs_grid(3)%rs_grid%r

    ! 7-point stencil central differences
    h(1:3) = 60.0_dp*f%pw_grid%dr(1:3)
    !$omp parallel do default(none) &
    !$omp             private(i,j,k) &
    !$omp             shared(drdx,drdy,drdz,h,lb,r,ub)
    DO k=lb(3),ub(3)
      DO j=lb(2),ub(2)
        DO i=lb(1),ub(1)
          drdx(i,j,k) = (r(i+3,j,k) - r(i-3,j,k) + 9.0_dp*(r(i-2,j,k) - r(i+2,j,k)) +&
                         45.0_dp*(r(i+1,j,k) - r(i-1,j,k)))/h(1)
          drdy(i,j,k) = (r(i,j+3,k) - r(i,j-3,k) + 9.0_dp*(r(i,j-2,k) - r(i,j+2,k)) +&
                         45.0_dp*(r(i,j+1,k) - r(i,j-1,k)))/h(2)
          drdz(i,j,k) = (r(i,j,k+3) - r(i,j,k-3) + 9.0_dp*(r(i,j,k-2) - r(i,j,k+2)) +&
                         45.0_dp*(r(i,j,k+1) - r(i,j,k-1)))/h(3)
        END DO
      END DO
    END DO
    !$omp end parallel do

    ! Cleanup
    DO i=1,3
       CALL rs_pw_transfer(drs_grid(i)%rs_grid, df(i)%pw,rs2pw)
       CALL rs_grid_release(drs_grid(i)%rs_grid)
    END DO

    CALL timestop(handle)

  END SUBROUTINE derive_fdm_cd7

! *****************************************************************************
!> \brief returns the global axes and the portion of the axes that are local to
!>        the current mpi rank
!> \param pw_grid plane wave grid
!> \param x_glbl x grid vetor of the simulation box
!> \param y_glbl y grid vetor of the simulation box
!> \param z_glbl z grid vetor of the simulation box
!> \param x_locl x grid vetor of the simulation box local to this process
!> \param y_locl y grid vetor of the simulation box local to this process
!> \param z_locl z grid vetor of the simulation box local to this process
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!>       07.2015 moved here from dirichlet_bc_utils.F [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE setup_grid_axes(pw_grid, x_glbl, y_glbl, z_glbl, x_locl, y_locl, z_locl)

    TYPE(pw_grid_type), INTENT(IN), POINTER  :: pw_grid
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(OUT)                            :: x_glbl, y_glbl, z_glbl, &
                                                x_locl, y_locl, z_locl

    CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_grid_axes', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: glb1, glb2, glb3, gub1, gub2, &
                                                gub3, handle, i, lb1, lb2, &
                                                lb3, ub1, ub2, ub3
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: gindx, gindy, gindz, lindx, &
                                                lindy, lindz
    INTEGER, DIMENSION(2, 3)                 :: bounds, bounds_local
    INTEGER, DIMENSION(3)                    :: npts, npts_local
    REAL(dp), DIMENSION(3)                   :: dr

    CALL timeset(routineN,handle)

    dr = pw_grid%dr
    bounds = pw_grid%bounds
    bounds_local = pw_grid%bounds_local
    npts = pw_grid%npts
    npts_local = pw_grid%npts_local

! local and global lower and upper bounds
    glb1 = bounds(1,1);          gub1 = bounds(2,1)
    glb2 = bounds(1,2);          gub2 = bounds(2,2)
    glb3 = bounds(1,3);          gub3 = bounds(2,3)
    lb1  = bounds_local(1,1);    ub1  = bounds_local(2,1)
    lb2  = bounds_local(1,2);    ub2  = bounds_local(2,2)
    lb3  = bounds_local(1,3);    ub3  = bounds_local(2,3)

    ALLOCATE(lindx( lb1 : ub1 ), lindy( lb2 : ub2 ), lindz( lb3 : ub3 ))
    ALLOCATE(gindx(glb1 : gub1), gindy(glb2 : gub2), gindz(glb3 : gub3))
    ALLOCATE(x_locl( lb1 : ub1 ), y_locl( lb2 : ub2 ), z_locl( lb3 : ub3 ))
    ALLOCATE(x_glbl(glb1 : gub1), y_glbl(glb2 : gub2), z_glbl(glb3 : gub3))

    gindx(:) = (/ (i, i=0, npts(1) - 1) /)
    gindy(:) = (/ (i, i=0, npts(2) - 1) /)
    gindz(:) = (/ (i, i=0, npts(3) - 1) /)
    lindx(:) = (/ (i, i=0, npts_local(1) - 1) /)
    lindy(:) = (/ (i, i=0, npts_local(2) - 1) /)
    lindz(:) = (/ (i, i=0, npts_local(3) - 1) /)

    x_glbl(:) = gindx * dr(1)
    y_glbl(:) = gindy * dr(2)
    z_glbl(:) = gindz * dr(3)

! map [0 .. (npts_local-1)] --> [lb .. ub]
    IF (lb1 .EQ. ub1) THEN
       lindx(:) = lb1
    ELSE
       lindx(:) = lindx(:) * ((ub1 - lb1) / (npts_local(1) - 1)) + lb1
    END IF
    IF (lb2 .EQ. ub2) THEN
       lindy(:) = lb2
    ELSE
       lindy(:) = lindy(:) * ((ub2 - lb2) / (npts_local(2) - 1)) + lb2
    END IF
    IF (lb3 .EQ. ub3) THEN
       lindz(:) = lb3
    ELSE
       lindz(:) = lindz(:) * ((ub3 - lb3) / (npts_local(3) - 1)) + lb3
    END IF

    x_locl(:) = x_glbl(lindx)
    y_locl(:) = y_glbl(lindy)
    z_locl(:) = z_glbl(lindz)

    CALL timestop(handle)

  END SUBROUTINE setup_grid_axes

! *****************************************************************************
!> \brief convolutes a function with a smoothing kernel K_\zeta
!>                         v * K_\zeta
!> K_\zeta is the standard mollifier defined as:
!>        K_\zeta(x) = \frac{1}{\zeta^3} K(\frac{x}{\zeta})
!> where
!>        K(x) = \kappa \exp (\frac{1}{|x|^2 - 1}),  if |x| <= 1
!>             = 0,                                  otherwise
!> \param pw_pool pool of pw grid
!> \param zeta paramter \zeta defining the width of the mollifier
!> \param x_glbl x grid vetor of the simulation box
!> \param y_glbl y grid vetor of the simulation box
!> \param z_glbl z grid vetor of the simulation box
!> \param pw_in the input function
!> \param pw_out the convoluted function
!> \par History
!>       10.2014 created [Hossein Bani-Hashemian]
!>       07.2015 moved here from ps_implicit_methods.F [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE pw_mollifier(pw_pool, zeta, x_glbl, y_glbl, z_glbl, pw_in, pw_out)

    TYPE(pw_pool_type), INTENT(IN), POINTER  :: pw_pool
    REAL(dp), INTENT(IN)                     :: zeta
    REAL(dp), ALLOCATABLE, DIMENSION(:), &
      INTENT(IN)                             :: x_glbl, y_glbl, z_glbl
    TYPE(pw_type), INTENT(IN), POINTER       :: pw_in
    TYPE(pw_type), INTENT(INOUT), POINTER    :: pw_out

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pw_mollifier', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, k, lb1, lb2, &
                                                lb3, ub1, ub2, ub3
    INTEGER, DIMENSION(2, 3)                 :: bounds, bounds_local
    REAL(dp)                                 :: normfact, xi, xmax, xmin, yj, &
                                                ymax, ymin, zk, zmax, zmin
    REAL(dp), DIMENSION(3, 3)                :: dh
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_type), POINTER                   :: G, G_gs, pw_in_gs, pw_out_gs

    CALL timeset(routineN,handle)

    pw_grid => pw_in%pw_grid
    dh = pw_grid%dh
    bounds_local = pw_grid%bounds_local
    bounds = pw_grid%bounds

    lb1 = bounds_local(1,1);     ub1 = bounds_local(2,1)
    lb2 = bounds_local(1,2);     ub2 = bounds_local(2,2)
    lb3 = bounds_local(1,3);     ub3 = bounds_local(2,3)

    CALL pw_pool_create_pw(pw_pool, G, use_data=REALDATA3D, in_space=REALSPACE)
    CALL pw_pool_create_pw(pw_pool, G_gs, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
    CALL pw_pool_create_pw(pw_pool, pw_in_gs, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
    CALL pw_pool_create_pw(pw_pool, pw_out_gs, use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

    CALL pw_zero(G)
    xmin = x_glbl(bounds(1,1));     xmax = x_glbl(bounds(2,1))
    ymin = y_glbl(bounds(1,2));     ymax = y_glbl(bounds(2,2))
    zmin = z_glbl(bounds(1,3));     zmax = z_glbl(bounds(2,3))

    DO k = lb3, ub3
      DO j = lb2, ub2
        DO i = lb1, ub1
          xi = x_glbl(i);     yj = y_glbl(j);     zk = z_glbl(k)
          IF      (vec_norm2((/(xi - xmin), (yj - ymin), (zk - zmin)/)) .LT. zeta - small_value) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmin), (yj - ymin), (zk - zmin)/)/zeta)**2 - 1))
          ELSE IF (vec_norm2((/(xi - xmax), (yj - ymax), (zk - zmax)/)) .LT. zeta - small_value) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmax), (yj - ymax), (zk - zmax)/)/zeta)**2 - 1))
          ELSE IF (vec_norm2((/(xi - xmin), (yj - ymax), (zk - zmax)/)) .LT. zeta - small_value) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmin), (yj - ymax), (zk - zmax)/)/zeta)**2 - 1))
          ELSE IF (vec_norm2((/(xi - xmax), (yj - ymin), (zk - zmax)/)) .LT. zeta - small_value) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmax), (yj - ymin), (zk - zmax)/)/zeta)**2 - 1))
          ELSE IF (vec_norm2((/(xi - xmax), (yj - ymax), (zk - zmin)/)) .LT. zeta - small_value) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmax), (yj - ymax), (zk - zmin)/)/zeta)**2 - 1))
          ELSE IF (vec_norm2((/(xi - xmin), (yj - ymin), (zk - zmax)/)) .LT. zeta - small_value) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmin), (yj - ymin), (zk - zmax)/)/zeta)**2 - 1))
          ELSE IF (vec_norm2((/(xi - xmin), (yj - ymax), (zk - zmin)/)) .LT. zeta - small_value) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmin), (yj - ymax), (zk - zmin)/)/zeta)**2 - 1))
          ELSE IF (vec_norm2((/(xi - xmax), (yj - ymin), (zk - zmin)/)) .LT. zeta - small_value) THEN
             G%cr3d(i,j,k) = EXP(1.0_dp / (vec_norm2((/(xi - xmax), (yj - ymin), (zk - zmin)/)/zeta)**2 - 1))
          END IF
        END DO
      END DO
    END DO
    CALL pw_scale(G, (1.0_dp/zeta)**3)
    normfact = pw_integrate_function(G)
    CALL pw_scale(G, 1.0_dp/normfact)

    CALL pw_transfer(G, G_gs)
    CALL pw_transfer(pw_in, pw_in_gs)
    pw_out_gs%cc = G_gs%cc * pw_in_gs%cc
    CALL pw_transfer(pw_out_gs, pw_out)

    ! multiply by the reciprocal of the forward Fourier transform normalization prefactor (here 1/N, by convention)
    CALL pw_scale(pw_out, REAL(pw_grid%ngpts,KIND=dp))
    ! from discrete convolution to continuous convolution
    CALL pw_scale(pw_out, pw_grid%dvol)

    DO k = lb3, ub3
      DO j = lb2, ub2
        DO i = lb1, ub1
          IF (ABS(pw_out%cr3d(i,j,k)) .LE. 1.0E-10_dp) pw_out%cr3d(i,j,k) = 0.0_dp
        END DO
      END DO
    END DO

    CALL pw_pool_give_back_pw(pw_pool, G)
    CALL pw_pool_give_back_pw(pw_pool, G_gs)
    CALL pw_pool_give_back_pw(pw_pool, pw_in_gs)
    CALL pw_pool_give_back_pw(pw_pool, pw_out_gs)
    CALL timestop(handle)

    CONTAINS
! *****************************************************************************
!> \brief computes the norm 2 of a vector in the three-dimensional real space
!> \param vec input vector
!> \retval vec_norm2 ...
! *****************************************************************************
       PURE FUNCTION vec_norm2(vec)
    REAL(dp), DIMENSION(3), INTENT(IN)       :: vec
    REAL(dp)                                 :: vec_norm2

          vec_norm2 = SQRT(SUM(vec**2))
       END FUNCTION vec_norm2
  END SUBROUTINE pw_mollifier

END MODULE rs_methods
