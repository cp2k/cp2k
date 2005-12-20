!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2003 CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/real_space_test [1.0] *
!!
!!   NAME
!!     real_space_test
!!
!!   FUNCTION
!!     Test of real space grids
!!
!!   AUTHOR
!!     JGH (06.06.2003)
!!
!!   MODIFICATION HISTORY
!!     JGH (19.06.2003) : moved get_my_tasks to realspace_grid_types
!!                        check on center of grid instead of corner
!!                        generalized first dimension of tasks array
!!
!!*****
!******************************************************************************

MODULE real_space_test

! *****************************************************************************

  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_cputime
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_environ,&
                                             mp_sum
  USE pw_grid_types,                   ONLY: FULLSPACE,&
                                             PW_MODE_LOCAL, &
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_create,&
                                             pw_grid_setup
  USE pw_types,                        ONLY: REALDATA3D,&
                                             pw_create,&
                                             pw_release,&
                                             pw_type
  USE realspace_grid_types,            ONLY: realspace_grid_type,&
                                             rs_get_my_tasks,&
                                             realspace_grid_p_type,&
                                             rs_grid_release,&
                                             rs_grid_create,&
                                             rs_grid_zero,&
                                             rs_pw_transfer,&
                                             rs2pw,&
                                             pw2rs
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_hinv,&
                                             cell_create,&
                                             cell_release
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE util,                            ONLY: get_limit

#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: rs_test
  CHARACTER(len=*), PARAMETER :: moduleN = 'real_space_test'

CONTAINS

!******************************************************************************
SUBROUTINE rs_test ( globenv, error )

    TYPE(global_environment_type), &
      pointer                             :: globenv
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_test', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: nparticle = 5000

    INTEGER                                  :: dir, i, ierr, iw, j, me, &
                                                memory1, memory2, memory3, n, &
                                                npme
    INTEGER, DIMENSION(2)                    :: bo
    INTEGER, DIMENSION(3)                    :: no, np, pt
    INTEGER, DIMENSION(3, nparticle)         :: refco, refpt
    INTEGER, DIMENSION(:, :), POINTER        :: tasks
    INTEGER, DIMENSION(nparticle)            :: bsize
    REAL(KIND=dp)                                :: energy1, energy2, energy3, &
                                                integral, t1, t2, t3, tend, &
                                                tstart
    REAL(KIND=dp), DIMENSION(3)                  :: s
    REAL(KIND=dp), DIMENSION(3, nparticle)       :: coord
    REAL(KIND=dp), DIMENSION(nparticle)          :: bval
    TYPE(cell_type), POINTER                 :: box
    TYPE(pw_grid_type),pointer               :: grid
    TYPE(pw_type), POINTER                   :: pw_a, pw_b, pw_c
    TYPE(realspace_grid_type), POINTER       :: rden

  iw =globenv % scr
  NULLIFY(pw_a,pw_b,pw_c)
  NULLIFY(box)
  CALL cell_create(box)
  IF ( globenv % ionode ) &
       WRITE(iw,'(/,A,/)') " Test of Real Space Grid Routines"

!..the unit cell
! box % hmat = RESHAPE ( (/10.0_dp,0.0_dp,0.0_dp,0.0_dp,10.0_dp,0.0_dp,&
!                          0.0_dp,0.0_dp,10.0_dp/), (/3,3/) )
  box % hmat(1,1) = 10.0_dp
  box % hmat(1,2) =  0.0_dp
  box % hmat(1,3) =  0.0_dp
  box % hmat(2,1) =  0.0_dp
  box % hmat(2,2) = 10.0_dp
  box % hmat(2,3) =  0.0_dp
  box % hmat(3,1) =  0.0_dp
  box % hmat(3,2) =  0.0_dp
  box % hmat(3,3) = 10.0_dp

  CALL get_hinv ( box )

  np ( : ) =  128

  IF ( globenv % ionode ) THEN
    CALL random_number ( coord )
    coord = matmul ( box % hmat, coord )
    CALL random_number ( bval )
    bsize = MAX ( NINT ( 32.0_dp*bval ), 1 )
    CALL random_number ( bval )
    bval = bval * 0.001_dp
  END IF
  call mp_bcast ( coord, globenv%source, globenv%group )
  call mp_bcast ( bsize, globenv%source, globenv%group )
  call mp_bcast ( bval , globenv%source, globenv%group )

  DO i = 1, nparticle
    s = MATMUL ( box % h_inv, coord ( :, i ) )
    s = s - NINT ( s )
    refpt ( :, i ) = NINT ( REAL ( np ( : ),KIND=dp) * s ( : ) )
    IF ( refpt ( 1, i ) <= 0 ) refpt ( 1, i ) = refpt ( 1, i ) + np ( 1 )
    IF ( refpt ( 2, i ) <= 0 ) refpt ( 2, i ) = refpt ( 2, i ) + np ( 2 )
    IF ( refpt ( 3, i ) <= 0 ) refpt ( 3, i ) = refpt ( 3, i ) + np ( 3 )
    IF ( refpt ( 1, i ) > np(1) ) refpt ( 1, i ) = refpt ( 1, i ) - np ( 1 )
    IF ( refpt ( 2, i ) > np(2) ) refpt ( 2, i ) = refpt ( 2, i ) - np ( 2 )
    IF ( refpt ( 3, i ) > np(3) ) refpt ( 3, i ) = refpt ( 3, i ) - np ( 3 )
  END DO

  CALL pw_grid_create ( grid, globenv%group )
  grid % grid_span = FULLSPACE
  grid % bounds ( 1, : ) = -np / 2
  grid % bounds ( 2, : ) = ( np - 1 ) / 2

  CALL pw_grid_setup ( box, grid, &
                       info = globenv % scr, rs_dims=(/globenv%num_pe,1/) )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TEST 1, replicated calculation

  IF ( globenv % ionode ) &
    WRITE(iw,'(///,A,/)') " Start Test 1: Replicated grid, replicated calculation"
  tstart = m_cputime ( )
  CALL rs_grid_create ( rden, grid, -1)
  CALL rs_grid_zero ( rden )
  DO i = 1, nparticle
     refco(:,i) = refpt(:,i) + rden%lb(:) - 1
  END DO
  DO i = 1, nparticle
    no ( : ) = bsize ( i )
    pt = refco(:,i) - rden%lb(:) + 1
    CALL put_on_grid ( rden%r, pt, no, bval(i), 0 )
  END DO
  CALL pw_create ( pw_a, grid, REALDATA3D )
  CALL rs_pw_get ( rden, pw_a )
  pw_a%cr3d = pw_a%cr3d**(1.0_dp/3.0_dp)
  CALL rs_pw_transfer ( rden, pw_a, pw2rs )
  energy1 = 0.0_dp
  DO i = 1, nparticle
    pt = refco(:,i) - rden%lb(:) + 1
    no ( : ) = bsize ( i )
    CALL int_on_grid ( rden%r, pt, no, bval(i), 0, integral )
    energy1 = energy1 + integral
  END DO
  memory1 = SIZE ( rden%r ) + SIZE ( pw_a%cr3d )
  CALL rs_grid_release ( rden )
  tend  = m_cputime ( )
  t1 = tend - tstart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TEST 2, replicated grid, parallel calculation

  IF ( globenv % ionode ) &
    WRITE(iw,'(/,A,/)') " Start Test 2: Replicated grid, parallel calculation"
  tstart = m_cputime ( )
  CALL rs_grid_create ( rden, grid, -1, error )
  CALL rs_grid_zero ( rden )
  DO i = 1, nparticle
     refco(:,i) = refpt(:,i) + rden%lb(:) - 1
  END DO
  IF ( rden%parallel ) THEN
    CALL mp_environ ( n, me, rden%group )
  ELSE
    n  = 1
    me = 0
  END IF
  bo = get_limit ( nparticle, n, me )
  DO i = bo(1), bo(2)
    pt = refco(:,i) - rden%lb(:) + 1
    no ( : ) = bsize ( i )
    CALL put_on_grid ( rden%r, pt, no, bval(i), 0 )
  END DO
  CALL pw_create ( pw_b, grid, REALDATA3D )
  CALL rs_pw_transfer ( rden, pw_b, rs2pw )
  pw_b%cr3d = pw_b%cr3d**(1.0_dp/3.0_dp)
  CALL rs_grid_zero ( rden )
  CALL rs_pw_transfer ( rden, pw_b, pw2rs )
  energy2 = 0.0_dp
  DO i = bo(1), bo(2)
    pt = refco(:,i) - rden%lb(:) + 1
    no ( : ) = bsize ( i )
    CALL int_on_grid ( rden%r, pt, no, bval(i), 0, integral )
    energy2 = energy2 + integral
  END DO
  IF ( rden%parallel ) CALL mp_sum ( energy2, rden%group )
  memory2 = SIZE ( rden%r ) + SIZE ( pw_b%cr3d )
  CALL rs_grid_release ( rden )
  tend  = m_cputime ( )
  t2 = tend - tstart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TEST 3, distributed grid, parallel calculation

  IF ( globenv % ionode ) &
    WRITE(iw,'(/,A,/)') " Start Test 3: Distributed grid, parallel calculation"
  tstart = m_cputime ( )
  CALL rs_grid_create ( rden, grid, MAXVAL ( bsize(1:nparticle) ), error )
  CALL rs_grid_zero ( rden )
  dir = rden%direction
  DO i = 1, nparticle
     refco(:,i) = refpt(:,i) + rden%lb(:) - 1
  END DO
  IF ( rden%parallel ) THEN
    CALL mp_environ ( n, me, rden%group )
  ELSE
    n  = 1
    me = 0
  END IF
  bo = get_limit ( nparticle, n, me )
  npme = bo(2) - bo(1) + 1
  ALLOCATE ( tasks ( 2, npme ), STAT=ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "rs_test", "tasks", 2*npme )
  DO i = bo(1), bo(2)
    j = i - bo(1) + 1
    tasks ( 1, j ) = i
    tasks ( 2, j ) = refco(dir,i)
  END DO
  CALL rs_get_my_tasks ( rden, tasks, npme )
  DO i = 1, npme
    j = tasks ( 1, i )
    pt = refco(:,j) - rden%lb_local(:) + 1
    no = bsize ( j )
    CALL put_on_grid ( rden%r, pt, no, bval(j), dir )
  END DO
  DEALLOCATE ( tasks, STAT=ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "rs_test", "tasks" )
  CALL pw_create ( pw_c, grid, REALDATA3D )
  CALL rs_pw_transfer ( rden, pw_c, rs2pw )
  pw_c%cr3d = pw_c%cr3d**(1.0_dp/3.0_dp)
  CALL rs_grid_zero ( rden )
  CALL rs_pw_transfer ( rden, pw_c, pw2rs )
  IF ( rden%parallel ) THEN
    CALL mp_environ ( n, me, rden%group )
  ELSE
    n  = 1
    me = 0
  END IF
  bo = get_limit ( nparticle, n, me )
  npme = bo(2) - bo(1) + 1
  ALLOCATE ( tasks ( 2, npme ), STAT=ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "rs_test", "tasks", 3*npme )
  DO i = bo(1), bo(2)
    j = i - bo(1) + 1
    tasks ( 1, j ) = i
    tasks ( 2, j ) = refco(dir,i)
  END DO
  CALL rs_get_my_tasks ( rden, tasks, npme )
  energy3 = 0.0_dp
  DO i = 1, npme
    j = tasks ( 1, i )
    pt = refco(:,j) - rden%lb_local(:) + 1
    no = bsize ( j )
    CALL int_on_grid ( rden%r, pt, no, bval(j), dir, integral )
    energy3 = energy3 + integral
  END DO
  DEALLOCATE ( tasks, STAT=ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "rs_test", "tasks" )
  IF ( rden%parallel ) CALL mp_sum ( energy3, rden%group )
  memory3 = SIZE ( rden%r ) + SIZE ( pw_c%cr3d )
  CALL rs_grid_release ( rden )
  tend  = m_cputime ( )
  t3 = tend - tstart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF ( globenv % ionode ) THEN
    WRITE(iw,'(///,A,/)') " Results of Tests"
    WRITE(iw,'(A,T61,F20.4)') " Timing Test 1 [seconds] ",t1
    WRITE(iw,'(A,T61,F20.4)') " Timing Test 2 [seconds] ",t2
    WRITE(iw,'(A,T61,F20.4)') " Timing Test 3 [seconds] ",t3
    WRITE(iw,'(/,A,T61,I20)') " Memory Test 1 [words] ",memory1
    WRITE(iw,'(A,T61,I20)') " Memory Test 2 [words] ",memory2
    WRITE(iw,'(A,T61,I20)') " Memory Test 3 [words] ",memory3
    WRITE(iw,'(/,A,T61,F20.12)') " Integral Test 1 ",energy1*grid%dvol
    WRITE(iw,'(A,T61,F20.12)') " Integral Test 2 ",energy2*grid%dvol
    WRITE(iw,'(A,T61,F20.12)') " Integral Test 3 ",energy3*grid%dvol
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL pw_release ( pw_a )
  CALL pw_release ( pw_b )
  CALL pw_release ( pw_c )
  call cell_release(box)
END SUBROUTINE rs_test

! *****************************************************************************

SUBROUTINE put_on_grid ( r, nc, np, val, dir )

    REAL(KIND=dp), DIMENSION(:, :, :)            :: r
    INTEGER, DIMENSION(3)                    :: nc, np
    REAL(KIND=dp)                                :: val
    INTEGER                                  :: dir

    INTEGER                                  :: i, ii, j, jj, k, kk, ns, nx, &
                                                ny, nz
    INTEGER, DIMENSION(3)                    :: nl, nu

  nx = SIZE ( r, 1 )
  ny = SIZE ( r, 2 )
  nz = SIZE ( r, 3 )
  nl = -np/2
  nu = np + nl - 1
  IF ( dir > 0 ) THEN
    ns = SIZE ( r, dir )
    IF ( nc(dir) + nl(dir) <= 0 .OR. nc(dir) + nu(dir) > ns ) stop "out of bounds"
  END IF

  DO i = nl(3), nu(3)
    ii = nc(3) + i
    IF (dir/=3) ii = MODULO ( nz + nc(3) + i - 1, nz ) + 1
    DO j = nu(2), nu(2)
      jj = nc(2) + j 
      IF (dir/=2) jj = MODULO ( ny + nc(2) + j - 1, ny ) + 1
      DO k = nl(1), nu(1)
        kk = nc(1) + k
        IF (dir/=1) kk = MODULO ( nx + nc(1) + k - 1, nx ) + 1
        r ( kk, jj, ii ) = r ( kk, jj, ii ) + val
      END DO
    END DO
  END DO

END SUBROUTINE put_on_grid

! *****************************************************************************

SUBROUTINE int_on_grid ( r, nc, np, val, dir, integral )

    REAL(KIND=dp), DIMENSION(:, :, :)            :: r
    INTEGER, DIMENSION(3)                    :: nc, np
    REAL(KIND=dp)                                :: val
    INTEGER                                  :: dir
    REAL(KIND=dp)                                :: integral

    INTEGER                                  :: i, ii, j, jj, k, kk, ns, nx, &
                                                ny, nz
    INTEGER, DIMENSION(3)                    :: nl, nu

  integral = 0.0_dp
  nx = SIZE ( r, 1 )
  ny = SIZE ( r, 2 )
  nz = SIZE ( r, 3 )
  nl = -np/2
  nu = np + nl - 1
  IF ( dir > 0 ) THEN
    ns = SIZE ( r, dir )
    IF ( nc(dir) + nl(dir) <= 0 .OR. nc(dir) + nu(dir) > ns ) stop "out of bounds"
  END IF

  DO i = nl(3), nu(3)
    ii = nc(3) + i
    IF (dir/=3) ii = MODULO ( nz + nc(3) + i - 1, nz ) + 1
    DO j = nl(2), nu(2)
      jj = nc(2) + j
      IF (dir/=2) jj = MODULO ( ny + nc(2) + j - 1, ny ) + 1
      DO k = nl(1), nu(1)
        kk = nc(1) + k
        IF (dir/=1) kk = MODULO ( nx + nc(1) + k - 1, nx ) + 1
        integral = integral + r ( kk, jj, ii ) * val
      END DO
    END DO
  END DO

END SUBROUTINE int_on_grid

! *****************************************************************************

SUBROUTINE rs_pw_get ( rs, pw )

    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(pw_type), POINTER                   :: pw

    INTEGER                                  :: group, ierr, ii, ip, ix, iy, &
                                                iz, mepos, nn, np
    INTEGER, DIMENSION(3)                    :: lb, ub
    INTEGER, DIMENSION(:, :), POINTER        :: pbo
    INTEGER, DIMENSION(:, :, :), POINTER     :: bo
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)     :: rlocal

!-----------------------------------------------------------------------------!

   IF ( rs % group_dim ( 1 ) /= 1 ) &
       CALL stop_program ( "rs_pw_get", "rs grid has to be replicated" )
   IF ( pw % in_use /= REALDATA3D ) &
       CALL stop_program ( "rs_pw_get", "rs grid has to be REALDATA3D" )

   IF ( pw % pw_grid % para % mode == PW_MODE_LOCAL ) THEN

     nn = SIZE ( rs%r )
     CALL dcopy ( nn, rs%r, 1, pw % cr3d, 1 )

   ELSE

     np = pw % pw_grid % para % group_size
     bo => pw % pw_grid % para % bo (1:2,1:3,0:np-1,1)
     pbo => pw % pw_grid % bounds
     group = pw % pw_grid % para % rs_group
     mepos = pw % pw_grid % para % rs_mpo
     ip = mepos + 1
     nn = PRODUCT ( bo(2,:,ip) - bo(1,:,ip) + 1 )
     ALLOCATE ( rlocal ( nn ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_get", "rlocal", nn )

     ii = 0
     lb = pbo(1,:)+bo(1,:,ip)-1
     ub = pbo(1,:)+bo(2,:,ip)-1
     DO iz = lb(3), ub(3)
       DO iy = lb(2), ub(2)
         DO ix = lb(1), ub(1)
            ii=ii+1
            rlocal(ii) = rs%r(ix,iy,iz)
         END DO
       END DO
     END DO
     CALL dcopy ( nn, rlocal, 1, pw % cr3d, 1 )

     DEALLOCATE ( rlocal, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer", "rlocal" )

   END IF

END SUBROUTINE rs_pw_get

! *****************************************************************************

END MODULE real_space_test

! *****************************************************************************
