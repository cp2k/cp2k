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
!!     none
!!
!!*****
!******************************************************************************

MODULE real_space_test

! *****************************************************************************

  USE kinds, ONLY: dbl
  USE simulation_cell, ONLY : cell_type, get_hinv, get_cell_param
  USE pw_grid_types, ONLY : pw_grid_type, FULLSPACE
  USE pw_grids, ONLY : pw_grid_setup, pw_grid_construct, pw_grid_destruct
  USE pw_types, ONLY : pw_type, pw_allocate, pw_deallocate, pw_zero, REALDATA3D
  USE global_types, ONLY : global_environment_type
  USE termination, ONLY : stop_memory, stop_program
  USE realspace_grid_types, ONLY : realspace_grid_type, rs_grid_allocate, &
       rs_grid_deallocate, rs_grid_setup, rs_pw_transfer, rs_grid_zero
  USE message_passing, ONLY : mp_sum, mp_max, mp_bcast, mp_environ, &
       mp_cart_sub, mp_comm_free, mp_shift
  USE util, ONLY : get_limit
  USE memory_utilities, ONLY: reallocate
  USE machine, ONLY : m_cputime

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: rs_test

! *****************************************************************************

CONTAINS

!******************************************************************************
SUBROUTINE rs_test ( globenv )

  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

  INTEGER :: iw, ierr, i, j, dir, sdir, pp, n, me, npme
  TYPE ( pw_grid_type ) :: grid
  TYPE ( cell_type ) :: box
  INTEGER, DIMENSION ( :, : ), POINTER :: tasks, tasks_local
  INTEGER, DIMENSION ( 3 ) :: np, no, pt
  INTEGER, DIMENSION ( 2 ) :: bo
  REAL ( dbl ), DIMENSION ( 3 ) :: s
  REAL ( dbl ) :: tstart, tend, t1, t2, t3
  REAL ( dbl ) :: energy1, energy2, energy3, integral
  INTEGER :: memory1, memory2, memory3

  TYPE ( pw_type ) :: pw_a, pw_b, pw_c
  TYPE ( realspace_grid_type ) :: rden

  INTEGER, PARAMETER :: nparticle = 25000
  INTEGER, DIMENSION ( nparticle ) :: bsize
  REAL ( dbl ), DIMENSION ( nparticle ) :: bval
  REAL ( dbl ), DIMENSION ( 3, nparticle ) :: coord
  INTEGER, DIMENSION ( 3, nparticle ) :: refpt, refco

  iw =globenv % scr

  IF ( globenv % ionode ) &
       WRITE(iw,'(/,A,/)') " Test of Real Space Grid Routines"

!..the unit cell
  box % hmat = RESHAPE ( (/10._dbl,0._dbl,0._dbl,0._dbl,10._dbl,0._dbl,&
                           0._dbl,0._dbl,10._dbl/), (/3,3/) )
  CALL get_hinv ( box )

  np ( : ) = 128

  IF ( globenv % ionode ) THEN
    CALL random_number ( coord )
    coord = matmul ( box % hmat, coord )
    CALL random_number ( bval )
    bsize = MAX ( NINT ( 32._dbl*bval ), 1 )
    CALL random_number ( bval )
    bval = bval * 0.001_dbl
  END IF
  call mp_bcast ( coord, globenv%source, globenv%group )
  call mp_bcast ( bsize, globenv%source, globenv%group )
  call mp_bcast ( bval , globenv%source, globenv%group )

  DO i = 1, nparticle
    s = MATMUL ( box % h_inv, coord ( :, i ) )
    s = s - NINT ( s )
    refpt ( :, i ) = NINT ( REAL ( np ( : ), dbl ) * s ( : ) ) - bsize(i)/2
    IF ( refpt ( 1, i ) <= 0 ) refpt ( 1, i ) = refpt ( 1, i ) + np ( 1 )
    IF ( refpt ( 2, i ) <= 0 ) refpt ( 2, i ) = refpt ( 2, i ) + np ( 2 )
    IF ( refpt ( 3, i ) <= 0 ) refpt ( 3, i ) = refpt ( 3, i ) + np ( 3 )
    IF ( refpt ( 1, i ) > np(1) ) refpt ( 1, i ) = refpt ( 1, i ) - np ( 1 )
    IF ( refpt ( 2, i ) > np(2) ) refpt ( 2, i ) = refpt ( 2, i ) - np ( 2 )
    IF ( refpt ( 3, i ) > np(3) ) refpt ( 3, i ) = refpt ( 3, i ) - np ( 3 )
  END DO

  call pw_grid_construct ( grid )
  grid % grid_span = FULLSPACE
  grid % para % rs_dims ( 1 ) = globenv % num_pe
  grid % para % rs_dims ( 2 ) = 1
  grid % bounds ( 1, : ) = -np / 2
  grid % bounds ( 2, : ) = ( np - 1 ) / 2

  CALL pw_grid_setup ( box, grid, pe_group = globenv % group, &
                       info = globenv % scr )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TEST 1, replicated calculation

  IF ( globenv % ionode ) &
    WRITE(iw,'(///,A,/)') " Start Test 1: Replicated grid, replicated calculation"
  tstart = m_cputime ( )
  CALL rs_grid_setup ( rden, grid, -1, iw )
  CALL rs_grid_allocate ( rden )
  CALL rs_grid_zero ( rden )
  DO i = 1, nparticle
     refco(:,i) = refpt(:,i) + rden%lb(:) - 1
  END DO
  DO i = 1, nparticle
    no ( : ) = bsize ( i )
    pt = refco(:,i) - rden%lb(:) + 1
    CALL put_on_grid ( rden%r, pt, no, bval(i), 0 )
  END DO
  CALL pw_allocate ( pw_a, grid, REALDATA3D )
  CALL rs_pw_get ( rden, pw_a )
  pw_a%cr3d = pw_a%cr3d**(1._dbl/3._dbl)
  CALL rs_pw_transfer ( rden, pw_a, "BACKWARD" )
  energy1 = 0._dbl
  DO i = 1, nparticle
    pt = refco(:,i) - rden%lb(:) + 1
    no ( : ) = bsize ( i )
    CALL int_on_grid ( rden%r, pt, no, bval(i), 0, integral )
    energy1 = energy1 + integral
  END DO
  memory1 = SIZE ( rden%r ) + SIZE ( pw_a%cr3d )
  CALL rs_grid_deallocate ( rden )
  tend  = m_cputime ( )
  t1 = tend - tstart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TEST 2, replicated grid, parallel calculation

  IF ( globenv % ionode ) &
    WRITE(iw,'(/,A,/)') " Start Test 2: Replicated grid, parallel calculation"
  tstart = m_cputime ( )
  CALL rs_grid_setup ( rden, grid, -1, iw )
  CALL rs_grid_allocate ( rden )
  CALL rs_grid_zero ( rden )
  DO i = 1, nparticle
     refco(:,i) = refpt(:,i) + rden%lb(:) - 1
  END DO
  CALL mp_environ ( n, me, rden%group )
  bo = get_limit ( nparticle, n, me )
  DO i = bo(1), bo(2)
    pt = refco(:,i) - rden%lb(:) + 1
    no ( : ) = bsize ( i )
    CALL put_on_grid ( rden%r, pt, no, bval(i), 0 )
  END DO
  CALL pw_allocate ( pw_b, grid, REALDATA3D )
  CALL rs_pw_transfer ( rden, pw_b, "FORWARD" )
  pw_b%cr3d = pw_b%cr3d**(1._dbl/3._dbl)
  CALL rs_grid_zero ( rden )
  CALL rs_pw_transfer ( rden, pw_b, "BACKWARD" )
  energy2 = 0._dbl
  DO i = bo(1), bo(2)
    pt = refco(:,i) - rden%lb(:) + 1
    no ( : ) = bsize ( i )
    CALL int_on_grid ( rden%r, pt, no, bval(i), 0, integral )
    energy2 = energy2 + integral
  END DO
  CALL mp_sum ( energy2, rden%group )
  memory2 = SIZE ( rden%r ) + SIZE ( pw_b%cr3d )
  CALL rs_grid_deallocate ( rden )
  tend  = m_cputime ( )
  t2 = tend - tstart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TEST 3, distributed grid, parallel calculation

  IF ( globenv % ionode ) &
    WRITE(iw,'(/,A,/)') " Start Test 3: Distributed grid, parallel calculation"
  tstart = m_cputime ( )
  CALL rs_grid_setup ( rden, grid, MAXVAL ( bsize ), iw )
  CALL rs_grid_allocate ( rden )
  CALL rs_grid_zero ( rden )
  dir = rden%direction
  DO i = 1, nparticle
     refco(:,i) = refpt(:,i) + rden%lb(:) - 1
  END DO
  IF ( dir > 0 ) THEN
    DO i = 1, nparticle
       IF ( refco(dir,i) >= rden%ub(dir) - rden%border + 1 ) &
         refco(dir,i) = refco(dir,i) - rden%npts(dir)
    END DO
  END IF
  CALL mp_environ ( n, me, rden%group )
  bo = get_limit ( nparticle, n, me )
  npme = bo(2) - bo(1) + 1
  ALLOCATE ( tasks ( 3, npme ), STAT=ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "rs_test", "tasks", 3*npme )
  DO i = bo(1), bo(2)
    j = i - bo(1) + 1
    tasks ( 1, j ) = i
    tasks ( 2, j ) = refco(dir,i)
    tasks ( 3, j ) = bsize(i)
  END DO
  ALLOCATE ( tasks_local ( 3, npme ), STAT=ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "rs_test", "tasks_local", 1 )
  CALL get_my_tasks ( rden, tasks, tasks_local, npme )
  DO i = 1, npme
    j = tasks_local ( 1, i )
    pt = refco(:,j) - rden%lb_local(:) + 1
    no = bsize ( j )
    CALL put_on_grid ( rden%r, pt, no, bval(j), dir )
  END DO
  DEALLOCATE ( tasks, STAT=ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "rs_test", "tasks" )
  DEALLOCATE ( tasks_local, STAT=ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "rs_test", "tasks_local" )
  CALL pw_allocate ( pw_c, grid, REALDATA3D )
  CALL rs_pw_transfer ( rden, pw_c, "FORWARD" )
  pw_c%cr3d = pw_c%cr3d**(1._dbl/3._dbl)
  CALL rs_grid_zero ( rden )
  CALL rs_pw_transfer ( rden, pw_c, "BACKWARD" )
  CALL mp_environ ( n, me, rden%group )
  bo = get_limit ( nparticle, n, me )
  npme = bo(2) - bo(1) + 1
  ALLOCATE ( tasks ( 3, npme ), STAT=ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "rs_test", "tasks", 3*npme )
  DO i = bo(1), bo(2)
    j = i - bo(1) + 1
    tasks ( 1, j ) = i
    tasks ( 2, j ) = refco(dir,i)
    tasks ( 3, j ) = bsize(i)
  END DO
  ALLOCATE ( tasks_local ( 3, npme ), STAT=ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "rs_test", "tasks_local", 1 )
  CALL get_my_tasks ( rden, tasks, tasks_local, npme )
  energy3 = 0._dbl
  DO i = 1, npme
    j = tasks_local ( 1, i )
    pt = refco(:,j) - rden%lb_local(:) + 1
    no = bsize ( j )
    CALL int_on_grid ( rden%r, pt, no, bval(j), dir, integral )
    energy3 = energy3 + integral
  END DO
  DEALLOCATE ( tasks, STAT=ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "rs_test", "tasks" )
  DEALLOCATE ( tasks_local, STAT=ierr )
  IF ( ierr /= 0 ) CALL stop_memory ( "rs_test", "tasks_local" )
  CALL mp_sum ( energy3, rden%group )
  memory3 = SIZE ( rden%r ) + SIZE ( pw_c%cr3d )
  CALL rs_grid_deallocate ( rden )
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

  CALL pw_deallocate ( pw_a )
  CALL pw_deallocate ( pw_b )
  CALL pw_deallocate ( pw_c )

END SUBROUTINE rs_test

! *****************************************************************************

SUBROUTINE put_on_grid ( r, nc, np, val, dir )

  REAL ( dbl ), DIMENSION (:,:,:) :: r
  INTEGER, DIMENSION ( 3 ) :: nc, np
  INTEGER :: dir
  REAL ( dbl ) :: val
 
  INTEGER :: i, j, k, ii, jj, kk, nx, ny, nz, ns

  nx = SIZE ( r, 1 )
  ny = SIZE ( r, 2 )
  nz = SIZE ( r, 3 )
  IF ( dir > 0 ) THEN
    ns = SIZE ( r, dir )
    IF ( nc(dir) <= 0 .OR. nc(dir) + np(dir) - 1 > ns ) stop "out of bounds"
  END IF

  DO i = 1, np(3)
    ii = nc(3) + i -1
    IF (dir/=3) ii = MOD ( nc(3) + i - 2, nz ) + 1
    DO j = 1, np(2)
      jj = nc(2) + j -1
      IF (dir/=2) jj = MOD ( nc(2) + j - 2, ny ) + 1
      DO k = 1, np(1)
        kk = nc(1) + k -1
        IF (dir/=1) kk = MOD ( nc(1) + k - 2, nx ) + 1
        r ( kk, jj, ii ) = r ( kk, jj, ii ) + val
      END DO
    END DO
  END DO

END SUBROUTINE put_on_grid

! *****************************************************************************

SUBROUTINE int_on_grid ( r, nc, np, val, dir, integral )

  REAL ( dbl ), DIMENSION (:,:,:) :: r
  INTEGER, DIMENSION ( 3 ) :: nc, np
  INTEGER :: dir
  REAL ( dbl ) :: val, integral

  INTEGER :: i, j, k, ii, jj, kk, nx, ny, nz, ns

  integral = 0._dbl
  nx = SIZE ( r, 1 )
  ny = SIZE ( r, 2 )
  nz = SIZE ( r, 3 )
  IF ( dir > 0 ) THEN
    ns = SIZE ( r, dir )
    IF ( nc(dir) <= 0 .OR. nc(dir) + np(dir) - 1 > ns ) stop "out of bounds"
  END IF

  DO i = 1, np(3)
    ii = nc(3) + i - 1
    IF (dir/=3) ii = MOD ( nc(3) + i - 2, nz ) + 1
    DO j = 1, np(2)
      jj = nc(2) + j - 1
      IF (dir/=2) jj = MOD ( nc(2) + j - 2, ny ) + 1
      DO k = 1, np(1)
        kk = nc(1) + k - 1
        IF (dir/=1) kk = MOD ( nc(1) + k - 2, nx ) + 1
        integral = integral + r ( kk, jj, ii ) * val
      END DO
    END DO
  END DO

END SUBROUTINE int_on_grid

! *****************************************************************************

SUBROUTINE rs_pw_get ( rs, pw )

!  Arguments
   TYPE ( realspace_grid_type ), INTENT ( INOUT ) :: rs
   TYPE ( pw_type ), TARGET, INTENT ( INOUT ) :: pw

!  Local
   INTEGER, DIMENSION ( :, : ), POINTER :: pbo
   INTEGER, DIMENSION ( 3 ) :: lb, ub, lc, uc
   INTEGER :: subgroup, source, dest, ierr, idir, nn, handle
   LOGICAL :: subdim ( 2 )
   REAL ( dbl ), DIMENSION ( :, :, : ), ALLOCATABLE :: buffer

   INTEGER :: np, ip, ix, iy, iz, ii, nma, mepos, group
   REAL ( dbl ), DIMENSION ( : ), ALLOCATABLE :: rlocal
   INTEGER, DIMENSION ( : ), ALLOCATABLE :: rcount
   INTEGER, DIMENSION (:,:,:), POINTER :: bo

!-----------------------------------------------------------------------------!
   IF ( rs % group_dim ( 1 ) /= 1 ) &
       CALL stop_program ( "rs_pw_get", "rs grid has to be replicated" )
   IF ( pw % in_use /= REALDATA3D ) &
       CALL stop_program ( "rs_pw_get", "rs grid has to be REALDATA3D" )

   IF ( pw % pw_grid % para % mode == 0 ) THEN

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

SUBROUTINE get_my_tasks ( rs, tasks, tasks_local, npme )

  TYPE ( realspace_grid_type ) :: rs
  INTEGER, DIMENSION ( :,: ), POINTER :: tasks, tasks_local
  INTEGER :: npme
  INTEGER :: subgroup, subsize, subpos
  INTEGER :: isend, ileft, i, j, lb, ub
  LOGICAL, DIMENSION ( 2 ) :: subdim
 
  IF ( .NOT. ASSOCIATED ( tasks ) ) &
        CALL stop_program ( "get_my_tasks", "tasks not associated" )
  IF ( .NOT. ASSOCIATED ( tasks_local ) ) &
        CALL stop_program ( "get_my_tasks", "tasks_local not associated" )

  IF ( rs%parallel .AND. rs%direction > 0 ) THEN

    npme = 0
    ! bounds of local grid
    lb = rs%lb_local ( rs%direction )
    ub = rs%ub_local ( rs%direction )

    ! We need subgroups along cartesian coordinate 1.
    subdim ( 1 ) = .TRUE.
    subdim ( 2 ) = .FALSE.
    CALL mp_cart_sub ( rs % group, subdim, subgroup )
    CALL mp_environ ( subsize, subpos, subgroup )

    ! keep tasks to be processed locally
    ! send remaining tasks to the next processor
    DO isend = 0, subsize - 1
      IF ( isend == 0 ) THEN
        ileft = SIZE ( tasks, 2)
      ELSE
        ! count left over tasks
        j = 0
        DO i = 1, ileft
          IF ( tasks(1,i) /= 0 ) THEN
            j = j + 1
            tasks(:,j) = tasks(:,i)
          END IF
        END DO
        tasks(:,j+1:ileft) = 0
        ileft = j
        CALL mp_max ( ileft, subgroup )
        IF ( ileft == 0 ) EXIT
        ! check if local send/receive array is big enough
        IF ( ileft > SIZE ( tasks, 2 ) ) THEN
          CALL reallocate ( tasks, 1, 3, 1, ileft )
        END IF
        ! send/receive tasks
        CALL mp_shift ( tasks ( 1:3, 1:ileft ), subgroup )
      END IF
      ! look for tasks do be done on this processor
      DO i = 1, ileft
        IF ( tasks(1,i) > 0 ) THEN
          IF ( tasks(2,i) >= lb .AND. tasks(2,i) + tasks(3,i) - 1 <= ub ) THEN
            ! found new local task
            npme = npme + 1
            IF ( npme > SIZE ( tasks_local, 2 ) ) THEN
              CALL reallocate ( tasks_local, 1, 3, 1, 2*npme )
            END IF
            tasks_local ( :, npme ) = tasks ( :, i )
            tasks ( :, i ) = 0
          END IF
        END IF
      END DO
    END DO
    ! are all tasks distributed?
    j = 0
    DO i = 1, SIZE ( tasks, 2 )
       IF ( tasks(1,i) /= 0 ) j = j + 1
    END DO
    IF ( j /= 0 ) CALL stop_program ( "get_my_tasks", "left over tasks" )

    ! Release the communicator
    CALL mp_comm_free ( subgroup )

  ELSE

    ! fully replicated grids, each processor can process all its tasks
    npme = SIZE ( tasks, 2)
    IF ( SIZE ( tasks_local, 2 ) < npme ) THEN
      CALL reallocate ( tasks_local, 1, 3, 1, npme )
    END IF
    tasks_local = tasks

  END IF

END SUBROUTINE get_my_tasks

! *****************************************************************************

END MODULE real_space_test

! *****************************************************************************
