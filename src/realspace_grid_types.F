!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****s* cp2k/realspace_grid_types [1.0] *
!!
!!   NAME
!!     realspace_grid_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH (18-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (22-May-2002) : New routine rs_grid_zero
!!     JGH (12-Jun-2002) : Bug fix for mpi groups
!!     JGH (19-Jun-2003) : Added routine for task distribution
!!     JGH (23-Nov-2003) : Added routine for task loop separation
!!
!!   NOTES
!!     Basic type for real space grid methods
!!
!!   SOURCE
!******************************************************************************

MODULE realspace_grid_types
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dbl
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_cart_create,&
                                             mp_cart_shift,&
                                             mp_cart_sub,&
                                             mp_comm_dup,&
                                             mp_comm_free,&
                                             mp_environ,&
                                             mp_max,&
                                             mp_min,&
                                             mp_sendrecv,&
                                             mp_shift,&
                                             mp_sum
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA1D,&
                                             REALDATA3D,&
                                             pw_type
  USE sparse_matrix_types,             ONLY: add_block_node,&
                                             get_block_node,&
                                             real_matrix_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: realspace_grid_type, realspace_grid_p_type
  PUBLIC :: rs_grid_allocate, rs_grid_deallocate, &
       rs_grid_setup, rs_pw_transfer, rs_grid_zero, rs_grid_to_cube, &
       rs_pw_to_cube, rs_grid_write
  PUBLIC :: rs_get_my_tasks, rs_get_loop_vars
  PUBLIC :: rs_grid_create, rs_grid_retain, rs_grid_release

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="realspace_grid_types"
  INTEGER, SAVE, PRIVATE :: last_rs_id=0
  INTEGER, SAVE, PRIVATE :: allocated_rs_grid_count=0

  TYPE realspace_grid_type
     INTEGER :: identifier                               ! tag of the pw_grid
     INTEGER :: id_nr                                    ! unique identifier of rs
     INTEGER :: ref_count                                ! reference count
     INTEGER :: ngpts                                    ! # grid points
     INTEGER, DIMENSION (3) :: npts                      ! # grid points per dimension
     INTEGER, DIMENSION (3) :: lb                        ! lower bounds
     INTEGER, DIMENSION (3) :: ub                        ! upper bounds
     INTEGER, DIMENSION (:), POINTER :: px,py,pz         ! index translators
     REAL ( dbl ), DIMENSION ( :, :, : ),POINTER :: r    ! the grid
     INTEGER :: border                                   ! border points
     INTEGER :: ngpts_local                              ! local dimensions
     INTEGER, DIMENSION (3) :: npts_local
     INTEGER, DIMENSION (3) :: lb_local
     INTEGER, DIMENSION (3) :: ub_local
     REAL(dbl ), DIMENSION(3) :: dr                      ! grid spacing
     LOGICAL :: parallel
     INTEGER :: group
     LOGICAL :: group_head
     INTEGER, DIMENSION (2) :: group_dim
     INTEGER, DIMENSION (2) :: group_coor
     INTEGER, DIMENSION (3) :: perd                      ! periodicity enforced
     INTEGER :: direction                                ! non-periodic direction 
  END TYPE realspace_grid_type

  TYPE realspace_grid_p_type
     TYPE(realspace_grid_type), POINTER :: rs_grid
  END TYPE realspace_grid_p_type

  INTERFACE rs_grid_allocate
     MODULE PROCEDURE rs_grid_allocate_1, rs_grid_allocate_n
  END INTERFACE

  INTERFACE rs_grid_deallocate
     MODULE PROCEDURE rs_grid_deallocate_1, rs_grid_deallocate_n
  END INTERFACE

  INTERFACE rs_grid_setup
     MODULE PROCEDURE rs_grid_setup_1, rs_grid_setup_n
  END INTERFACE

  !*****
  !******************************************************************************

CONTAINS

  !******************************************************************************
!!****** realspace_grid_types/rs_grid_allocate [1.0] *
!!
!!   NAME
!!     rs_grid_allocate
!!
!!   FUNCTION
!!     Allocate real space grids
!!
!!   AUTHOR
!!     JGH (18-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
  !******************************************************************************

  SUBROUTINE rs_grid_allocate_1 ( rs )

    !  Arguments
    TYPE ( realspace_grid_type ), TARGET, INTENT (INOUT) :: rs

    !  Local
    INTEGER :: ierr,handle
    INTEGER, DIMENSION (:), POINTER :: lb, ub

    !-----------------------------------------------------------------------------!

    CALL timeset("rs_grid_allocate",'I',' ',handle)

    allocated_rs_grid_count = allocated_rs_grid_count + 1
    ! write(6,*) "allocated_rs_grid_count ",allocated_rs_grid_count

    lb => rs % lb_local
    ub => rs % ub_local

    ALLOCATE ( rs % r (lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "rs_grid_allocate","rs % r", &
         rs % ngpts_local )
    ALLOCATE ( rs % px ( rs % npts ( 1 ) ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "rs_grid_allocate","rs % px", &
         rs % npts ( 1 ) )
    ALLOCATE ( rs % py ( rs % npts ( 2 ) ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "rs_grid_allocate","rs % py", &
         rs % npts ( 2 ) )
    ALLOCATE ( rs % pz ( rs % npts ( 3 ) ), STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "rs_grid_allocate","rs % pz", &
         rs % npts ( 3 ) )

    CALL timestop(0.0_dbl,handle)
  END SUBROUTINE rs_grid_allocate_1

  !******************************************************************************

  SUBROUTINE rs_grid_allocate_n ( rs )

    !  Arguments
    TYPE ( realspace_grid_type ), DIMENSION ( : ), INTENT (INOUT) :: rs

    !  Local
    INTEGER :: n, i

    !-----------------------------------------------------------------------------!

    n = SIZE ( rs )

    DO i = 1, n

       CALL rs_grid_allocate_1 ( rs ( i ) )

    END DO

  END SUBROUTINE rs_grid_allocate_n

  !*****
  !******************************************************************************
!!****** realspace_grid_types/rs_grid_deallocate [1.0] *
!!
!!   NAME
!!     rs_grid_deallocate
!!
!!   FUNCTION
!!     Deallocate real space grids
!!
!!   AUTHOR
!!     JGH (18-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
  !******************************************************************************

  SUBROUTINE rs_grid_deallocate_1 ( rs )

    !  Arguments
    TYPE ( realspace_grid_type ), INTENT (INOUT) :: rs

    !  Local
    INTEGER :: ierr,handle

    !-----------------------------------------------------------------------------!
    CALL timeset("rs_grid_dealloc",'I',' ',handle)
    allocated_rs_grid_count=allocated_rs_grid_count-1
    DEALLOCATE ( rs % r, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "rs_grid_deallocate","rs % r" )
    DEALLOCATE ( rs % px , STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "rs_grid_deallocate","rs % px" )
    DEALLOCATE ( rs % py , STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "rs_grid_deallocate","rs % py" )
    DEALLOCATE ( rs % pz , STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "rs_grid_deallocate","rs % pz" )
    rs % identifier = 0
    IF ( rs % parallel ) THEN
       ! release the group communicator
       CALL mp_comm_free ( rs % group )
    END IF

    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE rs_grid_deallocate_1

  !******************************************************************************

  SUBROUTINE rs_grid_deallocate_n ( rs )

    !  Arguments
    TYPE ( realspace_grid_type ), DIMENSION ( : ), INTENT (INOUT) :: rs

    !  Local
    INTEGER :: n, i, ierr

    !-----------------------------------------------------------------------------!

    n = SIZE ( rs )

    DO i = 1, n
       CALL rs_grid_deallocate_1 ( rs ( i ) )
    END DO

  END SUBROUTINE rs_grid_deallocate_n

  !*****
  !******************************************************************************
!!****** realspace_grid_types/rs_grid_setup [1.0] *
!!
!!   NAME
!!     rs_grid_setup
!!
!!   FUNCTION
!!     Determine the setup of real space grids
!!
!!   AUTHOR
!!     JGH (18-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (08-Jun-2003) : nsmax <= 0 indicates fully replicated grid
!!
!!   SOURCE
  !******************************************************************************

  SUBROUTINE rs_grid_setup_1 ( rs, pw_grid, nsmax, iounit )

    !  Arguments
    TYPE ( realspace_grid_type ), INTENT (INOUT) :: rs
    TYPE ( pw_grid_type ), INTENT ( IN ) :: pw_grid
    INTEGER, INTENT ( IN ) :: nsmax 
    INTEGER, OPTIONAL, INTENT ( IN ) :: iounit

    !  Local
    INTEGER :: dir, nmin, ngmax, ip, n_slices, pos ( 2 ), lb ( 2 ), nn, i
    REAL ( dbl ) :: pp ( 3 )

    !-----------------------------------------------------------------------------!
    rs % dr = pw_grid%dr
    rs % identifier = pw_grid % identifier
    last_rs_id = last_rs_id+1
    rs % id_nr = last_rs_id
    rs % ref_count = 1
    IF ( pw_grid % para % mode == 0 ) THEN
       ! The corresponding group has dimension 1 
       ! All operations will be done localy
       rs % npts = pw_grid % npts
       rs % ngpts = PRODUCT ( rs % npts )
       rs % lb = pw_grid % bounds ( 1, : )
       rs % ub = pw_grid % bounds ( 2, : )
       rs % perd = 1
       rs % direction = 0
       rs % border = 0
       rs % npts_local = pw_grid % npts
       rs % ngpts_local = PRODUCT ( rs % npts )
       rs % lb_local = pw_grid % bounds ( 1, : )
       rs % ub_local = pw_grid % bounds ( 2, : )
       rs % parallel = .FALSE.
       rs % group = 0
       rs % group_head = .TRUE.
       rs % group_dim = 1
       rs % group_coor = 0
    ELSE
       IF ( .NOT. ALL ( pw_grid % bounds ( 1, 2:3 ) == &
            pw_grid % bounds_local ( 1, 2:3 ) .AND. &
            pw_grid % bounds ( 2, 2:3 ) == &
            pw_grid % bounds_local ( 2, 2:3 ) ) ) THEN
          CALL stop_program ( "rs_grid_setup_1", &
               "This pw type not supported" )
       END IF
       ! x is the main distribution direction of pw grids in real space
       ! we distribute real space grids in z direction
       dir = 3
       rs % direction = 3
       IF ( nsmax <= 0 ) THEN
          n_slices = 1
       ELSE
          ! the minimum slice thickness is half of the small box size
          ! this way there is only overlap with direct neighbors
          nmin = ( nsmax + 1 ) / 2 
          ! maximum number of slices
          ngmax = pw_grid % npts ( dir ) / nmin
          ! now we reduce this number until we find a divisor of the total number
          ! of processors
          n_slices = 1
          DO ip = ngmax, 1, -1
             IF ( MOD ( pw_grid % para % group_size, ip ) == 0 ) THEN
                n_slices = ip
                EXIT
             END IF
          END DO
       END IF
       ! global grid dimensions are still the same
       rs % npts = pw_grid % npts
       rs % ngpts = PRODUCT ( rs % npts )
       rs % lb = pw_grid % bounds ( 1, : )
       rs % ub = pw_grid % bounds ( 2, : )
       rs % group_dim ( 1 ) = n_slices
       rs % group_dim ( 2 ) = pw_grid % para % group_size / n_slices
       IF ( n_slices == 1 ) THEN
          ! CASE 1 : only one slice: we do not need overlapping regions and special
          !          recombination of the total density
          rs % direction = 0
          rs % perd = 1
          rs % border = 0
          rs % npts_local = pw_grid % npts
          rs % ngpts_local = PRODUCT ( rs % npts )
          rs % lb_local = pw_grid % bounds ( 1, : )
          rs % ub_local = pw_grid % bounds ( 2, : )
          rs % parallel = .TRUE.
          CALL mp_comm_dup( pw_grid % para % rs_group , rs % group )
          rs % group_head = pw_grid % para % group_head
          rs % group_coor ( 1 ) = 0
          rs % group_coor ( 2 ) = pw_grid % para % rs_mpo
       ELSE
          ! CASE 2 : general case
          ! periodicity is no longer enforced along direction dir
          rs % perd = 1
          rs % perd ( dir ) = 0
          ! we keep a border of nmin points
          rs % border = nmin
          ! we are going parallel on the real space grid
          rs % parallel = .TRUE.
          ! the new cartesian group
          CALL mp_cart_create ( pw_grid % para % group, 2, &
               rs % group_dim, pos, rs % group )
          rs % group_head = ALL ( pos == 0 )
          rs % group_coor = pos
          ! local dimensions of the grid
          rs % npts_local = pw_grid % npts
          lb = get_limit ( pw_grid % npts ( dir ), rs % group_dim ( 1 ), pos ( 1 ) )
          rs % npts_local ( dir ) = 2 * rs % border + ( lb ( 2 ) - lb ( 1 ) + 1 )
          rs % ngpts_local = PRODUCT ( rs % npts_local )
          rs % lb_local = pw_grid % bounds ( 1, : )
          rs % ub_local = pw_grid % bounds ( 2, : )
          rs % lb_local ( dir ) = lb ( 1 ) + pw_grid % bounds ( 1, dir ) - rs % border - 1
          rs % ub_local ( dir ) = lb ( 2 ) + pw_grid % bounds ( 1, dir ) + rs % border - 1
       END IF

    END IF

    IF ( PRESENT ( iounit ) ) THEN
       IF ( iounit >= 0 ) THEN
          IF ( rs % parallel ) THEN
             IF ( rs % group_head ) THEN
                WRITE ( iounit, '(/,A,T71,I10)' ) &
                     " RS_GRID: Information for grid number ", rs % identifier
                DO i = 1, 3
                   WRITE ( iounit, '(A,I3,T30,2I8,T62,A,T71,I10)' ) " RS_GRID:   Bounds ", &
                        i, rs % lb ( i ), rs % ub ( i ), "Points:", rs % npts ( i )
                END DO
                IF ( rs % group_dim ( 1 ) == 1 ) THEN
                   WRITE ( iounit, '(A)' ) " RS_GRID: Real space fully replicated"
                   WRITE ( iounit, '(A,T71,I10)' ) &
                        " RS_GRID: Group size ", rs % group_dim ( 2 )
                ELSE
                   WRITE ( iounit, '(A,T71,I3,A)' ) &
                        " RS_GRID: Real space distribution over ", rs % group_dim ( 1 ), " groups"
                   WRITE ( iounit, '(A,T71,I10)' ) &
                        " RS_GRID: Group size ", rs % group_dim ( 2 )
                   WRITE ( iounit, '(A,T71,I10)' ) &
                        " RS_GRID: Real space distribution along direction ", rs % direction
                   WRITE ( iounit, '(A,T71,I10)' ) &
                        " RS_GRID: Border size ", rs % border
                END IF
             END IF
             nn = rs % npts_local ( dir )
             CALL mp_sum ( nn, rs % group )
             pp ( 1 ) = REAL ( nn, dbl ) / REAL ( PRODUCT ( rs % group_dim ), dbl )
             nn = rs % npts_local ( dir )
             CALL mp_max ( nn, rs % group )
             pp ( 2 ) = REAL ( nn, dbl )
             nn = rs % npts_local ( dir )
             CALL mp_min ( nn, rs % group )
             pp ( 3 ) = REAL ( nn, dbl )
             IF ( rs % group_head ) THEN
                WRITE ( iounit, '(A,T48,A)' ) " RS_GRID:   Distribution", &
                     "  Average         Max         Min"
                WRITE ( iounit, '(A,T45,F12.1,2I12)' ) " RS_GRID:   Planes   ", &
                     pp ( 1 ), NINT ( pp ( 2 ) ), NINT ( pp ( 3 ) )
                WRITE ( iounit, '(/)' )
             END IF
          ELSE
             WRITE ( iounit, '(/,A,T71,I10)' ) &
                  " RS_GRID: Information for grid number ", rs % identifier
             DO i = 1, 3
                WRITE ( iounit, '(A,I3,T30,2I8,T62,A,T71,I10)' ) " RS_GRID:   Bounds ", &
                     i, rs % lb ( i ), rs % ub ( i ), "Points:", rs % npts ( i )
             END DO
             WRITE ( iounit, '(/)' )
          END IF
       END IF
    END IF

  END SUBROUTINE rs_grid_setup_1

  !******************************************************************************

  SUBROUTINE rs_grid_setup_n ( rs, pw_grid, nsmax, iounit )

    !  Arguments
    TYPE ( realspace_grid_type ), DIMENSION ( : ), INTENT (INOUT) :: rs
    TYPE ( pw_grid_type ), INTENT ( IN ) :: pw_grid
    INTEGER, INTENT ( IN ) :: nsmax 
    INTEGER, OPTIONAL, INTENT ( IN ) :: iounit

    !  Local
    INTEGER :: i

    !-----------------------------------------------------------------------------!

    CALL rs_grid_setup_1 ( rs ( 1 ), pw_grid, nsmax, iounit )

    DO i = 2, SIZE ( rs )

       CALL rs_grid_setup_1 ( rs ( i ), pw_grid, nsmax )

    END DO

  END SUBROUTINE rs_grid_setup_n

  !*****
  !******************************************************************************
!!****** realspace_grid_types/rs_pw_transfer [1.0] *
!!
!!   NAME
!!     rs_pw_transfer
!!
!!   SYNOPSIS
!!     Subroutine rs_pw_transfer(rs, pw, dir)
!!       Type(realspace_grid_type), Intent (INOUT):: rs
!!       Type(pw_type), Target, Intent (INOUT):: pw
!!       Character(Len=*), Intent (IN):: dir
!!     End Subroutine rs_pw_transfer
!!
!!   FUNCTION
!!     Copy a function from/to a PW grid type to/from a real
!!     space type 
!!
!!   AUTHOR
!!     JGH (18-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (15-Feb-2003) reduced additional memory usage
!!     Joost VandeVondele (Sep-2003) moved from sum/bcast to shift
!!
!!   SOURCE
  !******************************************************************************

  SUBROUTINE rs_pw_transfer ( rs, pw, dir )

    !  Arguments
    TYPE ( realspace_grid_type ), INTENT ( INOUT ) :: rs
    TYPE ( pw_type ), TARGET, INTENT ( INOUT ) :: pw
    CHARACTER ( LEN = * ), INTENT ( IN ) :: dir

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
    CALL timeset("rs_pw_transfer",'I',' ',handle)
    IF ( rs % identifier /= pw % pw_grid % identifier ) THEN
       CALL stop_program ( "rs_pw_transfer", &
            "Real space grid and pw type not compatible" )
    END IF

    IF ( pw % in_use == REALDATA1D .OR. &
         pw % in_use == COMPLEXDATA1D ) THEN
       CALL stop_program ( "rs_pw_transfer", "PW type not compatible" )
    END IF

    IF (  rs % parallel ) THEN

       IF ( rs % group_dim ( 1 ) == 1 ) THEN
        
          np = pw % pw_grid % para % group_size
          bo => pw % pw_grid % para % bo (1:2,1:3,0:np-1,1)
          pbo => pw % pw_grid % bounds
          group = pw % pw_grid % para % rs_group
          mepos = pw % pw_grid % para % rs_mpo 
          ALLOCATE ( rcount ( 0 : np - 1 ), STAT = ierr )
          IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer", "rcount", 2*np )
          DO ip = 1, np
             rcount ( ip-1 ) = PRODUCT ( bo(2,:,ip) - bo(1,:,ip) + 1 )
          END DO
          nma = MAXVAL ( rcount ( 0 : np - 1 ) )
          ALLOCATE ( rlocal ( nma ), STAT = ierr )
          IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer", "rlocal", nma )

          IF ( dir == "FORWARD" ) THEN
             rlocal=0.0_dbl
             DO ip = 1, np
                lb = pbo(1,:)+bo(1,:,MODULO(mepos-ip,np)+1)-1
                ub = pbo(1,:)+bo(2,:,MODULO(mepos-ip,np)+1)-1
                ii = 0
                DO iz = lb(3), ub(3)
                   DO iy = lb(2), ub(2)
                      DO ix = lb(1), ub(1)
                         ii=ii+1
                         rlocal(ii) = rlocal(ii)+rs%r(ix,iy,iz)
                      END DO
                   END DO
                END DO
                IF ( ip .EQ. np ) EXIT
                CALL mp_shift(rlocal,group,1)
             ENDDO
             nn = rcount(mepos)
             IF ( pw % in_use == REALDATA3D ) THEN
                CALL dcopy ( nn, rlocal, 1, pw % cr3d, 1 )
             ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
                CALL copy_rc ( nn, rlocal, pw % cc3d )
             ELSE
                CALL stop_program ( "rs_pw_transfer", "PW type not compatible" )
             END IF
          ELSEIF ( dir == "BACKWARD" ) THEN
             nn = rcount ( mepos )
             IF ( pw % in_use == REALDATA3D ) THEN
                CALL dcopy ( nn, pw % cr3d, 1, rlocal, 1 )
             ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
                CALL dcopy ( nn, pw % cc3d, 2, rlocal, 1 )
             ELSE
                CALL stop_program ( "rs_pw_transfer", "PW type not compatible" )
             END IF
             lb = pbo(1,:)+bo(1,:,mepos+1)-1
             ub = pbo(1,:)+bo(2,:,mepos+1)-1
             ii = 0
             DO iz = lb(3), ub(3)
                DO iy = lb(2), ub(2)
                   DO ix = lb(1), ub(1)
                      ii=ii+1
                      rs%r(ix,iy,iz) = rlocal(ii)
                   END DO
                END DO
             END DO
             DO ip = 1, np-1

                CALL mp_shift(rlocal,group,1)

                lb = pbo(1,:)+bo(1,:,MODULO(mepos-ip,np)+1)-1
                ub = pbo(1,:)+bo(2,:,MODULO(mepos-ip,np)+1)-1
                ii = 0
                DO iz = lb(3), ub(3)
                   DO iy = lb(2), ub(2)
                      DO ix = lb(1), ub(1)
                         ii=ii+1
                         rs%r(ix,iy,iz) = rlocal(ii)
                      END DO
                   END DO
                END DO
             ENDDO
          ELSE

             CALL stop_program ( "rs_pw_transfer", &
                  "Parameter dir ="//dir//" not allowed" )

          END IF

          DEALLOCATE ( rcount, STAT = ierr )
          IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer", "rcount" )
          DEALLOCATE ( rlocal, STAT = ierr )
          IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer", "rlocal" )

       ELSE
          !Both grids are distributed, we have to do some reshuffling of data

          IF ( dir == "FORWARD" ) THEN

             ! The border data has to be send/received from the neighbors
             ! First we calculate the source and destination processes for the shift
             ! The first shift is "downwards"
             CALL mp_cart_shift ( rs % group, 0, -1, source, dest )
             idir = rs % direction
             lb ( : ) = rs % lb_local ( : )
             ub ( : ) = rs % ub_local ( : )
             ub ( idir ) = lb ( idir ) + rs % border - 1 
             ! Allocate a scratch array to receive the data
             nn = PRODUCT ( ub - lb + 1 )
             ALLOCATE ( buffer ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ), STAT=ierr )
             IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","buffer",nn)
             CALL mp_sendrecv ( rs % r ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ), &
                  dest, buffer, source, rs % group )
             lb ( : ) = rs % lb_local ( : )
             ub ( : ) = rs % ub_local ( : )
             ub ( idir ) = ub ( idir ) - rs % border
             lb ( idir ) = ub ( idir ) - rs % border + 1
             ! Sum the data in the RS Grid
             rs % r ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ) = &
                  rs % r ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ) + buffer ( :, :, : )
             ! Now for the other direction
             CALL mp_cart_shift ( rs % group, 0, 1, source, dest )
             lb ( : ) = rs % lb_local ( : )
             ub ( : ) = rs % ub_local ( : )
             lb ( idir ) = ub ( idir ) - rs % border + 1 
             CALL mp_sendrecv ( rs % r ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ), &
                  dest, buffer, source, rs % group )
             lb ( : ) = rs % lb_local ( : )
             ub ( : ) = rs % ub_local ( : )
             lb ( idir ) = lb ( idir ) + rs % border
             ub ( idir ) = lb ( idir ) + rs % border - 1
             ! Sum the data in the RS Grid
             rs % r ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ) = &
                  rs % r ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ) + buffer ( :, :, : )
             ! Get rid of the scratch space
             DEALLOCATE ( buffer, STAT=ierr )
             IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","buffer")

             ! We have to sum up the data within the RS groups
             ! We generate a subgroup that covers all processors with the same RS grid
             subdim ( 1 ) = .FALSE.
             subdim ( 2 ) = .TRUE.
             CALL mp_cart_sub ( rs % group, subdim, subgroup )
             ! Now we do the sum (exclude border areas)
             lb ( : ) = rs % lb_local ( : )
             ub ( : ) = rs % ub_local ( : )
             lb ( idir ) = lb ( idir ) + rs % border
             ub ( idir ) = ub ( idir ) - rs % border
             CALL mp_sum ( rs % r ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ), subgroup )
             ! And we release the communicator again
             CALL mp_comm_free ( subgroup )

             ! This is the real redistribution
             CALL send_forward ( rs, pw )

          ELSEIF ( dir == "BACKWARD" ) THEN

             ! This is the real redistribution
             CALL send_backward ( rs, pw )
             ! Redistribution of the border area
             CALL mp_cart_shift ( rs % group, 0, -1, source, dest )
             idir = rs % direction
             lb ( : ) = rs % lb_local ( : )
             ub ( : ) = rs % ub_local ( : )
             lb ( idir ) = lb ( idir ) + rs % border
             ub ( idir ) = lb ( idir ) + rs % border - 1
             lc ( : ) = rs % lb_local ( : )
             uc ( : ) = rs % ub_local ( : )
             lc ( idir ) = uc ( idir ) - rs % border + 1
             nn = PRODUCT ( uc - lc + 1 )
             ALLOCATE ( buffer ( lc(1):uc(1), lc(2):uc(2), lc(3):uc(3) ), STAT=ierr )
             IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","buffer",nn)
             CALL mp_sendrecv ( rs % r ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ), &
                  dest, buffer, source, rs % group )
             rs % r ( lc(1):uc(1), lc(2):uc(2), lc(3):uc(3) ) = buffer ( :, :, : )
             ! Now for the other direction
             CALL mp_cart_shift ( rs % group, 0, 1, source, dest )
             lb ( : ) = rs % lb_local ( : )
             ub ( : ) = rs % ub_local ( : )
             ub ( idir ) = ub ( idir ) - rs % border
             lb ( idir ) = ub ( idir ) - rs % border + 1
             lc ( : ) = rs % lb_local ( : )
             uc ( : ) = rs % ub_local ( : )
             uc ( idir ) = lc ( idir ) + rs % border - 1
             CALL mp_sendrecv ( rs % r ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ), &
                  dest, buffer, source, rs % group )
             rs % r ( lc(1):uc(1), lc(2):uc(2), lc(3):uc(3) ) = buffer ( :, :, : )
             DEALLOCATE ( buffer, STAT=ierr )
             IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","buffer")

          ELSE
             CALL stop_program ( "rs_pw_transfer", &
                  "Parameter dir ="//dir//" not allowed" )
          END IF

       END IF

    ELSE

       ! non parallel case
       nn = SIZE ( rs % r )
       IF ( dir == "FORWARD" ) THEN
          IF ( pw % in_use == REALDATA3D ) THEN
             CALL dcopy ( nn, rs % r, 1, pw % cr3d, 1 )
          ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
             CALL copy_rc ( nn, rs % r, pw % cc3d )
          ELSE
             CALL stop_program ( "rs_pw_transfer", "Pw type not compatible" )
          END IF
       ELSEIF ( dir == "BACKWARD" ) THEN
          IF ( pw % in_use == REALDATA3D ) THEN
             CALL dcopy ( nn, pw % cr3d, 1, rs % r, 1 )
          ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
             CALL copy_cr ( nn, pw % cc3d, rs % r )
          ELSE
             CALL stop_program ( "rs_pw_transfer", "Pw type not compatible" )
          END IF
       ELSE
          CALL stop_program ( "rs_pw_transfer", &
               "Parameter dir ="//dir//" not allowed" )
       END IF

    END IF
    CALL timestop(0.0_dbl,handle)

  END SUBROUTINE rs_pw_transfer

!******************************************************************************

SUBROUTINE send_forward ( rs, pw )

!  Arguments
   TYPE ( realspace_grid_type ), INTENT ( INOUT ) :: rs
   TYPE ( pw_type ), INTENT ( INOUT ) :: pw

!  Local
   LOGICAL, DIMENSION ( 2 ) :: subdim
   INTEGER :: subgroup, ierr, nn, idir, il, iu, source
   INTEGER, DIMENSION ( 3 ) :: lb, ub, lbp, ubp, lz, uz
   REAL ( dbl ), DIMENSION ( :, :, : ), ALLOCATABLE :: buffer

!-----------------------------------------------------------------------------!

   lbp ( : ) = pw % pw_grid % bounds_local ( 1, : )
   ubp ( : ) = pw % pw_grid % bounds_local ( 2, : )
   ! We only have to do the communication within the subgroups along cartesian
   ! coordinates 1. These processors have the complete data available
   subdim ( 1 ) = .TRUE.
   subdim ( 2 ) = .FALSE.
   CALL mp_cart_sub ( rs % group, subdim, subgroup )
   idir = rs % direction
   DO source = 0, rs % group_dim ( 1 ) - 1
     ! First send information on bounds
     IF ( source == rs % group_coor ( 1 ) ) THEN
       lb = rs % lb_local
       ub = rs % ub_local
       lb ( idir ) = lb ( idir ) + rs % border
       ub ( idir ) = ub ( idir ) - rs % border
     END IF
     CALL mp_bcast ( lb, source, subgroup )
     CALL mp_bcast ( ub, source, subgroup )
     ! allocate scratch array to hold the data
     nn = PRODUCT ( ub - lb + 1 )
     ALLOCATE ( buffer ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "send_forward", "buffer", nn )
     ! send the data to all other processors 
     IF ( source == rs % group_coor ( 1 ) ) THEN
       buffer ( :,:,: ) = rs % r ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )
     END IF
     CALL mp_bcast ( buffer, source, subgroup )
     ! Pick the data needed on this processor
     ! We only have to consider the special direction, on the other direction
     ! we know that all the data is in buffer 
     lz = lbp
     uz = ubp
     lz ( idir ) = MAX ( lb ( idir ), lbp ( idir ) )
     uz ( idir ) = MIN ( ub ( idir ), ubp ( idir ) )
     IF ( pw % in_use == REALDATA3D ) THEN
       pw % cr3d ( lz(1):uz(1), lz(2):uz(2), lz(3):uz(3) ) = &
          buffer ( lz(1):uz(1), lz(2):uz(2), lz(3):uz(3) )
     ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
       pw % cc3d ( lz(1):uz(1), lz(2):uz(2), lz(3):uz(3) ) = &
          CMPLX ( buffer ( lz(1):uz(1), lz(2):uz(2), lz(3):uz(3) ), KIND = dbl )
     ELSE
       CALL stop_program ( "rs_pw_transfer", "Pw type not compatible" )
     END IF
     ! get rid of the scratch array
     DEALLOCATE ( buffer, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "send_forward", "buffer" )
   END DO
   ! Release the communicator
   CALL mp_comm_free ( subgroup )

END SUBROUTINE send_forward

!******************************************************************************

SUBROUTINE send_backward ( rs, pw )

!  Arguments
   TYPE ( realspace_grid_type ), INTENT ( INOUT ) :: rs
   TYPE ( pw_type ), INTENT ( INOUT ) :: pw

!  Local
   LOGICAL, DIMENSION ( 2 ) :: subdim
   INTEGER :: subgroup, ierr, nn, idir, il, iu, source
   INTEGER, DIMENSION ( 3 ) :: lb, ub, lbp, ubp, lz, uz
   REAL ( dbl ), DIMENSION ( :, :, : ), ALLOCATABLE :: buffer

!-----------------------------------------------------------------------------!

   ! initialize the rs grid to zero
   CALL rs_grid_zero( rs )
   lbp ( : ) = pw % pw_grid % bounds_local ( 1, : )
   ubp ( : ) = pw % pw_grid % bounds_local ( 2, : )
   ! We only have to do the communication within the subgroups along cartesian
   ! coordinates 1. These processors have the complete data available
   subdim ( 1 ) = .TRUE.
   subdim ( 2 ) = .FALSE.
   CALL mp_cart_sub ( rs % group, subdim, subgroup )
   idir = rs % direction
   DO source = 0, rs % group_dim ( 1 ) - 1
     ! First send information on bounds
     IF ( source == rs % group_coor ( 1 ) ) THEN
       lb = rs % lb_local
       ub = rs % ub_local
       lb ( idir ) = lb ( idir ) + rs % border
       ub ( idir ) = ub ( idir ) - rs % border
     END IF
     CALL mp_bcast ( lb, source, subgroup )
     CALL mp_bcast ( ub, source, subgroup )
     ! allocate scratch array to hold the data
     nn = PRODUCT ( ub - lb + 1 )
     ALLOCATE ( buffer ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ), STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "send_backward", "buffer", nn )
     buffer = 0._dbl
     ! Pick the data needed on this processor
     ! We only have to consider the special direction
     lz = lbp
     uz = ubp
     lz ( idir ) = MAX ( lb ( idir ), lbp ( idir ) )
     uz ( idir ) = MIN ( ub ( idir ), ubp ( idir ) )
     IF ( pw % in_use == REALDATA3D ) THEN
       buffer ( lz(1):uz(1), lz(2):uz(2), lz(3):uz(3) ) = &
          pw % cr3d ( lz(1):uz(1), lz(2):uz(2), lz(3):uz(3) )
     ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
       buffer ( lz(1):uz(1), lz(2):uz(2), lz(3):uz(3) ) = &
          REAL ( pw % cc3d ( lz(1):uz(1), lz(2):uz(2), lz(3):uz(3) ), dbl )
     ELSE
       CALL stop_program ( "rs_pw_transfer", "Pw type not compatible" )
     END IF
     ! collect the data from all other processors
     CALL mp_sum ( buffer, source, subgroup )
     IF ( source == rs % group_coor ( 1 ) ) THEN
       rs % r ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ) = &
         rs % r ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) ) + buffer ( :, :, : )
     END IF
     ! get rid of the scratch array
     DEALLOCATE ( buffer, STAT = ierr )
     IF ( ierr /= 0 ) CALL stop_memory ( "send_forward", "buffer" )
   END DO
   ! Release the communicator
   CALL mp_comm_free ( subgroup )
   ! We have to sum up the data in the second direction too
   subdim ( 1 ) = .FALSE.
   subdim ( 2 ) = .TRUE.
   CALL mp_cart_sub ( rs % group, subdim, subgroup )
   CALL mp_sum ( rs % r, subgroup )
   CALL mp_comm_free ( subgroup )

END SUBROUTINE send_backward

!*****
!******************************************************************************
!!****** realspace_grid_types/rs_grid_zero [1.0] *
!!
!!   NAME
!!     rs_grid_zero
!!
!!   SYNOPSIS
!!     Subroutine rs_grid_zero(rs)
!!       Type(realspace_grid_type), Intent (INOUT):: rs
!!     End Subroutine rs_grid_zero
!!
!!   FUNCTION
!!     Initialize grid to zero
!!
!!   AUTHOR
!!     JGH (23-Mar-2002)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!******************************************************************************

SUBROUTINE rs_grid_zero ( rs )

!  Arguments
   TYPE ( realspace_grid_type ), INTENT (INOUT) :: rs

!  Local
   INTEGER :: n,handle

!-----------------------------------------------------------------------------!

   CALL timeset("rs_grid_zero",'I',' ',handle)
   n = SIZE ( rs % r )
   CALL dcopy ( n, 0._dbl, 0, rs % r, 1 )
   CALL timestop(0.0_dbl,handle)

END SUBROUTINE rs_grid_zero
!******************************************************************************
SUBROUTINE rs_grid_to_cube ( rs, iunit, ionode, title )
!  Arguments
   TYPE ( realspace_grid_type ), INTENT (IN) :: rs
   INTEGER, INTENT(IN) :: iunit
   LOGICAL, INTENT(IN) :: ionode
   CHARACTER(*), OPTIONAL :: title
   INTEGER :: I1,I2,I3
   INTEGER :: S1,S2,S3
   INTEGER :: L1,L2,L3,U1,U2,U3,N1,N2,N3
!-----------------------------------------------------------------------------!
   ! not working for distributed grids !
   IF (rs%ngpts.ne.rs%ngpts_local) CALL stop_program("rs_grid_to_cube","NYI")
   IF (ionode) THEN
     ! this format seems to work for e.g. molekel and gOpenmol
     WRITE(iunit,'(a11)') "-Quickstep-"
     IF (PRESENT(title)) THEN
        WRITE(iunit,*) title
     ELSE
        WRITE(iunit,*) "No Title"
     ENDIF
     ! put 0.0, 0.0 ,0.0 (approx) in the center of the cube
     S1=rs%npts(1)/2
     S2=rs%npts(2)/2
     S3=rs%npts(3)/2
     WRITE(iunit,'(I5,3f12.6)') 0,0.0_dbl-S1*rs%dr(1),0.0_dbl-S2*rs%dr(2),0.0_dbl-S3*rs%dr(3)
     WRITE(iunit,'(I5,3f12.6)') rs%npts(1),rs%dr(1),0.0_dbl,0.0_dbl
     WRITE(iunit,'(I5,3f12.6)') rs%npts(2),0.0_dbl,rs%dr(2),0.0_dbl
     WRITE(iunit,'(I5,3f12.6)') rs%npts(3),0.0_dbl,0.0_dbl,rs%dr(3)
     ! write(iunit,'(I5,4f12.6)') 1, 0.0,0.0,0.0,0.0
     L1=LBOUND(rs%r,1)
     L2=LBOUND(rs%r,2)
     L3=LBOUND(rs%r,3)
     U1=UBOUND(rs%r,1)
     U2=UBOUND(rs%r,2)
     U3=UBOUND(rs%r,3)
     DO I1=L1,U1
       DO I2=L2,U2
        WRITE(iunit,'(6E13.5)') (rs%r(MODULO(I1-S1-L1,U1-L1+1)+L1,MODULO(I2-S2-L2,U2-L2+1)+L2, &
                                      MODULO(I3-S3-L3,U3-L3+1)+L3),I3=L3,U3)
       ENDDO
     ENDDO
   ENDIF
END SUBROUTINE rs_grid_to_cube
!******************************************************************************

!******************************************************************************
SUBROUTINE rs_pw_to_cube ( pw, iunit, ionode, title )
 !  Arguments
   TYPE ( pw_type), INTENT ( INOUT )        :: pw
   INTEGER, INTENT(IN)                      :: iunit
   LOGICAL, INTENT(IN)                      :: ionode
   CHARACTER(*), INTENT(IN), OPTIONAL       :: title
 !-----------------------------------------------------------------------------!
   TYPE ( realspace_grid_type )   :: rs
   INTEGER                        :: nsmax
 !-----------------------------------------------------------------------------!
   nsmax=-1
   CALL rs_grid_setup(rs,pw%pw_grid,nsmax)
   CALL rs_grid_allocate(rs)
   CALL rs_pw_transfer(rs,pw,"BACKWARD")
   IF (PRESENT(title)) THEN
      CALL rs_grid_to_cube(rs,iunit,ionode,title)
   ELSE
      CALL rs_grid_to_cube(rs,iunit,ionode)
   ENDIF
   CALL rs_grid_deallocate(rs)
END SUBROUTINE rs_pw_to_cube
!******************************************************************************

!!****f* realspace_grid_types/rs_grid_write [1.0] *
!!
!!   NAME
!!     rs_grid_write
!!
!!   SYNOPSIS
!!     Subroutine rs_grid_write(rs, unit_nr, error)
!!       Type(realspace_grid_type), Pointer:: rs
!!       Integer, Intent (IN):: unit_nr
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine rs_grid_write
!!
!!   FUNCTION
!!     writes out information about the given rs grid
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rs: the realspace grid to output
!!     - unit_nr: the unit number to output to
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     03.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE rs_grid_write(rs,unit_nr,error)
  TYPE(realspace_grid_type), POINTER :: rs
  INTEGER, INTENT(in) :: unit_nr
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='rs_grid_write',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  
  CALL cp_unimplemented_error(fromWhere=routineP, message="at "//&
       CPSourceFileRef,&
       error=error)
END SUBROUTINE rs_grid_write
!***************************************************************************

!!****f* realspace_grid_types/rs_grid_create [1.0] *
!!
!!   NAME
!!     rs_grid_create
!!
!!   SYNOPSIS
!!     Subroutine rs_grid_create(rs_grid, pw_grid, nsmax, error)
!!       Type(realspace_grid_type), Pointer:: rs_grid
!!       Type(pw_grid_type), Pointer:: pw_grid
!!       Integer, Optional, Intent (IN):: nsmax
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine rs_grid_create
!!
!!   FUNCTION
!!     allocates and initialize an rs grid
!!
!!   NOTES
!!     nsmax removed
!!
!!   ARGUMENTS
!!     - rs_grid: the rs grid that is created
!!     - pw_grid: the corresponding pw_grid
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     03.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE rs_grid_create(rs_grid, pw_grid, nsmax, error)
  TYPE(realspace_grid_type), POINTER :: rs_grid
  TYPE ( pw_grid_type ), POINTER :: pw_grid
  INTEGER, OPTIONAL, INTENT(in) :: nsmax
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='rs_grid_create',&
        routineP=moduleN//':'//routineN
  TYPE(cp_logger_type), POINTER :: logger
  INTEGER :: stat, iounit
  INTEGER :: my_nsmax

  failure=.FALSE.
  logger => cp_error_get_logger(error)
  my_nsmax = -1
  IF (PRESENT(nsmax)) my_nsmax=nsmax
  
  CPPrecondition(ASSOCIATED(pw_grid),cp_failure_level,routineP,error,failure)
  
  ALLOCATE(rs_grid,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
     IF (logger%print_keys%pw_grid_information) THEN
        CALL rs_grid_setup(rs_grid, pw_grid, my_nsmax, &
             iounit=cp_logger_get_default_unit_nr(logger,local=.FALSE.))
     ELSE
        CALL rs_grid_setup(rs_grid, pw_grid, my_nsmax)
     END IF
     CALL rs_grid_allocate(rs_grid)
  END IF
END SUBROUTINE rs_grid_create
!***************************************************************************

!!****f* realspace_grid_types/rs_grid_retain [1.0] *
!!
!!   NAME
!!     rs_grid_retain
!!
!!   SYNOPSIS
!!     Subroutine rs_grid_retain(rs_grid, error)
!!       Type(realspace_grid_type), Pointer:: rs_grid
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine rs_grid_retain
!!
!!   FUNCTION
!!     retains the given rs grid (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rs_grid: the grid to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     03.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE rs_grid_retain(rs_grid, error)
  TYPE(realspace_grid_type), POINTER :: rs_grid
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='rs_grid_retain',&
        routineP=moduleN//':'//routineN

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(rs_grid),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPreconditionNoFail(rs_grid%ref_count>0,cp_failure_level,routineP,error)
     rs_grid%ref_count=rs_grid%ref_count+1
  END IF
END SUBROUTINE rs_grid_retain
!***************************************************************************

!!****f* realspace_grid_types/rs_grid_release [1.0] *
!!
!!   NAME
!!     rs_grid_release
!!
!!   SYNOPSIS
!!     Subroutine rs_grid_release(rs_grid, error)
!!       Type(realspace_grid_type), Pointer:: rs_grid
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine rs_grid_release
!!
!!   FUNCTION
!!     releases the given rs grid (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rs_grid: the rs grid to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     03.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE rs_grid_release(rs_grid,error)
  TYPE(realspace_grid_type), POINTER :: rs_grid
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='rs_grid_release',&
        routineP=moduleN//':'//routineN
  INTEGER :: stat

  failure=.FALSE.
  
  IF (ASSOCIATED(rs_grid)) THEN
     CPPreconditionNoFail(rs_grid%ref_count>0,cp_failure_level,routineP,error)
     rs_grid%ref_count=rs_grid%ref_count-1
     IF (rs_grid%ref_count==0) THEN
        rs_grid%ref_count=1
        CALL rs_grid_deallocate(rs_grid)
        rs_grid%ref_count=0
        DEALLOCATE(rs_grid, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(rs_grid)
END SUBROUTINE rs_grid_release

!***************************************************************************
!!****f* realspace_grid_types/rs_get_my_tasks [1.0] *
!!
!!   NAME
!!     rs_get_my_tasks
!!
!!   SYNOPSIS
!!     Subroutine rs_get_my_tasks(rs, tasks, npme, ival, rval, pmat, pcor)
!!       Type(realspace_grid_type), Intent (IN):: rs
!!       Integer, Dimension(:,:), Pointer:: tasks
!!       Integer, Intent (OUT):: npme
!!       Integer, Dimension(:,:), Pointer, Optional:: ival
!!       Real(Kind=dbl), Dimension(:,:), Pointer, Optional:: rval
!!       Type(real_matrix_type), Pointer, Optional:: pmat
!!       Integer, Dimension(:,:), Pointer, Optional:: pcor
!!     End Subroutine rs_get_my_tasks
!!
!!   FUNCTION
!!     Assembles tasks to be performed on local grid
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - rs: the grid 
!!     - tasks: the task set generate on this processor
!!     - tasks_local: the task set to be processed localy
!!     - npme: Number of tasks for local processing
!!
!!   AUTHOR
!!     JGH (19.06.2003)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** **********************************************************************

SUBROUTINE rs_get_my_tasks ( rs, tasks, npme, ival, rval, pmat, pcor )

  TYPE ( realspace_grid_type ), INTENT ( IN ) :: rs
  INTEGER, DIMENSION ( :, : ), POINTER :: tasks
  INTEGER, INTENT ( OUT ) :: npme
  INTEGER, DIMENSION ( :, : ), POINTER, OPTIONAL :: ival
  REAL(dbl), DIMENSION ( :, : ), POINTER, OPTIONAL :: rval
  TYPE ( real_matrix_type ), POINTER, OPTIONAL :: pmat
  INTEGER, DIMENSION ( :, : ), POINTER, OPTIONAL :: pcor

  INTEGER, DIMENSION ( :, : ), POINTER :: tlist, ilist, plist, pls
  REAL(dbl), DIMENSION ( :, : ), POINTER :: rlist
  REAL(dbl), DIMENSION (:,:), POINTER :: p_block
  REAL(dbl), DIMENSION (:), POINTER :: ppack
  INTEGER :: subgroup, subsize, subpos
  INTEGER :: isend, ileft, i, j, k, lb, ub, n1, stat, npmax, ic, acol, arow
  INTEGER :: icmax, plength, nppp, handle
  LOGICAL, DIMENSION ( 2 ) :: subdim
  LOGICAL :: matrix
 
  CALL timeset("rs_get_my_tasks",'I',' ',handle)

  IF ( .NOT. ASSOCIATED ( tasks ) ) &
        CALL stop_program ( "get_my_tasks", "tasks not associated" )

  IF ( PRESENT ( pmat ) .AND. PRESENT ( pcor ) ) THEN
    matrix = .TRUE.
  ELSE
    matrix = .FALSE.
  END IF

  ! get number of tasks available locally
  npme = SIZE ( tasks, 2 )
  DO i = 1, SIZE ( tasks, 2)
    IF ( tasks ( 1, i ) <= 0 ) THEN
      npme = i - 1
      EXIT
    END IF
  END DO

  IF ( rs%parallel .AND. rs%direction > 0 ) THEN

    ! bounds of local grid
    lb = rs%lb_local ( rs%direction ) + rs%border
    ub = rs%ub_local ( rs%direction ) - rs%border

    ! We need subgroups along cartesian coordinate 1.
    subdim ( 1 ) = .TRUE.
    subdim ( 2 ) = .FALSE.
    CALL mp_cart_sub ( rs % group, subdim, subgroup )
    CALL mp_environ ( subsize, subpos, subgroup )

    ! determine maximum number of tasks 
    npmax = npme
    CALL mp_max ( npmax, subgroup )
    nppp = npmax

    ! allocate local arrays
    ALLOCATE(tlist(2,npmax),STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","tlist",2*npmax)
    tlist(1:2,1:npme)=tasks(1:2,1:npme)
    tlist(1:2,npme+1:npmax)=0
    IF ( PRESENT ( ival ) ) THEN
      n1 = SIZE ( ival, 1 )
      ALLOCATE(ilist(n1,npmax),STAT=stat)
      IF (stat/=0) CALL stop_memory("get_my_tasks","ilist",n1*npmax)
      ilist(1:n1,1:npme)=ival(1:n1,1:npme)
      ilist(1:n1,npme+1:npmax)=0
      ival=0
    END IF
    IF ( PRESENT ( rval ) ) THEN
      n1 = SIZE ( rval, 1 )
      ALLOCATE(rlist(n1,npmax),STAT=stat)
      IF (stat/=0) CALL stop_memory("get_my_tasks","rlist",n1*npmax)
      rlist(1:n1,1:npme)=rval(1:n1,1:npme)
      rlist(1:n1,npme+1:npmax)=0._dbl
      rval=0._dbl
    END IF
    IF ( matrix ) THEN
      n1 = SIZE ( pcor, 1 )
      ALLOCATE(plist(n1,npmax),STAT=stat)
      IF (stat/=0) CALL stop_memory("get_my_tasks","plist",n1*npmax)
      plist(1:n1,1:npme)=pcor(1:n1,1:npme)
      plist(1:n1,npme+1:npmax)=0
      ALLOCATE(pls(n1,npmax),STAT=stat)
      IF (stat/=0) CALL stop_memory("get_my_tasks","pls",n1*npmax)
      pls = 0
      ALLOCATE(ppack(npmax),STAT=stat)
      IF (stat/=0) CALL stop_memory("get_my_tasks","ppack",npmax)
    END IF

    ! keep tasks to be processed locally
    ! send remaining tasks to the next processor
    npme = 0
    DO isend = 0, subsize - 1
      IF ( isend == 0 ) THEN
        ileft = SIZE ( tlist, 2)
        DO i = 1, SIZE ( tlist, 2)
          IF ( tlist ( 1, i ) <= 0 ) THEN
            ileft = i - 1
            EXIT
          END IF
        END DO
        npmax = ileft
        CALL mp_max ( npmax, subgroup )
      ELSE
        ! count left over tasks
        j = 0
        DO i = 1, npmax
          IF ( tlist ( 1, i ) > 0 ) THEN
            j = j + 1
            tlist(:,j) = tlist(:,i)
            IF ( PRESENT ( ival ) ) ilist(:,j) = ilist(:,i)
            IF ( PRESENT ( rval ) ) rlist(:,j) = rlist(:,i)
            IF ( matrix ) plist(:,j) = plist(:,i)
          END IF
        END DO
        tlist(:,j+1:npmax) = 0
        ileft = j
        npmax = ileft
        CALL mp_max ( npmax, subgroup )
        IF ( npmax > nppp ) &
          CALL stop_program ( "get_my_tasks","npmax too large")
        IF ( npmax == 0 ) EXIT
        ! send/receive tasks
        CALL mp_shift ( tlist ( 1:2, 1:npmax ), subgroup )
        IF ( PRESENT ( ival ) ) THEN
           n1 = SIZE ( ilist, 1 )
           CALL mp_shift ( ilist ( 1:n1, 1:npmax ), subgroup )
        END IF
        IF ( PRESENT ( rval ) ) THEN
           n1 = SIZE ( rlist, 1 )
           CALL mp_shift ( rlist ( 1:n1, 1:npmax ), subgroup )
        END IF
        IF ( matrix ) THEN
           n1 = SIZE ( plist, 1 )
           CALL mp_shift ( plist ( 1:n1, 1:npmax ), subgroup )
        END IF
      END IF
      ! look for tasks do be done on this processor
      ic = 0
      DO i = 1, npmax
        IF ( tlist(1,i) > 0 ) THEN
          IF ( tlist(2,i) >= lb .AND. tlist(2,i) <= ub ) THEN
            ! found new local task
            npme = npme + 1
            IF ( npme > SIZE ( tasks, 2 ) ) &
               CALL reallocate ( tasks, 1, 2, 1, 2*npme )
            tasks ( :, npme ) = tlist ( :, i )
            tlist ( :, i ) = 0
            IF ( PRESENT ( ival ) ) THEN
              IF ( npme > SIZE ( ival, 2 ) ) THEN
                n1 = SIZE ( ival, 1 )
                CALL reallocate ( ival, 1, n1, 1, 2*npme )
              END IF
              ival ( :, npme ) = ilist ( :, i )
            END IF
            IF ( PRESENT ( rval ) ) THEN
              IF ( npme > SIZE ( rval, 2 ) ) THEN
                n1 = SIZE ( rval, 1 )
                CALL reallocate ( rval, 1, n1, 1, 2*npme )
              END IF
              rval ( :, npme ) = rlist ( :, i )
            END IF
            IF ( matrix ) THEN
              ic = ic + 1
              IF ( ic > SIZE ( pcor, 2 ) ) THEN
                n1 = SIZE ( pcor, 1 )
                CALL reallocate ( pcor, 1, n1, 1, 2*ic )
              END IF
              pcor ( :, ic ) = plist ( :, i )
            END IF
          END IF
        END IF
      END DO
      ! now distribute matrix blocks
      IF ( matrix .AND. isend /= 0 ) THEN
        ! eliminate double entries
        IF ( ic > 1 ) THEN
          j=1
          DO i = 2, ic
            DO k = j, 1, -1
              IF ( pcor(1,i)==pcor(1,k) .AND. pcor(2,i)==pcor(2,k) ) EXIT
              IF ( k == 1 ) THEN
                j = j + 1 
                pcor(:,j)=pcor(:,i)
              END IF
            END DO
          END DO
          ic = j
        END IF
        ! eliminate blocks that are already local
        j = 0
        DO i = 1, ic
          IF ( pcor(1,i) <= pcor(2,i) ) THEN
             arow = pcor(1,i)
             acol = pcor(2,i)
          ELSE
             arow = pcor(2,i)
             acol = pcor(1,i)
          END IF
          CALL get_block_node(matrix=pmat,&
                              block_row=arow,&
                              block_col=acol,&
                              BLOCK=p_block)
          IF ( .NOT. ASSOCIATED ( p_block ) ) THEN
            j = j + 1
            pcor(:,j) = pcor(:,i)
          END IF
        END DO
        ic = j
        icmax = ic
        CALL mp_max ( icmax, subgroup )
        ! are there missing blocks?
        IF ( icmax > 0 ) THEN
          ! send block coordinates back to original pe
          pls(:,1:ic) = pcor(:,1:ic)
          pls(:,ic+1:) = 0
          CALL mp_shift ( pls ( 1:2, 1:icmax ), subgroup, -isend )
          ! pack the data to be sent
          CALL pack_matrix ( pmat, pls, ppack, plength, subgroup )
          CALL mp_shift ( ppack ( 1:plength ), subgroup, isend )
          CALL mp_shift ( pls ( 1:2, 1:icmax ), subgroup, isend )
          CALL unpack_matrix ( pmat, pcor, pls, ic, ppack )
        END IF
      END IF
    END DO
    ! are all tasks distributed?
    IF ( ANY ( tlist(1,:) > 0 ) ) &
          CALL stop_program ( "get_my_tasks", "left over tasks" )

    DEALLOCATE(tlist,STAT=stat)
    IF (stat/=0) CALL stop_memory("get_my_tasks","tlist")
    IF ( PRESENT ( ival ) ) THEN
      DEALLOCATE(ilist,STAT=stat)
      IF (stat/=0) CALL stop_memory("get_my_tasks","ilist")
    END IF
    IF ( PRESENT ( rval ) ) THEN
      DEALLOCATE(rlist,STAT=stat)
      IF (stat/=0) CALL stop_memory("get_my_tasks","rlist")
    END IF
    IF ( matrix ) THEN
      DEALLOCATE(plist,STAT=stat)
      IF (stat/=0) CALL stop_memory("get_my_tasks","plist")
      DEALLOCATE(pls,STAT=stat)
      IF (stat/=0) CALL stop_memory("get_my_tasks","pls")
      DEALLOCATE(ppack,STAT=stat)
      IF (stat/=0) CALL stop_memory("get_my_tasks","ppack")
    END IF

    ! Release the communicator
    CALL mp_comm_free ( subgroup )

  ELSE

    ! fully replicated grids, each processor can process all its tasks

  END IF

  CALL timestop(0.0_dbl,handle)

END SUBROUTINE rs_get_my_tasks

!***************************************************************************
SUBROUTINE pack_matrix ( pmat, pls, ppack , plength, group )

  TYPE ( real_matrix_type ), POINTER, OPTIONAL :: pmat
  INTEGER, DIMENSION ( :, : ), POINTER         :: pls
  REAL(dbl), DIMENSION ( : ), POINTER          :: ppack
  INTEGER, INTENT ( OUT )                      :: plength
  INTEGER, INTENT ( IN )                       :: group

  REAL(dbl), DIMENSION (:,:), POINTER :: p_block
  INTEGER :: i, j, k, arow, acol, nc, nr, n

  plength = 0
  DO i = 1, SIZE ( pls, 2 )
    IF ( pls ( 1, i ) == 0 ) EXIT
    IF ( pls(1,i) <= pls(2,i) ) THEN
       arow = pls(1,i)
       acol = pls(2,i)
    ELSE
       arow = pls(2,i)
       acol = pls(1,i)
    END IF
    CALL get_block_node(matrix=pmat,&
                        block_row=arow,&
                        block_col=acol,&
                        BLOCK=p_block)
    IF ( .NOT. ASSOCIATED ( p_block ) ) &
       CALL stop_program ( "pack_matrix", "Matrix block not found" )
    nr = SIZE ( p_block, 1 )
    nc = SIZE ( p_block, 2 )
    pls ( 1, i ) = nr
    pls ( 2, i ) = nc
    n = nc * nr
    IF ( plength + n > SIZE ( ppack ) ) THEN
      CALL reallocate ( ppack, 1, plength+5*n )
    END IF
    DO j = 1, nc
      DO k = 1, nr
        plength = plength + 1
        ppack(plength) = p_block(k,j)
      END DO
    END DO
  END DO
  CALL mp_max ( plength, group )
  IF ( plength > SIZE ( ppack ) ) THEN
    CALL reallocate ( ppack, 1, plength )
  END IF

END SUBROUTINE pack_matrix 

SUBROUTINE unpack_matrix ( pmat, pcor, pls, ic, ppack )

  TYPE ( real_matrix_type ), POINTER, OPTIONAL :: pmat
  INTEGER, DIMENSION ( :, : ), POINTER         :: pcor
  INTEGER, DIMENSION ( :, : ), POINTER         :: pls
  INTEGER, INTENT ( IN )                       :: ic
  REAL(dbl), DIMENSION ( : ), POINTER          :: ppack

  REAL(dbl), DIMENSION (:,:), POINTER :: p_block
  INTEGER :: i, j, k, nc, nr, pl, arow, acol, stat

  pl = 0
  DO i = 1, ic
    IF ( pcor(1,i) <= pcor(2,i) ) THEN
       arow = pcor(1,i)
       acol = pcor(2,i)
    ELSE
       arow = pcor(2,i)
       acol = pcor(1,i)
    END IF
    nr = pls(1,i)
    nc = pls(2,i)
    NULLIFY ( p_block )
    CALL add_block_node ( pmat, arow, acol, p_block )
    DO j = 1, nc
      DO k = 1, nr
        pl = pl + 1
        p_block(k,j) = ppack(pl)
      END DO
    END DO

  END DO

END SUBROUTINE unpack_matrix

!***************************************************************************

SUBROUTINE rs_get_loop_vars ( npme, ival, natom_pairs, asets, atasks )

     INTEGER, INTENT(IN) :: npme
     INTEGER, INTENT(OUT) :: natom_pairs
     INTEGER, DIMENSION(6,npme), INTENT(IN) :: ival
     INTEGER, DIMENSION(:,:), POINTER :: asets, atasks

     INTEGER :: itask, iatom, jatom, iatom_old, jatom_old, &
                nset_pairs, iset, jset, iset_old, jset_old

     IF(SIZE(asets,2) < npme) CALL reallocate(asets,1,2,1,npme)
     IF(SIZE(atasks,2) < npme) CALL reallocate(atasks,1,2,1,npme)
     natom_pairs = 0
     nset_pairs = 0
     iatom_old = 0
     jatom_old = 0
     DO itask = 1, npme
        iatom = ival(1,itask)
        jatom = ival(2,itask)
        iset = ival(3,itask)
        jset = ival(4,itask)
        IF ( iatom /= iatom_old .OR. jatom /= jatom_old ) THEN
          IF(natom_pairs>0) asets(2,natom_pairs) = nset_pairs
          natom_pairs = natom_pairs + 1
          iatom_old = iatom
          jatom_old = jatom
          IF(nset_pairs>0) atasks(2,nset_pairs) = itask - 1
          nset_pairs = nset_pairs + 1
          asets(1,natom_pairs) = nset_pairs
          atasks(1,nset_pairs) = itask
          iset_old = iset
          jset_old = jset
        ELSE IF ( iset /= iset_old .OR. jset /= jset_old ) THEN
          atasks(2,nset_pairs) = itask - 1
          nset_pairs = nset_pairs + 1
          atasks(1,nset_pairs) = itask
          iset_old = iset
          jset_old = jset
        END IF
     END DO
     IF(natom_pairs>0) asets(2,natom_pairs) = nset_pairs
     IF(nset_pairs>0) atasks(2,nset_pairs) = npme

END SUBROUTINE rs_get_loop_vars

!***************************************************************************

END MODULE realspace_grid_types

!******************************************************************************
