!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2001 - 2003  CP2K developers group                          !
!-----------------------------------------------------------------------------!


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
  USE atomic_kind_types,               ONLY: get_atomic_kind
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE memory_utilities,                ONLY: reallocate
  USE message_passing,                 ONLY: &
       mp_bcast, mp_cart_create, mp_cart_shift, mp_cart_sub, mp_comm_dup, &
       mp_comm_free, mp_environ, mp_isendrecv, mp_max, mp_min, mp_sendrecv, &
       mp_shift, mp_sum, mp_waitall
  USE particle_list_types,             ONLY: particle_list_type
  USE pw_grid_types,                   ONLY: PW_MODE_LOCAL,&
                                             pw_grid_type
  USE pw_types,                        ONLY: COMPLEXDATA3D,&
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
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: realspace_grid_type,&
            realspace_grid_p_type

  PUBLIC :: rs_pw_transfer,&
            rs_grid_zero,&
            rs_grid_to_cube, &
            rs_pw_to_cube,&
            rs_get_my_tasks,&
            rs_get_loop_vars,&
            rs_grid_create,&
            rs_grid_retain,&
            rs_grid_release,&
            rs_grid_p_create,&
            rs_grid_p_release

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN="realspace_grid_types"
  INTEGER, SAVE, PRIVATE :: last_rs_id=0
  INTEGER, SAVE, PRIVATE :: allocated_rs_grid_count=0
  INTEGER, PARAMETER, PUBLIC :: rs2pw=11,pw2rs=12

  TYPE realspace_grid_type
     INTEGER :: grid_id                               ! tag of the pw_grid
     INTEGER :: id_nr                                    ! unique identifier of rs
     INTEGER :: ref_count                                ! reference count
     INTEGER :: ngpts                                    ! # grid points
     INTEGER, DIMENSION (3) :: npts                      ! # grid points per dimension
     INTEGER, DIMENSION (3) :: lb                        ! lower bounds
     INTEGER, DIMENSION (3) :: ub                        ! upper bounds
     INTEGER, DIMENSION (:), POINTER :: px,py,pz         ! index translators
     REAL(KIND=dp), DIMENSION ( :, :, : ),POINTER :: r    ! the grid
     INTEGER :: border                                   ! border points
     INTEGER :: ngpts_local                              ! local dimensions
     INTEGER, DIMENSION (3) :: npts_local
     INTEGER, DIMENSION (3) :: lb_local
     INTEGER, DIMENSION (3) :: ub_local
     REAL(KIND=dp), DIMENSION(3) :: dr                      ! grid spacing
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

  INTERFACE rs_grid_deallocate
     MODULE PROCEDURE rs_grid_deallocate_1, rs_grid_deallocate_n
  END INTERFACE

  INTERFACE rs_grid_setup
     MODULE PROCEDURE rs_grid_setup_1, rs_grid_setup_n
  END INTERFACE


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

  SUBROUTINE rs_grid_allocate_1 ( rs, error )

    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "rs_grid_allocate_1", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr
    INTEGER, DIMENSION(:), POINTER           :: lb, ub

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

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE rs_grid_allocate_1

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

  SUBROUTINE rs_grid_deallocate_1 ( rs, error )

    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "rs_grid_deallocate_1", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ierr, stat
    LOGICAL                                  :: failure

    failure = .FALSE.
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
    rs % grid_id = 0
    IF ( rs % parallel ) THEN
       ! release the group communicator
       CALL mp_comm_free ( rs % group )
    END IF
 
    CALL timestop(0.0_dp,handle)    
  END SUBROUTINE rs_grid_deallocate_1

  !******************************************************************************

  SUBROUTINE rs_grid_deallocate_n ( rs, error )

    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "rs_grid_deallocate_n", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, n
    LOGICAL                                  :: failure

    failure = .FALSE.

    n = SIZE ( rs )
    DO i = 1, n
       CALL rs_grid_deallocate_1 ( rs ( i )%rs_grid, error )
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

  SUBROUTINE rs_grid_setup_1 ( rs, pw_grid, nsmax, should_output, iounit, error)
    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(pw_grid_type), INTENT(IN)           :: pw_grid
    INTEGER, INTENT(IN)                      :: nsmax
    LOGICAL, INTENT(IN)                      :: should_output
    INTEGER, INTENT(IN)                      :: iounit
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "rs_grid_setup_1", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dir, i, ip, lb( 2 ), &
                                                n_slices, ngmax, nmin, nn, &
                                                pos( 2 )
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: pp( 3 )

    failure = .FALSE.

    rs % dr = pw_grid%dr
    rs % grid_id = pw_grid % id_nr
    last_rs_id = last_rs_id+1
    rs % id_nr = last_rs_id
    rs % ref_count = 1
    IF ( pw_grid % para % mode == PW_MODE_LOCAL ) THEN
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
          rs % group_head = ALL ( pos == 0 ) ! for output would be needed to have this pw_grid % para % group_head
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

    IF ( should_output ) THEN
          IF ( rs % parallel ) THEN
             IF ( iounit>0) THEN
                WRITE ( iounit, '(/,A,T71,I10)' ) &
                     " RS_GRID: Information for grid number ", rs % grid_id 
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
             pp ( 1 ) = REAL ( nn, KIND=dp ) / REAL ( PRODUCT ( rs % group_dim ), KIND=dp )
             nn = rs % npts_local ( dir )
             CALL mp_max ( nn, rs % group )
             pp ( 2 ) = REAL ( nn, KIND=dp )
             nn = rs % npts_local ( dir )
             CALL mp_min ( nn, rs % group )
             pp ( 3 ) = REAL ( nn, KIND=dp )
             IF ( iounit > 0) THEN
                WRITE ( iounit, '(A,T48,A)' ) " RS_GRID:   Distribution", &
                     "  Average         Max         Min"
                WRITE ( iounit, '(A,T45,F12.1,2I12)' ) " RS_GRID:   Planes   ", &
                     pp ( 1 ), NINT ( pp ( 2 ) ), NINT ( pp ( 3 ) )
                WRITE ( iounit, '(/)' )
             END IF
          ELSE
             IF (iounit>0) THEN
               WRITE ( iounit, '(/,A,T71,I10)' ) &
                    " RS_GRID: Information for grid number ", rs % grid_id
               DO i = 1, 3
                  WRITE ( iounit, '(A,I3,T30,2I8,T62,A,T71,I10)' ) " RS_GRID:   Bounds ", &
                       i, rs % lb ( i ), rs % ub ( i ), "Points:", rs % npts ( i )
               END DO
               WRITE ( iounit, '(/)' )
             ENDIF
          END IF
    END IF
    CALL rs_grid_allocate_1 ( rs, error )

  END SUBROUTINE rs_grid_setup_1

  !******************************************************************************

  SUBROUTINE rs_grid_setup_n ( rs, pw_grid, nsmax, should_output, iounit, error)

    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs
    TYPE(pw_grid_type), INTENT(IN)           :: pw_grid
    INTEGER, INTENT(IN)                      :: nsmax
    LOGICAL, INTENT(IN)                      :: should_output
    INTEGER, INTENT(IN)                      :: iounit
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = "rs_grid_setup_n", &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

    failure = .FALSE.

    CALL rs_grid_setup_1 ( rs ( 1 ) % rs_grid , pw_grid, nsmax, should_output, iounit )
    DO i = 2, SIZE ( rs )
       CALL rs_grid_setup_1 ( rs ( i ) % rs_grid , pw_grid, nsmax, .FALSE., 0 )
    END DO

  END SUBROUTINE rs_grid_setup_n

!!****** realspace_grid_types/rs_pw_transfer [1.0] *
!!
!!   NAME
!!     rs_pw_transfer
!!
!!   FUNCTION
!!     Copy a function from/to a PW grid type to/from a real
!!     space type 
!!     dir is the named constant rs2pw or pw2rs
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

    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(pw_type), POINTER                   :: pw
    INTEGER, INTENT(IN)                      :: dir

    CHARACTER(LEN=5)                         :: dirname
    INTEGER                                  :: handle, handle2, nn

!-----------------------------------------------------------------------------!

    CALL timeset("rs_pw_transfer_all",'I',' ',handle2)
    IF (dir.EQ.rs2pw) dirname="RS2PW"
    IF (dir.EQ.pw2rs) dirname="PW2RS"
    CALL timeset("rs_pw_transfer_"//TRIM(dirname)//"_"// TRIM(ADJUSTL(&
              cp_to_string(CEILING(pw%pw_grid%cutoff/10)*10))),'I',' ',handle)

    IF ( rs % grid_id /= pw % pw_grid % id_nr ) &
       CALL stop_program ( "rs_pw_transfer", "different rs and pw indentitfiers")
    IF ( (pw%in_use .NE. REALDATA3D) .AND. (pw%in_use .NE. COMPLEXDATA3D)) &
       CALL stop_program ( "rs_pw_transfer", "Need REALDATA3D or COMPLEXDATA3D" )
    IF ( (dir.NE.rs2pw) .AND. (dir.NE.pw2rs) ) &
       CALL stop_program ( "rs_pw_transfer", "direction must be rs2pw or pw2rs")

    IF (  rs % parallel ) THEN
       IF ( rs % group_dim ( 1 ) == 1 ) THEN
          CALL rs_pw_transfer_replicated(rs,pw,dir)
       ELSE
          CALL rs_pw_transfer_distributed(rs,pw,dir)
       END IF
    ELSE ! treat simple serial case locally
       nn = SIZE ( rs % r )
       IF ( dir == rs2pw ) THEN
          IF ( pw % in_use == REALDATA3D ) THEN
             CALL dcopy ( nn, rs % r, 1, pw % cr3d, 1 )
          ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
             CALL copy_rc ( nn, rs % r, pw % cc3d )
          ELSE
             CALL stop_program ( "rs_pw_transfer", "Pw type not compatible" )
          END IF
       ELSE
          IF ( pw % in_use == REALDATA3D ) THEN
             CALL dcopy ( nn, pw % cr3d, 1, rs % r, 1 )
          ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
             CALL copy_cr ( nn, pw % cc3d, rs % r )
          ELSE
             CALL stop_program ( "rs_pw_transfer", "Pw type not compatible" )
          END IF
       END IF
    END IF
    CALL timestop(0.0_dp,handle)
    CALL timestop(0.0_dp,handle2)
  END SUBROUTINE rs_pw_transfer

  ! treat the two different parallel cases here
  SUBROUTINE rs_pw_transfer_replicated(rs,pw,dir)
    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(pw_type), POINTER                   :: pw
    INTEGER, INTENT(IN)                      :: dir

    INTEGER                                  :: dest, group, ierr, ii, ip, &
                                                ix, iy, iz, mepos, nma, nn, &
                                                np, req(2), s(3), source
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcount
    INTEGER, DIMENSION(3)                    :: lb, ub
    INTEGER, DIMENSION(:, :), POINTER        :: pbo
    INTEGER, DIMENSION(:, :, :), POINTER     :: bo
    REAL(KIND=dp), DIMENSION(:), POINTER     :: recvbuf, sendbuf, swapptr
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: grid

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
    ALLOCATE(sendbuf(nma),recvbuf(nma), STAT=ierr)
    IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer", "sendbuf/recvbuf", 2*nma )
    grid=>rs%r

    IF ( dir == rs2pw ) THEN
       dest  =MODULO(mepos+1,np)
       source=MODULO(mepos-1,np)
       sendbuf=0.0_dp

       DO ip = 1, np

          lb = pbo(1,:)+bo(1,:,MODULO(mepos-ip,np)+1)-1
          ub = pbo(1,:)+bo(2,:,MODULO(mepos-ip,np)+1)-1
          ! this loop takes about the same time as the message passing call
          ! notice that the range of ix is only a small fraction of the first index of grid
          ! therefore it seems faster to have the second index as the innermost loop
          ! if this runs on many cpus
          ! tested on itanium, pentium4, opteron, ultrasparc...
          s=ub-lb+1
          DO iz = lb(3), ub(3)
           DO ix = lb(1), ub(1)
             ii= (iz-lb(3))*s(1)*s(2)+(ix-lb(1))+1
             DO iy = lb(2), ub(2) 
                sendbuf(ii) = sendbuf(ii) + grid(ix,iy,iz)
                ii=ii+s(1)
            END DO
           END DO
          END DO
          IF ( ip .EQ. np ) EXIT
          CALL mp_isendrecv(sendbuf,dest,recvbuf,source, &
                                   group,req(1),req(2),13) 
          CALL mp_waitall(req)
          swapptr=>sendbuf
          sendbuf=>recvbuf
          recvbuf=>swapptr
       ENDDO
       nn = rcount(mepos)
       IF ( pw % in_use == REALDATA3D ) THEN
          CALL dcopy ( nn, sendbuf, 1, pw % cr3d, 1 )
       ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
          CALL copy_rc ( nn, sendbuf, pw % cc3d )
       ELSE
          CALL stop_program ( "rs_pw_transfer", "PW type not compatible" )
       END IF
    ELSE
       nn = rcount ( mepos )
       IF ( pw % in_use == REALDATA3D ) THEN
          CALL dcopy ( nn, pw % cr3d, 1, sendbuf, 1 )
       ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
          CALL dcopy ( nn, pw % cc3d, 2, sendbuf, 1 )
       ELSE
          CALL stop_program ( "rs_pw_transfer", "PW type not compatible" )
       END IF

       dest  =MODULO(mepos+1,np)
       source=MODULO(mepos-1,np)

       DO ip = 0, np-1
          ! we must shift the buffer only np-1 times around
          IF (ip .NE. np-1) THEN
              CALL mp_isendrecv(sendbuf,dest,recvbuf,source, &
                                group,req(1),req(2),13) 
          ENDIF
          lb = pbo(1,:)+bo(1,:,MODULO(mepos-ip,np)+1)-1
          ub = pbo(1,:)+bo(2,:,MODULO(mepos-ip,np)+1)-1
          ii = 0
          ! this loop takes about the same time as the message passing call
          DO iz = lb(3), ub(3)
             DO iy = lb(2), ub(2)
                DO ix = lb(1), ub(1)
                   ii=ii+1
                   grid(ix,iy,iz) = sendbuf(ii)
                END DO
             END DO
          END DO
          IF (ip .NE. np-1) THEN
              CALL mp_waitall(req)
          ENDIF
          swapptr=>sendbuf
          sendbuf=>recvbuf
          recvbuf=>swapptr
       ENDDO

    END IF

    DEALLOCATE ( rcount, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer", "rcount" )
    DEALLOCATE ( sendbuf, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer", "sendbuf" )
    DEALLOCATE ( recvbuf, STAT = ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer", "recvbuf" )
  END SUBROUTINE rs_pw_transfer_replicated

  SUBROUTINE rs_pw_transfer_distributed(rs,pw,dir)
    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(pw_type), POINTER                   :: pw
    INTEGER, INTENT(IN)                      :: dir

    INTEGER                                  :: dest, group, idir, ierr, nn, &
                                                source, subgroup
    INTEGER, DIMENSION(3)                    :: lb, lc, ub, uc
    LOGICAL                                  :: subdim( 2 )
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: buffer

!Both grids are distributed, we have to do some reshuffling of data

    IF ( dir == rs2pw ) THEN

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

    ELSE

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

    END IF
  END SUBROUTINE rs_pw_transfer_distributed

!******************************************************************************
SUBROUTINE send_forward ( rs, pw )

    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(pw_type), INTENT(INOUT)             :: pw

    INTEGER                                  :: idir, ierr, nn, source, &
                                                subgroup
    INTEGER, DIMENSION(3)                    :: lb, lbp, lz, ub, ubp, uz
    LOGICAL, DIMENSION(2)                    :: subdim
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: buffer

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
          CMPLX ( buffer ( lz(1):uz(1), lz(2):uz(2), lz(3):uz(3) ), KIND = dp )
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

    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(pw_type), INTENT(INOUT)             :: pw

    INTEGER                                  :: idir, ierr, nn, source, &
                                                subgroup
    INTEGER, DIMENSION(3)                    :: lb, lbp, lz, ub, ubp, uz
    LOGICAL, DIMENSION(2)                    :: subdim
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: buffer

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
     buffer = 0.0_dp
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
          REAL ( pw % cc3d ( lz(1):uz(1), lz(2):uz(2), lz(3):uz(3) ), KIND=dp )
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

    TYPE(realspace_grid_type), POINTER       :: rs

    INTEGER                                  :: handle, n

!-----------------------------------------------------------------------------!

   CALL timeset("rs_grid_zero",'I',' ',handle)
   n = SIZE ( rs % r )
   CALL dcopy ( n, 0.0_dp, 0, rs % r, 1 )
   CALL timestop(0.0_dp,handle)

END SUBROUTINE rs_grid_zero
!******************************************************************************
SUBROUTINE rs_grid_to_cube ( rs, iunit, ionode, title, particles, stride, error )
    TYPE(realspace_grid_type), POINTER       :: rs
    INTEGER, INTENT(IN)                      :: iunit
    LOGICAL, INTENT(IN)                      :: ionode
    CHARACTER(*), OPTIONAL                   :: title
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: particles
    INTEGER, INTENT(in), OPTIONAL            :: stride
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = "rs_grid_to_cube", &
      routineP = moduleN//":"//routineN

    INTEGER                                  :: I1, I2, I3, iat, L1, L2, L3, &
                                                my_stride, np, U1, U2, U3, z
    LOGICAL                                  :: failure
    TYPE(particle_list_type), POINTER        :: my_particles

!-----------------------------------------------------------------------------!
! not working for distributed grids !

   failure=.FALSE.
   my_stride=1
   IF (PRESENT(stride)) my_stride=stride
   CPPrecondition(my_stride>0,cp_failure_level,routineP,error,failure) 
   IF (rs%ngpts.ne.rs%ngpts_local) CALL stop_program("rs_grid_to_cube","NYI")
   IF (ionode) THEN
     ! this format seems to work for e.g. molekel and gOpenmol
     WRITE(iunit,'(a11)') "-Quickstep-"
     IF (PRESENT(title)) THEN
        WRITE(iunit,*) TRIM(title)
     ELSE
        WRITE(iunit,*) "No Title"
     ENDIF
     
     np=0
     NULLIFY(my_particles)
     IF (PRESENT(particles)) my_particles=>particles
     IF (ASSOCIATED(my_particles)) np=my_particles%n_els

     WRITE(iunit,'(I5,3f12.6)') np,0.0_dp,0._dp,0._dp !start of cube

     WRITE(iunit,'(I5,3f12.6)') (rs%npts(1)+my_stride-1)/my_stride,&
          rs%dr(1)*REAL(my_stride,dp),0.0_dp,0.0_dp
     WRITE(iunit,'(I5,3f12.6)') (rs%npts(2)+my_stride-1)/my_stride,&
          0.0_dp,rs%dr(2)*REAL(my_stride,dp),0.0_dp
     WRITE(iunit,'(I5,3f12.6)') (rs%npts(3)+my_stride-1)/my_stride,&
          0.0_dp,0.0_dp,rs%dr(3)*REAL(my_stride,dp)

     IF (ASSOCIATED(my_particles)) THEN
        DO iat=1,np
           CALL get_atomic_kind(my_particles%els(iat)%atomic_kind,z=z)

           WRITE(iunit,'(I5,4f12.6)') z,0._dp,my_particles%els(iat)%r
        END DO
     END IF

     ! write(iunit,'(I5,4f12.6)') 1, 0.0,0.0,0.0,0.0
     L1=LBOUND(rs%r,1)
     L2=LBOUND(rs%r,2)
     L3=LBOUND(rs%r,3)
     U1=UBOUND(rs%r,1)
     U2=UBOUND(rs%r,2)
     U3=UBOUND(rs%r,3)
     DO I1=L1,U1,my_stride
       DO I2=L2,U2,my_stride
        WRITE(iunit,'(6E13.5)') (rs%r(MODULO(I1-L1,U1-L1+1)+L1,MODULO(I2-L2,U2-L2+1)+L2, &
                                      MODULO(I3-L3,U3-L3+1)+L3),I3=L3,U3,my_stride)
       ENDDO
     ENDDO
   ENDIF
END SUBROUTINE rs_grid_to_cube
!******************************************************************************

!******************************************************************************
! optionally pass a preallocated rs grid
!******************************************************************************
SUBROUTINE rs_pw_to_cube ( pw, iunit, ionode, title, rs, particles, stride,&
     force_env_section, error )
    TYPE(pw_type), POINTER                   :: pw
    INTEGER, INTENT(IN)                      :: iunit
    LOGICAL, INTENT(IN)                      :: ionode
    CHARACTER(*), INTENT(IN), OPTIONAL       :: title
    TYPE(realspace_grid_type), OPTIONAL, &
      POINTER                                :: rs
    TYPE(particle_list_type), OPTIONAL, &
      POINTER                                :: particles
    INTEGER, INTENT(in), OPTIONAL            :: stride
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: handle, nsmax
    TYPE(realspace_grid_type), POINTER       :: rs_local

!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

   CALL timeset("rs_pw_to_cube",handle)
   IF (PRESENT(rs)) THEN
     rs_local=>rs
   ELSE
     nsmax=-1
     NULLIFY(rs_local)
     CALL rs_grid_create(rs_local,pw%pw_grid,nsmax,force_env_section=force_env_section)
   ENDIF

   CALL rs_pw_transfer(rs_local,pw,pw2rs)
   IF (PRESENT(particles)) THEN
      CALL rs_grid_to_cube(rs_local,iunit,ionode,title,particles=particles,&
           stride=stride,error=error)
   ELSE
      CALL rs_grid_to_cube(rs_local,iunit,ionode,title,stride=stride, error=error)
   END IF

   IF (.NOT. PRESENT(rs)) THEN
       CALL rs_grid_release(rs_local)
   ENDIF
   CALL timestop(handle)

END SUBROUTINE rs_pw_to_cube
!******************************************************************************


!!****f* realspace_grid_types/rs_grid_create [1.0] *
!!
!!   NAME
!!     rs_grid_create
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
SUBROUTINE rs_grid_create(rs_grid, pw_grid, nsmax, force_env_section, error)
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    TYPE(pw_grid_type), POINTER              :: pw_grid
    INTEGER, INTENT(in), OPTIONAL            :: nsmax
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_grid_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: my_nsmax, output_unit, stat
    LOGICAL                                  :: failure, should_output
    TYPE(cp_logger_type), POINTER            :: logger

  failure=.FALSE.
  logger => cp_error_get_logger(error)
  my_nsmax = -1
  IF (PRESENT(nsmax)) my_nsmax=nsmax
  
  CPPrecondition(ASSOCIATED(pw_grid),cp_failure_level,routineP,error,failure)
  
  ALLOCATE(rs_grid,stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

  IF (.NOT. failure) THEN
     output_unit = cp_print_key_unit_nr(logger,force_env_section,"PRINT%RS_GRID_INFORMATION",&
          extension=".Log",error=error)     
     should_output=BTEST(cp_print_key_should_output(logger%iter_info,&
          force_env_section,"PRINT%RS_GRID_INFORMATION",error=error),cp_p_file)
     CALL rs_grid_setup(rs_grid, pw_grid, my_nsmax,should_output=should_output, iounit=output_unit)
     CALL cp_print_key_finished_output(output_unit,logger,force_env_section,&
          "PRINT%RS_GRID_INFORMATION",error=error)  
  END IF
END SUBROUTINE rs_grid_create
!***************************************************************************

!!****f* realspace_grid_types/rs_grid_retain [1.0] *
!!
!!   NAME
!!     rs_grid_retain
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
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_grid_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

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
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_grid_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(rs_grid)) THEN
     CPPreconditionNoFail(rs_grid%ref_count>0,cp_failure_level,routineP,error)
     rs_grid%ref_count=rs_grid%ref_count-1
     IF (rs_grid%ref_count==0) THEN
        rs_grid%ref_count=1
        CALL rs_grid_deallocate(rs_grid)
        DEALLOCATE(rs_grid, stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
     END IF
  END IF
  NULLIFY(rs_grid)
END SUBROUTINE rs_grid_release

!!****f* realspace_grid_types/rs_grid_p_create [1.0] *
!!
!!   NAME
!!     rs_grid_p_create
!!
!!   FUNCTION
!!     allocates and initialize an rs_grid_p_type
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
!!     teo
!!
!!   MODIFICATION HISTORY
!!     12.2005 created [teo]
!!
!!*** **********************************************************************
SUBROUTINE rs_grid_p_create(rs_p_grid, ndim, pw_grid, nsmax, force_env_section, error)
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_p_grid
    INTEGER, INTENT(in), OPTIONAL            :: ndim
    TYPE(pw_grid_type), POINTER              :: pw_grid
    INTEGER, INTENT(in), OPTIONAL            :: nsmax
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_grid_p_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, my_nsmax, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  my_nsmax = -1
  IF (PRESENT(nsmax)) my_nsmax=nsmax
  
  CPPrecondition(ASSOCIATED(pw_grid),cp_failure_level,routineP,error,failure)
  
  ALLOCATE(rs_p_grid(ndim),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  DO j = 1, ndim
     NULLIFY(rs_p_grid(j)%rs_grid)
     CALL rs_grid_create(rs_p_grid(j)%rs_grid, pw_grid,my_nsmax, force_env_section)
  END DO

END SUBROUTINE rs_grid_p_create
!***************************************************************************

!!****f* realspace_grid_types/rs_grid_p_release [1.0] *
!!
!!   NAME
!!     rs_grid_p_release
!!
!!   FUNCTION
!!     releases the given rs_grid_p_type
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
!!     teo
!!
!!   MODIFICATION HISTORY
!!     12.2005 created [teo]
!!
!!*** **********************************************************************
SUBROUTINE rs_grid_p_release(rs_grid,error)
    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: rs_grid
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_grid_p_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: j, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(rs_grid)) THEN
     DO j = 1, SIZE(rs_grid)        
        CALL rs_grid_release(rs_grid(j)%rs_grid,error)
     END DO
     DEALLOCATE(rs_grid, stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
  NULLIFY(rs_grid)
END SUBROUTINE rs_grid_p_release


!***************************************************************************
!!****f* realspace_grid_types/rs_get_my_tasks [1.0] *
!!
!!   NAME
!!     rs_get_my_tasks
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

SUBROUTINE rs_get_my_tasks ( rs, tasks, npme, ival, rval, pmat,pmat2,pmat3, pcor ,symmetric)

    TYPE(realspace_grid_type), POINTER       :: rs
    INTEGER, DIMENSION(:, :), POINTER        :: tasks
    INTEGER, INTENT(OUT)                     :: npme
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: ival
    REAL(KIND=dp), DIMENSION(:, :), &
      OPTIONAL, POINTER                      :: rval
    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: pmat, pmat2, pmat3
    INTEGER, DIMENSION(:, :), OPTIONAL, &
      POINTER                                :: pcor
    LOGICAL, INTENT(IN), OPTIONAL            :: symmetric

    INTEGER :: acol, arow, handle, i, ic, icmax, ileft, isend, j, k, lb, n1, &
      npmax, nppp, plength, plength2, plength3, stat, subgroup, subpos, &
      subsize, ub
    INTEGER, DIMENSION(:, :), POINTER        :: ilist, plist, pls, pls2, &
                                                pls3, tlist
    LOGICAL                                  :: matrix, my_symmetric
    LOGICAL, DIMENSION(2)                    :: subdim
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ppack, ppack2, ppack3
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block, rlist

  CALL timeset("rs_get_my_tasks",'I',' ',handle)

  IF ( .NOT. ASSOCIATED ( tasks ) ) &
        CALL stop_program ( "get_my_tasks", "tasks not associated" )

  IF ( PRESENT ( pmat ) .AND. PRESENT ( pcor ) ) THEN
    matrix = .TRUE.
  ELSE
    matrix = .FALSE.
  END IF
  IF(PRESENT(symmetric) )THEN
    my_symmetric = symmetric
  ELSE
    my_symmetric = .TRUE.
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
      rlist(1:n1,npme+1:npmax)=0.0_dp
      rval=0.0_dp
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
      IF(PRESENT(pmat2))THEN
        ALLOCATE(pls2(n1,npmax),STAT=stat)
        IF (stat/=0) CALL stop_memory("get_my_tasks","pls",n1*npmax)
        pls2= 0
        ALLOCATE(ppack2(npmax),STAT=stat)
        IF (stat/=0) CALL stop_memory("get_my_tasks","ppack",npmax)
      END IF
      IF(PRESENT(pmat3))THEN
        ALLOCATE(pls3(n1,npmax),STAT=stat)
        IF (stat/=0) CALL stop_memory("get_my_tasks","pls",n1*npmax)
        pls3= 0
        ALLOCATE(ppack3(npmax),STAT=stat)
        IF (stat/=0) CALL stop_memory("get_my_tasks","ppack",npmax)
      END IF
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
          IF(my_symmetric) THEN
            IF ( pcor(1,i) <= pcor(2,i) ) THEN
               arow = pcor(1,i)
               acol = pcor(2,i)
            ELSE
               arow = pcor(2,i)
               acol = pcor(1,i)
            END IF
          ELSE
            arow = pcor(1,i)
            acol = pcor(2,i)
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
          CALL pack_matrix ( pmat, pls, ppack, plength, subgroup ,symmetric=my_symmetric)
          CALL mp_shift ( ppack ( 1:plength ), subgroup, isend )
          CALL mp_shift ( pls ( 1:2, 1:icmax ), subgroup, isend )
          CALL unpack_matrix ( pmat, pcor, pls, ic, ppack , symmetric=my_symmetric)
          IF(PRESENT(pmat2)) THEN
            pls2(:,1:ic) = pcor(:,1:ic)
            pls2(:,ic+1:) = 0
            CALL mp_shift ( pls2 ( 1:2, 1:icmax ), subgroup, -isend )
          ! pack the data to be sent
            CALL pack_matrix ( pmat2, pls2, ppack2, plength2, subgroup,symmetric=my_symmetric )
            CALL mp_shift ( ppack2 ( 1:plength2 ), subgroup, isend )
            CALL mp_shift ( pls2 ( 1:2, 1:icmax ), subgroup, isend )
            CALL unpack_matrix ( pmat2, pcor, pls2, ic, ppack2,symmetric=my_symmetric )
          END IF
          IF(PRESENT(pmat3)) THEN
            pls3(:,1:ic) = pcor(:,1:ic)
            pls3(:,ic+1:) = 0
            CALL mp_shift ( pls3 ( 1:2, 1:icmax ), subgroup, -isend )
          ! pack the data to be sent
            CALL pack_matrix ( pmat3, pls3, ppack3, plength3, subgroup,symmetric=my_symmetric )
            CALL mp_shift ( ppack3 ( 1:plength3 ), subgroup, isend )
            CALL mp_shift ( pls3 ( 1:2, 1:icmax ), subgroup, isend )
            CALL unpack_matrix ( pmat3, pcor, pls3, ic, ppack3,symmetric=my_symmetric )
          END IF
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
      IF(PRESENT(pmat2)) THEN
        DEALLOCATE(pls2,STAT=stat)
        IF (stat/=0) CALL stop_memory("get_my_tasks","pls")
        DEALLOCATE(ppack2,STAT=stat)
        IF (stat/=0) CALL stop_memory("get_my_tasks","ppack")
      END IF
      IF(PRESENT(pmat3)) THEN
        DEALLOCATE(pls3,STAT=stat)
        IF (stat/=0) CALL stop_memory("get_my_tasks","pls")
        DEALLOCATE(ppack3,STAT=stat)
        IF (stat/=0) CALL stop_memory("get_my_tasks","ppack")
      END IF
    END IF

    ! Release the communicator
    CALL mp_comm_free ( subgroup )

  ELSE

    ! fully replicated grids, each processor can process all its tasks

  END IF

  CALL timestop(0.0_dp,handle)

END SUBROUTINE rs_get_my_tasks

!***************************************************************************
SUBROUTINE pack_matrix ( pmat, pls, ppack , plength, group , symmetric)

    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: pmat
    INTEGER, DIMENSION(:, :), POINTER        :: pls
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ppack
    INTEGER, INTENT(OUT)                     :: plength
    INTEGER, INTENT(IN)                      :: group
    LOGICAL, INTENT(IN), OPTIONAL            :: symmetric

    INTEGER                                  :: acol, arow, i, j, k, n, nc, nr
    LOGICAL                                  :: my_symmetric
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block

  my_symmetric = .TRUE.
  IF(PRESENT(symmetric)) my_symmetric = symmetric
  plength = 0
  DO i = 1, SIZE ( pls, 2 )
    IF ( pls ( 1, i ) == 0 ) EXIT
    IF(my_symmetric) THEN
      IF ( pls(1,i) <= pls(2,i) ) THEN
         arow = pls(1,i)
         acol = pls(2,i)
      ELSE
         arow = pls(2,i)
         acol = pls(1,i)
      END IF
    ELSE
       arow = pls(1,i)
       acol = pls(2,i)
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

SUBROUTINE unpack_matrix ( pmat, pcor, pls, ic, ppack , symmetric)

    TYPE(real_matrix_type), OPTIONAL, &
      POINTER                                :: pmat
    INTEGER, DIMENSION(:, :), POINTER        :: pcor, pls
    INTEGER, INTENT(IN)                      :: ic
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ppack
    LOGICAL, INTENT(IN), OPTIONAL            :: symmetric

    INTEGER                                  :: acol, arow, i, j, k, nc, nr, &
                                                pl
    LOGICAL                                  :: my_symmetric
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: p_block

  my_symmetric = .TRUE.
  IF(PRESENT(symmetric)) my_symmetric = symmetric
  pl = 0
  DO i = 1, ic
    IF(my_symmetric) THEN
      IF ( pcor(1,i) <= pcor(2,i) ) THEN
         arow = pcor(1,i)
         acol = pcor(2,i)
      ELSE
         arow = pcor(2,i)
         acol = pcor(1,i)
      END IF
    ELSE
       arow = pcor(1,i)
       acol = pcor(2,i)
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

    INTEGER, INTENT(IN)                      :: npme
    INTEGER, DIMENSION(6, npme), INTENT(IN)  :: ival
    INTEGER, INTENT(OUT)                     :: natom_pairs
    INTEGER, DIMENSION(:, :), POINTER        :: asets, atasks

    INTEGER                                  :: iatom, iatom_old, iset, &
                                                iset_old, itask, jatom, &
                                                jatom_old, jset, jset_old, &
                                                nset_pairs

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

END MODULE realspace_grid_types

!******************************************************************************
