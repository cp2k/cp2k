!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
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
!*****
!******************************************************************************

MODULE realspace_grid_types
  USE input_constants,                 ONLY: rsgrid_replicated,rsgrid_automatic,rsgrid_distributed
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_get
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: dp
  USE mathlib,                         ONLY: det_3x3
  USE message_passing,                 ONLY: &
       mp_alltoall, mp_cart_coords, mp_cart_create, mp_cart_shift, &
       mp_comm_dup, mp_comm_free, mp_environ, mp_isendrecv, mp_max, mp_min, &
       mp_sendrecv, mp_sum, mp_sync, mp_waitall
  USE pw_grid_types,                   ONLY: PW_MODE_LOCAL,&
                                             pw_grid_type
  USE pw_types,                        ONLY: COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             pw_integrate_function,&
                                             pw_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE util,                            ONLY: get_limit
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: realspace_grid_type,&
            realspace_grid_p_type,&
            realspace_grid_input_type

  PUBLIC :: rs_pw_transfer,&
            rs_grid_zero,&
            rs_grid_set_box,&
            rs_grid_create,&
            rs_grid_retain,&
            rs_grid_release,&
            rs_grid_print, &
            rs_find_node,&
            init_input_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'realspace_grid_types'
  INTEGER, SAVE, PRIVATE :: last_rs_id=0
  INTEGER, SAVE, PRIVATE :: allocated_rs_grid_count=0
  INTEGER, PARAMETER, PUBLIC :: rs2pw=11,pw2rs=12

  TYPE realspace_grid_input_type
     INTEGER       :: distribution_type
     INTEGER       :: distribution_layout(3)
     REAL(KIND=dp) :: memory_factor
     INTEGER       :: nsmax
  END TYPE realspace_grid_input_type

  TYPE realspace_grid_type
     INTEGER :: grid_id                                  ! tag of the pw_grid
     INTEGER :: id_nr                                    ! unique identifier of rs
     INTEGER :: ref_count                                ! reference count

     REAL(KIND=dp), DIMENSION ( :, :, : ),POINTER :: r   ! the grid

     INTEGER :: ngpts                                    ! # grid points
     INTEGER, DIMENSION (3) :: npts                      ! # grid points per dimension
     INTEGER, DIMENSION (3) :: lb                        ! lower bounds
     INTEGER, DIMENSION (3) :: ub                        ! upper bounds

     INTEGER, DIMENSION (:), POINTER :: px,py,pz         ! index translators

     INTEGER :: border                                   ! border points
     INTEGER :: ngpts_local                              ! local dimensions
     INTEGER, DIMENSION (3) :: npts_local
     INTEGER, DIMENSION (3) :: lb_local
     INTEGER, DIMENSION (3) :: ub_local
     INTEGER, DIMENSION (3) :: lb_real                   ! lower bounds of the real local data
     INTEGER, DIMENSION (3) :: ub_real                   ! upper bounds of the real local data

     INTEGER, DIMENSION (3) :: perd                      ! periodicity enforced
     REAL(KIND=dp), DIMENSION(3,3) :: dh                 ! incremental grid matrix
     REAL(KIND=dp), DIMENSION(3,3) :: dh_inv             ! inverse incremental grid matrix
     LOGICAL :: orthorhombic                             ! grid symmetry

     LOGICAL :: parallel                                 ! whether the corresponding pw grid is distributed
     LOGICAL :: distributed                              ! whether the rs grid is distributed
     INTEGER :: group
     INTEGER :: my_pos
     LOGICAL :: group_head
     INTEGER :: group_size
     INTEGER, DIMENSION (3) :: group_dim
     INTEGER, DIMENSION (3) :: group_coor
     INTEGER, DIMENSION (3) :: neighbours

  END TYPE realspace_grid_type

  TYPE realspace_grid_p_type
     TYPE(realspace_grid_type), POINTER :: rs_grid
  END TYPE realspace_grid_p_type

!-----------------------------------------------------------------------------!

CONTAINS

!!****f* realspace_grid_types/init_input_type *
!!
!!   NAME
!!     init_input_type
!!
!!   FUNCTION
!!     parses an input section to assign the proper values to the input type
!!
!!   NOTES
!!     if rs_grid_section is not present we setup for an replicated setup
!!
!!   INPUTS
!!    -
!!    -
!!
!!   MODIFICATION HISTORY
!!     01.2008 created [Joost VandeVondele]
!!
!!*** **********************************************************************
SUBROUTINE init_input_type(input_settings,nsmax,rs_grid_section,ilevel,error)
   TYPE(realspace_grid_input_type), INTENT(OUT) :: input_settings
   INTEGER, INTENT(IN)                          :: nsmax
   TYPE(section_vals_type), POINTER, OPTIONAL   :: rs_grid_section
   INTEGER, INTENT(IN)                          :: ilevel
   TYPE(cp_error_type), INTENT(inout)           :: error

   INTEGER,POINTER, DIMENSION(:)                :: tmp
   INTEGER                                      :: isection, nsection, max_distributed_level

   IF (PRESENT(rs_grid_section)) THEN
       input_settings%nsmax=nsmax
       ! we use the section corresponding to the level, or the largest available one
       ! i.e. the last section defines all following ones
       CALL section_vals_get(rs_grid_section,n_repetition=nsection,error=error)
       isection=MAX(1,MIN(ilevel,nsection))
       CALL section_vals_val_get(rs_grid_section,"DISTRIBUTION_TYPE",&
                i_rep_section=isection,&
                i_val=input_settings%distribution_type,error=error)
       CALL section_vals_val_get(rs_grid_section,"DISTRIBUTION_LAYOUT",&
                i_rep_section=isection,&
                i_vals=tmp,error=error)
       input_settings%distribution_layout=tmp
       CALL section_vals_val_get(rs_grid_section,"MEMORY_FACTOR",&
                i_rep_section=isection,&
                r_val=input_settings%memory_factor,error=error)
       CALL section_vals_val_get(rs_grid_section,"MAX_DISTRIBUTED_LEVEL",&
                i_rep_section=isection,&
                i_val=max_distributed_level,error=error)

       ! multigrids that are to coarse are not distributed in the automatic scheme
       IF (input_settings%distribution_type == rsgrid_automatic) THEN
          IF (ilevel>max_distributed_level) THEN
              input_settings%distribution_type=rsgrid_replicated
          ENDIF
       ENDIF
   ELSE
       input_settings%nsmax=-1
       input_settings%distribution_type=rsgrid_replicated
   ENDIF
END SUBROUTINE init_input_type
!******************************************************************************
!!****** realspace_grid_types/rs_grid_create [1.0] *
!!
!!   NAME
!!     rs_grid_create
!!
!!   FUNCTION
!!     Determine the setup of real space grids
!!
!!   AUTHOR
!!     JGH (18-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     JGH (08-Jun-2003) : nsmax <= 0 indicates fully replicated grid
!*****
!******************************************************************************

  SUBROUTINE rs_grid_create ( rs, pw_grid, input_settings, error)
    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(realspace_grid_input_type), INTENT(IN) :: input_settings
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rs_grid_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dir, i, j, k, l, lb( 2 ), &
                                                n_slices(3), neighbours(3), &
                                                nmin, pos( 3 ), n_slices_tmp(3), nsmax, stat, handle
    LOGICAL                                  :: failure, overlap
    REAL(KIND=dp)                            :: ratio, ratio_best, volume, &
                                                volume_dist

    CALL timeset(routineN,handle)

    failure = .FALSE.

    CPPrecondition(ASSOCIATED(pw_grid),cp_failure_level,routineP,error,failure)

    ALLOCATE(rs,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL mp_sync(pw_grid % para % group)

    rs % dh = pw_grid%dh
    rs % dh_inv = pw_grid%dh_inv
    rs % orthorhombic = pw_grid%orthorhombic
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
       rs % border = 0
       rs % npts_local = pw_grid % npts
       rs % ngpts_local = PRODUCT ( rs % npts )
       rs % lb_local = pw_grid % bounds ( 1, : )
       rs % ub_local = pw_grid % bounds ( 2, : )
       rs % parallel = .FALSE.
       rs % distributed = .FALSE.
       rs % group = -1
       rs % group_size = 1
       rs % group_head = .TRUE.
       rs % group_dim = 1
       rs % group_coor = 0
       rs % my_pos = 0
    ELSE
       IF ( .NOT. ALL ( pw_grid % bounds ( 1, 2:3 ) == &
            pw_grid % bounds_local ( 1, 2:3 ) .AND. &
            pw_grid % bounds ( 2, 2:3 ) == &
            pw_grid % bounds_local ( 2, 2:3 ) ) ) THEN
          CALL stop_program ( routineN, &
               "This pw type not supported" )
       END IF
 
       ! group size of rs grid
       ! global grid dimensions are still the same
       rs % group_size = pw_grid % para % group_size
       rs % npts = pw_grid % npts
       rs % ngpts = PRODUCT ( rs % npts )
       rs % lb = pw_grid % bounds ( 1, : )
       rs % ub = pw_grid % bounds ( 2, : )

       ! this is the eventual border size
       nmin = ( input_settings % nsmax + 1 ) / 2

       ! x is the main distribution direction of pw grids in real space

       IF ( input_settings % distribution_type == rsgrid_replicated ) THEN

          n_slices = 1

       ELSE
          n_slices = 1
          ratio_best=-HUGE(ratio_best)

          ! don't allow distributions with more processors than real grid points
          DO k=1, MIN(rs % npts (3), rs % group_size)
          DO j=1, MIN(rs % npts (2), rs % group_size)
             i= MIN (rs % npts (1), rs % group_size/(j*k) )
             n_slices_tmp=(/i,j,k/)
     
             ! we don't match the actual number of CPUs
             IF (PRODUCT(n_slices_tmp) .NE. rs % group_size) CYCLE

             ! we see if there has been a input constraint
             ! i.e. if the layout is not -1 we need to fullfil it
             IF (.NOT. ALL(PACK(n_slices_tmp==input_settings % distribution_layout,&
                                (/-1,-1,-1/)/=input_settings % distribution_layout )&
                          )) CYCLE

             ! we currently can not work with a grid that has borders that overlaps with the local data of the grid itself
             overlap=.FALSE.
             DO dir=1,3
                IF (n_slices_tmp(dir)>1) THEN
                   neighbours(dir) = HUGE(0)
                   DO l=0,n_slices_tmp(dir)-1
                      lb = get_limit ( rs % npts ( dir ),  n_slices_tmp( dir ), l )
                      neighbours(dir) = MIN(lb (2) -lb (1) + 1, neighbours(dir) )
                   ENDDO
                   rs % neighbours(dir) =  (nmin + neighbours(dir) - 1)/neighbours(dir)
                   IF( rs % neighbours(dir) .GE. n_slices_tmp (dir) ) overlap=.TRUE.
                ELSE
                   neighbours(dir) = 0
                ENDIF
             ENDDO
             ! XXXXXXX I think the overlap criterium / neighbour calculation is too conservative
             ! write(6,*) n_slices_tmp,rs % neighbours, overlap
             IF (overlap) CYCLE

             ! a heuristic optimisation to reduce the memory usage
             ! we go for the smallest local to real volume
             ! volume of the box without the wings / volume of the box with the wings
             ratio = PRODUCT(REAL(rs % npts,KIND=dp) / n_slices_tmp ) /  & 
                     PRODUCT(REAL(rs % npts,KIND=dp) / n_slices_tmp  +   &
                       MERGE((/0,0,0/),2*(/nmin,nmin,nmin/), n_slices_tmp ==(/1,1,1/)))
             IF (ratio>ratio_best) THEN
                ratio_best=ratio
                n_slices=n_slices_tmp
             ENDIF

          ENDDO
          ENDDO


          ! if automatic we can still decide this is a replicated grid
          ! if the memory gain (or the gain is messages) is too small.
          IF (input_settings % distribution_type == rsgrid_automatic) THEN
             volume=PRODUCT(REAL(rs % npts,KIND=dp)) 
             volume_dist=PRODUCT(REAL(rs % npts,KIND=dp) / n_slices  +   &
                          MERGE((/0,0,0/),2*(/nmin,nmin,nmin/),n_slices==(/1,1,1/)))
             IF (volume<volume_dist*input_settings%memory_factor) THEN
                 n_slices=1
             ENDIF
          ENDIF
          
       END IF

       rs % group_dim (:) = n_slices(:)

       IF ( ALL (n_slices == 1 ) ) THEN
          ! CASE 1 : only one slice: we do not need overlapping regions and special
          !          recombination of the total density
          rs % perd = 1
          rs % border = 0
          rs % npts_local = pw_grid % npts
          rs % ngpts_local = PRODUCT ( rs % npts )
          rs % lb_local = rs % lb
          rs % ub_local = rs % ub
          rs % lb_real = rs % lb
          rs % ub_real = rs % ub
          rs % parallel = .TRUE.
          rs % distributed = .FALSE.
          CALL mp_comm_dup( pw_grid % para % rs_group , rs % group )
          CALL mp_environ(rs % group_size, rs % my_pos, rs % group )
          rs % group_head = pw_grid % para % group_head
          rs % group_coor ( : ) = 0
       ELSE
          ! CASE 2 : general case
          ! periodicity is no longer enforced arbritary  directions
          rs % perd = 1
          DO dir = 1,3
             IF (n_slices(dir).GT.1) THEN
                rs % perd ( dir ) = 0
             ENDIF
          ENDDO
          ! we keep a border of nmin points
          rs % border = nmin
          ! we are going parallel on the real space grid
          rs % parallel = .TRUE.
          rs % distributed = .TRUE.
          ! the new cartesian group
          CALL mp_cart_create ( pw_grid % para % group, 3, &
               rs % group_dim, pos, rs % group )
          CALL mp_environ(rs % group_size, rs % my_pos, rs % group )
          ! for output would be needed to have this pw_grid % para % group_head
          rs % group_head = ALL ( pos == 0 ) 
          rs % group_coor = pos
          ! local dimensions of the grid
          rs % npts_local = pw_grid % npts
          rs % lb_local = rs % lb
          rs % ub_local = rs % ub
          rs % lb_real = rs % lb
          rs % ub_real = rs % ub
          DO dir = 1,3
             rs % neighbours(dir) = 0
             IF (n_slices(dir).GT.1) THEN             
                neighbours(dir) = HUGE(0)
                DO l=0,n_slices(dir)-1
                   lb = get_limit ( rs % npts ( dir ),  n_slices( dir ), l )
                   neighbours(dir) = MIN(lb (2) -lb (1) + 1, neighbours(dir) )
                ENDDO
                rs % neighbours(dir) =  (nmin + neighbours(dir) - 1)/neighbours(dir)

                lb = get_limit ( rs % npts ( dir ), rs % group_dim ( dir ), pos ( dir ) )
                rs % npts_local ( dir ) = 2 * rs % border + ( lb ( 2 ) - lb ( 1 ) + 1 )
                rs % lb_local ( dir ) = lb ( 1 ) + rs % lb ( dir ) - rs % border - 1
                rs % ub_local ( dir ) = lb ( 2 ) + rs % lb ( dir ) + rs % border - 1
                rs % lb_real ( dir ) = rs % lb_local(dir) + rs % border 
                rs % ub_real ( dir ) = rs % ub_local(dir) - rs % border                
             ENDIF
          ENDDO
          rs % ngpts_local = PRODUCT ( rs % npts_local )
       END IF

    END IF


    allocated_rs_grid_count = allocated_rs_grid_count + 1

    ALLOCATE ( rs % r (rs % lb_local(1):rs % ub_local(1), &
                       rs % lb_local(2):rs % ub_local(2), &
                       rs % lb_local(3):rs % ub_local(3)), STAT = stat )
    IF ( stat /= 0 ) CALL stop_memory ( routineN,"rs % r", rs % ngpts_local )
    ALLOCATE ( rs % px ( rs % npts ( 1 ) ), STAT = stat )
    IF ( stat /= 0 ) CALL stop_memory ( routineN,"rs % px", rs % npts ( 1 ) )
    ALLOCATE ( rs % py ( rs % npts ( 2 ) ), STAT = stat )
    IF ( stat /= 0 ) CALL stop_memory ( routineN, "rs % py", rs % npts ( 2 ) )
    ALLOCATE ( rs % pz ( rs % npts ( 3 ) ), STAT = stat )
    IF ( stat /= 0 ) CALL stop_memory ( routineN,"rs % pz", rs % npts ( 3 ) )

    CALL timestop(handle)

  END SUBROUTINE rs_grid_create

!******************************************************************************
!!****** realspace_grid_types/rs_grid_print [1.0] *
!!
!!   NAME
!!     rs_grid_print
!!
!!   FUNCTION
!!     Print information on grids to output
!!
!!   AUTHOR
!!     JGH (17-May-2007)
!!
!!   MODIFICATION HISTORY
!*****
!******************************************************************************

  SUBROUTINE rs_grid_print ( rs, iounit, error)
    TYPE(realspace_grid_type), POINTER       :: rs
    INTEGER, INTENT(in)                      :: iounit
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rs_grid_print', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dir, i, nn
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: pp( 3 )

    failure = .FALSE.

    IF ( rs % parallel ) THEN
       IF ( iounit>0) THEN
          WRITE ( iounit, '(/,A,T71,I10)' ) &
               " RS_GRID: Information for grid number ", rs % grid_id
          DO i = 1, 3
             WRITE ( iounit, '(A,I3,T30,2I8,T62,A,T71,I10)' ) " RS_GRID:   Bounds ", &
                  i, rs % lb ( i ), rs % ub ( i ), "Points:", rs % npts ( i )
          END DO
          IF ( .NOT. rs % distributed ) THEN
             WRITE ( iounit, '(A)' ) " RS_GRID: Real space fully replicated"
             WRITE ( iounit, '(A,T71,I10)' ) &
                  " RS_GRID: Group size ", rs % group_dim (2)
          ELSE
             DO dir = 1,3
                IF ( rs % perd (dir) /= 1 ) THEN 
                   WRITE ( iounit, '(A,T71,I3,A)' ) &
                        " RS_GRID: Real space distribution over ", rs % group_dim ( dir ), " groups"
                   WRITE ( iounit, '(A,T71,I10)' ) &
                        " RS_GRID: Real space distribution along direction ", dir
                   WRITE ( iounit, '(A,T71,I10)' ) &
                        " RS_GRID: Border size ", rs % border
                END IF
             END DO
          END IF
       END IF
       IF ( rs %distributed ) THEN
          DO dir = 1 ,3
             IF ( rs % perd (dir) /= 1 ) THEN
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
                END IF
             END IF
          END DO
!          WRITE ( iounit, '(/)' )
       END IF
    ELSE
       IF (iounit>0) THEN
         WRITE ( iounit, '(/,A,T71,I10)' ) &
              " RS_GRID: Information for grid number ", rs % grid_id
         DO i = 1, 3
            WRITE ( iounit, '(A,I3,T30,2I8,T62,A,T71,I10)' ) " RS_GRID:   Bounds ", &
                 i, rs % lb ( i ), rs % ub ( i ), "Points:", rs % npts ( i )
         END DO
!         WRITE ( iounit, '(/)' )
       ENDIF
    END IF

  END SUBROUTINE rs_grid_print

!******************************************************************************
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
!*****
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
       IF ( .NOT. rs % distributed ) THEN
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

!!****f* realspace_grid_types/rs_pw_transfer_replicated *
!!
!!   NAME
!!     rs_pw_transfer_replicated
!!
!!   FUNCTION
!!     transfer from a realspace grid to a planewave grid or the other way around
!!
!!   NOTES
!!     rs2pw sums all data on the rs grid into the respective pw grid
!!     pw2rs will scatter all pw data to the rs grids
!!
!!   INPUTS
!!     - dir == rs2pw or dir == pw2rs
!!
!!   MODIFICATION HISTORY
!!
!!*** **********************************************************************
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

!!****f* realspace_grid_types/rs_pw_transfer_distributed *
!!
!!   NAME
!!      rs_pw_transfer_distributed
!!
!!   FUNCTION
!!      does the rs2pw and pw2rs transfer in the case where the rs grid is 
!!      distributed (3D domain decomposition)
!!
!!   NOTES
!!      the transfer is a two step procedure. For example, for the rs2pw transfer:
!!     
!!      1) Halo-exchange in 3D so that the local part of the rs_grid contains the full data
!!      2) an alltoall communication to redistribute the local rs_grid to the local pw_grid
!!
!!      the halo exchange is most expensive on a large number of CPUs. Particular in this halo
!!      exchange is that the border region is rather large (e.g. 20 points) and that it might overlap
!!      with the central domain of several CPUs (i.e. next nearest neighbors)
!!
!!   INPUTS
!!    -
!!    -
!!
!!   MODIFICATION HISTORY
!!     12.2007 created [Matt Watkins]
!!
!!*** **********************************************************************
  SUBROUTINE rs_pw_transfer_distributed(rs,pw,dir)
    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(pw_type), POINTER                   :: pw
    INTEGER, INTENT(IN)                      :: dir

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_pw_transfer_distributed', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=200)                       :: error_string
    INTEGER                                  :: dest, group, i, idir, ierr, &
                                                j, k, my_pw_rank, my_rs_rank, &
                                                n_shifts, nn, position, &
                                                source, x, y, z
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: dshifts, rank_pw2rs, &
                                                recv_disps, recv_sizes, &
                                                send_disps, send_sizes, &
                                                ushifts
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: bounds, recv_tasks, send_tasks
    INTEGER, DIMENSION(2)                    :: neighbours, pos
    INTEGER, DIMENSION(3) :: coords, lb_recv, lb_recv_new, lb_send, &
      lb_send_new, ub_recv, ub_recv_new, ub_send, ub_send_new
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: pw_sum, rs_sum
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_buf, send_buf
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: recv_buf_3d, send_buf_3d
    TYPE(cp_error_type)                      :: error

! Both grids are distributed, we have to do some reshuffling of data
! first map pw ranks (corresponding to the cartesian communicator) to rs ranks (and their cart. comm)
! this could become part of the rs type

    ALLOCATE(rank_pw2rs(0:rs%group_size-1), STAT=ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","buffer",nn)
    rank_pw2rs=0
    rank_pw2rs(pw%pw_grid%para%rs_mpo)=rs%my_pos
    CALL mp_sum(rank_pw2rs,rs%group)

    IF ( dir == rs2pw ) THEN

       ! safety check, to be removed once we're absolute sure the routine is correct
       rs_sum=accurate_sum(rs%r)*det_3x3(rs%dh)
       CALL mp_sum(rs_sum,rs%group)

       DO idir = 1,3

          IF ( rs % perd (idir) .NE. 1) THEN

             ALLOCATE ( dshifts ( 0:rs % neighbours (idir ) ), STAT=ierr )
             ALLOCATE ( ushifts ( 0:rs % neighbours (idir ) ), STAT=ierr )
             
             ushifts = 0
             dshifts = 0

             ! check that we don't try to send data to ourself
             DO n_shifts = 1, MIN(rs % neighbours(idir),rs%group_dim(idir)-1)

                ! need to take into account the possible varying widths of neighbouring cells
                ! offset_up and offset_down hold the real size of the neighbouring cells
                position = MODULO ( rs % group_coor (idir) - n_shifts, rs % group_dim (idir))
                neighbours = get_limit ( rs % npts ( idir ), rs % group_dim ( idir ), position )
                dshifts(n_shifts) = dshifts(n_shifts-1) + ( neighbours (2) - neighbours (1) + 1 ) 

                position = MODULO ( rs % group_coor (idir) + n_shifts, rs % group_dim (idir))
                neighbours = get_limit ( rs % npts ( idir ), rs % group_dim ( idir ), position )
                ushifts(n_shifts) = ushifts(n_shifts-1) + ( neighbours (2) - neighbours (1) + 1 )

                ! The border data has to be send/received from the neighbours
                ! First we calculate the source and destination processes for the shift
                ! The first shift is "downwards"

                CALL mp_cart_shift ( rs % group, idir -1, -1 * n_shifts, source, dest )

                lb_send ( : ) = rs % lb_local ( : )
                ub_send ( : ) = rs % ub_local ( : )
                ub_send ( idir ) = lb_send ( idir ) + rs % border - 1 - dshifts(n_shifts-1)

                lb_recv ( : ) = rs % lb_local ( : )
                ub_recv ( : ) = rs % ub_local ( : )
                ub_recv ( idir ) = ub_recv ( idir ) - rs % border
                lb_recv ( idir ) = ub_recv ( idir ) - rs % border + 1 + ushifts(n_shifts-1)

                lb_send_new ( : ) = rs % lb_local ( : )
                lb_recv_new ( : ) = rs % lb_local ( : )
                ub_recv_new ( : ) = rs % ub_local ( : )
                ub_send_new ( : ) = rs % ub_local ( : )

                IF(dshifts (n_shifts - 1) .LE. rs % border ) THEN
                   ub_send_new ( idir ) = lb_send_new ( idir ) + rs % border  - 1 - dshifts(n_shifts-1)
                   lb_send_new ( idir ) = MAX( lb_send_new ( idir ),&
                        lb_send_new ( idir ) + rs % border  - dshifts ( n_shifts ) )
                   
                   ub_recv_new ( idir ) = ub_recv_new ( idir ) - rs % border
                   lb_recv_new ( idir ) = MAX(lb_recv_new(idir) + rs % border,&
                        ub_recv_new ( idir ) - rs % border + 1 + ushifts(n_shifts-1) )
                ELSE
                   lb_send_new( idir ) = 0
                   ub_send_new( idir ) = lb_send_new( idir ) - 1
                   lb_recv_new( idir ) = 0
                   ub_recv_new( idir ) = lb_recv_new( idir ) - 1
                ENDIF

                lb_send = lb_send_new
                ub_send = ub_send_new
                lb_recv = lb_recv_new
                ub_recv = ub_recv_new

                ! Allocate a scratch arrays to receive the data
                nn = PRODUCT ( ub_send - lb_send + 1 )
                ALLOCATE ( send_buf_3d ( lb_send(1):ub_send(1), lb_send(2):ub_send(2), lb_send(3):ub_send(3) ), STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_buf_3d",nn)
                nn = PRODUCT ( ub_recv - lb_recv + 1 )
                ALLOCATE ( recv_buf_3d ( lb_recv(1):ub_recv(1), lb_recv(2):ub_recv(2), lb_recv(3):ub_recv(3) ), STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_buf_3d",nn)
                
                send_buf_3d = rs % r ( lb_send(1):ub_send(1), lb_send(2):ub_send(2), lb_send(3):ub_send(3) )

                CALL mp_sendrecv ( send_buf_3d, dest, recv_buf_3d, source, rs % group )

                ! only some procs may need later shifts
                IF(ub_recv(idir).GE.lb_recv(idir))THEN

                   ! Sum the data in the RS Grid
                   rs % r ( lb_recv(1):ub_recv(1), lb_recv(2):ub_recv(2), lb_recv(3):ub_recv(3) ) = &
                        rs % r ( lb_recv(1):ub_recv(1), lb_recv(2):ub_recv(2), lb_recv(3):ub_recv(3) ) + recv_buf_3d ( :, :, : )
                END IF

                DEALLOCATE ( send_buf_3d, STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_buf_3d")
                DEALLOCATE ( recv_buf_3d, STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_buf_3d")

                ! Now for the other direction

                CALL mp_cart_shift ( rs % group, idir -1 , n_shifts, source, dest )

                lb_send ( : ) = rs % lb_local ( : )
                ub_send ( : ) = rs % ub_local ( : )
                lb_send ( idir ) = ub_send ( idir ) - rs % border + 1 + ushifts(n_shifts-1)

                lb_recv ( : ) = rs % lb_local ( : )
                ub_recv ( : ) = rs % ub_local ( : )
                lb_recv ( idir ) = lb_recv ( idir ) + rs % border
                ub_recv ( idir ) = lb_recv ( idir ) + rs % border -1 - dshifts(n_shifts-1)

                lb_send_new ( : ) = rs % lb_local ( : )
                lb_recv_new ( : ) = rs % lb_local ( : )
                ub_recv_new ( : ) = rs % ub_local ( : )
                ub_send_new ( : ) = rs % ub_local ( : )

                IF (ushifts(n_shifts - 1) .LE. rs % border ) THEN

                   lb_send_new ( idir ) = ub_send_new ( idir ) - rs % border + 1 + ushifts(n_shifts-1)
                   ub_send_new ( idir ) = MIN( ub_send_new( idir ),&
                        ub_send_new ( idir ) - rs % border  + ushifts( n_shifts ) )

                   lb_recv_new ( idir ) = lb_recv_new ( idir ) + rs % border
                   ub_recv_new ( idir ) = MIN(ub_recv_new (idir)- rs%border,&
                        lb_recv_new ( idir ) + rs % border - 1 - dshifts(n_shifts-1))
                ELSE
                   lb_send_new( idir ) = 0
                   ub_send_new( idir ) = lb_send_new( idir ) - 1
                   lb_recv_new( idir ) = 0
                   ub_recv_new( idir ) = lb_recv_new( idir ) - 1
                ENDIF

                lb_send = lb_send_new
                ub_send = ub_send_new
                lb_recv = lb_recv_new
                ub_recv = ub_recv_new

                nn = PRODUCT ( ub_send - lb_send + 1 )
                ALLOCATE ( send_buf_3d ( lb_send(1):ub_send(1), lb_send(2):ub_send(2), lb_send(3):ub_send(3) ), STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_buf_3d",nn)
                nn = PRODUCT ( ub_recv - lb_recv + 1 )
                ALLOCATE ( recv_buf_3d ( lb_recv(1):ub_recv(1), lb_recv(2):ub_recv(2), lb_recv(3):ub_recv(3) ), STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_buf_3d",nn)
                
                send_buf_3d = rs % r ( lb_send(1):ub_send(1), lb_send(2):ub_send(2), lb_send(3):ub_send(3) )

                CALL mp_sendrecv ( send_buf_3d, dest, recv_buf_3d, source, rs % group )

                ! only some procs may need later shifts
                IF(ub_recv(idir).GE.lb_recv(idir))THEN

                   ! Sum the data in the RS Grid
                   rs % r ( lb_recv(1):ub_recv(1), lb_recv(2):ub_recv(2), lb_recv(3):ub_recv(3) ) = &
                        rs % r ( lb_recv(1):ub_recv(1), lb_recv(2):ub_recv(2), lb_recv(3):ub_recv(3) ) + recv_buf_3d ( :, :, : )
                END IF

                DEALLOCATE ( send_buf_3d, STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_buf_3d")
                DEALLOCATE ( recv_buf_3d, STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_buf_3d")

             END DO

             DEALLOCATE ( dshifts, STAT=ierr )
             IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","dshifts")
             DEALLOCATE ( ushifts, STAT=ierr )
             IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","ushifts")

          END IF
       END DO

       ! This is the real redistribution

       ALLOCATE ( bounds ( 0:pw % pw_grid % para % group_size - 1, 1:2 ), STAT=ierr )
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","bounds",nn)

       ! work out the pw grid points each proc holds
       DO i = 0, pw % pw_grid % para % group_size - 1
          bounds ( i , : ) = get_limit ( pw % pw_grid % npts ( 1 ), pw % pw_grid % para % rs_dims ( 1 ), i )
          bounds ( i , : ) = bounds ( i , : ) - pw % pw_grid % npts (1) / 2 - 1
       ENDDO

       ALLOCATE ( send_tasks ( 0:pw % pw_grid % para % group_size -1,1:6 ), STAT=ierr )
       ALLOCATE ( send_sizes ( 0:pw % pw_grid % para % group_size -1 ), STAT=ierr )
       ALLOCATE ( send_disps ( 0:pw % pw_grid % para % group_size -1 ), STAT=ierr )

       ALLOCATE ( recv_tasks ( 0:pw % pw_grid % para % group_size -1,1:6 ), STAT=ierr )
       ALLOCATE ( recv_sizes ( 0:pw % pw_grid % para % group_size -1 ), STAT=ierr )
       ALLOCATE ( recv_disps ( 0:pw % pw_grid % para % group_size -1 ), STAT=ierr )

       send_tasks = 0
       send_tasks(:,1)=1
       send_tasks(:,2)=0
       send_tasks(:,3)=1
       send_tasks(:,4)=0
       send_tasks(:,5)=1
       send_tasks(:,6)=0
       send_sizes = 0

       recv_tasks = 0
       recv_tasks(:,1)=1
       recv_tasks(:,2)=0
       send_tasks(:,3)=1
       send_tasks(:,4)=0
       send_tasks(:,5)=1
       send_tasks(:,6)=0
       recv_sizes = 0

       send_disps = 0
       recv_disps = 0

       my_rs_rank = rs % my_pos
       my_pw_rank = pw % pw_grid % para % rs_mpo

       ! find the processors that should hold our data
       ! should be part of the rs grid type
       DO i = 0, rs % group_size -1

          CALL mp_cart_coords(rs % group, i, coords)
          !calculate the real rs grid points on each processor
          DO idir = 1,3
             pos (:) = get_limit ( rs %  npts ( idir ), rs % group_dim ( idir ), coords(idir) )
             pos (:) = pos (:) - rs %  npts (idir) / 2 - 1
             lb_send(idir) = pos(1)
             ub_send(idir) = pos(2)                
          ENDDO

          IF(i==my_rs_rank)THEN
             lb_recv(:) = lb_send(:)
             ub_recv(:) = ub_send(:)
          ENDIF

          DO j = 0, pw % pw_grid % para % group_size - 1
          
             IF (lb_send(1) .GT.bounds(j,2)) CYCLE
             IF (ub_send(1) .LT.bounds(j,1)) CYCLE 

             IF(i==my_rs_rank) THEN
                send_tasks(j,1)= MAX(lb_send(1),bounds(j,1))
                send_tasks(j,2)= MIN(ub_send(1),bounds(j,2))
                send_tasks(j,3)= lb_send(2)
                send_tasks(j,4)= ub_send(2)
                send_tasks(j,5)= lb_send(3)
                send_tasks(j,6)= ub_send(3)
                send_sizes(j)  = (send_tasks( j ,2)-send_tasks( j ,1)+1)* &
                     (send_tasks( j ,4)-send_tasks( j ,3)+1)*(send_tasks( j ,6)-send_tasks( j ,5)+1)
             ENDIF

             IF(j==my_rs_rank) THEN
                recv_tasks(i,1)= MAX(lb_send(1),bounds(j,1))
                recv_tasks(i,2)= MIN(ub_send(1),bounds(j,2))
                recv_tasks(i,3)= lb_send(2)
                recv_tasks(i,4)= ub_send(2)
                recv_tasks(i,5)= lb_send(3)
                recv_tasks(i,6)= ub_send(3)
                recv_sizes(i)  = (recv_tasks(i,2)-recv_tasks( i ,1)+1)* &
                     (recv_tasks( i ,4)-recv_tasks( i ,3)+1)*(recv_tasks( i ,6)-recv_tasks( i ,5)+1)
             ENDIF
          ENDDO
       ENDDO

       send_disps(0)=0
       recv_disps(0)=0
       DO i = 1, pw % pw_grid % para % group_size - 1
          send_disps(i) = send_disps(i-1) + send_sizes(i-1)
          recv_disps(i) = recv_disps(i-1) + recv_sizes(i-1)
       ENDDO

       CPPrecondition(SUM(send_sizes)==PRODUCT(ub_recv - lb_recv + 1),cp_failure_level,routineP,error,failure)

       ALLOCATE ( send_buf ( SUM(send_sizes) ), STAT=ierr )
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_buf",nn)
       ALLOCATE ( recv_buf ( SUM(recv_sizes) ), STAT=ierr )
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_buf",nn)
       
       send_buf = 0
       recv_buf = 0

       ! do packing
       j=0
       DO i = 0, rs % group_size - 1
          k=0
          DO z = send_tasks(i,5) , send_tasks(i,6)
             DO y = send_tasks(i,3) , send_tasks(i,4)
                DO x = send_tasks(i,1) , send_tasks(i,2)
                   j=j+1
                   k=k+1
                   send_buf(j)= rs % r (x,y,z)
                ENDDO
             ENDDO
          ENDDO
          CPPrecondition(send_sizes(i)==k,cp_failure_level,routineP,error,failure)
       ENDDO       

       ! do the communication

       CALL mp_alltoall(send_buf, send_sizes, send_disps, recv_buf, recv_sizes, recv_disps, rs % group)

       ! do unpacking
       j=0
       DO i = 0, rs % group_size - 1
          k=0
          DO z = recv_tasks(i,5) , recv_tasks(i,6)
             DO y = recv_tasks(i,3) , recv_tasks(i,4)
                DO x = recv_tasks(i,1) , recv_tasks(i,2)
                   j=j+1
                   k=k+1
                   pw % cr3d (x,y,z) = recv_buf(j) 
                ENDDO
             ENDDO
          ENDDO
          CPPrecondition(recv_sizes(i)==k,cp_failure_level,routineP,error,failure)
       ENDDO

       DEALLOCATE ( send_buf, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_buf")
       DEALLOCATE ( recv_buf, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_buf")

       DEALLOCATE ( send_tasks, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_tasks")
       DEALLOCATE ( send_sizes, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_sizes")
       DEALLOCATE ( send_disps, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_disps")
       DEALLOCATE ( recv_tasks, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_tasks")
       DEALLOCATE ( recv_sizes, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_sizes")
       DEALLOCATE ( recv_disps, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_disps")

       ! safety check, to be removed once we're absolute sure the routine is correct
       pw_sum=pw_integrate_function(pw)
       IF (ABS(pw_sum-rs_sum)/MAX(1.0_dp,ABS(pw_sum),ABS(rs_sum))>EPSILON(rs_sum)*100) THEN
           WRITE(error_string,'(A,6(1X,I4.4),3F20.16)') "rs_pw_transfer_distributed", rs %  npts, rs % group_dim,&
            pw_sum,rs_sum,ABS(pw_sum-rs_sum)
           CALL stop_program(error_string, &
                "Please report this bug... quick workaround: use DISTRIBUTION_TYPE REPLICATED")
       ENDIF

    ELSE

       ! pw to rs transfer

       CALL rs_grid_zero( rs )

       ! This is the real redistribution

       ALLOCATE ( bounds ( 0:pw % pw_grid % para % group_size - 1, 1:2 ), STAT=ierr )
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","bounds",nn)

       ! work out the pw grid points each proc holds
       DO i = 0, pw % pw_grid % para % group_size - 1
          bounds ( i , : ) = get_limit ( pw % pw_grid % npts ( 1 ), pw % pw_grid % para % rs_dims ( 1 ), i )
          bounds ( i , : ) = bounds ( i , : ) - pw % pw_grid % npts (1) / 2 - 1
       ENDDO

       ALLOCATE ( send_tasks ( 0:pw % pw_grid % para % group_size -1,1:6 ), STAT=ierr )
       ALLOCATE ( send_sizes ( 0:pw % pw_grid % para % group_size -1 ), STAT=ierr )
       ALLOCATE ( send_disps ( 0:pw % pw_grid % para % group_size -1 ), STAT=ierr )

       ALLOCATE ( recv_tasks ( 0:pw % pw_grid % para % group_size -1,1:6 ), STAT=ierr )
       ALLOCATE ( recv_sizes ( 0:pw % pw_grid % para % group_size -1 ), STAT=ierr )
       ALLOCATE ( recv_disps ( 0:pw % pw_grid % para % group_size -1 ), STAT=ierr )

       send_tasks = 0
       send_tasks(:,1)=1
       send_tasks(:,2)=0
       send_tasks(:,3)=1
       send_tasks(:,4)=0
       send_tasks(:,5)=1
       send_tasks(:,6)=0
       send_sizes = 0

       recv_tasks = 0
       recv_tasks(:,1)=1
       recv_tasks(:,2)=0
       send_tasks(:,3)=1
       send_tasks(:,4)=0
       send_tasks(:,5)=1
       send_tasks(:,6)=0
       recv_sizes = 0

       my_rs_rank = rs % my_pos
       my_pw_rank = pw % pw_grid % para % rs_mpo

       ! find the processors that should hold our data
       ! should be part of the rs grid type

       DO i = 0, pw % pw_grid % para % group_size -1

          CALL mp_cart_coords(rs % group, i, coords)
          !calculate the real rs grid points on each processor
          DO idir = 1,3
             pos (:) = get_limit ( rs %  npts ( idir ), rs % group_dim ( idir ), coords(idir) )
             pos (:) = pos (:) - rs %  npts (idir) / 2 - 1
             lb_send(idir) = pos(1)
             ub_send(idir) = pos(2)                
          ENDDO

          IF(i==my_rs_rank)THEN
             lb_recv(:) = lb_send(:)
             ub_recv(:) = ub_send(:)
          ENDIF

          DO j = 0, pw % pw_grid % para % group_size - 1
          
             IF (lb_send(1) .GT.bounds(j,2)) CYCLE
             IF (ub_send(1) .LT.bounds(j,1)) CYCLE 

             ! this is the reverse of rs2pw: what were the sends are now the recvs
             IF(i==my_rs_rank) THEN
                   recv_tasks(j,1)= MAX(lb_send(1),bounds(j,1))
                   recv_tasks(j,2)= MIN(ub_send(1),bounds(j,2))
                   recv_tasks(j,3)= lb_send(2)
                   recv_tasks(j,4)= ub_send(2)
                   recv_tasks(j,5)= lb_send(3)
                   recv_tasks(j,6)= ub_send(3)
                   recv_sizes(j)  = (recv_tasks( j ,2)-recv_tasks( j ,1)+1)* &
                        (recv_tasks( j ,4)-recv_tasks( j ,3)+1)*(recv_tasks( j ,6)-recv_tasks( j ,5)+1)
             ENDIF

             IF(j==my_rs_rank) THEN
                   send_tasks(i,1)= MAX(lb_send(1),bounds(j,1))
                   send_tasks(i,2)= MIN(ub_send(1),bounds(j,2))
                   send_tasks(i,3)= lb_send(2)
                   send_tasks(i,4)= ub_send(2)
                   send_tasks(i,5)= lb_send(3)
                   send_tasks(i,6)= ub_send(3)
                   send_sizes(i)  = (send_tasks(i,2)-send_tasks( i ,1)+1)* &
                        (send_tasks( i ,4)-send_tasks( i ,3)+1)*(send_tasks( i ,6)-send_tasks( i ,5)+1)
             ENDIF
          ENDDO
       ENDDO

       send_disps(0)=0
       recv_disps(0)=0
       DO i = 1, pw % pw_grid % para % group_size - 1
          send_disps(i) = send_disps(i-1) + send_sizes(i-1)
          recv_disps(i) = recv_disps(i-1) + recv_sizes(i-1)
       ENDDO

       CPPrecondition(SUM(recv_sizes)==PRODUCT(ub_recv - lb_recv + 1),cp_failure_level,routineP,error,failure)

       ALLOCATE ( send_buf ( SUM(send_sizes) ), STAT=ierr )
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_buf",nn)
       ALLOCATE ( recv_buf ( SUM(recv_sizes) ), STAT=ierr )
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_buf",nn)
       
       send_buf = 0
       recv_buf = 0

       ! do packing
       j=0
       DO i = 0, rs % group_size - 1
          k=0
          DO z = send_tasks(i,5) , send_tasks(i,6)
             DO y = send_tasks(i,3) , send_tasks(i,4)
                DO x = send_tasks(i,1) , send_tasks(i,2)
                   j=j+1
                   k=k+1
                   send_buf(j)= pw % cr3d (x,y,z)
                ENDDO
             ENDDO
          ENDDO
          CPPrecondition(send_sizes(i)==k,cp_failure_level,routineP,error,failure)
       ENDDO       

       ! do the communication

       CALL mp_alltoall(send_buf, send_sizes, send_disps, recv_buf, recv_sizes, recv_disps, rs % group)

       ! do unpacking
       j=0
       DO i = 0, rs % group_size - 1
          k=0
          DO z = recv_tasks(i,5) , recv_tasks(i,6)
             DO y = recv_tasks(i,3) , recv_tasks(i,4)
                DO x = recv_tasks(i,1) , recv_tasks(i,2)
                   j=j+1
                   k=k+1
                   rs % r (x,y,z) = recv_buf(j) 
                ENDDO
             ENDDO
          ENDDO
          CPPrecondition(recv_sizes(i)==k,cp_failure_level,routineP,error,failure)
       ENDDO

 
       DEALLOCATE ( send_buf, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_buf")
       DEALLOCATE ( recv_buf, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_buf")

       DEALLOCATE ( send_tasks, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_tasks")
       DEALLOCATE ( send_sizes, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_sizes")
       DEALLOCATE ( send_disps, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_disps")
       DEALLOCATE ( recv_tasks, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_tasks")
       DEALLOCATE ( recv_sizes, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_sizes")
       DEALLOCATE ( recv_disps, STAT=ierr)
       IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_disps")

       ! now pass wings around

       DO idir = 1,3

          IF ( rs % perd (idir) /= 1) THEN

             ALLOCATE ( dshifts ( 0:rs % neighbours (idir ) ), STAT=ierr )
             ALLOCATE ( ushifts ( 0:rs % neighbours (idir ) ), STAT=ierr )

             ushifts = 0
             dshifts = 0

             DO n_shifts = 1, rs % neighbours(idir)

                ! need to take into account the possible varying widths of neighbouring cells
                ! ushifts and dshifts hold the real size of the neighbouring cells

                position = MODULO ( rs % group_coor (idir) - n_shifts, rs % group_dim (idir))
                neighbours = get_limit ( rs % npts ( idir ), rs % group_dim ( idir ), position )
                dshifts(n_shifts) = dshifts(n_shifts-1) + ( neighbours (2) - neighbours (1) + 1 ) 

                position = MODULO ( rs % group_coor (idir) + n_shifts, rs % group_dim (idir))
                neighbours = get_limit ( rs % npts ( idir ), rs % group_dim ( idir ), position )
                ushifts(n_shifts) = ushifts(n_shifts-1) + ( neighbours (2) - neighbours (1) + 1 )


                ! The border data has to be send/received from the neighbors
                ! First we calculate the source and destination processes for the shift
                ! The first shift is "downwards"

                CALL mp_cart_shift ( rs % group, idir -1, -1 * n_shifts, source, dest )

                   lb_send ( : ) = rs % lb_local ( : )
                   ub_send ( : ) = rs % ub_local ( : )
                   lb_recv ( : ) = rs % lb_local ( : )
                   ub_recv ( : ) = rs % ub_local ( : )

                IF( dshifts ( n_shifts - 1 ) .LE. rs % border ) THEN
                   lb_send ( idir ) = lb_send ( idir ) + rs % border 
                   ub_send ( idir ) = MIN (ub_send (idir) - rs% border,&
                        lb_send ( idir ) + rs % border - 1 - dshifts (n_shifts - 1))

                   lb_recv ( idir ) = ub_recv ( idir ) - rs % border + 1 + ushifts ( n_shifts - 1 ) 
                   ub_recv ( idir ) = MIN (ub_recv ( idir ),&
                        ub_recv ( idir ) - rs%border + ushifts ( n_shifts )) 
                ELSE
                   lb_send( idir ) = 0
                   ub_send( idir ) = lb_send( idir ) - 1
                   lb_recv( idir ) = 0
                   ub_recv( idir ) = lb_recv( idir ) - 1
                ENDIF

                nn = PRODUCT ( ub_send - lb_send + 1 )
                ALLOCATE ( send_buf_3d ( lb_send(1):ub_send(1), lb_send(2):ub_send(2), lb_send(3):ub_send(3) ), STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_buf_3d",nn)
                nn = PRODUCT ( ub_recv - lb_recv + 1 )
                ALLOCATE ( recv_buf_3d ( lb_recv(1):ub_recv(1), lb_recv(2):ub_recv(2), lb_recv(3):ub_recv(3) ), STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_buf_3d",nn)
                
                send_buf_3d = rs % r ( lb_send(1):ub_send(1), lb_send(2):ub_send(2), lb_send(3):ub_send(3) )

                CALL mp_sendrecv ( send_buf_3d, dest, recv_buf_3d, source, rs % group )

                ! only some procs may need later shifts
                IF(ub_recv(idir).GE.lb_recv(idir))THEN

                   ! Sum the data in the RS Grid
                   rs % r ( lb_recv(1):ub_recv(1), lb_recv(2):ub_recv(2), lb_recv(3):ub_recv(3) ) =  recv_buf_3d ( :, :, : )

                END IF

                DEALLOCATE ( send_buf_3d, STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_buf_3d")
                DEALLOCATE ( recv_buf_3d, STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_buf_3d")


                ! Now for the other direction

                CALL mp_cart_shift ( rs % group, idir -1 , n_shifts, source, dest )

                   lb_send ( : ) = rs % lb_local ( : )
                   ub_send ( : ) = rs % ub_local ( : )
                   lb_recv ( : ) = rs % lb_local ( : )
                   ub_recv ( : ) = rs % ub_local ( : )

                IF( ushifts ( n_shifts - 1 ) .LE. rs % border ) THEN
                   ub_send ( idir ) = ub_send ( idir ) - rs % border 
                   lb_send ( idir ) = MAX (lb_send( idir) + rs % border,&
                        ub_send ( idir ) - rs % border + 1 + ushifts (n_shifts - 1) )

                   ub_recv ( idir ) = lb_recv ( idir ) + rs % border - 1 - dshifts ( n_shifts - 1 ) 
                   lb_recv ( idir ) = MAX (lb_recv ( idir ),&
                        lb_recv ( idir ) + rs % border - dshifts ( n_shifts )) 
                ELSE
                   lb_send( idir ) = 0
                   ub_send( idir ) = lb_send( idir ) - 1
                   lb_recv( idir ) = 0
                   ub_recv( idir ) = lb_recv( idir ) - 1
                ENDIF

                nn = PRODUCT ( ub_send - lb_send + 1 )
                ALLOCATE ( send_buf_3d ( lb_send(1):ub_send(1), lb_send(2):ub_send(2), lb_send(3):ub_send(3) ), STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_buf_3d",nn)
                nn = PRODUCT ( ub_recv - lb_recv + 1 )
                ALLOCATE ( recv_buf_3d ( lb_recv(1):ub_recv(1), lb_recv(2):ub_recv(2), lb_recv(3):ub_recv(3) ), STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_buf_3d",nn)

                send_buf_3d = rs % r ( lb_send(1):ub_send(1), lb_send(2):ub_send(2), lb_send(3):ub_send(3) )

                   CALL mp_sendrecv ( send_buf_3d, dest, recv_buf_3d, source, rs % group )

                ! only some procs may need later shifts
                IF(ub_recv(idir).GE.lb_recv(idir))THEN

                   ! Add the data to the RS Grid
                   rs % r ( lb_recv(1):ub_recv(1), lb_recv(2):ub_recv(2), lb_recv(3):ub_recv(3) ) =  recv_buf_3d ( :, :, : )
                END IF

                DEALLOCATE ( send_buf_3d, STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","send_buf_3d")
                DEALLOCATE ( recv_buf_3d, STAT=ierr )
                IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","recv_buf_3d")
                   
             END DO

             DEALLOCATE ( ushifts, STAT=ierr)
             DEALLOCATE ( dshifts, STAT=ierr)

          END IF
       END DO       
    END IF

    DEALLOCATE(rank_pw2rs, STAT=ierr )
    IF ( ierr /= 0 ) CALL stop_memory ( "rs_pw_transfer","rank_pw2rs")


  END SUBROUTINE rs_pw_transfer_distributed

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
!*****
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
!!****** realspace_grid_types/rs_grid_set_box [1.0] *
!!
!!   NAME
!!     rs_grid_set_box
!!
!!   FUNCTION
!!     Set box matrix info for real space grid
!!     This is needed for variable cell simulations
!!
!!   AUTHOR
!!     JGH (15-May-2007)
!!
!!   MODIFICATION HISTORY
!!     none
!*****
!******************************************************************************

SUBROUTINE rs_grid_set_box ( pw_grid, rs, error )

    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_grid_set_box', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, id, n
    LOGICAL                                  :: failure

!-----------------------------------------------------------------------------!

    CPPrecondition(ASSOCIATED(pw_grid),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rs),cp_failure_level,routineP,error,failure)
    CPPrecondition(rs%grid_id==pw_grid%id_nr,cp_failure_level,routineP,error,failure)
    rs % dh = pw_grid%dh
    rs % dh_inv = pw_grid%dh_inv

END SUBROUTINE rs_grid_set_box

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
!*****
!!**************************************************************************

SUBROUTINE rs_grid_retain(rs_grid, error)
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    TYPE(cp_error_type), INTENT(inout)       :: error

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
!*****
!!**************************************************************************

SUBROUTINE rs_grid_release(rs_grid,error)
    TYPE(realspace_grid_type), POINTER       :: rs_grid
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_grid_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(rs_grid)) THEN
       CPPreconditionNoFail(rs_grid%ref_count>0,cp_failure_level,routineP,error)
       rs_grid%ref_count=rs_grid%ref_count-1
       IF (rs_grid%ref_count==0) THEN

         allocated_rs_grid_count=allocated_rs_grid_count-1
         DEALLOCATE ( rs_grid % r, STAT = stat )
         IF ( stat /= 0 ) CALL stop_memory ( "rs_grid_release","rs % r" )
         DEALLOCATE ( rs_grid % px , STAT = stat )
         IF ( stat /= 0 ) CALL stop_memory ( "rs_grid_release","rs % px" )
         DEALLOCATE ( rs_grid % py , STAT = stat )
         IF ( stat /= 0 ) CALL stop_memory ( "rs_grid_release","rs % py" )
         DEALLOCATE ( rs_grid % pz , STAT = stat )
         IF ( stat /= 0 ) CALL stop_memory ( "rs_grid_release","rs % pz" )
         IF ( rs_grid % parallel ) THEN
            ! release the group communicator
            CALL mp_comm_free ( rs_grid % group )
         END IF

         DEALLOCATE(rs_grid, stat=stat)
         CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF
    NULLIFY(rs_grid)
END SUBROUTINE rs_grid_release

!***************************************************************************
!!****f* realspace_grid_types/rs_find_node  
!!
!!   NAME
!!     rs_find_node
!!
!!   FUNCTION
!!     determines the rank of the processor who's real rs grid contains point
!!
!!   AUTHOR
!!     MattW
!!
!!   MODIFICATION HISTORY
!!     11.2007 created [MattW]
!*****
!!**************************************************************************

SUBROUTINE rs_find_node(rs,point,dest)

    TYPE(realspace_grid_type), POINTER       :: rs
    INTEGER, DIMENSION(3), INTENT(IN)        :: point
    INTEGER, INTENT(OUT)                     :: dest

    CHARACTER(len=*), PARAMETER :: routineN = 'rs_find_node', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dir, i
    INTEGER, DIMENSION(2)                    :: pos
    INTEGER, DIMENSION(3)                    :: coords, lb, ub
    TYPE(cp_error_type)                      :: error

!  CPPrecondition(rs%distributed,cp_failure_level,routineP,error,failure)

  dest = -1

  DO i = 0, rs % group_size - 1
     
     CALL mp_cart_coords(rs % group, i, coords)
     !calculate the real rs grid points on each processor
     DO dir = 1,3
        pos (:) = get_limit ( rs %  npts ( dir ), rs % group_dim ( dir ), coords(dir) )
        pos (:) = pos (:) - rs %  npts (dir) / 2 - 1
        lb(dir) = pos(1)
        ub(dir) = pos(2)                
     ENDDO

     ! check whether the point is local to proc i
     IF ( ALL (point.GE.lb) .AND. ALL (point.LE.ub ) ) THEN
        dest = i
        EXIT
     ENDIF

  ENDDO

  CPPostconditionNoFail(dest.GE.0 .AND. dest.LE.rs % group_size,cp_warning_level,routineP,error)

END SUBROUTINE rs_find_node

END MODULE realspace_grid_types

!******************************************************************************
