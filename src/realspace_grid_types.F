!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
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
!!     none
!!
!!   NOTES
!!     Basic type for real space grid methods
!!
!!   SOURCE
!******************************************************************************

MODULE realspace_grid_types
  
  USE kinds, ONLY : dbl
  USE message_passing, ONLY : mp_cart_create, mp_comm_free, mp_comm_dup, &
      mp_sum, mp_max, mp_min, mp_cart_sub, mp_cart_shift, mp_sendrecv, &
      mp_bcast
  USE pw_types, ONLY : pw_type, REALDATA1D, COMPLEXDATA1D, &
      REALDATA3D, COMPLEXDATA3D
  USE pw_grid_types, ONLY : pw_grid_type
  USE simulation_cell, ONLY : cell_type
  USE termination, ONLY : stop_memory, stop_program
  USE util, ONLY : get_limit
  USE timings, ONLy : timeset,timestop
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: realspace_grid_type, rs_grid_allocate, rs_grid_deallocate, &
            rs_grid_setup, rs_pw_transfer
  
  TYPE realspace_grid_type
     INTEGER :: identifier                               ! tag of this grid
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

   lb => rs % lb_local
   ub => rs % ub_local

   ALLOCATE ( rs % r (lb(1):ub(1),lb(2):ub(2),lb(3):ub(3)), STAT = ierr )
   IF ( ierr /= 0 ) CALL stop_memory ( "rs_grid_deallocate","rs % r", &
                                       rs % ngpts_local )
   ALLOCATE ( rs % px ( rs % npts ( 1 ) ), STAT = ierr )
   IF ( ierr /= 0 ) CALL stop_memory ( "rs_grid_deallocate","rs % px", &
                                       rs % npts ( 1 ) )
   ALLOCATE ( rs % py ( rs % npts ( 2 ) ), STAT = ierr )
   IF ( ierr /= 0 ) CALL stop_memory ( "rs_grid_deallocate","rs % py", &
                                       rs % npts ( 2 ) )
   ALLOCATE ( rs % pz ( rs % npts ( 3 ) ), STAT = ierr )
   IF ( ierr /= 0 ) CALL stop_memory ( "rs_grid_deallocate","rs % pz", &
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
   call timeset("rs_grid_dealloc",'I',' ',handle)
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
   call timestop(0.0_dbl,handle)

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
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE rs_grid_setup_1 ( rs, pw_grid, nsmax, iounit )

!  Arguments
   TYPE ( realspace_grid_type ), INTENT (INOUT) :: rs
   TYPE ( pw_grid_type ), INTENT ( IN ) :: pw_grid
   INTEGER, DIMENSION ( 3 ), INTENT ( IN ) :: nsmax 
   INTEGER, OPTIONAL, INTENT ( IN ) :: iounit

!  Local
   INTEGER :: dir, nmin, ngmax, ip, n_slices, pos ( 2 ), lb ( 2 ), nn, i
   REAL ( dbl ) :: pp ( 3 )

!-----------------------------------------------------------------------------!
   rs % dr = pw_grid%dr
   rs % identifier = pw_grid % identifier
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
     ! x is the main distribution direction of pw grids in real space
     ! we distribute real space grids accordingly
     dir = 1
     rs % direction = 1
     ! the minimum slice thickness is half of the small box size
     ! this way there is only overlap with direct neighbors
     nmin = ( nsmax ( dir ) + 1 ) / 2 
     ! maximum number of slices
     ngmax = pw_grid % npts ( dir ) / nmin
     ! now we reduce this number until we find a divisor of the total number
     ! of processors
     !
!cdeb
!the code works currently only for this setting
     ngmax=1
!cdeb
     n_slices = 1
     DO ip = ngmax, 1, -1
       IF ( MOD ( pw_grid % para % group_size, ip ) == 0 ) THEN
         n_slices = ip
         EXIT
       END IF
     END DO
     ! global grid dimensions are still the same
     rs % npts = pw_grid % npts
     rs % ngpts = PRODUCT ( rs % npts )
     rs % lb = pw_grid % bounds ( 1, : )
     rs % ub = pw_grid % bounds ( 2, : )
     rs % group_dim ( 1 ) = n_slices
     rs % group_dim ( 2 ) = pw_grid % para % group_size / n_slices
     IF ( n_slices == 1 ) THEN
       ! CASE 1 : only one slice: we do not need overlapping reagions and special
       !          recombination of the total density
       rs % perd = 1
       rs % border = 0
       rs % npts_local = pw_grid % npts
       rs % ngpts_local = PRODUCT ( rs % npts )
       rs % lb_local = pw_grid % bounds ( 1, : )
       rs % ub_local = pw_grid % bounds ( 2, : )
       rs % parallel = .TRUE.
       CALL mp_comm_dup( pw_grid % para % group , rs % group )
       rs % group_head = pw_grid % para % group_head
       rs % group_coor ( 1 ) = 0
       rs % group_coor ( 2 ) = pw_grid % para % my_pos
     ELSE
       ! CASE 2 : general case
       ! periodicity is no longer enforced along direction dir
       rs % perd = 1
       rs % perd ( dir ) = 0
       ! we keep a border on nmin points
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
             WRITE ( iounit, '(A,T72,I3,A)' ) &
             " RS_GRID: Real space distribution over ", rs % group_dim ( 1 ), " group"
           ELSE
             WRITE ( iounit, '(A,T71,I3,A)' ) &
             " RS_GRID: Real space distribution over ", rs % group_dim ( 1 ), " groups"
           END IF
           WRITE ( iounit, '(A,T71,I10)' ) &
             " RS_GRID: Group size ", rs % group_dim ( 2 )
           WRITE ( iounit, '(A,T71,I10)' ) &
             " RS_GRID: Real space distribution along direction ", dir
           WRITE ( iounit, '(A,T71,I10)' ) &
             " RS_GRID: Border size ", rs % border
         END IF
         nn = rs % npts_local ( dir )
         CALL mp_sum ( nn, rs % group )
         pp ( 1 ) = REAL ( nn, dbl ) / REAL ( pw_grid % para % group_size, dbl )
         nn = rs % npts_local ( dir )
         CALL mp_max ( nn, pw_grid % para % group )
         pp ( 2 ) = REAL ( nn, dbl )
         nn = rs % npts_local ( dir )
         CALL mp_min ( nn, pw_grid % para % group )
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
   INTEGER, DIMENSION ( 3 ), INTENT ( IN ) :: nsmax 
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
!!   FUNCTION
!!     Copy a function from/to a PW grid type to/from a real
!!     space type 
!!
!!   AUTHOR
!!     JGH (18-Mar-2001)
!!
!!   MODIFICATION HISTORY
!!     none
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

!-----------------------------------------------------------------------------!
   call timeset("rs_pw_transfer",'I',' ',handle)
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
       !Only one type of grid, this is easy

       pbo => pw % pw_grid % bounds_local
       IF ( dir == "FORWARD" ) THEN

         !We first sum the full array, this could be optimised 
         CALL mp_sum ( rs % r, rs % group )

         !Now copy the local section from the RS Grid to the PW Grid
         IF ( pw % in_use == REALDATA3D ) THEN
           pw % cr3d ( pbo(1,1):pbo(2,1), pbo(1,2):pbo(2,2), &
                       pbo(1,3):pbo(2,3) ) = &
              rs % r ( pbo(1,1):pbo(2,1), pbo(1,2):pbo(2,2), &
                       pbo(1,3):pbo(2,3) )
         ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
           pw % cc3d ( pbo(1,1):pbo(2,1), pbo(1,2):pbo(2,2), &
                       pbo(1,3):pbo(2,3) ) = &
                CMPLX ( rs % r ( pbo(1,1):pbo(2,1), pbo(1,2):pbo(2,2), &
                       pbo(1,3):pbo(2,3) ), KIND = dbl )
         ELSE
           CALL stop_program ( "rs_pw_transfer", "Pw type not compatible" )
         END IF
 
       ELSEIF ( dir == "BACKWARD" ) THEN

         rs % r = 0._dbl
         IF ( pw % in_use == REALDATA3D ) THEN
           rs % r ( pbo(1,1):pbo(2,1), pbo(1,2):pbo(2,2), &
                    pbo(1,3):pbo(2,3) ) = &
                pw % cr3d ( pbo(1,1):pbo(2,1), pbo(1,2):pbo(2,2), &
                            pbo(1,3):pbo(2,3) )
         ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
           rs % r ( pbo(1,1):pbo(2,1), pbo(1,2):pbo(2,2), &
                    pbo(1,3):pbo(2,3) ) = &
                REAL ( pw % cr3d ( pbo(1,1):pbo(2,1), pbo(1,2):pbo(2,2), &
                            pbo(1,3):pbo(2,3) ), dbl )
         ELSE
           CALL stop_program ( "rs_pw_transfer", "Pw type not compatible" )
         END IF

         !We sum the full array, this could be optimised 
         CALL mp_sum ( rs % r, rs % group )

       ELSE
         CALL stop_program ( "rs_pw_transfer", &
              "Parameter dir ="//dir//" not allowed" )
       END IF

     ELSE
       !Both grids are distributed, we have to do some reshuffling of data

       IF ( dir == "FORWARD" ) THEN

         ! First we have to sum up the data within the RS groups
         ! we generate a subgroup that covers all processors with the same RS grid
         subdim ( 1 ) = .FALSE.
         subdim ( 2 ) = .TRUE.
         CALL mp_cart_sub ( rs % group, subdim, subgroup )
         ! Now we do the sum
         CALL mp_sum ( rs % r, subgroup )
         ! And we release the communicator again
         CALL mp_comm_free ( subgroup )
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
     IF ( dir == "FORWARD" ) THEN
       IF ( pw % in_use == REALDATA3D ) THEN
         pw % cr3d = rs % r
       ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
         pw % cc3d = CMPLX ( rs % r, KIND = dbl )
       ELSE
         CALL stop_program ( "rs_pw_transfer", "Pw type not compatible" )
       END IF
     ELSEIF ( dir == "BACKWARD" ) THEN
       IF ( pw % in_use == REALDATA3D ) THEN
         rs % r = pw % cr3d
       ELSEIF ( pw % in_use == COMPLEXDATA3D ) THEN
         rs % r = REAL ( pw % cc3d, dbl )
       ELSE
         CALL stop_program ( "rs_pw_transfer", "Pw type not compatible" )
       END IF
     ELSE
       CALL stop_program ( "rs_pw_transfer", &
            "Parameter dir ="//dir//" not allowed" )
     END IF

   END IF
   call timestop(0.0_dbl,handle)

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
!  This code is for MPI only, and uses a lot of redunant communication
!  Optimal code will use one-sided-communication features
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
!  This code is for MPI only, and uses a lot of redunant communication
!  Optimal code will use one-sided-communication features
!-----------------------------------------------------------------------------!

   ! initialize the rs grid to zero
   rs % r = 0._dbl
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

END MODULE realspace_grid_types

!******************************************************************************
