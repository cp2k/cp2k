!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief the type I Discrete Cosine Transfrom (DCT-I)
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
MODULE dct

  USE fast,                            ONLY: copy_cr
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_cart_rank,&
                                             mp_irecv,&
                                             mp_isend,&
                                             mp_request_null,&
                                             mp_wait,&
                                             mp_waitall
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_create,&
                                             pw_grid_setup
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA1D,&
                                             RECIPROCALSPACE,&
                                             pw_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dct'

  TYPE :: dct_type
    INTEGER, DIMENSION(4)        :: dests_expand
    INTEGER, DIMENSION(4)        :: srcs_expand
    INTEGER, DIMENSION(4)        :: flipg_stat
    INTEGER, DIMENSION(4)        :: dests_shrink
    INTEGER                      :: srcs_shrink
    INTEGER, DIMENSION(2,3,4)    :: recv_msgs_bnds
    INTEGER, DIMENSION(2,3)      :: dct_bounds
    INTEGER, DIMENSION(2,3)      :: dct_bounds_local
    INTEGER, DIMENSION(2,3)      :: bounds_shftd
    INTEGER, DIMENSION(2,3)      :: bounds_local_shftd
  END TYPE dct_type

  TYPE dct_msg_type
      PRIVATE
      REAL(dp), DIMENSION(:, :, :), POINTER :: msg => NULL()
  END TYPE dct_msg_type

  PUBLIC dct_type, &
         dct_type_init, &
         setup_dct_pw_grids, &
         pw_shrink, &
         pw_expand, &
         pw_adjust_bounds_local

  INTEGER, PARAMETER, PRIVATE    :: NOT_FLIPPED = 0, &
                                    UD_FLIPPED  = 1, &
                                    LR_FLIPPED  = 2, &
                                    BF_FLIPPED  = 3, &
                                    ROTATED     = 4

  INTEGER, PARAMETER, PRIVATE    :: black = 0, &
                                    red   = 1 

  INTEGER, PARAMETER, PRIVATE    :: nomsg_received = -1, &
                                    zsmsg_received =  0, &
                                    msg_received   =  1
CONTAINS

! *****************************************************************************
!> \brief  Initializes a dct_type
!> \param pw_grid the original plane wave grid
!> \param dct_env dct_type to be initialized
!> \param error cp2k error
!> \par History
!>       08.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dct_type_init(pw_grid, dct_env, error)

    TYPE(pw_grid_type), INTENT(IN), POINTER  :: pw_grid
    TYPE(dct_type), INTENT(INOUT)            :: dct_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dct_type_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    CALL get_dests_srcs_pid(pw_grid, &
              dct_env%dests_expand, dct_env%srcs_expand, dct_env%flipg_stat, &
              dct_env%dests_shrink, dct_env%srcs_shrink, error)
    CALL expansion_bounds(pw_grid, &
              dct_env%srcs_expand, dct_env%flipg_stat, &
              dct_env%bounds_shftd, dct_env%bounds_local_shftd, &
              dct_env%recv_msgs_bnds, dct_env%dct_bounds, &
              dct_env%dct_bounds_local, error)

    CALL timestop(handle)

  END SUBROUTINE dct_type_init

! *****************************************************************************
!> \brief   sets up an extended pw_grid for Discrete Cosine Transform (DCT) 
!>          calculations 
!> \param pw_grid the original plane wave grid
!> \param cell_hmat cell hmat
!> \param dct_pw_grid DCT plane-wave grid
!> \param dct_aux_pw_grid auxiliary DCT plane-wave grid (i.e an expanded pw grid
!>                        whose local bounds are decided internally by cp2k)
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE setup_dct_pw_grids(pw_grid, cell_hmat, dct_pw_grid, dct_aux_pw_grid, error)

    TYPE(pw_grid_type), INTENT(IN), POINTER  :: pw_grid
    REAL(dp), DIMENSION(3, 3), INTENT(IN)    :: cell_hmat
    TYPE(pw_grid_type), INTENT(INOUT), &
      POINTER                                :: dct_pw_grid, dct_aux_pw_grid
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'setup_dct_pw_grids', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blocked, handle, srcs_shrink
    INTEGER, DIMENSION(2, 3)                 :: bounds_local_new, &
                                                bounds_local_shftd, &
                                                bounds_new, bounds_shftd
    INTEGER, DIMENSION(2, 3, 4)              :: recv_msgs_bnds
    INTEGER, DIMENSION(4)                    :: dests_expand, dests_shrink, &
                                                flipg_stat, srcs_expand
    REAL(KIND=dp), DIMENSION(3, 3)           :: hmat2

    CALL timeset(routineN,handle)

    CALL get_dests_srcs_pid(pw_grid, dests_expand, srcs_expand, flipg_stat, dests_shrink, srcs_shrink, error)
    CALL expansion_bounds(pw_grid, srcs_expand, flipg_stat, bounds_shftd, bounds_local_shftd, &
                          recv_msgs_bnds, bounds_new, bounds_local_new, error)
    CALL pw_grid_create(dct_pw_grid, pw_grid%para%rs_group, local=.FALSE., error=error)
    CALL pw_grid_create(dct_aux_pw_grid, pw_grid%para%rs_group, local=.FALSE., error=error)

    hmat2 = 0.0_dp
    hmat2(1,1) = 2*(cell_hmat(1,1) - pw_grid%dr(1))
    hmat2(2,2) = 2*(cell_hmat(2,2) - pw_grid%dr(2))
    hmat2(3,3) = 2*(cell_hmat(3,3) - pw_grid%dr(3))

    ! uses bounds_local_new that is 2*n-2 in size....this is only rarely fft-able by fftsg, and needs fftw3,
    ! where it might use sizes that are not declared available in fft_get_radix.

    IF (pw_grid%para%blocked) THEN
       blocked = 1
    ELSE IF (pw_grid%para%ray_distribution) THEN
       blocked = 0
    END IF

    CALL pw_grid_setup(hmat2, dct_pw_grid, &
                       bounds=bounds_new, &
                       rs_dims = pw_grid%para%rs_dims, &
                       blocked = blocked, &
                       bounds_local=bounds_local_new, &
                       error=error)

    CALL pw_grid_setup(hmat2, dct_aux_pw_grid, &
                       bounds=bounds_new, &
                       rs_dims = pw_grid%para%rs_dims, &
                       blocked = blocked, &
                       error=error)

    CALL timestop(handle)

  END SUBROUTINE setup_dct_pw_grids

! *****************************************************************************
!> \brief Finds the process ids for mpi_isend destiations and mpi_irecv sources
!>   for expanding and shrinking a pw_type data
!> \param pw_grid the original plane wave grid
!> \param dests_expand list of the destination processes (pw_expand)
!> \param srcs_expand list of the source processes (pw_expand)
!> \param flipg_stat flipping status for the received data chunks (pw_expand)
!> \param dests_shrink list of the destination processes (pw_shrink)
!> \param srcs_shrink list of the source proceses (pw_shrink)
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE get_dests_srcs_pid(pw_grid, dests_expand, srcs_expand, flipg_stat, &
                                dests_shrink, srcs_shrink, error)

    TYPE(pw_grid_type), INTENT(IN), POINTER  :: pw_grid
    INTEGER, DIMENSION(4), INTENT(OUT)       :: dests_expand, srcs_expand, &
                                                flipg_stat, dests_shrink
    INTEGER, INTENT(OUT)                     :: srcs_shrink
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'get_dests_srcs_pid', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: group_size, handle, i, j, k, &
                                                rs_dim1, rs_dim2, rs_group, &
                                                rs_mpo, tmp_size1, tmp_size2
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: tmp1_arr, tmp2_arr
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: dests_shrink_all, src_pos1, &
                                                src_pos2, srcs_expand_all
    INTEGER, DIMENSION(2)                    :: rs_dims, rs_pos
    INTEGER, DIMENSION(2, 4)                 :: srcs_coord

    CALL timeset(routineN,handle)

! example: 3x4 process grid
! rs_dim1 = 3 -->  src_pos1 = [1  2; 3 -3; -2 -1]'
! rs_dim2 = 4 -->  src_pos2 = [1  2; 3  4; -4 -3; -2 -1]' 
! => (1,1) receives from (1,1) ; (1,2) ; (2,1) ; (2,2) | flipping status 0 0 0 0
!    (2,4) receives from (3,2) ; (3,1) ; (3,2) ; (3,1) | flipping status 2 2 4 4 
! and so on ...
! or equivalently 
! =>   0   receives from 0 ; 1 ; 4 ; 5 | flipping status 0 0 0 0 
!      7   receives from 9 ; 8 ; 9 ; 8 | flipping status 2 2 4 4
! schematically :
! ij : coordinates in a 2D process grid (starting from 1 just for demonstration) 
! () : to be flipped from left to right
! [] : to be flipped from up to down
! {} : to be rotated 180 degrees
!       11   12  |  13   14  | (14) (13) | (12) (11)
!       21   22  |  23   24  | (24) (23) | (22) (21)
!      ---------------------------------------------
!       31   32  |  33   34  | (34) (33) | (32) (31)
!      [31] [32] | [33] [34] | {34} {33} | {32} {31}
!      ---------------------------------------------
!      [21] [22] | [23] [24] | {24} {23} | {22} {21}
!      [11] [12] | [13] [14] | {14} {13} | {12} {11}

    rs_group = pw_grid%para%rs_group
    rs_mpo = pw_grid%para%rs_mpo
    group_size = pw_grid%para%group_size
    rs_dims = pw_grid%para%rs_dims

    rs_pos = pw_grid%para%rs_pos
    rs_dim1 = rs_dims(1);     rs_dim2 = rs_dims(2)

    IF (MOD(rs_dim1,2) .EQ. 0) THEN
       tmp_size1 = rs_dim1
    ELSE
       tmp_size1 = rs_dim1 + 1
    END IF

    ALLOCATE(tmp1_arr(tmp_size1), src_pos1(2, 0:rs_dim1-1))

    IF (MOD(rs_dim1,2) .EQ. 0) THEN
       tmp1_arr(:) = (/ (i, i=1, rs_dim1) /)
       src_pos1(:,:) = RESHAPE((/ tmp1_arr, - tmp1_arr(tmp_size1:1:-1) /), (/2, rs_dim1/))
    ELSE
       tmp1_arr(:) = (/ (i, i=1, rs_dim1), - rs_dim1 /)
       src_pos1(:,:) = RESHAPE((/ tmp1_arr, - tmp1_arr(tmp_size1-2:1:-1) /), (/2, rs_dim1/))
    END IF
 
    IF (MOD(rs_dim2,2) .EQ. 0) THEN 
       tmp_size2 = rs_dim2
    ELSE
       tmp_size2 = rs_dim2 + 1
    END IF

    ALLOCATE(tmp2_arr(tmp_size2), src_pos2(2, 0:rs_dim2-1))

    IF (MOD(rs_dim2,2) .EQ. 0) THEN
       tmp2_arr(:) = (/ (i, i=1, rs_dim2) /)
       src_pos2(:,:) = RESHAPE((/ tmp2_arr, - tmp2_arr(tmp_size2:1:-1) /), (/2, rs_dim2/))
    ELSE
       tmp2_arr(:) = (/ (i, i=1, rs_dim2), - rs_dim2 /)
       src_pos2(:,:) = RESHAPE((/ tmp2_arr, - tmp2_arr(tmp_size2-2:1:-1) /), (/2, rs_dim2/))
    END IF

! default flipping status
    flipg_stat = NOT_FLIPPED

    srcs_coord(:,1) = (/src_pos1(1,rs_pos(1)) , src_pos2(1,rs_pos(2))/)
    srcs_coord(:,2) = (/src_pos1(1,rs_pos(1)) , src_pos2(2,rs_pos(2))/)
    srcs_coord(:,3) = (/src_pos1(2,rs_pos(1)) , src_pos2(1,rs_pos(2))/)
    srcs_coord(:,4) = (/src_pos1(2,rs_pos(1)) , src_pos2(2,rs_pos(2))/)

    DO k = 1, 4
! convert srcs_coord to pid
       CALL mp_cart_rank(pw_grid%para%rs_group, ABS(srcs_coord(:,k)) - 1, srcs_expand(k))
! find out the flipping status 
       IF ((srcs_coord(1,k) .GT. 0) .AND. (srcs_coord(2,k) .GT. 0)) THEN
          flipg_stat(k) = NOT_FLIPPED
       ELSE IF ((srcs_coord(1,k) .LT. 0) .AND. (srcs_coord(2,k) .GT. 0)) THEN
          flipg_stat(k) = UD_FLIPPED
       ELSE IF ((srcs_coord(1,k) .GT. 0) .AND. (srcs_coord(2,k) .LT. 0)) THEN
          flipg_stat(k) = LR_FLIPPED
       ELSE 
          flipg_stat(k) = ROTATED
       END IF
    END DO

! let all the nodes know about each others srcs_expand list 
    ALLOCATE(srcs_expand_all(4, group_size))
    CALL mp_allgather(srcs_expand, srcs_expand_all, rs_group)
! now scan the srcs_expand_all list and check if I am on the srcs_expand list of the other nodes
! if that the case then I am obliged to send data to those nodes (the nodes are on my dests_expand list) 
    k = 1
    DO i = 1, group_size
       DO j = 1, 4 
          IF (srcs_expand_all(j,i) .EQ.rs_mpo) THEN 
             dests_expand(k) = i - 1
             k = k + 1
          END IF 
       END DO 
    END DO

! find srcs and dests for the reverse procedure :
! initialize dests_shrink and srcs_shrink with invalid process id 
    dests_shrink = -1
    srcs_shrink = -1
! scan the flipping status of the data that I am supposed to receive  
! if the flipping status for a process is NOT_FLIPPED that means I will have to resend 
! data to that process in the reverse procedure (the process is on my dests_shrink list)
    DO i = 1, 4
       IF (flipg_stat(i) .EQ. NOT_FLIPPED) dests_shrink(i) = srcs_expand(i)
    END DO

! let all the nodes know about each others dests_shrink list 
    ALLOCATE(dests_shrink_all(4,group_size))
    CALL mp_allgather(dests_shrink, dests_shrink_all, rs_group)
! now scan the dests_shrink_all list and check if I am on the rev_dest list of any other node 
! if that is the case then I'll receive data from that node (the node is on my rev_src list)
! note that in the shrinking procedure I will receive from only one node  
    DO i = 1, group_size
       DO j = 1, 4
          IF (dests_shrink_all(j,i) .EQ. rs_mpo) THEN 
             srcs_shrink = i - 1
             EXIT
          END IF
       END DO
    END DO
 
    DEALLOCATE(tmp1_arr, tmp2_arr, src_pos1, src_pos2, srcs_expand_all, dests_shrink_all)
    CALL timestop(handle)

  END SUBROUTINE get_dests_srcs_pid

! *****************************************************************************
!> \brief expands a pw_type data to an evenly symmetric pw_type data that is 8 times 
!>   larger than the original one: 
!>   the even symmetry for a 1D sequence of length n is defined as:
!>      1 2 3 ... n-2 n-1 n --> 1 2 3 ... n-2 n-1 n n-1 n-2 ... 3 2 
!>   and is generalized to 3D by applying the same rule in all three directions 
!>
!> \param recv_msgs_bnds bounds of the messages to be received
!> \param dests_expand list of the destination processes
!> \param srcs_expand list of the source processes
!> \param flipg_stat flipping status for the received data chunks
!> \param bounds_shftd bounds of the original grid shifted to have g0 in the middle of the cell
!> \param pw_in the original plane wave data
!> \param pw_expanded the pw data after expansion
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE pw_expand(recv_msgs_bnds, dests_expand, srcs_expand, flipg_stat, bounds_shftd, &
                       pw_in, pw_expanded, error)

    INTEGER, DIMENSION(2, 3, 4), INTENT(IN)  :: recv_msgs_bnds
    INTEGER, DIMENSION(4), INTENT(IN)        :: dests_expand, srcs_expand, &
                                                flipg_stat
    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bounds_shftd
    TYPE(pw_type), INTENT(IN), POINTER       :: pw_in
    TYPE(pw_type), INTENT(INOUT), POINTER    :: pw_expanded
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pw_expand', &
      routineP = moduleN//':'//routineN

    INTEGER :: group_size, handle, i, ind, lb1, lb1_loc, lb1_new, lb2, &
      lb2_loc, lb2_new, lb3, lb3_loc, lb3_new, loc, rs_group, rs_mpo, ub1, &
      ub1_loc, ub1_new, ub2, ub2_loc, ub2_new, ub3, ub3_loc, ub3_new
    INTEGER, DIMENSION(2, 3)                 :: bounds_local_new
    INTEGER, DIMENSION(2, 3, 4)              :: pcs_bnds
    INTEGER, DIMENSION(4)                    :: dest_hist, recv_reqs, &
                                                send_reqs, src_hist
    REAL(dp), DIMENSION(:, :, :), POINTER    :: catd, catd_flipdbf, &
                                                cr3d_xpndd, send_msg
    TYPE(dct_msg_type), DIMENSION(4)         :: pcs, recv_msgs
    TYPE(pw_grid_type), POINTER              :: pw_grid

    CALL timeset(routineN,handle)

    pw_grid => pw_in%pw_grid
    rs_group = pw_grid%para%rs_group
    rs_mpo = pw_grid%para%my_pos
    group_size = pw_grid%para%group_size
   
    bounds_local_new = pw_expanded%pw_grid%bounds_local

    send_reqs = mp_request_null
    recv_reqs = mp_request_null

    send_msg => pw_in%cr3d

    src_hist = -1   ! keeps the history of sources 
    dest_hist = -1  ! keeps the history of destinations 

    DO i = 1, 4
! no need to send to myself or to the destination that I have already sent to
       IF ((dests_expand(i) .NE. rs_mpo) .AND. .NOT. ANY(dest_hist .EQ. dests_expand(i))) THEN
          CALL mp_isend(send_msg, dests_expand(i), rs_group, send_reqs(i))
       END IF
       dest_hist(i) = dests_expand(i)
    END DO

    DO i = 1, 4
       lb1 = recv_msgs_bnds(1,1,i)
       ub1 = recv_msgs_bnds(2,1,i)
       lb2 = recv_msgs_bnds(1,2,i)
       ub2 = recv_msgs_bnds(2,2,i)
       lb3 = recv_msgs_bnds(1,3,i)
       ub3 = recv_msgs_bnds(2,3,i)
! no need to receive from myself
       IF (srcs_expand(i) .EQ. rs_mpo) THEN
          ALLOCATE(recv_msgs(i)%msg(lb1:ub1, lb2:ub2, lb3:ub3))
          recv_msgs(i)%msg = send_msg
! if I have already received data from the source, just use the one from the last time
       ELSE IF (ANY(src_hist .EQ. srcs_expand(i))) THEN
          loc = MINLOC(ABS(src_hist - srcs_expand(i)), 1) 
          lb1_loc = recv_msgs_bnds(1,1,loc)
          ub1_loc = recv_msgs_bnds(2,1,loc)
          lb2_loc = recv_msgs_bnds(1,2,loc)
          ub2_loc = recv_msgs_bnds(2,2,loc)
          lb3_loc = recv_msgs_bnds(1,3,loc)
          ub3_loc = recv_msgs_bnds(2,3,loc)
          ALLOCATE(recv_msgs(i)%msg(lb1_loc:ub1_loc, lb2_loc:ub2_loc, lb3_loc:ub3_loc))
          recv_msgs(i)%msg = recv_msgs(loc)%msg
       ELSE
          ALLOCATE(recv_msgs(i)%msg(lb1:ub1, lb2:ub2, lb3:ub3))
          CALL mp_irecv(recv_msgs(i)%msg, srcs_expand(i), rs_group, recv_reqs(i)) 
          CALL mp_wait(recv_reqs(i))
       END IF
       src_hist(i) = srcs_expand(i)
    END DO
! cleanup mpi_request to prevent memory leak 
    CALL mp_waitall(send_reqs(:))

! flip the received data according on the flipping status
    DO i = 1, 4
       SELECT CASE (flipg_stat(i))
          CASE(NOT_FLIPPED)
              lb1 = recv_msgs_bnds(1,1,i)
              ub1 = recv_msgs_bnds(2,1,i)
              lb2 = recv_msgs_bnds(1,2,i)
              ub2 = recv_msgs_bnds(2,2,i)
              lb3 = recv_msgs_bnds(1,3,i)
              ub3 = recv_msgs_bnds(2,3,i)
              ALLOCATE(pcs(i)%msg(lb1:ub1, lb2:ub2, lb3:ub3))
              pcs(i)%msg = recv_msgs(i)%msg
          CASE(UD_FLIPPED)
              CALL flipud(recv_msgs(i)%msg, pcs(i)%msg, bounds_shftd, error)
          CASE(LR_FLIPPED)
              CALL fliplr(recv_msgs(i)%msg, pcs(i)%msg, bounds_shftd, error)
          CASE(BF_FLIPPED)
              CALL flipbf(recv_msgs(i)%msg, pcs(i)%msg, bounds_shftd, error)
          CASE(ROTATED)
              CALL rot180(recv_msgs(i)%msg, pcs(i)%msg, bounds_shftd, error)
       END SELECT
    END DO
! concatenate the received (flipped) data store the result as catd
! need the bounds of the four pieces for concatenation
    DO i = 1, 4
       pcs_bnds(1,1,i) = LBOUND(pcs(i)%msg,1) 
       pcs_bnds(2,1,i) = UBOUND(pcs(i)%msg,1)
       pcs_bnds(1,2,i) = LBOUND(pcs(i)%msg,2) 
       pcs_bnds(2,2,i) = UBOUND(pcs(i)%msg,2)
       pcs_bnds(1,3,i) = LBOUND(pcs(i)%msg,3) 
       pcs_bnds(2,3,i) = UBOUND(pcs(i)%msg,3)
    END DO

    lb1_new = bounds_local_new(1,1);    ub1_new = bounds_local_new(2,1)
    lb2_new = bounds_local_new(1,2);    ub2_new = bounds_local_new(2,2)
    lb3_new = bounds_local_new(1,3);    ub3_new = bounds_local_new(2,3)
    ind = INT(0.5*(ub3_new+lb3_new+1))
    ALLOCATE(catd(lb1_new:ub1_new, lb2_new:ub2_new, lb3_new:ind))
    DO i = 1, 4
       catd(pcs_bnds(1,1,i):pcs_bnds(2,1,i), &
            pcs_bnds(1,2,i):pcs_bnds(2,2,i), &
            pcs_bnds(1,3,i):pcs_bnds(2,3,i)) = pcs(i)%msg
    END DO

! flip catd from back to front 
    CALL flipbf(catd, catd_flipdbf, bounds_shftd, error)
! concatenate catd and catd_flipdbf to get cr3d_xpndd
    ALLOCATE(cr3d_xpndd(lb1_new:ub1_new, lb2_new:ub2_new, lb3_new:ub3_new))
    cr3d_xpndd(:,:,lb3_new:ind) = catd
    cr3d_xpndd(:,:,ind+1:ub3_new) = catd_flipdbf

    pw_expanded%cr3d = cr3d_xpndd

    DO i = 1, 4
       DEALLOCATE(pcs(i)%msg)
       DEALLOCATE(recv_msgs(i)%msg)
    END DO
    DEALLOCATE(catd, catd_flipdbf, cr3d_xpndd)

    CALL timestop(handle)

  END SUBROUTINE pw_expand

! *****************************************************************************
!> \brief shrinks an evenly symmetric pw_type data to a pw_type data that is 8 
!>        times smaller (the reverse procedure of pw_expand).
!>
!> \param dests_shrink list of the destination processes
!> \param srcs_shrink list of the source proceses
!> \param bounds_local_shftd local bounds of the original grid after shifting
!> \param pw_in the original plane wave data
!> \param pw_shrinked the shrinked plane wave data
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE pw_shrink(dests_shrink, srcs_shrink, bounds_local_shftd, &
                       pw_in, pw_shrinked, error)

    INTEGER, DIMENSION(4), INTENT(IN)        :: dests_shrink
    INTEGER, INTENT(IN)                      :: srcs_shrink
    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bounds_local_shftd
    TYPE(pw_type), INTENT(IN), POINTER       :: pw_in
    TYPE(pw_type), INTENT(INOUT), POINTER    :: pw_shrinked
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pw_shrink', &
      routineP = moduleN//':'//routineN

    COMPLEX(dp), DIMENSION(:, :, :), POINTER :: cc3d
    INTEGER :: group_size, handle, i, in_space, in_use, lb1_orig, lb1_xpnd, &
      lb2_orig, lb2_xpnd, lb3_orig, lb3_xpnd, recv_req, rs_group, rs_mpo, &
      send_lb1, send_lb2, send_lb3, send_req, send_ub1, send_ub2, send_ub3, &
      ub1_orig, ub1_xpnd, ub2_orig, ub2_xpnd, ub3_orig, ub3_xpnd
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: bounds_local_all
    INTEGER, DIMENSION(2, 3)                 :: bounds_local_xpnd
    REAL(dp), DIMENSION(:, :, :), POINTER    :: cr3d, send_crmsg
    TYPE(pw_grid_type), POINTER              :: pw_grid_orig

    CALL timeset(routineN,handle)

    pw_grid_orig => pw_shrinked%pw_grid 
    rs_group = pw_grid_orig%para%rs_group
    rs_mpo = pw_grid_orig%para%my_pos
    group_size = pw_grid_orig%para%group_size
    bounds_local_xpnd = pw_in%pw_grid%bounds_local
    in_space = pw_in%in_space
    in_use = pw_in%in_use

! cosine transform is a real transform. The cosine transfrom of a 3D data must be real and 3D.
    NULLIFY(cr3d, cc3d)
    lb1_xpnd = bounds_local_xpnd(1,1)
    ub1_xpnd = bounds_local_xpnd(2,1)
    lb2_xpnd = bounds_local_xpnd(1,2)
    ub2_xpnd = bounds_local_xpnd(2,2)
    lb3_xpnd = bounds_local_xpnd(1,3)
    ub3_xpnd = bounds_local_xpnd(2,3)
    IF ((in_use .EQ. REALDATA1D) .OR. (in_use .EQ. COMPLEXDATA1D)) THEN
       IF (in_space .EQ. RECIPROCALSPACE) THEN
          ALLOCATE(cr3d(lb1_xpnd:ub1_xpnd, lb2_xpnd:ub2_xpnd, lb3_xpnd:ub3_xpnd))
          ALLOCATE(cc3d(lb1_xpnd:ub1_xpnd, lb2_xpnd:ub2_xpnd, lb3_xpnd:ub3_xpnd))
          cc3d = RESHAPE(pw_in%cc,(/ub1_xpnd-lb1_xpnd+1, ub2_xpnd-lb2_xpnd+1, ub3_xpnd-lb3_xpnd+1/))
          CALL copy_cr(cc3d,cr3d)
          DEALLOCATE(cc3d)
       ELSE
          ALLOCATE(cr3d(lb1_xpnd:ub1_xpnd, lb2_xpnd:ub2_xpnd, lb3_xpnd:ub3_xpnd))
          cr3d = RESHAPE(pw_in%cr,(/ub1_xpnd-lb1_xpnd+1, ub2_xpnd-lb2_xpnd+1, ub3_xpnd-lb3_xpnd+1/))
       END IF
    ELSE
       IF (in_space .EQ. RECIPROCALSPACE) THEN
          CALL copy_cr(pw_in%cc3d,cr3d)
       ELSE
          ALLOCATE(cr3d(lb1_xpnd:ub1_xpnd, lb2_xpnd:ub2_xpnd, lb3_xpnd:ub3_xpnd))
          cr3d = pw_in%cr3d
       END IF
    END IF

! let all the nodes know about each others shifted local bounds
    ALLOCATE(bounds_local_all(2,3,group_size))
    CALL mp_allgather(bounds_local_shftd, bounds_local_all, rs_group) 

    DO i = 1, 4
! no need to send to myself or to an invalid destination (pid = -1)
       IF ((dests_shrink(i) .NE. rs_mpo) .AND. (dests_shrink(i) .NE. -1)) THEN
          send_lb1 = bounds_local_all(1,1,dests_shrink(i)+1)
          send_ub1 = bounds_local_all(2,1,dests_shrink(i)+1)
          send_lb2 = bounds_local_all(1,2,dests_shrink(i)+1)
          send_ub2 = bounds_local_all(2,2,dests_shrink(i)+1)
          send_lb3 = bounds_local_all(1,3,dests_shrink(i)+1)
          send_ub3 = bounds_local_all(2,3,dests_shrink(i)+1)

          ALLOCATE(send_crmsg(send_lb1:send_ub1, send_lb2:send_ub2, send_lb3:send_ub3))
          send_crmsg = cr3d(send_lb1:send_ub1, send_lb2:send_ub2, send_lb3:send_ub3)
          CALL mp_isend(send_crmsg, dests_shrink(i), rs_group, send_req)
          CALL mp_wait(send_req)
          DEALLOCATE(send_crmsg)
       END IF
    END DO

    lb1_orig = bounds_local_shftd(1,1)
    ub1_orig = bounds_local_shftd(2,1)
    lb2_orig = bounds_local_shftd(1,2)
    ub2_orig = bounds_local_shftd(2,2)
    lb3_orig = bounds_local_shftd(1,3)
    ub3_orig = bounds_local_shftd(2,3)

! no need to receive from myself
    IF (srcs_shrink .EQ. rs_mpo) THEN
       pw_shrinked%cr3d = cr3d(lb1_orig:ub1_orig, lb2_orig:ub2_orig, lb3_orig:ub3_orig)
    ELSE IF (srcs_shrink .EQ. -1) THEN
! the source is invalid ... do nothing
    ELSE
       CALL mp_irecv(pw_shrinked%cr3d, srcs_shrink, rs_group, recv_req) 
       CALL mp_wait(recv_req)
    END IF

    DEALLOCATE(bounds_local_all)
    DEALLOCATE(cr3d)
    CALL timestop(handle)

  END SUBROUTINE pw_shrink

! *****************************************************************************
!> \brief   flips a 3d (real dp) array up to down (the way needed to expand data 
!>          as explained in the description of the afore-defined subroutines)
!> \param cr3d_in input array
!> \param cr3d_out output array
!> \param bounds the bounds of the output array
!> \param error cp2k error
!> 
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE flipud(cr3d_in, cr3d_out, bounds, error)

    REAL(dp), DIMENSION(:, :, :), &
      INTENT(IN), POINTER                    :: cr3d_in
    REAL(dp), DIMENSION(:, :, :), &
      INTENT(OUT), POINTER                   :: cr3d_out
    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bounds
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'flipud', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, lb1, lb1_new, lb2, &
                                                lb2_new, lb3, lb3_new, ub1, &
                                                ub1_new, ub2, ub2_new, ub3, &
                                                ub3_new
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: indx
    INTEGER, DIMENSION(2, 3)                 :: bndsl, bndsl_new

    CALL timeset(routineN,handle)

    lb1 = LBOUND(cr3d_in,1);     ub1 = UBOUND(cr3d_in,1)
    lb2 = LBOUND(cr3d_in,2);     ub2 = UBOUND(cr3d_in,2)
    lb3 = LBOUND(cr3d_in,3);     ub3 = UBOUND(cr3d_in,3)

    bndsl = RESHAPE((/lb1, ub1, lb2, ub2, lb3, ub3/), (/2, 3/))
    bndsl_new = flipud_bounds_local(bndsl, bounds)

    lb1_new = bndsl_new(1,1);    ub1_new = bndsl_new(2,1)
    lb2_new = bndsl_new(1,2);    ub2_new = bndsl_new(2,2)
    lb3_new = bndsl_new(1,3);    ub3_new = bndsl_new(2,3)

    ALLOCATE(indx(ub1_new-lb1_new+1))
    ALLOCATE(cr3d_out(lb1_new : ub1_new, lb2_new : ub2_new, lb3_new : ub3_new))
    cr3d_out = 0.0_dp

    indx(:) = (/ (i, i=2*bounds(2,1) - lb1_new, 2*bounds(2,1) - ub1_new, -1) /)

    cr3d_out(lb1_new:ub1_new, lb2_new:ub2_new, lb3_new:ub3_new) = cr3d_in(indx, :, :)
  
    DEALLOCATE(indx)

    CALL timestop(handle)

  END SUBROUTINE flipud

! *****************************************************************************
!> \brief   flips a 3d (real dp) array left to right (the way needed to expand data 
!>          as explained in the description of the afore-defined subroutines)
!> \param cr3d_in input array
!> \param cr3d_out output array
!> \param bounds the bounds of the output array
!> \param error cp2k error
!> 
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE fliplr(cr3d_in, cr3d_out, bounds, error)

    REAL(dp), DIMENSION(:, :, :), &
      INTENT(IN), POINTER                    :: cr3d_in
    REAL(dp), DIMENSION(:, :, :), &
      INTENT(OUT), POINTER                   :: cr3d_out
    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bounds
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fliplr', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, lb1, lb1_new, lb2, &
                                                lb2_new, lb3, lb3_new, ub1, &
                                                ub1_new, ub2, ub2_new, ub3, &
                                                ub3_new
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: indy
    INTEGER, DIMENSION(2, 3)                 :: bndsl, bndsl_new

    CALL timeset(routineN,handle)

    lb1 = LBOUND(cr3d_in,1);     ub1 = UBOUND(cr3d_in,1)
    lb2 = LBOUND(cr3d_in,2);     ub2 = UBOUND(cr3d_in,2)
    lb3 = LBOUND(cr3d_in,3);     ub3 = UBOUND(cr3d_in,3)

    bndsl = RESHAPE((/lb1, ub1, lb2, ub2, lb3, ub3/), (/2, 3/))
    bndsl_new = fliplr_bounds_local(bndsl, bounds)

    lb1_new = bndsl_new(1,1);    ub1_new = bndsl_new(2,1)
    lb2_new = bndsl_new(1,2);    ub2_new = bndsl_new(2,2)
    lb3_new = bndsl_new(1,3);    ub3_new = bndsl_new(2,3)

    ALLOCATE(indy(ub2_new-lb2_new+1))
    ALLOCATE(cr3d_out(lb1_new : ub1_new, lb2_new : ub2_new, lb3_new : ub3_new))
    cr3d_out = 0.0_dp

    indy(:) = (/ (i, i=2*bounds(2,2) - lb2_new, 2*bounds(2,2) - ub2_new, -1) /)
    cr3d_out(lb1_new:ub1_new, lb2_new:ub2_new, lb3_new:ub3_new) = cr3d_in(:, indy, :)

    DEALLOCATE(indy)

    CALL timestop(handle)

  END SUBROUTINE fliplr

! *****************************************************************************
!> \brief   flips a 3d (real dp) array back to front (the way needed to expand data 
!>          as explained in the description of the afore-defined subroutines)
!> \param cr3d_in input array
!> \param cr3d_out output array
!> \param bounds the bounds of the output array
!> \param error cp2k error
!> 
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE flipbf(cr3d_in, cr3d_out, bounds, error)

    REAL(dp), DIMENSION(:, :, :), &
      INTENT(IN), POINTER                    :: cr3d_in
    REAL(dp), DIMENSION(:, :, :), &
      INTENT(OUT), POINTER                   :: cr3d_out
    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bounds
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'flipbf', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, lb1, lb1_new, lb2, &
                                                lb2_new, lb3, lb3_new, ub1, &
                                                ub1_new, ub2, ub2_new, ub3, &
                                                ub3_new
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: indz
    INTEGER, DIMENSION(2, 3)                 :: bndsl, bndsl_new

    CALL timeset(routineN,handle)

    lb1 = LBOUND(cr3d_in,1);     ub1 = UBOUND(cr3d_in,1)
    lb2 = LBOUND(cr3d_in,2);     ub2 = UBOUND(cr3d_in,2)
    lb3 = LBOUND(cr3d_in,3);     ub3 = UBOUND(cr3d_in,3)

    bndsl = RESHAPE((/lb1, ub1, lb2, ub2, lb3, ub3/), (/2, 3/))
    bndsl_new = flipbf_bounds_local(bndsl, bounds)

    lb1_new = bndsl_new(1,1);    ub1_new = bndsl_new(2,1)
    lb2_new = bndsl_new(1,2);    ub2_new = bndsl_new(2,2)
    lb3_new = bndsl_new(1,3);    ub3_new = bndsl_new(2,3)

    ALLOCATE(indz(ub3_new-lb3_new+1))
    ALLOCATE(cr3d_out(lb1_new : ub1_new, lb2_new : ub2_new, lb3_new : ub3_new))
    cr3d_out = 0.0_dp

    indz(:) = (/ (i, i=2*bounds(2,3) - lb3_new, 2*bounds(2,3) - ub3_new, -1) /)
    cr3d_out(lb1_new:ub1_new, lb2_new:ub2_new, lb3_new:ub3_new) = cr3d_in(:, :, indz)

    DEALLOCATE(indz)

    CALL timestop(handle)

  END SUBROUTINE flipbf

! *****************************************************************************
!> \brief   rotates a 3d (real dp) array by 180 degrees (the way needed to expand data 
!>          as explained in the description of the afore-defined subroutines)
!> \param cr3d_in input array
!> \param cr3d_out output array
!> \param bounds the bounds of the output array
!> \param error cp2k error
!> 
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE rot180(cr3d_in, cr3d_out, bounds, error)

    REAL(dp), DIMENSION(:, :, :), &
      INTENT(IN), POINTER                    :: cr3d_in
    REAL(dp), DIMENSION(:, :, :), &
      INTENT(OUT), POINTER                   :: cr3d_out
    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bounds
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rot180', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, lb1, lb1_new, lb2, &
                                                lb2_new, lb3, lb3_new, ub1, &
                                                ub1_new, ub2, ub2_new, ub3, &
                                                ub3_new
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: indx, indy
    INTEGER, DIMENSION(2, 3)                 :: bndsl, bndsl_new

    CALL timeset(routineN,handle)

    lb1 = LBOUND(cr3d_in,1);     ub1 = UBOUND(cr3d_in,1)
    lb2 = LBOUND(cr3d_in,2);     ub2 = UBOUND(cr3d_in,2)
    lb3 = LBOUND(cr3d_in,3);     ub3 = UBOUND(cr3d_in,3)

    bndsl = RESHAPE((/lb1, ub1, lb2, ub2, lb3, ub3/), (/2, 3/))
    bndsl_new = rot180_bounds_local(bndsl, bounds)

    lb1_new = bndsl_new(1,1);    ub1_new = bndsl_new(2,1)
    lb2_new = bndsl_new(1,2);    ub2_new = bndsl_new(2,2)
    lb3_new = bndsl_new(1,3);    ub3_new = bndsl_new(2,3)

    ALLOCATE(indx(ub1_new-lb1_new+1), indy(ub2_new-lb2_new+1))
    ALLOCATE(cr3d_out(lb1_new : ub1_new, lb2_new : ub2_new, lb3_new : ub3_new))
    cr3d_out = 0.0_dp

    indx(:) = (/ (i, i=2*bounds(2,1) - lb1_new, 2*bounds(2,1) - ub1_new, -1) /)
    indy(:) = (/ (i, i=2*bounds(2,2) - lb2_new, 2*bounds(2,2) - ub2_new, -1) /)

    cr3d_out(lb1_new:ub1_new, lb2_new:ub2_new, lb3_new:ub3_new) = cr3d_in(indx, indy, :)

    DEALLOCATE(indx, indy)

    CALL timestop(handle)

  END SUBROUTINE rot180

! *****************************************************************************
!> \brief   calculates the global and local bounds of the expanded data 
!> \param pw_grid original plane-wave grid
!> \param srcs_expand list of the source processes (pw_expand) 
!> \param flipg_stat flipping status for the received data chunks (pw_expand) 
!> \param bounds_shftd bounds of the original grid shifted to have g0 in the middle of the cell
!> \param bounds_local_shftd local bounds of the original grid after shifting
!> \param recv_msgs_bnds bounds of the messages to be received (pw_expand) 
!> \param bounds_new new global lower and upper bounds
!> \param bounds_local_new new local lower and upper bounds
!> \param error cp2k error
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE expansion_bounds(pw_grid, srcs_expand, flipg_stat, &
                              bounds_shftd, bounds_local_shftd, &
                              recv_msgs_bnds, bounds_new, bounds_local_new, error)

    TYPE(pw_grid_type), INTENT(IN), POINTER  :: pw_grid
    INTEGER, DIMENSION(4), INTENT(IN)        :: srcs_expand, flipg_stat
    INTEGER, DIMENSION(2, 3), INTENT(OUT)    :: bounds_shftd, &
                                                bounds_local_shftd
    INTEGER, DIMENSION(2, 3, 4), INTENT(OUT) :: recv_msgs_bnds
    INTEGER, DIMENSION(2, 3), INTENT(OUT)    :: bounds_new, bounds_local_new
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'expansion_bounds', &
      routineP = moduleN//':'//routineN

    INTEGER :: group_size, handle, i, lb1_new, lb2_new, lb3_new, loc, &
      rs_group, rs_mpo, ub1_new, ub2_new, ub3_new
    INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: bounds_local_all, &
                                                bounds_local_new_all
    INTEGER, DIMENSION(2, 3)                 :: bounds, bounds_local
    INTEGER, DIMENSION(2, 3, 4)              :: pcs_bnds
    INTEGER, DIMENSION(3)                    :: npts_new, shift
    INTEGER, DIMENSION(4)                    :: src_hist

    CALL timeset(routineN,handle)

    rs_group = pw_grid%para%rs_group
    rs_mpo = pw_grid%para%my_pos
    group_size = pw_grid%para%group_size
    bounds = pw_grid%bounds
    bounds_local = pw_grid%bounds_local

    ! Note that this is not easily FFT-able ... needed anyway, so link in FFTW.
    npts_new = 2 * pw_grid%npts - 2
    shift = - npts_new / 2
    shift = shift - bounds(1,:)
    bounds_shftd(:,1) = bounds(:,1) + shift(1) 
    bounds_shftd(:,2) = bounds(:,2) + shift(2) 
    bounds_shftd(:,3) = bounds(:,3) + shift(3)
    bounds_local_shftd(:,1) = bounds_local(:,1) + shift(1) 
    bounds_local_shftd(:,2) = bounds_local(:,2) + shift(2) 
    bounds_local_shftd(:,3) = bounds_local(:,3) + shift(3)

! let all the nodes know about each others local shifted bounds
    ALLOCATE(bounds_local_all(2,3,group_size))
    CALL mp_allgather(bounds_local_shftd, bounds_local_all, rs_group) 

    src_hist = -1   ! keeps the history of sources 

    DO i = 1, 4
! no need to receive from myself
       IF (srcs_expand(i) .EQ. rs_mpo) THEN
          recv_msgs_bnds(1,1,i) = bounds_local_shftd(1,1)
          recv_msgs_bnds(2,1,i) = bounds_local_shftd(2,1)
          recv_msgs_bnds(1,2,i) = bounds_local_shftd(1,2)
          recv_msgs_bnds(2,2,i) = bounds_local_shftd(2,2)
          recv_msgs_bnds(1,3,i) = bounds_local_shftd(1,3)
          recv_msgs_bnds(2,3,i) = bounds_local_shftd(2,3)
! if I have already received data from the source, just use the one from the last time
       ELSE IF (ANY(src_hist .EQ. srcs_expand(i))) THEN 
          loc = MINLOC(ABS(src_hist - srcs_expand(i)), 1)
          recv_msgs_bnds(1,1,i) = bounds_local_all(1,1,srcs_expand(loc)+1)
          recv_msgs_bnds(2,1,i) = bounds_local_all(2,1,srcs_expand(loc)+1)
          recv_msgs_bnds(1,2,i) = bounds_local_all(1,2,srcs_expand(loc)+1)
          recv_msgs_bnds(2,2,i) = bounds_local_all(2,2,srcs_expand(loc)+1)
          recv_msgs_bnds(1,3,i) = bounds_local_all(1,3,srcs_expand(loc)+1)
          recv_msgs_bnds(2,3,i) = bounds_local_all(2,3,srcs_expand(loc)+1)
       ELSE
          recv_msgs_bnds(1,1,i) = bounds_local_all(1,1,srcs_expand(i)+1)
          recv_msgs_bnds(2,1,i) = bounds_local_all(2,1,srcs_expand(i)+1)
          recv_msgs_bnds(1,2,i) = bounds_local_all(1,2,srcs_expand(i)+1)
          recv_msgs_bnds(2,2,i) = bounds_local_all(2,2,srcs_expand(i)+1)
          recv_msgs_bnds(1,3,i) = bounds_local_all(1,3,srcs_expand(i)+1)
          recv_msgs_bnds(2,3,i) = bounds_local_all(2,3,srcs_expand(i)+1)
       END IF
       src_hist(i) = srcs_expand(i)
    END DO

! flip the received data based on the flipping status
    DO i = 1, 4
       SELECT CASE (flipg_stat(i))
          CASE(NOT_FLIPPED)
              pcs_bnds(:,:,i) = recv_msgs_bnds(:,:,i)
          CASE(UD_FLIPPED)
              pcs_bnds(:,:,i) = flipud_bounds_local(recv_msgs_bnds(:,:,i), bounds_shftd)
          CASE(LR_FLIPPED)
              pcs_bnds(:,:,i) = fliplr_bounds_local(recv_msgs_bnds(:,:,i), bounds_shftd)
          CASE(BF_FLIPPED)
              pcs_bnds(:,:,i) = flipbf_bounds_local(recv_msgs_bnds(:,:,i), bounds_shftd)
          CASE(ROTATED)
              pcs_bnds(:,:,i) = rot180_bounds_local(recv_msgs_bnds(:,:,i), bounds_shftd)
       END SELECT
    END DO

    lb1_new = MINVAL(pcs_bnds(1,1,:));    ub1_new = MAXVAL(pcs_bnds(2,1,:))
    lb2_new = MINVAL(pcs_bnds(1,2,:));    ub2_new = MAXVAL(pcs_bnds(2,2,:))
    lb3_new = MINVAL(pcs_bnds(1,3,:));    ub3_new = MAXVAL(pcs_bnds(2,3,:))

! calculate the new local and global bounds
    bounds_local_new(1,1) = MINVAL(pcs_bnds(1,1,:))    
    bounds_local_new(2,1) = MAXVAL(pcs_bnds(2,1,:))
    bounds_local_new(1,2) = MINVAL(pcs_bnds(1,2,:))    
    bounds_local_new(2,2) = MAXVAL(pcs_bnds(2,2,:))
    bounds_local_new(1,3) = MINVAL(pcs_bnds(1,3,:)) 
    bounds_local_new(2,3) = 2*MAXVAL(pcs_bnds(2,3,:)) - bounds_local_new(1,3) - 1

    ALLOCATE(bounds_local_new_all(2,3,group_size))
    CALL mp_allgather(bounds_local_new, bounds_local_new_all, rs_group)
    bounds_new(1,1) = MINVAL(bounds_local_new_all(1,1,:)) 
    bounds_new(2,1) = MAXVAL(bounds_local_new_all(2,1,:)) 
    bounds_new(1,2) = MINVAL(bounds_local_new_all(1,2,:)) 
    bounds_new(2,2) = MAXVAL(bounds_local_new_all(2,2,:)) 
    bounds_new(1,3) = MINVAL(bounds_local_new_all(1,3,:)) 
    bounds_new(2,3) = MAXVAL(bounds_local_new_all(2,3,:)) 

    DEALLOCATE(bounds_local_all, bounds_local_new_all)

    CALL timestop(handle)

  END SUBROUTINE expansion_bounds

! *****************************************************************************
!> \brief   precalculates the local bounds of a 3d array after applying flipud
!> \param bndsl_in current local lower and upper bounds
!> \param bounds global lower and upper bounds
!> \retval bndsl_out new local lower and upper bounds
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  FUNCTION flipud_bounds_local(bndsl_in, bounds) RESULT(bndsl_out)

    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bndsl_in, bounds
    INTEGER, DIMENSION(2, 3)                 :: bndsl_out

    CHARACTER(LEN=*), PARAMETER :: routineN = 'flipud_bounds_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    bndsl_out(1,1) = 2*bounds(2,1) - bndsl_in(2,1)
    bndsl_out(2,1) = 2*bounds(2,1) - bndsl_in(1,1)
    IF (bndsl_out(1,1) .EQ. bounds(2,1)) bndsl_out(1,1) = bndsl_out(1,1) + 1
    IF (bndsl_out(2,1) .EQ. 2*bounds(2,1) - bounds(1,1)) bndsl_out(2,1) = bndsl_out(2,1) - 1 

    bndsl_out(1,2) = bndsl_in(1,2)    
    bndsl_out(2,2) = bndsl_in(2,2)

    bndsl_out(1,3) = bndsl_in(1,3)    
    bndsl_out(2,3) = bndsl_in(2,3)

    CALL timestop(handle)

  END FUNCTION flipud_bounds_local

! *****************************************************************************
!> \brief   precalculates the local bounds of a 3d array after applying fliplr
!> \param bndsl_in current local lower and upper bounds
!> \param bounds global lower and upper bounds
!> \retval bndsl_out new local lower and upper bounds
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  FUNCTION fliplr_bounds_local(bndsl_in, bounds) RESULT(bndsl_out)

    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bndsl_in, bounds
    INTEGER, DIMENSION(2, 3)                 :: bndsl_out

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fliplr_bounds_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    bndsl_out(1,1) = bndsl_in(1,1)    
    bndsl_out(2,1) = bndsl_in(2,1)

    bndsl_out(1,2) = 2*bounds(2,2) - bndsl_in(2,2)
    bndsl_out(2,2) = 2*bounds(2,2) - bndsl_in(1,2)
    IF (bndsl_out(1,2) .EQ. bounds(2,1)) bndsl_out(1,2) = bndsl_out(1,2) + 1
    IF (bndsl_out(2,2) .EQ. 2*bounds(2,2) - bounds(1,2)) bndsl_out(2,2) = bndsl_out(2,2) - 1 

    bndsl_out(1,3) = bndsl_in(1,3)    
    bndsl_out(2,3) = bndsl_in(2,3)

    CALL timestop(handle)

  END FUNCTION fliplr_bounds_local

! *****************************************************************************
!> \brief   precalculates the local bounds of a 3d array after applying flipbf
!> \param bndsl_in current local lower and upper bounds
!> \param bounds global lower and upper bounds
!> \retval bndsl_out new local lower and upper bounds
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  FUNCTION flipbf_bounds_local(bndsl_in, bounds) RESULT(bndsl_out)

    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bndsl_in, bounds
    INTEGER, DIMENSION(2, 3)                 :: bndsl_out

    CHARACTER(LEN=*), PARAMETER :: routineN = 'flipbf_bounds_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    bndsl_out(1,1) = bndsl_in(1,1)
    bndsl_out(2,1) = bndsl_in(2,1)

    bndsl_out(1,2) = bndsl_in(1,2)
    bndsl_out(2,2) = bndsl_in(2,2)

    bndsl_out(1,3) = 2*bounds(2,3) - bndsl_in(2,3)
    bndsl_out(2,3) = 2*bounds(2,3) - bndsl_in(1,3)
    IF (bndsl_out(1,3) .EQ. bounds(2,3)) bndsl_out(1,3) = bndsl_out(1,3) + 1
    IF (bndsl_out(2,3) .EQ. 2*bounds(2,3) - bounds(1,3)) bndsl_out(2,3) = bndsl_out(2,3) - 1

    CALL timestop(handle)

  END FUNCTION flipbf_bounds_local

! *****************************************************************************
!> \brief   precalculates the local bounds of a 3d array after applying rot180
!> \param bndsl_in current local lower and upper bounds
!> \param bounds global lower and upper bounds
!> \retval bndsl_out new local lower and upper bounds
!> \par History
!>       07.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  FUNCTION rot180_bounds_local(bndsl_in, bounds) RESULT(bndsl_out)

    INTEGER, DIMENSION(2, 3), INTENT(IN)     :: bndsl_in, bounds
    INTEGER, DIMENSION(2, 3)                 :: bndsl_out

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rot180_bounds_local', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    bndsl_out(1,1) = 2*bounds(2,1) - bndsl_in(2,1)
    bndsl_out(2,1) = 2*bounds(2,1) - bndsl_in(1,1)
    IF (bndsl_out(1,1) .EQ. bounds(2,1)) bndsl_out(1,1) = bndsl_out(1,1) + 1
    IF (bndsl_out(2,1) .EQ. 2*bounds(2,1) - bounds(1,1)) bndsl_out(2,1) = bndsl_out(2,1) - 1

    bndsl_out(1,2) = 2*bounds(2,2) - bndsl_in(2,2)
    bndsl_out(2,2) = 2*bounds(2,2) - bndsl_in(1,2)
    IF (bndsl_out(1,2) .EQ. bounds(2,2)) bndsl_out(1,2) = bndsl_out(1,2) + 1
    IF (bndsl_out(2,2) .EQ. 2*bounds(2,2) - bounds(1,2)) bndsl_out(2,2) = bndsl_out(2,2) - 1

    bndsl_out(1,3) = bndsl_in(1,3)
    bndsl_out(2,3) = bndsl_in(2,3)

    CALL timestop(handle)

  END FUNCTION rot180_bounds_local

! *****************************************************************************
!> \brief   Copies pw_old to pw_new that has different local bounds. 
!> \param pw_old input pw data
!> \param pw_new output pw data
!> \param error cp2k error
!> \par History
!>       12.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE pw_adjust_bounds_local(pw_old, pw_new, error)

    TYPE(pw_type), INTENT(IN), POINTER       :: pw_old
    TYPE(pw_type), INTENT(INOUT), POINTER    :: pw_new
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'pw_adjust_bounds_local', &
      routineP = moduleN//':'//routineN

    INTEGER :: enbr_pid, handle, I1, I2, J1, J2, lb1_new, lb1_old, lb1max, &
      lb1min, lb2_new, lb2_old, lb2max, lb2min, nnbr_pid, rs_group, snbr_pid, &
      ub1_new, ub1_old, ub1max, ub1min, ub2_new, ub2_old, ub2max, ub2min, &
      wnbr_pid
    INTEGER, DIMENSION(2)                    :: color, dim1swp_stat, &
                                                dim2swp_stat, edim2swp_stat, &
                                                rs_pos, wdim2swp_stat
    INTEGER, DIMENSION(2, 3)                 :: bounds, bounds_local_new, &
                                                bounds_local_old
    REAL(dp), DIMENSION(:, :, :), POINTER    :: e_rmsg, n_rmsg, ne_rmsg, &
                                                nw_rmsg, s_rmsg, se_rmsg, &
                                                sw_rmsg, w_rmsg
    TYPE(pw_grid_type), POINTER              :: pw_grid_new, pw_grid_old

    CALL timeset(routineN,handle)

    pw_grid_old => pw_old%pw_grid
    pw_grid_new => pw_new%pw_grid

    rs_pos = pw_grid_old%para%rs_pos
    rs_group = pw_grid_old%para%rs_group
    bounds_local_old = pw_grid_old%bounds_local
    bounds_local_new = pw_grid_new%bounds_local
    bounds = pw_grid_old%bounds

    lb1_old = bounds_local_old(1,1);    ub1_old = bounds_local_old(2,1)
    lb2_old = bounds_local_old(1,2);    ub2_old = bounds_local_old(2,2)
    lb1_new = bounds_local_new(1,1);    ub1_new = bounds_local_new(2,1)
    lb2_new = bounds_local_new(1,2);    ub2_new = bounds_local_new(2,2)

    lb1min = MIN(lb1_old, lb1_new);     lb1max = MAX(lb1_old, lb1_new)     
    lb2min = MIN(lb2_old, lb2_new);     lb2max = MAX(lb2_old, lb2_new)     
    ub1min = MIN(ub1_old, ub1_new);     ub1max = MAX(ub1_old, ub1_new)     
    ub2min = MIN(ub2_old, ub2_new);     ub2max = MAX(ub2_old, ub2_new)     
    ! copy data within the shared bounds
    pw_new%cr3d(lb1max:ub1min, lb2max:ub2min, :) = pw_old%cr3d(lb1max:ub1min, lb2max:ub2min, :)

    IF (lb1_old .NE. bounds(1,1)) THEN  ! west neighbor exists
       CALL mp_cart_rank(rs_group, (/rs_pos(1)-1, rs_pos(2)/), wnbr_pid) 
    ELSE
       wnbr_pid = -1
    END IF
    IF (ub1_old .NE. bounds(2,1)) THEN  ! east neighbor exists
       CALL mp_cart_rank(rs_group, (/rs_pos(1)+1, rs_pos(2)/), enbr_pid) 
    ELSE
       enbr_pid = -1
    END IF
    IF (lb2_old .NE. bounds(1,2)) THEN  ! south neighbor exists
       CALL mp_cart_rank(rs_group, (/rs_pos(1), rs_pos(2)+1/), snbr_pid) 
    ELSE
       snbr_pid = -1
    END IF
    IF (ub2_old .NE. bounds(2,2)) THEN  ! north neighbor exists
       CALL mp_cart_rank(rs_group, (/rs_pos(1), rs_pos(2)-1/), nnbr_pid) 
    ELSE
       nnbr_pid = -1
    END IF

    ! set the colors for red-black algorithm 
    color = (/MOD(rs_pos(1),2), MOD(rs_pos(2),2)/) 

    ! sweep edge data
    CALL dim1_sweep(pw_old%cr3d, bounds_local_new(:,1), wnbr_pid, enbr_pid, &
                    rs_group, color, e_rmsg, w_rmsg, dim1swp_stat, error)
    CALL dim2_sweep(pw_old%cr3d, bounds_local_new(:,2), nnbr_pid, snbr_pid, &
                    rs_group, color, n_rmsg, s_rmsg, dim2swp_stat, error)
    ! sweep corner data
    IF (dim1swp_stat(1) .NE. nomsg_received) THEN
       CALL dim2_sweep(w_rmsg, bounds_local_new(:,2), nnbr_pid, snbr_pid, &
                       rs_group, color, nw_rmsg, sw_rmsg, wdim2swp_stat, error)
    ELSE
       wdim2swp_stat = nomsg_received
    END IF
    IF (dim1swp_stat(2) .NE. nomsg_received) THEN
       CALL dim2_sweep(e_rmsg, bounds_local_new(:,2), nnbr_pid, snbr_pid, &
                       rs_group, color, ne_rmsg, se_rmsg, edim2swp_stat, error)
    ELSE
       edim2swp_stat = nomsg_received
    END IF

    ! copy edge data
    IF (dim1swp_stat(1) .EQ. msg_received) THEN 
       I1 = lb1_new;       I2 = lb1_old-1
       J1 = lb2_new;       J2 = ub2_new
       pw_new%cr3d(I1:I2, J1:J2, :) = w_rmsg(I1:I2, J1:J2, :)
    END IF
    IF (dim1swp_stat(2) .EQ. msg_received) THEN
       I1 = ub1_old+1;     I2 = ub1_new
       J1 = lb2_new;       J2 = ub2_new
       pw_new%cr3d(I1:I2, J1:J2, :) = e_rmsg(I1:I2, J1:J2, :)
    END IF
    IF (dim2swp_stat(1) .EQ. msg_received) THEN
       I1 = lb1_new;       I2 = ub1_new
       J1 = lb2_new;       J2 = lb2_old-1
       pw_new%cr3d(I1:I2, J1:J2, :) = s_rmsg(I1:I2, J1:J2, :)
    END IF
    IF (dim2swp_stat(2) .EQ. msg_received) THEN
       I1 = lb1_new;       I2 = ub1_new
       J1 = ub2_old+1;     J2 = ub2_new
       pw_new%cr3d(I1:I2, J1:J2, :) = n_rmsg(I1:I2, J1:J2, :)
    END IF

    ! copy corner data
    IF (edim2swp_stat(1) .EQ. msg_received) THEN
       I1 = ub1_old+1;     I2 = ub1_new
       J1 = lb2_new;       J2 = lb2_old-1
       pw_new%cr3d(I1:I2, J1:J2, :) = se_rmsg(I1:I2, J1:J2, :)
    END IF
    IF (edim2swp_stat(2) .EQ. msg_received) THEN 
       I1 = ub1_old+1;     I2 = ub1_new
       J1 = ub2_old+1;     J2 = ub2_new
       pw_new%cr3d(I1:I2, J1:J2, :) = ne_rmsg(I1:I2, J1:J2, :)
    END IF
    IF (wdim2swp_stat(1) .EQ. msg_received) THEN
       I1 = lb1_new;       I2 = lb1_old-1
       J1 = lb2_new;       J2 = lb2_old-1
       pw_new%cr3d(I1:I2, J1:J2, :) = sw_rmsg(I1:I2, J1:J2, :)
    END IF
    IF (wdim2swp_stat(2) .EQ. msg_received) THEN
       I1 = lb1_new;       I2 = lb1_old-1
       J1 = ub2_old+1;     J2 = ub2_new
       pw_new%cr3d(I1:I2, J1:J2, :) = nw_rmsg(I1:I2, J1:J2, :)
    END IF

    ! deallocate received data after using them
    IF (dim1swp_stat(1)  .NE. nomsg_received) DEALLOCATE(w_rmsg)
    IF (dim1swp_stat(2)  .NE. nomsg_received) DEALLOCATE(e_rmsg)
    IF (dim2swp_stat(1)  .NE. nomsg_received) DEALLOCATE(s_rmsg)
    IF (dim2swp_stat(2)  .NE. nomsg_received) DEALLOCATE(n_rmsg)
    IF (wdim2swp_stat(1) .NE. nomsg_received) DEALLOCATE(sw_rmsg)
    IF (wdim2swp_stat(2) .NE. nomsg_received) DEALLOCATE(nw_rmsg)
    IF (edim2swp_stat(1) .NE. nomsg_received) DEALLOCATE(se_rmsg)
    IF (edim2swp_stat(2) .NE. nomsg_received) DEALLOCATE(ne_rmsg)

    CALL timestop(handle)

  END SUBROUTINE pw_adjust_bounds_local

! *****************************************************************************
!> \brief sweeps data along the first dimension, in both directions. Neighboring 
!>   processes exchange boundary data using a 2D red-black (chess-board) algorithm
!> \param arr3d_old input 3d array
!> \param dim1blcl_new new local bounds along the first dimension
!> \param wnbr_pid process id of the west neighbor
!> \param enbr_pid process id of the east neighbor
!> \param mp_group mpi group id
!> \param color process coloring
!> \param e_rmsg received message from the east neighbor
!> \param w_rmsg received message from the west neighbor
!> \param sweep_stat sweeping status
!> \param error cp2k error
!> \par History
!>       12.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dim1_sweep(arr3d_old, dim1blcl_new, wnbr_pid, enbr_pid, &
                        mp_group, color, e_rmsg, w_rmsg, sweep_stat, error)

    REAL(dp), DIMENSION(:, :, :), &
      INTENT(IN), POINTER                    :: arr3d_old
    INTEGER, DIMENSION(2), INTENT(IN)        :: dim1blcl_new
    INTEGER, INTENT(IN)                      :: wnbr_pid, enbr_pid, mp_group
    INTEGER, DIMENSION(2), INTENT(IN)        :: color
    REAL(dp), DIMENSION(:, :, :), &
      INTENT(OUT), POINTER                   :: e_rmsg, w_rmsg
    INTEGER, DIMENSION(2), INTENT(OUT)       :: sweep_stat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dim1_sweep', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ew_rreq, ew_sreq, handle, &
                                                lb1_new, lb1_old, lb2, lb3, &
                                                ub1_new, ub1_old, ub2, ub3, &
                                                we_rreq, we_sreq
    LOGICAL                                  :: enbr_exists, should_recve, &
                                                should_recvw, should_sende, &
                                                should_sendw, wnbr_exists
    REAL(dp), DIMENSION(:, :, :), POINTER    :: e_smsg, w_smsg

    CALL timeset(routineN,handle)

    lb1_old = LBOUND(arr3d_old,1);      ub1_old = UBOUND(arr3d_old,1)
    lb1_new = dim1blcl_new(1);          ub1_new = dim1blcl_new(2)
    lb2 = LBOUND(arr3d_old,2);          ub2 = UBOUND(arr3d_old,2)
    lb3 = LBOUND(arr3d_old,3);          ub3 = UBOUND(arr3d_old,3)

    enbr_exists = enbr_pid .GE. 0;      wnbr_exists = wnbr_pid .GE. 0

    should_sendw = lb1_old .LT. lb1_new;     should_sende = ub1_old .GT. ub1_new
    should_recvw = lb1_new .LT. lb1_old;     should_recve = ub1_new .GT. ub1_old 

    ! initiate send buffers
    IF (should_sendw) THEN
       ALLOCATE(w_smsg(lb1_old:lb1_new-1 , lb2:ub2 , lb3:ub3))
       w_smsg  = arr3d_old(lb1_old:lb1_new-1 , : , :)
    END IF 
    IF (should_sende) THEN
       ALLOCATE(e_smsg(ub1_new+1:ub1_old , lb2:ub2 , lb3:ub3))
       e_smsg  = arr3d_old(ub1_new+1:ub1_old , : , :)
    END IF 

    ! allocate receive buffers
    IF (should_recvw) THEN
       ALLOCATE(w_rmsg(lb1_new:lb1_old-1 , lb2:ub2 , lb3:ub3))
    END IF
    IF (should_recve) THEN
       ALLOCATE(e_rmsg(ub1_old+1:ub1_new , lb2:ub2 , lb3:ub3))
    END IF

    ! WE , first black processes then red processes
    IF ((color(1) .EQ. black) .AND. enbr_exists) THEN 
       ! black sends to the east
       IF (should_sende) THEN
          CALL mp_isend(e_smsg, enbr_pid, mp_group, we_sreq)
          CALL mp_wait(we_sreq)
          DEALLOCATE(e_smsg)
       END IF
       IF (should_recve) THEN
       ! black receives from the east
          CALL mp_irecv(e_rmsg, enbr_pid, mp_group, ew_rreq)
          CALL mp_wait(ew_rreq)
       END IF
    ELSE IF ((color(1) .EQ. red) .AND. wnbr_exists) THEN
       ! red receives from the west
       IF (should_recvw) THEN
          CALL mp_irecv(w_rmsg, wnbr_pid, mp_group, we_rreq)
          CALL mp_wait(we_rreq)
       END IF
       ! red sends to the west
       IF (should_sendw) THEN
          CALL mp_isend(w_smsg, wnbr_pid, mp_group, ew_sreq)
          CALL mp_wait(ew_sreq)
          DEALLOCATE(w_smsg)
       END IF
    END IF

    ! EW
    IF ((color(1) .EQ. red) .AND. enbr_exists) THEN 
       ! red sends to the east
       IF (should_sende) THEN
          CALL mp_isend(e_smsg, enbr_pid, mp_group, we_sreq)
          CALL mp_wait(we_sreq)
          DEALLOCATE(e_smsg)
       END IF
       IF (should_recve) THEN
       ! red receives from the east
          CALL mp_irecv(e_rmsg, enbr_pid, mp_group, ew_rreq)
          CALL mp_wait(ew_rreq)
       END IF
    ELSE IF ((color(1) .EQ. black) .AND. wnbr_exists) THEN
       ! black receives from the west
       IF (should_recvw) THEN
          CALL mp_irecv(w_rmsg, wnbr_pid, mp_group, we_rreq)
          CALL mp_wait(we_rreq)
       END IF
       ! black sends to the west
       IF (should_sendw) THEN
          CALL mp_isend(w_smsg, wnbr_pid, mp_group, ew_sreq)
          CALL mp_wait(ew_sreq)
          DEALLOCATE(w_smsg)
       END IF
    END IF
    
    IF (.NOT. should_recvw) THEN
       sweep_stat(1) = nomsg_received
    ELSE IF (should_recvw .AND. (SIZE(w_rmsg) .EQ. 0)) THEN
       sweep_stat(1) = zsmsg_received
    ELSE IF (should_recvw .AND. (SIZE(w_rmsg) .GT. 0)) THEN
       sweep_stat(1) = msg_received
    END IF

    IF (.NOT. should_recve) THEN
       sweep_stat(2) = nomsg_received
    ELSE IF (should_recve .AND. (SIZE(e_rmsg) .EQ. 0)) THEN
       sweep_stat(2) = zsmsg_received
    ELSE IF (should_recve .AND. (SIZE(e_rmsg) .GT. 0)) THEN
       sweep_stat(2) = msg_received
    END IF

    CALL timestop(handle)

  END SUBROUTINE dim1_sweep

! *****************************************************************************
!> \brief sweeps data along the second dimension, in both directions. Neighboring 
!>   processes exchange boundary data using a 2D red-black (chess-board) algorithm
!> \param arr3d_old input 3d array
!> \param dim2blcl_new new local bound along the second dimension
!> \param nnbr_pid process id of the north neighbor
!> \param snbr_pid process id of the south neighbor
!> \param mp_group mpi group id
!> \param color process coloring
!> \param n_rmsg received message from the north neighbor
!> \param s_rmsg received message from the south neighbor
!> \param sweep_stat data sweeping status
!> \param error cp2k error
!> \par History
!>       12.2014 created [Hossein Bani-Hashemian]
!> \author Mohammad Hossein Bani-Hashemian
! *****************************************************************************
  SUBROUTINE dim2_sweep(arr3d_old, dim2blcl_new, nnbr_pid, snbr_pid, &
                        mp_group, color, n_rmsg, s_rmsg, sweep_stat, error)

    REAL(dp), DIMENSION(:, :, :), &
      INTENT(IN), POINTER                    :: arr3d_old
    INTEGER, DIMENSION(2), INTENT(IN)        :: dim2blcl_new
    INTEGER, INTENT(IN)                      :: nnbr_pid, snbr_pid, mp_group
    INTEGER, DIMENSION(2), INTENT(IN)        :: color
    REAL(dp), DIMENSION(:, :, :), &
      INTENT(OUT), POINTER                   :: n_rmsg, s_rmsg
    INTEGER, DIMENSION(2), INTENT(OUT)       :: sweep_stat
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dim2_sweep', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: black = 0, red = 1 

    INTEGER                                  :: handle, lb1, lb2_new, &
                                                lb2_old, lb3, ns_rreq, &
                                                ns_sreq, sn_rreq, sn_sreq, &
                                                ub1, ub2_new, ub2_old, ub3
    LOGICAL                                  :: nnbr_exists, should_recvn, &
                                                should_recvs, should_sendn, &
                                                should_sends, snbr_exists
    REAL(dp), DIMENSION(:, :, :), POINTER    :: n_smsg, s_smsg

    CALL timeset(routineN,handle)

    lb1 = LBOUND(arr3d_old,1);          ub1 = UBOUND(arr3d_old,1)
    lb2_old = LBOUND(arr3d_old,2);      ub2_old = UBOUND(arr3d_old,2)
    lb2_new = dim2blcl_new(1);          ub2_new = dim2blcl_new(2)
    lb3 = LBOUND(arr3d_old,3);          ub3 = UBOUND(arr3d_old,3)

    nnbr_exists = nnbr_pid .GE. 0;      snbr_exists = snbr_pid .GE. 0

    should_sends = lb2_old .LT. lb2_new;     should_sendn = ub2_old .GT. ub2_new
    should_recvs = lb2_new .LT. lb2_old;     should_recvn = ub2_new .GT. ub2_old 

    ! initiate send buffers
    IF (should_sends) THEN
       ALLOCATE(s_smsg(lb1:ub1 , lb2_old:lb2_new-1 , lb3:ub3))
       s_smsg  = arr3d_old(: , lb2_old:lb2_new-1 , :)
    END IF 
    IF (should_sendn) THEN
       ALLOCATE(n_smsg(lb1:ub1 , ub2_new+1:ub2_old , lb3:ub3))
       n_smsg  = arr3d_old(: , ub2_new+1:ub2_old , :)
    END IF 

    ! allocate receive buffers
    IF (should_recvs) THEN
       ALLOCATE(s_rmsg(lb1:ub1 , lb2_new:lb2_old-1 , lb3:ub3))
    END IF
    IF (should_recvn) THEN
       ALLOCATE(n_rmsg(lb1:ub1 , ub2_old+1:ub2_new , lb3:ub3))
    END IF

    ! SN , first black processes then red processes
    IF ((color(2) .EQ. black) .AND. nnbr_exists) THEN 
       ! black sends to the north
       IF (should_sendn) THEN
          CALL mp_isend(n_smsg, nnbr_pid, mp_group, sn_sreq)
          CALL mp_wait(sn_sreq)
          DEALLOCATE(n_smsg)
       END IF
       IF (should_recvn) THEN
       ! black receives from the north
          CALL mp_irecv(n_rmsg, nnbr_pid, mp_group, ns_rreq)
          CALL mp_wait(ns_rreq)
       END IF
    ELSE IF ((color(2) .EQ. red) .AND. snbr_exists) THEN
       ! red receives from the south
       IF (should_recvs) THEN
          CALL mp_irecv(s_rmsg, snbr_pid, mp_group, sn_rreq)
          CALL mp_wait(sn_rreq)
       END IF
       ! red sends to the south
       IF (should_sends) THEN
          CALL mp_isend(s_smsg, snbr_pid, mp_group, ns_sreq)
          CALL mp_wait(ns_sreq)
          DEALLOCATE(s_smsg)
       END IF
    END IF

    ! NS
    IF ((color(2) .EQ. red) .AND. nnbr_exists) THEN 
       ! red sends to the north
       IF (should_sendn) THEN
          CALL mp_isend(n_smsg, nnbr_pid, mp_group, sn_sreq)
          CALL mp_wait(sn_sreq)
          DEALLOCATE(n_smsg)
       END IF
       IF (should_recvn) THEN
       ! red receives from the north
          CALL mp_irecv(n_rmsg, nnbr_pid, mp_group, ns_rreq)
          CALL mp_wait(ns_rreq)
       END IF
    ELSE IF ((color(2) .EQ. black) .AND. snbr_exists) THEN
       ! black receives from the south
       IF (should_recvs) THEN
          CALL mp_irecv(s_rmsg, snbr_pid, mp_group, sn_rreq)
          CALL mp_wait(sn_rreq)
       END IF
       ! black sends to the south
       IF (should_sends) THEN
          CALL mp_isend(s_smsg, snbr_pid, mp_group, ns_sreq)
          CALL mp_wait(ns_sreq)
          DEALLOCATE(s_smsg)
       END IF
    END IF

    IF (.NOT. should_recvs) THEN
       sweep_stat(1) = nomsg_received
    ELSE IF (should_recvs .AND. (SIZE(s_rmsg) .EQ. 0)) THEN
       sweep_stat(1) = zsmsg_received
    ELSE IF (should_recvs .AND. (SIZE(s_rmsg) .GT. 0)) THEN
       sweep_stat(1) = msg_received
    END IF

    IF (.NOT. should_recvn) THEN
       sweep_stat(2) = nomsg_received
    ELSE IF (should_recvn .AND. (SIZE(n_rmsg) .EQ. 0)) THEN
       sweep_stat(2) = zsmsg_received
    ELSE IF (should_recvn .AND. (SIZE(n_rmsg) .GT. 0)) THEN
       sweep_stat(2) = msg_received
    END IF

    CALL timestop(handle)

  END SUBROUTINE dim2_sweep

END MODULE dct

