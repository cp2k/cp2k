!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

MODULE qs_fb_atomic_matrix_methods

  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_col_block_sizes,&
                                             cp_dbcsr_get_block_p,&
                                             cp_dbcsr_get_info,&
                                             cp_dbcsr_get_stored_coordinates,&
                                             cp_dbcsr_row_block_sizes,&
                                             cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE kinds,                           ONLY: dp,&
                                             int_8
  USE message_passing,                 ONLY: mp_alltoall
  USE qs_fb_atomic_halo_types,         ONLY: fb_atomic_halo_atom_global2halo,&
                                             fb_atomic_halo_get,&
                                             fb_atomic_halo_has_data,&
                                             fb_atomic_halo_list_get,&
                                             fb_atomic_halo_list_obj,&
                                             fb_atomic_halo_obj
  USE qs_fb_com_tasks_types,           ONLY: &
       TASK_COST, TASK_DEST, TASK_N_RECORDS, TASK_PAIR, TASK_SRC, &
       fb_com_atom_pairs_calc_buffer_sizes, fb_com_atom_pairs_create, &
       fb_com_atom_pairs_decode, fb_com_atom_pairs_get, &
       fb_com_atom_pairs_has_data, fb_com_atom_pairs_init, &
       fb_com_atom_pairs_nullify, fb_com_atom_pairs_obj, &
       fb_com_atom_pairs_release, fb_com_tasks_build_atom_pairs, &
       fb_com_tasks_create, fb_com_tasks_decode_pair, &
       fb_com_tasks_encode_pair, fb_com_tasks_get, fb_com_tasks_nullify, &
       fb_com_tasks_obj, fb_com_tasks_release, fb_com_tasks_set, &
       fb_com_tasks_transpose_dest_src
  USE qs_fb_matrix_data_types,         ONLY: fb_matrix_data_get,&
                                             fb_matrix_data_has_data,&
                                             fb_matrix_data_obj
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fb_atomic_matrix_methods'

  PUBLIC :: fb_atmatrix_calc_size,&
            fb_atmatrix_construct,&
            fb_atmatrix_construct_2,&
            fb_atmatrix_generate_com_pairs_2


CONTAINS


! **********************************************************************
!> \brief Calculates the atomic matrix size from a given DBCSR matrix
!>        and atomic halo. It also calculates the first row (col) or the
!>        row (col) atomic blocks in the atomic matrix
!> \param dbcsr_mat : pointer to the DBCSR matrix the atomic matrix is
!>                    to be constructed from
!> \param atomic_halo : the atomic halo used for defining the atomic
!>                      matrix from the DBCSR matrix
!> \param nrows : outputs total number of rows in the atomic matrix
!> \param ncols : outputs total number of cols in the atomic matrix
!> \param blk_row_start : first row in each atomic blk row in the
!>                        atomic matrix
!> \param blk_col_start : first col in each atomic blk col in the
!>                        atomic matrix
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atmatrix_calc_size(dbcsr_mat, &
                                   atomic_halo, &
                                   nrows, &
                                   ncols, &
                                   blk_row_start, &
                                   blk_col_start)
    TYPE(cp_dbcsr_type), POINTER             :: dbcsr_mat
    TYPE(fb_atomic_halo_obj), INTENT(IN)     :: atomic_halo
    INTEGER, INTENT(OUT)                     :: nrows, ncols
    INTEGER, DIMENSION(:), INTENT(OUT)       :: blk_row_start, blk_col_start

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_atmatrix_calc_size', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ii, natoms_in_halo
    INTEGER, DIMENSION(:), POINTER           :: col_block_size_data, &
                                                halo_atoms, &
                                                row_block_size_data
    LOGICAL                                  :: check_ok, failure

    failure = .FALSE.

    NULLIFY(halo_atoms, row_block_size_data, col_block_size_data)

    row_block_size_data => cp_dbcsr_row_block_sizes(dbcsr_mat)
    col_block_size_data => cp_dbcsr_col_block_sizes(dbcsr_mat)
    CALL fb_atomic_halo_get(atomic_halo=atomic_halo, &
                            natoms=natoms_in_halo, &
                            halo_atoms=halo_atoms)
    check_ok = SIZE(blk_row_start) .GE. (natoms_in_halo + 1)
    CPASSERT(check_ok)
    check_ok = SIZE(blk_col_start) .GE. (natoms_in_halo + 1)
    CPASSERT(check_ok)
    blk_row_start = 0
    blk_col_start = 0
    nrows = 0
    ncols = 0
    DO ii = 1, natoms_in_halo
       blk_row_start(ii) = nrows + 1
       blk_col_start(ii) = ncols + 1
       nrows = nrows + row_block_size_data(halo_atoms(ii))
       ncols = ncols + col_block_size_data(halo_atoms(ii))
    END DO
    blk_row_start(natoms_in_halo + 1) = nrows + 1
    blk_col_start(natoms_in_halo + 1) = ncols + 1
  END SUBROUTINE fb_atmatrix_calc_size


! ****************************************************************************
!> \brief Constructs atomic matrix for filter basis method from a given
!>        DBCSR matrix and a set of atomic send and recv pairs
!>        corresponding to the matrix blocks that needs to be included
!>        in the atomic matrix. This version is for when we do MPI
!>        communications at every step, for each atomic matrix.
!> \param dbcsr_mat : the DBCSR matrix the atomic matrix is to be
!>                    constructed from
!> \param atomic_halo : the atomic halo conrresponding to this atomic
!>                      matrix
!> \param para_env : cp2k parallel environment
!> \param atomic_matrix : the atomic matrix to be constructed, it should
!>                        have already been allocated prior entering
!>                        this subroutine
!> \param blk_row_start : first row in each atomic blk row in the
!>                        atomic matrix
!> \param blk_col_start : first col in each atomic blk col in the
!>                        atomic matrix
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atmatrix_construct(dbcsr_mat, &
                                   atomic_halo, &
                                   para_env, &
                                   atomic_matrix, &
                                   blk_row_start, &
                                   blk_col_start)
    TYPE(cp_dbcsr_type), POINTER             :: dbcsr_mat
    TYPE(fb_atomic_halo_obj), INTENT(IN)     :: atomic_halo
    TYPE(cp_para_env_type), POINTER          :: para_env
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: atomic_matrix
    INTEGER, DIMENSION(:), INTENT(IN)        :: blk_row_start, blk_col_start

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_atmatrix_construct', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iatom, iatom_in_halo, ii, ind, ipair, ipe, jatom, &
      jatom_in_halo, jj, ncols_blk, npairs_recv, npairs_send, nrows_blk, &
      numprocs, pe, recv_encode, send_encode, stat
    INTEGER(KIND=int_8), DIMENSION(:), &
      POINTER                                :: pairs_recv, pairs_send
    INTEGER, ALLOCATABLE, DIMENSION(:) :: recv_disps, recv_pair_count, &
      recv_pair_disps, recv_sizes, send_disps, send_pair_count, &
      send_pair_disps, send_sizes
    INTEGER, DIMENSION(:), POINTER           :: col_block_size_data, &
                                                row_block_size_data
    LOGICAL                                  :: failure, found
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: recv_buf, send_buf
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: mat_block
    TYPE(fb_com_atom_pairs_obj)              :: atom_pairs_recv, &
                                                atom_pairs_send

    CALL timeset(routineN, handle)

    failure = .FALSE.

    NULLIFY(pairs_send, pairs_recv, mat_block, &
            row_block_size_data, col_block_size_data)
    CALL fb_com_atom_pairs_nullify(atom_pairs_send)
    CALL fb_com_atom_pairs_nullify(atom_pairs_recv)

    ! initialise atomic matrix
    IF (SIZE(atomic_matrix,1) > 0 .AND. SIZE(atomic_matrix,2) > 0) THEN
       atomic_matrix = 0.0_dp
    END IF

    ! generate send and receiv atomic pairs
    CALL fb_com_atom_pairs_create(atom_pairs_send)
    CALL fb_com_atom_pairs_create(atom_pairs_recv)
    CALL fb_atmatrix_generate_com_pairs(dbcsr_mat, &
                                        atomic_halo, &
                                        para_env, &
                                        atom_pairs_send, &
                                        atom_pairs_recv)

    ! get com pair informations
    CALL fb_com_atom_pairs_get(atom_pairs=atom_pairs_send, &
                               pairs=pairs_send, &
                               npairs=npairs_send, &
                               natoms_encode=send_encode)
    CALL fb_com_atom_pairs_get(atom_pairs=atom_pairs_recv, &
                               pairs=pairs_recv, &
                               npairs=npairs_recv, &
                               natoms_encode=recv_encode)

    ! get para_env info
    numprocs = para_env%num_pe

    ! get dbcsr row and col block sizes
    row_block_size_data => cp_dbcsr_row_block_sizes(dbcsr_mat)
    col_block_size_data => cp_dbcsr_col_block_sizes(dbcsr_mat)

    ! allocate temporary arrays for send
    ALLOCATE(send_sizes(numprocs), STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(send_disps(numprocs), STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(send_pair_count(numprocs), STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(send_pair_disps(numprocs), STAT=stat)
    CPASSERT(stat==0)

    ! setup send buffer sizes
    CALL fb_com_atom_pairs_calc_buffer_sizes(atom_pairs_send, &
                                             numprocs, &
                                             row_block_size_data, &
                                             col_block_size_data, &
                                             send_sizes, &
                                             send_disps, &
                                             send_pair_count, &
                                             send_pair_disps)
    ! allocate send buffer
    ALLOCATE(send_buf(SUM(send_sizes)), STAT=stat)
    CPASSERT(stat==0)

    ! allocate temporary arrays for recv
    ALLOCATE(recv_sizes(numprocs), STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(recv_disps(numprocs), STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(recv_pair_count(numprocs), STAT=stat)
    CPASSERT(stat==0)
    ALLOCATE(recv_pair_disps(numprocs), STAT=stat)
    CPASSERT(stat==0)

    ! setup recv buffer sizes
    CALL fb_com_atom_pairs_calc_buffer_sizes(atom_pairs_recv, &
                                             numprocs, &
                                             row_block_size_data, &
                                             col_block_size_data, &
                                             recv_sizes, &
                                             recv_disps, &
                                             recv_pair_count, &
                                             recv_pair_disps)
    ! allocate recv buffer
    ALLOCATE(recv_buf(SUM(recv_sizes)), STAT=stat)
    CPASSERT(stat==0)
    ! do packing
    DO ipe = 1, numprocs
       ! need to reuse send_sizes as an accumulative displacement, so recalculate
       send_sizes(ipe) = 0
       DO ipair = 1, send_pair_count(ipe)
          CALL fb_com_atom_pairs_decode(pairs_send(send_pair_disps(ipe) + ipair), &
                                        pe, iatom, jatom, send_encode)
          nrows_blk = row_block_size_data(iatom)
          ncols_blk = col_block_size_data(jatom)
          CALL cp_dbcsr_get_block_p(matrix=dbcsr_mat, &
                                    row=iatom, col=jatom, block=mat_block, &
                                    found=found)
          IF (.NOT. found) THEN
             CALL stop_program(routineN, moduleN, __LINE__, "Matrix block not found")
          ELSE
             ! we have found the matrix block
             DO jj = 1, ncols_blk
                DO ii = 1, nrows_blk
                   ! column major format in blocks
                   ind = send_disps(ipe) + send_sizes(ipe) + ii + (jj-1)*nrows_blk
                   send_buf(ind) = mat_block(ii,jj)
                END DO ! ii
             END DO ! jj
             send_sizes(ipe) = send_sizes(ipe) + nrows_blk*ncols_blk
          END IF
       END DO ! ipair
    END DO ! ipe

    ! do communication
    CALL mp_alltoall(send_buf, send_sizes, send_disps, &
                     recv_buf, recv_sizes, recv_disps, &
                     para_env%group)

    ! cleanup temporary arrays no longer needed
    DEALLOCATE(send_buf, STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(send_sizes, STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(send_disps, STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(send_pair_count, STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(send_pair_disps, STAT=stat)
    CPASSERT(stat==0)

    ! do unpacking
    DO ipe = 1, numprocs
       recv_sizes(ipe) = 0
       DO ipair = 1, recv_pair_count(ipe)
          CALL fb_com_atom_pairs_decode(pairs_recv(recv_pair_disps(ipe) + ipair), &
                                        pe, iatom, jatom, recv_encode)
          ! nrows_blk = last_row(iatom) - first_row(iatom) + 1
          ! ncols_blk = last_col(jatom) - first_col(jatom) + 1
          nrows_blk = row_block_size_data(iatom)
          ncols_blk = col_block_size_data(jatom)
          ! get the corresponding atom indices in halo
          ! the atoms from the recv_pairs should be in the atomic_halo, because
          ! the recv_pairs are the matrix blocks requested by the local proc for
          ! this particular atomic_halo
          CALL fb_atomic_halo_atom_global2halo(atomic_halo, &
                                               iatom, iatom_in_halo, &
                                               found)
          CPASSERT(found)
          CALL fb_atomic_halo_atom_global2halo(atomic_halo, &
                                               jatom, jatom_in_halo, &
                                               found)
          CPASSERT(found)
          ! put block into the full conventional matrix
          DO jj = 1, ncols_blk
             DO ii = 1, nrows_blk
                ! column major format in blocks
                ind = recv_disps(ipe) + recv_sizes(ipe) + ii + (jj-1)*nrows_blk
                atomic_matrix(blk_row_start(iatom_in_halo) + ii - 1, &
                              blk_col_start(jatom_in_halo) + jj - 1) = recv_buf(ind)

             END DO ! ii
          END DO ! jj
          recv_sizes(ipe) = recv_sizes(ipe) + nrows_blk*ncols_blk
       END DO ! ipair
    END DO ! ipe

    ! cleanup rest of the temporary arrays
    DEALLOCATE(recv_buf, STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(recv_sizes, STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(recv_disps, STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(recv_pair_count, STAT=stat)
    CPASSERT(stat==0)
    DEALLOCATE(recv_pair_disps, STAT=stat)
    CPASSERT(stat==0)
    CALL fb_com_atom_pairs_release(atom_pairs_send)
    CALL fb_com_atom_pairs_release(atom_pairs_recv)

    CALL timestop(handle)

  END SUBROUTINE fb_atmatrix_construct


! ****************************************************************************
!> \brief Constructs atomic matrix for filter basis method from a given
!>        DBCSR matrix and a set of atomic send and recv pairs
!>        corresponding to the matrix blocks that needs to be included
!>        in the atomic matrix. This version is for when we do MPI
!>        communications collectively in one go at the beginning.
!> \param matrix_storage : data storing the relevant DBCSR matrix blocks
!>                         needed for constructing the atomic matrix
!> \param atomic_halo : the atomic halo conrresponding to this atomic
!>                      matrix
!> \param atomic_matrix : the atomic matrix to be constructed, it should
!>                        have already been allocated prior entering
!>                        this subroutine
!> \param blk_row_start : first row in each atomic blk row in the
!>                        atomic matrix
!> \param blk_col_start : first col in each atomic blk col in the
!>                        atomic matrix
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atmatrix_construct_2(matrix_storage, &
                                     atomic_halo, &
                                     atomic_matrix, &
                                     blk_row_start, &
                                     blk_col_start)
    TYPE(fb_matrix_data_obj), INTENT(IN)     :: matrix_storage
    TYPE(fb_atomic_halo_obj), INTENT(IN)     :: atomic_halo
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: atomic_matrix
    INTEGER, DIMENSION(:), INTENT(IN)        :: blk_row_start, blk_col_start

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fb_atmatrix_construct_2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iatom, iatom_global, &
                                                icol, ii, irow, jatom, &
                                                jatom_global, jj, &
                                                natoms_in_halo
    INTEGER, DIMENSION(:), POINTER           :: halo_atoms
    LOGICAL                                  :: check_ok, failure, found
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: blk_p

    CALL timeset(routineN, handle)

    failure = .FALSE.
    check_ok = fb_matrix_data_has_data(matrix_storage)
    CPASSERT(check_ok)
    check_ok = fb_atomic_halo_has_data(atomic_halo)
    CPASSERT(check_ok)

    NULLIFY(halo_atoms, blk_p)

    ! initialise atomic matrix
    IF (SIZE(atomic_matrix,1) > 0 .AND. SIZE(atomic_matrix,2) > 0) THEN
       atomic_matrix = 0.0_dp
    END IF

    ! get atomic halo information
    CALL fb_atomic_halo_get(atomic_halo=atomic_halo, &
                            natoms=natoms_in_halo, &
                            halo_atoms=halo_atoms)

    ! construct atomic matrix using data from matrix_storage
    DO iatom = 1, natoms_in_halo
       iatom_global = halo_atoms(iatom)
       DO jatom = 1, natoms_in_halo
          jatom_global = halo_atoms(jatom)
          CALL fb_matrix_data_get(matrix_storage, &
                                  iatom_global, &
                                  jatom_global, &
                                  blk_p, &
                                  found)
          ! copy data to atomic_matrix if found
          IF (found) THEN
             DO jj = 1, SIZE(blk_p, 2)
                icol = blk_col_start(jatom) + jj - 1
                DO ii = 1, SIZE(blk_p, 1)
                   irow = blk_row_start(iatom) + ii - 1
                   atomic_matrix(irow,icol) = blk_p(ii,jj)
                END DO ! ii
             END DO ! jj
          END IF
       END DO ! jatom
    END DO ! iatom

    CALL timestop(handle)

  END SUBROUTINE fb_atmatrix_construct_2


! ****************************************************************************
!> \brief generate list of blocks (atom pairs) of a DBCSR matrix to be
!>        sent and recived in order to construct an atomic matrix
!>        corresponding to a given atomic halo. This version is for the case
!>        when we do MPI communications at each step, for each atomic matrix.
!> \param dbcsr_mat : The DBCSR matrix the atom blocks come from
!> \param atomic_halo : the atomic halo used to construct the atomic
!>                      matrix
!> \param para_env : cp2k parallel environment
!> \param atom_pairs_send : list of atom blocks from local DBCSR matrix
!>                          data to be sent
!> \param atom_pairs_recv : list of atom blocks from remote DBCSR matrix
!>                          data to be recveived
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atmatrix_generate_com_pairs(dbcsr_mat, &
                                            atomic_halo, &
                                            para_env, &
                                            atom_pairs_send, &
                                            atom_pairs_recv)
    TYPE(cp_dbcsr_type), POINTER             :: dbcsr_mat
    TYPE(fb_atomic_halo_obj), INTENT(IN)     :: atomic_halo
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(fb_com_atom_pairs_obj), &
      INTENT(INOUT)                          :: atom_pairs_send, &
                                                atom_pairs_recv

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'fb_atmatrix_generate_com_pairs', &
      routineP = moduleN//':'//routineN

    INTEGER :: counter, dest, handle, iatom, iatom_global, itask, jatom, &
      jatom_global, natoms_in_halo, nblkrows_total, nencode, ntasks_recv, &
      ntasks_send, src, stat
    INTEGER(KIND=int_8)                      :: pair
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks_recv, tasks_send
    INTEGER, DIMENSION(:), POINTER           :: halo_atoms
    LOGICAL                                  :: failure, found, transpose
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: mat_block
    TYPE(fb_com_tasks_obj)                   :: com_tasks_recv, com_tasks_send

    CALL timeset(routineN, handle)

    failure = .FALSE.

    NULLIFY(halo_atoms, tasks_send, tasks_recv)
    CALL fb_com_tasks_nullify(com_tasks_send)
    CALL fb_com_tasks_nullify(com_tasks_recv)

    ! initialise atom_pairs_send and atom_pairs_receive
    IF (fb_com_atom_pairs_has_data(atom_pairs_send)) THEN
       CALL fb_com_atom_pairs_init(atom_pairs_send)
    ELSE
       CALL fb_com_atom_pairs_create(atom_pairs_send)
    END IF
    IF (fb_com_atom_pairs_has_data(atom_pairs_recv)) THEN
       CALL fb_com_atom_pairs_init(atom_pairs_recv)
    ELSE
       CALL fb_com_atom_pairs_create(atom_pairs_recv)
    END IF

    ! get atomic halo information
    CALL fb_atomic_halo_get(atomic_halo=atomic_halo, &
                            natoms=natoms_in_halo, &
                            halo_atoms=halo_atoms)

    ! get the total number of atoms, we can obtain this directly
    ! from the global block row dimension of the dbcsr matrix
    CALL cp_dbcsr_get_info(matrix=dbcsr_mat, &
                           nblkrows_total=nblkrows_total)

    ! destination proc is always the local processor
    dest = para_env%mepos  ! my MPI rank

    ! generate recv task list (tasks_recv)

    ! a recv task corresponds to the copying or transfering of a
    ! matrix block in the part of the DBCSR matrix owned by the src
    ! proc to this proc in order to construct the atomic matrix
    ! corresponding to the given atomic halo. As an upper-bound, the
    ! number of matrix blocks requred do not exceed natoms_in_halo**2
    ntasks_recv = natoms_in_halo*natoms_in_halo

    ALLOCATE(tasks_recv(TASK_N_RECORDS,ntasks_recv), STAT=stat)
    CPASSERT(stat==0)

    ! now that tasks_recv has been allocated, generate the tasks
    itask = 1
    DO iatom = 1, natoms_in_halo
       iatom_global = halo_atoms(iatom)
       DO jatom = 1, natoms_in_halo
          jatom_global = halo_atoms(jatom)
          ! find the source proc that supposed to own the block
          ! (iatom_global, jatom_global)
          transpose = .FALSE.
          CALL cp_dbcsr_get_stored_coordinates(dbcsr_mat, &
                                               iatom_global, &
                                               jatom_global, &
                                               transpose=transpose, &
                                               processor=src)
          ! we must encode the global atom indices rather the halo atomic indices
          ! in each task, because halo atomic indices are local to each halo, and
          ! each processor is working on a different halo local to them. So
          ! one processor would not have the information about the halo on another
          ! processor, rendering the halo atomic indices rather useless outside
          ! the local processor.
          tasks_recv(TASK_DEST,itask) = dest
          tasks_recv(TASK_SRC,itask) = src

          CALL fb_com_tasks_encode_pair(tasks_recv(TASK_PAIR,itask), &
                                        iatom_global, jatom_global, &
                                        nblkrows_total)
          ! calculation of cost not implemented at the moment
          tasks_recv(TASK_COST,itask) = 0
          itask = itask + 1
       END DO ! jatom
    END DO ! iatom

    ! create tasks
    CALL fb_com_tasks_create(com_tasks_recv)
    CALL fb_com_tasks_create(com_tasks_send)

    CALL fb_com_tasks_set(com_tasks=com_tasks_recv, &
                          task_dim=TASK_N_RECORDS, &
                          ntasks=ntasks_recv, &
                          nencode=nblkrows_total, &
                          tasks=tasks_recv)

    ! genearte the send task list (tasks_send) from the recv task list
    CALL fb_com_tasks_transpose_dest_src(com_tasks_recv, ">", com_tasks_send, &
                                         para_env)


    CALL fb_com_tasks_get(com_tasks=com_tasks_send, &
                          ntasks=ntasks_send, &
                          tasks=tasks_send, &
                          nencode=nencode)

    ! because the atomic_halos and the neighbor_list_set used to
    ! generate the sparse structure of the DBCSR matrix do not
    ! necessarily have to coincide, we must check of the blocks in
    ! tasks_send (these should be local to the processor) do indeed
    ! exist in the DBCSR matrix, if not, then we need to prune these
    ! out of the task list

    counter = 0
    DO itask = 1, ntasks_send
       pair = tasks_send(TASK_PAIR,itask)
       CALL fb_com_tasks_decode_pair(pair, iatom_global, jatom_global, nencode)
       ! check if block exists in DBCSR matrix
       CALL cp_dbcsr_get_block_p(matrix=dbcsr_mat, &
                                 row=iatom_global, col=jatom_global, block=mat_block, &
                                 found=found)
       IF (found) THEN
          counter = counter + 1
          ! we can do this here, because essencially we are inspecting
          ! the send tasks one by one, and then omit ones which the
          ! block is not found in the DBCSR matrix. itask is always
          ! .GE. counter
          tasks_send(1:TASK_N_RECORDS,counter) = tasks_send(1:TASK_N_RECORDS,itask)
       END IF
    END DO
    ! the new send task list should have size counter. counter
    ! .LE. the old ntasks_send, thus the task list does not really
    ! need to be reallocated (as it is just a temporary array), and
    ! the useful data will cutoff at counter, and the rest of the
    ! array will just be garbage
    ntasks_send = counter

    ! tasks_send is set through the pointer already
    CALL fb_com_tasks_set(com_tasks=com_tasks_send, &
                          ntasks=ntasks_send)

    ! now, re-distribute the new send tasks list to other processors
    ! to build the updated recv tasks list
    CALL fb_com_tasks_transpose_dest_src(com_tasks_recv, "<", com_tasks_send, &
                                         para_env)

    ! task lists are now complete, now construct the atom_pairs_send
    ! and atom_pairs_recv from the tasks lists
    CALL fb_com_tasks_build_atom_pairs(com_tasks=com_tasks_send, &
                                       atom_pairs=atom_pairs_send, &
                                       natoms_encode=nencode, &
                                       send_or_recv="send", &
                                       symmetric=.FALSE.)
    CALL fb_com_tasks_build_atom_pairs(com_tasks=com_tasks_recv, &
                                       atom_pairs=atom_pairs_recv, &
                                       natoms_encode=nencode, &
                                       send_or_recv="recv", &
                                       symmetric=.FALSE.)

    ! cleanup
    CALL fb_com_tasks_release(com_tasks_recv)
    CALL fb_com_tasks_release(com_tasks_send)

    CALL timestop(handle)

  END SUBROUTINE fb_atmatrix_generate_com_pairs


! ****************************************************************************
!> \brief generate list of blocks (atom pairs) of a DBCSR matrix to be
!>        sent and recived in order to construct all local atomic matrices
!>        corresponding to the atomic halos. This version is for the case
!>        when we do MPI communications collectively in one go at the
!>        beginning.
!> \param dbcsr_mat : The DBCSR matrix the atom blocks come from
!> \param atomic_halos : the list of all atomic halos local to the process
!> \param para_env : cp2k parallel environment
!> \param atom_pairs_send : list of atom blocks from local DBCSR matrix
!>                          data to be sent
!> \param atom_pairs_recv : list of atom blocks from remote DBCSR matrix
!>                          data to be recveived
!> \author Lianheng Tong (LT) lianheng.tong@kcl.ac.uk
! *****************************************************************************
  SUBROUTINE fb_atmatrix_generate_com_pairs_2(dbcsr_mat, &
                                              atomic_halos, &
                                              para_env, &
                                              atom_pairs_send, &
                                              atom_pairs_recv)
    TYPE(cp_dbcsr_type), POINTER             :: dbcsr_mat
    TYPE(fb_atomic_halo_list_obj), &
      INTENT(IN)                             :: atomic_halos
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(fb_com_atom_pairs_obj), &
      INTENT(INOUT)                          :: atom_pairs_send, &
                                                atom_pairs_recv

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'fb_atmatrix_generate_com_pairs_2', &
      routineP = moduleN//':'//routineN

    INTEGER :: counter, dest, handle, iatom, iatom_global, ihalo, itask, &
      jatom, jatom_global, natoms_in_halo, nblkrows_total, nencode, nhalos, &
      ntasks_recv, ntasks_send, src, stat
    INTEGER(KIND=int_8)                      :: pair
    INTEGER(KIND=int_8), DIMENSION(:, :), &
      POINTER                                :: tasks_recv, tasks_send
    INTEGER, DIMENSION(:), POINTER           :: halo_atoms
    LOGICAL                                  :: failure, found, transpose
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: mat_block
    TYPE(fb_atomic_halo_obj), DIMENSION(:), &
      POINTER                                :: halos
    TYPE(fb_com_tasks_obj)                   :: com_tasks_recv, com_tasks_send

    CALL timeset(routineN, handle)

    failure = .FALSE.
    
    NULLIFY(halo_atoms, tasks_send, tasks_recv)
    CALL fb_com_tasks_nullify(com_tasks_send)
    CALL fb_com_tasks_nullify(com_tasks_recv)

    ! initialise atom_pairs_send and atom_pairs_receive
    IF (fb_com_atom_pairs_has_data(atom_pairs_send)) THEN
       CALL fb_com_atom_pairs_init(atom_pairs_send)
    ELSE
       CALL fb_com_atom_pairs_create(atom_pairs_send)
    END IF
    IF (fb_com_atom_pairs_has_data(atom_pairs_recv)) THEN
       CALL fb_com_atom_pairs_init(atom_pairs_recv)
    ELSE
       CALL fb_com_atom_pairs_create(atom_pairs_recv)
    END IF

    ! get atomic halo list information
    CALL fb_atomic_halo_list_get(atomic_halos=atomic_halos, &
                                 nhalos=nhalos, &
                                 halos=halos)
    ! get the total number of atoms, we can obtain this directly
    ! from the global block row dimension of the dbcsr matrix
    CALL cp_dbcsr_get_info(matrix=dbcsr_mat, &
                           nblkrows_total=nblkrows_total)

    ! estimate the maximum number of blocks to be received
    ntasks_recv = 0
    DO ihalo = 1, nhalos
       CALL fb_atomic_halo_get(atomic_halo=halos(ihalo), &
                               natoms=natoms_in_halo)
       ntasks_recv = ntasks_recv + natoms_in_halo*natoms_in_halo
    END DO
    ALLOCATE(tasks_recv(TASK_N_RECORDS,ntasks_recv), STAT=stat)
    CPASSERT(stat==0)

    ! now that tasks_recv has been allocated, generate the tasks

    ! destination proc is always the local process
    dest = para_env%mepos
    itask = 1
    DO ihalo = 1, nhalos
       CALL fb_atomic_halo_get(atomic_halo=halos(ihalo), &
                               natoms=natoms_in_halo, &
                               halo_atoms=halo_atoms)
       DO iatom = 1, natoms_in_halo
          iatom_global = halo_atoms(iatom)
          DO jatom = 1, natoms_in_halo
             jatom_global = halo_atoms(jatom)
             ! find the source proc that supposed to own the block
             ! (iatom_global, jatom_global)
             transpose = .FALSE.
             CALL cp_dbcsr_get_stored_coordinates(dbcsr_mat, &
                                                  iatom_global, &
                                                  jatom_global, &
                                                  transpose=transpose, &
                                                  processor=src)
             ! we must encode the global atom indices rather the halo
             ! atomic indices in each task, because halo atomic indices
             ! are local to each halo, and each processor is working on a
             ! different halo local to them. So one processor would not
             ! have the information about the halo on another processor,
             ! rendering the halo atomic indices rather useless outside
             ! the local processor.
             tasks_recv(TASK_DEST,itask) = dest
             tasks_recv(TASK_SRC,itask) = src
             CALL fb_com_tasks_encode_pair(tasks_recv(TASK_PAIR,itask), &
                                           iatom_global, jatom_global, &
                                           nblkrows_total)
             ! calculation of cost not implemented at the moment
             tasks_recv(TASK_COST,itask) = 0
             itask = itask + 1
          END DO ! jatom
       END DO ! iatom
    END DO ! ihalo

    ! create tasks
    CALL fb_com_tasks_create(com_tasks_recv)
    CALL fb_com_tasks_create(com_tasks_send)

    CALL fb_com_tasks_set(com_tasks=com_tasks_recv, &
                          task_dim=TASK_N_RECORDS, &
                          ntasks=ntasks_recv, &
                          nencode=nblkrows_total, &
                          tasks=tasks_recv)

    ! genearte the send task list (tasks_send) from the recv task list
    CALL fb_com_tasks_transpose_dest_src(com_tasks_recv, ">", com_tasks_send, &
                                         para_env)

    CALL fb_com_tasks_get(com_tasks=com_tasks_send, &
                          ntasks=ntasks_send, &
                          tasks=tasks_send, &
                          nencode=nencode)

    ! because the atomic_halos and the neighbor_list_set used to
    ! generate the sparse structure of the DBCSR matrix do not
    ! necessarily have to coincide, we must check of the blocks in
    ! tasks_send (these should be local to the processor) do indeed
    ! exist in the DBCSR matrix, if not, then we need to prune these
    ! out of the task list

    counter = 0
    DO itask = 1, ntasks_send
       pair = tasks_send(TASK_PAIR,itask)
       CALL fb_com_tasks_decode_pair(pair, iatom_global, jatom_global, nencode)
       ! check if block exists in DBCSR matrix
       CALL cp_dbcsr_get_block_p(matrix=dbcsr_mat, row=iatom_global, &
                                 col=jatom_global, block=mat_block, &
                                 found=found)
       IF (found) THEN
          counter = counter + 1
          ! we can do this here, because essencially we are inspecting
          ! the send tasks one by one, and then omit ones which the
          ! block is not found in the DBCSR matrix. itask is always
          ! .GE. counter
          tasks_send(1:TASK_N_RECORDS,counter) = tasks_send(1:TASK_N_RECORDS,itask)
       END IF
    END DO
    ! the new send task list should have size counter. counter
    ! .LE. the old ntasks_send, thus the task list does not really
    ! need to be reallocated (as it is just a temporary array), and
    ! the useful data will cutoff at counter, and the rest of the
    ! array will just be garbage
    ntasks_send = counter

    ! tasks_send is set through the pointer already
    CALL fb_com_tasks_set(com_tasks=com_tasks_send, &
                          ntasks=ntasks_send)

    ! now, re-distribute the new send tasks list to other processors
    ! to build the updated recv tasks list
    CALL fb_com_tasks_transpose_dest_src(com_tasks_recv, "<", com_tasks_send, &
                                         para_env)

    ! task lists are now complete, now construct the atom_pairs_send
    ! and atom_pairs_recv from the tasks lists
    CALL fb_com_tasks_build_atom_pairs(com_tasks=com_tasks_send, &
                                       atom_pairs=atom_pairs_send, &
                                       natoms_encode=nencode, &
                                       send_or_recv="send", &
                                       symmetric=.FALSE.)
    CALL fb_com_tasks_build_atom_pairs(com_tasks=com_tasks_recv, &
                                       atom_pairs=atom_pairs_recv, &
                                       natoms_encode=nencode, &
                                       send_or_recv="recv", &
                                       symmetric=.FALSE.)

    ! cleanup
    CALL fb_com_tasks_release(com_tasks_recv)
    CALL fb_com_tasks_release(com_tasks_send)

    CALL timestop(handle)

  END SUBROUTINE fb_atmatrix_generate_com_pairs_2


END MODULE qs_fb_atomic_matrix_methods
