!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to convert sparse matrices between DBCSR (distributed-blocks compressed sparse rows)
!>        and SIESTA (distributed compressed sparse columns) formats.
!> \author Sergey Chulkov
!> \author Christian Ahart
!> \author Clotilde Cucinotta
! **************************************************************************************************
MODULE smeagol_matrix_utils
   USE cell_types, ONLY: cell_type, &
                         real_to_scaled, &
                         scaled_to_real
   USE cp_dbcsr_api, ONLY: dbcsr_get_block_p, &
                           dbcsr_get_info, &
                           dbcsr_p_type, &
                           dbcsr_set
   USE kinds, ONLY: dp, &
                    dp_size, &
                    int_8
   USE message_passing, ONLY: mp_para_env_type, &
                              mp_request_type, &
                              mp_waitall
   USE negf_matrix_utils, ONLY: get_index_by_cell
#if defined(__SMEAGOL)
   USE parallel, ONLY: GetNodeOrbs, &
                       GlobalToLocalOrb, &
                       LocalToGlobalOrb, &
                       WhichNodeOrb
#endif
   USE particle_types, ONLY: particle_type
   USE qs_neighbor_list_types, ONLY: get_iterator_info, &
                                     get_neighbor_list_set_p, &
                                     neighbor_list_iterate, &
                                     neighbor_list_iterator_create, &
                                     neighbor_list_iterator_p_type, &
                                     neighbor_list_iterator_release, &
                                     neighbor_list_set_p_type
   USE qs_subsys_types, ONLY: qs_subsys_get, &
                              qs_subsys_type
   USE util, ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'smeagol_matrix_utils'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

   INTEGER, PARAMETER, PRIVATE          :: neighbor_list_iatom_index = 1
   INTEGER, PARAMETER, PRIVATE          :: neighbor_list_jatom_index = 2
   INTEGER, PARAMETER, PRIVATE          :: neighbor_list_dbcsr_image_index = 3
   INTEGER, PARAMETER, PRIVATE          :: neighbor_list_siesta_image_index = 4
   INTEGER, PARAMETER, PRIVATE          :: neighbor_list_siesta_transp_image_index = 5
   INTEGER, PARAMETER, PRIVATE          :: neighbor_list_dim1 = neighbor_list_siesta_transp_image_index

   PUBLIC :: siesta_distrib_csc_struct_type
   PUBLIC :: siesta_struct_create, siesta_struct_release
   PUBLIC :: convert_dbcsr_to_distributed_siesta, convert_distributed_siesta_to_dbcsr

   PRIVATE :: get_negf_cell_ijk, index_in_canonical_enumeration, number_from_canonical_enumeration, pbc_0_1
   PRIVATE :: get_number_of_mpi_sendrecv_requests, assign_nonzero_elements_to_requests

   !> number of DBCSR matrix elements to receive from a given rank
   INTEGER, PARAMETER, PRIVATE          :: nelements_dbcsr_recv = 1
   !> number of DBCSR matrix elements to send to a given rank
   INTEGER, PARAMETER, PRIVATE          :: nelements_dbcsr_send = 2
   INTEGER, PARAMETER, PRIVATE          :: nelements_dbcsr_dim2 = nelements_dbcsr_send

   INTEGER(kind=int_8), PARAMETER, PRIVATE :: max_mpi_packet_size_bytes = 134217728 ! 128 MiB (to limit memory usage for matrix redistribution)
   INTEGER(kind=int_8), PARAMETER, PRIVATE :: max_mpi_packet_size_dp = max_mpi_packet_size_bytes/INT(dp_size, kind=int_8)

   ! a portable way to determine the upper bound for tag value is to call
   ! MPI_COMM_GET_ATTR(comm, MPI_TAG_UB, max_mpi_rank, flag, ierror).
   ! The MPI specification guarantees a value of 32767
   INTEGER, PARAMETER, PRIVATE          :: max_mpi_rank = 32767

! **************************************************************************************************
!> \brief Sparsity pattern of replicated SIESTA compressed sparse column (CSC) matrices
! **************************************************************************************************
   TYPE siesta_distrib_csc_struct_type
      !> gather all non-zero matrix elements on the given MPI process.
      !> Distribute the elements across MPI processes if gather_root < 0.
      !> The matrix elements should be located on I/O process in case of bulk transport,
      !> and should be distributed in case of SMEAGOL calculation.
      INTEGER                                            :: gather_root = 0
      !> Based of full (.FALSE.) or upper-triangular (.TRUE.) DBCSR matrix.
      !> It is used in CPASSERTs to allow access to lower-triangular matrix elements of non-symmetric DBCSR matrices
      !> In case there is no bugs in this module, these CPASSERTs should newer trigger.
      !> Therefore the 'symmetric' variable alongside with relevant CPASSERT calls are excessive and can be removed.
      LOGICAL                                            :: symmetric = .TRUE.
      !> number of neighbour list nodes for each MPI process (0:num_pe-1).
      !> If do_merge == .TRUE., nodes for different cell images along k cell vector are merged into one node
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nnodes_per_proc

      !> replicated neighbour list (1:neighbor_list_dim1, 1:SUM(nnodes_per_proc)).
      !> Neighbour list nodes are ordered according to their MPI ranks.
      !> Thus, the first nnodes_per_proc(0) nodes are stored on MPI rank 0,
      !> the next nnodes_per_proc(1) nodes reside on MPI rank 1, etc
      !> Nodes for cell images along transport direction are merged into one node.
      !> The number of non-zero DBCSR matrix blocks and their DBCSR cell image indices
      !> are stored into 'n_dbcsr_cell_images_to_merge' and 'dbcsr_cell_image_to_merge' arrays
      INTEGER, ALLOCATABLE, DIMENSION(:, :)               :: nl_repl

      !> number of DBCSR images for each local merged neighbour list node (1:nnodes_per_proc(para_env%mepos))
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: n_dbcsr_cell_images_to_merge
      !> list of DBCSR image indices to merge; (1:SUM(n_dbcsr_cell_images_to_merge))
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dbcsr_cell_image_to_merge

      !> number of DBCSR non-zero matrix elements that should be received/sent from each MPI rank
      !> (0:num_pe-1, 1:nelements_dbcsr_dim2)
      INTEGER(kind=int_8), ALLOCATABLE, DIMENSION(:, :)   :: nelements_per_proc

      !> number of non-zero matrix elements local to this MPI rank.
      INTEGER(kind=int_8)                                :: n_nonzero_elements = 0_int_8
      INTEGER                                            :: nrows = 0, ncols = 0
      !> Number of non-zero matrix elements (columns) on each row
      !> n_nonzero_cols(1:nrows); same as 'numh' in SMEAGOL code.
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: n_nonzero_cols
      !> offset of the first non-zero matrix elements on each row.
      !> column_offset(1:nrows); same as 'listhptr' in SMEAGOL code.
      !> It should be declared as INTEGER(kind=int_8), but SMEAGOL expects it to be INTEGER
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: row_offset
      !> column index of each non-zero matrix element.
      !> col_index(1:n_nonzero_elements); same as 'listh' in SMEAGOL code
      !> col index of the first non-zero matrix elements of irow row is col_index(row_offset(irow)+1)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: col_index
      !> index of the non-zero matrix element in a communication buffer for each rank
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: packed_index
      !> R_atom_row - R_atom_col
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)         :: xij
      !> equivalent atomic orbitals
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: indxuo
      !> atomic index on which the orbital is centred
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iaorb
      !> coordinates of all atoms in the supercell
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)         :: xa
   END TYPE siesta_distrib_csc_struct_type

CONTAINS

! **************************************************************************************************
!> \brief Map non-zero matrix blocks between sparse matrices in DBCSR and SIESTA formats.
!> \param siesta_struct      structure that stores metadata (sparsity pattern) of sparse SIESTA matrices
!> \param matrix_dbcsr_kp    DBCSR matrices for each cell image
!> \param subsys             QuickStep molecular system
!> \param cell_to_index      array to convert 3-D cell indices to 1-D DBCSR image indices
!> \param sab_nl             pair-wise neighbour list
!> \param para_env           MPI parallel environment
!> \param max_ij_cell_image  largest index of cell images along i and j cell vectors (e.g. (2,0) in
!>                           case of 5 cell images (0,0), (1,0), (-1,0), (2,0), and (-2,0))
!> \param do_merge           merge DBCSR images along transport direction (k cell vector)
!> \param gather_root        distribute non-zero matrix elements of SIESTA matrices across all
!>                           parallel processes (-1), or gather them on the given MPI rank (>= 0).
! **************************************************************************************************
   SUBROUTINE siesta_struct_create(siesta_struct, matrix_dbcsr_kp, subsys, cell_to_index, &
                                   sab_nl, para_env, max_ij_cell_image, do_merge, gather_root)
      TYPE(siesta_distrib_csc_struct_type), &
         INTENT(inout)                                   :: siesta_struct
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(in)       :: matrix_dbcsr_kp
      TYPE(qs_subsys_type), POINTER                      :: subsys
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, DIMENSION(2), INTENT(inout)               :: max_ij_cell_image
      LOGICAL, INTENT(in)                                :: do_merge
      INTEGER, INTENT(in)                                :: gather_root

      CHARACTER(len=*), PARAMETER :: routineN = 'siesta_struct_create'

      CHARACTER(len=20)                                  :: str_nelem, str_nelem_max
      INTEGER :: handle, iatom, icol, icol_blk, icol_local, image, image_j, image_k, irow, &
         irow_local, natoms, ncells_siesta_total, ncols_blk, ncols_total, nrows_local, &
         nrows_total, offset
      INTEGER(kind=int_8)                                :: n_nonzero_elements_local
      INTEGER, DIMENSION(3)                              :: max_ijk_cell_image, ncells_siesta
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_offset, col_blk_size
      LOGICAL                                            :: do_distribute, is_root_rank
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: particle_coords
      REAL(kind=dp), DIMENSION(3)                        :: real_cell_shift, scaled_cell_shift
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)
      do_distribute = gather_root < 0
      is_root_rank = gather_root == para_env%mepos

      ! here row_blk_offset / col_blk_offset are global indices of the first row / column of a given non-zero block.
      ! They are not offsets (index-1) but the actual indices.
      CALL dbcsr_get_info(matrix=matrix_dbcsr_kp(1)%matrix, &
                          nfullrows_total=nrows_total, nfullcols_total=ncols_total, &
                          nblkcols_total=ncols_blk, col_blk_size=col_blk_size, col_blk_offset=col_blk_offset)
      IF (debug_this_module) THEN
         CPASSERT(nrows_total == ncols_total)
         CPASSERT(gather_root < para_env%num_pe)
      END IF

      siesta_struct%gather_root = gather_root
      CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl, symmetric=siesta_struct%symmetric)

      ! apply periodic boundary conditions to atomic coordinates
      CALL qs_subsys_get(subsys, cell=cell, particle_set=particle_set, nparticle=natoms)
      ALLOCATE (particle_coords(3, natoms))
      DO iatom = 1, natoms
         CALL pbc_0_1(particle_coords(1:3, iatom), particle_set(iatom)%r(1:3), cell)
      END DO

      ! Note: in case we would like to limit the number of cell images along transport direction (k cell vector)
      !       by enabling 'BulkTransvCellSizeZ' keyword, we need to pass max_ijk_cell_image(1:3) vector
      !       to the subroutine instead of the reduced vector max_ij_cell_image(1:2)
      max_ijk_cell_image(1:2) = max_ij_cell_image(1:2)

      ! determine the actual number of cell images along k cell vector. Unless the third element is also passed
      ! via subroutine arguments, an extra MPI_Allreduce operation is needed each time we call
      ! replicate_neighbour_list() / get_nnodes_local().
      max_ijk_cell_image(3) = -1
      IF (.NOT. do_merge) max_ijk_cell_image(3) = 1 ! bulk-transport calculation expects exactly 3 cell images along transport direction

      ! replicate pair-wise neighbour list. Identical non-zero matrix blocks from cell image along transport direction
      ! are grouped together if do_merge == .TRUE.
      ALLOCATE (siesta_struct%nnodes_per_proc(0:para_env%num_pe - 1))
      CALL replicate_neighbour_list(siesta_struct%nl_repl, &
                                    siesta_struct%n_dbcsr_cell_images_to_merge, &
                                    siesta_struct%dbcsr_cell_image_to_merge, &
                                    siesta_struct%nnodes_per_proc, &
                                    max_ijk_cell_image, sab_nl, para_env, particle_coords, cell, cell_to_index, do_merge)
      max_ij_cell_image(1:2) = max_ijk_cell_image(1:2)

      ! count number of non-zero matrix elements that need to be send to and received from other parallel processes
      ALLOCATE (siesta_struct%nelements_per_proc(0:para_env%num_pe - 1, nelements_dbcsr_dim2))
      CALL count_remote_dbcsr_elements(siesta_struct%nelements_per_proc, siesta_struct%nnodes_per_proc, &
                                       siesta_struct%nl_repl, matrix_dbcsr_kp, siesta_struct%symmetric, para_env, gather_root)

      ! number of SIESTA non-zero matrix elements that are going to be stored on this parallel process
      n_nonzero_elements_local = SUM(siesta_struct%nelements_per_proc(:, nelements_dbcsr_recv))
      siesta_struct%n_nonzero_elements = n_nonzero_elements_local

      ! as SMEAGOL uses 32-bits integers, the number of non-zero matrix elements is limited by 2^31 per MPI rank.
      ! Abort CP2K if we are about to exceed this limit.
      IF (n_nonzero_elements_local > INT(HUGE(0), kind=int_8)) THEN
         WRITE (str_nelem, '(I0)') n_nonzero_elements_local
         WRITE (str_nelem_max, '(I0)') HUGE(0)
         CALL cp_abort(__LOCATION__, &
                       "The number of non-zero matrix elements per MPI process "//TRIM(str_nelem)// &
                       " cannot exceed "//TRIM(str_nelem_max)// &
                       ". Please increase the number of MPI processes to satisfy this SMEAGOL limitation.")
      END IF

      ! in case there is no nonzero matrix element stored on this process, allocate arrays with one element to avoid SEGFAULT
      IF (n_nonzero_elements_local == 0) n_nonzero_elements_local = 1

      ! number of SIESTA-matrix rows local to the given parallel process
      IF (do_distribute) THEN
#if defined(__SMEAGOL)
         CALL GetNodeOrbs(nrows_total, para_env%mepos, para_env%num_pe, nrows_local)
#else
         CALL cp_abort(__LOCATION__, &
                       "CP2K was compiled with no SMEAGOL support.")
#endif
      ELSE
         IF (is_root_rank) THEN
            nrows_local = nrows_total
         ELSE
            nrows_local = 0
         END IF
      END IF

      ! number of cell images along each cell vector. It is 2*m+1, as SIESTA images are ordered as 0, 1, -1, ..., m, -m
      ncells_siesta(1:3) = 2*max_ijk_cell_image(1:3) + 1
      ! in case of merged cell images along the transport direction, there will be just 1 'merged' cell image along it
      IF (do_merge) ncells_siesta(3) = 1

      ncells_siesta_total = ncells_siesta(1)*ncells_siesta(2)*ncells_siesta(3)

      ! number of rows local to the given parallel process. Rows are distributed in a block-cyclic manner
      siesta_struct%nrows = nrows_local

      ! number of columns of the matrix in its dense form. SIESTA uses 1-D (rows) block-cyclic distribution.
      ! All non-zero matrix elements on a given row are stored on the same parallel process
      siesta_struct%ncols = nrows_total*ncells_siesta_total

      ! allocate at least one array element to avoid SIGFAULT when passing unallocated arrays to subroutines
      IF (nrows_local == 0) nrows_local = 1

      ALLOCATE (siesta_struct%n_nonzero_cols(nrows_local))
      ALLOCATE (siesta_struct%row_offset(nrows_local))
      ALLOCATE (siesta_struct%col_index(n_nonzero_elements_local))
      ALLOCATE (siesta_struct%packed_index(n_nonzero_elements_local))

      ! restore the actual number of local rows
      nrows_local = siesta_struct%nrows

      ! get number of non-zero matrix element on each local row (n_nonzero_cols),
      ! offset of the first non-zero matrix element for each local row (row_offset),
      ! global column indices of all local non-zero matrix elements (col_index), and
      ! the indices of all local non-zero matrix elements in the communication buffer (packed_index)
      CALL get_nonzero_element_indices(siesta_struct%n_nonzero_cols, siesta_struct%row_offset, &
                                       siesta_struct%col_index, siesta_struct%packed_index, &
                                       siesta_struct%nl_repl, matrix_dbcsr_kp, &
                                       siesta_struct%symmetric, para_env, gather_root)

      ! indices of equivalent atomic orbitals
      ALLOCATE (siesta_struct%indxuo(siesta_struct%ncols))
      DO icol = 1, ncols_total
         siesta_struct%indxuo(icol) = icol
      END DO
      DO image = 2, ncells_siesta_total
         siesta_struct%indxuo((image - 1)*ncols_total + 1:image*ncols_total) = siesta_struct%indxuo(1:ncols_total)
      END DO

      ! particle index on which the orbital is centred
      ALLOCATE (siesta_struct%iaorb(siesta_struct%ncols))
      DO icol_blk = 1, ncols_blk
         ! col_blk_offset() is not an offset but the index of the first atomic orbital in the column block
         siesta_struct%iaorb(col_blk_offset(icol_blk):col_blk_offset(icol_blk) + col_blk_size(icol_blk) - 1) = icol_blk
      END DO
      DO image = 2, ncells_siesta_total
        siesta_struct%iaorb((image - 1)*ncols_total + 1:image*ncols_total) = siesta_struct%iaorb(1:ncols_total) + (image - 1)*natoms
      END DO

      ! coordinates of all particles in each cell images
      ALLOCATE (siesta_struct%xa(3, natoms*ncells_siesta_total))
      DO image_k = 1, ncells_siesta(3)
         !icell_siesta(3) = image_k
         scaled_cell_shift(3) = REAL(number_from_canonical_enumeration(image_k), kind=dp) ! SIESTA -> actual cell index
         DO image_j = 1, ncells_siesta(2)
            scaled_cell_shift(2) = REAL(number_from_canonical_enumeration(image_j), kind=dp)
            DO image = 1, ncells_siesta(1)
               scaled_cell_shift(1) = REAL(number_from_canonical_enumeration(image), kind=dp)
               CALL scaled_to_real(real_cell_shift, scaled_cell_shift, cell)
               offset = (((image_k - 1)*ncells_siesta(2) + image_j - 1)*ncells_siesta(1) + image - 1)*natoms
               DO iatom = 1, natoms
                  siesta_struct%xa(1:3, offset + iatom) = particle_set(iatom)%r(1:3) + real_cell_shift(1:3)
               END DO
            END DO
         END DO
      END DO

      ! inter-atomic distance
      ALLOCATE (siesta_struct%xij(3, n_nonzero_elements_local))
      DO irow_local = 1, nrows_local
         IF (do_distribute) THEN
#if defined(__SMEAGOL)
            CALL LocalToGlobalOrb(irow_local, para_env%mepos, para_env%num_pe, irow)
#else
            CALL cp_abort(__LOCATION__, &
                          "CP2K was compiled with no SMEAGOL support.")
#endif
         ELSE
            irow = irow_local
            IF (debug_this_module) THEN
               CPASSERT(is_root_rank)
            END IF
         END IF
         offset = siesta_struct%row_offset(irow_local)
         DO icol_local = offset + 1, offset + siesta_struct%n_nonzero_cols(irow_local)
            icol = siesta_struct%col_index(icol_local)
            siesta_struct%xij(1:3, icol_local) = siesta_struct%xa(1:3, siesta_struct%iaorb(icol)) - &
                                                 siesta_struct%xa(1:3, siesta_struct%iaorb(irow))
         END DO
      END DO

      DEALLOCATE (particle_coords)

      CALL timestop(handle)
   END SUBROUTINE siesta_struct_create

! **************************************************************************************************
!> \brief Release a SIESTA matrix structure
!> \param siesta_struct      structure to release
! **************************************************************************************************
   SUBROUTINE siesta_struct_release(siesta_struct)
      TYPE(siesta_distrib_csc_struct_type), &
         INTENT(inout)                                   :: siesta_struct

      CHARACTER(len=*), PARAMETER :: routineN = 'siesta_struct_release'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      siesta_struct%gather_root = -1

      IF (ALLOCATED(siesta_struct%nnodes_per_proc)) DEALLOCATE (siesta_struct%nnodes_per_proc)
      IF (ALLOCATED(siesta_struct%nl_repl)) DEALLOCATE (siesta_struct%nl_repl)
      IF (ALLOCATED(siesta_struct%n_dbcsr_cell_images_to_merge)) DEALLOCATE (siesta_struct%n_dbcsr_cell_images_to_merge)
      IF (ALLOCATED(siesta_struct%dbcsr_cell_image_to_merge)) DEALLOCATE (siesta_struct%dbcsr_cell_image_to_merge)
      IF (ALLOCATED(siesta_struct%nelements_per_proc)) DEALLOCATE (siesta_struct%nelements_per_proc)

      siesta_struct%n_nonzero_elements = 0
      siesta_struct%nrows = 0
      siesta_struct%ncols = 0

      IF (ALLOCATED(siesta_struct%n_nonzero_cols)) DEALLOCATE (siesta_struct%n_nonzero_cols)
      IF (ALLOCATED(siesta_struct%row_offset)) DEALLOCATE (siesta_struct%row_offset)
      IF (ALLOCATED(siesta_struct%col_index)) DEALLOCATE (siesta_struct%col_index)
      IF (ALLOCATED(siesta_struct%packed_index)) DEALLOCATE (siesta_struct%packed_index)

      IF (ALLOCATED(siesta_struct%xij)) DEALLOCATE (siesta_struct%xij)
      IF (ALLOCATED(siesta_struct%indxuo)) DEALLOCATE (siesta_struct%indxuo)
      IF (ALLOCATED(siesta_struct%iaorb)) DEALLOCATE (siesta_struct%iaorb)
      IF (ALLOCATED(siesta_struct%xa)) DEALLOCATE (siesta_struct%xa)

      CALL timestop(handle)
   END SUBROUTINE siesta_struct_release

! **************************************************************************************************
!> \brief Convert matrix from DBCSR to sparse SIESTA format.
!> \param matrix_siesta      matrix in SIESTA format [out]
!> \param matrix_dbcsr_kp    DBCSR matrix [in]
!> \param siesta_struct      structure to map matrix blocks between formats
!> \param para_env           MPI parallel environment
! **************************************************************************************************
   SUBROUTINE convert_dbcsr_to_distributed_siesta(matrix_siesta, matrix_dbcsr_kp, siesta_struct, para_env)
      REAL(kind=dp), DIMENSION(:), INTENT(out)           :: matrix_siesta
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(in)       :: matrix_dbcsr_kp
      TYPE(siesta_distrib_csc_struct_type), INTENT(in)   :: siesta_struct
      TYPE(mp_para_env_type), INTENT(in), POINTER        :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'convert_dbcsr_to_distributed_siesta'

      INTEGER :: first_col_minus_one, first_row_minus_one, handle, icol_blk, icol_local, &
         image_dbcsr, image_ind, image_ind_offset, image_siesta, image_siesta_transp, inode, &
         inode_proc, iproc, irequest, irow_blk, irow_local, irow_proc, mepos, n_image_ind, &
         ncols_blk, ncols_local, nnodes_proc, node_offset, nprocs, nrequests_recv, &
         nrequests_total, nrows_blk, nrows_local
      INTEGER(kind=int_8)                                :: n_nonzero_elements_dbcsr, &
                                                            n_nonzero_elements_siesta, &
                                                            offset_recv_mepos, offset_send_mepos
      INTEGER(kind=int_8), ALLOCATABLE, DIMENSION(:)     :: n_packed_elements_per_proc, &
                                                            nelements_per_request, &
                                                            offset_per_proc, offset_per_request
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: next_nonzero_element_offset, peer_rank, &
                                                            request_tag
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_offset, col_blk_size, &
                                                            row_blk_offset, row_blk_size
      LOGICAL                                            :: do_distribute, found, is_root_rank, &
                                                            symmetric
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: recv_buffer, reorder_recv_buffer, &
                                                            send_buffer
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: sm_block, sm_block_merged
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: requests

      CALL timeset(routineN, handle)
      matrix_siesta(:) = 0.0_dp

      mepos = para_env%mepos
      nprocs = para_env%num_pe
      do_distribute = siesta_struct%gather_root < 0
      is_root_rank = siesta_struct%gather_root == mepos

      CALL dbcsr_get_info(matrix=matrix_dbcsr_kp(1)%matrix, &
                          nblkrows_total=nrows_blk, nblkcols_total=ncols_blk, &
                          row_blk_size=row_blk_size, col_blk_size=col_blk_size, &
                          row_blk_offset=row_blk_offset, col_blk_offset=col_blk_offset)
      symmetric = siesta_struct%symmetric

      ! number of locally stored SIESTA non-zero matrix elements
      n_nonzero_elements_siesta = SUM(siesta_struct%nelements_per_proc(:, nelements_dbcsr_recv))
      ! number of locally stored DBCSR non-zero matrix elements
      n_nonzero_elements_dbcsr = SUM(siesta_struct%nelements_per_proc(:, nelements_dbcsr_send))

      ! number of concurrent MPI isend / irecv operations
      nrequests_recv = get_number_of_mpi_sendrecv_requests(mepos, siesta_struct%nelements_per_proc(:, nelements_dbcsr_recv), &
                                                           max_mpi_packet_size_dp)
      nrequests_total = get_number_of_mpi_sendrecv_requests(mepos, siesta_struct%nelements_per_proc(:, nelements_dbcsr_send), &
                                                            max_mpi_packet_size_dp) + nrequests_recv

      IF (nrequests_total > 0) THEN
         ! allocate MPI-related arrays. request_tag is not actually needed, as MPI standard guarantees the order of
         ! peer-to-peer messages with the same tag between same processes
         ALLOCATE (requests(nrequests_total))
         ALLOCATE (peer_rank(nrequests_total), request_tag(nrequests_total))
         ALLOCATE (offset_per_request(nrequests_total), nelements_per_request(nrequests_total))
         !requests(:) = mp_request_null

         ! split large messages into a number of smaller messages. It is not really needed
         ! unless we are going to send > 2^31 matrix elements per MPI request
         IF (nrequests_recv > 0) THEN
            CALL assign_nonzero_elements_to_requests(offset_per_request(1:nrequests_recv), &
                                                     nelements_per_request(1:nrequests_recv), &
                                                     peer_rank(1:nrequests_recv), &
                                                     request_tag(1:nrequests_recv), &
                                                     mepos, &
                                                     siesta_struct%nelements_per_proc(:, nelements_dbcsr_recv), &
                                                     max_mpi_packet_size_dp)
         END IF
         IF (nrequests_total > nrequests_recv) THEN
            CALL assign_nonzero_elements_to_requests(offset_per_request(nrequests_recv + 1:nrequests_total), &
                                                     nelements_per_request(nrequests_recv + 1:nrequests_total), &
                                                     peer_rank(nrequests_recv + 1:nrequests_total), &
                                                     request_tag(nrequests_recv + 1:nrequests_total), &
                                                     mepos, &
                                                     siesta_struct%nelements_per_proc(:, nelements_dbcsr_send), &
                                                     max_mpi_packet_size_dp)
         END IF
      END IF

      ! point-to-point recv/send can be replaced with alltoallv, if data to distribute per rank is < 2^31 elements (should be OK due to SMEAGOL limitation)
      ! in principle, it is possible to overcome this limit by using derived datatypes (which will require additional wrapper functions, indeed).
      !
      ! pre-post non-blocking receive operations
      IF (n_nonzero_elements_siesta > 0) THEN
         ALLOCATE (recv_buffer(n_nonzero_elements_siesta))
      END IF
      DO irequest = 1, nrequests_recv
         CALL para_env%irecv(recv_buffer(offset_per_request(irequest) + 1: &
                                         offset_per_request(irequest) + nelements_per_request(irequest)), &
                             peer_rank(irequest), requests(irequest), request_tag(irequest))
      END DO

      ! pack local DBCSR non-zero matrix elements ordering by their target parallel process
      ALLOCATE (offset_per_proc(0:nprocs - 1), n_packed_elements_per_proc(0:nprocs - 1))
      offset_per_proc(0) = 0
      DO iproc = 1, nprocs - 1
         offset_per_proc(iproc) = offset_per_proc(iproc - 1) + siesta_struct%nelements_per_proc(iproc - 1, nelements_dbcsr_send)
      END DO
      n_packed_elements_per_proc(:) = 0

      ! number of local neighbour-list nodes and offset of the first local neighbour-list node
      nnodes_proc = siesta_struct%nnodes_per_proc(mepos)
      !node_offset = SUM(siesta_struct%nnodes_per_proc(0:mepos)) - siesta_struct%nnodes_per_proc(mepos)
      node_offset = SUM(siesta_struct%nnodes_per_proc(0:mepos)) - nnodes_proc

      ! if do_distribute == .FALSE., send all matrix elements to MPI process with rank gather_root
      ! in case of do_distribute == .TRUE., iproc is determined by calling WhichNodeOrb()
      iproc = siesta_struct%gather_root

      IF (n_nonzero_elements_dbcsr > 0) THEN
         ALLOCATE (send_buffer(n_nonzero_elements_dbcsr))
         send_buffer(:) = 0.0_dp

         ! iterate over locally-stored DBCSR matrix blocks.
         ! inode_proc is the target parallel process (where data are going to be sent)
         image_ind_offset = 0
         DO inode_proc = 1, nnodes_proc
            n_image_ind = siesta_struct%n_dbcsr_cell_images_to_merge(inode_proc)
            IF (n_image_ind > 0) THEN
               inode = node_offset + inode_proc

               irow_blk = siesta_struct%nl_repl(neighbor_list_iatom_index, inode)
               icol_blk = siesta_struct%nl_repl(neighbor_list_jatom_index, inode)
               CPASSERT(irow_blk <= icol_blk .OR. (.NOT. symmetric))
               image_siesta = siesta_struct%nl_repl(neighbor_list_siesta_image_index, inode)
               image_siesta_transp = siesta_struct%nl_repl(neighbor_list_siesta_transp_image_index, inode)

               nrows_local = row_blk_size(irow_blk)
               ncols_local = col_blk_size(icol_blk)
               first_row_minus_one = row_blk_offset(irow_blk) - 1
               first_col_minus_one = col_blk_offset(icol_blk) - 1

               ! merging cell images along transport direction
               IF (n_image_ind == 1) THEN
                  ! the most common case. Nothing to merge, so there is no need to allocate memory for a merged block
                  image_dbcsr = siesta_struct%dbcsr_cell_image_to_merge(image_ind_offset + 1)
                  CALL dbcsr_get_block_p(matrix=matrix_dbcsr_kp(image_dbcsr)%matrix, &
                                         row=irow_blk, col=icol_blk, block=sm_block_merged, found=found)
                  CPASSERT(found)
               ELSE ! n_image_ind > 1
                  ALLOCATE (sm_block_merged(nrows_local, ncols_local))

                  DO image_ind = 1, n_image_ind
                     image_dbcsr = siesta_struct%dbcsr_cell_image_to_merge(image_ind + image_ind_offset)

                     CALL dbcsr_get_block_p(matrix=matrix_dbcsr_kp(image_dbcsr)%matrix, &
                                            row=irow_blk, col=icol_blk, block=sm_block, found=found)
                     CPASSERT(found)
                     sm_block_merged(1:nrows_local, 1:ncols_local) = sm_block_merged(1:nrows_local, 1:ncols_local) + &
                                                                     sm_block(1:nrows_local, 1:ncols_local)
                  END DO
               END IF

               ! pack matrix elements for the 'normal' SIESTA matrix block
               IF (image_siesta > 0) THEN
                  DO irow_local = 1, nrows_local
                     IF (do_distribute) THEN
#if defined(__SMEAGOL)
                        CALL WhichNodeOrb(irow_local + first_row_minus_one, nprocs, iproc)
#else
                        CALL cp_abort(__LOCATION__, &
                                      "CP2K was compiled with no SMEAGOL support.")
#endif
                     END IF

                     ! CPASSERT
                     IF (debug_this_module) THEN
                        CPASSERT(iproc >= 0 .AND. iproc < nprocs)
                        IF (n_packed_elements_per_proc(iproc) + ncols_local > &
                            siesta_struct%nelements_per_proc(iproc, nelements_dbcsr_send)) THEN
                           CALL cp__a(__SHORT_FILE__, __LINE__)
                        END IF
                     END IF

                     offset_send_mepos = offset_per_proc(iproc) + n_packed_elements_per_proc(iproc)
                     send_buffer(offset_send_mepos + 1:offset_send_mepos + ncols_local) = sm_block_merged(irow_local, 1:ncols_local)

                     n_packed_elements_per_proc(iproc) = n_packed_elements_per_proc(iproc) + ncols_local
                  END DO
               END IF

               ! pack matrix elements of the transposed SIESTA matrix block
               IF (image_siesta_transp > 0) THEN
                  DO icol_local = 1, ncols_local
                     IF (do_distribute) THEN
#if defined(__SMEAGOL)
                        CALL WhichNodeOrb(icol_local + first_col_minus_one, nprocs, iproc) ! iproc_orb
#else
                        CALL cp_abort(__LOCATION__, &
                                      "CP2K was compiled with no SMEAGOL support.")
#endif
                     END IF

                     ! CPASSERT
                     IF (debug_this_module) THEN
                        CPASSERT(iproc >= 0 .AND. iproc < nprocs)
                        IF (n_packed_elements_per_proc(iproc) + nrows_local > &
                            siesta_struct%nelements_per_proc(iproc, nelements_dbcsr_send)) THEN
                           CALL cp__a(__SHORT_FILE__, __LINE__)
                        END IF
                     END IF

                     offset_send_mepos = offset_per_proc(iproc) + n_packed_elements_per_proc(iproc)
                     send_buffer(offset_send_mepos + 1:offset_send_mepos + nrows_local) = sm_block_merged(1:nrows_local, icol_local)

                     n_packed_elements_per_proc(iproc) = n_packed_elements_per_proc(iproc) + nrows_local
                  END DO
               END IF

               IF (n_image_ind > 1) THEN
                  DEALLOCATE (sm_block_merged)
               END IF
            END IF

            image_ind_offset = image_ind_offset + siesta_struct%n_dbcsr_cell_images_to_merge(inode_proc)
         END DO

         IF (debug_this_module) THEN
            DO iproc = 0, nprocs - 1
               IF (n_packed_elements_per_proc(iproc) /= siesta_struct%nelements_per_proc(iproc, nelements_dbcsr_send)) THEN
                  CALL cp__a(__SHORT_FILE__, __LINE__)
               END IF
            END DO
         END IF

         ! send packed data to other parallel processes
         DO irequest = nrequests_recv + 1, nrequests_total
            CALL para_env%isend(send_buffer(offset_per_request(irequest) + 1: &
                                            offset_per_request(irequest) + nelements_per_request(irequest)), &
                                peer_rank(irequest), requests(irequest), request_tag(irequest))
         END DO

         ! copy data locally that stay on the same process.
         IF (mepos > 0) THEN
            offset_recv_mepos = SUM(siesta_struct%nelements_per_proc(0:mepos - 1, nelements_dbcsr_recv))
         ELSE
            offset_recv_mepos = 0
         END IF
         offset_send_mepos = offset_per_proc(mepos)

         IF (debug_this_module) THEN
            IF (n_packed_elements_per_proc(mepos) /= siesta_struct%nelements_per_proc(mepos, nelements_dbcsr_recv)) THEN
               CALL cp__a(__SHORT_FILE__, __LINE__)
            END IF
         END IF

         IF (n_packed_elements_per_proc(mepos) > 0) THEN
            recv_buffer(offset_recv_mepos + 1:offset_recv_mepos + n_packed_elements_per_proc(mepos)) = &
               send_buffer(offset_send_mepos + 1:offset_send_mepos + n_packed_elements_per_proc(mepos))
         END IF
      END IF

      IF (nrequests_total > 0) THEN
         ! wait for pending isend/irecv requests
         CALL mp_waitall(requests)
         DEALLOCATE (nelements_per_request, offset_per_request, peer_rank, requests, request_tag)
      END IF

      ! release send buffers
      IF (ALLOCATED(send_buffer)) DEALLOCATE (send_buffer)
      DEALLOCATE (offset_per_proc, n_packed_elements_per_proc)

      ! non-zero matrix elements in 'recv_buffer' array are grouped by their source MPI rank, local row index, and column index (in this order).
      ! Reorder the matrix elements ('reorder_recv_buffer') so they are grouped by their local row index, source MPI rank, and column index.
      ! (column indices are in the ascending order within each (row index, source MPI rank) block).
      ! The array 'packed_index' allows mapping matrix element between these intermediate order and SIESTA order : local row index, column index
      IF (n_nonzero_elements_siesta > 0) THEN
         ALLOCATE (reorder_recv_buffer(n_nonzero_elements_siesta))
         ALLOCATE (next_nonzero_element_offset(siesta_struct%nrows))
         next_nonzero_element_offset(:) = 0
         offset_recv_mepos = 0

         DO inode = 1, SIZE(siesta_struct%nl_repl, 2)
            irow_blk = siesta_struct%nl_repl(neighbor_list_iatom_index, inode)
            icol_blk = siesta_struct%nl_repl(neighbor_list_jatom_index, inode)
            CPASSERT(irow_blk <= icol_blk .OR. (.NOT. symmetric))
            image_siesta = siesta_struct%nl_repl(neighbor_list_siesta_image_index, inode)
            image_siesta_transp = siesta_struct%nl_repl(neighbor_list_siesta_transp_image_index, inode)

            nrows_local = row_blk_size(irow_blk)
            ncols_local = col_blk_size(icol_blk)
            first_row_minus_one = row_blk_offset(irow_blk) - 1
            first_col_minus_one = col_blk_offset(icol_blk) - 1

            ! normal block
            IF (image_siesta > 0) THEN
               DO irow_local = 1, nrows_local
                  IF (do_distribute) THEN
#if defined(__SMEAGOL)
                     CALL GlobalToLocalOrb(irow_local + first_row_minus_one, mepos, nprocs, irow_proc)
#else
                     CALL cp_abort(__LOCATION__, &
                                   "CP2K was compiled with no SMEAGOL support.")
#endif
                  ELSE
                     IF (is_root_rank) THEN
                        irow_proc = irow_local + first_row_minus_one
                     ELSE
                        irow_proc = 0
                     END IF
                  END IF
                  IF (irow_proc > 0) THEN
                     offset_send_mepos = siesta_struct%row_offset(irow_proc) + next_nonzero_element_offset(irow_proc)
                     reorder_recv_buffer(offset_send_mepos + 1:offset_send_mepos + ncols_local) = &
                        recv_buffer(offset_recv_mepos + 1:offset_recv_mepos + ncols_local)
                     offset_recv_mepos = offset_recv_mepos + ncols_local
                     next_nonzero_element_offset(irow_proc) = next_nonzero_element_offset(irow_proc) + ncols_local
                  END IF
               END DO
            END IF

            ! transposed block
            IF (image_siesta_transp > 0) THEN
               DO icol_local = 1, ncols_local
                  IF (do_distribute) THEN
#if defined(__SMEAGOL)
                     CALL GlobalToLocalOrb(icol_local + first_col_minus_one, mepos, nprocs, irow_proc)
#else
                     CALL cp_abort(__LOCATION__, &
                                   "CP2K was compiled with no SMEAGOL support.")
#endif
                  ELSE
                     IF (is_root_rank) THEN
                        irow_proc = icol_local + first_col_minus_one
                     ELSE
                        irow_proc = 0
                     END IF
                  END IF
                  IF (irow_proc > 0) THEN
                     offset_send_mepos = siesta_struct%row_offset(irow_proc) + next_nonzero_element_offset(irow_proc)
                     reorder_recv_buffer(offset_send_mepos + 1:offset_send_mepos + nrows_local) = &
                        recv_buffer(offset_recv_mepos + 1:offset_recv_mepos + nrows_local)
                     offset_recv_mepos = offset_recv_mepos + nrows_local
                     next_nonzero_element_offset(irow_proc) = next_nonzero_element_offset(irow_proc) + nrows_local
                  END IF
               END DO
            END IF
         END DO

         IF (debug_this_module) THEN
            DO irow_local = 1, siesta_struct%nrows
               IF (siesta_struct%n_nonzero_cols(irow_local) /= next_nonzero_element_offset(irow_local)) THEN
                  CALL cp__a(__SHORT_FILE__, __LINE__)
               END IF
            END DO
         END IF

         DEALLOCATE (next_nonzero_element_offset)
         DEALLOCATE (recv_buffer)

         ! Map non-zero matrix element between the intermediate order and SIESTA order
         DO irow_local = 1, siesta_struct%nrows
            offset_recv_mepos = siesta_struct%row_offset(irow_local)
            DO icol_local = 1, siesta_struct%n_nonzero_cols(irow_local)
               matrix_siesta(offset_recv_mepos + icol_local) = &
                  reorder_recv_buffer(offset_recv_mepos + siesta_struct%packed_index(offset_recv_mepos + icol_local))
            END DO
         END DO
         DEALLOCATE (reorder_recv_buffer)
      END IF

      CALL timestop(handle)
   END SUBROUTINE convert_dbcsr_to_distributed_siesta

! **************************************************************************************************
!> \brief Convert matrix from DBCSR to sparse SIESTA format.
!> \param matrix_dbcsr_kp    DBCSR matrix [out]. The matrix is declared as INTENT(in) as pointers to
!>                 dbcsr matrices remain intact. However we have intention to update matrix elements
!> \param matrix_siesta      matrix in SIESTA format [in]
!> \param siesta_struct      structure to map matrix blocks between formats
!> \param para_env           MPI parallel environment
! **************************************************************************************************
   SUBROUTINE convert_distributed_siesta_to_dbcsr(matrix_dbcsr_kp, matrix_siesta, siesta_struct, para_env)
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(in)       :: matrix_dbcsr_kp
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: matrix_siesta
      TYPE(siesta_distrib_csc_struct_type), INTENT(in)   :: siesta_struct
      TYPE(mp_para_env_type), INTENT(in), POINTER        :: para_env

      CHARACTER(len=*), PARAMETER :: routineN = 'convert_distributed_siesta_to_dbcsr'

      INTEGER :: first_col_minus_one, first_row_minus_one, handle, icol_blk, icol_local, &
         image_dbcsr, image_siesta, image_siesta_transp, inode, inode_proc, iproc, irequest, &
         irow_blk, irow_local, irow_proc, mepos, n_image_ind, ncols_blk, ncols_local, nnodes_proc, &
         node_offset, nprocs, nrequests_recv, nrequests_total, nrows_blk, nrows_local
      INTEGER(kind=int_8)                                :: n_nonzero_elements_dbcsr, &
                                                            n_nonzero_elements_siesta, &
                                                            offset_recv_mepos, offset_send_mepos
      INTEGER(kind=int_8), ALLOCATABLE, DIMENSION(:)     :: n_packed_elements_per_proc, &
                                                            nelements_per_request, &
                                                            offset_per_proc, offset_per_request
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: next_nonzero_element_offset, peer_rank, &
                                                            request_tag
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_offset, col_blk_size, &
                                                            row_blk_offset, row_blk_size
      LOGICAL                                            :: do_distribute, found, is_root_rank, &
                                                            symmetric
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: recv_buffer, reorder_send_buffer, &
                                                            send_buffer
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: sm_block
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: requests

      CALL timeset(routineN, handle)
      DO image_dbcsr = 1, SIZE(matrix_dbcsr_kp)
         CALL dbcsr_set(matrix_dbcsr_kp(image_dbcsr)%matrix, 0.0_dp)
      END DO

      mepos = para_env%mepos
      nprocs = para_env%num_pe
      do_distribute = siesta_struct%gather_root < 0
      is_root_rank = siesta_struct%gather_root == mepos

      CALL dbcsr_get_info(matrix=matrix_dbcsr_kp(1)%matrix, &
                          nblkrows_total=nrows_blk, nblkcols_total=ncols_blk, &
                          row_blk_size=row_blk_size, col_blk_size=col_blk_size, &
                          row_blk_offset=row_blk_offset, col_blk_offset=col_blk_offset)
      symmetric = siesta_struct%symmetric

      n_nonzero_elements_siesta = SUM(siesta_struct%nelements_per_proc(:, nelements_dbcsr_recv))
      n_nonzero_elements_dbcsr = SUM(siesta_struct%nelements_per_proc(:, nelements_dbcsr_send))

      nrequests_recv = get_number_of_mpi_sendrecv_requests(mepos, siesta_struct%nelements_per_proc(:, nelements_dbcsr_send), &
                                                           max_mpi_packet_size_dp)
      nrequests_total = get_number_of_mpi_sendrecv_requests(mepos, siesta_struct%nelements_per_proc(:, nelements_dbcsr_recv), &
                                                            max_mpi_packet_size_dp) + nrequests_recv
      IF (nrequests_total > 0) THEN
         ALLOCATE (requests(nrequests_total))
         ALLOCATE (peer_rank(nrequests_total), request_tag(nrequests_total))
         ALLOCATE (offset_per_request(nrequests_total), nelements_per_request(nrequests_total))
         !requests(:) = mp_request_null
         IF (nrequests_recv > 0) THEN
            CALL assign_nonzero_elements_to_requests(offset_per_request(1:nrequests_recv), &
                                                     nelements_per_request(1:nrequests_recv), &
                                                     peer_rank(1:nrequests_recv), &
                                                     request_tag(1:nrequests_recv), &
                                                     mepos, &
                                                     siesta_struct%nelements_per_proc(:, nelements_dbcsr_send), &
                                                     max_mpi_packet_size_dp)
         END IF
         IF (nrequests_total > nrequests_recv) THEN
            CALL assign_nonzero_elements_to_requests(offset_per_request(nrequests_recv + 1:nrequests_total), &
                                                     nelements_per_request(nrequests_recv + 1:nrequests_total), &
                                                     peer_rank(nrequests_recv + 1:nrequests_total), &
                                                     request_tag(nrequests_recv + 1:nrequests_total), &
                                                     mepos, &
                                                     siesta_struct%nelements_per_proc(:, nelements_dbcsr_recv), &
                                                     max_mpi_packet_size_dp)
         END IF
      END IF

      IF (n_nonzero_elements_dbcsr > 0) THEN
         ALLOCATE (recv_buffer(n_nonzero_elements_dbcsr))
      END IF
      DO irequest = 1, nrequests_recv
         CALL para_env%irecv(recv_buffer(offset_per_request(irequest) + 1: &
                                         offset_per_request(irequest) + nelements_per_request(irequest)), &
                             peer_rank(irequest), requests(irequest), request_tag(irequest))
      END DO

      ALLOCATE (offset_per_proc(0:nprocs - 1), n_packed_elements_per_proc(0:nprocs - 1))
      offset_per_proc(0) = 0
      DO iproc = 1, nprocs - 1
         offset_per_proc(iproc) = offset_per_proc(iproc - 1) + siesta_struct%nelements_per_proc(iproc - 1, nelements_dbcsr_send)
      END DO
      n_packed_elements_per_proc(:) = 0

      IF (mepos > 0) THEN
         node_offset = SUM(siesta_struct%nnodes_per_proc(0:mepos - 1))
      ELSE
         node_offset = 0
      END IF
      nnodes_proc = siesta_struct%nnodes_per_proc(mepos)

      IF (n_nonzero_elements_siesta > 0) THEN
         ALLOCATE (send_buffer(n_nonzero_elements_siesta))

         ALLOCATE (reorder_send_buffer(n_nonzero_elements_siesta))
         DO irow_local = 1, siesta_struct%nrows
            offset_send_mepos = siesta_struct%row_offset(irow_local)
            DO icol_local = 1, siesta_struct%n_nonzero_cols(irow_local)
               reorder_send_buffer(offset_send_mepos + siesta_struct%packed_index(offset_send_mepos + icol_local)) = &
                  matrix_siesta(offset_send_mepos + icol_local)
            END DO
         END DO

         ALLOCATE (next_nonzero_element_offset(siesta_struct%nrows))
         next_nonzero_element_offset(:) = 0
         offset_send_mepos = 0

         DO inode = 1, SIZE(siesta_struct%nl_repl, 2)
            irow_blk = siesta_struct%nl_repl(neighbor_list_iatom_index, inode)
            icol_blk = siesta_struct%nl_repl(neighbor_list_jatom_index, inode)
            CPASSERT(irow_blk <= icol_blk .OR. (.NOT. symmetric))
            image_siesta = siesta_struct%nl_repl(neighbor_list_siesta_image_index, inode)
            image_siesta_transp = siesta_struct%nl_repl(neighbor_list_siesta_transp_image_index, inode)

            nrows_local = row_blk_size(irow_blk)
            ncols_local = col_blk_size(icol_blk)
            first_row_minus_one = row_blk_offset(irow_blk) - 1
            first_col_minus_one = col_blk_offset(icol_blk) - 1

            IF (image_siesta > 0) THEN
               DO irow_local = 1, nrows_local
                  IF (do_distribute) THEN
#if defined(__SMEAGOL)
                     CALL GlobalToLocalOrb(irow_local + first_row_minus_one, mepos, nprocs, irow_proc)
#else
                     CALL cp_abort(__LOCATION__, &
                                   "CP2K was compiled with no SMEAGOL support.")
#endif
                  ELSE
                     IF (is_root_rank) THEN
                        irow_proc = irow_local + first_row_minus_one
                     ELSE
                        irow_proc = 0
                     END IF
                  END IF
                  IF (irow_proc > 0) THEN
                     offset_recv_mepos = siesta_struct%row_offset(irow_proc) + next_nonzero_element_offset(irow_proc)
                     send_buffer(offset_send_mepos + 1:offset_send_mepos + ncols_local) = &
                        reorder_send_buffer(offset_recv_mepos + 1:offset_recv_mepos + ncols_local)
                     offset_send_mepos = offset_send_mepos + ncols_local
                     next_nonzero_element_offset(irow_proc) = next_nonzero_element_offset(irow_proc) + ncols_local
                  END IF
               END DO
            END IF

            ! transposed block
            IF (image_siesta_transp > 0) THEN
               DO icol_local = 1, ncols_local
                  IF (do_distribute) THEN
#if defined(__SMEAGOL)
                     CALL GlobalToLocalOrb(icol_local + first_col_minus_one, mepos, nprocs, irow_proc)
#else
                     CALL cp_abort(__LOCATION__, &
                                   "CP2K was compiled with no SMEAGOL support.")
#endif
                  ELSE
                     IF (is_root_rank) THEN
                        irow_proc = icol_local + first_col_minus_one
                     ELSE
                        irow_proc = 0
                     END IF
                  END IF
                  IF (irow_proc > 0) THEN
                     offset_recv_mepos = siesta_struct%row_offset(irow_proc) + next_nonzero_element_offset(irow_proc)
                     send_buffer(offset_send_mepos + 1:offset_send_mepos + nrows_local) = &
                        reorder_send_buffer(offset_recv_mepos + 1:offset_recv_mepos + nrows_local)
                     offset_send_mepos = offset_send_mepos + nrows_local
                     next_nonzero_element_offset(irow_proc) = next_nonzero_element_offset(irow_proc) + nrows_local
                  END IF
               END DO
            END IF
         END DO

         IF (debug_this_module) THEN
            DO irow_local = 1, siesta_struct%nrows
               IF (siesta_struct%n_nonzero_cols(irow_local) /= next_nonzero_element_offset(irow_local)) THEN
                  CALL cp__a(__SHORT_FILE__, __LINE__)
               END IF
            END DO
         END IF

         DEALLOCATE (next_nonzero_element_offset)
         DEALLOCATE (reorder_send_buffer)

         DO irequest = nrequests_recv + 1, nrequests_total
            CALL para_env%isend(send_buffer(offset_per_request(irequest) + 1: &
                                            offset_per_request(irequest) + nelements_per_request(irequest)), &
                                peer_rank(irequest), requests(irequest), request_tag(irequest))
         END DO

         ! copy data locally that stay on the same process.
         IF (mepos > 0) THEN
            offset_send_mepos = SUM(siesta_struct%nelements_per_proc(0:mepos - 1, nelements_dbcsr_recv))
         ELSE
            offset_send_mepos = 0
         END IF
         offset_recv_mepos = offset_per_proc(mepos)

         IF (debug_this_module) THEN
            IF (siesta_struct%nelements_per_proc(mepos, nelements_dbcsr_recv) /= &
                siesta_struct%nelements_per_proc(mepos, nelements_dbcsr_send)) THEN
               CALL cp__a(__SHORT_FILE__, __LINE__)
            END IF
         END IF

         IF (siesta_struct%nelements_per_proc(mepos, nelements_dbcsr_send) > 0) THEN
            recv_buffer(offset_recv_mepos + 1:offset_recv_mepos + siesta_struct%nelements_per_proc(mepos, nelements_dbcsr_send)) = &
               send_buffer(offset_send_mepos + 1:offset_send_mepos + siesta_struct%nelements_per_proc(mepos, nelements_dbcsr_send))
         END IF
      END IF

      IF (nrequests_total > 0) THEN
         ! wait for pending isend/irecv requests
         CALL mp_waitall(requests)
         DEALLOCATE (nelements_per_request, offset_per_request, peer_rank, requests, request_tag)
      END IF

      IF (ALLOCATED(send_buffer)) DEALLOCATE (send_buffer)

      iproc = siesta_struct%gather_root ! if do_distribute == .FALSE., collect matrix elements from MPI process with rank gather_root
      IF (n_nonzero_elements_dbcsr > 0) THEN
         DO inode_proc = 1, nnodes_proc
            n_image_ind = siesta_struct%n_dbcsr_cell_images_to_merge(inode_proc)
            IF (n_image_ind > 0) THEN
               inode = node_offset + inode_proc

               irow_blk = siesta_struct%nl_repl(neighbor_list_iatom_index, inode)
               icol_blk = siesta_struct%nl_repl(neighbor_list_jatom_index, inode)
               image_dbcsr = siesta_struct%nl_repl(neighbor_list_dbcsr_image_index, inode)
               CPASSERT(irow_blk <= icol_blk .OR. (.NOT. symmetric))
               image_siesta = siesta_struct%nl_repl(neighbor_list_siesta_image_index, inode)
               image_siesta_transp = siesta_struct%nl_repl(neighbor_list_siesta_transp_image_index, inode)

               nrows_local = row_blk_size(irow_blk)
               ncols_local = col_blk_size(icol_blk)
               first_row_minus_one = row_blk_offset(irow_blk) - 1
               first_col_minus_one = col_blk_offset(icol_blk) - 1

               CALL dbcsr_get_block_p(matrix=matrix_dbcsr_kp(image_dbcsr)%matrix, &
                                      row=irow_blk, col=icol_blk, block=sm_block, found=found)
               CPASSERT(found)

               IF (image_siesta > 0) THEN
                  DO irow_local = 1, nrows_local
                     IF (do_distribute) THEN
#if defined(__SMEAGOL)
                        CALL WhichNodeOrb(irow_local + first_row_minus_one, nprocs, iproc) ! iproc_orb
#else
                        CALL cp_abort(__LOCATION__, &
                                      "CP2K was compiled with no SMEAGOL support.")
#endif
                     END IF
                     ! CPASSERT
                     IF (debug_this_module) THEN
                        CPASSERT(iproc >= 0 .AND. iproc < nprocs)
                        IF (n_packed_elements_per_proc(iproc) + ncols_local > &
                            siesta_struct%nelements_per_proc(iproc, nelements_dbcsr_send)) THEN
                           CALL cp__a(__SHORT_FILE__, __LINE__)
                        END IF
                     END IF

                     offset_recv_mepos = offset_per_proc(iproc) + n_packed_elements_per_proc(iproc)
                     sm_block(irow_local, 1:ncols_local) = recv_buffer(offset_recv_mepos + 1:offset_recv_mepos + ncols_local)

                     n_packed_elements_per_proc(iproc) = n_packed_elements_per_proc(iproc) + ncols_local
                  END DO
               END IF

               ! transposed block
               IF (image_siesta_transp > 0) THEN
                  DO icol_local = 1, ncols_local
                     IF (do_distribute) THEN
#if defined(__SMEAGOL)
                        CALL WhichNodeOrb(icol_local + first_col_minus_one, nprocs, iproc) ! iproc_orb
#else
                        CALL cp_abort(__LOCATION__, &
                                      "CP2K was compiled with no SMEAGOL support.")
#endif
                     END IF
                     ! CPASSERT
                     IF (debug_this_module) THEN
                        CPASSERT(iproc >= 0 .AND. iproc < nprocs)
                        IF (n_packed_elements_per_proc(iproc) + nrows_local > &
                            siesta_struct%nelements_per_proc(iproc, nelements_dbcsr_send)) THEN
                           CALL cp__a(__SHORT_FILE__, __LINE__)
                        END IF
                     END IF

                     offset_recv_mepos = offset_per_proc(iproc) + n_packed_elements_per_proc(iproc)
                     sm_block(1:nrows_local, icol_local) = recv_buffer(offset_recv_mepos + 1:offset_recv_mepos + nrows_local)

                     n_packed_elements_per_proc(iproc) = n_packed_elements_per_proc(iproc) + nrows_local
                  END DO
               END IF
            END IF
         END DO

         IF (debug_this_module) THEN
            DO iproc = 0, nprocs - 1
               IF (n_packed_elements_per_proc(iproc) /= siesta_struct%nelements_per_proc(iproc, nelements_dbcsr_send)) THEN
                  CALL cp__a(__SHORT_FILE__, __LINE__)
               END IF
            END DO
         END IF

         DEALLOCATE (recv_buffer)
      END IF

      DEALLOCATE (offset_per_proc, n_packed_elements_per_proc)

      CALL timestop(handle)
   END SUBROUTINE convert_distributed_siesta_to_dbcsr

   ! *** PRIVATE SUBROUTINES ***

! **************************************************************************************************
!> \brief Computes number of neighbour-list nodes on the current parallel process.
!> \param nnodes_local              number of nodes [out]
!> \param max_ijk_cell_image_local  largest index of cell images along i, j and k cell vectors
!>                                  on this parallel process [out]
!> \param max_ijk_cell_image        largest index of cell images along i, j and k cell vectors [inout]
!> \param sab_nl                    pair-wise neighbour list [in]
!> \param para_env                  MPI parallel environment [in]
!> \param particle_coords           list of atomic coordinates subject to periodic boundary conditions [in]
!> \param cell                      simulation unit cell [in]
! **************************************************************************************************
   SUBROUTINE get_nnodes_local(nnodes_local, max_ijk_cell_image_local, max_ijk_cell_image, sab_nl, para_env, particle_coords, cell)
      INTEGER, INTENT(out)                               :: nnodes_local
      INTEGER, DIMENSION(3), INTENT(out)                 :: max_ijk_cell_image_local
      INTEGER, DIMENSION(3), INTENT(inout)               :: max_ijk_cell_image
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(in), POINTER                             :: sab_nl
      TYPE(mp_para_env_type), INTENT(in), POINTER        :: para_env
      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :), &
         INTENT(in)                                      :: particle_coords
      TYPE(cell_type), INTENT(in), POINTER               :: cell

      CHARACTER(len=*), PARAMETER                        :: routineN = 'get_nnodes_local'

      INTEGER                                            :: handle, iatom, icoord, jatom
      INTEGER, DIMENSION(3)                              :: cell_ijk, max_ijk_cell_image_tmp
      LOGICAL                                            :: update_ncells
      REAL(kind=dp), DIMENSION(3)                        :: r_ij
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      CALL timeset(routineN, handle)

      update_ncells = .FALSE.
      DO icoord = 1, 3 ! x, y, z
         IF (max_ijk_cell_image(icoord) >= 0) THEN
            max_ijk_cell_image_tmp(icoord) = max_ijk_cell_image(icoord)
         ELSE
            max_ijk_cell_image_tmp(icoord) = HUGE(max_ijk_cell_image_tmp(icoord))
            update_ncells = .TRUE.
         END IF
      END DO

      nnodes_local = 0
      max_ijk_cell_image_local(:) = 0

      CALL neighbor_list_iterator_create(nl_iterator, sab_nl)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, r=r_ij)
         CALL get_negf_cell_ijk(cell_ijk, r_ij, r_i=particle_coords(1:3, iatom), r_j=particle_coords(1:3, jatom), cell=cell)
         cell_ijk(1:3) = ABS(cell_ijk(1:3))

         IF (cell_ijk(1) <= max_ijk_cell_image_tmp(1) .AND. cell_ijk(2) <= max_ijk_cell_image_tmp(2) .AND. &
             cell_ijk(3) <= max_ijk_cell_image_tmp(3)) THEN
            nnodes_local = nnodes_local + 1
            max_ijk_cell_image_local(1:3) = MAX(max_ijk_cell_image_local(1:3), cell_ijk(1:3))
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      IF (update_ncells) THEN
         max_ijk_cell_image_tmp(1:3) = max_ijk_cell_image_local(1:3)
         CALL para_env%max(max_ijk_cell_image_tmp)
         DO icoord = 1, 3
            IF (max_ijk_cell_image(icoord) < 0) THEN
               max_ijk_cell_image(icoord) = max_ijk_cell_image_tmp(icoord)
            END IF
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE get_nnodes_local

! **************************************************************************************************
!> \brief Construct list of neighbour-list's nodes on the current parallel process.
!> \param nl_local                  non-merged local neighbour-list's nodes [out]
!> \param max_ijk_cell_image_local  largest index of cell images along i, j and k cell vectors
!>                                  on this parallel process [in]
!> \param max_ijk_cell_image        largest index of cell images along i, j and k cell vectors [in]
!> \param sab_nl                    pair-wise neighbour list [in]
!> \param particle_coords           list of atomic coordinates subject to periodic boundary conditions [in]
!> \param cell                      simulation unit cell [in]
!> \param cell_to_index             array to convert 3-D cell indices to 1-D DBCSR image indices
!> \param do_merge                  merge cell images along transport direction [in]
!> \param node_merged_indices       nodes-related indices. Nodes with identical indices will be merged [out]
!> \param k_cells                   list of cell image indices along transport direction. Nodes to
!>                be merged have the same 'node_merged_indices' but different 'k_cells' indices [out]
! **************************************************************************************************
   SUBROUTINE get_nl_nodes_local(nl_local, max_ijk_cell_image_local, max_ijk_cell_image, sab_nl, particle_coords, &
                                 cell, cell_to_index, do_merge, node_merged_indices, k_cells)
      INTEGER, DIMENSION(:, :), INTENT(out)              :: nl_local
      INTEGER, DIMENSION(3), INTENT(in)                  :: max_ijk_cell_image_local, &
                                                            max_ijk_cell_image
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(in), POINTER                             :: sab_nl
      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :), &
         INTENT(in)                                      :: particle_coords
      TYPE(cell_type), INTENT(in), POINTER               :: cell
      INTEGER, DIMENSION(:, :, :), INTENT(in), POINTER   :: cell_to_index
      LOGICAL, INTENT(in)                                :: do_merge
      INTEGER(kind=int_8), DIMENSION(:), INTENT(out)     :: node_merged_indices
      INTEGER, DIMENSION(:), INTENT(out)                 :: k_cells

      CHARACTER(len=*), PARAMETER :: routineN = 'get_nl_nodes_local'

      INTEGER                                            :: handle, iatom, icol_blk, image, inode, &
                                                            irow_blk, jatom, natoms
      INTEGER(kind=8), DIMENSION(2)                      :: ncells_siesta_local
      INTEGER, DIMENSION(2)                              :: ncells_siesta
      INTEGER, DIMENSION(3)                              :: cell_ijk_abs, cell_ijk_dbcsr, &
                                                            cell_ijk_siesta
      LOGICAL                                            :: do_symmetric
      REAL(kind=dp), DIMENSION(3)                        :: r_ij
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

      CALL timeset(routineN, handle)
      ! natoms only used to compute a merged 1D index of each DBCSR block
      natoms = SIZE(particle_coords, 2)
      CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl, symmetric=do_symmetric)
      ncells_siesta(1:2) = 2*max_ijk_cell_image(1:2) + 1

      ncells_siesta_local(1:2) = INT(2*max_ijk_cell_image_local(1:2) + 1, kind=int_8)

      inode = 0
      CALL neighbor_list_iterator_create(nl_iterator, sab_nl)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, cell=cell_ijk_dbcsr, r=r_ij)
         CALL get_negf_cell_ijk(cell_ijk_abs, r_ij, r_i=particle_coords(1:3, iatom), r_j=particle_coords(1:3, jatom), cell=cell)

         IF (ABS(cell_ijk_abs(1)) <= max_ijk_cell_image(1) .AND. ABS(cell_ijk_abs(2)) <= max_ijk_cell_image(2) .AND. &
             ABS(cell_ijk_abs(3)) <= max_ijk_cell_image(3)) THEN

            inode = inode + 1

            image = get_index_by_cell(cell_ijk_dbcsr, cell_to_index)
            CPASSERT(image > 0)
            nl_local(neighbor_list_dbcsr_image_index, inode) = image

            IF (do_symmetric .AND. iatom > jatom) THEN
               irow_blk = jatom
               icol_blk = iatom
               cell_ijk_abs(1:3) = -cell_ijk_abs(1:3)
            ELSE
               irow_blk = iatom
               icol_blk = jatom
            END IF

            nl_local(neighbor_list_iatom_index, inode) = irow_blk
            nl_local(neighbor_list_jatom_index, inode) = icol_blk

            cell_ijk_siesta(1:3) = index_in_canonical_enumeration(cell_ijk_abs(1:3)) ! absolute -> SIESTA

            IF (do_merge) THEN
               node_merged_indices(inode) = (((cell_ijk_siesta(2) - 1)*ncells_siesta_local(1) + cell_ijk_siesta(1) - 1)* &
                                     INT(natoms, kind=int_8) + icol_blk - 1)*INT(natoms, kind=int_8) + INT(irow_blk - 1, kind=int_8)
               image = cell_ijk_siesta(1) + ncells_siesta(1)*(cell_ijk_siesta(2) - 1)
            ELSE
               node_merged_indices(inode) = ((((cell_ijk_siesta(3) - 1)*ncells_siesta_local(2) + &
                                               cell_ijk_siesta(2) - 1)*ncells_siesta_local(1) + cell_ijk_siesta(1) - 1)* &
                                     INT(natoms, kind=int_8) + icol_blk - 1)*INT(natoms, kind=int_8) + INT(irow_blk - 1, kind=int_8)
               image = cell_ijk_siesta(1) + ncells_siesta(1)*(cell_ijk_siesta(2) - 1 + ncells_siesta(2)*(cell_ijk_siesta(3) - 1))
            END IF
            k_cells(inode) = cell_ijk_siesta(3)
            nl_local(neighbor_list_siesta_image_index, inode) = image

            IF (do_symmetric .AND. irow_blk /= icol_blk) THEN
               cell_ijk_abs(1:3) = -cell_ijk_abs(1:3)
               cell_ijk_siesta(1:3) = index_in_canonical_enumeration(cell_ijk_abs(1:3)) ! absolute -> SIESTA
               IF (do_merge) cell_ijk_siesta(3) = 1
               nl_local(neighbor_list_siesta_transp_image_index, inode) = &
                  cell_ijk_siesta(1) + ncells_siesta(1)*(cell_ijk_siesta(2) - 1 + ncells_siesta(2)*(cell_ijk_siesta(3) - 1))
            ELSE
               nl_local(neighbor_list_siesta_transp_image_index, inode) = 0
            END IF
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      IF (debug_this_module) THEN
         CPASSERT(SIZE(nl_local, 2) == inode)
      END IF

      CALL timestop(handle)
   END SUBROUTINE get_nl_nodes_local

! **************************************************************************************************
!> \brief Replicate (and optionally merge) pair-wise neighbour list.
!> \param repl_nl                   replicated neighbour list. It needs to be deallocated elsewhere [allocated]
!> \param n_dbcsr_cell_images_to_merge  number of merged blocks per neighbour-list node [allocated]
!> \param dbcsr_cell_image_to_merge list of DBCSR image indices to merge [allocated]
!> \param nnodes_per_proc           number of merged nodes on each parallel processes [out]
!> \param max_ijk_cell_image        largest index of cell images along i, j and k cell vectors [inout]
!> \param sab_nl                    pair-wise neighbour list [in]
!> \param para_env                  MPI parallel environment [in]
!> \param particle_coords           list of atomic coordinates subject to periodic boundary conditions [in]
!> \param cell                      simulation unit cell [in]
!> \param cell_to_index             array to convert 3-D cell indices to 1-D DBCSR image indices [in]
!> \param do_merge                  merge cell images along transport direction [in]
! **************************************************************************************************
   SUBROUTINE replicate_neighbour_list(repl_nl, n_dbcsr_cell_images_to_merge, dbcsr_cell_image_to_merge, &
                                       nnodes_per_proc, max_ijk_cell_image, sab_nl, para_env, particle_coords, &
                                       cell, cell_to_index, do_merge)
      INTEGER, ALLOCATABLE, DIMENSION(:, :), &
         INTENT(inout)                                   :: repl_nl
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(inout)  :: n_dbcsr_cell_images_to_merge, &
                                                            dbcsr_cell_image_to_merge
      INTEGER, DIMENSION(0:), INTENT(out)                :: nnodes_per_proc
      INTEGER, DIMENSION(3), INTENT(inout)               :: max_ijk_cell_image
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         INTENT(in), POINTER                             :: sab_nl
      TYPE(mp_para_env_type), INTENT(in), POINTER        :: para_env
      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :), &
         INTENT(in)                                      :: particle_coords
      TYPE(cell_type), INTENT(in), POINTER               :: cell
      INTEGER, DIMENSION(:, :, :), INTENT(in), POINTER   :: cell_to_index
      LOGICAL, INTENT(in)                                :: do_merge

      CHARACTER(len=*), PARAMETER :: routineN = 'replicate_neighbour_list'

      INTEGER                                            :: handle, inode, iproc, kcell_closest, &
                                                            nnodes_local, nnodes_merged, &
                                                            nnodes_repl, offset_inode
      INTEGER(kind=int_8), ALLOCATABLE, DIMENSION(:)     :: node_merged_indices
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: inodes_orig, k_cells
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: nl_local
      INTEGER, DIMENSION(3)                              :: max_ijk_cell_image_local

      CALL timeset(routineN, handle)
      CPASSERT(.NOT. ALLOCATED(repl_nl))
      CPASSERT(.NOT. ALLOCATED(n_dbcsr_cell_images_to_merge))
      CPASSERT(.NOT. ALLOCATED(dbcsr_cell_image_to_merge))

      CALL get_nnodes_local(nnodes_local, max_ijk_cell_image_local, max_ijk_cell_image, sab_nl, para_env, particle_coords, cell)

      nnodes_per_proc(:) = 0

      IF (nnodes_local > 0) THEN
         ALLOCATE (nl_local(neighbor_list_dim1, nnodes_local))
         ALLOCATE (node_merged_indices(nnodes_local))
         ALLOCATE (k_cells(nnodes_local))
         CALL get_nl_nodes_local(nl_local, max_ijk_cell_image_local, max_ijk_cell_image, sab_nl, particle_coords, cell, &
                                 cell_to_index, do_merge, node_merged_indices, k_cells)

         ALLOCATE (inodes_orig(nnodes_local))
         CALL sort(node_merged_indices, nnodes_local, inodes_orig)

         nnodes_merged = 1
         DO inode = 2, nnodes_local
            IF (node_merged_indices(inode) > node_merged_indices(inode - 1)) nnodes_merged = nnodes_merged + 1
         END DO
      ELSE
         nnodes_merged = 0
      END IF

      nnodes_per_proc(para_env%mepos) = nnodes_merged
      CALL para_env%sum(nnodes_per_proc)

      nnodes_repl = SUM(nnodes_per_proc(:))
      ALLOCATE (repl_nl(neighbor_list_dim1, nnodes_repl))

      IF (nnodes_local > 0) THEN
         IF (para_env%mepos > 0) THEN
            offset_inode = SUM(nnodes_per_proc(0:para_env%mepos - 1))
         ELSE
            offset_inode = 0
         END IF

         ALLOCATE (n_dbcsr_cell_images_to_merge(nnodes_merged))
         ALLOCATE (dbcsr_cell_image_to_merge(nnodes_local))
         n_dbcsr_cell_images_to_merge(:) = 0

         nnodes_merged = 1 !offset_inode + 1
         repl_nl(:, offset_inode + 1) = nl_local(:, inodes_orig(1))
         n_dbcsr_cell_images_to_merge(1) = 1
         dbcsr_cell_image_to_merge(1) = nl_local(neighbor_list_dbcsr_image_index, inodes_orig(1))
         kcell_closest = k_cells(inodes_orig(1))
         DO inode = 2, nnodes_local
            IF (node_merged_indices(inode) > node_merged_indices(inode - 1)) THEN
               nnodes_merged = nnodes_merged + 1
               repl_nl(:, offset_inode + nnodes_merged) = nl_local(:, inodes_orig(inode))
               !n_dbcsr_cell_images_to_merge(nnodes_merged) = 1
               kcell_closest = k_cells(inodes_orig(inode))
            ELSE
               IF (ABS(k_cells(inodes_orig(inode))) < ABS(kcell_closest) .OR. &
                   (ABS(k_cells(inodes_orig(inode))) == ABS(kcell_closest) .AND. kcell_closest < 0)) THEN
                  repl_nl(:, offset_inode + nnodes_merged) = nl_local(:, inodes_orig(inode))
                  kcell_closest = k_cells(inodes_orig(inode))
               END IF
            END IF
            dbcsr_cell_image_to_merge(inode) = nl_local(neighbor_list_dbcsr_image_index, inodes_orig(inode))
            n_dbcsr_cell_images_to_merge(nnodes_merged) = n_dbcsr_cell_images_to_merge(nnodes_merged) + 1
         END DO

         IF (debug_this_module) THEN
            CPASSERT(SUM(n_dbcsr_cell_images_to_merge) == nnodes_local)
         END IF

         DEALLOCATE (inodes_orig)
         DEALLOCATE (node_merged_indices, k_cells)
         DEALLOCATE (nl_local)
      END IF

      IF (para_env%num_pe > 1) THEN
         offset_inode = 0
         DO iproc = 0, para_env%num_pe - 1
            IF (nnodes_per_proc(iproc) > 0) THEN
               CALL para_env%bcast(repl_nl(:, offset_inode + 1:offset_inode + nnodes_per_proc(iproc)), iproc)
               offset_inode = offset_inode + nnodes_per_proc(iproc)
            END IF
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE replicate_neighbour_list

! **************************************************************************************************
!> \brief Count number of DBCSR matrix elements that should be received from (*,nelements_dbcsr_recv)
!>        and send to (*,nelements_dbcsr_send) each parallel process.
!> \param nelements_per_proc number of non-zero matrix elements for each MPI process
!> \param nnodes_per_proc    number of non-zero DBCSR matrix blocks (neighbour-list nodes)
!> \param nl_repl            replicated neighbour list
!> \param matrix_dbcsr_kp    DBCSR matrix
!> \param symmetric          whether the DBCSR matrix is a symmetric one
!> \param para_env           parallel environment
!> \param gather_root        if >=0, gather all non-zero matrix element on the MPI process with
!>                           gather_root rank (useful for bulk transport calculation).
!>                           If <0, distribute non-zero matrix element across all MPI processes
! **************************************************************************************************
   SUBROUTINE count_remote_dbcsr_elements(nelements_per_proc, nnodes_per_proc, nl_repl, matrix_dbcsr_kp, &
                                          symmetric, para_env, gather_root)
      INTEGER(kind=int_8), DIMENSION(0:, :), INTENT(out) :: nelements_per_proc
      INTEGER, DIMENSION(0:), INTENT(in)                 :: nnodes_per_proc
      INTEGER, DIMENSION(:, :), INTENT(in)               :: nl_repl
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(in)       :: matrix_dbcsr_kp
      LOGICAL, INTENT(in)                                :: symmetric
      TYPE(mp_para_env_type), INTENT(in), POINTER        :: para_env
      INTEGER, INTENT(in)                                :: gather_root

      CHARACTER(len=*), PARAMETER :: routineN = 'count_remote_dbcsr_elements'

      INTEGER :: first_row_minus_one, handle, icol_blk, image, image_transp, inode, inode_proc, &
         iproc, iproc_orb, irow_blk, irow_local, mepos, ncols_blk, ncols_local, nnodes_proc, &
         nprocs, nrows_blk, nrows_local, offset_inode
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_offset, col_blk_size, &
                                                            row_blk_offset, row_blk_size
      LOGICAL                                            :: do_distribute

      CALL timeset(routineN, handle)
      nelements_per_proc(:, :) = 0
      mepos = para_env%mepos
      nprocs = para_env%num_pe
      do_distribute = gather_root < 0
      IF (debug_this_module) THEN
         CPASSERT(SIZE(nnodes_per_proc) == nprocs)
      END IF

      CALL dbcsr_get_info(matrix=matrix_dbcsr_kp(1)%matrix, &
                          nblkrows_total=nrows_blk, nblkcols_total=ncols_blk, &
                          row_blk_size=row_blk_size, col_blk_size=col_blk_size, &
                          row_blk_offset=row_blk_offset, col_blk_offset=col_blk_offset)

      offset_inode = 0
      iproc_orb = gather_root ! if do_distribute == .FALSE., send all matrix elements to MPI process with rank gather_root
      DO iproc = LBOUND(nnodes_per_proc, 1), UBOUND(nnodes_per_proc, 1)
         nnodes_proc = nnodes_per_proc(iproc)
         DO inode_proc = 1, nnodes_proc
            inode = inode_proc + offset_inode

            irow_blk = nl_repl(neighbor_list_iatom_index, inode)
            icol_blk = nl_repl(neighbor_list_jatom_index, inode)
            CPASSERT(irow_blk <= icol_blk .OR. (.NOT. symmetric))
            image = nl_repl(neighbor_list_siesta_image_index, inode)
            image_transp = nl_repl(neighbor_list_siesta_transp_image_index, inode)

            IF (image > 0) THEN
               nrows_local = row_blk_size(irow_blk)
               first_row_minus_one = row_blk_offset(irow_blk) - 1
               ncols_local = col_blk_size(icol_blk)
               DO irow_local = 1, nrows_local
                  IF (do_distribute) THEN
#if defined(__SMEAGOL)
                     CALL WhichNodeOrb(irow_local + first_row_minus_one, nprocs, iproc_orb)
#else
                     CALL cp_abort(__LOCATION__, &
                                   "CP2K was compiled with no SMEAGOL support.")
#endif
                  END IF
                  IF (iproc_orb == mepos) THEN
                     nelements_per_proc(iproc, nelements_dbcsr_recv) = nelements_per_proc(iproc, nelements_dbcsr_recv) + &
                                                                       ncols_local
                  END IF

                  IF (iproc == mepos) THEN
                     nelements_per_proc(iproc_orb, nelements_dbcsr_send) = nelements_per_proc(iproc_orb, nelements_dbcsr_send) + &
                                                                           ncols_local
                  END IF
               END DO
            END IF

            ! transposed block
            IF (image_transp > 0) THEN
               nrows_local = col_blk_size(icol_blk)
               first_row_minus_one = col_blk_offset(icol_blk) - 1
               ncols_local = row_blk_size(irow_blk)
               DO irow_local = 1, nrows_local
                  IF (do_distribute) THEN
#if defined(__SMEAGOL)
                     CALL WhichNodeOrb(irow_local + first_row_minus_one, nprocs, iproc_orb)
#else
                     CALL cp_abort(__LOCATION__, &
                                   "CP2K was compiled with no SMEAGOL support.")
#endif
                  END IF
                  IF (iproc_orb == mepos) THEN
                     nelements_per_proc(iproc, nelements_dbcsr_recv) = nelements_per_proc(iproc, nelements_dbcsr_recv) + &
                                                                       ncols_local
                  END IF

                  IF (iproc == mepos) THEN
                     nelements_per_proc(iproc_orb, nelements_dbcsr_send) = nelements_per_proc(iproc_orb, nelements_dbcsr_send) + &
                                                                           ncols_local
                  END IF
               END DO
            END IF
         END DO
         offset_inode = offset_inode + nnodes_proc
      END DO
      CALL timestop(handle)
   END SUBROUTINE count_remote_dbcsr_elements

! **************************************************************************************************
!> \brief Construct list of non-zero matrix elements' indices in SIESTA format.
!> \param n_nonzero_cols  number of non-zero matrix elements on each matrix row local to the current
!>                        MPI process
!> \param row_offset      offset of the first non-zero matrix elements for each locally-stores row
!> \param col_index       sorted list of column indices of non-zero matrix element
!> \param packed_index    original order of non-sorted column indices
!> \param nl_repl         replicated neighbour list
!> \param matrix_dbcsr_kp DBCSR matrix
!> \param symmetric       whether the DBCSR matrix is a symmetric one
!> \param para_env        parallel environment
!> \param gather_root     if >=0, gather all non-zero matrix element on the MPI process with
!>                        gather_root rank (useful for bulk transport calculation).
!>                        If <0, distribute non-zero matrix element across all MPI processes
! **************************************************************************************************
   SUBROUTINE get_nonzero_element_indices(n_nonzero_cols, row_offset, col_index, packed_index, &
                                          nl_repl, matrix_dbcsr_kp, symmetric, para_env, gather_root)
      INTEGER, DIMENSION(:), INTENT(out)                 :: n_nonzero_cols, row_offset, col_index, &
                                                            packed_index
      INTEGER, DIMENSION(:, :), INTENT(in)               :: nl_repl
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(in)       :: matrix_dbcsr_kp
      LOGICAL, INTENT(in)                                :: symmetric
      TYPE(mp_para_env_type), INTENT(in), POINTER        :: para_env
      INTEGER, INTENT(in)                                :: gather_root

      CHARACTER(len=*), PARAMETER :: routineN = 'get_nonzero_element_indices'

      INTEGER :: first_col_minus_one, first_row_minus_one, handle, icol_blk, icol_local, &
         icol_offset, image, image_transp, inode, irow_blk, irow_local, irow_proc, mepos, &
         ncols_blk, ncols_local, ncols_total, nnodes, nprocs, nrows_blk, nrows_local, nrows_total
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_offset, col_blk_size, &
                                                            row_blk_offset, row_blk_size
      LOGICAL                                            :: do_distribute, is_root_rank

      CALL timeset(routineN, handle)
      n_nonzero_cols(:) = 0
      mepos = para_env%mepos
      nprocs = para_env%num_pe
      do_distribute = gather_root < 0
      is_root_rank = gather_root == mepos

      CALL dbcsr_get_info(matrix=matrix_dbcsr_kp(1)%matrix, &
                          nblkrows_total=nrows_blk, nblkcols_total=ncols_blk, &
                          nfullrows_total=nrows_total, nfullcols_total=ncols_total, &
                          row_blk_size=row_blk_size, col_blk_size=col_blk_size, &
                          row_blk_offset=row_blk_offset, col_blk_offset=col_blk_offset)

      nnodes = SIZE(nl_repl, 2)
      DO inode = 1, nnodes
         irow_blk = nl_repl(neighbor_list_iatom_index, inode)
         icol_blk = nl_repl(neighbor_list_jatom_index, inode)
         CPASSERT(irow_blk <= icol_blk .OR. (.NOT. symmetric))
         image = nl_repl(neighbor_list_siesta_image_index, inode)
         image_transp = nl_repl(neighbor_list_siesta_transp_image_index, inode)

         IF (image > 0) THEN
            nrows_local = row_blk_size(irow_blk)
            first_row_minus_one = row_blk_offset(irow_blk) - 1
            ncols_local = col_blk_size(icol_blk)
            DO irow_local = 1, nrows_local
               IF (do_distribute) THEN
#if defined(__SMEAGOL)
                  CALL GlobalToLocalOrb(irow_local + first_row_minus_one, mepos, nprocs, irow_proc)
#else
                  CALL cp_abort(__LOCATION__, &
                                "CP2K was compiled with no SMEAGOL support.")
#endif
               ELSE
                  IF (is_root_rank) THEN
                     irow_proc = irow_local + first_row_minus_one
                  ELSE
                     irow_proc = 0
                  END IF
               END IF
               IF (irow_proc > 0) THEN
                  n_nonzero_cols(irow_proc) = n_nonzero_cols(irow_proc) + ncols_local
               END IF
            END DO
         END IF

         ! transposed block
         IF (image_transp > 0) THEN
            nrows_local = col_blk_size(icol_blk)
            first_row_minus_one = col_blk_offset(icol_blk) - 1
            ncols_local = row_blk_size(irow_blk)
            DO irow_local = 1, nrows_local
               IF (do_distribute) THEN
#if defined(__SMEAGOL)
                  CALL GlobalToLocalOrb(irow_local + first_row_minus_one, mepos, nprocs, irow_proc)
#else
                  CALL cp_abort(__LOCATION__, &
                                "CP2K was compiled with no SMEAGOL support.")
#endif
               ELSE
                  IF (is_root_rank) THEN
                     irow_proc = irow_local + first_row_minus_one
                  ELSE
                     irow_proc = 0
                  END IF
               END IF
               IF (irow_proc > 0) THEN
                  n_nonzero_cols(irow_proc) = n_nonzero_cols(irow_proc) + ncols_local
               END IF
            END DO
         END IF
      END DO

      row_offset(1) = 0
      DO irow_local = 1, SIZE(n_nonzero_cols) - 1
         row_offset(irow_local + 1) = row_offset(irow_local) + n_nonzero_cols(irow_local)
      END DO

      n_nonzero_cols(:) = 0
      col_index(:) = 0
      DO inode = 1, nnodes
         irow_blk = nl_repl(neighbor_list_iatom_index, inode)
         icol_blk = nl_repl(neighbor_list_jatom_index, inode)
         CPASSERT(irow_blk <= icol_blk .OR. (.NOT. symmetric))
         image = nl_repl(neighbor_list_siesta_image_index, inode)
         image_transp = nl_repl(neighbor_list_siesta_transp_image_index, inode)

         IF (image > 0) THEN
            nrows_local = row_blk_size(irow_blk)
            first_row_minus_one = row_blk_offset(irow_blk) - 1
            ncols_local = col_blk_size(icol_blk)
            first_col_minus_one = col_blk_offset(icol_blk) + (image - 1)*ncols_total - 1
            DO irow_local = 1, nrows_local
               IF (do_distribute) THEN
#if defined(__SMEAGOL)
                  CALL GlobalToLocalOrb(irow_local + first_row_minus_one, mepos, nprocs, irow_proc)
#else
                  CALL cp_abort(__LOCATION__, &
                                "CP2K was compiled with no SMEAGOL support.")
#endif
               ELSE
                  IF (is_root_rank) THEN
                     irow_proc = irow_local + first_row_minus_one
                  ELSE
                     irow_proc = 0
                  END IF
               END IF
               IF (irow_proc > 0) THEN
                  icol_offset = row_offset(irow_proc) + n_nonzero_cols(irow_proc)
                  DO icol_local = 1, ncols_local
                     col_index(icol_offset + icol_local) = first_col_minus_one + icol_local
                  END DO
                  n_nonzero_cols(irow_proc) = n_nonzero_cols(irow_proc) + ncols_local
               END IF
            END DO
         END IF

         ! transposed block
         IF (image_transp > 0) THEN
            nrows_local = col_blk_size(icol_blk)
            first_row_minus_one = col_blk_offset(icol_blk) - 1
            ncols_local = row_blk_size(irow_blk)
            first_col_minus_one = row_blk_offset(irow_blk) + (image_transp - 1)*nrows_total - 1
            DO irow_local = 1, nrows_local
               IF (do_distribute) THEN
#if defined(__SMEAGOL)
                  CALL GlobalToLocalOrb(irow_local + first_row_minus_one, mepos, nprocs, irow_proc)
#else
                  CALL cp_abort(__LOCATION__, &
                                "CP2K was compiled with no SMEAGOL support.")
#endif
               ELSE
                  IF (is_root_rank) THEN
                     irow_proc = irow_local + first_row_minus_one
                  ELSE
                     irow_proc = 0
                  END IF
               END IF
               IF (irow_proc > 0) THEN
                  icol_offset = row_offset(irow_proc) + n_nonzero_cols(irow_proc)
                  DO icol_local = 1, ncols_local
                     col_index(icol_offset + icol_local) = first_col_minus_one + icol_local
                  END DO
                  n_nonzero_cols(irow_proc) = n_nonzero_cols(irow_proc) + ncols_local
               END IF
            END DO
         END IF
      END DO

      IF (SIZE(n_nonzero_cols) > 0) THEN
         DO irow_local = 1, SIZE(n_nonzero_cols)
            CALL sort(col_index(row_offset(irow_local) + 1:row_offset(irow_local) + n_nonzero_cols(irow_local)), &
                      n_nonzero_cols(irow_local), &
                      packed_index(row_offset(irow_local) + 1:row_offset(irow_local) + n_nonzero_cols(irow_local)))
         END DO
      END IF

      CALL timestop(handle)
   END SUBROUTINE get_nonzero_element_indices

! **************************************************************************************************
!> \brief Get absolute i, j, and k indices of cell image for the given DBCSR matrix block.
!>        Effective cell image recorded in the neighbour-list (where the matrix block is actually
!>        stored) depends on atomic coordinates can be significantly different.
!> \param cell_ijk      array with 3 indices along the cell's vectors
!> \param r_ij          actual interatomic distance (vector R_j - r_i), where R_j is the coordinates
!>                      of the j-th atom in the supercell
!> \param r_i           coordinates of the i-th atom in the primary unit cell
!> \param r_j           coordinates of the j-th atom in the primary unit cell
!> \param cell          unit cell
! **************************************************************************************************
   SUBROUTINE get_negf_cell_ijk(cell_ijk, r_ij, r_i, r_j, cell)
      INTEGER, DIMENSION(3), INTENT(out)                 :: cell_ijk
      REAL(kind=dp), DIMENSION(3), INTENT(in)            :: r_ij, r_i, r_j
      TYPE(cell_type), INTENT(in), POINTER               :: cell

      REAL(kind=dp), DIMENSION(3)                        :: coords_scaled, r

      r(:) = r_ij(:) + r_i(:) - r_j(:)
      CALL real_to_scaled(coords_scaled, r, cell)
      cell_ijk(:) = NINT(coords_scaled(:))
   END SUBROUTINE get_negf_cell_ijk

! **************************************************************************************************
!> \brief Return the index of an integer number in the sequence 0, 1, -1, ..., n, -n, ...
!>        (canonical enumeration of integers, oeis.org/A001057).
!> \param  inum     integer number [in]
!> \return index of 'inum' in A001057
!> \note Cell images in SMEAGOL / SIESTA are ordered according to A001057. Therefore this
!>       function converts the absolute index of a cell replica along some (x/y/z) dimension
!>       into its corresponding SIESTA's index.
! **************************************************************************************************
   ELEMENTAL FUNCTION index_in_canonical_enumeration(inum) RESULT(ind)
      INTEGER, INTENT(in)                                :: inum
      INTEGER                                            :: ind

      INTEGER                                            :: inum_abs, is_non_positive

      inum_abs = ABS(inum)
      !IF (inum <= 0) THEN; is_non_positive = 1; ELSE; is_non_positive = 0; END IF
      is_non_positive = MERGE(1, 0, inum <= 0)

      ! inum =  0 -> inum_abs = 0, is_non_positive = 1 -> ind = 1
      ! inum =  1 -> inum_abs = 1, is_non_positive = 0 -> ind = 2
      ! inum = -1 -> inum_abs = 1, is_non_positive = 1 -> ind = 3
      ind = 2*inum_abs + is_non_positive
   END FUNCTION index_in_canonical_enumeration

! **************************************************************************************************
!> \brief Return an integer number according to its index in the sequence 0, 1, -1, ..., n, -n, ...
!>        (canonical enumeration of integers, oeis.org/A001057)
!> \param  ind       index in A001057 starting from 1
!> \return integer number according to its position 'ind' in A001057
!> \note Cell images in SMEAGOL / SIESTA are ordered according to A001057. Therefore this
!>       function converts SIESTA's index of a cell replica along some (x/y/z) dimension
!>       into the corresponding absolute index.
! **************************************************************************************************
   ELEMENTAL FUNCTION number_from_canonical_enumeration(ind) RESULT(inum)
      INTEGER, INTENT(in)                                :: ind
      INTEGER                                            :: inum

      ! ind < 1 is invalid
      ! ind = 1 -> SIGN(0, -1) =  0
      ! ind = 2 -> SIGN(1,  0) =  1
      ! ind = 3 -> SIGN(1, -1) = -1
      inum = SIGN(ind/2, -MOD(ind, 2))
   END FUNCTION number_from_canonical_enumeration

! **************************************************************************************************
!> \brief   Apply periodic boundary conditions defined by a simulation cell to a position
!>          vector r. Similar to pbc1 from cell_types.F but returns unscaled coordinates from
!>          the scaled range [0, 1) instead of [-0.5, 0.5)
!> \param r_pbc   position vector subject to the periodic boundary conditions [out]
!> \param r       initial position vector [in]
!> \param cell    simulation unit cell [in]
! **************************************************************************************************
   PURE SUBROUTINE pbc_0_1(r_pbc, r, cell)
      REAL(KIND=dp), DIMENSION(3), INTENT(out)           :: r_pbc
      REAL(KIND=dp), DIMENSION(3), INTENT(in)            :: r
      TYPE(cell_type), INTENT(in), POINTER               :: cell

      REAL(KIND=dp), DIMENSION(3)                        :: s

      IF (cell%orthorhombic) THEN
         r_pbc(1) = r(1) - cell%hmat(1, 1)*cell%perd(1)*REAL(FLOOR(cell%h_inv(1, 1)*r(1)), dp)
         r_pbc(2) = r(2) - cell%hmat(2, 2)*cell%perd(2)*REAL(FLOOR(cell%h_inv(2, 2)*r(2)), dp)
         r_pbc(3) = r(3) - cell%hmat(3, 3)*cell%perd(3)*REAL(FLOOR(cell%h_inv(3, 3)*r(3)), dp)
      ELSE
         s(1) = cell%h_inv(1, 1)*r(1) + cell%h_inv(1, 2)*r(2) + cell%h_inv(1, 3)*r(3)
         s(2) = cell%h_inv(2, 1)*r(1) + cell%h_inv(2, 2)*r(2) + cell%h_inv(2, 3)*r(3)
         s(3) = cell%h_inv(3, 1)*r(1) + cell%h_inv(3, 2)*r(2) + cell%h_inv(3, 3)*r(3)
         s(1) = s(1) - cell%perd(1)*REAL(FLOOR(s(1)), dp)
         s(2) = s(2) - cell%perd(2)*REAL(FLOOR(s(2)), dp)
         s(3) = s(3) - cell%perd(3)*REAL(FLOOR(s(3)), dp)
         r_pbc(1) = cell%hmat(1, 1)*s(1) + cell%hmat(1, 2)*s(2) + cell%hmat(1, 3)*s(3)
         r_pbc(2) = cell%hmat(2, 1)*s(1) + cell%hmat(2, 2)*s(2) + cell%hmat(2, 3)*s(3)
         r_pbc(3) = cell%hmat(3, 1)*s(1) + cell%hmat(3, 2)*s(2) + cell%hmat(3, 3)*s(3)
      END IF
   END SUBROUTINE pbc_0_1

! **************************************************************************************************
!> \brief Computes the number of send requests from this MPI process to all the other processes.
!>        Alternatively computes the number of recv requests per MPI process that the given process
!>        expects.
!> \param mepos                    MPI rank of the given process
!> \param nelements_per_proc       number of element to send / receive
!> \param max_nelements_per_packet maximum number of elements per single MPI request
!> \return number of MPI requests
! **************************************************************************************************
   PURE FUNCTION get_number_of_mpi_sendrecv_requests(mepos, nelements_per_proc, max_nelements_per_packet) RESULT(nrequests)
      INTEGER, INTENT(in)                                :: mepos
      INTEGER(kind=int_8), DIMENSION(0:), INTENT(in)     :: nelements_per_proc
      INTEGER(kind=int_8), INTENT(in)                    :: max_nelements_per_packet
      INTEGER                                            :: nrequests

      INTEGER                                            :: iproc

      nrequests = 0
      DO iproc = LBOUND(nelements_per_proc, 1), UBOUND(nelements_per_proc, 1)
         ! there is no need to send data to the same MPI process
         IF (iproc /= mepos) THEN
            nrequests = nrequests + INT(nelements_per_proc(iproc)/max_nelements_per_packet)
            IF (MOD(nelements_per_proc(iproc), max_nelements_per_packet) > 0) &
               nrequests = nrequests + 1
         END IF
      END DO
   END FUNCTION get_number_of_mpi_sendrecv_requests

! **************************************************************************************************
!> \brief Map non-zero matrix elements on to MPI requests.
!> \param element_offset           offset (index-1) of the first element [out]
!> \param nelements_per_request    number of element for each request [out]
!> \param peer_rank                rank of a peering MPI process
!> \param tag                      MPI tag
!> \param mepos                    MPI rank of a given MPI process
!> \param nelements_per_proc       number of element to send / receive by the current MPI process
!> \param max_nelements_per_packet maximum number of elements per single MPI request
! **************************************************************************************************
   SUBROUTINE assign_nonzero_elements_to_requests(element_offset, nelements_per_request, peer_rank, tag, &
                                                  mepos, nelements_per_proc, max_nelements_per_packet)
      INTEGER(kind=int_8), DIMENSION(:), INTENT(out)     :: element_offset, nelements_per_request
      INTEGER, DIMENSION(:), INTENT(out)                 :: peer_rank, tag
      INTEGER, INTENT(in)                                :: mepos
      INTEGER(kind=int_8), DIMENSION(0:), INTENT(in)     :: nelements_per_proc
      INTEGER(kind=int_8), INTENT(in)                    :: max_nelements_per_packet

      INTEGER                                            :: iproc, irequest, nrequests, &
                                                            request_offset
      INTEGER(kind=int_8)                                :: element_offset_tmp, nelements

      request_offset = 0
      element_offset_tmp = 0
      DO iproc = LBOUND(nelements_per_proc, 1), UBOUND(nelements_per_proc, 1)
         IF (iproc /= mepos) THEN
            nrequests = INT(nelements_per_proc(iproc)/max_nelements_per_packet)
            IF (MOD(nelements_per_proc(iproc), max_nelements_per_packet) > 0) nrequests = nrequests + 1
            CPASSERT(nrequests <= max_mpi_rank + 1)
            IF (nrequests > 0) THEN
               nelements = nelements_per_proc(iproc)/nrequests
               IF (nelements_per_proc(iproc) - nelements*nrequests > 0) nelements = nelements + 1
               CPASSERT(nelements <= max_nelements_per_packet)

               DO irequest = 1, nrequests
                  element_offset(request_offset + irequest) = (irequest - 1)*nelements + element_offset_tmp
                  IF (irequest < nrequests) THEN
                     nelements_per_request(request_offset + irequest) = nelements
                  ELSE
                     nelements_per_request(request_offset + irequest) = nelements_per_proc(iproc) - nelements*(nrequests - 1)
                  END IF
                  peer_rank(request_offset + irequest) = iproc
                  tag(request_offset + irequest) = irequest - 1
               END DO
            END IF
            request_offset = request_offset + nrequests
         END IF
         element_offset_tmp = element_offset_tmp + nelements_per_proc(iproc)
      END DO

      IF (debug_this_module) THEN
         CPASSERT(SIZE(element_offset) == request_offset)
         CPASSERT(SIZE(nelements_per_request) == request_offset)
         CPASSERT(SIZE(peer_rank) == request_offset)
         CPASSERT(SIZE(tag) == request_offset)
      END IF
   END SUBROUTINE assign_nonzero_elements_to_requests

END MODULE smeagol_matrix_utils
