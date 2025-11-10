!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Helper routines to manipulate with matrices.
! **************************************************************************************************
MODULE negf_matrix_utils
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_finalize, &
        dbcsr_get_block_p, dbcsr_init_p, dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_set, &
        dbcsr_type, dbcsr_type_no_symmetry
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_type,&
                                              mp_request_type
   USE negf_alloc_types,                ONLY: negf_allocatable_rvector
   USE negf_atom_map,                   ONLY: negf_atom_map_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_matrix_utils'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

   PUBLIC :: number_of_atomic_orbitals, negf_copy_fm_submat_to_dbcsr, negf_copy_sym_dbcsr_to_fm_submat
   PUBLIC :: negf_copy_contact_matrix, negf_reference_contact_matrix
   PUBLIC :: invert_cell_to_index, get_index_by_cell

CONTAINS

! **************************************************************************************************
!> \brief Compute the number of atomic orbitals of the given set of atoms.
!> \param subsys    QuickStep subsystem
!> \param atom_list list of selected atom; when absent all the atoms are taken into account
!> \return number of atomic orbitals
!> \par History
!>   * 02.2017 created [Sergey Chulkov]
! **************************************************************************************************
   FUNCTION number_of_atomic_orbitals(subsys, atom_list) RESULT(nao)
      TYPE(qs_subsys_type), POINTER                      :: subsys
      INTEGER, DIMENSION(:), INTENT(in), OPTIONAL        :: atom_list
      INTEGER                                            :: nao

      INTEGER                                            :: iatom, natoms
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nsgfs
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL qs_subsys_get(subsys, particle_set=particle_set, qs_kind_set=qs_kind_set)
      ALLOCATE (nsgfs(SIZE(particle_set)))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=nsgfs)

      IF (PRESENT(atom_list)) THEN
         natoms = SIZE(atom_list)
         nao = 0

         DO iatom = 1, natoms
            nao = nao + nsgfs(atom_list(iatom))
         END DO
      ELSE
         nao = SUM(nsgfs)
      END IF

      DEALLOCATE (nsgfs)
   END FUNCTION number_of_atomic_orbitals

! **************************************************************************************************
!> \brief Populate relevant blocks of the DBCSR matrix using data from a ScaLAPACK matrix.
!>        Irrelevant blocks of the DBCSR matrix are kept untouched.
!> \param fm              dense matrix to copy
!> \param matrix          DBCSR matrix (modified on exit)
!> \param atomlist_row    set of atomic indices along the 1st (row) dimension
!> \param atomlist_col    set of atomic indices along the 2nd (column) dimension
!> \param subsys          subsystem environment
!> \par History
!>   * 02.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE negf_copy_fm_submat_to_dbcsr(fm, matrix, atomlist_row, atomlist_col, subsys)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm
      TYPE(dbcsr_type), POINTER                          :: matrix
      INTEGER, DIMENSION(:), INTENT(in)                  :: atomlist_row, atomlist_col
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_copy_fm_submat_to_dbcsr'

      INTEGER :: first_sgf_col, first_sgf_row, handle, iatom_col, iatom_row, icol, irow, &
         natoms_col, natoms_row, ncols, nparticles, nrows
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nsgfs
      LOGICAL                                            :: found
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: fm_block
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: sm_block
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(matrix))
      CPASSERT(ASSOCIATED(subsys))

      CALL cp_fm_get_info(fm, nrow_global=nrows, ncol_global=ncols)

      CALL qs_subsys_get(subsys, particle_set=particle_set, qs_kind_set=qs_kind_set)

      natoms_row = SIZE(atomlist_row)
      natoms_col = SIZE(atomlist_col)
      nparticles = SIZE(particle_set)

      ALLOCATE (nsgfs(nparticles))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=nsgfs)

      ALLOCATE (fm_block(nrows, ncols))
      CALL cp_fm_get_submatrix(fm, fm_block)

      first_sgf_col = 1
      DO iatom_col = 1, natoms_col
         first_sgf_row = 1
         DO iatom_row = 1, natoms_row
            CALL dbcsr_get_block_p(matrix=matrix, row=atomlist_row(iatom_row), col=atomlist_col(iatom_col), &
                                   block=sm_block, found=found)
            IF (found) THEN
               ! the following LAPACK call violates the coding convention
               !CALL dlacpy('F', nsgfs(atomlist_row(iatom_row)), nsgfs(atomlist_col(iatom_col)), &
               !            fm_block(first_sgf_row, first_sgf_col), SIZE(fm_block, 1), sm_block(1, 1), SIZE(sm_block, 1))
               nrows = nsgfs(atomlist_row(iatom_row))
               ncols = nsgfs(atomlist_col(iatom_col))
               DO icol = 1, ncols
                  DO irow = 1, nrows
                     sm_block(irow, icol) = fm_block(first_sgf_row + irow - 1, first_sgf_col + icol - 1)
                  END DO
               END DO
            END IF

            first_sgf_row = first_sgf_row + nsgfs(atomlist_row(iatom_row))
         END DO
         first_sgf_col = first_sgf_col + nsgfs(atomlist_col(iatom_col))
      END DO

      DEALLOCATE (fm_block)
      DEALLOCATE (nsgfs)

      CALL timestop(handle)
   END SUBROUTINE negf_copy_fm_submat_to_dbcsr

! **************************************************************************************************
!> \brief Extract part of the DBCSR matrix based on selected atoms and copy it into a dense matrix.
!> \param matrix          DBCSR matrix
!> \param fm              dense matrix (created and initialised on exit)
!> \param atomlist_row    set of atomic indices along the 1st (row) dimension
!> \param atomlist_col    set of atomic indices along the 2nd (column) dimension
!> \param subsys          subsystem environment
!> \param mpi_comm_global MPI communicator which was used to distribute blocks of the DBCSR matrix.
!>                        If missed, assume that both DBCSR and ScaLapack matrices are distributed
!>                        across the same set of processors
!> \param do_upper_diag   initialise upper-triangular part of the dense matrix as well as diagonal elements
!> \param do_lower        initialise lower-triangular part of the dense matrix
!> \par History
!>   * 02.2017 created [Sergey Chulkov]
!> \note A naive implementation that copies relevant local DBCSR blocks into a 2-D matrix,
!>       performs collective summation, and then distributes the result. This approach seems to be
!>       optimal when processors are arranged into several independent MPI subgroups due to the fact
!>       that every subgroup automatically holds the copy of the dense matrix at the end, so
!>       we can avoid the final replication stage.
! **************************************************************************************************
   SUBROUTINE negf_copy_sym_dbcsr_to_fm_submat(matrix, fm, atomlist_row, atomlist_col, subsys, &
                                               mpi_comm_global, do_upper_diag, do_lower)
      TYPE(dbcsr_type), POINTER                          :: matrix
      TYPE(cp_fm_type), INTENT(IN)                       :: fm
      INTEGER, DIMENSION(:), INTENT(in)                  :: atomlist_row, atomlist_col
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CLASS(mp_comm_type), INTENT(in)                     :: mpi_comm_global
      LOGICAL, INTENT(in)                                :: do_upper_diag, do_lower

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_copy_sym_dbcsr_to_fm_submat'

      INTEGER :: handle, iatom_col, iatom_row, icol, irow, natoms_col, natoms_row, ncols_fm, &
                 nparticles, nrows_fm, offset_sgf_col, offset_sgf_row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nsgfs
      LOGICAL                                            :: found
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: r2d
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: sm_block
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(matrix))
      CPASSERT(ASSOCIATED(subsys))

      CALL qs_subsys_get(subsys, particle_set=particle_set, qs_kind_set=qs_kind_set)

      natoms_row = SIZE(atomlist_row)
      natoms_col = SIZE(atomlist_col)
      nparticles = SIZE(particle_set)

      ALLOCATE (nsgfs(nparticles))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=nsgfs)

      CALL cp_fm_get_info(fm, nrow_global=nrows_fm, ncol_global=ncols_fm, para_env=para_env)

      IF (debug_this_module) THEN
         CPASSERT(SUM(nsgfs(atomlist_row(:))) == nrows_fm)
         CPASSERT(SUM(nsgfs(atomlist_col(:))) == ncols_fm)
      END IF

      ALLOCATE (r2d(nrows_fm, ncols_fm))
      r2d(:, :) = 0.0_dp

      offset_sgf_col = 0
      DO iatom_col = 1, natoms_col
         offset_sgf_row = 0

         DO iatom_row = 1, natoms_row
            IF (atomlist_row(iatom_row) <= atomlist_col(iatom_col)) THEN
               IF (do_upper_diag) THEN
                  CALL dbcsr_get_block_p(matrix=matrix, row=atomlist_row(iatom_row), col=atomlist_col(iatom_col), &
                                         block=sm_block, found=found)
               END IF
            ELSE
               IF (do_lower) THEN
                  CALL dbcsr_get_block_p(matrix=matrix, row=atomlist_col(iatom_col), col=atomlist_row(iatom_row), &
                                         block=sm_block, found=found)
               END IF
            END IF

            IF (found) THEN
               IF (atomlist_row(iatom_row) <= atomlist_col(iatom_col)) THEN
                  IF (do_upper_diag) THEN
                     DO icol = nsgfs(atomlist_col(iatom_col)), 1, -1
                        DO irow = nsgfs(atomlist_row(iatom_row)), 1, -1
                           r2d(offset_sgf_row + irow, offset_sgf_col + icol) = sm_block(irow, icol)
                        END DO
                     END DO
                  END IF
               ELSE
                  IF (do_lower) THEN
                     DO icol = nsgfs(atomlist_col(iatom_col)), 1, -1
                        DO irow = nsgfs(atomlist_row(iatom_row)), 1, -1
                           r2d(offset_sgf_row + irow, offset_sgf_col + icol) = sm_block(icol, irow)
                        END DO
                     END DO
                  END IF
               END IF
            END IF

            offset_sgf_row = offset_sgf_row + nsgfs(atomlist_row(iatom_row))
         END DO
         offset_sgf_col = offset_sgf_col + nsgfs(atomlist_col(iatom_col))
      END DO

      CALL mpi_comm_global%sum(r2d)

      CALL cp_fm_set_submatrix(fm, r2d)

      DEALLOCATE (r2d)
      DEALLOCATE (nsgfs)

      CALL timestop(handle)
   END SUBROUTINE negf_copy_sym_dbcsr_to_fm_submat

! **************************************************************************************************
!> \brief Driver routine to extract diagonal and off-diagonal blocks from a symmetric DBCSR matrix.
!> \param fm_cell0        extracted diagonal matrix block
!> \param fm_cell1        extracted off-diagonal matrix block
!> \param direction_axis  axis towards the secondary unit cell
!> \param matrix_kp       set of DBCSR matrices
!> \param atom_list0      list of atoms which belong to the primary contact unit cell
!> \param atom_list1      list of atoms which belong to the secondary contact unit cell
!> \param subsys          QuickStep subsystem
!> \param mpi_comm_global global MPI communicator
!> \param kpoints ...
!> \par History
!>   * 10.2017 created [Sergey Chulkov]
!>   * 10.2025 The subroutine is essentially modified. [Dmitry Ryndyk]
! **************************************************************************************************
   SUBROUTINE negf_copy_contact_matrix(fm_cell0, fm_cell1, direction_axis, matrix_kp, &
                                       atom_list0, atom_list1, subsys, mpi_comm_global, kpoints)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_cell0, fm_cell1
      INTEGER, INTENT(in)                                :: direction_axis
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(in), &
         POINTER                                         :: matrix_kp
      INTEGER, DIMENSION(:), INTENT(in)                  :: atom_list0, atom_list1
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CLASS(mp_comm_type), INTENT(in)                    :: mpi_comm_global
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_copy_contact_matrix'

      INTEGER                                            :: direction_axis_abs, handle, rep, ncell, ic
      TYPE(dbcsr_p_type), ALLOCATABLE, DIMENSION(:)      :: matrix_cells_raw
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_nosym
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: i_to_c
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: c_to_i

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(subsys))

      direction_axis_abs = ABS(direction_axis)

      CALL desymmetrize_matrix(matrix_kp, mat_nosym, c_to_i, i_to_c, kpoints)
      ncell = SIZE(i_to_c, 2) ! update the number of cells

      ! 0 -- primary unit cell;
      ! +- 1 -- upper- and lower-diagonal matrices for neighbor-cell matrix elements;
      ! +- 2 -- for control
      ALLOCATE (matrix_cells_raw(-2:2))
      DO rep = -2, 2
         NULLIFY (matrix_cells_raw(rep)%matrix)
         CALL dbcsr_init_p(matrix_cells_raw(rep)%matrix)
         CALL dbcsr_copy(matrix_cells_raw(rep)%matrix, mat_nosym(1)%matrix)
         CALL dbcsr_set(matrix_cells_raw(rep)%matrix, 0.0_dp)
      END DO

      DO ic = 1, ncell
         rep = i_to_c(direction_axis_abs, ic)
         IF (ABS(rep) <= 2) &
            CALL dbcsr_add(matrix_cells_raw(rep)%matrix, mat_nosym(ic)%matrix, 1.0_dp, 1.0_dp)
      END DO

      IF (direction_axis >= 0) THEN

         CALL negf_copy_sym_dbcsr_to_fm_submat(matrix_cells_raw(1)%matrix, fm_cell1, atom_list0, atom_list1, &
                                               subsys, mpi_comm_global, do_upper_diag=.TRUE., do_lower=.FALSE.)
         CALL negf_copy_sym_dbcsr_to_fm_submat(matrix_cells_raw(-1)%matrix, fm_cell0, atom_list0, atom_list1, &
                                               subsys, mpi_comm_global, do_upper_diag=.FALSE., do_lower=.TRUE.)

      ELSE

         CALL negf_copy_sym_dbcsr_to_fm_submat(matrix_cells_raw(1)%matrix, fm_cell1, atom_list0, atom_list1, &
                                               subsys, mpi_comm_global, do_upper_diag=.FALSE., do_lower=.TRUE.)
         CALL negf_copy_sym_dbcsr_to_fm_submat(matrix_cells_raw(-1)%matrix, fm_cell0, atom_list0, atom_list1, &
                                               subsys, mpi_comm_global, do_upper_diag=.TRUE., do_lower=.FALSE.)

      END IF
      CALL cp_fm_scale_and_add(1.0_dp, fm_cell1, 1.0_dp, fm_cell0)

      ! symmetric matrix fm_cell0
      CALL negf_copy_sym_dbcsr_to_fm_submat(matrix_cells_raw(0)%matrix, fm_cell0, atom_list0, atom_list0, &
                                            subsys, mpi_comm_global, do_upper_diag=.TRUE., do_lower=.TRUE.)

      ! clean up
      DEALLOCATE (c_to_i, i_to_c)
      DO ic = 1, ncell
         CALL dbcsr_release(mat_nosym(ic)%matrix)
      END DO
      CALL dbcsr_deallocate_matrix_set(mat_nosym)
      DO rep = -2, 2
         CALL dbcsr_deallocate_matrix(matrix_cells_raw(rep)%matrix)
      END DO
      DEALLOCATE (matrix_cells_raw)

      CALL timestop(handle)
   END SUBROUTINE negf_copy_contact_matrix

! **************************************************************************************************
!> \brief Extract part of the DBCSR matrix based on selected atoms and copy it into another DBCSR
!>        matrix.
!> \param matrix_contact  extracted DBCSR matrix
!> \param matrix_device   original DBCSR matrix
!> \param atom_list       list of selected atoms
!> \param atom_map        atomic map between device and contact force environments
!> \param para_env        parallel environment
! **************************************************************************************************
   SUBROUTINE negf_reference_contact_matrix(matrix_contact, matrix_device, atom_list, atom_map, para_env)
      TYPE(dbcsr_type), POINTER                          :: matrix_contact, matrix_device
      INTEGER, DIMENSION(:), INTENT(in)                  :: atom_list
      TYPE(negf_atom_map_type), DIMENSION(:), INTENT(in) :: atom_map
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_reference_contact_matrix'

      INTEGER                                            :: handle, i1, i2, iatom_col, iatom_row, &
                                                            icol, iproc, irow, max_atom, &
                                                            mepos_plus1, n1, n2, natoms, offset
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: recv_nelems, send_nelems
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: rank_contact, rank_device
      LOGICAL                                            :: found, transp
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: rblock
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: recv_handlers, send_handlers
      TYPE(negf_allocatable_rvector), ALLOCATABLE, &
         DIMENSION(:)                                    :: recv_packed_blocks, send_packed_blocks

      CALL timeset(routineN, handle)
      mepos_plus1 = para_env%mepos + 1

      natoms = SIZE(atom_list)
      max_atom = 0
      DO iatom_row = 1, natoms
         IF (atom_map(iatom_row)%iatom > max_atom) max_atom = atom_map(iatom_row)%iatom
      END DO

      ! find out which block goes to which node
      ALLOCATE (rank_contact(max_atom, max_atom))
      ALLOCATE (rank_device(max_atom, max_atom))

      rank_contact(:, :) = 0
      rank_device(:, :) = 0

      DO iatom_col = 1, natoms
         DO iatom_row = 1, iatom_col
            IF (atom_map(iatom_row)%iatom <= atom_map(iatom_col)%iatom) THEN
               icol = atom_map(iatom_col)%iatom
               irow = atom_map(iatom_row)%iatom
            ELSE
               icol = atom_map(iatom_row)%iatom
               irow = atom_map(iatom_col)%iatom
            END IF

            CALL dbcsr_get_block_p(matrix=matrix_device, &
                                   row=atom_list(iatom_row), col=atom_list(iatom_col), &
                                   block=rblock, found=found)
            IF (found) rank_device(irow, icol) = mepos_plus1

            CALL dbcsr_get_block_p(matrix=matrix_contact, row=irow, col=icol, block=rblock, found=found)
            IF (found) rank_contact(irow, icol) = mepos_plus1
         END DO
      END DO

      CALL para_env%sum(rank_device)
      CALL para_env%sum(rank_contact)

      ! compute number of packed matrix elements to send to / receive from each processor
      ALLOCATE (recv_nelems(para_env%num_pe))
      ALLOCATE (send_nelems(para_env%num_pe))
      recv_nelems(:) = 0
      send_nelems(:) = 0

      DO iatom_col = 1, natoms
         DO iatom_row = 1, iatom_col
            IF (atom_map(iatom_row)%iatom <= atom_map(iatom_col)%iatom) THEN
               icol = atom_map(iatom_col)%iatom
               irow = atom_map(iatom_row)%iatom
            ELSE
               icol = atom_map(iatom_row)%iatom
               irow = atom_map(iatom_col)%iatom
            END IF

            CALL dbcsr_get_block_p(matrix=matrix_device, &
                                   row=atom_list(iatom_row), col=atom_list(iatom_col), &
                                   block=rblock, found=found)
            IF (found) THEN
               iproc = rank_contact(irow, icol)
               IF (iproc > 0) &
                  send_nelems(iproc) = send_nelems(iproc) + SIZE(rblock)
            END IF

            CALL dbcsr_get_block_p(matrix=matrix_contact, row=irow, col=icol, block=rblock, found=found)
            IF (found) THEN
               iproc = rank_device(irow, icol)
               IF (iproc > 0) &
                  recv_nelems(iproc) = recv_nelems(iproc) + SIZE(rblock)
            END IF
         END DO
      END DO

      ! pack blocks
      ALLOCATE (recv_packed_blocks(para_env%num_pe))
      DO iproc = 1, para_env%num_pe
         IF (iproc /= mepos_plus1 .AND. recv_nelems(iproc) > 0) &
            ALLOCATE (recv_packed_blocks(iproc)%vector(recv_nelems(iproc)))
      END DO

      ALLOCATE (send_packed_blocks(para_env%num_pe))
      DO iproc = 1, para_env%num_pe
         IF (send_nelems(iproc) > 0) &
            ALLOCATE (send_packed_blocks(iproc)%vector(send_nelems(iproc)))
      END DO

      send_nelems(:) = 0
      DO iatom_col = 1, natoms
         DO iatom_row = 1, iatom_col
            IF (atom_map(iatom_row)%iatom <= atom_map(iatom_col)%iatom) THEN
               icol = atom_map(iatom_col)%iatom
               irow = atom_map(iatom_row)%iatom
               transp = .FALSE.
            ELSE
               icol = atom_map(iatom_row)%iatom
               irow = atom_map(iatom_col)%iatom
               transp = .TRUE.
            END IF

            iproc = rank_contact(irow, icol)
            IF (iproc > 0) THEN
               CALL dbcsr_get_block_p(matrix=matrix_device, &
                                      row=atom_list(iatom_row), col=atom_list(iatom_col), &
                                      block=rblock, found=found)
               IF (found) THEN
                  offset = send_nelems(iproc)
                  n1 = SIZE(rblock, 1)
                  n2 = SIZE(rblock, 2)

                  IF (transp) THEN
                     DO i1 = 1, n1
                        DO i2 = 1, n2
                           send_packed_blocks(iproc)%vector(offset + i2) = rblock(i1, i2)
                        END DO
                        offset = offset + n2
                     END DO
                  ELSE
                     DO i2 = 1, n2
                        DO i1 = 1, n1
                           send_packed_blocks(iproc)%vector(offset + i1) = rblock(i1, i2)
                        END DO
                        offset = offset + n1
                     END DO
                  END IF

                  send_nelems(iproc) = offset
               END IF
            END IF
         END DO
      END DO

      ! send blocks
      ALLOCATE (recv_handlers(para_env%num_pe), send_handlers(para_env%num_pe))

      DO iproc = 1, para_env%num_pe
         IF (iproc /= mepos_plus1 .AND. send_nelems(iproc) > 0) THEN
            CALL para_env%isend(send_packed_blocks(iproc)%vector, iproc - 1, send_handlers(iproc), 1)
         END IF
      END DO

      ! receive blocks
      DO iproc = 1, para_env%num_pe
         IF (iproc /= mepos_plus1) THEN
            IF (recv_nelems(iproc) > 0) THEN
               CALL para_env%irecv(recv_packed_blocks(iproc)%vector, iproc - 1, recv_handlers(iproc), 1)
            END IF
         ELSE
            IF (ALLOCATED(send_packed_blocks(iproc)%vector)) &
               CALL MOVE_ALLOC(send_packed_blocks(iproc)%vector, recv_packed_blocks(iproc)%vector)
         END IF
      END DO

      ! unpack blocks
      DO iproc = 1, para_env%num_pe
         IF (iproc /= mepos_plus1 .AND. recv_nelems(iproc) > 0) &
            CALL recv_handlers(iproc)%wait()
      END DO

      recv_nelems(:) = 0
      DO iatom_col = 1, natoms
         DO iatom_row = 1, iatom_col
            IF (atom_map(iatom_row)%iatom <= atom_map(iatom_col)%iatom) THEN
               icol = atom_map(iatom_col)%iatom
               irow = atom_map(iatom_row)%iatom
            ELSE
               icol = atom_map(iatom_row)%iatom
               irow = atom_map(iatom_col)%iatom
            END IF

            iproc = rank_device(irow, icol)
            IF (iproc > 0) THEN
               CALL dbcsr_get_block_p(matrix=matrix_contact, row=irow, col=icol, block=rblock, found=found)

               IF (found) THEN
                  offset = recv_nelems(iproc)
                  n1 = SIZE(rblock, 1)
                  n2 = SIZE(rblock, 2)

                  DO i2 = 1, n2
                     DO i1 = 1, n1
                        rblock(i1, i2) = recv_packed_blocks(iproc)%vector(offset + i1)
                     END DO
                     offset = offset + n1
                  END DO

                  recv_nelems(iproc) = offset
               END IF
            END IF
         END DO
      END DO

      DO iproc = 1, para_env%num_pe
         IF (iproc /= mepos_plus1 .AND. send_nelems(iproc) > 0) &
            CALL send_handlers(iproc)%wait()
      END DO

      ! release memory
      DEALLOCATE (recv_handlers, send_handlers)

      DO iproc = para_env%num_pe, 1, -1
         IF (ALLOCATED(send_packed_blocks(iproc)%vector)) &
            DEALLOCATE (send_packed_blocks(iproc)%vector)
      END DO
      DEALLOCATE (send_packed_blocks)

      DO iproc = para_env%num_pe, 1, -1
         IF (ALLOCATED(recv_packed_blocks(iproc)%vector)) &
            DEALLOCATE (recv_packed_blocks(iproc)%vector)
      END DO
      DEALLOCATE (recv_packed_blocks)

      DEALLOCATE (rank_contact, rank_device)
      CALL timestop(handle)
   END SUBROUTINE negf_reference_contact_matrix

! **************************************************************************************************
!> \brief Invert cell_to_index mapping between unit cells and DBCSR matrix images.
!> \param cell_to_index  mapping: unit_cell -> image_index
!> \param nimages        number of images
!> \param index_to_cell  inverted mapping: image_index -> unit_cell
!> \par History
!>   * 10.2017 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE invert_cell_to_index(cell_to_index, nimages, index_to_cell)
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER, INTENT(in)                                :: nimages
      INTEGER, DIMENSION(3, nimages), INTENT(out)        :: index_to_cell

      CHARACTER(LEN=*), PARAMETER :: routineN = 'invert_cell_to_index'

      INTEGER                                            :: handle, i1, i2, i3, image
      INTEGER, DIMENSION(3)                              :: lbounds, ubounds

      CALL timeset(routineN, handle)

      index_to_cell(:, :) = 0
      lbounds = LBOUND(cell_to_index)
      ubounds = UBOUND(cell_to_index)

      DO i3 = lbounds(3), ubounds(3) ! z
         DO i2 = lbounds(2), ubounds(2) ! y
            DO i1 = lbounds(1), ubounds(1) ! x
               image = cell_to_index(i1, i2, i3)
               IF (image > 0 .AND. image <= nimages) THEN
                  index_to_cell(1, image) = i1
                  index_to_cell(2, image) = i2
                  index_to_cell(3, image) = i3
               END IF
            END DO
         END DO
      END DO

      CALL timestop(handle)
   END SUBROUTINE invert_cell_to_index

! **************************************************************************************************
!> \brief Helper routine to obtain index of a DBCSR matrix image by its unit cell replica.
!>        Can be used with any usin cell.
!> \param cell           indices of the unit cell
!> \param cell_to_index  mapping: unit_cell -> image_index
!> \return DBCSR matrix images
!>                       (0 means there are no non-zero matrix elements in the image)
!> \par History
!>   * 10.2017 created [Sergey Chulkov]
! **************************************************************************************************
   PURE FUNCTION get_index_by_cell(cell, cell_to_index) RESULT(image)
      INTEGER, DIMENSION(3), INTENT(in)                  :: cell
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER                                            :: image

      IF (LBOUND(cell_to_index, 1) <= cell(1) .AND. UBOUND(cell_to_index, 1) >= cell(1) .AND. &
          LBOUND(cell_to_index, 2) <= cell(2) .AND. UBOUND(cell_to_index, 2) >= cell(2) .AND. &
          LBOUND(cell_to_index, 3) <= cell(3) .AND. UBOUND(cell_to_index, 3) >= cell(3)) THEN

         image = cell_to_index(cell(1), cell(2), cell(3))
      ELSE
         image = 0
      END IF
   END FUNCTION get_index_by_cell

! **************************************************************************************************
!> \brief Desymmetrizes the KS or S matrices for one of spin components
!> \param mat Hamiltonian or overlap matrices
!> \param mat_nosym Desymmetrized Hamiltonian or overlap matrices
!> \param cell_to_index Mapping of cell indices to linear RS indices
!> \param index_to_cell Mapping of linear RS indices to cell indices
!> \param kpoints Kpoint environment
!> \par History
!>   * 05.2020 created [Fabian Ducry]
!>   * 11.2025 Modified for one spin component. [Dmitry Ryndyk]
!> \author Fabian Ducry
! **************************************************************************************************
   SUBROUTINE desymmetrize_matrix(mat, mat_nosym, cell_to_index, index_to_cell, kpoints)
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         POINTER                                         :: mat
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT), &
         POINTER                                         :: mat_nosym
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: cell_to_index
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: index_to_cell
      TYPE(kpoint_type), INTENT(IN), POINTER             :: kpoints

      CHARACTER(len=*), PARAMETER :: routineN = 'desymmetrize_matrix'

      INTEGER                                            :: handle, iatom, ic, icn, icol, irow, &
                                                            jatom, ncell, nomirror, nx, ny, nz
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:, :), POINTER                  :: i2c
      INTEGER, DIMENSION(:, :, :), POINTER               :: c2i
      LOGICAL                                            :: found, lwtr
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CALL timeset(routineN, handle)

      i2c => kpoints%index_to_cell
      c2i => kpoints%cell_to_index

      ncell = SIZE(i2c, 2)

      nx = MAX(ABS(LBOUND(c2i, 1)), ABS(UBOUND(c2i, 1)))
      ny = MAX(ABS(LBOUND(c2i, 2)), ABS(UBOUND(c2i, 3)))
      nz = MAX(ABS(LBOUND(c2i, 3)), ABS(UBOUND(c2i, 3)))
      ALLOCATE (cell_to_index(-nx:nx, -ny:ny, -nz:nz))
      cell_to_index(LBOUND(c2i, 1):UBOUND(c2i, 1), &
                    LBOUND(c2i, 2):UBOUND(c2i, 2), &
                    LBOUND(c2i, 3):UBOUND(c2i, 3)) = c2i

      ! identify cells with no mirror img
      nomirror = 0
      DO ic = 1, ncell
         cell = i2c(:, ic)
         IF (cell_to_index(-cell(1), -cell(2), -cell(3)) == 0) &
            nomirror = nomirror + 1
      END DO

      ! create the mirror imgs
      ALLOCATE (index_to_cell(3, ncell + nomirror))
      index_to_cell(:, 1:ncell) = i2c

      nomirror = 0 ! count the imgs without mirror
      DO ic = 1, ncell
         cell = index_to_cell(:, ic)
         IF (cell_to_index(-cell(1), -cell(2), -cell(3)) == 0) THEN
            nomirror = nomirror + 1
            index_to_cell(:, ncell + nomirror) = -cell
            cell_to_index(-cell(1), -cell(2), -cell(3)) = ncell + nomirror
         END IF
      END DO
      ncell = ncell + nomirror

      CALL get_kpoint_info(kpoints, sab_nl=sab_nl)
      ! allocate the nonsymmetric matrices
      NULLIFY (mat_nosym)
      CALL dbcsr_allocate_matrix_set(mat_nosym, ncell)
      DO ic = 1, ncell
         ALLOCATE (mat_nosym(ic)%matrix)
         CALL dbcsr_create(matrix=mat_nosym(ic)%matrix, &
                           template=mat(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_alloc_block_from_nbl(mat_nosym(ic)%matrix, &
                                            sab_nl, desymmetrize=.TRUE.)
         CALL dbcsr_set(mat_nosym(ic)%matrix, 0.0_dp)
      END DO

      ! desymmetrize the matrix for real space printing
      CALL neighbor_list_iterator_create(nl_iterator, sab_nl)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, cell=cell)

         ic = cell_to_index(cell(1), cell(2), cell(3))
         icn = cell_to_index(-cell(1), -cell(2), -cell(3))
         CPASSERT(icn > 0)

         irow = iatom
         icol = jatom
         lwtr = .FALSE.
         ! always copy from the top
         IF (iatom > jatom) THEN
            irow = jatom
            icol = iatom
            lwtr = .TRUE.
         END IF

         CALL dbcsr_get_block_p(matrix=mat(ic)%matrix, &
                                row=irow, col=icol, block=block, found=found)
         CPASSERT(found)

         ! copy to M(R) at (iatom,jatom)
         ! copy to M(-R) at (jatom,iatom)
         IF (lwtr) THEN
            CALL dbcsr_put_block(matrix=mat_nosym(ic)%matrix, &
                                 row=iatom, col=jatom, block=TRANSPOSE(block))
            CALL dbcsr_put_block(matrix=mat_nosym(icn)%matrix, &
                                 row=jatom, col=iatom, block=block)
         ELSE
            CALL dbcsr_put_block(matrix=mat_nosym(ic)%matrix, &
                                 row=iatom, col=jatom, block=block)
            CALL dbcsr_put_block(matrix=mat_nosym(icn)%matrix, &
                                 row=jatom, col=iatom, block=TRANSPOSE(block))
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO ic = 1, ncell
         CALL dbcsr_finalize(mat_nosym(ic)%matrix)
      END DO

      CALL timestop(handle)

   END SUBROUTINE desymmetrize_matrix

END MODULE negf_matrix_utils
