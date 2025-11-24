!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Implements transformations from k-space to R-space for Fortran array matrices
!> \note This code is less performant/more memory consuming than the methods in kpoint_methods.F
!>       Use only when transformations are not the computational bottleneck.
!> \par History
!>       2025.11 created [Stepan Marek]
!> \author Stepan Marek
! **************************************************************************************************
MODULE kpoint_k_r_trafo_simple
   USE kinds,                             ONLY: dp
   USE mathconstants,                     ONLY: twopi
   USE cp_dbcsr_api,                      ONLY: dbcsr_p_type, &
                                                dbcsr_type, &
                                                dbcsr_get_info, &
                                                dbcsr_get_readonly_block_p, &
                                                dbcsr_type_no_symmetry
   USE message_passing,                   ONLY: mp_comm_type
   USE qs_neighbor_list_types,            ONLY: neighbor_list_set_p_type, &
                                                neighbor_list_iterator_p_type, &
                                                get_iterator_info, &
                                                neighbor_list_iterate, &
                                                neighbor_list_iterator_create, &
                                                neighbor_list_iterator_release
   USE kpoint_types,                      ONLY: kpoint_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_k_r_trafo_simple'

   PUBLIC :: replicate_rs_matrices, &
      rs_to_kp

CONTAINS
! **************************************************************************************************
!> \brief Convert dbcsr matrices representing operators in real-space image cells to arrays
!> \param rs_dbcsr_in Array of dbcsr matrices
!> \param kpoint_in The kpoint environment of the source matrix (providing neighbor list and cell_to_index)
!> \param rs_array_out Multidimensional array - matrices are duplicated on each MPI rank
!> \param cell_to_index_out Cell to index array for the destination array
!> \author Stepan Marek
!> \date 11.2025
! **************************************************************************************************
   SUBROUTINE replicate_rs_matrices(rs_dbcsr_in, kpoint_in, rs_array_out, cell_to_index_out)
      ! dimension 1 : spin, dimension 2 : image cell index
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rs_dbcsr_in
      TYPE(kpoint_type), POINTER                         :: kpoint_in
      ! dimension 1 : spin/derivs, dimension 2 : image cell index, dimension 3 : rows, dimension 4 : cols
      REAL(kind=dp), DIMENSION(:,:,:,:)                  :: rs_array_out
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index_out
      CHARACTER(len=*), PARAMETER                        :: routineN = "replicate_rs_matrices"
      INTEGER                                            :: n_spin, ispin, &
                                                            handle, &
                                                            !i, j, &
                                                            !row_offset, col_offset, &
                                                            iatom, jatom, &
                                                            src_index!src_index, &
                                                            !col_size, row_size
      ! INTEGER, DIMENSION(:), POINTER                     :: row_offsets, col_offsets
      TYPE(mp_comm_type)                                 :: group
      ! REAL(kind=dp), DIMENSION(:,:), POINTER             :: block
      INTEGER, DIMENSION(3)                              :: cell
      ! LOGICAL                                            :: found
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_kp_src
      CHARACTER                                          :: matrix_sym

      CALL timeset(routineN, handle)

      IF (SIZE(rs_dbcsr_in, 2) < 1) THEN
         CALL cp_abort(__LOCATION__, "No source image cells provided!")
      END IF
      ! Start by constructing the cell_to_index_src
      sab_kp_src => kpoint_in%sab_nl
      ! NOTE : The first index in matrix_s_kp is not spin index, but number of derivatives.
      !        But, for matrix_ks_kp, this is indeed the spin index.
      n_spin = SIZE(rs_dbcsr_in, 1)
      CALL dbcsr_get_info(rs_dbcsr_in(1, 1)%matrix, group=group, matrix_type=matrix_sym)
      DO ispin=1,n_spin
         CALL neighbor_list_iterator_create(iterator, sab_kp_src)
         DO WHILE (neighbor_list_iterate(iterator) == 0)
            CALL get_iterator_info(iterator, cell=cell, iatom=iatom, jatom=jatom)
            src_index = kpoint_in%cell_to_index(cell(1), cell(2), cell(3))
            IF (src_index == 0) THEN
               CALL cp_abort(__LOCATION__, "Image not found in the source array.")
            END IF
            ! NOTE : Expect only specific symmetry storage relevant for kpoint calculations
            IF (matrix_sym == dbcsr_type_no_symmetry) THEN
               CALL write_block_no_sym(iatom, jatom, cell, rs_dbcsr_in(ispin, src_index)%matrix, &
                                       rs_array_out(ispin, :, :, :), cell_to_index_out)
            ELSE
               CALL write_block_symmetric(iatom, jatom, cell, rs_dbcsr_in(ispin, src_index)%matrix, &
                                          rs_array_out(ispin, :, :, :), cell_to_index_out)
            END IF
         END DO
         CALL neighbor_list_iterator_release(iterator)
      END DO
      CALL group%sum(rs_array_out(:,:,:,:))
      CALL timestop(handle)
   END SUBROUTINE replicate_rs_matrices
! **************************************************************************************************
!> \brief Write a single block from the dbcsr matrix to correct place, with assumed symmetric dbcsr
!> \param iatom first atom index
!> \param jatom second atom index
!> \param cell Current cell (of second atom)
!> \param matrix_in DBCSR matrix input, symmetric assumed (A_(mu nu)^(R) = A_(nu mu)^(-R))
!> \param array_out Multidimensional array - dim 1 : rs_index, dim 2 : rows, dim 3 : cols
!> \param cell_to_index_out Mapping of cell coords to rs_indices
!> \author Stepan Marek
!> \date 11.2025
! **************************************************************************************************
   SUBROUTINE write_block_symmetric(iatom, jatom, cell, matrix_in, array_out, cell_to_index_out)
      INTEGER, INTENT(IN)                                :: iatom, &
                                                            jatom
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_in
      REAL(kind=dp), DIMENSION(:,:,:), INTENT(INOUT)     :: array_out
      INTEGER, DIMENSION(:,:,:), POINTER, INTENT(IN)     :: cell_to_index_out
      REAL(kind=dp), DIMENSION(:,:), POINTER             :: block
      LOGICAL                                            :: found
      INTEGER, DIMENSION(:), POINTER                     :: row_offsets, &
                                                            col_offsets
      INTEGER                                            :: row_size, &
                                                            col_size, &
                                                            dest_index, &
                                                            dest_index_t, &
                                                            row_offset, &
                                                            col_offset, &
                                                            i, j, i_g, j_g

      CALL dbcsr_get_info(matrix_in, row_blk_offset=row_offsets, col_blk_offset=col_offsets)
      IF (iatom > jatom) THEN
         CALL dbcsr_get_readonly_block_p(matrix_in, row=jatom, col=iatom, block=block, &
                                row_size=row_size, col_size=col_size, found=found)
         IF (.NOT. found) RETURN
         ! Block found, prepare for write
         dest_index = cell_to_index_out(-cell(1), -cell(2), -cell(3))
         IF (dest_index == 0) CPABORT("Mirror image index not present.")
         dest_index_t = cell_to_index_out(cell(1), cell(2), cell(3))
         IF (dest_index_t == 0) CPABORT("Image index not present.")
         row_offset = row_offsets(jatom)
         col_offset = col_offsets(iatom)
      ELSE
         CALL dbcsr_get_readonly_block_p(matrix_in, row=iatom, col=jatom, block=block, &
                                row_size=row_size, col_size=col_size, found=found)
         IF (.NOT. found) RETURN
         ! Block found, prepare for write
         dest_index = cell_to_index_out(cell(1), cell(2), cell(3))
         IF (dest_index == 0) CPABORT("Image index not present.")
         dest_index_t = cell_to_index_out(-cell(1), -cell(2), -cell(3))
         IF (dest_index_t == 0) CPABORT("Mirror image index not present.")
         row_offset = row_offsets(iatom)
         col_offset = col_offsets(jatom)
      END IF
      ! Do the write
      !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i,j,i_g,j_g) &
      !$OMP SHARED(row_size, col_size, row_offset, col_offset, array_out, &
      !$OMP dest_index, dest_index_t, block, iatom, jatom)
      DO i=1,row_size
         i_g = i + row_offset - 1
         DO j=1, col_size
            j_g = j + col_offset - 1
            array_out(dest_index, i_g, j_g) = block(i, j)
            IF (iatom /= jatom) array_out(dest_index_t, j_g, i_g) = block(i, j)
         END DO
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE write_block_symmetric
! **************************************************************************************************
!> \brief Write a single block from the dbcsr matrix to correct place, assuming no symmetry dbcsr
!> \param iatom first atom index
!> \param jatom second atom index
!> \param cell Current cell (of second atom)
!> \param matrix_in DBCSR matrix input, symmetric assumed (A_(mu nu)^(R) = A_(nu mu)^(-R))
!> \param array_out Multidimensional array - dim 1 : rs_index, dim 2 : rows, dim 3 : cols
!> \param cell_to_index_out Mapping of cell coords to rs_indices
!> \author Stepan Marek
!> \date 11.2025
! **************************************************************************************************
   SUBROUTINE write_block_no_sym(iatom, jatom, cell, matrix_in, array_out, cell_to_index_out)
      INTEGER, INTENT(IN)                                :: iatom, &
                                                            jatom
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_in
      REAL(kind=dp), DIMENSION(:,:,:), INTENT(INOUT)     :: array_out
      INTEGER, DIMENSION(:,:,:), POINTER, INTENT(IN)     :: cell_to_index_out
      REAL(kind=dp), DIMENSION(:,:), POINTER             :: block
      LOGICAL                                            :: found
      INTEGER, DIMENSION(:), POINTER                     :: row_offsets, &
                                                            col_offsets
      INTEGER                                            :: row_size, &
                                                            col_size, &
                                                            dest_index, &
                                                            row_offset, &
                                                            col_offset, &
                                                            i, j, i_g, j_g

      dest_index = cell_to_index_out(cell(1), cell(2), cell(3))
      IF (dest_index == 0) CPABORT("Image index not present.")
      CALL dbcsr_get_info(matrix_in, row_blk_offset=row_offsets, col_blk_offset=col_offsets)
      row_offset = row_offsets(iatom)
      col_offset = col_offsets(jatom)
      CALL dbcsr_get_readonly_block_p(matrix_in, row=iatom, col=jatom, block=block, found=found, &
                                      row_size=row_size, col_size=col_size)
      !$OMP PARALLEL DO DEFAULT(none) PRIVATE(i, j, i_g, j_g) &
      !$OMP SHARED(row_size, col_size, row_offset, col_offset, dest_index, block, array_out)
      DO i=1, row_size
         i_g = i + row_offset - 1
         DO j=1, col_size
            j_g = j + col_offset - 1
            array_out(dest_index, i_g, j_g) = block(i, j)
         END DO
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE write_block_no_sym
! **************************************************************************************************
!> \brief Integrate RS matrices (stored as Fortran array) into a kpoint matrix at given kp
!> \param ks_real Array of resulting transformed 
!> \param rs_real Multidimensional array of real parts of the matrix
!> \param array Multidimensional array - matrices are duplicated on each MPI rank
!> \param index_to_cell Gets the image cell coordinates from the rs_dbcsr index
!> \param cell_to_index Gets the rs_dbcsr index from image cell coordinates
!> \param deriv_direction Derivative direction - vector to take scalar product with the gradient
!> \author Stepan Marek
!> \date 11.2025
! **************************************************************************************************
   SUBROUTINE rs_to_kp(rs_real, ks_complex, index_to_cell, xkp, deriv_direction, hmat)
      ! index 1 : rows, index 2 : columns of the matrix, Index 3 : image cell index
      REAL(kind=dp), DIMENSION(:,:,:), INTENT(IN)        :: rs_real
      ! Index 1 : rows, Index 2 : columns of the matrix (AO or MO indices)
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(OUT)      :: ks_complex
      ! Index 1 : direction, index 2 : image cell index
      INTEGER, DIMENSION(:,:), POINTER                   :: index_to_cell
      REAL(kind=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, INTENT(IN), OPTIONAL                      :: deriv_direction
      REAL(kind=dp), DIMENSION(3,3), INTENT(IN), OPTIONAL:: hmat
      CHARACTER(len=*), PARAMETER                        :: routineN="rs_to_kp"
      INTEGER                                            :: handle, i, &
                                                            n_images
                                                            ! row_size, &
                                                            ! col_size, &
                                                            ! n_images, &
                                                            ! i, j, k
      ! REAL(kind=dp)                                      :: cosfactor, &
      !                                                       sinfactor
      ! LOGICAL                                            :: derivs_required

      CALL timeset(routineN, handle)
      ! Get matrix constants
      n_images = SIZE(rs_real, 1)
      ! Get the required derivatives for the deriv direction
      ! derivs_required = .FALSE.
      ! IF (PRESENT(deriv_direction)) THEN
      !    derivs_required = (deriv_direction /= 0)
      !    IF (.NOT. PRESENT(hmat)) THEN
      !       CALL cp_abort(__LOCATION__, "derivative requested but h matrix not provided")
      !    END IF
      ! END IF
      ! Now, iterate over realspace and build the sum
      ! TODO : OMP
      ks_complex(:,:) = CMPLX(0.0, 0.0, kind=dp)
      DO i=1,n_images
         CALL add_rs_to_kp(ks_complex, rs_real, xkp, i, index_to_cell, deriv_direction, hmat)
         ! ! Construct the phase and factors
         ! deriv_factor = 1.0_dp
         ! IF (derivs_required) THEN
         !    cell_vectors = MATMUL(hmat, REAL(index_to_cell(1:3, i), dp))
         !    deriv_factor = cell_vectors(deriv_direction)
         ! END IF
         ! cosfactor = COS(twopi * SUM(xkp * index_to_cell(:,i)))
         ! sinfactor = SIN(twopi * SUM(xkp * index_to_cell(:,i)))
         ! DO j=1,row_size
         !    DO k=1,col_size
         !       ks_real(j,k) = ks_real(j,k) + deriv_factor * cosfactor * rs_real(i, j, k)
         !       ks_imag(j,k) = ks_imag(j,k) + deriv_factor * sinfactor * rs_real(i, j, k)
         !    END DO
         ! END DO
      END DO
      CALL timestop(handle)
   END SUBROUTINE rs_to_kp
   ! Extra routine, reusing the same code for fm-type transform
   SUBROUTINE add_rs_to_kp(ks_array_out, rs_array_in, xkp, imindex, index_to_cell, deriv_direction, hmat)
      COMPLEX(kind=dp), DIMENSION(:,:), INTENT(INOUT)    :: ks_array_out
      REAL(kind=dp), DIMENSION(:,:,:), INTENT(IN)        :: rs_array_in
      REAL(kind=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, INTENT(IN)                                :: imindex
      INTEGER, DIMENSION(:,:), INTENT(IN)                :: index_to_cell
      INTEGER, INTENT(IN), OPTIONAL                      :: deriv_direction
      REAL(kind=dp), DIMENSION(3,3), &
         INTENT(IN), OPTIONAL                            :: hmat
      REAL(kind=dp)                                      :: deriv_factor
      REAL(kind=dp), DIMENSION(3)                        :: cell_vector
      COMPLEX(kind=dp)                                   :: phase_factor
      INTEGER                                            :: i, j

      IF (PRESENT(deriv_direction) .AND. (.NOT. PRESENT(hmat))) THEN
         CALL cp_abort(__LOCATION__, "Deriv. direction given but no hmat provided")
      END IF

      deriv_factor = 1.0_dp
      IF (PRESENT(deriv_direction)) THEN
         cell_vector = MATMUL(hmat, index_to_cell(1:3,imindex))
         deriv_factor = cell_vector(deriv_direction)
      END IF

      phase_factor = CMPLX(COS(twopi*SUM(xkp(:) * index_to_cell(:, imindex))), &
                           SIN(twopi*SUM(xkp(:) * index_to_cell(:, imindex))), kind=dp)

      DO i=1, SIZE(ks_array_out, 1)
         DO j=1, SIZE(ks_array_out, 2)
            ks_array_out(i,j) = ks_array_out(i, j) + &
                                deriv_factor * phase_factor * rs_array_in(i, j, imindex)
         END DO
      END DO
   END SUBROUTINE add_rs_to_kp
END MODULE kpoint_k_r_trafo_simple
