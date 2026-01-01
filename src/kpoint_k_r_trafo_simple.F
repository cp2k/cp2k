!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Implements transformations from k-space to R-space for Fortran array matrices
!> \note This code is less performant/more memory consuming than the methods in kpoint_methods.F
!>       Use only when transformations are not the computational bottleneck.
!> \par History
!>       11.2025 created [Stepan Marek]
!> \author Stepan Marek
! **************************************************************************************************
MODULE kpoint_k_r_trafo_simple
   USE cp_cfm_types,                    ONLY: cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_info,&
                                              dbcsr_get_readonly_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE mathconstants,                   ONLY: twopi,&
                                              z_zero
   USE message_passing,                 ONLY: mp_comm_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_k_r_trafo_simple'

   PUBLIC :: replicate_rs_matrices, &
             rs_to_kp, &
             fm_rs_to_kp, &
             add_kp_to_all_rs, &
             fm_add_kp_to_all_rs, &
             add_kp_to_rs

CONTAINS
! **************************************************************************************************
!> \brief Convert dbcsr matrices representing operators in real-space image cells to arrays
!> \param rs_dbcsr_in Array of dbcsr matrices
!> \param kpoint_in The kpoint environment of the source matrix (providing neighbor list and cell_to_index)
!> \param rs_array_out Multidimensional array - matrices are duplicated on each MPI rank
!>                     dim  1 : spin, dim 2 : rows, dim 3 : cols, dim 4 : image index
!> \param cell_to_index_out Cell to index array for the destination array
!> \date 11.2025
!> \author Stepan Marek
! **************************************************************************************************
   SUBROUTINE replicate_rs_matrices(rs_dbcsr_in, kpoint_in, rs_array_out, cell_to_index_out)
      ! dimension 1 : spin, dimension 2 : image cell index
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rs_dbcsr_in
      TYPE(kpoint_type), POINTER                         :: kpoint_in
      REAL(kind=dp), DIMENSION(:, :, :, :)               :: rs_array_out
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index_out

      CHARACTER(len=*), PARAMETER :: routineN = 'replicate_rs_matrices'

      CHARACTER                                          :: matrix_sym
      INTEGER                                            :: handle, iatom, ispin, jatom, n_spin, &
                                                            src_index
      INTEGER, DIMENSION(3)                              :: cell
      TYPE(mp_comm_type)                                 :: group
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_kp_src

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
      DO ispin = 1, n_spin
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
      CALL group%sum(rs_array_out(:, :, :, :))
      CALL timestop(handle)
   END SUBROUTINE replicate_rs_matrices
! **************************************************************************************************
!> \brief Write a single block from the dbcsr matrix to correct place, with assumed symmetric dbcsr
!> \param iatom first atom index
!> \param jatom second atom index
!> \param cell Current cell (of second atom)
!> \param matrix_in DBCSR matrix input, symmetric assumed (A_(mu nu)^(R) = A_(nu mu)^(-R))
!> \param array_out Multidimensional array - dim 1 : rows, dim 2 : cols, dim 3 : rs_index
!> \param cell_to_index_out Mapping of cell coords to rs_indices
!> \date 11.2025
!> \author Stepan Marek
! **************************************************************************************************
   SUBROUTINE write_block_symmetric(iatom, jatom, cell, matrix_in, array_out, cell_to_index_out)
      INTEGER, INTENT(IN)                                :: iatom, jatom
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_in
      REAL(kind=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: array_out
      INTEGER, DIMENSION(:, :, :), INTENT(IN), POINTER   :: cell_to_index_out

      INTEGER                                            :: col_offset, col_size, dest_index, &
                                                            dest_index_t, i, i_g, j, j_g, &
                                                            row_offset, row_size
      INTEGER, DIMENSION(:), POINTER                     :: col_offsets, row_offsets
      LOGICAL                                            :: found
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: block

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
      !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,i_g,j_g) &
      !$OMP SHARED(row_size, col_size, row_offset, col_offset, array_out, &
      !$OMP dest_index, dest_index_t, block, iatom, jatom)
      DO i = 1, row_size
         i_g = i + row_offset - 1
         DO j = 1, col_size
            j_g = j + col_offset - 1
            array_out(i_g, j_g, dest_index) = block(i, j)
            IF (iatom /= jatom) array_out(j_g, i_g, dest_index_t) = block(i, j)
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
!> \param array_out Multidimensional array - dim 1 : rows, dim 2 : cols, dim 3 : rs_index
!> \param cell_to_index_out Mapping of cell coords to rs_indices
!> \date 11.2025
!> \author Stepan Marek
! **************************************************************************************************
   SUBROUTINE write_block_no_sym(iatom, jatom, cell, matrix_in, array_out, cell_to_index_out)
      INTEGER, INTENT(IN)                                :: iatom, jatom
      INTEGER, DIMENSION(3), INTENT(IN)                  :: cell
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_in
      REAL(kind=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: array_out
      INTEGER, DIMENSION(:, :, :), INTENT(IN), POINTER   :: cell_to_index_out

      INTEGER                                            :: col_offset, col_size, dest_index, i, &
                                                            i_g, j, j_g, row_offset, row_size
      INTEGER, DIMENSION(:), POINTER                     :: col_offsets, row_offsets
      LOGICAL                                            :: found
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: block

      dest_index = cell_to_index_out(cell(1), cell(2), cell(3))
      IF (dest_index == 0) CPABORT("Image index not present.")
      CALL dbcsr_get_info(matrix_in, row_blk_offset=row_offsets, col_blk_offset=col_offsets)
      row_offset = row_offsets(iatom)
      col_offset = col_offsets(jatom)
      CALL dbcsr_get_readonly_block_p(matrix_in, row=iatom, col=jatom, block=block, found=found, &
                                      row_size=row_size, col_size=col_size)
      !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, i_g, j_g) &
      !$OMP SHARED(row_size, col_size, row_offset, col_offset, dest_index, block, array_out)
      DO i = 1, row_size
         i_g = i + row_offset - 1
         DO j = 1, col_size
            j_g = j + col_offset - 1
            array_out(i_g, j_g, dest_index) = block(i, j)
         END DO
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE write_block_no_sym
! **************************************************************************************************
!> \brief Integrate RS matrices (stored as Fortran array) into a kpoint matrix at given kp
!> \param rs_real Multidimensional array of real parts of the matrix, dim 1 : rows, dim 2 : cols, dim 3 : image index
!> \param ks_complex Target complex k-space matrix
!> \param index_to_cell Gets the image cell coordinates from the rs_dbcsr index
!> \param xkp Single kpoint where the transformation is evalued at
!> \param deriv_direction Derivative direction - indicates along which cartesian direction to take the derivative
!> \param hmat Cell size matrix, required for the derivative
!> \date 11.2025
!> \author Stepan Marek
! **************************************************************************************************
   SUBROUTINE rs_to_kp(rs_real, ks_complex, index_to_cell, xkp, deriv_direction, hmat)
      REAL(kind=dp), DIMENSION(:, :, :), INTENT(IN)      :: rs_real
      COMPLEX(kind=dp), DIMENSION(:, :), INTENT(OUT)     :: ks_complex
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      REAL(kind=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, INTENT(IN), OPTIONAL                      :: deriv_direction
      REAL(kind=dp), DIMENSION(3, 3), INTENT(IN), &
         OPTIONAL                                        :: hmat

      CHARACTER(len=*), PARAMETER                        :: routineN = 'rs_to_kp'

      INTEGER                                            :: handle, i, n_images

      CALL timeset(routineN, handle)
      ! Get matrix constants
      n_images = SIZE(rs_real, 3)
      ! Get the required derivatives for the deriv direction
      IF (PRESENT(deriv_direction)) THEN
         IF (.NOT. PRESENT(hmat) .AND. deriv_direction /= 0) THEN
            CALL cp_abort(__LOCATION__, "Derivative requested but h matrix not provided!")
         END IF
      END IF
      ! Now, iterate over realspace and build the sum
      ks_complex(:, :) = CMPLX(0.0, 0.0, kind=dp)
      DO i = 1, n_images
         CALL add_rs_to_kp(ks_complex, rs_real(:, :, i), xkp, i, index_to_cell, deriv_direction, hmat)
      END DO
      CALL timestop(handle)
   END SUBROUTINE rs_to_kp
! **************************************************************************************************
!> \brief Transforms array of fm RS matrices into cfm k-space matrix, at given kpoint index
!> \param cfm_kp The resulting complex fm matrix, containing the k-space matrix elements
!> \param fm_rs The real space matrix array
!> \param kpoints Structure containing information about number, positions, weights etc. of kpoints
!> \param ikp Index of the target kpoint at which the transformation is evaluated
!> \par History
!>    05.2024 Created [Jan Wilhelm]
!>    11.2025 Moved to general file [Stepan Marek]
! **************************************************************************************************
   SUBROUTINE fm_rs_to_kp(cfm_kp, fm_rs, kpoints, ikp)
      TYPE(cp_cfm_type)                                  :: cfm_kp
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_rs
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fm_rs_to_kp'

      INTEGER                                            :: handle, img, nimages, nimages_fm_rs

      CALL timeset(routineN, handle)

      nimages = SIZE(kpoints%index_to_cell, 2)
      nimages_fm_rs = SIZE(fm_rs)

      IF (nimages /= nimages_fm_rs) CALL cp_abort(__LOCATION__, "Index to cell and provided fm "// &
                                                  "array are inconsistent.")

      cfm_kp%local_data(:, :) = z_zero
      DO img = 1, nimages

         CALL add_rs_to_kp(cfm_kp%local_data, fm_rs(img)%local_data, kpoints%xkp(1:3, ikp), &
                           img, kpoints%index_to_cell)

      END DO

      CALL timestop(handle)
   END SUBROUTINE fm_rs_to_kp
! **************************************************************************************************
!> \brief Adds a single RS array with correct phase factor to the target k-space argument
!> \param ks_array_out The array representing (part) of the complex k-space matrix
!> \param rs_array_in The RS input array
!> \param xkp The kpoint coordinates, in units of unit cell
!> \param imindex The RS index, used for index to cell determination
!> \param index_to_cell Getting cell coordinates from the imindex
!> \param deriv_direction The cartesian direction of the derivative
!> \param hmat The unit cell dimensions
!> \par History
!>    05.2024 Created [Jan Wilhelm]
!>    11.2025 Added k-derivative option [Shridhar Shanbhag]
!>    11.2025 Moved to general file [Stepan Marek]
!> \note Always produces non-symmetric matrix, given a non-symmetric RS array
! **************************************************************************************************
   SUBROUTINE add_rs_to_kp(ks_array_out, rs_array_in, xkp, imindex, index_to_cell, deriv_direction, hmat)
      COMPLEX(kind=dp), DIMENSION(:, :), INTENT(INOUT)   :: ks_array_out
      REAL(kind=dp), DIMENSION(:, :), INTENT(IN)         :: rs_array_in
      REAL(kind=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, INTENT(IN)                                :: imindex
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: index_to_cell
      INTEGER, INTENT(IN), OPTIONAL                      :: deriv_direction
      REAL(kind=dp), DIMENSION(3, 3), INTENT(IN), &
         OPTIONAL                                        :: hmat

      COMPLEX(kind=dp)                                   :: deriv_factor, phase_factor
      INTEGER                                            :: i, j
      REAL(kind=dp), DIMENSION(3)                        :: cell_vector

      IF (PRESENT(deriv_direction) .AND. (.NOT. PRESENT(hmat))) THEN
         CALL cp_abort(__LOCATION__, "Deriv. direction given but no hmat provided")
      END IF

      deriv_factor = CMPLX(1.0, 0.0, kind=dp)
      IF (PRESENT(deriv_direction)) THEN
         cell_vector = MATMUL(hmat, index_to_cell(1:3, imindex))
         ! Multiply by i R(direction)
         deriv_factor = CMPLX(0.0, cell_vector(deriv_direction), kind=dp)
      END IF

      phase_factor = CMPLX(COS(twopi*SUM(xkp(:)*index_to_cell(:, imindex))), &
                           SIN(twopi*SUM(xkp(:)*index_to_cell(:, imindex))), kind=dp)

      DO i = 1, SIZE(ks_array_out, 1)
         DO j = 1, SIZE(ks_array_out, 2)
            ks_array_out(i, j) = ks_array_out(i, j) + &
                                 deriv_factor*phase_factor*rs_array_in(i, j)
         END DO
      END DO
   END SUBROUTINE add_rs_to_kp
! **************************************************************************************************
!> \brief Adds given kpoint matrix to all rs matrices
!> \param array_kp Input k-space matrix array
!> \param array_rs Output rs arrays
!> \param kpoints Structure containing information about number, positions, weights etc. of kpoints
!> \param ikp Index of the target kpoint at which the transformation is evaluated
!> \param index_to_cell_ext External supplied index_to_cell
!> \par History
!>    05.2024 Created [Jan Wilhelm]
!>    11.2025 Moved to general file [Stepan Marek]
!> \note Part of transform in parallelism over k-points (hence why not for a single rs cell)
! **************************************************************************************************
   SUBROUTINE add_kp_to_all_rs(array_kp, array_rs, kpoints, ikp, index_to_cell_ext)
      COMPLEX(kind=dp), DIMENSION(:, :)                  :: array_kp
      REAL(kind=dp), DIMENSION(:, :, :)                  :: array_rs
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp
      INTEGER, DIMENSION(:, :), OPTIONAL, POINTER        :: index_to_cell_ext

      CHARACTER(len=*), PARAMETER                        :: routineN = 'add_kp_to_all_rs'

      INTEGER                                            :: handle, img
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell

      index_to_cell => kpoints%index_to_cell
      IF (PRESENT(index_to_cell_ext)) index_to_cell => index_to_cell_ext

      CALL timeset(routineN, handle)

      IF (SIZE(array_rs, 3) /= SIZE(index_to_cell, 2)) CALL cp_abort(__LOCATION__, &
                                                                  "The provided RS array and cell_to_index array are inconsistent.")

      DO img = 1, SIZE(array_rs, 3)
         cell(:) = index_to_cell(:, img)
         CALL add_kp_to_rs(array_kp, array_rs(:, :, img), cell, kpoints, ikp)
      END DO

      CALL timestop(handle)
   END SUBROUTINE add_kp_to_all_rs
! **************************************************************************************************
!> \brief Adds given kpoint matrix to a single rs matrix
!> \param cfm_kp The input complex fm matrix, containing the k-space matrix elements
!> \param fm_rs The real space matrix array
!> \param kpoints Structure containing information about number, positions, weights etc. of kpoints
!> \param ikp Index of the target kpoint at which the transformation is evaluated
!> \par History
!>    05.2024 Created [Jan Wilhelm]
!>    11.2025 Moved to general file [Stepan Marek]
!> \note Part of transform in parallelism over k-points (hence why not for a single rs cell)
! **************************************************************************************************
   SUBROUTINE fm_add_kp_to_all_rs(cfm_kp, fm_rs, kpoints, ikp)
      TYPE(cp_cfm_type)                                  :: cfm_kp
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_rs
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp

      CHARACTER(len=*), PARAMETER :: routineN = 'fm_add_kp_to_all_rs'

      INTEGER                                            :: handle, img
      INTEGER, DIMENSION(3)                              :: cell

      CALL timeset(routineN, handle)

      IF (SIZE(fm_rs, 1) /= SIZE(kpoints%index_to_cell, 2)) CALL cp_abort(__LOCATION__, &
                                                                  "The provided RS array and cell_to_index array are inconsistent.")

      DO img = 1, SIZE(fm_rs, 1)
         cell(:) = kpoints%index_to_cell(:, img)
         CALL add_kp_to_rs(cfm_kp%local_data, fm_rs(img)%local_data, cell, kpoints, ikp)
      END DO
      CALL timestop(handle)
   END SUBROUTINE fm_add_kp_to_all_rs
! **************************************************************************************************
!> \brief Adds given kpoint matrix to a single rs matrix
!> \param array_kp The resulting complex fm matrix, containing the k-space matrix elements
!> \param array_rs The real space matrix array
!> \param cell The image cell of the target RS matrix
!> \param kpoints Structure containing information about number, positions, weights etc. of kpoints
!> \param ikp Index of the target kpoint at which the transformation is evaluated
!> \par History
!>    05.2024 Created [Jan Wilhelm]
!>    11.2025 Moved to general file [Stepan Marek]
!> \note Part of transform in parallelism over k-points (hence why not for a single rs cell)
! **************************************************************************************************
   SUBROUTINE add_kp_to_rs(array_kp, array_rs, cell, kpoints, ikp)
      COMPLEX(kind=dp), DIMENSION(:, :)                  :: array_kp
      REAL(kind=dp), DIMENSION(:, :)                     :: array_rs
      INTEGER, DIMENSION(3)                              :: cell
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp

      REAL(kind=dp)                                      :: phase

      ! Determine the phase
      phase = twopi*SUM(kpoints%xkp(:, ikp)*cell(:))
      ! The array is real - it is assumed that imaginary parts cancel each other, so they are not handled
      ! Then e^(-i phi) S(k) = Re(S(k)) cos(phi) - i sin (phi) i Im(S(k)) =
      ! = Re(S(k)) cos(phi) + Im(S(k)) sin(phi)
      array_rs(:, :) = array_rs(:, :) + (REAL(array_kp(:, :))*COS(phase) + &
                                         AIMAG(array_kp(:, :))*SIN(phase))*kpoints%wkp(ikp)
   END SUBROUTINE add_kp_to_rs
! **************************************************************************************************
!> \brief Calculates the inverse transform, assuming all kpoints are present on the given rank
!> \param array_kp The input complex fm matrix, containing the k-space matrix elements
!> \param array_rs The output real space matrix array
!> \param cell The image cell of the target RS matrix
!> \param kpoints Structure containing information about number, positions, weights etc. of kpoints
!> \par History
!>    05.2024 Created [Jan Wilhelm]
!>    11.2025 Moved to general file [Stepan Marek]
!> \note Part of transform in parallelism over k-points
! **************************************************************************************************
   SUBROUTINE kp_to_rs(array_kp, array_rs, cell, kpoints)
      COMPLEX(kind=dp), DIMENSION(:, :, :)               :: array_kp
      REAL(kind=dp), DIMENSION(:, :)                     :: array_rs
      INTEGER, DIMENSION(3)                              :: cell
      TYPE(kpoint_type), POINTER                         :: kpoints

      CHARACTER(len=*), PARAMETER                        :: routineN = 'kp_to_rs'

      INTEGER                                            :: handle, ikp

      IF (kpoints%nkp /= SIZE(array_kp, 3)) CALL cp_abort(__LOCATION__, &
                                                          "Provided kpoints and array_kp are inconsistent.")

      CALL timeset(routineN, handle)

      array_rs(:, :) = 0.0_dp
      DO ikp = 1, kpoints%nkp
         CALL add_kp_to_rs(array_kp(:, :, ikp), array_rs, cell, kpoints, ikp)
      END DO

      CALL timestop(handle)
   END SUBROUTINE kp_to_rs
END MODULE kpoint_k_r_trafo_simple
