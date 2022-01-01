!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to compute the Coulomb integral V_(alpha beta)(k) for a k-point k using lattice
!>        summation in real space. These integrals are e.g. needed in periodic RI for RPA, GW
!> \par History
!>       2018.05 created [Jan Wilhelm]
!> \author Jan Wilhelm
! **************************************************************************************************
MODULE kpoint_coulomb_2c
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE dbcsr_api,                       ONLY: &
        dbcsr_create, dbcsr_init_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_release_p, dbcsr_reserve_all_blocks, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE generic_shg_integrals,           ONLY: int_operators_r12_ab_shg
   USE generic_shg_integrals_init,      ONLY: contraction_matrix_shg
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE mathconstants,                   ONLY: twopi
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_coulomb_2c'

   PUBLIC :: build_2c_coulomb_matrix_kp

! **************************************************************************************************

   TYPE two_d_util_type
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)  :: block
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param matrix_v_kp ...
!> \param kpoints ...
!> \param basis_type ...
!> \param cell ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
!> \param size_lattice_sum ...
!> \param operator_type ...
! **************************************************************************************************
   SUBROUTINE build_2c_coulomb_matrix_kp(matrix_v_kp, kpoints, basis_type, cell, particle_set, qs_kind_set, &
                                         atomic_kind_set, size_lattice_sum, operator_type)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_v_kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      INTEGER                                            :: size_lattice_sum, operator_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_2c_coulomb_matrix_kp'

      INTEGER                                            :: handle, total_periodicity
      TYPE(dbcsr_type), POINTER                          :: matrix_v_L_tmp

      CALL timeset(routineN, handle)

      CALL check_periodicity(cell, kpoints, total_periodicity)

      CALL allocate_tmp(matrix_v_L_tmp, matrix_v_kp)

      CALL lattice_sum(matrix_v_kp, kpoints, basis_type, cell, particle_set, &
                       qs_kind_set, atomic_kind_set, size_lattice_sum, matrix_v_L_tmp, &
                       operator_type)

      CALL deallocate_tmp(matrix_v_L_tmp)

      CALL timestop(handle)

   END SUBROUTINE build_2c_coulomb_matrix_kp

! **************************************************************************************************
!> \brief ...
!> \param matrix_v_kp ...
!> \param kpoints ...
!> \param basis_type ...
!> \param cell ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
!> \param size_lattice_sum ...
!> \param matrix_v_L_tmp ...
!> \param operator_type ...
! **************************************************************************************************
   SUBROUTINE lattice_sum(matrix_v_kp, kpoints, basis_type, cell, particle_set, &
                          qs_kind_set, atomic_kind_set, size_lattice_sum, matrix_v_L_tmp, &
                          operator_type)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_v_kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      INTEGER                                            :: size_lattice_sum
      TYPE(dbcsr_type), POINTER                          :: matrix_v_L_tmp
      INTEGER                                            :: operator_type

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'lattice_sum'

      INTEGER :: factor, handle, i_block, i_x, i_x_inner, i_x_outer, ik, j_y, j_y_inner, &
         j_y_outer, k_z, k_z_inner, k_z_outer, nkp, x_max, x_min, y_max, y_min, z_max, z_min
      INTEGER, DIMENSION(3)                              :: nkp_grid, periodic
      REAL(KIND=dp)                                      :: coskl, sinkl
      REAL(KIND=dp), DIMENSION(3)                        :: vec_L, vec_s
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      TYPE(two_d_util_type), ALLOCATABLE, DIMENSION(:)   :: blocks_of_v_L, blocks_of_v_L_store
      TYPE(two_d_util_type), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: blocks_of_v_kp

      CALL timeset(routineN, handle)

      CALL get_kpoint_info(kpoints, nkp_grid=nkp_grid, nkp=nkp)
      CALL get_cell(cell=cell, h=hmat, periodic=periodic)

      IF (MODULO(nkp_grid(1), 2) == 1) THEN
         factor = 3**(size_lattice_sum - 1)
      ELSE IF (MODULO(nkp_grid(1), 2) == 0) THEN
         factor = 2**(size_lattice_sum - 1)
      END IF

      IF (MODULO(nkp_grid(1), 2) == 1) THEN
         x_min = -(factor*nkp_grid(1) - 1)/2
         x_max = (factor*nkp_grid(1) - 1)/2
      ELSE IF (MODULO(nkp_grid(1), 2) == 0) THEN
         x_min = -factor*nkp_grid(1)/2
         x_max = factor*nkp_grid(1)/2 - 1
      END IF
      IF (periodic(1) == 0) THEN
         x_min = 0
         x_max = 0
      END IF

      IF (MODULO(nkp_grid(2), 2) == 1) THEN
         y_min = -(factor*nkp_grid(2) - 1)/2
         y_max = (factor*nkp_grid(2) - 1)/2
      ELSE IF (MODULO(nkp_grid(2), 2) == 0) THEN
         y_min = -factor*nkp_grid(2)/2
         y_max = factor*nkp_grid(2)/2 - 1
      END IF
      IF (periodic(2) == 0) THEN
         y_min = 0
         y_max = 0
      END IF

      IF (MODULO(nkp_grid(3), 2) == 1) THEN
         z_min = -(factor*nkp_grid(3) - 1)/2
         z_max = (factor*nkp_grid(3) - 1)/2
      ELSE IF (MODULO(nkp_grid(3), 2) == 0) THEN
         z_min = -factor*nkp_grid(3)/2
         z_max = factor*nkp_grid(3)/2 - 1
      END IF
      IF (periodic(3) == 0) THEN
         z_min = 0
         z_max = 0
      END IF

      CALL allocate_blocks_of_v_kp(blocks_of_v_kp, matrix_v_kp)
      CALL allocate_blocks_of_v_L(blocks_of_v_L, matrix_v_L_tmp)
      CALL allocate_blocks_of_v_L(blocks_of_v_L_store, matrix_v_L_tmp)

      DO i_x_inner = 0, 2*nkp_grid(1) - 1
         DO j_y_inner = 0, 2*nkp_grid(2) - 1
            DO k_z_inner = 0, 2*nkp_grid(3) - 1

               DO i_x_outer = x_min, x_max + nkp_grid(1), 2*nkp_grid(1)
                  DO j_y_outer = y_min, y_max + nkp_grid(2), 2*nkp_grid(2)
                     DO k_z_outer = z_min, z_max + nkp_grid(3), 2*nkp_grid(3)

                        i_x = i_x_inner + i_x_outer
                        j_y = j_y_inner + j_y_outer
                        k_z = k_z_inner + k_z_outer

                        IF (i_x > x_max .OR. i_x < x_min .OR. &
                            j_y > y_max .OR. j_y < y_min .OR. &
                            k_z > z_max .OR. k_z < z_min) CYCLE

                        vec_s = [REAL(i_x, dp), REAL(j_y, dp), REAL(k_z, dp)]

                        vec_L = MATMUL(hmat, vec_s)

                        ! Compute (P 0 | Q vec_L) and store it in matrix_v_L_tmp
                        CALL compute_v_transl(matrix_v_L_tmp, blocks_of_v_L, vec_L, particle_set, &
                                              qs_kind_set, atomic_kind_set, basis_type, cell, &
                                              operator_type)

                        DO i_block = 1, SIZE(blocks_of_v_L)
                           blocks_of_v_L_store(i_block)%block(:, :) = blocks_of_v_L_store(i_block)%block(:, :) &
                                                                      + blocks_of_v_L(i_block)%block(:, :)
                        END DO

                     END DO
                  END DO
               END DO

               ! add exp(iq*vec_L) * (P 0 | Q vec_L) to V_PQ(q)
               DO ik = 1, nkp

                  ! coskl and sinkl are identical for all i_x_outer, j_y_outer, k_z_outer
                  coskl = COS(twopi*DOT_PRODUCT(vec_s(1:3), kpoints%xkp(1:3, ik)))
                  sinkl = SIN(twopi*DOT_PRODUCT(vec_s(1:3), kpoints%xkp(1:3, ik)))

                  DO i_block = 1, SIZE(blocks_of_v_L)

                     blocks_of_v_kp(ik, 1, i_block)%block(:, :) = blocks_of_v_kp(ik, 1, i_block)%block(:, :) &
                                                                  + coskl*blocks_of_v_L_store(i_block)%block(:, :)
                     blocks_of_v_kp(ik, 2, i_block)%block(:, :) = blocks_of_v_kp(ik, 2, i_block)%block(:, :) &
                                                                  + sinkl*blocks_of_v_L_store(i_block)%block(:, :)

                  END DO

               END DO

               DO i_block = 1, SIZE(blocks_of_v_L)

                  blocks_of_v_L_store(i_block)%block(:, :) = 0.0_dp

               END DO

            END DO
         END DO
      END DO

      CALL set_blocks_to_matrix_v_kp(matrix_v_kp, blocks_of_v_kp)

      CALL deallocate_blocks_of_v_kp(blocks_of_v_kp)
      CALL deallocate_blocks_of_v_L(blocks_of_v_L)
      CALL deallocate_blocks_of_v_L(blocks_of_v_L_store)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix_v_kp ...
!> \param blocks_of_v_kp ...
! **************************************************************************************************
   SUBROUTINE set_blocks_to_matrix_v_kp(matrix_v_kp, blocks_of_v_kp)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_v_kp
      TYPE(two_d_util_type), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: blocks_of_v_kp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'set_blocks_to_matrix_v_kp'

      INTEGER                                            :: col, handle, i_block, i_real_im, ik, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      DO ik = 1, SIZE(matrix_v_kp, 1)

         DO i_real_im = 1, 2

            i_block = 1

            CALL dbcsr_iterator_start(iter, matrix_v_kp(ik, i_real_im)%matrix)

            DO WHILE (dbcsr_iterator_blocks_left(iter))

               CALL dbcsr_iterator_next_block(iter, row, col, data_block)

               data_block(:, :) = blocks_of_v_kp(ik, i_real_im, i_block)%block(:, :)

               i_block = i_block + 1

            END DO

            CALL dbcsr_iterator_stop(iter)

         END DO

      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix_v_L_tmp ...
!> \param blocks_of_v_L ...
!> \param vec_L ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
!> \param basis_type ...
!> \param cell ...
!> \param operator_type ...
! **************************************************************************************************
   SUBROUTINE compute_v_transl(matrix_v_L_tmp, blocks_of_v_L, vec_L, particle_set, &
                               qs_kind_set, atomic_kind_set, basis_type, cell, operator_type)

      TYPE(dbcsr_type), POINTER                          :: matrix_v_L_tmp
      TYPE(two_d_util_type), ALLOCATABLE, DIMENSION(:)   :: blocks_of_v_L
      REAL(KIND=dp), DIMENSION(3)                        :: vec_L
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type
      TYPE(cell_type), POINTER                           :: cell
      INTEGER                                            :: operator_type

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_v_transl'

      INTEGER                                            :: col, handle, i_block, kind_a, kind_b, &
                                                            natom, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      REAL(dp), DIMENSION(3)                             :: ra, rab_L, rb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: contr_a, contr_b
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_a, basis_set_b

      CALL timeset(routineN, handle)

      NULLIFY (basis_set_a, basis_set_b, data_block)

      natom = SIZE(particle_set)

      ALLOCATE (kind_of(natom))

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      CALL dbcsr_set(matrix_v_L_tmp, 0.0_dp)

      i_block = 1

      CALL dbcsr_iterator_start(iter, matrix_v_L_tmp)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block)

         kind_a = kind_of(row)
         kind_b = kind_of(col)

         CALL get_qs_kind(qs_kind=qs_kind_set(kind_a), basis_set=basis_set_a, basis_type=basis_type)
         CALL get_qs_kind(qs_kind=qs_kind_set(kind_b), basis_set=basis_set_b, basis_type=basis_type)

         ra(1:3) = pbc(particle_set(row)%r(1:3), cell)
         rb(1:3) = pbc(particle_set(col)%r(1:3), cell)

         rab_L(1:3) = rb(1:3) - ra(1:3) + vec_L(1:3)

         CALL contraction_matrix_shg(basis_set_a, contr_a)
         CALL contraction_matrix_shg(basis_set_b, contr_b)

         blocks_of_v_L(i_block)%block = 0.0_dp

         CALL int_operators_r12_ab_shg(operator_type, blocks_of_v_L(i_block)%block, rab=rab_L, &
                                       fba=basis_set_a, fbb=basis_set_b, scona_shg=contr_a, sconb_shg=contr_b, &
                                       calculate_forces=.FALSE.)

         i_block = i_block + 1

         DEALLOCATE (contr_a, contr_b)

      END DO

      CALL dbcsr_iterator_stop(iter)

      DEALLOCATE (kind_of)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param blocks_of_v_kp ...
! **************************************************************************************************
   SUBROUTINE deallocate_blocks_of_v_kp(blocks_of_v_kp)
      TYPE(two_d_util_type), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: blocks_of_v_kp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_blocks_of_v_kp'

      INTEGER                                            :: handle, i_block, i_real_img, ik

      CALL timeset(routineN, handle)

      DO ik = 1, SIZE(blocks_of_v_kp, 1)
         DO i_real_img = 1, SIZE(blocks_of_v_kp, 2)
            DO i_block = 1, SIZE(blocks_of_v_kp, 3)
               DEALLOCATE (blocks_of_v_kp(ik, i_real_img, i_block)%block)
            END DO
         END DO
      END DO

      DEALLOCATE (blocks_of_v_kp)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param blocks_of_v_L ...
! **************************************************************************************************
   SUBROUTINE deallocate_blocks_of_v_L(blocks_of_v_L)
      TYPE(two_d_util_type), ALLOCATABLE, DIMENSION(:)   :: blocks_of_v_L

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_blocks_of_v_L'

      INTEGER                                            :: handle, i_block

      CALL timeset(routineN, handle)

      DO i_block = 1, SIZE(blocks_of_v_L, 1)
         DEALLOCATE (blocks_of_v_L(i_block)%block)
      END DO

      DEALLOCATE (blocks_of_v_L)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param blocks_of_v_L ...
!> \param matrix_v_L_tmp ...
! **************************************************************************************************
   SUBROUTINE allocate_blocks_of_v_L(blocks_of_v_L, matrix_v_L_tmp)
      TYPE(two_d_util_type), ALLOCATABLE, DIMENSION(:)   :: blocks_of_v_L
      TYPE(dbcsr_type), POINTER                          :: matrix_v_L_tmp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_blocks_of_v_L'

      INTEGER                                            :: col, handle, i_block, nblocks, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      nblocks = 0

      CALL dbcsr_iterator_start(iter, matrix_v_L_tmp)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block)

         nblocks = nblocks + 1

      END DO

      CALL dbcsr_iterator_stop(iter)

      ALLOCATE (blocks_of_v_L(nblocks))

      i_block = 1

      CALL dbcsr_iterator_start(iter, matrix_v_L_tmp)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block)

         ALLOCATE (blocks_of_v_L(i_block)%block(SIZE(data_block, 1), SIZE(data_block, 2)))
         blocks_of_v_L(i_block)%block = 0.0_dp

         i_block = i_block + 1

      END DO

      CALL dbcsr_iterator_stop(iter)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param blocks_of_v_kp ...
!> \param matrix_v_kp ...
! **************************************************************************************************
   SUBROUTINE allocate_blocks_of_v_kp(blocks_of_v_kp, matrix_v_kp)
      TYPE(two_d_util_type), ALLOCATABLE, &
         DIMENSION(:, :, :)                              :: blocks_of_v_kp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_v_kp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_blocks_of_v_kp'

      INTEGER                                            :: col, handle, i_block, i_real_img, ik, &
                                                            nblocks, row
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL timeset(routineN, handle)

      nblocks = 0

      CALL dbcsr_iterator_start(iter, matrix_v_kp(1, 1)%matrix)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block)

         nblocks = nblocks + 1

      END DO

      CALL dbcsr_iterator_stop(iter)

      ALLOCATE (blocks_of_v_kp(SIZE(matrix_v_kp, 1), SIZE(matrix_v_kp, 2), nblocks))

      DO ik = 1, SIZE(matrix_v_kp, 1)

         DO i_real_img = 1, SIZE(matrix_v_kp, 2)

            i_block = 1

            CALL dbcsr_iterator_start(iter, matrix_v_kp(ik, i_real_img)%matrix)

            DO WHILE (dbcsr_iterator_blocks_left(iter))

               CALL dbcsr_iterator_next_block(iter, row, col, data_block)

               ALLOCATE (blocks_of_v_kp(ik, i_real_img, i_block)%block(SIZE(data_block, 1), SIZE(data_block, 2)))
               blocks_of_v_kp(ik, i_real_img, i_block)%block = 0.0_dp

               i_block = i_block + 1

            END DO

            CALL dbcsr_iterator_stop(iter)

         END DO

      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param cell ...
!> \param kpoints ...
!> \param total_periodicity ...
! **************************************************************************************************
   SUBROUTINE check_periodicity(cell, kpoints, total_periodicity)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: total_periodicity

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'check_periodicity'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(3)                              :: nkp_grid, periodic

      CALL timeset(routineN, handle)

      CALL get_cell(cell=cell, periodic=periodic)
      CALL get_kpoint_info(kpoints, nkp_grid=nkp_grid)

      IF (periodic(1) == 0) THEN
         CPASSERT(nkp_grid(1) == 1)
      END IF
      IF (periodic(2) == 0) THEN
         CPASSERT(nkp_grid(2) == 1)
      END IF
      IF (periodic(3) == 0) THEN
         CPASSERT(nkp_grid(3) == 1)
      END IF

      total_periodicity = SUM(periodic)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix_v_L_tmp ...
!> \param matrix_v_kp ...
! **************************************************************************************************
   SUBROUTINE allocate_tmp(matrix_v_L_tmp, matrix_v_kp)

      TYPE(dbcsr_type), POINTER                          :: matrix_v_L_tmp
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_v_kp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'allocate_tmp'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      NULLIFY (matrix_v_L_tmp)
      CALL dbcsr_init_p(matrix_v_L_tmp)
      CALL dbcsr_create(matrix=matrix_v_L_tmp, &
                        template=matrix_v_kp(1, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_reserve_all_blocks(matrix_v_L_tmp)
      CALL dbcsr_set(matrix_v_L_tmp, 0.0_dp)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix_v_L_tmp ...
! **************************************************************************************************
   SUBROUTINE deallocate_tmp(matrix_v_L_tmp)

      TYPE(dbcsr_type), POINTER                          :: matrix_v_L_tmp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'deallocate_tmp'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL dbcsr_release_p(matrix_v_L_tmp)

      CALL timestop(handle)

   END SUBROUTINE

END MODULE kpoint_coulomb_2c
