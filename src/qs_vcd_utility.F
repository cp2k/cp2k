!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Helper functions for the vcd code
!> \author E. Ditler 2020
! *****************************************************************************
MODULE qs_vcd_utility
! the ones above are needed for sure, below who knows
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE pw_methods,                      ONLY: pw_transfer
   USE pw_types,                        ONLY: pw_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
#include "base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_vcd_utility'

   PUBLIC :: build_d_matrix
   PUBLIC :: write_cpfm, write_cpfm3, write_cpfm9
   PUBLIC :: write_dbcsr, write_dbcsr3, write_dbcsr9
   PUBLIC :: write_pw_grid, write_pw_g_grid, write_pw_grid_coordinates
   PUBLIC :: trace_occupied
   PUBLIC :: read_matrix, read_and_add_matrix, read_pw_grid, read_and_add_pw_grid
   PUBLIC :: fm_clear_lower
   PUBLIC :: cp_fm_checksum

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \return ...
! **************************************************************************************************
   FUNCTION cp_fm_checksum(matrix) RESULT(checksum)
      TYPE(cp_fm_type), POINTER                          :: matrix
      REAL(KIND=dp)                                      :: checksum

      INTEGER                                            :: i
      REAL(KIND=dp)                                      :: row_sum

      DO i = 1, UBOUND(matrix%local_data, 1)
         row_sum = REAL(DOT_PRODUCT(matrix%local_data(i, :), matrix%local_data(i, :)), dp)
         checksum = checksum + row_sum
      END DO

   END FUNCTION cp_fm_checksum

! **************************************************************************************************
!> \brief ...
!> \param filename ...
!> \param nao ...
!> \param matrix ...
!> \param nmo ...
! **************************************************************************************************
   SUBROUTINE read_matrix(filename, nao, matrix, nmo)
      CHARACTER(len=*)                                   :: filename
      INTEGER                                            :: nao
      TYPE(cp_fm_type), POINTER                          :: matrix
      INTEGER, OPTIONAL                                  :: nmo

      INTEGER                                            :: i, j, l, n_elements
      REAL(dp)                                           :: temp

      n_elements = nao*nao
      IF (PRESENT(nmo)) n_elements = nao*nmo

      matrix%local_data(:, :) = HUGE(1._dp)
      PRINT *, 'reading ', TRIM(filename)
      PRINT *, 'reading ', n_elements, nao
      OPEN (unit=12312, file=TRIM(filename), status='old')
      DO l = 1, n_elements
         READ (12312, *) i, j, temp
         matrix%local_data(i, j) = temp
      END DO
      CLOSE (12312)

   END SUBROUTINE read_matrix

! **************************************************************************************************
!> \brief ...
!> \param filename ...
!> \param nao ...
!> \param matrix ...
!> \param nmo ...
! **************************************************************************************************
   SUBROUTINE read_and_add_matrix(filename, nao, matrix, nmo)
      CHARACTER(len=*)                                   :: filename
      INTEGER                                            :: nao
      TYPE(cp_fm_type), POINTER                          :: matrix
      INTEGER, OPTIONAL                                  :: nmo

      INTEGER                                            :: i, j, l, n_elements
      REAL(dp)                                           :: temp

      n_elements = nao*nao
      IF (PRESENT(nmo)) n_elements = nao*nmo

      PRINT *, 'reading ', TRIM(filename)
      PRINT *, 'reading ', n_elements, nao
      OPEN (unit=12312, file=TRIM(filename), status='old')
      DO l = 1, n_elements
         READ (12312, *) i, j, temp
         matrix%local_data(i, j) = matrix%local_data(i, j) + temp
      END DO
      CLOSE (12312)

   END SUBROUTINE read_and_add_matrix

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
! **************************************************************************************************
   SUBROUTINE fm_clear_lower(matrix)
      TYPE(cp_fm_type), POINTER                          :: matrix

      INTEGER                                            :: i, j

      DO i = 1, UBOUND(matrix%local_data, 1)
      DO j = 1, UBOUND(matrix%local_data, 2)
         IF (j < i) THEN
            matrix%local_data(i, j) = 0._dp
         END IF
      END DO
      END DO

   END SUBROUTINE fm_clear_lower

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param occupied ...
!> \param trace ...
! **************************************************************************************************
   SUBROUTINE trace_occupied(matrix, occupied, trace)
      TYPE(cp_fm_type), INTENT(in)                       :: matrix
      INTEGER, INTENT(in)                                :: occupied
      REAL(dp), INTENT(out)                              :: trace

      INTEGER                                            :: i

      trace = 0._dp
      DO i = 1, occupied
         trace = trace + matrix%local_data(i, i)
      END DO

   END SUBROUTINE trace_occupied

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param name ...
! **************************************************************************************************
   SUBROUTINE write_pw_grid_coordinates(pw, name)
      TYPE(pw_type), POINTER                             :: pw
      CHARACTER(len=*)                                   :: name

      INTEGER                                            :: i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(dp)                                           :: x, y, z
      REAL(dp), DIMENSION(3)                             :: roffset

      lb(:) = pw%pw_grid%bounds_local(1, 1:3)
      ub(:) = pw%pw_grid%bounds_local(2, 1:3)
      roffset = REAL(MODULO(pw%pw_grid%npts, 2), dp)*pw%pw_grid%dr/2._dp

      ! We write to a file:
      OPEN (12345, file=name)

      DO i = lb(1), ub(1)
         DO j = lb(2), ub(2)
            DO k = lb(3), ub(3)
               x = pw%pw_grid%dr(1)*(i - lb(1)) + roffset(1)
               y = pw%pw_grid%dr(2)*(j - lb(2)) + roffset(2)
               z = pw%pw_grid%dr(3)*(k - lb(3)) + roffset(3)

               WRITE (12345, *) x, y, z
            END DO
         END DO
      END DO

      CLOSE (12345)
   END SUBROUTINE write_pw_grid_coordinates

! **************************************************************************************************
!> \brief ...
!> \param filename ...
!> \param pw ...
! **************************************************************************************************
   SUBROUTINE read_pw_grid(filename, pw)
      CHARACTER(len=*)                                   :: filename
      TYPE(pw_type), POINTER                             :: pw

      INTEGER                                            :: i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub

      lb(:) = pw%pw_grid%bounds_local(1, 1:3)
      ub(:) = pw%pw_grid%bounds_local(2, 1:3)

      PRINT *, 'reading ', TRIM(filename)
      OPEN (unit=12312, file=TRIM(filename), status='old')
      DO i = lb(1), ub(1)
         DO j = lb(2), ub(2)
            DO k = lb(3), ub(3)
               READ (12312, *) pw%cr3d(i, j, k)
            END DO
         END DO
      END DO
      CLOSE (12312)
   END SUBROUTINE read_pw_grid

! **************************************************************************************************
!> \brief ...
!> \param filename ...
!> \param pw ...
! **************************************************************************************************
   SUBROUTINE read_and_add_pw_grid(filename, pw)
      CHARACTER(len=*)                                   :: filename
      TYPE(pw_type), POINTER                             :: pw

      INTEGER                                            :: i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(dp)                                           :: temp

      lb(:) = pw%pw_grid%bounds_local(1, 1:3)
      ub(:) = pw%pw_grid%bounds_local(2, 1:3)

      PRINT *, 'reading ', TRIM(filename)
      OPEN (unit=12312, file=TRIM(filename), status='old')
      DO i = lb(1), ub(1)
         DO j = lb(2), ub(2)
            DO k = lb(3), ub(3)
               READ (12312, *) temp
               pw%cr3d(i, j, k) = pw%cr3d(i, j, k) + temp
            END DO
         END DO
      END DO
      CLOSE (12312)
   END SUBROUTINE read_and_add_pw_grid

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param name ...
!> \param scale ...
! **************************************************************************************************
   SUBROUTINE write_pw_grid(pw, name, scale)
      TYPE(pw_type), POINTER                             :: pw
      CHARACTER(len=*)                                   :: name
      REAL(dp), OPTIONAL                                 :: scale

      INTEGER                                            :: i, j, k
      INTEGER, DIMENSION(3)                              :: lb, ub
      REAL(dp)                                           :: my_scale
      REAL(dp), DIMENSION(3)                             :: roffset

      my_scale = 1._dp
      IF (PRESENT(scale)) my_scale = scale

      lb(:) = pw%pw_grid%bounds_local(1, 1:3)
      ub(:) = pw%pw_grid%bounds_local(2, 1:3)
      roffset = -REAL(MODULO(pw%pw_grid%npts, 2), dp)*pw%pw_grid%dr/2.0_dp
      ! We write to a file:
      OPEN (12345, file=TRIM(name))

      ! Now we can loop over the 3 dimensions + spin:
      DO i = lb(1), ub(1)
         DO j = lb(2), ub(2)
            DO k = lb(3), ub(3)
               WRITE (12345, *) pw%cr3d(i, j, k)*my_scale
            END DO
         END DO
      END DO

      CLOSE (12345)
   END SUBROUTINE write_pw_grid

! **************************************************************************************************
!> \brief ...
!> \param pw ...
!> \param pw_transfer_to ...
!> \param name ...
! **************************************************************************************************
   SUBROUTINE write_pw_g_grid(pw, pw_transfer_to, name)
      TYPE(pw_type), POINTER                             :: pw, pw_transfer_to
      CHARACTER(len=*)                                   :: name

      INTEGER                                            :: i, j, k

      CALL pw_transfer(pw, pw_transfer_to)
      ! We write to a file:
      OPEN (12345, file=name)

      ! Now we can loop over the 3 dimensions + spin:
      DO i = LBOUND(pw_transfer_to%cc3d, 1), UBOUND(pw_transfer_to%cc3d, 1)
      DO j = LBOUND(pw_transfer_to%cc3d, 2), UBOUND(pw_transfer_to%cc3d, 2)
      DO k = LBOUND(pw_transfer_to%cc3d, 3), UBOUND(pw_transfer_to%cc3d, 3)
         ! WRITE(12345, *) i, j, k, real(real(pw_transfer_to%cc3d(i, j, k))), aimag(pw_transfer_to%cc3d(i, j, k))
      END DO
      END DO
      END DO

      CLOSE (12345)

   END SUBROUTINE write_pw_g_grid

! **************************************************************************************************
!> \brief Generates the transformation matrix from hessian in cartesian into
!>      internal coordinates (based on Gram-Schmidt orthogonalization)
!> \param mat ...
!> \param dof ...
!> \param dout ...
!> \param natoms ...
!> \author Teodoro Laino 08.2006
! **************************************************************************************************
   SUBROUTINE build_d_matrix(mat, dof, dout, natoms)
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: mat
      INTEGER, INTENT(in)                                :: dof
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: dout
      INTEGER, INTENT(in)                                :: natoms

      CHARACTER(len=*), PARAMETER :: routineN = 'build_d_matrix', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ifound, iseq, j, nvib
      REAL(kind=dp)                                      :: norm
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: work
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: d

      CALL timeset(routinen, handle)
      ! generate the missing vectors of the orthogonal basis set
      nvib = 3*natoms - dof
      ALLOCATE (work(3*natoms))
      ALLOCATE (d(3*natoms, 3*natoms))
      ALLOCATE (dout(3*natoms, nvib))
      ! check first orthogonality in the first element of the basis set
      DO i = 1, dof
         d(:, i) = mat(:, i)
         DO j = i + 1, dof
            norm = DOT_PRODUCT(mat(:, i), mat(:, j))
            CPASSERT(ABS(norm) < 1.0E4_dp*EPSILON(0.0_dp))
         END DO
      END DO
      ! generate the nvib orthogonal vectors
      iseq = 0
      ifound = 0
      DO WHILE (ifound /= nvib)
         iseq = iseq + 1
         CPASSERT(iseq <= 3*natoms)
         work = 0.0_dp
         work(iseq) = 1.0_dp
         ! gram schmidt orthogonalization
         DO i = 1, dof + ifound
            norm = DOT_PRODUCT(work, d(:, i))
            work(:) = work - norm*d(:, i)
         END DO
         ! check norm of the new generated vector
         norm = SQRT(DOT_PRODUCT(work, work))
         IF (norm >= 1e9_dp*EPSILON(0.0_dp)) THEN
            ! accept new vector
            ifound = ifound + 1
            d(:, dof + ifound) = work/norm
         END IF
      END DO
      CPASSERT(dof + ifound == 3*natoms)
      dout = d(:, dof + 1:)
      DEALLOCATE (work)
      DEALLOCATE (d)
      DEALLOCATE (mat)
      CALL timestop(handle)
   END SUBROUTINE build_d_matrix

! Routines to write matrices out
! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param name ...
!> \param iatom ...
!> \param idir ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE write_cpfm(matrix, name, iatom, idir, unit_nr)
      TYPE(cp_fm_type), POINTER                          :: matrix
      CHARACTER(len=*)                                   :: name
      INTEGER                                            :: iatom, idir
      INTEGER, OPTIONAL                                  :: unit_nr

      INTEGER                                            :: ii, jj

      INTEGER:: myunit_nr

      myunit_nr = 11111
      IF (PRESENT(unit_nr)) THEN
         myunit_nr = unit_nr
      END IF
      WRITE (myunit_nr, *) 'start '//name, iatom, idir
      DO ii = LBOUND(matrix%local_data, 1), UBOUND(matrix%local_data, 1)
      DO jj = LBOUND(matrix%local_data, 2), UBOUND(matrix%local_data, 2)
         WRITE (myunit_nr, "(I4,I4,F25.10)") ii, jj, matrix%local_data(ii, jj)
      END DO
      END DO
      WRITE (myunit_nr, *) 'end '//name
   END SUBROUTINE write_cpfm

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param name ...
!> \param iatom ...
!> \param idir ...
! **************************************************************************************************
   SUBROUTINE write_cpfm3(matrix, name, iatom, idir)
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: matrix
      CHARACTER(len=*)                                   :: name
      INTEGER                                            :: iatom, idir

      INTEGER                                            :: ii, jj, kk

      WRITE (11111, *) 'start '//name, iatom, idir
      DO kk = 1, 3
         WRITE (11111, *) 'dir kk', kk
         DO ii = LBOUND(matrix(kk)%matrix%local_data, 1), UBOUND(matrix(kk)%matrix%local_data, 1)
         DO jj = LBOUND(matrix(kk)%matrix%local_data, 2), UBOUND(matrix(kk)%matrix%local_data, 2)
            WRITE (11111, "(I4,I4,F25.10)") ii, jj, matrix(kk)%matrix%local_data(ii, jj)
         END DO
         END DO
      END DO
      WRITE (11111, *) 'end '//name
   END SUBROUTINE write_cpfm3

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param name ...
!> \param iatom ...
!> \param idir ...
! **************************************************************************************************
   SUBROUTINE write_cpfm9(matrix, name, iatom, idir)
      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: matrix
      CHARACTER(len=*)                                   :: name
      INTEGER                                            :: iatom, idir

      INTEGER                                            :: ii, jj, kk, ll

      WRITE (11111, *) 'start '//name, iatom, idir
      DO kk = 1, 3
      DO ll = 1, 3
         WRITE (11111, *) 'dir kk ll', kk, ll
         DO ii = LBOUND(matrix(kk, ll)%matrix%local_data, 1), UBOUND(matrix(kk, ll)%matrix%local_data, 1)
         DO jj = LBOUND(matrix(kk, ll)%matrix%local_data, 2), UBOUND(matrix(kk, ll)%matrix%local_data, 2)
            WRITE (11111, "(I4,I4,F25.10)") ii, jj, matrix(kk, ll)%matrix%local_data(ii, jj)
         END DO
         END DO
      END DO
      END DO
      WRITE (11111, *) 'end '//name
   END SUBROUTINE write_cpfm9

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param name ...
!> \param iatom ...
!> \param idir ...
!> \param qs_env ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE write_dbcsr(matrix, name, iatom, idir, qs_env, unit_nr)
      TYPE(dbcsr_type), POINTER                          :: matrix
      CHARACTER(len=*)                                   :: name
      INTEGER                                            :: iatom, idir
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, OPTIONAL                                  :: unit_nr

      INTEGER                                            :: ii, jj, myunit_nr, nao, nmo
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fmatrix, mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      myunit_nr = 11111
      IF (PRESENT(unit_nr)) myunit_nr = unit_nr

      CALL get_qs_env(qs_env=qs_env, para_env=para_env, mos=mos)
      CALL get_mo_set(mo_set=mos(1), mo_coeff=mo_coeff, nao=nao, nmo=nmo)

      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, &
                               para_env=para_env, context=mo_coeff%matrix_struct%context)
      ALLOCATE (fmatrix)
      CALL cp_fm_create(fmatrix, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      CALL copy_dbcsr_to_fm(matrix, fmatrix)
      WRITE (myunit_nr, *) 'start '//name, iatom, idir
      DO ii = LBOUND(fmatrix%local_data, 1), UBOUND(fmatrix%local_data, 1)
      DO jj = LBOUND(fmatrix%local_data, 2), UBOUND(fmatrix%local_data, 2)
         WRITE (myunit_nr, "(I4,I4,F25.10)") ii, jj, fmatrix%local_data(ii, jj)
      END DO
      END DO
      WRITE (myunit_nr, *) 'end '//name

      CALL cp_fm_release(fmatrix)
      DEALLOCATE (fmatrix)
   END SUBROUTINE write_dbcsr

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param name ...
!> \param iatom ...
!> \param idir ...
!> \param qs_env ...
!> \param scale ...
! **************************************************************************************************
   SUBROUTINE write_dbcsr3(matrix, name, iatom, idir, qs_env, scale)
      TYPE(dbcsr_p_type), DIMENSION(:)                   :: matrix
      CHARACTER(len=*)                                   :: name
      INTEGER                                            :: iatom, idir
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), OPTIONAL                                 :: scale

      INTEGER                                            :: ii, jj, kk, nao, nmo
      REAL(dp)                                           :: factor
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fmatrix, mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL get_qs_env(qs_env=qs_env, para_env=para_env, mos=mos)
      CALL get_mo_set(mo_set=mos(1), mo_coeff=mo_coeff, nao=nao, nmo=nmo)

      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, &
                               para_env=para_env, context=mo_coeff%matrix_struct%context)
      CALL cp_fm_create(fmatrix, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      factor = 1._dp
      IF (PRESENT(scale)) factor = scale

      WRITE (11111, *) 'start '//name, iatom, idir
      DO kk = 1, 3
         WRITE (11111, *) 'dir kk', kk
         IF (ASSOCIATED(matrix(kk)%matrix)) THEN
            CALL copy_dbcsr_to_fm(matrix(kk)%matrix, fmatrix)
            DO ii = LBOUND(fmatrix%local_data, 1), UBOUND(fmatrix%local_data, 1)
            DO jj = LBOUND(fmatrix%local_data, 2), UBOUND(fmatrix%local_data, 2)
               WRITE (11111, "(I4,I4,F25.10)") ii, jj, fmatrix%local_data(ii, jj)*factor
               ! write (11111, "(I4,I4,F15.8)") ii, jj, fmatrix%local_data(ii,jj)
            END DO
            END DO
         END IF
      END DO
      WRITE (11111, *) 'end '//name
      CALL cp_fm_release(fmatrix)
   END SUBROUTINE write_dbcsr3

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param name ...
!> \param iatom ...
!> \param idir ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE write_dbcsr9(matrix, name, iatom, idir, qs_env)
      TYPE(dbcsr_p_type), DIMENSION(:, :)                :: matrix
      CHARACTER(len=*)                                   :: name
      INTEGER                                            :: iatom, idir
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: ii, jj, kk, ll, nao, nmo
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: fmatrix, mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL get_qs_env(qs_env=qs_env, para_env=para_env, mos=mos)
      CALL get_mo_set(mo_set=mos(1), mo_coeff=mo_coeff, nao=nao, nmo=nmo)

      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, &
                               para_env=para_env, context=mo_coeff%matrix_struct%context)
      CALL cp_fm_create(fmatrix, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      WRITE (11111, *) 'start '//name, iatom, idir
      DO kk = 1, 3
      DO ll = 1, 3
         WRITE (11111, *) 'dir kk ll', kk, ll
         CALL copy_dbcsr_to_fm(matrix(kk, ll)%matrix, fmatrix)
         DO ii = LBOUND(fmatrix%local_data, 1), UBOUND(fmatrix%local_data, 1)
         DO jj = LBOUND(fmatrix%local_data, 2), UBOUND(fmatrix%local_data, 2)
            WRITE (11111, "(I4,I4,F15.8)") ii, jj, fmatrix%local_data(ii, jj)
         END DO !jj
         END DO !ii
      END DO !ll
      END DO !kk
      WRITE (11111, *) 'end '//name
      CALL cp_fm_release(fmatrix)
   END SUBROUTINE write_dbcsr9

END MODULE qs_vcd_utility
