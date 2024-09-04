!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Contains utility routines for the active space module
!> \par History
!>      04.2023 created [SB]
!> \author SB
! **************************************************************************************************
MODULE qs_active_space_utils

   USE cp_dbcsr_api,                    ONLY: dbcsr_csr_type
   USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                              cp_fm_get_info,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE qs_active_space_types,           ONLY: csr_idx_from_combined,&
                                              csr_idx_to_combined,&
                                              eri_type,&
                                              get_irange_csr
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_active_space_utils'

   PUBLIC :: subspace_matrix_to_array, eri_to_array

CONTAINS

! **************************************************************************************************
!> \brief Copy a (square portion) of a `cp_fm_type` matrix to a standard 1D Fortran array
!> \param source_matrix the matrix from where the data is taken
!> \param target_array the array were the data is copied to
!> \param row_index a list containing the row subspace indices
!> \param col_index a list containing the column subspace indices
! **************************************************************************************************
   SUBROUTINE subspace_matrix_to_array(source_matrix, target_array, row_index, col_index)
      TYPE(cp_fm_type), INTENT(IN)                       :: source_matrix
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: target_array
      INTEGER, DIMENSION(:), INTENT(IN)                  :: row_index, col_index

      INTEGER                                            :: i, i_sub, j, j_sub, max_col, max_row, &
                                                            ncols, nrows
      REAL(KIND=dp)                                      :: mval

      CALL cp_fm_get_info(source_matrix, nrow_global=max_row, ncol_global=max_col)
      nrows = SIZE(row_index)
      ncols = SIZE(col_index)

      CPASSERT(MAXVAL(row_index) <= max_row)
      CPASSERT(MAXVAL(col_index) <= max_col)
      CPASSERT(MINVAL(row_index) > 0)
      CPASSERT(MINVAL(col_index) > 0)
      CPASSERT(nrows <= max_row)
      CPASSERT(ncols <= max_col)

      CPASSERT(SIZE(target_array) == nrows*ncols)

      DO j = 1, ncols
         j_sub = col_index(j)
         DO i = 1, nrows
            i_sub = row_index(i)
            CALL cp_fm_get_element(source_matrix, i_sub, j_sub, mval)
            target_array(i + (j - 1)*nrows) = mval
         END DO
      END DO
   END SUBROUTINE subspace_matrix_to_array

! **************************************************************************************************
!> \brief Copy the eri tensor for spins isp1 and isp2 to a standard 1D Fortran array
!> \param eri_env the eri environment
!> \param array the 1D Fortran array where the eri are copied to
!> \param active_orbitals a list containing the active orbitals indices
!> \param spin1 the spin of the bra
!> \param spin2 the spin of the ket
! **************************************************************************************************
   SUBROUTINE eri_to_array(eri_env, array, active_orbitals, spin1, spin2)
      TYPE(eri_type), INTENT(IN)                         :: eri_env
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: array
      INTEGER, DIMENSION(:, :), INTENT(IN)               :: active_orbitals
      INTEGER, INTENT(IN)                                :: spin1, spin2

      INTEGER                                            :: i, i1, i12, i12l, i2, i3, i34, i34l, i4, &
                                                            ijkl, ijlk, irptr, j, jikl, jilk, k, &
                                                            klij, klji, l, lkij, lkji, nindex, &
                                                            nmo_active, nmo_max
      INTEGER, DIMENSION(2)                              :: irange
      REAL(KIND=dp)                                      :: erival
      TYPE(dbcsr_csr_type), POINTER                      :: eri
      TYPE(mp_comm_type)                                 :: mp_group

      nmo_active = SIZE(active_orbitals, 1)
      nmo_max = eri_env%norb
      nindex = (nmo_max*(nmo_max + 1))/2
      IF (spin1 == 1 .AND. spin2 == 1) THEN
         eri => eri_env%eri(1)%csr_mat
      ELSE IF ((spin1 == 1 .AND. spin2 == 2) .OR. (spin1 == 2 .AND. spin2 == 1)) THEN
         eri => eri_env%eri(2)%csr_mat
      ELSE
         eri => eri_env%eri(3)%csr_mat
      END IF

      CALL mp_group%set_handle(eri%mp_group%get_handle())
      irange = get_irange_csr(nindex, mp_group)

      array = 0.0_dp

      DO i = 1, nmo_active
         i1 = active_orbitals(i, spin1)
         DO j = i, nmo_active
            i2 = active_orbitals(j, spin1)
            i12 = csr_idx_to_combined(i1, i2, nmo_max)
            IF (i12 >= irange(1) .AND. i12 <= irange(2)) THEN
               i12l = i12 - irange(1) + 1
               irptr = eri%rowptr_local(i12l) - 1
               DO i34l = 1, eri%nzerow_local(i12l)
                  i34 = eri%colind_local(irptr + i34l)
                  CALL csr_idx_from_combined(i34, nmo_max, i3, i4)
! The FINDLOC intrinsic function of the Fortran 2008 standard is only available since GCC 9
! That is why we use a custom-made implementation of this function for this compiler
#if __GNUC__ < 9
                  k = cp_findloc(active_orbitals(:, spin2), i3)
                  l = cp_findloc(active_orbitals(:, spin2), i4)
#else
                  k = FINDLOC(active_orbitals(:, spin2), i3, dim=1)
                  l = FINDLOC(active_orbitals(:, spin2), i4, dim=1)
#endif
                  erival = eri%nzval_local%r_dp(irptr + i34l)

                  ! 8-fold permutational symmetry
                  ijkl = i + (j - 1)*nmo_active + (k - 1)*nmo_active**2 + (l - 1)*nmo_active**3
                  jikl = j + (i - 1)*nmo_active + (k - 1)*nmo_active**2 + (l - 1)*nmo_active**3
                  ijlk = i + (j - 1)*nmo_active + (l - 1)*nmo_active**2 + (k - 1)*nmo_active**3
                  jilk = j + (i - 1)*nmo_active + (l - 1)*nmo_active**2 + (k - 1)*nmo_active**3
                  array(ijkl) = erival
                  array(jikl) = erival
                  array(ijlk) = erival
                  array(jilk) = erival
                  IF (spin1 == spin2) THEN
                     klij = k + (l - 1)*nmo_active + (i - 1)*nmo_active**2 + (j - 1)*nmo_active**3
                     lkij = l + (k - 1)*nmo_active + (i - 1)*nmo_active**2 + (j - 1)*nmo_active**3
                     klji = k + (l - 1)*nmo_active + (j - 1)*nmo_active**2 + (i - 1)*nmo_active**3
                     lkji = l + (k - 1)*nmo_active + (j - 1)*nmo_active**2 + (i - 1)*nmo_active**3
                     array(klij) = erival
                     array(lkij) = erival
                     array(klji) = erival
                     array(lkji) = erival
                  END IF
               END DO
            END IF
         END DO
      END DO
      CALL mp_group%sum(array)

   END SUBROUTINE eri_to_array

#if __GNUC__ < 9
! **************************************************************************************************
!> \brief This function implements the FINDLOC function of the Fortran 2008 standard for the case needed above
!>        To be removed as soon GCC 8 is dropped.
!> \param array ...
!> \param value ...
!> \return ...
! **************************************************************************************************
   PURE INTEGER FUNCTION cp_findloc(array, value) RESULT(loc)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: array
      INTEGER, INTENT(IN)                                :: value

      INTEGER                                            :: idx

      loc = 0

      DO idx = 1, SIZE(array)
      IF (array(idx) == value) THEN
         loc = idx
         RETURN
      END IF
      END DO

   END FUNCTION cp_findloc
#endif

END MODULE qs_active_space_utils
