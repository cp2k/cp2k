!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE cp_fm_dlaf_api

   USE cp_fm_basic_linalg, ONLY: cp_fm_uplo_to_full
   USE cp_fm_types, ONLY: cp_fm_type
   USE kinds, ONLY: dp
#include "../base/base_uses.f90"

#if defined(__DLAF)
   USE dlaf_fortran, ONLY: dlaf_pdpotrf, &
                           dlaf_pdsyevd, &
                           dlaf_pdsygvd, &
                           dlaf_pspotrf
#endif

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_fm_dlaf_api'

   PUBLIC :: cp_pdpotrf_dlaf, cp_pspotrf_dlaf
   PUBLIC :: cp_fm_diag_dlaf, cp_fm_diag_gen_dlaf

CONTAINS

!***************************************************************************************************
!> \brief Cholesky factorization using DLA-Future
!> \param uplo ...
!> \param n Matrix size
!> \param a Local matrix
!> \param ia Row index of first row (has to be 1)
!> \param ja Col index of first column ()
!> \param desca ScaLAPACK matrix descriptor
!> \param info 0 if factorization completed normally
!> \author Rocco Meli
!> \author Mikael Simberg
!> \author Mathieu Taillefumier
! **************************************************************************************************
   SUBROUTINE cp_pdpotrf_dlaf(uplo, n, a, ia, ja, desca, info)
      CHARACTER, INTENT(IN)                              :: uplo
      INTEGER, INTENT(IN)                                :: n
      REAL(KIND=dp), DIMENSION(:, :), TARGET             :: a
      INTEGER, INTENT(IN)                                :: ia, ja
      INTEGER, DIMENSION(9)                              :: desca
      INTEGER, TARGET                                    :: info

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_pdpotrf_dlaf'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
#if defined(__DLAF)
      CALL dlaf_pdpotrf(uplo, n, a, ia, ja, desca, info)
#else
      MARK_USED(uplo)
      MARK_USED(n)
      MARK_USED(a)
      MARK_USED(ia)
      MARK_USED(ja)
      MARK_USED(desca)
      MARK_USED(info)
      CPABORT("CP2K compiled without the DLA-Future library.")
#endif
      CALL timestop(handle)
   END SUBROUTINE cp_pdpotrf_dlaf

!***************************************************************************************************
!> \brief Cholesky factorization using DLA-Future
!> \param uplo ...
!> \param n Matrix size
!> \param a Local matrix
!> \param ia Row index of first row (has to be 1)
!> \param ja Col index of first column ()
!> \param desca ScaLAPACK matrix descriptor
!> \param info 0 if factorization completed normally
!> \author Rocco Meli
!> \author Mikael Simberg
!> \author Mathieu Taillefumier
! **************************************************************************************************
   SUBROUTINE cp_pspotrf_dlaf(uplo, n, a, ia, ja, desca, info)
      CHARACTER, INTENT(IN)                              :: uplo
      INTEGER, INTENT(IN)                                :: n
      REAL, DIMENSION(:, :), TARGET                      :: a
      INTEGER, INTENT(IN)                                :: ia, ja
      INTEGER, DIMENSION(9)                              :: desca
      INTEGER, TARGET                                    :: info

      CHARACTER(len=*), PARAMETER                        :: routineN = 'cp_pspotrf_dlaf'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

#if defined(__DLAF)
      CALL dlaf_pspotrf(uplo, n, a, ia, ja, desca, info)
#else
      MARK_USED(uplo)
      MARK_USED(n)
      MARK_USED(a)
      MARK_USED(ia)
      MARK_USED(ja)
      MARK_USED(desca)
      MARK_USED(info)
      CPABORT("CP2K compiled without the DLA-Future library.")
#endif
      CALL timestop(handle)
   END SUBROUTINE cp_pspotrf_dlaf

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
! **************************************************************************************************
   SUBROUTINE cp_fm_diag_dlaf(matrix, eigenvectors, eigenvalues)

      TYPE(cp_fm_type), INTENT(IN)                       :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'cp_fm_diag_dlaf'

      INTEGER                                            :: handle, n, nmo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET   :: eig

      CALL timeset(routineN, handle)

      n = matrix%matrix_struct%nrow_global
      ALLOCATE (eig(n))

      CALL cp_fm_diag_dlaf_base(matrix, eigenvectors, eig)

      nmo = SIZE(eigenvalues, 1)
      IF (nmo > n) THEN
         eigenvalues(1:n) = eig(1:n)
      ELSE
         eigenvalues(1:nmo) = eig(1:nmo)
      END IF

      DEALLOCATE (eig)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_diag_dlaf

!***************************************************************************************************
!> \brief DLA-Future eigensolver
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \author Rocco Meli
! **************************************************************************************************
   SUBROUTINE cp_fm_diag_dlaf_base(matrix, eigenvectors, eigenvalues)
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix, eigenvectors
      REAL(kind=dp), DIMENSION(:), INTENT(OUT), TARGET   :: eigenvalues

      CHARACTER(len=*), PARAMETER :: dlaf_name = 'pdsyevd_dlaf', routineN = 'cp_fm_diag_dlaf_base'
      CHARACTER, PARAMETER                               :: uplo = 'L'

      CHARACTER(LEN=100)                                 :: message
      INTEGER                                            :: dlaf_handle, handle, n
      INTEGER, DIMENSION(9)                              :: desca, descz
      INTEGER, TARGET                                    :: info
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: a, z

      CALL timeset(routineN, handle)

#if defined(__DLAF)
      ! DLAF needs the lower triangular part
      ! Use eigenvectors matrix as workspace
      CALL cp_fm_uplo_to_full(matrix, eigenvectors)

      n = matrix%matrix_struct%nrow_global

      a => matrix%local_data
      z => eigenvectors%local_data

      desca(:) = matrix%matrix_struct%descriptor(:)
      descz(:) = eigenvectors%matrix_struct%descriptor(:)

      info = -1
      CALL timeset(dlaf_name, dlaf_handle)

      CALL dlaf_pdsyevd(uplo, n, a, 1, 1, desca, eigenvalues, z, 1, 1, descz, info)

      CALL timestop(dlaf_handle)

      IF (info /= 0) THEN
         WRITE (message, "(A,I0,A)") "ERROR in DLAF_PDSYEVD: Eigensolver failed (INFO = ", info, ")"
         CPABORT(TRIM(message))
      END IF
#else
      MARK_USED(a)
      MARK_USED(z)
      MARK_USED(desca)
      MARK_USED(descz)
      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(eigenvalues)
      MARK_USED(uplo)
      MARK_USED(n)
      MARK_USED(info)
      MARK_USED(dlaf_handle)
      MARK_USED(dlaf_name)
      MARK_USED(message)
      CPABORT("CP2K compiled without DLA-Future-Fortran library.")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_diag_dlaf_base

! **************************************************************************************************
!> \brief ...
!> \param a_matrix ...
!> \param b_matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \author Rocco Meli
! **************************************************************************************************
   SUBROUTINE cp_fm_diag_gen_dlaf(a_matrix, b_matrix, eigenvectors, eigenvalues)

      TYPE(cp_fm_type), INTENT(IN)                       :: a_matrix, b_matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues

      CHARACTER(LEN=*), PARAMETER :: routineN = 'cp_fm_diag_gen_dlaf'

      INTEGER                                            :: handle, n, nmo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET   :: eig

      CALL timeset(routineN, handle)

      n = a_matrix%matrix_struct%nrow_global
      ALLOCATE (eig(n))

      CALL cp_fm_diag_gen_dlaf_base(a_matrix, b_matrix, eigenvectors, eig)

      nmo = SIZE(eigenvalues, 1)
      IF (nmo > n) THEN
         eigenvalues(1:n) = eig(1:n)
      ELSE
         eigenvalues(1:nmo) = eig(1:nmo)
      END IF

      DEALLOCATE (eig)

      CALL timestop(handle)

   END SUBROUTINE cp_fm_diag_gen_dlaf

!***************************************************************************************************
!> \brief DLA-Future generalized eigensolver
!> \param a_matrix ...
!> \param b_matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
!> \author Rocco Meli
! **************************************************************************************************
   SUBROUTINE cp_fm_diag_gen_dlaf_base(a_matrix, b_matrix, eigenvectors, eigenvalues)
      TYPE(cp_fm_type), INTENT(IN)                       :: a_matrix, b_matrix, eigenvectors
      REAL(kind=dp), DIMENSION(:), INTENT(OUT), TARGET   :: eigenvalues

      CHARACTER(len=*), PARAMETER :: dlaf_name = 'pdsyevd_dlaf', &
         routineN = 'cp_fm_diag_gen_dlaf_base'
      CHARACTER, PARAMETER                               :: uplo = 'L'

      CHARACTER(LEN=100)                                 :: message
      INTEGER                                            :: dlaf_handle, handle, n
      INTEGER, DIMENSION(9)                              :: desca, descb, descz
      INTEGER, TARGET                                    :: info
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: a, b, z

      CALL timeset(routineN, handle)

#if defined(__DLAF)
      ! DLAF needs the lower triangular part
      ! Use eigenvectors matrix as workspace
      CALL cp_fm_uplo_to_full(a_matrix, eigenvectors)
      CALL cp_fm_uplo_to_full(b_matrix, eigenvectors)

      n = a_matrix%matrix_struct%nrow_global

      a => a_matrix%local_data
      b => b_matrix%local_data
      z => eigenvectors%local_data

      desca(:) = a_matrix%matrix_struct%descriptor(:)
      descb(:) = b_matrix%matrix_struct%descriptor(:)
      descz(:) = eigenvectors%matrix_struct%descriptor(:)

      info = -1
      CALL timeset(dlaf_name, dlaf_handle)

      CALL dlaf_pdsygvd(uplo, n, a, 1, 1, desca, b, 1, 1, descb, eigenvalues, z, 1, 1, descz, info)

      CALL timestop(dlaf_handle)

      IF (info /= 0) THEN
         WRITE (message, "(A,I0,A)") "ERROR in DLAF_PDSYGVD: Generalized Eigensolver failed (INFO = ", info, ")"
         CPABORT(TRIM(message))
      END IF
#else
      MARK_USED(a)
      MARK_USED(b)
      MARK_USED(z)
      MARK_USED(desca)
      MARK_USED(descb)
      MARK_USED(descz)
      MARK_USED(a_matrix)
      MARK_USED(b_matrix)
      MARK_USED(eigenvectors)
      MARK_USED(eigenvalues)
      MARK_USED(uplo)
      MARK_USED(n)
      MARK_USED(info)
      MARK_USED(dlaf_handle)
      MARK_USED(dlaf_name)
      MARK_USED(message)
      CPABORT("CP2K compiled without DLA-Future-Fortran library.")
#endif

      CALL timestop(handle)

   END SUBROUTINE cp_fm_diag_gen_dlaf_base

END MODULE cp_fm_dlaf_api
