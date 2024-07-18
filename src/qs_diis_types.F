!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief buffer for the diis of the scf
!> \par History
!>      02.2003 revamped [fawzi]
!> \author Matthias Krack
! **************************************************************************************************
MODULE qs_diis_types
   USE cp_cfm_types,                    ONLY: cp_cfm_release,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_release
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_diis_types'

   PUBLIC :: qs_diis_buffer_type
   PUBLIC :: qs_diis_b_release
   PUBLIC :: qs_diis_buffer_type_sparse, &
             qs_diis_b_release_sparse, &
             qs_diis_buffer_type_kp, &
             qs_diis_b_release_kp

! **************************************************************************************************
!> \brief keeps a buffer with the previous values of s,p,k
!> \par History
!>      02.2003 rewamped [fawzi]
!> \author Matthias Krack
! **************************************************************************************************
   TYPE qs_diis_buffer_type
      INTEGER                                          :: nbuffer = -1
      INTEGER                                          :: ncall = -1
      TYPE(cp_fm_type), DIMENSION(:, :), POINTER       :: error => NULL()
      TYPE(cp_fm_type), DIMENSION(:, :), POINTER       :: param => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER          :: b_matrix => NULL()
   END TYPE qs_diis_buffer_type

! **************************************************************************************************
!> \brief build array of pointers to diis buffers
!> \param diis_buffer the diis buffer pointer
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   TYPE qs_diis_buffer_p_type
      TYPE(qs_diis_buffer_type), POINTER :: diis_buffer => NULL()
   END TYPE qs_diis_buffer_p_type

! **************************************************************************************************
!> \brief build array of pointers to diis buffers for sparse matrix case
!> \param diis_buffer the diis buffer pointer
!> \par History
!>      10.2014 Modified from non-sparse case by Fredy W. Aquino
!> \author fwaq
! **************************************************************************************************
   TYPE qs_diis_buffer_type_sparse
      INTEGER                                          :: nbuffer = -1
      INTEGER                                          :: ncall = -1
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER     :: error => NULL()
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER     :: param => NULL()
      REAL(KIND=dp), DIMENSION(:, :), POINTER          :: b_matrix => NULL()
   END TYPE qs_diis_buffer_type_sparse

   TYPE qs_diis_buffer_p_type_sparse
      TYPE(qs_diis_buffer_type_sparse), POINTER        :: diis_buffer => NULL()
   END TYPE qs_diis_buffer_p_type_sparse

! **************************************************************************************************
!> \brief build array of pointers to diis buffers in the k-point (complex full matrices) case
!> \note we keep the overlap matrix and 1 parameter matrix per k-point, in the corresponding kp subgroup
!> \par History
!>      04.2023 created
!> \author Augustin Bussy
! **************************************************************************************************
   TYPE qs_diis_buffer_type_kp
      INTEGER                                          :: nbuffer = -1
      INTEGER                                          :: ncall = -1
      TYPE(cp_cfm_type), DIMENSION(:), POINTER         :: smat => NULL()
      TYPE(cp_cfm_type), DIMENSION(:, :, :), POINTER   :: param => NULL()
      TYPE(cp_cfm_type), DIMENSION(:, :, :), POINTER   :: error => NULL()
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER       :: b_matrix => NULL()
   END TYPE qs_diis_buffer_type_kp

   TYPE qs_diis_buffer_p_type_kp
      TYPE(qs_diis_buffer_type_kp), POINTER            :: diis_buffer => NULL()
   END TYPE qs_diis_buffer_p_type_kp

CONTAINS

! **************************************************************************************************
!> \brief releases the given diis buffer (see doc/ReferenceCounting.html)
!> \param diis_buffer the buffer to release
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE qs_diis_b_release(diis_buffer)
      TYPE(qs_diis_buffer_type), INTENT(INOUT)           :: diis_buffer

      IF (ASSOCIATED(diis_buffer%b_matrix)) THEN
         DEALLOCATE (diis_buffer%b_matrix)
      END IF
      CALL cp_fm_release(diis_buffer%error)
      CALL cp_fm_release(diis_buffer%param)

   END SUBROUTINE qs_diis_b_release

! **************************************************************************************************
!> \brief releases the given diis buffer (see doc/ReferenceCounting.html)
!> \param diis_buffer the buffer to release
!> \par History
!>      10-11-14 created [FA] modified from qs_diis_b_release
!> \author Fredy W. Aquino
! **************************************************************************************************
   SUBROUTINE qs_diis_b_release_sparse(diis_buffer)

      TYPE(qs_diis_buffer_type_sparse), INTENT(INOUT)    :: diis_buffer

      INTEGER                                            :: i, j

      IF (ASSOCIATED(diis_buffer%b_matrix)) THEN
         DEALLOCATE (diis_buffer%b_matrix)
      END IF
      IF (ASSOCIATED(diis_buffer%error)) THEN
         DO j = 1, SIZE(diis_buffer%error, 2)
            DO i = 1, SIZE(diis_buffer%error, 1)
               CALL dbcsr_release(diis_buffer%error(i, j)%matrix)
               DEALLOCATE (diis_buffer%error(i, j)%matrix)
            END DO
         END DO
         DEALLOCATE (diis_buffer%error)
      END IF
      IF (ASSOCIATED(diis_buffer%param)) THEN
         DO j = 1, SIZE(diis_buffer%param, 2)
            DO i = 1, SIZE(diis_buffer%param, 1)
               CALL dbcsr_release(diis_buffer%param(i, j)%matrix)
               DEALLOCATE (diis_buffer%param(i, j)%matrix)
            END DO
         END DO
         DEALLOCATE (diis_buffer%param)
      END IF
   END SUBROUTINE qs_diis_b_release_sparse

! **************************************************************************************************
!> \brief releases the given diis KP buffer
!> \param diis_buffer the buffer to release
! **************************************************************************************************
   SUBROUTINE qs_diis_b_release_kp(diis_buffer)
      TYPE(qs_diis_buffer_type_kp), INTENT(INOUT)        :: diis_buffer

      INTEGER                                            :: i, j, k

      IF (ASSOCIATED(diis_buffer%b_matrix)) THEN
         DEALLOCATE (diis_buffer%b_matrix)
      END IF
      IF (ASSOCIATED(diis_buffer%smat)) THEN
         DO i = 1, SIZE(diis_buffer%smat)
            CALL cp_cfm_release(diis_buffer%smat(i))
         END DO
         DEALLOCATE (diis_buffer%smat)
      END IF
      IF (ASSOCIATED(diis_buffer%error)) THEN
         DO k = 1, SIZE(diis_buffer%error, 3)
            DO j = 1, SIZE(diis_buffer%error, 2)
               DO i = 1, SIZE(diis_buffer%error, 1)
                  CALL cp_cfm_release(diis_buffer%error(i, j, k))
               END DO
            END DO
         END DO
         DEALLOCATE (diis_buffer%error)
      END IF
      IF (ASSOCIATED(diis_buffer%param)) THEN
         DO k = 1, SIZE(diis_buffer%param, 3)
            DO j = 1, SIZE(diis_buffer%param, 2)
               DO i = 1, SIZE(diis_buffer%param, 1)
                  CALL cp_cfm_release(diis_buffer%param(i, j, k))
               END DO
            END DO
         END DO
         DEALLOCATE (diis_buffer%param)
      END IF
   END SUBROUTINE qs_diis_b_release_kp
END MODULE qs_diis_types
