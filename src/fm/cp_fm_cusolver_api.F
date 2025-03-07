!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper for cuSOLVERMp
!> \author Ole Schuett
! **************************************************************************************************
MODULE cp_fm_cusolver_api
   USE ISO_C_BINDING,                   ONLY: C_DOUBLE,&
                                              C_INT
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: cp_fm_diag_cusolver
   PUBLIC :: cp_fm_general_cusolver

CONTAINS

! **************************************************************************************************
!> \brief Driver routine to diagonalize a FM matrix with the cuSOLVERMp library.
!> \param matrix the matrix that is diagonalized
!> \param eigenvectors eigenvectors of the input matrix
!> \param eigenvalues eigenvalues of the input matrix
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE cp_fm_diag_cusolver(matrix, eigenvectors, eigenvalues)
      TYPE(cp_fm_type), INTENT(IN)                       :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_diag_cusolver'

      INTEGER                                            :: handle, n, nmo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues_buffer
      TYPE(cp_blacs_env_type), POINTER                   :: context
      INTERFACE
         SUBROUTINE cp_fm_diag_cusolver_c(fortran_comm, matrix_desc, &
                                          nprow, npcol, myprow, mypcol, &
                                          n, matrix, eigenvectors, eigenvalues) &
            BIND(C, name="cp_fm_diag_cusolver")
            IMPORT :: C_INT, C_DOUBLE
            INTEGER(kind=C_INT), VALUE                :: fortran_comm
            INTEGER(kind=C_INT), DIMENSION(*)         :: matrix_desc
            INTEGER(kind=C_INT), VALUE                :: nprow
            INTEGER(kind=C_INT), VALUE                :: npcol
            INTEGER(kind=C_INT), VALUE                :: myprow
            INTEGER(kind=C_INT), VALUE                :: mypcol
            INTEGER(kind=C_INT), VALUE                :: n
            REAL(kind=C_DOUBLE), DIMENSION(*)         :: matrix
            REAL(kind=C_DOUBLE), DIMENSION(*)         :: eigenvectors
            REAL(kind=C_DOUBLE), DIMENSION(*)         :: eigenvalues
         END SUBROUTINE cp_fm_diag_cusolver_c
      END INTERFACE

      CALL timeset(routineN, handle)

#if defined(__CUSOLVERMP)
      n = matrix%matrix_struct%nrow_global
      context => matrix%matrix_struct%context

      ! The passed eigenvalues array might be smaller than n.
      ALLOCATE (eigenvalues_buffer(n))

      CALL cp_fm_diag_cusolver_c( &
         fortran_comm=matrix%matrix_struct%para_env%get_handle(), &
         matrix_desc=matrix%matrix_struct%descriptor, &
         nprow=context%num_pe(1), &
         npcol=context%num_pe(2), &
         myprow=context%mepos(1), &
         mypcol=context%mepos(2), &
         n=matrix%matrix_struct%nrow_global, &
         matrix=matrix%local_data, &
         eigenvectors=eigenvectors%local_data, &
         eigenvalues=eigenvalues_buffer)

      nmo = SIZE(eigenvalues)
      eigenvalues(1:nmo) = eigenvalues_buffer(1:nmo)

#else
      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      eigenvalues = 0.0_dp
      MARK_USED(n)
      MARK_USED(nmo)
      MARK_USED(eigenvalues_buffer)
      MARK_USED(context)
      CPABORT("CP2K compiled without the cuSOLVERMp library.")
#endif

      CALL timestop(handle)
   END SUBROUTINE cp_fm_diag_cusolver

! **************************************************************************************************
!> \brief Driver routine to solve generalized eigenvalue problem A*x = lambda*B*x with cuSOLVERMp.
!> \param aMatrix the first matrix for the generalized eigenvalue problem
!> \param bMatrix the second matrix for the generalized eigenvalue problem
!> \param eigenvectors eigenvectors of the input matrix
!> \param eigenvalues eigenvalues of the input matrix
! **************************************************************************************************
   SUBROUTINE cp_fm_general_cusolver(aMatrix, bMatrix, eigenvectors, eigenvalues)
      USE ISO_C_BINDING, ONLY: C_INT, C_DOUBLE
      TYPE(cp_fm_type), INTENT(IN)                       :: aMatrix, bMatrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_fm_general_cusolver'

      INTEGER(kind=C_INT)                                :: handle, n, nmo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues_buffer
      TYPE(cp_blacs_env_type), POINTER                   :: context
      INTERFACE
         SUBROUTINE cp_fm_general_cusolver_c(fortran_comm, a_matrix_desc, b_matrix_desc, &
                                             nprow, npcol, myprow, mypcol, &
                                             n, aMatrix, bMatrix, eigenvectors, eigenvalues) &
            BIND(C, name="cp_fm_diag_cusolver_sygvd")
            IMPORT :: C_INT, C_DOUBLE
            INTEGER(kind=C_INT), VALUE                :: fortran_comm
            INTEGER(kind=C_INT), DIMENSION(*)         :: a_matrix_desc, b_matrix_desc
            INTEGER(kind=C_INT), VALUE                :: nprow
            INTEGER(kind=C_INT), VALUE                :: npcol
            INTEGER(kind=C_INT), VALUE                :: myprow
            INTEGER(kind=C_INT), VALUE                :: mypcol
            INTEGER(kind=C_INT), VALUE                :: n
            REAL(kind=C_DOUBLE), DIMENSION(*)         :: aMatrix
            REAL(kind=C_DOUBLE), DIMENSION(*)         :: bMatrix
            REAL(kind=C_DOUBLE), DIMENSION(*)         :: eigenvectors
            REAL(kind=C_DOUBLE), DIMENSION(*)         :: eigenvalues
         END SUBROUTINE cp_fm_general_cusolver_c
      END INTERFACE

      CALL timeset(routineN, handle)

#if defined(__CUSOLVERMP)
      n = INT(aMatrix%matrix_struct%nrow_global, C_INT)
      context => aMatrix%matrix_struct%context

      ! Allocate eigenvalues_buffer
      ALLOCATE (eigenvalues_buffer(n))

      CALL cp_fm_general_cusolver_c( &
         fortran_comm=INT(aMatrix%matrix_struct%para_env%get_handle(), C_INT), &
         a_matrix_desc=INT(aMatrix%matrix_struct%descriptor, C_INT), &
         b_matrix_desc=INT(bMatrix%matrix_struct%descriptor, C_INT), &
         nprow=INT(context%num_pe(1), C_INT), &
         npcol=INT(context%num_pe(2), C_INT), &
         myprow=INT(context%mepos(1), C_INT), &
         mypcol=INT(context%mepos(2), C_INT), &
         n=n, &
         aMatrix=aMatrix%local_data, &
         bMatrix=bMatrix%local_data, &
         eigenvectors=eigenvectors%local_data, &
         eigenvalues=eigenvalues_buffer)

      nmo = SIZE(eigenvalues)
      eigenvalues(1:nmo) = eigenvalues_buffer(1:nmo)

      DEALLOCATE (eigenvalues_buffer)
#else
      MARK_USED(aMatrix)
      MARK_USED(bMatrix)
      MARK_USED(eigenvectors)
      eigenvalues = 0.0_dp
      MARK_USED(n)
      MARK_USED(nmo)
      MARK_USED(eigenvalues_buffer)
      MARK_USED(context)
      CPABORT("CP2K compiled without the cuSOLVERMp library.")
#endif

      CALL timestop(handle)
   END SUBROUTINE cp_fm_general_cusolver

END MODULE cp_fm_cusolver_api
