!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Wrapper for cuSOLVERMp
!> \author Ole Schuett
! **************************************************************************************************
MODULE cp_fm_cusolver_api
   USE ISO_C_BINDING,                   ONLY: C_DOUBLE,&
                                              C_DOUBLE_COMPLEX,&
                                              C_INT,&
                                              C_LOC,&
                                              C_PTR
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_types,                    ONLY: cp_cfm_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: cp_fm_diag_cusolver
   PUBLIC :: cp_cfm_diag_cusolver
   PUBLIC :: cp_cfm_general_cusolver
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

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      CALL standard_cusolver(matrix%matrix_struct, C_LOC(matrix%local_data(1, 1)), &
                             C_LOC(eigenvectors%local_data(1, 1)), SIZE(eigenvectors%local_data, 1), eigenvalues, 0_C_INT)
      CALL timestop(handle)
   END SUBROUTINE cp_fm_diag_cusolver

! **************************************************************************************************
!> \brief Diagonalize a complex Hermitian matrix with cuSOLVERMp.
!> \param matrix Hermitian matrix (upper triangle)
!> \param eigenvectors eigenvectors of the input matrix
!> \param eigenvalues real eigenvalues in ascending order
! **************************************************************************************************
   SUBROUTINE cp_cfm_diag_cusolver(matrix, eigenvectors, eigenvalues)
      TYPE(cp_cfm_type), INTENT(IN)                      :: matrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_diag_cusolver'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      CALL standard_cusolver(matrix%matrix_struct, C_LOC(matrix%local_data(1, 1)), &
                             C_LOC(eigenvectors%local_data(1, 1)), SIZE(eigenvectors%local_data, 1), eigenvalues, 1_C_INT)
      CALL timestop(handle)
   END SUBROUTINE cp_cfm_diag_cusolver

! **************************************************************************************************
!> \brief Call the type-generic cuSOLVERMp standard eigensolver on a BLACS grid.
!> \param layout input matrix distribution
!> \param matrix local input data
!> \param eigenvectors local output data with the same block distribution
!> \param ldz local leading dimension of eigenvectors
!> \param eigenvalues requested prefix of the spectrum
!> \param is_complex nonzero for double complex matrices
! **************************************************************************************************
   SUBROUTINE standard_cusolver(layout, matrix, eigenvectors, ldz, eigenvalues, is_complex)
      TYPE(cp_fm_struct_type), INTENT(IN)                :: layout
      TYPE(C_PTR), INTENT(IN)                            :: matrix, eigenvectors
      INTEGER, INTENT(IN)                                :: ldz
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues
      INTEGER(KIND=C_INT), INTENT(IN)                    :: is_complex

#if defined(__CUSOLVERMP)
      INTERFACE
         SUBROUTINE cp_fm_diag_cusolver_c(cusolver_context, fortran_comm, matrix_desc, &
                                          nprow, npcol, myprow, mypcol, &
                                          matrix, eigenvectors, ldz, eigenvalues, neig, is_complex) &
            BIND(C, name="cp_fm_diag_cusolver")
            IMPORT :: C_INT, C_DOUBLE, C_PTR
            TYPE(C_PTR), INTENT(INOUT)                 :: cusolver_context
            INTEGER(kind=C_INT), VALUE                :: fortran_comm
            INTEGER(kind=C_INT), DIMENSION(*)         :: matrix_desc
            INTEGER(kind=C_INT), VALUE                :: nprow
            INTEGER(kind=C_INT), VALUE                :: npcol
            INTEGER(kind=C_INT), VALUE                :: myprow
            INTEGER(kind=C_INT), VALUE                :: mypcol
            TYPE(C_PTR), VALUE                        :: matrix, eigenvectors
            INTEGER(kind=C_INT), VALUE                :: ldz, neig, is_complex
            REAL(kind=C_DOUBLE), DIMENSION(*)         :: eigenvalues
         END SUBROUTINE cp_fm_diag_cusolver_c
      END INTERFACE

      CALL cp_fm_diag_cusolver_c( &
         cusolver_context=layout%context%cusolver_context, fortran_comm=INT(layout%para_env%get_handle(), C_INT), &
         matrix_desc=INT(layout%descriptor, C_INT), nprow=INT(layout%context%num_pe(1), C_INT), &
         npcol=INT(layout%context%num_pe(2), C_INT), myprow=INT(layout%context%mepos(1), C_INT), &
         mypcol=INT(layout%context%mepos(2), C_INT), matrix=matrix, eigenvectors=eigenvectors, ldz=INT(ldz, C_INT), &
         eigenvalues=eigenvalues, neig=INT(MIN(SIZE(eigenvalues), layout%nrow_global), C_INT), is_complex=is_complex)
#else
      MARK_USED(layout)
      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(ldz)
      MARK_USED(is_complex)
      eigenvalues = 0.0_dp
      CPABORT("CP2K compiled without the cuSOLVERMp library.")
#endif
   END SUBROUTINE standard_cusolver

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

      INTEGER(kind=C_INT)                                :: handle, n
      TYPE(cp_blacs_env_type), POINTER                   :: context
      INTERFACE
         SUBROUTINE cp_fm_general_cusolver_c(cusolver_context, fortran_comm, a_matrix_desc, b_matrix_desc, &
                                             nprow, npcol, myprow, mypcol, &
                                             n, nmo, aMatrix, bMatrix, eigenvectors, eigenvalues) &
            BIND(C, name="cp_fm_diag_cusolver_sygvd")
            IMPORT :: C_INT, C_DOUBLE, C_PTR
            TYPE(C_PTR), INTENT(INOUT)                 :: cusolver_context
            INTEGER(kind=C_INT), VALUE                :: fortran_comm
            INTEGER(kind=C_INT), DIMENSION(*)         :: a_matrix_desc, b_matrix_desc
            INTEGER(kind=C_INT), VALUE                :: nprow
            INTEGER(kind=C_INT), VALUE                :: npcol
            INTEGER(kind=C_INT), VALUE                :: myprow
            INTEGER(kind=C_INT), VALUE                :: mypcol
            INTEGER(kind=C_INT), VALUE                :: n
            INTEGER(kind=C_INT), VALUE                :: nmo
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

      CALL cp_fm_general_cusolver_c( &
         cusolver_context=context%cusolver_context, &
         fortran_comm=INT(aMatrix%matrix_struct%para_env%get_handle(), C_INT), &
         a_matrix_desc=INT(aMatrix%matrix_struct%descriptor, C_INT), &
         b_matrix_desc=INT(bMatrix%matrix_struct%descriptor, C_INT), &
         nprow=INT(context%num_pe(1), C_INT), &
         npcol=INT(context%num_pe(2), C_INT), &
         myprow=INT(context%mepos(1), C_INT), &
         mypcol=INT(context%mepos(2), C_INT), &
         n=n, &
         nmo=INT(SIZE(eigenvalues), C_INT), &
         aMatrix=aMatrix%local_data, &
         bMatrix=bMatrix%local_data, &
         eigenvectors=eigenvectors%local_data, &
         eigenvalues=eigenvalues)
#else
      MARK_USED(aMatrix)
      MARK_USED(bMatrix)
      MARK_USED(eigenvectors)
      eigenvalues = 0.0_dp
      MARK_USED(n)
      MARK_USED(context)
      CPABORT("CP2K compiled without the cuSOLVERMp library.")
#endif

      CALL timestop(handle)
   END SUBROUTINE cp_fm_general_cusolver

! **************************************************************************************************
!> \brief Driver routine to solve generalized complex eigenvalue problem A*x = lambda*B*x with
!>        cuSOLVERMp.
!> \param aMatrix the first matrix for the generalized eigenvalue problem
!> \param bMatrix the second matrix for the generalized eigenvalue problem
!> \param eigenvectors eigenvectors of the input matrix
!> \param eigenvalues eigenvalues of the input matrix
! **************************************************************************************************
   SUBROUTINE cp_cfm_general_cusolver(aMatrix, bMatrix, eigenvectors, eigenvalues)
      USE ISO_C_BINDING, ONLY: C_DOUBLE, C_DOUBLE_COMPLEX, C_INT
      TYPE(cp_cfm_type), INTENT(IN)                      :: aMatrix, bMatrix, eigenvectors
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues

      CHARACTER(len=*), PARAMETER :: routineN = 'cp_cfm_general_cusolver'

      INTEGER(kind=C_INT)                                :: handle, n
      TYPE(cp_blacs_env_type), POINTER                   :: context
      INTERFACE
         SUBROUTINE cp_cfm_general_cusolver_c(cusolver_context, fortran_comm, a_matrix_desc, b_matrix_desc, &
                                              nprow, npcol, myprow, mypcol, &
                                              n, nmo, aMatrix, bMatrix, eigenvectors, eigenvalues) &
            BIND(C, name="cp_cfm_diag_cusolver_hegvd")
            IMPORT :: C_DOUBLE, C_DOUBLE_COMPLEX, C_INT, C_PTR
            TYPE(C_PTR), INTENT(INOUT)                 :: cusolver_context
            INTEGER(kind=C_INT), VALUE                :: fortran_comm
            INTEGER(kind=C_INT), DIMENSION(*)         :: a_matrix_desc, b_matrix_desc
            INTEGER(kind=C_INT), VALUE                :: nprow
            INTEGER(kind=C_INT), VALUE                :: npcol
            INTEGER(kind=C_INT), VALUE                :: myprow
            INTEGER(kind=C_INT), VALUE                :: mypcol
            INTEGER(kind=C_INT), VALUE                :: n
            INTEGER(kind=C_INT), VALUE                :: nmo
            COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(*) :: aMatrix
            COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(*) :: bMatrix
            COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(*) :: eigenvectors
            REAL(kind=C_DOUBLE), DIMENSION(*)         :: eigenvalues
         END SUBROUTINE cp_cfm_general_cusolver_c
      END INTERFACE

      CALL timeset(routineN, handle)

#if defined(__CUSOLVERMP)
      n = INT(aMatrix%matrix_struct%nrow_global, C_INT)
      context => aMatrix%matrix_struct%context

      CALL cp_cfm_general_cusolver_c( &
         cusolver_context=context%cusolver_context, &
         fortran_comm=INT(aMatrix%matrix_struct%para_env%get_handle(), C_INT), &
         a_matrix_desc=INT(aMatrix%matrix_struct%descriptor, C_INT), &
         b_matrix_desc=INT(bMatrix%matrix_struct%descriptor, C_INT), &
         nprow=INT(context%num_pe(1), C_INT), &
         npcol=INT(context%num_pe(2), C_INT), &
         myprow=INT(context%mepos(1), C_INT), &
         mypcol=INT(context%mepos(2), C_INT), &
         n=n, &
         nmo=INT(SIZE(eigenvalues), C_INT), &
         aMatrix=aMatrix%local_data, &
         bMatrix=bMatrix%local_data, &
         eigenvectors=eigenvectors%local_data, &
         eigenvalues=eigenvalues)
#else
      MARK_USED(aMatrix)
      MARK_USED(bMatrix)
      MARK_USED(eigenvectors)
      eigenvalues = 0.0_dp
      MARK_USED(n)
      MARK_USED(context)
      CPABORT("CP2K compiled without the cuSOLVERMp library.")
#endif

      CALL timestop(handle)
   END SUBROUTINE cp_cfm_general_cusolver

END MODULE cp_fm_cusolver_api
