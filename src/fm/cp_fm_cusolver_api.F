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
                                              C_INT,&
                                              C_LOC,&
                                              C_NULL_PTR,&
                                              C_PTR
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

#if defined(__CUSOLVERMP)
   INTERFACE
      FUNCTION uses_dense(n, nprow, npcol) BIND(C, name="cp_fm_cusolver_uses_dense")
         IMPORT :: C_INT
         INTEGER(KIND=C_INT), VALUE :: n, nprow, npcol
         INTEGER(KIND=C_INT) :: uses_dense
      END FUNCTION uses_dense
   END INTERFACE
#endif

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
      CALL solve_cusolver(matrix%matrix_struct, C_LOC(matrix%local_data), &
                          C_LOC(eigenvectors%local_data), SIZE(eigenvectors%local_data, 1), eigenvalues, 0_C_INT, &
                          C_NULL_PTR, -1)
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
      CALL solve_cusolver(matrix%matrix_struct, C_LOC(matrix%local_data), &
                          C_LOC(eigenvectors%local_data), SIZE(eigenvectors%local_data, 1), eigenvalues, 1_C_INT, &
                          C_NULL_PTR, -1)
      CALL timestop(handle)
   END SUBROUTINE cp_cfm_diag_cusolver

! **************************************************************************************************
!> \brief Call the real or complex cuSOLVER eigensolver on a BLACS grid.
!> \param layout input matrix distribution
!> \param matrix local input data
!> \param eigenvectors local output data with the same block distribution
!> \param ldz local leading dimension of eigenvectors
!> \param eigenvalues requested prefix of the spectrum
!> \param overlap local overlap matrix, or a null pointer for a standard solve
!> \param ldb local leading dimension of the overlap, or -1 for a standard solve
!> \param is_complex nonzero for double complex matrices
! **************************************************************************************************
   SUBROUTINE solve_cusolver(layout, matrix, eigenvectors, ldz, eigenvalues, is_complex, overlap, ldb)
      TYPE(cp_fm_struct_type), INTENT(IN)                :: layout
      TYPE(C_PTR), INTENT(IN)                            :: matrix, eigenvectors, overlap
      INTEGER, INTENT(IN)                                :: ldz, ldb
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: eigenvalues
      INTEGER(KIND=C_INT), INTENT(IN)                    :: is_complex

#if defined(__CUSOLVERMP)
      INTERFACE
         SUBROUTINE cp_fm_diag_cusolver_c(cusolver_context, fortran_comm, matrix_desc, &
                                          nprow, npcol, myprow, mypcol, &
                                          matrix, overlap, ldb, eigenvectors, ldz, eigenvalues, neig, is_complex) &
            BIND(C, name="cp_fm_cusolver_solve")
            IMPORT :: C_INT, C_DOUBLE, C_PTR
            TYPE(C_PTR), INTENT(INOUT)                 :: cusolver_context
            INTEGER(kind=C_INT), VALUE                :: fortran_comm
            INTEGER(kind=C_INT), DIMENSION(*)         :: matrix_desc
            INTEGER(kind=C_INT), VALUE                :: nprow
            INTEGER(kind=C_INT), VALUE                :: npcol
            INTEGER(kind=C_INT), VALUE                :: myprow
            INTEGER(kind=C_INT), VALUE                :: mypcol
            TYPE(C_PTR), VALUE                        :: matrix, overlap, eigenvectors
            INTEGER(kind=C_INT), VALUE                :: ldb, ldz, neig, is_complex
            REAL(kind=C_DOUBLE), DIMENSION(*)         :: eigenvalues
         END SUBROUTINE cp_fm_diag_cusolver_c
      END INTERFACE

      CALL cp_fm_diag_cusolver_c( &
         cusolver_context=layout%context%cusolver_context, fortran_comm=INT(layout%para_env%get_handle(), C_INT), &
         matrix_desc=INT(layout%descriptor, C_INT), nprow=INT(layout%context%num_pe(1), C_INT), &
         npcol=INT(layout%context%num_pe(2), C_INT), myprow=INT(layout%context%mepos(1), C_INT), &
         mypcol=INT(layout%context%mepos(2), C_INT), matrix=matrix, overlap=overlap, ldb=INT(ldb, C_INT), &
         eigenvectors=eigenvectors, ldz=INT(ldz, C_INT), &
         eigenvalues=eigenvalues, neig=INT(MIN(SIZE(eigenvalues), layout%nrow_global), C_INT), is_complex=is_complex)
#else
      MARK_USED(layout)
      MARK_USED(matrix)
      MARK_USED(eigenvectors)
      MARK_USED(overlap)
      MARK_USED(ldb)
      MARK_USED(ldz)
      MARK_USED(is_complex)
      eigenvalues = 0.0_dp
      CPABORT("CP2K compiled without the cuSOLVERMp library.")
#endif
   END SUBROUTINE solve_cusolver

   #:for fm, scalar, complex_flag in [("fm", "REAL", 0), ("cfm", "COMPLEX", 1)]
! **************************************************************************************************
!> \brief Solve A*x = lambda*B*x with cuSOLVER.
!> \param aMatrix symmetric or Hermitian matrix (upper triangle)
!> \param bMatrix positive definite metric (upper triangle)
!> \param eigenvectors eigenvectors of the input matrix
!> \param eigenvalues requested prefix of the spectrum
! **************************************************************************************************
      SUBROUTINE cp_${fm}$_general_cusolver(aMatrix, bMatrix, eigenvectors, eigenvalues)
         TYPE(cp_${fm}$_type), INTENT(IN) :: aMatrix, bMatrix, eigenvectors
         REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: eigenvalues

         CHARACTER(len=*), PARAMETER :: routineN = 'cp_${fm}$_general_cusolver'
         INTEGER :: handle, ldb
         TYPE(C_PTR) :: b
#if defined(__CUSOLVERMP)
         ${scalar}$(KIND=dp), ALLOCATABLE, TARGET, DIMENSION(:, :) :: lower
         INTEGER :: n
#endif

         CALL timeset(routineN, handle)
         b = C_LOC(bMatrix%local_data)
         ldb = SIZE(bMatrix%local_data, 1)
#if defined(__CUSOLVERMP)
         n = aMatrix%matrix_struct%nrow_global
         ASSOCIATE (grid => aMatrix%matrix_struct%context%num_pe)
            IF (uses_dense(INT(n, C_INT), INT(grid(1), C_INT), INT(grid(2), C_INT)) == 0) THEN
               ! Mp type-1 eigenvectors require lower A and B. Reuse A's undefined lower
               ! triangle and a single scratch tile in A's distribution for B.
               ALLOCATE (lower, MOLD=aMatrix%local_data)
               CALL lower_triangle(aMatrix)
               ASSOCIATE (desc => aMatrix%matrix_struct%descriptor)
                  CALL p${'z' if complex_flag else 'd'}$lacpy('L', n, n, lower, 1, 1, desc, aMatrix%local_data, 1, 1, desc)
               END ASSOCIATE
               CALL lower_triangle(bMatrix)
               b = C_LOC(lower)
               ldb = SIZE(lower, 1)
            END IF
         END ASSOCIATE
#endif
         CALL solve_cusolver(aMatrix%matrix_struct, C_LOC(aMatrix%local_data), C_LOC(eigenvectors%local_data), &
                             SIZE(eigenvectors%local_data, 1), eigenvalues, ${complex_flag}$_C_INT, b, ldb)
         CALL timestop(handle)
#if defined(__CUSOLVERMP)

      CONTAINS

         SUBROUTINE lower_triangle(matrix)
            TYPE(cp_${fm}$_type), INTENT(IN) :: matrix

            ASSOCIATE (src => matrix%matrix_struct%descriptor, dst => aMatrix%matrix_struct%descriptor)
               #:if complex_flag
                  CALL pztradd('L', 'C', n, n, (1.0_dp, 0.0_dp), matrix%local_data, 1, 1, src, &
                               (0.0_dp, 0.0_dp), lower, 1, 1, dst)
               #:else
                  CALL pdtradd('L', 'T', n, n, 1.0_dp, matrix%local_data, 1, 1, src, &
                               0.0_dp, lower, 1, 1, dst)
               #:endif
            END ASSOCIATE
         END SUBROUTINE lower_triangle
#endif
      END SUBROUTINE cp_${fm}$_general_cusolver
   #:endfor

END MODULE cp_fm_cusolver_api
