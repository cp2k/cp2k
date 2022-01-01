!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Set of routines to:
!>        Contract integrals over primitive Gaussians
!>        Decontract (density) matrices
!>        Trace matrices to get forces
!>        Block copy and add matrices
!> \par History
!>      Replace dgemm by MATMUL: Massive speedups in openMP loops (JGH, 12.2019)
!> \author JGH (01.07.2014)
! **************************************************************************************************
MODULE ai_contraction

   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ai_contraction'

   PRIVATE

   PUBLIC :: contraction, decontraction, block_add, force_trace

   INTERFACE contraction
      MODULE PROCEDURE contraction_ab, contraction_abc
   END INTERFACE

   INTERFACE decontraction
      MODULE PROCEDURE decontraction_ab
   END INTERFACE

   INTERFACE force_trace
      MODULE PROCEDURE force_trace_ab
   END INTERFACE

   INTERFACE block_add
      MODULE PROCEDURE block_add_ab
   END INTERFACE

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Applying the contraction coefficients to a set of two-center primitive
!>        integrals
!>        QAB <- CA(T) * SAB * CB
!>        QAB is optionally scaled with "fscale"
!>        Variable "trans" requests the output to be QAB(T)
!>        If only one of the transformation matrix is given, only a half
!>        transformation is done
!>        Active dimensions are: QAB(ma,mb), SAB(na,nb)
!> \param sab     Input matrix, dimension(:,:)
!> \param qab     Output matrix, dimension(:,:)
!> \param ca      Left transformation matrix, optional
!> \param na      First dimension of ca, optional
!> \param ma      Second dimension of ca, optional
!> \param cb      Right transformation matrix, optional
!> \param nb      First dimension of cb, optional
!> \param mb      Second dimension of cb, optional
!> \param fscale  Optional scaling of output
!> \param trans   Optional transposition of output
! **************************************************************************************************
   SUBROUTINE contraction_ab(sab, qab, ca, na, ma, cb, nb, mb, fscale, trans)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: qab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: ca
      INTEGER, INTENT(IN), OPTIONAL                      :: na, ma
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: cb
      INTEGER, INTENT(IN), OPTIONAL                      :: nb, mb
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: fscale
      LOGICAL, INTENT(IN), OPTIONAL                      :: trans

      INTEGER                                            :: lda, ldb, ldq, lds, ldw, mal, mbl, nal, &
                                                            nbl
      LOGICAL                                            :: my_trans
      REAL(KIND=dp)                                      :: fs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work

! Should output matrix be transposed?

      IF (PRESENT(trans)) THEN
         my_trans = trans
      ELSE
         my_trans = .FALSE.
      END IF

      ! Scaling of output matrix
      IF (PRESENT(fscale)) THEN
         fs = fscale
      ELSE
         fs = 1.0_dp
      END IF

      ! Active matrix size
      IF (PRESENT(ca)) THEN
         IF (PRESENT(na)) THEN
            nal = na
         ELSE
            nal = SIZE(ca, 1)
         END IF
         IF (PRESENT(ma)) THEN
            mal = ma
         ELSE
            mal = SIZE(ca, 2)
         END IF
         lda = SIZE(ca, 1)
      END IF
      IF (PRESENT(cb)) THEN
         IF (PRESENT(nb)) THEN
            nbl = nb
         ELSE
            nbl = SIZE(cb, 1)
         END IF
         IF (PRESENT(mb)) THEN
            mbl = mb
         ELSE
            mbl = SIZE(cb, 2)
         END IF
         ldb = SIZE(cb, 1)
      END IF

      lds = SIZE(sab, 1)
      ldq = SIZE(qab, 1)

      IF (PRESENT(ca) .AND. PRESENT(cb)) THEN
         ! Full transform
         ALLOCATE (work(nal, mbl))
         ldw = nal
!dg      CALL dgemm("N", "N", nal, mbl, nbl, 1.0_dp, sab(1, 1), lds, cb(1, 1), ldb, 0.0_dp, work(1, 1), ldw)
         work(1:nal, 1:mbl) = MATMUL(sab(1:nal, 1:nbl), cb(1:nbl, 1:mbl))
         IF (my_trans) THEN
!dg         CALL dgemm("T", "N", mbl, mal, nal, fs, work(1, 1), ldw, ca(1, 1), lda, 0.0_dp, qab(1, 1), ldq)
            qab(1:mbl, 1:mal) = fs*MATMUL(TRANSPOSE(work(1:nal, 1:mbl)), ca(1:nal, 1:mal))
         ELSE
!dg         CALL dgemm("T", "N", mal, mbl, nal, fs, ca(1, 1), lda, work(1, 1), ldw, 0.0_dp, qab(1, 1), ldq)
            qab(1:mal, 1:mbl) = fs*MATMUL(TRANSPOSE(ca(1:nal, 1:mal)), work(1:nal, 1:mbl))
         END IF
         DEALLOCATE (work)
      ELSE IF (PRESENT(ca)) THEN
         IF (PRESENT(nb)) THEN
            nbl = nb
         ELSE
            nbl = SIZE(sab, 2)
         END IF
         IF (my_trans) THEN
!dg         CALL dgemm("T", "N", nbl, mal, nal, fs, sab(1, 1), lds, ca(1, 1), lda, 0.0_dp, qab(1, 1), ldq)
            qab(1:nbl, 1:mal) = fs*MATMUL(TRANSPOSE(sab(1:nal, 1:nbl)), ca(1:nal, 1:mal))
         ELSE
!dg         CALL dgemm("T", "N", mal, nbl, nal, fs, ca(1, 1), lda, sab(1, 1), lds, 0.0_dp, qab(1, 1), ldq)
            qab(1:mal, 1:nbl) = fs*MATMUL(TRANSPOSE(ca(1:nal, 1:mal)), sab(1:nal, 1:nbl))
         END IF
      ELSE IF (PRESENT(cb)) THEN
         IF (PRESENT(na)) THEN
            nal = na
         ELSE
            nal = SIZE(sab, 1)
         END IF
         IF (my_trans) THEN
!dg         CALL dgemm("N", "N", nal, mbl, nbl, fs, sab(1, 1), lds, cb(1, 1), ldb, 0.0_dp, qab, ldq)
            qab(1:nal, 1:mbl) = fs*MATMUL(sab(1:nal, 1:nbl), cb(1:nbl, 1:mbl))
         ELSE
!dg         CALL dgemm("T", "T", mbl, nal, nbl, fs, cb(1, 1), ldb, sab(1, 1), lds, 0.0_dp, qab, ldq)
            qab(1:mbl, 1:nal) = fs*MATMUL(TRANSPOSE(cb(1:nbl, 1:mbl)), TRANSPOSE(sab(1:nal, 1:nbl)))
         END IF
      ELSE
         ! Copy of arrays is not covered here
         CPABORT("Copy of arrays is not covered here")
      END IF

   END SUBROUTINE contraction_ab

! **************************************************************************************************
!> \brief Applying the contraction coefficients to a tripple set integrals
!>        QABC <- CA(T) * SABC * CB * CC
!>        If only one or two of the transformation matrices are given, only a
!>        part transformation is done
!> \param sabc    Input matrix, dimension(:,:)
!> \param qabc    Output matrix, dimension(:,:)
!> \param ca      Transformation matrix (index 1), optional
!> \param na      First dimension of ca, optional
!> \param ma      Second dimension of ca, optional
!> \param cb      Transformation matrix (index 2), optional
!> \param nb      First dimension of cb, optional
!> \param mb      Second dimension of cb, optional
!> \param cc      Transformation matrix (index 3), optional
!> \param nc      First dimension of cc, optional
!> \param mc      Second dimension of cc, optional
! **************************************************************************************************
   SUBROUTINE contraction_abc(sabc, qabc, ca, na, ma, cb, nb, mb, cc, nc, mc)

      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: sabc
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: qabc
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: ca
      INTEGER, INTENT(IN), OPTIONAL                      :: na, ma
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: cb
      INTEGER, INTENT(IN), OPTIONAL                      :: nb, mb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: cc
      INTEGER, INTENT(IN), OPTIONAL                      :: nc, mc

      INTEGER                                            :: lda, ldb, ldc, mal, mbl, mcl, nal, nbl, &
                                                            ncl
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: work1, work2, work3, work4

! Active matrix size

      IF (PRESENT(ca)) THEN
         IF (PRESENT(na)) THEN
            nal = na
         ELSE
            nal = SIZE(ca, 1)
         END IF
         IF (PRESENT(ma)) THEN
            mal = ma
         ELSE
            mal = SIZE(ca, 2)
         END IF
         lda = SIZE(ca, 1)
      END IF
      IF (PRESENT(cb)) THEN
         IF (PRESENT(nb)) THEN
            nbl = nb
         ELSE
            nbl = SIZE(cb, 1)
         END IF
         IF (PRESENT(mb)) THEN
            mbl = mb
         ELSE
            mbl = SIZE(cb, 2)
         END IF
         ldb = SIZE(cb, 1)
      END IF
      IF (PRESENT(cc)) THEN
         IF (PRESENT(nc)) THEN
            ncl = nc
         ELSE
            ncl = SIZE(cc, 1)
         END IF
         IF (PRESENT(mc)) THEN
            mcl = mc
         ELSE
            mcl = SIZE(cc, 2)
         END IF
         ldc = SIZE(cc, 1)
      END IF

      IF (PRESENT(ca) .AND. PRESENT(cb) .AND. PRESENT(cc)) THEN
         ! Full transform
         ALLOCATE (work1(nal, nbl, ncl))
         ! make sure that we have contiguous memory, needed for transpose algorithm
         work1(1:nal, 1:nbl, 1:ncl) = sabc(1:nal, 1:nbl, 1:ncl)
         !
         ALLOCATE (work2(nbl, ncl, mal))
         CALL dgemm("T", "N", nbl*ncl, mal, nal, 1.0_dp, work1(1, 1, 1), nal, ca(1, 1), lda, &
                    0.0_dp, work2(1, 1, 1), nbl*ncl)
         !
         ALLOCATE (work3(ncl, mal, mbl))
         CALL dgemm("T", "N", ncl*mal, mbl, nbl, 1.0_dp, work2(1, 1, 1), nbl, cb(1, 1), ldb, &
                    0.0_dp, work3(1, 1, 1), ncl*mal)
         !
         ALLOCATE (work4(mal, mbl, mcl))
         CALL dgemm("T", "N", mal*mbl, mcl, ncl, 1.0_dp, work3(1, 1, 1), ncl, cc(1, 1), ldc, &
                    0.0_dp, work4(1, 1, 1), mal*mbl)
         !
         qabc(1:mal, 1:mbl, 1:mcl) = work4(1:mal, 1:mbl, 1:mcl)
         !
         DEALLOCATE (work1, work2, work3, work4)
         !
      ELSE IF (PRESENT(ca) .AND. PRESENT(cb)) THEN
         CPABORT("Not implemented")
      ELSE IF (PRESENT(ca) .AND. PRESENT(cc)) THEN
         CPABORT("Not implemented")
      ELSE IF (PRESENT(cb) .AND. PRESENT(cc)) THEN
         CPABORT("Not implemented")
      ELSE IF (PRESENT(ca)) THEN
         CPABORT("Not implemented")
      ELSE IF (PRESENT(cb)) THEN
         CPABORT("Not implemented")
      ELSE IF (PRESENT(cc)) THEN
         CPABORT("Not implemented")
      ELSE
         ! Copy of arrays is not covered here
         CPABORT("Copy of arrays is not covered here")
      END IF

   END SUBROUTINE contraction_abc

! **************************************************************************************************
!> \brief Applying the de-contraction coefficients to a matrix
!>        QAB <- CA * SAB * CB(T)
!>        Variable "trans" requests the input matrix to be SAB(T)
!>        Active dimensions are: QAB(na,nb), SAB(ma,mb)
!> \param sab     Input matrix, dimension(:,:)
!> \param qab     Output matrix, dimension(:,:)
!> \param ca      Left transformation matrix
!> \param na      First dimension of ca
!> \param ma      Second dimension of ca
!> \param cb      Right transformation matrix
!> \param nb      First dimension of cb
!> \param mb      Second dimension of cb
!> \param trans   Optional transposition of input matrix
! **************************************************************************************************
   SUBROUTINE decontraction_ab(sab, qab, ca, na, ma, cb, nb, mb, trans)

      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: sab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: qab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: ca
      INTEGER, INTENT(IN)                                :: na, ma
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: cb
      INTEGER, INTENT(IN)                                :: nb, mb
      LOGICAL, INTENT(IN), OPTIONAL                      :: trans

      INTEGER                                            :: lda, ldb, ldq, lds, ldw
      LOGICAL                                            :: my_trans
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: work

      ! Should input matrix be transposed?
      IF (PRESENT(trans)) THEN
         my_trans = trans
      ELSE
         my_trans = .FALSE.
      END IF

      lds = SIZE(sab, 1)
      ldq = SIZE(qab, 1)
      lda = SIZE(ca, 1)
      ldb = SIZE(cb, 1)

      ALLOCATE (work(na, mb))
      ldw = na

      IF (my_trans) THEN
!dg      CALL dgemm("N", "T", na, mb, ma, 1.0_dp, ca, lda, sab, lds, 0.0_dp, work, ldw)
         work(1:na, 1:mb) = MATMUL(ca(1:na, 1:ma), TRANSPOSE(sab(1:mb, 1:ma)))
      ELSE
!dg      CALL dgemm("N", "N", na, mb, ma, 1.0_dp, ca, lda, sab, lds, 0.0_dp, work, ldw)
         work(1:na, 1:mb) = MATMUL(ca(1:na, 1:ma), sab(1:ma, 1:mb))
      END IF
!dg   CALL dgemm("N", "T", na, nb, mb, 1.0_dp, work, ldw, cb, ldb, 0.0_dp, qab, ldq)
      qab(1:na, 1:nb) = MATMUL(work(1:na, 1:mb), TRANSPOSE(cb(1:nb, 1:mb)))

      DEALLOCATE (work)

   END SUBROUTINE decontraction_ab

! **************************************************************************************************
!> \brief Routine to trace a series of matrices with another matrix
!>        Calculate forces of type f(:) = Trace(Pab*Sab(:))
!> \param force   Vector to hold output forces
!> \param sab     Input vector of matrices, dimension (:,:,:)
!> \param pab     Input matrix
!> \param na      Active first dimension
!> \param nb      Active second dimension
!> \param m       Number of matrices to be traced
!> \param trans   Matrices are transposed (Sab and Pab)
! **************************************************************************************************
   SUBROUTINE force_trace_ab(force, sab, pab, na, nb, m, trans)

      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: force
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: sab
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: pab
      INTEGER, INTENT(IN)                                :: na, nb, m
      LOGICAL, INTENT(IN), OPTIONAL                      :: trans

      INTEGER                                            :: i
      LOGICAL                                            :: my_trans

      CPASSERT(m <= SIZE(SAB, 3))
      CPASSERT(m <= SIZE(force, 1))

      ! are matrices transposed?
      IF (PRESENT(trans)) THEN
         my_trans = trans
      ELSE
         my_trans = .FALSE.
      END IF

      DO i = 1, m
         IF (my_trans) THEN
            force(i) = SUM(sab(1:nb, 1:na, i)*pab(1:nb, 1:na))
         ELSE
            force(i) = SUM(sab(1:na, 1:nb, i)*pab(1:na, 1:nb))
         END IF
      END DO

   END SUBROUTINE force_trace_ab

! **************************************************************************************************
!> \brief Copy a block out of a matrix and add it to another matrix
!>        SAB = SAB + QAB  or  QAB = QAB + SAB
!>        QAB(ia:,ib:) and SAB(1:,1:)
!> \param dir    "IN" and "OUT" defines direction of copy
!> \param sab    Matrix input for "IN", output for "OUT"
!> \param na     first dimension of matrix to copy
!> \param nb     second dimension of matrix to copy
!> \param qab    Matrix output for "IN", input for "OUT"
!>               Use subblock of this matrix
!> \param ia     Starting index in qab first dimension
!> \param ib     Starting index in qab second dimension
!> \param trans  Matrices (qab and sab) are transposed
! **************************************************************************************************
   SUBROUTINE block_add_ab(dir, sab, na, nb, qab, ia, ib, trans)

      CHARACTER(LEN=*), INTENT(IN)                       :: dir
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: sab
      INTEGER, INTENT(IN)                                :: na, nb
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: qab
      INTEGER, INTENT(IN)                                :: ia, ib
      LOGICAL, INTENT(IN), OPTIONAL                      :: trans

      INTEGER                                            :: ja, jb
      LOGICAL                                            :: my_trans

      IF (PRESENT(trans)) THEN
         my_trans = trans
      ELSE
         my_trans = .FALSE.
      END IF

      IF (dir == "IN" .OR. dir == "in") THEN
         !  QAB(block) <= SAB
         ja = ia + na - 1
         jb = ib + nb - 1
         IF (my_trans) THEN
            qab(ib:jb, ia:ja) = qab(ib:jb, ia:ja) + sab(1:nb, 1:na)
         ELSE
            qab(ia:ja, ib:jb) = qab(ia:ja, ib:jb) + sab(1:na, 1:nb)
         END IF
      ELSEIF (dir == "OUT" .OR. dir == "out") THEN
         !  SAB <= QAB(block)
         ja = ia + na - 1
         jb = ib + nb - 1
         IF (my_trans) THEN
            sab(1:nb, 1:na) = sab(1:nb, 1:na) + qab(ib:jb, ia:ja)
         ELSE
            sab(1:na, 1:nb) = sab(1:na, 1:nb) + qab(ia:ja, ib:jb)
         END IF
      ELSE
         CPABORT("")
      END IF

   END SUBROUTINE block_add_ab
! **************************************************************************************************

END MODULE ai_contraction
