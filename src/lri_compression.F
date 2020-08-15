!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief integral compression (fix point accuracy)
!> \par History
!>      created JGH [11.2017]
!> \authors JGH
! **************************************************************************************************
MODULE lri_compression
   USE kinds,                           ONLY: dp,&
                                              sp
   USE lri_environment_types,           ONLY: carray,&
                                              int_container
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! **************************************************************************************************

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'lri_compression'

   PUBLIC :: lri_comp, lri_decomp_i, lri_cont_mem

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param aval ...
!> \param amax ...
!> \param cont ...
! **************************************************************************************************
   SUBROUTINE lri_comp(aval, amax, cont)
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(IN)      :: aval
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: amax
      TYPE(int_container), INTENT(INOUT)                 :: cont

      INTEGER                                            :: i, ia, ib, ii, na, nb, nc, nn
      REAL(KIND=dp)                                      :: xm
      TYPE(carray), POINTER                              :: ca

      IF (ASSOCIATED(cont%ca)) THEN
         DO i = 1, SIZE(cont%ca)
            IF (ASSOCIATED(cont%ca(i)%cdp)) DEALLOCATE (cont%ca(i)%cdp)
            IF (ASSOCIATED(cont%ca(i)%csp)) DEALLOCATE (cont%ca(i)%csp)
            IF (ASSOCIATED(cont%ca(i)%cip)) DEALLOCATE (cont%ca(i)%cip)
         END DO
      END IF

      na = SIZE(aval, 1)
      nb = SIZE(aval, 2)
      nc = SIZE(aval, 3)
      nn = na*nb
      cont%na = na
      cont%nb = nb
      cont%nc = nc

      IF (nc > 0) THEN
         ALLOCATE (cont%ca(nc))
         DO i = 1, nc
            ca => cont%ca(i)
            NULLIFY (ca%cdp, ca%csp, ca%cip)
            xm = MAXVAL(ABS(aval(:, :, i)))
            IF (xm >= 1.0e-05_dp) THEN
               ca%compression = 1
               ALLOCATE (ca%cdp(nn))
               ii = 0
               DO ib = 1, nb
                  DO ia = 1, na
                     ii = ii + 1
                     ca%cdp(ii) = aval(ia, ib, i)
                  END DO
               END DO
            ELSE IF (xm >= 1.0e-10_dp) THEN
               ca%compression = 2
               ALLOCATE (ca%csp(nn))
               ii = 0
               DO ib = 1, nb
                  DO ia = 1, na
                     ii = ii + 1
                     ca%csp(ii) = REAL(aval(ia, ib, i), KIND=sp)
                  END DO
               END DO
            ELSE
               ca%compression = 0
            END IF
            amax(i) = xm
         END DO
      END IF

   END SUBROUTINE lri_comp

! **************************************************************************************************
!> \brief ...
!> \param cont ...
!> \return ...
! **************************************************************************************************
   FUNCTION lri_cont_mem(cont) RESULT(cmem)
      TYPE(int_container), INTENT(IN)                    :: cont
      REAL(KIND=dp)                                      :: cmem

      INTEGER                                            :: i

      cmem = 0.0_dp
      IF (ASSOCIATED(cont%ca)) THEN
         DO i = 1, SIZE(cont%ca)
            IF (ASSOCIATED(cont%ca(i)%cdp)) THEN
               cmem = cmem + SIZE(cont%ca(i)%cdp)
            END IF
            IF (ASSOCIATED(cont%ca(i)%csp)) THEN
               cmem = cmem + 0.5_dp*SIZE(cont%ca(i)%csp)
            END IF
            IF (ASSOCIATED(cont%ca(i)%cip)) THEN
               cmem = cmem + SIZE(cont%ca(i)%cip)
            END IF
         END DO
      END IF

   END FUNCTION lri_cont_mem
! **************************************************************************************************
!> \brief ...
!> \param aval ...
!> \param cont ...
!> \param ival ...
! **************************************************************************************************
   SUBROUTINE lri_decomp_i(aval, cont, ival)
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: aval
      TYPE(int_container), INTENT(INOUT)                 :: cont
      INTEGER                                            :: ival

      INTEGER                                            :: ia, ib, ii, na, nb, nn
      TYPE(carray), POINTER                              :: ca

      na = SIZE(aval, 1)
      nb = SIZE(aval, 2)
      nn = na*nb
      CPASSERT(na == cont%na)
      CPASSERT(nb == cont%nb)
      CPASSERT(ival <= cont%nc)

      ca => cont%ca(ival)
      !
      SELECT CASE (ca%compression)
      CASE (0)
         aval(1:na, 1:nb) = 0.0_dp
      CASE (1)
         ii = 0
         DO ib = 1, nb
            DO ia = 1, na
               ii = ii + 1
               aval(ia, ib) = ca%cdp(ii)
            END DO
         END DO
      CASE (2)
         ii = 0
         DO ib = 1, nb
            DO ia = 1, na
               ii = ii + 1
               aval(ia, ib) = REAL(ca%csp(ii), KIND=dp)
            END DO
         END DO
      CASE DEFAULT
         CPABORT("lri_decomp_i: compression label invalid")
      END SELECT

   END SUBROUTINE lri_decomp_i

END MODULE lri_compression

