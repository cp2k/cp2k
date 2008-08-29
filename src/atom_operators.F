!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Calculate the atomic operator matrices
!> \author  jgh
!> \date    03.03.2008
!> \version 1.0
!>
! *****************************************************************************
MODULE atom_operators

  USE ai_onecenter,                    ONLY: sg_coulomb,&
                                             sg_exchange,&
                                             sg_kinetic,&
                                             sg_nuclear,&
                                             sg_overlap
  USE atom_types,                      ONLY: CGTO_BASIS,&
                                             GTO_BASIS,&
                                             NUM_BASIS,&
                                             STO_BASIS,&
                                             atom_basis_type,&
                                             atom_integrals
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_operators'

  PUBLIC :: atom_int_setup, atom_int_release

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE atom_int_setup(integrals,basis,eri_coulomb,eri_exchange,all_nu,error)
    TYPE(atom_integrals), INTENT(INOUT)      :: integrals
    TYPE(atom_basis_type), INTENT(INOUT)     :: basis
    LOGICAL, INTENT(IN), OPTIONAL            :: eri_coulomb, eri_exchange, &
                                                all_nu
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_int_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, l, l1, l2, ll, n, n1, &
       n2, nn1, nn2, nu, nx, m, m1, m2, mm1, mm2
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE :: omat
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: eri

    failure = .FALSE.

    IF ( integrals%status == 0 ) THEN
      n = MAXVAL(basis%nbas)
      integrals%n = basis%nbas

      IF ( PRESENT(eri_coulomb) ) THEN
        integrals%eri_coulomb = eri_coulomb
      ELSE
        integrals%eri_coulomb = .FALSE.
      END IF
      IF ( PRESENT(eri_exchange) ) THEN
        integrals%eri_exchange = eri_exchange
      ELSE
        integrals%eri_exchange = .FALSE.
      END IF
      IF ( PRESENT(all_nu) ) THEN
        integrals%all_nu = all_nu
      ELSE
        integrals%all_nu = .FALSE.
      END IF

      NULLIFY ( integrals%ovlp, integrals%kin, integrals%core, &
                integrals%proj, integrals%conf )
      DO ll=1,SIZE(integrals%ceri)
        NULLIFY ( integrals%ceri(ll)%int, integrals%eeri(ll)%int )
      END DO
    
      ALLOCATE (integrals%ovlp(n,n,0:3),STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      integrals%ovlp = 0._dp

      ALLOCATE (integrals%kin(n,n,0:3),STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      integrals%kin = 0._dp

      ALLOCATE (integrals%core(n,n,0:3),STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      integrals%core = 0._dp
      integrals%status = 1

      SELECT CASE (basis%basis_type)
        CASE DEFAULT
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        CASE (GTO_BASIS)
          DO l=0,3
            n = integrals%n(l)
            CALL  sg_overlap ( integrals%ovlp(1:n,1:n,l), l, basis%am(1:n,l), basis%am(1:n,l) )
            CALL  sg_kinetic ( integrals%kin(1:n,1:n,l), l, basis%am(1:n,l), basis%am(1:n,l) )
            CALL  sg_nuclear ( integrals%core(1:n,1:n,l), l, basis%am(1:n,l), basis%am(1:n,l) )
          END DO
          IF ( integrals%eri_coulomb ) THEN
            ll = 0
            DO l1=0,3
              n1 = integrals%n(l1)
              nn1 = (n1*(n1+1))/2
              DO l2=0,l1
                n2 = integrals%n(l2)
                nn2 = (n2*(n2+1))/2
                IF ( integrals%all_nu ) THEN
                  nx = MIN(2*l1,2*l2)
                ELSE
                  nx = 0
                END IF
                DO nu = 0, nx, 2
                  ll = ll + 1
                  ALLOCATE (integrals%ceri(ll)%int(nn1,nn2),STAT=ierr)
                  CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
                  integrals%ceri(ll)%int = 0._dp
                  eri => integrals%ceri(ll)%int
                  CALL sg_coulomb ( eri, nu, basis%am(1:n1,l1), l1, basis%am(1:n2,l2), l2 )
                END DO
              END DO
            END DO
          END IF
          IF ( integrals%eri_exchange ) THEN
            ll = 0
            DO l1=0,3
              n1 = integrals%n(l1)
              nn1 = (n1*(n1+1))/2
              DO l2=0,l1
                n2 = integrals%n(l2)
                nn2 = (n2*(n2+1))/2
                DO nu = ABS(l1-l2),l1+l2,2
                  ll = ll + 1
                  ALLOCATE (integrals%eeri(ll)%int(nn1,nn2),STAT=ierr)
                  CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
                  integrals%eeri(ll)%int = 0._dp
                  eri => integrals%eeri(ll)%int
                  CALL sg_exchange ( eri, nu, basis%am(1:n1,l1), l1, basis%am(1:n2,l2), l2 )
                END DO
              END DO
            END DO
          END IF
        CASE (CGTO_BASIS)
          DO l=0,3
            n = integrals%n(l)
            m = basis%nprim(l)
            ALLOCATE (omat(m,m),STAT=ierr)
            CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

            CALL  sg_overlap ( omat(1:m,1:m), l, basis%am(1:m,l), basis%am(1:m,l) )
            CALL contract2(integrals%ovlp(1:n,1:n,l),omat,basis%cm(1:m,1:n,l), error)
            CALL  sg_kinetic ( omat(1:m,1:m), l, basis%am(1:m,l), basis%am(1:m,l) )
            CALL contract2(integrals%kin(1:n,1:n,l),omat,basis%cm(1:m,1:n,l), error)
            CALL  sg_nuclear ( omat(1:m,1:m), l, basis%am(1:m,l), basis%am(1:m,l) )
            CALL contract2(integrals%core(1:n,1:n,l),omat,basis%cm(1:m,1:n,l), error)

            DEALLOCATE (omat,STAT=ierr)
            CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          END DO
          IF ( integrals%eri_coulomb ) THEN
            ll = 0
            DO l1=0,3
              n1 = integrals%n(l1)
              nn1 = (n1*(n1+1))/2
              m1 = basis%nprim(l1)
              mm1 = (m1*(m1+1))/2
              DO l2=0,l1
                n2 = integrals%n(l2)
                nn2 = (n2*(n2+1))/2
                m2 = basis%nprim(l2)
                mm2 = (m2*(m2+1))/2
                IF ( integrals%all_nu ) THEN
                  nx = MIN(2*l1,2*l2)
                ELSE
                  nx = 0
                END IF
                DO nu = 0, nx, 2
                  ll = ll + 1
                  ALLOCATE (integrals%ceri(ll)%int(nn1,nn2),STAT=ierr)
                  CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
                  integrals%ceri(ll)%int = 0._dp
                  ALLOCATE (omat(mm1,mm2),STAT=ierr)
                  CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

                  eri => integrals%ceri(ll)%int
                  CALL sg_coulomb ( omat, nu, basis%am(1:m1,l1), l1, basis%am(1:m2,l2), l2 )
                  CALL contract4 ( eri, omat, basis%cm(1:m1,1:n1,l1), basis%cm(1:m2,1:n2,l2), error )

                  DEALLOCATE (omat,STAT=ierr)
                  CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
                END DO
              END DO
            END DO
          END IF
          IF ( integrals%eri_exchange ) THEN
            ll = 0
            DO l1=0,3
              n1 = integrals%n(l1)
              nn1 = (n1*(n1+1))/2
              m1 = basis%nprim(l1)
              mm1 = (m1*(m1+1))/2
              DO l2=0,l1
                n2 = integrals%n(l2)
                nn2 = (n2*(n2+1))/2
                m2 = basis%nprim(l2)
                mm2 = (m2*(m2+1))/2
                DO nu = ABS(l1-l2),l1+l2,2
                  ll = ll + 1
                  ALLOCATE (integrals%eeri(ll)%int(nn1,nn2),STAT=ierr)
                  CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
                  integrals%eeri(ll)%int = 0._dp
                  ALLOCATE (omat(mm1,mm2),STAT=ierr)
                  CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

                  eri => integrals%eeri(ll)%int
                  CALL sg_exchange ( omat, nu, basis%am(1:m1,l1), l1, basis%am(1:m2,l2), l2 )
                  CALL contract4 ( eri, omat, basis%cm(1:m1,1:n1,l1), basis%cm(1:m2,1:n2,l2), error )

                  DEALLOCATE (omat,STAT=ierr)
                  CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
                END DO
              END DO
            END DO
          END IF
        CASE (STO_BASIS)
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        CASE (NUM_BASIS)
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT

    END IF

  END SUBROUTINE atom_int_setup
! *****************************************************************************
  SUBROUTINE atom_int_release(integrals,error)
    TYPE(atom_integrals), INTENT(INOUT)      :: integrals
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_int_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, ll, n
    LOGICAL                                  :: failure

    failure = .FALSE.

    integrals%n = 0

    IF ( ASSOCIATED(integrals%ovlp) ) THEN
      DEALLOCATE (integrals%ovlp,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF
    IF ( ASSOCIATED(integrals%kin) ) THEN
      DEALLOCATE (integrals%kin,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF
    IF ( ASSOCIATED(integrals%core) ) THEN
      DEALLOCATE (integrals%core,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF
    IF ( ASSOCIATED(integrals%proj) ) THEN
      DEALLOCATE (integrals%proj,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF
    IF ( ASSOCIATED(integrals%conf) ) THEN
      DEALLOCATE (integrals%conf,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF
    DO ll=1,SIZE(integrals%ceri)
      IF ( ASSOCIATED(integrals%ceri(ll)%int) ) THEN
        DEALLOCATE (integrals%ceri(ll)%int,STAT=ierr)
        CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      END IF
      IF ( ASSOCIATED(integrals%eeri(ll)%int) ) THEN
        DEALLOCATE (integrals%eeri(ll)%int,STAT=ierr)
        CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      END IF
    END DO

  END SUBROUTINE atom_int_release
! *****************************************************************************
  SUBROUTINE contract2 ( int, omat, cm, error )
    REAL(dp), DIMENSION(:,:), INTENT(INOUT)  :: int
    REAL(dp), DIMENSION(:,:), INTENT(IN)     :: omat, cm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'contract2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: n, m

    n = SIZE(int,1)
    m = SIZE(omat,1)

    int(1:n,1:n) = MATMUL(TRANSPOSE(cm(1:m,1:n)),MATMUL(omat(1:m,1:m),cm(1:m,1:n)))

  END SUBROUTINE contract2 
! *****************************************************************************
  SUBROUTINE contract4 ( eri, omat, cm1, cm2, error )
    REAL(dp), DIMENSION(:,:), INTENT(INOUT)  :: eri
    REAL(dp), DIMENSION(:,:), INTENT(IN)     :: omat, cm1, cm2
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'contract4', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure=.FALSE.
    INTEGER                                  :: n1, n2, m1, m2, nn1, nn2, mm1, mm2
    INTEGER                                  :: i1, i2, ierr
    REAL(dp), DIMENSION(:,:), ALLOCATABLE    :: amat, atran, bmat, btran, hint

    m1 = SIZE(cm1,1)
    n1 = SIZE(cm1,2)
    m2 = SIZE(cm2,1)
    n2 = SIZE(cm2,2)
    nn1 = SIZE(eri,1)
    nn2 = SIZE(eri,2)
    mm1 = SIZE(omat,1)
    mm2 = SIZE(omat,2)

    ALLOCATE(amat(m1,m1),atran(n1,n1),bmat(m2,m2),btran(n2,n2),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    ALLOCATE(hint(mm1,nn2),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    DO i1=1,mm1
      CALL iunpack(bmat(1:m2,1:m2),omat(i1,1:mm2),m2)
      CALL contract2( btran(1:n2,1:n2), bmat(1:m2,1:m2), cm2(1:m2,1:n2), error )
      CALL ipack(btran(1:n2,1:n2),hint(i1,1:nn2),n2)
    END DO

    DO i2=1,nn2
      CALL iunpack(amat(1:m1,1:m1),hint(1:mm1,i2),m1)
      CALL contract2( atran(1:n1,1:n1), amat(1:m1,1:m1), cm1(1:m1,1:n1), error )
      CALL ipack(atran(1:n1,1:n1),eri(1:nn1,i2),n1)
    END DO

    DEALLOCATE(amat,atran,bmat,btran,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    DEALLOCATE(hint,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE contract4
! *****************************************************************************
  SUBROUTINE ipack(mat,vec,n)
    REAL(dp), DIMENSION(:,:), INTENT(IN)     :: mat
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: vec
    INTEGER, INTENT(IN)                      :: n

    INTEGER                                  :: i, j, ij

    ij = 0
    DO i=1,n
      DO j=i,n
        ij = ij + 1
        vec(ij) = mat(i,j)
      END DO
    END DO

  END SUBROUTINE ipack

  SUBROUTINE iunpack(mat,vec,n)
    REAL(dp), DIMENSION(:,:), INTENT(INOUT)  :: mat
    REAL(dp), DIMENSION(:), INTENT(IN)       :: vec
    INTEGER, INTENT(IN)                      :: n

    INTEGER                                  :: i, j, ij

    ij = 0
    DO i=1,n
      DO j=i,n
        ij = ij + 1
        mat(i,j) = vec(ij)
        mat(j,i) = vec(ij)
      END DO
    END DO

  END SUBROUTINE iunpack
! *****************************************************************************

END MODULE atom_operators
