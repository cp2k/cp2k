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

  USE ai_onecenter,                    ONLY: &
       sg_conf, sg_coulomb, sg_erf, sg_exchange, sg_gpot, sg_kinetic, &
       sg_nuclear, sg_overlap, sg_proj_ol
  USE atom_types,                      ONLY: &
       CGTO_BASIS, GTH_PSEUDO, GTO_BASIS, NO_PSEUDO, NUM_BASIS, STO_BASIS, &
       atom_basis_type, atom_integrals, atom_potential_type
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_operators'

  PUBLIC :: atom_int_setup, atom_ppint_setup, atom_int_release, atom_ppint_release

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE atom_int_setup(integrals,basis,potential,&
                            eri_coulomb,eri_exchange,all_nu,error)
    TYPE(atom_integrals), INTENT(INOUT)      :: integrals
    TYPE(atom_basis_type), INTENT(INOUT)     :: basis
    TYPE(atom_potential_type), INTENT(IN), &
      OPTIONAL                               :: potential
    LOGICAL, INTENT(IN), OPTIONAL            :: eri_coulomb, eri_exchange, &
                                                all_nu
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_int_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, l, l1, l2, ll, m, m1, &
                                                m2, mm1, mm2, n, n1, n2, nn1, &
                                                nn2, nu, nx
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: omat
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

      NULLIFY ( integrals%ovlp, integrals%kin, integrals%core, integrals%conf )
      DO ll=1,SIZE(integrals%ceri)
        NULLIFY ( integrals%ceri(ll)%int, integrals%eeri(ll)%int )
      END DO
    
      ALLOCATE (integrals%ovlp(n,n,0:3),STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      integrals%ovlp = 0._dp

      ALLOCATE (integrals%kin(n,n,0:3),STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      integrals%kin = 0._dp

      integrals%status = 1

      IF ( PRESENT(potential) ) THEN
        IF ( potential%confinement ) THEN
          ALLOCATE (integrals%conf(n,n,0:3),STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
        END IF
      END IF

      SELECT CASE (basis%basis_type)
        CASE DEFAULT
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        CASE (GTO_BASIS)
          DO l=0,3
            n = integrals%n(l)
            CALL  sg_overlap ( integrals%ovlp(1:n,1:n,l), l, basis%am(1:n,l), basis%am(1:n,l) )
            CALL  sg_kinetic ( integrals%kin(1:n,1:n,l), l, basis%am(1:n,l), basis%am(1:n,l) )
            IF ( PRESENT(potential) ) THEN
              IF ( potential%confinement ) THEN
                CALL  sg_conf ( integrals%conf(1:n,1:n,l), potential%rcon, potential%ncon/2, &
                                l, basis%am(1:n,l), basis%am(1:n,l) )
              END IF
            END IF
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
            CALL contract2(integrals%ovlp(1:n,1:n,l),omat(1:m,1:m),basis%cm(1:m,1:n,l), error)
            CALL  sg_kinetic ( omat(1:m,1:m), l, basis%am(1:m,l), basis%am(1:m,l) )
            CALL contract2(integrals%kin(1:n,1:n,l),omat(1:m,1:m),basis%cm(1:m,1:n,l), error)
            IF ( PRESENT(potential) ) THEN
              IF ( potential%confinement ) THEN
                CALL  sg_conf ( omat(1:m,1:m), potential%rcon, potential%ncon/2, &
                                l, basis%am(1:m,l), basis%am(1:m,l) )
                CALL contract2(integrals%conf(1:n,1:n,l),omat(1:m,1:m),basis%cm(1:m,1:n,l), error)
              END IF
            END IF
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
  SUBROUTINE atom_ppint_setup(integrals,basis,potential,error)
    TYPE(atom_integrals), INTENT(INOUT)      :: integrals
    TYPE(atom_basis_type), INTENT(INOUT)     :: basis
    TYPE(atom_potential_type), INTENT(IN)    :: potential
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_ppint_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, k, l, m, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: alpha
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: xmat
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: omat, spmat

    failure = .FALSE.

    IF ( integrals%ppstat == 0 ) THEN
      n = MAXVAL(basis%nbas)
      integrals%n = basis%nbas

      NULLIFY ( integrals%core, integrals%hnl )
    
      ALLOCATE (integrals%hnl(n,n,0:3),STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      integrals%hnl = 0._dp

      ALLOCATE (integrals%core(n,n,0:3),STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      integrals%core = 0._dp

      integrals%ppstat = 1

      SELECT CASE (basis%basis_type)
        CASE DEFAULT
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        CASE (GTO_BASIS)

          SELECT CASE (potential%ppot_type)
            CASE (NO_PSEUDO)
              DO l=0,3
                n = integrals%n(l)
                CALL sg_nuclear ( integrals%core(1:n,1:n,l), l, basis%am(1:n,l), basis%am(1:n,l) )
              END DO
            CASE (GTH_PSEUDO)
              alpha = 1._dp/potential%gth_pot%rc/SQRT(2._dp)
              DO l=0,3
                n = integrals%n(l)
                ALLOCATE (omat(n,n),spmat(n,5),STAT=ierr)
                CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

                omat = 0._dp
                CALL sg_erf ( omat(1:n,1:n), l, alpha, basis%am(1:n,l), basis%am(1:n,l) )
                integrals%core(1:n,1:n,l) = -potential%gth_pot%zion*omat(1:n,1:n)
                DO i=1,potential%gth_pot%ncl
                  omat = 0._dp
                  CALL sg_gpot ( omat(1:n,1:n), i-1, potential%gth_pot%rc, l, basis%am(1:n,l), basis%am(1:n,l) )
                  integrals%core(1:n,1:n,l) = integrals%core(1:n,1:n,l) + &
                     potential%gth_pot%cl(i)*omat(1:n,1:n)
                END DO

                spmat = 0._dp
                m = potential%gth_pot%nl(l)
                DO i=1,m
                  CALL sg_proj_ol ( spmat(1:n,i), l, basis%am(1:n,l), i-1, potential%gth_pot%rcnl(l) )
                END DO
                integrals%hnl(1:n,1:n,l) = MATMUL(spmat(1:n,1:m),&
                    MATMUL(potential%gth_pot%hnl(1:m,1:m,l),TRANSPOSE(spmat(1:n,1:m))))

                DEALLOCATE (omat,spmat,STAT=ierr)
                CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
              END DO
            CASE DEFAULT
              CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT

        CASE (CGTO_BASIS)

          SELECT CASE (potential%ppot_type)
            CASE (NO_PSEUDO)
              DO l=0,3
                n = integrals%n(l)
                m = basis%nprim(l)
                ALLOCATE (omat(m,m),STAT=ierr)
                CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

                CALL sg_nuclear ( omat(1:m,1:m), l, basis%am(1:m,l), basis%am(1:m,l) )
                CALL contract2(integrals%core(1:n,1:n,l),omat(1:m,1:m),basis%cm(1:m,1:n,l), error)

                DEALLOCATE (omat,STAT=ierr)
                CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
              END DO
            CASE (GTH_PSEUDO)
              alpha = 1._dp/potential%gth_pot%rc/SQRT(2._dp)
              DO l=0,3
                n = integrals%n(l)
                m = basis%nprim(l)
                ALLOCATE (omat(m,m),spmat(n,n),xmat(m),STAT=ierr)
                CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

                omat = 0._dp
                CALL sg_erf ( omat(1:m,1:m), l, alpha, basis%am(1:m,l), basis%am(1:m,l) )
                omat(1:m,1:m) = -potential%gth_pot%zion*omat(1:m,1:m)
                CALL contract2(integrals%core(1:n,1:n,l),omat(1:m,1:m),basis%cm(1:m,1:n,l), error)
                DO i=1,potential%gth_pot%ncl
                  omat = 0._dp
                  CALL sg_gpot ( omat(1:m,1:m), i-1, potential%gth_pot%rc, l, basis%am(1:m,l), basis%am(1:m,l) )
                  omat(1:m,1:m) = potential%gth_pot%cl(i)*omat(1:m,1:m)
                  CALL contract2add(integrals%core(1:n,1:n,l),omat(1:m,1:m),basis%cm(1:m,1:n,l), error)
                END DO

                spmat = 0._dp
                k = potential%gth_pot%nl(l)
                DO i=1,k
                  CALL sg_proj_ol ( xmat(1:m), l, basis%am(1:m,l), i-1, potential%gth_pot%rcnl(l) )
                  spmat(1:n,i) = MATMUL(TRANSPOSE(basis%cm(1:m,1:n,l)),xmat(1:m))
                END DO
                integrals%hnl(1:n,1:n,l) = MATMUL(spmat(1:n,1:k),&
                    MATMUL(potential%gth_pot%hnl(1:k,1:k,l),TRANSPOSE(spmat(1:n,1:k))))

                DEALLOCATE (omat,spmat,xmat,STAT=ierr)
                CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
              END DO
            CASE DEFAULT
              CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT

        CASE (STO_BASIS)
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        CASE (NUM_BASIS)
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT

    END IF

  END SUBROUTINE atom_ppint_setup
! *****************************************************************************
  SUBROUTINE atom_int_release(integrals,error)
    TYPE(atom_integrals), INTENT(INOUT)      :: integrals
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_int_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, ll
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF ( ASSOCIATED(integrals%ovlp) ) THEN
      DEALLOCATE (integrals%ovlp,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF
    IF ( ASSOCIATED(integrals%kin) ) THEN
      DEALLOCATE (integrals%kin,STAT=ierr)
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

    integrals%status = 0

  END SUBROUTINE atom_int_release
! *****************************************************************************
  SUBROUTINE atom_ppint_release(integrals,error)
    TYPE(atom_integrals), INTENT(INOUT)      :: integrals
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_ppint_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF ( ASSOCIATED(integrals%hnl) ) THEN
      DEALLOCATE (integrals%hnl,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF
    IF ( ASSOCIATED(integrals%core) ) THEN
      DEALLOCATE (integrals%core,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF

    integrals%ppstat = 0

  END SUBROUTINE atom_ppint_release
! *****************************************************************************
  SUBROUTINE contract2 ( int, omat, cm, error )
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: int
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: omat, cm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'contract2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: m, n

    n = SIZE(int,1)
    m = SIZE(omat,1)

    INT(1:n,1:n) = MATMUL(TRANSPOSE(cm(1:m,1:n)),MATMUL(omat(1:m,1:m),cm(1:m,1:n)))

  END SUBROUTINE contract2 
! *****************************************************************************
  SUBROUTINE contract2add ( int, omat, cm, error )
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: int
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: omat, cm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'contract2add', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: m, n

    n = SIZE(int,1)
    m = SIZE(omat,1)

    INT(1:n,1:n) = INT(1:n,1:n) + MATMUL(TRANSPOSE(cm(1:m,1:n)),MATMUL(omat(1:m,1:m),cm(1:m,1:n)))

  END SUBROUTINE contract2add
! *****************************************************************************
  SUBROUTINE contract4 ( eri, omat, cm1, cm2, error )
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: eri
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: omat, cm1, cm2
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'contract4', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i2, ierr, m1, m2, mm1, &
                                                mm2, n1, n2, nn1, nn2
    LOGICAL                                  :: failure = .FALSE.
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: amat, atran, bmat, btran, hint

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
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: mat
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: vec
    INTEGER, INTENT(IN)                      :: n

    INTEGER                                  :: i, ij, j

    ij = 0
    DO i=1,n
      DO j=i,n
        ij = ij + 1
        vec(ij) = mat(i,j)
      END DO
    END DO

  END SUBROUTINE ipack

  SUBROUTINE iunpack(mat,vec,n)
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: mat
    REAL(dp), DIMENSION(:), INTENT(IN)       :: vec
    INTEGER, INTENT(IN)                      :: n

    INTEGER                                  :: i, ij, j

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
