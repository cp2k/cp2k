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
                                                n2, nn1, nn2, nu, nx
    LOGICAL                                  :: failure
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
            CALL  sg_overlap ( integrals%ovlp(1:n,1:n,l), l, basis%al, basis%al )
            CALL  sg_kinetic ( integrals%kin(1:n,1:n,l), l, basis%al, basis%al )
            CALL  sg_nuclear ( integrals%core(1:n,1:n,l), l, basis%al, basis%al )
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
                  CALL sg_coulomb ( eri, nu, basis%al, l1, basis%al, l2 )
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
                  CALL sg_exchange ( eri, nu, basis%al, l1, basis%al, l2 )
                END DO
              END DO
            END DO
          END IF
        CASE (CGTO_BASIS)
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
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

END MODULE atom_operators
