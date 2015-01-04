!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
       sg_coulomb, sg_erf, sg_erfc, sg_exchange, sg_gpot, sg_kinetic, &
       sg_kinnuc, sg_nuclear, sg_overlap, sg_proj_ol, sto_kinetic, &
       sto_nuclear, sto_overlap
  USE atom_types,                      ONLY: &
       CGTO_BASIS, GTH_PSEUDO, GTO_BASIS, NO_PSEUDO, NUM_BASIS, STO_BASIS, &
       atom_basis_type, atom_integrals, atom_potential_type, atom_state
  USE atom_utils,                      ONLY: atom_solve,&
                                             coulomb_potential_numeric,&
                                             integrate_grid,&
                                             numpot_matrix,&
                                             slater_density,&
                                             wigner_slater_functional
  USE dkh_main,                        ONLY: dkh_atom_transformation
  USE erf_fn,                          ONLY: erfc
  USE input_constants,                 ONLY: &
       do_dkh0_atom, do_dkh1_atom, do_dkh2_atom, do_dkh3_atom, do_dkh4_atom, &
       do_dkh5_atom, do_nonrel_atom, do_sczoramp_atom, do_zoramp_atom
  USE kinds,                           ONLY: dp
  USE lapack,                          ONLY: lapack_sgesv
  USE mathconstants,                   ONLY: gamma1,&
                                             sqrt2
  USE mathlib,                         ONLY: jacobi
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: a_fine
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_operators'

  PUBLIC :: atom_int_setup, atom_ppint_setup, atom_int_release, atom_ppint_release
  PUBLIC :: atom_relint_setup, atom_relint_release

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param integrals ...
!> \param basis ...
!> \param potential ...
!> \param eri_coulomb ...
!> \param eri_exchange ...
!> \param all_nu ...
!> \param error ...
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

    INTEGER                                  :: handle, i, ierr, ii, info, &
                                                ipiv(1000), l, l1, l2, ll, &
                                                lwork, m, m1, m2, mm1, mm2, &
                                                n, n1, n2, nn1, nn2, nu, nx
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: rc, sc
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cpot, w, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: omat, vmat
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: eri

    failure = .FALSE.

    CALL timeset(routineN,handle)

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
          integrals%conf = 0._dp
          m = basis%grid%nr
          ALLOCATE (cpot(1:m),STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          rc = potential%rcon
          sc = potential%scon
          cpot(1:m) = (basis%grid%rad(1:m)/rc)**sc
          CALL numpot_matrix(integrals%conf,cpot,basis,0,error)
          DEALLOCATE (cpot,STAT=ierr)
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
            IF (n>0 .AND. m>0) THEN
              ALLOCATE (omat(m,m),STAT=ierr)
              CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

              CALL  sg_overlap ( omat(1:m,1:m), l, basis%am(1:m,l), basis%am(1:m,l) )
              CALL contract2(integrals%ovlp(1:n,1:n,l),omat(1:m,1:m),basis%cm(1:m,1:n,l), error)
              CALL  sg_kinetic ( omat(1:m,1:m), l, basis%am(1:m,l), basis%am(1:m,l) )
              CALL contract2(integrals%kin(1:n,1:n,l),omat(1:m,1:m),basis%cm(1:m,1:n,l), error)
              DEALLOCATE (omat,STAT=ierr)
              CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
            END IF
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
          DO l=0,3
            n = integrals%n(l)
            CALL sto_overlap ( integrals%ovlp(1:n,1:n,l), basis%ns(1:n,l), basis%as(1:n,l), &
                               basis%ns(1:n,l), basis%as(1:n,l) )
            CALL sto_kinetic ( integrals%kin(1:n,1:n,l), l, basis%ns(1:n,l), basis%as(1:n,l), &
                               basis%ns(1:n,l), basis%as(1:n,l) )
          END DO
          CPAssert(.NOT.integrals%eri_coulomb,cp_failure_level,routineP,error,failure)
          CPAssert(.NOT.integrals%eri_exchange,cp_failure_level,routineP,error,failure)
        CASE (NUM_BASIS)
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT

      ! setup transformation matrix to get an orthogonal basis, remove linear dependencies
      NULLIFY(integrals%utrans,integrals%uptrans)
      n = MAXVAL(basis%nbas)
      ALLOCATE (integrals%utrans(n,n,0:3),integrals%uptrans(n,n,0:3),STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      integrals%utrans = 0._dp
      integrals%uptrans = 0._dp
      integrals%nne = integrals%n
      lwork = 10*n
      ALLOCATE(omat(n,n),w(n),vmat(n,n),work(lwork),STAT=ierr)
      CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
      DO l = 0, 3
        n = integrals%n(l)
        IF ( n > 0 ) THEN
          omat(1:n,1:n) = integrals%ovlp(1:n,1:n,l)
          CALL jacobi(omat(1:n,1:n),w(1:n),vmat(1:n,1:n))
          omat(1:n,1:n) = vmat(1:n,1:n)
          ii = 0
          DO i=1,n
            IF (w(i) > basis%eps_eig) THEN
              ii = ii + 1
              integrals%utrans(1:n,ii,l) = omat(1:n,i)/SQRT(w(i))
            END IF
          END DO
          integrals%nne(l) = ii
          omat(1:ii,1:ii)=MATMUL(TRANSPOSE(integrals%utrans(1:n,1:ii,l)),integrals%utrans(1:n,1:ii,l))
          DO i=1,ii
            integrals%uptrans(i,i,l)=1._dp
          ENDDO
          CALL lapack_sgesv ( ii, ii, omat(1:ii,1:ii), ii, ipiv, integrals%uptrans(1:ii,1:ii,l), ii, info )
          CPPostcondition(info==0,cp_failure_level,routineP,error,failure)
        END IF
      END DO
      DEALLOCATE(omat,vmat,w,work,STAT=ierr)
      CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

    END IF

    CALL timestop(handle)

  END SUBROUTINE atom_int_setup
! *****************************************************************************
!> \brief ...
!> \param integrals ...
!> \param basis ...
!> \param potential ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_ppint_setup(integrals,basis,potential,error)
    TYPE(atom_integrals), INTENT(INOUT)      :: integrals
    TYPE(atom_basis_type), INTENT(INOUT)     :: basis
    TYPE(atom_potential_type), INTENT(IN)    :: potential
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_ppint_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ierr, ii, j, k, l, &
                                                m, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: al, alpha, rc
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cpot, xmat
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: omat, spmat
    REAL(KIND=dp), DIMENSION(:), POINTER     :: rad

    failure = .FALSE.

    CALL timeset(routineN,handle)

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

      ALLOCATE (integrals%clsd(n,n,0:3),STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
      integrals%clsd = 0._dp

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
                IF (potential%gth_pot%lpotextended) THEN
                  DO k=1,potential%gth_pot%nexp_lpot
                    DO i=1,potential%gth_pot%nct_lpot(k)
                      omat = 0._dp
                      CALL sg_gpot ( omat(1:n,1:n), i-1, potential%gth_pot%alpha_lpot(k), l, &
                                     basis%am(1:n,l), basis%am(1:n,l) )
                      integrals%core(1:n,1:n,l) = integrals%core(1:n,1:n,l) + &
                         potential%gth_pot%cval_lpot(i,k)*omat(1:n,1:n)
                    END DO
                  END DO
                END IF
                IF (potential%gth_pot%lsdpot) THEN
                  DO k=1,potential%gth_pot%nexp_lsd
                    DO i=1,potential%gth_pot%nct_lsd(k)
                      omat = 0._dp
                      CALL sg_gpot ( omat(1:n,1:n), i-1, potential%gth_pot%alpha_lsd(k), l, &
                                     basis%am(1:n,l), basis%am(1:n,l) )
                      integrals%clsd(1:n,1:n,l) = integrals%clsd(1:n,1:n,l) + &
                         potential%gth_pot%cval_lsd(i,k)*omat(1:n,1:n)
                    END DO
                  END DO
                END IF

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
                IF(n>0 .AND. m>0) THEN
                  ALLOCATE (omat(m,m),spmat(n,5),xmat(m),STAT=ierr)
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
                  IF (potential%gth_pot%lpotextended) THEN
                    DO k=1,potential%gth_pot%nexp_lpot
                      DO i=1,potential%gth_pot%nct_lpot(k)
                        omat = 0._dp
                        CALL sg_gpot ( omat(1:m,1:m), i-1, potential%gth_pot%alpha_lpot(k), l, &
                                       basis%am(1:m,l), basis%am(1:m,l) )
                        omat(1:m,1:m) = potential%gth_pot%cval_lpot(i,k)*omat(1:m,1:m)
                        CALL contract2add(integrals%core(1:n,1:n,l),omat(1:m,1:m),basis%cm(1:m,1:n,l), error)
                      END DO
                    END DO
                  END IF
                  IF (potential%gth_pot%lsdpot) THEN
                    DO k=1,potential%gth_pot%nexp_lsd
                      DO i=1,potential%gth_pot%nct_lsd(k)
                        omat = 0._dp
                        CALL sg_gpot ( omat(1:m,1:m), i-1, potential%gth_pot%alpha_lsd(k), l, &
                                       basis%am(1:m,l), basis%am(1:m,l) )
                        omat(1:m,1:m) = potential%gth_pot%cval_lsd(i,k)*omat(1:m,1:m)
                        CALL contract2add(integrals%clsd(1:n,1:n,l),omat(1:m,1:m),basis%cm(1:m,1:n,l), error)
                      END DO
                    END DO
                  END IF

                  spmat = 0._dp
                  k = potential%gth_pot%nl(l)
                  DO i=1,k
                    CALL sg_proj_ol ( xmat(1:m), l, basis%am(1:m,l), i-1, potential%gth_pot%rcnl(l) )
                    spmat(1:n,i) = MATMUL(TRANSPOSE(basis%cm(1:m,1:n,l)),xmat(1:m))
                  END DO
                  IF(k>0) THEN
                    integrals%hnl(1:n,1:n,l) = MATMUL(spmat(1:n,1:k),&
                      MATMUL(potential%gth_pot%hnl(1:k,1:k,l),TRANSPOSE(spmat(1:n,1:k))))
                  END IF

                  DEALLOCATE (omat,spmat,xmat,STAT=ierr)
                  CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
                END IF
              END DO
            CASE DEFAULT
              CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT

        CASE (STO_BASIS)

          SELECT CASE (potential%ppot_type)
            CASE (NO_PSEUDO)
              DO l=0,3
                n = integrals%n(l)
                CALL sto_nuclear ( integrals%core(1:n,1:n,l), basis%ns(1:n,l), basis%as(1:n,l),&
                                   basis%ns(1:n,l), basis%as(1:n,l) )
              END DO
            CASE (GTH_PSEUDO)
              rad => basis%grid%rad
              m = basis%grid%nr
              ALLOCATE (cpot(1:m),STAT=ierr)
              CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
              rc = potential%gth_pot%rc
              alpha = 1._dp/rc/SQRT(2._dp)

              ! local pseudopotential, we use erf = 1/r - erfc
              integrals%core = 0._dp
              DO i=1,m
                cpot(i) = potential%gth_pot%zion*erfc(alpha*rad(i))/rad(i)
              END DO
              DO i=1,potential%gth_pot%ncl
                ii = 2*(i-1)
                cpot(1:m) = cpot(1:m) + potential%gth_pot%cl(i) * (rad/rc)**ii * EXP(-0.5_dp*(rad/rc)**2)
              END DO
              IF (potential%gth_pot%lpotextended) THEN
                DO k=1,potential%gth_pot%nexp_lpot
                  al = potential%gth_pot%alpha_lpot(k)
                  DO i=1,potential%gth_pot%nct_lpot(k)
                    ii = 2*(i-1)
                    cpot(1:m) = cpot(1:m) + potential%gth_pot%cval_lpot(i,k)*(rad/al)**ii * EXP(-0.5_dp*(rad/al)**2)
                  END DO
                END DO
              END IF
              CALL numpot_matrix(integrals%core,cpot,basis,0,error)
              DO l=0,3
                n = integrals%n(l)
                ALLOCATE (omat(n,n),STAT=ierr)
                CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
                omat = 0._dp
                CALL sto_nuclear ( omat(1:n,1:n), basis%ns(1:n,l), basis%as(1:n,l),&
                                   basis%ns(1:n,l), basis%as(1:n,l) )
                integrals%core(1:n,1:n,l) = integrals%core(1:n,1:n,l) - potential%gth_pot%zion*omat(1:n,1:n)
                DEALLOCATE (omat,STAT=ierr)
                CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
              END DO

              IF (potential%gth_pot%lsdpot) THEN
                cpot = 0._dp
                DO k=1,potential%gth_pot%nexp_lsd
                  al = potential%gth_pot%alpha_lsd(k)
                  DO i=1,potential%gth_pot%nct_lsd(k)
                    ii = 2*(i-1)
                    cpot(:) = cpot + potential%gth_pot%cval_lsd(i,k)*(rad/al)**ii * EXP(-0.5_dp*(rad/al)**2)
                  END DO
                END DO
                CALL numpot_matrix(integrals%clsd,cpot,basis,0,error)
              END IF

              DO l=0,3
                n = integrals%n(l)
                ! non local pseudopotential
                ALLOCATE (spmat(n,5),STAT=ierr)
                CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
                spmat = 0._dp
                k = potential%gth_pot%nl(l)
                DO i=1,k
                  rc = potential%gth_pot%rcnl(l)
                  cpot(:) = sqrt2/SQRT(gamma1(l+2*i-1))*rad**(l+2*i-2) * EXP(-0.5_dp*(rad/rc)**2) / rc**(l+2*i-0.5_dp)
                  DO j=1,basis%nbas(l)
                    spmat(j,i) = integrate_grid(cpot,basis%bf(:,j,l),basis%grid)
                  END DO
                END DO
                integrals%hnl(1:n,1:n,l) = MATMUL(spmat(1:n,1:k),&
                    MATMUL(potential%gth_pot%hnl(1:k,1:k,l),TRANSPOSE(spmat(1:n,1:k))))

                DEALLOCATE (spmat,STAT=ierr)
                CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

              END DO

              DEALLOCATE (cpot,STAT=ierr)
              CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

            CASE DEFAULT
              CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT

        CASE (NUM_BASIS)
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
      END SELECT

    END IF

    CALL timestop(handle)

  END SUBROUTINE atom_ppint_setup
! *****************************************************************************
!> \brief ...
!> \param integrals ...
!> \param basis ...
!> \param reltyp ...
!> \param zcore ...
!> \param alpha ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_relint_setup(integrals,basis,reltyp,zcore,alpha,error)
    TYPE(atom_integrals), INTENT(INOUT)      :: integrals
    TYPE(atom_basis_type), INTENT(INOUT)     :: basis
    INTEGER, INTENT(IN)                      :: reltyp
    REAL(dp), INTENT(IN)                     :: zcore
    REAL(dp), INTENT(IN), OPTIONAL           :: alpha
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_relint_setup', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: dkhorder, handle, i, ierr, &
                                                k1, k2, l, m, n, nl
    LOGICAL                                  :: failure
    REAL(dp)                                 :: ascal, cspeed
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: cpot
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: modpot
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: ener, sps
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: hmat, pvp, sp, tp, vp, wfn

    failure = .FALSE.

    CALL timeset(routineN,handle)

    cspeed = 1._dp/a_fine

    SELECT CASE (reltyp)
      CASE DEFAULT
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_nonrel_atom,do_zoramp_atom,do_sczoramp_atom)
        dkhorder = -1
      CASE (do_dkh0_atom)
        dkhorder = 0
      CASE (do_dkh1_atom)
        dkhorder = 1
      CASE (do_dkh2_atom)
        dkhorder = 2
      CASE (do_dkh3_atom)
        dkhorder = 3
      CASE (do_dkh4_atom)
        dkhorder = 4
      CASE (do_dkh5_atom)
        dkhorder = 5
    END SELECT

    SELECT CASE (reltyp)
      CASE DEFAULT
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_nonrel_atom)
        ! nothing to do
        NULLIFY (integrals%tzora,integrals%hdkh)
      CASE (do_zoramp_atom, do_sczoramp_atom)

        NULLIFY (integrals%hdkh)

        IF ( integrals%zorastat == 0 ) THEN
          n = MAXVAL(basis%nbas)
          ALLOCATE (integrals%tzora(n,n,0:3),STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          integrals%tzora = 0._dp
          m = basis%grid%nr
          ALLOCATE (modpot(1:m),cpot(1:m),STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          CALL calculate_model_potential(modpot,basis%grid,zcore,error)
          ! Zora potential
          cpot(1:m) = modpot(1:m)/(4._dp*cspeed*cspeed - 2._dp*modpot(1:m))
          cpot(1:m) = cpot(1:m)/basis%grid%rad2(1:m)
          CALL numpot_matrix(integrals%tzora,cpot,basis,0,error)
          DO l=0,3
             nl = basis%nbas(l)
             integrals%tzora(1:nl,1:nl,l) = REAL(l*(l+1),dp) * integrals%tzora(1:nl,1:nl,l)
          END DO
          cpot(1:m) = cpot(1:m)*basis%grid%rad2(1:m)
          CALL numpot_matrix(integrals%tzora,cpot,basis,2,error)
          !
          ! scaled ZORA
          IF( reltyp == do_sczoramp_atom ) THEN
             ALLOCATE (hmat(n,n,0:3),wfn(n,n,0:3),ener(n,0:3),pvp(n,n,0:3),sps(n,n),STAT=ierr)
             CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
             hmat(:,:,:) = integrals%kin + integrals%tzora
             ! model potential
             CALL numpot_matrix(hmat,modpot,basis,0,error)
             ! eigenvalues and eigenvectors
             CALL atom_solve(hmat,integrals%utrans,wfn,ener,basis%nbas,integrals%nne,3,error)
             ! relativistic kinetic energy
             cpot(1:m) = cspeed*cspeed/(2._dp*cspeed*cspeed - modpot(1:m))**2
             cpot(1:m) = cpot(1:m)/basis%grid%rad2(1:m)
             pvp = 0.0_dp
             CALL numpot_matrix(pvp,cpot,basis,0,error)
             DO l=0,3
                nl = basis%nbas(l)
                pvp(1:nl,1:nl,l) = REAL(l*(l+1),dp) * pvp(1:nl,1:nl,l)
             END DO
             cpot(1:m) = cpot(1:m)*basis%grid%rad2(1:m)
             CALL numpot_matrix(pvp,cpot,basis,2,error)
             ! calculate psi*pvp*psi and the scaled orbital energies
             ! actually, we directly calculate the energy difference
             DO l=0,3
                nl = basis%nbas(l)
                DO i=1,integrals%nne(l)
                   IF (ener(i,l) < 0._dp) THEN
                      ascal = SUM(wfn(1:nl,i,l)*MATMUL(pvp(1:nl,1:nl,l),wfn(1:nl,i,l)))
                      ener(i,l) = ener(i,l) * ascal/(1.0_dp + ascal)
                   ELSE
                      ener(i,l) = 0.0_dp
                   END IF
                END DO
             END DO
             ! correction term is calculated as a projector
             hmat = 0.0_dp
             DO l=0,3
                nl = basis%nbas(l)
                DO i=1,integrals%nne(l)
                   DO k1=1,nl
                      DO k2=1,nl
                         hmat(k1,k2,l) = hmat(k1,k2,l) + ener(i,l)*wfn(k1,i,l)*wfn(k2,i,l)
                      END DO
                   END DO
                END DO
                ! transform with overlap matrix
                sps(1:nl,1:nl) = MATMUL(integrals%ovlp(1:nl,1:nl,l),&
                                 MATMUL(hmat(1:nl,1:nl,l),integrals%ovlp(1:nl,1:nl,l)) ) 
                ! add scaling correction to tzora
                integrals%tzora(1:nl,1:nl,l) = integrals%tzora(1:nl,1:nl,l) - sps(1:nl,1:nl)
             END DO

             DEALLOCATE (hmat,wfn,ener,pvp,sps,STAT=ierr)
             CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          END IF
          !
          DEALLOCATE (modpot,cpot,STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

          integrals%zorastat = 1

        END IF

      CASE (do_dkh0_atom, do_dkh1_atom, do_dkh2_atom, do_dkh3_atom, do_dkh4_atom, do_dkh5_atom)

        NULLIFY (integrals%tzora)

        IF ( integrals%dkhstat == 0 ) THEN
          n = MAXVAL(basis%nbas)
          ALLOCATE (integrals%hdkh(n,n,0:3),STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          integrals%hdkh = 0._dp

          m = MAXVAL(basis%nprim)
          ALLOCATE (tp(m,m,0:3),sp(m,m,0:3),vp(m,m,0:3),pvp(m,m,0:3),STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
          tp = 0._dp
          sp = 0._dp
          vp = 0._dp
          pvp = 0._dp

          SELECT CASE (basis%basis_type)
            CASE DEFAULT
              CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
            CASE (GTO_BASIS, CGTO_BASIS)

              DO l=0,3
                m = basis%nprim(l)
                IF ( m > 0 ) THEN
                  CALL sg_kinetic ( tp(1:m,1:m,l), l, basis%am(1:m,l), basis%am(1:m,l) )
                  CALL sg_overlap ( sp(1:m,1:m,l), l, basis%am(1:m,l), basis%am(1:m,l) )
                  IF(PRESENT(alpha)) THEN
                    CALL sg_erfc ( vp(1:m,1:m,l), l, alpha, basis%am(1:m,l), basis%am(1:m,l) )
                  ELSE
                    CALL sg_nuclear ( vp(1:m,1:m,l), l, basis%am(1:m,l), basis%am(1:m,l) )
                  END IF
                  CALL sg_kinnuc ( pvp(1:m,1:m,l), l, basis%am(1:m,l), basis%am(1:m,l) )
                  vp(1:m,1:m,l) = -zcore*vp(1:m,1:m,l)
                  pvp(1:m,1:m,l) = -zcore*pvp(1:m,1:m,l)
                END IF
              END DO

            CASE (STO_BASIS)
              CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
            CASE (NUM_BASIS)
              CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          END SELECT

          CALL dkh_integrals(integrals,basis,dkhorder,sp,tp,vp,pvp,error)

          integrals%dkhstat = 1

          DEALLOCATE (tp,sp,vp,pvp,STAT=ierr)
          CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

        ELSE
          CPPostcondition(ASSOCIATED(integrals%hdkh), cp_failure_level, routineP, error, failure)
        END IF

    END SELECT

    CALL timestop(handle)

  END SUBROUTINE atom_relint_setup
! *****************************************************************************
!> \brief ...
!> \param integrals ...
!> \param basis ...
!> \param order ...
!> \param sp ...
!> \param tp ...
!> \param vp ...
!> \param pvp ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE dkh_integrals(integrals,basis,order,sp,tp,vp,pvp,error)
    TYPE(atom_integrals), INTENT(INOUT)      :: integrals
    TYPE(atom_basis_type), INTENT(INOUT)     :: basis
    INTEGER, INTENT(IN)                      :: order
    REAL(dp), DIMENSION(:, :, 0:)            :: sp, tp, vp, pvp
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dkh_integrals', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: l, m, n
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:, :, :), POINTER    :: hdkh

    failure = .FALSE.
    CPPrecondition(order>=0, cp_failure_level, routineP, error, failure)

    hdkh => integrals%hdkh

    DO l=0,3
      n = integrals%n(l)
      m = basis%nprim(l)
      IF ( n > 0 ) THEN
        CALL dkh_atom_transformation(sp(1:m,1:m,l),vp(1:m,1:m,l),tp(1:m,1:m,l),pvp(1:m,1:m,l),m,order)
        SELECT CASE (basis%basis_type)
          CASE DEFAULT
            CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          CASE (GTO_BASIS)
            CPAssert(n==m,cp_failure_level,routineP,error,failure)
            integrals%hdkh(1:n,1:n,l) = tp(1:n,1:n,l) + vp(1:n,1:n,l)
          CASE (CGTO_BASIS)
            CALL contract2(integrals%hdkh(1:n,1:n,l),tp(1:m,1:m,l),basis%cm(1:m,1:n,l), error)
            CALL contract2add(integrals%hdkh(1:n,1:n,l),vp(1:m,1:m,l),basis%cm(1:m,1:n,l), error)
          CASE (STO_BASIS)
            CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
          CASE (NUM_BASIS)
            CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
        END SELECT
      ELSE
        integrals%hdkh(1:n,1:n,l) = 0._dp
      END IF
    END DO

  END SUBROUTINE dkh_integrals
! *****************************************************************************
!> \brief ...
!> \param integrals ...
!> \param error ...
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
    IF ( ASSOCIATED(integrals%utrans) ) THEN
      DEALLOCATE (integrals%utrans,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF
    IF ( ASSOCIATED(integrals%uptrans) ) THEN
      DEALLOCATE (integrals%uptrans,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF

    integrals%status = 0

  END SUBROUTINE atom_int_release
! *****************************************************************************
!> \brief ...
!> \param integrals ...
!> \param error ...
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
    IF ( ASSOCIATED(integrals%clsd) ) THEN
      DEALLOCATE (integrals%clsd,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF

    integrals%ppstat = 0

  END SUBROUTINE atom_ppint_release
! *****************************************************************************
!> \brief ...
!> \param integrals ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_relint_release(integrals,error)
    TYPE(atom_integrals), INTENT(INOUT)      :: integrals
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_relint_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF ( ASSOCIATED(integrals%tzora) ) THEN
      DEALLOCATE (integrals%tzora,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF
    IF ( ASSOCIATED(integrals%hdkh) ) THEN
      DEALLOCATE (integrals%hdkh,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF

    integrals%zorastat = 0
    integrals%dkhstat = 0

  END SUBROUTINE atom_relint_release
! *****************************************************************************
!> \brief ...
!> \param modpot ...
!> \param grid ...
!> \param zcore ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_model_potential(modpot,grid,zcore,error)
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: modpot
    TYPE(grid_atom_type), INTENT(IN)         :: grid
    REAL(dp), INTENT(IN)                     :: zcore
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_model_potential', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, ii, l, ll, n, z
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: pot, rho
    TYPE(atom_state)                         :: state

    n = SIZE(modpot)
    ALLOCATE(rho(n),pot(n),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    ! fill default occupation
    state%core = 0._dp
    state%occ = 0._dp
    state%multiplicity = -1
    z = NINT(zcore)
    DO l=0,3
      IF ( ptable(z)%e_conv(l) /= 0 ) THEN
        state%maxl_occ = l
        ll = 2*(2*l+1)
        DO i=1,SIZE(state%occ,2)
          ii = ptable(z)%e_conv(l) - (i-1)*ll
          IF ( ii <= ll ) THEN
            state%occ(l,i) = ii
            EXIT
          ELSE
            state%occ(l,i) = ll
          END IF
        END DO
      END IF
    END DO

    modpot = -zcore/grid%rad(:)

    ! Coulomb potential
    CALL slater_density(rho,pot,NINT(zcore),state,grid,error)
    CALL coulomb_potential_numeric(pot,rho,grid,error)
    modpot = modpot + pot

    ! XC potential
    CALL wigner_slater_functional(rho,pot,error)
    modpot = modpot + pot

    DEALLOCATE(rho,pot,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE calculate_model_potential
! *****************************************************************************
!> \brief ...
!> \param int ...
!> \param omat ...
!> \param cm ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE contract2 ( int, omat, cm, error )
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: int
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: omat, cm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'contract2', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, m, n

    CALL timeset(routineN,handle)

    n = SIZE(int,1)
    m = SIZE(omat,1)

    INT(1:n,1:n) = MATMUL(TRANSPOSE(cm(1:m,1:n)),MATMUL(omat(1:m,1:m),cm(1:m,1:n)))

    CALL timestop(handle)

  END SUBROUTINE contract2
! *****************************************************************************
!> \brief ...
!> \param int ...
!> \param omat ...
!> \param cm ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE contract2add ( int, omat, cm, error )
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: int
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: omat, cm
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'contract2add', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, m, n

    CALL timeset(routineN,handle)

    n = SIZE(int,1)
    m = SIZE(omat,1)

    INT(1:n,1:n) = INT(1:n,1:n) + MATMUL(TRANSPOSE(cm(1:m,1:n)),MATMUL(omat(1:m,1:m),cm(1:m,1:n)))

    CALL timestop(handle)

  END SUBROUTINE contract2add
! *****************************************************************************
!> \brief ...
!> \param eri ...
!> \param omat ...
!> \param cm1 ...
!> \param cm2 ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE contract4 ( eri, omat, cm1, cm2, error )
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: eri
    REAL(dp), DIMENSION(:, :), INTENT(IN)    :: omat, cm1, cm2
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'contract4', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i1, i2, ierr, m1, m2, &
                                                mm1, mm2, n1, n2, nn1, nn2
    LOGICAL                                  :: failure
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: amat, atran, bmat, btran, hint

    failure = .FALSE.
    CALL timeset(routineN,handle)

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

    CALL timestop(handle)

  END SUBROUTINE contract4
! *****************************************************************************
!> \brief ...
!> \param mat ...
!> \param vec ...
!> \param n ...
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

! *****************************************************************************
!> \brief ...
!> \param mat ...
!> \param vec ...
!> \param n ...
! *****************************************************************************
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
