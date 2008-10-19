!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Some basic routines for atomic calculations
!> \author  jgh
!> \date    01.04.2008
!> \version 1.0
!>
! *****************************************************************************
MODULE atom_utils
  USE atom_types,                      ONLY: CGTO_BASIS,&
                                             GTO_BASIS,&
                                             atom_basis_type,&
                                             eri
  USE erf_fn,                          ONLY: erf
  USE f77_blas
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE lapack
  USE mathconstants,                   ONLY: fac,&
                                             fourpi,&
                                             maxfac,&
                                             rootpi
  USE periodic_table,                  ONLY: nelem,&
                                             ptable
  USE qs_grid_atom,                    ONLY: grid_atom_type
  USE string_utilities,                ONLY: uppercase
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_utils'

  PUBLIC :: atom_set_occupation, get_maxl_occ, get_maxn_occ
  PUBLIC :: atom_denmat, atom_density
  PUBLIC :: prune_grid, integrate_grid, atom_trace, atom_solve
  PUBLIC :: coulomb_potential_numeric, coulomb_potential_analytic
  PUBLIC :: numpot_matrix, ceri_contract, eeri_contract, err_matrix

!-----------------------------------------------------------------------------!

  INTERFACE integrate_grid
     MODULE PROCEDURE integrate_grid_function1, &
                      integrate_grid_function2, &
                      integrate_grid_function3
  END INTERFACE


! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE atom_set_occupation(ostring,occupation,wfnocc,error)
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: ostring
    REAL(Kind=dp), DIMENSION(0:3, 10)        :: occupation, wfnocc
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_set_occupation', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=2)                         :: elem
    CHARACTER(LEN=default_string_length)     :: pstring
    INTEGER                                  :: i, i1, i2, ielem, is, jd, jf, &
                                                jp, js, k, l, n
    LOGICAL                                  :: failure
    REAL(Kind=dp)                            :: e0, el, oo

    failure = .FALSE.

    occupation = 0._dp

    CPPrecondition(ASSOCIATED(ostring),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ostring) > 0,cp_failure_level,routineP,error,failure)

    is = 1
    IF ( INDEX(ostring(1),"CORE") /= 0 ) is = 2   !Pseudopotential detected
    IF ( INDEX(ostring(1),"none") /= 0 ) is = 2   !no electrons, used with CORE
    IF ( INDEX(ostring(1),"[") /= 0 ) THEN
      ! core occupation from element [XX]
      i1 = INDEX(ostring(1),"[")
      i2 = INDEX(ostring(1),"]")
      CPAssert((i2-i1-1 > 0).AND.(i2-i1-1 < 3),cp_failure_level,routineP,error,failure)
      elem = ostring(1)(i1+1:i2-1)
      ielem=0
      DO k=1,nelem
        IF(elem==ptable(k)%symbol) THEN
          ielem=k
          EXIT
        END IF
      END DO
      CPAssert(ielem /= 0,cp_failure_level,routineP,error,failure)
      DO l=0,3
        el = 2._dp * (2._dp*REAL(l,dp) + 1._dp)
        e0 = ptable(ielem)%e_conv(l)
        DO k=1,10
          occupation(l,k) = MIN(el,e0)
          e0 = e0 - el
          IF ( e0 <= 0._dp ) EXIT
        END DO
      END DO
      is = 2
    END IF
    DO i=is,SIZE(ostring)
      pstring = ostring(i)
      CALL uppercase(pstring)
      js = INDEX(pstring,"S")
      jp = INDEX(pstring,"P")
      jd = INDEX(pstring,"D")
      jf = INDEX(pstring,"F")
      CPAssert(js+jp+jd+jf > 0,cp_failure_level,routineP,error,failure)
      IF ( js > 0 ) THEN
        CPAssert(jp+jd+jf == 0,cp_failure_level,routineP,error,failure)
        READ(pstring(1:js-1),*) n
        READ(pstring(js+1:),*) oo
        CPAssert(n > 0,cp_failure_level,routineP,error,failure)
        CPAssert(oo >= 0._dp,cp_failure_level,routineP,error,failure)
        CPAssert(occupation(0,n)==0,cp_failure_level,routineP,error,failure)
        occupation(0,n) = oo
      END IF
      IF ( jp > 0 ) THEN
        CPAssert(js+jd+jf == 0,cp_failure_level,routineP,error,failure)
        READ(pstring(1:jp-1),*) n
        READ(pstring(jp+1:),*) oo
        CPAssert(n > 1,cp_failure_level,routineP,error,failure)
        CPAssert(oo >= 0._dp,cp_failure_level,routineP,error,failure)
        CPAssert(occupation(1,n-1)==0,cp_failure_level,routineP,error,failure)
        occupation(1,n-1) = oo
      END IF
      IF ( jd > 0 ) THEN
        CPAssert(js+jp+jf == 0,cp_failure_level,routineP,error,failure)
        READ(pstring(1:jd-1),*) n
        READ(pstring(jd+1:),*) oo
        CPAssert(n > 2,cp_failure_level,routineP,error,failure)
        CPAssert(oo >= 0._dp,cp_failure_level,routineP,error,failure)
        CPAssert(occupation(2,n-2)==0,cp_failure_level,routineP,error,failure)
        occupation(2,n-2) = oo
      END IF
      IF ( jf > 0 ) THEN
        CPAssert(js+jp+jd == 0,cp_failure_level,routineP,error,failure)
        READ(pstring(1:jf-1),*) n
        READ(pstring(jf+1:),*) oo
        CPAssert(n > 3,cp_failure_level,routineP,error,failure)
        CPAssert(oo >= 0._dp,cp_failure_level,routineP,error,failure)
        CPAssert(occupation(3,n-3)==0,cp_failure_level,routineP,error,failure)
        occupation(3,n-3) = oo
      END IF
      
    END DO

    wfnocc = 0._dp
    DO l=0,3
      k = 0
      DO i=1,10
        IF ( occupation(l,i) /= 0._dp ) THEN
          k = k + 1
          wfnocc(l,k) = occupation(l,i)
        END IF
      END DO
    END DO

  END SUBROUTINE atom_set_occupation

  FUNCTION get_maxl_occ(occupation) RESULT(maxl)
    REAL(Kind=dp), DIMENSION(0:3, 10)        :: occupation
    INTEGER                                  :: maxl

    INTEGER                                  :: l

    maxl = 0
    DO l=0,3
      IF(SUM(occupation(l,:)) /= 0._dp) maxl=l
    END DO

  END FUNCTION get_maxl_occ

  FUNCTION get_maxn_occ(occupation) RESULT(maxn)
    REAL(Kind=dp), DIMENSION(0:3, 10)        :: occupation
    INTEGER, DIMENSION(0:3)                  :: maxn

    INTEGER                                  :: k, l

    maxn = 0
    DO l=0,3
      DO k=1,10
        IF(occupation(l,k) /= 0._dp) maxn(l)=maxn(l)+1
      END DO
    END DO

  END FUNCTION get_maxn_occ

! *****************************************************************************
  SUBROUTINE atom_denmat(pmat,wfn,nbas,occ,maxl,maxn,error)
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(INOUT)                          :: pmat
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: wfn
    INTEGER, DIMENSION(0:3), INTENT(IN)      :: nbas
    REAL(KIND=dp), DIMENSION(0:, :), &
      INTENT(IN)                             :: occ
    INTEGER, INTENT(IN)                      :: maxl
    INTEGER, DIMENSION(0:3), INTENT(IN)      :: maxn
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_denmat', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, l

    pmat = 0._dp
    DO l=0,maxl
      DO i=1,maxn(l)
        DO k=1,nbas(l)
          DO j=1,nbas(l)
            pmat(j,k,l) =  pmat(j,k,l) + occ(l,i)*wfn(j,i,l)*wfn(k,i,l)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE atom_denmat
! *****************************************************************************
  SUBROUTINE atom_density(density,pmat,basis,maxl,typ,rr,error)
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: density
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: pmat
    TYPE(atom_basis_type), INTENT(IN)        :: basis
    INTEGER, INTENT(IN)                      :: maxl
    CHARACTER(LEN=*), OPTIONAL               :: typ
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(IN), OPTIONAL                   :: rr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_density', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=3)                         :: my_typ = "RHO"
    INTEGER                                  :: i, j, l, n
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: ff

    IF(PRESENT(typ)) my_typ = typ(1:3)
    IF(my_typ == "KIN") THEN
      CPPrecondition(PRESENT(rr),cp_failure_level,routineP,error,failure)
    END IF

    density=0._dp
    DO l=0,maxl
      n=basis%nbas(l)
      DO i=1,n
        DO j=i,n
          ff=pmat(i,j,l)
          IF(i/=j) ff=2._dp*pmat(i,j,l)
          IF(my_typ == "RHO") THEN
            density(:) = density(:) + ff*basis%bf(:,i,l)*basis%bf(:,j,l)
          ELSE IF(my_typ == "DER") THEN
            density(:) = density(:) + ff*basis%dbf(:,i,l)*basis%bf(:,j,l) &
                                    + ff*basis%bf(:,i,l)*basis%dbf(:,j,l)
          ELSE IF(my_typ == "KIN") THEN
            density(:) = density(:) + 0.5_dp*ff*(&
              basis%dbf(:,i,l)*basis%dbf(:,j,l) + &
              REAL(l*(l+1),dp)*basis%bf(:,i,l)*basis%bf(:,j,l)/rr(:))
          ELSE
            CPPrecondition(.FALSE.,cp_failure_level,routineP,error,failure)
          END IF
        END DO
      END DO
    END DO
    ! this factor from the product of two spherical harmonics
    density=density/fourpi

  END SUBROUTINE atom_density
! *****************************************************************************
  SUBROUTINE atom_solve(hmat,umat,orb,ener,nb,nv,maxl,error)
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: hmat, umat
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(INOUT)                          :: orb
    REAL(KIND=dp), DIMENSION(:, 0:), &
      INTENT(INOUT)                          :: ener
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nb, nv
    INTEGER, INTENT(IN)                      :: maxl
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_solve', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, info, l, lwork, m, n
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: w, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: a

    CPPrecondition(ALL(n >= nv),cp_failure_level,routineP,error,failure)

    orb = 0._dp
    DO l=0,maxl
      n = nb(l)
      m = nv(l)
      IF ( n > 0 .AND. m > 0 ) THEN
        lwork = 10*m
        ALLOCATE(a(n,n),w(n),work(lwork),STAT=ierr)
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
        a(1:m,1:m) = MATMUL(TRANSPOSE(umat(1:n,1:m,l)),MATMUL(hmat(1:n,1:n,l),umat(1:n,1:m,l)))
        CALL lapack_ssyev ( "V", "U", m, a(1:m,1:m), m, w(1:m), work, lwork, info )
        a(1:n,1:m) = MATMUL(umat(1:n,1:m,l),a(1:m,1:m))

        m = MIN(m,SIZE(orb,2))
        orb(1:n,1:m,l) = a(1:n,1:m)
        ener(1:m,l) = w(1:m)

        DEALLOCATE(a,w,work,STAT=ierr)
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
      END IF
    END DO

  END SUBROUTINE atom_solve
! *****************************************************************************

  FUNCTION prune_grid ( fun, deps ) RESULT(nc)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fun
    REAL(KIND=dp), INTENT(IN), OPTIONAL      :: deps
    INTEGER                                  :: nc

    INTEGER                                  :: i, nr
    REAL(KIND=dp)                            :: meps = 1.e-14_dp

    IF ( PRESENT(deps) ) meps=deps

    nr = SIZE(fun)
    nc = 0
    DO i=nr,1,-1
      IF (ABS(fun(i)) > meps) THEN
        nc = i
        EXIT
      END IF
    END DO

  END FUNCTION prune_grid
! *****************************************************************************
  FUNCTION integrate_grid_function1 ( fun, grid ) RESULT(integral)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fun
    TYPE(grid_atom_type)                     :: grid
    REAL(KIND=dp)                            :: integral

    INTEGER                                  :: nc

    nc = SIZE(fun)
    integral = SUM ( fun(1:nc) * grid%wr(1:nc) )
    
  END FUNCTION integrate_grid_function1

  FUNCTION integrate_grid_function2 ( fun1, fun2, grid ) RESULT(integral)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fun1, fun2
    TYPE(grid_atom_type)                     :: grid
    REAL(KIND=dp)                            :: integral

    INTEGER                                  :: nc

    nc = MIN(SIZE(fun1),SIZE(fun2))
    integral = SUM ( fun1(1:nc) * fun2(1:nc) * grid%wr(1:nc) )

  END FUNCTION integrate_grid_function2

  FUNCTION integrate_grid_function3 ( fun1, fun2, fun3, grid ) RESULT(integral)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: fun1, fun2, fun3
    TYPE(grid_atom_type)                     :: grid
    REAL(KIND=dp)                            :: integral

    INTEGER                                  :: nc

    nc = MIN(SIZE(fun1),SIZE(fun2),SIZE(fun3))
    integral = SUM ( fun1(1:nc) * fun2(1:nc) * fun3(1:nc) * grid%wr(1:nc) )

  END FUNCTION integrate_grid_function3
! *****************************************************************************
  SUBROUTINE coulomb_potential_numeric ( cpot, density, grid, error )
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: cpot
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: density
    TYPE(grid_atom_type), INTENT(IN)         :: grid
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'coulomb_potential_numeric', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, nc
    LOGICAL                                  :: failure = .FALSE.
    REAL(dp)                                 :: int1, int2
    REAL(dp), DIMENSION(:), POINTER          :: r, wr

    nc = MIN(SIZE(cpot),SIZE(density))
    r => grid%rad
    wr => grid%wr

    int1 = fourpi * integrate_grid ( density, grid )
    int2 = 0._dp
    cpot(nc:) = int1/r(nc:)

    ! test that grid is decreasing
    CPPostcondition(r(1)>r(nc),cp_failure_level,routineP,error,failure)
    DO i = 1, nc
      cpot(i) = int1/r(i) + int2
      int1 = int1 - fourpi * density(i)*wr(i)
      int2 = int2 + fourpi * density(i)*wr(i)/r(i)
    END DO

  END SUBROUTINE coulomb_potential_numeric

! *****************************************************************************

  SUBROUTINE coulomb_potential_analytic ( cpot, pmat, basis, grid, maxl, error )
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(INOUT)                          :: cpot
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: pmat
    TYPE(atom_basis_type), INTENT(IN)        :: basis
    TYPE(grid_atom_type)                     :: grid
    INTEGER, INTENT(IN)                      :: maxl
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'coulomb_potential_analytic', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, j, k, l, m, n
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: a, b, ff, oab, sab
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: erfa, expa, z
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: unp

    m = SIZE(cpot)
    ALLOCATE ( erfa(1:m), expa(1:m), z(1:m), STAT=ierr )
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

    cpot = 0._dp

    DO l = 0, maxl
      IF ( MAXVAL(ABS(pmat(:,:,l))) < 1.e-14_dp ) CYCLE
      SELECT CASE ( basis%basis_type )
        CASE DEFAULT
          CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
        CASE ( GTO_BASIS )
          DO i = 1, basis%nbas(l)
            DO j = i, basis%nbas(l)
              IF ( ABS(pmat(i,j,l)) < 1.e-14_dp ) CYCLE
              ff = pmat(i,j,l)
              IF (i /= j) ff = 2._dp*ff
              a = basis%am(i,l)
              b = basis%am(j,l)
              sab = SQRT(a+b)
              oab = rootpi/(a+b)**(l+1.5_dp)*ff
              z(:) = sab*grid%rad(:)
              DO k=1,SIZE(erfa)
                erfa(k)= oab*erf(z(k))/grid%rad(k)
              END DO
              expa(:)= EXP(-z(:)**2)*ff/(a+b)**(l+1)
              SELECT CASE (l)
                CASE DEFAULT
                  CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                CASE (0)
                  cpot(:) = cpot(:) + 0.25_dp*erfa(:)
                CASE (1)
                  cpot(:) = cpot(:) + 0.375_dp*erfa(:)-0.25_dp*expa(:)
                CASE (2)
                  cpot(:) = cpot(:) + 0.9375_dp*erfa(:)-expa(:)*(0.875_dp+0.25_dp*z(:)**2)
                CASE (3)
                  cpot(:) = cpot(:) + 3.28125_dp*erfa(:)-expa(:)*(3.5625_dp+1.375_dp*z(:)**2+0.25*z(:)**4)
              END SELECT
            END DO
          END DO
        CASE ( CGTO_BASIS )
          n = basis%nprim(l)
          m = basis%nbas(l)
          ALLOCATE ( unp(n,n), STAT=ierr )
          CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

          unp(1:n,1:n) = MATMUL ( MATMUL(basis%cm(1:n,1:m,l),pmat(1:m,1:m,l)),&
                                  TRANSPOSE(basis%cm(1:n,1:m,l)))
          DO i = 1, basis%nprim(l)
            DO j = i, basis%nprim(l)
              IF ( ABS(unp(i,j)) < 1.e-14_dp ) CYCLE
              ff = unp(i,j)
              IF (i /= j) ff = 2._dp*ff
              a = basis%am(i,l)
              b = basis%am(j,l)
              sab = SQRT(a+b)
              oab = rootpi/(a+b)**(l+1.5_dp)*ff
              z(:) = sab*grid%rad(:)
              DO k=1,SIZE(erfa)
                erfa(k)= oab*erf(z(k))/grid%rad(k)
              END DO
              expa(:)= EXP(-z(:)**2)*ff/(a+b)**(l+1)
              SELECT CASE (l)
                CASE DEFAULT
                  CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
                CASE (0)
                  cpot(:) = cpot(:) + 0.25_dp*erfa(:)
                CASE (1)
                  cpot(:) = cpot(:) + 0.375_dp*erfa(:)-0.25_dp*expa(:)
                CASE (2)
                  cpot(:) = cpot(:) + 0.9375_dp*erfa(:)-expa(:)*(0.875_dp+0.25_dp*z(:)**2)
                CASE (3)
                  cpot(:) = cpot(:) + 3.28125_dp*erfa(:)-expa(:)*(3.5625_dp+1.375_dp*z(:)**2+0.25*z(:)**4)
              END SELECT
            END DO
          END DO

          DEALLOCATE ( unp, STAT=ierr )
          CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
      END SELECT
    END DO
    DEALLOCATE ( erfa, expa, z, STAT=ierr )
    CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

  END SUBROUTINE coulomb_potential_analytic

! *****************************************************************************

  FUNCTION atom_trace(opmat,pmat,error) RESULT(trace)
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: opmat, pmat
    TYPE(cp_error_type), INTENT(INOUT)       :: error
    REAL(KIND=dp)                            :: trace

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_trace', &
      routineP = moduleN//':'//routineN

    trace = SUM(opmat(:,:,:)*pmat(:,:,:))

  END FUNCTION atom_trace

! *****************************************************************************
  SUBROUTINE numpot_matrix(imat,cpot,basis,derivatives,error)
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(INOUT)                          :: imat
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: cpot
    TYPE(atom_basis_type), INTENT(INOUT)     :: basis
    INTEGER, INTENT(IN)                      :: derivatives
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'numpot_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, l, n
    LOGICAL                                  :: failure

    failure = .FALSE.

    SELECT CASE (derivatives)
      CASE (0)
        DO l=0,3
          n=basis%nbas(l)
          DO i=1,n
            DO j=i,n
              imat(i,j,l) = imat(i,j,l) + &
                  integrate_grid ( cpot, basis%bf(:,i,l), basis%bf(:,j,l), basis%grid )
              imat(j,i,l) = imat(i,j,l)
            END DO
          END DO
        END DO
      CASE (1)
        DO l=0,3
          n=basis%nbas(l)
          DO i=1,n
            DO j=i,n
              imat(i,j,l) = imat(i,j,l) + &
                  integrate_grid ( cpot, basis%dbf(:,i,l), basis%bf(:,j,l), basis%grid )
              imat(i,j,l) = imat(i,j,l) + &
                  integrate_grid ( cpot, basis%bf(:,i,l), basis%dbf(:,j,l), basis%grid )
              imat(j,i,l) = imat(i,j,l)
            END DO
          END DO
        END DO
      CASE (2)
        DO l=0,3
          n=basis%nbas(l)
          DO i=1,n
            DO j=i,n
              imat(i,j,l) = imat(i,j,l) + &
                 integrate_grid ( cpot, basis%dbf(:,i,l), basis%dbf(:,j,l), basis%grid )
              imat(j,i,l) = imat(i,j,l)
            END DO
          END DO
        END DO
      CASE DEFAULT
        CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

  END SUBROUTINE numpot_matrix
! *****************************************************************************
  SUBROUTINE ceri_contract(jmat,erint,pmat,nsize,all_nu,error)
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(INOUT)                          :: jmat
    TYPE(eri), DIMENSION(:), INTENT(IN)      :: erint
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: pmat
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nsize
    LOGICAL, INTENT(IN), OPTIONAL            :: all_nu
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'ceri_contract', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i2, ij1, ij2, j1, j2, l1, &
                                                l2, ll, n1, n2
    LOGICAL                                  :: failure, have_all_nu
    REAL(KIND=dp)                            :: eint, f1, f2

    failure = .FALSE.
    IF ( PRESENT(all_nu) ) THEN
      have_all_nu = all_nu
    ELSE
      have_all_nu = .FALSE.
    END IF

    jmat(:,:,:) = 0._dp
    ll = 0
    DO l1=0,3
      n1 = nsize(l1)
      DO l2=0,l1
        n2 = nsize(l2)
        ll = ll + 1
        ij1=0
        DO i1=1,n1
          DO j1=i1,n1
            ij1=ij1+1
            f1 = 1._dp
            IF(i1/=j1) f1 = 2._dp
            ij2=0
            DO i2=1,n2
              DO j2=i2,n2
                ij2=ij2+1
                f2 = 1._dp
                IF(i2/=j2) f2 = 2._dp
                eint = erint(ll)%int(ij1,ij2)
                IF ( l1==l2 ) THEN
                  jmat(i1,j1,l1) = jmat(i1,j1,l1) + f2*pmat(i2,j2,l2)*eint
                ELSE
                  jmat(i1,j1,l1) = jmat(i1,j1,l1) + f2*pmat(i2,j2,l2)*eint
                  jmat(i2,j2,l2) = jmat(i2,j2,l2) + f1*pmat(i1,j1,l1)*eint
                END IF
              END DO
            END DO
          END DO
        END DO
        IF ( have_all_nu ) THEN
          ! skip integral blocks with nu/=0
          ll = ll + l2
        END IF
      END DO
    END DO
    DO l1=0,3
      n1 = nsize(l1)
      DO i1=1,n1
        DO j1=i1,n1
          jmat(j1,i1,l1) = jmat(i1,j1,l1)
        END DO
      END DO
    END DO

  END SUBROUTINE ceri_contract
! *****************************************************************************
  SUBROUTINE eeri_contract(kmat,erint,pmat,nsize,error)
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(INOUT)                          :: kmat
    TYPE(eri), DIMENSION(:), INTENT(IN)      :: erint
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: pmat
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nsize
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eeri_contract', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i1, i2, ij1, ij2, j1, j2, l1, &
                                                l2, lh, ll, n1, n2, nu
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: almn, eint, f1, f2
    REAL(KIND=dp), DIMENSION(0:maxfac)       :: arho

    failure = .FALSE.

    arho = 0._dp
    DO ll=0,maxfac,2
      lh=ll/2
      arho(ll)=fac(ll)/fac(lh)**2
    END DO

    kmat(:,:,:) = 0._dp
    ll = 0
    DO l1=0,3
      n1 = nsize(l1)
      DO l2=0,l1
        n2 = nsize(l2)
        DO nu = ABS(l1-l2),l1+l2,2
          almn = arho(-l1+l2+nu)*arho(l1-l2+nu)*arho(l1+l2-nu)/(REAL(l1+l2+nu+1,dp)*arho(l1+l2+nu))
          almn = -0.5_dp*almn
          ll = ll + 1
          ij1=0
          DO i1=1,n1
            DO j1=i1,n1
              ij1=ij1+1
              f1 = 1._dp
              IF(i1/=j1) f1 = 2._dp
              ij2=0
              DO i2=1,n2
                DO j2=i2,n2
                  ij2=ij2+1
                  f2 = 1._dp
                  IF(i2/=j2) f2 = 2._dp
                  eint = erint(ll)%int(ij1,ij2)
                  IF ( l1==l2 ) THEN
                    kmat(i1,j1,l1) = kmat(i1,j1,l1) + f2*almn*pmat(i2,j2,l2)*eint
                  ELSE
                    kmat(i1,j1,l1) = kmat(i1,j1,l1) + f2*almn*pmat(i2,j2,l2)*eint
                    kmat(i2,j2,l2) = kmat(i2,j2,l2) + f1*almn*pmat(i1,j1,l1)*eint
                  END IF
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
    DO l1=0,3
      n1 = nsize(l1)
      DO i1=1,n1
        DO j1=i1,n1
          kmat(j1,i1,l1) = kmat(i1,j1,l1)
        END DO
      END DO
    END DO

  END SUBROUTINE eeri_contract
! *****************************************************************************
  SUBROUTINE err_matrix(emat,demax,kmat,pmat,umat,upmat,nval,nbs,error)
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(OUT)                            :: emat
    REAL(KIND=dp), INTENT(OUT)               :: demax
    REAL(KIND=dp), DIMENSION(:, :, 0:), &
      INTENT(IN)                             :: kmat, pmat, umat, upmat
    INTEGER, DIMENSION(0:), INTENT(IN)       :: nval, nbs
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'err_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, l, m, n
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: tkmat, tpmat

    emat = 0._dp
    DO l=0,3
      n = nval(l)
      m = nbs(l)
      IF (m > 0) THEN
        ALLOCATE(tkmat(1:m,1:m),tpmat(1:m,1:m),STAT=ierr)
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
        tkmat = 0._dp
        tpmat = 0._dp
        tkmat(1:m,1:m) = MATMUL(TRANSPOSE(umat(1:n,1:m,l)),MATMUL(kmat(1:n,1:n,l),umat(1:n,1:m,l)))
        tpmat(1:m,1:m) = MATMUL(TRANSPOSE(umat(1:n,1:m,l)),MATMUL(pmat(1:n,1:n,l),umat(1:n,1:m,l)))
        tpmat(1:m,1:m) = MATMUL(upmat(1:m,1:m,l),MATMUL(tpmat(1:m,1:m),upmat(1:m,1:m,l)))

        emat(1:m,1:m,l) = MATMUL(tkmat(1:m,1:m),tpmat(1:m,1:m)) - MATMUL(tpmat(1:m,1:m),tkmat(1:m,1:m))

        DEALLOCATE(tkmat,tpmat,STAT=ierr)
        CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
      END IF
    END DO
    demax = MAXVAL(ABS(emat))

  END SUBROUTINE err_matrix
! *****************************************************************************

END MODULE atom_utils
