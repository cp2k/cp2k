!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that fit parameters for /from atomic calculations
! *****************************************************************************
MODULE atom_fit

  USE atom_types,                      ONLY: atom_type,&
                                             create_opgrid,&
                                             opgrid_type,&
                                             release_opgrid
  USE atom_utils,                      ONLY: atom_denmat,&
                                             atom_density,&
                                             integrate_grid
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE lapack,                          ONLY: lapack_sgesv
  USE mathconstants,                   ONLY: fourpi,&
                                             pi
  USE powell,                          ONLY: opt_state_type,&
                                             powell_optimize
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_fit'

  PUBLIC :: atom_fit_density

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE atom_fit_density (atom,atom_section,num_gto,iunit,error)
    TYPE(atom_type), POINTER                 :: atom
    TYPE(section_vals_type), POINTER         :: atom_section
    INTEGER, INTENT(IN)                      :: num_gto, iunit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_fit_density', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, n10
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: co
    REAL(KIND=dp), DIMENSION(2)              :: x
    TYPE(opgrid_type), POINTER               :: density
    TYPE(opt_state_type)                     :: ostate

    ALLOCATE(co(num_gto),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    co = 0._dp
    NULLIFY(density)
    CALL create_opgrid(density,atom%basis%grid,error)
    CALL atom_denmat(atom%orbitals%pmat,atom%orbitals%wfn,atom%basis%nbas,atom%state%occupation,&
                     atom%state%maxl_occ,atom%state%maxn_occ,error)
    CALL atom_density(density%op,atom%orbitals%pmat,atom%basis,atom%state%maxl_occ,&
                      typ="RHO",error=error)
    density%op = fourpi*density%op

    ostate%nvar = 2
    x(1) = 0.10_dp     !starting point of geometric series
    x(2) = 2.00_dp      !factor of series
    ostate%rhobeg = 0.005_dp
    ostate%rhoend = 0.000000001_dp
    ostate%maxfun = 5000
    ostate%iprint = 1
    ostate%unit  = iunit

    ostate%state = 0
    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(/," POWELL| Start optimization procedure")') 
    END IF
    n10 = ostate%maxfun/10

    DO

      IF ( ostate%state == 2 ) THEN
        CALL density_fit (density,atom,num_gto,x(1),x(2),co,ostate%f,error)
      END IF

      IF ( ostate%state == -1 ) EXIT

      CALL powell_optimize (ostate%nvar, x, ostate)

      IF ( MOD(ostate%nf,n10) == 0 .AND. iunit > 0 ) THEN
        WRITE(iunit,'(/," POWELL| Reached",i4,"% of maximal function calls")') ostate%nf/n10*10
      END IF

    END DO

    ostate%state = 8
    CALL powell_optimize (ostate%nvar, x, ostate)

    CALL release_opgrid(density,error)

    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(" POWELL| Number of function evaluations",T71,I10)') ostate%nf
      WRITE(iunit,'(" POWELL| Final value of function",T61,G20.10)') ostate%f
      WRITE(iunit,'(" Optimized representation of density using a Geometrical Gaussian basis")')
      WRITE(iunit,'(A,I3,/,T10,A,F10.6,T48,A,F10.6)') " Number of Gaussians:",num_gto,&
             "Starting exponent:",x(1),"Proportionality factor:",x(2)
      WRITE(iunit,'(A)') " Expansion coefficients"
      WRITE(iunit,'(4F20.10)') co(1:num_gto)
    END IF

    DEALLOCATE(co,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE atom_fit_density
! *****************************************************************************
  SUBROUTINE density_fit (density,atom,n,aval,cval,co,aerr,error)
    TYPE(opgrid_type), POINTER               :: density
    TYPE(atom_type), POINTER                 :: atom
    INTEGER, INTENT(IN)                      :: n
    REAL(dp), INTENT(IN)                     :: aval, cval
    REAL(dp), DIMENSION(:), INTENT(INOUT)    :: co
    REAL(dp), INTENT(OUT)                    :: aerr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'density_fit', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, info, ip, &
                                                ipiv(1000), iq, k, nr
    LOGICAL                                  :: failure = .FALSE.
    REAL(dp)                                 :: a, rk, zval
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: den, pe, uval
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: bf, smat, tval

! create basis

    ALLOCATE(pe(n),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    nr = atom%basis%grid%nr
    ALLOCATE (bf(nr,n),den(nr),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    bf =  0._dp

    DO i=1,n
      pe(i) = aval * cval**i
      a = pe(i)
      DO k=1,nr
        rk  = atom%basis%grid%rad(k)
        bf(k,i) = EXP(-a*rk**2)
      END DO
    END DO

    ! total charge
    zval = integrate_grid(density%op,atom%basis%grid)

    ! allocate vectors and matrices for overlaps
    ALLOCATE(tval(n+1,1),uval(n),smat(n+1,n+1),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    DO i=1,n
      uval(i) = (pi/pe(i))**1.5_dp
      tval(i,1) = integrate_grid(density%op,bf(:,i),atom%basis%grid)
    END DO
    tval(n+1,1) = zval

    DO iq = 1, n
      DO ip = 1, n
        smat(ip,iq) = (pi/(pe(ip)+pe(iq)))**1.5_dp
      END DO
    END DO
    smat(1:n,n+1) = uval(1:n)
    smat(n+1,1:n) = uval(1:n)
    smat(n+1,n+1) = 0._dp

    CALL lapack_sgesv ( n+1, 1, smat, n+1, ipiv, tval, n+1, info )
    CPPostcondition(info==0, cp_failure_level, routineP, error, failure)
    co(1:n) = tval(1:n,1)

    ! calculate density
    den = 0._dp
    DO i=1,n
      den(:) = den(:) + co(i)*bf(:,i)
    END DO
    den = den * fourpi
    den(:) = (den(:)-density%op(:))**2
    aerr = SQRT(integrate_grid(den,atom%basis%grid))

    DEALLOCATE(pe,bf,den,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    DEALLOCATE(tval,uval,smat,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE density_fit
! *****************************************************************************

END MODULE atom_fit
