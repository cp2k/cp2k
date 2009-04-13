!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief routines that fit parameters for /from atomic calculations
! *****************************************************************************
MODULE atom_fit
  USE atom_electronic_structure,       ONLY: calculate_atom
  USE atom_operators,                  ONLY: atom_int_release,&
                                             atom_int_setup,&
                                             atom_ppint_release,&
                                             atom_ppint_setup,&
                                             atom_relint_release,&
                                             atom_relint_setup
  USE atom_output,                     ONLY: atom_print_basis
  USE atom_types,                      ONLY: &
       GTO_BASIS, atom_basis_type, atom_integrals, atom_p_type, &
       atom_potential_type, atom_type, create_opgrid, opgrid_type, &
       release_opgrid, set_atom
  USE atom_utils,                      ONLY: atom_denmat,&
                                             atom_density,&
                                             integrate_grid
  USE f77_blas
  USE input_constants,                 ONLY: do_analytic
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
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

  PUBLIC :: atom_fit_density, atom_fit_basis, atom_fit_pseudo

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE atom_fit_density (atom,num_gto,iunit,powell_section,error)
    TYPE(atom_type), POINTER                 :: atom
    INTEGER, INTENT(IN)                      :: num_gto, iunit
    TYPE(section_vals_type), POINTER         :: powell_section
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

    ostate%nf = 0
    ostate%nvar = 2
    x(1) = 0.10_dp      !starting point of geometric series
    x(2) = 2.00_dp      !factor of series
    CALL section_vals_val_get(powell_section,"ACCURACY", r_val=ostate%rhoend, error=error)
    CALL section_vals_val_get(powell_section,"STEP_SIZE", r_val=ostate%rhobeg, error=error)
    CALL section_vals_val_get(powell_section,"MAX_FUN", i_val=ostate%maxfun, error=error)
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
        WRITE(iunit,'(" POWELL| Reached",i4,"% of maximal function calls")') &
              INT(REAL(ostate%nf,dp)/REAL(ostate%maxfun,dp)*100._dp)
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
  SUBROUTINE atom_fit_basis (atom_info,basis,pptype,iunit,powell_section,error)
    TYPE(atom_p_type), DIMENSION(:, :), &
      POINTER                                :: atom_info
    TYPE(atom_basis_type), POINTER           :: basis
    LOGICAL, INTENT(IN)                      :: pptype
    INTEGER, INTENT(IN)                      :: iunit
    TYPE(section_vals_type), POINTER         :: powell_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_fit_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, k, l, ll, n10, nr
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: xtob
    LOGICAL                                  :: failure = .FALSE., mult
    REAL(KIND=dp)                            :: al, ear, fopt, rk
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: x
    TYPE(opt_state_type)                     :: ostate

    SELECT CASE (basis%basis_type)
       CASE DEFAULT
         CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       CASE (GTO_BASIS)
         IF ( basis%geometrical ) THEN
           ostate%nvar = 2
           ALLOCATE(x(2),STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
           x(1) = SQRT(basis%aval)
           x(2) = SQRT(basis%cval)
         ELSE
           ll = MAXVAL(basis%nprim(:))
           ALLOCATE(xtob(ll,0:3),STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
           xtob = 0
           ll = SUM(basis%nprim(:))
           ALLOCATE(x(ll),STAT=ierr)
           CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
           x = 0._dp
           ll = 0
           DO l=0,3
             DO i=1,basis%nprim(l)
               mult = .FALSE.
               DO k=1,ll
                 IF ( ABS(basis%am(i,l)-x(k)) < 1.e-6_dp ) THEN
                   mult=.TRUE.
                   xtob(i,l) = k
                 END IF
               END DO
               IF (.NOT. mult) THEN
                 ll = ll + 1
                 x(ll) = basis%am(i,l)
                 xtob(i,l) = ll
               END IF
             END DO
           END DO
           ostate%nvar = ll
           DO i=1,ostate%nvar
             x(i) = SQRT(LOG(1._dp+x(i)))
           END DO
         END IF
    END SELECT

    CALL section_vals_val_get(powell_section,"ACCURACY", r_val=ostate%rhoend, error=error)
    CALL section_vals_val_get(powell_section,"STEP_SIZE", r_val=ostate%rhobeg, error=error)
    CALL section_vals_val_get(powell_section,"MAX_FUN", i_val=ostate%maxfun, error=error)

    ostate%nf = 0
    ostate%iprint = 1
    ostate%unit  = iunit

    ostate%state = 0
    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(/," POWELL| Start optimization procedure")')
      WRITE(iunit,'(/," POWELL| Total number of parameters in optimization",T71,I10)') ostate%nvar
    END IF
    n10 = MAX(ostate%maxfun/100,1)

    fopt = HUGE(0._dp)

    DO

      IF ( ostate%state == 2 ) THEN
         SELECT CASE (basis%basis_type)
            CASE DEFAULT
              CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
            CASE (GTO_BASIS)
              IF ( basis%geometrical ) THEN
                basis%am = 0._dp
                DO l=0,3
                  DO i=1,basis%nbas(l)
                    ll = i - 1 + basis%start(l)
                    basis%am(i,l) = x(1)*x(1) * (x(2)*x(2))**(ll)
                  END DO
                END DO
                basis%aval = x(1)*x(1)
                basis%cval = x(2)*x(2)
              ELSE
                DO l=0,3
                  DO i=1,basis%nprim(l)
                    al = x(xtob(i,l))**2
                    basis%am(i,l) = EXP(al) - 1._dp
                  END DO
                END DO
              END IF
              basis%bf =  0._dp
              basis%dbf = 0._dp
              nr = basis%grid%nr
              DO l=0,3
                DO i=1,basis%nbas(l)
                  al  = basis%am(i,l)
                  DO k=1,nr
                    rk  = basis%grid%rad(k)
                    ear = EXP(-al*basis%grid%rad(k)**2)
                    basis%bf(k,i,l) = rk**l * ear
                    basis%dbf(k,i,l) = ( REAL(l,dp)*rk**(l-1) - 2._dp*al*rk**(l+1) ) * ear
                  END DO
                END DO
              END DO
         END SELECT
         CALL basis_fit (atom_info,basis,pptype,ostate%f,error)
         fopt = MIN(fopt,ostate%f)
      END IF

      IF ( ostate%state == -1 ) EXIT

      CALL powell_optimize (ostate%nvar, x, ostate)

      IF ( ostate%nf == 2 .AND. iunit > 0 ) THEN
        WRITE(iunit,'(" POWELL| Inital value of function",T61,F20.10)') ostate%f
      END IF
      IF ( MOD(ostate%nf,n10) == 0 .AND. iunit > 0 ) THEN
        WRITE(iunit,'(" POWELL| Reached",i4,"% of maximal function calls",T61,F20.10)') &
              INT(REAL(ostate%nf,dp)/REAL(ostate%maxfun,dp)*100._dp), fopt
      END IF

    END DO

    ostate%state = 8
    CALL powell_optimize (ostate%nvar, x, ostate)

    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(" POWELL| Number of function evaluations",T71,I10)') ostate%nf
      WRITE(iunit,'(" POWELL| Final value of function",T61,F20.10)') ostate%f
      CALL atom_print_basis(basis,iunit," Optimized Basis",error)
    END IF

    DEALLOCATE(x,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    IF ( ALLOCATED(xtob) ) THEN
      DEALLOCATE(xtob,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF

  END SUBROUTINE atom_fit_basis
! *****************************************************************************
  SUBROUTINE atom_fit_pseudo (atom_info,atom_refs,ppot,iunit,&
                              pseudo_section,powell_section,error)
    TYPE(atom_p_type), DIMENSION(:, :), &
      POINTER                                :: atom_info, atom_refs
    TYPE(atom_potential_type), POINTER       :: ppot
    INTEGER, INTENT(IN)                      :: iunit
    TYPE(section_vals_type), POINTER         :: pseudo_section, powell_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_fit_pseudo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, n10
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: xtob
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: fopt
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: x
    TYPE(opt_state_type)                     :: ostate

    CALL section_vals_val_get(powell_section,"ACCURACY", r_val=ostate%rhoend, error=error)
    CALL section_vals_val_get(powell_section,"STEP_SIZE", r_val=ostate%rhobeg, error=error)
    CALL section_vals_val_get(powell_section,"MAX_FUN", i_val=ostate%maxfun, error=error)

    ostate%nf = 0
    ostate%nvar = 0
    ALLOCATE(x(ostate%nvar),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    ostate%iprint = 1
    ostate%unit  = iunit

    ostate%state = 0
    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(/," POWELL| Start optimization procedure")')
      WRITE(iunit,'(/," POWELL| Total number of parameters in optimization",T71,I10)') ostate%nvar
    END IF
    n10 = MAX(ostate%maxfun/100,1)

    fopt = HUGE(0._dp)

    DO

      IF ( ostate%state == 2 ) THEN
!deb         CALL pseudo_fit (atom_info,basis,pptype,ostate%f,error)
         fopt = MIN(fopt,ostate%f)
      END IF

      IF ( ostate%state == -1 ) EXIT

      CALL powell_optimize (ostate%nvar, x, ostate)

      IF ( ostate%nf == 2 .AND. iunit > 0 ) THEN
        WRITE(iunit,'(" POWELL| Inital value of function",T61,F20.10)') ostate%f
      END IF
      IF ( MOD(ostate%nf,n10) == 0 .AND. iunit > 0 ) THEN
        WRITE(iunit,'(" POWELL| Reached",i4,"% of maximal function calls",T61,F20.10)') &
              INT(REAL(ostate%nf,dp)/REAL(ostate%maxfun,dp)*100._dp), fopt
      END IF

    END DO

    ostate%state = 8
    CALL powell_optimize (ostate%nvar, x, ostate)

    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(" POWELL| Number of function evaluations",T71,I10)') ostate%nf
      WRITE(iunit,'(" POWELL| Final value of function",T61,F20.10)') ostate%f
    END IF

    DEALLOCATE(x,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    IF ( ALLOCATED(xtob) ) THEN
      DEALLOCATE(xtob,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF

  END SUBROUTINE atom_fit_pseudo
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
  SUBROUTINE basis_fit (atom_info,basis,pptype,afun,error)
    TYPE(atom_p_type), DIMENSION(:, :), &
      POINTER                                :: atom_info
    TYPE(atom_basis_type), POINTER           :: basis
    LOGICAL, INTENT(IN)                      :: pptype
    REAL(dp), INTENT(OUT)                    :: afun
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'basis_fit', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: do_eric, do_erie, ierr, im, &
                                                in, iw, nm, nn, reltyp, zval
    LOGICAL                                  :: eri_c, eri_e, &
                                                failure = .FALSE.
    TYPE(atom_integrals), POINTER            :: atint
    TYPE(atom_potential_type), POINTER       :: pot
    TYPE(atom_type), POINTER                 :: atom

! integrals

    ALLOCATE(atint,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    nn = SIZE(atom_info,1)
    nm = SIZE(atom_info,2)

    ! calculate integrals
    NULLIFY(pot)
    zval = 0
    eri_c = .FALSE.
    eri_e = .FALSE.
    DO in=1,nn
      DO im=1,nm
        atom => atom_info(in,im)%atom
        IF ( pptype .EQV. atom%pp_calc ) THEN
          pot => atom%potential
          zval = atom_info(in,im)%atom%z
          reltyp = atom_info(in,im)%atom%relativistic
          do_eric = atom_info(in,im)%atom%coulomb_integral_type
          do_erie = atom_info(in,im)%atom%exchange_integral_type
          IF(do_eric==do_analytic) eri_c = .TRUE.
          IF(do_erie==do_analytic) eri_e = .TRUE.
          EXIT
        END IF
      END DO
      IF(ASSOCIATED(pot)) EXIT
    END DO
    ! general integrals
    CALL atom_int_setup(atint,basis,potential=pot,eri_coulomb=eri_c,eri_exchange=eri_e,error=error)
    ! potential
    CALL atom_ppint_setup(atint,basis,potential=pot,error=error)
    IF ( pptype ) THEN
      NULLIFY(atint%tzora,atint%hdkh)
    ELSE
      ! relativistic correction terms
      CALL atom_relint_setup(atint,basis,reltyp,zcore=REAL(zval,dp),error=error)
    END IF

    iw = 0
    afun = 0._dp

    DO in=1,nn
      DO im=1,nm
        atom => atom_info(in,im)%atom
        IF ( pptype .EQV. atom%pp_calc ) THEN
          CALL set_atom(atom,basis=basis,error=error)
          CALL set_atom(atom,integrals=atint,error=error)
          CALL calculate_atom(atom,iw,error)

          afun =afun + atom%energy%etot
        END IF
      END DO
    END DO

    CALL atom_int_release(atint,error)
    CALL atom_ppint_release(atint,error)
    CALL atom_relint_release(atint,error)

    DEALLOCATE(atint,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE basis_fit
! *****************************************************************************

END MODULE atom_fit
