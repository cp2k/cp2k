!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
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
  USE atom_output,                     ONLY: atom_print_basis,&
                                             atom_print_basis_file,&
                                             atom_write_pseudo_param
  USE atom_types,                      ONLY: &
       GTO_BASIS, STO_BASIS, atom_basis_type, atom_gthpot_type, &
       atom_integrals, atom_p_type, atom_potential_type, atom_type, &
       create_opgrid, opgrid_type, release_opgrid, set_atom
  USE atom_utils,                      ONLY: &
       atom_consistent_method, atom_denmat, atom_density, &
       atom_orbital_charge, atom_orbital_max, atom_orbital_nodes, atom_wfnr0, &
       get_maxn_occ, integrate_grid
  USE cp_files,                        ONLY: close_file,&
                                             open_file
  USE input_constants,                 ONLY: do_analytic
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE lapack,                          ONLY: lapack_sgesv
  USE mathconstants,                   ONLY: fac,&
                                             fourpi,&
                                             pi
  USE periodic_table,                  ONLY: ptable
  USE physcon,                         ONLY: bohr,&
                                             evolt
  USE powell,                          ONLY: opt_state_type,&
                                             powell_optimize
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_fit'

  PUBLIC :: atom_fit_density, atom_fit_basis, atom_fit_pseudo, atom_fit_kgpot

! *****************************************************************************

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param atom ...
!> \param num_gto ...
!> \param norder ...
!> \param iunit ...
!> \param powell_section ...
!> \param results ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_fit_density (atom,num_gto,norder,iunit,powell_section,results,error)
    TYPE(atom_type), POINTER                 :: atom
    INTEGER, INTENT(IN)                      :: num_gto, norder, iunit
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: powell_section
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: results
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_fit_density', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ierr, n10
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: co
    REAL(KIND=dp), DIMENSION(2)              :: x
    TYPE(opgrid_type), POINTER               :: density
    TYPE(opt_state_type)                     :: ostate

    failure = .FALSE.
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
    IF (norder /= 0) THEN
      density%op = density%op*atom%basis%grid%rad**norder
    END IF

    ostate%nf = 0
    ostate%nvar = 2
    x(1) = 0.10_dp      !starting point of geometric series
    x(2) = 2.00_dp      !factor of series
    IF(PRESENT(powell_section)) THEN
       CALL section_vals_val_get(powell_section,"ACCURACY", r_val=ostate%rhoend, error=error)
       CALL section_vals_val_get(powell_section,"STEP_SIZE", r_val=ostate%rhobeg, error=error)
       CALL section_vals_val_get(powell_section,"MAX_FUN", i_val=ostate%maxfun, error=error)
    ELSE
       ostate%rhoend = 1.e-8_dp
       ostate%rhobeg = 5.e-2_dp
       ostate%maxfun = 1000
    END IF
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
      WRITE(iunit,'(" POWELL| Final value of function",T61,G20.10)') ostate%fopt
      WRITE(iunit,'(" Optimized representation of density using a Geometrical Gaussian basis")')
      WRITE(iunit,'(A,I3,/,T10,A,F10.6,T48,A,F10.6)') " Number of Gaussians:",num_gto,&
             "Starting exponent:",x(1),"Proportionality factor:",x(2)
      WRITE(iunit,'(A)') " Expansion coefficients"
      WRITE(iunit,'(4F20.10)') co(1:num_gto)
    END IF

    IF(PRESENT(results)) THEN
      CPPrecondition(SIZE(results)>=num_gto+2, cp_failure_level, routineP, error, failure)
      results(1) = x(1)
      results(2) = x(2)
      results(3:2+num_gto) = co(1:num_gto)
    END IF

    DEALLOCATE(co,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE atom_fit_density
! *****************************************************************************
!> \brief ...
!> \param atom_info ...
!> \param basis ...
!> \param pptype ...
!> \param iunit ...
!> \param powell_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_fit_basis (atom_info,basis,pptype,iunit,powell_section,&
                             error)
    TYPE(atom_p_type), DIMENSION(:, :), &
      POINTER                                :: atom_info
    TYPE(atom_basis_type), POINTER           :: basis
    LOGICAL, INTENT(IN)                      :: pptype
    INTEGER, INTENT(IN)                      :: iunit
    TYPE(section_vals_type), POINTER         :: powell_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_fit_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, j, k, l, ll, m, n, &
                                                n10, nl, nr
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: xtob
    LOGICAL                                  :: explicit, failure, mult
    REAL(KIND=dp)                            :: al, ear, fopt, pf, rk
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: x
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: wem
    REAL(KIND=dp), DIMENSION(:), POINTER     :: w
    TYPE(opt_state_type)                     :: ostate

    failure = .FALSE.

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
       CASE (STO_BASIS)
         ll = MAXVAL(basis%nbas(:))
         ALLOCATE(xtob(ll,0:3),STAT=ierr)
         CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
         xtob = 0
         ll = SUM(basis%nbas(:))
         ALLOCATE(x(ll),STAT=ierr)
         CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
         x = 0._dp
         ll = 0
         DO l=0,3
           DO i=1,basis%nbas(l)
             ll = ll + 1
             x(ll) = basis%as(i,l)
             xtob(i,l) = ll
           END DO
         END DO
         ostate%nvar = ll
         DO i=1,ostate%nvar
           x(i) = SQRT(LOG(1._dp+x(i)))
         END DO
    END SELECT

    CALL section_vals_val_get(powell_section,"ACCURACY", r_val=ostate%rhoend, error=error)
    CALL section_vals_val_get(powell_section,"STEP_SIZE", r_val=ostate%rhobeg, error=error)
    CALL section_vals_val_get(powell_section,"MAX_FUN", i_val=ostate%maxfun, error=error)

    n=SIZE(atom_info,1)
    m=SIZE(atom_info,2)
    ALLOCATE(wem(n,m),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    wem = 1._dp
    CALL section_vals_val_get(powell_section,"WEIGHT_ELECTRON_CONFIGURATION", explicit=explicit, error=error)
    IF(explicit) THEN
      CALL section_vals_val_get(powell_section,"WEIGHT_ELECTRON_CONFIGURATION", r_vals=w, error=error)
      DO i=1,MIN(SIZE(w),n)
        wem(i,:)=w(i)*wem(i,:)
      END DO
    END IF
    CALL section_vals_val_get(powell_section,"WEIGHT_METHOD", explicit=explicit, error=error)
    IF(explicit) THEN
      CALL section_vals_val_get(powell_section,"WEIGHT_METHOD", r_vals=w, error=error)
      DO i=1,MIN(SIZE(w),m)
        wem(:,i)=w(i)*wem(:,i)
      END DO
    END IF

    DO i=1,SIZE(atom_info,1)
      DO j=1,SIZE(atom_info,2)
        atom_info(i,j)%atom%weight = wem(i,j)
      END DO
    END DO
    DEALLOCATE(wem,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

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
                  al = basis%am(i,l)
                  DO k=1,nr
                    rk  = basis%grid%rad(k)
                    ear = EXP(-al*basis%grid%rad(k)**2)
                    basis%bf(k,i,l) = rk**l * ear
                    basis%dbf(k,i,l) = ( REAL(l,dp)*rk**(l-1) - 2._dp*al*rk**(l+1) ) * ear
                  END DO
                END DO
              END DO
            CASE (STO_BASIS)
              DO l=0,3
                DO i=1,basis%nbas(l)
                  al = x(xtob(i,l))**2
                  basis%as(i,l) = EXP(al) - 1._dp
                END DO
              END DO
              basis%bf =  0._dp
              basis%dbf = 0._dp
              nr = basis%grid%nr
              DO l=0,3
                DO i=1,basis%nbas(l)
                  al  = basis%as(i,l)
                  nl  = basis%ns(i,l)
                  pf  = (2._dp*al)**nl * SQRT(2._dp*al/fac(2*nl))
                  DO k=1,nr
                    rk  = basis%grid%rad(k)
                    ear = EXP(-al*basis%grid%rad(k))
                    basis%bf(k,i,l) = pf * rk**(nl-1) * ear
                    basis%dbf(k,i,l) = pf * ( REAL(nl-1,dp)/rk - al ) * rk**(nl-1) * ear
                  END DO
                END DO
              END DO
         END SELECT
         CALL basis_fit (atom_info,basis,pptype,ostate%f,0,error)
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
       WRITE(iunit,'(" POWELL| Final value of function",T61,F20.10)') ostate%fopt
       ! x->basis
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
          CASE (STO_BASIS)
            DO l=0,3
              DO i=1,basis%nprim(l)
                al = x(xtob(i,l))**2
                basis%as(i,l) = EXP(al) - 1._dp
              END DO
            END DO
            basis%bf =  0._dp
            basis%dbf = 0._dp
            nr = basis%grid%nr
            DO l=0,3
              DO i=1,basis%nbas(l)
                al  = basis%as(i,l)
                nl  = basis%ns(i,l)
                pf  = (2._dp*al)**nl * SQRT(2._dp*al/fac(2*nl))
                DO k=1,nr
                  rk  = basis%grid%rad(k)
                  ear = EXP(-al*basis%grid%rad(k))
                  basis%bf(k,i,l) = pf * rk**(nl-1) * ear
                  basis%dbf(k,i,l) = pf * ( REAL(nl-1,dp)/rk - al ) * rk**(nl-1) * ear
                END DO
              END DO
            END DO
       END SELECT
       CALL atom_print_basis(basis,iunit," Optimized Basis",error)
       CALL atom_print_basis_file(basis,error)
    END IF

    DEALLOCATE(x,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    IF ( ALLOCATED(xtob) ) THEN
      DEALLOCATE(xtob,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF

  END SUBROUTINE atom_fit_basis
! *****************************************************************************
!> \brief ...
!> \param atom_info ...
!> \param atom_refs ...
!> \param ppot ...
!> \param iunit ...
!> \param powell_section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_fit_pseudo (atom_info,atom_refs,ppot,iunit,powell_section,error)
    TYPE(atom_p_type), DIMENSION(:, :), &
      POINTER                                :: atom_info, atom_refs
    TYPE(atom_potential_type), POINTER       :: ppot
    INTEGER, INTENT(IN)                      :: iunit
    TYPE(section_vals_type), POINTER         :: powell_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_fit_pseudo', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, iw, j, k, l, m, n, &
                                                n10, np
    INTEGER, ALLOCATABLE, DIMENSION(:, :)    :: xtob
    INTEGER, DIMENSION(0:3)                  :: ncore
    LOGICAL                                  :: explicit, failure
    REAL(KIND=dp) :: charge, deig, drho, eig, fopt, oc, pchg, peig, pv, rcm, &
      rcov, rmax, t_semi, t_valence, t_virt, w_node, w_psir0, w_semi, &
      w_valence, w_virt
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: x, xi
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: wem
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :, :)               :: pval
    REAL(KIND=dp), DIMENSION(:), POINTER     :: w
    TYPE(atom_type), POINTER                 :: atom
    TYPE(opt_state_type)                     :: ostate

    failure = .FALSE.
! weights for the optimization

    CALL section_vals_val_get(powell_section,"ACCURACY", r_val=ostate%rhoend, error=error)
    CALL section_vals_val_get(powell_section,"STEP_SIZE", r_val=ostate%rhobeg, error=error)
    CALL section_vals_val_get(powell_section,"MAX_FUN", i_val=ostate%maxfun, error=error)

    CALL section_vals_val_get(powell_section,"WEIGHT_POT_VALENCE", r_val=w_valence, error=error)
    CALL section_vals_val_get(powell_section,"WEIGHT_POT_VIRTUAL", r_val=w_virt, error=error)
    CALL section_vals_val_get(powell_section,"WEIGHT_POT_SEMICORE", r_val=w_semi, error=error)
    CALL section_vals_val_get(powell_section,"WEIGHT_POT_NODE", r_val=w_node, error=error)

    CALL section_vals_val_get(powell_section,"WEIGHT_PSIR0", r_val=w_psir0, error=error)
    CALL section_vals_val_get(powell_section,"RCOV_MULTIPLICATION", r_val=rcm, error=error)

    CALL section_vals_val_get(powell_section,"TARGET_POT_VALENCE", r_val=t_valence, error=error)
    CALL section_vals_val_get(powell_section,"TARGET_POT_VIRTUAL", r_val=t_virt, error=error)
    CALL section_vals_val_get(powell_section,"TARGET_POT_SEMICORE", r_val=t_semi, error=error)

    n=SIZE(atom_info,1)
    m=SIZE(atom_info,2)
    ALLOCATE(wem(n,m),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    wem = 1._dp
    ALLOCATE(pval(4,10,0:3,m,n),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    CALL section_vals_val_get(powell_section,"WEIGHT_ELECTRON_CONFIGURATION", explicit=explicit, error=error)
    IF(explicit) THEN
      CALL section_vals_val_get(powell_section,"WEIGHT_ELECTRON_CONFIGURATION", r_vals=w, error=error)
      DO i=1,MIN(SIZE(w),n)
        wem(i,:)=w(i)*wem(i,:)
      END DO
    END IF
    CALL section_vals_val_get(powell_section,"WEIGHT_METHOD", explicit=explicit, error=error)
    IF(explicit) THEN
      CALL section_vals_val_get(powell_section,"WEIGHT_METHOD", r_vals=w, error=error)
      DO i=1,MIN(SIZE(w),m)
        wem(:,i)=w(i)*wem(:,i)
      END DO
    END IF

    CALL open_file(file_name="POWELL_RESULT",file_status="UNKNOWN",file_action="WRITE",unit_number=iw)

    ALLOCATE(xi(200),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    !decide here what to optimize
    CALL get_pseudo_param(xi,ostate%nvar,ppot%gth_pot,error)
    ALLOCATE(x(ostate%nvar),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    x(1:ostate%nvar) = xi(1:ostate%nvar)
    DEALLOCATE(xi,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    ostate%nf = 0
    ostate%iprint = 1
    ostate%unit  = iunit

    ostate%state = 0
    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(/," POWELL| Start optimization procedure")')
      WRITE(iunit,'(/," POWELL| Total number of parameters in optimization",T71,I10)') ostate%nvar
    END IF
    n10 = MAX(ostate%maxfun/100,1)

    rcov = ptable(atom_refs(1,1)%atom%z)%covalent_radius*bohr * rcm
    !set all reference values
    DO i=1,SIZE(atom_info,1)
      DO j=1,SIZE(atom_info,2)
        atom => atom_info(i,j)%atom
        IF(atom_consistent_method(atom%method_type,atom%state%multiplicity)) THEN
          IF(atom%state%multiplicity==-1) THEN
            ! no spin polarization
            atom%weight = wem(i,j)
            ncore = get_maxn_occ(atom_info(i,j)%atom%state%core)
            DO l=0,3
              np = atom%state%maxn_calc(l)
              DO k=1,np
                CALL atom_orbital_max(rmax,atom_refs(i,j)%atom%orbitals%wfn(:,ncore(l)+k,l),&
                                      rcov,l,atom_refs(i,j)%atom%basis,error)
                atom%orbitals%rcmax(k,l,1) = MAX(rcov,rmax)
                CALL atom_orbital_charge(charge,atom_refs(i,j)%atom%orbitals%wfn(:,ncore(l)+k,l),&
                                         atom%orbitals%rcmax(k,l,1),l,atom_refs(i,j)%atom%basis,error)
                atom%orbitals%refene(k,l,1) = atom_refs(i,j)%atom%orbitals%ener(ncore(l)+k,l)
                atom%orbitals%refchg(k,l,1) = charge
                IF ( k > atom%state%maxn_occ(l) ) THEN
                  IF ( k <= atom%state%maxn_occ(l)+1 ) THEN
                    atom%orbitals%wrefene(k,l,1) = w_virt
                    atom%orbitals%wrefchg(k,l,1) = w_virt/100._dp
                  ELSE
                    atom%orbitals%wrefene(k,l,1) = 0._dp
                    atom%orbitals%wrefchg(k,l,1) = 0._dp
                  END IF
                  atom%orbitals%crefene(k,l,1) = t_virt
                ELSEIF ( k < atom%state%maxn_occ(l) ) THEN
                  atom%orbitals%wrefene(k,l,1) = w_semi
                  atom%orbitals%wrefchg(k,l,1) = w_semi/100._dp
                  atom%orbitals%crefene(k,l,1) = t_semi
                ELSE
                  atom%orbitals%wrefene(k,l,1) = w_valence
                  atom%orbitals%wrefchg(k,l,1) = w_valence/100._dp
                  atom%orbitals%crefene(k,l,1) = t_valence
                  IF ( l==0 ) THEN
                    atom%orbitals%wpsir0(k,1) = w_psir0
                  END IF
                END IF
              END DO
              DO k=1,np
                atom%orbitals%refnod(k,l,1) = REAL(k-1,KIND=dp)
                ! we only enforce 0-nodes for the first state
                IF ( k==1 .AND. atom%state%occupation(l,k)/=0._dp ) atom%orbitals%wrefnod(k,l,1) = w_node
              END DO
            END DO
          ELSE
            ! spin polarization
            atom%weight = wem(i,j)
            ncore = get_maxn_occ(atom_info(i,j)%atom%state%core)
            DO l=0,3
              np = atom%state%maxn_calc(l)
              DO k=1,np
                CALL atom_orbital_max(rmax,atom_refs(i,j)%atom%orbitals%wfna(:,ncore(l)+k,l),&
                                      rcov,l,atom_refs(i,j)%atom%basis,error)
                atom%orbitals%rcmax(k,l,1) = MAX(rcov,rmax)
                CALL atom_orbital_max(rmax,atom_refs(i,j)%atom%orbitals%wfnb(:,ncore(l)+k,l),&
                                      rcov,l,atom_refs(i,j)%atom%basis,error)
                atom%orbitals%rcmax(k,l,2) = MAX(rcov,rmax)
                CALL atom_orbital_charge(charge,atom_refs(i,j)%atom%orbitals%wfna(:,ncore(l)+k,l),&
                                         atom%orbitals%rcmax(k,l,1),l,atom_refs(i,j)%atom%basis,error)
                atom%orbitals%refene(k,l,1) = atom_refs(i,j)%atom%orbitals%enera(ncore(l)+k,l)
                atom%orbitals%refchg(k,l,1) = charge
                CALL atom_orbital_charge(charge,atom_refs(i,j)%atom%orbitals%wfnb(:,ncore(l)+k,l),&
                                         atom%orbitals%rcmax(k,l,1),l,atom_refs(i,j)%atom%basis,error)
                atom%orbitals%refene(k,l,2) = atom_refs(i,j)%atom%orbitals%enerb(ncore(l)+k,l)
                atom%orbitals%refchg(k,l,2) = charge
                ! the following assignments could be further specialized
                IF ( k > atom%state%maxn_occ(l) ) THEN
                  IF ( k <= atom%state%maxn_occ(l)+1 ) THEN
                    atom%orbitals%wrefene(k,l,1:2) = w_virt
                    atom%orbitals%wrefchg(k,l,1:2) = w_virt/100._dp
                  ELSE
                    atom%orbitals%wrefene(k,l,1:2) = 0._dp
                    atom%orbitals%wrefchg(k,l,1:2) = 0._dp
                  END IF
                  atom%orbitals%crefene(k,l,1:2) = t_virt
                ELSEIF ( k < atom%state%maxn_occ(l) ) THEN
                  atom%orbitals%wrefene(k,l,1:2) = w_semi
                  atom%orbitals%wrefchg(k,l,1:2) = w_semi/100._dp
                  atom%orbitals%crefene(k,l,1:2) = t_semi
                ELSE
                  atom%orbitals%wrefene(k,l,1:2) = w_valence
                  atom%orbitals%wrefchg(k,l,1:2) = w_valence/100._dp
                  atom%orbitals%crefene(k,l,1:2) = t_valence
                  IF ( l==0 ) THEN
                    atom%orbitals%wpsir0(k,1:2) = w_psir0
                  END IF
                END IF
              END DO
              DO k=1,np
                atom%orbitals%refnod(k,l,1:2) = REAL(k-1,KIND=dp)
                ! we only enforce 0-nodes for the first state
                IF ( k==1 .AND. atom%state%occa(l,k)/=0._dp ) atom%orbitals%wrefnod(k,l,1) = w_node
                IF ( k==1 .AND. atom%state%occb(l,k)/=0._dp ) atom%orbitals%wrefnod(k,l,2) = w_node
              END DO
            END DO
          END IF
          CALL calculate_atom(atom,0,error=error)
        END IF
      END DO
    END DO

    DEALLOCATE(wem,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    WRITE(iunit,'(/," POWELL| Initial errors of target values")')
    CALL put_pseudo_param(x,ppot%gth_pot,error)
    CALL pseudo_fit(atom_info,ppot,ostate%f,pval,error)
    DO i=1,SIZE(atom_info,1)
      DO j=1,SIZE(atom_info,2)
        atom => atom_info(i,j)%atom
        IF(atom_consistent_method(atom%method_type,atom%state%multiplicity)) THEN
          WRITE(iunit,'(/," Reference configuration  ",T31,i5,T50," Method number ",T76,i5)') i,j
          IF(atom%state%multiplicity==-1) THEN
            ! no spin polarization
            WRITE(iunit,'("    L    N    Occupation      Eigenvalue [eV]           dE [eV]          dCharge ")')
            DO l=0,3
              np = atom%state%maxn_calc(l)
              IF (np > 0) THEN
                DO k=1,np
                  oc = atom%state%occupation(l,k)
                  eig = atom%orbitals%ener(k,l)
                  deig = eig - atom%orbitals%refene(k,l,1)
                  peig = pval(1,k,l,j,i)/ostate%f * 100._dp
                  CALL atom_orbital_charge(charge,atom%orbitals%wfn(:,k,l),atom%orbitals%rcmax(k,l,1),l,atom%basis,error)
                  drho = charge - atom%orbitals%refchg(k,l,1)
                  pchg = pval(2,k,l,j,i)/ostate%f * 100._dp
                  WRITE(iunit,'(I5,I5,F14.2,F21.10,F14.6,"[",I2,"]",F13.6,"[",I2,"]")') &
                              l,k,oc,eig*evolt,deig*evolt,NINT(peig),drho,NINT(pchg)
                END DO
              END IF
            END DO
          ELSE
            ! spin polarization
            WRITE(iunit,'("    L    N   Spin Occupation    Eigenvalue [eV]          dE [eV]         dCharge ")')
            DO l=0,3
              np = atom%state%maxn_calc(l)
              IF (np > 0) THEN
                DO k=1,np
                  oc = atom%state%occa(l,k)
                  eig = atom%orbitals%enera(k,l)
                  deig = eig - atom%orbitals%refene(k,l,1)
                  peig = pval(1,k,l,j,i)/ostate%f * 100._dp
                  CALL atom_orbital_charge(charge,atom%orbitals%wfna(:,k,l),atom%orbitals%rcmax(k,l,1),l,atom%basis,error)
                  drho = charge - atom%orbitals%refchg(k,l,1)
                  pchg = pval(2,k,l,j,i)/ostate%f * 100._dp
                  WRITE(iunit,'(I5,I5,2X,A5,F11.2,F19.10,F13.6,"[",I2,"]",F12.6,"[",I2,"]")') &
                              l,k,"alpha",oc,eig*evolt,deig*evolt,NINT(peig),drho,NINT(pchg)
                  oc = atom%state%occb(l,k)
                  eig = atom%orbitals%enerb(k,l)
                  deig = eig - atom%orbitals%refene(k,l,2)
                  peig = pval(3,k,l,j,i)/ostate%f * 100._dp
                  CALL atom_orbital_charge(charge,atom%orbitals%wfnb(:,k,l),atom%orbitals%rcmax(k,l,2),l,atom%basis,error)
                  drho = charge - atom%orbitals%refchg(k,l,2)
                  pchg = pval(4,k,l,j,i)/ostate%f * 100._dp
                  WRITE(iunit,'(I5,I5,2X,A5,F11.2,F19.10,F13.6,"[",I2,"]",F12.6,"[",I2,"]")') &
                              l,k," beta",oc,eig*evolt,deig*evolt,NINT(peig),drho,NINT(pchg)
                END DO
              END IF
            END DO
          END IF
        END IF
      END DO
      WRITE(iunit,*)
    END DO

    fopt = HUGE(0._dp)

    DO

      IF ( ostate%state == 2 ) THEN
         CALL put_pseudo_param(x,ppot%gth_pot,error)
         CALL pseudo_fit (atom_info,ppot,ostate%f,pval,error)
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
        CALL put_pseudo_param(ostate%xopt,ppot%gth_pot,error)
        CALL atom_write_pseudo_param(ppot%gth_pot,error=error)
      END IF

      WRITE(iw,*) ostate%nf,ostate%f,x(1:ostate%nvar)

    END DO

    ostate%state = 8
    CALL powell_optimize (ostate%nvar, x, ostate)
    CALL put_pseudo_param(x,ppot%gth_pot,error)
    CALL atom_write_pseudo_param(ppot%gth_pot,error=error)

    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(" POWELL| Number of function evaluations",T71,I10)') ostate%nf
      WRITE(iunit,'(" POWELL| Final value of function",T61,F20.10)') ostate%fopt

      CALL put_pseudo_param(x,ppot%gth_pot,error)
      CALL pseudo_fit(atom_info,ppot,ostate%f,pval,error)

      WRITE(iunit,'(/," POWELL| Final errors of target values")')
      DO i=1,SIZE(atom_info,1)
        DO j=1,SIZE(atom_info,2)
          atom => atom_info(i,j)%atom
          IF(atom_consistent_method(atom%method_type,atom%state%multiplicity)) THEN
            WRITE(iunit,'(/," Reference configuration  ",T31,i5,T50," Method number ",T76,i5)') i,j
            IF(atom%state%multiplicity==-1) THEN
              !no spin polarization
              WRITE(iunit,'("    L    N    Occupation      Eigenvalue [eV]           dE [eV]          dCharge ")')
              DO l=0,3
                np = atom%state%maxn_calc(l)
                IF (np > 0) THEN
                  DO k=1,np
                    oc = atom%state%occupation(l,k)
                    eig = atom%orbitals%ener(k,l)
                    deig = eig - atom%orbitals%refene(k,l,1)
                    peig = pval(1,k,l,j,i)/ostate%f * 100._dp
                    CALL atom_orbital_charge(charge,atom%orbitals%wfn(:,k,l),atom%orbitals%rcmax(k,l,1),l,atom%basis,error)
                    drho = charge - atom%orbitals%refchg(k,l,1)
                    pchg = pval(2,k,l,j,i)/ostate%f * 100._dp
                    WRITE(iunit,'(I5,I5,F14.2,F21.10,F14.6,"[",I2,"]",F13.6,"[",I2,"]")') &
                                l,k,oc,eig*evolt,deig*evolt,NINT(peig),drho,NINT(pchg)
                  END DO
                END IF
              END DO
              np = atom%state%maxn_calc(0)
              DO k=1,np
                 CALL atom_wfnr0(pv,atom%orbitals%wfn(:,k,0),atom%basis,error)
                 pchg = atom%weight*atom%orbitals%wpsir0(k,1)*pv*pv/ostate%f * 100._dp
                 WRITE(iunit,'("    s-states"," N=",I5,T40,"Wavefunction at r=0:",T64,F13.6,"[",I2,"]")') k,pv,NINT(pchg)
              END DO
            ELSE
              !spin polarization
              WRITE(iunit,'("    L    N   Spin Occupation     Eigenvalue [eV]          dE [eV]        dCharge ")')
              DO l=0,3
                np = atom%state%maxn_calc(l)
                IF (np > 0) THEN
                  DO k=1,np
                    oc = atom%state%occa(l,k)
                    eig = atom%orbitals%enera(k,l)
                    deig = eig - atom%orbitals%refene(k,l,1)
                    peig = pval(1,k,l,j,i)/ostate%f * 100._dp
                    CALL atom_orbital_charge(charge,atom%orbitals%wfna(:,k,l),atom%orbitals%rcmax(k,l,1),l,atom%basis,error)
                    drho = charge - atom%orbitals%refchg(k,l,1)
                    pchg = pval(2,k,l,j,i)/ostate%f * 100._dp
                    WRITE(iunit,'(I5,I5,A,F11.2,F20.10,F13.6,"[",I2,"]",F11.6,"[",I2,"]")') &
                                l,k,"  alpha",oc,eig*evolt,deig*evolt,NINT(peig),drho,NINT(pchg)
                    oc = atom%state%occb(l,k)
                    eig = atom%orbitals%enerb(k,l)
                    deig = eig - atom%orbitals%refene(k,l,2)
                    peig = pval(3,k,l,j,i)/ostate%f * 100._dp
                    CALL atom_orbital_charge(charge,atom%orbitals%wfnb(:,k,l),atom%orbitals%rcmax(k,l,2),l,atom%basis,error)
                    drho = charge - atom%orbitals%refchg(k,l,2)
                    pchg = pval(4,k,l,j,i)/ostate%f * 100._dp
                    WRITE(iunit,'(I5,I5,A,F11.2,F20.10,F13.6,"[",I2,"]",F11.6,"[",I2,"]")') &
                                l,k,"   beta",oc,eig*evolt,deig*evolt,NINT(peig),drho,NINT(pchg)
                  END DO
                END IF
              END DO
              np = atom%state%maxn_calc(0)
              DO k=1,np
                 CALL atom_wfnr0(pv,atom%orbitals%wfna(:,k,0),atom%basis,error)
                 pchg = atom%weight*atom%orbitals%wpsir0(k,1)*pv*pv/ostate%f * 100._dp
                 WRITE(iunit,'("    s-states"," N=",I5,T35,"Alpha Wavefunction at r=0:",T64,F13.6,"[",I2,"]")') k,pv,NINT(pchg)
                 CALL atom_wfnr0(pv,atom%orbitals%wfnb(:,k,0),atom%basis,error)
                 pchg = atom%weight*atom%orbitals%wpsir0(k,2)*pv*pv/ostate%f * 100._dp
                 WRITE(iunit,'("    s-states"," N=",I5,T36,"Beta Wavefunction at r=0:",T64,F13.6,"[",I2,"]")') k,pv,NINT(pchg)
              END DO
            END IF
          END IF
        END DO
      END DO
    END IF

    DEALLOCATE(x,pval,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    IF ( ALLOCATED(xtob) ) THEN
      DEALLOCATE(xtob,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    END IF

    CALL close_file(unit_number=iw)

  END SUBROUTINE atom_fit_pseudo
! *****************************************************************************
!> \brief ...
!> \param density ...
!> \param atom ...
!> \param n ...
!> \param aval ...
!> \param cval ...
!> \param co ...
!> \param aerr ...
!> \param error ...
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
    LOGICAL                                  :: failure
    REAL(dp)                                 :: a, rk, zval
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: den, pe, uval
    REAL(dp), ALLOCATABLE, DIMENSION(:, :)   :: bf, smat, tval

    failure = .FALSE.

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
    den(:) = 0._dp
    DO i=1,n
      den(:) = den(:) + co(i)*bf(:,i)
    END DO
    den(:) = den(:) * fourpi
    den(:) = (den(:)-density%op(:))**2
    aerr = SQRT(integrate_grid(den,atom%basis%grid))

    DEALLOCATE(pe,bf,den,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    DEALLOCATE(tval,uval,smat,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE density_fit
! *****************************************************************************
!> \brief ...
!> \param atom_info ...
!> \param basis ...
!> \param pptype ...
!> \param afun ...
!> \param iw ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE basis_fit (atom_info,basis,pptype,afun,iw,error)
    TYPE(atom_p_type), DIMENSION(:, :), &
      POINTER                                :: atom_info
    TYPE(atom_basis_type), POINTER           :: basis
    LOGICAL, INTENT(IN)                      :: pptype
    REAL(dp), INTENT(OUT)                    :: afun
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'basis_fit', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: do_eric, do_erie, ierr, im, &
                                                in, nm, nn, reltyp, zval
    LOGICAL                                  :: eri_c, eri_e, failure
    TYPE(atom_integrals), POINTER            :: atint
    TYPE(atom_potential_type), POINTER       :: pot
    TYPE(atom_type), POINTER                 :: atom

    failure = .FALSE.

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

    afun = 0._dp

    DO in=1,nn
      DO im=1,nm
        atom => atom_info(in,im)%atom
        IF(atom_consistent_method(atom%method_type,atom%state%multiplicity)) THEN
          IF ( pptype .EQV. atom%pp_calc ) THEN
            CALL set_atom(atom,basis=basis,error=error)
            CALL set_atom(atom,integrals=atint,error=error)
            CALL calculate_atom(atom,iw,error=error)
            afun =afun + atom%energy%etot*atom%weight
          END IF
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
!> \brief ...
!> \param atom_info ...
!> \param ppot ...
!> \param afun ...
!> \param pval ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE pseudo_fit (atom_info,ppot,afun,pval,error)
    TYPE(atom_p_type), DIMENSION(:, :), &
      POINTER                                :: atom_info
    TYPE(atom_potential_type), POINTER       :: ppot
    REAL(dp), INTENT(OUT)                    :: afun
    REAL(dp), DIMENSION(:, :, 0:, :, :), &
      INTENT(OUT)                            :: pval
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pseudo_fit', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, j, k, l, n, node
    REAL(KIND=dp)                            :: charge, de, fde, pv, rcov, &
                                                rcov1, rcov2, tv
    TYPE(atom_integrals), POINTER            :: pp_int
    TYPE(atom_type), POINTER                 :: atom

    afun = 0._dp
    pval = 0._dp

    pp_int => atom_info(1,1)%atom%integrals
    CALL atom_ppint_release(pp_int,error)
    CALL atom_ppint_setup(pp_int,atom_info(1,1)%atom%basis,potential=ppot,error=error)

    DO i=1,SIZE(atom_info,1)
      DO j=1,SIZE(atom_info,2)
        atom => atom_info(i,j)%atom
        IF(atom_consistent_method(atom%method_type,atom%state%multiplicity)) THEN
          CALL set_atom(atom,integrals=pp_int,potential=ppot,error=error)
          CALL calculate_atom(atom,0,noguess=.TRUE.,error=error)
          DO l=0,atom%state%maxl_calc
            n = atom%state%maxn_calc(l)
            DO k=1,n
              IF(atom%state%multiplicity==-1) THEN
                !no spin polarization
                rcov = atom%orbitals%rcmax(k,l,1)
                tv = atom%orbitals%crefene(k,l,1)
                de = ABS(atom%orbitals%ener(k,l)-atom%orbitals%refene(k,l,1))
                fde = 1._dp/(1._dp+EXP(10._dp*(1._dp-de/tv)))
                pv = atom%weight*atom%orbitals%wrefene(k,l,1)*de*de*fde
                afun = afun + pv
                pval(1,k,l,j,i) = pv
                IF(atom%orbitals%wrefchg(k,l,1) > 0._dp) THEN
                  CALL atom_orbital_charge(charge,atom%orbitals%wfn(:,k,l),rcov,l,atom%basis,error)
                  pv = 0.01_dp*atom%weight*atom%orbitals%wrefchg(k,l,1)*(charge-atom%orbitals%refchg(k,l,1))**2
                  afun = afun + pv
                  pval(2,k,l,j,i) = pv
                END IF
                IF(atom%orbitals%wrefnod(k,l,1) > 0._dp) THEN
                  CALL atom_orbital_nodes(node,atom%orbitals%wfn(:,k,l),2._dp*rcov,l,atom%basis,error)
                  afun = afun + atom%weight*atom%orbitals%wrefnod(k,l,1)*ABS(REAL(node,dp)-atom%orbitals%refnod(k,l,1))
                END IF
                IF ( l==0 ) THEN
                  IF(atom%orbitals%wpsir0(k,1) > 0._dp) THEN
                    CALL atom_wfnr0(pv,atom%orbitals%wfn(:,k,l),atom%basis,error)
                    pv = atom%weight*atom%orbitals%wpsir0(k,1)*pv*pv
                    afun = afun + pv
                  END IF
                END IF
              ELSE
                !spin polarization
                rcov1 = atom%orbitals%rcmax(k,l,1)
                rcov2 = atom%orbitals%rcmax(k,l,2)
                tv = atom%orbitals%crefene(k,l,1)
                de = ABS(atom%orbitals%enera(k,l)-atom%orbitals%refene(k,l,1))
                fde = 1._dp/(1._dp+EXP(10._dp*(1._dp-de/tv)))
                pv = atom%weight*atom%orbitals%wrefene(k,l,1)*de*de*fde
                afun = afun + pv
                pval(1,k,l,j,i) = pv
                tv = atom%orbitals%crefene(k,l,2)
                de = ABS(atom%orbitals%enerb(k,l)-atom%orbitals%refene(k,l,2))
                fde = 1._dp/(1._dp+EXP(10._dp*(1._dp-de/tv)))
                pv = atom%weight*atom%orbitals%wrefene(k,l,2)*de*de*fde
                afun = afun + pv
                pval(3,k,l,j,i) = pv
                IF(atom%orbitals%wrefchg(k,l,1) > 0._dp) THEN
                  CALL atom_orbital_charge(charge,atom%orbitals%wfna(:,k,l),rcov1,l,atom%basis,error)
                  pv = 0.01_dp*atom%weight*atom%orbitals%wrefchg(k,l,1)*(charge-atom%orbitals%refchg(k,l,1))**2
                  afun = afun + pv
                  pval(2,k,l,j,i) = pv
                END IF
                IF(atom%orbitals%wrefchg(k,l,2) > 0._dp) THEN
                  CALL atom_orbital_charge(charge,atom%orbitals%wfnb(:,k,l),rcov2,l,atom%basis,error)
                  pv = 0.01_dp*atom%weight*atom%orbitals%wrefchg(k,l,2)*(charge-atom%orbitals%refchg(k,l,2))**2
                  afun = afun + pv
                  pval(4,k,l,j,i) = pv
                END IF
                IF(atom%orbitals%wrefnod(k,l,1) > 0._dp) THEN
                  CALL atom_orbital_nodes(node,atom%orbitals%wfna(:,k,l),2._dp*rcov1,l,atom%basis,error)
                  afun = afun + atom%weight*atom%orbitals%wrefnod(k,l,1)*ABS(REAL(node,dp)-atom%orbitals%refnod(k,l,1))
                END IF
                IF(atom%orbitals%wrefnod(k,l,2) > 0._dp) THEN
                  CALL atom_orbital_nodes(node,atom%orbitals%wfnb(:,k,l),2._dp*rcov2,l,atom%basis,error)
                  afun = afun + atom%weight*atom%orbitals%wrefnod(k,l,2)*ABS(REAL(node,dp)-atom%orbitals%refnod(k,l,2))
                END IF
                IF ( l==0 ) THEN
                  IF(atom%orbitals%wpsir0(k,1) > 0._dp) THEN
                    CALL atom_wfnr0(pv,atom%orbitals%wfna(:,k,l),atom%basis,error)
                    pv = atom%weight*atom%orbitals%wpsir0(k,1)*pv*pv
                    afun = afun + pv
                  END IF
                  IF(atom%orbitals%wpsir0(k,2) > 0._dp) THEN
                    CALL atom_wfnr0(pv,atom%orbitals%wfnb(:,k,l),atom%basis,error)
                    pv = atom%weight*atom%orbitals%wpsir0(k,2)*pv*pv
                    afun = afun + pv
                  END IF
                END IF
              ENDIF
            END DO
          END DO
        END IF
      END DO
    END DO

  END SUBROUTINE pseudo_fit
! *****************************************************************************
!> \brief ...
!> \param pvec ...
!> \param nval ...
!> \param gthpot ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE get_pseudo_param (pvec,nval,gthpot,error)
    REAL(KIND=dp), DIMENSION(:), INTENT(out) :: pvec
    INTEGER, INTENT(out)                     :: nval
    TYPE(atom_gthpot_type), INTENT(in)       :: gthpot
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'get_pseudo_param', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ival, j, l, n

    IF(gthpot%lsdpot) THEN
      pvec = 0
      ival = 0
      DO j=1,gthpot%nexp_lsd
        ival=ival+1
        pvec(ival) = rcpro(-1,gthpot%alpha_lsd(j))
        DO i=1,gthpot%nct_lsd(j)
          ival=ival+1
          pvec(ival) = gthpot%cval_lsd(i,j)
        END DO
      END DO
    ELSE
      pvec = 0
      ival = 1
      pvec(ival) = rcpro(-1,gthpot%rc)
      DO i=1,gthpot%ncl
        ival=ival+1
        pvec(ival) = gthpot%cl(i)
      END DO
      IF(gthpot%lpotextended) THEN
        DO j=1,gthpot%nexp_lpot
          ival=ival+1
          pvec(ival) = rcpro(-1,gthpot%alpha_lpot(j))
          DO i=1,gthpot%nct_lpot(j)
            ival=ival+1
            pvec(ival) = gthpot%cval_lpot(i,j)
          END DO
        END DO
      END IF
      IF(gthpot%nlcc) THEN
        DO j=1,gthpot%nexp_nlcc
          ival=ival+1
          pvec(ival) = rcpro(-1,gthpot%alpha_nlcc(j))
          DO i=1,gthpot%nct_nlcc(j)
            ival=ival+1
            pvec(ival) = gthpot%cval_nlcc(i,j)
          END DO
        END DO
      END IF
      DO l=0,3
        IF(gthpot%nl(l) > 0) THEN
          ival = ival+1
          pvec(ival) = rcpro(-1,gthpot%rcnl(l))
        END IF
      END DO
      DO l=0,3
        IF(gthpot%nl(l) > 0) THEN
          n = gthpot%nl(l)
          DO i=1,n
            DO j=i,n
              ival = ival+1
              pvec(ival) = gthpot%hnl(i,j,l)
            END DO
          END DO
        END IF
      END DO
    END IF
    nval = ival

  END SUBROUTINE get_pseudo_param

! *****************************************************************************
!> \brief ...
!> \param pvec ...
!> \param gthpot ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE put_pseudo_param (pvec,gthpot,error)
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: pvec
    TYPE(atom_gthpot_type), INTENT(inout)    :: gthpot
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'put_pseudo_param', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ival, j, l, n

    IF(gthpot%lsdpot) THEN
      ival=0
      DO j=1,gthpot%nexp_lsd
        ival=ival+1
        gthpot%alpha_lsd(j)=rcpro(1,pvec(ival))
        DO i=1,gthpot%nct_lsd(j)
          ival=ival+1
          gthpot%cval_lsd(i,j) = pvec(ival)
        END DO
      END DO
    ELSE
      ival = 1
      gthpot%rc = rcpro(1,pvec(ival))
      DO i=1,gthpot%ncl
        ival=ival+1
        gthpot%cl(i) = pvec(ival)
      END DO
      IF(gthpot%lpotextended) THEN
        DO j=1,gthpot%nexp_lpot
          ival=ival+1
          gthpot%alpha_lpot(j)=rcpro(1,pvec(ival))
          DO i=1,gthpot%nct_lpot(j)
            ival=ival+1
            gthpot%cval_lpot(i,j) = pvec(ival)
          END DO
        END DO
      END IF
      IF(gthpot%nlcc) THEN
        DO j=1,gthpot%nexp_nlcc
          ival=ival+1
          gthpot%alpha_nlcc(j)=rcpro(1,pvec(ival))
          DO i=1,gthpot%nct_nlcc(j)
            ival=ival+1
            gthpot%cval_nlcc(i,j) = pvec(ival)
          END DO
        END DO
      END IF
      DO l=0,3
        IF(gthpot%nl(l) > 0) THEN
          ival = ival+1
          gthpot%rcnl(l) = rcpro(1,pvec(ival))
        END IF
      END DO
      DO l=0,3
        IF(gthpot%nl(l) > 0) THEN
          n = gthpot%nl(l)
          DO i=1,n
            DO j=i,n
              ival = ival+1
              gthpot%hnl(i,j,l) = pvec(ival)
            END DO
          END DO
        END IF
      END DO
   END IF

  END SUBROUTINE put_pseudo_param

! *****************************************************************************
!> \brief ...
!> \param id ...
!> \param xval ...
!> \retval yval ...
! *****************************************************************************
  FUNCTION rcpro(id,xval) RESULT(yval)
    INTEGER, INTENT(IN)                      :: id
    REAL(dp), INTENT(IN)                     :: xval
    REAL(dp)                                 :: yval

    REAL(dp)                                 :: x1, x2

    IF ( id==1 ) THEN
      yval = TANH(0.1_dp*xval)**2
    ELSE IF ( id==-1 ) THEN
      x1 = SQRT(xval)
      IF ( x1 > 1._dp ) STOP "rcpro"
      x2 = 0.5_dp * LOG((1._dp+x1)/(1._dp-x1))
      yval = x2/0.1_dp
    ELSE
      STOP "rcpro"
    END IF
  END FUNCTION rcpro

! *****************************************************************************
!> \brief ...
!> \param atom ...
!> \param num_gau ...
!> \param num_pol ...
!> \param iunit ...
!> \param powell_section ...
!> \param results ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE atom_fit_kgpot (atom,num_gau,num_pol,iunit,powell_section,results,error)
    TYPE(atom_type), POINTER                 :: atom
    INTEGER, INTENT(IN)                      :: num_gau, num_pol, iunit
    TYPE(section_vals_type), OPTIONAL, &
      POINTER                                :: powell_section
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: results
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_fit_kgpot', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: t23 = 2._dp/3._dp

    INTEGER                                  :: i, ierr, ig, ip, iw, j, n10
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: x
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: co
    TYPE(opgrid_type), POINTER               :: density
    TYPE(opt_state_type)                     :: ostate

    failure = .FALSE.
    ! at least one parameter to be optimized
    CPPrecondition(num_pol*num_gau > 0, cp_failure_level, routineP, error, failure)

    ALLOCATE(co(num_pol+1,num_gau),x(num_pol*num_gau+num_gau),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    co = 0._dp

    ! calculate density
    NULLIFY(density)
    CALL create_opgrid(density,atom%basis%grid,error)
    CALL atom_denmat(atom%orbitals%pmat,atom%orbitals%wfn,atom%basis%nbas,atom%state%occupation,&
                     atom%state%maxl_occ,atom%state%maxn_occ,error)
    CALL atom_density(density%op,atom%orbitals%pmat,atom%basis,atom%state%maxl_occ,typ="RHO",error=error)
    ! target functional 
    density%op = t23*(0.3_dp*(3.0_dp*pi*pi)**t23)*density%op**t23

    ! initiallize parameter
    ostate%nf = 0
    ostate%nvar = num_pol*num_gau + num_gau
    DO i=1,num_gau
       co(1,i) = 0.5_dp + REAL(i-1,KIND=dp)
       co(2,i) = 1.0_dp
       DO j=3,num_pol+1
          co(j,i) = 0.1_dp
       END DO
    END DO
   CALL putvar(x,co,num_pol,num_gau)

    IF(PRESENT(powell_section)) THEN
       CALL section_vals_val_get(powell_section,"ACCURACY", r_val=ostate%rhoend, error=error)
       CALL section_vals_val_get(powell_section,"STEP_SIZE", r_val=ostate%rhobeg, error=error)
       CALL section_vals_val_get(powell_section,"MAX_FUN", i_val=ostate%maxfun, error=error)
    ELSE
       ostate%rhoend = 1.e-8_dp
       ostate%rhobeg = 5.e-2_dp
       ostate%maxfun = 1000
    END IF
    ostate%iprint = 1
    ostate%unit  = iunit

    ostate%state = 0
    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(/," ",13("*")," Approximated Nonadditive Kinetic Energy Functional ",14("*"))')
      WRITE(iunit,'(" POWELL| Start optimization procedure")')
    END IF
    n10 = 50

    DO

      IF ( ostate%state == 2 ) THEN
         CALL getvar(x,co,num_pol,num_gau)
         CALL kgpot_fit (density,num_gau,num_pol,co,ostate%f,error)
      END IF

      IF ( ostate%state == -1 ) EXIT

      CALL powell_optimize (ostate%nvar, x, ostate)

      IF ( MOD(ostate%nf,n10) == 0 .AND. iunit > 0 ) THEN
        WRITE(iunit,'(" POWELL| Reached",i4,"% of maximal function calls",T66,G15.7)') &
              INT(REAL(ostate%nf,dp)/REAL(ostate%maxfun,dp)*100._dp),ostate%fopt
      END IF

    END DO

    ostate%state = 8
    CALL powell_optimize (ostate%nvar, x, ostate)
    CALL getvar(x,co,num_pol,num_gau)

    CALL release_opgrid(density,error)

    IF ( iunit > 0 ) THEN
      WRITE(iunit,'(" POWELL| Number of function evaluations",T71,I10)') ostate%nf
      WRITE(iunit,'(" POWELL| Final value of function",T61,G20.10)') ostate%fopt
      WRITE(iunit,'(" Optimized local potential of approximated nonadditive kinetic energy functional")')
      DO ig=1,num_gau
         WRITE(iunit,'(I2,T15,"Gaussian polynomial expansion",T66,"Rc=",F12.4)') ig,co(1,ig)
         WRITE(iunit,'(T15,"Coefficients",T33,4F12.4)') (co(1+ip,ig),ip=1,num_pol)
      END DO
    END IF

    CALL open_file(file_name="FIT_RESULT",file_status="UNKNOWN",file_action="WRITE",unit_number=iw)
    WRITE(iw,*) ptable(atom%z)%symbol
    WRITE(iw,*) num_gau, num_pol
    DO ig=1,num_gau
      WRITE(iw,'(T10,F12.4,6X,4F12.4)') (co(ip,ig),ip=1,num_pol+1)
    END DO
    CALL close_file(unit_number=iw)

    IF(PRESENT(results)) THEN
      CPPrecondition(SIZE(results)>=SIZE(x), cp_failure_level, routineP, error, failure)
      results = x
    END IF

    DEALLOCATE(co,x,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE atom_fit_kgpot

! *****************************************************************************
!> \brief ...
!> \param kgpot ...
!> \param ng ...
!> \param np ...
!> \param cval ...
!> \param aerr ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE kgpot_fit (kgpot,ng,np,cval,aerr,error)
    TYPE(opgrid_type), POINTER               :: kgpot
    INTEGER, INTENT(IN)                      :: ng, np
    REAL(dp), DIMENSION(:, :), INTENT(INOUT) :: cval
    REAL(dp), INTENT(OUT)                    :: aerr
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kgpot_fit', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ierr, ig, ip, n
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: pc, rc
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: pval

    n = kgpot%grid%nr
    ALLOCATE(pval(n),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    pval = 0.0_dp
    DO i=1,n
       DO ig=1,ng
          rc = kgpot%grid%rad(i)/cval(1,ig)
          pc = 0.0_dp
          DO ip=1,np
             pc = pc + cval(ip+1,ig)*rc**(2*ip-2)
          END DO
          pval(i) = pval(i) + pc*EXP(-0.5_dp*rc*rc)
       END DO
    END DO
    pval(1:n) = (pval(1:n) - kgpot%op(1:n))**2
    aerr = fourpi*SUM(pval(1:n)*kgpot%grid%wr(1:n))

    DEALLOCATE(pval,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE kgpot_fit 

! *****************************************************************************
!> \brief ...
!> \param xvar ...
!> \param cvar ...
!> \param np ...
!> \param ng ...
! *****************************************************************************
  SUBROUTINE getvar(xvar,cvar,np,ng)
    REAL(KIND=dp), DIMENSION(:), INTENT(in)  :: xvar
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(inout)                          :: cvar
    INTEGER                                  :: np, ng

    INTEGER                                  :: ig, ii, ip

    ii = 0
    DO ig=1,ng
       ii=ii+1
       cvar(1,ig) = xvar(ii)
       DO ip=1,np
          ii=ii+1
          cvar(ip+1,ig) = xvar(ii)**2
       END DO
    END DO

  END SUBROUTINE getvar

! *****************************************************************************
!> \brief ...
!> \param xvar ...
!> \param cvar ...
!> \param np ...
!> \param ng ...
! *****************************************************************************
  SUBROUTINE putvar(xvar,cvar,np,ng)
    REAL(KIND=dp), DIMENSION(:), &
      INTENT(inout)                          :: xvar
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(in)                             :: cvar
    INTEGER                                  :: np, ng

    INTEGER                                  :: ig, ii, ip

    ii = 0
    DO ig=1,ng
       ii=ii+1
       xvar(ii) = cvar(1,ig)
       DO ip=1,np
          ii=ii+1
          xvar(ii) = SQRT(ABS(cvar(ip+1,ig)))
       END DO
    END DO

  END SUBROUTINE putvar
! *****************************************************************************

END MODULE atom_fit
