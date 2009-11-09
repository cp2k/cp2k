!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE atom_energy
  USE atom_electronic_structure,       ONLY: calculate_atom
  USE atom_fit,                        ONLY: atom_fit_density
  USE atom_operators,                  ONLY: atom_int_release,&
                                             atom_int_setup,&
                                             atom_ppint_release,&
                                             atom_ppint_setup,&
                                             atom_relint_release,&
                                             atom_relint_setup
  USE atom_output,                     ONLY: atom_print_basis,&
                                             atom_print_info,&
                                             atom_print_method,&
                                             atom_print_potential
  USE atom_types,                      ONLY: &
       atom_basis_type, atom_integrals, atom_optimization_type, &
       atom_orbitals, atom_p_type, atom_potential_type, atom_state, &
       atom_type, create_atom_orbs, create_atom_type, init_atom_basis, &
       init_atom_potential, read_atom_opt_section, release_atom_basis, &
       release_atom_potential, release_atom_type, set_atom
  USE atom_utils,                      ONLY: atom_consistent_method,&
                                             atom_set_occupation,&
                                             get_maxl_occ,&
                                             get_maxn_occ
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE f77_blas
  USE input_constants,                 ONLY: do_analytic
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE lapack
  USE periodic_table,                  ONLY: nelem,&
                                             ptable
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: atom_energy_opt

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_energy'

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE atom_energy_opt(atom_section,error)
    TYPE(section_vals_type), POINTER         :: atom_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_energy_opt', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: elem
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: tmpstringlist
    INTEGER :: do_eric, do_erie, handle, i, ierr, im, in, iw, k, maxl, mb, &
      method, mo, multiplicity, n_meth, n_rep, nder, num_gto, reltyp, zcore, &
      zval, zz
    INTEGER, DIMENSION(0:3)                  :: maxn
    INTEGER, DIMENSION(:), POINTER           :: cn
    LOGICAL                                  :: eri_c, eri_e, failure, &
                                                had_ae, had_pp, pp_calc
    REAL(KIND=dp)                            :: delta
    REAL(KIND=dp), DIMENSION(0:3, 10)        :: pocc
    TYPE(atom_basis_type), POINTER           :: ae_basis, pp_basis
    TYPE(atom_integrals), POINTER            :: ae_int, pp_int
    TYPE(atom_optimization_type)             :: optimization
    TYPE(atom_orbitals), POINTER             :: orbitals
    TYPE(atom_p_type), DIMENSION(:, :), &
      POINTER                                :: atom_info
    TYPE(atom_potential_type), POINTER       :: ae_pot, p_pot
    TYPE(atom_state), POINTER                :: state
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER :: basis_section, method_section, &
      opt_section, potential_section, powell_section, xc_section

    failure = .FALSE.

    CALL timeset(routineN,handle)

    ! What atom do we calculate
    CALL section_vals_val_get(atom_section,"ATOMIC_NUMBER", i_val=zval, error=error)
    CALL section_vals_val_get(atom_section,"ELEMENT", c_val=elem, error=error)
    zz = 0
    DO i=1,nelem
      IF ( ptable(i)%symbol == elem ) THEN
        zz = i
        EXIT
      END IF
    END DO
    IF ( zz /= 1 ) zval = zz

    ! read and set up inofrmation on the basis sets
    ALLOCATE(ae_basis,pp_basis,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    basis_section => section_vals_get_subs_vals(atom_section,"AE_BASIS",error=error)
    NULLIFY(ae_basis%grid)
    CALL init_atom_basis(ae_basis,basis_section,zval,"AE",error)
    NULLIFY(pp_basis%grid)
    basis_section => section_vals_get_subs_vals(atom_section,"PP_BASIS",error=error)
    CALL init_atom_basis(pp_basis,basis_section,zval,"PP",error)

    ! print general and basis set information
    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%PROGRAM_BANNER",extension=".log",error=error)
    IF(iw > 0) CALL atom_print_info(zval,"Atomic Energy Calculation",iw,error)
    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%PROGRAM_BANNER",error=error)
    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%BASIS_SET",extension=".log",error=error)
    IF(iw > 0) THEN
      CALL atom_print_basis(ae_basis,iw," All Electron Basis",error)
      CALL atom_print_basis(pp_basis,iw," Pseudopotential Basis",error)
    END IF
    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%BASIS_SET",error=error)

    ! read and setup information on the pseudopotential
    NULLIFY(potential_section)
    potential_section => section_vals_get_subs_vals(atom_section,"POTENTIAL",error=error)
    ALLOCATE(ae_pot,p_pot,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    CALL init_atom_potential(p_pot,potential_section,zval,error)
    CALL init_atom_potential(ae_pot,potential_section,-1,error)

    ! if the ERI's are calculated analytically, we have to precalculate them
    eri_c = .FALSE.
    CALL section_vals_val_get(atom_section,"COULOMB_INTEGRALS", i_val=do_eric, error=error)
    IF(do_eric==do_analytic) eri_c = .TRUE.
    eri_e = .FALSE.
    CALL section_vals_val_get(atom_section,"EXCHANGE_INTEGRALS", i_val=do_erie, error=error)
    IF(do_erie==do_analytic) eri_e = .TRUE.

    ! information on the states to be calculated
    CALL section_vals_val_get(atom_section,"MAX_ANGULAR_MOMENTUM", i_val=maxl, error=error)
    maxn=0
    CALL section_vals_val_get(atom_section,"CALCULATE_STATES", i_vals=cn, error=error)
    DO in = 1, MIN(SIZE(cn),4)
      maxn(in-1) = cn(in)
    END DO
    DO in = 0, 3
      maxn(in) = MIN(maxn(in),ae_basis%nbas(in))
      maxn(in) = MIN(maxn(in),pp_basis%nbas(in))
    END DO

    ! read optimization section
    opt_section => section_vals_get_subs_vals(atom_section,"OPTIMIZATION",error=error)
    CALL read_atom_opt_section(optimization,opt_section,error)

    had_ae = .FALSE.
    had_pp = .FALSE.

    ! Check for the total number of electron configurations to be calculated
    CALL section_vals_val_get(atom_section,"ELECTRON_CONFIGURATION", n_rep_val=n_rep, error=error)
    ! Check for the total number of method types to be calculated
    method_section => section_vals_get_subs_vals(atom_section,"METHOD",error=error)
    CALL section_vals_get(method_section,n_repetition=n_meth,error=error)
 
    ! integrals
    ALLOCATE(ae_int, pp_int,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    ALLOCATE(atom_info(n_rep,n_meth),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    DO in = 1, n_rep
      DO im = 1, n_meth

        NULLIFY(atom_info(in,im)%atom)
        CALL create_atom_type(atom_info(in,im)%atom,error)

        atom_info(in,im)%atom%optimization = optimization 

        atom_info(in,im)%atom%z = zval
        xc_section => section_vals_get_subs_vals(method_section,"XC",i_rep_section=im,error=error)
        atom_info(in,im)%atom%xc_section => xc_section

        ALLOCATE(state,STAT=ierr)
        CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

        ! get the electronic configuration
        CALL section_vals_val_get(atom_section,"ELECTRON_CONFIGURATION", i_rep_val=in,&
                                  c_vals=tmpstringlist, error=error)

        ! set occupations 
        CALL atom_set_occupation(tmpstringlist,state%occ,state%occupation,state%multiplicity,error)
        state%maxl_occ = get_maxl_occ(state%occ)
        state%maxn_occ = get_maxn_occ(state%occ)

        ! set number of states to be calculated
        state%maxl_calc = MAX(maxl,state%maxl_occ)
        state%maxl_calc = MIN(3,state%maxl_calc)
        state%maxn_calc = 0
        DO k=0,state%maxl_calc
          state%maxn_calc(k) = MAX(maxn(k),state%maxn_occ(k))
        END DO

        ! is there a pseudo potential
        pp_calc = INDEX(tmpstringlist(1),"CORE") /= 0
        IF ( pp_calc ) THEN
          ! get and set the core occupations 
          CALL section_vals_val_get(atom_section,"CORE", c_vals=tmpstringlist, error=error)
          CALL atom_set_occupation(tmpstringlist,state%core,pocc,error=error)
          zcore = zval - SUM(state%core)
          CALL set_atom(atom_info(in,im)%atom,zcore=zcore,pp_calc=.TRUE.,error=error)
        ELSE
          state%core=0._dp
          CALL set_atom(atom_info(in,im)%atom,zcore=zval,pp_calc=.FALSE.,error=error)
        END IF

        CALL section_vals_val_get(method_section,"METHOD_TYPE",i_val=method,i_rep_section=im,error=error)
        CALL section_vals_val_get(method_section,"RELATIVISTIC",i_val=reltyp,i_rep_section=im,error=error)
        CALL set_atom(atom_info(in,im)%atom,method_type=method,relativistic=reltyp,error=error)

        IF(atom_consistent_method(method,state%multiplicity)) THEN
          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%METHOD_INFO",extension=".log",error=error)
          CALL atom_print_method(atom_info(in,im)%atom,iw,error)
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%METHOD_INFO",error=error)
  
          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%POTENTIAL",extension=".log",error=error)
          IF ( pp_calc ) THEN
            IF(iw > 0) CALL atom_print_potential(p_pot,iw,error)
          ELSE
            IF(iw > 0) CALL atom_print_potential(ae_pot,iw,error)
          END IF
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%POTENTIAL",error=error)
        END IF

        ! calculate integrals
        IF ( pp_calc ) THEN
          ! general integrals
          CALL atom_int_setup(pp_int,pp_basis,&
                              potential=p_pot,eri_coulomb=eri_c,eri_exchange=eri_e,error=error)
          ! potential
          CALL atom_ppint_setup(pp_int,pp_basis,potential=p_pot,error=error)
          !
          NULLIFY(pp_int%tzora,pp_int%hdkh)
          !
          CALL set_atom(atom_info(in,im)%atom,basis=pp_basis,integrals=pp_int,potential=p_pot,error=error)
          state%maxn_calc(:) = MIN( state%maxn_calc(:), pp_basis%nbas(:) )
          CPPostcondition(ALL(state%maxn_calc(:) >= state%maxn_occ), cp_failure_level, routineP, error, failure)
          had_pp = .TRUE.
        ELSE
          ! general integrals
          CALL atom_int_setup(ae_int,ae_basis,potential=ae_pot,&
                              eri_coulomb=eri_c,eri_exchange=eri_e,error=error)
          ! potential
          CALL atom_ppint_setup(ae_int,ae_basis,potential=ae_pot,error=error)
          ! relativistic correction terms
          CALL atom_relint_setup(ae_int,ae_basis,reltyp,zcore=REAL(zval,dp),error=error)
          !
          CALL set_atom(atom_info(in,im)%atom,basis=ae_basis,integrals=ae_int,potential=ae_pot,error=error)
          state%maxn_calc(:) = MIN( state%maxn_calc(:), ae_basis%nbas(:) )
          CPPostcondition(ALL(state%maxn_calc(:) >= state%maxn_occ), cp_failure_level, routineP, error, failure)
          had_ae = .TRUE.
        END IF

        CALL set_atom(atom_info(in,im)%atom,state=state,error=error)

        CALL set_atom(atom_info(in,im)%atom,coulomb_integral_type=do_eric,&
                      exchange_integral_type=do_erie,error=error)

        NULLIFY(orbitals)
        mo = MAXVAL(state%maxn_calc)
        mb = MAXVAL(atom_info(in,im)%atom%basis%nbas)
        CALL create_atom_orbs(orbitals,mb,mo,error)
        CALL set_atom(atom_info(in,im)%atom,orbitals=orbitals,error=error)

        IF(atom_consistent_method(method,state%multiplicity)) THEN
          !Calculate the electronic structure
          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%SCF_INFO",extension=".log",error=error)
          CALL calculate_atom(atom_info(in,im)%atom,iw,error=error)
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%SCF_INFO",error=error)

          ! perform a fit of the total electronic density
          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%FIT_DENSITY",extension=".log",error=error)
          IF (iw>0) THEN
            CALL section_vals_val_get(atom_section,"PRINT%FIT_DENSITY%NUM_GTO", i_val=num_gto, error=error)
            powell_section => section_vals_get_subs_vals(atom_section,"POWELL",error=error)
            CALL atom_fit_density(atom_info(in,im)%atom,num_gto,iw,powell_section,error)
          END IF
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%FIT_DENSITY",error=error)

          ! generate a response basis
          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%RESPONSE_BASIS",extension=".log",error=error)
          IF (iw>0) THEN
            CALL section_vals_val_get(atom_section,"PRINT%RESPONSE_BASIS%DELTA_CHARGE", r_val=delta, error=error)
            CALL section_vals_val_get(atom_section,"PRINT%RESPONSE_BASIS%DERIVATIVES", i_val=nder, error=error)
            CALL atom_response_basis(atom_info(in,im)%atom,delta,nder,iw,error)
          END IF
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%RESPONSE_BASIS",error=error)
        END IF

      END DO
    END DO

    ! clean up
    IF ( had_ae ) THEN
      CALL atom_int_release(ae_int,error)
      CALL atom_ppint_release(ae_int,error)
      CALL atom_relint_release(ae_int,error)
    END IF
    IF ( had_pp ) THEN
      CALL atom_int_release(pp_int,error)
      CALL atom_ppint_release(pp_int,error)
      CALL atom_relint_release(pp_int,error)
    END IF
    CALL release_atom_basis(ae_basis,error)
    CALL release_atom_basis(pp_basis,error)

    CALL release_atom_potential(p_pot,error)
    CALL release_atom_potential(ae_pot,error)

    DO in = 1, n_rep
      DO im = 1, n_meth
        CALL release_atom_type(atom_info(in,im)%atom,error)
      END DO
    END DO
    DEALLOCATE(atom_info,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    DEALLOCATE(ae_pot,p_pot,ae_basis,pp_basis,ae_int,pp_int,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    CALL timestop(handle)

  END SUBROUTINE atom_energy_opt

! *****************************************************************************
  SUBROUTINE atom_response_basis(atom,delta,nder,iw,error)
 
    TYPE(atom_type), POINTER                 :: atom
    REAL(KIND=dp), INTENT(IN)                :: delta
    INTEGER, INTENT(IN)                      :: nder, iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_response_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ider, ierr, k, l, lhomo, &
                                                lwork, m, n, nhomo, s1, s2
    LOGICAL                                  :: failure = .FALSE.
    REAL(KIND=dp)                            :: dene, emax, fhomo
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: w, work
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: amat
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: rbasis
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :, :)                  :: wfn
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: ovlp
    TYPE(atom_state), POINTER                :: state

    WRITE(iw,'(/," ",79("*"),/,T34,A,/," ",79("*"))') "RESPONSE BASIS"

    state => atom%state
    ovlp => atom%integrals%ovlp

    ! find HOMO
    lhomo = -1
    nhomo = -1
    emax=-HUGE(1._dp)
    DO l=0,state%maxl_occ
      DO i=1,state%maxn_occ(l)
        IF (atom%orbitals%ener(i,l) > emax) THEN
          lhomo=l
          nhomo=i
          emax=atom%orbitals%ener(i,l)
          fhomo=state%occupation(l,i)
        END IF
      END DO
    END DO

    s1=SIZE(atom%orbitals%wfn,1)
    s2=SIZE(atom%orbitals%wfn,2)
    ALLOCATE(wfn(s1,s2,0:3,-nder:nder),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)
    s2=MAXVAL(state%maxn_occ)+nder
    ALLOCATE(rbasis(s1,s2,0:3),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    DO ider=-nder,nder
      dene = REAL(ider,KIND=dp)*delta
      CPPostcondition(fhomo>ABS(dene), cp_failure_level, routineP, error, failure)
      state%occupation(lhomo,nhomo)=fhomo+dene
      CALL calculate_atom(atom,iw=0,noguess=.TRUE.,error=error)
      wfn(:,:,:,ider) = atom%orbitals%wfn
      state%occupation(lhomo,nhomo)=fhomo
    END DO

    DO l=0,state%maxl_occ
      ! occupied states
      DO i=1,state%maxn_occ(l)
        rbasis(:,i,l) = wfn(:,i,l,0)
      END DO
      ! differentiation
      DO ider=1,nder
        i=state%maxn_occ(l)
        SELECT CASE (ider)
          CASE (1)
            rbasis(:,i+1,l) = 0.5_dp*(wfn(:,i,l,1) - wfn(:,i,l,-1))/delta
          CASE (2)
            rbasis(:,i+2,l) = 0.25_dp*(wfn(:,i,l,2) - 2._dp*wfn(:,i,l,0) + wfn(:,i,l,-2))/delta**2
          CASE (3)
            rbasis(:,i+3,l) = 0.125_dp*(wfn(:,i,l,3) - 3._dp*wfn(:,i,l,1) &
                                      + 3._dp*wfn(:,i,l,-1) - wfn(:,i,l,-3))/delta**3
          CASE DEFAULT
            CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
        END SELECT
      END DO
      ! orthogonalization
      n = state%maxn_occ(l)+nder
      m = atom%basis%nbas(l)
      lwork = n*n
      ALLOCATE(amat(n,n),w(n),work(lwork),STAT=ierr)
      CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
      amat(1:n,1:n) = MATMUL(TRANSPOSE(rbasis(1:m,1:n,l)),MATMUL(ovlp(1:m,1:m,l),rbasis(1:m,1:n,l)))
      CALL lapack_ssyev ( "V", "U", n, amat(1:n,1:n), n, w(1:n), work, lwork, ierr )
      CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
      DO i=1,n
        amat(1:n,i)=amat(1:n,i)/SQRT(w(i))
      END DO
      rbasis(1:m,1:n,l) = MATMUL(rbasis(1:m,1:n,l),amat(1:n,1:n))
      DEALLOCATE(amat,w,work,STAT=ierr)
      CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)

      WRITE(iw,'(/,A,T30,I3)') " Angular momentum :",l
      WRITE(iw,'(A)') " Basis             Coefficients"
      DO i=1,m
        WRITE(iw,'(I6,4X,15ES20.6)') i,(rbasis(i,k,l),k=1,n) 
      END DO

    END DO

    DEALLOCATE(wfn,rbasis,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    WRITE(iw,'(" ",79("*"))')

  END SUBROUTINE atom_response_basis
! *****************************************************************************

END MODULE atom_energy
