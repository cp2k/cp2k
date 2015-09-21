!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE atom_basis
  USE atom_fit,                        ONLY: atom_fit_basis
  USE atom_output,                     ONLY: atom_print_basis,&
                                             atom_print_info,&
                                             atom_print_method,&
                                             atom_print_potential
  USE atom_types,                      ONLY: &
       atom_basis_type, atom_integrals, atom_optimization_type, &
       atom_orbitals, atom_p_type, atom_potential_type, atom_state, &
       create_atom_orbs, create_atom_type, init_atom_basis, &
       init_atom_potential, read_atom_opt_section, release_atom_basis, &
       release_atom_potential, release_atom_type, set_atom
  USE atom_utils,                      ONLY: atom_consistent_method,&
                                             atom_set_occupation,&
                                             get_maxl_occ,&
                                             get_maxn_occ
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE input_constants,                 ONLY: do_analytic
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE periodic_table,                  ONLY: nelem,&
                                             ptable
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: atom_basis_opt

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_basis'

! *****************************************************************************

CONTAINS

! *****************************************************************************

! *****************************************************************************
!> \brief ...
!> \param atom_section ...
! *****************************************************************************
  SUBROUTINE atom_basis_opt(atom_section)
    TYPE(section_vals_type), POINTER         :: atom_section

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_basis_opt', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: elem
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: tmpstringlist
    INTEGER                                  :: do_eric, do_erie, handle, i, &
                                                im, in, iunit, iw, k, maxl, &
                                                mb, method, mo, n_meth, &
                                                n_rep, reltyp, zcore, zval, zz
    INTEGER, DIMENSION(0:3)                  :: maxn
    INTEGER, DIMENSION(:), POINTER           :: cn
    LOGICAL                                  :: eri_c, eri_e, failure, &
                                                had_ae, had_pp, pp_calc
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
    CALL section_vals_val_get(atom_section,"ATOMIC_NUMBER", i_val=zval)
    CALL section_vals_val_get(atom_section,"ELEMENT", c_val=elem)
    zz = 0
    DO i=1,nelem
      IF ( ptable(i)%symbol == elem ) THEN
        zz = i
        EXIT
      END IF
    END DO
    IF ( zz /= 1 ) zval = zz

    ! read and set up inofrmation on the basis sets
    ALLOCATE(ae_basis,pp_basis)
    basis_section => section_vals_get_subs_vals(atom_section,"AE_BASIS")
    NULLIFY(ae_basis%grid)
    CALL init_atom_basis(ae_basis,basis_section,zval,"AE")
    NULLIFY(pp_basis%grid)
    basis_section => section_vals_get_subs_vals(atom_section,"PP_BASIS")
    CALL init_atom_basis(pp_basis,basis_section,zval,"PP")

    ! print general and basis set information
    logger => cp_get_default_logger()
    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%PROGRAM_BANNER",extension=".log")
    IF(iw > 0) CALL atom_print_info(zval,"Atomic Basis Optimization",iw)
    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%PROGRAM_BANNER")

    ! read and setup information on the pseudopotential
    NULLIFY(potential_section)
    potential_section => section_vals_get_subs_vals(atom_section,"POTENTIAL")
    ALLOCATE(ae_pot,p_pot)
    CALL init_atom_potential(p_pot,potential_section,zval)
    CALL init_atom_potential(ae_pot,potential_section,-1)

    ! if the ERI's are calculated analytically, we have to precalculate them
    eri_c = .FALSE.
    CALL section_vals_val_get(atom_section,"COULOMB_INTEGRALS", i_val=do_eric)
    IF(do_eric==do_analytic) eri_c = .TRUE.
    eri_e = .FALSE.
    CALL section_vals_val_get(atom_section,"EXCHANGE_INTEGRALS", i_val=do_erie)
    IF(do_erie==do_analytic) eri_e = .TRUE.

    ! information on the states to be calculated
    CALL section_vals_val_get(atom_section,"MAX_ANGULAR_MOMENTUM", i_val=maxl)
    maxn=0
    CALL section_vals_val_get(atom_section,"CALCULATE_STATES", i_vals=cn)
    DO in = 1, MIN(SIZE(cn),4)
      maxn(in-1) = cn(in)
    END DO
    DO in = 0, 3
      maxn(in) = MIN(maxn(in),ae_basis%nbas(in))
      maxn(in) = MIN(maxn(in),pp_basis%nbas(in))
    END DO

    ! read optimization section
    opt_section => section_vals_get_subs_vals(atom_section,"OPTIMIZATION")
    CALL read_atom_opt_section(optimization,opt_section)

    had_ae = .FALSE.
    had_pp = .FALSE.

    ! Check for the total number of electron configurations to be calculated
    CALL section_vals_val_get(atom_section,"ELECTRON_CONFIGURATION", n_rep_val=n_rep)
    ! Check for the total number of method types to be calculated
    method_section => section_vals_get_subs_vals(atom_section,"METHOD")
    CALL section_vals_get(method_section,n_repetition=n_meth)

    ! integrals
    ALLOCATE(ae_int, pp_int)

    ALLOCATE(atom_info(n_rep,n_meth))

    DO in = 1, n_rep
      DO im = 1, n_meth

        NULLIFY(atom_info(in,im)%atom)
        CALL create_atom_type(atom_info(in,im)%atom)

        atom_info(in,im)%atom%optimization = optimization

        atom_info(in,im)%atom%z = zval
        xc_section => section_vals_get_subs_vals(method_section,"XC",i_rep_section=im)
        atom_info(in,im)%atom%xc_section => xc_section

        ALLOCATE(state)

        ! get the electronic configuration
        CALL section_vals_val_get(atom_section,"ELECTRON_CONFIGURATION", i_rep_val=in,&
                                  c_vals=tmpstringlist)

        ! set occupations
        CALL atom_set_occupation(tmpstringlist,state%occ,state%occupation,state%multiplicity)
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
        pp_calc = ANY(INDEX(tmpstringlist(1:),"CORE") /= 0)
        IF ( pp_calc ) THEN
          ! get and set the core occupations
          CALL section_vals_val_get(atom_section,"CORE", c_vals=tmpstringlist)
          CALL atom_set_occupation(tmpstringlist,state%core,pocc)
          zcore = zval - NINT(SUM(state%core))
          CALL set_atom(atom_info(in,im)%atom,zcore=zcore,pp_calc=.TRUE.)
          had_pp = .TRUE.
          CALL set_atom(atom_info(in,im)%atom,basis=pp_basis,potential=p_pot)
          state%maxn_calc(:) = MIN( state%maxn_calc(:), pp_basis%nbas(:) )
          CPASSERT(ALL(state%maxn_calc(:) >= state%maxn_occ))
        ELSE
          state%core=0._dp
          CALL set_atom(atom_info(in,im)%atom,zcore=zval,pp_calc=.FALSE.)
          had_ae = .TRUE.
          CALL set_atom(atom_info(in,im)%atom,basis=ae_basis,potential=ae_pot)
          state%maxn_calc(:) = MIN( state%maxn_calc(:), ae_basis%nbas(:) )
          CPASSERT(ALL(state%maxn_calc(:) >= state%maxn_occ))
        END IF

        CALL section_vals_val_get(method_section,"METHOD_TYPE",i_val=method,i_rep_val=im)
        CALL section_vals_val_get(method_section,"RELATIVISTIC",i_val=reltyp,i_rep_section=im)
        CALL set_atom(atom_info(in,im)%atom,method_type=method,relativistic=reltyp)
        CALL set_atom(atom_info(in,im)%atom,state=state)
        CALL set_atom(atom_info(in,im)%atom,coulomb_integral_type=do_eric,&
                      exchange_integral_type=do_erie)

        IF (atom_consistent_method(method,state%multiplicity)) THEN
          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%METHOD_INFO",extension=".log")
          CALL atom_print_method(atom_info(in,im)%atom,iw)
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%METHOD_INFO")
          iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%POTENTIAL",extension=".log")
          IF (pp_calc) THEN
            IF (iw > 0) CALL atom_print_potential(p_pot,iw)
          ELSE
            IF (iw > 0) CALL atom_print_potential(ae_pot,iw)
          END IF
          CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%POTENTIAL")
        ELSE
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,routineP,&
                         "METHOD_TYPE and MULTIPLICITY are incompatible",&
                         only_ionode=.TRUE.)
        END IF

        NULLIFY(orbitals)
        mo = MAXVAL(state%maxn_calc)
        mb = MAXVAL(atom_info(in,im)%atom%basis%nbas)
        CALL create_atom_orbs(orbitals,mb,mo)
        CALL set_atom(atom_info(in,im)%atom,orbitals=orbitals)

      END DO
    END DO

    ! Start the Optimization
    powell_section => section_vals_get_subs_vals(atom_section,"POWELL")
    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%SCF_INFO",extension=".log")
    iunit = cp_print_key_unit_nr(logger,atom_section,"PRINT%FIT_BASIS",extension=".log")
    IF ( had_ae ) THEN
      pp_calc = .FALSE.
      CALL atom_fit_basis(atom_info,ae_basis,pp_calc,iunit,powell_section)
    END IF
    IF ( had_pp ) THEN
      pp_calc = .TRUE.
      CALL atom_fit_basis(atom_info,pp_basis,pp_calc,iunit,powell_section)
    END IF
    CALL cp_print_key_finished_output(iunit,logger,atom_section,"PRINT%FIT_BASIS")
    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%SCF_INFO")
    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%BASIS_SET",extension=".log")
    IF(iw > 0) THEN
      CALL atom_print_basis(ae_basis,iw," All Electron Basis")
      CALL atom_print_basis(pp_basis,iw," Pseudopotential Basis")
    END IF
    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%BASIS_SET")

    CALL release_atom_basis(ae_basis)
    CALL release_atom_basis(pp_basis)

    CALL release_atom_potential(p_pot)
    CALL release_atom_potential(ae_pot)

    DO in = 1, n_rep
      DO im = 1, n_meth
        CALL release_atom_type(atom_info(in,im)%atom)
      END DO
    END DO
    DEALLOCATE(atom_info)

    DEALLOCATE(ae_pot,p_pot,ae_basis,pp_basis,ae_int,pp_int)

    CALL timestop(handle)

  END SUBROUTINE atom_basis_opt

! *****************************************************************************

END MODULE atom_basis
