!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE atom
  USE atom_fit,                        ONLY: atom_fit_density
  USE atom_operators,                  ONLY: atom_int_release,&
                                             atom_int_setup,&
                                             atom_ppint_release,&
                                             atom_ppint_setup,&
                                             atom_relint_release,&
                                             atom_relint_setup
  USE atom_optimization,               ONLY: atom_history_init,&
                                             atom_history_release,&
                                             atom_history_type,&
                                             atom_history_update,&
                                             atom_opt
  USE atom_output,                     ONLY: atom_print_basis,&
                                             atom_print_energies,&
                                             atom_print_info,&
                                             atom_print_iteration,&
                                             atom_print_method,&
                                             atom_print_potential,&
                                             atom_print_state
  USE atom_types,                      ONLY: &
       GTH_PSEUDO, NO_PSEUDO, atom_basis_type, atom_integrals, atom_orbitals, &
       atom_p_type, atom_potential_type, atom_state, atom_type, &
       create_atom_orbs, create_atom_type, create_opgrid, create_opmat, &
       init_atom_basis, init_atom_potential, opgrid_type, opmat_type, &
       release_atom_basis, release_atom_potential, release_atom_type, &
       release_opgrid, release_opmat, set_atom
  USE atom_utils,                      ONLY: &
       atom_denmat, atom_density, atom_set_occupation, atom_solve, &
       atom_trace, ceri_contract, coulomb_potential_analytic, &
       coulomb_potential_numeric, eeri_contract, err_matrix, get_maxl_occ, &
       get_maxn_occ, numpot_matrix
  USE atom_xc,                         ONLY: calculate_atom_vxc
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE header,                          ONLY: atom_footer,&
                                             atom_header
  USE input_constants,                 ONLY: &
       atom_energy_run, atom_no_run, do_analytic, do_dkh0_atom, do_dkh1_atom, &
       do_dkh2_atom, do_dkh3_atom, do_dkh4_atom, do_dkh5_atom, &
       do_nonrel_atom, do_numeric, do_rhf_atom, do_rks_atom, do_rohf_atom, &
       do_semi_analytic, do_uhf_atom, do_uks_atom, do_zoramp_atom
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
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: atom_code

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom'

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE atom_code(root_section,para_env,error)
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_code', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, iw, run_type_id
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: atom_section

    failure=.FALSE.

    CALL timeset(routineN,handle)

    logger => cp_error_get_logger(error)
    NULLIFY(atom_section)
    atom_section => section_vals_get_subs_vals(root_section,"ATOM",error=error)

    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%PROGRAM_BANNER",extension=".log",error=error)
    IF (iw>0) CALL atom_header(iw)
    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%PROGRAM_BANNER",error=error)

    CALL atom_test(atom_section,iw,error)

    CALL section_vals_val_get(atom_section,"RUN_TYPE",i_val=run_type_id,error=error)
    SELECT CASE (run_type_id)
      CASE (atom_no_run)
        ! do (almost) nothing
      CASE (atom_energy_run)
        CALL atom_energy(atom_section,para_env,error)
      CASE default
        CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%PROGRAM_BANNER",extension=".log",error=error)
    IF (iw>0) CALL atom_footer(iw)
    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%PROGRAM_BANNER",error=error)

    CALL timestop(handle)

  END SUBROUTINE atom_code

! *****************************************************************************
  SUBROUTINE atom_test(atom_section,iw,error)
    TYPE(section_vals_type), POINTER         :: atom_section
    INTEGER, INTENT(IN)                      :: iw
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_test', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: elem
    INTEGER                                  :: i, ie, z
    LOGICAL                                  :: failure = .FALSE.

    CALL section_vals_val_get(atom_section,"ATOMIC_NUMBER", i_val=z, error=error)
    CALL section_vals_val_get(atom_section,"ELEMENT", c_val=elem, error=error)

    ie = 0
    DO i=1,nelem
      IF ( ptable(i)%symbol == elem ) THEN
        ie = i
        EXIT
      END IF
    END DO
    IF (ie /= z) THEN
      IF ( ie /= 1 .AND. z /= 1 ) THEN
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      END IF
    END IF

  END SUBROUTINE atom_test
! *****************************************************************************

  SUBROUTINE atom_energy(atom_section,para_env,error)
    TYPE(section_vals_type), POINTER         :: atom_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_energy', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=2)                         :: elem
    CHARACTER(LEN=default_string_length), &
      DIMENSION(:), POINTER                  :: tmpstringlist
    INTEGER                                  :: do_eri, handle, i, ierr, im, &
                                                in, iw, k, maxl, mb, mo, &
                                                n_meth, n_rep, num_gto, &
                                                zcore, zval, zz
    INTEGER, DIMENSION(0:3)                  :: maxn
    INTEGER, DIMENSION(:), POINTER           :: cn
    LOGICAL                                  :: eri_c, eri_e, failure, &
                                                had_ae, had_pp, pp_calc
    TYPE(atom_basis_type), TARGET            :: ae_basis, pp_basis
    TYPE(atom_integrals), TARGET             :: ae_int, pp_int
    TYPE(atom_orbitals), POINTER             :: orbitals
    TYPE(atom_p_type), DIMENSION(:, :), &
      POINTER                                :: atom_info
    TYPE(atom_potential_type), TARGET        :: ae_pot, p_pot
    TYPE(atom_state), POINTER                :: state
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: basis_section, &
                                                method_section, &
                                                potential_section

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
    CALL init_atom_potential(p_pot,potential_section,zval,error)
    CALL init_atom_potential(ae_pot,potential_section,-1,error)

    ! if the ERI's are calculated analytically, we have to precalculate them
    eri_c = .FALSE.
    CALL section_vals_val_get(atom_section,"COULOMB_INTEGRALS", i_val=do_eri, error=error)
    IF(do_eri==do_analytic) eri_c = .TRUE.
    eri_e = .FALSE.
    CALL section_vals_val_get(atom_section,"EXCHANGE_INTEGRALS", i_val=do_eri, error=error)
    IF(do_eri==do_analytic) eri_e = .TRUE.

    ! information on the states to be calculated
    CALL section_vals_val_get(atom_section,"MAX_ANGULAR_MOMENTUM", i_val=maxl, error=error)
    maxn=0
    CALL section_vals_val_get(atom_section,"CALCULATE_STATES", i_vals=cn, error=error)
    DO in = 1, MIN(SIZE(cn),4)
      maxn(in-1) = cn(in)
    END DO

    had_ae = .FALSE.
    had_pp = .FALSE.

    ! Check for the total number of electron configurations to be calculated
    CALL section_vals_val_get(atom_section,"ELECTRON_CONFIGURATION", n_rep_val=n_rep, error=error)
    ! Check for the total number of method types to be calculated
    method_section => section_vals_get_subs_vals(atom_section,"METHOD",error=error)
    CALL section_vals_get(method_section,n_repetition=n_meth,error=error)

    ALLOCATE(atom_info(n_rep,n_meth),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    DO in = 1, n_rep
      DO im = 1, n_meth

        NULLIFY(atom_info(in,im)%atom)
        CALL create_atom_type(atom_info(in,im)%atom,error)

        atom_info(in,im)%atom%method_section => method_section
        atom_info(in,im)%atom%method_section_number = im
        atom_info(in,im)%atom%z = zval

        ALLOCATE(state,STAT=ierr)
        CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

        ! get the electronic configuration
        CALL section_vals_val_get(atom_section,"ELECTRON_CONFIGURATION", i_rep_val=in,&
                                  c_vals=tmpstringlist, error=error)

        ! set occupations 
        CALL atom_set_occupation(tmpstringlist,state%occ,state%occupation,error)
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
          CALL atom_set_occupation(tmpstringlist,state%core,state%occupation,error)
          zcore = zval - SUM(state%core)
          CALL set_atom(atom_info(in,im)%atom,zcore=zcore,pp_calc=.TRUE.,error=error)
        ELSE
          state%core=0._dp
          CALL set_atom(atom_info(in,im)%atom,zcore=zval,pp_calc=.FALSE.,error=error)
        END IF

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
          CALL atom_relint_setup(ae_int,ae_basis,method_section,im,zcore=REAL(zval,dp),error=error)
          !
          CALL set_atom(atom_info(in,im)%atom,basis=ae_basis,integrals=ae_int,potential=ae_pot,error=error)
          state%maxn_calc(:) = MIN( state%maxn_calc(:), ae_basis%nbas(:) )
          CPPostcondition(ALL(state%maxn_calc(:) >= state%maxn_occ), cp_failure_level, routineP, error, failure)
          had_ae = .TRUE.
        END IF

        CALL set_atom(atom_info(in,im)%atom,state=state,error=error)

        NULLIFY(orbitals)
        mo = MAXVAL(state%maxn_calc)
        mb = MAXVAL(atom_info(in,im)%atom%basis%nbas)
        CALL create_atom_orbs(orbitals,mb,mo,error)
        CALL set_atom(atom_info(in,im)%atom,orbitals=orbitals,error=error)

        !Calculate the electronic structure
        CALL atom_electronic_structure(atom_info(in,im)%atom,atom_section,error)

        ! perform a fit of the total electronic density
        iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%FIT_DENSITY",extension=".log",error=error)
        IF (iw>0) THEN
          CALL section_vals_val_get(atom_section,"PRINT%FIT_DENSITY%NUM_GTO", i_val=num_gto, error=error)
          CALL atom_fit_density(atom_info(in,im)%atom,atom_section,num_gto,iw,error)
        END IF
        CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%FIT_DENSITY",error=error)

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

    CALL timestop(handle)

  END SUBROUTINE atom_energy

! *****************************************************************************
  SUBROUTINE atom_electronic_structure(atom,atom_section,error)
    TYPE(atom_type), POINTER                 :: atom
    TYPE(section_vals_type), POINTER         :: atom_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_electronic_structure', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: do_eri, handle, i, im, iter, &
                                                iw, l, method, reltyp
    LOGICAL                                  :: do_hfx, failure, need_x, &
                                                need_xc
    REAL(KIND=dp)                            :: deps, eps_scf, hf_frac
    TYPE(atom_history_type)                  :: history
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(opgrid_type), POINTER               :: cpot, density
    TYPE(opmat_type), POINTER                :: fmat, hcore, jmat, kmat, xcmat
    TYPE(section_vals_type), POINTER         :: hfx_sections, method_section, &
                                                opt_section, xc_section

    failure = .FALSE.

    CALL timeset(routineN,handle)

    method_section => atom%method_section
    im = atom%method_section_number
    CALL section_vals_val_get(method_section,"METHOD_TYPE", i_val=method, i_rep_val=im, error=error)
    xc_section => section_vals_get_subs_vals(method_section,"XC",i_rep_section=im,error=error)
    hfx_sections => section_vals_get_subs_vals(xc_section,"HF",error=error)
    CALL section_vals_get(hfx_sections,explicit=do_hfx,error=error)
    IF ( do_hfx ) THEN
      CALL section_vals_val_get(hfx_sections,"FRACTION", r_val=hf_frac, error=error)
    ELSE
      hf_frac = 0._dp
    END IF
!not yet installed
!   adiabatic_rescaling_section => section_vals_get_subs_vals(input,"DFT%XC%ADIABATIC_RESCALING",error=error)
!   CALL section_vals_get(adiabatic_rescaling_section,explicit=do_adiabatic_rescaling,error=error)

    opt_section => section_vals_get_subs_vals(atom_section,"OPTIMIZATION",error=error)

    SELECT CASE (method)
      CASE DEFAULT
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_rks_atom)
        need_x = do_hfx
        need_xc = .TRUE.
      CASE (do_uks_atom)
        need_x = do_hfx
        need_xc = .TRUE.
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_rhf_atom)
        need_x = .TRUE.
        need_xc = .FALSE.
        hf_frac = 1._dp
      CASE (do_uhf_atom)
        need_x = .TRUE.
        need_xc = .FALSE.
        hf_frac = 1._dp
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_rohf_atom)
        need_x = .TRUE.
        need_xc = .FALSE.
        hf_frac = 1._dp
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
    END SELECT

    ! check for relativistic method
    CALL section_vals_val_get(method_section,"RELATIVISTIC",i_val=reltyp,i_rep_section=im,error=error)

    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%SCF_INFO",extension=".log",error=error)
    IF (iw>0) THEN
      CALL atom_print_state(atom%state,iw,error)
    END IF

    eps_scf = 1.e-8_dp

    NULLIFY(hcore)
    CALL create_opmat(hcore,atom%basis%nbas,error)
    ! Pseudopotentials
    SELECT CASE (atom%potential%ppot_type)
      CASE DEFAULT
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (NO_PSEUDO)
        SELECT CASE (reltyp)
          CASE DEFAULT
            CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
          CASE (do_nonrel_atom)
            hcore%op = atom%integrals%kin - atom%zcore*atom%integrals%core
          CASE (do_zoramp_atom)
            hcore%op = atom%integrals%kin + atom%integrals%tzora - atom%zcore*atom%integrals%core
          CASE (do_dkh0_atom)
            hcore%op = atom%integrals%hdkh - atom%zcore*atom%integrals%core
          CASE (do_dkh1_atom,do_dkh2_atom,do_dkh3_atom,do_dkh4_atom,do_dkh5_atom)
            hcore%op = atom%integrals%hdkh
        END SELECT
      CASE (GTH_PSEUDO)
        hcore%op = atom%integrals%kin + atom%integrals%core + atom%integrals%hnl
    END SELECT
    ! add confinement potential (not included in relativistic transformations)
    IF ( atom%potential%confinement ) THEN
      hcore%op = hcore%op + atom%potential%acon * atom%integrals%conf
    END IF

    ! initial guess
    CALL atom_solve(hcore%op,atom%integrals%utrans,atom%orbitals%wfn,atom%orbitals%ener,&
                    atom%basis%nbas,atom%integrals%nne,atom%state%maxl_calc,error)
    CALL atom_denmat(atom%orbitals%pmat,atom%orbitals%wfn,atom%basis%nbas,atom%state%occupation,&
                     atom%state%maxl_occ,atom%state%maxn_occ,error)

    NULLIFY(fmat,jmat,kmat,xcmat)
    CALL create_opmat(fmat,atom%basis%nbas,error)
    CALL create_opmat(jmat,atom%basis%nbas,error)
    CALL create_opmat(kmat,atom%basis%nbas,error)
    CALL create_opmat(xcmat,atom%basis%nbas,error)

    NULLIFY(density,cpot)
    CALL create_opgrid(density,atom%basis%grid,error)
    CALL create_opgrid(cpot,atom%basis%grid,error)

    ! wavefunction history
    NULLIFY(history%dmat,history%hmat)
    CALL atom_history_init (history,opt_section,fmat%op,error)

    iter = 0
    DO            !SCF Loop

      ! Kinetic energy
      atom%energy%ekin = atom_trace(atom%integrals%kin,atom%orbitals%pmat,error)

      ! Band energy
      atom%energy%eband = 0._dp
      DO l=0,3
        DO i=1,MIN(SIZE(atom%state%occupation,2),SIZE(atom%orbitals%ener,1))
          atom%energy%eband = atom%energy%eband + atom%orbitals%ener(i,l)*atom%state%occupation(l,i)
        END DO
      END DO

      ! Pseudopotential energy
      SELECT CASE (atom%potential%ppot_type)
        CASE DEFAULT
          CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
        CASE (NO_PSEUDO)
          atom%energy%eploc = 0._dp
          atom%energy%epnl = 0._dp
        CASE (GTH_PSEUDO)
          atom%energy%eploc = atom_trace(atom%integrals%core,atom%orbitals%pmat,error)
          atom%energy%epnl = atom_trace(atom%integrals%hnl,atom%orbitals%pmat,error)
      END SELECT
      atom%energy%epseudo = atom%energy%eploc + atom%energy%epnl
  
      ! Core energy
      atom%energy%ecore = atom_trace(hcore%op,atom%orbitals%pmat,error)

      ! Confinement energy
      IF ( atom%potential%confinement ) THEN
        atom%energy%econfinement = atom_trace(atom%integrals%conf,atom%orbitals%pmat,error)
      ELSE
        atom%energy%econfinement = 0._dp
      END IF

      ! Hartree Term
      jmat%op = 0._dp
      CALL section_vals_val_get(atom_section,"COULOMB_INTEGRALS", i_val=do_eri, error=error)
      SELECT CASE (do_eri)
        CASE DEFAULT
          CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
        CASE (do_analytic)
          CALL ceri_contract(jmat%op,atom%integrals%ceri,atom%orbitals%pmat,atom%integrals%n,error=error)
        CASE (do_semi_analytic)
          CALL coulomb_potential_analytic(cpot%op,atom%orbitals%pmat,atom%basis,atom%basis%grid,&
                 atom%state%maxl_occ,error)
          CALL numpot_matrix(jmat%op,cpot%op,atom%basis,0,error)
        CASE (do_numeric)
          CALL atom_density(density%op,atom%orbitals%pmat,atom%basis,atom%state%maxl_occ,typ="RHO",error=error)
          CALL coulomb_potential_numeric(cpot%op,density%op,density%grid,error)
          CALL numpot_matrix(jmat%op,cpot%op,atom%basis,0,error)
      END SELECT
      atom%energy%ecoulomb = 0.5_dp * atom_trace(jmat%op,atom%orbitals%pmat,error)

      ! Exchange Term
      IF (need_x) THEN
        kmat%op = 0._dp
        CALL section_vals_val_get(atom_section,"EXCHANGE_INTEGRALS", i_val=do_eri, error=error)
        SELECT CASE (do_eri)
          CASE DEFAULT
            CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
          CASE (do_analytic)
            CALL eeri_contract(kmat%op,atom%integrals%eeri,atom%orbitals%pmat,atom%integrals%n,error=error)
          CASE (do_semi_analytic)
            WRITE(6,*) "Exchange semi-analytic"
          CASE (do_numeric)
            WRITE(6,*) "Exchange numeric"
        END SELECT
        atom%energy%eexchange = hf_frac * 0.5_dp * atom_trace(kmat%op,atom%orbitals%pmat,error)
        kmat%op = hf_frac*kmat%op
      ELSE
        kmat%op = 0._dp
        atom%energy%eexchange = 0._dp
      END IF

      ! XC
      IF (need_xc) THEN
        xcmat%op = 0._dp
        CALL calculate_atom_vxc(xcmat,atom,method_section,error)
      ELSE
        xcmat%op = 0._dp
        atom%energy%exc = 0._dp
      END IF
     
      ! Total energy
      atom%energy%etot = atom%energy%ecore + atom%energy%ecoulomb + atom%energy%eexchange + atom%energy%exc

      ! Total HF/KS matrix
      fmat%op = hcore%op + jmat%op + kmat%op + xcmat%op

      ! calculate error matrix
      CALL err_matrix(jmat%op,deps,fmat%op,atom%orbitals%pmat,atom%integrals%utrans,&
                      atom%integrals%uptrans,atom%basis%nbas,atom%integrals%nne,error)

      iter = iter + 1

      CALL atom_print_iteration(iter,deps,atom%energy%etot,iw,error)

      IF ( deps < eps_scf ) EXIT

      ! update history container and extrapolate KS matrix
      CALL atom_history_update (history,fmat%op,jmat%op,error)
      CALL atom_opt (fmat%op,history,deps,error)

      ! Solve HF/KS equations
      CALL atom_solve(fmat%op,atom%integrals%utrans,atom%orbitals%wfn,atom%orbitals%ener,&
                      atom%basis%nbas,atom%integrals%nne,atom%state%maxl_calc,error)
      CALL atom_denmat(atom%orbitals%pmat,atom%orbitals%wfn,atom%basis%nbas,atom%state%occupation,&
                       atom%state%maxl_occ,atom%state%maxn_occ,error)

    END DO        !SCF Loop

    CALL atom_print_energies(atom,iw,error)

    CALL atom_history_release(history,error)

    CALL release_opmat(fmat,error)
    CALL release_opmat(jmat,error)
    CALL release_opmat(kmat,error)
    CALL release_opmat(xcmat,error)
    CALL release_opmat(hcore,error)

    CALL release_opgrid(density,error)
    CALL release_opgrid(cpot,error)

    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%SCF_INFO",error=error)

    CALL timestop(handle)

  END SUBROUTINE atom_electronic_structure

! *****************************************************************************

END MODULE atom
