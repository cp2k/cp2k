!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE atom
  USE atom_operators,                  ONLY: atom_int_release,&
                                             atom_int_setup
  USE atom_optimization,               ONLY: atom_history_init,&
                                             atom_history_release,&
                                             atom_history_type,&
                                             atom_history_update,&
                                             atom_opt
  USE atom_output,                     ONLY: atom_print_energies,&
                                             atom_print_iteration,&
                                             atom_print_state
  USE atom_types,                      ONLY: &
       atom_basis_type, atom_integrals, atom_orbitals, atom_p_type, &
       atom_state, atom_type, create_atom_orbs, create_atom_type, &
       create_opgrid, create_opmat, init_atom_basis, opgrid_type, opmat_type, &
       release_atom_basis, release_atom_type, release_opgrid, release_opmat, &
       set_atom
  USE atom_utils,                      ONLY: &
       atom_denmat, atom_density, atom_set_occupation, atom_solve, &
       atom_trace, ceri_contract, coulomb_potential_analytic, &
       coulomb_potential_numeric, eeri_contract, err_matrix, get_maxl_occ, &
       get_maxn_occ, integrate_grid, numpot_matrix
  USE atom_xc,                         ONLY: calculate_atom_vxc
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE header,                          ONLY: atom_footer,&
                                             atom_header
  USE input_constants,                 ONLY: &
       atom_energy_run, atom_no_run, do_analytic, do_numeric, do_rhf_atom, &
       do_rks_atom, do_rohf_atom, do_semi_analytic, do_uhf_atom, do_uks_atom
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE mathconstants,                   ONLY: fourpi
  USE periodic_table,                  ONLY: nelem,&
                                             ptable
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

    INTEGER                                  :: iw, run_type_id
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(section_vals_type), POINTER         :: atom_section

    failure=.FALSE.

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
    INTEGER                                  :: do_eri, i, ierr, in, k, maxl, &
                                                mb, mo, n_rep, zval, zz
    INTEGER, DIMENSION(0:3)                  :: maxn
    INTEGER, DIMENSION(:), POINTER           :: cn
    LOGICAL                                  :: eri_c, eri_e, failure, &
                                                had_ae, had_pp, pp_calc
    TYPE(atom_basis_type), TARGET            :: ae_basis, pp_basis
    TYPE(atom_integrals), TARGET             :: ae_int, pp_int
    TYPE(atom_orbitals), POINTER             :: orbitals
    TYPE(atom_p_type), DIMENSION(:), POINTER :: atom_info
    TYPE(atom_state), POINTER                :: state
    TYPE(section_vals_type), POINTER         :: basis_section

    failure = .FALSE.

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

    basis_section => section_vals_get_subs_vals(atom_section,"AE_BASIS",error=error)
    NULLIFY(ae_basis%grid)
    CALL init_atom_basis(ae_basis,basis_section,zval,error)
    NULLIFY(pp_basis%grid)
    basis_section => section_vals_get_subs_vals(atom_section,"PP_BASIS",error=error)
    CALL init_atom_basis(pp_basis,basis_section,zval,error)

    eri_c = .FALSE.
    CALL section_vals_val_get(atom_section,"COULOMB_INTEGRALS", i_val=do_eri, error=error)
    IF(do_eri==do_analytic) eri_c = .TRUE.
    eri_e = .FALSE.
    CALL section_vals_val_get(atom_section,"EXCHANGE_INTEGRALS", i_val=do_eri, error=error)
    IF(do_eri==do_analytic) eri_e = .TRUE.

    CALL section_vals_val_get(atom_section,"MAX_ANGULAR_MOMENTUM", i_val=maxl, error=error)
    maxn=0
    CALL section_vals_val_get(atom_section,"CALCULATE_STATES", i_vals=cn, error=error)
    DO in = 1, MIN(SIZE(cn),4)
      maxn(in-1) = cn(in)
    END DO

    ae_int%status = 0
    pp_int%status = 0
    had_ae = .FALSE.
    had_pp = .FALSE.

    ! Check for the total number of electron configurations to be calculated
    CALL section_vals_val_get(atom_section,"ELECTRON_CONFIGURATION", n_rep_val=n_rep, error=error)

    ALLOCATE(atom_info(n_rep),STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

    DO in = 1, n_rep

      NULLIFY(atom_info(in)%atom)
      CALL create_atom_type(atom_info(in)%atom,error)
      atom_info(in)%atom%z = zval

      ALLOCATE(state,STAT=ierr)
      CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

      ! get the electronic configuration
      CALL section_vals_val_get(atom_section,"ELECTRON_CONFIGURATION", i_rep_val=in,&
                                c_vals=tmpstringlist, error=error)

      ! set occupations 
      CALL atom_set_occupation(tmpstringlist,state%occ,error)
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
        CALL atom_set_occupation(tmpstringlist,state%core,error)
        atom_info(in)%atom%zcore = zval - SUM(state%core)
      ELSE
        state%core=0._dp
        atom_info(in)%atom%zcore = zval
      END IF

      ! calculate integrals
      IF ( pp_calc ) THEN
        CALL atom_int_setup(pp_int,pp_basis,eri_coulomb=eri_c,eri_exchange=eri_e,error=error)
        CALL set_atom(atom_info(in)%atom,basis=pp_basis,integrals=pp_int,error=error)
        state%maxn_calc(:) = MIN( state%maxn_calc(:), pp_basis%nbas(:) )
        CPPostcondition(ALL(state%maxn_calc(:) >= state%maxn_occ), cp_failure_level, routineP, error, failure)
        had_pp = .TRUE.
      ELSE
        CALL atom_int_setup(ae_int,ae_basis,eri_coulomb=eri_c,eri_exchange=eri_e,error=error)
        CALL set_atom(atom_info(in)%atom,basis=ae_basis,integrals=ae_int,error=error)
        state%maxn_calc(:) = MIN( state%maxn_calc(:), ae_basis%nbas(:) )
        CPPostcondition(ALL(state%maxn_calc(:) >= state%maxn_occ), cp_failure_level, routineP, error, failure)
        had_ae = .TRUE.
      END IF

      CALL set_atom(atom_info(in)%atom,state=state,error=error)

      NULLIFY(orbitals)
      mo = MAXVAL(state%maxn_calc)
      mb = MAXVAL(atom_info(in)%atom%basis%nbas)
      CALL create_atom_orbs(orbitals,mb,mo,error)
      CALL set_atom(atom_info(in)%atom,orbitals=orbitals,error=error)

      !Calculate the electronic structure
      CALL atom_electronic_structure(atom_info(in)%atom,atom_section,error)

    END DO

    ! clean up
    IF ( had_ae ) THEN
      CALL atom_int_release(ae_int,error)
    END IF
    IF ( had_pp ) THEN
      CALL atom_int_release(pp_int,error)
    END IF
    CALL release_atom_basis(ae_basis,error)
    CALL release_atom_basis(pp_basis,error)

    DO in = 1, n_rep
      CALL release_atom_type(atom_info(in)%atom,error)
    END DO
    DEALLOCATE(atom_info,STAT=ierr)
    CPPostcondition(ierr==0, cp_failure_level, routineP, error, failure)

  END SUBROUTINE atom_energy

! *****************************************************************************
  SUBROUTINE atom_electronic_structure(atom,atom_section,error)
    TYPE(atom_type), POINTER                 :: atom
    TYPE(section_vals_type), POINTER         :: atom_section
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'atom_electronic_structure', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: do_eri, i, iter, iw, l, method
    LOGICAL                                  :: failure, need_x, need_xc
    REAL(KIND=dp)                            :: deps, eps_scf
    TYPE(atom_history_type)                  :: history
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(opgrid_type), POINTER               :: cpot, density
    TYPE(opmat_type), POINTER                :: fmat, hcore, jmat, kmat, xcmat
    TYPE(section_vals_type), POINTER         :: method_section, opt_section

    failure = .FALSE.

    method_section => section_vals_get_subs_vals(atom_section,"METHOD",error=error)
    opt_section => section_vals_get_subs_vals(atom_section,"OPTIMIZATION",error=error)
    CALL section_vals_val_get(method_section,"METHOD_TYPE", i_val=method, error=error)

    SELECT CASE (method)
      CASE DEFAULT
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_rks_atom)
        need_x = .FALSE.
        need_xc = .TRUE.
      CASE (do_uks_atom)
        need_x = .FALSE.
        need_xc = .TRUE.
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_rhf_atom)
        need_x = .TRUE.
        need_xc = .FALSE.
      CASE (do_uhf_atom)
        need_x = .TRUE.
        need_xc = .FALSE.
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_rohf_atom)
        need_x = .TRUE.
        need_xc = .FALSE.
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
    END SELECT

    logger => cp_error_get_logger(error)
    iw = cp_print_key_unit_nr(logger,atom_section,"PRINT%SCF_INFO",extension=".log",error=error)
    IF (iw>0) THEN
      CALL atom_print_state(atom%state,iw,error)
    END IF

    eps_scf = 1.e-8_dp

    NULLIFY(hcore)
    CALL create_opmat(hcore,atom%basis%nbas,error)
    hcore%op = atom%integrals%kin - atom%zcore*atom%integrals%core
    ! Pseudopotentials
!   SELECT CASE (atom%ppot_type)
!     CASE DEFAULT
!       CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
!     CASE
!   END SELECT

    ! initial guess
    CALL atom_solve(hcore%op,atom%integrals%ovlp,atom%orbitals%wfn,atom%orbitals%ener,&
                    atom%basis%nbas,atom%state%maxl_calc,error)
    CALL atom_denmat(atom%orbitals%pmat,atom%orbitals%wfn,atom%basis%nbas,atom%state%occ,&
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
        DO i=1,MIN(SIZE(atom%state%occ,2),SIZE(atom%orbitals%ener,1))
          atom%energy%eband = atom%energy%eband + atom%orbitals%ener(i,l)*atom%state%occ(l,i)
        END DO
      END DO

      ! Pseudopotential energy
      atom%energy%eploc = 0._dp
      atom%energy%epnl = 0._dp
!     SELECT CASE (atom%ppot_type)
!       CASE DEFAULT
!         CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
!       CASE
!     END SELECT
      atom%energy%epseudo = atom%energy%eploc + atom%energy%epnl
  
      ! Core energy
      atom%energy%ecore = atom_trace(hcore%op,atom%orbitals%pmat,error)

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
          CALL atom_density(density%op,atom%orbitals%pmat,atom%basis,atom%state%maxl_occ,"RHO",error)
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
        atom%energy%eexchange = 0.5_dp * atom_trace(kmat%op,atom%orbitals%pmat,error)
      ELSE
        kmat%op = 0._dp
        atom%energy%eexchange = 0._dp
      END IF

      ! XC
      IF (need_xc) THEN
        xcmat%op = 0._dp
        CALL calculate_atom_vxc(xcmat,atom,atom_section,error)
      ELSE
        xcmat%op = 0._dp
        atom%energy%exc = 0._dp
      END IF
     
      ! Total energy
      atom%energy%etot = atom%energy%ecore + atom%energy%ecoulomb + atom%energy%eexchange + atom%energy%exc

      ! Total HF/KS matrix
      fmat%op = hcore%op + jmat%op + kmat%op + xcmat%op

      ! calculate error matrix
      CALL err_matrix(jmat%op,deps,fmat%op,atom%orbitals%pmat,atom%integrals%ovlp,atom%basis%nbas,error)

      iter = iter + 1

      CALL atom_print_iteration(iter,deps,atom%energy%etot,iw,error)

      IF ( deps < eps_scf ) EXIT

      ! update history container and extrapolate KS matrix
      CALL atom_history_update (history,fmat%op,jmat%op,error)
      CALL atom_opt (fmat%op,history,deps,error)

      ! Solve HF/KS equations
      CALL atom_solve(fmat%op,atom%integrals%ovlp,atom%orbitals%wfn,atom%orbitals%ener,&
                      atom%basis%nbas,atom%state%maxl_calc,error)
      CALL atom_denmat(atom%orbitals%pmat,atom%orbitals%wfn,atom%basis%nbas,atom%state%occ,&
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

    CALL cp_print_key_finished_output(iw,logger,atom_section,"PRINT%PROGRAM_BANNER",error=error)

  END SUBROUTINE atom_electronic_structure

! *****************************************************************************

END MODULE atom
