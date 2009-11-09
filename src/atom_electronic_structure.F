!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
MODULE atom_electronic_structure
  USE atom_optimization,               ONLY: atom_history_init,&
                                             atom_history_release,&
                                             atom_history_type,&
                                             atom_history_update,&
                                             atom_opt
  USE atom_output,                     ONLY: atom_print_energies,&
                                             atom_print_iteration,&
                                             atom_print_state
  USE atom_types,                      ONLY: &
       GTH_PSEUDO, NO_PSEUDO, atom_type, create_opgrid, create_opmat, &
       opgrid_type, opmat_type, release_opgrid, release_opmat
  USE atom_utils,                      ONLY: &
       atom_denmat, atom_density, atom_solve, atom_trace, ceri_contract, &
       coulomb_potential_analytic, coulomb_potential_numeric, eeri_contract, &
       err_matrix, exchange_numeric, exchange_semi_analytic, numpot_matrix, &
       slater_density, wigner_slater_functional
  USE atom_xc,                         ONLY: calculate_atom_vxc
  USE f77_blas
  USE input_constants,                 ONLY: &
       do_analytic, do_dkh0_atom, do_dkh1_atom, do_dkh2_atom, do_dkh3_atom, &
       do_dkh4_atom, do_dkh5_atom, do_nonrel_atom, do_numeric, do_rhf_atom, &
       do_rks_atom, do_rohf_atom, do_semi_analytic, do_uhf_atom, do_uks_atom, &
       do_zoramp_atom
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: calculate_atom

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_electronic_structure'

! *****************************************************************************

CONTAINS

! *****************************************************************************
  SUBROUTINE calculate_atom(atom,iw,noguess,error)
    TYPE(atom_type), POINTER                 :: atom
    INTEGER, INTENT(IN)                      :: iw
    LOGICAL, INTENT(IN), OPTIONAL            :: noguess
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_atom', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, method
    LOGICAL                                  :: failure

    failure = .FALSE.

    CALL timeset(routineN,handle)

    method   = atom%method_type

    SELECT CASE (method)
      CASE DEFAULT
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_rks_atom)
        CALL calculate_atom_restricted(atom,iw,noguess,error)
      CASE (do_uks_atom)
        CALL calculate_atom_unrestricted(atom,iw,noguess,error)
      CASE (do_rhf_atom)
        CALL calculate_atom_restricted(atom,iw,noguess,error)
      CASE (do_uhf_atom)
        CALL calculate_atom_unrestricted(atom,iw,noguess,error)
      CASE (do_rohf_atom)
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE calculate_atom

! *****************************************************************************
  SUBROUTINE calculate_atom_restricted(atom,iw,noguess,error)
    TYPE(atom_type), POINTER                 :: atom
    INTEGER, INTENT(IN)                      :: iw
    LOGICAL, INTENT(IN), OPTIONAL            :: noguess
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_atom_restricted', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iter, l, max_iter, &
                                                method, reltyp
    LOGICAL                                  :: do_hfx, doguess, failure, &
                                                need_x, need_xc
    REAL(KIND=dp)                            :: deps, eps_scf, hf_frac
    TYPE(atom_history_type)                  :: history
    TYPE(opgrid_type), POINTER               :: cpot, density
    TYPE(opmat_type), POINTER                :: fmat, hcore, jmat, kmat, xcmat
    TYPE(section_vals_type), POINTER         :: hfx_sections, xc_section

    failure = .FALSE.

    CALL timeset(routineN,handle)

    IF (PRESENT(noguess)) THEN
      doguess=.NOT.noguess
    ELSE
      doguess=.TRUE.
    END IF

    hf_frac = 0._dp
    IF (ASSOCIATED(atom%xc_section)) THEN
      xc_section => atom%xc_section
      hfx_sections => section_vals_get_subs_vals(xc_section,"HF",error=error)
      CALL section_vals_get(hfx_sections,explicit=do_hfx,error=error)
      IF ( do_hfx ) THEN
        CALL section_vals_val_get(hfx_sections,"FRACTION", r_val=hf_frac, error=error)
      END IF
    ELSE
      NULLIFY(xc_section)
      do_hfx = .FALSE.
    END IF

    method   = atom%method_type
    max_iter = atom%optimization%max_iter
    eps_scf  = atom%optimization%eps_scf 

    SELECT CASE (method)
      CASE DEFAULT
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_rks_atom)
        need_x = do_hfx
        need_xc = .TRUE.
      CASE (do_uks_atom)
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_rhf_atom)
        need_x = .TRUE.
        need_xc = .FALSE.
        hf_frac = 1._dp
      CASE (do_uhf_atom)
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_rohf_atom)
        need_x = .TRUE.
        need_xc = .FALSE.
        hf_frac = 1._dp
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
    END SELECT

    ! check for relativistic method
    reltyp = atom%relativistic

    IF (iw>0) CALL atom_print_state(atom%state,iw,error)

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
          CASE (do_dkh0_atom, do_dkh1_atom,do_dkh2_atom,do_dkh3_atom,do_dkh4_atom,do_dkh5_atom)
            hcore%op = atom%integrals%hdkh
        END SELECT
      CASE (GTH_PSEUDO)
        hcore%op = atom%integrals%kin + atom%integrals%core + atom%integrals%hnl
    END SELECT
    ! add confinement potential (not included in relativistic transformations)
    IF ( atom%potential%confinement ) THEN
      hcore%op = hcore%op + atom%potential%acon * atom%integrals%conf
    END IF

    NULLIFY(fmat,jmat,kmat,xcmat)
    CALL create_opmat(fmat,atom%basis%nbas,error)
    CALL create_opmat(jmat,atom%basis%nbas,error)
    CALL create_opmat(kmat,atom%basis%nbas,error)
    CALL create_opmat(xcmat,atom%basis%nbas,error)

    NULLIFY(density,cpot)
    CALL create_opgrid(density,atom%basis%grid,error)
    CALL create_opgrid(cpot,atom%basis%grid,error)

    IF (doguess) THEN
      ! initial guess
      CALL slater_density(density%op,atom%z,atom%state,atom%basis%grid,error)
      CALL coulomb_potential_numeric(cpot%op,density%op,density%grid,error)
      CALL numpot_matrix(jmat%op,cpot%op,atom%basis,0,error)
      CALL wigner_slater_functional(density%op,cpot%op,error)
      CALL numpot_matrix(xcmat%op,cpot%op,atom%basis,0,error)
      fmat%op = hcore%op + jmat%op + xcmat%op
      CALL atom_solve(fmat%op,atom%integrals%utrans,atom%orbitals%wfn,atom%orbitals%ener,&
                      atom%basis%nbas,atom%integrals%nne,atom%state%maxl_calc,error)
    END IF
    CALL atom_denmat(atom%orbitals%pmat,atom%orbitals%wfn,atom%basis%nbas,atom%state%occupation,&
                     atom%state%maxl_occ,atom%state%maxn_occ,error)

    ! wavefunction history
    NULLIFY(history%dmat,history%hmat)
    CALL atom_history_init (history,atom%optimization,fmat%op,error)

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
      SELECT CASE (atom%coulomb_integral_type)
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
        SELECT CASE (atom%exchange_integral_type)
          CASE DEFAULT
            CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
          CASE (do_analytic)
            CALL eeri_contract(kmat%op,atom%integrals%eeri,atom%orbitals%pmat,atom%integrals%n,error=error)
          CASE (do_semi_analytic)
            CALL exchange_semi_analytic(kmat%op,atom%state,atom%orbitals%wfn,atom%basis,error=error)
          CASE (do_numeric)
            CALL exchange_numeric(kmat%op,atom%state,atom%orbitals%wfn,atom%basis,error=error)
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
        CALL calculate_atom_vxc(xcmat,atom,xc_section,error)
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

      IF ( iw > 0 ) THEN
        CALL atom_print_iteration(iter,deps,atom%energy%etot,iw,error)
      END IF

      IF ( deps < eps_scf ) EXIT
      IF ( iter >= max_iter ) THEN
        IF ( iw > 0 ) THEN
          WRITE(iw,"(A)") " No convergence within maximum number of iterations "
        END IF
        EXIT
      END IF

      ! update history container and extrapolate KS matrix
      CALL atom_history_update (history,fmat%op,jmat%op,error)
      CALL atom_opt (fmat%op,history,deps,error)

      ! Solve HF/KS equations
      CALL atom_solve(fmat%op,atom%integrals%utrans,atom%orbitals%wfn,atom%orbitals%ener,&
                      atom%basis%nbas,atom%integrals%nne,atom%state%maxl_calc,error)
      CALL atom_denmat(atom%orbitals%pmat,atom%orbitals%wfn,atom%basis%nbas,atom%state%occupation,&
                       atom%state%maxl_occ,atom%state%maxn_occ,error)

    END DO        !SCF Loop

    IF ( iw > 0 ) THEN
      CALL atom_print_energies(atom,iw,error)
    END IF

    CALL atom_history_release(history,error)

    CALL release_opmat(fmat,error)
    CALL release_opmat(jmat,error)
    CALL release_opmat(kmat,error)
    CALL release_opmat(xcmat,error)
    CALL release_opmat(hcore,error)

    CALL release_opgrid(density,error)
    CALL release_opgrid(cpot,error)

    CALL timestop(handle)

  END SUBROUTINE calculate_atom_restricted

! *****************************************************************************
  SUBROUTINE calculate_atom_unrestricted(atom,iw,noguess,error)
    TYPE(atom_type), POINTER                 :: atom
    INTEGER, INTENT(IN)                      :: iw
    LOGICAL, INTENT(IN), OPTIONAL            :: noguess
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_atom_unrestricted', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iter, k, l, &
                                                max_iter, method, reltyp
    LOGICAL                                  :: do_hfx, doguess, failure, &
                                                need_x, need_xc
    REAL(KIND=dp)                            :: deps, eps_scf, hf_frac, ne, nm
    TYPE(atom_history_type)                  :: history
    TYPE(opgrid_type), POINTER               :: cpot, density
    TYPE(opmat_type), POINTER                :: fmat, hcore, jmat, kmat, xcmat
    TYPE(section_vals_type), POINTER         :: hfx_sections, xc_section

    failure = .FALSE.

    CALL timeset(routineN,handle)

    IF (PRESENT(noguess)) THEN
      doguess=.NOT.noguess
    ELSE
      doguess=.TRUE.
    END IF

    hf_frac = 0._dp
    IF (ASSOCIATED(atom%xc_section)) THEN
      xc_section => atom%xc_section
      hfx_sections => section_vals_get_subs_vals(xc_section,"HF",error=error)
      CALL section_vals_get(hfx_sections,explicit=do_hfx,error=error)
      IF ( do_hfx ) THEN
        CALL section_vals_val_get(hfx_sections,"FRACTION", r_val=hf_frac, error=error)
      END IF
    ELSE
      NULLIFY(xc_section)
      do_hfx = .FALSE.
    END IF

    method   = atom%method_type
    max_iter = atom%optimization%max_iter
    eps_scf  = atom%optimization%eps_scf 

    SELECT CASE (method)
      CASE DEFAULT
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_rks_atom)
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_uks_atom)
        need_x = do_hfx
        need_xc = .TRUE.
      CASE (do_rhf_atom)
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
      CASE (do_uhf_atom)
        need_x = .TRUE.
        need_xc = .FALSE.
        hf_frac = 1._dp
      CASE (do_rohf_atom)
        need_x = .TRUE.
        need_xc = .FALSE.
        hf_frac = 1._dp
        CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
    END SELECT

    ! set alpha and beta occupations
    atom%state%occa = 0._dp
    atom%state%occb = 0._dp
    DO l=0,3
      nm=REAL((2*l+1),KIND=dp)
      DO k=1,10
        ne = atom%state%occupation(k,l)
        IF(ne==0._dp) THEN  !empty shell
          EXIT  !assume there are no holes
        ELSEIF(ne==2._dp*nm) THEN  !closed shell
          atom%state%occa(k,l)=nm
          atom%state%occb(k,l)=nm
        ELSEIF(atom%state%multiplicity==-2) THEN !High spin case
          atom%state%occa(k,l)=MIN(ne,nm)
          atom%state%occb(k,l)=MAX(0._dp,ne-nm)
        ELSE
          atom%state%occa(k,l)=0.5_dp*(ne+atom%state%multiplicity-1._dp)
          atom%state%occb(k,l)=ne-atom%state%occa(k,l)
        END IF
      END DO
    END DO
    ! check for relativistic method
    reltyp = atom%relativistic

    IF (iw>0) CALL atom_print_state(atom%state,iw,error)

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
          CASE (do_dkh0_atom, do_dkh1_atom,do_dkh2_atom,do_dkh3_atom,do_dkh4_atom,do_dkh5_atom)
            hcore%op = atom%integrals%hdkh
        END SELECT
      CASE (GTH_PSEUDO)
        hcore%op = atom%integrals%kin + atom%integrals%core + atom%integrals%hnl
    END SELECT
    ! add confinement potential (not included in relativistic transformations)
    IF ( atom%potential%confinement ) THEN
      hcore%op = hcore%op + atom%potential%acon * atom%integrals%conf
    END IF

!uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu
    NULLIFY(fmat,jmat,kmat,xcmat)
    CALL create_opmat(fmat,atom%basis%nbas,error)
    CALL create_opmat(jmat,atom%basis%nbas,error)
    CALL create_opmat(kmat,atom%basis%nbas,error)
    CALL create_opmat(xcmat,atom%basis%nbas,error)

    NULLIFY(density,cpot)
    CALL create_opgrid(density,atom%basis%grid,error)
    CALL create_opgrid(cpot,atom%basis%grid,error)

    IF (doguess) THEN
      ! initial guess
      CALL slater_density(density%op,atom%z,atom%state,atom%basis%grid,error)
      CALL coulomb_potential_numeric(cpot%op,density%op,density%grid,error)
      CALL numpot_matrix(jmat%op,cpot%op,atom%basis,0,error)
      CALL wigner_slater_functional(density%op,cpot%op,error)
      CALL numpot_matrix(xcmat%op,cpot%op,atom%basis,0,error)
      fmat%op = hcore%op + jmat%op + xcmat%op
      CALL atom_solve(fmat%op,atom%integrals%utrans,atom%orbitals%wfn,atom%orbitals%ener,&
                      atom%basis%nbas,atom%integrals%nne,atom%state%maxl_calc,error)
    END IF
    CALL atom_denmat(atom%orbitals%pmat,atom%orbitals%wfn,atom%basis%nbas,atom%state%occupation,&
                     atom%state%maxl_occ,atom%state%maxn_occ,error)

    ! wavefunction history
    NULLIFY(history%dmat,history%hmat)
    CALL atom_history_init (history,atom%optimization,fmat%op,error)

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
      SELECT CASE (atom%coulomb_integral_type)
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
        SELECT CASE (atom%exchange_integral_type)
          CASE DEFAULT
            CPPostcondition(.FALSE., cp_failure_level, routineP, error, failure)
          CASE (do_analytic)
            CALL eeri_contract(kmat%op,atom%integrals%eeri,atom%orbitals%pmat,atom%integrals%n,error=error)
          CASE (do_semi_analytic)
            CALL exchange_semi_analytic(kmat%op,atom%state,atom%orbitals%wfn,atom%basis,error=error)
          CASE (do_numeric)
            CALL exchange_numeric(kmat%op,atom%state,atom%orbitals%wfn,atom%basis,error=error)
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
        CALL calculate_atom_vxc(xcmat,atom,xc_section,error)
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

      IF ( iw > 0 ) THEN
        CALL atom_print_iteration(iter,deps,atom%energy%etot,iw,error)
      END IF

      IF ( deps < eps_scf ) EXIT
      IF ( iter >= max_iter ) THEN
        IF ( iw > 0 ) THEN
          WRITE(iw,"(A)") " No convergence within maximum number of iterations "
        END IF
        EXIT
      END IF

      ! update history container and extrapolate KS matrix
      CALL atom_history_update (history,fmat%op,jmat%op,error)
      CALL atom_opt (fmat%op,history,deps,error)

      ! Solve HF/KS equations
      CALL atom_solve(fmat%op,atom%integrals%utrans,atom%orbitals%wfn,atom%orbitals%ener,&
                      atom%basis%nbas,atom%integrals%nne,atom%state%maxl_calc,error)
      CALL atom_denmat(atom%orbitals%pmat,atom%orbitals%wfn,atom%basis%nbas,atom%state%occupation,&
                       atom%state%maxl_occ,atom%state%maxn_occ,error)

    END DO        !SCF Loop

    IF ( iw > 0 ) THEN
      CALL atom_print_energies(atom,iw,error)
    END IF

    CALL atom_history_release(history,error)

    CALL release_opmat(fmat,error)
    CALL release_opmat(jmat,error)
    CALL release_opmat(kmat,error)
    CALL release_opmat(xcmat,error)
    CALL release_opmat(hcore,error)

    CALL release_opgrid(density,error)
    CALL release_opgrid(cpot,error)

    CALL timestop(handle)

  END SUBROUTINE calculate_atom_unrestricted

! *****************************************************************************

END MODULE atom_electronic_structure
