!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE atom_electronic_structure
   USE atom_optimization,               ONLY: atom_history_init,&
                                              atom_history_release,&
                                              atom_history_type,&
                                              atom_history_update,&
                                              atom_opt_fmat
   USE atom_output,                     ONLY: atom_print_energies,&
                                              atom_print_iteration,&
                                              atom_print_state,&
                                              atom_print_zmp_iteration
   USE atom_types,                      ONLY: &
        atom_type, create_opgrid, create_opmat, ecp_pseudo, gth_pseudo, lmat, no_pseudo, &
        opgrid_type, opmat_type, release_opgrid, release_opmat, set_atom, setup_hf_section, &
        sgp_pseudo, upf_pseudo
   USE atom_utils,                      ONLY: &
        atom_denmat, atom_density, atom_read_external_density, atom_read_external_vxc, &
        atom_read_zmp_restart, atom_solve, atom_trace, ceri_contract, coulomb_potential_analytic, &
        coulomb_potential_numeric, eeri_contract, err_matrix, exchange_numeric, &
        exchange_semi_analytic, numpot_matrix, slater_density, wigner_slater_functional
   USE atom_xc,                         ONLY: calculate_atom_ext_vxc,&
                                              calculate_atom_vxc_lda,&
                                              calculate_atom_vxc_lsd,&
                                              calculate_atom_zmp
   USE input_constants,                 ONLY: &
        do_analytic, do_dkh0_atom, do_dkh1_atom, do_dkh2_atom, do_dkh3_atom, do_nonrel_atom, &
        do_numeric, do_rhf_atom, do_rks_atom, do_rohf_atom, do_sczoramp_atom, do_semi_analytic, &
        do_uhf_atom, do_uks_atom, do_zoramp_atom
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC  :: calculate_atom

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom_electronic_structure'

CONTAINS

! **************************************************************************************************
!> \brief General routine to perform electronic structure atomic calculations.
!> \param atom       information about the atomic kind
!> \param iw         output file unit
!> \param noguess    skip initial guess
!> \param converged  whether SCF iterations have been converged
!> \par History
!>    * 06.2017 disable XC input options [Juerg Hutter]
!>    * 11.2009 created from the subroutine calculate_atom() [Juerg Hutter]
! **************************************************************************************************
   SUBROUTINE calculate_atom(atom, iw, noguess, converged)
      TYPE(atom_type), POINTER                           :: atom
      INTEGER, INTENT(IN)                                :: iw
      LOGICAL, INTENT(IN), OPTIONAL                      :: noguess
      LOGICAL, INTENT(OUT), OPTIONAL                     :: converged

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_atom', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, method
      LOGICAL                                            :: explicit
      TYPE(section_vals_type), POINTER                   :: sub_section

      CALL timeset(routineN, handle)

      ! test for not supported methods and functionals
      IF (ASSOCIATED(atom%xc_section)) THEN
         !
         sub_section => section_vals_get_subs_vals(atom%xc_section, "ADIABATIC_RESCALING")
         CALL section_vals_get(sub_section, explicit=explicit)
         IF (explicit) CALL cp_abort(__LOCATION__, "ADIABATIC_RESCALING not supported in ATOM code")
         !
         sub_section => section_vals_get_subs_vals(atom%xc_section, "VDW_POTENTIAL")
         CALL section_vals_get(sub_section, explicit=explicit)
         IF (explicit) CALL cp_abort(__LOCATION__, "VDW_POTENTIAL not supported in ATOM code")
         !
         sub_section => section_vals_get_subs_vals(atom%xc_section, "XC_POTENTIAL")
         CALL section_vals_get(sub_section, explicit=explicit)
         IF (explicit) CALL cp_abort(__LOCATION__, "XC_POTENTIAL not supported in ATOM code")
         !
         sub_section => section_vals_get_subs_vals(atom%xc_section, "WF_CORRELATION")
         CALL section_vals_get(sub_section, explicit=explicit)
         IF (explicit) CALL cp_abort(__LOCATION__, "WF_CORRELATION methods not supported in ATOM code")
         !
      END IF

      method = atom%method_type
      SELECT CASE (method)
      CASE (do_rks_atom, do_rhf_atom)
         CALL calculate_atom_restricted(atom, iw, noguess, converged)
      CASE (do_uks_atom, do_uhf_atom)
         CALL calculate_atom_unrestricted(atom, iw, noguess, converged)
      CASE (do_rohf_atom)
         CPABORT("")
      CASE DEFAULT
         CPABORT("")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE calculate_atom

! **************************************************************************************************
!> \brief Perform restricted (closed shell) electronic structure atomic calculations.
!> \param atom       information about the atomic kind
!> \param iw         output file unit
!> \param noguess    skip initial guess
!> \param converged  whether SCF iterations have been converged
!> \par History
!>    * 02.2016 support UPF files and ECP potentials [Juerg Hutter]
!>    * 09.2015 Polarized Atomic Orbitals (PAO) method [Ole Schuett]
!>    * 11.2013 Zhao, Morrison, and Parr (ZMP) potential [Daniele Varsano]
!>    * 07.2013 scaled ZORA with model potential [Juerg Hutter]
!>    * 11.2009 split into three subroutines calculate_atom(), calculate_atom_restricted(),
!>              and calculate_atom_unrestricted() [Juerg Hutter]
!>    * 12.2008 refactored and renamed to calculate_atom() [Juerg Hutter]
!>    * 08.2008 created as atom_electronic_structure() [Juerg Hutter]
! **************************************************************************************************
   SUBROUTINE calculate_atom_restricted(atom, iw, noguess, converged)
      TYPE(atom_type), POINTER                           :: atom
      INTEGER, INTENT(IN)                                :: iw
      LOGICAL, INTENT(IN), OPTIONAL                      :: noguess
      LOGICAL, INTENT(OUT), OPTIONAL                     :: converged

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_atom_restricted', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iter, l, max_iter, method, &
                                                            reltyp
      LOGICAL                                            :: do_hfx, doguess, is_converged, need_vxc, &
                                                            need_x, need_xc, need_zmp
      REAL(KIND=dp)                                      :: deps, eps_scf, hf_frac
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ext_density, ext_vxc, tmp_dens
      TYPE(atom_history_type)                            :: history
      TYPE(opgrid_type), POINTER                         :: cpot, density
      TYPE(opmat_type), POINTER                          :: fmat, hcore, jmat, kmat, xcmat
      TYPE(section_vals_type), POINTER                   :: xc_section

      CALL timeset(routineN, handle)

      IF (PRESENT(noguess)) THEN
         doguess = .NOT. noguess
      ELSE
         doguess = .TRUE.
      END IF

      CALL setup_hf_section(hf_frac, do_hfx, atom, xc_section, atom%exchange_integral_type)

      method = atom%method_type
      max_iter = atom%optimization%max_iter
      eps_scf = atom%optimization%eps_scf

      SELECT CASE (method)
      CASE DEFAULT
         CPABORT("")
      CASE (do_rks_atom)
         need_x = do_hfx
         need_xc = .TRUE.
      CASE (do_uks_atom)
         CPABORT("")
      CASE (do_rhf_atom)
         need_x = .TRUE.
         need_xc = .FALSE.
         hf_frac = 1._dp
      CASE (do_uhf_atom)
         CPABORT("")
      CASE (do_rohf_atom)
         need_x = .TRUE.
         need_xc = .FALSE.
         hf_frac = 1._dp
         CPABORT("")
      END SELECT

      ! ZMP starting to read external density for the zmp calculation
      need_zmp = atom%do_zmp
      IF (need_zmp) THEN
         need_x = .FALSE.
         need_xc = .FALSE.
         ALLOCATE (ext_density(atom%basis%grid%nr))
         ext_density = 0._dp
         CALL atom_read_external_density(ext_density, atom, iw)
      ENDIF

      ! ZMP starting to read external vxc for the zmp calculation
      need_vxc = atom%read_vxc

      IF (need_vxc) THEN
         need_x = .FALSE.
         need_xc = .FALSE.
         need_zmp = .FALSE.
         ALLOCATE (ext_vxc(atom%basis%grid%nr))
         ext_vxc = 0._dp
         CALL atom_read_external_vxc(ext_vxc, atom, iw)
      ENDIF

      ! check for relativistic method
      reltyp = atom%relativistic

      IF (iw > 0) CALL atom_print_state(atom%state, iw)

      NULLIFY (hcore)
      CALL create_opmat(hcore, atom%basis%nbas)
      ! Pseudopotentials
      SELECT CASE (atom%potential%ppot_type)
      CASE DEFAULT
         CPABORT("")
      CASE (NO_PSEUDO)
         SELECT CASE (reltyp)
         CASE DEFAULT
            CPABORT("")
         CASE (do_nonrel_atom)
            hcore%op = atom%integrals%kin - atom%zcore*atom%integrals%core
         CASE (do_zoramp_atom, do_sczoramp_atom)
            hcore%op = atom%integrals%kin + atom%integrals%tzora - atom%zcore*atom%integrals%core
         CASE (do_dkh0_atom, do_dkh1_atom, do_dkh2_atom, do_dkh3_atom)
            hcore%op = atom%integrals%hdkh
         END SELECT
      CASE (gth_pseudo, upf_pseudo, sgp_pseudo, ecp_pseudo)
         hcore%op = atom%integrals%kin + atom%integrals%core + atom%integrals%hnl
      END SELECT
      ! add confinement potential (not included in relativistic transformations)
      IF (atom%potential%confinement) THEN
         hcore%op = hcore%op + atom%potential%acon*atom%integrals%conf
      END IF

      NULLIFY (fmat, jmat, kmat, xcmat)
      CALL create_opmat(fmat, atom%basis%nbas)
      CALL create_opmat(jmat, atom%basis%nbas)
      CALL create_opmat(kmat, atom%basis%nbas)
      CALL create_opmat(xcmat, atom%basis%nbas)

      NULLIFY (density, cpot)
      CALL create_opgrid(density, atom%basis%grid)
      CALL create_opgrid(cpot, atom%basis%grid)

      ! ZMP reading the file to restart
      IF (atom%doread) CALL atom_read_zmp_restart(atom, doguess, iw)

      IF (doguess) THEN
         ! initial guess
         ALLOCATE (tmp_dens(SIZE(density%op)))
         CALL slater_density(density%op, tmp_dens, atom%z, atom%state, atom%basis%grid)
         density%op = density%op + tmp_dens
         DEALLOCATE (tmp_dens)
         CALL coulomb_potential_numeric(cpot%op, density%op, density%grid)
         CALL numpot_matrix(jmat%op, cpot%op, atom%basis, 0)
         CALL wigner_slater_functional(density%op, cpot%op)
         CALL numpot_matrix(xcmat%op, cpot%op, atom%basis, 0)
         fmat%op = hcore%op + jmat%op + xcmat%op
         CALL atom_solve(fmat%op, atom%integrals%utrans, atom%orbitals%wfn, atom%orbitals%ener, &
                         atom%basis%nbas, atom%integrals%nne, atom%state%maxl_calc)
      END IF
      CALL atom_denmat(atom%orbitals%pmat, atom%orbitals%wfn, atom%basis%nbas, atom%state%occupation, &
                       atom%state%maxl_occ, atom%state%maxn_occ)

      ! wavefunction history
      NULLIFY (history%dmat, history%hmat)
      CALL atom_history_init(history, atom%optimization, fmat%op)
      is_converged = .FALSE.
      iter = 0
      DO !SCF Loop

         ! Kinetic energy
         atom%energy%ekin = atom_trace(atom%integrals%kin, atom%orbitals%pmat)

         ! Band energy
         atom%energy%eband = 0._dp
         DO l = 0, lmat
            DO i = 1, MIN(SIZE(atom%state%occupation, 2), SIZE(atom%orbitals%ener, 1))
               atom%energy%eband = atom%energy%eband + atom%orbitals%ener(i, l)*atom%state%occupation(l, i)
            END DO
         END DO

         ! Pseudopotential energy
         SELECT CASE (atom%potential%ppot_type)
         CASE DEFAULT
            CPABORT("")
         CASE (no_pseudo)
            atom%energy%eploc = 0._dp
            atom%energy%epnl = 0._dp
         CASE (gth_pseudo, upf_pseudo, sgp_pseudo, ecp_pseudo)
            atom%energy%eploc = atom_trace(atom%integrals%core, atom%orbitals%pmat)
            atom%energy%epnl = atom_trace(atom%integrals%hnl, atom%orbitals%pmat)
         END SELECT
         atom%energy%epseudo = atom%energy%eploc + atom%energy%epnl

         ! Core energy
         atom%energy%ecore = atom_trace(hcore%op, atom%orbitals%pmat)

         ! Confinement energy
         IF (atom%potential%confinement) THEN
            atom%energy%econfinement = atom_trace(atom%integrals%conf, atom%orbitals%pmat)
         ELSE
            atom%energy%econfinement = 0._dp
         END IF

         ! Hartree Term
         jmat%op = 0._dp
         SELECT CASE (atom%coulomb_integral_type)
         CASE DEFAULT
            CPABORT("")
         CASE (do_analytic)
            CALL ceri_contract(jmat%op, atom%integrals%ceri, atom%orbitals%pmat, atom%integrals%n)
         CASE (do_semi_analytic)
            CALL coulomb_potential_analytic(cpot%op, atom%orbitals%pmat, atom%basis, atom%basis%grid, &
                                            atom%state%maxl_occ)
            CALL numpot_matrix(jmat%op, cpot%op, atom%basis, 0)
         CASE (do_numeric)
            CALL atom_density(density%op, atom%orbitals%pmat, atom%basis, atom%state%maxl_occ, typ="RHO")
            CALL coulomb_potential_numeric(cpot%op, density%op, density%grid)
            CALL numpot_matrix(jmat%op, cpot%op, atom%basis, 0)
         END SELECT
         atom%energy%ecoulomb = 0.5_dp*atom_trace(jmat%op, atom%orbitals%pmat)

         ! Exchange Term
         IF (need_x) THEN
            kmat%op = 0._dp
            SELECT CASE (atom%exchange_integral_type)
            CASE DEFAULT
               CPABORT("")
            CASE (do_analytic)
               CALL eeri_contract(kmat%op, atom%integrals%eeri, atom%orbitals%pmat, atom%integrals%n)
            CASE (do_semi_analytic)
               CALL exchange_semi_analytic(kmat%op, atom%state, atom%state%occupation, atom%orbitals%wfn, atom%basis, atom%hfx_pot)
            CASE (do_numeric)
               CALL exchange_numeric(kmat%op, atom%state, atom%state%occupation, atom%orbitals%wfn, atom%basis, atom%hfx_pot)
            END SELECT
            atom%energy%eexchange = hf_frac*0.5_dp*atom_trace(kmat%op, atom%orbitals%pmat)
            kmat%op = hf_frac*kmat%op
         ELSE
            kmat%op = 0._dp
            atom%energy%eexchange = 0._dp
         END IF

         ! XC
         IF (need_xc) THEN
            xcmat%op = 0._dp
            CALL calculate_atom_vxc_lda(xcmat, atom, xc_section)
            ! ZMP added options for the zmp calculations, building external density and vxc potential
         ELSEIF (need_zmp) THEN
            xcmat%op = 0._dp
            CALL calculate_atom_zmp(ext_density=ext_density, atom=atom, lprint=.FALSE., xcmat=xcmat)
         ELSEIF (need_vxc) THEN
            xcmat%op = 0._dp
            CALL calculate_atom_ext_vxc(vxc=ext_vxc, atom=atom, lprint=.FALSE., xcmat=xcmat)
         ELSE
            xcmat%op = 0._dp
            atom%energy%exc = 0._dp
         END IF

         ! Zero this contribution
         atom%energy%elsd = 0._dp

         ! Total energy
         atom%energy%etot = atom%energy%ecore + atom%energy%ecoulomb + atom%energy%eexchange + atom%energy%exc

         ! Potential energy
         atom%energy%epot = atom%energy%etot - atom%energy%ekin

         ! Total HF/KS matrix
         fmat%op = hcore%op + jmat%op + kmat%op + xcmat%op

         ! calculate error matrix
         CALL err_matrix(jmat%op, deps, fmat%op, atom%orbitals%pmat, atom%integrals%utrans, &
                         atom%integrals%uptrans, atom%basis%nbas, atom%integrals%nne)

         iter = iter + 1

         IF (iw > 0) THEN
            IF (need_zmp) THEN
               CALL atom_print_zmp_iteration(iter, deps, atom, iw)
            ELSE
               CALL atom_print_iteration(iter, deps, atom%energy%etot, iw)
            ENDIF
         END IF

         IF (deps < eps_scf) THEN
            is_converged = .TRUE.
            EXIT
         END IF
         IF (iter >= max_iter) THEN
            IF (iw > 0) THEN
               WRITE (iw, "(A)") " No convergence within maximum number of iterations "
            END IF
            EXIT
         END IF

         ! update history container and extrapolate KS matrix
         CALL atom_history_update(history, atom%orbitals%pmat, fmat%op, jmat%op, atom%energy%etot, deps)
         CALL atom_opt_fmat(fmat%op, history, deps)

         ! Solve HF/KS equations
         CALL atom_solve(fmat%op, atom%integrals%utrans, atom%orbitals%wfn, atom%orbitals%ener, &
                         atom%basis%nbas, atom%integrals%nne, atom%state%maxl_calc)
         CALL atom_denmat(atom%orbitals%pmat, atom%orbitals%wfn, atom%basis%nbas, atom%state%occupation, &
                          atom%state%maxl_occ, atom%state%maxn_occ)

      END DO !SCF Loop

      IF (iw > 0) THEN
         CALL atom_print_energies(atom, iw)
      END IF

      CALL atom_history_release(history)

      ! conserve fmat within atom_type
      CALL set_atom(atom, fmat=fmat)
      CALL release_opmat(jmat)
      CALL release_opmat(kmat)
      CALL release_opmat(xcmat)
      CALL release_opmat(hcore)

      CALL release_opgrid(density)
      CALL release_opgrid(cpot)

      ! ZMP deallocating ext_density ext_vxc
      IF (need_zmp) DEALLOCATE (ext_density)
      IF (need_vxc) DEALLOCATE (ext_vxc)

      IF (PRESENT(converged)) THEN
         converged = is_converged
      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_atom_restricted

! **************************************************************************************************
!> \brief Perform unrestricted (spin-polarised) electronic structure atomic calculations.
!> \param atom       information about the atomic kind
!> \param iw         output file unit
!> \param noguess    skip initial guess
!> \param converged  whether SCF iterations have been converged
!> \par History
!>    * identical to the subroutine calculate_atom_restricted()
! **************************************************************************************************
   SUBROUTINE calculate_atom_unrestricted(atom, iw, noguess, converged)
      TYPE(atom_type), POINTER                           :: atom
      INTEGER, INTENT(IN)                                :: iw
      LOGICAL, INTENT(IN), OPTIONAL                      :: noguess
      LOGICAL, INTENT(OUT), OPTIONAL                     :: converged

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_atom_unrestricted', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, iter, k, l, max_iter, method, &
                                                            reltyp
      LOGICAL                                            :: do_hfx, doguess, is_converged, lsdpot, &
                                                            need_x, need_xc
      REAL(KIND=dp)                                      :: deps, depsa, depsb, eps_scf, hf_frac, &
                                                            ne, nm
      TYPE(atom_history_type)                            :: historya, historyb
      TYPE(opgrid_type), POINTER                         :: cpot, density, rhoa, rhob
      TYPE(opmat_type), POINTER                          :: fmata, fmatb, hcore, hlsd, jmat, kmata, &
                                                            kmatb, xcmata, xcmatb
      TYPE(section_vals_type), POINTER                   :: xc_section

      CALL timeset(routineN, handle)

      IF (PRESENT(noguess)) THEN
         doguess = .NOT. noguess
      ELSE
         doguess = .TRUE.
      END IF

      CALL setup_hf_section(hf_frac, do_hfx, atom, xc_section, atom%exchange_integral_type)

      method = atom%method_type
      max_iter = atom%optimization%max_iter
      eps_scf = atom%optimization%eps_scf

      SELECT CASE (method)
      CASE DEFAULT
         CPABORT("")
      CASE (do_rks_atom)
         CPABORT("")
      CASE (do_uks_atom)
         need_x = do_hfx
         need_xc = .TRUE.
      CASE (do_rhf_atom)
         CPABORT("")
      CASE (do_uhf_atom)
         need_x = .TRUE.
         need_xc = .FALSE.
         hf_frac = 1._dp
      CASE (do_rohf_atom)
         need_x = .TRUE.
         need_xc = .FALSE.
         hf_frac = 1._dp
         CPABORT("")
      END SELECT

      ! set alpha and beta occupations
      IF (SUM(ABS(atom%state%occa) + ABS(atom%state%occb)) == 0.0_dp) THEN
         DO l = 0, 3
            nm = REAL((2*l + 1), KIND=dp)
            DO k = 1, 10
               ne = atom%state%occupation(l, k)
               IF (ne == 0._dp) THEN !empty shell
                  EXIT !assume there are no holes
               ELSEIF (ne == 2._dp*nm) THEN !closed shell
                  atom%state%occa(l, k) = nm
                  atom%state%occb(l, k) = nm
               ELSEIF (atom%state%multiplicity == -2) THEN !High spin case
                  atom%state%occa(l, k) = MIN(ne, nm)
                  atom%state%occb(l, k) = MAX(0._dp, ne - nm)
               ELSE
                  atom%state%occa(l, k) = 0.5_dp*(ne + atom%state%multiplicity - 1._dp)
                  atom%state%occb(l, k) = ne - atom%state%occa(l, k)
               END IF
            END DO
         END DO
      END IF
      ! check for relativistic method
      reltyp = atom%relativistic

      IF (iw > 0) CALL atom_print_state(atom%state, iw)

      NULLIFY (hcore, hlsd)
      CALL create_opmat(hcore, atom%basis%nbas)
      CALL create_opmat(hlsd, atom%basis%nbas)
      hlsd%op = 0._dp
      ! Pseudopotentials
      lsdpot = .FALSE.
      SELECT CASE (atom%potential%ppot_type)
      CASE DEFAULT
         CPABORT("")
      CASE (no_pseudo)
         SELECT CASE (reltyp)
         CASE DEFAULT
            CPABORT("")
         CASE (do_nonrel_atom)
            hcore%op = atom%integrals%kin - atom%zcore*atom%integrals%core
         CASE (do_zoramp_atom, do_sczoramp_atom)
            hcore%op = atom%integrals%kin + atom%integrals%tzora - atom%zcore*atom%integrals%core
         CASE (do_dkh0_atom, do_dkh1_atom, do_dkh2_atom, do_dkh3_atom)
            hcore%op = atom%integrals%hdkh
         END SELECT
      CASE (gth_pseudo)
         hcore%op = atom%integrals%kin + atom%integrals%core + atom%integrals%hnl
         IF (atom%potential%gth_pot%lsdpot) THEN
            lsdpot = .TRUE.
            hlsd%op = atom%integrals%clsd
         END IF
      CASE (upf_pseudo, sgp_pseudo, ecp_pseudo)
         hcore%op = atom%integrals%kin + atom%integrals%core + atom%integrals%hnl
      END SELECT
      ! add confinement potential (not included in relativistic transformations)
      IF (atom%potential%confinement) THEN
         hcore%op = hcore%op + atom%potential%acon*atom%integrals%conf
      END IF

      NULLIFY (fmata, fmatb, jmat, kmata, kmatb, xcmata, xcmatb)
      CALL create_opmat(fmata, atom%basis%nbas)
      CALL create_opmat(fmatb, atom%basis%nbas)
      CALL create_opmat(jmat, atom%basis%nbas)
      CALL create_opmat(kmata, atom%basis%nbas)
      CALL create_opmat(kmatb, atom%basis%nbas)
      CALL create_opmat(xcmata, atom%basis%nbas)
      CALL create_opmat(xcmatb, atom%basis%nbas)

      NULLIFY (density, rhoa, rhob, cpot)
      CALL create_opgrid(density, atom%basis%grid)
      CALL create_opgrid(rhoa, atom%basis%grid)
      CALL create_opgrid(rhob, atom%basis%grid)
      CALL create_opgrid(cpot, atom%basis%grid)

      IF (doguess) THEN
         ! initial guess
         CALL slater_density(rhoa%op, rhob%op, atom%z, atom%state, atom%basis%grid)
         density%op = rhoa%op + rhob%op
         CALL coulomb_potential_numeric(cpot%op, density%op, density%grid)
         CALL numpot_matrix(jmat%op, cpot%op, atom%basis, 0)
         ! alpha spin
         density%op = 2._dp*rhoa%op
         CALL wigner_slater_functional(density%op, cpot%op)
         CALL numpot_matrix(xcmata%op, cpot%op, atom%basis, 0)
         fmata%op = hcore%op + hlsd%op + jmat%op + xcmata%op
         CALL atom_solve(fmata%op, atom%integrals%utrans, atom%orbitals%wfna, atom%orbitals%enera, &
                         atom%basis%nbas, atom%integrals%nne, atom%state%maxl_calc)
         ! beta spin
         density%op = 2._dp*rhob%op
         CALL wigner_slater_functional(density%op, cpot%op)
         CALL numpot_matrix(xcmatb%op, cpot%op, atom%basis, 0)
         fmatb%op = hcore%op - hlsd%op + jmat%op + xcmatb%op
         CALL atom_solve(fmatb%op, atom%integrals%utrans, atom%orbitals%wfnb, atom%orbitals%enerb, &
                         atom%basis%nbas, atom%integrals%nne, atom%state%maxl_calc)
      END IF
      CALL atom_denmat(atom%orbitals%pmata, atom%orbitals%wfna, atom%basis%nbas, atom%state%occa, &
                       atom%state%maxl_occ, atom%state%maxn_occ)
      CALL atom_denmat(atom%orbitals%pmatb, atom%orbitals%wfnb, atom%basis%nbas, atom%state%occb, &
                       atom%state%maxl_occ, atom%state%maxn_occ)
      atom%orbitals%pmat = atom%orbitals%pmata + atom%orbitals%pmatb

      ! wavefunction history
      NULLIFY (historya%dmat, historya%hmat)
      CALL atom_history_init(historya, atom%optimization, fmata%op)
      NULLIFY (historyb%dmat, historyb%hmat)
      CALL atom_history_init(historyb, atom%optimization, fmatb%op)

      is_converged = .FALSE.
      iter = 0
      DO !SCF Loop

         ! Kinetic energy
         atom%energy%ekin = atom_trace(atom%integrals%kin, atom%orbitals%pmat)

         ! Band energy
         atom%energy%eband = 0._dp
         DO l = 0, 3
            DO i = 1, MIN(SIZE(atom%state%occupation, 2), SIZE(atom%orbitals%ener, 1))
               atom%energy%eband = atom%energy%eband + atom%orbitals%enera(i, l)*atom%state%occa(l, i)
               atom%energy%eband = atom%energy%eband + atom%orbitals%enerb(i, l)*atom%state%occb(l, i)
            END DO
         END DO

         ! Pseudopotential energy
         SELECT CASE (atom%potential%ppot_type)
         CASE DEFAULT
            CPABORT("")
         CASE (no_pseudo)
            atom%energy%eploc = 0._dp
            atom%energy%epnl = 0._dp
         CASE (gth_pseudo, upf_pseudo, sgp_pseudo, ecp_pseudo)
            atom%energy%eploc = atom_trace(atom%integrals%core, atom%orbitals%pmat)
            atom%energy%epnl = atom_trace(atom%integrals%hnl, atom%orbitals%pmat)
         END SELECT
         atom%energy%epseudo = atom%energy%eploc + atom%energy%epnl

         ! Core energy
         atom%energy%ecore = atom_trace(hcore%op, atom%orbitals%pmat)

         ! Confinement energy
         IF (atom%potential%confinement) THEN
            atom%energy%econfinement = atom_trace(atom%integrals%conf, atom%orbitals%pmat)
         ELSE
            atom%energy%econfinement = 0._dp
         END IF

         ! Hartree Term
         jmat%op = 0._dp
         SELECT CASE (atom%coulomb_integral_type)
         CASE DEFAULT
            CPABORT("")
         CASE (do_analytic)
            CALL ceri_contract(jmat%op, atom%integrals%ceri, atom%orbitals%pmat, atom%integrals%n)
         CASE (do_semi_analytic)
            CALL coulomb_potential_analytic(cpot%op, atom%orbitals%pmat, atom%basis, atom%basis%grid, &
                                            atom%state%maxl_occ)
            CALL numpot_matrix(jmat%op, cpot%op, atom%basis, 0)
         CASE (do_numeric)
            CALL atom_density(density%op, atom%orbitals%pmat, atom%basis, atom%state%maxl_occ, typ="RHO")
            CALL coulomb_potential_numeric(cpot%op, density%op, density%grid)
            CALL numpot_matrix(jmat%op, cpot%op, atom%basis, 0)
         END SELECT
         atom%energy%ecoulomb = 0.5_dp*atom_trace(jmat%op, atom%orbitals%pmat)

         ! Exchange Term
         IF (need_x) THEN
            kmata%op = 0._dp
            kmatb%op = 0._dp
            SELECT CASE (atom%exchange_integral_type)
            CASE DEFAULT
               CPABORT("")
            CASE (do_analytic)
               CALL eeri_contract(kmata%op, atom%integrals%eeri, atom%orbitals%pmata, atom%integrals%n)
               CALL eeri_contract(kmatb%op, atom%integrals%eeri, atom%orbitals%pmatb, atom%integrals%n)
            CASE (do_semi_analytic)
               CALL exchange_semi_analytic(kmata%op, atom%state, atom%state%occa, atom%orbitals%wfna, atom%basis, atom%hfx_pot)
               CALL exchange_semi_analytic(kmatb%op, atom%state, atom%state%occb, atom%orbitals%wfnb, atom%basis, atom%hfx_pot)
            CASE (do_numeric)
               CALL exchange_numeric(kmata%op, atom%state, atom%state%occa, atom%orbitals%wfna, atom%basis, atom%hfx_pot)
               CALL exchange_numeric(kmatb%op, atom%state, atom%state%occb, atom%orbitals%wfnb, atom%basis, atom%hfx_pot)
            END SELECT
            atom%energy%eexchange = hf_frac*(atom_trace(kmata%op, atom%orbitals%pmata) + &
                                             atom_trace(kmatb%op, atom%orbitals%pmatb))
            kmata%op = 2._dp*hf_frac*kmata%op
            kmatb%op = 2._dp*hf_frac*kmatb%op
         ELSE
            kmata%op = 0._dp
            kmatb%op = 0._dp
            atom%energy%eexchange = 0._dp
         END IF

         ! XC
         IF (need_xc) THEN
            xcmata%op = 0._dp
            xcmatb%op = 0._dp
            CALL calculate_atom_vxc_lsd(xcmata, xcmatb, atom, xc_section)
         ELSE
            xcmata%op = 0._dp
            xcmatb%op = 0._dp
            atom%energy%exc = 0._dp
         END IF

         IF (lsdpot) THEN
            atom%energy%elsd = atom_trace(hlsd%op, atom%orbitals%pmata) - &
                               atom_trace(hlsd%op, atom%orbitals%pmatb)
            atom%energy%epseudo = atom%energy%epseudo + atom%energy%elsd
            atom%energy%ecore = atom%energy%ecore + atom%energy%elsd
         ELSE
            atom%energy%elsd = 0._dp
         END IF

         ! Total energy
         atom%energy%etot = atom%energy%ecore + atom%energy%ecoulomb + atom%energy%eexchange + atom%energy%exc

         ! Potential energy
         atom%energy%epot = atom%energy%etot - atom%energy%ekin

         ! Total HF/KS matrix
         fmata%op = hcore%op + hlsd%op + jmat%op + kmata%op + xcmata%op
         fmatb%op = hcore%op - hlsd%op + jmat%op + kmatb%op + xcmatb%op

         ! calculate error matrix
         CALL err_matrix(xcmata%op, depsa, fmata%op, atom%orbitals%pmata, atom%integrals%utrans, &
                         atom%integrals%uptrans, atom%basis%nbas, atom%integrals%nne)
         CALL err_matrix(xcmatb%op, depsb, fmatb%op, atom%orbitals%pmatb, atom%integrals%utrans, &
                         atom%integrals%uptrans, atom%basis%nbas, atom%integrals%nne)
         deps = 2._dp*MAX(depsa, depsb)

         iter = iter + 1

         IF (iw > 0) THEN
            CALL atom_print_iteration(iter, deps, atom%energy%etot, iw)
         END IF

         IF (deps < eps_scf) THEN
            is_converged = .TRUE.
            EXIT
         END IF
         IF (iter >= max_iter) THEN
            IF (iw > 0) THEN
               WRITE (iw, "(A)") " No convergence within maximum number of iterations "
            END IF
            EXIT
         END IF

         ! update history container and extrapolate KS matrix
         CALL atom_history_update(historya, atom%orbitals%pmata, fmata%op, xcmata%op, atom%energy%etot, deps)
         CALL atom_history_update(historyb, atom%orbitals%pmatb, fmatb%op, xcmatb%op, atom%energy%etot, deps)
         CALL atom_opt_fmat(fmata%op, historya, depsa)
         CALL atom_opt_fmat(fmatb%op, historyb, depsb)

         ! Solve HF/KS equations
         CALL atom_solve(fmata%op, atom%integrals%utrans, atom%orbitals%wfna, atom%orbitals%enera, &
                         atom%basis%nbas, atom%integrals%nne, atom%state%maxl_calc)
         CALL atom_denmat(atom%orbitals%pmata, atom%orbitals%wfna, atom%basis%nbas, atom%state%occa, &
                          atom%state%maxl_occ, atom%state%maxn_occ)
         CALL atom_solve(fmatb%op, atom%integrals%utrans, atom%orbitals%wfnb, atom%orbitals%enerb, &
                         atom%basis%nbas, atom%integrals%nne, atom%state%maxl_calc)
         CALL atom_denmat(atom%orbitals%pmatb, atom%orbitals%wfnb, atom%basis%nbas, atom%state%occb, &
                          atom%state%maxl_occ, atom%state%maxn_occ)
         atom%orbitals%pmat = atom%orbitals%pmata + atom%orbitals%pmatb

      END DO !SCF Loop

      IF (iw > 0) THEN
         CALL atom_print_energies(atom, iw)
      END IF

      CALL atom_history_release(historya)
      CALL atom_history_release(historyb)

      CALL release_opgrid(density)
      CALL release_opgrid(rhoa)
      CALL release_opgrid(rhob)
      CALL release_opgrid(cpot)

      ! conserve fmata as fmat within atom_type.
      CALL set_atom(atom, fmat=fmata)
      CALL release_opmat(fmatb)
      CALL release_opmat(jmat)
      CALL release_opmat(kmata)
      CALL release_opmat(kmatb)
      CALL release_opmat(xcmata)
      CALL release_opmat(xcmatb)
      CALL release_opmat(hlsd)
      CALL release_opmat(hcore)

      IF (PRESENT(converged)) THEN
         converged = is_converged
      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_atom_unrestricted

END MODULE atom_electronic_structure
