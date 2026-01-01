!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of the core Hamiltonian integral matrix <a|H|b> over
!>      Cartesian Gaussian-type functions.
!>
!>      Nuclear potential energy:
!>
!>      a) Allelectron calculation:
!>
!>                          erfc(r)
!>         <a|V|b> = -Z*<a|---------|b>
!>                             r
!>
!>                          1 - erf(r)
!>                 = -Z*<a|------------|b>
!>                              r
!>
!>                           1           erf(r)
!>                 = -Z*(<a|---|b> - <a|--------|b>)
!>                           r             r
!>
!>                           1
!>                 = -Z*(<a|---|b> - N*<ab||c>)
!>                           r
!>
!>                      -Z
!>                 = <a|---|b> + Z*N*<ab||c>
!>                       r
!>                   \_______/       \_____/
!>                       |              |
!>                    nuclear        coulomb
!>
!>      b) Pseudopotential calculation (Goedecker, Teter and Hutter; GTH):
!>
!>         <a|V|b> = <a|(V(local) + V(non-local))|b>
!>
!>                 = <a|(V(local)|b> + <a|V(non-local))|b>
!>
!>         <a|V(local)|b> = <a|-Z(eff)*erf(SQRT(2)*alpha*r)/r +
!>                             (C1 + C2*(alpha*r)**2 + C3*(alpha*r)**4 +
!>                              C4*(alpha*r)**6)*exp(-(alpha*r)**2/2))|b>
!>
!>         <a|V(non-local)|b> = <a|p(l,i)>*h(i,j)*<p(l,j)|b>
!> \par Literature
!>      S. Goedecker, M. Teter and J. Hutter, Phys. Rev. B 54, 1703 (1996)
!>      C. Hartwigsen, S. Goedecker and J. Hutter, Phys. Rev. B 58, 3641 (1998)
!>      M. Krack and M. Parrinello, Phys. Chem. Chem. Phys. 2, 2105 (2000)
!>      S. Obara and A. Saika, J. Chem. Phys. 84, 3963 (1986)
!> \par History
!>      - Joost VandeVondele (April 2003) : added LSD forces
!>      - Non-redundant calculation of the non-local part of the GTH PP
!>        (22.05.2003,MK)
!>      - New parallelization scheme (27.06.2003,MK)
!>      - OpenMP version (07.12.2003,JGH)
!>      - Binary search loop for VPPNL operators (09.01.2004,JGH,MK)
!>      - Refactoring of pseudopotential and nuclear attraction integrals (25.02.2009,JGH)
!>      - General refactoring (01.10.2010,JGH)
!>      - Refactoring related to the new kinetic energy and overlap routines (07.2014,JGH)
!>      - k-point functionality (07.2015,JGH)
!>      - Refactor for PP functionality (08.2025,JGH)
!> \author Matthias Krack (14.09.2000,21.03.02)
! **************************************************************************************************
MODULE qs_core_matrices
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cell_types,                      ONLY: cell_type
   USE core_ae,                         ONLY: build_core_ae
   USE core_ppl,                        ONLY: build_core_ppl
   USE core_ppnl,                       ONLY: build_core_ppnl
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_iterator_blocks_left,&
                                              dbcsr_iterator_next_block,&
                                              dbcsr_iterator_start,&
                                              dbcsr_iterator_stop,&
                                              dbcsr_iterator_type,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE ec_env_types,                    ONLY: energy_correction_type
   USE input_constants,                 ONLY: do_ppl_analytic,&
                                              rel_none,&
                                              rel_trans_atom
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE lri_environment_types,           ONLY: lri_environment_type
   USE mathlib,                         ONLY: det_3x3
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: pascal
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_linres_types,                 ONLY: dcdr_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE virial_methods,                  ONLY: one_third_sum_diag
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_core_matrices'

   PUBLIC :: core_matrices, kinetic_energy_matrix

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix_h ...
!> \param matrix_p ...
!> \param calculate_forces ...
!> \param nder ...
!> \param ec_env ...
!> \param dcdr_env ...
!> \param ec_env_matrices ...
!> \param basis_type ...
!> \param debug_forces ...
!> \param debug_stress ...
!> \param atcore ...
! **************************************************************************************************
   SUBROUTINE core_matrices(qs_env, matrix_h, matrix_p, calculate_forces, nder, ec_env, dcdr_env, &
                            ec_env_matrices, basis_type, debug_forces, debug_stress, atcore)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER, INTENT(IN)                                :: nder
      TYPE(energy_correction_type), OPTIONAL, POINTER    :: ec_env
      TYPE(dcdr_env_type), OPTIONAL                      :: dcdr_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: ec_env_matrices
      CHARACTER(LEN=*), OPTIONAL                         :: basis_type
      LOGICAL, INTENT(IN), OPTIONAL                      :: debug_forces, debug_stress
      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: atcore

      CHARACTER(LEN=default_string_length)               :: my_basis_type
      INTEGER                                            :: iounit, natom, nimages
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: all_present, debfor, debstr, my_gt_nl, &
                                                            ppl_present, ppnl_present, use_lrigpw, &
                                                            use_virial
      REAL(KIND=dp)                                      :: eps_ppnl, fconv
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: stdeb
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: deltaR
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_hloc, matrix_ploc
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(lri_environment_type), POINTER                :: lri_env
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ae, sac_ppl, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial

      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      END IF

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, para_env=para_env)

      ! basis type
      IF (PRESENT(basis_type)) THEN
         my_basis_type = basis_type
      ELSE
         my_basis_type = "ORB"
      END IF

      IF (PRESENT(dcdr_env) .AND. PRESENT(ec_env)) THEN
         CPABORT("core_matrices: conflicting options")
      END IF

      ! check size of atcore array
      IF (PRESENT(atcore)) THEN
         CALL get_qs_env(qs_env=qs_env, natom=natom)
         CPASSERT(SIZE(atcore) >= natom)
      END IF

      ! check whether a gauge transformed version of the non-local potential part has to be used
      my_gt_nl = .FALSE.
      IF (qs_env%run_rtp) THEN
         CPASSERT(ASSOCIATED(dft_control%rtp_control))
         IF (dft_control%rtp_control%velocity_gauge) THEN
            my_gt_nl = dft_control%rtp_control%nl_gauge_transform
         END IF
      END IF

      ! prepare for k-points
      nimages = dft_control%nimages
      NULLIFY (cell_to_index)
      IF (nimages > 1) THEN
         CALL get_qs_env(qs_env=qs_env, ks_env=ks_env)
         CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
         dft_control%qs_control%do_ppl_method = do_ppl_analytic
      END IF

      ! Possible dc/dr terms
      IF (PRESENT(dcdr_env)) THEN
         deltaR => dcdr_env%delta_basis_function
         matrix_hloc(1:3, 1:1) => dcdr_env%matrix_hc(1:3)
      ELSE
         NULLIFY (deltaR)
         matrix_hloc => matrix_h
      END IF
      matrix_ploc => matrix_p

      ! force analytic ppl calculation for GAPW methods
      IF (dft_control%qs_control%gapw .OR. dft_control%qs_control%gapw_xc) THEN
         dft_control%qs_control%do_ppl_method = do_ppl_analytic
      END IF

      ! force
      NULLIFY (force)
      IF (calculate_forces) CALL get_qs_env(qs_env=qs_env, force=force)
      ! virial
      CALL get_qs_env(qs_env=qs_env, virial=virial)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      ! lri/gpw
      use_lrigpw = .FALSE.

      !debug
      debfor = .FALSE.
      IF (PRESENT(debug_forces)) debfor = debug_forces
      debfor = debfor .AND. calculate_forces
      debstr = .FALSE.
      IF (PRESENT(debug_stress)) debstr = debug_stress
      debstr = debstr .AND. use_virial
      IF (debstr) THEN
         CALL get_qs_env(qs_env=qs_env, cell=cell)
         fconv = 1.0E-9_dp*pascal/cell%deth
      END IF

      NULLIFY (qs_kind_set, atomic_kind_set, particle_set)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set)

      NULLIFY (sab_orb, sac_ae, sac_ppl, sap_ppnl)
      IF (PRESENT(ec_env)) THEN
         sab_orb => ec_env%sab_orb
         sac_ae => ec_env%sac_ae
         sac_ppl => ec_env%sac_ppl
         sap_ppnl => ec_env%sap_ppnl
      ELSE
         CALL get_qs_env(qs_env=qs_env, &
                         sab_orb=sab_orb, &
                         sac_ae=sac_ae, &
                         sac_ppl=sac_ppl, &
                         sap_ppnl=sap_ppnl)
      END IF

      IF (PRESENT(ec_env) .AND. PRESENT(ec_env_matrices)) THEN
         IF (ec_env_matrices) THEN
            matrix_hloc => ec_env%matrix_h
            matrix_ploc => ec_env%matrix_p
         END IF
      END IF

      ! *** compute the nuclear attraction contribution to the core hamiltonian ***
      all_present = ASSOCIATED(sac_ae)
      IF (all_present) THEN
         IF (PRESENT(dcdr_env)) THEN
            CPABORT("ECP/AE functionality for dcdr missing")
         END IF
         IF (debfor) fodeb(1:3) = force(1)%all_potential(1:3, 1)
         IF (debstr) stdeb = virial%pv_ppl
         CALL build_core_ae(matrix_hloc, matrix_ploc, force, virial, calculate_forces, use_virial, nder, &
                            qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ae, &
                            nimages, cell_to_index, my_basis_type, atcore=atcore)
         IF (debfor) THEN
            fodeb(1:3) = force(1)%all_potential(1:3, 1) - fodeb(1:3)
            CALL para_env%sum(fodeb)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*dHae    ", fodeb
         END IF
         IF (debstr) THEN
            stdeb = fconv*(virial%pv_ppl - stdeb)
            CALL para_env%sum(stdeb)
            IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
               'STRESS| P*dHae    ', one_third_sum_diag(stdeb), det_3x3(stdeb)
         END IF
      END IF

      ! *** compute the ppl contribution to the core hamiltonian ***
      ppl_present = ASSOCIATED(sac_ppl)
      IF (ppl_present) THEN
         IF (dft_control%qs_control%do_ppl_method == do_ppl_analytic) THEN
            CALL get_qs_env(qs_env, lri_env=lri_env)
            IF (ASSOCIATED(lri_env)) THEN
               use_lrigpw = (dft_control%qs_control%lrigpw .AND. lri_env%ppl_ri)
            END IF
            IF (use_lrigpw) THEN
               IF (lri_env%exact_1c_terms) THEN
                  CPABORT("not implemented")
               END IF
            ELSE
               IF (debfor) fodeb(1:3) = force(1)%gth_ppl(1:3, 1)
               IF (debstr) stdeb = virial%pv_ppl
               CALL build_core_ppl(matrix_hloc, matrix_ploc, force, &
                                   virial, calculate_forces, use_virial, nder, &
                                   qs_kind_set, atomic_kind_set, particle_set, &
                                   sab_orb, sac_ppl, nimages, cell_to_index, my_basis_type, &
                                   deltaR=deltaR, atcore=atcore)
               IF (debfor) THEN
                  fodeb(1:3) = force(1)%gth_ppl(1:3, 1) - fodeb(1:3)
                  CALL para_env%sum(fodeb)
                  IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*dHppl   ", fodeb
               END IF
               IF (debstr) THEN
                  stdeb = fconv*(virial%pv_ppl - stdeb)
                  CALL para_env%sum(stdeb)
                  IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                     'STRESS| P*dHppl   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
               END IF
            END IF
         END IF
      END IF

      ! *** compute the ppnl contribution to the core hamiltonian ***
      eps_ppnl = dft_control%qs_control%eps_ppnl
      ppnl_present = ASSOCIATED(sap_ppnl)
      IF (ppnl_present) THEN
         IF (PRESENT(dcdr_env)) THEN
            matrix_hloc(1:3, 1:1) => dcdr_env%matrix_ppnl_1(1:3)
         END IF
         IF (.NOT. my_gt_nl) THEN
            IF (debfor) fodeb(1:3) = force(1)%gth_ppnl(1:3, 1)
            IF (debstr) stdeb = virial%pv_ppnl
            CALL build_core_ppnl(matrix_hloc, matrix_ploc, force, &
                                 virial, calculate_forces, use_virial, nder, &
                                 qs_kind_set, atomic_kind_set, particle_set, &
                                 sab_orb, sap_ppnl, eps_ppnl, nimages, cell_to_index, my_basis_type, &
                                 deltaR=deltaR, atcore=atcore)
            IF (debfor) THEN
               fodeb(1:3) = force(1)%gth_ppnl(1:3, 1) - fodeb(1:3)
               CALL para_env%sum(fodeb)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*dHppnl  ", fodeb
            END IF
            IF (debstr) THEN
               stdeb = fconv*(virial%pv_ppnl - stdeb)
               CALL para_env%sum(stdeb)
               IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
                  'STRESS| P*dHppnl   ', one_third_sum_diag(stdeb), det_3x3(stdeb)
            END IF
         END IF
      END IF

   END SUBROUTINE core_matrices

! **************************************************************************************************
!> \brief Calculate kinetic energy matrix and possible relativistic correction
!> \param qs_env ...
!> \param matrixkp_t ...
!> \param matrix_t ...
!> \param matrix_p ...
!> \param matrix_name ...
!> \param calculate_forces ...
!> \param nderivative ...
!> \param sab_orb ...
!> \param eps_filter ...
!> \param basis_type ...
!> \param debug_forces ...
!> \param debug_stress ...
!> \author Creation (21.08.2025,JGH)
! **************************************************************************************************
   SUBROUTINE kinetic_energy_matrix(qs_env, matrixkp_t, matrix_t, matrix_p, &
                                    matrix_name, calculate_forces, nderivative, &
                                    sab_orb, eps_filter, basis_type, debug_forces, debug_stress)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: matrixkp_t
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: matrix_t
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: matrix_p
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: matrix_name
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER, INTENT(IN), OPTIONAL                      :: nderivative
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_orb
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps_filter
      CHARACTER(LEN=*), OPTIONAL                         :: basis_type
      LOGICAL, INTENT(IN), OPTIONAL                      :: debug_forces, debug_stress

      CHARACTER(LEN=default_string_length)               :: my_basis_type
      INTEGER                                            :: ic, img, iounit, nimages
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: debfor, debstr, use_virial
      REAL(KIND=dp)                                      :: fconv
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      REAL(KIND=dp), DIMENSION(3, 3)                     :: stdeb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(virial_type), POINTER                         :: virial

      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      END IF

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, para_env=para_env)
      ! is this a orbital-free method calculation
      IF (dft_control%qs_control%ofgpw) RETURN

      ! Matrix images (kp)
      nimages = dft_control%nimages

      ! basis type
      IF (PRESENT(basis_type)) THEN
         my_basis_type = basis_type
      ELSE
         my_basis_type = "ORB"
      END IF

      debfor = .FALSE.
      IF (PRESENT(debug_forces)) debfor = debug_forces
      debfor = debfor .AND. calculate_forces

      CALL get_qs_env(qs_env=qs_env, virial=virial)
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      debstr = .FALSE.
      IF (PRESENT(debug_stress)) debstr = debug_stress
      debstr = debstr .AND. use_virial
      IF (debstr) THEN
         CALL get_qs_env(qs_env=qs_env, cell=cell)
         fconv = 1.0E-9_dp*pascal/cell%deth
      END IF

      NULLIFY (ks_env)
      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env)
      NULLIFY (sab_nl)
      IF (PRESENT(sab_orb)) THEN
         sab_nl => sab_orb
      ELSE
         CALL get_qs_env(qs_env=qs_env, sab_orb=sab_nl)
      END IF

      IF (calculate_forces) THEN
         IF (SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimages
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            END DO
         END IF
      END IF

      IF (debfor) THEN
         CALL get_qs_env(qs_env=qs_env, force=force)
         fodeb(1:3) = force(1)%kinetic(1:3, 1)
      END IF
      IF (debstr) THEN
         stdeb = virial%pv_ekinetic
      END IF

      ! T matrix
      CALL build_kinetic_matrix(ks_env, matrixkp_t=matrixkp_t, matrix_t=matrix_t, &
                                matrix_name=matrix_name, &
                                basis_type=my_basis_type, &
                                sab_nl=sab_nl, &
                                calculate_forces=calculate_forces, &
                                matrixkp_p=matrix_p, &
                                nderivative=nderivative, &
                                eps_filter=eps_filter)

      IF (debfor) THEN
         fodeb(1:3) = force(1)%kinetic(1:3, 1) - fodeb(1:3)
         CALL para_env%sum(fodeb)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: P*dT    ", fodeb
      END IF
      IF (debstr) THEN
         stdeb = fconv*(virial%pv_ekinetic - stdeb)
         CALL para_env%sum(stdeb)
         IF (iounit > 0) WRITE (UNIT=iounit, FMT="(T2,A,T41,2(1X,ES19.11))") &
            'STRESS| P*dT', one_third_sum_diag(stdeb), det_3x3(stdeb)
      END IF

      IF (calculate_forces) THEN
         IF (SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimages
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=-1.0_dp)
            END DO
         END IF
      END IF

      ! relativistic atomic correction to kinetic energy
      IF (qs_env%rel_control%rel_method /= rel_none) THEN
         IF (qs_env%rel_control%rel_transformation == rel_trans_atom) THEN
            CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set)
            IF (nimages == 1) THEN
               ic = 1
            ELSE
               CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
               CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
               ic = cell_to_index(0, 0, 0)
            END IF
            IF (my_basis_type /= "ORB") THEN
               CPABORT("Basis mismatch for relativistic correction")
            END IF
            IF (PRESENT(matrixkp_t)) THEN
               CALL build_atomic_relmat(matrixkp_t(1, ic)%matrix, atomic_kind_set, qs_kind_set)
            ELSEIF (PRESENT(matrix_t)) THEN
               CALL build_atomic_relmat(matrix_t(1)%matrix, atomic_kind_set, qs_kind_set)
            END IF
         ELSE
            CPABORT("Relativistic corrections of this type are currently not implemented")
         END IF
      END IF

   END SUBROUTINE kinetic_energy_matrix

! **************************************************************************************************
!> \brief Adds atomic blocks of relativistic correction for the kinetic energy
!> \param matrix_t ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
! **************************************************************************************************
   SUBROUTINE build_atomic_relmat(matrix_t, atomic_kind_set, qs_kind_set)
      TYPE(dbcsr_type), POINTER                          :: matrix_t
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      INTEGER                                            :: iatom, ikind, jatom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: reltmat, tblock
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      CALL dbcsr_iterator_start(iter, matrix_t)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, iatom, jatom, tblock)
         IF (iatom == jatom) THEN
            ikind = kind_of(iatom)
            CALL get_qs_kind(qs_kind_set(ikind), reltmat=reltmat)
            IF (ASSOCIATED(reltmat)) tblock = tblock + reltmat
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

   END SUBROUTINE build_atomic_relmat

END MODULE qs_core_matrices
