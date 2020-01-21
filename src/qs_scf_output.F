!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_scf_output
   USE admm_types,                      ONLY: admm_type
   USE admm_utils,                      ONLY: admm_correct_for_eigenvalues,&
                                              admm_uncorrect_for_eigenvalues
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_init_random,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_type
   USE input_constants,                 ONLY: &
        becke_cutoff_element, becke_cutoff_global, cdft_alpha_constraint, cdft_beta_constraint, &
        cdft_charge_constraint, cdft_magnetization_constraint, ot_precond_full_all, &
        outer_scf_becke_constraint, outer_scf_hirshfeld_constraint, outer_scf_optimizer_bisect, &
        outer_scf_optimizer_broyden, outer_scf_optimizer_diis, outer_scf_optimizer_newton, &
        outer_scf_optimizer_newton_ls, outer_scf_optimizer_sd, outer_scf_optimizer_secant, &
        radius_covalent, radius_default, radius_single, radius_user, radius_vdw, &
        shape_function_density, shape_function_gaussian
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kahan_sum,                       ONLY: accurate_sum
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE machine,                         ONLY: m_flush
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt,&
                                              kcalmol
   USE preconditioner_types,            ONLY: preconditioner_type
   USE ps_implicit_types,               ONLY: MIXED_BC,&
                                              MIXED_PERIODIC_BC,&
                                              NEUMANN_BC,&
                                              PERIODIC_BC
   USE pw_env_types,                    ONLY: pw_env_type
   USE pw_poisson_types,                ONLY: pw_poisson_implicit
   USE qmmm_image_charge,               ONLY: print_image_coefficients
   USE qs_cdft_opt_types,               ONLY: cdft_opt_type_write
   USE qs_cdft_types,                   ONLY: cdft_control_type
   USE qs_charges_types,                ONLY: qs_charges_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_io,                        ONLY: write_mo_set_to_output_unit
   USE qs_mo_methods,                   ONLY: calculate_magnitude,&
                                              calculate_orthonormality,&
                                              calculate_subspace_eigenvalues
   USE qs_mo_occupation,                ONLY: set_mo_occupation
   USE qs_mo_types,                     ONLY: allocate_mo_set,&
                                              deallocate_mo_set,&
                                              get_mo_set,&
                                              init_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_ot_eigensolver,               ONLY: ot_eigensolver
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: ot_method_nr,&
                                              qs_scf_env_type,&
                                              special_diag_method_nr
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_output'

   PUBLIC :: qs_scf_loop_info, &
             qs_scf_print_summary, &
             qs_scf_loop_print, &
             qs_scf_outer_loop_info, &
             qs_scf_initial_info, &
             qs_scf_write_mos, &
             qs_scf_cdft_info, &
             qs_scf_cdft_initial_info, &
             qs_scf_cdft_constraint_info

CONTAINS

! **************************************************************************************************
!> \brief writes a summary of information after scf
!> \param output_unit ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE qs_scf_print_summary(output_unit, qs_env)
      INTEGER, INTENT(IN)                                :: output_unit
      TYPE(qs_environment_type), POINTER                 :: qs_env

      INTEGER                                            :: nelectron_total
      LOGICAL                                            :: gapw, gapw_xc, qmmm
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_charges_type), POINTER                     :: qs_charges
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      NULLIFY (rho, energy, dft_control, scf_env, qs_charges)
      CALL get_qs_env(qs_env=qs_env, rho=rho, energy=energy, dft_control=dft_control, &
                      scf_env=scf_env, qs_charges=qs_charges)

      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc
      qmmm = qs_env%qmmm
      nelectron_total = scf_env%nelectron

      CALL qs_scf_print_scf_summary(output_unit, rho, qs_charges, energy, nelectron_total, &
                                    dft_control, qmmm, qs_env, gapw, gapw_xc)

   END SUBROUTINE qs_scf_print_summary

! **************************************************************************************************
!> \brief writes basic information at the beginning of an scf run
!> \param output_unit ...
!> \param mos ...
!> \param dft_control ...
! **************************************************************************************************
   SUBROUTINE qs_scf_initial_info(output_unit, mos, dft_control)
      INTEGER                                            :: output_unit
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(dft_control_type), POINTER                    :: dft_control

      INTEGER                                            :: homo, ispin, nao, nelectron_spin, nmo

      IF (output_unit > 0) THEN
         DO ispin = 1, dft_control%nspins
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
                            homo=homo, &
                            nelectron=nelectron_spin, &
                            nao=nao, &
                            nmo=nmo)
            IF (dft_control%nspins > 1) THEN
               WRITE (UNIT=output_unit, FMT="(/,T2,A,I2)") "Spin", ispin
            END IF
            WRITE (UNIT=output_unit, FMT="(/,(T2,A,T71,I10))") &
               "Number of electrons:", nelectron_spin, &
               "Number of occupied orbitals:", homo, &
               "Number of molecular orbitals:", nmo
         END DO
         WRITE (UNIT=output_unit, FMT="(/,T2,A,T71,I10)") &
            "Number of orbital functions:", nao
      END IF

   END SUBROUTINE qs_scf_initial_info

! **************************************************************************************************
!> \brief Write the MO eigenvector, eigenvalues, and occupation numbers to the output unit
!> \param qs_env ...
!> \param scf_env ...
!> \param final_mos ...
!> \par History
!>      - Revise MO printout to enable eigenvalues with OT (05.05.2021, MK)
! **************************************************************************************************
   SUBROUTINE qs_scf_write_mos(qs_env, scf_env, final_mos)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      LOGICAL, INTENT(IN)                                :: final_mos

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'qs_scf_write_mos'

      CHARACTER(LEN=2)                                   :: solver_method
      CHARACTER(LEN=3*default_string_length)             :: message
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: tmpstringlist
      INTEGER                                            :: handle, ikp, ispin, iw, kpoint, nao, &
                                                            nelectron, nkp, nmo, nmo_tmp, nspin
      INTEGER, DIMENSION(2)                              :: nmos_occ
      INTEGER, DIMENSION(:), POINTER                     :: mo_index_range
      LOGICAL                                            :: do_kpoints, print_eigvals, &
                                                            print_eigvecs, print_mo_info, &
                                                            print_occup, print_occup_stats
      REAL(KIND=dp)                                      :: flexible_electron_count, maxocc, n_el_f, &
                                                            occup_stats_occ_threshold
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues, mo_eigenvalues_tmp
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: mo_coeff, mo_coeff_tmp
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks, s
      TYPE(dbcsr_type), POINTER                          :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(mo_set_type), POINTER                         :: mo_set_tmp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(preconditioner_type), POINTER                 :: local_preconditioner
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))

      ! Retrieve the required information for the requested print output
      CALL get_qs_env(qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      blacs_env=blacs_env, &
                      dft_control=dft_control, &
                      do_kpoints=do_kpoints, &
                      input=input, &
                      qs_kind_set=qs_kind_set, &
                      para_env=para_env, &
                      particle_set=particle_set, &
                      scf_control=scf_control)

      ! Quick return, if no printout of MO information is requested
      dft_section => section_vals_get_subs_vals(input, "DFT")
      CALL section_vals_val_get(dft_section, "PRINT%MO%EIGENVALUES", l_val=print_eigvals)
      CALL section_vals_val_get(dft_section, "PRINT%MO%EIGENVECTORS", l_val=print_eigvecs)
      CALL section_vals_val_get(dft_section, "PRINT%MO%OCCUPATION_NUMBERS", l_val=print_occup)
      CALL section_vals_val_get(dft_section, "PRINT%MO%OCCUPATION_NUMBERS_STATS", c_vals=tmpstringlist)

      print_occup_stats = .FALSE.
      occup_stats_occ_threshold = 1e-6_dp
      IF (SIZE(tmpstringlist) > 0) THEN  ! the lone_keyword_c_vals doesn't work as advertised, handle it manually
         print_occup_stats = .TRUE.
         IF (LEN_TRIM(tmpstringlist(1)) > 0) &
            READ (tmpstringlist(1), *) print_occup_stats
      END IF
      IF (SIZE(tmpstringlist) > 1) &
         READ (tmpstringlist(2), *) occup_stats_occ_threshold

      logger => cp_get_default_logger()
      print_mo_info = (cp_print_key_should_output(logger%iter_info, dft_section, "PRINT%MO") /= 0)

      IF ((.NOT. print_mo_info) .OR. (.NOT. (print_eigvals .OR. print_eigvecs .OR. print_occup .OR. print_occup_stats))) THEN
         CALL timestop(handle)
         RETURN
      END IF

      NULLIFY (fm_struct_tmp)
      NULLIFY (mo_coeff)
      NULLIFY (mo_coeff_tmp)
      NULLIFY (mo_set_tmp)

      nspin = dft_control%nspins
      nmos_occ = 0

      ! Check, if we have k points
      IF (do_kpoints) THEN
         CALL get_qs_env(qs_env, kpoints=kpoints)
         nkp = SIZE(kpoints%kp_env)
      ELSE
         CALL get_qs_env(qs_env, matrix_ks=ks, matrix_s=s, mos=mos)
         CPASSERT(ASSOCIATED(ks))
         CPASSERT(ASSOCIATED(s))
         nkp = 1
      END IF

      DO ikp = 1, nkp

         IF (do_kpoints) THEN
            mos => kpoints%kp_env(ikp)%kpoint_env%mos(1, :)
            kpoint = ikp
         ELSE
            CALL get_qs_env(qs_env, matrix_ks=ks, mos=mos)
            kpoint = 0 ! Gamma point only
         END IF
         CPASSERT(ASSOCIATED(mos))

         ! Prepare MO information for printout
         DO ispin = 1, nspin

            ! Calculate MO eigenvalues and eigenvector when OT is used
            IF (scf_env%method == ot_method_nr) THEN
               solver_method = "OT"

               IF (do_kpoints) THEN
                  CPABORT("The OT method is not implemented for k points")
               END IF

               matrix_ks => ks(ispin)%matrix
               matrix_s => s(1)%matrix

               ! With ADMM, we have to modify the Kohn-Sham matrix
               IF (dft_control%do_admm) THEN
                  CALL get_qs_env(qs_env, admm_env=admm_env)
                  CALL admm_correct_for_eigenvalues(ispin, admm_env, matrix_ks)
               END IF

               CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
                               mo_coeff=mo_coeff, &
                               eigenvalues=mo_eigenvalues, &
                               maxocc=maxocc, &
                               nelectron=nelectron, &
                               n_el_f=n_el_f, &
                               nao=nao, &
                               nmo=nmo, &
                               flexible_electron_count=flexible_electron_count)

               mo_index_range => section_get_ivals(dft_section, "PRINT%MO%MO_INDEX_RANGE")
               CPASSERT(ASSOCIATED(mo_index_range))

               ! Retrieve the index of the last MO for which a printout is requested
               nmo_tmp = mo_index_range(2)
               IF (nmo_tmp <= nmo) THEN
                  nmo_tmp = nmo
               ELSE
                  IF (.NOT. final_mos) THEN
                     message = "The MO information for unoccupied MOs is only calculated after "// &
                               "SCF convergence is achieved when the orbital transformation (OT) "// &
                               "method is used"
                     CPWARN(TRIM(message))
                     nmo_tmp = nmo
                  END IF
               END IF

               ! Create temporary MO set for printout
               CALL cp_fm_struct_create(fm_struct_tmp, &
                                        context=blacs_env, &
                                        para_env=para_env, &
                                        nrow_global=nao, &
                                        ncol_global=nmo_tmp)
               CALL allocate_mo_set(mo_set=mo_set_tmp, &
                                    nao=nao, &
                                    nmo=nmo_tmp, &
                                    nelectron=nelectron, &
                                    n_el_f=n_el_f, &
                                    maxocc=maxocc, &
                                    flexible_electron_count=flexible_electron_count)
               CALL init_mo_set(mo_set=mo_set_tmp, &
                                fm_struct=fm_struct_tmp, &
                                name="Temporary MO set for printout")
               CALL cp_fm_struct_release(fm_struct_tmp)

               mo_coeff_tmp => mo_set_tmp%mo_coeff
               mo_eigenvalues_tmp => mo_set_tmp%eigenvalues

               ! Calculate MOs information for the request MO index range
               IF (final_mos) THEN
                  ! Prepare printout of additional unoccupied MOs
                  IF (nmo_tmp > nmo) THEN
                     CALL cp_fm_init_random(mo_set_tmp%mo_coeff, nmo_tmp)
                     ! The FULL_ALL preconditioner makes not much sense for the unoccupied orbitals
                     NULLIFY (local_preconditioner)
                     IF (ASSOCIATED(scf_env%ot_preconditioner)) THEN
                        local_preconditioner => scf_env%ot_preconditioner(1)%preconditioner
                        IF (local_preconditioner%in_use == ot_precond_full_all) THEN
                           NULLIFY (local_preconditioner)
                        END IF
                     END IF
                     CALL ot_eigensolver(matrix_h=matrix_ks, &
                                         matrix_s=matrix_s, &
                                         matrix_c_fm=mo_coeff_tmp, &
                                         matrix_orthogonal_space_fm=mo_coeff, &
                                         eps_gradient=scf_control%eps_lumos, &
                                         preconditioner=local_preconditioner, &
                                         iter_max=scf_control%max_iter_lumos, &
                                         size_ortho_space=nmo)
                  ELSE
                     CALL cp_fm_to_fm(mo_coeff, mo_coeff_tmp, nmo)
                  END IF
               ELSE
                  CALL cp_fm_to_fm(mo_coeff, mo_coeff_tmp, nmo)
               END IF ! final MOs

               CALL calculate_subspace_eigenvalues(orbitals=mo_coeff_tmp, &
                                                   ks_matrix=matrix_ks, &
                                                   evals_arg=mo_eigenvalues_tmp, &
                                                   do_rotation=.TRUE.)
               CALL set_mo_occupation(mo_set=mo_set_tmp)

               ! With ADMM, we have to undo the modification of the Kohn-Sham matrix
               IF (dft_control%do_admm) THEN
                  CALL admm_uncorrect_for_eigenvalues(ispin, admm_env, matrix_ks)
               END IF

            ELSE

               solver_method = "TD"
               mo_set_tmp => mos(ispin)%mo_set

            END IF ! OT is used

            ! Print MO information
            IF (nspin > 1) THEN
               SELECT CASE (ispin)
               CASE (1)
                  CALL write_mo_set_to_output_unit(mo_set_tmp, atomic_kind_set, qs_kind_set, particle_set, &
                                                   dft_section, 4, kpoint, final_mos=final_mos, spin="ALPHA", &
                                                   solver_method=solver_method)
               CASE (2)
                  CALL write_mo_set_to_output_unit(mo_set_tmp, atomic_kind_set, qs_kind_set, particle_set, &
                                                   dft_section, 4, kpoint, final_mos=final_mos, spin="BETA", &
                                                   solver_method=solver_method)
               CASE DEFAULT
                  CPABORT("Invalid spin")
               END SELECT
            ELSE
               CALL write_mo_set_to_output_unit(mo_set_tmp, atomic_kind_set, qs_kind_set, particle_set, &
                                                dft_section, 4, kpoint, final_mos=final_mos, &
                                                solver_method=solver_method)
            END IF

            nmos_occ(ispin) = MAX(nmos_occ(ispin), COUNT(mo_set_tmp%occupation_numbers > occup_stats_occ_threshold))

            ! Deallocate temporary objects needed for OT
            IF (scf_env%method == ot_method_nr) THEN
               CALL deallocate_mo_set(mo_set_tmp)
               NULLIFY (matrix_ks)
               NULLIFY (matrix_s)
            END IF
            NULLIFY (mo_set_tmp)

         END DO ! ispin

      END DO ! k point loop
      IF (print_mo_info .AND. print_occup_stats) THEN
         iw = cp_print_key_unit_nr(logger, dft_section, "PRINT%MO", &
                                   ignore_should_output=print_mo_info, &
                                   extension=".MOLog")
         IF (iw > 0) THEN
            IF (SIZE(mos) > 1) THEN
               WRITE (UNIT=iw, FMT="(A,I4)") " MO| Total occupied (ALPHA):", nmos_occ(1)
               WRITE (UNIT=iw, FMT="(A,I4)") " MO| Total occupied (BETA): ", nmos_occ(2)
            ELSE
               WRITE (UNIT=iw, FMT="(A,I4)") " MO| Total occupied: ", nmos_occ(1)
            END IF
            WRITE (UNIT=iw, FMT="(A)") ""
         END IF
         CALL cp_print_key_finished_output(iw, logger, dft_section, "PRINT%MO", &
                                           ignore_should_output=print_mo_info)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_write_mos

! **************************************************************************************************
!> \brief writes basic information obtained in a scf outer loop step
!> \param output_unit ...
!> \param scf_control ...
!> \param scf_env ...
!> \param energy ...
!> \param total_steps ...
!> \param should_stop ...
!> \param outer_loop_converged ...
! **************************************************************************************************
   SUBROUTINE qs_scf_outer_loop_info(output_unit, scf_control, scf_env, &
                                     energy, total_steps, should_stop, outer_loop_converged)
      INTEGER                                            :: output_unit
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_energy_type), POINTER                      :: energy
      INTEGER                                            :: total_steps
      LOGICAL, INTENT(IN)                                :: should_stop, outer_loop_converged

      REAL(KIND=dp)                                      :: outer_loop_eps

      outer_loop_eps = SQRT(MAXVAL(scf_env%outer_scf%gradient(:, scf_env%outer_scf%iter_count)**2))
      IF (output_unit > 0) WRITE (output_unit, '(/,T3,A,I4,A,E10.2,A,F22.10)') &
         "outer SCF iter = ", scf_env%outer_scf%iter_count, &
         " RMS gradient = ", outer_loop_eps, " energy =", energy%total

      IF (outer_loop_converged) THEN
         IF (output_unit > 0) WRITE (output_unit, '(T3,A,I4,A,I4,A,/)') &
            "outer SCF loop converged in", scf_env%outer_scf%iter_count, &
            " iterations or ", total_steps, " steps"
      ELSE IF (scf_env%outer_scf%iter_count > scf_control%outer_scf%max_scf &
               .OR. should_stop) THEN
         IF (output_unit > 0) WRITE (output_unit, '(T3,A,I4,A,I4,A,/)') &
            "outer SCF loop FAILED to converge after ", &
            scf_env%outer_scf%iter_count, " iterations or ", total_steps, " steps"
      END IF

   END SUBROUTINE qs_scf_outer_loop_info

! **************************************************************************************************
!> \brief writes basic information obtained in a scf step
!> \param scf_env ...
!> \param output_unit ...
!> \param just_energy ...
!> \param t1 ...
!> \param t2 ...
!> \param energy ...
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_info(scf_env, output_unit, just_energy, t1, t2, energy)

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      INTEGER                                            :: output_unit
      LOGICAL                                            :: just_energy
      REAL(KIND=dp)                                      :: t1, t2
      TYPE(qs_energy_type), POINTER                      :: energy

      IF ((output_unit > 0) .AND. scf_env%print_iter_line) THEN
         IF (just_energy) THEN
            WRITE (UNIT=output_unit, &
                   FMT="(T2,I5,1X,A,T20,E8.2,1X,F6.1,16X,F20.10)") &
               scf_env%iter_count, TRIM(scf_env%iter_method), scf_env%iter_param, &
               t2 - t1, energy%total
         ELSE
            IF ((ABS(scf_env%iter_delta) < 1.0E-8_dp) .OR. &
                (ABS(scf_env%iter_delta) >= 1.0E5_dp)) THEN
               WRITE (UNIT=output_unit, &
                      FMT="(T2,I5,1X,A,T20,E8.2,1X,F6.1,1X,ES14.4,1X,F20.10,1X,ES9.2)") &
                  scf_env%iter_count, TRIM(scf_env%iter_method), scf_env%iter_param, &
                  t2 - t1, scf_env%iter_delta, energy%total, energy%total - energy%tot_old
            ELSE
               WRITE (UNIT=output_unit, &
                      FMT="(T2,I5,1X,A,T20,E8.2,1X,F6.1,1X,F14.8,1X,F20.10,1X,ES9.2)") &
                  scf_env%iter_count, TRIM(scf_env%iter_method), scf_env%iter_param, &
                  t2 - t1, scf_env%iter_delta, energy%total, energy%total - energy%tot_old
            END IF
         END IF
      END IF

   END SUBROUTINE qs_scf_loop_info

! **************************************************************************************************
!> \brief writes rather detailed summary of densities and energies
!>      after the SCF
!> \param output_unit ...
!> \param rho ...
!> \param qs_charges ...
!> \param energy ...
!> \param nelectron_total ...
!> \param dft_control ...
!> \param qmmm ...
!> \param qs_env ...
!> \param gapw ...
!> \param gapw_xc ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
!>      10.2019 print dipole moment [SGh]
! **************************************************************************************************
   SUBROUTINE qs_scf_print_scf_summary(output_unit, rho, qs_charges, energy, nelectron_total, &
                                       dft_control, qmmm, qs_env, gapw, gapw_xc)
      INTEGER, INTENT(IN)                                :: output_unit
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_charges_type), POINTER                     :: qs_charges
      TYPE(qs_energy_type), POINTER                      :: energy
      INTEGER, INTENT(IN)                                :: nelectron_total
      TYPE(dft_control_type), POINTER                    :: dft_control
      LOGICAL, INTENT(IN)                                :: qmmm
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: gapw, gapw_xc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_print_scf_summary'

      INTEGER                                            :: bc, handle, ispin, psolver
      REAL(kind=dp)                                      :: exc1_energy, exc_energy, &
                                                            implicit_ps_ehartree, tot1_h, tot1_s
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      TYPE(pw_env_type), POINTER                         :: pw_env

      NULLIFY (tot_rho_r, pw_env)
      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      psolver = pw_env%poisson_env%parameters%solver

      IF (output_unit > 0) THEN
         CALL qs_rho_get(rho, tot_rho_r=tot_rho_r)
         IF (.NOT. (dft_control%qs_control%semi_empirical .OR. &
                    dft_control%qs_control%xtb .OR. &
                    dft_control%qs_control%dftb)) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T41,2F20.10))") &
               "Electronic density on regular grids: ", &
               accurate_sum(tot_rho_r), &
               accurate_sum(tot_rho_r) + nelectron_total, &
               "Core density on regular grids:", &
               qs_charges%total_rho_core_rspace, &
               qs_charges%total_rho_core_rspace - REAL(nelectron_total + dft_control%charge, dp)

            IF (dft_control%correct_surf_dip) THEN
               WRITE (UNIT=output_unit, FMT="((T3,A,/,T3,A,T41,F20.10))") &
                  "Total dipole moment perpendicular to ", &
                  "the slab [electrons-Angstroem]: ", &
                  qs_env%surface_dipole_moment
            END IF

            IF (gapw) THEN
               tot1_h = qs_charges%total_rho1_hard(1)
               tot1_s = qs_charges%total_rho1_soft(1)
               DO ispin = 2, dft_control%nspins
                  tot1_h = tot1_h + qs_charges%total_rho1_hard(ispin)
                  tot1_s = tot1_s + qs_charges%total_rho1_soft(ispin)
               END DO
               WRITE (UNIT=output_unit, FMT="((T3,A,T41,2F20.10))") &
                  "Hard and soft densities (Lebedev):", &
                  tot1_h, tot1_s
               WRITE (UNIT=output_unit, FMT="(T3,A,T41,F20.10)") &
                  "Total Rho_soft + Rho1_hard - Rho1_soft (r-space): ", &
                  accurate_sum(tot_rho_r) + tot1_h - tot1_s, &
                  "Total charge density (r-space):      ", &
                  accurate_sum(tot_rho_r) + tot1_h - tot1_s &
                  + qs_charges%total_rho_core_rspace, &
                  "Total Rho_soft + Rho0_soft (g-space):", &
                  qs_charges%total_rho_gspace
            ELSE
               WRITE (UNIT=output_unit, FMT="(T3,A,T41,F20.10)") &
                  "Total charge density on r-space grids:     ", &
                  accurate_sum(tot_rho_r) + &
                  qs_charges%total_rho_core_rspace, &
                  "Total charge density g-space grids:     ", &
                  qs_charges%total_rho_gspace
            END IF
         END IF
         IF (dft_control%qs_control%semi_empirical) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
               "Core-core repulsion energy [eV]:               ", energy%core_overlap*evolt, &
               "Core Hamiltonian energy [eV]:                  ", energy%core*evolt, &
               "Two-electron integral energy [eV]:             ", energy%hartree*evolt, &
               "Electronic energy [eV]:                        ", &
               (energy%core + 0.5_dp*energy%hartree)*evolt
            IF (energy%dispersion /= 0.0_dp) &
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "Dispersion energy [eV]:                     ", energy%dispersion*evolt
         ELSEIF (dft_control%qs_control%dftb) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
               "Core Hamiltonian energy:                       ", energy%core, &
               "Repulsive potential energy:                    ", energy%repulsive, &
               "Electronic energy:                             ", energy%hartree, &
               "Dispersion energy:                             ", energy%dispersion
            IF (energy%dftb3 /= 0.0_dp) &
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "DFTB3 3rd order energy:                     ", energy%dftb3
            IF (energy%efield /= 0.0_dp) &
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "Electric field interaction energy:          ", energy%efield
         ELSEIF (dft_control%qs_control%xtb) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
               "Core Hamiltonian energy:                       ", energy%core, &
               "Repulsive potential energy:                    ", energy%repulsive, &
               "Electronic energy:                             ", energy%hartree, &
               "DFTB3 3rd order energy:                        ", energy%dftb3, &
               "Dispersion energy:                             ", energy%dispersion
            IF (dft_control%qs_control%xtb_control%xb_interaction) &
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "Correction for halogen bonding:                ", energy%xtb_xb_inter
            IF (dft_control%qs_control%xtb_control%do_nonbonded) &
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "Correction for nonbonded interactions:         ", energy%xtb_nonbonded
            IF (energy%efield /= 0.0_dp) &
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "Electric field interaction energy:          ", energy%efield
         ELSE
            IF (dft_control%do_admm) THEN
               exc_energy = energy%exc + energy%exc_aux_fit
               IF (gapw .OR. gapw_xc) exc1_energy = energy%exc1 + energy%exc1_aux_fit
            ELSE
               exc_energy = energy%exc
               IF (gapw .OR. gapw_xc) exc1_energy = energy%exc1
            END IF

            IF (psolver .EQ. pw_poisson_implicit) THEN
               implicit_ps_ehartree = pw_env%poisson_env%implicit_env%ehartree
               bc = pw_env%poisson_env%parameters%ps_implicit_params%boundary_condition
               SELECT CASE (bc)
               CASE (MIXED_PERIODIC_BC, MIXED_BC)
                  WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
                     "Overlap energy of the core charge distribution:", energy%core_overlap, &
                     "Self energy of the core charge distribution:   ", energy%core_self, &
                     "Core Hamiltonian energy:                       ", energy%core, &
                     "Hartree energy:                                ", implicit_ps_ehartree, &
                     "Electric enthalpy:                             ", energy%hartree, &
                     "Exchange-correlation energy:                   ", exc_energy
               CASE (PERIODIC_BC, NEUMANN_BC)
                  WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
                     "Overlap energy of the core charge distribution:", energy%core_overlap, &
                     "Self energy of the core charge distribution:   ", energy%core_self, &
                     "Core Hamiltonian energy:                       ", energy%core, &
                     "Hartree energy:                                ", energy%hartree, &
                     "Exchange-correlation energy:                   ", exc_energy
               END SELECT
            ELSE
               WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
                  "Overlap energy of the core charge distribution:", energy%core_overlap, &
                  "Self energy of the core charge distribution:   ", energy%core_self, &
                  "Core Hamiltonian energy:                       ", energy%core, &
                  "Hartree energy:                                ", energy%hartree, &
                  "Exchange-correlation energy:                   ", exc_energy
            END IF
            IF (energy%e_hartree /= 0.0_dp) &
               WRITE (UNIT=output_unit, FMT="(T3,A,/,T3,A,T56,F25.14)") &
               "Coulomb Electron-Electron Interaction Energy ", &
               "- Already included in the total Hartree term ", energy%e_hartree
            IF (energy%ex /= 0.0_dp) &
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "Hartree-Fock Exchange energy:                  ", energy%ex
            IF (energy%dispersion /= 0.0_dp) &
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "Dispersion energy:                             ", energy%dispersion
            IF (energy%gcp /= 0.0_dp) &
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "gCP energy:                                    ", energy%gcp
            IF (gapw) THEN
               WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
                  "GAPW| Exc from hard and soft atomic rho1:      ", exc1_energy, &
                  "GAPW| local Eh = 1 center integrals:           ", energy%hartree_1c
            END IF
            IF (gapw_xc) THEN
               WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
                  "GAPW_XC| Exc from hard and soft atomic rho1:      ", exc1_energy
            END IF
         END IF
         IF (dft_control%smear) THEN
            WRITE (UNIT=output_unit, FMT="((T3,A,T56,F25.14))") &
               "Electronic entropic energy:", energy%kTS
            WRITE (UNIT=output_unit, FMT="((T3,A,T56,F25.14))") &
               "Fermi energy:", energy%efermi
         END IF
         IF (dft_control%dft_plus_u) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
               "DFT+U energy:", energy%dft_plus_u
         END IF
         IF (dft_control%do_sccs) THEN
            WRITE (UNIT=output_unit, FMT="(/,T3,A,T56,F25.14)") &
               "SCCS| Hartree energy of solute and solvent [Hartree]", energy%sccs_hartree
            WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14,/,T3,A,T61,F20.3)") &
               "SCCS| Polarisation energy                  [Hartree]", energy%sccs_pol, &
               "SCCS|                                      [kcal/mol]", &
               cp_unit_from_cp2k(energy%sccs_pol, "kcalmol")
         END IF
         IF (qmmm) THEN
            WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "QM/MM Electrostatic energy:                    ", energy%qmmm_el
            IF (qs_env%qmmm_env_qm%image_charge) THEN
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
                  "QM/MM image charge energy:                ", energy%image_charge
            END IF
         END IF
         IF (dft_control%qs_control%mulliken_restraint) THEN
            WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "Mulliken restraint energy: ", energy%mulliken
         END IF
         IF (qs_env%excited_state) THEN
            IF (energy%excited_state /= 0.0_dp) &
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "Excited State energy: ", energy%excited_state
         END IF
         IF (dft_control%qs_control%semi_empirical) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
               "Total energy [eV]:                             ", energy%total*evolt
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
               "Atomic reference energy [eV]:                  ", energy%core_self*evolt, &
               "Heat of formation [kcal/mol]:                  ", &
               (energy%total + energy%core_self)*kcalmol
         ELSE
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
               "Total energy:                                  ", energy%total
         END IF
         IF (qmmm) THEN
            IF (qs_env%qmmm_env_qm%image_charge) THEN
               CALL print_image_coefficients(qs_env%image_coeff, qs_env)
            END IF
         END IF
         CALL m_flush(output_unit)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_print_scf_summary

! **************************************************************************************************
!> \brief collects the 'heavy duty' printing tasks out of the SCF loop
!> \param qs_env ...
!> \param scf_env ...
!> \param para_env ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_print(qs_env, scf_env, para_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'qs_scf_loop_print'

      INTEGER                                            :: after, handle, ic, ispin, iw
      LOGICAL                                            :: do_kpoints, omit_headers
      REAL(KIND=dp)                                      :: mo_mag_max, mo_mag_min, orthonormality
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_p, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

      logger => cp_get_default_logger()
      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, input=input, dft_control=dft_control, &
                      do_kpoints=do_kpoints)

      dft_section => section_vals_get_subs_vals(input, "DFT")
      scf_section => section_vals_get_subs_vals(dft_section, "SCF")

      CALL section_vals_val_get(input, "DFT%PRINT%AO_MATRICES%OMIT_HEADERS", l_val=omit_headers)
      DO ispin = 1, dft_control%nspins

         IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                              dft_section, "PRINT%AO_MATRICES/DENSITY"), cp_p_file)) THEN
            CALL get_qs_env(qs_env, rho=rho)
            CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
            iw = cp_print_key_unit_nr(logger, dft_section, "PRINT%AO_MATRICES/DENSITY", &
                                      extension=".Log")
            CALL section_vals_val_get(dft_section, "PRINT%AO_MATRICES%NDIGITS", i_val=after)
            after = MIN(MAX(after, 1), 16)
            DO ic = 1, SIZE(matrix_p, 2)
               CALL cp_dbcsr_write_sparse_matrix(matrix_p(ispin, ic)%matrix, 4, after, qs_env, para_env, &
                                                 output_unit=iw, omit_headers=omit_headers)
            END DO
            CALL cp_print_key_finished_output(iw, logger, dft_section, &
                                              "PRINT%AO_MATRICES/DENSITY")
         END IF

         IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                              dft_section, "PRINT%AO_MATRICES/KOHN_SHAM_MATRIX"), cp_p_file)) THEN
            iw = cp_print_key_unit_nr(logger, dft_section, "PRINT%AO_MATRICES/KOHN_SHAM_MATRIX", &
                                      extension=".Log")
            CALL section_vals_val_get(dft_section, "PRINT%AO_MATRICES%NDIGITS", i_val=after)
            after = MIN(MAX(after, 1), 16)
            CALL get_qs_env(qs_env=qs_env, matrix_ks_kp=matrix_ks)
            DO ic = 1, SIZE(matrix_ks, 2)
               IF (dft_control%qs_control%semi_empirical) THEN
                  CALL cp_dbcsr_write_sparse_matrix(matrix_ks(ispin, ic)%matrix, 4, after, qs_env, para_env, &
                                                    scale=evolt, output_unit=iw, omit_headers=omit_headers)
               ELSE
                  CALL cp_dbcsr_write_sparse_matrix(matrix_ks(ispin, ic)%matrix, 4, after, qs_env, para_env, &
                                                    output_unit=iw, omit_headers=omit_headers)
               END IF
            END DO
            CALL cp_print_key_finished_output(iw, logger, dft_section, &
                                              "PRINT%AO_MATRICES/KOHN_SHAM_MATRIX")
         END IF

      END DO

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           scf_section, "PRINT%MO_ORTHONORMALITY"), cp_p_file)) THEN
         IF (do_kpoints) THEN
            iw = cp_print_key_unit_nr(logger, scf_section, "PRINT%MO_ORTHONORMALITY", &
                                      extension=".scfLog")
            IF (iw > 0) THEN
               WRITE (iw, '(T8,A)') &
                  " K-points: Maximum deviation from MO S-orthonormality not determined"
            END IF
            CALL cp_print_key_finished_output(iw, logger, scf_section, &
                                              "PRINT%MO_ORTHONORMALITY")
         ELSE
            CALL get_qs_env(qs_env, mos=mos)
            IF (scf_env%method == special_diag_method_nr) THEN
               CALL calculate_orthonormality(orthonormality, mos)
            ELSE
               CALL get_qs_env(qs_env=qs_env, matrix_s_kp=matrix_s)
               CALL calculate_orthonormality(orthonormality, mos, matrix_s(1, 1)%matrix)
            END IF
            iw = cp_print_key_unit_nr(logger, scf_section, "PRINT%MO_ORTHONORMALITY", &
                                      extension=".scfLog")
            IF (iw > 0) THEN
               WRITE (iw, '(T8,A,T61,E20.4)') &
                  " Maximum deviation from MO S-orthonormality", orthonormality
            END IF
            CALL cp_print_key_finished_output(iw, logger, scf_section, &
                                              "PRINT%MO_ORTHONORMALITY")
         END IF
      END IF
      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           scf_section, "PRINT%MO_MAGNITUDE"), cp_p_file)) THEN
         IF (do_kpoints) THEN
            iw = cp_print_key_unit_nr(logger, scf_section, "PRINT%MO_MAGNITUDE", &
                                      extension=".scfLog")
            IF (iw > 0) THEN
               WRITE (iw, '(T8,A)') &
                  " K-points: Minimum/Maximum MO magnitude not determined"
            END IF
            CALL cp_print_key_finished_output(iw, logger, scf_section, &
                                              "PRINT%MO_MAGNITUDE")
         ELSE
            CALL get_qs_env(qs_env, mos=mos)
            CALL calculate_magnitude(mos, mo_mag_min, mo_mag_max)
            iw = cp_print_key_unit_nr(logger, scf_section, "PRINT%MO_MAGNITUDE", &
                                      extension=".scfLog")
            IF (iw > 0) THEN
               WRITE (iw, '(T8,A,T41,2E20.4)') &
                  " Minimum/Maximum MO magnitude ", mo_mag_min, mo_mag_max
            END IF
            CALL cp_print_key_finished_output(iw, logger, scf_section, &
                                              "PRINT%MO_MAGNITUDE")
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_loop_print

! **************************************************************************************************
!> \brief writes CDFT constraint information and optionally CDFT scf loop info
!> \param output_unit where to write the information
!> \param scf_control settings of the SCF loop
!> \param scf_env the env which holds convergence data
!> \param cdft_control the env which holds information about the constraint
!> \param energy the total energy
!> \param total_steps the total number of performed SCF iterations
!> \param should_stop if the calculation should stop
!> \param outer_loop_converged logical which determines if the CDFT SCF loop converged
!> \param cdft_loop logical which determines a CDFT SCF loop is active
!> \par History
!>      12.2015 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE qs_scf_cdft_info(output_unit, scf_control, scf_env, cdft_control, &
                               energy, total_steps, should_stop, outer_loop_converged, &
                               cdft_loop)
      INTEGER                                            :: output_unit
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(cdft_control_type), POINTER                   :: cdft_control
      TYPE(qs_energy_type), POINTER                      :: energy
      INTEGER                                            :: total_steps
      LOGICAL, INTENT(IN)                                :: should_stop, outer_loop_converged, &
                                                            cdft_loop

      REAL(KIND=dp)                                      :: outer_loop_eps

      IF (cdft_loop) THEN
         outer_loop_eps = SQRT(MAXVAL(scf_env%outer_scf%gradient(:, scf_env%outer_scf%iter_count)**2))
         IF (output_unit > 0) WRITE (output_unit, '(/,T3,A,I4,A,E10.2,A,F22.10)') &
            "CDFT SCF iter =  ", scf_env%outer_scf%iter_count, &
            " RMS gradient = ", outer_loop_eps, " energy =", energy%total
         IF (outer_loop_converged) THEN
            IF (output_unit > 0) WRITE (output_unit, '(T3,A,I4,A,I4,A,/)') &
               "CDFT SCF loop converged in", scf_env%outer_scf%iter_count, &
               " iterations or ", total_steps, " steps"
         END IF
         IF ((scf_env%outer_scf%iter_count > scf_control%outer_scf%max_scf .OR. should_stop) &
             .AND. .NOT. outer_loop_converged) THEN
            IF (output_unit > 0) WRITE (output_unit, '(T3,A,I4,A,I4,A,/)') &
               "CDFT SCF loop FAILED to converge after ", &
               scf_env%outer_scf%iter_count, " iterations or ", total_steps, " steps"
         END IF
      END IF
      CALL qs_scf_cdft_constraint_info(output_unit, cdft_control)

   END SUBROUTINE qs_scf_cdft_info

! **************************************************************************************************
!> \brief writes information about the CDFT env
!> \param output_unit where to write the information
!> \param cdft_control the CDFT env that stores information about the constraint calculation
!> \par History
!>      12.2015 created [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE qs_scf_cdft_initial_info(output_unit, cdft_control)
      INTEGER                                            :: output_unit
      TYPE(cdft_control_type), POINTER                   :: cdft_control

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(/,A)') &
            "  ---------------------------------- CDFT --------------------------------------"
         WRITE (output_unit, '(A)') &
            "  Optimizing a density constraint in an external SCF loop "
         WRITE (output_unit, '(A)') "  "
         SELECT CASE (cdft_control%type)
         CASE (outer_scf_hirshfeld_constraint)
            WRITE (output_unit, '(A)') "  Type of constraint:     Hirshfeld"
         CASE (outer_scf_becke_constraint)
            WRITE (output_unit, '(A)') "  Type of constraint:         Becke"
         END SELECT
         WRITE (output_unit, '(A,I8)') "  Number of constraints:   ", SIZE(cdft_control%group)
         WRITE (output_unit, '(A,L8)') "  Using fragment densities:", cdft_control%fragment_density
         WRITE (output_unit, '(A)') "  "
         IF (cdft_control%atomic_charges) WRITE (output_unit, '(A,/)') "  Calculating atomic CDFT charges"
         SELECT CASE (cdft_control%constraint_control%optimizer)
         CASE (outer_scf_optimizer_sd)
            WRITE (output_unit, '(A)') &
               "  Minimizer               : SD                  : steepest descent"
         CASE (outer_scf_optimizer_diis)
            WRITE (output_unit, '(A)') &
               "  Minimizer               : DIIS                : direct inversion"
            WRITE (output_unit, '(A)') &
               "                                                       in the iterative subspace"
            WRITE (output_unit, '(A,I3,A)') &
               "                                                  using ", &
               cdft_control%constraint_control%diis_buffer_length, " DIIS vectors"
         CASE (outer_scf_optimizer_bisect)
            WRITE (output_unit, '(A)') &
               "  Minimizer               : BISECT              : gradient bisection"
            WRITE (output_unit, '(A,I3)') &
               "                                                  using a trust count of", &
               cdft_control%constraint_control%bisect_trust_count
         CASE (outer_scf_optimizer_broyden, outer_scf_optimizer_newton, &
               outer_scf_optimizer_newton_ls)
            CALL cdft_opt_type_write(cdft_control%constraint_control%cdft_opt_control, &
                                     cdft_control%constraint_control%optimizer, output_unit)
         CASE (outer_scf_optimizer_secant)
            WRITE (output_unit, '(A)') "  Minimizer               : Secant"
         CASE DEFAULT
            CPABORT("")
         END SELECT
         WRITE (output_unit, '(/,A,L7)') &
            "  Reusing OT preconditioner: ", cdft_control%reuse_precond
         IF (cdft_control%reuse_precond) THEN
            WRITE (output_unit, '(A,I3,A,I3,A)') &
               "       using old preconditioner for up to ", &
               cdft_control%max_reuse, " subsequent CDFT SCF"
            WRITE (output_unit, '(A,I3,A,I3,A)') &
               "       iterations if the relevant loop converged in less than ", &
               cdft_control%precond_freq, " steps"
         END IF
         SELECT CASE (cdft_control%type)
         CASE (outer_scf_hirshfeld_constraint)
            WRITE (output_unit, '(/,A)') "  Hirshfeld constraint settings"
            WRITE (output_unit, '(A)') "  "
            SELECT CASE (cdft_control%hirshfeld_control%shape_function)
            CASE (shape_function_gaussian)
               WRITE (output_unit, '(A, A8)') &
                  "  Shape function type:     ", "Gaussian"
               WRITE (output_unit, '(A)', ADVANCE='NO') &
                  "  Type of Gaussian:   "
               SELECT CASE (cdft_control%hirshfeld_control%gaussian_shape)
               CASE (radius_default)
                  WRITE (output_unit, '(A13)') "Default"
               CASE (radius_covalent)
                  WRITE (output_unit, '(A13)') "Covalent"
               CASE (radius_single)
                  WRITE (output_unit, '(A13)') "Fixed radius"
               CASE (radius_vdw)
                  WRITE (output_unit, '(A13)') "Van der Waals"
               CASE (radius_user)
                  WRITE (output_unit, '(A13)') "User-defined"

               END SELECT
            CASE (shape_function_density)
               WRITE (output_unit, '(A, A8)') &
                  "  Shape function type:     ", "Density"
            END SELECT
         CASE (outer_scf_becke_constraint)
            WRITE (output_unit, '(/, A)') "  Becke constraint settings"
            WRITE (output_unit, '(A)') "  "
            SELECT CASE (cdft_control%becke_control%cutoff_type)
            CASE (becke_cutoff_global)
               WRITE (output_unit, '(A,F8.3,A)') &
                  "  Cutoff for partitioning :", cp_unit_from_cp2k(cdft_control%becke_control%rglobal, &
                                                                   "angstrom"), " angstrom"
            CASE (becke_cutoff_element)
               WRITE (output_unit, '(A)') &
                  "  Using element specific cutoffs for partitioning"
            END SELECT
            WRITE (output_unit, '(A,L7)') &
               "  Skipping distant gpoints: ", cdft_control%becke_control%should_skip
            WRITE (output_unit, '(A,L7)') &
               "  Precompute gradients    : ", cdft_control%becke_control%in_memory
            WRITE (output_unit, '(A)') "  "
            IF (cdft_control%becke_control%adjust) &
               WRITE (output_unit, '(A)') &
               "  Using atomic radii to generate a heteronuclear charge partitioning"
            WRITE (output_unit, '(A)') "  "
            IF (.NOT. cdft_control%becke_control%cavity_confine) THEN
               WRITE (output_unit, '(A)') &
                  "  No confinement is active"
            ELSE
               WRITE (output_unit, '(A)') "  Confinement using a Gaussian shaped cavity is active"
               SELECT CASE (cdft_control%becke_control%cavity_shape)
               CASE (radius_single)
                  WRITE (output_unit, '(A,F8.4, A)') &
                     "  Type of Gaussian        : Fixed radius: ", &
                     cp_unit_from_cp2k(cdft_control%becke_control%rcavity, "angstrom"), " angstrom"
               CASE (radius_covalent)
                  WRITE (output_unit, '(A)') &
                     "  Type of Gaussian        : Covalent radius "
               CASE (radius_vdw)
                  WRITE (output_unit, '(A)') &
                     "  Type of Gaussian        : vdW radius "
               CASE (radius_user)
                  WRITE (output_unit, '(A)') &
                     "  Type of Gaussian        : User radius "
               END SELECT
               WRITE (output_unit, '(A,ES12.4)') &
                  "  Cavity threshold        : ", cdft_control%becke_control%eps_cavity
            END IF
         END SELECT
         WRITE (output_unit, '(/,A)') &
            "  ---------------------------------- CDFT --------------------------------------"
      END IF

   END SUBROUTINE qs_scf_cdft_initial_info

! **************************************************************************************************
!> \brief writes CDFT constraint information
!> \param output_unit where to write the information
!> \param cdft_control the env which holds information about the constraint
!> \par History
!>      08.2018 separated from qs_scf_cdft_info to make code callable elsewhere  [Nico Holmberg]
! **************************************************************************************************
   SUBROUTINE qs_scf_cdft_constraint_info(output_unit, cdft_control)
      INTEGER                                            :: output_unit
      TYPE(cdft_control_type), POINTER                   :: cdft_control

      INTEGER                                            :: igroup

      IF (output_unit > 0) THEN
         SELECT CASE (cdft_control%type)
         CASE (outer_scf_hirshfeld_constraint)
            WRITE (output_unit, '(/,T3,A,T60)') &
               '------------------- Hirshfeld constraint information -------------------'
         CASE (outer_scf_becke_constraint)
            WRITE (output_unit, '(/,T3,A,T60)') &
               '--------------------- Becke constraint information ---------------------'
         CASE DEFAULT
            CPABORT("Unknown CDFT constraint.")
         END SELECT
         DO igroup = 1, SIZE(cdft_control%target)
            IF (igroup > 1) WRITE (output_unit, '(T3,A)') ' '
            WRITE (output_unit, '(T3,A,T54,(3X,I18))') &
               'Atomic group                :', igroup
            SELECT CASE (cdft_control%group(igroup)%constraint_type)
            CASE (cdft_charge_constraint)
               IF (cdft_control%group(igroup)%is_fragment_constraint) THEN
                  WRITE (output_unit, '(T3,A,T42,A)') &
                     'Type of constraint          :', ADJUSTR('Charge density constraint (frag.)')
               ELSE
                  WRITE (output_unit, '(T3,A,T50,A)') &
                     'Type of constraint          :', ADJUSTR('Charge density constraint')
               END IF
            CASE (cdft_magnetization_constraint)
               IF (cdft_control%group(igroup)%is_fragment_constraint) THEN
                  WRITE (output_unit, '(T3,A,T35,A)') &
                     'Type of constraint          :', ADJUSTR('Magnetization density constraint (frag.)')
               ELSE
                  WRITE (output_unit, '(T3,A,T43,A)') &
                     'Type of constraint          :', ADJUSTR('Magnetization density constraint')
               END IF
            CASE (cdft_alpha_constraint)
               IF (cdft_control%group(igroup)%is_fragment_constraint) THEN
                  WRITE (output_unit, '(T3,A,T38,A)') &
                     'Type of constraint          :', ADJUSTR('Alpha spin density constraint (frag.)')
               ELSE
                  WRITE (output_unit, '(T3,A,T46,A)') &
                     'Type of constraint          :', ADJUSTR('Alpha spin density constraint')
               END IF
            CASE (cdft_beta_constraint)
               IF (cdft_control%group(igroup)%is_fragment_constraint) THEN
                  WRITE (output_unit, '(T3,A,T39,A)') &
                     'Type of constraint          :', ADJUSTR('Beta spin density constraint (frag.)')
               ELSE
                  WRITE (output_unit, '(T3,A,T47,A)') &
                     'Type of constraint          :', ADJUSTR('Beta spin density constraint')
               END IF
            CASE DEFAULT
               CPABORT("Unknown constraint type.")
            END SELECT
            WRITE (output_unit, '(T3,A,T54,(3X,F18.12))') &
               'Target value of constraint  :', cdft_control%target(igroup)
            WRITE (output_unit, '(T3,A,T54,(3X,F18.12))') &
               'Current value of constraint :', cdft_control%value(igroup)
            WRITE (output_unit, '(T3,A,T59,(3X,ES13.3))') &
               'Deviation from target       :', cdft_control%value(igroup) - cdft_control%target(igroup)
            WRITE (output_unit, '(T3,A,T54,(3X,F18.12))') &
               'Strength of constraint      :', cdft_control%strength(igroup)
         END DO
         WRITE (output_unit, '(T3,A)') &
            '------------------------------------------------------------------------'
      END IF

   END SUBROUTINE qs_scf_cdft_constraint_info

END MODULE qs_scf_output
