!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_restart
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_trace
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pool_create_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_read_unformatted,&
                                              cp_fm_release,&
                                              cp_fm_type,&
                                              cp_fm_write_formatted,&
                                              cp_fm_write_info,&
                                              cp_fm_write_unformatted
   USE cp_log_handling,                 ONLY: cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: tddfpt_ground_state_mos
   USE string_utilities,                ONLY: integer_to_string
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_restart'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   ! number of first derivative components (3: d/dx, d/dy, d/dz)
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC :: tddfpt_write_restart, tddfpt_read_restart, tddfpt_write_newtonx_output, tddfpt_check_orthonormality

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Write Ritz vectors to a binary restart file.
!> \param evects               vectors to store
!> \param evals                TDDFPT eigenvalues
!> \param gs_mos               structure that holds ground state occupied and virtual
!>                             molecular orbitals
!> \param logger               a logger object
!> \param tddfpt_print_section TDDFPT%PRINT input section
!> \par History
!>    * 08.2016 created [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_write_restart(evects, evals, gs_mos, logger, tddfpt_print_section)
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(in)      :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_write_restart'

      INTEGER                                            :: handle, ispin, istate, nao, nspins, &
                                                            nstates, ounit
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ

      IF (BTEST(cp_print_key_should_output(logger%iter_info, tddfpt_print_section, "RESTART"), cp_p_file)) THEN
         CALL timeset(routineN, handle)

         nspins = SIZE(evects, 1)
         nstates = SIZE(evects, 2)

         IF (debug_this_module) THEN
            CPASSERT(SIZE(evals) == nstates)
            CPASSERT(nspins > 0)
            CPASSERT(nstates > 0)
         END IF

         CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=nao)

         DO ispin = 1, nspins
            nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
         END DO

         ounit = cp_print_key_unit_nr(logger, tddfpt_print_section, "RESTART", &
                                      extension=".tdwfn", file_status="REPLACE", file_action="WRITE", &
                                      do_backup=.TRUE., file_form="UNFORMATTED")

         IF (ounit > 0) THEN
            WRITE (ounit) nstates, nspins, nao
            WRITE (ounit) nmo_occ(1:nspins)
            WRITE (ounit) evals
         END IF

         DO istate = 1, nstates
            DO ispin = 1, nspins
               ! TDDFPT wave function is actually stored as a linear combination of virtual MOs
               ! that replaces the corresponding deoccupied MO. Unfortunately, the phase
               ! of the occupied MOs varies depending on the eigensolver used as well as
               ! how eigenvectors are distributed across computational cores. The phase is important
               ! because TDDFPT wave functions are used to compute a response electron density
               ! \rho^{-} = 1/2 * [C_{0} * evect^T + evect * C_{0}^{-}], where C_{0} is the expansion
               ! coefficients of the reference ground-state wave function. To make the restart file
               ! transferable, TDDFPT wave functions are stored in assumption that all ground state
               ! MOs have a positive phase.
               CALL cp_fm_column_scale(evects(ispin, istate), gs_mos(ispin)%phases_occ)

               CALL cp_fm_write_unformatted(evects(ispin, istate), ounit)

               CALL cp_fm_column_scale(evects(ispin, istate), gs_mos(ispin)%phases_occ)
            END DO
         END DO

         CALL cp_print_key_finished_output(ounit, logger, tddfpt_print_section, "RESTART")

         CALL timestop(handle)
      END IF

   END SUBROUTINE tddfpt_write_restart

! **************************************************************************************************
!> \brief Initialise initial guess vectors by reading (un-normalised) Ritz vectors
!>        from a binary restart file.
!> \param evects               vectors to initialise (initialised on exit)
!> \param evals                TDDFPT eigenvalues (initialised on exit)
!> \param gs_mos               structure that holds ground state occupied and virtual
!>                             molecular orbitals
!> \param logger               a logger object
!> \param tddfpt_section       TDDFPT input section
!> \param tddfpt_print_section TDDFPT%PRINT input section
!> \param fm_pool_ao_mo_occ    pools of dense matrices with shape [nao x nmo_occ(spin)]
!> \param blacs_env_global     BLACS parallel environment involving all the processor
!> \return the number of excited states found in the restart file
!> \par History
!>    * 08.2016 created [Sergey Chulkov]
! **************************************************************************************************
   FUNCTION tddfpt_read_restart(evects, evals, gs_mos, logger, tddfpt_section, tddfpt_print_section, &
                                fm_pool_ao_mo_occ, blacs_env_global) RESULT(nstates_read)
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(inout)   :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(out)           :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: tddfpt_section, tddfpt_print_section
      TYPE(cp_fm_pool_p_type), DIMENSION(:), INTENT(in)  :: fm_pool_ao_mo_occ
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      INTEGER                                            :: nstates_read

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_read_restart'

      CHARACTER(len=20)                                  :: read_str, ref_str
      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: handle, ispin, istate, iunit, n_rep_val, &
                                                            nao, nao_read, nspins, nspins_read, &
                                                            nstates
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ, nmo_occ_read
      LOGICAL                                            :: file_exists
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals_read
      TYPE(mp_para_env_type), POINTER                    :: para_env_global
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(tddfpt_section))

      ! generate restart file name
      CALL section_vals_val_get(tddfpt_section, "WFN_RESTART_FILE_NAME", n_rep_val=n_rep_val)
      IF (n_rep_val > 0) THEN
         CALL section_vals_val_get(tddfpt_section, "WFN_RESTART_FILE_NAME", c_val=filename)
      ELSE
         print_key => section_vals_get_subs_vals(tddfpt_print_section, "RESTART")
         filename = cp_print_key_generate_filename(logger, print_key, &
                                                   extension=".tdwfn", my_local=.FALSE.)
      END IF

      CALL blacs_env_global%get(para_env=para_env_global)

      IF (para_env_global%is_source()) THEN
         INQUIRE (FILE=filename, exist=file_exists)

         IF (.NOT. file_exists) THEN
            nstates_read = 0
            CALL para_env_global%bcast(nstates_read)

            CALL cp_warn(__LOCATION__, &
                         "User requested to restart the TDDFPT wave functions from the file '"//TRIM(filename)// &
                         "' which does not exist. Guess wave functions will be constructed using Kohn-Sham orbitals.")
            CALL timestop(handle)
            RETURN
         END IF

         CALL open_file(file_name=filename, file_action="READ", file_form="UNFORMATTED", &
                        file_status="OLD", unit_number=iunit)
      END IF

      nspins = SIZE(evects, 1)
      nstates = SIZE(evects, 2)
      CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=nao)

      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
      END DO

      IF (para_env_global%is_source()) THEN
         READ (iunit) nstates_read, nspins_read, nao_read

         IF (nspins_read /= nspins) THEN
            CALL integer_to_string(nspins, ref_str)
            CALL integer_to_string(nspins_read, read_str)
            CALL cp_abort(__LOCATION__, &
                          "Restarted TDDFPT wave function contains incompatible number of spin components ("// &
                          TRIM(read_str)//" instead of "//TRIM(ref_str)//").")
         END IF

         IF (nao_read /= nao) THEN
            CALL integer_to_string(nao, ref_str)
            CALL integer_to_string(nao_read, read_str)
            CALL cp_abort(__LOCATION__, &
                          "Incompatible number of atomic orbitals ("//TRIM(read_str)//" instead of "//TRIM(ref_str)//").")
         END IF

         READ (iunit) nmo_occ_read(1:nspins)

         DO ispin = 1, nspins
            IF (nmo_occ_read(ispin) /= nmo_occ(ispin)) THEN
               CALL cp_abort(__LOCATION__, &
                             "Incompatible number of electrons and/or multiplicity.")
            END IF
         END DO

         IF (nstates_read /= nstates) THEN
            CALL integer_to_string(nstates, ref_str)
            CALL integer_to_string(nstates_read, read_str)
            CALL cp_warn(__LOCATION__, &
                         "TDDFPT restart file contains "//TRIM(read_str)// &
                         " wave function(s) however "//TRIM(ref_str)// &
                         " excited states were requested.")
         END IF
      END IF
      CALL para_env_global%bcast(nstates_read)

      ! exit if restart file does not exist
      IF (nstates_read <= 0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      IF (para_env_global%is_source()) THEN
         ALLOCATE (evals_read(nstates_read))
         READ (iunit) evals_read
         IF (nstates_read <= nstates) THEN
            evals(1:nstates_read) = evals_read(1:nstates_read)
         ELSE
            evals(1:nstates) = evals_read(1:nstates)
         END IF
         DEALLOCATE (evals_read)
      END IF
      CALL para_env_global%bcast(evals)

      DO istate = 1, nstates_read
         DO ispin = 1, nspins
            IF (istate <= nstates) THEN
               CALL fm_pool_create_fm(fm_pool_ao_mo_occ(ispin)%pool, evects(ispin, istate))

               CALL cp_fm_read_unformatted(evects(ispin, istate), iunit)

               CALL cp_fm_column_scale(evects(ispin, istate), gs_mos(ispin)%phases_occ)
            END IF
         END DO
      END DO

      IF (para_env_global%is_source()) &
         CALL close_file(unit_number=iunit)

      CALL timestop(handle)

   END FUNCTION tddfpt_read_restart
! **************************************************************************************************
!> \brief Write Ritz vectors to a binary restart file.
!> \param evects               vectors to store
!> \param evals                TDDFPT eigenvalues
!> \param gs_mos               structure that holds ground state occupied and virtual
!>                             molecular orbitals
!> \param logger               a logger object
!> \param tddfpt_print_section TDDFPT%PRINT input section
!> \param matrix_s ...
!> \param S_evects ...
!> \param sub_env ...
! **************************************************************************************************
   SUBROUTINE tddfpt_write_newtonx_output(evects, evals, gs_mos, logger, tddfpt_print_section, &
                                          matrix_s, S_evects, sub_env)

      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(in)      :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(cp_logger_type), INTENT(in), POINTER          :: logger
      TYPE(section_vals_type), INTENT(in), POINTER       :: tddfpt_print_section
      TYPE(dbcsr_type), INTENT(in), POINTER              :: matrix_s
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(in)      :: S_evects
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_write_newtonx_output'

      INTEGER                                            :: handle, iocc, ispin, istate, ivirt, nao, &
                                                            nspins, nstates, ounit
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ, nmo_virt
      LOGICAL                                            :: print_phases, print_virtuals, &
                                                            scale_with_phases
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: phase_evects
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: evects_mo

      IF (BTEST(cp_print_key_should_output(logger%iter_info, tddfpt_print_section, "NAMD_PRINT"), cp_p_file)) THEN
         CALL timeset(routineN, handle)
         CALL section_vals_val_get(tddfpt_print_section, "NAMD_PRINT%PRINT_VIRTUALS", l_val=print_virtuals)
         CALL section_vals_val_get(tddfpt_print_section, "NAMD_PRINT%PRINT_PHASES", l_val=print_phases)
         CALL section_vals_val_get(tddfpt_print_section, "NAMD_PRINT%SCALE_WITH_PHASES", l_val=scale_with_phases)

         nspins = SIZE(evects, 1)
         nstates = SIZE(evects, 2)

         IF (debug_this_module) THEN
            CPASSERT(SIZE(evals) == nstates)
            CPASSERT(nspins > 0)
            CPASSERT(nstates > 0)
         END IF

         CALL cp_fm_get_info(gs_mos(1)%mos_occ, nrow_global=nao)

         IF (sub_env%is_split) THEN
            CALL cp_abort(__LOCATION__, "NEWTONX interface print not possible when states"// &
                          " are distributed to different CPU pools.")
         END IF

         ounit = cp_print_key_unit_nr(logger, tddfpt_print_section, "NAMD_PRINT", &
                                      extension=".inp", file_form="FORMATTED", file_action="WRITE", file_status="REPLACE")
         IF (debug_this_module) CALL tddfpt_check_orthonormality(evects, ounit, S_evects, matrix_s)

         ! print eigenvectors
         IF (print_virtuals) THEN
            ALLOCATE (evects_mo(nspins, nstates))
            DO istate = 1, nstates
               DO ispin = 1, nspins

                  ! transform eigenvectors
                  NULLIFY (fmstruct)
                  nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
                  nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
                  CALL cp_fm_struct_create(fmstruct, para_env=sub_env%para_env, &
                                           context=sub_env%blacs_env, &
                                           nrow_global=nmo_virt(ispin), ncol_global=nmo_occ(ispin))
                  CALL cp_fm_create(evects_mo(ispin, istate), fmstruct)
                  CALL cp_fm_struct_release(fmstruct)
                  CALL cp_dbcsr_sm_fm_multiply(matrix_s, evects(ispin, istate), S_evects(ispin, istate), &
                                               ncol=nmo_occ(ispin), alpha=1.0_dp, beta=0.0_dp)
               END DO
            END DO
            DO istate = 1, nstates
               DO ispin = 1, nspins
                  CALL parallel_gemm("T", "N", &
                                     nmo_virt(ispin), &
                                     nmo_occ(ispin), &
                                     nao, &
                                     1.0_dp, &
                                     gs_mos(ispin)%mos_virt, &
                                     S_evects(ispin, istate), & !this also needs to be orthogonalized
                                     0.0_dp, &
                                     evects_mo(ispin, istate))
               END DO
            END DO
         END IF

         DO istate = 1, nstates
            DO ispin = 1, nspins

               IF (.NOT. print_virtuals) THEN
                  CALL cp_fm_column_scale(evects(ispin, istate), gs_mos(ispin)%phases_occ)
                  IF (ounit > 0) THEN
                     WRITE (ounit, "(/,A)") "ES EIGENVECTORS SIZE"
                     CALL cp_fm_write_info(evects(ispin, istate), ounit)
                  END IF
                  CALL cp_fm_write_formatted(evects(ispin, istate), ounit, "ES EIGENVECTORS")
               ELSE
                  CALL cp_fm_column_scale(evects_mo(ispin, istate), gs_mos(ispin)%phases_occ)
                  IF (ounit > 0) THEN
                     WRITE (ounit, "(/,A)") "ES EIGENVECTORS SIZE"
                     CALL cp_fm_write_info(evects_mo(ispin, istate), ounit)
                  END IF
                  CALL cp_fm_write_formatted(evects_mo(ispin, istate), ounit, "ES EIGENVECTORS")
               END IF

               ! compute and print phase of eigenvectors
               nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
               ALLOCATE (phase_evects(nmo_occ(ispin)))
               IF (print_virtuals) THEN
                  CALL compute_phase_eigenvectors(evects_mo(ispin, istate), phase_evects, sub_env)
               ELSE
                  CALL compute_phase_eigenvectors(evects(ispin, istate), phase_evects, sub_env)
               END IF
               IF (ounit > 0) THEN
                  WRITE (ounit, "(/,A,/)") "PHASES ES EIGENVECTORS"
                  DO iocc = 1, nmo_occ(ispin)
                     WRITE (ounit, "(F20.14)") phase_evects(iocc)
                  END DO
               END IF
               DEALLOCATE (phase_evects)

            END DO
         END DO

         IF (print_virtuals) THEN
            CALL cp_fm_release(evects_mo)
         END IF

         DO ispin = 1, nspins
            IF (ounit > 0) THEN
               WRITE (ounit, "(/,A)") "OCCUPIED MOS SIZE"
               CALL cp_fm_write_info(gs_mos(ispin)%mos_occ, ounit)
            END IF
            CALL cp_fm_write_formatted(gs_mos(ispin)%mos_occ, ounit, "OCCUPIED MO COEFFICIENTS")
         END DO

         IF (ounit > 0) THEN
            WRITE (ounit, "(A)") "OCCUPIED MO EIGENVALUES"
            DO ispin = 1, nspins
               nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
               DO iocc = 1, nmo_occ(ispin)
                  WRITE (ounit, "(F20.14)") gs_mos(ispin)%evals_occ(iocc)
               END DO
            END DO
         END IF
!
         IF (print_virtuals) THEN
            DO ispin = 1, nspins
               IF (ounit > 0) THEN
                  WRITE (ounit, "(/,A)") "VIRTUAL MOS SIZE"
                  CALL cp_fm_write_info(gs_mos(ispin)%mos_virt, ounit)
               END IF
               CALL cp_fm_write_formatted(gs_mos(ispin)%mos_virt, ounit, "VIRTUAL MO COEFFICIENTS")
            END DO

            IF (ounit > 0) THEN
               WRITE (ounit, "(A)") "VIRTUAL MO EIGENVALUES"
               DO ispin = 1, nspins
                  nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
                  DO ivirt = 1, nmo_virt(ispin)
                     WRITE (ounit, "(F20.14)") gs_mos(ispin)%evals_virt(ivirt)
                  END DO
               END DO
            END IF
         END IF

         ! print phases of molecular orbitals

         IF (print_phases) THEN
            IF (ounit > 0) THEN
               WRITE (ounit, "(A)") "PHASES OCCUPIED ORBITALS"
               DO ispin = 1, nspins
                  DO iocc = 1, nmo_occ(ispin)
                     WRITE (ounit, "(F20.14)") gs_mos(ispin)%phases_occ(iocc)
                  END DO
               END DO
               IF (print_virtuals) THEN
                  WRITE (ounit, "(A)") "PHASES VIRTUAL ORBITALS"
                  DO ispin = 1, nspins
                     DO ivirt = 1, nmo_virt(ispin)
                        WRITE (ounit, "(F20.14)") gs_mos(ispin)%phases_virt(ivirt)
                     END DO
                  END DO
               END IF
            END IF
         END IF

         CALL cp_print_key_finished_output(ounit, logger, tddfpt_print_section, "NAMD_PRINT")

         CALL timestop(handle)
      END IF

   END SUBROUTINE tddfpt_write_newtonx_output
! **************************************************************************************************
!> \brief ...
!> \param evects ...
!> \param ounit ...
!> \param S_evects ...
!> \param matrix_s ...
! **************************************************************************************************
   SUBROUTINE tddfpt_check_orthonormality(evects, ounit, S_evects, matrix_s)

      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(in)      :: evects
      INTEGER, INTENT(in)                                :: ounit
      TYPE(cp_fm_type), DIMENSION(:, :), INTENT(in)      :: S_evects
      TYPE(dbcsr_type), INTENT(in), POINTER              :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_check_orthonormality'

      INTEGER                                            :: handle, ispin, ivect, jvect, nspins, &
                                                            nvects_total
      INTEGER, DIMENSION(maxspins)                       :: nactive
      REAL(kind=dp)                                      :: norm
      REAL(kind=dp), DIMENSION(maxspins)                 :: weights

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects_total = SIZE(evects, 2)

      IF (debug_this_module) THEN
         CPASSERT(SIZE(S_evects, 1) == nspins)
         CPASSERT(SIZE(S_evects, 2) == nvects_total)
      END IF

      DO ispin = 1, nspins
         CALL cp_fm_get_info(matrix=evects(ispin, 1), ncol_global=nactive(ispin))
      END DO

      DO jvect = 1, nvects_total
         ! <psi1_i | psi1_j>
         DO ivect = 1, jvect - 1
            CALL cp_fm_trace(evects(:, jvect), S_evects(:, ivect), weights(1:nspins), accurate=.FALSE.)
            norm = SUM(weights(1:nspins))

            DO ispin = 1, nspins
               CALL cp_fm_scale_and_add(1.0_dp, evects(ispin, jvect), -norm, evects(ispin, ivect))
            END DO
         END DO

         ! <psi1_j | psi1_j>
         DO ispin = 1, nspins
            CALL cp_dbcsr_sm_fm_multiply(matrix_s, evects(ispin, jvect), S_evects(ispin, jvect), &
                                         ncol=nactive(ispin), alpha=1.0_dp, beta=0.0_dp)
         END DO

         CALL cp_fm_trace(evects(:, jvect), S_evects(:, jvect), weights(1:nspins), accurate=.FALSE.)

         norm = SUM(weights(1:nspins))
         norm = 1.0_dp/SQRT(norm)

         IF ((ounit > 0) .AND. debug_this_module) WRITE (ounit, '(A,F10.8)') "norm", norm

      END DO

      CALL timestop(handle)

   END SUBROUTINE tddfpt_check_orthonormality
! **************************************************************************************************
!> \brief ...
!> \param evects ...
!> \param phase_evects ...
!> \param sub_env ...
! **************************************************************************************************
   SUBROUTINE compute_phase_eigenvectors(evects, phase_evects, sub_env)

      ! copied from parts of tddgpt_init_ground_state_mos by S. Chulkov

      TYPE(cp_fm_type), INTENT(in)                       :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(out)           :: phase_evects
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_phase_eigenvectors'
      REAL(kind=dp), PARAMETER                           :: eps_dp = EPSILON(0.0_dp)

      INTEGER :: handle, icol_global, icol_local, irow_global, irow_local, ncol_global, &
         ncol_local, nrow_global, nrow_local, sign_int
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: minrow_neg_array, minrow_pos_array, &
                                                            sum_sign_array
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(kind=dp)                                      :: element
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: my_block

      CALL timeset(routineN, handle)

      ! compute and print the phase of excited-state eigenvectors:
      CALL cp_fm_get_info(evects, nrow_global=nrow_global, ncol_global=ncol_global, &
                          nrow_local=nrow_local, ncol_local=ncol_local, local_data=my_block, &
                          row_indices=row_indices, col_indices=col_indices) ! nrow_global either nao or nocc

      ALLOCATE (minrow_neg_array(ncol_global), minrow_pos_array(ncol_global), sum_sign_array(ncol_global))
      minrow_neg_array(:) = nrow_global
      minrow_pos_array(:) = nrow_global
      sum_sign_array(:) = 0

      DO icol_local = 1, ncol_local
         icol_global = col_indices(icol_local)

         DO irow_local = 1, nrow_local
            irow_global = row_indices(irow_local)

            element = my_block(irow_local, icol_local)

            sign_int = 0
            IF (element >= eps_dp) THEN
               sign_int = 1
            ELSE IF (element <= -eps_dp) THEN
               sign_int = -1
            END IF

            sum_sign_array(icol_global) = sum_sign_array(icol_global) + sign_int

            IF (sign_int > 0) THEN
               IF (minrow_pos_array(icol_global) > irow_global) &
                  minrow_pos_array(icol_global) = irow_global
            ELSE IF (sign_int < 0) THEN
               IF (minrow_neg_array(icol_global) > irow_global) &
                  minrow_neg_array(icol_global) = irow_global
            END IF

         END DO
      END DO

      CALL sub_env%para_env%sum(sum_sign_array)
      CALL sub_env%para_env%min(minrow_neg_array)
      CALL sub_env%para_env%min(minrow_pos_array)

      DO icol_global = 1, ncol_global

         IF (sum_sign_array(icol_global) > 0) THEN
            ! most of the expansion coefficients are positive => MO's phase = +1
            phase_evects(icol_global) = 1.0_dp
         ELSE IF (sum_sign_array(icol_global) < 0) THEN
            ! most of the expansion coefficients are negative => MO's phase = -1
            phase_evects(icol_global) = -1.0_dp
         ELSE
            ! equal number of positive and negative expansion coefficients
            IF (minrow_pos_array(icol_global) <= minrow_neg_array(icol_global)) THEN
               ! the first positive expansion coefficient has a lower index then
               ! the first negative expansion coefficient; MO's phase = +1
               phase_evects(icol_global) = 1.0_dp
            ELSE
               ! MO's phase = -1
               phase_evects(icol_global) = -1.0_dp
            END IF
         END IF

      END DO

      DEALLOCATE (minrow_neg_array, minrow_pos_array, sum_sign_array)

      CALL timestop(handle)

   END SUBROUTINE compute_phase_eigenvectors

END MODULE qs_tddfpt2_restart
