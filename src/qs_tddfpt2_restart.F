!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_restart
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type,&
                                              get_blacs_info
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pool_create_fm
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_read_unformatted,&
                                              cp_fm_write_unformatted
   USE cp_log_handling,                 ONLY: cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE message_passing,                 ONLY: mp_bcast
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

   PUBLIC :: tddfpt_write_restart, tddfpt_read_restart

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
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
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
               CALL cp_fm_column_scale(evects(ispin, istate)%matrix, gs_mos(ispin)%phases_occ)

               CALL cp_fm_write_unformatted(evects(ispin, istate)%matrix, ounit)

               CALL cp_fm_column_scale(evects(ispin, istate)%matrix, gs_mos(ispin)%phases_occ)
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
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
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
      TYPE(cp_para_env_type), POINTER                    :: para_env_global
      TYPE(section_vals_type), POINTER                   :: print_key

      CALL timeset(routineN, handle)

      ! generate restart file name
      CALL section_vals_val_get(tddfpt_section, "WFN_RESTART_FILE_NAME", n_rep_val=n_rep_val)
      IF (n_rep_val > 0) THEN
         CALL section_vals_val_get(tddfpt_section, "WFN_RESTART_FILE_NAME", c_val=filename)
      ELSE
         print_key => section_vals_get_subs_vals(tddfpt_print_section, "RESTART")
         filename = cp_print_key_generate_filename(logger, print_key, &
                                                   extension=".tdwfn", my_local=.FALSE.)
      END IF

      CALL get_blacs_info(blacs_env_global, para_env=para_env_global)

      IF (para_env_global%ionode) THEN
         INQUIRE (FILE=filename, exist=file_exists)

         IF (.NOT. file_exists) THEN
            nstates_read = 0
            CALL mp_bcast(nstates_read, para_env_global%source, para_env_global%group)

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

      IF (para_env_global%ionode) THEN
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
      CALL mp_bcast(nstates_read, para_env_global%source, para_env_global%group)

      ! exit if restart file does not exist
      IF (nstates_read <= 0) THEN
         CALL timestop(handle)
         RETURN
      END IF

      IF (para_env_global%ionode) THEN
         ALLOCATE (evals_read(nstates_read))
         READ (iunit) evals_read
         IF (nstates_read <= nstates) THEN
            evals(1:nstates_read) = evals_read(1:nstates_read)
         ELSE
            evals(1:nstates) = evals_read(1:nstates)
         END IF
         DEALLOCATE (evals_read)
      END IF
      CALL mp_bcast(evals, para_env_global%source, para_env_global%group)

      DO istate = 1, nstates_read
         DO ispin = 1, nspins
            IF (istate <= nstates) THEN
               CALL fm_pool_create_fm(fm_pool_ao_mo_occ(ispin)%pool, evects(ispin, istate)%matrix)

               CALL cp_fm_read_unformatted(evects(ispin, istate)%matrix, iunit)

               CALL cp_fm_column_scale(evects(ispin, istate)%matrix, gs_mos(ispin)%phases_occ)
            END IF
         END DO
      END DO

      IF (para_env_global%ionode) &
         CALL close_file(unit_number=iunit)

      CALL timestop(handle)

   END FUNCTION tddfpt_read_restart

END MODULE qs_tddfpt2_restart
