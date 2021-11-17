!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate the derivatives of the MO coefficients wrt nuclear coordinates
!> \author Sandra Luber, Edward Ditler
! **************************************************************************************************

MODULE qs_dcdr_utils
!#include "./common/cp_common_uses.f90"
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_p_type, cp_fm_release, &
        cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_result_methods,               ONLY: get_results
   USE cp_result_types,                 ONLY: cp_result_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_init_p, dbcsr_p_type, dbcsr_set, &
        dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE input_constants,                 ONLY: current_orb_center_wannier,&
                                              use_mom_ref_user
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_bcast
   USE molecule_types,                  ONLY: molecule_type
   USE moments_utils,                   ONLY: get_reference_point
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_linres_types,                 ONLY: dcdr_env_type,&
                                              linres_control_type
   USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                              localized_wfn_control_type,&
                                              qs_loc_env_new_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_matrix
   USE string_utilities,                ONLY: xstring
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: dcdr_env_cleanup, dcdr_env_init, dcdr_print, &
             get_loc_setting, shift_wannier_into_cell, &
             dcdr_write_restart, dcdr_read_restart, &
             multiply_localization

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dcdr_utils'

   REAL(dp), DIMENSION(3, 3, 3), PARAMETER  :: Levi_Civita = RESHAPE((/ &
                                                          0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                                                          0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, &
                                                        0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/), &
                                                                     (/3, 3, 3/))

CONTAINS

! **************************************************************************************************
!> \brief Multiply (ao_matrix @ mo_coeff) and store the column icenter in res
!> \param ao_matrix ...
!> \param mo_coeff ...
!> \param work Working space
!> \param nmo ...
!> \param icenter ...
!> \param res ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE multiply_localization(ao_matrix, mo_coeff, work, nmo, icenter, res)
      TYPE(dbcsr_type), INTENT(IN), POINTER              :: ao_matrix
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: mo_coeff, work
      INTEGER, INTENT(IN)                                :: nmo, icenter
      TYPE(cp_fm_type), INTENT(IN), POINTER              :: res

      CHARACTER(LEN=*), PARAMETER :: routineN = 'multiply_localization', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! Multiply by the MO coefficients
      CALL cp_dbcsr_sm_fm_multiply(ao_matrix, mo_coeff, work, ncol=nmo)

      ! Only keep the icenter-th column
      CALL cp_fm_to_fm(work, res, 1, icenter, icenter)

      ! Reset the matrices
      CALL cp_fm_set_all(work, 0.0_dp)

      CALL timestop(handle)
   END SUBROUTINE multiply_localization

! **************************************************************************************************
!> \brief Copied from linres_read_restart
!> \param qs_env ...
!> \param linres_section ...
!> \param vec ...
!> \param lambda ...
!> \param beta ...
!> \param tag ...
!> \note Adapted from linres_read_restart (ED)
!>       Would be nice not to crash but to start from zero if the present file doesn't match.
! **************************************************************************************************
   SUBROUTINE dcdr_read_restart(qs_env, linres_section, vec, lambda, beta, tag)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: linres_section
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: vec
      INTEGER, INTENT(IN)                                :: lambda, beta
      CHARACTER(LEN=*)                                   :: tag

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'dcdr_read_restart'

      CHARACTER(LEN=default_path_length)                 :: filename
      CHARACTER(LEN=default_string_length)               :: my_middle
      INTEGER :: beta_tmp, group, handle, i, i_block, ia, ie, iostat, iounit, ispin, j, &
         lambda_tmp, max_block, n_rep_val, nao, nao_tmp, nmo, nmo_tmp, nspins, nspins_tmp, &
         rst_unit, source
      LOGICAL                                            :: file_exists
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vecbuffer
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(section_vals_type), POINTER                   :: print_key

      file_exists = .FALSE.

      CALL timeset(routineN, handle)

      NULLIFY (mos, para_env, logger, print_key, vecbuffer)
      logger => cp_get_default_logger()

      iounit = cp_print_key_unit_nr(logger, linres_section, &
                                    "PRINT%PROGRAM_RUN_INFO", extension=".Log")

      CALL get_qs_env(qs_env=qs_env, &
                      para_env=para_env, &
                      mos=mos)

      nspins = SIZE(mos)
      group = para_env%group
      source = para_env%source !ionode???

      rst_unit = -1
      IF (para_env%ionode) THEN
         CALL section_vals_val_get(linres_section, "WFN_RESTART_FILE_NAME", &
                                   n_rep_val=n_rep_val)

         CALL XSTRING(tag, ia, ie)
         my_middle = "RESTART-"//tag(ia:ie)//TRIM("-")//TRIM(ADJUSTL(cp_to_string(beta))) &
                     //TRIM("-")//TRIM(ADJUSTL(cp_to_string(lambda)))

         IF (n_rep_val > 0) THEN
            CALL section_vals_val_get(linres_section, "WFN_RESTART_FILE_NAME", c_val=filename)
            CALL xstring(filename, ia, ie)
            filename = filename(ia:ie)//TRIM(my_middle)//".lr"
         ELSE
            ! try to read from the filename that is generated automatically from the printkey
            print_key => section_vals_get_subs_vals(linres_section, "PRINT%RESTART")
            filename = cp_print_key_generate_filename(logger, print_key, &
                                                      extension=".lr", middle_name=TRIM(my_middle), my_local=.FALSE.)
         END IF
         INQUIRE (FILE=filename, exist=file_exists)
         !
         ! open file
         IF (file_exists) THEN
            CALL open_file(file_name=TRIM(filename), &
                           file_action="READ", &
                           file_form="UNFORMATTED", &
                           file_position="REWIND", &
                           file_status="OLD", &
                           unit_number=rst_unit)

            IF (iounit > 0) WRITE (iounit, "(T2,A)") &
               "LINRES| Reading response wavefunctions from the restart file <"//TRIM(ADJUSTL(filename))//">"
         ELSE
            IF (iounit > 0) WRITE (iounit, "(T2,A)") &
               "LINRES| Restart file  <"//TRIM(ADJUSTL(filename))//"> not found"
         END IF
      END IF

      CALL mp_bcast(file_exists, source, group)

      IF (file_exists) THEN

         CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, nrow_global=nao, ncol_block=max_block)

         ALLOCATE (vecbuffer(nao, max_block))
         !
         ! read headers
         IF (rst_unit > 0) READ (rst_unit, IOSTAT=iostat) lambda_tmp, beta_tmp, nspins_tmp, nao_tmp
         CALL mp_bcast(iostat, source, group)

         CALL mp_bcast(beta_tmp, source, group)
         CALL mp_bcast(lambda_tmp, source, group)
         CALL mp_bcast(nspins_tmp, source, group)
         CALL mp_bcast(nao_tmp, source, group)

         ! check that the number nao, nmo and nspins are
         ! the same as in the current mos
         IF (nspins_tmp .NE. nspins) THEN
            CPABORT("nspins not consistent")
         END IF
         IF (nao_tmp .NE. nao) CPABORT("nao not consistent")
         ! check that it's the right file
         ! the same as in the current mos
         IF (lambda_tmp .NE. lambda) CPABORT("lambda not consistent")
         IF (beta_tmp .NE. beta) CPABORT("beta not consistent")
         !
         DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin)%mo_set, mo_coeff=mo_coeff)
            CALL cp_fm_get_info(mo_coeff, ncol_global=nmo)
            !
            IF (rst_unit > 0) READ (rst_unit) nmo_tmp
            CALL mp_bcast(nmo_tmp, source, group)
            IF (nmo_tmp .NE. nmo) CPABORT("nmo not consistent")
            !
            ! read the response
            DO i = 1, nmo, MAX(max_block, 1)
               i_block = MIN(max_block, nmo - i + 1)
               DO j = 1, i_block
                  IF (rst_unit > 0) READ (rst_unit) vecbuffer(1:nao, j)
               END DO
               CALL mp_bcast(vecbuffer, source, group)
               CALL cp_fm_set_submatrix(vec(ispin)%matrix, vecbuffer, 1, i, nao, i_block)
            END DO
         END DO

         IF (iostat /= 0) THEN
            IF (iounit > 0) WRITE (iounit, "(T2,A)") &
               "LINRES| Restart file <"//TRIM(ADJUSTL(filename))//"> not found"
         END IF

         DEALLOCATE (vecbuffer)

      END IF

      IF (para_env%ionode) THEN
         IF (file_exists) CALL close_file(unit_number=rst_unit)
      END IF

      CALL timestop(handle)

   END SUBROUTINE dcdr_read_restart

! **************************************************************************************************
!> \brief Copied from linres_write_restart
!> \param qs_env ...
!> \param linres_section ...
!> \param vec ...
!> \param lambda ...
!> \param beta ...
!> \param tag ...
!> \note Adapted from linres_read_restart (ED)
!>       Would be nice not to crash but to start from zero if the present file doesn't match.
! **************************************************************************************************
   SUBROUTINE dcdr_write_restart(qs_env, linres_section, vec, lambda, beta, tag)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: linres_section
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: vec
      INTEGER, INTENT(IN)                                :: lambda, beta
      CHARACTER(LEN=*)                                   :: tag

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dcdr_write_restart'

      CHARACTER(LEN=default_path_length)                 :: filename
      CHARACTER(LEN=default_string_length)               :: my_middle, my_pos, my_status
      INTEGER                                            :: handle, i, i_block, ia, ie, iounit, &
                                                            ispin, j, max_block, nao, nmo, nspins, &
                                                            rst_unit
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: vecbuffer
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(section_vals_type), POINTER                   :: print_key

      NULLIFY (logger, mo_coeff, mos, para_env, print_key, vecbuffer)

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()

      IF (BTEST(cp_print_key_should_output(logger%iter_info, linres_section, "PRINT%RESTART", &
                                           used_print_key=print_key), &
                cp_p_file)) THEN

         iounit = cp_print_key_unit_nr(logger, linres_section, &
                                       "PRINT%PROGRAM_RUN_INFO", extension=".Log")

         CALL get_qs_env(qs_env=qs_env, &
                         mos=mos, &
                         para_env=para_env)

         nspins = SIZE(mos)

         my_status = "REPLACE"
         my_pos = "REWIND"
         CALL XSTRING(tag, ia, ie)
         my_middle = "RESTART-"//tag(ia:ie)//TRIM("-")//TRIM(ADJUSTL(cp_to_string(beta))) &
                     //TRIM("-")//TRIM(ADJUSTL(cp_to_string(lambda)))
         rst_unit = cp_print_key_unit_nr(logger, linres_section, "PRINT%RESTART", &
                                         extension=".lr", middle_name=TRIM(my_middle), file_status=TRIM(my_status), &
                                         file_position=TRIM(my_pos), file_action="WRITE", file_form="UNFORMATTED")

         filename = cp_print_key_generate_filename(logger, print_key, &
                                                   extension=".lr", middle_name=TRIM(my_middle), my_local=.FALSE.)

         IF (iounit > 0) THEN
            WRITE (UNIT=iounit, FMT="(T2,A)") &
               "LINRES| Writing response functions to the restart file <"//TRIM(ADJUSTL(filename))//">"
         END IF

         !
         ! write data to file
         ! use the scalapack block size as a default for buffering columns
         CALL get_mo_set(mos(1)%mo_set, mo_coeff=mo_coeff)
         CALL cp_fm_get_info(mo_coeff, nrow_global=nao, ncol_block=max_block)
         ALLOCATE (vecbuffer(nao, max_block))

         IF (rst_unit > 0) WRITE (rst_unit) lambda, beta, nspins, nao

         DO ispin = 1, nspins
            CALL cp_fm_get_info(vec(ispin)%matrix, ncol_global=nmo)

            IF (rst_unit > 0) WRITE (rst_unit) nmo

            DO i = 1, nmo, MAX(max_block, 1)
               i_block = MIN(max_block, nmo - i + 1)
               CALL cp_fm_get_submatrix(vec(ispin)%matrix, vecbuffer, 1, i, nao, i_block)
               ! doing this in one write would increase efficiency, but breaks RESTART compatibility.
               ! to old ones, and in cases where max_block is different between runs, as might happen during
               ! restarts with a different number of CPUs
               DO j = 1, i_block
                  IF (rst_unit > 0) WRITE (rst_unit) vecbuffer(1:nao, j)
               END DO
            END DO
         END DO

         DEALLOCATE (vecbuffer)

         CALL cp_print_key_finished_output(rst_unit, logger, linres_section, &
                                           "PRINT%RESTART")
      END IF

      CALL timestop(handle)

   END SUBROUTINE dcdr_write_restart

! **************************************************************************************************
!> \brief Print the APT and sum rules
!> \param dcdr_env ...
!> \param qs_env ...
!> \author Edward Ditler
! **************************************************************************************************
   SUBROUTINE dcdr_print(dcdr_env, qs_env)
      TYPE(dcdr_env_type)                                :: dcdr_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'dcdr_print', routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_string_length)               :: description
      INTEGER                                            :: alpha, beta, delta, gamma, i, k, l, &
                                                            lambda, natom, nsubset, output_unit
      REAL(dp), DIMENSION(:, :, :), POINTER              :: apt_el_dcdr, apt_nuc_dcdr, apt_total_dcdr
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: apt_center_dcdr, apt_subset_dcdr
      REAL(kind=dp), DIMENSION(3, 3)                     :: sum_rule_0, sum_rule_1, sum_rule_2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_result_type), POINTER                      :: results
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: dcdr_section

      NULLIFY (logger)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      dcdr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES%DCDR")

      NULLIFY (particle_set)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, molecule_set=molecule_set)
      natom = SIZE(particle_set)
      nsubset = SIZE(molecule_set)

      apt_el_dcdr => dcdr_env%apt_el_dcdr
      apt_nuc_dcdr => dcdr_env%apt_nuc_dcdr
      apt_total_dcdr => dcdr_env%apt_total_dcdr
      apt_subset_dcdr => dcdr_env%apt_el_dcdr_per_subset
      apt_center_dcdr => dcdr_env%apt_el_dcdr_per_center

      IF (dcdr_env%localized_psi0) THEN
         IF (dcdr_env%output_unit > 0) WRITE (dcdr_env%output_unit, "(A)") 'APT | Write the final apt matrix per atom per subset'
         DO k = 1, natom
            DO l = 1, nsubset
               IF (dcdr_env%output_unit > 0) WRITE (dcdr_env%output_unit, *) 'APT | Subset', l
               DO i = 1, 3
                  IF (dcdr_env%output_unit > 0) WRITE (dcdr_env%output_unit, "(A,I3,F15.6,F15.6,F15.6)") &
                     'APT | apt_subset ', i, apt_subset_dcdr(i, :, k, l)
               END DO
            END DO
         END DO
      END IF

      IF (dcdr_env%output_unit > 0) WRITE (dcdr_env%output_unit, "(A)") &
         'APT | Write the final apt matrix per atom (Position perturbation)'
      DO l = 1, natom
         IF (dcdr_env%output_unit > 0) WRITE (dcdr_env%output_unit, "(A,I3,A,F15.6)") &
            'APT | Atom', l, ' - GAPT ', &
            (apt_total_dcdr(1, 1, l) &
             + apt_total_dcdr(2, 2, l) &
             + apt_total_dcdr(3, 3, l))/3._dp
         DO i = 1, 3
            IF (dcdr_env%output_unit > 0) WRITE (dcdr_env%output_unit, "(A,F15.6,F15.6,F15.6)") "APT | ", apt_total_dcdr(i, :, l)
         END DO
      END DO

      IF (dcdr_env%output_unit > 0) WRITE (dcdr_env%output_unit, "(A)") 'APT | Write the total apt matrix'
      DO i = 1, 3
         IF (dcdr_env%output_unit > 0) WRITE (dcdr_env%output_unit, &
                                              "(A,F15.6,F15.6,F15.6)") "APT | ", SUM(apt_total_dcdr(i, :, :), dim=2)
      END DO
      IF (dcdr_env%output_unit > 0) WRITE (dcdr_env%output_unit, "(A)") 'APT | End Write the final apt matrix'

      ! Get the dipole
      CALL get_qs_env(qs_env, results=results)
      description = "[DIPOLE]"
      CALL get_results(results=results, description=description, values=dcdr_env%dipole_pos(1:3))

      ! Sum rules [for all alpha, beta]
      sum_rule_0 = 0._dp
      sum_rule_1 = 0._dp
      sum_rule_2 = 0._dp

      DO alpha = 1, 3
         DO beta = 1, 3
            ! 0: sum_lambda apt(alpha, beta, lambda)
            DO lambda = 1, natom
               sum_rule_0(alpha, beta) = sum_rule_0(alpha, beta) &
                                         + apt_total_dcdr(alpha, beta, lambda)
            END DO

            ! 1: sum_gamma epsilon_(alpha beta gamma) mu_gamma
            DO gamma = 1, 3
               sum_rule_1(alpha, beta) = sum_rule_1(alpha, beta) &
                                         + Levi_Civita(alpha, beta, gamma)*dcdr_env%dipole_pos(gamma)
            END DO

            ! 2: sum_(lambda gamma delta) R^lambda_gamma apt(delta, alpha, lambda)
            DO lambda = 1, natom
               DO gamma = 1, 3
                  DO delta = 1, 3
                     sum_rule_2(alpha, beta) = sum_rule_2(alpha, beta) &
                                               + Levi_Civita(beta, gamma, delta) &
                                               *particle_set(lambda)%r(gamma) &
                                               *apt_total_dcdr(delta, alpha, lambda)
                  END DO
               END DO
            END DO

         END DO ! beta
      END DO   ! alpha

      IF (dcdr_env%output_unit > 0) WRITE (dcdr_env%output_unit, "(A)") "APT | Position perturbation sum rules"
      IF (dcdr_env%output_unit > 0) WRITE (dcdr_env%output_unit, "(A,T18,A,T34,A,T49,A)") &
         "APT |", "Total APT", "Dipole", "R * APT"
      DO alpha = 1, 3
         DO beta = 1, 3
            IF (dcdr_env%output_unit > 0) WRITE (dcdr_env%output_unit, &
                                                 "(A,I3,I3,F15.6,F15.6,F15.6)") &
               "APT | ", &
               alpha, beta, &
               sum_rule_0(alpha, beta), &
               sum_rule_1(alpha, beta), &
               sum_rule_2(alpha, beta)
         END DO
      END DO

   END SUBROUTINE dcdr_print

! **************************************************************************************************
!> \brief ...
!> \param r ...
!> \param cell ...
!> \param r_shifted ...
! **************************************************************************************************
   SUBROUTINE shift_wannier_into_cell(r, cell, r_shifted)
      REAL(dp), DIMENSION(3), INTENT(in)                 :: r
      TYPE(cell_type), INTENT(in), POINTER               :: cell
      REAL(dp), DIMENSION(3), INTENT(out)                :: r_shifted

      INTEGER                                            :: i
      REAL(kind=dp), DIMENSION(3)                        :: abc

      ! Only orthorombic cell for now
      CALL get_cell(cell, abc=abc)

      DO i = 1, 3
         IF (r(i) < 0._dp) THEN
            r_shifted(i) = r(i) + abc(i)
         ELSE IF (r(i) > abc(i)) THEN
            r_shifted(i) = r(i) - abc(i)
         ELSE
            r_shifted(i) = r(i)
         END IF
      END DO
   END SUBROUTINE shift_wannier_into_cell

! **************************************************************************************************
!> \brief ...
!> \param dcdr_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_loc_setting(dcdr_env, qs_env)
      TYPE(dcdr_env_type)                                :: dcdr_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_loc_setting', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, is, istate, max_states, nmo, &
                                                            nmoloc, nstate, nstate_list(2)
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: state_list
      REAL(dp), DIMENSION(:, :), POINTER                 :: center_array
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(localized_wfn_control_type), POINTER          :: localized_wfn_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_loc_env_new_type), POINTER                 :: qs_loc_env
      TYPE(section_vals_type), POINTER                   :: dcdr_section

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, &
                      linres_control=linres_control, &
                      mos=mos)

      ! Some checks
      max_states = 0
      CALL get_mo_set(mo_set=mos(1)%mo_set, nmo=nmo)
      max_states = MAX(max_states, nmo)

      ! check that the number of localized states is equal to the number of states
      nmoloc = SIZE(linres_control%qs_loc_env%localized_wfn_control%centers_set(1)%array, 2)
      IF (nmoloc .NE. nmo) THEN
         CPABORT("The number of localized functions is not equal to the number of states.")
      END IF

      ! which center for the orbitals shall we use
      dcdr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES%DCDR")
      CALL section_vals_val_get(dcdr_section, "ORBITAL_CENTER", i_val=dcdr_env%orb_center)
      SELECT CASE (dcdr_env%orb_center)
      CASE (current_orb_center_wannier)
         dcdr_env%orb_center_name = "WANNIER"
      CASE DEFAULT
         CPABORT(" ")
      END SELECT

      qs_loc_env => linres_control%qs_loc_env
      CALL get_qs_loc_env(qs_loc_env, localized_wfn_control=localized_wfn_control)

      ALLOCATE (dcdr_env%centers_set(1))
      ALLOCATE (dcdr_env%center_list(1))
      ALLOCATE (state_list(max_states, 1))
      state_list(:, :) = HUGE(0)
      nstate_list(:) = HUGE(0)

      ! Build the state_list
      center_array => localized_wfn_control%centers_set(1)%array
      nstate = 0
      DO istate = 1, SIZE(center_array, 2)
         nstate = nstate + 1
         state_list(nstate, 1) = istate
      END DO
      nstate_list(1) = nstate

      ! clustering the states
      nstate = nstate_list(1)
      dcdr_env%nstates(1) = nstate

      ALLOCATE (dcdr_env%center_list(1)%array(2, nstate + 1))
      ALLOCATE (dcdr_env%centers_set(1)%array(3, nstate))
      dcdr_env%center_list(1)%array(:, :) = HUGE(0)
      dcdr_env%centers_set(1)%array(:, :) = HUGE(0.0_dp)

      center_array => localized_wfn_control%centers_set(1)%array

      ! point to the psi0 centers
      SELECT CASE (dcdr_env%orb_center)
      CASE (current_orb_center_wannier)
         ! use the wannier center as -center-
         dcdr_env%nbr_center(1) = nstate
         DO is = 1, nstate
            istate = state_list(is, 1)
            dcdr_env%centers_set(1)%array(1:3, is) = center_array(1:3, istate)
            dcdr_env%center_list(1)%array(1, is) = is
            dcdr_env%center_list(1)%array(2, is) = istate
         END DO
         dcdr_env%center_list(1)%array(1, nstate + 1) = nstate + 1

      CASE DEFAULT
         CPABORT("Unknown orbital center...")
      END SELECT

      ! This doesn't do anything right?
      ! ! hack to make centers_set%array closest to chosen reference
      ! nsubset = SIZE(molecule_set)
      ! natoms = SIZE(particle_set)
      ! ALLOCATE (map_rc(3, natoms))
      ! ALLOCATE (mapm_rc(3, nsubset))

      ! ! Get reference point for this molecule
      ! NULLIFY (ref_point)
      ! reference = section_get_ival(section_vals=qs_env%input, &
      !                              keyword_name="DFT%PRINT%MOMENTS%REFERENCE")
      ! CALL section_vals_val_get(qs_env%input, "DFT%PRINT%MOMENTS%REF_POINT", r_vals=ref_point)

      ! DO i = 1, nsubset
      !    first_atom = molecule_set(i)%first_atom
      !    last_atom = molecule_set(i)%last_atom

      !    CALL get_reference_point(rcc, qs_env=qs_env, reference=reference, &
      !                             ref_point=ref_point, &
      !                             ifirst=first_atom, ilast=last_atom)
      !    mapm_rc(:, i) = rcc
      !    DO j = first_atom, last_atom
      !       map_rc(:, j) = rcc
      !    END DO
      ! END DO

      ! DEALLOCATE (state_list)
      ! DEALLOCATE (map_rc)
      ! DEALLOCATE (mapm_rc)

      CALL timestop(handle)
   END SUBROUTINE get_loc_setting

! **************************************************************************************************
!> \brief Initialize the dcdr environment
!> \param dcdr_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE dcdr_env_init(dcdr_env, qs_env)
      TYPE(dcdr_env_type)                                :: dcdr_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dcdr_env_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, homo, i, isize, j, jg, n_rep, &
                                                            nao, natom, nmo, nspins, nsubset, &
                                                            output_unit, reference, unit_number
      INTEGER, DIMENSION(:), POINTER                     :: tmplist
      LOGICAL                                            :: explicit
      REAL(KIND=dp), DIMENSION(:), POINTER               :: ref_point
      TYPE(cp_fm_type), POINTER                          :: buf, mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_all, sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(section_vals_type), POINTER                   :: dcdr_section, loc_section, lr_section

      CALL timeset(routineN, handle)

      ! Set up the logger
      NULLIFY (logger, loc_section, dcdr_section, lr_section)
      logger => cp_get_default_logger()
      loc_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES%LOCALIZE")
      dcdr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES%DCDR")
      dcdr_env%output_unit = cp_print_key_unit_nr(logger, dcdr_section, "PRINT%APT", &
                                                  extension=".data", middle_name="dcdr", log_filename=.FALSE., &
                                                  file_position="REWIND", file_status="REPLACE")

      lr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES")
      output_unit = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", &
                                         extension=".linresLog")
      unit_number = cp_print_key_unit_nr(logger, lr_section, "PRINT%PROGRAM_RUN_INFO", extension=".linresLog")

      IF (output_unit > 0) THEN
         WRITE (output_unit, "(/,T20,A,/)") "*** Start DCDR calculation ***"
      END IF

      NULLIFY (ks_env, dft_control, sab_orb, sab_all, particle_set, molecule_set, matrix_s, matrix_ks, mos, para_env)
      CALL get_qs_env(qs_env=qs_env, &
                      ks_env=ks_env, &
                      dft_control=dft_control, &
                      sab_orb=sab_orb, &
                      sab_all=sab_all, &
                      particle_set=particle_set, &
                      molecule_set=molecule_set, &
                      matrix_s=matrix_s, &
                      matrix_ks=matrix_ks, &
                      mos=mos, &
                      para_env=para_env)

      natom = SIZE(particle_set)
      nsubset = SIZE(molecule_set)
      nspins = dft_control%nspins
      IF (nspins == 2) THEN
         CPABORT("NYE")
      END IF

      CALL build_overlap_matrix(ks_env, matrix_s=dcdr_env%matrix_s, &
                                matrix_name="OVERLAP MATRIX", &
                                nderivative=1, &
                                basis_type_a="ORB", &
                                basis_type_b="ORB", &
                                sab_nl=sab_orb)

      CALL build_kinetic_matrix(ks_env, matrix_t=dcdr_env%matrix_t, &
                                matrix_name="KINETIC ENERGY MATRIX", &
                                basis_type="ORB", &
                                sab_nl=sab_orb, nderivative=1, &
                                eps_filter=dft_control%qs_control%eps_filter_matrix)

      ! Get inputs
      CALL section_vals_val_get(dcdr_section, "DISTRIBUTED_ORIGIN", l_val=dcdr_env%distributed_origin)
      CALL section_vals_val_get(loc_section, "_SECTION_PARAMETERS_", l_val=dcdr_env%localized_psi0)
      CALL section_vals_val_get(lr_section, "RESTART", l_val=qs_env%linres_control%linres_restart)
      dcdr_env%ref_point = 0._dp

      ! List of atoms
      NULLIFY (tmplist)
      isize = 0
      CALL section_vals_val_get(dcdr_section, "LIST_OF_ATOMS", n_rep_val=n_rep)
      IF (n_rep == 0) THEN
         ALLOCATE (dcdr_env%list_of_atoms(natom))
         DO jg = 1, natom
            dcdr_env%list_of_atoms(jg) = jg
         END DO
      ELSE
         DO jg = 1, n_rep
            ALLOCATE (dcdr_env%list_of_atoms(isize))
            CALL section_vals_val_get(dcdr_section, "LIST_OF_ATOMS", i_rep_val=jg, i_vals=tmplist)
            CALL reallocate(dcdr_env%list_of_atoms, 1, isize + SIZE(tmplist))
            dcdr_env%list_of_atoms(isize + 1:isize + SIZE(tmplist)) = tmplist
            isize = SIZE(dcdr_env%list_of_atoms)
         END DO
      END IF

      ! Reference point
      IF (dcdr_env%localized_psi0) THEN
         ! Get the Wannier localized wave functions and centers
         CALL get_loc_setting(dcdr_env, qs_env)
      ELSE
         ! Get the reference point from the input
         CALL section_vals_val_get(dcdr_section, "REFERENCE", i_val=reference)
         CALL section_vals_val_get(dcdr_section, "REFERENCE_POINT", explicit=explicit)
         IF (explicit) THEN
            CALL section_vals_val_get(dcdr_section, "REFERENCE_POINT", r_vals=ref_point)
         ELSE
            IF (reference == use_mom_ref_user) &
               CPABORT("User-defined reference point should be given explicitly")
         END IF

         CALL get_reference_point(rpoint=dcdr_env%ref_point, qs_env=qs_env, &
                                  reference=reference, &
                                  ref_point=ref_point)
      END IF

      ! Helper matrix structs
      NULLIFY (dcdr_env%aoao_fm_struct, &
               dcdr_env%momo_fm_struct, &
               dcdr_env%likemos_fm_struct, &
               dcdr_env%homohomo_fm_struct)
      CALL get_mo_set(mo_set=mos(1)%mo_set, mo_coeff=mo_coeff, &
                      nao=nao, nmo=nmo, homo=homo)
      CALL cp_fm_struct_create(dcdr_env%aoao_fm_struct, nrow_global=nao, &
                               ncol_global=nao, para_env=para_env, &
                               context=mo_coeff%matrix_struct%context)
      CALL cp_fm_struct_create(dcdr_env%momo_fm_struct, nrow_global=nmo, &
                               ncol_global=nmo, para_env=para_env, &
                               context=mo_coeff%matrix_struct%context)
      CALL cp_fm_struct_create(dcdr_env%homohomo_fm_struct, nrow_global=homo, &
                               ncol_global=homo, para_env=para_env, &
                               context=mo_coeff%matrix_struct%context)
      CALL cp_fm_struct_create(dcdr_env%likemos_fm_struct, &
                               template_fmstruct=mo_coeff%matrix_struct)
      dcdr_env%nao = nao
      dcdr_env%nmo = nmo

      ! Fields of reals
      ALLOCATE (dcdr_env%deltaR(3, natom))
      ALLOCATE (dcdr_env%delta_basis_function(3, natom))
      ALLOCATE (dcdr_env%apt_el_dcdr(3, 3, natom))
      ALLOCATE (dcdr_env%apt_nuc_dcdr(3, 3, natom))
      ALLOCATE (dcdr_env%apt_total_dcdr(3, 3, natom))

      dcdr_env%apt_el_dcdr = 0._dp
      dcdr_env%apt_nuc_dcdr = 0._dp
      dcdr_env%apt_total_dcdr = 0._dp

      dcdr_env%deltaR = 0.0_dp
      dcdr_env%delta_basis_function = 0._dp

      ! Localization
      IF (dcdr_env%localized_psi0) THEN
         ALLOCATE (dcdr_env%apt_el_dcdr_per_center(3, 3, natom, dcdr_env%nbr_center(1)))
         ALLOCATE (dcdr_env%apt_el_dcdr_per_subset(3, 3, natom, nsubset))
         ALLOCATE (dcdr_env%apt_subset(3, 3, natom, nsubset))
         dcdr_env%apt_el_dcdr_per_center = 0._dp
         dcdr_env%apt_el_dcdr_per_subset = 0._dp
         dcdr_env%apt_subset = 0.0_dp
      END IF

      ! Full matrices
      ALLOCATE (dcdr_env%mo_coeff(nspins))
      ALLOCATE (dcdr_env%dCR(nspins))
      ALLOCATE (dcdr_env%dCR_prime(nspins))
      ALLOCATE (dcdr_env%chc(nspins))
      ALLOCATE (dcdr_env%op_dR(nspins))

      NULLIFY (dcdr_env%mo_coeff(1)%matrix)
      NULLIFY (dcdr_env%dCR(1)%matrix)
      NULLIFY (dcdr_env%dCR_prime(1)%matrix)
      NULLIFY (dcdr_env%chc(1)%matrix)
      NULLIFY (dcdr_env%op_dR(1)%matrix)

      CALL cp_fm_create(dcdr_env%dCR(1)%matrix, dcdr_env%likemos_fm_struct)
      CALL cp_fm_create(dcdr_env%dCR_prime(1)%matrix, dcdr_env%likemos_fm_struct)
      CALL cp_fm_create(dcdr_env%mo_coeff(1)%matrix, dcdr_env%likemos_fm_struct)
      CALL cp_fm_create(dcdr_env%chc(1)%matrix, dcdr_env%momo_fm_struct)
      CALL cp_fm_create(dcdr_env%op_dR(1)%matrix, dcdr_env%likemos_fm_struct)

      CALL cp_fm_to_fm(mo_coeff, dcdr_env%mo_coeff(1)%matrix)

      ! DBCSR matrices
      NULLIFY (dcdr_env%hamiltonian1)
      NULLIFY (dcdr_env%moments)
      NULLIFY (dcdr_env%matrix_difdip)
      NULLIFY (dcdr_env%matrix_core_charge_1)
      NULLIFY (dcdr_env%matrix_nosym_temp)
      NULLIFY (dcdr_env%matrix_s1)
      NULLIFY (dcdr_env%matrix_t1)
      NULLIFY (dcdr_env%matrix_apply_op_constant)
      NULLIFY (dcdr_env%matrix_d_vhxc_dR)
      NULLIFY (dcdr_env%matrix_vhxc_perturbed_basis)
      NULLIFY (dcdr_env%matrix_hc)
      NULLIFY (dcdr_env%matrix_ppnl_1)
      CALL dbcsr_allocate_matrix_set(dcdr_env%hamiltonian1, nspins)
      CALL dbcsr_allocate_matrix_set(dcdr_env%moments, 3)
      CALL dbcsr_allocate_matrix_set(dcdr_env%matrix_difdip, 3, 3)
      CALL dbcsr_allocate_matrix_set(dcdr_env%matrix_core_charge_1, 3)
      CALL dbcsr_allocate_matrix_set(dcdr_env%matrix_nosym_temp, 3)
      CALL dbcsr_allocate_matrix_set(dcdr_env%matrix_s1, 4)
      CALL dbcsr_allocate_matrix_set(dcdr_env%matrix_t1, 4)
      CALL dbcsr_allocate_matrix_set(dcdr_env%matrix_apply_op_constant, 1)
      CALL dbcsr_allocate_matrix_set(dcdr_env%matrix_d_vhxc_dR, 3)
      CALL dbcsr_allocate_matrix_set(dcdr_env%matrix_vhxc_perturbed_basis, nspins, 6)
      CALL dbcsr_allocate_matrix_set(dcdr_env%matrix_hc, 3)
      CALL dbcsr_allocate_matrix_set(dcdr_env%matrix_ppnl_1, 3)

      CALL dbcsr_init_p(dcdr_env%perturbed_dm_correction)
      CALL dbcsr_copy(dcdr_env%perturbed_dm_correction, matrix_ks(1)%matrix)

      ! temporary no_symmetry matrix:
      DO i = 1, 3
         CALL dbcsr_init_p(dcdr_env%matrix_nosym_temp(i)%matrix)
         CALL dbcsr_create(dcdr_env%matrix_nosym_temp(i)%matrix, template=matrix_ks(1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_alloc_block_from_nbl(dcdr_env%matrix_nosym_temp(i)%matrix, sab_all)
         CALL dbcsr_set(dcdr_env%matrix_nosym_temp(i)%matrix, 0._dp)

      END DO

      ! moments carry the result of build_local_moment_matrix
      DO i = 1, 3
         CALL dbcsr_init_p(dcdr_env%moments(i)%matrix)
         CALL dbcsr_copy(dcdr_env%moments(i)%matrix, matrix_ks(1)%matrix, "dcdr_env%moments")
         CALL dbcsr_set(dcdr_env%moments(i)%matrix, 0.0_dp)
      END DO
      CALL build_local_moment_matrix(qs_env, dcdr_env%moments, 1, ref_point=[0._dp, 0._dp, 0._dp])

      DO i = 1, 3
         DO j = 1, 3
            CALL dbcsr_init_p(dcdr_env%matrix_difdip(i, j)%matrix)
            CALL dbcsr_copy(dcdr_env%matrix_difdip(i, j)%matrix, dcdr_env%matrix_nosym_temp(1)%matrix)
            CALL dbcsr_set(dcdr_env%matrix_difdip(i, j)%matrix, 0.0_dp)
         END DO
      END DO

      CALL dbcsr_init_p(dcdr_env%hamiltonian1(1)%matrix)

      ! overlap/kinetic matrix: s(1)    normal overlap matrix;
      !                         s(2:4)  derivatives wrt. nuclear coordinates
      CALL dbcsr_init_p(dcdr_env%matrix_s1(1)%matrix)
      CALL dbcsr_init_p(dcdr_env%matrix_t1(1)%matrix)

      CALL dbcsr_copy(dcdr_env%matrix_s1(1)%matrix, matrix_s(1)%matrix)
      CALL dbcsr_copy(dcdr_env%matrix_t1(1)%matrix, dcdr_env%matrix_t(1)%matrix)

      DO i = 2, 4
         CALL dbcsr_init_p(dcdr_env%matrix_s1(i)%matrix)
         CALL dbcsr_copy(dcdr_env%matrix_s1(i)%matrix, dcdr_env%matrix_nosym_temp(1)%matrix)

         CALL dbcsr_init_p(dcdr_env%matrix_t1(i)%matrix)
         CALL dbcsr_copy(dcdr_env%matrix_t1(i)%matrix, dcdr_env%matrix_nosym_temp(1)%matrix)
      END DO

      CALL dbcsr_init_p(dcdr_env%matrix_apply_op_constant(1)%matrix)
      CALL dbcsr_copy(dcdr_env%matrix_apply_op_constant(1)%matrix, dcdr_env%matrix_s1(1)%matrix)

      ! j=1...3: derivative wrt nucleus A, 4...6: wrt nucleus B
      DO j = 1, 6
         CALL dbcsr_init_p(dcdr_env%matrix_vhxc_perturbed_basis(1, j)%matrix)
         CALL dbcsr_copy(dcdr_env%matrix_vhxc_perturbed_basis(1, j)%matrix, dcdr_env%matrix_s1(1)%matrix)
      END DO

      DO i = 1, 3
         CALL dbcsr_init_p(dcdr_env%matrix_hc(i)%matrix)
         CALL dbcsr_create(dcdr_env%matrix_hc(i)%matrix, template=matrix_ks(1)%matrix, &
                           matrix_type=dbcsr_type_symmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(dcdr_env%matrix_hc(i)%matrix, sab_orb)
         CALL dbcsr_set(dcdr_env%matrix_hc(i)%matrix, 0.0_dp)
      END DO

      DO i = 1, 3
         CALL dbcsr_init_p(dcdr_env%matrix_ppnl_1(i)%matrix)
         CALL dbcsr_create(dcdr_env%matrix_ppnl_1(i)%matrix, template=matrix_ks(1)%matrix, &
                           matrix_type=dbcsr_type_symmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(dcdr_env%matrix_ppnl_1(i)%matrix, sab_orb)
         CALL dbcsr_set(dcdr_env%matrix_ppnl_1(i)%matrix, 0.0_dp)
      END DO

      DO i = 1, 3
         CALL dbcsr_init_p(dcdr_env%matrix_d_vhxc_dR(i)%matrix)
         CALL dbcsr_copy(dcdr_env%matrix_d_vhxc_dR(i)%matrix, dcdr_env%matrix_s1(1)%matrix)

         CALL dbcsr_init_p(dcdr_env%matrix_core_charge_1(i)%matrix)
         CALL dbcsr_copy(dcdr_env%matrix_core_charge_1(i)%matrix, dcdr_env%matrix_s1(1)%matrix)
         CALL dbcsr_set(dcdr_env%matrix_core_charge_1(i)%matrix, 0.0_dp)
      END DO

      ! CHC
      CALL get_mo_set(mo_set=mos(1)%mo_set, mo_coeff=mo_coeff, nao=nao, nmo=nmo)
      CALL cp_fm_create(buf, dcdr_env%likemos_fm_struct)

      CALL cp_dbcsr_sm_fm_multiply(matrix_ks(1)%matrix, mo_coeff, buf, nmo)
      ! chc = mo * matrix_ks * mo
      CALL cp_fm_get_info(buf)
      CALL cp_fm_get_info(dcdr_env%chc(1)%matrix)
      CALL cp_gemm('T', 'N', nmo, nmo, nao, &
                   1.0_dp, mo_coeff, buf, &
                   0.0_dp, dcdr_env%chc(1)%matrix)

      CALL cp_fm_release(buf)

      CALL cp_print_key_finished_output(output_unit, logger, lr_section, &
                                        "PRINT%PROGRAM_RUN_INFO")

      IF (dcdr_env%output_unit > 0) &
         WRITE (dcdr_env%output_unit, *) 'List of atoms: ', dcdr_env%list_of_atoms

      CALL timestop(handle)
   END SUBROUTINE dcdr_env_init

! **************************************************************************************************
!> \brief Deallocate the dcdr environment
!> \param qs_env ...
!> \param dcdr_env ...
! **************************************************************************************************
   SUBROUTINE dcdr_env_cleanup(qs_env, dcdr_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dcdr_env_type)                                :: dcdr_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'dcdr_env_cleanup', &
         routineP = moduleN//':'//routineN

      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: dcdr_section

      ! Destroy the logger
      logger => cp_get_default_logger()
      dcdr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES%DCDR")
      CALL cp_print_key_finished_output(dcdr_env%output_unit, logger, dcdr_section, "PRINT%APT")

      DEALLOCATE (dcdr_env%list_of_atoms)

      CALL cp_fm_struct_release(dcdr_env%aoao_fm_struct)
      CALL cp_fm_struct_release(dcdr_env%momo_fm_struct)
      CALL cp_fm_struct_release(dcdr_env%homohomo_fm_struct)
      CALL cp_fm_struct_release(dcdr_env%likemos_fm_struct)

      DEALLOCATE (dcdr_env%deltar)
      DEALLOCATE (dcdr_env%delta_basis_function)

      IF (dcdr_env%localized_psi0) THEN
         ! DEALLOCATE (dcdr_env%psi0_order)
         DEALLOCATE (dcdr_env%centers_set(1)%array)
         DEALLOCATE (dcdr_env%center_list(1)%array)
         DEALLOCATE (dcdr_env%centers_set)
         DEALLOCATE (dcdr_env%center_list)
         DEALLOCATE (dcdr_env%apt_subset)
      END IF

      DEALLOCATE (dcdr_env%apt_el_dcdr)
      DEALLOCATE (dcdr_env%apt_nuc_dcdr)
      DEALLOCATE (dcdr_env%apt_total_dcdr)
      IF (dcdr_env%localized_psi0) THEN
         DEALLOCATE (dcdr_env%apt_el_dcdr_per_center)
         DEALLOCATE (dcdr_env%apt_el_dcdr_per_subset)
      END IF

      ! Full matrices
      CALL cp_fm_release(dcdr_env%dCR(1)%matrix)
      CALL cp_fm_release(dcdr_env%dCR_prime(1)%matrix)
      CALL cp_fm_release(dcdr_env%mo_coeff(1)%matrix)
      CALL cp_fm_release(dcdr_env%chc(1)%matrix)
      CALL cp_fm_release(dcdr_env%op_dR(1)%matrix)

      DEALLOCATE (dcdr_env%op_dR)
      DEALLOCATE (dcdr_env%dCR)
      DEALLOCATE (dcdr_env%dCR_prime)
      DEALLOCATE (dcdr_env%mo_coeff)
      DEALLOCATE (dcdr_env%chc)

      ! DBCSR matrices
      CALL dbcsr_deallocate_matrix(dcdr_env%perturbed_dm_correction)
      CALL dbcsr_deallocate_matrix_set(dcdr_env%hamiltonian1)
      CALL dbcsr_deallocate_matrix_set(dcdr_env%moments)
      CALL dbcsr_deallocate_matrix_set(dcdr_env%matrix_difdip)
      CALL dbcsr_deallocate_matrix_set(dcdr_env%matrix_core_charge_1)
      CALL dbcsr_deallocate_matrix_set(dcdr_env%matrix_nosym_temp)
      CALL dbcsr_deallocate_matrix_set(dcdr_env%matrix_s)
      CALL dbcsr_deallocate_matrix_set(dcdr_env%matrix_t)
      CALL dbcsr_deallocate_matrix_set(dcdr_env%matrix_s1)
      CALL dbcsr_deallocate_matrix_set(dcdr_env%matrix_t1)
      CALL dbcsr_deallocate_matrix_set(dcdr_env%matrix_apply_op_constant)
      CALL dbcsr_deallocate_matrix_set(dcdr_env%matrix_d_vhxc_dR)
      CALL dbcsr_deallocate_matrix_set(dcdr_env%matrix_vhxc_perturbed_basis)
      CALL dbcsr_deallocate_matrix_set(dcdr_env%matrix_hc)
      CALL dbcsr_deallocate_matrix_set(dcdr_env%matrix_ppnl_1)


   END SUBROUTINE dcdr_env_cleanup

END MODULE qs_dcdr_utils
