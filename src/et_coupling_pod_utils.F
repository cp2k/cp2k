!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief different utilities used by POD method routines
!> \author Z. Futera (03.2019)
! **************************************************************************************************
MODULE et_coupling_pod_utils

   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_basic_linalg,              ONLY: cp_fm_syrk,&
                                              cp_fm_upper_to_full
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_generate_filename
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE et_coupling_pod_save,            ONLY: pod_save_transf_mat
   USE et_coupling_pod_types,           ONLY: et_cpl
   USE input_section_types,             ONLY: section_get_lval,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling_pod_utils'

   PUBLIC :: pod_check_restart, &
             pod_get_transf_mat

CONTAINS

! **************************************************************************************************
!> \brief calculated transformation matrix for basis-set orthogonalization (S^{-1/2})
!> \param output_unit unit number of the open output stream
!> \param qs_env QuickStep environment containing all system data
!> \param ec electronic coupling data structure
!> \param mat_w working matrix of the same dimension
!> \param n_ao number of atomic orbitLals (dimension of S matrix)
!> \author Z. Futera (02.2017)
! **************************************************************************************************
   SUBROUTINE pod_get_transf_mat(output_unit, qs_env, ec, mat_w, n_ao)

      ! Routine arguments
      INTEGER                                            :: output_unit
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(et_cpl), POINTER                              :: ec
      TYPE(cp_fm_type), POINTER                          :: mat_w
      INTEGER                                            :: n_ao

      ! Routine name for debug purposes
      CHARACTER(len=*), PARAMETER :: routineN = 'pod_get_transf_mat', &
                                     routineP = moduleN//':'//routineN

      ! Local variables
      INTEGER                                            :: i
      INTEGER                                            :: n_deps
      REAL(KIND=dp)                                      :: f
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: et, ei
      TYPE(cp_fm_type), POINTER                          :: mat_u
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_s
      TYPE(scf_control_type), POINTER                    :: scf_cntrl
#if defined(__SCALAPACK)
      INTEGER, EXTERNAL                                  :: indxg2l, indxg2p

      INTEGER                                            :: j, ic, ir, ic_loc, ir_loc
      INTEGER                                            :: n_b_row, n_b_col
      INTEGER                                            :: n_p_row, n_p_col
      INTEGER                                            :: p_row_0, p_col_0
      INTEGER                                            :: p_row_i, p_col_i

      ! Number of rows/columns in Scalapack blocks
      n_b_row = mat_w%matrix_struct%nrow_block
      n_b_col = mat_w%matrix_struct%ncol_block
      ! Number of processors
      n_p_row = mat_w%matrix_struct%context%num_pe(1)
      n_p_col = mat_w%matrix_struct%context%num_pe(2)
      ! Position of the first processor
      p_row_0 = mat_w%matrix_struct%first_p_pos(1)
      p_col_0 = mat_w%matrix_struct%first_p_pos(2)
      ! Position of the actual processor
      p_row_i = mat_w%matrix_struct%context%mepos(1)
      p_col_i = mat_w%matrix_struct%context%mepos(2)
#endif

      ! Initialization
      NULLIFY (mat_u)
      NULLIFY (mat_s)
      NULLIFY (scf_cntrl)

      ALLOCATE (et(n_ao))
      ALLOCATE (ei(n_ao))
      et(:) = 0.0_dp
      ei(:) = 0.0_dp

      ! Diagonalization of S
      CALL get_qs_env(qs_env, matrix_s=mat_s)
      CALL cp_fm_create(matrix=mat_u, matrix_struct=mat_w%matrix_struct, &
                        name='OVERLAP MATRIX EIGENVECTORS')
      CALL copy_dbcsr_to_fm(mat_s(1)%matrix, ec%tmatf)
      CALL choose_eigv_solver(ec%tmatf, mat_w, et)
      CALL cp_fm_to_fm(ec%tmatf, ec%tmatr)
      ei(:) = et(:)

      ! Eigenvector quenching
      CALL get_qs_env(qs_env, scf_control=scf_cntrl)

      n_deps = 0
      DO i = 1, n_ao
         IF (et(i) < scf_cntrl%eps_eigval) THEN
#if defined(__SCALAPACK)
            ic = indxg2p(i, n_b_col, p_col_i, p_col_0, n_p_col)
            IF (p_col_i == ic) THEN
               ic_loc = indxg2l(i, n_b_col, p_col_i, p_col_0, n_p_col)
               DO j = 1, n_ao
                  ir = indxg2p(j, n_b_row, p_row_i, p_row_0, n_p_row)
                  IF (p_row_i == ir) THEN
                     ir_loc = indxg2l(j, n_b_row, p_row_i, p_row_0, n_p_row)
                     mat_w%local_data(ir_loc, ic_loc) = 0.0_dp
                  END IF
               END DO
            END IF
#else
            mat_w%local_data(1:n_ao, i) = 0.0_dp
#endif
            n_deps = n_deps + 1
         END IF
      END DO

      IF (n_deps /= 0) THEN
         CALL cp_warn(__LOCATION__, &
                      "Overlap matrix exhibits linear dependencies. At least some "// &
                      "eigenvalues have been quenched.")
      END IF

      CALL cp_fm_to_fm(mat_w, mat_u)

      ! Forward ransformations S^{-1/2}
      IF (output_unit > 0) &
         WRITE (output_unit, '(/,T3,A)') &
         'Constructing forward basis-set transformation matrix'

      DO i = 1, n_ao
         IF (et(i) >= scf_cntrl%eps_eigval) THEN
            f = et(i)**(-0.25_dp)
#if defined(__SCALAPACK)
            ic = indxg2p(i, n_b_col, p_col_i, p_col_0, n_p_col)
            IF (p_col_i == ic) THEN
               ic_loc = indxg2l(i, n_b_col, p_col_i, p_col_0, n_p_col)
               DO j = 1, n_ao
                  ir = indxg2p(j, n_b_row, p_row_i, p_row_0, n_p_row)
                  IF (p_row_i == ir) THEN
                     ir_loc = indxg2l(j, n_b_row, p_row_i, p_row_0, n_p_row)
                     mat_w%local_data(ir_loc, ic_loc) = f*mat_w%local_data(ir_loc, ic_loc)
                  END IF
               END DO
            END IF
#else
            mat_w%local_data(1:n_ao, i) = f*mat_w%local_data(1:n_ao, i)
#endif
         END IF
      END DO

      CALL cp_fm_syrk("U", "N", n_ao, 1.0_dp, mat_w, 1, 1, 0.0_dp, ec%tmatf)
      CALL cp_fm_upper_to_full(ec%tmatf, mat_w)

      ! Save the forward-transformation matrix to file
      CALL pod_save_transf_mat(qs_env, ec, 'F')

      ! Reverse transformation S^{+1/2}
      IF (output_unit > 0) &
         WRITE (output_unit, '(T3,A)') &
         'Constructing reverse basis-set transformation matrix'

      CALL cp_fm_to_fm(mat_u, mat_w)
      DO i = 1, n_ao
         IF (ei(i) >= scf_cntrl%eps_eigval) THEN
            f = ei(i)**(+0.25_dp)
#if defined(__SCALAPACK)
            ic = indxg2p(i, n_b_col, p_col_i, p_col_0, n_p_col)
            IF (p_col_i == ic) THEN
               ic_loc = indxg2l(i, n_b_col, p_col_i, p_col_0, n_p_col)
               DO j = 1, n_ao
                  ir = indxg2p(j, n_b_row, p_row_i, p_row_0, n_p_row)
                  IF (p_row_i == ir) THEN
                     ir_loc = indxg2l(j, n_b_row, p_row_i, p_row_0, n_p_row)
                     mat_w%local_data(ir_loc, ic_loc) = f*mat_w%local_data(ir_loc, ic_loc)
                  END IF
               END DO
            END IF
#else
            mat_w%local_data(1:n_ao, i) = f*mat_w%local_data(1:n_ao, i)
#endif
         END IF
      END DO

      CALL cp_fm_syrk("U", "N", n_ao, 1.0_dp, mat_w, 1, 1, 0.0_dp, ec%tmatr)
      CALL cp_fm_upper_to_full(ec%tmatr, mat_w)

      ! Save the backward-transformation matrix to file
      CALL pod_save_transf_mat(qs_env, ec, 'B')

      ! Clean memory
      CALL cp_fm_release(matrix=mat_u)
      DEALLOCATE (et)
      DEALLOCATE (ei)

   END SUBROUTINE pod_get_transf_mat

! **************************************************************************************************
!> \brief check whether the POD data are read from the restart files
!> \param logger output stream data
!> \param ec_sec POD input-file section
!> \param file_tmf name of the POD restart file with forward transformation matrix
!> \param file_tmr name of the POD restart file with reverse transformation matrix
!> \param file_hab name of the POD restart file with Hamiltonian
!> \param file_wfn name of the POD restart file with MOs
!> \return ...
!> \author Z. Futera (03.2019)
! **************************************************************************************************
   FUNCTION pod_check_restart(logger, ec_sec, file_tmf, file_tmr, file_hab, file_wfn) RESULT(restart)

      ! Routine arguments
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: ec_sec
      CHARACTER(LEN=default_path_length), OPTIONAL       :: file_tmf, file_tmr, file_hab, file_wfn
      LOGICAL                                            :: restart

      CHARACTER(len=*), PARAMETER :: routineN = 'pod_check_restart', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=default_path_length)                 :: filename
      LOGICAL                                            :: ok
      TYPE(section_vals_type), POINTER                   :: print_key

! Returning value
! Routine name for debug purposes
! Local variables

      ! restart
      restart = section_get_lval(ec_sec, 'RESTART')

      IF (restart) THEN

         ! Forward transformation matrix
         IF (PRESENT(file_tmf)) THEN
            print_key => section_vals_get_subs_vals(ec_sec, &
                                                    'PRINT%TRANSF_MATRIX_F')
            filename = cp_print_key_generate_filename(logger, print_key, &
                                                      middle_name='RESTART', extension='.tmf', my_local=.FALSE.)
            INQUIRE (FILE=filename, exist=ok)
            IF (.NOT. ok) THEN
               CPABORT('restart file "'//TRIM(ADJUSTL(filename))//'" not found')
            ELSE
               file_tmf = filename
            END IF
         END IF

         ! Reverse transformation matrix
         IF (PRESENT(file_tmf)) THEN
            print_key => section_vals_get_subs_vals(ec_sec, &
                                                    'PRINT%TRANSF_MATRIX_R')
            filename = cp_print_key_generate_filename(logger, print_key, &
                                                      middle_name='RESTART', extension='.tmf', my_local=.FALSE.)
            INQUIRE (FILE=filename, exist=ok)
            IF (.NOT. ok) THEN
               CPABORT('restart file "'//TRIM(ADJUSTL(filename))//'" not found')
            ELSE
               file_tmr = filename
            END IF
         END IF

         ! Hamiltonian restart file
         IF (PRESENT(file_hab)) THEN
            print_key => section_vals_get_subs_vals(ec_sec, &
                                                    'PRINT%HAMILTONIAN')
            filename = cp_print_key_generate_filename(logger, print_key, &
                                                      middle_name='RESTART', extension='.hab', my_local=.FALSE.)
            INQUIRE (FILE=filename, exist=ok)
            IF (.NOT. ok) THEN
               CPABORT('restart file "'//TRIM(ADJUSTL(filename))//'" not found')
            ELSE
               file_hab = filename
            END IF
         END IF

         ! MO restart file
         IF (PRESENT(file_wfn)) THEN
            print_key => section_vals_get_subs_vals(ec_sec, &
                                                    'PRINT%WAVEFUNCTION')
            filename = cp_print_key_generate_filename(logger, print_key, &
                                                      middle_name='RESTART', extension='.wfn', my_local=.FALSE.)
            INQUIRE (FILE=filename, exist=ok)
            IF (.NOT. ok) THEN
               CPABORT('restart file "'//TRIM(ADJUSTL(filename))//'" not found')
            ELSE
               file_wfn = filename
            END IF
         END IF

      END IF

   END FUNCTION pod_check_restart

END MODULE et_coupling_pod_utils
