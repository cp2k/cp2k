!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Functions to print the KS and S matrix in the CSR format to file
!> \par History
!>      Started as a copy from the relevant part of qs_scf_post_gpw
!> \author Fabian Ducry (05.2020)
! **************************************************************************************************
MODULE qs_scf_csr_write
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_convert_dbcsr_to_csr, dbcsr_copy, dbcsr_create, &
        dbcsr_csr_create_from_dbcsr, dbcsr_csr_dbcsr_blkrow_dist, dbcsr_csr_destroy, &
        dbcsr_csr_type, dbcsr_csr_write, dbcsr_desymmetrize, dbcsr_has_symmetry, dbcsr_p_type, &
        dbcsr_release, dbcsr_set, dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_complex_8, &
        dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE kpoint_methods,                  ONLY: rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   ! Global parameters
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_csr_write'
   PUBLIC :: write_ks_matrix_csr, &
             write_s_matrix_csr

! **************************************************************************************************

CONTAINS

!**************************************************************************************************
!> \brief writing the KS matrix in csr format into a file
!> \param qs_env qs environment
!> \param input the input
!> \par History
!>       Moved to module qs_scf_csr_write (05.2020)
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE write_ks_matrix_csr(qs_env, input)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input

      CHARACTER(len=*), PARAMETER :: routineN = 'write_ks_matrix_csr', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, output_unit
      LOGICAL                                            :: do_kpoints, do_ks_csr_write
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks_kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(section_vals_type), POINTER                   :: dft_section

      CALL timeset(routineN, handle)

      NULLIFY (dft_section)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      dft_section => section_vals_get_subs_vals(input, "DFT")
      do_ks_csr_write = BTEST(cp_print_key_should_output(logger%iter_info, dft_section, &
                                                         "PRINT%KS_CSR_WRITE"), cp_p_file)

      IF (do_ks_csr_write) THEN
         CALL get_qs_env(qs_env=qs_env, do_kpoints=do_kpoints)

         IF (do_kpoints) THEN
            CALL get_qs_env(qs_env=qs_env, kpoints=kpoints, matrix_ks_kp=matrix_ks_kp)
            CALL write_matrix_csr(dft_section, matrix_to_print_kp=matrix_ks_kp, &
                                  kpoints=kpoints, prefix="KS")
         ELSE
            CALL get_qs_env(qs_env=qs_env, matrix_ks=matrix_ks)
            CALL write_matrix_csr(dft_section, matrix_to_print=matrix_ks, &
                                  prefix="KS")
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_ks_matrix_csr

!**************************************************************************************************
!> \brief writing the overlap matrix in csr format into a file
!> \param qs_env qs environment
!> \param input the input
!> \par History
!>      Moved to module qs_scf_csr_write
!> \author Mohammad Hossein Bani-Hashemian
! **************************************************************************************************
   SUBROUTINE write_s_matrix_csr(qs_env, input)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input

      CHARACTER(len=*), PARAMETER :: routineN = 'write_s_matrix_csr', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, output_unit
      LOGICAL                                            :: do_kpoints, do_s_csr_write
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(section_vals_type), POINTER                   :: dft_section

      CALL timeset(routineN, handle)

      NULLIFY (dft_section)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      dft_section => section_vals_get_subs_vals(input, "DFT")
      do_s_csr_write = BTEST(cp_print_key_should_output(logger%iter_info, dft_section, &
                                                        "PRINT%S_CSR_WRITE"), cp_p_file)

      IF (do_s_csr_write) THEN
         CALL get_qs_env(qs_env=qs_env, do_kpoints=do_kpoints)

         IF (do_kpoints) THEN
            CALL get_qs_env(qs_env=qs_env, kpoints=kpoints, matrix_s_kp=matrix_s_kp)
            CALL write_matrix_csr(dft_section, matrix_to_print_kp=matrix_s_kp, &
                                  kpoints=kpoints, prefix="S")
         ELSE
            CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
            CALL write_matrix_csr(dft_section, matrix_to_print=matrix_s, &
                                  prefix="S")
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_s_matrix_csr

! **************************************************************************************************
!> \brief helper function to print the KS or S matrix to file
!> \param dft_section the dft_section
!> \param matrix_to_print Hamiltonian or overlap matrix for Gamma point calculations
!> \param matrix_to_print_kp Hamiltonian or overlap matrix for k-point calculations
!> \param kpoints the kpoints
!> \param prefix string to distinguish between KS and S matrix
!> \par History
!>       Moved most of the code from write_ks_matrix_csr and write_s_matrix_csr
! **************************************************************************************************
   SUBROUTINE write_matrix_csr(dft_section, matrix_to_print, matrix_to_print_kp, kpoints, prefix)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: dft_section
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(IN), &
         OPTIONAL, POINTER                               :: matrix_to_print
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL, POINTER                               :: matrix_to_print_kp
      TYPE(kpoint_type), INTENT(IN), OPTIONAL, POINTER   :: kpoints
      CHARACTER(*), INTENT(in)                           :: prefix

      CHARACTER(len=*), PARAMETER :: routineN = 'write_matrix_csr', &
         routineP = moduleN//':'//routineN
      COMPLEX(KIND=dp), PARAMETER :: cone = CMPLX(1.0_dp, 0.0_dp, KIND=dp), &
         ione = CMPLX(0.0_dp, 1.0_dp, KIND=dp)

      CHARACTER(LEN=default_path_length)                 :: file_name, fileformat, subs_string
      INTEGER                                            :: handle, igroup, ik, ikp, ispin, kplocal, &
                                                            nkp_groups, output_unit, unit_nr
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: bin, uptr, use_real_wfn
      REAL(KIND=dp)                                      :: thld
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_csr_type)                               :: mat_csr
      TYPE(dbcsr_type)                                   :: matrix_nosym
      TYPE(dbcsr_type), POINTER                          :: cmatrix, imatrix, imatrix_nosym, &
                                                            rmatrix, rmatrix_nosym, tmatrix
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      subs_string = "PRINT%"//prefix//"_CSR_WRITE"

      CALL section_vals_val_get(dft_section, subs_string//"%THRESHOLD", r_val=thld)
      CALL section_vals_val_get(dft_section, subs_string//"%UPPER_TRIANGULAR", l_val=uptr)
      CALL section_vals_val_get(dft_section, subs_string//"%BINARY", l_val=bin)

      IF (bin) THEN
         fileformat = "UNFORMATTED"
      ELSE
         fileformat = "FORMATTED"
      END IF

      IF (PRESENT(kpoints)) THEN
         !  Calculate the Hamiltonian at the k-points
         NULLIFY (sab_nl)
         CALL get_kpoint_info(kpoints, xkp=xkp, use_real_wfn=use_real_wfn, kp_range=kp_range, &
                              nkp_groups=nkp_groups, kp_dist=kp_dist, sab_nl=sab_nl, &
                              cell_to_index=cell_to_index)

         ALLOCATE (rmatrix)
         CALL dbcsr_create(rmatrix, template=matrix_to_print_kp(1, 1)%matrix, &
                           matrix_type=dbcsr_type_symmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_nl)

         IF (.NOT. use_real_wfn) THEN
            ! Allocate temporary variables
            ALLOCATE (rmatrix_nosym, imatrix, imatrix_nosym, cmatrix, tmatrix)
            CALL dbcsr_create(rmatrix_nosym, template=matrix_to_print_kp(1, 1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(imatrix, template=matrix_to_print_kp(1, 1)%matrix, &
                              matrix_type=dbcsr_type_antisymmetric)
            CALL dbcsr_create(imatrix_nosym, template=matrix_to_print_kp(1, 1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(cmatrix, template=matrix_to_print_kp(1, 1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry, &
                              data_type=dbcsr_type_complex_8)
            CALL dbcsr_create(tmatrix, template=matrix_to_print_kp(1, 1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry, &
                              data_type=dbcsr_type_complex_8)
            CALL cp_dbcsr_alloc_block_from_nbl(rmatrix_nosym, sab_nl)
            CALL cp_dbcsr_alloc_block_from_nbl(imatrix, sab_nl)
            CALL cp_dbcsr_alloc_block_from_nbl(imatrix_nosym, sab_nl)
            CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_nl)
            CALL cp_dbcsr_alloc_block_from_nbl(tmatrix, sab_nl)
         END IF

         kplocal = kp_range(2) - kp_range(1) + 1
         DO ikp = 1, kplocal
            DO ispin = 1, SIZE(matrix_to_print_kp, 1)
               DO igroup = 1, nkp_groups
                  ! number of current kpoint
                  ik = kp_dist(1, igroup) + ikp - 1
                  CALL dbcsr_set(rmatrix, 0.0_dp)
                  IF (use_real_wfn) THEN
                     ! FT of KS matrix
                     CALL rskp_transform(rmatrix=rmatrix, rsmat=matrix_to_print_kp, ispin=ispin, &
                                         xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_nl)
                     ! Convert to desymmetrized csr matrix
                     CALL dbcsr_desymmetrize(rmatrix, matrix_nosym)
                     CALL dbcsr_csr_create_from_dbcsr(matrix_nosym, mat_csr, dbcsr_csr_dbcsr_blkrow_dist)
                     CALL dbcsr_convert_dbcsr_to_csr(matrix_nosym, mat_csr)
                     CALL dbcsr_release(matrix_nosym)
                  ELSE
                     ! FT of KS matrix
                     CALL dbcsr_set(imatrix, 0.0_dp)
                     CALL rskp_transform(rmatrix=rmatrix, cmatrix=imatrix, rsmat=matrix_to_print_kp, ispin=ispin, &
                                         xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_nl)

                     ! Desymmetrize and sum the real and imaginary part into
                     ! cmatrix
                     CALL dbcsr_desymmetrize(rmatrix, rmatrix_nosym)
                     CALL dbcsr_desymmetrize(imatrix, imatrix_nosym)
                     CALL dbcsr_copy(cmatrix, rmatrix_nosym)
                     CALL dbcsr_copy(tmatrix, imatrix_nosym)
                     CALL dbcsr_add(cmatrix, tmatrix, cone, ione)
                     ! Convert to csr
                     CALL dbcsr_csr_create_from_dbcsr(cmatrix, mat_csr, dbcsr_csr_dbcsr_blkrow_dist)
                     CALL dbcsr_convert_dbcsr_to_csr(cmatrix, mat_csr)
                  END IF
                  ! Write to file
                  WRITE (file_name, '(2(A,I0))') prefix//"_SPIN_", ispin, "_K_", ik
                  unit_nr = cp_print_key_unit_nr(logger, dft_section, subs_string, &
                                                 extension=".csr", middle_name=TRIM(file_name), &
                                                 file_status="REPLACE", file_form=fileformat)
                  CALL dbcsr_csr_write(mat_csr, unit_nr, upper_triangle=uptr, threshold=thld, binary=bin)

                  CALL cp_print_key_finished_output(unit_nr, logger, dft_section, subs_string)

                  CALL dbcsr_csr_destroy(mat_csr)
               END DO
            END DO
         END DO
         CALL dbcsr_release(rmatrix)
         DEALLOCATE (rmatrix)
         IF (.NOT. use_real_wfn) THEN
            CALL dbcsr_release(rmatrix_nosym)
            CALL dbcsr_release(imatrix)
            CALL dbcsr_release(imatrix_nosym)
            CALL dbcsr_release(cmatrix)
            CALL dbcsr_release(tmatrix)
            DEALLOCATE (rmatrix_nosym, imatrix, imatrix_nosym, cmatrix, tmatrix)
         END IF
      ELSE
         DO ispin = 1, SIZE(matrix_to_print)
            ! Desymmetrize
            IF (dbcsr_has_symmetry(matrix_to_print(ispin)%matrix)) THEN
               CALL dbcsr_desymmetrize(matrix_to_print(ispin)%matrix, matrix_nosym)
            ELSE
               CALL dbcsr_copy(matrix_nosym, matrix_to_print(ispin)%matrix)
            END IF
            ! Convert dbcsr to csr
            CALL dbcsr_csr_create_from_dbcsr(matrix_nosym, mat_csr, dbcsr_csr_dbcsr_blkrow_dist)
            CALL dbcsr_convert_dbcsr_to_csr(matrix_nosym, mat_csr)
            ! Write to file
            WRITE (file_name, '(A,I0)') prefix//"_SPIN_", ispin
            unit_nr = cp_print_key_unit_nr(logger, dft_section, subs_string, &
                                           extension=".csr", middle_name=TRIM(file_name), &
                                           file_status="REPLACE", file_form=fileformat)
            CALL dbcsr_csr_write(mat_csr, unit_nr, upper_triangle=uptr, threshold=thld, binary=bin)

            CALL cp_print_key_finished_output(unit_nr, logger, dft_section, subs_string)

            CALL dbcsr_csr_destroy(mat_csr)
            CALL dbcsr_release(matrix_nosym)
         END DO
      END IF
      CALL timestop(handle)

   END SUBROUTINE write_matrix_csr

END MODULE qs_scf_csr_write
