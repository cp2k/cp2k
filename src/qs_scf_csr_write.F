!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Functions to print the KS and S matrix in the CSR format to file
!> \par History
!>      Started as a copy from the relevant part of qs_scf_post_gpw
!> \author Fabian Ducry (05.2020)
! **************************************************************************************************
MODULE qs_scf_csr_write
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
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
        dbcsr_csr_type, dbcsr_csr_write, dbcsr_desymmetrize, dbcsr_finalize, dbcsr_get_block_p, &
        dbcsr_has_symmetry, dbcsr_p_type, dbcsr_put_block, dbcsr_release, dbcsr_set, dbcsr_type, &
        dbcsr_type_antisymmetric, dbcsr_type_complex_8, dbcsr_type_no_symmetry, dbcsr_type_real_8, &
        dbcsr_type_symmetric
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
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              get_neighbor_list_set_p,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
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

      CHARACTER(len=*), PARAMETER :: routineN = 'write_ks_matrix_csr'

      INTEGER                                            :: handle, output_unit
      LOGICAL                                            :: do_kpoints, do_ks_csr_write, real_space
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks
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
         CALL get_qs_env(qs_env=qs_env, kpoints=kpoints, matrix_ks_kp=matrix_ks, do_kpoints=do_kpoints)
         CALL section_vals_val_get(dft_section, "PRINT%KS_CSR_WRITE%REAL_SPACE", &
                                   l_val=real_space)

         IF (do_kpoints .AND. .NOT. real_space) THEN
            CALL write_matrix_kp_csr(mat=matrix_ks, dft_section=dft_section, &
                                     kpoints=kpoints, prefix="KS")
         ELSE
            CALL write_matrix_csr(dft_section, mat=matrix_ks, kpoints=kpoints, do_kpoints=do_kpoints, &
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

      CHARACTER(len=*), PARAMETER :: routineN = 'write_s_matrix_csr'

      INTEGER                                            :: handle, output_unit
      LOGICAL                                            :: do_kpoints, do_s_csr_write, real_space
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
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
         CALL get_qs_env(qs_env=qs_env, kpoints=kpoints, matrix_s_kp=matrix_s, do_kpoints=do_kpoints)
         CALL section_vals_val_get(dft_section, "PRINT%S_CSR_WRITE%REAL_SPACE", &
                                   l_val=real_space)

         IF (do_kpoints .AND. .NOT. real_space) THEN
            CALL write_matrix_kp_csr(mat=matrix_s, dft_section=dft_section, &
                                     kpoints=kpoints, prefix="S")
         ELSE
            CALL write_matrix_csr(dft_section, mat=matrix_s, kpoints=kpoints, do_kpoints=do_kpoints, &
                                  prefix="S")
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE write_s_matrix_csr

! **************************************************************************************************
!> \brief helper function to print the real space representation of KS or S matrix to file
!> \param dft_section the dft_section
!> \param mat Hamiltonian or overlap matrix
!> \param kpoints Kpoint environment
!> \param prefix string to distinguish between KS and S matrix
!> \param do_kpoints Whether it is a gamma-point or k-point calculation
!> \par History
!>       Moved most of the code from write_ks_matrix_csr and write_s_matrix_csr
!>       Removed the code for printing k-point dependent matrices and added
!>              printing of the real space representation
! **************************************************************************************************
   SUBROUTINE write_matrix_csr(dft_section, mat, kpoints, prefix, do_kpoints)
      TYPE(section_vals_type), INTENT(IN), POINTER       :: dft_section
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(IN), &
         POINTER                                         :: mat
      TYPE(kpoint_type), INTENT(IN), POINTER             :: kpoints
      CHARACTER(*), INTENT(in)                           :: prefix
      LOGICAL, INTENT(IN)                                :: do_kpoints

      CHARACTER(len=*), PARAMETER                        :: routineN = 'write_matrix_csr'

      CHARACTER(LEN=default_path_length)                 :: file_name, fileformat, subs_string
      INTEGER                                            :: handle, ic, ispin, ncell, nspin, &
                                                            output_unit, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cell_to_index
      INTEGER, DIMENSION(:, :), POINTER                  :: i2c
      INTEGER, DIMENSION(:, :, :), POINTER               :: c2i
      LOGICAL                                            :: bin, do_symmetric, uptr
      REAL(KIND=dp)                                      :: thld
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_csr_type)                               :: mat_csr
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_nosym
      TYPE(dbcsr_type)                                   :: matrix_nosym
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

      nspin = SIZE(mat, 1)
      ncell = SIZE(mat, 2)

      IF (do_kpoints) THEN

         i2c => kpoints%index_to_cell
         c2i => kpoints%cell_to_index

         NULLIFY (sab_nl)
         CALL get_kpoint_info(kpoints, sab_nl=sab_nl)
         CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl, symmetric=do_symmetric)

         ! desymmetrize the KS or S matrices if necessary
         IF (do_symmetric) THEN
            CALL desymmetrize_rs_matrix(mat, mat_nosym, cell_to_index, index_to_cell, kpoints)
            ncell = SIZE(index_to_cell, 2) ! update the number of cells
         ELSE
            ALLOCATE (cell_to_index(LBOUND(c2i, 1):UBOUND(c2i, 1), &
                                    LBOUND(c2i, 2):UBOUND(c2i, 2), &
                                    LBOUND(c2i, 3):UBOUND(c2i, 3)))
            cell_to_index(LBOUND(c2i, 1):UBOUND(c2i, 1), &
                          LBOUND(c2i, 2):UBOUND(c2i, 2), &
                          LBOUND(c2i, 3):UBOUND(c2i, 3)) = c2i

            ALLOCATE (index_to_cell(3, ncell))
            index_to_cell(1:3, 1:ncell) = i2c

            mat_nosym => mat
         END IF

         ! print the index to cell mapping to the output
         IF (output_unit > 0) THEN
            WRITE (output_unit, "(/,A15,T15,I4,A)") prefix//" CSR write| ", &
               ncell, " periodic images"
            WRITE (output_unit, "(T7,A,T17,A,T24,A,T31,A)") "Number", "X", "Y", "Z"
            DO ic = 1, ncell
               WRITE (output_unit, "(T8,I3,T15,I3,T22,I3,T29,I3)") ic, index_to_cell(:, ic)
            END DO
         END IF
      END IF

      ! write the csr file(s)
      DO ispin = 1, nspin
         DO ic = 1, ncell
            IF (do_kpoints) THEN
               CALL dbcsr_copy(matrix_nosym, mat_nosym(ispin, ic)%matrix)
               WRITE (file_name, '(2(A,I0))') prefix//"_SPIN_", ispin, "_R_", ic
            ELSE
               IF (dbcsr_has_symmetry(mat(ispin, ic)%matrix)) THEN
                  CALL dbcsr_desymmetrize(mat(ispin, ic)%matrix, matrix_nosym)
               ELSE
                  CALL dbcsr_copy(matrix_nosym, mat(ispin, ic)%matrix)
               END IF
               WRITE (file_name, '(A,I0)') prefix//"_SPIN_", ispin
            END IF
            ! Convert dbcsr to csr
            CALL dbcsr_csr_create_from_dbcsr(matrix_nosym, &
                                             mat_csr, dbcsr_csr_dbcsr_blkrow_dist)
            CALL dbcsr_convert_dbcsr_to_csr(matrix_nosym, mat_csr)
            ! Write to file
            unit_nr = cp_print_key_unit_nr(logger, dft_section, subs_string, &
                                           extension=".csr", middle_name=TRIM(file_name), &
                                           file_status="REPLACE", file_form=fileformat)
            CALL dbcsr_csr_write(mat_csr, unit_nr, upper_triangle=uptr, threshold=thld, binary=bin)

            CALL cp_print_key_finished_output(unit_nr, logger, dft_section, subs_string)
            CALL dbcsr_csr_destroy(mat_csr)
            CALL dbcsr_release(matrix_nosym)
         END DO
      END DO

      ! clean up
      IF (do_kpoints) THEN
         DEALLOCATE (cell_to_index, index_to_cell)
         IF (do_symmetric) THEN
            DO ispin = 1, nspin
               DO ic = 1, ncell
                  CALL dbcsr_release(mat_nosym(ispin, ic)%matrix)
               END DO
            END DO
            CALL dbcsr_deallocate_matrix_set(mat_nosym)
         END IF
      END IF
      CALL timestop(handle)

   END SUBROUTINE write_matrix_csr

! **************************************************************************************************
!> \brief helper function to print the k-dependent KS or S matrix to file
!> \param mat Hamiltonian or overlap matrix for k-point calculations
!> \param dft_section the dft_section
!> \param kpoints Kpoint environment
!> \param prefix string to distinguish between KS and S matrix
!> \par History
!>       Moved the code from write_matrix_csr to write_matrix_kp_csr
!> \author Fabian Ducry
! **************************************************************************************************
   SUBROUTINE write_matrix_kp_csr(mat, dft_section, kpoints, prefix)
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(IN), &
         POINTER                                         :: mat
      TYPE(section_vals_type), INTENT(IN), POINTER       :: dft_section
      TYPE(kpoint_type), INTENT(IN), POINTER             :: kpoints
      CHARACTER(*), INTENT(in)                           :: prefix

      CHARACTER(len=*), PARAMETER :: routineN = 'write_matrix_kp_csr'
      COMPLEX(KIND=dp), PARAMETER :: cone = CMPLX(1.0_dp, 0.0_dp, KIND=dp), &
         ione = CMPLX(0.0_dp, 1.0_dp, KIND=dp)

      CHARACTER(LEN=default_path_length)                 :: file_name, fileformat, subs_string
      INTEGER                                            :: handle, igroup, ik, ikp, ispin, kplocal, &
                                                            nkp_groups, output_unit, unit_nr
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
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

      NULLIFY (sab_nl)

      !  Calculate the Hamiltonian at the k-points
      CALL get_kpoint_info(kpoints, xkp=xkp, use_real_wfn=use_real_wfn, kp_range=kp_range, &
                           nkp_groups=nkp_groups, kp_dist=kp_dist, sab_nl=sab_nl)

      ALLOCATE (rmatrix)
      CALL dbcsr_create(rmatrix, template=mat(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_nl)

      IF (.NOT. use_real_wfn) THEN
         ! Allocate temporary variables
         ALLOCATE (rmatrix_nosym, imatrix, imatrix_nosym, cmatrix, tmatrix)
         CALL dbcsr_create(rmatrix_nosym, template=mat(1, 1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(imatrix, template=mat(1, 1)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric)
         CALL dbcsr_create(imatrix_nosym, template=mat(1, 1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(cmatrix, template=mat(1, 1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry, &
                           data_type=dbcsr_type_complex_8)
         CALL dbcsr_create(tmatrix, template=mat(1, 1)%matrix, &
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
         DO ispin = 1, SIZE(mat, 1)
            DO igroup = 1, nkp_groups
               ! number of current kpoint
               ik = kp_dist(1, igroup) + ikp - 1
               CALL dbcsr_set(rmatrix, 0.0_dp)
               IF (use_real_wfn) THEN
                  ! FT of the matrix
                  CALL rskp_transform(rmatrix=rmatrix, rsmat=mat, ispin=ispin, &
                                      xkp=xkp(1:3, ik), cell_to_index=kpoints%cell_to_index, sab_nl=sab_nl)
                  ! Convert to desymmetrized csr matrix
                  CALL dbcsr_desymmetrize(rmatrix, matrix_nosym)
                  CALL dbcsr_csr_create_from_dbcsr(matrix_nosym, mat_csr, dbcsr_csr_dbcsr_blkrow_dist)
                  CALL dbcsr_convert_dbcsr_to_csr(matrix_nosym, mat_csr)
                  CALL dbcsr_release(matrix_nosym)
               ELSE
                  ! FT of the matrix
                  CALL dbcsr_set(imatrix, 0.0_dp)
                  CALL rskp_transform(rmatrix=rmatrix, cmatrix=imatrix, rsmat=mat, ispin=ispin, &
                                      xkp=xkp(1:3, ik), cell_to_index=kpoints%cell_to_index, sab_nl=sab_nl)

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
      CALL timestop(handle)

   END SUBROUTINE write_matrix_kp_csr

! **************************************************************************************************
!> \brief Desymmetrizes the KS or S matrices which are stored in symmetric !matrices
!> \param mat Hamiltonian or overlap matrices
!> \param mat_nosym Desymmetrized Hamiltonian or overlap matrices
!> \param cell_to_index Mapping of cell indices to linear RS indices
!> \param index_to_cell Mapping of linear RS indices to cell indices
!> \param kpoints Kpoint environment
!> \author Fabian Ducry
! **************************************************************************************************
   SUBROUTINE desymmetrize_rs_matrix(mat, mat_nosym, cell_to_index, index_to_cell, kpoints)
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(IN), &
         POINTER                                         :: mat
      TYPE(dbcsr_p_type), DIMENSION(:, :), &
         INTENT(INOUT), POINTER                          :: mat_nosym
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: cell_to_index
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: index_to_cell
      TYPE(kpoint_type), INTENT(IN), POINTER             :: kpoints

      CHARACTER(len=*), PARAMETER :: routineN = 'desymmetrize_rs_matrix'

      INTEGER                                            :: handle, iatom, ic, icn, icol, irow, &
                                                            ispin, jatom, ncell, nomirror, nspin, &
                                                            nx, ny, nz
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:, :), POINTER                  :: i2c
      INTEGER, DIMENSION(:, :, :), POINTER               :: c2i
      LOGICAL                                            :: found, lwtr
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CALL timeset(routineN, handle)

      i2c => kpoints%index_to_cell
      c2i => kpoints%cell_to_index

      ncell = SIZE(i2c, 2)
      nspin = SIZE(mat, 1)

      nx = MAX(ABS(LBOUND(c2i, 1)), ABS(UBOUND(c2i, 1)))
      ny = MAX(ABS(LBOUND(c2i, 2)), ABS(UBOUND(c2i, 3)))
      nz = MAX(ABS(LBOUND(c2i, 3)), ABS(UBOUND(c2i, 3)))
      ALLOCATE (cell_to_index(-nx:nx, -ny:ny, -nz:nz))
      cell_to_index(LBOUND(c2i, 1):UBOUND(c2i, 1), &
                    LBOUND(c2i, 2):UBOUND(c2i, 2), &
                    LBOUND(c2i, 3):UBOUND(c2i, 3)) = c2i

      ! identify cells with no mirror img
      nomirror = 0
      DO ic = 1, ncell
         cell = i2c(:, ic)
         IF (cell_to_index(-cell(1), -cell(2), -cell(3)) == 0) &
            nomirror = nomirror + 1
      END DO

      ! create the mirror imgs
      ALLOCATE (index_to_cell(3, ncell + nomirror))
      index_to_cell(:, 1:ncell) = i2c

      nomirror = 0 ! count the imgs without mirror
      DO ic = 1, ncell
         cell = index_to_cell(:, ic)
         IF (cell_to_index(-cell(1), -cell(2), -cell(3)) == 0) THEN
            nomirror = nomirror + 1
            index_to_cell(:, ncell + nomirror) = -cell
            cell_to_index(-cell(1), -cell(2), -cell(3)) = ncell + nomirror
         END IF
      END DO
      ncell = ncell + nomirror

      CALL get_kpoint_info(kpoints, sab_nl=sab_nl)
      ! allocate the nonsymmetric matrices
      NULLIFY (mat_nosym)
      CALL dbcsr_allocate_matrix_set(mat_nosym, nspin, ncell)
      DO ispin = 1, nspin
         DO ic = 1, ncell
            ALLOCATE (mat_nosym(ispin, ic)%matrix)
            CALL dbcsr_create(matrix=mat_nosym(ispin, ic)%matrix, &
                              template=mat(1, 1)%matrix, &
                              matrix_type=dbcsr_type_no_symmetry, &
                              data_type=dbcsr_type_real_8)
            CALL cp_dbcsr_alloc_block_from_nbl(mat_nosym(ispin, ic)%matrix, &
                                               sab_nl, desymmetrize=.TRUE.)
            CALL dbcsr_set(mat_nosym(ispin, ic)%matrix, 0.0_dp)
         END DO
      END DO

      DO ispin = 1, nspin
         ! desymmetrize the matrix for real space printing
         CALL neighbor_list_iterator_create(nl_iterator, sab_nl)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, cell=cell)

            ic = cell_to_index(cell(1), cell(2), cell(3))
            icn = cell_to_index(-cell(1), -cell(2), -cell(3))
            CPASSERT(icn > 0)

            irow = iatom
            icol = jatom
            lwtr = .FALSE.
            ! always copy from the top
            IF (iatom > jatom) THEN
               irow = jatom
               icol = iatom
               lwtr = .TRUE.
            END IF

            CALL dbcsr_get_block_p(matrix=mat(ispin, ic)%matrix, &
                                   row=irow, col=icol, block=block, found=found)
            CPASSERT(found)

            ! copy to M(R) at (iatom,jatom)
            ! copy to M(-R) at (jatom,iatom)
            IF (lwtr) THEN
               CALL dbcsr_put_block(matrix=mat_nosym(ispin, ic)%matrix, &
                                    row=iatom, col=jatom, block=TRANSPOSE(block))
               CALL dbcsr_put_block(matrix=mat_nosym(ispin, icn)%matrix, &
                                    row=jatom, col=iatom, block=block)
            ELSE
               CALL dbcsr_put_block(matrix=mat_nosym(ispin, ic)%matrix, &
                                    row=iatom, col=jatom, block=block)
               CALL dbcsr_put_block(matrix=mat_nosym(ispin, icn)%matrix, &
                                    row=jatom, col=iatom, block=TRANSPOSE(block))
            END IF
         END DO
         CALL neighbor_list_iterator_release(nl_iterator)
      END DO

      DO ispin = 1, nspin
         DO ic = 1, ncell
            CALL dbcsr_finalize(mat_nosym(ispin, ic)%matrix)
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE desymmetrize_rs_matrix

END MODULE qs_scf_csr_write
