!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE xc_gauxc_functional
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_create,&
                                              dbcsr_finalize,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_release
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE input_constants,                 ONLY: xc_vdw_fun_none
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE iso_c_binding,                   ONLY: c_double
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE particle_types,                  ONLY: particle_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE input_constants,                 ONLY: xc_vdw_fun_none
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_gauxc_functional'

   PUBLIC :: apply_gauxc

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_mat ...
!> \param dense_mat ...
! **************************************************************************************************
   SUBROUTINE dbcsr_to_dense(dbcsr_mat, dense_mat)
      USE cp_dbcsr_api, ONLY: dbcsr_get_info, dbcsr_iterator_type, &
                              dbcsr_iterator_start, dbcsr_iterator_next_block, dbcsr_iterator_blocks_left, &
                              dbcsr_iterator_stop
      TYPE(dbcsr_p_type), INTENT(IN)                     :: dbcsr_mat
      REAL(c_double), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: dense_mat

      INTEGER                                            :: col, i, icol, irow, j, nblkcols_total, &
                                                            nblkrows_total, ncols, nrows, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: c_offset, r_offset
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      REAL(c_double), POINTER                            :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_get_info(dbcsr_mat%matrix, &
                          row_blk_size=row_blk_size, &
                          col_blk_size=col_blk_size, &
                          nblkrows_total=nblkrows_total, &
                          nblkcols_total=nblkcols_total, &
                          nfullrows_total=nrows, &
                          nfullcols_total=ncols)

      IF (.NOT. ALLOCATED(dense_mat)) THEN
         ALLOCATE (dense_mat(nrows, ncols))
         dense_mat = 0.0_c_double
      END IF

      ALLOCATE (r_offset(nblkrows_total), c_offset(nblkcols_total))

      r_offset(1) = 1
      DO row = 2, nblkrows_total
         r_offset(row) = r_offset(row - 1) + row_blk_size(row - 1)
      END DO
      c_offset(1) = 1
      DO col = 2, nblkcols_total
         c_offset(col) = c_offset(col - 1) + col_blk_size(col - 1)
      END DO

      CALL dbcsr_iterator_start(iter, dbcsr_mat%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, block)
         DO j = 1, SIZE(block, 2)
            DO i = 1, SIZE(block, 1)
               dense_mat(r_offset(irow) + i - 1, c_offset(icol) + j - 1) = &
                  dense_mat(r_offset(irow) + i - 1, c_offset(icol) + j - 1) + block(i, j)
            END DO
         END DO
      END DO
      CALL dbcsr_iterator_stop(iter)

      DEALLOCATE (r_offset, c_offset)

   END SUBROUTINE dbcsr_to_dense

   ! **************************************************************************************************
   !> \brief Convert a dense matrix to dbcsr format using a template matrix
   !> \param dense_mat Input dense matrix
   !> \param template_dbcsr Template dbcsr matrix for structure
   !> \param dbcsr_mat Output dbcsr matrix
   ! **************************************************************************************************
! **************************************************************************************************
!> \brief ...
!> \param dense_mat ...
!> \param template_dbcsr ...
!> \return ...
! **************************************************************************************************
   FUNCTION dense_to_dbcsr(dense_mat, template_dbcsr) RESULT(dbcsr_mat)
      USE cp_dbcsr_api, ONLY: dbcsr_get_info, dbcsr_iterator_type, &
                              dbcsr_iterator_start, dbcsr_iterator_next_block, &
                              dbcsr_iterator_blocks_left, dbcsr_iterator_stop, &
                              dbcsr_put_block, dbcsr_finalize, dbcsr_work_create, &
                              dbcsr_init_p, dbcsr_create, dbcsr_copy, dbcsr_release
      REAL(c_double), DIMENSION(:, :), INTENT(IN)        :: dense_mat
      TYPE(dbcsr_p_type), INTENT(IN)                     :: template_dbcsr
      TYPE(dbcsr_p_type)                                 :: dbcsr_mat

      INTEGER                                            :: col, icol, irow, nblkcols_total, &
                                                            nblkrows_total, ncols, nrows, row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: c_offset, r_offset
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      REAL(c_double), POINTER                            :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_get_info(template_dbcsr%matrix, &
                          row_blk_size=row_blk_size, &
                          col_blk_size=col_blk_size, &
                          nblkrows_total=nblkrows_total, &
                          nblkcols_total=nblkcols_total, &
                          nfullrows_total=nrows, &
                          nfullcols_total=ncols)

      CALL dbcsr_init_p(dbcsr_mat%matrix)
      CALL dbcsr_create(dbcsr_mat%matrix, template=template_dbcsr%matrix, name="VXC from GauXC")
      CALL dbcsr_work_create(dbcsr_mat%matrix, work_mutable=.TRUE.)
      CALL dbcsr_finalize(dbcsr_mat%matrix)

      ALLOCATE (r_offset(nblkrows_total), c_offset(nblkcols_total))

      r_offset(1) = 1
      DO row = 2, nblkrows_total
         r_offset(row) = r_offset(row - 1) + row_blk_size(row - 1)
      END DO
      c_offset(1) = 1
      DO col = 2, nblkcols_total
         c_offset(col) = c_offset(col - 1) + col_blk_size(col - 1)
      END DO

      CALL dbcsr_iterator_start(iter, template_dbcsr%matrix)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, irow, icol, block)
         CALL dbcsr_put_block(dbcsr_mat%matrix, irow, icol, &
                              dense_mat(r_offset(irow):r_offset(irow) + SIZE(block, 1) - 1, &
                                        c_offset(icol):c_offset(icol) + SIZE(block, 2) - 1))
      END DO
      CALL dbcsr_iterator_stop(iter)

      DEALLOCATE (r_offset, c_offset)

   END FUNCTION dense_to_dbcsr

! **************************************************************************************************
!> \brief ...
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_gauxc_functional(xc_section) RESULT(gauxc_functional_section)
      TYPE(section_vals_type), INTENT(in), POINTER       :: xc_section
      TYPE(section_vals_type), POINTER                   :: gauxc_functional_section

      INTEGER                                            :: ifun
      TYPE(section_vals_type), POINTER                   :: functionals, xc_fun

      NULLIFY (gauxc_functional_section)

      functionals => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      IF (.NOT. ASSOCIATED(functionals)) THEN
         CPABORT("XC_FUNCTIONAL section not found")
      END IF

      ifun = 0
      DO
         ifun = ifun + 1
         xc_fun => section_vals_get_subs_vals2(functionals, i_section=ifun)
         IF (.NOT. ASSOCIATED(xc_fun)) EXIT
         IF (xc_fun%section%name /= "SKALA" .OR. ifun > 1) THEN
            CPABORT("SKALA functional is mutually exclusive with any other functional.")
         END IF
         gauxc_functional_section => xc_fun
      END DO

      IF (.NOT. ASSOCIATED(gauxc_functional_section)) THEN
         CPABORT("No XC functional found in XC_FUNCTIONAL section")
      END IF
   END FUNCTION get_gauxc_functional

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param xc_section ...
! **************************************************************************************************
   SUBROUTINE apply_gauxc(qs_env, xc_section)
      TYPE(qs_environment_type), INTENT(in), POINTER     :: qs_env
      TYPE(section_vals_type), INTENT(in), POINTER       :: xc_section

      CHARACTER(len=*), PARAMETER :: abort_message = &
         "GauXC does not support dispersion corrections. "// &
         "Please disable the VDW_POTENTIAL section or use the legacy XC (USE_GAUXC = .FALSE.)."

      CHARACTER(len=default_string_length)               :: xc_fun_name
      INTEGER                                            :: img, ispin, natom, nimages, nspins
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: density_scalar, density_zeta
      TYPE(cp_gauxc_basisset_type)                       :: gauxc_basis
      TYPE(cp_gauxc_grid_type)                           :: gauxc_grid_result
      TYPE(cp_gauxc_integrator_type)                     :: gauxc_integrator_result
      TYPE(cp_gauxc_molecule_type)                       :: gauxc_mol
      TYPE(cp_gauxc_status_type)                         :: gauxc_status
      TYPE(cp_gauxc_xc_type)                             :: gauxc_xc_result
      TYPE(dbcsr_p_type)                                 :: vxc_zeta_dbcsr
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_vxc
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: gauxc_functional_section

      NULLIFY ( &
         dft_control, &
         energy, &
         ks_env, &
         matrix_vxc, &
         particle_set, &
         qs_kind_set, &
         rho, &
         rho_ao)

      CALL get_qs_env( &
         qs_env, &
         dft_control=dft_control, &
         energy=energy, &
         ks_env=ks_env, &
         matrix_vxc=matrix_vxc, &
         natom=natom, &
         particle_set=particle_set, &
         qs_kind_set=qs_kind_set, &
         rho=rho)

      CALL qs_rho_get( &
         rho, &
         rho_ao_kp=rho_ao)

      nimages = dft_control%nimages
      nspins = dft_control%nspins

      IF (ASSOCIATED(qs_env%dispersion_env)) THEN
         IF (qs_env%dispersion_env%type /= xc_vdw_fun_none) THEN
            CPABORT(abort_message)
         END IF
      END IF
      NULLIFY (vxc_zeta_dbcsr%matrix)

      gauxc_functional_section => get_gauxc_functional(xc_section)
      CALL section_vals_val_get(gauxc_functional_section, "FUNCTIONAL", c_val=xc_fun_name)

      gauxc_mol = gauxc_create_molecule( &
                  particle_set, &
                  natom, &
                  gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      gauxc_basis = gauxc_create_basisset( &
                    qs_kind_set, &
                    particle_set, &
                    natom, &
                    gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      gauxc_grid_result = gauxc_create_grid( &
                          gauxc_mol, &
                          gauxc_basis, &
                          gauxc_functional_section, &
                          gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      gauxc_integrator_result = gauxc_create_integrator( &
                                TRIM(xc_fun_name), &
                                gauxc_grid_result, &
                                gauxc_functional_section, &
                                nspins, &
                                gauxc_status)
      CALL gauxc_check_status(gauxc_status)

      IF (qs_env%run_rtp) THEN
         CPABORT("GAUXC XC energy currently does not support real-time propagation")
      END IF

      DO img = 1, nimages
         IF (nspins == 1) THEN
            CALL dbcsr_to_dense(rho_ao(1, img), density_scalar)
            gauxc_xc_result = gauxc_compute_xc( &
                              gauxc_integrator_result, &
                              density_scalar, &
                              nspins=nspins, &
                              status=gauxc_status)
            CALL gauxc_check_status(gauxc_status)
         ELSE
            CPASSERT(nspins == 2)
            CALL dbcsr_to_dense(rho_ao(1, img), density_scalar)
            CALL dbcsr_to_dense(rho_ao(2, img), density_zeta)
            ! Do NOT reorder the following lines!
            density_scalar(:, :) = density_scalar(:, :) + density_zeta(:, :)
            ! Factor two because the next line is evaluated after the above line.
            ! We need to subtract density_zeta once to undo the above line and
            ! a second time because that is what UKS requires.
            ! This style lowers memory footprint.
            density_zeta(:, :) = density_scalar(:, :) - 2.0_dp*density_zeta(:, :)
            gauxc_xc_result = gauxc_compute_xc( &
                              gauxc_integrator_result, &
                              density_scalar, &
                              density_zeta, &
                              nspins, &
                              gauxc_status)
            CALL gauxc_check_status(gauxc_status)
         END IF

         energy%exc = gauxc_xc_result%exc

         IF (img == 1) THEN
            IF (ASSOCIATED(matrix_vxc)) CALL dbcsr_deallocate_matrix_set(matrix_vxc)
            CALL dbcsr_allocate_matrix_set(matrix_vxc, nspins)
         END IF

         IF (nspins == 1) THEN
            matrix_vxc(1) = dense_to_dbcsr(gauxc_xc_result%vxc_scalar, rho_ao(1, img))
         ELSE
            CPASSERT(nspins == 2)
            vxc_zeta_dbcsr = dense_to_dbcsr(gauxc_xc_result%vxc_zeta, rho_ao(1, img))
            DO ispin = 1, 2
               matrix_vxc(ispin) = dense_to_dbcsr(gauxc_xc_result%vxc_scalar, rho_ao(ispin, img))
               CALL dbcsr_add(matrix_vxc(ispin)%matrix, vxc_zeta_dbcsr%matrix, 1.0_dp, 1.0_dp - (ispin - 1)*2.0_dp)
            END DO
            CALL dbcsr_release(vxc_zeta_dbcsr%matrix)
            DEALLOCATE (vxc_zeta_dbcsr%matrix)
         END IF

         IF (img == 1) THEN
            DO ispin = 1, nspins
               CALL dbcsr_finalize(matrix_vxc(ispin)%matrix)
            END DO
         END IF
      END DO

      CALL set_ks_env(ks_env, matrix_vxc=matrix_vxc)

      CALL gauxc_destroy_integrator(gauxc_integrator_result, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      CALL gauxc_destroy_grid(gauxc_grid_result, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      CALL gauxc_destroy_basisset(gauxc_basis, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      CALL gauxc_destroy_molecule(gauxc_mol, gauxc_status)
      CALL gauxc_check_status(gauxc_status)

      DEALLOCATE (density_scalar)
      IF (ALLOCATED(density_zeta)) DEALLOCATE (density_zeta)

   END SUBROUTINE apply_gauxc

END MODULE xc_gauxc_functional
