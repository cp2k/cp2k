!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE xc_gauxc_functional
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_create,&
                                              dbcsr_finalize,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_release
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE external_potential_types,        ONLY: gth_potential_type,&
                                              sgp_potential_type
   USE input_constants,                 ONLY: xc_vdw_fun_nonloc
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_get_subs_vals2,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE iso_c_binding,                   ONLY: c_double
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_comm_self,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE string_utilities,                ONLY: uppercase
   USE xc_gauxc_interface,              ONLY: &
        cp_gauxc_basisset_type, cp_gauxc_grid_type, cp_gauxc_integrator_type, &
        cp_gauxc_molecule_type, cp_gauxc_status_type, cp_gauxc_xc_gradient_type, cp_gauxc_xc_type, &
        gauxc_check_status, gauxc_compute_xc, gauxc_compute_xc_gradient, gauxc_create_basisset, &
        gauxc_create_grid, gauxc_create_integrator, gauxc_create_molecule, gauxc_destroy_basisset, &
        gauxc_destroy_grid, gauxc_destroy_integrator, gauxc_destroy_molecule, &
        gauxc_write_basisset_hdf5, gauxc_write_molecule_hdf5
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'xc_gauxc_functional'

   PUBLIC :: apply_gauxc, skala_info, xc_section_uses_gauxc

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param dbcsr_mat ...
!> \param dense_mat ...
! **************************************************************************************************
   SUBROUTINE dbcsr_to_dense(dbcsr_mat, dense_mat)
      USE cp_dbcsr_api, ONLY: dbcsr_get_info, dbcsr_get_matrix_type, dbcsr_iterator_type, &
                              dbcsr_iterator_start, dbcsr_iterator_next_block, dbcsr_iterator_blocks_left, &
                              dbcsr_iterator_stop, dbcsr_type_antisymmetric, dbcsr_type_symmetric
      TYPE(dbcsr_p_type), INTENT(IN)                     :: dbcsr_mat
      REAL(c_double), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: dense_mat

      CHARACTER                                          :: matrix_type
      INTEGER                                            :: col, col_end, col_start, icol, irow, &
                                                            nblkcols_total, nblkrows_total, ncols, &
                                                            nrows, row, row_end, row_start
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: c_offset, r_offset
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      LOGICAL                                            :: transposed
      REAL(c_double), POINTER                            :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_get_info(dbcsr_mat%matrix, &
                          row_blk_size=row_blk_size, &
                          col_blk_size=col_blk_size, &
                          nblkrows_total=nblkrows_total, &
                          nblkcols_total=nblkcols_total, &
                          nfullrows_total=nrows, &
                          nfullcols_total=ncols)
      matrix_type = dbcsr_get_matrix_type(dbcsr_mat%matrix)

      IF (.NOT. ALLOCATED(dense_mat)) THEN
         ALLOCATE (dense_mat(nrows, ncols))
      ELSE IF (.NOT. ALL(SHAPE(dense_mat) == [nrows, ncols])) THEN
         DEALLOCATE (dense_mat)
         ALLOCATE (dense_mat(nrows, ncols))
      ELSE
         CPASSERT(ALL(SHAPE(dense_mat) == [nrows, ncols]))
      END IF
      dense_mat = 0._dp

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
         CALL dbcsr_iterator_next_block(iter, irow, icol, block, transposed=transposed)
         row_start = r_offset(irow)
         row_end = row_start + row_blk_size(irow) - 1
         col_start = c_offset(icol)
         col_end = col_start + col_blk_size(icol) - 1
         IF (transposed) THEN
            dense_mat(row_start:row_end, col_start:col_end) = TRANSPOSE(block)
            IF (irow /= icol) THEN
               IF (matrix_type == dbcsr_type_symmetric) THEN
                  dense_mat(col_start:col_end, row_start:row_end) = block
               ELSE IF (matrix_type == dbcsr_type_antisymmetric) THEN
                  dense_mat(col_start:col_end, row_start:row_end) = -block
               END IF
            END IF
         ELSE
            dense_mat(row_start:row_end, col_start:col_end) = block
            IF (irow /= icol) THEN
               IF (matrix_type == dbcsr_type_symmetric) THEN
                  dense_mat(col_start:col_end, row_start:row_end) = TRANSPOSE(block)
               ELSE IF (matrix_type == dbcsr_type_antisymmetric) THEN
                  dense_mat(col_start:col_end, row_start:row_end) = -TRANSPOSE(block)
               END IF
            END IF
         END IF
      END DO
      CALL dbcsr_iterator_stop(iter)

      DEALLOCATE (r_offset, c_offset)

   END SUBROUTINE dbcsr_to_dense

! ******, ***********************************************************************************
!> \brief Convert a dense symmetric matrix to a DBCSR matrix with full upper block structure.
!>        This creates all upper-triangular blocks, not just those present in a template.
!>        This is needed because GauXC computes VXC for the full dense density matrix.
!> \param dense_mat Input dense matrix
!> \param template_dbcsr Template DBCSR matrix for distribution and block sizes
!> \return dbcsr_mat Output DBCSR matrix with full upper block structure
! **************************************************************************************************
   FUNCTION dense_to_dbcsr(dense_mat, template_dbcsr) RESULT(dbcsr_mat)
      USE cp_dbcsr_api, ONLY: &
         dbcsr_create, &
         dbcsr_distribution_get, &
         dbcsr_distribution_type, &
         dbcsr_finalize, &
         dbcsr_get_info, &
         dbcsr_get_stored_coordinates, &
         dbcsr_init_p, &
         dbcsr_put_block, &
         dbcsr_release, &
         dbcsr_type_symmetric, &
         dbcsr_work_create
      REAL(c_double), DIMENSION(:, :), INTENT(IN)        :: dense_mat
      TYPE(dbcsr_p_type), INTENT(IN)                     :: template_dbcsr
      TYPE(dbcsr_p_type)                                 :: dbcsr_mat

      INTEGER                                            :: col, icol, irow, mynode, nblkcols_total, &
                                                            nblkrows_total, ncols, nrows, owner, &
                                                            row
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: c_offset, r_offset
      INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, row_blk_size
      TYPE(dbcsr_distribution_type)                      :: dist

      CALL dbcsr_get_info(template_dbcsr%matrix, &
                          row_blk_size=row_blk_size, &
                          col_blk_size=col_blk_size, &
                          nblkrows_total=nblkrows_total, &
                          nblkcols_total=nblkcols_total, &
                          nfullrows_total=nrows, &
                          nfullcols_total=ncols, &
                          distribution=dist)
      CALL dbcsr_distribution_get(dist, mynode=mynode)

      CPASSERT(nrows == SIZE(dense_mat, 1))
      CPASSERT(ncols == SIZE(dense_mat, 2))

      CALL dbcsr_init_p(dbcsr_mat%matrix)
      CALL dbcsr_create(dbcsr_mat%matrix, &
                        template=template_dbcsr%matrix, &
                        name="VXC from GauXC (dense)", &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_work_create(dbcsr_mat%matrix, work_mutable=.TRUE.)

      ALLOCATE (r_offset(nblkrows_total), c_offset(nblkcols_total))

      r_offset(1) = 1
      DO row = 2, nblkrows_total
         r_offset(row) = r_offset(row - 1) + row_blk_size(row - 1)
      END DO
      c_offset(1) = 1
      DO col = 2, nblkcols_total
         c_offset(col) = c_offset(col - 1) + col_blk_size(col - 1)
      END DO

      DO irow = 1, nblkrows_total
         DO icol = 1, nblkcols_total
            IF (irow > icol) CYCLE
            CALL dbcsr_get_stored_coordinates(dbcsr_mat%matrix, irow, icol, owner)
            IF (owner /= mynode) CYCLE
            CALL dbcsr_put_block(dbcsr_mat%matrix, irow, icol, &
                                 0.5_dp*( &
                                 dense_mat(r_offset(irow):r_offset(irow) + row_blk_size(irow) - 1, &
                                           c_offset(icol):c_offset(icol) + col_blk_size(icol) - 1) + &
                                 TRANSPOSE(dense_mat(r_offset(icol):r_offset(icol) + row_blk_size(icol) - 1, &
                                                     c_offset(irow):c_offset(irow) + col_blk_size(irow) - 1))))
         END DO
      END DO

      CALL dbcsr_finalize(dbcsr_mat%matrix)

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
         IF (xc_fun%section%name /= "GAUXC" .OR. ifun > 1) THEN
            CPABORT("GauXC functionals are mutually exclusive with any other functional.")
         END IF
         gauxc_functional_section => xc_fun
      END DO

      IF (.NOT. ASSOCIATED(gauxc_functional_section)) THEN
         CPABORT("No XC functional found in XC_FUNCTIONAL section")
      END IF
   END FUNCTION get_gauxc_functional

! **************************************************************************************************
!> \brief ...
!> \param xc_section ...
!> \return ...
! **************************************************************************************************
   FUNCTION xc_section_uses_gauxc(xc_section) RESULT(uses_gauxc)
      TYPE(section_vals_type), INTENT(in), POINTER       :: xc_section
      LOGICAL                                            :: uses_gauxc

      INTEGER                                            :: ifun
      TYPE(section_vals_type), POINTER                   :: functionals, xc_fun

      uses_gauxc = .FALSE.
      IF (.NOT. ASSOCIATED(xc_section)) RETURN

      functionals => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      IF (.NOT. ASSOCIATED(functionals)) RETURN

      ifun = 0
      DO
         ifun = ifun + 1
         xc_fun => section_vals_get_subs_vals2(functionals, i_section=ifun)
         IF (.NOT. ASSOCIATED(xc_fun)) EXIT
         IF (xc_fun%section%name == "GAUXC") THEN
            uses_gauxc = .TRUE.
            EXIT
         END IF
      END DO

   END FUNCTION xc_section_uses_gauxc

! **************************************************************************************************
!> \brief Reject unsupported pseudopotential variants in GauXC GAPW mode
!> \param qs_kind_set ...
! **************************************************************************************************
   SUBROUTINE ensure_gauxc_gapw_all_electron(qs_kind_set)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      INTEGER                                            :: ikind
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(sgp_potential_type), POINTER                  :: sgp_potential

      CPASSERT(ASSOCIATED(qs_kind_set))

      DO ikind = 1, SIZE(qs_kind_set)
         NULLIFY (gth_potential, sgp_potential)
         CALL get_qs_kind(qs_kind_set(ikind), &
                          gth_potential=gth_potential, &
                          sgp_potential=sgp_potential)
         IF (ASSOCIATED(gth_potential) .OR. ASSOCIATED(sgp_potential)) THEN
            CALL cp_abort(__LOCATION__, &
                          "GauXC with METHOD GAPW currently supports all-electron potentials only. "// &
                          "Use POTENTIAL ALL for GAPW validation or METHOD GPW with pseudopotentials.")
         END IF
      END DO

   END SUBROUTINE ensure_gauxc_gapw_all_electron

! **************************************************************************************************
!> \brief adds a replicated GauXC energy gradient to the local CP2K force accumulator
!> \param exc_grad ...
!> \param force ...
!> \param atomic_kind_set ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE add_gauxc_gradient_to_force(exc_grad, force, atomic_kind_set, para_env)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: exc_grad
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: ia, iatom, ikind, natom_kind
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind

      CPASSERT(ASSOCIATED(force))
      CPASSERT(ASSOCIATED(atomic_kind_set))

      IF (para_env%mepos /= 0) RETURN

      DO ikind = 1, SIZE(atomic_kind_set, 1)
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind=atomic_kind, natom=natom_kind)
         DO ia = 1, natom_kind
            iatom = atomic_kind%atom_list(ia)
            force(ikind)%rho_elec(:, ia) = force(ikind)%rho_elec(:, ia) + &
                                           exc_grad(3*iatom - 2:3*iatom)
         END DO
      END DO

   END SUBROUTINE add_gauxc_gradient_to_force

! **************************************************************************************************
!> \brief compute a GauXC XC energy for diagnostic finite differences
!> \param particle_set_eval ...
!> \param qs_kind_set ...
!> \param density_scalar ...
!> \param nspins ...
!> \param model_name ...
!> \param xc_fun_name ...
!> \param grid_type ...
!> \param radial_quadrature ...
!> \param pruning_scheme ...
!> \param lb_exec_space ...
!> \param int_exec_space ...
!> \param batch_size ...
!> \param exc ...
!> \param density_zeta ...
! **************************************************************************************************
   SUBROUTINE gauxc_xc_energy_for_particles( &
      particle_set_eval, qs_kind_set, density_scalar, nspins, model_name, &
      xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
      int_exec_space, batch_size, exc, density_zeta)
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set_eval
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: density_scalar
      INTEGER, INTENT(IN)                                :: nspins
      CHARACTER(len=*), INTENT(IN)                       :: model_name, xc_fun_name, grid_type, &
                                                            radial_quadrature, pruning_scheme, &
                                                            lb_exec_space, int_exec_space
      INTEGER, INTENT(IN)                                :: batch_size
      REAL(KIND=dp), INTENT(OUT)                         :: exc
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: density_zeta

      TYPE(cp_gauxc_basisset_type)                       :: gauxc_basis_fd
      TYPE(cp_gauxc_grid_type)                           :: gauxc_grid_fd
      TYPE(cp_gauxc_integrator_type)                     :: gauxc_integrator_fd
      TYPE(cp_gauxc_molecule_type)                       :: gauxc_mol_fd
      TYPE(cp_gauxc_status_type)                         :: gauxc_status
      TYPE(cp_gauxc_xc_type)                             :: gauxc_xc_result

      gauxc_mol_fd = gauxc_create_molecule(particle_set_eval, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      gauxc_basis_fd = gauxc_create_basisset(qs_kind_set, particle_set_eval, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      gauxc_grid_fd = gauxc_create_grid( &
                      gauxc_mol_fd, &
                      gauxc_basis_fd, &
                      grid_type, &
                      radial_quadrature, &
                      pruning_scheme, &
                      lb_exec_space, &
                      batch_size, &
                      gauxc_status, &
                      mpi_comm=mp_comm_self%get_handle())
      CALL gauxc_check_status(gauxc_status)
      gauxc_integrator_fd = gauxc_create_integrator( &
                            TRIM(xc_fun_name), &
                            gauxc_grid_fd, &
                            int_exec_space, &
                            nspins, &
                            gauxc_status)
      CALL gauxc_check_status(gauxc_status)

      IF (nspins == 1) THEN
         gauxc_xc_result = gauxc_compute_xc( &
                           gauxc_integrator_fd, &
                           density_scalar, &
                           nspins=nspins, &
                           status=gauxc_status, &
                           model=TRIM(model_name))
      ELSE
         CPASSERT(nspins == 2)
         CPASSERT(PRESENT(density_zeta))
         gauxc_xc_result = gauxc_compute_xc( &
                           gauxc_integrator_fd, &
                           density_scalar, &
                           density_zeta, &
                           nspins, &
                           gauxc_status, &
                           model=TRIM(model_name))
      END IF
      CALL gauxc_check_status(gauxc_status)
      exc = gauxc_xc_result%exc

      IF (ALLOCATED(gauxc_xc_result%vxc_scalar)) DEALLOCATE (gauxc_xc_result%vxc_scalar)
      IF (ALLOCATED(gauxc_xc_result%vxc_zeta)) DEALLOCATE (gauxc_xc_result%vxc_zeta)

      CALL gauxc_destroy_integrator(gauxc_integrator_fd, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      CALL gauxc_destroy_grid(gauxc_grid_fd, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      CALL gauxc_destroy_basisset(gauxc_basis_fd, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      CALL gauxc_destroy_molecule(gauxc_mol_fd, gauxc_status)
      CALL gauxc_check_status(gauxc_status)

   END SUBROUTINE gauxc_xc_energy_for_particles

! **************************************************************************************************
!> \brief compute a finite-difference GauXC XC nuclear gradient at fixed density
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param density_scalar ...
!> \param nspins ...
!> \param model_name ...
!> \param xc_fun_name ...
!> \param grid_type ...
!> \param radial_quadrature ...
!> \param pruning_scheme ...
!> \param lb_exec_space ...
!> \param int_exec_space ...
!> \param batch_size ...
!> \param dx ...
!> \param para_env ...
!> \param exc_grad ...
!> \param density_zeta ...
! **************************************************************************************************
   SUBROUTINE gauxc_xc_gradient_fd( &
      particle_set, qs_kind_set, density_scalar, nspins, model_name, &
      xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
      int_exec_space, batch_size, dx, para_env, exc_grad, density_zeta)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: density_scalar
      INTEGER, INTENT(IN)                                :: nspins
      CHARACTER(len=*), INTENT(IN)                       :: model_name, xc_fun_name, grid_type, &
                                                            radial_quadrature, pruning_scheme, &
                                                            lb_exec_space, int_exec_space
      INTEGER, INTENT(IN)                                :: batch_size
      REAL(KIND=dp), INTENT(IN)                          :: dx
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: exc_grad
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: density_zeta

      INTEGER                                            :: iatom, idir
      REAL(KIND=dp)                                      :: xc_minus, xc_plus
      TYPE(particle_type), ALLOCATABLE, DIMENSION(:)     :: particle_set_minus, particle_set_plus

      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(dx > 0.0_dp)

      ALLOCATE (exc_grad(3*SIZE(particle_set)))
      exc_grad = 0.0_dp

      IF (para_env%mepos == 0) THEN
         ALLOCATE (particle_set_minus(SIZE(particle_set)), particle_set_plus(SIZE(particle_set)))

         DO iatom = 1, SIZE(particle_set)
            DO idir = 1, 3
               particle_set_minus = particle_set
               particle_set_plus = particle_set
               particle_set_minus(iatom)%r(idir) = particle_set_minus(iatom)%r(idir) - dx
               particle_set_plus(iatom)%r(idir) = particle_set_plus(iatom)%r(idir) + dx
               IF (PRESENT(density_zeta)) THEN
                  CALL gauxc_xc_energy_for_particles( &
                     particle_set_plus, qs_kind_set, density_scalar, nspins, model_name, &
                     xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
                     int_exec_space, batch_size, xc_plus, density_zeta=density_zeta)
                  CALL gauxc_xc_energy_for_particles( &
                     particle_set_minus, qs_kind_set, density_scalar, nspins, model_name, &
                     xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
                     int_exec_space, batch_size, xc_minus, density_zeta=density_zeta)
               ELSE
                  CALL gauxc_xc_energy_for_particles( &
                     particle_set_plus, qs_kind_set, density_scalar, nspins, model_name, &
                     xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
                     int_exec_space, batch_size, xc_plus)
                  CALL gauxc_xc_energy_for_particles( &
                     particle_set_minus, qs_kind_set, density_scalar, nspins, model_name, &
                     xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
                     int_exec_space, batch_size, xc_minus)
               END IF
               exc_grad(3*iatom - 3 + idir) = (xc_plus - xc_minus)/(2.0_dp*dx)
            END DO
         END DO

         DEALLOCATE (particle_set_minus, particle_set_plus)
      END IF

      CALL para_env%bcast(exc_grad, 0)

   END SUBROUTINE gauxc_xc_gradient_fd

! **************************************************************************************************
!> \brief finite-difference check of the molecular GauXC XC virial diagnostic
!> \param exc_grad ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param density_scalar ...
!> \param nspins ...
!> \param model_name ...
!> \param xc_fun_name ...
!> \param grid_type ...
!> \param radial_quadrature ...
!> \param pruning_scheme ...
!> \param lb_exec_space ...
!> \param int_exec_space ...
!> \param batch_size ...
!> \param dx ...
!> \param para_env ...
!> \param density_zeta ...
! **************************************************************************************************
   SUBROUTINE debug_gauxc_molecular_virial( &
      exc_grad, particle_set, qs_kind_set, density_scalar, nspins, model_name, &
      xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
      int_exec_space, batch_size, dx, para_env, density_zeta)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: exc_grad
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: density_scalar
      INTEGER, INTENT(IN)                                :: nspins
      CHARACTER(len=*), INTENT(IN)                       :: model_name, xc_fun_name, grid_type, &
                                                            radial_quadrature, pruning_scheme, &
                                                            lb_exec_space, int_exec_space
      INTEGER, INTENT(IN)                                :: batch_size
      REAL(KIND=dp), INTENT(IN)                          :: dx
      TYPE(mp_para_env_type), POINTER                    :: para_env
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: density_zeta

      INTEGER                                            :: iatom, iw
      REAL(KIND=dp)                                      :: analytic_trace, diff_trace, &
                                                            numerical_trace, xc_minus, xc_plus
      REAL(KIND=dp), DIMENSION(3)                        :: center, displacement, grad
      TYPE(particle_type), ALLOCATABLE, DIMENSION(:)     :: particle_set_minus, particle_set_plus

      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(SIZE(exc_grad) == 3*SIZE(particle_set))

      IF (para_env%mepos /= 0) RETURN

      center = 0.0_dp
      DO iatom = 1, SIZE(particle_set)
         center = center + particle_set(iatom)%r
      END DO
      center = center/REAL(SIZE(particle_set), dp)

      ALLOCATE (particle_set_minus(SIZE(particle_set)), particle_set_plus(SIZE(particle_set)))
      particle_set_minus = particle_set
      particle_set_plus = particle_set

      analytic_trace = 0.0_dp
      DO iatom = 1, SIZE(particle_set)
         grad = exc_grad(3*iatom - 2:3*iatom)
         displacement = particle_set(iatom)%r - center
         analytic_trace = analytic_trace + DOT_PRODUCT(grad, displacement)
         particle_set_minus(iatom)%r = center + (1.0_dp - dx)*displacement
         particle_set_plus(iatom)%r = center + (1.0_dp + dx)*displacement
      END DO
      analytic_trace = analytic_trace/3.0_dp

      IF (PRESENT(density_zeta)) THEN
         CALL gauxc_xc_energy_for_particles( &
            particle_set_plus, qs_kind_set, density_scalar, nspins, model_name, &
            xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
            int_exec_space, batch_size, xc_plus, density_zeta=density_zeta)
         CALL gauxc_xc_energy_for_particles( &
            particle_set_minus, qs_kind_set, density_scalar, nspins, model_name, &
            xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
            int_exec_space, batch_size, xc_minus, density_zeta=density_zeta)
      ELSE
         CALL gauxc_xc_energy_for_particles( &
            particle_set_plus, qs_kind_set, density_scalar, nspins, model_name, &
            xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
            int_exec_space, batch_size, xc_plus)
         CALL gauxc_xc_energy_for_particles( &
            particle_set_minus, qs_kind_set, density_scalar, nspins, model_name, &
            xc_fun_name, grid_type, radial_quadrature, pruning_scheme, lb_exec_space, &
            int_exec_space, batch_size, xc_minus)
      END IF

      numerical_trace = (xc_plus - xc_minus)/(2.0_dp*dx)/3.0_dp
      diff_trace = analytic_trace - numerical_trace

      iw = cp_logger_get_default_io_unit()
      IF (iw > 0) THEN
         WRITE (UNIT=iw, FMT="(/,T2,A,1X,ES11.4)") &
            "GAUXC| Molecular XC virial finite-difference dx", dx
         WRITE (UNIT=iw, FMT="(T2,A,3(1X,ES19.11))") &
            "GAUXC| Molecular XC virial FD 1/3 Trace", &
            analytic_trace, numerical_trace, diff_trace
      END IF

      DEALLOCATE (particle_set_minus, particle_set_plus)

   END SUBROUTINE debug_gauxc_molecular_virial

! **************************************************************************************************
!> \brief prints a force-based molecular XC virial diagnostic from GauXC gradients
!> \param exc_grad ...
!> \param particle_set ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE print_gauxc_molecular_virial(exc_grad, particle_set, para_env)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: exc_grad
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CHARACTER(len=1), DIMENSION(3), PARAMETER          :: label = ["x", "y", "z"]

      INTEGER                                            :: i, iatom, iw, j
      REAL(KIND=dp), DIMENSION(3)                        :: center, displacement, grad, grad_sum
      REAL(KIND=dp), DIMENSION(3, 3)                     :: molecular_virial

      CPASSERT(ASSOCIATED(particle_set))
      CPASSERT(SIZE(exc_grad) == 3*SIZE(particle_set))

      IF (para_env%mepos /= 0) RETURN

      center = 0.0_dp
      DO iatom = 1, SIZE(particle_set)
         center = center + particle_set(iatom)%r
      END DO
      center = center/REAL(SIZE(particle_set), dp)

      grad_sum = 0.0_dp
      molecular_virial = 0.0_dp
      DO iatom = 1, SIZE(particle_set)
         grad = exc_grad(3*iatom - 2:3*iatom)
         displacement = particle_set(iatom)%r - center
         grad_sum = grad_sum + grad
         DO i = 1, 3
            DO j = 1, 3
               molecular_virial(i, j) = molecular_virial(i, j) + grad(i)*displacement(j)
            END DO
         END DO
      END DO

      iw = cp_logger_get_default_io_unit()
      IF (iw <= 0) RETURN

      WRITE (UNIT=iw, FMT="(/,T2,A)") &
         "GAUXC| Molecular XC gradient virial diagnostic [a.u.]"
      WRITE (UNIT=iw, FMT="(T2,A,T20,A,T40,A,T60,A)") "GAUXC|", "x", "y", "z"
      DO i = 1, 3
         WRITE (UNIT=iw, FMT="(T2,A,1X,A1,3(1X,ES19.11))") &
            "GAUXC|", label(i), molecular_virial(i, :)
      END DO
      WRITE (UNIT=iw, FMT="(T2,A,1X,ES19.11)") &
         "GAUXC| Molecular XC gradient virial 1/3 Trace", &
         (molecular_virial(1, 1) + molecular_virial(2, 2) + molecular_virial(3, 3))/3.0_dp
      WRITE (UNIT=iw, FMT="(T2,A,3(1X,ES19.11))") &
         "GAUXC| Molecular XC gradient sum", grad_sum
      WRITE (UNIT=iw, FMT="(T2,A)") &
         "GAUXC| Diagnostic only; this is not an analytical periodic stress tensor."

   END SUBROUTINE print_gauxc_molecular_virial

! **************************************************************************************************
!> \brief Return information about the Skala functional
!> \param functional section containing the SKALA subsection
!> \param lsd if you are using lsd or lda
!> \param reference the reference to the article where the functional is explained
!> \param shortform the short definition of the functional
!> \param needs the flags corresponding to the inputs needed by this
!>        functional are set to true (the flags not needed aren't touched)
!> \param max_deriv the maximal derivative available
! **************************************************************************************************
   SUBROUTINE skala_info(functional, lsd, reference, shortform, needs, max_deriv)
      TYPE(section_vals_type), POINTER                   :: functional
      LOGICAL, INTENT(in)                                :: lsd
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL            :: reference, shortform
      TYPE(xc_rho_cflags_type), INTENT(inout), OPTIONAL  :: needs
      INTEGER, INTENT(out), OPTIONAL                     :: max_deriv

      CHARACTER(len=default_path_length)                 :: model_key, model_name
      CHARACTER(len=default_string_length)               :: xc_fun_name

      CALL section_vals_val_get(functional, "FUNCTIONAL", c_val=xc_fun_name)
      CALL section_vals_val_get(functional, "MODEL", c_val=model_name)
      model_key = ADJUSTL(model_name)
      CALL uppercase(model_key)

      IF (PRESENT(reference)) THEN
         IF (TRIM(model_key) == "NONE" .OR. TRIM(model_key) == "") THEN
            reference = "Functional computed by GauXC (underlying: "//TRIM(xc_fun_name)//")"
         ELSE
            reference = "Functional computed by GauXC OneDFT model "//TRIM(model_name)
         END IF
      END IF
      IF (PRESENT(shortform)) THEN
         IF (TRIM(model_key) == "NONE" .OR. TRIM(model_key) == "") THEN
            shortform = "GAUXC ("//TRIM(xc_fun_name)//")"
         ELSE
            shortform = "GAUXC OneDFT"
         END IF
      END IF
      IF (PRESENT(needs)) THEN
         needs%rho = .TRUE.
         IF (lsd) THEN
            needs%rho_spin = .TRUE.
         END IF
      END IF
      IF (PRESENT(max_deriv)) max_deriv = 1

   END SUBROUTINE skala_info

! GauXC uses replicated dense density and VXC matrices. The DBCSR density matrix
! is distributed over MPI ranks, so apply_gauxc allreduces the dense copy before
! passing it to GauXC.

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param xc_section ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE apply_gauxc(qs_env, xc_section, calculate_forces)
      TYPE(qs_environment_type), INTENT(in), POINTER     :: qs_env
      TYPE(section_vals_type), INTENT(in), POINTER       :: xc_section
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(len=*), PARAMETER :: gapw_xc_abort_message = &
         "GauXC with METHOD GAPW_XC is not supported yet. "// &
         "The GAPW_XC one-center XC correction needs a dedicated GauXC design.", &
         nonlocal_vdw_abort_message = &
         "GauXC does not support non-local VDW_POTENTIAL corrections. "// &
         "Use an additive PAIR_POTENTIAL dispersion correction or disable GauXC."
      REAL(KIND=dp), PARAMETER :: gapw_fd_gradient_dx = 1.0E-4_dp

      CHARACTER(len=default_path_length)                 :: model_key, model_name, output_path
      CHARACTER(len=default_string_length) :: grid_key, grid_type, int_exec_space, lb_exec_space, &
         pruning_key, pruning_scheme, radial_quadrature, xc_fun_name
      INTEGER                                            :: batch_size, env_status, img, ispin, &
                                                            natom, nimages, nspins
      LOGICAL :: grid_explicit, hdf5_output, molecular_virial, molecular_virial_debug, &
         pruning_explicit, use_fd_gradient, use_gradient_self_runtime, use_onedft, &
         use_self_runtime, use_skala_model, write_hdf5_output
      REAL(KIND=dp)                                      :: molecular_virial_debug_dx
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: density_scalar, density_zeta
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_gauxc_basisset_type)                       :: gauxc_basis
      TYPE(cp_gauxc_grid_type)                           :: gauxc_gradient_grid_result, &
                                                            gauxc_grid_result
      TYPE(cp_gauxc_integrator_type) :: gauxc_gradient_integrator_result, gauxc_integrator_result
      TYPE(cp_gauxc_molecule_type)                       :: gauxc_mol
      TYPE(cp_gauxc_status_type)                         :: gauxc_status
      TYPE(cp_gauxc_xc_gradient_type)                    :: exc_grad
      TYPE(cp_gauxc_xc_type)                             :: gauxc_xc_result
      TYPE(dbcsr_p_type)                                 :: vxc_zeta_tmp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_vxc
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_use, rho_xc
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: gauxc_functional_section

      NULLIFY ( &
         atomic_kind_set, &
         dft_control, &
         energy, &
         force, &
         ks_env, &
         matrix_vxc, &
         para_env, &
         particle_set, &
         qs_kind_set, &
         rho, &
         rho_use, &
         rho_xc, &
         rho_ao, &
         scf_env)

      CALL get_qs_env( &
         qs_env, &
         dft_control=dft_control, &
         energy=energy, &
         ks_env=ks_env, &
         matrix_vxc=matrix_vxc, &
         natom=natom, &
         atomic_kind_set=atomic_kind_set, &
         force=force, &
         para_env=para_env, &
         particle_set=particle_set, &
         qs_kind_set=qs_kind_set, &
         rho=rho, &
         rho_xc=rho_xc, &
         scf_env=scf_env)

      IF (dft_control%qs_control%gapw_xc) THEN
         CPABORT(gapw_xc_abort_message)
      END IF
      IF (dft_control%qs_control%gapw) THEN
         CALL ensure_gauxc_gapw_all_electron(qs_kind_set)
      END IF
      CPASSERT(ASSOCIATED(rho))
      rho_use => rho
      CALL qs_rho_get( &
         rho_use, &
         rho_ao_kp=rho_ao)

      nimages = dft_control%nimages
      nspins = dft_control%nspins

      IF (ASSOCIATED(qs_env%dispersion_env)) THEN
         IF (qs_env%dispersion_env%type == xc_vdw_fun_nonloc) THEN
            CPABORT(nonlocal_vdw_abort_message)
         END IF
      END IF
      NULLIFY (vxc_zeta_tmp%matrix)

      gauxc_functional_section => get_gauxc_functional(xc_section)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "FUNCTIONAL", &
         c_val=xc_fun_name)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "MODEL", &
         c_val=model_name)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "GRID", &
         c_val=grid_type, &
         explicit=grid_explicit)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "RADIAL_QUADRATURE", &
         c_val=radial_quadrature)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "PRUNING_SCHEME", &
         c_val=pruning_scheme, &
         explicit=pruning_explicit)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "BATCH_SIZE", &
         i_val=batch_size)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "MOLECULAR_VIRIAL", &
         l_val=molecular_virial)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "MOLECULAR_VIRIAL_DEBUG", &
         l_val=molecular_virial_debug)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "MOLECULAR_VIRIAL_DEBUG_DX", &
         r_val=molecular_virial_debug_dx)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "LB_EXECUTION_SPACE", &
         c_val=lb_exec_space)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "INT_EXECUTION_SPACE", &
         c_val=int_exec_space)
      CALL section_vals_val_get( &
         gauxc_functional_section, &
         "OUTPUT_PATH", &
         c_val=output_path)

      model_key = ADJUSTL(model_name)
      CALL uppercase(model_key)
      use_onedft = (TRIM(model_key) /= "" .AND. TRIM(model_key) /= "NONE")
      use_skala_model = (INDEX(TRIM(model_key), "SKALA") > 0)
      IF (molecular_virial_debug) THEN
         IF (molecular_virial_debug_dx <= 0.0_dp) THEN
            CALL cp_abort(__LOCATION__, &
                          "GauXC MOLECULAR_VIRIAL_DEBUG_DX must be positive.")
         END IF
         molecular_virial = .TRUE.
      END IF
      IF (dft_control%qs_control%gapw_xc .AND. use_onedft) THEN
         CALL cp_abort(__LOCATION__, &
                       "GauXC OneDFT/SKALA with METHOD GAPW_XC is not implemented. "// &
                       "The GAPW_XC one-center XC correction must be evaluated by the same GauXC model.")
      END IF
      IF (use_onedft) THEN
         IF (.NOT. grid_explicit) grid_type = "SUPERFINE"
         IF (.NOT. pruning_explicit) pruning_scheme = "UNPRUNED"

         grid_key = ADJUSTL(grid_type)
         pruning_key = ADJUSTL(pruning_scheme)
         CALL uppercase(grid_key)
         CALL uppercase(pruning_key)
         IF (use_skala_model .AND. (calculate_forces .OR. molecular_virial) .AND. &
             (TRIM(grid_key) /= "SUPERFINE" .OR. TRIM(pruning_key) /= "UNPRUNED")) THEN
            CALL cp_warn( &
               __LOCATION__, &
               "GauXC OneDFT/SKALA nuclear gradients are sensitive to the GauXC molecular grid. "// &
               "Use GRID SUPERFINE and PRUNING_SCHEME UNPRUNED for quantitative force checks.")
         END IF
         IF (TRIM(model_key) == "SKALA") THEN
            CALL GET_ENVIRONMENT_VARIABLE("GAUXC_SKALA_MODEL", model_name, STATUS=env_status)
            IF (env_status /= 0 .OR. LEN_TRIM(model_name) == 0) THEN
               CPABORT("MODEL SKALA requires the GAUXC_SKALA_MODEL environment variable")
            END IF
         END IF
      END IF
      use_self_runtime = use_skala_model .AND. para_env%num_pe > 1
      use_gradient_self_runtime = calculate_forces .AND. use_onedft .AND. &
                                  para_env%num_pe > 1 .AND. .NOT. use_self_runtime

      gauxc_mol = gauxc_create_molecule( &
                  particle_set, &
                  gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      gauxc_basis = gauxc_create_basisset( &
                    qs_kind_set, &
                    particle_set, &
                    gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      hdf5_output = (TRIM(output_path) /= "")
      write_hdf5_output = hdf5_output .AND. para_env%mepos == 0
      IF (write_hdf5_output .AND. ASSOCIATED(scf_env)) THEN
         write_hdf5_output = scf_env%iter_count == 1
      END IF
      IF (write_hdf5_output) THEN
         CALL gauxc_write_molecule_hdf5( &
            gauxc_mol, &
            output_path, &
            "molecule.h5", &
            "molecule", &
            gauxc_status)
         CALL gauxc_check_status(gauxc_status)
         CALL gauxc_write_basisset_hdf5( &
            gauxc_basis, &
            output_path, &
            "basisset.h5", &
            "basisset", &
            gauxc_status)
         CALL gauxc_check_status(gauxc_status)
      END IF
      use_fd_gradient = dft_control%qs_control%gapw .AND. calculate_forces .AND. &
                        (use_skala_model .OR. gauxc_basis%max_l > 3)
      IF (use_fd_gradient .AND. para_env%mepos == 0) THEN
         CALL cp_warn( &
            __LOCATION__, &
            "Using finite-difference GauXC XC gradients for METHOD GAPW with SKALA or "// &
            "all-electron basis functions beyond f shells.  The upstream analytical GauXC "// &
            "gradient path is not yet reliable for this case.")
      END IF
      IF (use_self_runtime) THEN
         ! SKALA currently needs a replicated molecular runtime for reproducible
         ! open-shell densities across CP2K MPI ranks.
         gauxc_grid_result = gauxc_create_grid( &
                             gauxc_mol, &
                             gauxc_basis, &
                             grid_type, &
                             radial_quadrature, &
                             pruning_scheme, &
                             lb_exec_space, &
                             batch_size, &
                             gauxc_status, &
                             mpi_comm=mp_comm_self%get_handle())
      ELSE
         ! Use the QS force-evaluation communicator rather than GauXC's global
         ! runtime.  In mixed CDFT the diabatic states can be built on disjoint
         ! MPI subgroups; using the global communicator would make GauXC wait
         ! for ranks that are working on another state.
         gauxc_grid_result = gauxc_create_grid( &
                             gauxc_mol, &
                             gauxc_basis, &
                             grid_type, &
                             radial_quadrature, &
                             pruning_scheme, &
                             lb_exec_space, &
                             batch_size, &
                             gauxc_status, &
                             mpi_comm=para_env%get_handle())
      END IF
      CALL gauxc_check_status(gauxc_status)
      gauxc_integrator_result = gauxc_create_integrator( &
                                TRIM(xc_fun_name), &
                                gauxc_grid_result, &
                                int_exec_space, &
                                nspins, &
                                gauxc_status)
      CALL gauxc_check_status(gauxc_status)

      IF (use_gradient_self_runtime) THEN
         ! Upstream GauXC does not yet support OneDFT/SKALA nuclear gradients
         ! on an MPI runtime.  Keep the energy/VXC path on the normal MPI
         ! runtime and use an isolated runtime only for the replicated gradient.
         gauxc_gradient_grid_result = gauxc_create_grid( &
                                      gauxc_mol, &
                                      gauxc_basis, &
                                      grid_type, &
                                      radial_quadrature, &
                                      pruning_scheme, &
                                      lb_exec_space, &
                                      batch_size, &
                                      gauxc_status, &
                                      mpi_comm=mp_comm_self%get_handle())
         CALL gauxc_check_status(gauxc_status)
         gauxc_gradient_integrator_result = gauxc_create_integrator( &
                                            TRIM(xc_fun_name), &
                                            gauxc_gradient_grid_result, &
                                            int_exec_space, &
                                            nspins, &
                                            gauxc_status)
         CALL gauxc_check_status(gauxc_status)
      END IF

      IF (qs_env%run_rtp) THEN
         CPABORT("GAUXC XC energy currently does not support real-time propagation")
      END IF

      energy%exc = 0

      IF (ASSOCIATED(matrix_vxc)) CALL dbcsr_deallocate_matrix_set(matrix_vxc)
      CALL dbcsr_allocate_matrix_set(matrix_vxc, nspins)

      DO img = 1, nimages
         IF (img > 1) THEN
            CPABORT("UNIMPLEMENTED: Handling nimg>1 in k-point integration")
         END IF
         CALL dbcsr_to_dense(rho_ao(1, img), density_scalar)
         CALL para_env%sum(density_scalar)
         IF (nspins == 1) THEN
            gauxc_xc_result = gauxc_compute_xc( &
                              gauxc_integrator_result, &
                              density_scalar, &
                              nspins=nspins, &
                              status=gauxc_status, &
                              model=TRIM(model_name))
            CALL gauxc_check_status(gauxc_status)
            IF (calculate_forces) THEN
               IF (use_fd_gradient) THEN
                  CALL gauxc_xc_gradient_fd( &
                     particle_set, qs_kind_set, density_scalar, nspins, model_name, &
                     xc_fun_name, grid_type, radial_quadrature, pruning_scheme, &
                     lb_exec_space, int_exec_space, batch_size, gapw_fd_gradient_dx, &
                     para_env, exc_grad%exc_grad)
               ELSE IF (use_gradient_self_runtime) THEN
                  exc_grad = gauxc_compute_xc_gradient( &
                             gauxc_gradient_integrator_result, &
                             density_scalar, &
                             nspins=nspins, &
                             natom=natom, &
                             status=gauxc_status, &
                             model=TRIM(model_name))
               ELSE
                  exc_grad = gauxc_compute_xc_gradient( &
                             gauxc_integrator_result, &
                             density_scalar, &
                             nspins=nspins, &
                             natom=natom, &
                             status=gauxc_status, &
                             model=TRIM(model_name))
               END IF
               CALL gauxc_check_status(gauxc_status)
               CALL add_gauxc_gradient_to_force( &
                  exc_grad%exc_grad, &
                  force, &
                  atomic_kind_set, &
                  para_env)
               IF (molecular_virial) THEN
                  CALL print_gauxc_molecular_virial(exc_grad%exc_grad, particle_set, para_env)
               END IF
               IF (molecular_virial_debug) THEN
                  CALL debug_gauxc_molecular_virial( &
                     exc_grad%exc_grad, particle_set, qs_kind_set, density_scalar, nspins, &
                     model_name, xc_fun_name, grid_type, radial_quadrature, pruning_scheme, &
                     lb_exec_space, int_exec_space, batch_size, molecular_virial_debug_dx, para_env)
               END IF
               DEALLOCATE (exc_grad%exc_grad)
            END IF
         ELSE
            CPASSERT(nspins == 2)
            ! In here:
            ! scalar <- rho_ao(1, :) + rho_ao(2, :)
            ! zeta   <- rho_ao(1, :) - rho_ao(2, :)
            CALL dbcsr_to_dense(rho_ao(2, img), density_zeta)
            CALL para_env%sum(density_zeta)
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
                              gauxc_status, &
                              model=TRIM(model_name))
            CALL gauxc_check_status(gauxc_status)
            IF (calculate_forces) THEN
               IF (use_fd_gradient) THEN
                  CALL gauxc_xc_gradient_fd( &
                     particle_set, qs_kind_set, density_scalar, nspins, model_name, &
                     xc_fun_name, grid_type, radial_quadrature, pruning_scheme, &
                     lb_exec_space, int_exec_space, batch_size, gapw_fd_gradient_dx, &
                     para_env, exc_grad%exc_grad, density_zeta=density_zeta)
               ELSE IF (use_gradient_self_runtime) THEN
                  exc_grad = gauxc_compute_xc_gradient( &
                             gauxc_gradient_integrator_result, &
                             density_scalar, &
                             density_zeta, &
                             nspins, &
                             natom, &
                             gauxc_status, &
                             model=TRIM(model_name))
               ELSE
                  exc_grad = gauxc_compute_xc_gradient( &
                             gauxc_integrator_result, &
                             density_scalar, &
                             density_zeta, &
                             nspins, &
                             natom, &
                             gauxc_status, &
                             model=TRIM(model_name))
               END IF
               CALL gauxc_check_status(gauxc_status)
               CALL add_gauxc_gradient_to_force( &
                  exc_grad%exc_grad, &
                  force, &
                  atomic_kind_set, &
                  para_env)
               IF (molecular_virial) THEN
                  CALL print_gauxc_molecular_virial(exc_grad%exc_grad, particle_set, para_env)
               END IF
               IF (molecular_virial_debug) THEN
                  CALL debug_gauxc_molecular_virial( &
                     exc_grad%exc_grad, particle_set, qs_kind_set, density_scalar, nspins, &
                     model_name, xc_fun_name, grid_type, radial_quadrature, pruning_scheme, &
                     lb_exec_space, int_exec_space, batch_size, molecular_virial_debug_dx, para_env, &
                     density_zeta=density_zeta)
               END IF
               DEALLOCATE (exc_grad%exc_grad)
            END IF
         END IF

         energy%exc = energy%exc + gauxc_xc_result%exc

         IF (nspins == 1) THEN
            IF (img == 1) THEN
               matrix_vxc(1) = dense_to_dbcsr(gauxc_xc_result%vxc_scalar, rho_ao(1, img))
            ELSE
               CPABORT("UNIMPLEMENTED: Handling multiple result matrices in k-point integration")
            END IF
         ELSE
            CPASSERT(nspins == 2)
            ! Transform derivatives from total/spin density back to alpha/beta channels.
            vxc_zeta_tmp = dense_to_dbcsr(gauxc_xc_result%vxc_zeta, rho_ao(1, img))
            IF (img == 1) THEN
               DO ispin = 1, 2
                  matrix_vxc(ispin) = dense_to_dbcsr(gauxc_xc_result%vxc_scalar, rho_ao(ispin, 1))
                  CALL dbcsr_add( &
                     matrix_vxc(ispin)%matrix, &
                     vxc_zeta_tmp%matrix, &
                     1.0_dp, &
                     ! 1.0 for ispin==1, -1.0 for ispin==2
                     1.0_dp - REAL(ispin - 1, dp)*2.0_dp)
               END DO
            ELSE
               CPABORT("UNIMPLEMENTED: Handling multiple result matrices in k-point integration")
            END IF
            CALL dbcsr_release(vxc_zeta_tmp%matrix)
            DEALLOCATE (vxc_zeta_tmp%matrix)
         END IF
      END DO

      DEALLOCATE (density_scalar)
      IF (ALLOCATED(density_zeta)) DEALLOCATE (density_zeta)
      DEALLOCATE (gauxc_xc_result%vxc_scalar)
      IF (ALLOCATED(gauxc_xc_result%vxc_zeta)) DEALLOCATE (gauxc_xc_result%vxc_zeta)

      CALL set_ks_env(ks_env, matrix_vxc=matrix_vxc)
      DO ispin = 1, nspins
         CALL dbcsr_finalize(matrix_vxc(ispin)%matrix)
      END DO

      IF (use_gradient_self_runtime) THEN
         CALL gauxc_destroy_integrator(gauxc_gradient_integrator_result, gauxc_status)
         CALL gauxc_check_status(gauxc_status)
         CALL gauxc_destroy_grid(gauxc_gradient_grid_result, gauxc_status)
         CALL gauxc_check_status(gauxc_status)
      END IF
      CALL gauxc_destroy_integrator(gauxc_integrator_result, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      CALL gauxc_destroy_grid(gauxc_grid_result, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      CALL gauxc_destroy_basisset(gauxc_basis, gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      CALL gauxc_destroy_molecule(gauxc_mol, gauxc_status)
      CALL gauxc_check_status(gauxc_status)

   END SUBROUTINE apply_gauxc

END MODULE xc_gauxc_functional
