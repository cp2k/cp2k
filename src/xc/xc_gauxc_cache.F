!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

#ifdef __GAUXC
#include "gauxc/gauxc_config.f"
#endif

MODULE xc_gauxc_cache

   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE iso_c_binding,                   ONLY: c_double
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length
   USE message_passing,                 ONLY: mp_comm_self,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE xc_gauxc_gradient_policy,        ONLY: gauxc_resolve_gradient_policy
   USE xc_gauxc_interface,              ONLY: &
        cp_gauxc_basisset_type, cp_gauxc_grid_type, cp_gauxc_integrator_type, &
        cp_gauxc_molecule_type, cp_gauxc_status_type, cp_gauxc_xc_gradient_type, cp_gauxc_xc_type, &
        gauxc_check_status, gauxc_create_basisset, gauxc_create_grid, gauxc_create_integrator, &
        gauxc_create_molecule, gauxc_destroy_basisset, gauxc_destroy_grid, &
        gauxc_destroy_integrator, gauxc_destroy_molecule
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE cp_gauxc_cache_params
      INTEGER                                      :: natom = -1
      INTEGER                                      :: nspins = -1
      INTEGER                                      :: batch_size = -1
      REAL(c_double)                               :: device_runtime_fill_fraction = 0.0_c_double
      CHARACTER(LEN=default_string_length)         :: xc_fun_name = ""
      CHARACTER(LEN=default_string_length)         :: grid_type = ""
      CHARACTER(LEN=default_string_length)         :: radial_quadrature = ""
      CHARACTER(LEN=default_string_length)         :: pruning_scheme = ""
      CHARACTER(LEN=default_string_length)         :: lb_exec_space = ""
      CHARACTER(LEN=default_string_length)         :: int_exec_space = ""
      CHARACTER(LEN=default_path_length)           :: model_eval_name = ""
      CHARACTER(len=default_string_length)         :: lwd_kernel = ""
      CHARACTER(len=default_string_length)         :: gradient_runtime = "AUTO"
      LOGICAL                                      :: need_xc_gradient = .FALSE.
      LOGICAL                                      :: gapw_method = .FALSE.
      LOGICAL                                      :: force_fd_gradient = .FALSE.
      LOGICAL                                      :: use_gradient_mpi_runtime = .FALSE.
      LOGICAL                                      :: use_self_runtime = .FALSE.
      LOGICAL                                      :: use_gradient_self_runtime = .FALSE.
      LOGICAL                                      :: use_fd_gradient = .FALSE.
      LOGICAL                                      :: use_gauxc_model = .FALSE.
      INTEGER                                      :: nimages = -1
      INTEGER                                      :: nkind = -1
      CHARACTER(LEN=default_string_length), ALLOCATABLE :: kind_basis_names(:)
   END TYPE cp_gauxc_cache_params

   TYPE dense_to_dbcsr_layout_type
      LOGICAL                                      :: valid = .FALSE.
      INTEGER                                      :: nrows = -1
      INTEGER                                      :: ncols = -1
      INTEGER                                      :: nblkrows = -1
      INTEGER                                      :: nblkcols = -1
      INTEGER                                      :: mynode = -1
      INTEGER                                      :: myprow = -1
      INTEGER                                      :: mypcol = -1
      INTEGER, ALLOCATABLE, DIMENSION(:)           :: row_blk_size
      INTEGER, ALLOCATABLE, DIMENSION(:)           :: col_blk_size
      INTEGER, ALLOCATABLE, DIMENSION(:)           :: row_dist
      INTEGER, ALLOCATABLE, DIMENSION(:)           :: col_dist
      INTEGER, ALLOCATABLE, DIMENSION(:)           :: r_offset
      INTEGER, ALLOCATABLE, DIMENSION(:)           :: c_offset
      INTEGER, ALLOCATABLE, DIMENSION(:)           :: blk_row
      INTEGER, ALLOCATABLE, DIMENSION(:)           :: blk_col
      REAL(c_double), ALLOCATABLE, DIMENSION(:, :) :: blk
   END TYPE dense_to_dbcsr_layout_type

   TYPE, extends(cp_gauxc_cache_params)            :: cp_gauxc_cache_type
      TYPE(cp_gauxc_molecule_type)                 :: molecule = cp_gauxc_molecule_type()
      TYPE(cp_gauxc_basisset_type)                 :: basisset = cp_gauxc_basisset_type()
      TYPE(cp_gauxc_grid_type)                     :: grid = cp_gauxc_grid_type()
      TYPE(cp_gauxc_grid_type)                     :: gradient_grid = cp_gauxc_grid_type()
      TYPE(cp_gauxc_integrator_type)               :: integrator = cp_gauxc_integrator_type()
      TYPE(cp_gauxc_integrator_type)               :: gradient_integrator = cp_gauxc_integrator_type()
      TYPE(cp_gauxc_xc_type)                       :: xc_result = cp_gauxc_xc_type()
      TYPE(cp_gauxc_xc_gradient_type)              :: xc_gradient = cp_gauxc_xc_gradient_type()
      TYPE(dense_to_dbcsr_layout_type)             :: dense_layout
      REAL(c_double), DIMENSION(:, :), ALLOCATABLE :: density_scalar, density_zeta
      INTEGER                                      :: mpi_comm = -1
      LOGICAL                                      :: is_init = .FALSE.
      LOGICAL                                      :: gradient_is_init = .FALSE.
      LOGICAL                                      :: gradient_notice_printed = .FALSE.
      LOGICAL                                      :: execution_spaces_logged = .FALSE.
      REAL(c_double), DIMENSION(:, :), ALLOCATABLE :: last_positions
   CONTAINS

      PROCEDURE, PUBLIC                            :: has_gradient_grid => gauxc_cache_type_has_gradient_grid
      PROCEDURE, PUBLIC                            :: max_l => gauxc_cache_type_max_l
   END TYPE cp_gauxc_cache_type

   PUBLIC :: &
      cp_gauxc_cache_type, &
      dense_to_dbcsr_layout_type, &
      gauxc_cache_init, &
      gauxc_cache_is_valid, &
      gauxc_cache_ensure_gradient, &
      gauxc_cache_release, cp_gauxc_cache_params

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_cache_type_has_gradient_grid(this) RESULT(res)
      class(cp_gauxc_cache_type) :: this
      LOGICAL :: res

      res = this%gradient_is_init
   END FUNCTION gauxc_cache_type_has_gradient_grid

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_cache_type_max_l(this) RESULT(res)
      class(cp_gauxc_cache_type) :: this
      INTEGER :: res

      res = this%basisset%max_l
   END FUNCTION gauxc_cache_type_max_l

! **************************************************************************************************
!> \brief ...
!> \param cache ...
!> \param params ...
!> \param para_env ...
!> \param particle_set ...
!> \param qs_kind_set ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_cache_init( &
      cache, &
      params, &
      para_env, &
      particle_set, &
      qs_kind_set, &
      status)

      TYPE(cp_gauxc_cache_type), INTENT(inout)           :: cache
      TYPE(cp_gauxc_cache_params), INTENT(inout)         :: params
      TYPE(mp_para_env_type), INTENT(in), POINTER        :: para_env
      TYPE(particle_type), DIMENSION(:), INTENT(in), &
         POINTER                                         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), INTENT(in), &
         POINTER                                         :: qs_kind_set
      TYPE(cp_gauxc_status_type), INTENT(inout)          :: status

#ifdef __GAUXC
      INTEGER                                            :: gradient_error, unit_nr

      IF (gauxc_cache_is_valid( &
          cache, &
          params, &
          para_env, &
          particle_set)) THEN
         params%use_fd_gradient = cache%use_fd_gradient
         params%use_gradient_mpi_runtime = cache%use_gradient_mpi_runtime
         params%use_gradient_self_runtime = cache%use_gradient_self_runtime
         status%status%code = 0
         RETURN
      END IF

      CALL gauxc_cache_release(cache)

      cache%molecule = gauxc_create_molecule( &
                       particle_set, &
                       status)
      CALL gauxc_check_status(status)

      cache%basisset = gauxc_create_basisset( &
                       qs_kind_set, &
                       particle_set, &
                       status)
      CALL gauxc_check_status(status)

      CALL gauxc_resolve_gradient_policy( &
         params%need_xc_gradient, params%gapw_method, cache%max_l(), params%use_gauxc_model, &
         para_env%num_pe, params%use_self_runtime, TRIM(params%gradient_runtime), &
         params%force_fd_gradient, params%use_fd_gradient, params%use_gradient_mpi_runtime, &
         params%use_gradient_self_runtime, gradient_error)
      SELECT CASE (gradient_error)
      CASE (1)
         CALL cp_abort(__LOCATION__, &
                       "MODEL_GRADIENT_RUNTIME MPI requires verified GauXC OneDFT gradient support (PR #222).")
      CASE (2)
         CALL cp_abort(__LOCATION__, &
                       "MODEL_GRADIENT_RUNTIME MPI requires an MPI energy/VXC runtime. "// &
                       "Use AUTO/SELF to retain the open-shell replicated-runtime safeguard.")
      CASE (3)
         CALL cp_abort(__LOCATION__, &
                       "MODEL_GRADIENT_RUNTIME MPI is incompatible with FORCE_NUMERICAL_GRADIENT.")
      CASE (4)
         CALL cp_abort(__LOCATION__, "Unknown GAUXC%MODEL_GRADIENT_RUNTIME value.")
      END SELECT

      IF (params%use_self_runtime) THEN
         ! open-shell densities across CP2K MPI ranks.
         cache%grid = gauxc_create_grid( &
                      cache%molecule, cache%basisset, params%grid_type, params%radial_quadrature, &
                      params%pruning_scheme, params%lb_exec_space, params%batch_size, &
                      params%device_runtime_fill_fraction, status, &
                      mpi_comm=mp_comm_self%get_handle(), force_new_runtime=.TRUE., &
                      int_exec_space=params%int_exec_space)
      ELSE
         ! Use the QS force-evaluation communicator rather than GauXC's global
         ! runtime.  In mixed CDFT the diabatic states can be built on disjoint
         ! MPI subgroups; using the global communicator would make GauXC wait
         ! for ranks that are working on another state.
         cache%grid = gauxc_create_grid( &
                      cache%molecule, cache%basisset, params%grid_type, params%radial_quadrature, &
                      params%pruning_scheme, params%lb_exec_space, params%batch_size, &
                      params%device_runtime_fill_fraction, status, &
                      mpi_comm=para_env%get_handle(), force_new_runtime=.FALSE., &
                      int_exec_space=params%int_exec_space)
      END IF
      CALL gauxc_check_status(status)

      cache%integrator = gauxc_create_integrator( &
                         TRIM(params%xc_fun_name), &
                         cache%grid, &
                         params%int_exec_space, &
                         params%lwd_kernel, &
                         params%nspins, &
                         status)
      CALL gauxc_check_status(status)

      IF (.NOT. cache%execution_spaces_logged .AND. para_env%mepos == 0) THEN
         unit_nr = cp_logger_get_default_io_unit()
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(A)') " GAUXC| Execution spaces: LB "// &
               TRIM(params%lb_exec_space)//", INT "//TRIM(params%int_exec_space)// &
               ", runtime "//MERGE("DEVICE", "HOST  ", cache%grid%uses_device_runtime)
         END IF
         cache%execution_spaces_logged = .TRUE.
      END IF

      cache%gradient_is_init = .FALSE.
      cache%gradient_notice_printed = .FALSE.

      cache%natom = params%natom
      cache%nspins = params%nspins
      cache%batch_size = params%batch_size
      cache%device_runtime_fill_fraction = params%device_runtime_fill_fraction
      cache%xc_fun_name = TRIM(params%xc_fun_name)
      cache%grid_type = TRIM(params%grid_type)
      cache%radial_quadrature = TRIM(params%radial_quadrature)
      cache%pruning_scheme = TRIM(params%pruning_scheme)
      cache%lb_exec_space = TRIM(params%lb_exec_space)
      cache%int_exec_space = TRIM(params%int_exec_space)
      cache%model_eval_name = TRIM(params%model_eval_name)
      cache%lwd_kernel = TRIM(params%lwd_kernel)
      cache%gradient_runtime = params%gradient_runtime
      cache%need_xc_gradient = params%need_xc_gradient
      cache%gapw_method = params%gapw_method
      cache%force_fd_gradient = params%force_fd_gradient
      cache%use_gauxc_model = params%use_gauxc_model
      cache%use_gradient_mpi_runtime = params%use_gradient_mpi_runtime
      cache%use_gradient_self_runtime = params%use_gradient_self_runtime
      cache%use_self_runtime = params%use_self_runtime
      cache%use_fd_gradient = params%use_fd_gradient
      cache%nimages = params%nimages
      cache%nkind = params%nkind
      IF (ALLOCATED(cache%kind_basis_names)) DEALLOCATE (cache%kind_basis_names)
      IF (ALLOCATED(params%kind_basis_names)) THEN
         ALLOCATE (cache%kind_basis_names(SIZE(params%kind_basis_names)))
         cache%kind_basis_names(:) = params%kind_basis_names
      END IF
      IF (ALLOCATED(cache%last_positions)) THEN
         DEALLOCATE (cache%last_positions)
      END IF
      ALLOCATE (cache%last_positions(3, params%natom))

      BLOCK
         INTEGER :: iatom
         DO iatom = 1, params%natom
            cache%last_positions(:, iatom) = particle_set(iatom)%r(:)
         END DO
      END BLOCK

      cache%mpi_comm = para_env%get_handle()
      cache%is_init = .TRUE.

#else
      MARK_USED(params)
      MARK_USED(para_env)
      MARK_USED(particle_set)
      MARK_USED(qs_kind_set)
      MARK_USED(status)
      cache%is_init = .FALSE.
#endif

   END SUBROUTINE gauxc_cache_init

! **************************************************************************************************
!> \brief Release all GauXC objects in a cache and deallocate arrays
!> \param cache ...
! **************************************************************************************************
   SUBROUTINE gauxc_cache_release(cache)
      TYPE(cp_gauxc_cache_type), INTENT(INOUT)              :: cache

#ifdef __GAUXC
      TYPE(cp_gauxc_status_type)                            :: status

      IF (cache%is_init) THEN
         IF (cache%has_gradient_grid()) THEN
            CALL gauxc_destroy_integrator(cache%gradient_integrator, status)
            CALL gauxc_check_status(status)
            CALL gauxc_destroy_grid(cache%gradient_grid, status)
            CALL gauxc_check_status(status)
         END IF
         CALL gauxc_destroy_integrator(cache%integrator, status)
         CALL gauxc_check_status(status)
         CALL gauxc_destroy_grid(cache%grid, status)
         CALL gauxc_check_status(status)
         CALL gauxc_destroy_basisset(cache%basisset, status)
         CALL gauxc_check_status(status)
         CALL gauxc_destroy_molecule(cache%molecule, status)
         CALL gauxc_check_status(status)
      END IF
#else
      MARK_USED(cache)
#endif

      IF (ALLOCATED(cache%last_positions)) DEALLOCATE (cache%last_positions)
      IF (ALLOCATED(cache%kind_basis_names)) DEALLOCATE (cache%kind_basis_names)
      IF (ALLOCATED(cache%xc_result%vxc_scalar)) DEALLOCATE (cache%xc_result%vxc_scalar)
      IF (ALLOCATED(cache%xc_result%vxc_zeta)) DEALLOCATE (cache%xc_result%vxc_zeta)
      IF (ALLOCATED(cache%xc_result%scratch_rks)) DEALLOCATE (cache%xc_result%scratch_rks)
      IF (ALLOCATED(cache%xc_result%scratch_zeta_zero)) DEALLOCATE (cache%xc_result%scratch_zeta_zero)
      IF (ALLOCATED(cache%xc_gradient%exc_grad)) DEALLOCATE (cache%xc_gradient%exc_grad)
      IF (ALLOCATED(cache%xc_gradient%scratch_rks)) DEALLOCATE (cache%xc_gradient%scratch_rks)
      IF (ALLOCATED(cache%xc_gradient%scratch_zeta_zero)) DEALLOCATE (cache%xc_gradient%scratch_zeta_zero)
      IF (ALLOCATED(cache%density_scalar)) DEALLOCATE (cache%density_scalar)
      IF (ALLOCATED(cache%density_zeta)) DEALLOCATE (cache%density_zeta)
      IF (ALLOCATED(cache%dense_layout%row_blk_size)) THEN
         DEALLOCATE (cache%dense_layout%row_blk_size, cache%dense_layout%col_blk_size, &
                     cache%dense_layout%row_dist, cache%dense_layout%col_dist, &
                     cache%dense_layout%r_offset, cache%dense_layout%c_offset, &
                     cache%dense_layout%blk_row, cache%dense_layout%blk_col, &
                     cache%dense_layout%blk)
      END IF
      cache%dense_layout%valid = .FALSE.
      cache%gradient_is_init = .FALSE.
      cache%gradient_notice_printed = .FALSE.
      cache%execution_spaces_logged = .FALSE.
      cache%is_init = .FALSE.

   END SUBROUTINE gauxc_cache_release

! **************************************************************************************************
!> \brief Build the separate gradient grid/integrator on first use and keep it cached.
!>
!> \param cache ...
!> \param params ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_cache_ensure_gradient(cache, params, status)
      TYPE(cp_gauxc_cache_type), INTENT(INOUT)           :: cache
      TYPE(cp_gauxc_cache_params), INTENT(INOUT)         :: params
      TYPE(cp_gauxc_status_type), INTENT(INOUT)          :: status

#ifdef __GAUXC
      INTEGER                                            :: error_code, unit_nr
      CHARACTER(LEN=6)                                   :: gradient_runtime_label

      IF (.NOT. params%need_xc_gradient) RETURN

      IF (.NOT. cache%gradient_notice_printed) THEN
         IF (params%use_fd_gradient .AND. .NOT. params%force_fd_gradient) THEN
            CALL cp_warn( &
               __LOCATION__, &
               "Using finite-difference GauXC XC gradients for METHOD GAPW/GAPW_XC with "// &
               "basis functions beyond f shells. Analytical high-l gradients have not been "// &
               "enabled for this GauXC installation and functional path.")
         END IF
         unit_nr = cp_logger_get_default_io_unit()
         IF (unit_nr > 0) THEN
            IF (params%use_gradient_self_runtime) THEN
               gradient_runtime_label = "SELF"
            ELSE IF (params%use_gradient_mpi_runtime) THEN
               gradient_runtime_label = "MPI"
            ELSE
               gradient_runtime_label = "SHARED"
            END IF
            WRITE (unit_nr, '(A,A)') " GAUXC| Molecular XC gradient method: ", &
               TRIM(MERGE("FD        ", "ANALYTICAL", params%use_fd_gradient))
            WRITE (unit_nr, '(A,A)') " GAUXC| Molecular XC gradient runtime: ", &
               TRIM(gradient_runtime_label)
         END IF
         cache%gradient_notice_printed = .TRUE.
      END IF

      IF (params%use_fd_gradient) RETURN
      IF (.NOT. params%use_gradient_self_runtime) RETURN
      IF (.NOT. cache%is_init) RETURN
      IF (cache%gradient_is_init) RETURN

      cache%gradient_grid = gauxc_create_grid( &
                            cache%molecule, cache%basisset, params%grid_type, &
                            params%radial_quadrature, params%pruning_scheme, &
                            params%lb_exec_space, params%batch_size, &
                            params%device_runtime_fill_fraction, status, &
                            mpi_comm=mp_comm_self%get_handle(), force_new_runtime=.TRUE., &
                            int_exec_space=params%int_exec_space)
      IF (status%status%code /= 0) THEN
         error_code = status%status%code
         CALL gauxc_destroy_grid(cache%gradient_grid, status)
         status%status%code = error_code
         RETURN
      END IF

      cache%gradient_integrator = gauxc_create_integrator( &
                                  TRIM(params%xc_fun_name), &
                                  cache%gradient_grid, &
                                  params%int_exec_space, &
                                  params%lwd_kernel, &
                                  params%nspins, &
                                  status)
      IF (status%status%code /= 0) THEN
         error_code = status%status%code
         CALL gauxc_destroy_integrator(cache%gradient_integrator, status)
         CALL gauxc_destroy_grid(cache%gradient_grid, status)
         status%status%code = error_code
         RETURN
      END IF
      cache%gradient_is_init = .TRUE.
#else
      MARK_USED(cache)
      MARK_USED(params)
      MARK_USED(status)
#endif
   END SUBROUTINE gauxc_cache_ensure_gradient

! **************************************************************************************************
!> \brief Check if the GauXC cache is valid for the given parameters
!> \param cache ...
!> \param params ...
!> \param para_env ...
!> \param particle_set ...
!> \return ...
!> \retval cache_valid ...
! **************************************************************************************************
   FUNCTION gauxc_cache_is_valid( &
      cache, &
      params, &
      para_env, &
      particle_set) &
      RESULT(cache_valid)

      TYPE(cp_gauxc_cache_type), INTENT(IN)                :: cache
      TYPE(cp_gauxc_cache_params), INTENT(in)              :: params
      TYPE(mp_para_env_type), INTENT(in), POINTER          :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER           :: particle_set
      LOGICAL                                              :: cache_valid

#ifdef __GAUXC
      INTEGER                                              :: iatom
      REAL(c_double), PARAMETER                            :: pos_tol = 1.0E-8_c_double

      cache_valid = .FALSE.
      IF (.NOT. cache%is_init) RETURN
      IF (cache%natom /= params%natom) RETURN
      IF (cache%nspins /= params%nspins) RETURN
      IF (cache%batch_size /= params%batch_size) RETURN
      IF (ABS(cache%device_runtime_fill_fraction - params%device_runtime_fill_fraction) > &
          1.0E-10_c_double) RETURN
      IF (TRIM(cache%xc_fun_name) /= TRIM(params%xc_fun_name)) RETURN
      IF (TRIM(cache%grid_type) /= TRIM(params%grid_type)) RETURN
      IF (TRIM(cache%radial_quadrature) /= TRIM(params%radial_quadrature)) RETURN
      IF (TRIM(cache%pruning_scheme) /= TRIM(params%pruning_scheme)) RETURN
      IF (TRIM(cache%lb_exec_space) /= TRIM(params%lb_exec_space)) RETURN
      IF (TRIM(cache%int_exec_space) /= TRIM(params%int_exec_space)) RETURN
      IF (TRIM(cache%model_eval_name) /= TRIM(params%model_eval_name)) RETURN
      IF (TRIM(cache%lwd_kernel) /= TRIM(params%lwd_kernel)) RETURN
      IF (cache%use_self_runtime .NEQV. params%use_self_runtime) RETURN

      IF (cache%need_xc_gradient .NEQV. params%need_xc_gradient) RETURN
      IF (cache%gradient_runtime /= params%gradient_runtime) RETURN
      IF (cache%gapw_method .NEQV. params%gapw_method) RETURN
      IF (cache%force_fd_gradient .NEQV. params%force_fd_gradient) RETURN
      IF (cache%use_gauxc_model .NEQV. params%use_gauxc_model) RETURN
      IF (cache%mpi_comm /= para_env%get_handle()) RETURN
      IF (cache%nimages /= params%nimages) RETURN
      IF (cache%nkind /= params%nkind) RETURN
      IF (ALLOCATED(cache%kind_basis_names) .NEQV. ALLOCATED(params%kind_basis_names)) RETURN
      IF (ALLOCATED(params%kind_basis_names)) THEN
         IF (SIZE(cache%kind_basis_names) /= SIZE(params%kind_basis_names)) RETURN
         IF (ANY(cache%kind_basis_names /= params%kind_basis_names)) RETURN
      END IF

      IF (.NOT. ALLOCATED(cache%last_positions)) RETURN
      IF (SIZE(cache%last_positions, 2) /= params%natom) RETURN
      DO iatom = 1, params%natom
         IF (ANY(ABS(cache%last_positions(:, iatom) - &
                     particle_set(iatom)%r(:)) > pos_tol)) RETURN
      END DO
      cache_valid = .TRUE.
#else
      MARK_USED(cache)
      MARK_USED(particle_set)
      MARK_USED(params)
      MARK_USED(para_env)
      cache_valid = .FALSE.
#endif
   END FUNCTION gauxc_cache_is_valid

END MODULE xc_gauxc_cache
