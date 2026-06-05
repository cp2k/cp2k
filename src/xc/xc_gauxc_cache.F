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

   USE iso_c_binding,                   ONLY: c_double
   USE kinds,                           ONLY: default_string_length
   USE message_passing,                 ONLY: mp_comm_self,&
                                              mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE xc_gauxc_interface,              ONLY: &
        cp_gauxc_basisset_type, cp_gauxc_grid_type, cp_gauxc_integrator_type, &
        cp_gauxc_molecule_type, cp_gauxc_status_type, gauxc_check_status, gauxc_create_basisset, &
        gauxc_create_grid, gauxc_create_integrator, gauxc_create_molecule, gauxc_destroy_basisset, &
        gauxc_destroy_grid, gauxc_destroy_integrator, gauxc_destroy_molecule
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
      CHARACTER(LEN=default_string_length)         :: model_name = ""
      LOGICAL                                      :: use_self_runtime = .FALSE.
      LOGICAL                                      :: use_gradient_self_runtime = .FALSE.
   END TYPE cp_gauxc_cache_params

   TYPE, extends(cp_gauxc_cache_params)            :: cp_gauxc_cache_type
      TYPE(cp_gauxc_molecule_type)                 :: molecule
      TYPE(cp_gauxc_basisset_type)                 :: basisset
      TYPE(cp_gauxc_grid_type)                     :: grid
      TYPE(cp_gauxc_grid_type)                     :: gradient_grid
      TYPE(cp_gauxc_integrator_type)               :: integrator
      TYPE(cp_gauxc_integrator_type)               :: gradient_integrator
      LOGICAL                                      :: is_init = .FALSE.
      LOGICAL                                      :: use_fd_gradient = .FALSE.
      REAL(c_double), DIMENSION(:, :), ALLOCATABLE :: last_positions

   CONTAINS

      PROCEDURE, PUBLIC                            :: needs_gradient_grid => gauxc_cache_type_needs_gradient_grid
      PROCEDURE, PUBLIC                            :: max_l => gauxc_cache_type_max_l
   END TYPE cp_gauxc_cache_type

   PUBLIC :: &
      cp_gauxc_cache_type, &
      gauxc_cache_init, &
      gauxc_cache_is_valid, &
      gauxc_cache_release, cp_gauxc_cache_params

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param this ...
!> \return ...
! **************************************************************************************************
   FUNCTION gauxc_cache_type_needs_gradient_grid(this) RESULT(res)
      class(cp_gauxc_cache_type) :: this
      LOGICAL :: res

      res = this%use_gradient_self_runtime
   END FUNCTION gauxc_cache_type_needs_gradient_grid

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
!> \param gapw ...
!> \param calculate_forces ...
!> \param use_skala_model ...
!> \param status ...
! **************************************************************************************************
   SUBROUTINE gauxc_cache_init( &
      cache, &
      params, &
      para_env, &
      particle_set, &
      qs_kind_set, &
      gapw, &
      calculate_forces, &
      use_skala_model, &
      status)

      TYPE(cp_gauxc_cache_type), INTENT(inout)           :: cache
      TYPE(cp_gauxc_cache_params), INTENT(in)            :: params
      TYPE(mp_para_env_type), INTENT(in), POINTER        :: para_env
      TYPE(particle_type), DIMENSION(:), INTENT(in), &
         POINTER                                         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), INTENT(in), &
         POINTER                                         :: qs_kind_set
      LOGICAL, INTENT(in)                                :: gapw, calculate_forces, use_skala_model
      TYPE(cp_gauxc_status_type), INTENT(inout)          :: status

      IF (gauxc_cache_is_valid( &
          cache, &
          params, &
          particle_set)) THEN
         status%status%code = 0
         RETURN
      END IF

      CALL gauxc_cache_release(cache)

      cache%molecule = gauxc_create_molecule( &
                       particle_set, &
                       status)
      IF (status%status%code /= 0) GO TO 100

      cache%basisset = gauxc_create_basisset( &
                       qs_kind_set, &
                       particle_set, &
                       status)
      IF (status%status%code /= 0) GO TO 200

      cache%use_fd_gradient = gapw .AND. calculate_forces .AND. &
                              (use_skala_model .OR. cache%basisset%max_l > 3)
      IF (cache%use_fd_gradient .AND. para_env%mepos == 0) THEN
         CALL cp_warn( &
            __LOCATION__, &
            "Using finite-difference GauXC XC gradients for METHOD GAPW with SKALA or "// &
            "all-electron basis functions beyond f shells.  The upstream analytical GauXC "// &
            "gradient path is not yet reliable for this case.")
      END IF

      IF (params%use_self_runtime) THEN
         cache%grid = gauxc_create_grid( &
                      cache%molecule, &
                      cache%basisset, &
                      params%grid_type, &
                      params%radial_quadrature, &
                      params%pruning_scheme, &
                      params%lb_exec_space, &
                      params%batch_size, &
                      params%device_runtime_fill_fraction, &
                      status, &
                      mpi_comm=mp_comm_self%get_handle())
      ELSE
         cache%grid = gauxc_create_grid( &
                      cache%molecule, &
                      cache%basisset, &
                      params%grid_type, &
                      params%radial_quadrature, &
                      params%pruning_scheme, &
                      params%lb_exec_space, &
                      params%batch_size, &
                      params%device_runtime_fill_fraction, &
                      status, &
                      mpi_comm=para_env%get_handle())
      END IF
      IF (status%status%code /= 0) GO TO 300

      cache%integrator = gauxc_create_integrator( &
                         TRIM(params%xc_fun_name), &
                         cache%grid, &
                         params%int_exec_space, &
                         params%nspins, &
                         status)
      IF (status%status%code /= 0) GO TO 400

      IF (params%use_gradient_self_runtime) THEN
         cache%gradient_grid = gauxc_create_grid( &
                               cache%molecule, &
                               cache%basisset, &
                               params%grid_type, &
                               params%radial_quadrature, &
                               params%pruning_scheme, &
                               params%lb_exec_space, &
                               params%batch_size, &
                               params%device_runtime_fill_fraction, &
                               status, &
                               mpi_comm=mp_comm_self%get_handle())
         IF (status%status%code /= 0) GO TO 500
         cache%gradient_integrator = gauxc_create_integrator( &
                                     TRIM(params%xc_fun_name), &
                                     cache%gradient_grid, &
                                     params%int_exec_space, &
                                     params%nspins, &
                                     status)
         IF (status%status%code /= 0) GO TO 600
      END IF

      cache%natom = params%natom
      cache%nspins = params%nspins
      cache%batch_size = params%batch_size
      cache%device_runtime_fill_fraction = REAL(params%device_runtime_fill_fraction, c_double)
      cache%xc_fun_name = TRIM(params%xc_fun_name)
      cache%grid_type = TRIM(params%grid_type)
      cache%radial_quadrature = TRIM(params%radial_quadrature)
      cache%pruning_scheme = TRIM(params%pruning_scheme)
      cache%lb_exec_space = TRIM(params%lb_exec_space)
      cache%int_exec_space = TRIM(params%int_exec_space)
      cache%model_name = TRIM(params%model_name)
      cache%use_gradient_self_runtime = params%use_gradient_self_runtime
      cache%use_self_runtime = params%use_self_runtime
      IF (ALLOCATED(cache%last_positions)) &
         DEALLOCATE (cache%last_positions)
      ALLOCATE (cache%last_positions(3, params%natom))

      BLOCK
         INTEGER :: iatom
         DO iatom = 1, params%natom
            cache%last_positions(:, iatom) = REAL(particle_set(iatom)%r(:), c_double)
         END DO
      END BLOCK

      cache%is_init = .TRUE.

      RETURN

      ! Cleanup section in the error case.

100   CONTINUE
      CALL gauxc_destroy_molecule(cache%molecule, status)
      CALL gauxc_check_status(status)
      GO TO 999

200   CONTINUE
      CALL gauxc_destroy_basisset(cache%basisset, status)
      CALL gauxc_check_status(status)
      GO TO 100

300   CONTINUE
      CALL gauxc_destroy_grid(cache%grid, status)
      CALL gauxc_check_status(status)
      GO TO 200

400   CONTINUE
      CALL gauxc_destroy_integrator(cache%integrator, status)
      CALL gauxc_check_status(status)
      GO TO 300

500   CONTINUE
      CALL gauxc_destroy_grid(cache%gradient_grid, status)
      CALL gauxc_check_status(status)
      GO TO 400

600   CONTINUE
      CALL gauxc_destroy_integrator(cache%gradient_integrator, status)
      CALL gauxc_check_status(status)

999   CONTINUE
      IF (ALLOCATED(cache%last_positions)) DEALLOCATE (cache%last_positions)
      cache%is_init = .FALSE.

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
         IF (cache%needs_gradient_grid()) THEN
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
      IF (ALLOCATED(cache%last_positions)) DEALLOCATE (cache%last_positions)
      cache%is_init = .FALSE.
#else
      MARK_USED(cache)
      IF (ALLOCATED(cache%last_positions)) DEALLOCATE (cache%last_positions)
      cache%is_init = .FALSE.
#endif
   END SUBROUTINE gauxc_cache_release

! **************************************************************************************************
!> \brief Check if the GauXC cache is valid for the given parameters
!> \param cache ...
!> \param params ...
!> \param particle_set ...
!> \return ...
!> \retval cache_valid ...
! **************************************************************************************************
   FUNCTION gauxc_cache_is_valid( &
      cache, &
      params, &
      particle_set) &
      RESULT(cache_valid)

      TYPE(cp_gauxc_cache_type), INTENT(IN)                :: cache
      TYPE(cp_gauxc_cache_params), INTENT(in) :: params
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
      IF (ABS(cache%device_runtime_fill_fraction - REAL(params%device_runtime_fill_fraction, c_double)) > &
          1.0E-10_c_double) RETURN
      IF (TRIM(cache%xc_fun_name) /= TRIM(params%xc_fun_name)) RETURN
      IF (TRIM(cache%grid_type) /= TRIM(params%grid_type)) RETURN
      IF (TRIM(cache%radial_quadrature) /= TRIM(params%radial_quadrature)) RETURN
      IF (TRIM(cache%pruning_scheme) /= TRIM(params%pruning_scheme)) RETURN
      IF (TRIM(cache%lb_exec_space) /= TRIM(params%lb_exec_space)) RETURN
      IF (TRIM(cache%int_exec_space) /= TRIM(params%int_exec_space)) RETURN
      IF (TRIM(cache%model_name) /= TRIM(params%model_name)) RETURN
      IF (cache%use_gradient_self_runtime .NEQV. params%use_gradient_self_runtime) RETURN
      IF (cache%use_self_runtime .NEQV. params%use_self_runtime) RETURN
      IF (.NOT. ALLOCATED(cache%last_positions)) RETURN
      IF (SIZE(cache%last_positions, 2) /= params%natom) RETURN
      DO iatom = 1, params%natom
         IF (ANY(ABS(cache%last_positions(:, iatom) - &
                     REAL(particle_set(iatom)%r(:), c_double)) > pos_tol)) RETURN
      END DO
      cache_valid = .TRUE.
#else
      MARK_USED(cache)
      MARK_USED(particle_set)
      MARK_USED(params)
      MARK_USED(use_gradient_self_runtime)
      MARK_USED(gapw)
      MARK_USED(calculate_forces)
      MARK_USED(use_skala_model)
      MARK_USED(use_self_runtime)
      MARK_USED(para_env)
      cache_valid = .FALSE.
#endif
   END FUNCTION gauxc_cache_is_valid

END MODULE xc_gauxc_cache
