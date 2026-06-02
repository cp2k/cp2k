!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

#ifdef __GAUXC
#include "gauxc/gauxc_config.f"
#endif

#define GAUXC_RETURN_IF_ERROR(status) IF (status%status%code /= 0) RETURN

MODULE xc_gauxc_cache

   USE iso_c_binding,                   ONLY: c_double
   USE kinds,                           ONLY: default_string_length,&
                                              dp
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

   TYPE cp_gauxc_cache_type
#ifdef __GAUXC
      TYPE(cp_gauxc_molecule_type)       :: molecule
      TYPE(cp_gauxc_basisset_type)       :: basisset
      TYPE(cp_gauxc_grid_type)           :: grid
      TYPE(cp_gauxc_grid_type)           :: gradient_grid
      TYPE(cp_gauxc_integrator_type)     :: integrator
      TYPE(cp_gauxc_integrator_type)     :: gradient_integrator
      INTEGER                             :: natom = -1
      INTEGER                             :: nspins = -1
      INTEGER                             :: batch_size = -1
      REAL(c_double)                      :: device_runtime_ff = 0.0_c_double
      CHARACTER(LEN=default_string_length) :: xc_fun_name = ""
      CHARACTER(LEN=default_string_length) :: grid_type = ""
      CHARACTER(LEN=default_string_length) :: radial_quad = ""
      CHARACTER(LEN=default_string_length) :: pruning_scheme = ""
      CHARACTER(LEN=default_string_length) :: lb_exec_space = ""
      CHARACTER(LEN=default_string_length) :: int_exec_space = ""
      CHARACTER(LEN=default_string_length) :: model_name = ""
      LOGICAL                             :: gradient_self_runtime = .FALSE.
      LOGICAL                             :: is_init = .FALSE. ! param
      REAL(c_double), DIMENSION(:, :), ALLOCATABLE :: last_positions ! param
      TYPE(cp_gauxc_status_type)         :: status
   CONTAINS
      PROCEDURE, PUBLIC :: needs_gradient_grid => gauxc_cache_type_needs_gradient_grid
      PROCEDURE, PUBLIC :: max_l => gauxc_cache_type_max_l
#endif
   END TYPE cp_gauxc_cache_type

   PUBLIC :: &
      cp_gauxc_cache_type, &
      gauxc_cache_init, &
      gauxc_cache_is_valid, &
      gauxc_cache_release

CONTAINS

#ifdef __GAUXC
   FUNCTION gauxc_cache_type_needs_gradient_grid(this) RESULT(res)
      class(cp_gauxc_cache_type) :: this
      LOGICAL :: res

      res = this%gradient_self_runtime
   END FUNCTION

   FUNCTION gauxc_cache_type_max_l(this) RESULT(res)
      class(cp_gauxc_cache_type) :: this
      INTEGER :: res

      res = this%basisset%max_l
   END FUNCTION
#endif

! **************************************************************************************************
!> \brief ...
!> \param ! ...
!> \param inout ...
!> \param params ...
! **************************************************************************************************
   SUBROUTINE gauxc_cache_init( &
      ! inout params
      cache, &
      gauxc_mol, & ! param
      gauxc_basis, & ! param
      gauxc_grid_result, & ! param
      gauxc_integrator_result, & ! param
      gauxc_gradient_grid_result, & ! param
      gauxc_gradient_integrator_result, & ! param
      use_fd_gradient, &
      ! input data
      para_env, &
      particle_set, &
      qs_kind_set, &
      ! params
      gapw, &
      calculate_forces, &
      natom, & ! param
      nspins, & ! param
      grid_type, & ! param
      xc_fun_name, & ! param
      int_exec_space, & ! param
      lb_exec_space, & ! param
      model_name, & ! param
      pruning_scheme, & ! param
      radial_quadrature, & ! param
      device_runtime_fill_fraction, & ! param
      use_skala_model, &
      use_gradient_mpi_runtime, &
      use_gradient_self_runtime, & ! param
      use_self_runtime, &
      batch_size)

      TYPE(cp_gauxc_cache_type), INTENT(inout)           :: cache
      TYPE(cp_gauxc_molecule_type), INTENT(inout)        :: gauxc_mol
      TYPE(cp_gauxc_basisset_type), INTENT(inout)        :: gauxc_basis
      TYPE(cp_gauxc_grid_type), INTENT(inout)            :: gauxc_grid_result
      TYPE(cp_gauxc_integrator_type), INTENT(inout)      :: gauxc_integrator_result
      TYPE(cp_gauxc_grid_type), INTENT(inout)            :: gauxc_gradient_grid_result
      TYPE(cp_gauxc_integrator_type), INTENT(inout) :: gauxc_gradient_integrator_result
      LOGICAL, INTENT(inout)                             :: use_fd_gradient
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      LOGICAL:: gapw, use_skala_model, calculate_forces, use_gradient_mpi_runtime, use_self_runtime, use_gradient_self_runtime
      CHARACTER(len=*) :: grid_type, int_exec_space, lb_exec_space, model_name, pruning_scheme, radial_quadrature, xc_fun_name
      REAL(KIND=dp)                                      :: device_runtime_fill_fraction
      INTEGER :: batch_size, natom, nspins

      TYPE(cp_gauxc_status_type)                         :: gauxc_status

      IF (gauxc_cache_is_valid( &
          cache, &
          particle_set, &
          natom, &
          nspins, &
          batch_size, &
          device_runtime_fill_fraction, &
          xc_fun_name, &
          grid_type, &
          radial_quadrature, &
          pruning_scheme, &
          lb_exec_space, &
          int_exec_space, &
          model_name, &
          use_gradient_self_runtime)) THEN
         gauxc_mol = cache%molecule
         gauxc_basis = cache%basisset
         gauxc_grid_result = cache%grid
         gauxc_integrator_result = cache%integrator
         IF (use_gradient_self_runtime) THEN
            gauxc_gradient_grid_result = cache%gradient_grid
            gauxc_gradient_integrator_result = cache%gradient_integrator
         END IF

         ! TODO remove duplication
         use_fd_gradient = gapw .AND. calculate_forces .AND. &
                           (use_skala_model .OR. gauxc_basis%max_l > 3)
         RETURN
      END IF

      CALL gauxc_cache_release(cache)

      gauxc_mol = gauxc_create_molecule( &
                  particle_set, &
                  gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      gauxc_basis = gauxc_create_basisset( &
                    qs_kind_set, &
                    particle_set, &
                    gauxc_status)
      CALL gauxc_check_status(gauxc_status)
      use_fd_gradient = gapw .AND. calculate_forces .AND. &
                        (use_skala_model .OR. gauxc_basis%max_l > 3)
      IF (use_gradient_mpi_runtime) THEN
         use_gradient_self_runtime = .FALSE.
      END IF
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
                             device_runtime_fill_fraction, &
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
                             device_runtime_fill_fraction, &
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
                                      device_runtime_fill_fraction, &
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

      cache%molecule = gauxc_mol
      cache%basisset = gauxc_basis
      cache%grid = gauxc_grid_result
      cache%integrator = gauxc_integrator_result
      IF (use_gradient_self_runtime) THEN
         cache%gradient_grid = gauxc_gradient_grid_result
         cache%gradient_integrator = gauxc_gradient_integrator_result
      END IF
      cache%natom = natom
      cache%nspins = nspins
      cache%batch_size = batch_size
      cache%device_runtime_ff = REAL(device_runtime_fill_fraction, c_double)
      cache%xc_fun_name = TRIM(xc_fun_name)
      cache%grid_type = TRIM(grid_type)
      cache%radial_quad = TRIM(radial_quadrature)
      cache%pruning_scheme = TRIM(pruning_scheme)
      cache%lb_exec_space = TRIM(lb_exec_space)
      cache%int_exec_space = TRIM(int_exec_space)
      cache%model_name = TRIM(model_name)
      cache%gradient_self_runtime = use_gradient_self_runtime
      IF (ALLOCATED(cache%last_positions)) &
         DEALLOCATE (cache%last_positions)
      ALLOCATE (cache%last_positions(3, natom))

      BLOCK
         INTEGER :: iatom
         DO iatom = 1, natom
            cache%last_positions(:, iatom) = REAL(particle_set(iatom)%r(:), c_double)
         END DO
      END BLOCK

      cache%is_init = .TRUE.

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
            CALL gauxc_destroy_grid(cache%gradient_grid, status)
         END IF
         CALL gauxc_destroy_integrator(cache%integrator, status)
         CALL gauxc_destroy_grid(cache%grid, status)
         CALL gauxc_destroy_basisset(cache%basisset, status)
         CALL gauxc_destroy_molecule(cache%molecule, status)
      END IF
      IF (ALLOCATED(cache%last_positions)) DEALLOCATE (cache%last_positions)
      cache%is_init = .FALSE.
#else
      MARK_USED(cache)
#endif
   END SUBROUTINE gauxc_cache_release

! **************************************************************************************************
!> \brief Check if the GauXC cache is valid for the given parameters
!> \param cache ...
!> \param particle_set ...
!> \param natom ...
!> \param nspins ...
!> \param batch_size ...
!> \param device_runtime_fill_fraction ...
!> \param xc_fun_name ...
!> \param grid_type ...
!> \param radial_quadrature ...
!> \param pruning_scheme ...
!> \param lb_exec_space ...
!> \param int_exec_space ...
!> \param model_name ...
!> \param use_gradient_self_runtime ...
!> \return ...
!> \retval cache_valid ...
! **************************************************************************************************
   FUNCTION gauxc_cache_is_valid( &
      cache, &
      particle_set, &
      natom, &
      nspins, &
      batch_size, &
      device_runtime_fill_fraction, &
      xc_fun_name, &
      grid_type, &
      radial_quadrature, &
      pruning_scheme, &
      lb_exec_space, &
      int_exec_space, &
      model_name, &
      use_gradient_self_runtime) &
      RESULT(cache_valid)

      TYPE(cp_gauxc_cache_type), INTENT(IN)                :: cache
      TYPE(particle_type), DIMENSION(:), POINTER           :: particle_set
      INTEGER, INTENT(IN)                                  :: natom, nspins, batch_size
      REAL(KIND=dp), INTENT(IN)                            :: device_runtime_fill_fraction
      CHARACTER(LEN=*), INTENT(IN)                         :: xc_fun_name, grid_type, &
                                                              radial_quadrature, pruning_scheme, &
                                                              lb_exec_space, int_exec_space, &
                                                              model_name
      LOGICAL, INTENT(IN)                                  :: use_gradient_self_runtime
      LOGICAL                                              :: cache_valid

#ifdef __GAUXC
      INTEGER                                              :: iatom
      REAL(c_double), PARAMETER                            :: pos_tol = 1.0E-8_c_double

      cache_valid = .FALSE.
      IF (.NOT. cache%is_init) RETURN
      IF (cache%natom /= natom) RETURN
      IF (cache%nspins /= nspins) RETURN
      IF (cache%batch_size /= batch_size) RETURN
      IF (ABS(cache%device_runtime_ff - REAL(device_runtime_fill_fraction, c_double)) > &
          1.0E-10_c_double) RETURN
      IF (TRIM(cache%xc_fun_name) /= TRIM(xc_fun_name)) RETURN
      IF (TRIM(cache%grid_type) /= TRIM(grid_type)) RETURN
      IF (TRIM(cache%radial_quad) /= TRIM(radial_quadrature)) RETURN
      IF (TRIM(cache%pruning_scheme) /= TRIM(pruning_scheme)) RETURN
      IF (TRIM(cache%lb_exec_space) /= TRIM(lb_exec_space)) RETURN
      IF (TRIM(cache%int_exec_space) /= TRIM(int_exec_space)) RETURN
      IF (TRIM(cache%model_name) /= TRIM(model_name)) RETURN
      IF (cache%gradient_self_runtime .NEQV. use_gradient_self_runtime) RETURN
      IF (.NOT. ALLOCATED(cache%last_positions)) RETURN
      IF (SIZE(cache%last_positions, 2) /= natom) RETURN
      DO iatom = 1, natom
         IF (ANY(ABS(cache%last_positions(:, iatom) - &
                     REAL(particle_set(iatom)%r(:), c_double)) > pos_tol)) RETURN
      END DO
      cache_valid = .TRUE.
#else
      MARK_USED(cache)
      MARK_USED(particle_set)
      MARK_USED(natom)
      MARK_USED(nspins)
      MARK_USED(batch_size)
      MARK_USED(device_runtime_fill_fraction)
      MARK_USED(xc_fun_name)
      MARK_USED(grid_type)
      MARK_USED(radial_quadrature)
      MARK_USED(pruning_scheme)
      MARK_USED(lb_exec_space)
      MARK_USED(int_exec_space)
      MARK_USED(model_name)
      MARK_USED(use_gradient_self_runtime)
      cache_valid = .FALSE.
#endif
   END FUNCTION gauxc_cache_is_valid

END MODULE xc_gauxc_cache
