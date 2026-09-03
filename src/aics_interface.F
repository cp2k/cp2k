!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   C-interoperable bridge between CP2K-owned AICS fields and the embedded solver.
!> \author  Ziwei Chai
! **************************************************************************************************
MODULE aics_interface
   USE iso_c_binding,                   ONLY: c_char,&
                                              c_double,&
                                              c_int,&
                                              c_null_char,&
                                              c_null_ptr,&
                                              c_ptr,&
                                              c_size_t

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: clear_aics_field_registrations, cp2k_aics_finalize_from_fortran, &
             cp2k_aics_get_electrostatic_potential, &
             cp2k_aics_get_fd_field_squared_sensitivity_parallel, &
             cp2k_aics_get_fd_field_squared_sensitivity_perpendicular, &
             cp2k_aics_get_grid_metadata, cp2k_aics_get_relative_permittivity_parallel, &
             cp2k_aics_get_relative_permittivity_perpendicular, &
             cp2k_aics_get_solute_charge_density, cp2k_aics_get_last_error_from_fortran, &
             cp2k_aics_cancel_initialize_from_fortran, cp2k_aics_initialize_from_fortran, &
             cp2k_aics_initialize_prepared_from_fortran, &
             cp2k_aics_scope_enter_from_fortran, cp2k_aics_scope_leave_from_fortran, &
             cp2k_python_finalize_from_fortran, &
             cp2k_aics_set_fd_field_squared_sensitivity_parallel, &
             cp2k_aics_set_fd_field_squared_sensitivity_perpendicular, &
             cp2k_aics_set_electrostatic_potential, cp2k_aics_set_grid_metadata, &
             cp2k_aics_set_relative_permittivity_parallel, &
             cp2k_aics_set_relative_permittivity_perpendicular, &
             cp2k_aics_set_solute_charge_density, cp2k_aics_solve_from_fortran

   INTERFACE
      FUNCTION cp2k_aics_initialize(comm, normal_output_enabled, poisson_solver, normal_axis, relative_tolerance, &
                                    absolute_tolerance, max_iterations) RESULT(status) &
         BIND(C, name="cp2k_aics_initialize")
         IMPORT :: c_double, c_int
         INTEGER(c_int), INTENT(IN), VALUE :: comm
         INTEGER(c_int), INTENT(IN), VALUE :: normal_output_enabled
         INTEGER(c_int), INTENT(IN), VALUE :: poisson_solver
         INTEGER(c_int), INTENT(IN), VALUE :: normal_axis
         REAL(c_double), INTENT(IN), VALUE :: relative_tolerance, absolute_tolerance
         INTEGER(c_int), INTENT(IN), VALUE :: max_iterations
         INTEGER(c_int) :: status
      END FUNCTION cp2k_aics_initialize

      FUNCTION cp2k_aics_initialize_prepared() RESULT(status) &
         BIND(C, name="cp2k_aics_initialize_prepared")
         IMPORT :: c_int
         INTEGER(c_int) :: status
      END FUNCTION cp2k_aics_initialize_prepared

      SUBROUTINE cp2k_aics_cancel_initialize() &
         BIND(C, name="cp2k_aics_cancel_initialize")
      END SUBROUTINE cp2k_aics_cancel_initialize

      FUNCTION cp2k_aics_solve() RESULT(status) &
         BIND(C, name="cp2k_aics_solve")
         IMPORT :: c_int
         INTEGER(c_int) :: status
      END FUNCTION cp2k_aics_solve

      FUNCTION cp2k_aics_finalize() RESULT(status) &
         BIND(C, name="cp2k_aics_finalize")
         IMPORT :: c_int
         INTEGER(c_int) :: status
      END FUNCTION cp2k_aics_finalize

      FUNCTION cp2k_aics_scope_enter() RESULT(status) &
         BIND(C, name="cp2k_aics_scope_enter")
         IMPORT :: c_int
         INTEGER(c_int) :: status
      END FUNCTION cp2k_aics_scope_enter

      FUNCTION cp2k_aics_scope_leave() RESULT(status) &
         BIND(C, name="cp2k_aics_scope_leave")
         IMPORT :: c_int
         INTEGER(c_int) :: status
      END FUNCTION cp2k_aics_scope_leave

      FUNCTION cp2k_python_finalize() RESULT(status) &
         BIND(C, name="cp2k_python_finalize")
         IMPORT :: c_int
         INTEGER(c_int) :: status
      END FUNCTION cp2k_python_finalize

      SUBROUTINE cp2k_aics_get_last_error(buffer, buffer_size) &
         BIND(C, name="cp2k_aics_get_last_error")
         IMPORT :: c_char, c_size_t
         CHARACTER(kind=c_char), DIMENSION(*), INTENT(out) :: buffer
         INTEGER(c_size_t), INTENT(IN), VALUE :: buffer_size
      END SUBROUTINE cp2k_aics_get_last_error
   END INTERFACE

   ! The local metadata contract is half-open: end = start + extent on every axis.
   ! Every registered array is C_NULL_PTR exactly when any local extent is zero;
   ! otherwise it addresses the first element of the contiguous Fortran-order block.
   TYPE(c_ptr), SAVE :: solute_charge_density_ptr = c_null_ptr, &
                        relative_permittivity_parallel_ptr = c_null_ptr, &
                        relative_permittivity_perpendicular_ptr = c_null_ptr, &
                        dielectric_field_squared_sensitivity_parallel_ptr = c_null_ptr, &
                        dielectric_field_squared_sensitivity_perpendicular_ptr = c_null_ptr, &
                        electrostatic_potential_ptr = c_null_ptr
   INTEGER(c_size_t), SAVE :: global_npts1 = 0, global_npts2 = 0, global_npts3 = 0, &
                              local_npts1 = 0, local_npts2 = 0, local_npts3 = 0, local_l1 = 0, &
                              local_l2 = 0, local_l3 = 0, local_u1 = 0, local_u2 = 0, local_u3 = 0
   REAL(c_double), SAVE :: grid_dr1 = 0.0_c_double, grid_dr2 = 0.0_c_double, grid_dr3 = 0.0_c_double

CONTAINS

! **************************************************************************************************
!> \brief   Clear every non-owning CP2K field registration and its local grid metadata.
!> \author  Ziwei Chai
!> \note    This routine is process-local, non-collective, and non-fallible. A zero-local
!>          registration during a solve remains distinguishable because it has nonzero global grid
!>          metadata and a legal empty local slice; the cleared state is all-zero metadata with
!>          C_NULL_PTR fields.
! **************************************************************************************************
   SUBROUTINE clear_aics_field_registrations()
      solute_charge_density_ptr = c_null_ptr
      relative_permittivity_parallel_ptr = c_null_ptr
      relative_permittivity_perpendicular_ptr = c_null_ptr
      dielectric_field_squared_sensitivity_parallel_ptr = c_null_ptr
      dielectric_field_squared_sensitivity_perpendicular_ptr = c_null_ptr
      electrostatic_potential_ptr = c_null_ptr
      global_npts1 = 0_c_size_t
      global_npts2 = 0_c_size_t
      global_npts3 = 0_c_size_t
      local_npts1 = 0_c_size_t
      local_npts2 = 0_c_size_t
      local_npts3 = 0_c_size_t
      local_l1 = 0_c_size_t
      local_l2 = 0_c_size_t
      local_l3 = 0_c_size_t
      local_u1 = 0_c_size_t
      local_u2 = 0_c_size_t
      local_u3 = 0_c_size_t
      grid_dr1 = 0.0_c_double
      grid_dr2 = 0.0_c_double
      grid_dr3 = 0.0_c_double
   END SUBROUTINE clear_aics_field_registrations

! **************************************************************************************************
!> \brief   Copy the process-local bridge failure diagnostic into a Fortran string.
!> \param   message   Bounded, synchronous copy of the most recent bridge error.
!> \author  Ziwei Chai
! **************************************************************************************************
   SUBROUTINE cp2k_aics_get_last_error_from_fortran(message)
      CHARACTER(len=*), INTENT(out)                      :: message

      CHARACTER(kind=c_char), DIMENSION(1024)            :: buffer
      INTEGER                                            :: i

      buffer(:) = c_null_char
      CALL cp2k_aics_get_last_error(buffer, INT(SIZE(buffer), kind=c_size_t))
      message = ""
      DO i = 1, MIN(SIZE(buffer), LEN(message))
         IF (buffer(i) == c_null_char) EXIT
         message(i:i) = buffer(i)
      END DO
   END SUBROUTINE cp2k_aics_get_last_error_from_fortran

! **************************************************************************************************
!> \brief   Borrow the local total solute charge-density array for a synchronous AICS solve.
!> \param   ptr   Pointer to a contiguous local real(c_double) array in Fortran order.
!> \author  Ziwei Chai
!> \note    The view is owned by CP2K and is valid only until the synchronous solve returns.
! **************************************************************************************************
   SUBROUTINE cp2k_aics_set_solute_charge_density(ptr) BIND(C, name="cp2k_aics_set_solute_charge_density")
      TYPE(c_ptr), INTENT(IN), VALUE                     :: ptr

      solute_charge_density_ptr = ptr
   END SUBROUTINE cp2k_aics_set_solute_charge_density

! **************************************************************************************************
!> \brief   Return the borrowed local solute charge-density pointer.
!> \return  Borrowed Fortran-order view, valid only during the synchronous AICS solve.
!> \author  Ziwei Chai
! **************************************************************************************************
   FUNCTION cp2k_aics_get_solute_charge_density() RESULT(ptr) BIND(C, name="cp2k_aics_get_solute_charge_density")
      TYPE(c_ptr)                                        :: ptr

      ptr = solute_charge_density_ptr
   END FUNCTION cp2k_aics_get_solute_charge_density

! **************************************************************************************************
!> \brief   Borrow the local relative-permittivity array for both in-plane Cartesian axes.
!> \param   ptr   Pointer to a contiguous local real(c_double) array in Fortran order.
!> \author  Ziwei Chai
!> \note    The non-owning view is valid only for the synchronous AICS solve.
! **************************************************************************************************
   SUBROUTINE cp2k_aics_set_relative_permittivity_parallel(ptr) BIND(C, name="cp2k_aics_set_relative_permittivity_parallel")
      TYPE(c_ptr), INTENT(IN), VALUE                     :: ptr

      relative_permittivity_parallel_ptr = ptr
   END SUBROUTINE cp2k_aics_set_relative_permittivity_parallel

! **************************************************************************************************
!> \brief   Return the borrowed in-plane relative-permittivity pointer.
!> \return  Borrowed Fortran-order view, valid only during the synchronous AICS solve.
!> \author  Ziwei Chai
! **************************************************************************************************
   FUNCTION cp2k_aics_get_relative_permittivity_parallel() RESULT(ptr) BIND(C, name="cp2k_aics_get_relative_permittivity_parallel")
      TYPE(c_ptr)                                        :: ptr

      ptr = relative_permittivity_parallel_ptr
   END FUNCTION cp2k_aics_get_relative_permittivity_parallel

! **************************************************************************************************
!> \brief   Borrow the local relative-permittivity array for the oriented normal axis.
!> \param   ptr   Pointer to a contiguous local real(c_double) array in Fortran order.
!> \author  Ziwei Chai
!> \note    The non-owning view is valid only for the synchronous AICS solve.
! **************************************************************************************************
   SUBROUTINE cp2k_aics_set_relative_permittivity_perpendicular( &
      ptr) BIND(C, name="cp2k_aics_set_relative_permittivity_perpendicular")
      TYPE(c_ptr), INTENT(IN), VALUE                     :: ptr

      relative_permittivity_perpendicular_ptr = ptr
   END SUBROUTINE cp2k_aics_set_relative_permittivity_perpendicular

! **************************************************************************************************
!> \brief   Return the borrowed oriented-normal relative-permittivity pointer.
!> \return  Borrowed Fortran-order view, valid only during the synchronous AICS solve.
!> \author  Ziwei Chai
! **************************************************************************************************
   FUNCTION cp2k_aics_get_relative_permittivity_perpendicular( &
      ) RESULT(ptr) BIND(C, name="cp2k_aics_get_relative_permittivity_perpendicular")
      TYPE(c_ptr)                                        :: ptr

      ptr = relative_permittivity_perpendicular_ptr
   END FUNCTION cp2k_aics_get_relative_permittivity_perpendicular

! **************************************************************************************************
!> \brief   Borrow the local electrostatic-potential output array.
!> \param   ptr   Pointer to a contiguous local real(c_double) array in Fortran order.
!> \author  Ziwei Chai
!> \note    The non-owning view is valid only for the synchronous AICS solve.
! **************************************************************************************************
   SUBROUTINE cp2k_aics_set_electrostatic_potential(ptr) BIND(C, name="cp2k_aics_set_electrostatic_potential")
      TYPE(c_ptr), INTENT(IN), VALUE                     :: ptr

      electrostatic_potential_ptr = ptr
   END SUBROUTINE cp2k_aics_set_electrostatic_potential

! **************************************************************************************************
!> \brief   Return the borrowed local electrostatic-potential pointer.
!> \return  Borrowed Fortran-order view, valid only during the synchronous AICS solve.
!> \author  Ziwei Chai
! **************************************************************************************************
   FUNCTION cp2k_aics_get_electrostatic_potential() RESULT(ptr) BIND(C, name="cp2k_aics_get_electrostatic_potential")
      TYPE(c_ptr)                                        :: ptr

      ptr = electrostatic_potential_ptr
   END FUNCTION cp2k_aics_get_electrostatic_potential

! **************************************************************************************************
!> \brief   Borrow the local parallel FD dielectric field-squared sensitivity output array.
!> \param   ptr   Pointer to a contiguous local real(c_double) array in Fortran order.
!> \author  Ziwei Chai
!> \note    CP2K owns the array; Python may write it only during the synchronous FD solve.
! **************************************************************************************************
   SUBROUTINE cp2k_aics_set_fd_field_squared_sensitivity_parallel( &
      ptr) BIND(C, name="cp2k_aics_set_fd_field_squared_sensitivity_parallel")
      TYPE(c_ptr), INTENT(IN), VALUE                     :: ptr

      dielectric_field_squared_sensitivity_parallel_ptr = ptr
   END SUBROUTINE cp2k_aics_set_fd_field_squared_sensitivity_parallel

! **************************************************************************************************
!> \brief   Return the borrowed parallel FD dielectric field-squared sensitivity pointer.
!> \return  CP2K-owned Fortran-order output view, valid only during the synchronous FD solve.
!> \author  Ziwei Chai
! **************************************************************************************************
   FUNCTION cp2k_aics_get_fd_field_squared_sensitivity_parallel( &
      ) RESULT(ptr) BIND(C, name="cp2k_aics_get_fd_field_squared_sensitivity_parallel")
      TYPE(c_ptr)                                        :: ptr

      ptr = dielectric_field_squared_sensitivity_parallel_ptr
   END FUNCTION cp2k_aics_get_fd_field_squared_sensitivity_parallel

! **************************************************************************************************
!> \brief   Borrow the local perpendicular FD dielectric field-squared sensitivity output array.
!> \param   ptr   Pointer to a contiguous local real(c_double) array in Fortran order.
!> \author  Ziwei Chai
!> \note    CP2K owns the array; Python may write it only during the synchronous FD solve.
! **************************************************************************************************
   SUBROUTINE cp2k_aics_set_fd_field_squared_sensitivity_perpendicular( &
      ptr) BIND(C, name="cp2k_aics_set_fd_field_squared_sensitivity_perpendicular")
      TYPE(c_ptr), INTENT(IN), VALUE                     :: ptr

      dielectric_field_squared_sensitivity_perpendicular_ptr = ptr
   END SUBROUTINE cp2k_aics_set_fd_field_squared_sensitivity_perpendicular

! **************************************************************************************************
!> \brief   Return the borrowed perpendicular FD dielectric field-squared sensitivity pointer.
!> \return  CP2K-owned Fortran-order output view, valid only during the synchronous FD solve.
!> \author  Ziwei Chai
! **************************************************************************************************
   FUNCTION cp2k_aics_get_fd_field_squared_sensitivity_perpendicular( &
      ) RESULT(ptr) BIND(C, name="cp2k_aics_get_fd_field_squared_sensitivity_perpendicular")
      TYPE(c_ptr)                                        :: ptr

      ptr = dielectric_field_squared_sensitivity_perpendicular_ptr
   END FUNCTION cp2k_aics_get_fd_field_squared_sensitivity_perpendicular

! **************************************************************************************************
!> \brief   Store global/local grid sizes, local half-open bounds, and grid spacings.
!> \param   npts1    Global Cartesian grid dimension along x.
!> \param   npts2    Global Cartesian grid dimension along y.
!> \param   npts3    Global Cartesian grid dimension along z.
!> \param   lnpts1   Local Cartesian grid dimension along x.
!> \param   lnpts2   Local Cartesian grid dimension along y.
!> \param   lnpts3   Local Cartesian grid dimension along z.
!> \param   l1       Zero-based local x start.
!> \param   l2       Zero-based local y start.
!> \param   l3       Zero-based local z start.
!> \param   u1       Zero-based exclusive local x end, equal to l1 + lnpts1.
!> \param   u2       Zero-based exclusive local y end, equal to l2 + lnpts2.
!> \param   u3       Zero-based exclusive local z end, equal to l3 + lnpts3.
!> \param   dr1      Cartesian x spacing in CP2K internal units.
!> \param   dr2      Cartesian y spacing in CP2K internal units.
!> \param   dr3      Cartesian z spacing in CP2K internal units.
!> \author  Ziwei Chai
! **************************************************************************************************
   SUBROUTINE cp2k_aics_set_grid_metadata( &
      npts1, npts2, npts3, lnpts1, lnpts2, lnpts3, l1, l2, l3, &
      u1, u2, u3, dr1, dr2, dr3) BIND(C, name="cp2k_aics_set_grid_metadata")
      INTEGER(c_size_t), INTENT(IN), VALUE               :: npts1, npts2, npts3, lnpts1, lnpts2, &
                                                            lnpts3, l1, l2, l3, u1, u2, u3
      REAL(c_double), INTENT(IN), VALUE                  :: dr1, dr2, dr3

      global_npts1 = npts1
      global_npts2 = npts2
      global_npts3 = npts3
      local_npts1 = lnpts1
      local_npts2 = lnpts2
      local_npts3 = lnpts3
      grid_dr1 = dr1
      grid_dr2 = dr2
      grid_dr3 = dr3
      local_l1 = l1
      local_l2 = l2
      local_l3 = l3
      local_u1 = u1
      local_u2 = u2
      local_u3 = u3
   END SUBROUTINE cp2k_aics_set_grid_metadata

! **************************************************************************************************
!> \brief   Return global/local grid sizes, local half-open bounds, and grid spacings.
!> \param   npts1    Global Cartesian grid dimension along x.
!> \param   npts2    Global Cartesian grid dimension along y.
!> \param   npts3    Global Cartesian grid dimension along z.
!> \param   lnpts1   Local Cartesian grid dimension along x.
!> \param   lnpts2   Local Cartesian grid dimension along y.
!> \param   lnpts3   Local Cartesian grid dimension along z.
!> \param   l1       Zero-based local x start.
!> \param   l2       Zero-based local y start.
!> \param   l3       Zero-based local z start.
!> \param   u1       Zero-based exclusive local x end, equal to l1 + lnpts1.
!> \param   u2       Zero-based exclusive local y end, equal to l2 + lnpts2.
!> \param   u3       Zero-based exclusive local z end, equal to l3 + lnpts3.
!> \param   dr1      Cartesian x spacing in CP2K internal units.
!> \param   dr2      Cartesian y spacing in CP2K internal units.
!> \param   dr3      Cartesian z spacing in CP2K internal units.
!> \author  Ziwei Chai
! **************************************************************************************************
   SUBROUTINE cp2k_aics_get_grid_metadata( &
      npts1, npts2, npts3, lnpts1, lnpts2, lnpts3, l1, l2, l3, &
      u1, u2, u3, dr1, dr2, dr3) BIND(C, name="cp2k_aics_get_grid_metadata")
      INTEGER(c_size_t), INTENT(out)                     :: npts1, npts2, npts3, lnpts1, lnpts2, &
                                                            lnpts3, l1, l2, l3, u1, u2, u3
      REAL(c_double), INTENT(out)                        :: dr1, dr2, dr3

      npts1 = global_npts1
      npts2 = global_npts2
      npts3 = global_npts3
      lnpts1 = local_npts1
      lnpts2 = local_npts2
      lnpts3 = local_npts3
      dr1 = grid_dr1
      dr2 = grid_dr2
      dr3 = grid_dr3
      l1 = local_l1
      l2 = local_l2
      l3 = local_l3
      u1 = local_u1
      u2 = local_u2
      u3 = local_u3
   END SUBROUTINE cp2k_aics_get_grid_metadata

! **************************************************************************************************
!> \brief   Prepare AICS rank-local runtime and import state before communicator duplication.
!> \param   comm                    Active AICS Fortran MPI communicator handle.
!> \param   normal_output_enabled   Integer flag enabling normal solver output.
!> \param   poisson_solver          Selected AICS Poisson solver.
!> \param   normal_axis             Cartesian oriented-normal axis index.
!> \param   relative_tolerance      Relative linear-solver tolerance.
!> \param   absolute_tolerance      Absolute linear-solver tolerance.
!> \param   max_iterations          Maximum linear-solver iteration count.
!> \param   status                  Status returned by embedded AICS context initialization.
!> \author  Ziwei Chai
! **************************************************************************************************
   SUBROUTINE cp2k_aics_initialize_from_fortran(comm, normal_output_enabled, poisson_solver, normal_axis, &
                                                relative_tolerance, &
                                                absolute_tolerance, max_iterations, status)
      INTEGER, INTENT(in)                                :: comm, normal_output_enabled, &
                                                            poisson_solver, normal_axis
      REAL(c_double), INTENT(in)                         :: relative_tolerance, absolute_tolerance
      INTEGER, INTENT(in)                                :: max_iterations
      INTEGER(c_int), INTENT(out)                        :: status

      status = cp2k_aics_initialize(INT(comm, KIND=c_int), &
                                    INT(normal_output_enabled, KIND=c_int), &
                                    INT(poisson_solver, KIND=c_int), INT(normal_axis, KIND=c_int), relative_tolerance, &
                                    absolute_tolerance, INT(max_iterations, KIND=c_int))
   END SUBROUTINE cp2k_aics_initialize_from_fortran

! **************************************************************************************************
!> \brief   Collectively initialize the prepared AICS context.
!> \param   status   Status returned by embedded AICS context initialization.
!> \author  Ziwei Chai
! **************************************************************************************************
   SUBROUTINE cp2k_aics_initialize_prepared_from_fortran(status)
      INTEGER(c_int), INTENT(out)                        :: status

      status = cp2k_aics_initialize_prepared()
   END SUBROUTINE cp2k_aics_initialize_prepared_from_fortran

! **************************************************************************************************
!> \brief   Discard a prepared rank-local AICS request without creating a context.
!> \author  Ziwei Chai
! **************************************************************************************************
   SUBROUTINE cp2k_aics_cancel_initialize_from_fortran()
      CALL cp2k_aics_cancel_initialize()
   END SUBROUTINE cp2k_aics_cancel_initialize_from_fortran

! **************************************************************************************************
!> \brief   Run the synchronous collective solve on the active AICS context.
!> \param   status   Status returned by the embedded AICS solve.
!> \author  Ziwei Chai
!> \note    All registered array pointers must remain valid until this routine returns.
! **************************************************************************************************
   SUBROUTINE cp2k_aics_solve_from_fortran(status)
      INTEGER(c_int), INTENT(out)                        :: status

      status = cp2k_aics_solve()
   END SUBROUTINE cp2k_aics_solve_from_fortran

! **************************************************************************************************
!> \brief   Collectively finalize the active AICS context.
!> \param   status   Status returned by AICS context finalization.
!> \author  Ziwei Chai
! **************************************************************************************************
   SUBROUTINE cp2k_aics_finalize_from_fortran(status)
      INTEGER(c_int), INTENT(out)                        :: status

      CALL clear_aics_field_registrations()
      status = cp2k_aics_finalize()
      CALL clear_aics_field_registrations()
   END SUBROUTINE cp2k_aics_finalize_from_fortran

! **************************************************************************************************
!> \brief   Collectively enter one CP2K AICS execution scope.
!> \param   status   Status returned by the AICS context-scope transition.
!> \author  Ziwei Chai
! **************************************************************************************************
   SUBROUTINE cp2k_aics_scope_enter_from_fortran(status)
      INTEGER(c_int), INTENT(out)                        :: status

      status = cp2k_aics_scope_enter()
   END SUBROUTINE cp2k_aics_scope_enter_from_fortran

! **************************************************************************************************
!> \brief   Collectively leave the current CP2K AICS execution scope.
!> \param   status   Status returned by the AICS context-scope transition.
!> \author  Ziwei Chai
! **************************************************************************************************
   SUBROUTINE cp2k_aics_scope_leave_from_fortran(status)
      INTEGER(c_int), INTENT(out)                        :: status

      status = cp2k_aics_scope_leave()
   END SUBROUTINE cp2k_aics_scope_leave_from_fortran

! **************************************************************************************************
!> \brief   Finalize the CP2K-owned Python interpreter at process-level CP2K teardown.
!> \param   status   Status returned by the embedded Python runtime finalization.
!> \author  Ziwei Chai
! **************************************************************************************************
   SUBROUTINE cp2k_python_finalize_from_fortran(status)
      INTEGER(c_int), INTENT(out)                        :: status

      status = cp2k_python_finalize()
   END SUBROUTINE cp2k_python_finalize_from_fortran

END MODULE aics_interface
