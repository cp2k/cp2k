!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Cell-gauge scalar/spinor reconstruction on independent complete property meshes.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE qs_property_symmetry
   USE cell_types,                      ONLY: cell_type
   USE input_cp2k_kpoints,              ONLY: use_k290_kpoint_symmetry,&
                                              use_spglib_kpoint_backend,&
                                              use_spglib_kpoint_symmetry
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_initialize
   USE kpoint_mo_symmetry_methods,      ONLY: kpoint_cell_gauge_action,&
                                              kpoint_same_periodic
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_sym_type,&
                                              kpoint_type
   USE particle_types,                  ONLY: particle_type
   USE qs_basis_rotation_methods,       ONLY: qs_basis_rotation
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE topology_symmetry,               ONLY: spin_rotation,&
                                              transform_spinor_frame
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_property_symmetry'
   PUBLIC :: property_symmetry_plan, property_symmetry_frame
CONTAINS
! **************************************************************************************************
!> \brief Construct a geometric orbit plan; each reconstructed frame still needs H/S validation.
!> \param qs_env Frozen property environment
!> \param k Full mesh in fractional reciprocal coordinates
!> \param grid Complete Gamma-centered property mesh dimensions
!> \param projection Cartesian periodic-space projector
!> \param backend Existing K290 or SPGLIB backend identifier
!> \param tolerance Geometric symmetry tolerance
!> \param points Owned irreducible k-point object; caller must release it
!> \param representative Source irreducible index for each full point
!> \param operation Source-to-target operation, zero identity or minus one time reversal
! **************************************************************************************************
   SUBROUTINE property_symmetry_plan(qs_env, k, grid, projection, backend, tolerance, points, representative, operation)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), INTENT(IN)                               :: k(:, :)
      INTEGER, INTENT(IN)                                :: grid(3)
      REAL(dp), INTENT(IN)                               :: projection(3, 3)
      INTEGER, INTENT(IN)                                :: backend
      REAL(dp), INTENT(IN)                               :: tolerance
      TYPE(kpoint_type), POINTER                         :: points
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: representative(:), operation(:)

      INTEGER                                            :: ik, ir, is, nk
      REAL(dp)                                           :: r(3, 3)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(kpoint_sym_type), POINTER                     :: sym
      TYPE(particle_type), POINTER                       :: particles(:)

      CALL get_qs_env(qs_env, cell=cell, particle_set=particles)
      nk = SIZE(k, 2)
      CALL kpoint_create(points)
      points%kp_scheme = 'MONKHORST-PACK'
      points%nkp_grid = grid
      points%gamma_centered = .TRUE.
      points%symmetry = .TRUE.
      points%symmetry_backend = backend
      points%symmetry_reduction_method = use_k290_kpoint_symmetry
      IF (backend == use_spglib_kpoint_backend) points%symmetry_reduction_method = use_spglib_kpoint_symmetry
      points%eps_geo = tolerance
      CALL kpoint_initialize(points, particles, cell)
      CALL qs_basis_rotation(qs_env, points)
      ALLOCATE (representative(nk), operation(nk))
      representative = 0
      operation = 0
      DO ik = 1, nk
         DO ir = 1, points%nkp
            IF (kpoint_same_periodic(k(:, ik), points%xkp(:, ir))) THEN
               representative(ik) = ir
               EXIT
            END IF
            IF (kpoint_same_periodic(k(:, ik), -points%xkp(:, ir))) THEN
               representative(ik) = ir
               operation(ik) = -1
               EXIT
            END IF
            sym => points%kp_sym(ir)%kpoint_sym
            IF (.NOT. sym%apply_symmetry) CYCLE
            DO is = 1, sym%nwred
               IF (.NOT. kpoint_same_periodic(k(:, ik), sym%xkp(:, is))) CYCLE
               r = sym%rot(:, :, is)
               IF (MAXVAL(ABS(MATMUL(r, MATMUL(projection, TRANSPOSE(r))) - projection)) > 1.e-9_dp) CYCLE
               representative(ik) = ir
               operation(ik) = is
               EXIT
            END DO
            IF (representative(ik) > 0) EXIT
         END DO
      END DO
      IF (ANY(representative == 0)) THEN
         CPABORT('Property symmetry orbit does not preserve the periodic subspace.')
      END IF
   END SUBROUTINE property_symmetry_plan
! **************************************************************************************************
!> \brief Reconstruct a target eigenframe with atomic image phases and physical SU(2)/Theta action.
!> \param points Property symmetry plan
!> \param representative Source index
!> \param operation Signed operation index
!> \param k Target point in the cell Fourier gauge
!> \param source Source eigenframe
!> \param soc Whether the frame is spin-major doubled
!> \param TARGET Reconstructed frame, not rediagonalized or gauge-aligned
!> \param rotation Cartesian spatial rotation for property tensors
!> \param antiunitary Whether physical time reversal was applied
!> \param status Zero on success; the caller must additionally check the target H/S eigensystem
! **************************************************************************************************
   SUBROUTINE property_symmetry_frame(points, representative, operation, k, source, soc, TARGET, &
                                      rotation, antiunitary, status)
      TYPE(kpoint_type), POINTER                         :: points
      INTEGER, INTENT(IN)                                :: representative, operation
      REAL(dp), INTENT(IN)                               :: k(3)
      COMPLEX(dp), INTENT(IN)                            :: source(:, :)
      LOGICAL, INTENT(IN)                                :: soc
      COMPLEX(dp), INTENT(OUT)                           :: target(:, :)
      REAL(dp), INTENT(OUT)                              :: rotation(3, 3)
      LOGICAL, INTENT(OUT)                               :: antiunitary
      INTEGER, INTENT(OUT)                               :: status

      COMPLEX(dp)                                        :: spin(2, 2)
      COMPLEX(dp), ALLOCATABLE                           :: orbital(:, :)
      INTEGER                                            :: n

      n = SIZE(source, 1)
      IF (soc) n = n/2
      ALLOCATE (orbital(n, n))
      CALL kpoint_cell_gauge_action(points, representative, operation, k, orbital, rotation, antiunitary, status)
      IF (status /= 0) RETURN
      IF (soc) THEN
         CALL spin_rotation(rotation, spin, 1.e-9_dp, status)
         IF (status /= 0) RETURN
         ! CP2K assembles SOC with (sigma_x,-sigma_y,sigma_z), the conjugate spin representation.
         spin = CONJG(spin)
         CALL transform_spinor_frame(orbital, spin, source, antiunitary, TARGET, 1.e-9_dp, status)
      ELSE IF (antiunitary) THEN
         TARGET = MATMUL(orbital, CONJG(source))
      ELSE
         TARGET = MATMUL(orbital, source)
      END IF
   END SUBROUTINE property_symmetry_frame
END MODULE qs_property_symmetry
